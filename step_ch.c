
#include "decs.h"
#include "metric.h"
#include "patchwork.h"


/******************************************************************************************
  step_ch(): 
 -----------
     -- top-level time-stepping driver routine; 
     -- responsible for handling the order in which the timestep, boundary conditions 
        and fixups are done;
******************************************************************************************/

void step_ch(void)
{
	int i, j, k, l, failure_exists_world;

	TRACE_BEG;
  
        
	/***  PREDICTOR STEP *******/
	advance( p, p, 0.5*dx[0], ph );


	if( I.i_belong_global_patch ) {/*CHANGE*/
		// Do it only for global patch, because bounds()
		// updates local patch's edge-zones, too, but that
		// needs to be followed by connect_patches() to
		// wash away the update of normal edge-zones that
		// are inside the global patch.
		// connect_patches should not be performed at
		// half-step.
		bounds( ph, 0 );
	} else {
		// Just do internal ghost zones update for local patch.
		// Here is ignoring updating local edge-zones outside
		// global patch, but I don't know how to update them
		// without updating those standard local ghost zones... // LATER
		bounds_mpi( ph );
		bounds_pflag_mpi( pflag );
	}


	/* full-step */
	advance( p, ph, dx[0], p );

	if( sync_patches ) {
		connect_patches( p ); // patchwork_gl.c

		MPI_Allreduce( &failure_exists,
			       &failure_exists_world,
			       1,
			       MPI_INT,
			       MPI_MAX,
			       MPI_COMM_WORLD );

		if( failure_exists_world ) {
			connect_patches( p ); // patchwork_gl.c
		}
	}


	return;
}





/******************************************************************************************
  advance(): 
 -----------
     -- handles details pertaining to how a time step is performed, e.g. reconstruction, 
         flux calculation, and EOS evaluation; 
     -- timestep here means integrating the primitives from t to (t+Dt);
     -- arguments: 
           p_i = prim. var's at t
           p_h = prim. var's to be used for flux calculation 
           p_f = prim. var's at t+Dt  (on exit)
******************************************************************************************/

void advance(
	double ****p_i,
	double ****p_h,
	double Dt,
	double ****p_f
	)
{
	int i, j, k, l, d;
	long cnt;
	long int ind_U;
	double dU[NP], U[NP];

	struct of_state  q;
	struct of_geom  *geom_i, *geom_h, *geom_f;
	struct of_coord *coords;

	TRACE_BEG;

	/****************************************************************
          Find the fluxes in all directions : 
	*****************************************************************/
	numerical_flux( p_L, p_R, F );

	LOOP {

		cnt = get_cnt( i, j, k, N1TOT, N2TOT, N3TOT );
		if( cover_flags[ cnt ] == 0 ) {

			get_geometry( i, j, k, CENT, n_beg, geom_i ); 

			/* calculate conserved quantities */
			if( n_substep == 0 ) {	// 1st time
				get_geometry( i, j, k, CENT, n_mid, geom_f );
				primtoflux( p_i[i][j][k], &q, 0, geom_i, U );

			}  else {		// 2nd time
				get_geometry( i, j, k, CENT, n_mid, geom_h );
			}

			/* evolve conserved quantities */
			PLOOP {
				U[l] += Dt*(
					- ( F[i+1][j  ][k  ][0][l] - F[i][j][k][0][l] ) * invdx[1] 
					- ( F[i  ][j+1][k  ][1][l] - F[i][j][k][1][l] ) * invdx[2] 
					- ( F[i  ][j  ][k+1][2][l] - F[i][j][k][2][l] ) * invdx[3] 
					+ dU[l] 
					);
			}

			/**************************************************************************
			  Find primitive var's from conserved var's, taking into special 
			  conditions, most fixups (e.g. floor and GAMMAMAX), inversion checks...
			**************************************************************************/
			recover_primitives( i, j, k, p_f[i][j][k], U, geom_f, geom_i, coords );
		}
	}


	TRACE_END;

	return;
}





/******************************************************************************************
  numerical_flux(): 
 ------------------
     -- Using the reconstructed prim. variables, calculates the numerical flux used 
        when integrating the equations of motion;
     -- Responsible for calculating wave speeds and such needed for the 
        flux formula that is used here;
     -- Currently uses HLLE or LF-like flux formula 
     -- Assumes that dir = 0,1,2 = x1,x2,x3
******************************************************************************************/

static void numerical_flux( 
		    double *****p_L,
		    double *****p_R,
		    double *****F
		    )
{
	unsigned int i, j, k, l, i_s, j_s, k_s, i_e, j_e, k_e, d, dim, pos, ind;
	long cnt, cnt_1;
	double Ul[NP], Ur[NP], Fl[NP], Fr[NP], btmp ;
	double *pl, *pr;
	double cminl, cminr, cmaxl, cmaxr, cmax, cmin, ctop;
	struct of_state ql, qr;
	struct of_geom *geom;

	TRACE_BEG;


	/* Reset max speeds for timestep calculation: */
	DLOOP1 { max_char_speed[i] = 0.; }
	dt_char_min = 1.e200;


	FACE_LOOP {

		i_s = N1S + idir[d] - 1;
		j_s = N2S + jdir[d] - 1;
		k_s = N3S + kdir[d] - 1;
		i_e = N1E+1;
		j_e = N2E+1;
		k_e = N3E+1;


		dim = d + 1;  // Position in a 4-vector to which this direction corresponds
		pos = d;      // Position of the flux in the cell, assumes FACE1-3 = 0-2

    		ind = 0;
		for( i=i_s; i<=i_e; i++ ) {
		for( j=j_s; j<=j_e; j++ ) {
		for( k=k_s; k<=k_e; k++ ) {

			cnt   = get_cnt( i        , j        , k        , N1TOT, N2TOT, N3TOT ); // misc.c
			cnt_1 = get_cnt( i-idir[d], j-jdir[d], k-kdir[d], N1TOT, N2TOT, N3TOT );

			if( ( cover_flags[ cnt ] == 0 ) ||
			    ( cover_flags[ cnt ] == 2 && cover_flags[ cnt_1 ] == 0 ) ) {

				get_geometry( i, j, k, pos, ncurr, geom );

				/* Make local copy for efficiency's sake:  */ 
				pl = p_L[i][j][k][d];
				pr = p_R[i][j][k][d];

				/* Calculate the fluxes and conserved variables */
				get_state( pl, geom, &ql );		// phys.c
				get_state( pr, geom, &qr );
				primtoflux( pl, &ql,   0, geom, Ul );	// phys.c
				primtoflux( pr, &qr,   0, geom, Ur );

				/* calculate maximum and minimum speeds */
				cmax = fabs( MAX( MAX( 0.,  cmaxl ),  cmaxr ) );

				/* Lax-Friedrichs flux */
				PLOOP {
					F[i][j][k][d][l] = 0.5*(Fl[l] + Fr[l] - ctop*(Ur[l] - Ul[l]));
				}
			}

		}}}  /* ijk loop */

	}   /* direction loop */ 


	TRACE_END;

	return;
}





/******************************************************************************************/
/******************************************************************************************
  recover_primitives():
 ------------------
     -- responsible for returning with a reasonable set of primitive variables; 
     -- first performs the inversion, then sees if they are valid (e.g. good
        exit status or that gamma is not too large)
     -- using gamma > GAMMAMAX2  for unphysicality test;
******************************************************************************************/

static void recover_primitives( int i, int j, int k, double *pf, double *U, 
			       struct of_geom *geom, struct of_geom *geom_old, struct of_coord *coords )
{
	ret = Utoprim_2d_fast( U, geom, pf, &gamma, &bsq ); // utoprim_2d_fast.c

	/* Retry with different inverter */
	if( ret ) {
		ret = Utoprim_1d( U, geom->gcov, geom->gcon, geom->g, pf, &gamma, &bsq ); // utoprim_1d.c
	}

	return;
}



