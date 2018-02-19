
#include "decs.h"
#include "patchwork.h"



static double rhomax=0., umax=0.;

/***********************************************************************************
  init(): 
  ---------
   -- primary initialization routine; 

***********************************************************************************/

void init( void )
{
	double ftmp;

	void init_data( void );


	/* Do some checks on the pre-compiler definitions : */
	ftmp = (double) N1TOT;
	ftmp *= N2TOT;
	ftmp *= N3TOT;
	ftmp *= NP;

	if( ftmp > ( (double) UINT_MAX ) ) {
		fprintf(stderr,"init_base(): Max. no. of elements too large !!! N = %g \n", ftmp);
		fflush(stderr);
		fail( FAIL_BASIC, 0 );
	}

	/* Set global grid parameters and static work arrays */
	init_base();

	calc_all_geom(); // metric.c

	record_distant_coord();	// patchwork_setup.c


	/* Set MHD grid functions */
	init_data();  


	return;
}





/***********************************************************************************
  init_base(): 
  ---------
   -- initializes grid structure and constant grid functions

***********************************************************************************/

void init_base( void )
{
	int i, j, k, l;
	double x[NDIM], xp[NDIM];
	double x1s, x1e, x2s, x2e, th_s, th_e;


	gam = 1.4;
	cour = 0.8;

	t = 0.;		// Initial time
	dx[0] = 1.e-4;	// Initial time step

	Rout = 1000.;
	Rin = 290;

	th_s = M_PI/10.;
	th_e = M_PI-M_PI/10.;

	x1s = log(Rin);
	x1e = log(Rout);
	x2s = th_s;
	x2e = th_e;

	// Initial patch center in global cartesian
	centx_cart[1] = 0.0;
	centx_cart[2] = 0.0;
	centx_cart[3] = 0.0;

	/* Length of each dimension in numerical coordinate */
	GridLength[0] = 2500.;
	GridLength[1] = x1e-x1s;
	GridLength[2] = x2e-x2s;
	GridLength[3] = 2.*M_PI;

	// Length of each side of grid of local patch in numerical coordinate
 	dx[1] = GridLength[1]/totalsize[1];
 	dx[2] = GridLength[2]/totalsize[2];
 	dx[3] = GridLength[3]/totalsize[3];
	DLOOP1 { invdx[i] = 1./dx[i] ; }
	dV = dx[1]*dx[2]*dx[3];	// Coordinate volume

	// Location of starting point of each side of local patch in numerical coordinate
	startx[0] = 0.;
	startx[1] = x1s;
	startx[2] = x2s;
	startx[3] = 0.;

	// Retreat starting point to add gohst zones
	for( i=1; i<NDIM; i++ ) {
		startx[i] -= NG*dx[i];
	}


        // Dumping fewquency
        DT_out[ OUT_ASCII ] = 50.0;
        DT_out[ OUT_HDF5  ] = 50.0;
	DT_out[ OUT_TRAJ  ] = 50.0;

        // Restart dumps frequency
	dt_restart = 2499.;


	return;
}




/***********************************************************************************
  init_data(): 
  ---------
   -- calculates the initial distribution of the MHD fields;
   -- driver routine for various initial data prescriptions;

***********************************************************************************/

void init_data( void )
{
	int i, j, k, l;
	long cnt;
	struct of_coord *coords;

	void init_prim( int i, int j, int k, double *pr );
	static void bounds_phys( double ****pr );


	/* Set the hydrodynamic quantities (rho,uu,v^i) */
	ALL_LOOP { // ALL_LOOP, just in case
		init_prim( i, j, k, p[i][j][k] );
	}


	/* Correct bad points and setup boundary values since we will require them for B^i */
	fixup( p );		// fixup.c

	connect_patches( p );	// patchwork_gl.c

	return;
}





/***********************************************************************************
  init_prim(): 
  ------------------
   -- Test 
***********************************************************************************/

void init_prim( int i, int j, int k, double *pr )
{
	int ii, jj, l;
	double vcon[NDIM], ucon[NDIM], ucon_cart[NDIM];
	struct of_coord *coords;
	struct of_geom geom_cart, *geom;
	static int first_call = 1;


	/* Set parameters used for all points */
	if( first_call ) {
		first_call = 0;
	}

	get_coord(    i, j, k, CENT, ncurr, coords ); // decs.h
	get_geometry( i, j, k, CENT, ncurr, geom   ); // decs.h
	get_geom_cart( coords->x,
		       coords->xp,
		       coords->xcart,
		       coords->xspher, &geom_cart ); // metric.c


	// transform ucon to numerical coordinate
	for( ii=0; ii< NDIM; ii++ ) { ucon[ii] = 0.; }
	for( ii=0; ii< NDIM; ii++ ) {
		for( jj=0; jj< NDIM; jj++ ) {

			ucon[ii] += ucon_cart[jj]*coords->dxp_dxcart[ii][jj];
	}}

	
	return;

}



