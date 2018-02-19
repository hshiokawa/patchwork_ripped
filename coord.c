
#include "decs.h"




/***************************************************************************
    coord():
    -------
       -- given the indices i,j,k and location in the cell, return with 
          the values of xp1,xp2,xp3 there;  
       -- the locations are defined by : 

            -----------------------       
            |                     |   xp2 ^    
            |                     |       |       
            |FACE1   CENT         |       |    
            |                     |       |    
            |CORN    FACE2        |       ---> 
            ----------------------          xp1
   

            -----------------------   
            |                     |   xp2 ^    
            |                     |       |    
            |FACE3   CENT         |       |    
            |                     |       |    
            |CORN    FACE2        |       ---> 
            ----------------------          xp3
        

***************************************************************************/

void coord(int i, int j, int k, int loc, double *xp)
{
	int ig, jg, kg;


	get_global_ijk( i, j, k, &ig, &jg, &kg ); // harm_mpi.c

	xp[0] = t; 

	switch( loc ) { 
  
		case FACE1 : 
			xp[1] = startx[1] +  ig       *dx[1];
			xp[2] = startx[2] + (jg + 0.5)*dx[2];
			xp[3] = startx[3] + (kg + 0.5)*dx[3];
			return ;

		case FACE2 : 
			xp[1] = startx[1] + (ig + 0.5)*dx[1];
			xp[2] = startx[2] +  jg       *dx[2];
			xp[3] = startx[3] + (kg + 0.5)*dx[3];
			return ;

		case FACE3 : 
			xp[1] = startx[1] + (ig + 0.5)*dx[1];
			xp[2] = startx[2] + (jg + 0.5)*dx[2];
			xp[3] = startx[3] +  kg       *dx[3];
			return ;

		case CENT :
			xp[1] = startx[1] + (ig + 0.5)*dx[1];
			xp[2] = startx[2] + (jg + 0.5)*dx[2];
			xp[3] = startx[3] + (kg + 0.5)*dx[3];
			return ;

		case CORN :
			xp[1] = startx[1] + ig*dx[1];
			xp[2] = startx[2] + jg*dx[2];
			xp[3] = startx[3] + kg*dx[3];
			return ;

		default : 
			fprintf(stderr,"coord(): Bad location value :  %d \n", loc); 
			fflush(stderr);
			fail(FAIL_BASIC,0);
			return ;
	}
}





/***********************************************************************
  coord_of_xp():
  --------------------
   -- calculates rest of of_coord structure from xp 
***********************************************************************/

void coord_of_xp( double *xp, struct of_coord *coords )
{
	int i, j;

	TRACE_BEG;


	coords->xp[0] = xp[0];
	coords->xp[1] = xp[1];
	coords->xp[2] = xp[2];
	coords->xp[3] = xp[3];


	x_of_xp( coords->x, coords->xp ); // coord.c

	xspher_of_xcart( coords->xspher, coords->xcart ); // coord.c

	dx_dxp_calc( coords->x, coords->xp, coords->dx_dxp ); // coord.c

	dxcart_dxp_calc( coords->xcart, coords->xp, coords->x, coords->dxcart_dxp ); // coord.c

	// assumes the only non-diagonal terms are dxdxp[1][2], dxdxp[2][1]
	coords->det_dx_dxp = det_dx_dxp_calc2_default( coords->dx_dxp );  // coord.c


	/* Determine the floor state: */
	#if( FLOOR_DIM == 0 ) 
	coords->rhoflr = RHOMIN;
	coords->uuflr  = UUMIN;
	#else
	coords->rhoflr = RHOMIN*pow(coords->xspher[RR], RHOPOWER);
	coords->uuflr  =  UUMIN*pow(coords->xspher[RR],  UUPOWER);
	#endif


	TRACE_END;

	return;
}





/***********************************************************************
  calc_all_coord():
  --------------------
   -- loads the array of of_coord structures;
***********************************************************************/

void calc_all_coord( int n, double t_now )
{
	int i, j, k, pos;
	struct of_coord *coords;

	TRACE_BEG;


	ALL_LOOP POSLOOP {
		get_coord( i, j, k, pos, n, coords );	// decs.h

		coords->i   = i  ;
		coords->j   = j  ;
		coords->k   = k  ;
		coords->pos = pos;

		coord( i, j, k, pos, coords->xp );	// coord.c
		coord_of_xp( coords->xp, coords );	// coord.c

	}


	TRACE_END;

	return;
}





/*============================================================*/
/* x and xp relations */

void x_of_xp( double *x, double *xp ) 
{

	#if( COORD_TYPE_CHOICE == COORD_IDENTITY )
	x[0] = xp[0]; 
	x[1] = xp[1]; 
	x[2] = xp[2]; 
	x[3] = xp[3]; 

	#elif( COORD_TYPE_CHOICE == COORD_LOG_R )
	x[0] = xp[0]; 
	x[1] = exp(xp[1]);
	x[2] = xp[2]; 
	x[3] = xp[3]; 

	#elif( COORD_TYPE_CHOICE == COORD_LOG_R_SINH_Z )
	x[0] = xp[0];
	x[1] = R0 + exp(xp[1]); 
	x[2] = Z0*sinh(xp[2]) + Zd;
	x[3] = xp[3];

	#else
	fprintf(stderr,"x_of_xp(): Invalid value of COORD_TYPE_CHOICE : %d \n", COORD_TYPE_CHOICE);
	fflush(stderr);
	fail( FAIL_BASIC, 0 );

	#endif

	return;
}





void xcart_of_xp( double *xcart, double *xp ) 
{
	double x[NDIM];

	x_of_xp( x, xp );	// coord.c
	xcart_of_x( xcart, x );	// coord.c

	return;
}




void dxcart_dxp_calc( double *xcart, double *xp, double *x, double dxcart_dxp[][NDIM] )
{
	int i, j, l;
	double dxcart_dx[NDIM][NDIM], dx_dxp[NDIM][NDIM], tmp;

	dxcart_dx_calc( xcart, x, dxcart_dx );	// coord.c
	dx_dxp_calc( x, xp, dx_dxp );		// coord.c

	for( i=0; i< NDIM; i++ ) {
		for( j=0; j< NDIM; j++ ) {
			tmp = 0.;
			for( l=0; l< NDIM; l++ ) {
				tmp += dxcart_dx[i][l]*dx_dxp[l][j];
			}
			dxcart_dxp[i][j] = tmp;
		}
	}

	return;
}





void dxcart_dxspher_calc( double *xcart, double *xspher, double dxcart_dxspher[][NDIM] )
{
	dxcart_dxspher[0][0] = 1.;
	dxcart_dxspher[0][1] = 0.;
	dxcart_dxspher[0][2] = 0.;
	dxcart_dxspher[0][3] = 0.;
	dxcart_dxspher[1][0] = 0.;
	dxcart_dxspher[1][1] =            cos(xspher[3])*sin(xspher[2]);
	dxcart_dxspher[1][2] =  xspher[1]*cos(xspher[2])*cos(xspher[3]);
	dxcart_dxspher[1][3] = -xspher[1]*sin(xspher[2])*sin(xspher[3]);
	dxcart_dxspher[2][0] = 0.;
	dxcart_dxspher[2][1] =           sin(xspher[2])*sin(xspher[3]);
	dxcart_dxspher[2][2] = xspher[1]*cos(xspher[2])*sin(xspher[3]);
	dxcart_dxspher[2][3] = xspher[1]*cos(xspher[3])*sin(xspher[2]);
	dxcart_dxspher[3][0] = 0.;
	dxcart_dxspher[3][1] =            cos(xspher[2]);
	dxcart_dxspher[3][2] = -xspher[1]*sin(xspher[2]);
	dxcart_dxspher[3][3] = 0.;
}





/* xspher and xp relations */
void dxspher_dxp_calc( double *xspher, double *xp, double *x, double *xcart, double dxspher_dxp[][NDIM] )
{
	int i, j, l;
	double dxspher_dxcart[NDIM][NDIM], dxcart_dxp[NDIM][NDIM], tmp;

	dxspher_dxcart_calc( xspher, xcart, dxspher_dxcart );	// coord.c
	dxcart_dxp_calc( xcart, xp, x, dxcart_dxp );		// coord.c

	for( i=0; i< NDIM; i++ ) {
		for( j=0; j< NDIM; j++ ) {
			tmp = 0.;
			for( l=0; l< NDIM; l++ ) {
				tmp += dxspher_dxcart[i][l]*dxcart_dxp[l][j];
			}
			dxspher_dxp[i][j] = tmp;
		}
	}

	return;
}





void dxp_dxspher_calc( double *xp, double *xspher, double *x, double *xcart, double dxp_dxspher[][NDIM] )
{
	int i, j, l;
	double dxp_dxcart[NDIM][NDIM], dxcart_dxspher[NDIM][NDIM], tmp;

	dxp_dxcart_calc( xp, xcart, x, dxp_dxcart );		// coord.c
	dxcart_dxspher_calc( xcart, xspher, dxcart_dxspher );	// coord.c

	for( i=0; i< NDIM; i++ ) {
		for( j=0; j< NDIM; j++ ) {
			tmp = 0.;
			for( l=0; l< NDIM; l++ ) {
				tmp += dxp_dxcart[i][l]*dxcart_dxspher[l][j];
			}
			dxp_dxspher[i][j] = tmp;
		}
	}

	return;
}





/* Rotation matrix that rotates x by th about z-axis. */
void rotate( double *xrot, double *x, double th )
{
        int i, j;
        double Rot[NDIM][NDIM];

        /* For rotating about z-axis */
	Rot[0][0] = 1.;
	Rot[0][1] = 0.;
	Rot[0][2] = 0.;
	Rot[0][3] = 0.;

	Rot[1][0] = 0.;
        Rot[1][1] =  cos(th);
        Rot[1][2] = -sin(th);
        Rot[1][3] = 0.;

	Rot[2][0] = 0.;
        Rot[2][1] =  sin(th);
        Rot[2][2] =  cos(th);
        Rot[2][3] = 0.;

	Rot[3][0] = 0.;
        Rot[3][1] = 0.;
        Rot[3][2] = 0.;
        Rot[3][3] = 1.;

	/* For rotating about y-axis */
        //Rot[1][1] =  cos(th);
        //Rot[1][2] = 0.;
        //Rot[1][3] =  sin(th);
        //Rot[2][1] = 0.;
        //Rot[2][2] = 1.;
        //Rot[2][3] = 0.;
        //Rot[3][1] = -sin(th);
        //Rot[3][2] = 0.;
        //Rot[3][3] =  cos(th);

	/* For rotating about x-axis */
        //Rot[1][1] = 1.;
        //Rot[1][2] = 0.;
        //Rot[1][3] = 0.;
        //Rot[2][1] = 0.;
        //Rot[2][2] =  cos(th);
        //Rot[2][3] = -sin(th);
        //Rot[3][1] = 0.;
        //Rot[3][2] =  sin(th);
        //Rot[3][3] =  cos(th);

	xrot[0] = x[0];
        for( i=1; i<NDIM; i++) { xrot[i] = 0.; }

        for( i=1; i<NDIM; i++) {
        for( j=1; j<NDIM; j++) {
                xrot[i] += Rot[i][j]*x[j];
        }}

	return;
}





// find ijk with respect to patch corner
void patch_ijk_of_xp(int *ip, int *jp, int *kp, double *xp)
{
	*ip = (int) ((xp[1] - startx[1]) / dx[1] + 1000) - 1000;
	*jp = (int) ((xp[2] - startx[2]) / dx[2] + 1000) - 1000;
	*kp = (int) ((xp[3] - startx[3]) / dx[3] + 1000) - 1000;

	return;
}



