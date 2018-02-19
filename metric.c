#include "decs.h"
#include "metric.h"
#include "mpi.h"



void gcov_func__ks_spherical( double *xspher, double gcov[][NDIM] );
void gcon_func__ks_spherical( double *xspher, double gcon[][NDIM] );
void conn_func__ks_spherical(   struct of_coord *coords, double connp[][NDIM][NDIM] );


/********************************************************************************/
/********************************************************************************
 calc_all_geom
 -------------
  -- Driver routine for calculating the metric, inverse metric and connection 
     coefficients 

********************************************************************************/

void calc_all_geom( void )
{
	int ii, jj, kk, ll, i, j, k, l, pos;
	struct of_geom *geom;
	struct of_coord *coords;

	double connp[NDIM][NDIM][NDIM], dxcart_dxp_dxp[NDIM][NDIM][NDIM], dxspher_dxp_dxp[NDIM][NDIM][NDIM];
	double dxspher_dxp[NDIM][NDIM], dxp_dxspher[NDIM][NDIM];


	/* Set certian coordinate constants : */
	GEOM_LOOP {
		get_geometry( i, j, k, pos, ncurr, geom   );	// decs.h
		get_coord(    i, j, k, pos, ncurr, coords );	// decs.h

		gcov_func__mink_numerical( coords->x, coords->xp, coords->xcart, coords->xspher, geom->gcov );	// metric.c
		geom->g = gdet_func( geom->gcov );	// metric.c
	}

	l = 0;
	CONN2_LOOP {

		get_geometry( i, j, k, CENT, ncurr, geom   );	// decs.h
		get_coord(    i, j, k, CENT, ncurr, coords );	// decs.h
		
		conn_func__finite_diff( coords->xp, geom->gcon, conn[l], gcov_func__mink_numerical );

		l++;
	}


	/* Calculate the lapse and shift from the just calculated metric :  */
	GEOM_LOOP {
		get_geometry( i, j, k, pos, ncurr, geom ); // decs.h

		geom->ncon[0] =  sqrt(-geom->gcon[0][0]) ; // 1/alpha
		geom->alpha   =  1./geom->ncon[0];
	}


	return;
}





void gcov_func__mink_numerical( double *x, double *xp, double *xcart, double *xspher, double gcov[][NDIM] )
{
	int i, j, ii, jj;
	double dxcart_dxp[NDIM][NDIM], gcov_temp[NDIM][NDIM];


	DLOOP2 { gcov[i][j] = mink[i][j]; }


	/* Transform to numerical coordinate */
	dxcart_dxp_calc( xcart, xp, x, dxcart_dxp );	// coord.c

	DLOOP2 { gcov_temp[i][j] = 0.; }
	DLOOP2 {
	        for( ii=0; ii< NDIM; ii++ ) {
	        for( jj=0; jj< NDIM; jj++ ) {

	                gcov_temp[i][j] += gcov[ii][jj]*dxcart_dxp[ii][i]*dxcart_dxp[jj][j];
	}}}

	DLOOP2 { gcov[i][j] = gcov_temp[i][j]; }

        return;
}




void gcov_func__ks_spherical( double *xspher, double gcov[][NDIM] )
{
	double r, th;
	double t10, t12, t13, t2, t3, t4, t5, t7, t9;


	/* Calculate the KS metric using derivation with MAPLE */
	r  = xspher[RR];
	th = xspher[TH];

	t4  = cos(th);
	t12 = sin(th);
	t2  = r*r;
	t3  = asq;
	t5  = t4*t4;
	t7  = t2+t3*t5;
	t9  = M*r/t7;
	t10 = 2.0*t9;

	gcov[0][0] = t10-1.0;
	gcov[0][1] = t10;
	gcov[0][2] = 0.0;
	t13 = t12*t12;
	gcov[0][3] = -2.0*t9*a*t13;
	gcov[1][0] = t10;
	gcov[1][1] = 1.0+t10;
	gcov[1][2] = 0.0;
	gcov[1][3] = -a*gcov[1][1]*t13;
	gcov[2][0] = 0.0;
	gcov[2][1] = 0.0;
	gcov[2][2] = t7;
	gcov[2][3] = 0.0;
	gcov[3][0] = gcov[0][3];
	gcov[3][1] = gcov[1][3];
	gcov[3][2] = 0.0;
	gcov[3][3] = t13*(t2+t3*(1.0+2.0*t9*t13));

	return;
}





void gcon_func__ks_spherical( double *xspher, double gcon[][NDIM] )
{
	double r, th;
	double t1, t10, t13, t14, t2, t3, t4, t5, t8;


	/* Calculate the KS metric using derivation with MAPLE */
	r  = xspher[RR];
	th = xspher[TH];

	t4 = cos( th );
	t13 = sin( th );
	t1 = M*r;
	t2 = r*r;
	t3 = asq;
	t5 = t4*t4;
	t8 = 1/(t2+t3*t5);
	t10 = 2.0*t1*t8;

	gcon[0][0] = -1.0-t10;
	gcon[0][1] = t10;
	gcon[0][2] = 0.0;
	gcon[0][3] = 0.0;
	gcon[1][0] = t10;
	gcon[1][1] = (t2-2.0*t1+t3)*t8;
	gcon[1][2] = 0.0;
	gcon[1][3] = a*t8;
	gcon[2][0] = 0.0;
	gcon[2][1] = 0.0;
	gcon[2][2] = t8;
	gcon[2][3] = 0.0;
	gcon[3][0] = 0.0;
	gcon[3][1] = gcon[1][3];
	gcon[3][2] = 0.0;
	t14 = t13*t13;
	gcon[3][3] = t8/t14;

	return;
}

