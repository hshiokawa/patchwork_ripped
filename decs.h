
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>



#define RUN_TAG   "TEST"
#define TRACE_CALLS (0) 
#define FAST_AND_FURIOUS    (0)    /* Whether to skip some MPI communications used to check consistency and MPI alive states that may or may not be needed */



//===================================================================
#define PATCH_TAG "ANIMALS"
#define PATCH_ID	(0)
#define PATCH_MOVE	(0) // whether this patch moves or not
#define HYDRO_ONLY	(0) // set to 1 if you want B^i = 0 and no evolution of B-field

#define GLOBAL_PATCH_ID (0) // PATCH_ID of global patch
#define N_PATCHES	(2) // number of patches

#define N1	(16)	/* number of grid in X1 direction of this cpu domain */
#define N2	(16)
#define N3	(16)
#define NCELLS	(N1*N2*N3)

#define Ncpu1	(3)	/* number of cpu in X1 direction of this patch */
#define Ncpu2	(3)
#define Ncpu3	(3)
#define Ncpu	(Ncpu1*Ncpu2*Ncpu3)


//===================================================================
/* Choose the type of metric/coordinates we will be using */
#define METRIC_MINK_CARTESIAN	   (0)	/* Flatspace in cartesian x,y,z coordinates             */
#define METRIC_MINK_CARTESIAN_TILT (1)	/* Flatspace in tilted-cartesian xt,yt,zt coordinates   */
#define METRIC_MINK_CYLINDRICAL	   (2)	/* Flatspace in cylindrical r,z,phi coordinates         */
#define METRIC_MINK_SPHERICAL	   (3)	/* Flatspace in spherical r,th,phi coordinates          */
#define METRIC_KS_CARTESIAN	   (4)	/* Kerr-Schild in tilted-cartesian xt,yt,zt coordinates */
#define METRIC_KS_CARTESIAN_TILT   (5)	/* Kerr-Schild in tilted_cartesian x,y,z coordinates    */
#define METRIC_KS_CYLINDRICAL	   (6)	/* Kerr-Schild in cylindrical r,z,phi coordinates       */
#define METRIC_KS_SPHERICAL	   (7)	/* Kerr-Schild in spherical r,th,phi coordinates        */

#define METRIC_TYPE_CHOICE ( METRIC_MINK_SPHERICAL )


//===================================================================
/* Numerical Flux method */
#define FLUX_LF      (0)
#define FLUX_HLLE    (1)

#define FLUX_TYPE_CHOICE  ( FLUX_LF )



//===================================================================
/* Selects the sequence of fixups, see ./info/fixups */
/* 3 -> like 2 but never changes the B-field components */
#define FIXUP_TREE  (3)   

                          

//===================================================================
/* Option to use the entropy evolution equation instead of total energy equation */ 
#define USE_ENTROPY_EQ   (1)       /* Whether to ever use entropy equation */
#define BETA_MIN         (1.e-2)   /* use EE fix when u < BETA_MIN * bsq */


//===================================================================
/* mnemonics for primitive vars */
#define RHO	(0)	/* rest-mass density */
#define UU	(1)	/* internal energy density */ 
#define U1	(2)	/* spatial velocity component */
#define U2	(3)	/* spatial velocity component */
#define U3	(4)	/* spatial velocity component */
#define B1	(5)	/* spatial magnetic field component */
#define B2	(6)	/* spatial magnetic field component */
#define B3 	(7)	/* spatial magnetic field component */
#define NP	(8)	/* number of EOM and # of primitive  variables */


/* Position types within a cell */
#define FACE1 (0)
#define FACE2 (1)
#define FACE3 (2)
#define CENT  (3)
#define CORN  (4)   
#define NPOS  (5)   /* Number of unique positions in a cell  */ 


//===================================================================
/* for-loop aliases */
#define LOOP       for(i=N1S;i<=N1E ;i++) for(j=N2S;j<=N2E;j++) for(k=N3S;k<=N3E;k++)
#define N1_LOOP    for(i=N1S;i<=N1E ;i++)
#define N2_LOOP    for(j=N2S;j<=N2E ;j++)
#define N3_LOOP    for(k=N3S;k<=N3E ;k++)
#define N1ALL_LOOP for(i=0  ;i<N1TOT;i++)
#define N2ALL_LOOP for(j=0  ;j<N2TOT;j++)
#define N3ALL_LOOP for(k=0  ;k<N3TOT;k++)
#define ALL_LOOP   for(i=0  ;i<N1TOT;i++) for(j=0  ;j<N2TOT;j++) for(k=0  ;k<N3TOT;k++)
#define SDLOOP1    for(i=1  ;i<NDIM ;i++)
#define POSLOOP    for(pos=0;pos<NPOS;pos++)
#define FACE_LOOP  for(d=0  ;d<(NDIM-1);d++)

/* For patchwork */
#define PATCH_LOOP for(i_patch=0; i_patch<N_PATCHES; i_patch++)


//===================================================================
/* Macros used for memory allocations */
#define ALLOC_ARRAY( _array, _npts ) {\
		if( ( (_array) = malloc( (_npts)*sizeof(*(_array)) ) ) == NULL ) {\
			fprintf(stdout,"%s(): Cannot allocate %s on line %d  of %s \n",__func__,#_array,__LINE__,__FILE__);\
			exit(0);\
		}\
		memset( (_array), 0, (_npts)*sizeof(*(_array)) );\
	}

#define ALLOC_2D_ARRAY(_array,_n1,_n2) { \
 int _i1; \
 ALLOC_ARRAY(_array,_n1); \
 for(_i1=0; _i1<_n1; _i1++) { ALLOC_ARRAY(_array[_i1],_n2); } \
 }

#define FREE(_p) { free(_p); _p = NULL; }

#define DEALLOC_ARRAY( _array, _npts ) { FREE( (_array) ); }

#define DEALLOC_2D_ARRAY(_array,_n1,_n2) { \
 int _i1;\
 for(_i1=0; _i1<_n1; _i1++) { DEALLOC_ARRAY(_array[_i1],_n2); }\
 DEALLOC_ARRAY(_array,_n1);\
 }




//==================================================================
/* Global arrays and variables */

struct of_geom {
	double gcon[NDIM][NDIM];
	double gcov[NDIM][NDIM];
	double g;
	double g_inv;
};

extern struct of_geom *geom_arr[N0_GEOM];

struct of_coord {
	int i, j, k, pos;
	double x[NDIM];
	double xp[NDIM];
	double dx_dxp[NDIM][NDIM];
	double det_dx_dxp;
	double xcart[NDIM];
};

extern struct of_coord *coord_arr[N0_COORD];


/* 123  MHD  functions per cell */
extern double    ****p		; /* space for primitive vars */
extern double   *****p_L	; /* Left state for flux at i,j,k */
extern double   *****p_R	; /* slopes */
extern double   *****F		; /* fluxes */

extern void check_boundary_pflag( int i, int j, int k );
extern void myexit( int ret );

/* patchwork_gl.c / patchwork_lg.c */
extern void connect_patches(	       double ****prim_arr			);
extern void connect_patches_direction( double ****prim_arr, int local_to_global );
extern void connect_patches_gl(	       double ****prim_arr, int i_patch		);
extern void connect_patches_lg(	       double ****prim_arr, int i_patch		);
