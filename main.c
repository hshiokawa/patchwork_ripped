#include "decs.h"
#include "defs.h"
#include "patchwork.h"
#include "patchwork_defs.h"

#include <time.h>


#define MAXFILENAME (1000)
#define CPUTXT ".%04d"


static void set_initial_globals( void );
static void set_arrays( void ) ;
static void global_init( void );



int main( int argc, char *argv[] )
{
	int i, j, k, l;

	time_t time1, time2;
	double nsecs;


	/* Initialize MPI structures and domain decomp for each patch.
	   Share property of all patches among all cpus.
	   Setup intra/inter-patcg data communication structures. */
	setup_mpi_patchwork( argc, argv );	// harm_mpi.c

	/* set various global var's based on init()/restart_init() */
	global_init();				// main.c


	/* Main routine */
        while( t < GridLength[0] ) {

		time2 = time(NULL);
		nsecs = difftime( time2, time1 );

		if( (myid == printer_pid) && sync_patches ) {
			fprintf(stdout,"timestep %10d %20g %16.6e  %16.6e  :  %26.16e %26.16e %26.16e \n",
						nstep, nsecs, t, dx[0], M_tot, E_tot, L_tot);
			fflush(stdout);
		}

		/* Step variables forward in time */
		step_ch(); // step_ch.c

		nstep++;

		/* Output events periodic in time */
		diag( OUT_EVERY );
		OUT_LOOP {
			if( (t + 1.e-10) >= T_out[i] ) { diag(i); }
		}

	}


	myexit(0);

}





/******************************************************************************************
  set_arrays(): 
 --------------
      -- Initializes all arrays to zero , and initializes recon_type[] to default 
         reconstruction method specified by RECON_TYPE_CHOICE ;
      -- should not require anything to be set (e.g. via init()) 
           but the pre-compiler macros ; 
******************************************************************************************/
void set_arrays( void ) 
{ 
	int i,j,k,d,l; 

	alloc_grid_arrays(); // in decomp.c

	t = 0.;
	t_old = 0.;

	ALL_LOOP PLOOP {       p[i][j][k][l] = 0.; }
	ALL_LOOP FACE_LOOP PLOOP  {   F[i][j][k][d][l] = 0.; }

	sprintf( DIR_out[ OUT_TRAJ    ], "/data/bh-home/hotaka/patchwork_jhu_08_sedov/execute/dumps");

}

