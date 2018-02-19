
#include "decs.h"
#include "harm_mpi.h"
#include "patchwork.h"
#include "patchwork_defs.h"


#define BASE_TAG_BOUNDS ( (PATCH_ID + 2)*10*MAXCPUS )


static int    *bc_pflag_recv[NDIM][2];  /* pointers to pflag mpi arrays used to receive BC data */
void setup_mpibc_info( void );
void setup_mpibc_info_sym( void );
void get_cpupos_pid( int patch_rank, int tmp_cpupos[NDIM] );
int init_MPI( int argc, char *argv[] );




/*********************************************************************************
  setup_mpi():
 -------------
     -- collects a given grid function to the master node from all 
        the cells from every domain; 
     -- this means that the master node must have enough memory to allocate 
        a double array with as many elements as there are total points in the 
        simulation (n_cells_global) 
     -- the collection is done using a simple loop through the cpus;
*********************************************************************************/
void setup_mpi_patchwork( int argc, char *argv[] )
{
	int i; 

	void print_pid_info( void );


	/* Initialize MPI run */
	init_MPI( argc, argv );		// harm_mpi.c

	/* Assign pids to patches. Setup patchwork structures. Fill structure "patches" (see patchwork.h). */
	if( myid == 0 ) { fprintf(stdout, " - identify_patches()\n"); fflush(stdout); }
	identify_patches( PATCH_TAG );	// patchwork_setup.c

	/* set printer_pid, io_pid, ...etc. & out_pid[] */
	set_pids();

	/* secondary initialization procedure, e.g., make MPI comm groups for each patch. */
	if( myid == printer_pid ) { fprintf(stdout, " - organize_patches()\n"); fflush(stdout); }
	organize_patches();		// patchwork_setup.c

	/* Setup the domain decomposition grid parameter */
	if( myid == printer_pid ) { fprintf(stdout, " - set_mpi_grid()\n"); fflush(stdout); }
	set_mpi_grid();			// harm_mpi.c

	/* Set boundary condition information (e.g. neighbors, etc.) */
	if( myid == printer_pid ) { fprintf(stdout, " - set_mpibc_info()\n"); fflush(stdout); }
	setup_mpibc_info();		// harm_mpi.c

	/* */
	//if( myid == printer_pid ) { fprintf(stdout, "set_mpi_misc()\n"); fflush(stdout); }
	// set_mpi_misc();		// harm_mpi.c, LATER

	/* find number of points that data exchange is needed for each proc. */
	if( myid == printer_pid ) { fprintf(stdout, " - setup_data_exchange()\n"); fflush(stdout); }
	setup_data_exchange();		// patchwork_setup.c

	/* Decide procs in other patches that are in charge of me in inter-patch commu and vice versa */
	if( myid == printer_pid ) { fprintf(stdout, " - setup_charge_relation()\n\n"); fflush(stdout); }
	setup_charge_relation();	// patchwork_setup.c

	/* Allocate variables used in intra-patch communications in connect_patches_*() */
	if( myid == printer_pid ) { fprintf(stdout, " - setup_intra_patch_comm()\n\n"); fflush(stdout); }
	setup_intra_patch_comm();	// patchwork_setup.c

	/* Print out established MPI/patchwork information. */
	print_pid_info();		// patchwork_setup.c


	return;
}





/*********************************************************************************/
/*********************************************************************************
  init_MPI():
  ----------
     -- initialize MPI run; 
     -- set global MPI parametesr, CPU's identity, etc. 
*********************************************************************************/
int init_MPI( int argc, char *argv[] )
{
 	int i;

	TRACE_BEG;


	//fprintf(stdout, "Begin: init_MPI\n");
	//fflush(stdout);

	/* not expecting MPI command line arguments so this is a formality; */
	//fprintf(stdout, "MPI_Init....\n"); fflush(stdout);
	MPI_Init( &argc, &argv );

	// get total number of cpus that will be used in this simulation
	//fprintf(stdout, "MPI_Comm_size....\n"); fflush(stdout);
	MPI_Comm_size( MPI_COMM_WORLD, &numprocs );

	// get this cpu's ID number in entire simulation group
	//fprintf(stdout, "MPI_Comm_rank....\n"); fflush(stdout);
	MPI_Comm_rank( MPI_COMM_WORLD, &myid );

	// get cpu name (e.g. hostname)
	//fprintf(stdout, "MPI_Get_processor_name...\n"); fflush(stdout);
	sprintf( myidtxt, CPUTXT, myid );
	MPI_Get_processor_name( processor_name, &procnamelen );


	if( MAXCPUS < numprocs ) {
		if( myid == printer_pid ) {
			fprintf(stdout, "Must increase MAXCPUS in global.h, %d is too many\n", numprocs);
			fflush(stdout);
		}
		fail( FAIL_MPI_BASIC, 0 ); 
	}

	if( myid == printer_pid ) { 
		fprintf(stdout, "\n\n\n - init_MPI(): Success\n");
		fflush(stdout);
	}


	TRACE_END;

	return (0);
}




/*********************************************************************************
  get_cpupos_pid():
  ----------
     -- Given the processor's local rank in patch,
        return the MPI position of the point's subgrid in the patch domain
     -- subdomains are distributed fastest in x1, and slowest in x3 since x1 is 
        always used 
*********************************************************************************/
void get_cpupos_pid( int tmp_patch_rank, int tmp_cpupos[NDIM] ) 
{ 
	TRACE_BEG;


	tmp_cpupos[0] = 0;
	tmp_cpupos[1] = ( tmp_patch_rank                       ) % ncpux[1] ;
	tmp_cpupos[2] = ( tmp_patch_rank /  ncpux[1]           ) % ncpux[2] ; 
	tmp_cpupos[3] = ( tmp_patch_rank / (ncpux[1]*ncpux[2]) )            ;


	TRACE_END;

	return;
}





// Take cpupos and return its rank in patch (patch_rank)
int cpupos_to_patch_rank( int *tmpcpupos ) 
{ 
	int i, tmp_patch_rank;

	TRACE_BEG;


	/* Check to see if this is a valid position: */
	for( i=1; i<NDIM; i++ ) { 
		if( (tmpcpupos[i] < 0) || (tmpcpupos[i] >= ncpux[i]) ) { 
			return( BAD_PID ); // general failure
		}
	}

	tmp_patch_rank =  tmpcpupos[1] + ncpux[1]*(tmpcpupos[2] + tmpcpupos[3]*ncpux[2]);


	TRACE_END;

	return( tmp_patch_rank );
}





/*********************************************************************************
  get_global_ijk():
  ----------------
     -- return with a cell's GLOBAL indices that correspond given 
        the LOCAL ones:
*********************************************************************************/

void get_global_ijk( int i, int j, int k, int  *i_glob, int *j_glob, int *k_glob )
{
	*i_glob = i + globalpos[1] ;
	*j_glob = j + globalpos[2] ;
	*k_glob = k + globalpos[3] ;

	return;
}





/*********************************************************************************
  mpi_global_max():
  ----------
     -- finds the maximum value among all domains of a scalar quantity;
     -- routine is called with a pointer to the domain's maximum value; 
     -- on exit, pointer is to the global maximum value;
*********************************************************************************/

void mpi_global_max( double *maxval )
{
	double sendbuf;

	TRACE_BEG;


	sendbuf = *maxval;
	exit_status();
	MPI_Allreduce( &sendbuf, maxval, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );


	TRACE_END;

	return;
}






/*********************************************************************************
  mpi_global_vmax():
  ----------
     -- finds the maximum value among all domains of a vector quantity;
     -- routine is called with a pointer to the domain's maximum value; 
     -- on exit, pointer is to the global maximum values;
*********************************************************************************/

void mpi_global_vmax( int n, double *maxval )
{
	int i;
	double *sendbuf;

	TRACE_BEG;


	sendbuf = (double *) calloc( n, sizeof(double) );
	if( sendbuf == NULL ) {
		fprintf(stderr,"mpi_global_vamx(): Cannot allocate sendbuf !! \n");
		fflush(stderr);
		fail( FAIL_MPI_BASIC, 0 );
	}

	exit_status();

	for( i=0; i< n; i++ ) { sendbuf[i] = maxval[i]; }

	MPI_Allreduce( sendbuf, maxval, n, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

	free( sendbuf );


	TRACE_END;

	return;
}





/*********************************************************************************
  myexit():
  ----------
     -- calls final MPI routines before exiting simulation;
*********************************************************************************/
void myexit( int ret )
{
	int i, j, k;

	extern void free_global_arrays( void );

	TRACE_BEG;


	/* Free dynamically allocated arrays : */
	for( i=1; i< NDIM; i++ ) for( j=0; j< 2; j++ ) { 
		if( bc_data_send[ i][j] != NULL ) { free( bc_data_send[ i][j] ); }
		if( bc_data_recv[ i][j] != NULL ) { free( bc_data_recv[ i][j] ); }
		if( bc_pflag_send[i][j] != NULL ) { free( bc_pflag_send[i][j] ); }
		if( bc_pflag_recv[i][j] != NULL ) { free( bc_pflag_recv[i][j] ); }
	}


	/* Finish up MPI  */
	if( ret ) { MPI_Abort( MPI_COMM_WORLD, ret ); }

	MPI_Barrier( MPI_COMM_WORLD ); /*  Required!  */
	MPI_Finalize();


	TRACE_END;

	exit( ret ) ; 

	return;
}
