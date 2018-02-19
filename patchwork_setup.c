
/* Written by Hotaka Shiokawa, modified on 02/19/2018 */

#include "decs.h"
#include "patchwork.h"


/* patchwork_seup.c:

	Set of routines that setup inter/intra patch communication info.
	Most of them are called in setup_mpi_patchwork() in harm_mpi.c in following
	order:

		- identify_patches(): called in setup_mpi_patchwork() in harm_mpi.c

			- Identify domain decomposition of each patch.
			- Identify which CPU belongs which patch.

		- organize_patches(): called in setup_mpi_patchwork() in harm_mpi.c

			- Defines intra-patch communicator "patch_mpi_comm" and local
			  rank in the patch "patch_rank".

		- n_pts_calc(): called in setup_mpi_patchwork() in harm_mpi.c

			- Count number of edge zones "n_pts" for each proc, and create
			  a list of n_pts for all procs.
			  n_pts is needed to decide inter-patch communication size

		- setup_charge_relation(): called in setup_mpi_patchwork() in harm_mpi.c

			- Setup inter-patch cpu-cpu relation for data exchange.
			- Each proc on local patch communicates with only one
			  specific CPU on global patch. Likewise, each proc on global
			  pach communicates with only one specific CPU per local patch.
			  setup_charge_relation() decides this inter-patch communication
			  partner for each CPU.

		- print_pid_info(): called in setup_mpi_patchwork() in harm_mpi.c

			- Print out all patchwork related information after all setup is done.

		- record_distant_coord(): called in init_base() in init.c

			- CPUs on local patch store grid coordinates of global patch
			  because global patch never moves. This can omit sending
			  list of coordinate from global CPU to local CPU at each time
			  step.


	Notes:
		- Each patch has one "decider" that deals with inter-patch information
		  exchange.

		- Lines that are filling components of "patches" are marked by *fil* in
		  front of the line.

		- Structure "I" defined in patchwork.h stores properties specific to
		  each CPU, and lines that are filling components of structure "I" are
		  marked by ** in front of the line.

		- In general, "pids" are for world rank and "rank" are for rank in local group
		  (e.g.I.patch_rank)
*/




/*****************************************************************************
  identify_patches(): (called in setup_mpi_patchwork() in harm_mpi.c)
  
	- Identify domain decomposition info of each patch.
	- Identify which CPUs belong which patch.
	- Share the identified properties among all CPUs involved in simulation.
	  For example, a cpu on patch-0 knows which CPUs belong to patch-2 at
	  the end of this routine.
	- Properties of each patch is broadcasted by "decider". Each patch has
	  one decider.
	- The identified properties are stored in structure "patches" defined
	  in patchwork.h.
	  At the end of this routine, you can access, say, total number of grids
	  in x-direction of patch-2 by patches[2].N1TOT

	Notes:
		- "*fil*" in front of a line marks where components of "patches"
		  are filled

******************************************************************************/

void identify_patches( const char *local_patch_name ) 
{
	int i, j, i_patch;
	int patch_name_length, patch_n_procs, patch_decider_pid;
	int *list_patch_n_procs;


	// Basic info of all patches, patches[N_PATCHES].
	// Every procs eventually get completed "patches"
	ALLOC_ARRAY( patches, N_PATCHES );

	// For convenience, point to memory space where properties of patch that this
	// proc belongs is going to be stored
	my_patch = &( patches[ PATCH_ID ] );


	I.i_belong_global_patch = 0;
	if( PATCH_ID == GLOBAL_PATCH_ID ) {
/**/		I.i_belong_global_patch = 1; // /**/ denotes structure "I" is filled
	}



	//----------------------------------------------------------//
	/* Fill basic components of my_patch */

/*fil*/	my_patch->id			= PATCH_ID;
/*fil*/	my_patch->n_procs		= Ncpu;  // Ncpu1*Ncpu2*Ncpu3 (total num of cpu on this patch)
/*fil*/	my_patch->N1tot			= N1TOT; // N1+2*NG (total num of grids on 1st axis of this cpu comain)
/*fil*/	my_patch->N2tot			= N2TOT;
/*fil*/	my_patch->N3tot			= N3TOT;
/*fil*/	my_patch->Ntot			= NTOT;	 // N1TOT*N2TOT*N3TOT;
/*fil*/	my_patch->move			= PATCH_MOVE; // if this patch moves or not
/*fil*/	my_patch->top_type_choice	= TOP_TYPE_CHOICE;   // coordinate topology info
/*fil*/	my_patch->coord_type_choice	= COORD_TYPE_CHOICE; // coordinate type

	patch_name_length = ((int) strlen( local_patch_name )) + 1;
/*fil*/	my_patch->name_length = patch_name_length; // name length of this patch

	ALLOC_ARRAY( my_patch->name, patch_name_length );
/*fil*/	strcpy( my_patch->name, local_patch_name ); // name of this patch



	//----------------------------------------------------------//
	/* Fill my_patch->pids, my_patch->decider_pid, and my_patch->patch_decider_rank */

	ALLOC_ARRAY( list_patch_n_procs, numprocs ); 

	// gather # of processors per patch from each processor and create
	// a sorted list "list_patch_n_procs" = e.g. [3,3,3,4,4,4,4,2,2]
	MPI_Allgather( &(my_patch->n_procs),
		       1,
		       MPI_INT,
		       //-----
		       list_patch_n_procs, // every procs has the list (Allgather!)
		       1,
		       MPI_INT,
		       MPI_COMM_WORLD );

	// Now, fill "pids" of all patches, not only my_patch, because info of all
	// patches are available, why not?
	i = 0;
	PATCH_LOOP {

		patch_n_procs = list_patch_n_procs[i]; // n_procs of patch_id = i_patch
		ALLOC_ARRAY( patches[i_patch].pids, patch_n_procs );

		for( j=0; j< patch_n_procs; j++) {
/*fil*/			patches[i_patch].pids[j] = i;
			i++;
		}

		// decider is proc of final rank of each patch (DO NOT CHANGE!!)
/*fil*/		patches[i_patch].decider_pid = i - 1;
		// therefore, rank of decider in local patch is simply n_procs - 1
/*fil*/		patches[i_patch].patch_decider_rank = patch_n_procs - 1;
	}


	I.i_am_decider = 0;
	if( myid == my_patch->decider_pid ) {
/**/		I.i_am_decider = 1;
	}

	DEALLOC_ARRAY( list_patch_n_procs, numprocs );



	//----------------------------------------------------------//
	/* Broadcast info of each patch from decider to everyone */

	PATCH_LOOP {
		patch_decider_pid = patches[ i_patch ].decider_pid;
		
		MPI_Bcast( &( patches[ i_patch ].id	 ),    1, MPI_INT, // id
			   patch_decider_pid, MPI_COMM_WORLD );
		MPI_Bcast( &( patches[ i_patch ].n_procs ),    1, MPI_INT, // n_procs
			   patch_decider_pid, MPI_COMM_WORLD );
		MPI_Bcast( &( patches[ i_patch ].N1tot ),      1, MPI_INT, // N1TOT
			   patch_decider_pid, MPI_COMM_WORLD );
		MPI_Bcast( &( patches[ i_patch ].N2tot ),      1, MPI_INT, // N2TOT
			   patch_decider_pid, MPI_COMM_WORLD );
		MPI_Bcast( &( patches[ i_patch ].N3tot ),      1, MPI_INT, // N3TOT
			   patch_decider_pid, MPI_COMM_WORLD );
		MPI_Bcast( &( patches[ i_patch ].Ntot ),       1, MPI_INT, // N1TOT*N2TOT*N3TOT
			   patch_decider_pid, MPI_COMM_WORLD );
		MPI_Bcast( &( patches[ i_patch ].move ),       1, MPI_INT, // PATCH_MOVE
			   patch_decider_pid, MPI_COMM_WORLD );
		MPI_Bcast( &( patches[ i_patch ].top_type_choice   ), 1, MPI_INT, // top_type_choice
			   patch_decider_pid, MPI_COMM_WORLD );
		MPI_Bcast( &( patches[ i_patch ].coord_type_choice ), 1, MPI_INT, // coord_type_choice
			   patch_decider_pid, MPI_COMM_WORLD );
		MPI_Bcast( &( patches[ i_patch ].name_length	   ), 1, MPI_INT, // name_length
			   patch_decider_pid, MPI_COMM_WORLD );

		patch_name_length = patches[i_patch].name_length;
		if( my_patch->id != i_patch ) { ALLOC_ARRAY( patches[i_patch].name, patch_name_length ); }

		MPI_Bcast( patches[i_patch].name, patch_name_length, MPI_CHAR,	// name
			   patch_decider_pid, MPI_COMM_WORLD );
	}


	MPI_Barrier( MPI_COMM_WORLD );

	return;
}





/*****************************************************************************
  organize_patches(): (called in setup_mpi_patchwork() in harm_mpi.c)

	- Create MPI communicators "patch_mpi_comm" and "decider_mpi_comm".

		patch_mpi_comm:
			- Communicator for intra-patch communication
			- "patch_rank" is rank of a cpu in the patch it belongs

		decider_mpi_comm:
			- Communicator for "deciders" that deal with inter-patch
			  informaion exchange (this communicator has not been used
			  in current version yet(?)).
			- "deciders" are already assigned in
			  identify_patches() (patchwork_setup.c).

******************************************************************************/

void organize_patches( void ) 
{
	int i_patch, num_decider, *decider_pids;



	MPI_Comm_group( MPI_COMM_WORLD, &world_group ); // get world_goup from communicator MPI_COMM_WORLD


	//----------------------------------------------------------//
	/* Create patch_mpi_comm */
	MPI_Group_incl( world_group, my_patch->n_procs, my_patch->pids, &patch_mpi_group );
	// patch_mpi_comm = communicator of local patch procs
	MPI_Comm_create( MPI_COMM_WORLD, patch_mpi_group, &patch_mpi_comm );
	// patch_rank = rank of this proc in local patch group
/**/	MPI_Comm_rank( patch_mpi_comm, &(I.patch_rank) );



	//----------------------------------------------------------//
	/* Create decider_mpi_comm */

	ALLOC_ARRAY( decider_pids, N_PATCHES );
	
	// list of world rank of decider pids
	PATCH_LOOP { decider_pids[ i_patch ] = patches[ i_patch ].decider_pid; }

	MPI_Group_incl( world_group, N_PATCHES, decider_pids, &decider_mpi_group );
	// decider_mpi_comm = communicator of all decider
	MPI_Comm_create( MPI_COMM_WORLD, decider_mpi_group, &decider_mpi_comm );
	// decider_rank = rank of decider in decider_mpi_group ( which is equal to id of patch of decider)
	if( I.i_am_decider ) { MPI_Comm_rank( decider_mpi_comm, &decider_rank ); }

	// set decider_master
	decider_master_rank = GLOBAL_PATCH_ID; // MUST be set to the one on global patch

	I.i_am_decider_master = 0;
	if( I.i_am_decider ) {
		if( decider_rank == decider_master_rank ) {
/**/			I.i_am_decider_master = 1;
		}
	}

	// some chekings
	if( I.i_am_decider ) {
		MPI_Comm_size( decider_mpi_comm, &num_decider ); // get the number of decider procs
		if( num_decider != N_PATCHES ) {
			fprintf( stderr, "Number of decider procs must be equal to number of patches.\n" );
		}
	}



	DEALLOC_ARRAY( decider_pids, N_PATCHES );

	MPI_Barrier( MPI_COMM_WORLD );

	return;
}





/*****************************************************************************
  setup_data_exchange(): (called in setup_mpi_patchwork() in harm_mpi.c)

	- Allocate cover flag arrays etc.

******************************************************************************/

void setup_data_exchange( void )
{
	MPI_Barrier( MPI_COMM_WORLD );
}





/*****************************************************************************
  setup_charge_relation(): (called in setup_mpi_patchwork() in harm_mpi.c)

	- Setup inter-patch cpu-cpu relation for data exchange.

		- Each proc on local patch communicates with only one
		  specific CPU on global patch. This CPU on global patch
		  receives list of coordinate of the local grid's ghost zones,
		  and gathers/manages required data to send back interpolated
		  data to the local proc.

		- Ghost zones of local patch is usually edge-grids on
		  edge-procs, but sometimes entire patch needs to be ghost
		  zones if, e.g. first time step of adding a new local patch
		  to global patch during a simulation

		- Likewise, each proc on global pach communicates with only
		  one specific CPU per local patch.

		- setup_charge_relation() assigns this specific partner for
		  inter patch communication

	- local_charges_pid[ i_patch ]: (often called "advisor")

		- WORLD RANK of CPU on patch of PATCH_ID = i_patch that "I" am
		  asking for interpolated data.

		- Allocated for total number of patch.
		  Say, total number of patch (N_PATCHES) is 3, and global patch
		  ID (GLOBAL_PATCH_ID) is 0 by default.
		  For a CPU on global patch, local_charges_pid is something
		  like:
			local_charges_pid[0] = -1
			local_charges_pid[1] = 10
			local_charges_pid[2] = 20.

		  This means that, if this global CPU wants to get
		  interpolated data from patch-1, it always contact CPU on
		  patch-1 whose world rank is 10, and for patch-2 contact
		  CPU on patch-2 of world rank = 20.
		  local_charges_pid[0] is set to -1 (invalid) since it doesn't
		  contact its own patch.

		  For a CPU on local patch, say patch-1, local_charges_pid is
		  something like:
			local_charges_pid[0] = 5
			local_charges_pid[1] = -1
			local_charges_pid[2] = -1.

		  This means that, if this CPU wants to get interpolated data
		  from global patch, it always ask a CPU on global patch of
		  world rank = 5 for it. No inter patch communication between
		  2 local patches is integrated yet, so
		  local_charges_pid[2] = -1 (invalid).

		- local_charges_pid are often called "advisor" because they
		  take care of your request.

	- distant_charges[ i_patch ]: (often called "students")

		- WORLD RANK of CPU on patch of PATCH_ID = i_patch that "I" am
		  providing interpolated data when asked.

		- Allocated for total number of patch.

		- Example:
			Total number of patch (N_PATCHES) is 3, and global
			patch ID (GLOBAL_PATCH_ID) is 0 by default.

			For a CPU on global patch, distant_charges[] is
			something like:

				distant_charges[0].n_procs = 0
						  .n_procs_max = 0
						  .pids = []
				distant_charges[1].n_procs = 2
						  .n_procs_max = 2
						  .pids = [13,16]
				distant_charges[2].n_procs = 3
						  .n_procs_max = 4
						  .pids = [24,28,32]

			This means that CPU of world rank 13 and 16 on patch-1
			and 24, 28 and 32 on patch-2 will ask the global CPU for
			interpolated data for their ghost zones (often edge zones).
			"n_procs_max" is max possible number of procs that the
			global CPU may need to take care for each local patch.
			This number is useful for setting upper limit of
			for-loops.

			For a CPU on local patch patch-1, distant_charges[] is
			something like:
				distant_charges[0].n_procs = 2
						  .n_procs_max = 3
						  .pids = [3,5]
				distant_charges[1].n_procs = 0
						  .n_procs_max = 0
						  .pids = []
				distant_charges[2].n_procs = 0
						  .n_procs_max = 0
						  .pids = []

			This means CPU of world rank 3 and 5 on global patch
			ask the local CPU for interpolated data for their
			"covered zones". Inter-patch communication between 2
			local patches is not integrated in the code yet.

		- CPUs of distant_charges are often called "students" because
		  they ask for care.

	Notes:
		- distant_charges[].*****coords is filled in
		  record_distant_coord() (patchwork_setup.c).

******************************************************************************/

void setup_charge_relation( void )
{
	int i, i_patch;
	int n_multiple, n_student, advisor, student, student_pid, n_pts;

	struct of_patch *patch_global, *patch_local;



	//----------------------------------------------------------//
	/* Assign local_charges_pid (advisor) to each proc.
	   local_charges_pid is world-rank of a proc in another patch
	   that helps communications of me and the other patch. */

	if( I.i_belong_global_patch ) { // global patch

		ALLOC_ARRAY( local_charges_pid, N_PATCHES );

		PATCH_LOOP {

			if( i_patch != GLOBAL_PATCH_ID ) {

				patch_local = &( patches[ i_patch ] );

				// advisor = patch_rank of cpu in a local patch that I ask
				//      for data interpolation or questions about the local patch
				advisor = I.patch_rank // Let's choose "my patch_rank" = "advisor's patch_rank".
					// But, if "my n_procs" > "local patch's n_procs", need some overlaps
					- (int) (I.patch_rank / patch_local->n_procs) * patch_local->n_procs;

				local_charges_pid[ i_patch ] = patch_local->pids[ advisor ]; // world_rank of advisor
			} else {
				local_charges_pid[ i_patch ] = -1; // never access
			}
		}
	}


	if( I.i_belong_global_patch == 0 /* local patch */) {

		ALLOC_ARRAY( local_charges_pid, N_PATCHES );

		PATCH_LOOP {

			if( i_patch == GLOBAL_PATCH_ID ) {

				patch_global = &( patches[ i_patch ] );

				// advisor = patch_rank of cpu in global patch that I ask
				//      for data interpolation or questions about the global patch
				advisor = I.patch_rank // Let's choose "my patch_rank" = "advisor's patch_rank".
					// But, if "my patch_n_procs" > "global patch's n_procs", need some overlaps
					- (int) (I.patch_rank / patch_global->n_procs) * patch_global->n_procs;

				local_charges_pid[ i_patch ] = patch_global->pids[ advisor ]; // world_rank of advisor
			} else {
				local_charges_pid[ i_patch ] = -1; // never access
			}
		}
	}



	//----------------------------------------------------------//
	/* Setup structure distant_charges (students) for each proc.
	   distant_charges is info of a proc in another patch
	   that I help communications of it and my patch. */

	ALLOC_ARRAY( distant_charges, N_PATCHES );


	if( I.i_belong_global_patch ) { // global patch

		// deal with local patches, one by one
		PATCH_LOOP {

			if( i_patch != GLOBAL_PATCH_ID ) {

				patch_local = &( patches[ i_patch ] );

				/* whose advisor of local patch am I ? */

				// if local patch has more procs than n_procs of global patch,
				// I may need to take care of multiple local cpus
				n_multiple = (int) (patch_local->n_procs / my_patch->n_procs);
								// note here my_patch = global patch

				distant_charges[ i_patch ].n_procs_max = n_multiple + 1; // possible max # of student I 'may'
										 	 // need to take care. At least 1.

				// check if this proc really need to take care of
				// n_procs_max local procs or not
				for( i=0; i<= n_multiple; i++ ) {

				        // student: patch_rank of local patch that I take care of
				        student = I.patch_rank + i*my_patch->n_procs;
				        n_student = i + 1;

				        if( student >= patch_local->n_procs ) {
						n_student--;
						break;
				        }
				} // obtained actual number of local procs I take care of, n_student

												 // Array of:
				ALLOC_ARRAY( distant_charges[ i_patch ].coord_list, n_student ); // list of ghost zones of each student
				ALLOC_ARRAY( distant_charges[ i_patch ].n_pts	  , n_student ); // # of ghost zones of each student
				ALLOC_ARRAY( distant_charges[ i_patch ].pids      , n_student ); // world-rank of each student

				for( i=0; i< n_student; i++ ) {

					student = I.patch_rank + i*my_patch->n_procs; // patch_rank of student
					student_pid = patch_local->pids[ student ];  // world rank of the student

					distant_charges[ i_patch ].pids[ i ] = student_pid;

					// number of ghost zones of the student
					//n_pts = list_n_pts[ student_pid ];
					//distant_charges[ i_patch ].n_pts[ i ] = n_pts; // TBD

					// allocate ghost zone coordinate storage.
					// allocate maximum possible number.
					//ALLOC_ARRAY( distant_charges[ i_patch ].coord_list[ i ], 3*patch_local->Ntot ); // 3 for i,j,k // TBD
				}

				distant_charges[ i_patch ].n_procs = n_student;

			} else {
				distant_charges[ i_patch ].n_procs     = 0;
				distant_charges[ i_patch ].n_procs_max = 0;
			}

		} // end PATCH_LOOP

	}



	if( I.i_belong_global_patch == 0 ) { // if I belong local patch

		PATCH_LOOP {

		} // end PATCH_LOOP

	}



	MPI_Barrier( MPI_COMM_WORLD );

	return;
}





/*****************************************************************************
  setup_intra_patch_comm(): (called in setup_mpi_patchwork() in harm_mpi.c)

	- Allocate intra-patch communication related variables in
	  connect_patches_*.

******************************************************************************/

void setup_intra_patch_comm( void )
{
	int i_patch, n_student_max;

	ALLOC_ARRAY( intra_patch_comm, N_PATCHES );

	PATCH_LOOP {
		n_student_max = distant_charges[ i_patch ].n_procs_max;

		ALLOC_2D_ARRAY( intra_patch_comm[ i_patch ].assign_n_pts_3, 	 n_student_max, Ncpu );
		ALLOC_2D_ARRAY( intra_patch_comm[ i_patch ].list_assign_n_pts_3, n_student_max, Ncpu*Ncpu );
		ALLOC_2D_ARRAY( intra_patch_comm[ i_patch ].assign_index_map,	 n_student_max, Ncpu );
		ALLOC_2D_ARRAY( intra_patch_comm[ i_patch ].assign_coord_list,	 n_student_max, Ncpu );
		ALLOC_2D_ARRAY( intra_patch_comm[ i_patch ].data_recv,		 n_student_max, Ncpu );
	}
}





/*****************************************************************************
  record_distant_coord: (called in init_base() in init.c)

	- CPUs store grid coordinates of their distant_charges (students) for
	  stationary patches (note global patch never moves). This can omit
	  sending list of coordinate from global CPU to local CPU at each time
	  step.

	- The coordinate list is stored in
	  distant_charges[student's PATCH_ID].*****coord.

	- coords[ # of student on student patch ]
		[ N1TOT (N1+2*NG) of student cpu ]
		[ N2TOT of student cpu ]
		[ N3TOT of student cpu ]
		[ NDIM (for 4D coordinate values) ]

	- This routine must be called 'after' grid coordinates are set which
	  usually happens in init_base() in init.c

******************************************************************************/

void record_distant_coord( void )
{
	int N1TOT_s, N2TOT_s, N3TOT_s, i, j, k, d, ip, i_patch;
	int n_student, student_patch_rank, tag_send, tag_recv;
	long cnt, n_pts;
	double *data_send, *data_recv;

	struct of_patch *patch_student;
	struct of_charges *student;
	struct of_coord *coords;



	//----------------------------------------------------------//
	/* Send grid coordinate of procs on global patch to their
	   advisors on local patch */

	if( I.i_belong_global_patch ) {

	} // end of global patch routine


	if( I.i_belong_global_patch == 0 ) { // local patch

		student = &( distant_charges[ GLOBAL_PATCH_ID ] );
		n_student = student->n_procs;

		if( n_student > 0 ) { // if I have student on global patch

			patch_student = &( patches[ GLOBAL_PATCH_ID ] ); // just for convenience

			N1TOT_s = patch_student->N1tot; // "_s" denotes "student (on global patch)"
			N2TOT_s = patch_student->N2tot;
			N3TOT_s = patch_student->N3tot;

			n_pts = patch_student->Ntot * NDIM; // = N1TOT_s * N2TOT_s * N3TOT_s * NDIM
			ALLOC_ARRAY( data_recv, n_pts );

			// allocate memory space to store grid coordinates of students
			ALLOC_5D_ARRAY( student->coord,
					n_student,
					N1TOT_s,
					N2TOT_s,
					N3TOT_s,
					NDIM );

			for( ip=0; ip< n_student; ip++ ) { // receive coordinate student by student

				tag_recv = patch_student->n_procs * my_patch->id + student->pids[ ip ];

				MPI_Recv( data_recv,
					  n_pts,
					  MPI_DOUBLE,
					  //-------
					  student->pids[ ip ],
					  tag_recv,
					  MPI_COMM_WORLD,
					  MPI_STATUS_IGNORE );

				// unpack the coordinate list
				cnt = 0;
				for( i=0; i< N1TOT_s; i++ ) {
				for( j=0; j< N2TOT_s; j++ ) {
				for( k=0; k< N3TOT_s; k++ ) {
				for( d=0; d< NDIM   ; d++ ) {

					student->coord[ ip ][i][j][k][d] = data_recv[ cnt ];

					cnt++;
				}}}}
			}

			DEALLOC_ARRAY( data_recv, n_pts );

		} // end of 'n_student >0' routine

	} // end of local patch routine



	//----------------------------------------------------------//
	/* Send grid coordinate of procs on local patch to their
	   advisors on global patch. It canges if local patch moves. */

	if( I.i_belong_global_patch == 0 ) { // local patch

	} // end of local patch routine


	if( I.i_belong_global_patch ) {

		PATCH_LOOP {

			if( i_patch != GLOBAL_PATCH_ID ) {

				student = &( distant_charges[ i_patch ] );
				n_student = student->n_procs;

				if( n_student > 0 ) { // if I have student on local patch of PATCH_ID=i_patch

					patch_student = &( patches[ i_patch ] ); // just for convenience

					N1TOT_s = patch_student->N1tot; // "_s" denotes "student (on local patch)"
					N2TOT_s = patch_student->N2tot;
					N3TOT_s = patch_student->N3tot;

					n_pts = patch_student->Ntot * NDIM; // = N1TOT_s * N2TOT_s * N3TOT_s * NDIM
					ALLOC_ARRAY( data_recv, n_pts );

					// allocate memory space to store grid coordinates of students
					ALLOC_5D_ARRAY( student->coord,
							n_student,
							N1TOT_s,
							N2TOT_s,
							N3TOT_s,
							NDIM );

					for( ip=0; ip< n_student; ip++ ) { // receive coordinate student by student

						student_patch_rank = student->pids[ ip ] - patch_student->pids[0];
						tag_recv = myid + numprocs*i_patch + student_patch_rank;

						MPI_Recv( data_recv,
							  n_pts,
							  MPI_DOUBLE,
							  //-------
							  student->pids[ ip ],
							  tag_recv,
							  MPI_COMM_WORLD,
							  MPI_STATUS_IGNORE );

						// unpack the coordinate list
						cnt = 0;
						for( i=0; i< N1TOT_s; i++ ) {
						for( j=0; j< N2TOT_s; j++ ) {
						for( k=0; k< N3TOT_s; k++ ) {
						for( d=0; d< NDIM   ; d++ ) {

							student->coord[ ip ][i][j][k][d] = data_recv[ cnt ];

							cnt++;
						}}}}
					}

					DEALLOC_ARRAY( data_recv, n_pts );

				} // end of 'n_student >0' routine

			} // end of if-not-global_patch

		} // end of PATCH_LOOP

	} // end of global patch routine

}


