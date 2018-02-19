
/* Written by Hotaka Shiokawa, modified on 02/19/2018 */

#include "decs.h"
#include "patchwork.h"
#include "harm_mpi.h"




/*****************************************************************************
  connect_patches_lg():
******************************************************************************/

void connect_patches_lg( double ****prim_arr, int i_patch )
{
	long cnt;
	int i, j, k, l, d, ii, jj, kk, ip, jp, kp, is, js, ks;
	int p_rank, index, cnt1, cnt2, cnt3, n_pts, n_pts_3;
	int advisor_pid, student_pid, student_n_pts, student_n_pts_3;;
	int tag_recv, tag_send, tag_recv_intra, tag_send_intra;
	int patch_rank_of_x, patch_cpupos[NDIM];
	int N1TOT_s, N2TOT_s, N3TOT_s, NTOT_s;

	struct of_charges *student;
	struct of_patch *patch_advisor;
	struct of_patch *patch_student;

	double xcart[NDIM], xp[NDIM];
	double ucon[NDIM], ucon_cart[NDIM];
	double *pl;
	struct of_coord *coords;
	struct of_geom *geom, geom_cart;

	int *assign_n_pts_3, *list_assign_n_pts_3, **assign_index_map;
	double **assign_coord_list, **data_recv;

	double *data, *(data_interp_send[Ncpu]), *(data_interp_recv[Ncpu]);

	MPI_Request request_recv[Ncpu];

	void interp_prims( double *xphys, double *prim_dest, double ****prim_arr );
	void edge_softener( double ****prim_arr, int i_patch );
	void fill_rim( double ****prim_arr );




	//----------------------------------------------------------//
	/* Put flags on grids in global patch that are covered by 
	   local patch */

	if( I.i_belong_global_patch == 0 ) { // local patch


		// students of a local proc always belong to global patch
		student	     = &( distant_charges[  GLOBAL_PATCH_ID ] );
		patch_student = &( patches[ GLOBAL_PATCH_ID ] );

		N1TOT_s = patch_student->N1tot; // '_s' denotes 'student'
		N2TOT_s = patch_student->N2tot;
		N3TOT_s = patch_student->N3tot;
		NTOT_s  = patch_student->Ntot;


		// for advisors who have multiple students, deal with students one by one
		for( i=0; i< student->n_procs_max; i++ ) {


			/* (1) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   			/* Find cover flags (cover_flags_temp[whatever][Global patch's NTOT])
			   and get coordinate of ghost zones on global patch
			   (student->coord_list[student #][cnt]) */

			if( i < student->n_procs ) {

				// Ghost zones are stationary for stationary patches, so,
				// below is needed only very first time or when patch is moving.
				// Dump files want covered grids to be filled, so, below is
				// also required when we create dump files 
				if( first_connect_lg || my_patch->move || all_ghost > 0 ) {

					// Obtain cover_flags_temp[whatever][Global patch's NTOT].
					// Using cover_flags_temp, fill list of ghost zones coordinate
					// (student->coord_list[student #][cnt])
					get_cover_flags_temp( patch_student, student, i, i_patch ); // patchwork_lg.c

				} // end of if( first_connect_lg || my_patch->move )

				student_n_pts   = student->n_pts[i];
				student_n_pts_3 = student_n_pts*3;

			} // end of 'if( i < student->n_procs )'




			/* (2) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   			/* Find procs that can deal with interpolation of each coordinate */




			/* (3) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			/* gather communication info to decider and distribute, so that everyone knows
			   who to communicate how much of data */

 			list_assign_n_pts_3 = intra_patch_comm[ i_patch ].list_assign_n_pts_3[ i ];

			if( first_connect_lg || my_patch->move || all_ghost > 0 ) {

				// MPI_Allgather works, but Gather-Bcast through decider is supposed to be faster (?)
				MPI_Gather( assign_n_pts_3,
					    Ncpu,
					    MPI_INT,
					    //-------
					    list_assign_n_pts_3,
					    Ncpu,
					    MPI_INT,
					    //-------
					    my_patch->patch_decider_rank,
					    patch_mpi_comm );

				MPI_Bcast( list_assign_n_pts_3,
					   Ncpu*Ncpu,
					   MPI_INT,
					   //-------
					   my_patch->patch_decider_rank,
					   patch_mpi_comm );
			}




			/* (4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			/* Send out/receive list of coordinate to each cpu for interpolation */

			data_recv = intra_patch_comm[ i_patch ].data_recv[ i ];

			if( first_connect_lg || my_patch->move || all_ghost > 0 ) {

				for( p_rank=0; p_rank< Ncpu; p_rank++ ) {

					n_pts_3 = list_assign_n_pts_3[ p_rank*Ncpu + I.patch_rank ];

					if( n_pts_3 > 0 ) {

						ALLOC_ARRAY( data_recv[ p_rank ], n_pts_3 );

						tag_recv_intra = p_rank + Ncpu*I.patch_rank + 2*numprocs;

						MPI_Irecv( data_recv[ p_rank ],
							   n_pts_3,
							   MPI_DOUBLE,
							   //-------
							   p_rank,
							   tag_recv_intra,
							   patch_mpi_comm,
							   &( request_recv[ p_rank ]) );
					}

				}

				for( p_rank=0; p_rank< Ncpu; p_rank++ ) {

					n_pts_3 = assign_n_pts_3[ p_rank ]; // # of data points to be sent to
								    // cpu of patch_rank = p_rank
					if( n_pts_3 > 0 ) {

						tag_send_intra = I.patch_rank + Ncpu*p_rank + 2*numprocs;

						MPI_Send( assign_coord_list[ p_rank ],
							  n_pts_3,
							  MPI_DOUBLE,
							  //-------
							  p_rank,
							  tag_send_intra,
							  patch_mpi_comm );
					}
				}



				for( p_rank=0; p_rank< Ncpu; p_rank++ ) {
					if( list_assign_n_pts_3[ p_rank*Ncpu + I.patch_rank ] ) {
						MPI_Wait( &(request_recv[ p_rank ]), MPI_STATUS_IGNORE );
					}
				}

			}




			/* (5) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			/* interpolation */



			/* (6) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			/* send/receive the interpolated data */

			for( p_rank=0; p_rank< Ncpu; p_rank++ ) {

				n_pts = assign_n_pts_3[ p_rank ]/3;

				if( n_pts > 0 ) {

					ALLOC_ARRAY( data_interp_recv[ p_rank ], n_pts*NP );

					tag_recv_intra = p_rank + Ncpu*I.patch_rank + numprocs*numprocs;

					MPI_Irecv( data_interp_recv[ p_rank ],
						   n_pts*NP,
						   MPI_DOUBLE,
						   //-------
						   p_rank,
						   tag_recv_intra,
						   patch_mpi_comm,
						   &( request_recv[ p_rank ] ) );
				}

			}

			for( p_rank=0; p_rank< Ncpu; p_rank++ ) {

				n_pts = list_assign_n_pts_3[ p_rank*Ncpu + I.patch_rank ]/3;
								   
				if( n_pts > 0 ) {

					tag_send_intra = I.patch_rank + Ncpu*p_rank + numprocs*numprocs;

					MPI_Send( data_interp_send[ p_rank ],
						  n_pts*NP,
						  MPI_DOUBLE,
						  //-------
						  p_rank,
						  tag_send_intra,
						  patch_mpi_comm );
				}
			}


			for( p_rank=0; p_rank< Ncpu; p_rank++ ) {
				if( assign_n_pts_3[ p_rank ] ) {
					MPI_Wait( &(request_recv[ p_rank ]), MPI_STATUS_IGNORE );
				}
			}




			/* (7) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			/* reconstruct data structure */



			/* (8) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			/* send interpolated data as well as patch motion info to my "student(s)" */



			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			/* deallocate arrays */

			if( i < student->n_procs ) {

				DEALLOC_ARRAY( data, student_n_pts*NP );

				if( my_patch->move || all_ghost == ALL_GHOST ) {
					for( p_rank=0; p_rank< Ncpu; p_rank++ ) {

						DEALLOC_ARRAY( assign_coord_list[ p_rank ], student_n_pts_3 );
						DEALLOC_ARRAY( assign_index_map[  p_rank ], student_n_pts   );
					}
				}
			}

			for( p_rank=0; p_rank< Ncpu; p_rank++ ) {
				n_pts_3 = assign_n_pts_3[ p_rank ];
				n_pts = n_pts_3/3;
				if( n_pts_3 > 0 ) {
					DEALLOC_ARRAY( data_interp_recv[ p_rank ], n_pts*NP );
				}

				n_pts_3 = list_assign_n_pts_3[ p_rank*Ncpu + I.patch_rank ];
				n_pts = n_pts_3/3;
				if( n_pts_3 > 0 ) {
					if( my_patch->move || all_ghost == ALL_GHOST ) {
						DEALLOC_ARRAY( data_recv[ p_rank ], n_pts_3 );
					}
					DEALLOC_ARRAY( data_interp_send[ p_rank ], n_pts*NP );
				}
			}


		} // end of student loop

	} // end of 'if I belong to local patch'





	if( I.i_belong_global_patch ) {

		// Hear from advisor
		patch_advisor = &( patches[i_patch] );
		advisor_pid = local_charges_pid[ i_patch ];

		tag_recv = (my_patch->n_procs + 3) * i_patch + myid;

		// Need only first time for stationary local patches.
		if( first_connect_lg || patch_advisor->move || all_ghost > 0 ) {

			// First, get number of data points to be communicated
			MPI_Recv( &(data_ex[ i_patch ].n_pts),
				  1,
				  MPI_INT,
				  //-------
				  advisor_pid,
				  tag_recv,
				  MPI_COMM_WORLD,
				  MPI_STATUS_IGNORE );

			MPI_Recv( cover_flags_temp[ i_patch ],
				  my_patch->Ntot,
				  MPI_INT,
				  //-------
				  advisor_pid,
				  tag_recv+2,
				  MPI_COMM_WORLD,
				  MPI_STATUS_IGNORE );
		}



		//----------------------------------------------------------//
		/* Receive Interpolated data, reconstruct data structure,
		   and transform vectors */

		// Allocate space for incoming interpolated data.
		// Even if local patch is stationary, always allocate and deallocate
		// this array since n_pts change if all_ghosh==ALL_GHOST
		ALLOC_ARRAY( data_ex[ i_patch ].data, data_ex[ i_patch ].n_pts*NP );

		MPI_Recv( data_ex[ i_patch ].data,    // interpolated data
			  data_ex[ i_patch ].n_pts*NP,
			  MPI_DOUBLE,
			  //-------
			  advisor_pid,
			  tag_recv+1,
			  MPI_COMM_WORLD,
			  MPI_STATUS_IGNORE );

		cnt1 = 0;
		cnt2 = 0;
		ALL_LOOP {

			if( (  cover_flags_temp[ i_patch ][ cnt1 ] > 0				    ) ||
			    // Note that cover_flags_temp[][] is NOT modified even if all_ghost==ALL_GHOST,
			    // so, grids with covered-but-not-ghost flag (flag<0) need to be filled, too.
			    ( (cover_flags_temp[ i_patch ][ cnt1 ] < 0) && (all_ghost == ALL_GHOST) )    ) {

				PLOOP {
					prim_arr[i][j][k][l] = data_ex[ i_patch ].data[ cnt2++ ];
				}

				/* Transform vector components */
				// calculate ucon in global-cartesian coordinate
				get_coord( i, j, k, CENT, ncurr, coords );	// decs.h
				get_geom_cart( coords->x,
					       coords->xp,
					       coords->xcart,
					       coords->xspher, &geom_cart );	// metric.c

				ucon_calc( prim_arr[i][j][k], &geom_cart, ucon_cart ); // phys.c

				// transform ucon to numerical coordinate
				for( ii=0; ii< NDIM; ii++ ) { ucon[ii] = 0.; }
				for( ii=0; ii< NDIM; ii++ ) {
					for( jj=0; jj< NDIM; jj++ ) {

						ucon[ii] += ucon_cart[jj]*coords->dxp_dxcart[ii][jj];
				}}

				// calculate prim in numerical coordinate
				get_geometry( i, j, k, CENT, ncurr, geom ); 
				ucon2pr( prim_arr[i][j][k], ucon, geom->gcon ); // phys.c
			}


			// Unused covered grids have their prim. values set to be -1
			// unless all_ghost==ALL_GHOST
			if( (cover_flags_temp[ i_patch ][ cnt1 ] < 0) && (all_ghost != ALL_GHOST) ) {
				PLOOP {
					prim_arr[i][j][k][l] = -1;
				}
			}

			cnt1++;

		} // end of ALL_LOOP

		DEALLOC_ARRAY( data_ex[ i_patch ].data, data_ex[ i_patch ].n_pts*NP );


		/* Average edge ghost zone values (grids with cover_flags=2) with
		   their neighbors to avoid drastic change that may happen when
		   local grids are much smaller than global grids. This operation
		   is effectively similar to larger area interpolation from local to
		   global. Just comment it our if you don't need. */
		edge_softener( prim_arr, i_patch );


	} // end of if( I.i_belong_global_patch)


}





void get_cover_flags_temp( struct of_patch *patch_student,
			   struct of_charges *student,
			   int student_tag,
			   int i_patch )
{
	static double *xcart_g, xp_ls[NDIM], xp_le[NDIM];
	double xp_l[NDIM];
	long cnt, cnt1, cnts;
	int i, g, is, js, ks;
	int N1TOT_s, N2TOT_s, N3TOT_s, NTOT_s;


	if( first_connect_lg ) {

		ALLOC_ARRAY( xcart_g, NDIM );

		for( i=1; i< NDIM; i++ ) {
			xp_ls[i] = startx[i] + NG*dx[i];	  // xp_l(ocal)s(tart)
			xp_le[i] = xp_ls[i] + totalsize[i]*dx[i]; // xp_l(ocal)e(nd)
		}

		// avoid any usage of ghost zones in interpolation
		for( i=1; i< NDIM; i++ ) {
			xp_ls[i] += 0.5*dx[i];
			xp_le[i] -= 0.5*dx[i];
		}
		// but not for periodic boundaries
		if( TOP_TYPE_CHOICE == TOP_SPHERICAL || TOP_TYPE_CHOICE == TOP_CYLINDRICAL ) {
			xp_ls[3] -= 0.5*dx[3];
			xp_le[3] += 0.5*dx[3];
		}
	}

	N1TOT_s = patch_student->N1tot; // '_s' denotes 'student'
	N2TOT_s = patch_student->N2tot;
	N3TOT_s = patch_student->N3tot;
	NTOT_s  = patch_student->Ntot;


	/* mark global zones covered by local patch */
	cnt = 0;
	for( is=0; is< N1TOT_s; is++) {
	for( js=0; js< N2TOT_s; js++) {
	for( ks=0; ks< N3TOT_s; ks++) {

		xcart_g = student->coord[ student_tag ][is][js][ks]; // x_g(lobal)
		xp_of_xcart( xp_l, xcart_g ); // xp_l(ocal) corresponding to xcart_g(lobal), coord.c

		cover_flags_temp[ i_patch ][ cnt ] = 0; // uncovered grids are flagged by 0
		cover_flags_temp[ i_patch ][ cnt ] = 0;

		if( (xp_l[1] > xp_ls[1]) && (xp_l[1] < xp_le[1]) ) {
		if( (xp_l[2] > xp_ls[2]) && (xp_l[2] < xp_le[2]) ) {
		if( (xp_l[3] > xp_ls[3]) && (xp_l[3] < xp_le[3]) ) {
			cover_flags_temp[ i_patch ][ cnt ] = -1; // covered grids are flagged by -1
			cover_flags_temp[ i_patch ][ cnt ] = -1;
		}}}

		cnt++;
	}}}


	/* mark ghost zones */
	cnt = 0;
	cnt1 = 0;
	student->n_pts[ student_tag ] = 0;
	for( is=0; is< N1TOT_s; is++ ) { // whole range required
	for( js=0; js< N2TOT_s; js++ ) {
	for( ks=0; ks< N3TOT_s; ks++ ) {

		cnt = get_cnt( is, js, ks, N1TOT_s, N2TOT_s, N3TOT_s ); // misc.c

		if( cover_flags_temp[ i_patch ][ cnt ] < 0 ) { // if this grid is covered

			// check if this grid is ghost
			for( g=1; g<= NG; g++ ) {

				if( is-g >= 0 ) {
					cnts = get_cnt( is-g, js, ks, N1TOT_s, N2TOT_s, N3TOT_s ); // check negative direc
					if( cover_flags_temp[ i_patch ][ cnts ] == 0 ) { // - if grid(i-g, j, k) is not covered
						if( g == 1 ) {
							cover_flags_temp[ i_patch ][ cnt ] = 2; // - mark grids at the edge of
												//   covered region by 2
							break;
						} else {
							cover_flags_temp[ i_patch ][ cnt ] = 1; // - +1 denotes ghost zone
						}
					}
				}

				if( is+g < N1TOT_s ) {
					cnts = get_cnt( is+g, js, ks, N1TOT_s, N2TOT_s, N3TOT_s ); // check positive direc
					if( cover_flags_temp[ i_patch ][ cnts ] == 0 ) {
						if( g == 1 ) {
							cover_flags_temp[ i_patch ][ cnt ] = 2;
							break;
						} else {
							cover_flags_temp[ i_patch ][ cnt ] = 1;
						}
					}
				}

				if( js-g >= 0 ) {
					cnts = get_cnt( is, js-g, ks, N1TOT_s, N2TOT_s, N3TOT_s );
					if( cover_flags_temp[ i_patch ][ cnts ] == 0 ) {
						if( g == 1 ) {
							cover_flags_temp[ i_patch ][ cnt ] = 2;
							break;
						} else {
							cover_flags_temp[ i_patch ][ cnt ] = 1;
						}
					}
				}

				if( js+g < N2TOT_s ) {
					cnts = get_cnt( is, js+g, ks, N1TOT_s, N2TOT_s, N3TOT_s );
					if( cover_flags_temp[ i_patch ][ cnts ] == 0 ) {
						if( g == 1 ) {
							cover_flags_temp[ i_patch ][ cnt ] = 2;
							break;
						} else {
							cover_flags_temp[ i_patch ][ cnt ] = 1;
						}
					}
				}

				if( ks-g >= 0 ) {
					cnts = get_cnt( is, js, ks-g, N1TOT_s, N2TOT_s, N3TOT_s );
					if( cover_flags_temp[ i_patch ][ cnts ] == 0 ) {
						if( g == 1 ) {
							cover_flags_temp[ i_patch ][ cnt ] = 2;
							break;
						} else {
							cover_flags_temp[ i_patch ][ cnt ] = 1;
						}
					}
				}

				if( ks+g < N3TOT_s ) {
					cnts = get_cnt( is, js, ks+g, N1TOT_s, N2TOT_s, N3TOT_s );
					if( cover_flags_temp[ i_patch ][ cnts ] == 0 ) {
						if( g == 1 ) {
							cover_flags_temp[ i_patch ][ cnt ] = 2;
							break;
						} else {
							cover_flags_temp[ i_patch ][ cnt ] = 1;
						}
					}
				}

			} // end of "for( g=1; g<= NG; g++ )"



			if( cover_flags_temp[ i_patch ][ cnt ] < 0 ) { // if this grid is not set to ghost yet

				// Even if a covered grid is not used as ghost zone in numerical integration,
				// it still may be accessed by local grids in the global->local interpolation routine.
				// Below checks "closeness" to local patch edge by checking grids in diagonal direction.
				// Complete checking requires whole range of -NG<=i<=NG, -NG<=j<=NG, and -NG<=k<=NG,
				// but diagonal direction is usually sufficient and warned in interp.c if it fails.

				for( g=1; g<= NG; g++ ) {

					if( is-g >= 0 && js-g >= 0 && ks-g >= 0 ) {
						cnts = get_cnt( is-g, js-g, ks-g, N1TOT_s, N2TOT_s, N3TOT_s ); 
						if( cover_flags_temp[ i_patch ][ cnts ] == 0 ) {
							cover_flags_temp[ i_patch ][ cnt ] = 3; // - +3 denotes zones that are covered,
												// not ghost, but may be accessed in interpolation
							break;
						}
					}

					if( is-g >=0 && js-g >=0 && ks+g <= N3TOT_s ) {
						cnts = get_cnt( is-g, js-g, ks+g, N1TOT_s, N2TOT_s, N3TOT_s );
						if( cover_flags_temp[ i_patch ][ cnts ] == 0 ) {
							cover_flags_temp[ i_patch ][ cnt ] = 3;
							break;
						}
					}

					if( is-g >= 0 && js+g <= N2TOT_s && ks-g >= 0 ) {
						cnts = get_cnt( is-g, js+g, ks-g, N1TOT_s, N2TOT_s, N3TOT_s );
						if( cover_flags_temp[ i_patch ][ cnts ] == 0 ) {
							cover_flags_temp[ i_patch ][ cnt ] = 3;
							break;
						}
					}

					if( is+g <= N1TOT_s && js-g >= 0 && ks-g >= 0 ) {
						cnts = get_cnt( is+g, js-g, ks-g, N1TOT_s, N2TOT_s, N3TOT_s );
						if( cover_flags_temp[ i_patch ][ cnts ] == 0 ) {
							cover_flags_temp[ i_patch ][ cnt ] = 3;
							break;
						}
					}

					if( is-g >=0 && js+g <= N2TOT_s && ks+g <= N3TOT_s ) {
						cnts = get_cnt( is-g, js+g, ks+g, N1TOT_s, N2TOT_s, N3TOT_s );
						if( cover_flags_temp[ i_patch ][ cnts ] == 0 ) {
							cover_flags_temp[ i_patch ][ cnt ] = 3;
							break;
						}
					}

					if( is+g <= N1TOT_s && js-g >= 0 && ks+g <= N3TOT_s ) {
						cnts = get_cnt( is+g, js-g, ks+g, N1TOT_s, N2TOT_s, N3TOT_s );
						if( cover_flags_temp[ i_patch ][ cnts ] == 0 ) {
							cover_flags_temp[ i_patch ][ cnt ] = 3;
							break;
						}
					}

					if( is+g <= N1TOT_s && js+g <= N2TOT_s && ks-g >= 0 ) {
						cnts = get_cnt( is+g, js+g, ks-g, N1TOT_s, N2TOT_s, N3TOT_s );
						if( cover_flags_temp[ i_patch ][ cnts ] == 0 ) {
							cover_flags_temp[ i_patch ][ cnt ] = 3;
							break;
						}
					}

					if( is+g <= N1TOT_s && js+g <= N2TOT_s && ks+g <=N3TOT_s ) {
						cnts = get_cnt( is+g, js+g, ks+g, N1TOT_s, N2TOT_s, N3TOT_s );
						if( cover_flags_temp[ i_patch ][ cnts ] == 0 ) {
							cover_flags_temp[ i_patch ][ cnt ] = 3;
							break;
						}
					}

				} // end of "for( g=1; g<= NG; g++ )"

			} // end of inner 'if cover_flags_temp < 0' (meaning if the grid is covered)


			if( cover_flags_temp[ i_patch ][ cnt ] > 0 || all_ghost == ALL_GHOST ) { // if ghost

				(student->n_pts[ student_tag ])++;

				// record coordinate of ghost zones
				student->coord_list[ student_tag ][ cnt1++ ] = student->coord[ student_tag ][is][js][ks][1];
				student->coord_list[ student_tag ][ cnt1++ ] = student->coord[ student_tag ][is][js][ks][2];
				student->coord_list[ student_tag ][ cnt1++ ] = student->coord[ student_tag ][is][js][ks][3];
			}

		} // end of outer 'if cover_flags_temp < 0' (meaning if the grid is covered)

	}}} // end of ijk loop for ghost zone marking

}





/*****************************************************************************
	Average edge ghost zone values (grids with cover_flags=2) with
	their neighbors to avoid drastic change that may happen when
	local grids are much smaller than global grids. This operation
	is, to some level, equivalet to interpolation using
	global-grid-size-scale separated points on local patch, although
	at least one of the global grids in this averaging process is
	physical grid i.e. not derived from local patch information.
*****************************************************************************/

void edge_softener( double ****prim_arr, int i_patch )
{
	long cnt;
	int i, j, k, l;


	cnt = 0;
	ALL_LOOP {

		if( cover_flags_temp[ i_patch ][ cnt ] == 2 ) {

			if( i>0       && j>0       && k>0       &&
			    i<N1TOT-1 && j<N2TOT-1 && k<N3TOT-1    ) {

				if(  prim_arr[i-1][j][k][RHO]*prim_arr[i+1][j][k][RHO]
				    *prim_arr[i][j-1][k][RHO]*prim_arr[i][j+1][k][RHO]
				    *prim_arr[i][j][k-1][RHO]*prim_arr[i][j][k+1][RHO] < 0. ) {

					fprintf(stderr, "WARNING: covered grid is accessed in patchwork_lg.c:\n");
					fprintf(stderr, "myid=%d, ijk=(%d,%d,%d)\n", myid, i, j, k);

					for( l=RHO; l<=UU; l++ ) {
						prim_arr_temp[i][j][k][l] = prim_arr[i][j][k][l];
					}
					
				} else {

					for( l=RHO; l<=UU; l++ ) {
						prim_arr_temp[i][j][k][l] =   6.*prim_arr[i  ][j  ][k  ][l]
									    +    prim_arr[i-1][j  ][k  ][l]
									    +    prim_arr[i+1][j  ][k  ][l]
									    +    prim_arr[i  ][j-1][k  ][l]
									    +    prim_arr[i  ][j+1][k  ][l]
									    +    prim_arr[i  ][j  ][k-1][l]
									    +    prim_arr[i  ][j  ][k+1][l];

						prim_arr_temp[i][j][k][l] /= 12.;
					}
				}
			}
		}

		cnt++;
	}

	cnt = 0;
	ALL_LOOP {

		if( cover_flags_temp[ i_patch ][ cnt ] == 2 ) {

			if( i>0       && j>0       && k>0       &&
			    i<N1TOT-1 && j<N2TOT-1 && k<N3TOT-1    ) {
				//PLOOP {
				for( l=RHO; l<=UU; l++ ) {
					prim_arr[i][j][k][l] = prim_arr_temp[i][j][k][l];
				}
			}
		}

		cnt++;
	}

}
