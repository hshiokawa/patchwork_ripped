
/* Written by Hotaka Shiokawa, modified on 02/19/2018 */

/*****************************************************************************
     P A T C H W O R K    H E A D E R     F I L E  
******************************************************************************/

#include "harm_mpi.h"


/********************************************************************************************
  PATCHWORK MACROS : 
*********************************************************************************************/
#define  SEAM_INTRAPATCH   (0)
#define  SEAM_STATIC       (1)
#define  SEAM_DYANMIC      (2)

#define ALL_GHOST (1)
#define FIX_GHOST (2)

#define GLOBAL_TO_LOCAL	(0)
#define LOCAL_TO_GLOBAL	(1)



/********************************************************************************************
  PATCHWORK STRUCTURES : 
*********************************************************************************************/

/*******************************************************************************************************
 Data structure used for communicating boundary data:  Different interpretations for different types:
     seam_type = SEAM_INTRAPATCH:  -- use bbox instead of coord_list;
     seam_type = SEAM_STATIC:      -- use coord_list and only communicate coord_list once; 
     seam_type = SEAM_DYNAMIC:     -- use coord_list and communicate coord_list always; 
******************************************************************************************************/



struct of_data_exchange {
	double   *coord_list;    /* list of coordinates at which data is communicated */
	double   *data;          /* data[n_pts*NP] interpolated data to be communicated */
	int      *index_map;     /* index_map[n_pts*3] are the local i,j,k indices of the cell within which we interpolate */
	int       n_pts;         /* number of points in coord_list[] and in data[] */

} ;



/*******************************************************************************************************
 Data structure used describing a patch to another patch: 
******************************************************************************************************/
struct of_patch {
	int	id;			/* unique id number for patch */
	int	n_procs;		/* number of processors within the patch  */
	int	N1tot;			/* N1TOT */
	int	N2tot;			/* N2TOT */
	int	N3tot;			/* N3TOT */
	long	Ntot;			/* N1TOT * N2TOT * N3TOT */
	int	*pids;			/* WORLD-level pids of all members of the patch */
	int	name_length;		/* length of tag name */ 
	char	*name;			/* name of the patch (i.e. RUN_TAG); */ 
	int	patch_decider_rank;	/* local group id of rank of decider */
} ;


/* Intra-patch communication related variables in connect_patches_*() */
struct of_intra_patch_comm {
	int **assign_n_pts_3;
	int **list_assign_n_pts_3;
	int ***assign_index_map;
	double ***assign_coord_list;
	double ***data_recv;
} ;

extern struct of_intra_patch_comm *intra_patch_comm; 



/* */
struct of_charges { // info of procs on other patches that I will take care / will be taken care
	int	n_procs;	/* number of distant-procs (procs in other patches) */
	int	n_procs_max;	/* maximum possible number of distant-procs */
	double	*****coord;	/* coordinate of each grid of distant-proc's domain */
	double	**coord_list;	/* coordinate of each ghost grid of distant-proc's domain */
} ;

extern struct of_charges *distant_charges;
extern int *local_charges_pid;



/********************************************************************************************
  PATCHWORK GLOBALS: 
*********************************************************************************************/
extern struct of_patch		*patches;	  /* array of patches */ 
extern struct of_patch		*my_patch;	  /* pointer to the patch I belong to */
extern struct of_data_exchange	*data_ex;	  /* array of data exchange strucures of all patches */

extern double ****prim_arr_temp;   /* temporal primitives on global patch for whatever useage */

extern int first_connect_gl;    /* identify if this is the first call of connect_patches_gl */
extern int first_connect_lg;    /* identify if this is the first call of connect_patches_lg */
extern int all_ghost;		/* boolean to make all covered grids on global patch be ghost zones */

extern MPI_Group  world_group;		/* group of all processors */
extern MPI_Comm   patch_mpi_comm;	/* MPI comm handle for local patch, used for intrapatch communications */
extern MPI_Group  decider_mpi_group;	/* MPI group handle for all the deciders  */


//===================================================================
/* patchwork_setup.c */
extern void identify_patches( const char *local_patch_name );
extern void setup_data_exchange( void );
extern void record_distant_coord( void );

/* patchwork_gl.c */
extern void get_coord_list( void );


