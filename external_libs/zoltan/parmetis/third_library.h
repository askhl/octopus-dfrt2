/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile: third_library.h,v $
 *    $Author: kddevin $
 *    $Date: 2008/09/15 23:00:11 $
 *    Revision: 1.3.6.1 $
 ****************************************************************************/


#ifndef __THIRD_LIBRARY_H
#define __THIRD_LIBRARY_H

#include <limits.h>
#include "zoltan_comm.h"
#include "third_library_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#define INT_EPSILON (1e-5) /* How close we need to be to say it's an integer */

/* Macro for error handling */
#define ZOLTAN_PARMETIS_ERROR(error,str) {ierr = error ; \
 ZOLTAN_PRINT_ERROR(zz->Proc, __func__, str) ; goto End ;}

#define ZOLTAN_THIRD_ERROR(error,str) { \
 ZOLTAN_PRINT_ERROR(zz->Proc, __func__ , str) ; return(error) ;}


/* Structure that defines a graph for third party libraries like ParMetis. */
typedef struct ZOLTAN_Third_Graph_ {
  int graph_type;                       /* Type of the graph */
  int check_graph;                      /* We have to check graph consistency */
  int final_output;                     /* Do some output computations */
  int showMoveVol;                      /* How works the final output */
  int id_known;                         /* Associated gids & lids are already known */
  int scatter;                          /* Graph has been scattered */
  int scatter_min;                      /* Minimum level of scatter */
  int get_data;                         /* Construct edge datas */
  int obj_wgt_dim;                      /* Number of weights by vertex */
  int edge_wgt_dim;                     /* Number of weights by edge */
  int num_obj;                          /* Number of vertices */
  int num_obj_orig;                     /* Number of vertices in original graph */
  int num_edges;                        /* Number of edges */
  indextype * vtxdist;                  /* How vertices are distributed */
  indextype * xadj;                     /* Indexes on adjency array */
  indextype * adjncy;                   /* adjency array (CSR) */
  weighttype * vwgt;                    /* Array of vertex weights */
  weighttype * ewgts;                   /* Array of edge weights */
  float * float_ewgts;
  int * adjproc;                        /* Array of procs ? */
  ZOLTAN_COMM_OBJ *comm_plan;           /* Communication plan used by scattering process */
} ZOLTAN_Third_Graph;

/* Structure that defines a geometry for third party libraries like ParMetis. */
typedef struct ZOLTAN_Third_Geom_ {
  int ndims;                            /* Number of dimensions */
  float *xyz;                           /* Coordinates */
} ZOLTAN_Third_Geom;

/* Structure that defines a partition for third party libraries like ParMetis. */
typedef struct ZOLTAN_Third_Part_ {
  indextype *part;                      /* Current partition */
  indextype *input_part;
  indextype *part_orig;
  float     *part_sizes;
  float     *input_part_sizes;
} ZOLTAN_Third_Part;

typedef struct ZOLTAN_Third_Vsize_ {
  int vsize_malloc;
  indextype *vsize;
  indextype *vsizeBACKUP;
} ZOLTAN_Third_Vsize;

/* Structure that defines an ordering output for third party libraries like ParMetis. */
typedef struct ZOLTAN_Output_Order_ {
  int num_part;
  int start_index;
  int *rank;            /* rank[i] is the rank of gids[i] */
  int *iperm;           /* inverse permutation of rank */
  ZOOS *order_opt;	/* ordering options */
  ZOS *order_info;	/* ordering info */
  indextype *sep_sizes;
} ZOLTAN_Output_Order;

/* Structure that describes a partitioning output in Zoltan lingo
   for third party libraries like ParMetis. */
typedef struct ZOLTAN_Output_Part_ {
  int  compute_only_part_changes;

  int  num_imp;
  ZOLTAN_ID_PTR* imp_gids;
  ZOLTAN_ID_PTR* imp_lids;
  int **imp_procs;
  int **imp_part;

  int  num_exp;
  ZOLTAN_ID_PTR *exp_gids;
  ZOLTAN_ID_PTR *exp_lids;
  int **exp_procs;
  int **exp_part;
} ZOLTAN_Output_Part;

/*****************************
 *  Functions which can be used to interface with third party libraries
 ****************************************************************************/

/* Construct graph and associated datas in order to call third library */
int Zoltan_Preprocess_Graph(
  ZZ *zz,                               /* Zoltan structure */
  ZOLTAN_ID_PTR *global_ids,
  ZOLTAN_ID_PTR *local_ids,
  ZOLTAN_Third_Graph *gr,              /* Graph for third part libs */
  ZOLTAN_Third_Geom  *geo,
  ZOLTAN_Third_Part  *prt,
  ZOLTAN_Third_Vsize *vsp
  );

/* Activate timers if requested */
int Zoltan_Preprocess_Timer(ZZ *zz, int *use_timer);

/* Display timing informations */
void Zoltan_Third_DisplayTime(ZZ* zz, double* times);

/* Free temporary datas */
void Zoltan_Third_Exit(ZOLTAN_Third_Graph *gr, ZOLTAN_Third_Geom *geo,
		       ZOLTAN_Third_Part *prt, ZOLTAN_Third_Vsize *vsp,
		       ZOLTAN_Output_Part *part, ZOLTAN_Output_Order *ord);

/* Do some postprocesing in order to exploit third party results in Zoltan */
int Zoltan_Postprocess_Graph(
  ZZ *zz,                               /* Zoltan structure */
  ZOLTAN_ID_PTR      global_ids,
  ZOLTAN_ID_PTR      local_ids,
  ZOLTAN_Third_Graph *gr,              /* Graph for third part libs */
  ZOLTAN_Third_Geom  *geo,
  ZOLTAN_Third_Part  *prt,
  ZOLTAN_Third_Vsize *vsp,
  ZOLTAN_Output_Order *ord,
  ZOLTAN_Output_Part  *part);

int
Zoltan_Postprocess_FinalOutput (ZZ* zz, ZOLTAN_Third_Graph *gr,
				ZOLTAN_Third_Part *prt, ZOLTAN_Third_Vsize *vsp,
				int use_timers, double itr);

#ifdef __cplusplus
}
#endif

#endif
