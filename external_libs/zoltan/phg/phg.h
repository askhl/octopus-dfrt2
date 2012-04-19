/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile: phg.h,v $
 *    $Author: kddevin $
 *    $Date: 2008/09/17 20:29:21 $
 *    Revision: 1.77.6.1 $
 ****************************************************************************/

#ifndef __ZOLTAN_PHG_H
#define __ZOLTAN_PHG_H

#include "params_const.h"
#include "zoltan_comm.h"
#include "zz_rand.h"
#include "phg_const.h"
#include "phg_util.h"
#include "phg_hypergraph.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif




/* Function types for options to hypergraph partitioning */
struct PHGPartParamsStruct;  /* Forward declaration */

typedef int ZOLTAN_PHG_MATCHING_FN(ZZ*, HGraph*, Matching, 
                                   struct PHGPartParamsStruct*);
typedef int ZOLTAN_PHG_COARSEPARTITION_FN(ZZ*, HGraph*, int, float *, Partition,
                                          struct PHGPartParamsStruct*);
typedef int ZOLTAN_PHG_REFINEMENT_FN(ZZ*, HGraph*, int, float *, Partition,
                                     struct PHGPartParamsStruct*, float);


    
/******************************************/
/* Parameters to the hypergraph functions */
/******************************************/
struct PHGPartParamsStruct {
  char hgraph_pkg[MAX_PARAM_STRING_LEN];/* Package (Zoltan,PaToH,ParKway) */
  char hgraph_method[MAX_PARAM_STRING_LEN]; /* Operation(partition,repart,refine)*/
  int useMultilevel;             /* Flag indicating whether to use multilevel
                                    method */
  float bal_tol;                 /* Balance tolerance in % of average */
  float bal_tol_adjustment;      /* balance tolerance adjustment;
                                    between [0-1.0]: it is % of bal_tol
                                    > 1.0: it is a multiplier to bal_tol/log(k) */
  float *part_sizes;             /* Pointer to part_sizes array passed to
                                    Zoltan_PHG; lists part_sizes (%ages) 
                                    for all requested partitions. */
  int kway;                      /* 1 -> direct kway, 0->recursive bisection */
  int redl;                      /* Reduction limit (constant). */
  int nCand;                     /* Maximum number of candidates per round;
                                    currently only used in agglomerative matching */
  int UseFixedVtx;               /* Flag indicating whether any vertices of
                                    hypergraph are fixed. */
  int UsePrefPart;               /* Flag indicating that the coarsening
                                    should be restricted using pref_part and
                                    coarse partitioning should take the part
                                    preferences into account. */
                                    

  char convert_str[MAX_PARAM_STRING_LEN];  /* Graph->HG conversion method. */
  char redm_str[MAX_PARAM_STRING_LEN];  /* Reduction method string. */
  char redm_fast[MAX_PARAM_STRING_LEN]; /* Fast reduction method string. */
  char redmo_str[MAX_PARAM_STRING_LEN]; /* Matching optimization string*/
    
  ZOLTAN_PHG_MATCHING_FN *matching;    /* Pointers to Matching function */
  int  match_array_type;         /* interpretation of match array:
                                    assuming (u,v, w) are matched
                                    0 -> old style link-list
                                         match[u] = v; match[v] = w;  match[w] = u;
                                    1 -> all vertices point the representative
                                         match[u] = match[v] = match[w] = u;
                                 */
    
  int edge_scaling;              /* type of hyperedge weight scaling */
  int vtx_scaling;               /* type of vertex scaling for inner product */
  int vtx_scal_size;                   /* size of vtx_scal array */
  float *vtx_scal;                     /* vtx scaling array */
  int LocalCoarsePartition;            /* 1 -> apply coarse partitioner locally;
                                          0 -> gather entire HG to each proc
                                          and apply coarse partitioner. */
  char coarsepartition_str[MAX_PARAM_STRING_LEN]; 
                                         /* Coarse partitioning string */
  ZOLTAN_PHG_COARSEPARTITION_FN *CoarsePartition;
                                         /* pointer to coarse partitioning fn */
  char refinement_str[MAX_PARAM_STRING_LEN]; /* Refinement string and */
  ZOLTAN_PHG_REFINEMENT_FN *Refinement;      /* pointer to refinement fn */

  int fm_loop_limit;    /* Number of FM loops (if the refinement is FM) */
  int fm_max_neg_move;  /* Maximum number of vertex moves with negative gain;
                           an early stop condition. <0 means try all moves */
  float refinement_quality; /* Adjustment to refinement parameter settings; 
                               1 means the defaults, higher means more ref. */
    
  int check_graph;      /* Flag indicating whether the input hypergraph should 
                         * be checked for errors. */
  int output_level;     /* Flag indicating amount of output from HG algorithms.
                         * See levels PHG_DEBUG_* below.  */
  int final_output;     /* Prints final timing and quality info at end of PHG
                           (regardless of value of output_level) */

    /* NOTE THAT this comm refers to "GLOBAL" comm structure
       (hence the name change: just to make sure it has not been used
        accidentally in an inccorrect way)
       comm of hg should be used in coarsening/initpart/refinement codes
       because of possible processor splits in recursive bisection
    */
  PHGComm globalcomm;
  int nProc_x_req;  /* user's request for nProc_x (PHG_NPROC_X)
                       -1 if not specificed/auto */
  int nProc_y_req;  /* user's request for nProc_y (PHG_NPROC_Y)
                       -1 if not specificed/auto */

  int num_coarse_iter;  /* Number of coarse partitions to try on each proc. */
  int visit_order;      /* Vertex visit order. */
  int use_timers;       /* Flag indicating whether to time the PHG code. */
  float EdgeSizeThreshold;  /* % of global vtxs beyond which an edge is 
                               considered to be dense. */
  int MatchEdgeSizeThreshold;  /* Edges with sizes bigger than this threshold
                               considered to be dense for matching and ignored. */
  float hybrid_keep_factor; /* h-ipm only: keep matches with i.p. values
                               greater than this factor times the mean */
  char parkway_serpart[MAX_PARAM_STRING_LEN];  /* SerialPartitioner for parKway. */
  int connectivity_cut;     /* 1: if cut_objective==connectivity,
                               0: if it is hyperedges */
  int add_obj_weight;       /* Calculated weight: unit vertex, non-zeroes,
                                                  or none */
  int edge_weight_op;   /* What to do when more than one process returns a
                        weight for the same edge: add, take max, flag error */
  int RandomizeInitDist;  /* Flag indicating whether to randomly distribute
                             vertices and edges passed as input to PHG. */
  float ProRedL;           /* V-cycle processor reduction limit in % of pins */
  float RepartMultiplier;  /* In PHG_REPART, multiply input edge weights by
                              RepartMultiplier, a parameter that should be 
                              proportional to the number of times the 
                              communication described by the hypergraph 
                              will be done before the next repartitioning.  
                              This cost allows us to consider migration cost
                              compared to application communication cost. */
  int patoh_alloc_pool0,   /* to adjust patoh's memory pre-allocation amount */
      patoh_alloc_pool1;    
};

typedef struct PHGPartParamsStruct PHGPartParams;

    
/*****************************/
/* Hypergraph output levels: */
/*****************************/
#define PHG_DEBUG_NONE 0
#define PHG_DEBUG_LIST 1
#define PHG_DEBUG_ALL  2
#define PHG_DEBUG_PRINT 3
#define PHG_DEBUG_PLOT 4


/*********************************************************/
/* Structure for storing timer indices from call to call */
/*********************************************************/
struct phg_timer_indices {
  int all;  /* All of Zoltan_PHG; this timer includes other timers and their
               synchronization time, so it will be high. */
  int build;  /* hypergraph build time */
  int setupvmap;
  int parkway;  /* ParKway time */
  int patoh;  /* PaToH time */
  int retlist;  /* Building return lists time */
  int finaloutput;  /* printing final output time */
  int match;   /* Matching time */
  int coarse;  /* Coarsening time */
  int refine; /* Refinement time */
  int coarsepart;  /* Coarse partitioning time */
  int project;    /* Project coarse-to-fine */
  int procred;    /* Processor reduction */
  int vcycle;    /* Vcycle time */
  int comerge;   /* Part of coarsening */
  int coshuffle;   /* Part of coarsening */
  int coremove;   /* Part of coarsening */
  int cotheend;   /* Part of coarsening */
  int matchstage[7];  /* Matching stages */
  int rdrdivide;  /* Rdivide time. */
  int rdbefore;   /* Part of Rdivide */
  int rdafter;    /* Part of Rdivide */
  int rdsplit;    /* Part of Rdivide */
  int rdredist;    /* Part of Rdivide */
  int rdsend;    /* Part of Rdivide */
  int rdwait;    /* Part of Rdivide */
  int rfrefine;    /* Refinement time */
  int rfpins;   /* Part of Refinement */
  int rfiso;   /* Part of Refinement */
  int rfgain;   /* Part of Refinement */
  int rfheap;   /* Part of Refinement */
  int rfpass;   /* Part of Refinement */
  int rfroll;   /* Part of Refinement */
  int rfnonroot;   /* Part of Refinement */
  int cpart;   /* Coarse partitioning time */
  int cpgather;   /* Part of Coarse Partitioning */
  int cprefine;   /* Part of Coarse Partitioning */
};

/**********************/
/* Matching functions */
/**********************/
int Zoltan_PHG_Matching (ZZ*, HGraph*, Matching, PHGPartParams*);
int Zoltan_PHG_Set_Matching_Fn (PHGPartParams*);
int Zoltan_PHG_Scale_Edges (ZZ*, HGraph*, float*, int);
int Zoltan_PHG_Scale_Vtx (ZZ*, HGraph*, PHGPartParams*);
int Zoltan_PHG_Vertex_Visit_Order (ZZ*, HGraph*, PHGPartParams*, int*);

/**************/
/* Coarsening */
/**************/
int Zoltan_PHG_Coarsening(ZZ*, HGraph*, Matching, HGraph*, int*, int*, int*,
 int**, struct Zoltan_Comm_Obj**, PHGPartParams*);
 
/*********************************/
/* Coarse Partitioning functions */
/*********************************/
extern int Zoltan_PHG_Gather_To_All_Procs(ZZ*, HGraph*, PHGPartParams*,
                                          PHGComm*, HGraph**);
extern int Zoltan_PHG_CoarsePartition(ZZ*, HGraph*, int, float *, Partition, 
                                      PHGPartParams*);
ZOLTAN_PHG_COARSEPARTITION_FN *Zoltan_PHG_Set_CoarsePartition_Fn(PHGPartParams*,
                                                                 int*);
/************************/
/* Refinement functions */ 
/************************/
int Zoltan_PHG_Refinement (ZZ*, HGraph*, int, float *, Partition, PHGPartParams*);
ZOLTAN_PHG_REFINEMENT_FN *Zoltan_PHG_Set_Refinement_Fn(char*);

/*******************/
/* 2D Distribution */
/*******************/
extern int Zoltan_PHG_Set_2D_Proc_Distrib(ZZ *, MPI_Comm, int, int, int, int,
                                          PHGComm *);

extern int Zoltan_PHG_Gno_To_Proc_Block(int gno, int*, int);
    
/*****************************/
/* Other Function Prototypes */
/*****************************/

extern int Zoltan_PHG_Initialize_Params(ZZ*, float *, PHGPartParams*);

extern void Zoltan_PHG_Free_Hypergraph_Data(ZHG *zoltan_hg);

extern int Zoltan_PHG_Fill_Hypergraph(ZZ*, ZHG*, PHGPartParams *, Partition*);
extern int Zoltan_PHG_Remove_Repart_Data(ZZ*, ZHG*, HGraph*, PHGPartParams*);

extern int Zoltan_PHG_rdivide (int,  int, Partition, ZZ *, HGraph *,
                               PHGPartParams *, int);
    
extern int Zoltan_PHG_Set_Part_Options(ZZ*, PHGPartParams*);
extern int Zoltan_PHG_Partition(ZZ*, HGraph*, int, float *, Partition, 
                                PHGPartParams*);
extern double Zoltan_PHG_Compute_NetCut(PHGComm*, HGraph*, Partition, int);
extern double Zoltan_PHG_Compute_ConCut(PHGComm*, HGraph*, Partition, int, 
                                        int*);    
extern int Zoltan_PHG_Removed_Cuts(ZZ *, ZHG *, double *);

extern double Zoltan_PHG_Compute_Balance(ZZ*, HGraph*, float *, int, 
                                         int, Partition);

extern int Zoltan_PHG_Build_Hypergraph(ZZ*, ZHG**, Partition*, PHGPartParams*);
extern void Zoltan_PHG_Plot(int, int, int, int*, int*, int*, char*);
extern void Zoltan_PHG_Plot_2D_Distrib(ZZ*, HGraph*);

extern int Zoltan_PHG_PaToH(ZZ *, HGraph *, int, int *, PHGPartParams*);    
extern int Zoltan_PHG_ParKway(ZZ *, HGraph *, int, Partition, PHGPartParams* );

/* Functions that are used in more than one PHG source file, 
   but not called from outside PHG. */
#include "zz_heap.h"  /* defines type HEAP */
int Zoltan_HG_move_vertex (HGraph *hg, int vertex, int sour, int dest,
    int *part, int **cut, double *gain, HEAP *heap);

int Zoltan_Convert_To_CSR( ZZ *zz, int num_pins, int *col_ptr,
    int *num_lists, ZOLTAN_ID_PTR *vtx_GID,
    int **row_ptr, ZOLTAN_ID_PTR *edg_GID);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __ZOLTAN_PHG_H */
