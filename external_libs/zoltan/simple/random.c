/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile: random.c,v $
 *    $Author: lafisk $
 *    $Date: 2007/05/21 21:12:40 $
 *    Revision: 1.4 $
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <memory.h>
#include "zz_const.h"
#include "zz_rand.h"
#include "params_const.h"
#include "all_allo_const.h"

/* Random partitioning, does not attempt balance! */
/* Each processor selects a random subset of objects (default is all)
   to give to random partitions (and processors). */

/*****************************************************************************/
/*  Parameters structure for Random method. */

static PARAM_VARS Random_params[] = {
                  { "RANDOM_MOVE_FRACTION", NULL, "DOUBLE", 0 },
                  { NULL, NULL, NULL, 0 } };
/*****************************************************************************/

int Zoltan_Random(
  ZZ *zz,                       /* The Zoltan structure.                     */
  float *part_sizes,            /* Input:  Array of size zz->LB.Num_Global_Parts
                                   * zz->Obj_Weight_Dim
                                   containing the percentage of work to be
                                   assigned to each partition.               */
  int *num_import,              /* Return -1. Random uses only export lists. */
  ZOLTAN_ID_PTR *import_global_ids, /* Not used. */
  ZOLTAN_ID_PTR *import_local_ids,  /* Not used. */
  int **import_procs,           /* Not used. */
  int **import_to_part,         /* Not used. */
  int *num_export,              /* Output: Number of objects to export. */
  ZOLTAN_ID_PTR *export_global_ids, /* Output: GIDs to export. */
  ZOLTAN_ID_PTR *export_local_ids,  /* Output: LIDs to export. */
  int **export_procs,           /* Output: Processsors to export to. */
  int **export_to_part          /* Output: Partitions to export to. */
)
{
  int ierr = ZOLTAN_OK;
  int i, count, num_obj;
  int max_export;
  double rand_frac = 1.0;       /* Default is to move all objects. */
  ZOLTAN_ID_PTR global_ids = NULL;
  ZOLTAN_ID_PTR local_ids = NULL; 
  int *parts = NULL;
  float *dummy = NULL;
  static char *yo = "Zoltan_Random";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* No import lists computed. */
  *num_import = -1;

  /* Get parameter values. */
  Zoltan_Bind_Param(Random_params, "RANDOM_MOVE_FRACTION", (void *) &rand_frac);
  Zoltan_Assign_Param_Vals(zz->Params, Random_params, zz->Debug_Level, 
                           zz->Proc, zz->Debug_Proc);

  /* Get list of local objects. */
  ierr = Zoltan_Get_Obj_List(zz, &num_obj, &global_ids, &local_ids, 0,
                             &dummy, &parts);

  /* Bound number of objects to export. */
  max_export = 1.5*rand_frac*num_obj;

  /* Allocate export lists. */
  *export_global_ids = *export_local_ids = NULL;
  *export_procs = *export_to_part = NULL;
  if (max_export > 0) {
    if (!Zoltan_Special_Malloc(zz, (void **)export_global_ids, max_export,
                               ZOLTAN_SPECIAL_MALLOC_GID)
     || !Zoltan_Special_Malloc(zz, (void **)export_local_ids, max_export,
                               ZOLTAN_SPECIAL_MALLOC_LID)
     || !Zoltan_Special_Malloc(zz, (void **)export_procs, max_export,
                               ZOLTAN_SPECIAL_MALLOC_INT)
     || !Zoltan_Special_Malloc(zz, (void **)export_to_part, max_export,
                               ZOLTAN_SPECIAL_MALLOC_INT)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }

  /* Randomly assign ids to procs. */
  count=0;
  for (i=0; i<num_obj; i++){
    /* Randomly select some objects to move (export) */
    if ((count<max_export) && (Zoltan_Rand(NULL)<rand_frac*ZOLTAN_RAND_MAX)){
      /* export_global_ids[count] = global_ids[i]; */
      ZOLTAN_SET_GID(zz, &((*export_global_ids)[count*zz->Num_GID]),
                     &global_ids[i*zz->Num_GID]);
      if (local_ids)
        /* export_local_ids[count] = local_ids[i]; */
        ZOLTAN_SET_LID(zz, &((*export_local_ids)[count*zz->Num_LID]),
                       &local_ids[i*zz->Num_LID]);
      /* Randomly pick new partition number. */
      (*export_to_part)[count] = Zoltan_Rand_InRange(NULL, zz->LB.Num_Global_Parts);
      /* Processor number is derived from partition number. */
      (*export_procs)[count] = Zoltan_LB_Part_To_Proc(zz, 
                     (*export_to_part)[count], &global_ids[i*zz->Num_GID]);

      /* printf("Debug: Export gid %u to part %d and proc %d.\n", (*export_global_ids)[count], (*export_to_part)[count], (*export_procs)[count]); */

      ++count;
    }
  }
  (*num_export) = count;

End:
  /* Free local memory, but not export lists. */
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);
  ZOLTAN_FREE(&parts);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

