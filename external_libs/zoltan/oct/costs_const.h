/*****************************************************************************
 * CVS File Information :
 *    $RCSfile: costs_const.h,v $
 *    $Author: kddevin $
 *    $Date: 2002/06/19 23:56:40 $
 *    Revision: 1.7 $
 ****************************************************************************/

#ifndef __COSTS_CONST_H
#define __COSTS_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


extern void  Zoltan_Oct_costs_free(OCT_Global_Info *OCT_info,pOctant octree); 
extern float Zoltan_Oct_costs_value(pOctant octant);
extern float Zoltan_Oct_costs_global_compute(OCT_Global_Info *OCT_info);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
