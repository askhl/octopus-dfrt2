/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * This software is distributed under the GNU Lesser General Public License. *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile: DD_Set_Neighbor_Hash_Fn1.c,v $
 *    $Author: rheaphy $
 *    $Date: 2007/09/12 17:53:45 $
 *    Revision: 1.9 $
 ****************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include "DD.h"



#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


static unsigned int dd_nh1 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int nproc) ;

static int max_gid ;
static int groupsize ;


/*************  Zoltan_DD_Set_Hash_Fn1() ***********************/
/*
**  These routines associate the first n=groupsize GIDs to proc 0, the
**  next n to proc 1, etc.  It assumes the GIDs are consecutive numbers.
**  It assumes that GIDs primarily stay near their original owner. The
**  GID length is assumed to be 1. GIDs outside of range are evenly
**  distributed among the processors via modulo(nproc).
*/


int Zoltan_DD_Set_Neighbor_Hash_Fn1 (
 Zoltan_DD_Directory *dd,          /* directory state information */
 int size)                         /* number of reserved GIDs per CPU */
   {
   char *yo = "Zoltan_DD_Set_Hash_Fn1";

   if (dd == NULL || size < 1)  {
      ZOLTAN_PRINT_ERROR (0, yo, "Invalid input argument");
      return ZOLTAN_FATAL;
   }

   groupsize   = size;
   dd->hash    = dd_nh1;
   dd->cleanup = NULL;                 /* no need to free anything */
   max_gid     = size * dd->nproc;     /* larger GIDs out of range */

   return ZOLTAN_OK;
   }



static unsigned int dd_nh1 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int nproc)
   {
   int id = (signed) *gid;
   return (id < max_gid) ? (id / groupsize) : (id % nproc);
   }

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
