/*****************************************************************************
 * CVS File Information :
 *    $RCSfile: migreg.h,v $
 *    $Author: lafisk $
 *    $Date: 2004/12/03 15:40:27 $
 *    Revision: 1.14 $
 ****************************************************************************/

#ifndef __MIGREG_H
#define __MIGREG_H

#include "octant_const.h"
#include "octupdate_const.h"
#include "oct_util_const.h"
#include "migreg_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


typedef struct
{
  pRegion region;
  int npid;
} Message;


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
