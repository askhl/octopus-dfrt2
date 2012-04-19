/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile: g2l_hash.h,v $
 *    $Author: uvcatal $
 *    $Date: 2006/06/19 21:58:54 $
 *    Revision: 1.3 $
 ****************************************************************************/
#ifndef _G2L_HASH_H_
#define _G2L_HASH_H_

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Structure used for hashing */
struct G2L_Hash_Node {
    int gno;           /* Global number */
    int lno;           /* Mapped id of gno*/
    struct G2L_Hash_Node * next;
};

typedef struct G2L_Hash_Node G2LHashNode;

struct G2L_Hash {
    int   size;
    int   lastlno;
    
    G2LHashNode **table;
    G2LHashNode *nodes;
};

typedef struct G2L_Hash G2LHash;

/* Returns the prime number closest to (and smaller than) stop */
int Zoltan_GenPrime(int stop, int *prime_num);
    

int Zoltan_G2LHash_Create(G2LHash *hash, int size);
int Zoltan_G2LHash_Destroy(G2LHash *hash);
int Zoltan_G2LHash_G2L(G2LHash *hash, int gno);
/*
  if gno exist it returns lno, if it does not exist,
  it inserts andr returns newly assigned lno */
int Zoltan_G2LHash_Insert(G2LHash *hash, int gno);
#define Zoltan_G2LHash_L2G(hash, lno) ((hash)->nodes[lno].gno)


/* Key&Value hash functions using same datastructure above
   the only difference will be the insert function */
#define KVHash  G2LHash

#define Zoltan_KVHash_Create(hash, size)   Zoltan_G2LHash_Create(hash, size)
#define Zoltan_KVHash_Destroy(hash)        Zoltan_G2LHash_Destroy(hash)


int Zoltan_KVHash_Insert(KVHash *hash, int key, int value);
#define Zoltan_KVHash_GetValue(hash, key)  Zoltan_G2LHash_G2L(hash, key)

    
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
    

#endif

