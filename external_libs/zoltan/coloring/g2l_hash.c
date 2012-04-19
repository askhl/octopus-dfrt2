/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile: g2l_hash.c,v $
 *    $Author: kddevin $
 *    $Date: 2008/09/16 22:05:16 $
 *    Revision: 1.6.8.1 $
 ****************************************************************************/
#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zoltan_util.h"
#include "coloring.h"    
#include "g2l_hash.h"



/*****************************************************************************/
/* Returns the prime number closest to (and smaller than) stop */
int Zoltan_GenPrime(int stop, int *prime_num)
{
    int nap;
    int num, c, i;
    int *prime;
    
    prime = (int *) ZOLTAN_MALLOC(((stop/2) > 2 ? (stop/2) : 2) * sizeof(int));
    if (!prime)
        return ZOLTAN_MEMERR;
    
    prime[0] = 2; prime[1] = 3;
    c = 2; /* initial primes */
    
    /* only have to check odd numbers */
    for (num=5; num < stop; num = num + 2) {
        nap = 0;  /* set not-a-prime false */        
        /* cycle through list of known primes */
        for (i=0; i < c; i++) { 
            /* check if a previous prime divides evenly */
            /* if so the number is not a prime */
            if ((num % prime[i]) == 0) {
                nap = 1;
                break;
            }            
            /* stop if prime squared is bigger than the number */
            if ((prime[i] * prime[i]) > num)
                break;
        }        
        /* if not-a-prime, then we found a prime */
        if (nap != 1) {
           /* add prime to list of known primes */
            prime[c] = num;
            c++;
        }
    }
    *prime_num = prime[c-1];

    ZOLTAN_FREE(&prime);
    return ZOLTAN_OK;
}



int Zoltan_G2LHash_Create(G2LHash *hash, int size)
{
    if (size == 0) /* to avoid memory allocation errors */
        size = 1;
    
    hash->table = NULL;
    hash->nodes = NULL;
    hash->size = size;
    hash->lastlno = 0;
    hash->table = (G2LHashNode **) ZOLTAN_CALLOC(size, sizeof(G2LHashNode *));
    hash->nodes = (G2LHashNode *) ZOLTAN_MALLOC(size * sizeof(G2LHashNode));
    if (!hash->table || !hash->nodes) {
        Zoltan_G2LHash_Destroy(hash);
        return ZOLTAN_MEMERR;
    }
    return ZOLTAN_OK;
}

int Zoltan_G2LHash_Destroy(G2LHash *hash)
{
    ZOLTAN_FREE(&hash->table);
    ZOLTAN_FREE(&hash->nodes);

    return ZOLTAN_OK;
}

int Zoltan_G2LHash_Insert(G2LHash *hash, int gno)
{
    int i, lno;
    G2LHashNode *ptr;

    i = Zoltan_Hash((ZOLTAN_ID_PTR) &gno, 1, (unsigned int) hash->size);
    for (ptr=hash->table[i]; ptr && ptr->gno!=gno; ptr = ptr->next);
    if (!ptr) {
        lno = hash->lastlno++;

        if (lno >= hash->size) {
            char st[2048];
            sprintf(st, "Hash is full! #entries=%d  size=%d", lno, hash->size);
            ZOLTAN_PRINT_ERROR(-1, "Zoltan_G2LHash_G2L", st);
            return -1;
        }
        
        ptr = &(hash->nodes[lno]);
        ptr->gno = gno;
        ptr->lno = lno;
        ptr->next = hash->table[i];
        hash->table[i] = ptr;
    } else
        lno = ptr->lno;

    return lno;
}

int Zoltan_G2LHash_G2L(G2LHash *hash, int gno)
{
    int i;
    G2LHashNode *ptr;

    i = Zoltan_Hash((ZOLTAN_ID_PTR) &gno, 1, (unsigned int) hash->size);
    for (ptr=hash->table[i]; ptr && ptr->gno!=gno; ptr = ptr->next);
    if (!ptr)
        return -1;
    else
        return ptr->lno;
}




int Zoltan_KVHash_Insert(KVHash *hash, int key, int value)
{
    int i, lno;
    
    G2LHashNode *ptr;

    i = Zoltan_Hash((ZOLTAN_ID_PTR) &key, 1, (unsigned int) hash->size);
    for (ptr=hash->table[i]; ptr && ptr->gno!=key; ptr = ptr->next);
    if (!ptr) {
        lno = hash->lastlno++;

        if (lno >= hash->size) {
            ZOLTAN_PRINT_ERROR(-1, "Zoltan_KVHash_Insert", "Hash is full!");
            return -1;
        }
        
        ptr = &(hash->nodes[lno]);
        ptr->gno = key;
        ptr->lno = value;
        ptr->next = hash->table[i];
        hash->table[i] = ptr;
    } else
        value = ptr->lno;

    return value;   
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
