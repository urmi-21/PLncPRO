/****************************************************************\
*  ESTate - Expressed sequence tag analysis tools etc.           *
*  Copyright (C) 1996-1999.  Guy St.C. Slater.                   *
*  All Rights Reserved.                                          *
*                                                                *
*  gslater@hgmp.mrc.ac.uk  http://www.hgmp.mrc.ac.uk/~gslater/   *
*                                                                *
*  This source code is distributed under the terms of the        *
*  GNU General Public License. See the file COPYING for details. *
*                                                                *
*  If you use this code, please keep this notice intact.         *
\****************************************************************/

#ifndef INCLUDED_INTDICT_C
#define INCLUDED_INTDICT_C

#include "intdict.h"
#include "../general/common.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>  /* FOR CHAR_BITS */

#define INTDICT_MAXDEPTH (sizeof(int)*CHAR_BIT)
#define INTDICTother(n) ((~(n))&1)

/* INTDICT_CHUNKSIZE MUST BE GREATER THAN INTDICT_MAXDEPTH 
*/
#define INTDICTNODE_CHUNKSIZE (1000 + INTDICT_MAXDEPTH + 1)

INTDICT *newINTDICT(){
    register INTDICT *id = NEW(INTDICT);
    id->bufleft = INTDICTNODE_CHUNKSIZE;
    id->buffer = malloc(sizeof(INTDICTNODE)*(id->bufleft-=2));
    id->buffer->p[0] = NULL;
    id->root = id->buffer+1;
    id->root->p[0] = id->root->p[1] = NULL;
    return id;
    }

void freeINTDICT(INTDICT *id){
    INTDICTNODE *curr = id->buffer, *next = curr->p[0];
    do {
        free(curr);
        curr = next;
        if(!curr)
            break;
        next = curr->p[0];
    } while(TRUE);
    free(id);
    return;
    }

static INTDICTNODE *newINTDICTNODES(INTDICT *id, int num){
    register INTDICTNODE *idn;
    if(id->bufleft < num){
        id->bufleft = INTDICTNODE_CHUNKSIZE;
        idn = id->buffer;
        id->buffer = malloc(sizeof(INTDICTNODE)*id->bufleft--);
        id->buffer->p[0] = idn;
        }
    idn = id->buffer+(INTDICTNODE_CHUNKSIZE-id->bufleft);
    id->bufleft -= num;
    return idn;
    }

static INTDICTNODE *makeINTDICTNODEtail(INTDICT *id, int k, void *v,
                                        int len){
    register int i = 0, z = len-1, s;
    register INTDICTNODE *n = newINTDICTNODES(id, len);
    do {
        s = (k>>=1)&1;
        n[i].p[INTDICTother(s)] = NULL;
        n[i].p[s] = &n[i+1];
        i++;
    } while(i < z);
    s = (k>>1)&1;
    n[i].v[INTDICTother(s)] = NULL;
    n[i].v[s] = v;
    return n;
    }

/* uniqAddINTDICT : ADD INT:VALUE (k,v) PAIR TO id. 
                  id UNCHANGED IF NAME PRESENT, RETURNS CURRENT VALUE.
                  OTHERWISE CHANGES id, RETURNS NULL.
*/
void *uniqAddINTDICT(INTDICT *id, int k, void *v){
    register INTDICTNODE **ptp = &id->root->p[k&1];
    register int i = INTDICT_MAXDEPTH-1;
    do {
        if(!*ptp){
            *ptp = makeINTDICTNODEtail(id, k, v, i);
            return NULL;
            }
        ptp = &(*ptp)->p[(k>>=1)&1];
    } while(--i);
    if(!*ptp){ /* ADD LATE */
        *ptp = v;
        return NULL;
        }
    return *ptp; /* ALREADY PRESENT */
    }

/* replaceAddINTDICT : ADD INT:VALUE PAIR TO id, REPLACING ANY
                       PREVIOUS VALUE, WHICH IS RETURNED.
*/
void *replaceAddINTDICT(INTDICT *id, int k, void *v){
    register INTDICTNODE **ptp = &id->root->p[k&1];
    register int i = INTDICT_MAXDEPTH-1;
    register void *prev;
    do {
        if(!*ptp){
            *ptp = makeINTDICTNODEtail(id, k, v, i);
            return NULL;
            }
        ptp = &(*ptp)->p[(k>>=1)&1];
    } while(--i);
    if(!*ptp){ /* ADD LATE */
        *ptp = v;
        return NULL;
        }
    prev = (void*)*ptp;
    *ptp = v; /* REPLACE VALUE */
    return prev; /* RETURN PREVIOUS VALUE */
    }

/* lookupINTDICT : RETURNS VALUE CORRESPONDING TO k 
                   OR NULL IF CANNOT BE FOUND IN id.
*/
void *lookupINTDICT(INTDICT *id, int k){
    register INTDICTNODE *curr = id->root->p[k&1];
    register int i = INTDICT_MAXDEPTH-2;
    do {
        if(!curr)
            return NULL; /* DID NOT FIND KEY */
        curr = curr->p[(k>>=1)&1];
    } while(--i);
    return curr->v[(k>>1)&1]; /* MAYBE FOUND KEY */
    }

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

int main(){
    register INTDICT *id = newINTDICT();

    uniqAddINTDICT(id, 10, "1010");
    uniqAddINTDICT(id,  9, "1001");
    uniqAddINTDICT(id,  5, "0101");
    uniqAddINTDICT(id,  1, "0001");
    uniqAddINTDICT(id,  7, "0111");
    uniqAddINTDICT(id, 12, "1100");
    uniqAddINTDICT(id,  2, "0010");
    uniqAddINTDICT(id, 13, "1101");
    uniqAddINTDICT(id, 11, "1011");
    uniqAddINTDICT(id,  8, "1000");
    uniqAddINTDICT(id,  0, "0000");
    uniqAddINTDICT(id,  3, "0011");
    uniqAddINTDICT(id, 15, "1111");
    uniqAddINTDICT(id,  6, "0110");
    uniqAddINTDICT(id,  4, "0100");
    uniqAddINTDICT(id, 14, "1110");

    uniqAddINTDICT(id, 214, "11010110");
    uniqAddINTDICT(id, 215, "11010111");
    uniqAddINTDICT(id,  86, "01010110");

    printf("Found [%s]\n", (char*)lookupINTDICT(id, 10) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id,  9) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id,  5) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id,  1) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id,  7) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id, 12) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id,  2) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id, 13) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id, 11) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id,  8) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id,  0) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id,  3) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id, 15) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id,  6) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id,  4) );

    printf("Found [%s]\n", (char*)lookupINTDICT(id, 214) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id, 215) );
    printf("Found [%s]\n", (char*)lookupINTDICT(id, 86) );

    freeINTDICT(id);
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_INTDICT_C */

/* TODO:
   Rewrite so that it uses only goes to the required depth.
   include lookupNearestAfter() and lookupNearestBefore();
*/

