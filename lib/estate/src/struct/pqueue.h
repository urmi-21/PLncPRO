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

/* Priority Queue. (Algorithm based on leftist trees.)  
   Guy St.C. Slater..  April 1997.

   Reference: 
       "The Art of Computer Programming", Donald E. Knuth.
        Volume 3, "Sorting and searching", Pages 150-152.
       (First Edition). 
*/

#ifndef INCLUDED_PQUEUE_H
#define INCLUDED_PQUEUE_H

#include <stdlib.h>
#include "../general/common.h"

typedef struct PQNODE {
             void *data;   /* DATA STORED IN TREE    */
              int  dist;   /* DISTANCE TO BLANK NODE */
    struct PQNODE *left;   /* LEFT BRANCH            */
    struct PQNODE *right;  /* RIGHT BRANCH           */
    } PQNODE;

typedef struct {
        PQNODE *root;    /* THE LEFTIST TREE                */
           int  total;   /* CURRENT NUMBER OF NODES         */
           int  max;     /* MAXIMUM NUMBER OF NODES TO KEEP */
    PQCOMPFUNC  comp;    /* POINTER TO COMPARISON FUNCTION  */
    } PQUEUE;

PQUEUE * newPQUEUE(int max, PQCOMPFUNC comp);
  void  joinPQUEUE(PQUEUE *a, PQUEUE *b);
  void  joinLimitedPQUEUE(PQUEUE *a, PQUEUE *b);
  void * popPQUEUE(PQUEUE *pq);
  void  pushPQUEUE(PQUEUE *pq, void *data);
  void *pushLimitedPQUEUE(PQUEUE *pq, void *data);
  void  freePQUEUE(PQUEUE *pq);
  void  freePQUEUEptr(PQUEUE *pq, FREEFUNC ff);
  void  bestFirstPopPQUEUE(PQUEUE *pq, PQPOPFUNC func, void *v);
#define  viewTopPQUEUE(pq) ( (pq)->root?(pq)->root->data:NULL)
#define   atlimitPQUEUE(pq) ( (pq)->total == (pq)->max)

#endif /* INCLUDED_PQUEUE_H */

