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

/* Pairing heaps.
   Guy St.C. Slater.. October 1998.

   Simple, efficient priority queues,
   supporting 'decrease key' operations.

   Reference:
       "The Pairing Heap: A New Form of Self-Adjusting Heap"
       Michael L. Fredman, Robert Sedgewick,
       Daniel D. Sleator, and Robert E. Tarjan.
       Algorithmica. 1: 111-129. (1986).
*/
#ifndef INCLUDED_PAIRHEAP_H
#define INCLUDED_PAIRHEAP_H

#include "../general/common.h"

typedef struct PAIRHEAPNODE {
                void  *data;
                void  *priority;
struct PAIRHEAPNODE   *left;
struct PAIRHEAPNODE   *next;
struct PAIRHEAPNODE   *prev;
    } PAIRHEAPNODE;

typedef struct {
    PAIRHEAPNODE  *root;         /* THE PAIRING HEAP          */
             int   total;        /* NUMBER OF MEMBERS         */
      PQCOMPFUNC   comp;         /* COMPARISON FUNCTION       */
    PAIRHEAPNODE **combine;      /* COMBINATION ARRAY         */
             int   combinealloc; /* ALLOCATED COMBINE SIZE    */
       FREEFUNC    freedata;     /* FUNCTION TO FREE DATA     */
       FREEFUNC    freepriority; /* FUNCTION TO FREE PRIORITY */
    } PAIRHEAP;

    PAIRHEAP *  newPAIRHEAP(PQCOMPFUNC comp,
                            FREEFUNC freedata, FREEFUNC freepri);
        void   freePAIRHEAP(PAIRHEAP *ph);
PAIRHEAPNODE * pushPAIRHEAP(PAIRHEAP *ph, void *data, void *pri);
        void  raisePAIRHEAP(PAIRHEAP *ph, PAIRHEAPNODE *acc, void *pri);
        void *  popPAIRHEAP(PAIRHEAP *ph);

#define occupiedPAIRHEAP(ph) ((ph)->total)
#define     compPAIRHEAP(ph) ((ph)->comp)
#define   toppriPAIRHEAP(ph) ((ph)->root->priority)
#define nodedataPAIRHEAPNODE(phn) ((phn)->data)
#define  nodepriPAIRHEAPNODE(phn) ((phn)->priority)

#endif /* INCLUDED_PAIRHEAP_H */

