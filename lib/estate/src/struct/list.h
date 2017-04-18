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

/* Simple queue and stack based linked list.
   Guy St.C. Slater.. April 1997.
*/
#ifndef INCLUDED_LIST_H
#define INCLUDED_LIST_H

typedef struct LISTNODE {
    struct LISTNODE *n; void *v;
} LISTNODE;

typedef struct { LISTNODE *a, *z; } LIST;

#define occupiedLIST(l) ((l)->a != (l)->z)

#include "../general/common.h"

LIST   *newLIST();
void  queueLIST(LIST *l, void *v);
void  stackLIST(LIST *l, void *v);
void  *popLIST(LIST *l);
void   freeLIST(LIST *l);
void   freeLISTptr(LIST *l, FREEFUNC ff);
LIST  *joinLIST(LIST *a, LIST *b);
void **toarrayLIST(LIST *l, int *c);
LIST  *arraytoLIST(void **v, int c);

#define viewEndLIST(list) ((list)->a->n->v)

/* MACRO VERSION OF queueLIST */
#define MqueueLIST(l,d) ((l)->z->n=NEW(LISTNODE), \
                         (l)->z->n->v=(d),        \
                         (l)->z=(l)->z->n)

typedef struct {
   LIST *list;
    int  count;
    } CLIST;

CLIST   *newCLIST();
void  queueCLIST(CLIST *cl, void *v);
void  stackCLIST(CLIST *cl, void *v);
void  *popCLIST(CLIST *cl);
void   freeCLIST(CLIST *cl);
CLIST  *joinCLIST(CLIST *a, CLIST *b);
void **toarrayCLIST(CLIST *cl, int *c);
CLIST  *arraytoCLIST(void **v, int c);

/* MACRO VERSION OF queueCLIST */
#define MqueueCLIST(cl,d) ((cl)->count++, \
                          MqueueLIST((cl)->list, (d)))

#define countCLIST(cl) ((cl)->count)
#define occupiedCLIST(cl) (occupiedLIST((cl)->list))

#endif /* INCLUDED_LIST_H */

