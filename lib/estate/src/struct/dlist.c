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

#ifndef INCLUDED_DLIST_C
#define INCLUDED_DLIST_C

#include "dlist.h"
#include "../general/common.h"

DLIST *newDLIST(){
    register DLIST *d = NEW(DLIST);
    d->a = &d->sa;
    d->z = &d->sz;
    d->a->n = d->z;
    d->z->p = d->a;
    return d;
    }

void freeDLIST(DLIST *d){
    while(occupiedDLIST(d))
        popDLISThead(d);
    free(d);
    return;
    }

void *popDLIST(DLIST *d, DLISTNODE *n){
    register void *v = n->v;
    n->p->n = n->n;
    n->n->p = n->p;
    free(n);
    return v;
    }

void *popDLISThead(DLIST *d){
    register DLISTNODE *t = d->a->n;
    register void *v = t->v;
    d->a->n = t->n;
    d->a->n->p = d->a;
    free(t);
    return v;
    }

void *popDLISTtail(DLIST *d){
    register DLISTNODE *t = d->z->p;
    register void *v = t->v;
    d->z->p = t->p;
    d->z->p->n = d->z;
    free(t);
    return v;
    }

DLISTNODE *pushDLISThead(DLIST *d, void *v){
    register DLISTNODE *t = NEW(DLISTNODE);
    t->v = v;
    t->p = d->a;
    t->n = d->a->n;
    d->a->n = t->n->p = t;
    return t;
    }

DLISTNODE *pushDLISTtail(DLIST *d, void *v){
    register DLISTNODE *t = NEW(DLISTNODE);
    t->v = v;
    t->n = d->z;
    t->p = d->z->p;
    d->z->p = t->p->n = t;
    return t;
    }

DLISTNODE *pushDLISTbefore(DLIST *d, DLISTNODE *n, void *v){
    register DLISTNODE *t = NEW(DLISTNODE);
    t->v = v;
    t->n = n;
    t->p = n->p;
    n->p = t->p->n = t;
    return t;
    }

DLISTNODE *pushDLISTafter(DLIST *d, DLISTNODE *n, void *v){
    register DLISTNODE *t = NEW(DLISTNODE);
    t->v = v;
    t->p = n;
    t->n = n->n;
    n->n = t->n->p = t;
    return t;
    }

void raiseDLIST(DLIST *d, DLISTNODE *n){ /* PUT TO START */
    n->p->n = n->n;
    n->n->p = n->p;
    n->p = d->a;
    n->n = d->a->n;
    d->a->n = n->n->p = n;
    return;
    }

void lowerDLIST(DLIST *d, DLISTNODE *n){ /* PUT TO END */
    n->p->n = n->n;
    n->n->p = n->p;
    n->n = d->z;
    n->p = d->z->p;
    d->z->p = n->p->n = n;
    return;
    }

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

#include <stdio.h>

int main(){
    register DLIST *d = newDLIST();
    pushDLISThead(d, "not this");
    pushDLISTtail(d, "nor this");
    popDLIST(d, pushDLISThead(d, "or this"));
    popDLISThead(d);
    popDLISTtail(d);
    pushDLISThead(d, "the");
    pushDLISTtail(d, "quick");
    lowerDLIST(d, pushDLISThead(d, "brown"));
    pushDLISTbefore(d, pushDLISTtail(d, "jumps"), "fox");
    pushDLISTafter(d, pushDLISTtail(d, "over"), "the");
    lowerDLIST(d, pushDLISThead(d, "lazy"));
    lowerDLIST(d, pushDLISThead(d, "dog."));
    while(occupiedDLIST(d))
        fprintf(stderr, "Popped [%s]\n", (char*)popDLISThead(d));
    freeDLIST(d);
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */

#endif /* INCLUDED_DLIST_C  */

