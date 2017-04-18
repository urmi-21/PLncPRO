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

#ifndef INCLUDED_LIST_C
#define INCLUDED_LIST_C

#include "list.h"

LIST *newLIST(){
    register LIST *l = NEW(LIST);
    l->a = l->z = NEW(LISTNODE);
    return l;
    }

void queueLIST(LIST *l, void *v){
    l->z->n = NEW(LISTNODE);
    l->z = l->z->n;
    l->z->v = v;
    return;
    }

void stackLIST(LIST *l, void *v){
    register LISTNODE *t = NEW(LISTNODE);
    l->a->v = v;
    t->n = l->a;
    l->a = t;
    return;
    }

void *popLIST(LIST *l){
    register void *v = l->a->n->v;
    register LISTNODE *t = l->a;
    l->a = t->n;
    free(t);
    return v;
    }

void freeLIST(LIST *l){
    while(occupiedLIST(l))
        popLIST(l);
    free(l->a);
    free(l);
    return;
    }

void freeLISTptr(LIST *l, FREEFUNC ff){
    while(occupiedLIST(l))
        ff(popLIST(l));
    free(l->a);
    free(l);
    return;
    }

LIST *joinLIST(LIST *l, LIST *r){
    l->z->n = r->a->n;
    l->z = r->z;
    free(r->a);
    free(r);
    return l;
    }

#define LIST_ARRAY_CHUNK_SIZE 100

void **toarrayLIST(LIST *l, int *c){
    register int alloc = LIST_ARRAY_CHUNK_SIZE;
    register void **v = malloc(sizeof(void*)*alloc);
    *c = 0;
    while(occupiedLIST(l)){
        v[(*c)++] = popLIST(l);
        if((*c) > alloc){
            alloc+=LIST_ARRAY_CHUNK_SIZE;
            v = realloc(v, sizeof(void*)*alloc);
            }
        }
    free(l->a);
    free(l); 
    return v;
    }

LIST *arraytoLIST(void **v, int c){
    register LIST *l = newLIST();
    while(c--)
        stackLIST(l, v[c]);
    free(v); 
    return l;
    }

/* START OF CLIST CODE */

CLIST *newCLIST(){
    register CLIST *cl = NEW(CLIST);
    cl->list = newLIST();
    cl->count = 0;
    return cl;
    }

void queueCLIST(CLIST *cl, void *v){
    cl->count++;
    queueLIST(cl->list, v);
    return;
    }

void stackCLIST(CLIST *cl, void *v){
    cl->count++;
    stackLIST(cl->list, v);
    return;
    }

void *popCLIST(CLIST *cl){
    cl->count--;
    return popLIST(cl->list);
    }

void freeCLIST(CLIST *cl){
    freeLIST(cl->list);
    free(cl);
    return;
    }

CLIST *joinCLIST(CLIST *a, CLIST *b){
    a->count += b->count;
    a->list = joinLIST(a->list, b->list);
    free(b);
    return a;
    }

void **toarrayCLIST(CLIST *cl, int *c){
    register void **v = toarrayLIST(cl->list, c); 
    free(cl);
    return v;
    }

CLIST *arraytoCLIST(void **v, int c){
    register CLIST *cl = NEW(CLIST);
    cl->list = arraytoLIST(v, c);
    cl->count = c;
    return cl;
    }


/* END OF CLIST CODE   */


#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

#include <stdio.h>

int main(){
    LIST *a = newLIST();
    LIST *b = newLIST();
    char **word;
    int count, i;
    queueLIST(a, "jumps");
    stackLIST(a, "fox");
    queueLIST(a, "over");
    stackLIST(a, "brown");
    queueLIST(a, "a");
    stackLIST(b, "quick");
    queueLIST(a, "lazy");
    stackLIST(b, "The");
    queueLIST(a, "dog");
    a = joinLIST(b,a);

    word = (char**)toarrayLIST(a, &count); 
    printf("get count [%d]\n", count);
    for(i = 0; i < count; i++)
        printf("%s ", word[i]);
    printf("\n"); 
    a = arraytoLIST((void**)word, count);
    while(occupiedLIST(a))
        printf("%s ", (char*)popLIST(a));
    printf("\n");
    freeLIST(a);
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 


#endif /* INCLUDED_LIST_C  */

