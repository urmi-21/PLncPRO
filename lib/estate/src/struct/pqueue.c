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

#ifndef INCLUDED_PQUEUE_C
#define INCLUDED_PQUEUE_C

#include "../struct/pqueue.h"
#include "../general/error.h"

ALLOW_ERROR_MESSAGES;

/* 
PROPERTIES OF THE LEFTIST TREE: 
------------------------------
[1] ANY KEY IS GREATER THAN ANY OF IT'S DESCENDANTS.
     -> THUS ROOT WILL HAVE THE LARGEST KEY.
[2] DIST IS DISTANCE TO NEXT FREE NODE   
     -> FOUND BY MIN(LEFT DIST, RIGHT DIST) + 1
[3] LEFT DIST IS ALWAYS MORE THAN RIGHT DIST      
     -> THUS FIND A FREE NODE BY MOVING TO THE RIGHT
*/

/* newPQUEUE : MAKE A NEW PRIORITY QUEUE, 
               WILL BE MAINTAINED WITH UP TO max NODES.
               COMPARISONS WILL BE MADE WITH THE FUNCTION comp
*/
PQUEUE *newPQUEUE(int max, PQCOMPFUNC comp){
    register PQUEUE *pq = NEW(PQUEUE);
    pq->total = 0;
    pq->root  = NULL;
    pq->max   = max; 
    pq->comp  = comp;
    return pq;
    }

#define getPQUEUEdist(node) ((node)?(node->dist):0)

/* updatePQNODE : UPDATE THE DISTANCE NODES.
                 (USED FOR TREE INSERTION AND JOINING)
*/
static void updatePQNODE(PQNODE *n){
    register PQNODE *swap;
    if(getPQUEUEdist(n->left) < getPQUEUEdist(n->right))
        Swap(n->left, n->right, swap);
    n->dist = getPQUEUEdist(n->right)+1;
    return;
    }

/* joinPQNODErecur : RECURSIVE ELEMENT OF joinPQNODE
*/
static PQNODE *joinPQNODErecur(PQCOMPFUNC comp, PQNODE *a, PQNODE *b){
    register PQNODE *x, *y;
    if(!a)return b;
    comp(a->data, b->data)?(x=a,y=b):(x=b,y=a);
    x->right = joinPQNODErecur(comp, x->right, y);
    updatePQNODE(x);
    return x;
    }

/* joinPQNODE : JOIN NODES a AND b.
*/
static PQNODE *joinPQNODE(PQCOMPFUNC comp, PQNODE *a, PQNODE *b){
    return (b)?(joinPQNODErecur(comp, a, b)):a;
    }

/*  joinPQUEUE : JOIN QUEUES a AND b.
*/
void joinPQUEUE(PQUEUE *a, PQUEUE *b){
    if(a->comp != b->comp)
        errmsg(ERROR_FATAL, "Joining priority queues with different"
                            " comparison functions");
    a->root = joinPQNODE(a->comp, a->root, b->root);
    a->total+=b->total;
    free(b);
    return;
    }

/* joinLimitedPQUEUE : JOIN QUEUES a AND b, KEEPING SIZE LARGEST LIMIT.
*/
void joinLimitedPQUEUE(PQUEUE *a, PQUEUE *b){
    register int newmax = Max(a->max, b->max);
    joinPQUEUE(a,b);
    a->max = newmax;
    while(a->total > a->max)
        popPQUEUE(a);
    return;
    }

/*    popPQUEUE : RETURN LOWEST PRIORITY data.
*/
void *popPQUEUE(PQUEUE *pq){
    register PQNODE *temp;
    register void *data;
    if(pq->total < 1) /* IF TREE EMPTY */
        return NULL;
    data  = pq->root->data;
    temp = joinPQNODE(pq->comp, pq->root->left, pq->root->right);
    free(pq->root);
    pq->root = temp;
    pq->total--;
    return data;
    }

/* pushPQNODE : ADD b TO a.
*/
static PQNODE *pushPQNODE(PQCOMPFUNC comp, PQNODE *a, PQNODE *b){
    if(!b) return a; 
    if(comp(b->data, a->data)){ 
        b->right = pushPQNODE(comp, a, b->right); 
        updatePQNODE(b);
        return b;
        }
    a->left = b;
    return a;
    }

/* pushPQUEUE : ADD data TO pq. 
*/
void pushPQUEUE(PQUEUE *pq, void *data){
    register PQNODE *pqn = NEW(PQNODE);
    pqn->data  = data;
    pqn->left = pqn->right = NULL; 
    pq->root  = pq->root?pushPQNODE(pq->comp, pqn, pq->root):pqn; 
    pq->total++;
    return;
    }

/* pushLimitedPQUEUE : ADD data WITH KEY TO pq.
                       KEEP QUEUE SIZE DOWN TO MAX. 
                       RETURNS THE DATA DISCARDED, 
                       OR NULL IF NO DATA DISCARDED.
*/
void *pushLimitedPQUEUE(PQUEUE *pq, void *data){
    register PQNODE *pqn;
    register void *prevdata;
    if(pq->total < pq->max){           /* QUEUE NOT FULL */
        pushPQUEUE(pq, data);
        return NULL;
        }
    if(pq->comp(data, pq->root->data)) /* NOT INSERTING */
        return data;
    pqn = pq->root;                    /* REPLACE INSERT */ 
    prevdata = pqn->data;
    pqn->data  = data;
    pq->root  = joinPQNODE(pq->comp, pqn->left, pqn->right); 
    pqn->left = pqn->right = NULL; 
    pq->root  = pushPQNODE(pq->comp, pqn, pq->root); 
    return prevdata;
    }

/* freePQUEUE : DESTROY AND FREE pq.
*/
void freePQUEUE(PQUEUE *pq){
    while(popPQUEUE(pq));
    free(pq); 
    return;
    }

/* freePQUEUEptr : DESTROY AND FREE pq, USING ff to free members.
*/
void freePQUEUEptr(PQUEUE *pq, FREEFUNC ff){
    register void *v;
    while((v = popPQUEUE(pq)))
        ff(v);
    free(pq); 
    return;
    }

/* bfpPQUEUErecur : RECURSIVE FUNCTION FOR bestFirstPopPQUEUE.
  
*/
static void bfpPQUEUErecur(PQUEUE *pq, PQPOPFUNC func, 
                           int count, int total, void *v){
   register void *data = popPQUEUE(pq);
   if(!data)
       return;
   bfpPQUEUErecur(pq, func, count+1, total, v);
   func(data, total-count, total, v);
   return;       
   }

/* bestFirstPopPQUEUE : APPLY func TO EVERY ELEMENT IN THE pq, 
                        IN REVERSE ORDER (HIGHEST PRIORITY FIRST).
                        THIS DESTROYS THE QUEUE.
*/
void bestFirstPopPQUEUE(PQUEUE *pq, PQPOPFUNC func, void *v){
    bfpPQUEUErecur(pq, func, 0, pq->total, v);
    free(pq);
    return;
    }

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

#include <stdio.h>
#include "pqueue.h"

typedef struct { char *word; int pri; } DATA;

#ifdef NOTNOW
static BOOLEAN lowcompare(void *low, void *high){ 
    register DATA *da = low, *db = high;
    return (da->pri < db->pri);
    }
#endif /* NOTNOW */

static BOOLEAN highcompare(void *low, void *high){ 
    register DATA *da = low, *db = high;
    return (da->pri > db->pri);
    }

static void printData(void *data, int count, int total, void *info){ 
    DATA *d = data;
    printf("[%d/%d] PRI:[%d] WORD[%s]\n", 
           count, total, d->pri, d->word); 
    return;
    }

int main(){
    DATA dat[9] = {
         { "jumps", 5 }, { "dog",   9 }, { "brown", 3 },  
         { "over",  6 }, { "fox",   4 }, { "quick", 2 },
         { "lazy",  8 }, { "the",   1 }, { "a",     7 } 
    };
    PQUEUE *pqa = newPQUEUE(2, highcompare);
    PQUEUE *pqb = newPQUEUE(2, highcompare);
    int i; 
    for(i = 0; i < 5; i++){
        pushPQUEUE(pqa, &dat[i]);
        printf("INSERT (A) [%d] [%s]\n", dat[i].pri, dat[i].word);
        }
    while(i < 9){ 
        pushPQUEUE(pqb, &dat[i]);
        printf("INSERT (B) [%d] [%s]\n", dat[i].pri, dat[i].word);
        i++;
        }
    joinPQUEUE(pqa, pqb);
    bestFirstPopPQUEUE(pqa, printData, (void*)0); 
    return 0;  
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_PQUEUE_C */

