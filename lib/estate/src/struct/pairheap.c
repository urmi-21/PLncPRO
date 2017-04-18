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

#ifndef INCLUDED_PAIRHEAP_C
#define INCLUDED_PAIRHEAP_C

#include <stdio.h>
#include <stdlib.h>
#include "../general/error.h"

ALLOW_ERROR_MESSAGES;
#include "../struct/pairheap.h"

#define PAIRHEAP_COMBINE_CHUNKSIZE 128

static void blankfreefunc(void *v){return;}

PAIRHEAP *newPAIRHEAP(PQCOMPFUNC comp, 
                      FREEFUNC freedata, FREEFUNC freepri){
    register PAIRHEAP *ph = NEW(PAIRHEAP);
    ph->root = NULL;
    ph->total = 0;
    ph->comp = comp;
    ph->combinealloc = PAIRHEAP_COMBINE_CHUNKSIZE;
    ph->combine = malloc(sizeof(PAIRHEAPNODE)*ph->combinealloc);
    ph->freedata = freedata?freedata:blankfreefunc;
    ph->freepriority = freepri?freedata:blankfreefunc;
    return ph;
    }

static void freePAIRHEAPNODE(PAIRHEAP *ph, PAIRHEAPNODE *phn){
    if(phn){
        freePAIRHEAPNODE(ph, phn->left);
        freePAIRHEAPNODE(ph, phn->next);
        ph->freedata(phn->data);
        ph->freepriority(phn->priority);
        free(phn);
        }
    return;
    }

void freePAIRHEAP(PAIRHEAP *ph){
    freePAIRHEAPNODE(ph, ph->root);
    free(ph->combine);
    free(ph);
    return;
    }

static PAIRHEAPNODE *orderPAIRHEAPNODE(PAIRHEAP *ph, PAIRHEAPNODE *a, 
                                      PAIRHEAPNODE *b){
    register PAIRHEAPNODE *tmp;
    if(ph->comp(a->priority, b->priority)){ /* ADD b BEFORE a */
        a->next = b->next;
        if(a->next)
            a->next->prev = a;
        Swap(a, b, tmp);
    } else  /* ADD a BEFORE b */
        b->prev = a->prev;
    a->prev = b;
    a->next = b->left;
    if(a->next)
        a->next->prev = a;
    b->left = a;
    return b;
    }

PAIRHEAPNODE *pushPAIRHEAP(PAIRHEAP *ph, void *data, void *pri){
    register PAIRHEAPNODE *phn = NEW(PAIRHEAPNODE);
    phn->data      = data;
    phn->priority  = pri;
    phn->left = phn->next = phn->prev = NULL;
    if(ph->root)
        ph->root = orderPAIRHEAPNODE(ph, ph->root, phn);
    else
        ph->root = phn;  /* IS FIRST NODE */
    ph->total++;
    return phn;
    }

void raisePAIRHEAP(PAIRHEAP *ph, PAIRHEAPNODE *acc, void *pri){
    if(ph->comp(acc->priority, pri))
        errmsg(ERROR_FATAL, "raisePAIRHEAP with higher priority");
    ph->freepriority(acc->priority); 
    acc->priority = pri;
    if(ph->root == acc) /* IF THE ROOT NODE */
        return;
    if(acc->next)
        acc->next->prev = acc->prev;
    if(acc->prev->left == acc)
        acc->prev->left = acc->next;
    else
        acc->prev->next = acc->next;
    acc->next = NULL;
    ph->root = orderPAIRHEAPNODE(ph, ph->root, acc);
    return;
    }

static PAIRHEAPNODE *combinePAIRHEAPNODE(PAIRHEAP *ph, 
                                        PAIRHEAPNODE *phn){
    register int i, count;
    if(!phn->next) /* SINGLE MEMBER */
        return phn;
    for(count = 0; phn; count++){
        if(ph->combinealloc <= count){
            ph->combinealloc += PAIRHEAP_COMBINE_CHUNKSIZE;
            ph->combine = realloc(ph->combine, ph->combinealloc);
            }
        ph->combine[count] = phn;
        phn->prev->next = NULL; /* BREAK LINKS */
        phn = phn->next;
        }
    count--;
    for(i = 0; i < count; i+=2)
        ph->combine[i] = orderPAIRHEAPNODE(ph, 
                              ph->combine[i], ph->combine[i+1]);
    if(!(count&1)) /* PICK UP SPARE IF ODD COUNT */
        ph->combine[i-2] = orderPAIRHEAPNODE(ph, 
                                ph->combine[i-2], ph->combine[i]);
    for(i-=2; i>=2; i-=2) /* TIDY HERE */
        ph->combine[i-2] = orderPAIRHEAPNODE(ph, 
                             ph->combine[i-2], ph->combine[i]);
    return *ph->combine;
    }

void *popPAIRHEAP(PAIRHEAP *ph){
    register PAIRHEAPNODE *phn = ph->root;
    register void *data; 
    if(!ph->root)
        errmsg(ERROR_FATAL, "Cannot pop from an empty pairheap");
    ph->freepriority(ph->root->priority);
    data = ph->root->data;
    if(ph->root->left)
        ph->root = combinePAIRHEAPNODE(ph, ph->root->left);
    free(phn);
    if(!--ph->total) /* SET ROOT AS NULL IF PAIRHEAP NOW EMPTY */
        ph->root = NULL;
    return data;
    }

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

#include <stdio.h>

#define MaxSize 500

typedef struct {
    char *word;
     int  priority1;
     int  priority2;
     int  priority3;
    } DATA;

#define DATA_TOTAL 9

static BOOLEAN complowint(void *low, void *high){ 
    return (low < high);
    }

int main(){
    register int i;
    register PAIRHEAP *ph = newPAIRHEAP(complowint, NULL, NULL);
    PAIRHEAPNODE *acc[DATA_TOTAL];
    DATA data[DATA_TOTAL] = {
        { "jumps", 39, 38, 32}, { "dog",   72, 61, 60},
        { "brown", 99, 40,  8}, { "over",  68, 35, 33},
        { "fox",   20, 19, 17}, { "quick", 99, 19,  3},
        { "lazy",  99, 67, 51}, { "the",   3,   2,  1},
        { "a",     82, 41, 34}
    };
    printf("Add words with initial priorities\n");
    for(i = 0; i < DATA_TOTAL; i++){
        printf("Adding [%s]\n", data[i].word);
        acc[i] = pushPAIRHEAP(ph, 
                     (void*)data[i].word, (void*)data[i].priority1);
        }
    printf("Set intermediate priorities\n");
    for(i = 0; i < DATA_TOTAL; i++){
        printf("Raising [%s] from %d to %d\n", 
              data[i].word, data[i].priority1, data[i].priority2);
        raisePAIRHEAP(ph, acc[i], (void*)data[i].priority2);
        }
    printf("Set final priorities\n");
    for(i = 0; i < DATA_TOTAL; i++){
        printf("Raising [%s] from %d to %d\n", 
              data[i].word, data[i].priority2, data[i].priority3);
        raisePAIRHEAP(ph, acc[i], (void*)data[i].priority3);
        }

    printf("Popping data\n");
    for(i = 0; i < DATA_TOTAL; i++){
        printf("Word: [%s]\n", (char*)popPAIRHEAP(ph));
        }
    freePAIRHEAP(ph);
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_PAIRHEAP_C */

