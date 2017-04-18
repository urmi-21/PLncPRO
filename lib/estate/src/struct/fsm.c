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

#ifndef INCLUDED_FSM_C
#define INCLUDED_FSM_C

#include "fsm.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>   /* FOR MEMSET */

#define FSM_MINIMUM_NODES 16

/* allocFSMNODE : ALLOCATE AT LEAST FSM_MINIMUM_NODES FOR f,
                  RETURNS THE DISTANCE WHICH f->root HAS MOVED.
*/
static long allocFSMNODE(FSM *f, long extra){
    register FSMNODE *rootprev = f->root;
    register long diff;
    register size_t fresh = f->width*Max(extra, FSM_MINIMUM_NODES);
    f->root = (FSMNODE*)realloc(f->root,
                      sizeof(FSMNODE)*(f->alloc+fresh));
    memset(f->root+f->alloc, 0, sizeof(FSMNODE)*fresh);
    if((diff = (f->root - rootprev))){
        register FSMNODE *n, *end;
        for(n = f->root, end = f->root+f->alloc; n < end; n++)
            if(n->next)
                n->next += diff;
        }
    f->alloc+=fresh;
    return diff;
    }

/* submitFSM : ADD WORD s OF length TO THE PRE-FSM TRIE, f.
   1: A maximum of one reallocation per insertion will be required.
*/
void submitFSM(FSM *f, unsigned char *s, long length){
    register FSMNODE **ptp, *n = f->root;
    register unsigned char *end = s+length-1;
    for(end = s+length-1; s < end; s++){
        if(!*(ptp = &n[f->index[*s]].next)){
            if(f->count+f->width >= f->alloc){
                for(n += allocFSMNODE(f, end-s); s < end; s++){
                    if(!*(ptp = &n[f->index[*s]].next))
                        *ptp = &f->root[f->count+=f->width];
                    n = *ptp;
                    }
                n[f->index[*s]].data.l++;
                return;
                }
            *ptp = &f->root[f->count+=f->width];
            }
        n = *ptp;
        }
    n[f->index[*s]].data.l++;
    return;
    }

void submitFSMptr(FSM *f, unsigned char *s, long length, void *v){
    register FSMNODE **ptp, *n = f->root;
    register unsigned char *end = s+length-1;
    for(end = s+length-1; s < end; s++){
        if(!*(ptp = &n[f->index[*s]].next)){
            if(f->count+f->width >= f->alloc){
                for(n += allocFSMNODE(f, end-s); s < end; s++){
                    if(!*(ptp = &n[f->index[*s]].next))
                        *ptp = &f->root[f->count+=f->width];
                    n = *ptp;
                    }
                n[f->index[*s]].data.n = v;
                return;
                }
            *ptp = &f->root[f->count+=f->width];
            }
        n = *ptp;
        }
    n[f->index[*s]].data.n = v;
    return;
    }

void submitFSMptrJOIN(FSM *f, unsigned char *s, long length,
                      void *v, JOINFUNC combine){
    register FSMNODE **ptp, *n = f->root;
    register unsigned char *end = s+length-1;
    for(end = s+length-1; s < end; s++){
        if(!*(ptp = &n[f->index[*s]].next)){
            if(f->count+f->width >= f->alloc){
                for(n += allocFSMNODE(f, end-s); s < end; s++){
                    if(!*(ptp = &n[f->index[*s]].next))
                        *ptp = &f->root[f->count+=f->width];
                    n = *ptp;
                    }
                if(n[f->index[*s]].data.n)
                    n[f->index[*s]].data.n
                        = combine(n[f->index[*s]].data.n, v);
                else
                    n[f->index[*s]].data.n = v;
                return;
                }
            *ptp = &f->root[f->count+=f->width];
            }
        n = *ptp;
        }
    if( n[f->index[*s]].data.n)
        n[f->index[*s]].data.n
               = combine(n[f->index[*s]].data.n, v);
    else
        n[f->index[*s]].data.n = v;
    return;
    }

/* compileFSM : CONVERTS THE PRE-FSM TRIE TO THE FSM.
   1: All position zero nodes must have no score and point to root.
   2: All unused nodes must me made to point to the node which
      would be reached by insertion of their longest suffix
      (this will always be above the current node, and it's location
      is independent of any visits to root on this path).
   3: The scores from subsequences must also be inherited.
   4: This is achieved with a level-order trie traversal.
   5: A queue is used to remove traversal recursion and backtracking.
   6: No space is required for the queue, as it is makes temporary
      use of the position zero nodes during compilation.
   7: This algorithm is linear with the number of states.
*/
void compileFSM(FSM *f){
    register int i;
    FSMNODE a;
    register FSMNODE *suffix, *prev;
    register FSMNODE *out = &a, *in = out;
    out->next = in;              /* INITIALISE QUEUE */
    f->root->next = f->root;
    for(i = 1; i < f->width; i++)
        if(f->root[i].next){
            in->data.n = f->root;
            in->next   = f->root[i].next;
            in         = f->root[i].next;
        } else
            f->root[i].next   = f->root;
    while(out != in){
        prev   = out;
        suffix = out->data.n;
        out    = out->next;
        for(i = 1; i < f->width; i++){
            out[i].data.l += suffix[i].data.l;  /* INHERIT */
            if(out[i].next){
                in->data.n = suffix[i].next;
                in->next   = out[i].next;
                in         = out[i].next;
            } else
                out[i].next = suffix[i].next;
            }
        prev->next = f->root;
        prev->data.n = NULL;
        }
    out->next = f->root;
    out->data.n = NULL;
    return;
    }

/* compileFSMptr : COMPILE FSM USING POINTERS AT THE NODES.
                   THUS ABANDON INHERITANCE, AND JUST REPORT
                   HIT FOR THE LONGER OF THE TWO WORDS.
                   ONLY THE INHERIT LINE IS REMOVED.
*/
void compileFSMptr(FSM *f){
    register int i;
    FSMNODE a;
    register FSMNODE *suffix, *prev;
    register FSMNODE *out = &a, *in = out;
    out->next = in;              /* INITIALISE QUEUE */
    f->root->next = f->root;
    for(i = 1; i < f->width; i++)
        if(f->root[i].next){
            in->data.n = f->root;
            in->next   = f->root[i].next;
            in         = f->root[i].next;
        } else
            f->root[i].next   = f->root;
    while(out != in){
        prev   = out;
        suffix = out->data.n;
        out    = out->next;
        for(i = 1; i < f->width; i++){
            if(out[i].next){
                in->data.n = suffix[i].next;
                in->next   = out[i].next;
                in         = out[i].next;
            } else
                out[i].next = suffix[i].next;
            }
        prev->next = f->root;
        prev->data.n = NULL;
        }
    out->next = f->root;
    out->data.n = NULL;
    return;
    }

/* initFSM : COMPLETE INDEX AND INITIAL FSM MEMORY ALLOCATION */
void initFSM(FSM *f){
    register unsigned char *p;
    for(p = f->index+256; p > f->index; p--)
        if(*p)
            *p = ++f->width;
    f->count = ++f->width; /* ADD ONE FOR NULL */
    f->alloc = f->width*FSM_MINIMUM_NODES;
    f->root = calloc(f->alloc, sizeof(FSMNODE));
    return;
    }

/* makeFSM : RETURNS AN FSM TO SEARCH FOR total words.
*/
FSM *makeFSM(char **words, long total){
    register unsigned char *p;
    register long i;
    FSM *f = BLANK(FSM);
    for(i = 0; i < total; i++) /* PREPARE INDEX */
        for(p = (unsigned char*)words[i]; *p; p++)
            f->index[*p] = 1;
    initFSM(f);
    for(i = 0; i < total; i++)
        submitFSM(f, (unsigned char*)words[i], strlen(words[i]));
    compileFSM(f);
    return f;
    }

/* makeFSMptr : RETURNS AN FSM TO SEARCH FOR total words.
*/
FSM *makeFSMptr(char **words, long total, void **ptr){
    register unsigned char *p;
    register long i;
    FSM *f = BLANK(FSM);
    for(i = 0; i < total; i++) /* PREPARE INDEX */
        for(p = (unsigned char*)words[i]; *p; p++)
            f->index[*p] = 1;
    initFSM(f);
    for(i = 0; i < total; i++)
        submitFSMptr(f, (unsigned char*)words[i],
                     strlen(words[i]), ptr[i]);
    compileFSMptr(f);
    return f;
    }

/* makeFSMptrJOIN : RETURNS AN FSM TO SEARCH FOR total words.
*/
FSM *makeFSMptrJOIN(char **words, long total, void **ptr,
                JOINFUNC combine){
    register unsigned char *p;
    register long i;
    FSM *f = BLANK(FSM);
    for(i = 0; i < total; i++) /* PREPARE INDEX */
        for(p = (unsigned char*)words[i]; *p; p++)
            f->index[*p] = 1;
    initFSM(f);
    for(i = 0; i < total; i++)
        submitFSMptrJOIN(f, (unsigned char*)words[i], strlen(words[i]),
                         ptr[i], combine);
    compileFSMptr(f);
    return f;
    }

/* newFSM : PREPARE FSM FOR ALPHABET alpha.
*/
FSM *newFSM(char *alpha){
    register unsigned char *p;
    FSM *f = BLANK(FSM);
    for(p = (unsigned char*)alpha; *p; p++)
        f->index[*p] = 1;
    initFSM(f);
    return f;
    }

/* wordFSM : CONSTRUCT A FSM WITH ALL WORDS OF LENGTH wlen OF str.
*/
FSM *wordFSM(unsigned char *str, long wlen){
    register unsigned char *p;
    register long len;
    FSM *f = BLANK(FSM);
    for(p = str; *p; p++)  /* PREPARE INDEX */
        f->index[*p] = 1;
    len = p-str;
    initFSM(f);
    for(p = str+len-wlen; str <= p; p--)
        submitFSM(f, p, wlen);
    compileFSM(f);
    return f;
    }

/* countFSM : NAVIGATE FSM WITH data, RETURN THE NUMBER OF WORDS
              IN THE FSM HIT BY DATA */
long countFSM(FSM *f, unsigned char *str){
    register FSMNODE *n = f->root;
    register long count = 0;
    register unsigned char c;
    do {
        count += n[c = f->index[*str]].data.l;
        n = n[c].next;
    } while(*++str);
    return count;
    }

#ifdef DEBUG
/* countFSM : NAVIGATE FSM WITH data, RETURN THE NUMBER OF WORDS
              IN THE FSM HIT BY DATA */
long countFSM_DEBUG(FSM *f, unsigned char *str){
    register FSMNODE *n = f->root;
    register long count = 0;
    register unsigned char c;
    printf("search [%c", *str);
    do {
        printf("%c", *str);
        c = f->index[*str];
        if(n[c].data.l)
            printf("+");
        count += n[c].data.l;
        n = n[c].next;
    } while(*++str);
    printf("]\n\n");
    return count;
    }
#endif /* DEBUG */

/* infoFSM : PRINTS OUT INFORMATION ABOUT FSM f.
*/
void infoFSM(FSM *f){
    printf("FSM information:\n"
           " Alphabet size   = %d\n"
           " Size of node    = %d\n"
           " Nodes allocated = %d\n"
           " Nodes used      = %d\n"
           " Nodes unused    = %d\n"
           " Bytes unused    = %d\n",
           f->width,
           sizeof(FSMNODE),
           f->alloc/f->width,
           f->count/f->width,
           f->alloc-f->count,
           (f->alloc-f->count)*sizeof(FSMNODE) );
    return;
    }

/* freeFSM : FREE ALL MEMORY USED BY f
*/
void freeFSM(FSM *f){
    free(f->root);
    free(f);
    return;
    }

/* FREE FSM INCLUDING STATE OBJECTS
*/
void freeFSMptr(FSM *f, FREEFUNC ff){
    register FSMNODE *n, *end;
    for(n = f->root, end = f->root+f->alloc; n < end; n++)
        if(n->data.n)
            ff(n->data.n);
    free(f->root);
    free(f);
    return;
    }



/* copyFSM : RETURNS A COPY OF f (DOES NOT COPY UNUSED NODES).
*/

FSM *copyFSM(FSM *f){
    register FSM *nf = NEW(FSM);
    register long diff;
    register size_t datasize;
    register FSMNODE *a, *z;
    memcpy(nf, f, sizeof(FSM));
    nf->alloc = nf->count;
    datasize = nf->alloc*sizeof(FSMNODE);
    nf->root = (FSMNODE*)malloc(datasize);
    memcpy(nf->root, f->root, datasize);
    diff = nf->root-f->root;
    for(a = nf->root, z = nf->root+nf->alloc; a < z; a++)
        if(a->next)
            a->next += diff;
    return nf;
    }

/* HISTOGRAM ROUTINES */

/* makeFSMHIST: MAKE THE HISTOGRAM OF WORDS (h) FOR FSM (f).
*/
FSMHIST *makeFSMHIST(FSM *f){
    register FSMHIST *h = NEW(FSMHIST);
    register FSMNODE *a, *z;
    register FSMHISTNODE *hp;
    register long wc;
    h->fsm = f;
    wc = 1; /* ADD ONE FOR NULL */
    for(a = f->root, z = a+f->count+f->width; a < z; a++)
        if(a->data.l)
            wc++;
    h->hist = (FSMHISTNODE*)malloc(wc*sizeof(FSMHISTNODE));
    for(a = f->root, hp = h->hist; a < z; a++)
        if(a->data.l){
            hp->id = 0; /* MARK ALL AS UNVISITED */
            hp->count = a->data.l;
            a->data.n = hp++;
            }
    return h;
    }

void freeFSMHIST(FSMHIST *h){
    freeFSM(h->fsm);
    free(h->hist);
    free(h);
    return;
    }

/* pairsFSMHIST : COUNTS NUMBER OF WORD PAIRS MADE FROM SEQS
*/
long pairsFSMHIST(FSMHIST *h, unsigned char *seq, long id){
    register long pairs = 0;
    register unsigned char c, *s = seq, *index = h->fsm->index;
    register FSMNODE *n = h->fsm->root;
    register FSMHISTNODE *d;
    do {
       c = index[*s];
       if((d = n[c].data.n)){    /* IF WORD OCCURS IN QUERY */
           if(d->id == id){ /* IF ALREADY VISITED THIS SEQUENCE */
               if(d->match++ < 0)
                   pairs++;
           } else {
               d->id = id;  /* MARK AS VISITED */
               d->match = 1-d->count; /* SUBTRACT FOR 1ST MATCH */
               pairs++;
               }
           }
       n = n[c].next;
    } while(*++s);
    return pairs;
    }

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

int main(){
    FSM *f, *t;
#define SET_JAPANCITY
#ifdef SET_JAPANCITY
#define TESTSETSIZE 46
    char *testset[TESTSETSIZE] = { /* JAPANESE CITIES */
        "AOMORI", "AKITA", "CHIBA", "FUKUSHIMA", "FUKUI", "GIFU",
        "HIROSHIMA", "HAKATA", "KAGOSHIMA", "KANAZAWA", "KUMAMOTO",
        "KYOTO", "KOCHI", "KOBE", "KOFO", "MATUYAMA", "MATUE",
        "MAEBASHI", "MIYAZAKI", "MITO", "MORIOKA", "NAGASAKI",
        "NAGANO", "NAGOYA", "NARA", "NIIGATA", "OKAYAMA", "OSAKA",
        "OOITA", "OTSU", "SHIZUOKA", "SAPPRO", "SAGA", "SENDAI",
        "TAKAMATU", "TOKUSHIMA", "TOKYO", "TOTTORI", "TOYAMA",
        "TSU", "USTUNOMIYA", "URAWA", "WAKAYAMA", "YAMAGUCHI",
        "YAMAGATA", "YOKOHAMA"}; /* FROM Aoe1989 */
    char *testsen = "testing TOKYO ...";
    f = makeFSM(testset, TESTSETSIZE);
#endif /* SET_JAPANCITY */

#ifdef SET_KEYWORDS
#define TESTSETSIZE 32
    char *testset[TESTSETSIZE] = { /* C KEYWORDS */
        "auto",     "break",    "case",     "char",
        "const",    "continue", "default",  "do",
        "double",   "else",     "enum",     "extern",
        "float",    "for",      "goto",     "if",
        "int",      "long",     "register", "return",
        "short",    "signed",   "sizeof",   "static",
        "struct",   "switch",   "typedef",  "union",
        "unsigned", "void",     "volatile", "while"  };

    char *testsen = "-unsigned--switchar---doublelse--sizeofor";
/*        TOTAL=9            2       1 1    1   1  1       1 1 */
    f = makeFSM(testset, TESTSETSIZE);
#endif /* SET_KEYWORDS */

    t = copyFSM(f);
    freeFSM(f);
    f = t;

#ifdef SET_SIMPLE_A
#define TESTSETSIZE 5
    char *testset[TESTSETSIZE] = {
        "start", "star", "tar", "tart", "art" };
    char *testsen = "-start-";
                  /*     23  : TOTAL=5 */
    f = makeFSM(testset, TESTSETSIZE);
#endif /* SET_SIMPLE_A */

#ifdef SET_SIMPLE_B
#define TESTSETSIZE 2
    char *testset[TESTSETSIZE] = { "start", "tar"};
    char *testsen = "-start-";
    f = makeFSM(testset, TESTSETSIZE);
#endif /* SET_SIMPLE_B */

#ifdef SET_WORDSET
#define TESTSETSIZE 3
    char *testset = "appleggplantomatorangetc";
    char *testsen = "----ggp-----mat---get---";
    f = wordFSM(testset, TESTSETSIZE);
#endif /* SET_WORDSET */
    infoFSM(f);
    printf("FSM count is [%ld]\n", countFSM(f, testsen) );
    freeFSM(f);
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */

#endif /* INCLUDED_FSM_C */

