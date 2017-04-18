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

#ifndef INCLUDED_METRIC_C
#define INCLUDED_METRIC_C

#include <stdio.h> 
#include <stdlib.h> 

#include "../general/error.h"
#include "../metric/metric.h"
#include "../dynamic/editdist.h"
#include "../struct/fsm.h"

ALLOW_ERROR_MESSAGES;

/* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */
/* FUNCITONS FOR EDIT DISTANCE METRIC                          */

static void *prepEDITdist(char *qy, METRICPARAM *param){
    return newEDITDIST(qy, strlen(qy));
    }

static int cacheEDITdist(void *qydat, char *seq, int len){
    return calcEDITDIST((EDITDIST*)qydat, seq, len);
    }

static int distEDITdist(void *qydat, long pos, FILE *fp){
    return calcEDITDISTonline((EDITDIST*)qydat, fp, pos);
    }

static void freeEDITdist(void *qydat){
    freeEDITDIST((EDITDIST*)qydat);
    return;     
    }

/* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

/* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */
/* FUNCTIONS FOR WORD BASED DISTANCE METRIC                    */

typedef struct {
    FSMHIST *hist;
       char *seq;
        int  len;
        int  mailboxctr;
        int  wordlen;
        int  qywords;
    } METRIC_WORDDIST_DATA;

static void *prepWORDdist(char *qy, METRICPARAM *param){
    register METRIC_WORDDIST_DATA *t = NEW(METRIC_WORDDIST_DATA);
    t->hist = makeFSMHIST(wordFSM((unsigned char*)qy,
                                 param->wd.wordlength));
    t->seq = qy;
    t->len = strlen(qy);
    t->wordlen = param->wd.wordlength;
    t->mailboxctr = 0;
    t->qywords = t->len-t->wordlen+1;
    return t;
    }

static int cacheWORDdist(void *qydat, char *seq, int len){
    register METRIC_WORDDIST_DATA *t = qydat;
    register int dist = t->qywords+(len-t->wordlen+1)
                      - (pairsFSMHIST(t->hist, (unsigned char*)seq,
                        t->mailboxctr)<<1); 
    t->mailboxctr++;
    return dist;
    }

static int distWORDdist(void *qydat, long pos, FILE *fp){
    register METRIC_WORDDIST_DATA *t = qydat;
    register int ch, dblen = 0, pairs = 0;
    register unsigned char c;
    register unsigned char *index = t->hist->fsm->index;
    register FSMHISTNODE *d;
    register FSMNODE *n = t->hist->fsm->root;
    fseek(fp, pos, SEEK_SET);
    t->mailboxctr++;
    do {
        switch(ch = getc(fp)){
            case 'A': case 'B': case 'C': case 'D': case 'E':
            case 'F': case 'G': case 'H': case 'I': case 'J':
            case 'K': case 'L': case 'M': case 'N': case 'O':
            case 'P': case 'Q': case 'R': case 'S': case 'T':
            case 'U': case 'V': case 'W': case 'X': case 'Y':
            case 'Z': /* VALID CHARACTER */
            case 'a': case 'b': case 'c': case 'd': case 'e':
            case 'f': case 'g': case 'h': case 'i': case 'j':
            case 'k': case 'l': case 'm': case 'n': case 'o':
            case 'p': case 'q': case 'r': case 's': case 't':
            case 'u': case 'v': case 'w': case 'x': case 'y':
            case 'z': /* VALID CHARACTER */
                c = index[ch];
                if((d = n[c].data.n)){
                    if(d->id == t->mailboxctr){ /* IF ALREADY VISITED */
                        if(d->match++ < 0)
                            pairs++;
                    } else {
                        d->id = t->mailboxctr; /* MARK AS VISITED */ 
                        d->match = 1-d->count; /* SUBSTACT FOR 1ST */
                        pairs++; 
                        }
                    }
                n = n[c].next; /* CHANGE STATE */
                dblen++;
                break;
            case '>': case EOF:
                 return (t->qywords+(dblen-t->wordlen+1))-(pairs<<1);
            }
    } while(TRUE);
    /* WILL ONLY RETURN AFTER "\n>" OR EOF */
    }

static void freeWORDdist(void *qydat){
    register METRIC_WORDDIST_DATA *t = qydat;
    freeFSMHIST(t->hist);
    free(t);
    return; 
    }

/* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

METRIC *newMETRIC(int type, METRICPARAM *param){
    METRIC *m = NEW(METRIC);
    switch(type){
        case METRIC_UNKNOWN:
            errmsg(ERROR_FATAL, "Unknown metric type\n");
            break;
        case METRIC_EDITDIST:
            m->prep  =  prepEDITdist;  
            m->dist  =  distEDITdist;  
            m->cache = cacheEDITdist;  
            m->free  =  freeEDITdist;  
            break;
        case METRIC_WORDDIST:
            m->prep  =  prepWORDdist;  
            m->dist  =  distWORDdist;  
            m->cache = cacheWORDdist;  
            m->free  =  freeWORDdist;  
            break;
        default:
            errmsg(ERROR_FATAL, "Invalid metric id [%d]");
            break;
        }
    m->type  = type;
    m->param = param;
    return m; 
    }

METRICPARAM *newMETRICPARAMwd(int wordlength){
    METRICPARAM *m = NEW(METRICPARAM);
    m->wd.wordlength = wordlength;
    return m;
    }

void freeMETRIC(METRIC *m){
    if(m->param)
         free(m->param);
    free(m);
    return;
    }

METRICPARAM *readMETRICPARAM(FILE *fp, char type){
    register int wordlen; 
    char line[100];
    if(type == METRIC_WORDDIST){
        fgets(line, 100, fp);
        wordlen = atoi(line);
        errmsg(ERROR_INFO, "Word length=%d", wordlen);
        return newMETRICPARAMwd(wordlen);
        }
    return NULL;
    }

void writeMETRICPARAM(FILE *fp, METRIC *metric){
    if(metric->type == METRIC_WORDDIST)
        fprintf(fp, "%d", metric->param->wd.wordlength);
    return;
    }

void writeMETRIC(FILE *fp, METRIC *metric){
    fprintf(fp, "%d\n", metric->type);
    writeMETRICPARAM(fp, metric);
    return;
    }

METRIC *readMETRIC(FILE *fp){
    char line[100];
    char *names[METRIC_TOTAL] = 
               {"Unknown", "Edit Distance", "Word Distance"};
    int type;
    fgets(line, 100, fp);
    type = atoi(line); 
    errmsg(ERROR_INFO, "Metric of type: [%s]\n", names[type]);
    return newMETRIC(type, readMETRICPARAM(fp, type));
    }

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

int main(){
    METRIC *a;
    a = newMETRIC(METRIC_EDITDIST, newMETRICPARAMwd(9));
    freeMETRIC(a);
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_METRIC_C */

