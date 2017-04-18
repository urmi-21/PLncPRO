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

#ifndef INCLUDED_AVASTORE_C
#define INCLUDED_AVASTORE_C

#include <stdlib.h>
#include "../general/common.h"
#include "../general/error.h"
#include "avastore.h"

ALLOW_ERROR_MESSAGES;

#define ALLVSALL_BUFFERCHUNKSIZE 100
#define ALLVSALL_MEMBERCHUNKSIZE 100

static void newALLVSALLcommon(ALLVSALL *ava){
    ava->bufferalloc = ALLVSALL_BUFFERCHUNKSIZE;
    ava->buffer = malloc(sizeof(int)*ALLVSALL_BUFFERCHUNKSIZE);
    return;
    }

ALLVSALL *newALLVSALL(char *path){
    register ALLVSALL *ava = BLANK(ALLVSALL);
    int space = 0;
    newALLVSALLcommon(ava);
    ava->fp = fopen(path, "wb");
    if(!ava->fp)
        errmsg(ERROR_FATAL, "Could not write output to [%s]", path);
    ava->result = malloc(sizeof(ALLVSALLWRESULT)
                         *ALLVSALL_MEMBERCHUNKSIZE);
    ava->head = calloc(ALLVSALL_MEMBERCHUNKSIZE, sizeof(ALLVSALLHEAD));
    /* LEAVE SPACE FOR MEMBERTOTAL AND SCORETOTAL */
    fwrite(&space, sizeof(int), 1, ava->fp);
    fwrite(&space, sizeof(int), 1, ava->fp);
    return ava;
    }

static BOOLEAN avapqcomp(void *low, void *high){
    return (((ALLVSALLHEAD*) low)->score >
            ((ALLVSALLHEAD*)high)->score);
    }

ALLVSALL *openALLVSALL(char *path){
    register ALLVSALL *ava = BLANK(ALLVSALL);
    int data[3];
    register int i;
    ava->fp = fopen(path, "rb");
    if(!ava->fp)
        errmsg(ERROR_FATAL, "Could not open [%s]", path);
    /* READ MEMBERTOTAL AND SCORE TOTAL */
    fread(&ava->membertotal, sizeof(int), 1, ava->fp);
    fread(&ava->scoretotal,  sizeof(int), 1, ava->fp);
    newALLVSALLcommon(ava);
    ava->head = malloc(sizeof(ALLVSALLHEAD)*ava->membertotal);
    /* READ THE ROOT INFO */
    ava->pq = newPQUEUE(0, avapqcomp);
    /* SEEK THE END - END DATA POSITION */
    fseek(ava->fp, -(ava->membertotal*3*sizeof(int)), SEEK_END);
    for(i = 0; i < ava->membertotal; i++){
        fread(data, sizeof(int), 3, ava->fp);
        ava->head[i].position = data[0];
        ava->head[i].score    = data[1];
        ava->head[i].fetch    = data[2];
        if(data[2]) /* IF ANY TO FETCH  */
            pushPQUEUE(ava->pq, &ava->head[i]);
        }
    return ava;
    }

int viewTopALLVSALL(ALLVSALL *ava){
    register ALLVSALLHEAD *h = viewTopPQUEUE(ava->pq);
    if(h)
       return h->score;
    return 0;
    }

void freeALLVSALL(ALLVSALL *ava){
    if(ava->result) /* ONLY FOR WRITE */
        free(ava->result);
    if(ava->pq) /* ONLY FOR READ */
        freePQUEUE(ava->pq);
    free(ava->head);
    fclose(ava->fp);
    free(ava);
    return;
    }

void infoALLVSALL(ALLVSALL *ava, FILE *fp){
    fprintf(fp, "AllvsALL members: %d\n"
                "          scores: %d\n",
                ava->membertotal, ava->scoretotal);
    return;
    }

/* compallvsallresult: COMPARISON FUNCTION FOR SORTING RESULTS.
*/
static int compallvsallresult(void *a, void *b){
    register ALLVSALLWRESULT *avara = a, *avarb = b;
    return (avara->score == avarb->score)
          ?(avara->seqid - avarb->seqid)
          :(avarb->score - avara->score);
    }

/* writesectiontobufferALLVSALL:
    CHANGE TO USE COMPRESSED BUFFER FORMAT HERE.
*/
static int writesectiontobufferALLVSALL(ALLVSALL *ava, int start){
    register int i, t, prev = ava->result[start].seqid;
    ava->buffer[0] = prev;
    ava->bufferused = 1;
    for(i = start+1; i < ava->currentcount; i++){
        if(ava->result[i].score != ava->result[start].score)
            break; /* STOP IF AT SECTION END */
        if(ava->bufferused >= ava->bufferalloc){
            ava->bufferalloc += ALLVSALL_BUFFERCHUNKSIZE;
            ava->buffer = realloc(ava->buffer,
                                  sizeof(int)*ava->bufferalloc);
            }
        t = ava->result[i].seqid;
        ava->buffer[ava->bufferused++] = t-prev;
        prev = t;
        }
    return i-start;
    }

/* writecurrentALLVSALL: WRITE RESULTS FOR CURRENT QUERY.
*/
static void writecurrentALLVSALL(ALLVSALL *ava){
    int ctr, fetch, space = 0;
    if(ava->currentcount > 1) /* SORT IF MORE THAN ONE RESULT */
        qsort(ava->result, ava->currentcount, sizeof(ALLVSALLWRESULT),
            (int(*)(const void *, const void *))compallvsallresult);
    /* WRITE FIRST SECTION TO BUFFER */
    ctr = writesectiontobufferALLVSALL(ava, 0);
    /* WRITE ROOT INFO TO ava->head */
    ava->head[ava->currentquery].position = ftell(ava->fp);
    ava->head[ava->currentquery].score = ava->result[0].score;
    ava->head[ava->currentquery].fetch = ctr;
    while(ctr != ava->currentcount){
        /* WRITE BUFFER TO FILE */
        fwrite(ava->buffer, sizeof(int), ava->bufferused, ava->fp);
        /* WRITE NEXT SECTION TO BUFFER */
        fetch = writesectiontobufferALLVSALL(ava, ctr);
        /* WRITE NEXT SECTION ROOT INFO TO FILE */
        fwrite(&ava->result[ctr].score, sizeof(int), 1, ava->fp);
        fwrite(&fetch, sizeof(int), 1, ava->fp);
        ctr += fetch;
        }
    fwrite(ava->buffer, sizeof(int), ava->bufferused, ava->fp);
    /* WRITE BLANK ROOT INFO */
    fwrite(&space, sizeof(int), 1, ava->fp);
    fwrite(&space, sizeof(int), 1, ava->fp);
    return;
    }

/* setqueryALLVSALL: SET CURRENT QUERY SEQUENCE.
*/
void setqueryALLVSALL(ALLVSALL *ava, int query){
    register int prev;
    if(ava->currentcount)
        writecurrentALLVSALL(ava);
    if(ava->membertotal <= query){
        ava->membertotal = query;
        /* printf("Chal with [%d] [%d]\n", */
                /* ava->membertotal, ava->memberalloc); */
        if(ava->membertotal >= ava->memberalloc){
            prev = ava->memberalloc;
            ava->memberalloc += ALLVSALL_MEMBERCHUNKSIZE;
            ava->head = realloc(ava->head,
                        (sizeof(ALLVSALLHEAD)*ava->memberalloc));
            memset(ava->head+prev, 0,
                        sizeof(ALLVSALLHEAD)*ALLVSALL_MEMBERCHUNKSIZE);
            ava->result = realloc(ava->result,
                        sizeof(ALLVSALLWRESULT)*ava->memberalloc);
            /* printf("Realloced to %d for result\n", */
                    /* ava->memberalloc); */
            }
    } else {
        if(ava->head[query].position)
            errmsg(ERROR_FATAL, "Results registered for %d twice",
                                 query);
        }
    ava->currentquery = query;
    ava->scoretotal += ava->currentcount;
    ava->currentcount = 0;
    return;
    }

/* finishALLVSALL: FINISH ADDITION OF SCORES TO ALLVSALL.
*/
void finishALLVSALL(ALLVSALL *ava){
    register int i;
    int data[3];
    if(ava->currentcount)
        writecurrentALLVSALL(ava);
    ava->scoretotal += ava->currentcount;
    for(i = 0; i < ava->membertotal; i++){
        data[0] = ava->head[i].position;
        data[1] = ava->head[i].score;
        data[2] = ava->head[i].fetch;
        fwrite(data, sizeof(int), 3, ava->fp);
        }
    rewind(ava->fp);
    fwrite(&ava->membertotal, sizeof(int), 1, ava->fp);
    fwrite(&ava->scoretotal,  sizeof(int), 1, ava->fp);
    return;
    }

/* addresultALLVSALL: ADD NEW RESULT
*/
void addresultALLVSALL(ALLVSALL *ava, int seqid, int score){
    /*
    register ALLVSALLWRESULT *avar = &ava->result[ava->currentcount++];
    */
    register ALLVSALLWRESULT *avar;
    /* printf("%p: make %d from %d %d\n", */
          /* ava, ava->currentcount+1, */
          /* ava->membertotal, ava->memberalloc); */
    avar = &ava->result[ava->currentcount++];
    avar->seqid = seqid;
    avar->score = score;
    return;
    }

/* fillbuffernextALLVSALL: READ TO FILL SCORE BUFFER
*/
static BOOLEAN fillbuffernextALLVSALL(ALLVSALL *ava){
    register ALLVSALLHEAD *avahead;
    avahead = popPQUEUE(ava->pq);   /* POP NODE */
    if(!avahead)
        return FALSE; /* END OF SEQUENCE */
    if((ava->bufferused+avahead->fetch) >= ava->bufferalloc){
        ava->bufferalloc += (avahead->fetch+ALLVSALL_BUFFERCHUNKSIZE);
        ava->buffer = realloc(ava->buffer,
                  sizeof(int)*(ava->bufferalloc));
        }
    ava->bufferpos = 0;
    ava->bufferused = avahead->fetch;
    fseek(ava->fp, avahead->position, SEEK_SET); /* USE POS */
    ava->currentscore = avahead->score; /* USE SCORE */
    ava->currentquery = avahead-ava->head;
    fread(ava->buffer, sizeof(int),  /* USE FETCH */
           avahead->fetch+2, ava->fp);
    if(ava->buffer[ava->bufferused+1]){ /* MORE TO FETCH */
        avahead->position = ftell(ava->fp); /* SET POS */
        avahead->score
            = ava->buffer[ava->bufferused]; /* SET SCORE */
        avahead->fetch
            = ava->buffer[ava->bufferused+1]; /* SET FETCH */
        pushPQUEUE(ava->pq, avahead); /* REPLACE NODE */
        }
    return TRUE;
    }

/* nextALLVSALL: FETCH NEXT HIGHEST SCORING RESULT.
                 RETURNS:
                       0 IF NO MORE SCORES.
                      -1 IF SCORE READ FROM DISK.
                       1 IF SCORE READ FROM BUFFER.
*/
int nextALLVSALL(ALLVSALL *ava, ALLVSALLRESULT *avaresult){
    register int retval = 1;
    if(ava->bufferpos == ava->bufferused){ /* IF BUFFER EMPTY */
        if(fillbuffernextALLVSALL(ava))
            retval = -1;
        else
            return 0;
        }
    avaresult->from  = ava->currentquery;
    avaresult->to    = ava->buffer[ava->bufferpos];
    avaresult->score = ava->currentscore;
    ava->buffer[++ava->bufferpos] += avaresult->to;
    return retval;
    }

/* nextALLVSALLpair: FETCH NEXT HIGHEST SCORING RESULT FROM a OR b.
                     RETURNS:
                         0 IF NO MORE SCORES
                         1 IF SCORE FROM A IS USED
                         2 IF SCORE FROM B IS USED
*/
int nextALLVSALLpair(ALLVSALL *a, ALLVSALL *b,
                     ALLVSALLRESULT *avaresult){
    register BOOLEAN amt = FALSE, bmt = FALSE;
    register ALLVSALL *ava;
    register int retval;
    if(a->bufferpos == a->bufferused) /* IF BUFFER EMPTY */
        if(!fillbuffernextALLVSALL(a))
            amt = TRUE;
    if(b->bufferpos == b->bufferused) /* IF BUFFER EMPTY */
        if(!fillbuffernextALLVSALL(b))
            bmt = TRUE;
    if(amt){
        if(bmt){ /* 0,0 */
            return 0; /* NOTHING LEFT */
        } else { /* 0,1 */
            ava = b;
            retval = 2;
            }
    } else {
        if(bmt){ /* 1,0 */
            ava = a;
            retval = 1;
        } else { /* 1,1 */
            if(a->currentscore > b->currentscore){
                ava = a;
                retval = 1;
            } else {
                ava = b;
                retval = 2;
                }
            }
        }
    avaresult->from  = ava->currentquery;
    avaresult->to    = ava->buffer[ava->bufferpos];
    avaresult->score = ava->currentscore;
    ava->buffer[++ava->bufferpos] += avaresult->to;
    return retval;
    }

#ifdef TEST_THIS_MODULE

BOOLEAN filtertestfunc(int id, void *info){
    return (id <= 5)?FALSE:TRUE;
    }

int main(){
    register char *path = "test.ava";
    register ALLVSALL *ava = newALLVSALL(path);
    ALLVSALLRESULT avaresult;

    setqueryALLVSALL(ava, 9);
    addresultALLVSALL(ava, 7, 11);
    addresultALLVSALL(ava, 3, 10);
    addresultALLVSALL(ava, 5, 10);
    addresultALLVSALL(ava, 8, 12);
    addresultALLVSALL(ava, 2, 12);
    addresultALLVSALL(ava, 6, 10);

    setqueryALLVSALL(ava, 6);
    addresultALLVSALL(ava, 2, 9);
    addresultALLVSALL(ava, 6,11);
    addresultALLVSALL(ava, 1, 9);
    addresultALLVSALL(ava, 3, 8);

    setqueryALLVSALL(ava, 5);
    addresultALLVSALL(ava, 3, 8);
    addresultALLVSALL(ava, 2,14);
    addresultALLVSALL(ava, 7, 3);
    addresultALLVSALL(ava, 9,10);
    addresultALLVSALL(ava, 8, 3);

    finishALLVSALL(ava);
    freeALLVSALL(ava);

    errmsg(ERROR_INFO, "Finished write");
    ava = openALLVSALL(path, filtertestfunc, NULL);
    while(nextALLVSALL(ava, &avaresult)){
        printf("Result %d vs %d score=%d\n",
               avaresult.from, avaresult.to, avaresult.score);
        }
    freeALLVSALL(ava);
    remove(path);
    return 0;
    }

#endif /* TEST_THIS_MODULE */

#endif /* INCLUDED_AVASTORE_C */

/*
Format Definition:
-----------------
[numseqs,resultstotal]
atrootposition(rootfetch(seqid,seqiddiff,...)
              ([nextscore,nextfetch],...))
foreachseq([rootposition,rootscore,rootfetch],...)
--
Todo:
    Make big/little-endian safe.
*/

