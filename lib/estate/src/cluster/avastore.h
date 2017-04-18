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

/* Code for out-of-memory storage, compression and sorting
   of the results of an all vs all comparison.

   Guy St.C. Slater.. December 1998.
*/
#ifndef INCLUDED_AVASTORE_H
#define INCLUDED_AVASTORE_H

#include <stdio.h>
#include "../general/common.h"
#include "../struct/pqueue.h"

typedef struct {
    int seqid;
    int score;
    } ALLVSALLWRESULT;

typedef struct {
    int position;
    int score;
    int fetch;
    } ALLVSALLHEAD;

typedef struct {
           FILE *fp;           /* ALLVSALL DATA FILE               */
            int  membertotal;  /* TOTAL NUMBER OF MEMBERS          */
            int  memberalloc;  /* MEMBERS ALLOCATED (READ)         */
            int  currentquery; /* CURRENT QUERY ID                 */
            int  currentcount; /* SCORES FOR CURRENT QUERY (WRITE) */
            int  currentscore; /* CURRENT SCORE (READ)             */
            int  scoretotal;   /* TOTAL NUMBER OF SCORES           */
ALLVSALLWRESULT *result;       /* FOR CURRENT QUERY (WRITE)        */
            int *buffer;       /* I/O BUFFER                       */
            int  bufferalloc;  /* I/O BUFFER ALLOCATED             */
            int  bufferused;   /* I/O BUFFER USED                  */
            int  bufferpos;    /* I/O BUFFER POSITION (READ)       */
   ALLVSALLHEAD *head;         /* HEAD OF PQUEUE DATA              */
         PQUEUE *pq;           /* PQUEUE OF SCORE HEADS (READ)     */
} ALLVSALL;

ALLVSALL *     newALLVSALL(char *path);
ALLVSALL *    openALLVSALL(char *path);
    void      freeALLVSALL(ALLVSALL *ava);
    void      infoALLVSALL(ALLVSALL *ava, FILE *fp);
    void  setqueryALLVSALL(ALLVSALL *ava, int query);
    void    finishALLVSALL(ALLVSALL *ava);
    void addresultALLVSALL(ALLVSALL *ava, int seqid, int score);
     int   viewTopALLVSALL(ALLVSALL *ava);

#define MaddresultALLVSALL(ava, id, sc)              \
    (ava)->result[(ava)->currentcount].seqid = (id), \
    (ava)->result[(ava)->currentcount++].score = (sc)

typedef struct {
    int from;
    int to;
    int score;
} ALLVSALLRESULT;

int nextALLVSALL(ALLVSALL *ava, ALLVSALLRESULT *avaresult);
int nextALLVSALLpair(ALLVSALL *a, ALLVSALL *b,
                     ALLVSALLRESULT *avaresult);

#endif /* INCLUDED_AVASTORE_H */

