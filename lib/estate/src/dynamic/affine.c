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

#ifndef INCLUDED_AFFINE_C
#define INCLUDED_AFFINE_C

#include "affine.h"
#include <stdlib.h>
#include <string.h> /* FOR memset()  */
#include <ctype.h>  /* FOR tolower() */
#include "../general/common.h"
#include "../general/error.h"

ALLOW_ERROR_MESSAGES;

/* IT'S PRETTY DIFFICULT TO WRITE THIS SORT OF CODE
   WITHOUT IT BECOMING EITHER VERBOSE OR SLOW.

   I'VE ABUSED THE PREPROCESSOR QUITE A LOT,
   SO USE PLENTY OF OPTIMISATION SO THE COMPILER CAN
   MAKE SOME SENSE OF IT.

   NEXT TIME, I'D SERIOUSLY CONSIDER USING DYNAMITE:
   http://www.sanger.ac.uk/Software/Dynamite/
*/

static void prepareAFFINE_CELL(AFFINE *affine){
    register int i;
    affine->cell = malloc(sizeof(AFFINE_CELL)*(affine->qylen));
    for(i = 0; i < affine->qylen; i++)
        affine->cell[i].matrix = &affine->submat->matrix
                   [global_submat_index[affine->qyseq[i]]][0];
    return;
    }

AFFINE *newAFFINE(unsigned char *qyseq,  int qylen,
                 int gapopen, int gapextend, SUBMAT *submat){
    register AFFINE *affine = NEW(AFFINE);
    affine->qyseq     = (unsigned char*)strdup((char*)qyseq);
    affine->qylen     = qylen;
    affine->submat    = submat;
    affine->gapopen   = gapopen;
    affine->gapextend = gapextend;
    prepareAFFINE_CELL(affine);
    if(gapopen >= 0)
        errmsg(ERROR_FATAL, "Gap open penalty should be negative");
    if(gapextend >= 0)
        errmsg(ERROR_FATAL, "Gap extend penalty should be negative");
    return affine;
    }

void freeAFFINE(AFFINE *affine){
    free(affine->qyseq);
    free(affine->cell);
    free(affine);
    return;
    }

#ifdef NOTNOW
static void printAFFINE(AFFINE *affine, char *msg, char row){
    register int i;
    printf("ROW:%s:", msg);
    for(i = 0; i < affine->qylen; i++)
        printf("[%3d:%3d]",
            (row=='a')?affine->cell[i].scorea:affine->cell[i].scoreb,
            affine->cell[i].extend);
    printf("\n");
    return;
    }
#endif /* NOTNOW */

#define AffineNullCommand (j=j)

#define ChallengeScore(target, candidate, qypos, dbpos) \
    Challenge(target, candidate);

#define AffineMainLoop(initextra, incrextra) \
    for(i = 0, (initextra); i < iz; i+=2, (incrextra))

#define ChallengeAffine(extend, candextend, score, candscore, \
                        qypos, dbpos, chext, chscore)         \
    (extend) += (candextend);                                 \
    Challenge ## chext (extend, candscore, qypos, dbpos);     \
    Challenge ## chscore (score, extend, qypos, dbpos);

#define ChallengeCell(srcrow, dstrow, qypos, dbpos,          \
                      chde, chds, chre, chrs)                \
/*DIAG*/    cell[(qypos)+1].score ## dstrow =                \
                     cell[(qypos)].score ## srcrow           \
                   + cell[(qypos)+1].matrix[submatpos];      \
            temp = cell[(qypos)].score ## srcrow + gapopen;  \
/*DOWN*/    ChallengeAffine(cell[(qypos)].extend, gapextend, \
                   cell[(qypos)].score ## dstrow, temp,      \
                   qypos, dbpos, chde, chds);                \
/*RIGHT*/   ChallengeAffine(rightextend, gapextend,          \
                   cell[(qypos)+1].score ## srcrow, temp,    \
                   qypos, dbpos, chre, chrs);

#define AffineMainInnerLoop(srcrow, dstrow, dbseqos, leftinit,     \
                            dbseq, initextra, loopextra,           \
                            chde, chds, chre, chrs)                \
        submatpos = global_submat_index[(dbseq)[i+(dbseqos)]];     \
        cell[0].score ## dstrow = (leftinit)                       \
                                + (cell[0].matrix[submatpos]);     \
        rightextend = gapopen + (leftinit);                        \
        Challenge ## chrs (cell[0].score ## srcrow, rightextend,   \
                           -1, i+(dbseqos)-1);                     \
        initextra ;                                                \
        for(j = 0; j < jz; j++){                                   \
            ChallengeCell(srcrow, dstrow, j, i+(dbseqos)-1,        \
                          chde, chds, chre, chrs);                 \
            loopextra ;                                            \
            }                                                      \
        ChallengeAffine(cell[j].extend, gapextend,                 \
                        cell[j].score ## dstrow,                   \
                        cell[j].score ## srcrow + gapopen,         \
                        j, i+(dbseqos)-1, chde, chds)

#define AffineLocalZero(score, row, pos)            \
    Challenge(cell[(pos)].score ## row, 0);         \
    else Challenge(score, cell[(pos)].score ## row)

#define AffineLastRow(row, leftinit, loopextra, chre, chrs)  \
    rightextend = gapopen + (leftinit);                      \
    Challenge ## chrs (cell[0].score ## row, rightextend,    \
                       -1, i+oddrow);                        \
    for(j = 0; j < jz; j++){                                 \
        ChallengeAffine(rightextend, gapextend,              \
                   cell[j+1].score ## row,                   \
                   cell[j].score ## row + gapopen,           \
                   j, i+oddrow, chre, chrs);                 \
        loopextra ;                                          \
        }

#define AffineInit(dyn, qystart, dbseq, dbseqlen,                   \
                   initextra, incrextra, topinit, chds)             \
    register int i, j, temp, submatpos, rightextend;                \
    register AFFINE_CELL *cell = (dyn)->cell+(qystart);             \
    register int iz = ((((dbseqlen)-1)|1)>>1)<<1,                   \
                 jz = (dyn)->qylen-1+(qystart);                     \
    register int gapopen = dyn->gapopen,                            \
                 gapextend = (dyn)->gapextend;                      \
    register int oddrow = (((dbseqlen)-1)&1);                       \
    submatpos = global_submat_index[(dbseq)[0]];                    \
    cell[0].scorea = cell[0].matrix[submatpos];                     \
    for(j = 0, (initextra); j < jz; j++, (incrextra)){              \
/*DIAG*/cell[j+1].scorea = (topinit) + cell[j+1].matrix[submatpos]; \
        cell[j].extend = gapopen+(topinit);                         \
        Challenge ## chds (cell[j].scorea, cell[j].extend, j, -1);  \
        }                                                           \
    cell[j].extend = gapopen + (topinit);                           \
    Challenge ## chds (cell[j].scorea, cell[j].extend, j, -1);

int globalAFFINEscore(AFFINE *affine,
                      unsigned char *dbseq, int dbseqlen){
    register int gec;
    AffineInit(affine, 0, dbseq, dbseqlen,
               (gec = 0), (gec += gapextend),
               gapopen+gec, Score);
    AffineMainLoop((gec = 0), (gec += (gapextend+gapextend))){
        AffineMainInnerLoop(a, b, 1, gapopen+gec, dbseq,
                            AffineNullCommand, AffineNullCommand,
                            Score, Score, Score, Score);
        AffineMainInnerLoop(b, a, 2, gapopen+gec+gapextend, dbseq,
                            AffineNullCommand, AffineNullCommand,
                            Score, Score, Score, Score);
        }
    if(oddrow){
        AffineMainInnerLoop(a, b, 1, gapopen+gec, dbseq,
                            AffineNullCommand, AffineNullCommand,
                            Score, Score, Score, Score);
        AffineLastRow(b, gapopen+gec, AffineNullCommand,
                            Score, Score);
        return cell[j].scoreb;
        }
    AffineLastRow(a, gapopen+gec, AffineNullCommand, Score, Score);
    return cell[j].scorea;
    }

int bestfitAFFINEscore(AFFINE *affine,
                       unsigned char *dbseq, int dbseqlen){
    register int gec;
    register int score;
    AffineInit(affine, 0, dbseq, dbseqlen,
               (gec = 0), (gec += gapextend),
               gapopen+gec, Score);
    score = cell[j].scorea;
    AffineMainLoop(AffineNullCommand, AffineNullCommand){
        AffineMainInnerLoop(a, b, 1, 0, dbseq,
                            AffineNullCommand, AffineNullCommand,
                            Score, Score, Score, Score);
        Challenge(score, cell[j].scorea);
        AffineMainInnerLoop(b, a, 2, 0, dbseq,
                            AffineNullCommand, AffineNullCommand,
                            Score, Score, Score, Score);
        Challenge(score, cell[j].scoreb);
        }
    if(oddrow){
        AffineMainInnerLoop(a, b, 1, 0, dbseq,
                            AffineNullCommand, AffineNullCommand,
                            Score, Score, Score, Score);
        Challenge(score, cell[j].scorea);
        AffineLastRow(b, 0, AffineNullCommand, Score, Score);
        Challenge(score, cell[j].scoreb);
        return score;
        }
    AffineLastRow(a, 0, AffineNullCommand, Score, Score);
    Challenge(score, cell[j].scorea);
    return score;
    }

int localAFFINEscore(AFFINE *affine,
                     unsigned char *dbseq, int dbseqlen){
    register int score = 0;
    AffineInit(affine, 0, dbseq, dbseqlen,
               AffineNullCommand, AffineNullCommand, 0, Score);
    AffineMainLoop(AffineNullCommand, AffineNullCommand){
        AffineMainInnerLoop(a, b, 1, 0, dbseq,
                            AffineLocalZero(score, a, 0),
                            AffineLocalZero(score, a, j+1),
                            Score, Score, Score, Score);
        AffineMainInnerLoop(b, a, 2, 0, dbseq,
                            AffineLocalZero(score, b, 0),
                            AffineLocalZero(score, b, j+1),
                            Score, Score, Score, Score);
        }
    if(oddrow){
        AffineMainInnerLoop(a, b, 1, 0, dbseq,
                            AffineLocalZero(score, a, 0),
                            AffineLocalZero(score, a, j+1),
                            Score, Score, Score, Score);
        AffineLastRow(b, 0, AffineLocalZero(score, b, j+1),
                            Score, Score);
        return score;
        }
    AffineLastRow(a, 0, AffineLocalZero(score, a, j+1), Score, Score);
    return score;
    }

int overlapAFFINEscore(AFFINE *affine,
                       unsigned char *dbseq, int dbseqlen){
    register int score = 0;
    AffineInit(affine, 0, dbseq, dbseqlen,
               AffineNullCommand, AffineNullCommand, 0, Score);
    AffineMainLoop(AffineNullCommand, AffineNullCommand){
        AffineMainInnerLoop(a, b, 1, 0, dbseq,
                            AffineNullCommand, AffineNullCommand,
                            Score, Score, Score, Score);
        Challenge(score, cell[j].scorea);
        AffineMainInnerLoop(b, a, 2, 0, dbseq,
                            AffineNullCommand, AffineNullCommand,
                            Score, Score, Score, Score);
        Challenge(score, cell[j].scoreb);
        }
    if(oddrow){
        AffineMainInnerLoop(a, b, 1, 0, dbseq,
                            AffineNullCommand, AffineNullCommand,
                            Score, Score, Score, Score);
        Challenge(score, cell[j].scorea);
        AffineLastRow(b, 0, Challenge(score, cell[j+1].scoreb),
                            Score, Score);
        return score;
        }
    AffineLastRow(a, 0, Challenge(score, cell[j+1].scorea),
                            Score, Score);
    return score;
    }

/* START OF ALIGNMENT GENERATION CODE */

static void **newMATRIX(int x, int y, size_t chunk){
    register int i, j;
    register char **m;
    register size_t s = sizeof(char*)*x, r = chunk*y, t = r*x;
    m = (char**)malloc( s + t );
    if(!m)
        errmsg(ERROR_FATAL, "Insufficient memory for matrix");
    *m = (char*)m + s;
    memset(*m, 0, t);
    for(i = 1, j = (y*chunk); i < x; i++, j+=(y*chunk))
        m[i] = *m+j;
    return (void**)m;
    }

#define AFFINE_ALIGN_END         (1<<0)
#define AFFINE_ALIGN_MATCH       (1<<1)
#define AFFINE_ALIGN_INSERT      (1<<2)
#define AFFINE_ALIGN_DELETE      (1<<3)
#define AFFINE_ALIGN_INSERT_OPEN (1<<4)
#define AFFINE_ALIGN_DELETE_OPEN (1<<5)

#define AFFINE_SET_ALIGN_AFFINE(qypos, dbpos, type) \
    matrix[(qypos)][(dbpos)] |= (type)

#define AFFINE_GET_ALIGN_AFFINE(qypos, dbpos, type) \
    matrix[(qypos)][(dbpos)] & (type)

#define AFFINE_UNSET_ALIGN_AFFINE(qypos, dbpos, type) \
    matrix[(qypos)][(dbpos)] &= (~(type))

#define ChallengeDownScore(target, candidate, qypos, dbpos) \
    if((target) <= (candidate)){                            \
        if((target) == (candidate)){                        \
            AFFINE_SET_ALIGN_AFFINE(qypos, (dbpos)+1,       \
                  AFFINE_ALIGN_MATCH|AFFINE_ALIGN_INSERT);  \
        } else {                                            \
            (target) = (candidate);                         \
            AFFINE_SET_ALIGN_AFFINE(qypos, (dbpos)+1,       \
                                      AFFINE_ALIGN_INSERT); \
            }                                               \
    } else {                                                \
        AFFINE_SET_ALIGN_AFFINE(qypos, (dbpos)+1,           \
                                AFFINE_ALIGN_MATCH);        \
        }

#define ChallengeDownExtend(target, candidate, qypos, dbpos) \
    if((target) <= (candidate)){                             \
        (target) = (candidate);                              \
        AFFINE_SET_ALIGN_AFFINE(qypos, (dbpos)+1,            \
                                AFFINE_ALIGN_INSERT_OPEN);   \
        }

#define ChallengeRightScore(target, candidate, qypos, dbpos) \
    if((target) <= (candidate)){                             \
        if((target) != (candidate)){                         \
            (target) = (candidate);                          \
            AFFINE_UNSET_ALIGN_AFFINE((qypos)+1, dbpos,      \
                  AFFINE_ALIGN_MATCH|AFFINE_ALIGN_INSERT);   \
            }                                                \
        AFFINE_SET_ALIGN_AFFINE((qypos)+1, dbpos,            \
                                AFFINE_ALIGN_DELETE);        \
        }

#define ChallengeRightExtend(target, candidate, qypos, dbpos) \
    if((target) <= (candidate)){                              \
        (target) = (candidate);                               \
        AFFINE_SET_ALIGN_AFFINE((qypos)+1, dbpos,             \
                  AFFINE_ALIGN_DELETE_OPEN);                  \
        }

static void reverseAFFINEtranscript(int *transcript, int length){
    register int *a, *z, swap;
    for(a = transcript, z = transcript+length-1; a < z; a++, z--)
        Swap(*a,*z,swap);
    return;
    }

static void ExtractBestTranscript(AFFINE *affine,
                                  AFFINEALIGN *affinealign,
                                  char **matrix){
    register int qypos = affinealign->qyalignlen-1,
                         dbpos = affinealign->dbalignlen-1;
    affinealign->transcript = malloc(sizeof(int)
                            *(affinealign->qyalignlen
                             +affinealign->dbalignlen));
    affinealign->transcriptlen = 0;
    while((qypos >= 0) && (dbpos >= 0)){
        if(AFFINE_GET_ALIGN_AFFINE(qypos, dbpos,
                                    AFFINE_ALIGN_MATCH)){
            affinealign->transcript[affinealign->transcriptlen++]
                     = AFFINE_ALIGN_MATCH;
            qypos--;
            dbpos--;
        } else if(AFFINE_GET_ALIGN_AFFINE(qypos, dbpos,
                                    AFFINE_ALIGN_INSERT)){
            do {
                affinealign->transcript[affinealign->transcriptlen++]
                               = AFFINE_ALIGN_INSERT;
                dbpos--;
                if(AFFINE_GET_ALIGN_AFFINE(qypos, dbpos+1,
                                     AFFINE_ALIGN_INSERT_OPEN))
                    break;
            } while(dbpos > 0);
        } else if(AFFINE_GET_ALIGN_AFFINE(qypos, dbpos,
                                    AFFINE_ALIGN_DELETE)){
            do {
                affinealign->transcript[affinealign->transcriptlen++]
                               = AFFINE_ALIGN_DELETE;
                qypos--;
                if(AFFINE_GET_ALIGN_AFFINE(qypos+1, dbpos,
                                     AFFINE_ALIGN_DELETE_OPEN))
                    break;
            } while(qypos > 0);
            }
        }
    while(qypos >= 0){
        affinealign->transcript[affinealign->transcriptlen++]
                      = AFFINE_ALIGN_DELETE;
        qypos--;
        }
    while(dbpos >= 0){
        affinealign->transcript[affinealign->transcriptlen++]
                  = AFFINE_ALIGN_INSERT;
        dbpos--;
        }
    affinealign->transcript[affinealign->transcriptlen]
                  = AFFINE_ALIGN_END;
    reverseAFFINEtranscript(affinealign->transcript,
                            affinealign->transcriptlen);
    return;
    }

#ifdef NOTNOW
static void printTraceBackMatrix(AFFINE *affine,
                                 AFFINEALIGN *affinealign,
                                 unsigned char *dbseq, int dbseqlen,
                                 char **matrix){
    register int i, j;
    printf("Affine Traceback matrix for [%s] vs [%s]\n",
            affine->qyseq, dbseq);
    for(i = 0; i < affinealign->dbalignlen; i++){
        for(j = 0; j < da->qyalignlen; j++){
            printf("[%c%c%c%c%c]",
                AFFINE_GET_ALIGN_AFFINE(j, i,
                AFFINE_ALIGN_INSERT_OPEN)? '^':' ',
                AFFINE_GET_ALIGN_AFFINE(j, i,
                    AFFINE_ALIGN_INSERT    )? '!':' ',
                AFFINE_GET_ALIGN_AFFINE(j, i,
                    AFFINE_ALIGN_MATCH     )?'\\':' ',
                AFFINE_GET_ALIGN_AFFINE(j, i,
                    AFFINE_ALIGN_DELETE_OPEN)? '<':' ',
                AFFINE_GET_ALIGN_AFFINE(j, i,
                    AFFINE_ALIGN_DELETE    )? '-':' '
            );
            }
        printf("\n");
        }
    return;
    }
#endif /* NOTNOW */

static int generateAffineTranscript(AFFINE *affine,
                                    AFFINEALIGN *affinealign,
                            unsigned char *dbseq, int dbseqlen){
    register char **matrix = (char**)newMATRIX(affinealign->qyalignlen,
                                               affinealign->dbalignlen,
                                               sizeof(char));
    register int gec;
    register int score;
    AffineInit(affine, affinealign->qyalignstart,
               dbseq+affinealign->dbalignstart, affinealign->dbalignlen,
               (gec = 0), (gec += gapextend), gapopen+gec,
               DownScore);
    AffineMainLoop((gec = 0), (gec += (gapextend+gapextend))){
        AffineMainInnerLoop(a, b, 1, gapopen+gec,
                            dbseq+affinealign->dbalignstart,
                            AffineNullCommand, AffineNullCommand,
                DownExtend, DownScore, RightExtend, RightScore);
        AffineMainInnerLoop(b, a, 2, gapopen+gec+gapextend,
                            dbseq+affinealign->dbalignstart,
                            AffineNullCommand, AffineNullCommand,
                DownExtend, DownScore, RightExtend, RightScore);
        }
    if(oddrow){
        AffineMainInnerLoop(a, b, 1, gapopen+gec,
                            dbseq+affinealign->dbalignstart,
                            AffineNullCommand, AffineNullCommand,
                DownExtend, DownScore, RightExtend, RightScore);
        AffineLastRow(b, gapopen+gec, AffineNullCommand,
                       RightExtend, RightScore);
        score = cell[j].scoreb;
    } else {
        AffineLastRow(a, gapopen+gec, AffineNullCommand,
                       RightExtend, RightScore);
        score = cell[j].scorea;
        }
    /* printTraceBackMatrix(affine, affinealign, */
                            /* dbseq+affinealign->dbalignstart, */
                         /* affinealign->dbalignlen, matrix); */
    /* printf("Generated transcript with score = %d\n", score); */
    ExtractBestTranscript(affine, affinealign, matrix);
    free(matrix);
    return score;
    }

static void checkTranscriptScore(AFFINE *affine,
                                 AFFINEALIGN *affinealign,
              unsigned char *dbseq, int  dbseqlen, int checkscore){
    register int score = 0, prev = AFFINE_ALIGN_END, *tsp;
    register unsigned char *qyptr = affine->qyseq
                                  + affinealign->qyalignstart,
                           *dbptr = dbseq
                                  + affinealign->dbalignstart;
    register unsigned char *smi = global_submat_index;
    for(tsp = affinealign->transcript; *tsp != AFFINE_ALIGN_END; tsp++){
        switch(*tsp){
            case AFFINE_ALIGN_MATCH:
                score
                  += affine->submat->matrix[smi[*qyptr]][smi[*dbptr]];
                qyptr++;
                dbptr++;
                break;
            case AFFINE_ALIGN_INSERT:
                score
                  += (*tsp == prev)?affine->gapextend:affine->gapopen;
                dbptr++;
                break;
            case AFFINE_ALIGN_DELETE:
                score
                  += (*tsp == prev)?affine->gapextend:affine->gapopen;
                qyptr++;
                break;
            case AFFINE_ALIGN_END:
                break;
            }
        prev = *tsp;
        }
    if(score != checkscore)
        errmsg(ERROR_INFO, "Score mismatch [%d] [%d]\n[%s]\n[%s]",
                  score, checkscore, affine->qyseq, dbseq);
    else /* REMOVE AND CHANGE ABOVE BACK TO FATAL */
        errmsg(ERROR_INFO, "Score checked as %d", score);
    return;
    }

AFFINEALIGN *newAFFINEALIGNglobal(AFFINE *affine,
                               unsigned char *dbseq, int dbseqlen){
    register AFFINEALIGN *affinealign = NEW(AFFINEALIGN);
    register int checkscore;
    affinealign->qyalignstart = 0;
    affinealign->qyalignlen   = affine->qylen;
    affinealign->dbalignstart = 0;
    affinealign->dbalignlen   = dbseqlen;
    checkscore = generateAffineTranscript(affine, affinealign,
                                          dbseq, dbseqlen);
    checkTranscriptScore(affine, affinealign,
                         dbseq, dbseqlen, checkscore);
    return affinealign;
    }

#define ChallengeBestPos(target, candidate, bestpos, currpos) \
    if((target) < (candidate)){                               \
        (target) = (candidate);                               \
        (bestpos) = (currpos);                                \
        }

/* ??? */
#define ChallengeDownExtendTrace(target, candidate, qypos, dbpos) \
    if((target) < (candidate)){                                   \
        (target) = (candidate);                                   \
        trace[qypos].extend = dbpos;                              \
        }

#define ChallengeDownScoreTrace(target, candidate, qypos, dbpos) \
    if((target) < (candidate)){                                  \
        (target) = (candidate);                                  \
        OFFSET_ITEM(int, traceos, trace) = dbpos;                \
        }

typedef struct {
    int scorea;
    int scoreb;
    int extend;
    } AFFINETRACE_DBPOS;

AFFINEALIGN *newAFFINEALIGNbestfit(AFFINE *affine,
                               unsigned char *dbseq, int dbseqlen){
    errmsg(ERROR_WARNING, "Not yet reimplemented bestfit alignments");
    errmsg(ERROR_WARNING, "Displaying Global alignment instead");
    return newAFFINEALIGNglobal(affine, dbseq, dbseqlen);
    }

#ifdef NOTYET
    register AFFINEALIGN *affinealign = NEW(AFFINEALIGN);
    register int score, checkscore, traceright, dbalignend;
    int gec;
    register int traceos,
              traceosa = offsetof(AFFINETRACE_DBPOS, scorea),
              traceosb = offsetof(AFFINETRACE_DBPOS, scoreb);
    register AFFINETRACE_DBPOS *trace
                     = calloc(affine->qylen, sizeof(AFFINETRACE_DBPOS));
    AffineInit(affine, 0, dbseq, dbseqlen,
               (gec = 0), (gec += gapextend),
               gapopen+gec, Score);
    affinealign->qyalignstart = 0;
    affinealign->qyalignlen   = affine->qylen;
    score = cell[j].scorea;
    affinealign->dbalignstart = 0;
    dbalignend = 0;
    AffineMainLoop(AffineNullCommand, AffineNullCommand){
        traceright = i+1;
        traceos = traceosa;
        AffineMainInnerLoop(a, b, 1, 0, dbseq,
                            AffineNullCommand, AffineNullCommand,
                            DownExtendTrace, Score, Score, Score);
        ChallengeBestPos(score, cell[j].scorea, dbalignend, i);
        traceright = i+2;
        traceos = traceosb;
        AffineMainInnerLoop(b, a, 2, 0, dbseq,
                            AffineNullCommand, AffineNullCommand,
                            Score, Score, Score, Score);
        ChallengeBestPos(score, cell[j].scoreb, dbalignend, i+1);
        }
    if(oddrow){
        AffineMainInnerLoop(a, b, 1, 0, dbseq,
                            AffineNullCommand, AffineNullCommand,
                            Score, Score, Score, Score);
        ChallengeBestPos(score, cell[j].scorea, dbalignend, i);
        AffineLastRow(b, 0, AffineNullCommand, Score, Score);
        ChallengeBestPos(score, cell[j].scoreb, dbalignend, i+1);
    } else {
        AffineLastRow(a, 0, AffineNullCommand, Score, Score);
        ChallengeBestPos(score, cell[j].scorea, dbalignend, i+1);
        }
    printf("dbalignstart = %d\n", affinealign->dbalignstart);
    printf("dbalignend   = %d\n", dbalignend);
    affinealign->dbalignlen = dbalignend-affinealign->dbalignstart;
    checkscore = generateAffineTranscript(affine, affinealign,
                                          dbseq, dbseqlen);
    checkTranscriptScore(affine, affinealign,
                         dbseq, dbseqlen, checkscore);
    free(trace);
    return affinealign;
    }
#endif /* NOTYET */

AFFINEALIGN *newAFFINEALIGNlocal(AFFINE *affine,
                            unsigned char *dbseq, int dbseqlen){
    errmsg(ERROR_WARNING, "Not yet reimplemented local alignments");
    errmsg(ERROR_WARNING, "Displaying global alignment instead");
    return newAFFINEALIGNglobal(affine, dbseq, dbseqlen);
    }

AFFINEALIGN *newAFFINEALIGNoverlap(AFFINE *affine,
                            unsigned char *dbseq, int dbseqlen){
    errmsg(ERROR_WARNING, "Not yet reimplemented overlap alignments");
    errmsg(ERROR_WARNING, "Displaying global alignment instead");
    return newAFFINEALIGNglobal(affine, dbseq, dbseqlen);
    }

void freeAFFINEALIGN(AFFINEALIGN *affinealign){
    free(affinealign->transcript);
    free(affinealign);
    return;
    }

/* END OF ALIGNMENT GENERATION CODE */

/* START OF ALIGNMENT DISPLAY CODE */

#define AFFINE_ALIGN_STYLE_MIDSEQ    (1<<0)
#define AFFINE_ALIGN_STYLE_UNALIGNED (1<<1)

#define AFFINE_UNALIGNED_MATCH  (1<<6)
#define AFFINE_UNALIGNED_DELETE (1<<7)
#define AFFINE_UNALIGNED_INSERT (1<<8)

static int *newUnalignedTranscript(AFFINEALIGN *affinealign,
                                   int qyseqlen, int dbseqlen){
    register int i, startmax, diff, qytoend, dbtoend,
                 endmin, endmax;
    register int *ts;
    startmax = Max(affinealign->qyalignstart,
                   affinealign->dbalignstart);
    qytoend = qyseqlen-affinealign->qyalignstart
                      -affinealign->qyalignlen;
    dbtoend = dbseqlen-affinealign->dbalignstart
                      -affinealign->dbalignlen;
    endmax = Max(qytoend, dbtoend);
    ts = malloc(sizeof(int)*(affinealign->transcriptlen
                             +startmax+endmax+1));
    i = 0;
    if(affinealign->qyalignstart > affinealign->dbalignstart){
        diff = affinealign->qyalignstart-affinealign->dbalignstart;
        while(i < diff)
            ts[i++] = AFFINE_UNALIGNED_DELETE;
    } else {
        diff = affinealign->dbalignstart-affinealign->qyalignstart;
        while(i < diff)
            ts[i++] = AFFINE_UNALIGNED_INSERT;
        }
    while(i < startmax)
        ts[i++] = AFFINE_UNALIGNED_MATCH;
    memcpy(ts+i, affinealign->transcript,
                 sizeof(int)*(affinealign->transcriptlen));
    i += affinealign->transcriptlen;
    endmin = Min(qytoend, dbtoend)+i;
    while(i < endmin)
        ts[i++] = AFFINE_UNALIGNED_MATCH;
    if(qytoend > dbtoend){
        diff = qytoend-dbtoend+i;
        while(i < diff)
            ts[i++] = AFFINE_UNALIGNED_DELETE;
    } else {
        diff = dbtoend-qytoend+i;
        while(i < diff)
            ts[i++] = AFFINE_UNALIGNED_INSERT;
        }
    ts[i] = AFFINE_ALIGN_END;
    if(startmax+qytoend+dbtoend)
        errmsg(ERROR_INFO, "Unaligned sequence displayed in lowercase");
    return ts;
    }

int displayAFFINEALIGN(AFFINE *affine, AFFINEALIGN *affinealign,
                       unsigned char *dbseq, int dbseqlen,
                       FILE *output, int style, int width){
    register int i, total = 0, submatscore, *tsaddr, *tsptr;
    register char **align = (char**)newMATRIX(3, width+1, sizeof(char));
    register int qyleft, qyright, dbleft, dbright;
    register unsigned char *qyptr, *dbptr;
    register unsigned char *smi = global_submat_index;
    align[0][width] = align[1][width] = align[2][width] = '\0';
    style |= AFFINE_ALIGN_STYLE_UNALIGNED; /* TEMP HERE */
    if(style & AFFINE_ALIGN_STYLE_UNALIGNED){
        tsptr = newUnalignedTranscript(affinealign,
                             affine->qylen, dbseqlen);
        qyleft = 0;
        dbleft = 0;
        qyptr  = affine->qyseq;
        dbptr  = dbseq;
    } else {
        tsptr = affinealign->transcript;
        qyleft = affinealign->qyalignstart;
        dbleft = affinealign->dbalignstart;
        qyptr  = affine->qyseq+qyleft;
        dbptr  = dbseq+dbleft;
        }
    tsaddr = tsptr;
    total += fprintf(output, "   Query Sequence Length %d\n"
                             "Database Sequence Length %d\n\n",
                             affine->qylen,
                             dbseqlen);
    while(*tsptr != AFFINE_ALIGN_END){
        for(i = 0; i < width; i++){
            switch(*tsptr++){
                case AFFINE_ALIGN_MATCH:
                    align[0][i] = *qyptr;
                    submatscore
                    = affine->submat->matrix[smi[*qyptr]][smi[*dbptr]];
                    if(style & AFFINE_ALIGN_STYLE_MIDSEQ){
                        align[1][i] = (submatscore > 0)
                                    ? ((*qyptr == *dbptr)
                                    ? *qyptr:'+'):' ';
                    } else { /* DEFAULT */
                        align[1][i] = (submatscore > 0)
                                    ? ((*qyptr == *dbptr)
                                    ?'|':':'):' ';
                        }
                    align[2][i] = *dbptr;
                    qyptr++;
                    dbptr++;
                    break;
                case AFFINE_ALIGN_INSERT:
                    align[0][i] = '-';
                    align[1][i] = ' ';
                    align[2][i] = *dbptr++;
                    break;
                case AFFINE_ALIGN_DELETE:
                    align[0][i] = *qyptr++;
                    align[1][i] = ' ';
                    align[2][i] = '-';
                    break;
                case AFFINE_ALIGN_END:
                    align[0][i] = align[1][i] = align[2][i] = '\0';
                    i = width;
                    tsptr--;
                    break;
                case AFFINE_UNALIGNED_MATCH:
                    align[0][i] = tolower(*qyptr);
                    align[1][i] = ' ';
                    align[2][i] = tolower(*dbptr);
                    qyptr++;
                    dbptr++;
                    break;
                case AFFINE_UNALIGNED_DELETE:
                    align[0][i] = tolower(*qyptr);
                    align[1][i] = ' ';
                    align[2][i] = ' ';
                    qyptr++;
                    break;
                case AFFINE_UNALIGNED_INSERT:
                    align[0][i] = ' ';
                    align[1][i] = ' ';
                    align[2][i] = tolower(*dbptr);
                    dbptr++;
                    break;
                }
            }
        qyright = qyptr-affine->qyseq+affinealign->qyalignstart;
        dbright = dbptr-dbseq+affinealign->dbalignstart;
        total += fprintf(output,
                 "%10d  %s  %d\n            %s\n%10d  %s  %d\n\n",
                  qyleft+1, align[0], qyright,
                          align[1],
                  dbleft+1, align[2], dbright);
        qyleft = qyright;
        dbleft = dbright;
        }
    free(align);
    if(style & AFFINE_ALIGN_STYLE_UNALIGNED)
        free(tsaddr);
    return total;
    }

/* END OF ALIGNMENT DISPLAY CODE */

/* setAFFINEfunc: GIVEN name [global|bestfit|local|overlap],
                  SETS dsf AND daf APPROPRIATELY,
                  RETURNS TRUE IF NAME IS VALID, OTHERWISE FALSE.
*/
 BOOLEAN setAFFINEfunc(char *name,
                      AFFINESCOREFUNC *asf, AFFINEALIGNFUNC *aaf){
    if(!name)
        return FALSE;
    switch(*name){
        case 'G': case 'g':
            if(strcasecmp(name, "global"))
                return FALSE;
            *asf = globalAFFINEscore;
            *aaf = newAFFINEALIGNglobal;
            return TRUE;
        case 'B': case 'b':
            if(strcasecmp(name, "bestfit"))
                return FALSE;
            *asf = bestfitAFFINEscore;
            *aaf = newAFFINEALIGNbestfit;
            return TRUE;
        case 'L': case 'l':
            if(strcasecmp(name, "local"))
                return FALSE;
            *asf = localAFFINEscore;
            *aaf = newAFFINEALIGNlocal;
            return TRUE;
        case 'O': case 'o':
            if(strcasecmp(name, "overlap"))
                return FALSE;
            *asf = overlapAFFINEscore;
            *aaf = newAFFINEALIGNoverlap;
            return TRUE;
        }
    return FALSE;
    }

#ifdef TEST_THIS_MODULE

int main(){
/*
    register char *a =
      "AKVAVLGASGGIGQPLSLLLKNSPLVSRLTLYDIAHTPGVAADLSHIETRATVKGYLGPE"
      "QLPDCLKGCDVVVIPAGVPRKPGMTRDDLFNTNATIVATLTAACAQHCPDAMICIISNPV"
      "NSTIPITAEVFKKHGVYNPNKIF";
    register char *b =
      "ATLKDKLIGHLATSQEPRSYNKITVVGVGAVGMACAISILMKDLADEVALVDVMEDKLKG"
      "EMMDLQHGSLFLHTAKIVSGKDYSVSAGSKLVVITAGARQQEGESRLNLVQRNVNIFKFI"
      "IPNIVKHSPDCIILVVSNPVDVLTYVAWKLSGLPMHRII";

*/
    /* register SUBMAT *s = newSUBMAT("blosum"); */
    /* register char *a = "YEM", *b = "YEEFM"; */
    /* register char *a = "YEEFM", *b = "YEM"; */
    /* register char *a = "ARDCQ", *b = "ANQ"; */
    /* register char *a = "YEEFM", *b = "YEM"; */
    /* register char *a = "QCDRA", *b = "QNA"; */

    register SUBMAT *s = newSUBMAT("nucleic");
    /* register char *a = "ACGT", *b = "AT"; */
    /* register char *a = "TTTT", *b = "AAAATTTTAAAA"; */
    register char *a = "CGCCCGGCCTCT", *b = "TTTATAGCCTCT";
    /* register char *a = "AAAAAAAAAAAA", *b = "AAAAAAAAAAAA"; */
    /* register char *b = "CGCCCGGCCTCT", *a = "TTTATAGCCTCT"; */
    /* register char *a = "AAAAAATTT", *b = "TTTCCCCCCCC"; */
/*
    "CGCCCGGCCTCT" q
    "------GCCTCT" d ge=1: 30-12- 5 = 13
                     ge=2: 30-12-10 =  8
*/
    /* register char *a = "CGCCCGGCCTCT", */
                  /* *b = "TTTATAGCCTCT"; */
    /* register char *a = "AC", *b = "TC"; */
    /* register char *a = "C", *b = "ACG"; */
    /* register char *a = "ACGT", *b = "TTTACGTAAA"; */
    /* register char *a = "ACG", *b = "C"; */
    /* register char *b = "ACGT", *a = "CT"; */
/*
   char *a =
     "MAYGWNGCGMGVQVNGSNGAIGLSSKYSRNTELRRVEDNDIYRLAKILDE"
     "NSCWRKLMSIIPKGMDVQACSGAGCLNFPAEIKKGFKYTAQDVFQIDEAA"
     "NRLPPDQSKSQMMIDEWKTSGKLNERPTVGVLLQLLVQAELFSAADFVAL"
     "DFLNESTPARPVDGPGALISLELLEEEMEVDNEGLSLKYQSSTATLGADA"
     "QGSVGLNLDNFEKDIVRRDKSVPQPSGNTPPIAPPRRQQRSTTNSNFATL"
     "TGTGTTSTTIPNVPNLTILNPSEQIQEPVLQPRPMNIPDLSILISNSGDL"
     "RATVSDNPSNRTSSTDPPNIPRITLLIDNSGDVNSRPNHAPAKASTATTP"
     "TASSNNLPMISALNISKGSKETLRPESRSSSSSLSKDDDDDNDGEEDGEE"
     "EYPDAFLPNLSNSEQQSSNNDSSLTTVTGTSGDNSFELTNDSSSTSNDDY"
     "ACNIPDLSELQQ";
*/
/* TUBE_DROME */
/*
   char *b =
     "MEMAKANGWAVVCTSTNTNTVPIYSKYTRCTELRRVDDNDIYKLATILDVNGCWRKLMSI"
     "IPKRLDAQACSAPGALNYQEIAAKVGLKYTAQQISLIDGTAERLTPGQSISQVMIDEWKT"
     "SGKLNERPTVGVLLQLLVHAEIYSAADFVALHFLNEPKPERPTDGPAAHISLDLCSEDLS"
     "EDMDVEDGASYQPNTSALNAAVEQARGTGMNLDYFDKHMVRRDKSVPQQLENGTSSTVPV"
     "PPPRAARSSRLLKATASNVAPTTASNAPSASNTANVPNLSILNASSKKLAASSEQTLQPQ"
     "NIPNLSILNGSSEAVLMATTSTTLDAGKSDNASNGRSSASTSTATIPNVPLITLLIENSS"
     "CEISDASDATQITSKSTATKTVPDMSTASYNNLPAISALNLNIASGAGELDGNGAKARGD"
     "NADNNSSGTNSLSNDDDEQKEDDDDDDDDDVVDVDDEEADVSLPNLSNSDHNDSSLTTVT"
     "CTSGENSFEFTNDSSSASNDDYTNNIPNLSELQQ";
*/
/* TUBE_DROVI */

    register AFFINE *affine = newAFFINE((unsigned char*)a, strlen(a),
                              -12, -2, s);
    register AFFINEALIGN *affinealign = newAFFINEALIGNglobal(affine,
                                (unsigned char*)b, strlen(b));
    register AFFINEALIGN *affinealign2 = newAFFINEALIGNglobal(affine,
                                (unsigned char*)b, strlen(b));

    printf(" global %d\n",  globalAFFINEscore(affine,
                                    (unsigned char*)b, strlen(b)));
    printf("bestfit %d\n", bestfitAFFINEscore(affine,
                                    (unsigned char*)b, strlen(b)));
    printf("  local %d\n",   localAFFINEscore(affine,
                                    (unsigned char*)b, strlen(b)));
    printf("overlap %d\n", overlapAFFINEscore(affine,
                                    (unsigned char*)b, strlen(b)));

    displayAFFINEALIGN(affine, affinealign,
                       (unsigned char*)b, strlen(b), stdout, 0, 60);
    /* displayAFFINEALIGN(affine, affinealign2, */
                       /* b, strlen(b), stdout, 0, 60); */
    freeAFFINEALIGN(affinealign);
    freeAFFINEALIGN(affinealign2);
    freeAFFINE(affine);
    freeSUBMAT(s);
    return 0;
    }

#endif /* TEST_THIS_MODULE */

#endif /* INCLUDED_AFFINE_C */

/*
TODO:
    - SEQUENCE NAME SUPPORT.
    - FIX ALIGNMENTS FOR NON-GLOBAL MODELS.
*/

