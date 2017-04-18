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

#ifndef INCLUDED_AFFINE_H
#define INCLUDED_AFFINE_H

#include <stdio.h>
#include "submat.h"
#include "../general/common.h"

typedef struct {
  signed  int  scorea;
  signed  int  scoreb;
  signed  int  extend;
  signed char *matrix;
} AFFINE_CELL;

typedef struct {
unsigned char *qyseq;
          int  qylen;
  AFFINE_CELL *cell;
       SUBMAT *submat;
          int  gapopen;
          int  gapextend;
} AFFINE;

AFFINE *newAFFINE(unsigned char *qyseq,  int qylen,
                  int gapopen, int gapextend, SUBMAT *submat);
void freeAFFINE(AFFINE *affine);

int  globalAFFINEscore(AFFINE *affine,
                       unsigned char *dbseq, int dbseqlen);
int bestfitAFFINEscore(AFFINE *affine,
                       unsigned char *dbseq, int dbseqlen);
int   localAFFINEscore(AFFINE *affine,
                       unsigned char *dbseq, int dbseqlen);
int overlapAFFINEscore(AFFINE *affine,
                       unsigned char *dbseq, int dbseqlen);

typedef struct {
    int  qyalignstart;
    int  qyalignlen;
    int  dbalignstart;
    int  dbalignlen;
    int *transcript;
    int  transcriptlen;
} AFFINEALIGN;

AFFINEALIGN *newAFFINEALIGNglobal(AFFINE *affine,
                                  unsigned char *dbseq, int dbseqlen);
AFFINEALIGN *newAFFINEALIGNbestfit(AFFINE *affine,
                                  unsigned char *dbseq, int dbseqlen);
AFFINEALIGN *newAFFINEALIGNlocal(AFFINE *affine,
                                  unsigned char *dbseq, int dbseqlen);
AFFINEALIGN *newAFFINEALIGNoverlap(AFFINE *affine,
                                  unsigned char *dbseq, int dbseqlen);

void freeAFFINEALIGN(AFFINEALIGN *affinealign);

int displayAFFINEALIGN(AFFINE *affine, AFFINEALIGN *affinealign,
                       unsigned char *dbseq, int dbseqlen,
                       FILE *output, int style, int width);

typedef int (*AFFINESCOREFUNC) (AFFINE *affine,
                                unsigned char *dbseq, int dbseqlen);
typedef AFFINEALIGN *(*AFFINEALIGNFUNC) (AFFINE *affine,
                                unsigned char *dbseq, int dbseqlen);
BOOLEAN setAFFINEfunc(char *name,
                      AFFINESCOREFUNC *asf, AFFINEALIGNFUNC *aaf);

#endif /* INCLUDED_AFFINE_H */

