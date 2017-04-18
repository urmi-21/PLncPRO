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

#ifndef INCLUDED_EDITDIST_H
#define INCLUDED_EDITDIST_H

#include <stdio.h>

typedef struct {
    char *query;
     int  length;
     int *vectorA;
     int *vectorB;
    } EDITDIST;

EDITDIST *newEDITDIST(char *qyseq, int qylen);
   void  freeEDITDIST(EDITDIST *ed);
    int  calcEDITDIST(EDITDIST *ed, char *dbseq, int dblen);
    int  calcEDITDISTonline(EDITDIST *ed, FILE *fp, long pos);

#endif /* INCLUDED_EDITDIST_H */

