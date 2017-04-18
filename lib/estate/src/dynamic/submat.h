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

#ifndef INCLUDED_SUBMAT_H
#define INCLUDED_SUBMAT_H

#define SUBMAT_ALPHASIZE 24

#include "../general/common.h"

typedef signed char SUBMATMATRIX[SUBMAT_ALPHASIZE][SUBMAT_ALPHASIZE];

typedef struct {
    SUBMATMATRIX matrix;
           char *path;
    } SUBMAT;

extern char global_submat_aa[SUBMAT_ALPHASIZE];
extern unsigned char global_submat_index[ALPHABETSIZE];

SUBMAT * newSUBMAT(char *path);
  void  freeSUBMAT(SUBMAT *s);

   int  selfSUBMAT(SUBMAT *s, unsigned char *seq);

/* IF PATH IS {blosum,pam,nucleic,edit},
   WILL USE A HARD-CODED MATRIX.
*/

#endif /* INCLUDED_SUBMAT_H */

