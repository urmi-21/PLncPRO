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

/* Word graph generating code.
   (Generate word match scores for all vs all comparison of sequences.)
   Guy St.C. Slater.  November 1998.
*/

#ifndef INCLUDED_WORDGRAPH_H
#define INCLUDED_WORDGRAPH_H

#include "../general/common.h"

int  fsmWORDGRAPH(char *fastapath, char *avapath,
                  int wordlen, int wordin);
int vfsmWORDGRAPH(char *fastapath, char *avapath,
                  int wordlen, int wordmin);

int   fsmWORDGRAPHrevcomp(char *fastapath, char *avapath,
                          char *avarcpath, int wordlen, int wordmin);
int  vfsmWORDGRAPHrevcomp(char *fastapath, char *avapath,
                          char *avarcpath, int wordlen, int wordmin);

#endif /* INCLUDED_WORDGRAPH_H */

