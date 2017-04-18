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

/* Fast Adaptive Dictionary Data Structure.
   Guy St.C. Slater. December 1996.  Version 2.0
*/

#ifndef INCLUDED_DICT_H
#define INCLUDED_DICT_H

#include "common.h"

typedef struct DICT { struct DICT *a, *d; char c; } DICT;

DICT *newDICT();
void  freeDICT(DICT *d);
void freeDICTptr(DICT *n, FREEFUNC ff);
void *uniqAddDICT(DICT *n, char *w, void *v);
void *uniqAddDICTlen(DICT *n, char *w, void *v, int len);
void *replaceAddDICT(DICT *n, char *w, void *v);
long  countDICT(DICT *n, char *w);
void *lookupDICT(DICT *n, char *w);
long  walkDICT(DICT *n, WALKFUNC f);
long  walkDICTinfo(DICT *n, WALKFUNCINFO f, void *info);
long  sortWalkDICT(DICT *n, WALKFUNC f);
DICT *joinDICT(DICT *a, DICT *b, JOINFUNC f);
void mergeAddDICT(DICT *n, char *w, void *v, JOINFUNC f);

#endif /* INCLUDED_DICT_H */

