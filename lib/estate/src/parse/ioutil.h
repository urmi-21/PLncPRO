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

/* ioutil: Input/Output Utilities
   Guy St.C. Slater.  December 1996.
*/

#ifndef INCLUDED_IOUTIL_H
#define INCLUDED_IOUTIL_H

#include "../general/common.h"

void  catFILE(FILE *in, FILE *out);
void  bitprint(FILE *fp, int i, int n);
char *readFILE(char *path, long *length);
 int  seekFILE(FILE *fp, char *s);
 int  seekFILEpipe(FILE *fp, char *s);
char *tmpFILEpath(char *base);
void reversestring(char *str, int len);

#endif /* INCLUDED_IOUTIL_H */

