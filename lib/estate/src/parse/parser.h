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

/* parser : Basic line based file parsing.
   Guy St.C. Slater.  May 1998.
*/

#ifndef INCLUDED_PARSER_H
#define INCLUDED_PARSER_H

#include <stdio.h>

#define PARSER_CHUNK_SIZE 100

typedef struct {
         FILE *fp;
unsigned char *line;
          int  linealloc;
unsigned char **word;
          int  wordalloc;
    } PARSER;

PARSER *newPARSER(FILE *fp);
  void freePARSER(PARSER *p);
   int linePARSER(PARSER *p);
   int wordPARSER(PARSER *p);
   int wordPARSERoffset(PARSER *p, int length, int start);
char **linePARSERfile(FILE *fp, int *total);

#endif /* INCLUDED_PARSER_H */


