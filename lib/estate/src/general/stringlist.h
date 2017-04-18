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

/* stringlist : Routines for a simple string list data structure 
   Guy St.C. Slater.  January 1997.
*/

#ifndef INCLUDED_STRINGLIST_H
#define INCLUDED_STRINGLIST_H

#include <stdio.h>

typedef struct { 
    char **v;  /* VECTOR    */
     int   c;  /* COUNT     */
     int   a;  /* ALLOCATED */
    } STRINGLIST;

STRINGLIST *  newSTRINGLIST();
      void   freeSTRINGLIST(STRINGLIST *s);
      void    addSTRINGLIST(STRINGLIST *s, char *line, int len);
STRINGLIST * copySTRINGLIST(STRINGLIST *s);
STRINGLIST * joinSTRINGLIST(STRINGLIST *a, STRINGLIST *b); 
       int  printSTRINGLIST(FILE *fp, STRINGLIST *s); 
STRINGLIST *   toSTRINGLIST(char *line, int len, char delim);

#endif /* INCLUDED_STRINGLIST_H */


