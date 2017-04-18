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

/* string: Basic String Data Structure 
   Guy St.Slater.  January 1997.
*/

#ifndef INCLUDED_STRING_H
#define INCLUDED_STRING_H

typedef struct {
    char *str;
     int  len;
     int  alloc;
    } STRING;

STRING *newSTRING();
  void  freeSTRING(STRING *s);
STRING *makeSTRING(char *string, int len);
STRING *copySTRING(STRING *a);
STRING *joinSTRING(STRING *a, STRING *b);
  void  linetoSTRING(STRING *s, char *line, int len);
  void  addSTRING(STRING *s, char ch);
STRING *dlinetoSTRING(STRING *s, char *line, int len, char delim);
  void  reverseSTRING(STRING *s);

#endif /* INCLUDED_STRING_H */

