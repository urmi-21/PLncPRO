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

/* Boyer-Moore Exact String Matching. 
   Guy St.C.Slater.  April 1997.
*/
#ifndef INCLUDED_BOYERMOORE_H
#define INCLUDED_BOYERMOORE_H

 int *BMindex(unsigned char *p, int len); 
unsigned char *BMsearch(int *index, unsigned char *p,
                                    unsigned char *a, int alen);

#endif /* INCLUDED_BOYERMOORE_H */

