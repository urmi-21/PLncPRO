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

/* Integer based dicionary lookup (binary trie).
   Guy St.C. Slater..  October 1998.
*/

#ifndef INCLUDED_INTDICT_H
#define INCLUDED_INTDICT_H

typedef union union_intdictnode {
    union union_intdictnode *p[2];
                       void *v[2];
    } INTDICTNODE;

typedef struct {
    INTDICTNODE *root;    /* ROOT OF THE TRIE            */
            int  total;   /* NUMBER OF RECORDS PRESENT   */
    INTDICTNODE *buffer;  /* INTDICTNODE BUFFER          */
            int  bufleft; /* AMMOUNT OF BUFFER REMAINING */
    } INTDICT;

INTDICT *        newINTDICT();
   void         freeINTDICT(INTDICT *id);
   void *    uniqAddINTDICT(INTDICT *id, int k, void *v);
   void * replaceAddINTDICT(INTDICT *id, int k, void *v);
   void *     lookupINTDICT(INTDICT *id, int k);

# endif /* INCLUDED_INTDICT_H */

