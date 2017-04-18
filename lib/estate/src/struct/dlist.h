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

/* Simple Double-Linked List.
   Guy St.C. Slater.. September 1999.
*/
#ifndef INCLUDED_DLIST_H
#define INCLUDED_DLIST_H

typedef struct DLISTNODE {
    struct DLISTNODE *n;
    struct DLISTNODE *p;
           void *v;
    } DLISTNODE;

typedef struct {
    DLISTNODE *a;
    DLISTNODE *z;
    DLISTNODE  sa;
    DLISTNODE  sz;
    } DLIST;

#define occupiedDLIST(d) ((d)->a->n != (d)->z)

    DLIST *  newDLIST();
     void   freeDLIST(DLIST *d);

     void *  popDLIST(DLIST *d, DLISTNODE *n);
     void *  popDLISThead(DLIST *d);
     void *  popDLISTtail(DLIST *d);

DLISTNODE * pushDLISThead(DLIST *d, void *v);
DLISTNODE * pushDLISTtail(DLIST *d, void *v);
DLISTNODE * pushDLISTbefore(DLIST *d, DLISTNODE *n, void *v);
DLISTNODE * pushDLISTafter(DLIST *d, DLISTNODE *n, void *v);

     void  raiseDLIST(DLIST *d, DLISTNODE *n);
     void  lowerDLIST(DLIST *d, DLISTNODE *n);

#endif /* INCLUDED_DLIST_H */

