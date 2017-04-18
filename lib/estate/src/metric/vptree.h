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

/* vptree.h : Vantage Point Tree Algorithm.
   Guy St.C. Slater.  August 1997.

   Reference:
   "Data Structures and Algorithms for Nearest
    Neighbour Search in General Metric Spaces"
   Proc 4th ACM-SIAM Symposium Discrete Algorithms.
   pp 311-321, 1993. Peter N. Yianilos.
*/

#ifndef INCLUDED_VPTREE_H
#define INCLUDED_VPTREE_H

#include "metric.h"

void   makeVPTREE(char *dbpath, char *idxpath, METRIC *metric);
void searchVPTREE(char *idxpath, char *query);

#endif /* INCLUDED_VPTREE_H */

