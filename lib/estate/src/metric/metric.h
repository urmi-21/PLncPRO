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

/* Metric Sequence Comparison Routines
   Guy St.C. Slater..  September 1997. Version 1.1
*/

#ifndef INCLUDED_METRIC_H
#define INCLUDED_METRIC_H

#include <stdio.h>

#define METRIC_TOTAL    3

#define METRIC_UNKNOWN  0
#define METRIC_EDITDIST 1
#define METRIC_WORDDIST 2

typedef union {
    struct { int wordlength; } wd;
    } METRICPARAM;

typedef struct {
   void *(*prep )(char *qy, METRICPARAM *param);    /* PREPARE QUERY */
    int  (*dist )(void *qydat, long pos, FILE *fp); /* GET DISTANCE  */
    int  (*cache)(void *qydat, char *seq, int len); /* GET DISTANCE  */
   void  (*free )(void *qydat);               /* FREE THE QUERY DATA */
METRICPARAM *param;
    char     type;
    } METRIC; /* THE SEQUENCE COMPARISON METRIC FUNCTIONS */

METRICPARAM  *newMETRICPARAMwd();
     METRIC  *newMETRIC(int type, METRICPARAM *param);
       void  freeMETRIC(METRIC *metric);

       void  writeMETRIC(FILE *fp, METRIC *metric);
     METRIC * readMETRIC(FILE *fp);

# endif /* INCLUDED_METRIC_H */

