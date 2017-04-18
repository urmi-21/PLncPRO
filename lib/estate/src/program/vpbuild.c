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

/* vpbuild : Index a database with a vptree
             for sublinear time searching.
   Ideally, Clustering O(N lg N), Searching O(lg N).
   Guy St.C. Slater.   August 1997.
*/

#include <stdio.h>
#include <stdlib.h>

#include "../general/common.h"
#include "../general/error.h"
#include "../general/arg.h"
#include "../metric/vptree.h"

ALLOW_ERROR_MESSAGES;

/* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

#define ARGUMENT_COUNT 4

int estateMAIN(){
    static char *dbpath = NULL;
    static char *idxpath = NULL;
    static char *metricname = NULL;
    static int  wordlen = 0;
    ARGUMENT argdata[ARGUMENT_COUNT] = {
  {ARG_STRING, 'd', "database", "fasta format DNA sequences",
   &dbpath, NULL, "VPBUILD_DATABASE", NULL, TRUE},
  {ARG_STRING, 'i', "index", "output path for index",
   &idxpath, NULL, "VPBUILD_INDEX", NULL, TRUE},
  {ARG_STRING, 'm', "metric", "metric name",
   &metricname, NULL, "VPBUILD_METRIC", "[edit,word]", TRUE},
  {ARG_INT, 'w', "wordlen", "word length",
   &wordlen, NULL, "VPBUILD_WORDLEN", NULL, FALSE}
  };
    METRIC *metric = (METRIC*)0;
    processARGUMENT(argdata, ARGUMENT_COUNT, global_emn,
                   "build vptree for sublinear time searching");
    if(!strcasecmp(metricname, "edit")){
        errmsg(ERROR_WARNING, "Using edit distance metric");
        metric = newMETRIC(METRIC_EDITDIST, NULL);
    } else if(!strcasecmp(metricname, "word")){
        errmsg(ERROR_WARNING, "Using word based distance metric");
        if(wordlen == 0)
            errmsg(ERROR_WARNING, "Using default wordlen=%d",
                    wordlen=9);
        metric = newMETRIC(METRIC_WORDDIST, newMETRICPARAMwd(wordlen));
        }
    makeVPTREE(dbpath, idxpath, metric);
    return 0;
    }


