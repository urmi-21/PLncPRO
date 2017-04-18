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

/* vpsearch: Search a database in sub-linear time,
             which has been indexed by vpbuild.
   Ideally, Clustering O(N lg N), Searching O(lg N).
   Guy St.C. Slater.   August 1997.
*/

#include <stdio.h>
#include <stdlib.h>

#include "../general/common.h"
#include "../general/error.h"
#include "../general/arg.h"
#include "../metric/metric.h"
#include "../metric/vptree.h"
#include "../sequence/sequtil.h"
#include "../parse/ioutil.h"

ALLOW_ERROR_MESSAGES;

/* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

#define ARGUMENT_COUNT 2

int estateMAIN(){
    static char *idxpath = NULL;
    static char *qypath = NULL;
    ARGUMENT argdata[ARGUMENT_COUNT] = {
   {ARG_STRING, 'i', "index", "vpbuild index path",
    &idxpath, NULL, "VPSEARCH_INDEX", NULL, TRUE},
   {ARG_STRING, 'q', "query", "query sequence path",
    &qypath, NULL, "VPSEARCH_QUERY", NULL, TRUE},
   };
    register char *qyseq;
    long qylen;
    processARGUMENT(argdata, ARGUMENT_COUNT, global_emn,
        "sublinear time database searching using a vptree");
    if(!(qyseq = (char*)readFILE(qypath, &qylen)))
        errmsg(ERROR_FATAL, "Could not open query file \"%s\"",
                             qypath);
    qylen = SEQUTILclean(qyseq);
    /* displayVPTREE(index); */
    searchVPTREE(idxpath, qyseq);
    free(qyseq);
    return 0;
    }

