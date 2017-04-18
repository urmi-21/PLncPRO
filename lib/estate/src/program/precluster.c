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

/* precluster : preclustering step for EST clustering.
   Guy St.C. Slater.  January 1999.  Version 2.0
*/

#include "../general/common.h"
#include "../general/error.h"
#include "../general/arg.h"

#include "../parse/rafasta.h"

#include "../cluster/wordgraph.h"
#include "../cluster/avastore.h"
#include "../cluster/clusterset.h"

ALLOW_ERROR_MESSAGES;

#define ARGUMENT_COUNT 7

int estateMAIN(){
    static BOOLEAN usevfsm = TRUE;
    static FILE *dbfp, *outfp, *outrcfp;
    static char *dbpath = NULL, *outpath = NULL, *outrcpath;
    static int wordlen = 9;
    static int wordmin = 24;
    static BOOLEAN revcomp = TRUE;
    ARGUMENT argdata[ARGUMENT_COUNT] = {
    { ARG_BOOLEAN, 'V', "vfsm", "use vfsm optimisation",
      &usevfsm, NULL, "PRECLUSTER_USE_VFSM", NULL, FALSE},
    { ARG_FILE,  'd', "database", "fasta database path",
      &dbfp, &dbpath, "PRECLUSTER_DATABASE_PATH", "r", TRUE},
    { ARG_INT,   'w', "wordlen", "word length",
      &wordlen, NULL, "PRECLUSTER_WORD_LENGTH", NULL, FALSE},
    { ARG_INT,   'm', "wordmin", "minimum word matches required",
      &wordmin, NULL, "PRECLUSTER_WORD_MINIMUM", NULL, FALSE},
    { ARG_FILE,  'A', "allvsall", "allvsall results path",
      &outfp, &outpath, "PRECLUSTER_AVARESULT_PATH", "wb", TRUE},
    { ARG_FILE,  'R', "allvsallrc", "allvsall revcomp results path",
      &outrcfp, &outrcpath, "PRECLUSTER_AVARESULT_REVCOMP_PATH",
      "wb", FALSE},
    { ARG_BOOLEAN,  'r', "revcomp", "use reverse complement",
      &revcomp, NULL, "PRECLUSTER_REVCOMP", NULL, FALSE}
    };
    register int total;
    dbfp = stdin;
    processARGUMENT(argdata, ARGUMENT_COUNT, global_emn,
                    "wordscore prestep to EST clustering");
    if(!usevfsm)
        errmsg(ERROR_WARNING, "Not using VFSM - this be quite slow");
    if(revcomp){
        if(!outrcpath)
            errmsg(ERROR_FATAL, "Output path required for revcomp");
        total = usevfsm
              ? vfsmWORDGRAPHrevcomp(dbpath, outpath,
                                    outrcpath, wordlen, wordmin)
              : fsmWORDGRAPHrevcomp(dbpath, outpath,
                                    outrcpath, wordlen, wordmin);
    } else {
        total = usevfsm
              ? vfsmWORDGRAPH(dbpath, outpath, wordlen, wordmin)
              : fsmWORDGRAPH(dbpath, outpath, wordlen, wordmin);
        }
    errmsg(ERROR_INFO, "Wrote a total of %d scores", total);
    return 0;
    }

