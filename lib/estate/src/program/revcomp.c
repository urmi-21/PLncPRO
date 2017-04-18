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

/* revcomp : generate reverse complement of nucleotide sequences.
*/

#include <ctype.h> /* for toupper() */

#include "../general/arg.h"
#include "../general/error.h"
#include "../general/common.h"
#include "../parse/readfasta.h"
#include "../sequence/sequtil.h"

ALLOW_ERROR_MESSAGES;

#define TRANSLATED_SEQ_CHUNKSIZE 100

static void reverseComplementDatabase(char *path){
    register READFASTA *rf = newREADFASTA(path);
    while(nextREADFASTA(rf)){
        SEQUTILrevcomp((unsigned char*)rf->seq, rf->seqlen);
        printf(">%s [revcomp]\n", rf->def);
        SEQUTILwriteFASTAblock(stdout, rf->seq, rf->seqlen);
        }
    freeREADFASTA(rf);
    return;
    }

#define ARGUMENT_COUNT 1

int estateMAIN(){
    static char *dbpath = NULL;
    ARGUMENT argdata[ARGUMENT_COUNT] = {
    { ARG_STRING,  'd', "database", "fasta format DNA sequences",
     &dbpath, NULL, "REVCOMP_DATABASE", NULL, TRUE},
    };
    processARGUMENT(argdata, ARGUMENT_COUNT, global_emn,
      "generate reverse complement of nucleotide sequences");
    reverseComplementDatabase(dbpath);
    return 0;
    }

