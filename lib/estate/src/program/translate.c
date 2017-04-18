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

/* translate : translate sequences in up to six frames.
*/

#include <ctype.h> /* for toupper() */

#include "../general/arg.h"
#include "../general/error.h"
#include "../general/common.h"
#include "../parse/readfasta.h"
#include "../sequence/ntaa.h"
#include "../sequence/sequtil.h"

ALLOW_ERROR_MESSAGES;

#define TRANSLATED_SEQ_CHUNKSIZE 100

static void sixFrameTranslate(char *path){
    register READFASTA *rf = newREADFASTA(path);
    register int aaseqalloc = TRANSLATED_SEQ_CHUNKSIZE;
    register char *aaseq = malloc(sizeof(char)*aaseqalloc);
    register int t, frame;
    NTAAInitialiseData();
    while(nextREADFASTA(rf)){
        t = (rf->seqlen/3)+2;
        if(t > aaseqalloc){
            aaseqalloc += TRANSLATED_SEQ_CHUNKSIZE;
            if(aaseqalloc < t)
                aaseqalloc = t+TRANSLATED_SEQ_CHUNKSIZE;
            aaseq = realloc(aaseq, sizeof(char)*aaseqalloc);
            }
        for(frame = 1; frame < 4; frame++){
            printf(">%s (Frame %d)\n", rf->def, frame);
            translateNTAAseq((unsigned char*)rf->seq, rf->seqlen,
                            frame, (unsigned char*)aaseq);
            SEQUTILwriteFASTAblock(stdout, aaseq, strlen(aaseq));
            }
        for(frame = -1; frame > -4; frame--){
            printf(">%s (Frame %d)\n", rf->def, frame);
            translateNTAAseq((unsigned char*)rf->seq, rf->seqlen,
                            frame, (unsigned char*)aaseq);
            SEQUTILwriteFASTAblock(stdout, aaseq, strlen(aaseq));
            }
        }
    freeREADFASTA(rf);
    free(aaseq);
    return;
    }

#define ARGUMENT_COUNT 1

int estateMAIN(){
    static char *dbpath = NULL;
    ARGUMENT argdata[ARGUMENT_COUNT] = {
    { ARG_STRING,  'd', "database", "fasta format DNA sequences",
     &dbpath, NULL, "TRANSLATE_DATABASE", NULL, TRUE},
    };
    processARGUMENT(argdata, ARGUMENT_COUNT, global_emn,
             "six frame DNA translation");
    NTAAInitialiseData();
    sixFrameTranslate(dbpath);
    return 0;
    }

