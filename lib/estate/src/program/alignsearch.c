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

/* alignsearch: Dynamic programming based database search.
   Guy St.C. Slater..  May 1999.
*/

#include "../general/arg.h"
#include "../general/error.h"
#include "../general/common.h"
#include "../general/stringlist.h"

ALLOW_ERROR_MESSAGES;

#include "../dynamic/submat.h"
#include "../dynamic/affine.h"

#include "../parse/readfasta.h"
#include "../parse/ioutil.h"
#include "../parse/parser.h"
#include "../parse/rafasta.h"

#include "../sequence/sequtil.h"
#include "../struct/pqueue.h"

typedef struct {
    long pos;
     int score;
    } SEARCHHIT;

static BOOLEAN searchhitcomp(void *low, void *high){
    register SEARCHHIT *ha = low, *hb = high;
    return (ha->score < hb->score);
    }

typedef struct {
         char *def;
unsigned char *seq;
          int  len;
    SEARCHHIT *hit;
} SEARCHRESULT;

static void showresults(PQUEUE *pq, FILE *dbfp, AFFINE *affine,
                 int hitcount, int alignmentcount,
                 int compcount, char *searchtype,
                 AFFINEALIGNFUNC aaf){
    register SEARCHRESULT *sr = malloc(sizeof(SEARCHRESULT)
                                     *(pq->total+1));
    register int i, total = pq->total, limit;
    register AFFINEALIGN *affinealign;
    /* COLLECT DATA (IN REVERSE ORDER) */
    for(i = total-1; i >= 0; i--){
        sr[i].hit = popPQUEUE(pq);
        sr[i].def = getRAFASTAdef(dbfp, sr[i].hit->pos);
        sr[i].seq = (unsigned char*)getRAFASTAseq(dbfp,
                                    sr[i].hit->pos, &sr[i].len);
        }
    limit = Min(total, hitcount);
    if((limit > 0) && (compcount > 1)){ /* DON'T SHOW WHEN PAIRWISE */
        errmsg(ERROR_INFO, "Showing hitlist\n");
        for(i = 0; i < limit; i++){
            printf("%3d %.40s length = %d  score = %d\n",
                i+1, sr[i].def, sr[i].len, sr[i].hit->score);
            }
        }
    limit = Min(total, alignmentcount);
    if(limit > 0){
        errmsg(ERROR_INFO, "Showing alignments\n");
        for(i = 0; i < limit; i++){
            printf("------------\n"
                   "Hit %d %.40s length = %d  score = %d\n",
                i+1, sr[i].def, sr[i].len, sr[i].hit->score);
            affinealign = aaf(affine, sr[i].seq, sr[i].len);
            displayAFFINEALIGN(affine, affinealign,
                               sr[i].seq, sr[i].len,
                               stdout, 0, 60);
            freeAFFINEALIGN(affinealign);
            }
        }
    for(i = 0; i < limit; i++){
        free(sr[i].def);
        free(sr[i].seq);
        free(sr[i].hit);
        }
    free(sr);
    return;
    }

static void searchdb(char *dbpath, char *qypath,
                 int hitcount, int alignmentcount, char *searchtype,
                 char *submatpath, int gapopen, int gapextend){
    long spare;
    register READFASTA *rf = newREADFASTA(dbpath);
    register char *qyseq = readFILE(qypath, &spare);
    register long qylen = SEQUTILclean(qyseq);
    register SUBMAT *submat = newSUBMAT(submatpath);
    register AFFINE *affine = newAFFINE((unsigned char*)qyseq, qylen,
                                        gapopen, gapextend, submat);
    register PQUEUE *pq = newPQUEUE(Max(hitcount,alignmentcount),
                                    searchhitcomp);
    register SEARCHHIT *sh = NEW(SEARCHHIT);
    register int tcompctr = 0, compctr = 0;
    AFFINESCOREFUNC scorefunc = NULL;
    AFFINEALIGNFUNC alignfunc = NULL;
    if(!setAFFINEfunc(searchtype, &scorefunc, &alignfunc))
        errmsg(ERROR_FATAL, "Unknown search type [%s]", searchtype);
    errmsg(ERROR_INFO, "Searching database");
    while(nextREADFASTA(rf)){
        sh->score = scorefunc(affine, (unsigned char*)rf->seq,
                                                      rf->seqlen);
        sh->pos = READFASTA2RAFASTApos(rf);
        /* printf("got [%s][%s][%d][%ld]\n", rf->def, rf->seq, */
                                         /* sh->score, sh->pos); */
        /* printRAFASTAdef(tfp, sh->pos, stdout); */
        sh = pushLimitedPQUEUE(pq, (void*)sh);
        if(!sh)/* IF sh IS NULL, HAS BEEN USED */
            sh = NEW(SEARCHHIT);
        if(++tcompctr == 1000){
            tcompctr = 0;
            compctr += 1000;
            fputc('.', stderr);
            }
        }
    showresults(pq, rf->fp, affine, hitcount, alignmentcount,
                rf->count, searchtype, alignfunc);
    freeREADFASTA(rf);
    freeSUBMAT(submat);
    freeAFFINE(affine);
    free(qyseq);
    freePQUEUE(pq);
    free(sh);
    return;
    }

#define ARGUMENT_COUNT 8

int estateMAIN(){
    static FILE *db = NULL;
    static FILE *qy = NULL;
    static int hitcount   = 10;
    static int alignmentcount = 10;
    static int gapopen = -12;
    static int gapextend = -2;
    static char *dbname, *qyname;
    static char *searchtype = "local";
    static char *submat = "nucleic";

    ARGUMENT argdata[ARGUMENT_COUNT] = {
    {ARG_FILE, 'd', "database", "fasta format database",
    &db, &dbname, "DPDBSEARCH_DATABASE_PATH", "r", TRUE},
    {ARG_FILE, 'q', "query", "fasta format query",
    &qy, &qyname, "DPDBSEARCH_QUERY_PATH", "r", TRUE},
    {ARG_INT,  'H', "hits", "number of hits to display",
    &hitcount, NULL, "DPDBSEARCH_HITS",  NULL, FALSE},
    {ARG_INT,  'A', "alignments", "number of alignments to display",
    &alignmentcount, NULL, "DPDBSEARCH_ALIGNMENTS",  NULL, FALSE},
    {ARG_INT,  'O', "gapopen", "gap open penalty",
    &gapopen, NULL, "DPDBSEARCH_GAPOPEN",  NULL, FALSE},
    {ARG_INT,  'E', "gapextend", "gap extend penalty",
    &gapextend, NULL, "DPDBSEARCH_GAPEXTEND",  NULL, FALSE},
    {ARG_STRING,  't', "type", "search type to perform",
    &searchtype, NULL, "DPDBSEARCH_TYPE",
                           "[global,bestfit,local,overlap]", FALSE},
    {ARG_STRING,  'm', "submat", "substitution matrix to use",
    &submat, NULL, "DPDBSEARCH_SUBMAT",  NULL, FALSE}
    };
    processARGUMENT(argdata, ARGUMENT_COUNT, global_emn,
              "dynamic programming based database searching");
    searchdb(dbname, qyname, hitcount, alignmentcount, searchtype,
             submat, gapopen, gapextend);
    return 0;
    }

/*
TODO:
   DONE:search type selection
   DONE:parameter selection
   DONE:submat selection
   o put back in the m_fork() and pthreading code.
   o semi-automatic submat selection.
   o Add timing statistics. ( Mcell updates/second )
   o Add stats, local stat balancing.
*/

