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

/* estcluster : Rapid and accurate EST clustering.
   Guy St.C. Slater.  January 1999.  Version 2.1
*/

#include <signal.h> /* FOR signal() */
#include <sys/types.h> /* FOR chmod() */
#include <sys/stat.h>  /* FOR chmod() */

#include "../general/common.h"
#include "../general/error.h"
#include "../general/arg.h"

#include "../parse/rafasta.h"

#include "../cluster/avastore.h"
#include "../cluster/clusterset.h"

#include "../dynamic/affine.h"
#include "../sequence/sequtil.h"
#include "../parse/readfasta.h"
#include "../struct/list.h"

ALLOW_ERROR_MESSAGES;

typedef struct {
    FILE *fp;
    long *index;
     int  compcount;
  SUBMAT *submat;
#define SHOW_PERFORMANCE
#ifdef SHOW_PERFORMANCE
  CLUSTERSET *cs;
         int  wordscore;
         int  lastchange;
#endif /* SHOW_PERFORMANCE */
         int  joinscore;
     BOOLEAN  thinkaloud;
     BOOLEAN  showalignments;
         int  rctype;
         int  gapopen;
         int  gapextend;
    } CLUSTERDATA;

int maxScoreLimit(SUBMAT *s, unsigned char *seq){
    register int score = 0, pos;
    while(*seq){
        pos = global_submat_index[*seq];
        if(s->matrix[pos][pos] > 0)
            score += s->matrix[pos][pos];
        seq++;
        }
    return score;
    }

static int compare_EST_sequences(void *info, int a, int b){
    register CLUSTERDATA *cd = (CLUSTERDATA*)info;
    register unsigned char *seqa, *seqb;
    int seqalen, seqblen;
    register AFFINE *af;
    register AFFINEALIGN *afa;
    register int score;

    seqa = (unsigned char*)getRAFASTAseq(cd->fp, cd->index[a],
                                         &seqalen);
    seqb = (unsigned char*)getRAFASTAseq(cd->fp, cd->index[b],
                                         &seqblen);

    af = newAFFINE(seqa, seqalen,
                   cd->gapopen, cd->gapextend,
                   cd->submat);
    if(cd->rctype == 2)
        SEQUTILrevcomp(seqb, seqblen);
    score = overlapAFFINEscore(af, seqb, seqblen);

    /* printf("Action: %s %d %d score=%d limit=%d\n", */
          /* (score >= cd->joinscore)?"join":"fail", */
           /* a, b, score, cd->joinscore); */
    if(!(++cd->compcount & 1023))
        fputc('.', stderr);
    if(cd->showalignments){
        if(score >= cd->joinscore){
            printf("Adding [%s] alignment of score [%d]\n",
                 (cd->rctype == 2)?"revcomp":"forward", score);
            afa = newAFFINEALIGNoverlap(af, seqb, seqblen);
            displayAFFINEALIGN(af, afa, seqb, seqblen,
                               stdout, 0, 60);
            freeAFFINEALIGN(afa);
            fflush(stdout);
            }
        }
    freeAFFINE(af);
    free(seqa);
    free(seqb);
#ifdef SHOW_PERFORMANCE
    if(score >= cd->joinscore)
        cd->lastchange = cd->compcount;
    if(cd->thinkaloud){
        /* seqa, seqb, clusters, singletons, dpscore, wordscore,
           compcount, lastchange */
        fprintf(stdout, "sa %d sb %d cl %d si %d ds %d "
                        "ws %d cc %d lc %d cs %d rc %d\n",
                a, b,
                cd->cs->clustertotal, cd->cs->singletontotal, score,
                cd->wordscore, cd->compcount, cd->lastchange,
                cd->compcount-cd->lastchange, cd->rctype);
        fflush(stdout);
        }
#endif /* SHOW_PERFORMANCE */
    return (score < cd->joinscore);
    }

/* ---------------------------------- */
/* START OF CLUSTERSET REPORTING CODE */

typedef struct {
    int currid;       /* CLUSTERID BEING REPORTED */
    int membercount;  /* MEMBERS SEEN IN CURRENT CLUSTER */
    int clustercount; /* CLUSTERS SEEN */
    CLUSTERDATA *cd;  /* CLUSTER DATA FOR FASTA HEADERS ETC */
    } CLUSTERSETREPORTSTATUS;

/* END OF CLUSTERSET REPORTING CODE */
/* -------------------------------- */

/* ---------------------------------- */
/* START OF CLUSTERSET REPORTING CODE */

typedef struct {
    CLUSTERSETREPORTSTATUS *crs;
                      char *outputdir; 
                      FILE *clusterfp;   /* CURRENTLY OPEN FILE */
                      FILE *singletonfp; /*      SINGLETON FILE */
    } CLUSTERSETREPORTSTATUSOUTPUT;

static void initreportclustersoutput(void *info, int membertotal,
                            int clustertotal, int singletontotal){
    register CLUSTERSETREPORTSTATUSOUTPUT *crso = info;
    register char *singletonpath;
    errmsg(ERROR_INFO, "Writing clusters to [%s]", crso->outputdir);
    crso->crs->clustercount = 0;
    singletonpath = malloc(sizeof(char)*strlen(crso->outputdir)+20);
    sprintf(singletonpath, "%s/singletons.fasta", crso->outputdir);
    crso->singletonfp = fopen(singletonpath, "w");
    if(!crso->singletonfp)
        errmsg(ERROR_FATAL, "Could not write singletons to [%s]",
                             singletonpath);
    free(singletonpath);
    return;
    }

static void singletonreportclustersoutput(void *info, int singletonid){
    register CLUSTERSETREPORTSTATUSOUTPUT *crso = info;
    register long pos = crso->crs->cd->index[singletonid-1];
    fprintf(crso->singletonfp, ">");
    printRAFASTAdef(crso->crs->cd->fp, pos, crso->singletonfp);
    printRAFASTAseq(crso->crs->cd->fp, pos, crso->singletonfp);
    return;
    }

static void startreportclustersoutput(void *info, int clusterid){
    register CLUSTERSETREPORTSTATUSOUTPUT *crso = info;
    register char *clusterpath;
    crso->crs->clustercount++;
    crso->crs->currid = clusterid;
    crso->crs->membercount  = 0;
    clusterpath = malloc(sizeof(char)*strlen(crso->outputdir)+100);
    sprintf(clusterpath, "%s/cluster%08d.fasta",
            crso->outputdir, crso->crs->clustercount);
    crso->clusterfp = fopen(clusterpath, "w");
    if(!crso->clusterfp)
        errmsg(ERROR_FATAL, "Could not write cluster [%d] to [%s]",
                             clusterid, clusterpath);
    free(clusterpath);
    fprintf(stderr, ".");
    return;
    }

static void finishreportclustersoutput(void *info){
    register CLUSTERSETREPORTSTATUSOUTPUT *crso = info;
    fclose(crso->clusterfp);
    return;
    }

static void memberreportclustersoutput(void *info, int memberid){
    register CLUSTERSETREPORTSTATUSOUTPUT *crso = info;
    register long pos = crso->crs->cd->index[memberid-1];
    fprintf(crso->clusterfp, ">");
    printRAFASTAdef(crso->crs->cd->fp, pos, crso->clusterfp);
    printRAFASTAseq(crso->crs->cd->fp, pos, crso->clusterfp);
    return;
    }

/* END OF CLUSTERSET REPORTING CODE */
/* -------------------------------- */

/* ----------------------------- */
/* START OF REPORT ON BREAK CODE */

static BOOLEAN local_reportonbreakflag = FALSE;

static void setreportonbreakflag(int sig){
    if(local_reportonbreakflag)
        exit(1); /* STOP IF TWICE */
    errmsg(ERROR_INFO, "Will stop at end of current wordscore");
    local_reportonbreakflag = TRUE;
    return;
    }

/* END OF REPORT ON BREAK CODE */
/* --------------------------- */

/* ------------------------------ */
/* START OF SINGLETON FILTER CODE */

#define ESTCLUSTERINDEX_CHUNKSIZE 1024

static long *indexFastaWithFilter(char *path,
                                  int *seqtotal, int joinscore,
                                  SUBMAT *submat){
    register READFASTA *rf = newREADFASTA(path);
    register int total = 0, idxalloc = ESTCLUSTERINDEX_CHUNKSIZE;
    register long *idx = malloc(sizeof(long)*idxalloc);
    /* idx = getRAFASTAindices(fopen(path, "r"), NULL); */
    for(*seqtotal = 0; nextREADFASTA(rf); (*seqtotal)++){
        if(*seqtotal >= idxalloc){
            idxalloc <<=1;
            idx = realloc(idx, sizeof(long)*idxalloc);
            }
        if(maxScoreLimit(submat, (unsigned char*)rf->seq) < joinscore){
            idx[*seqtotal] = -READFASTA2RAFASTApos(rf);
            total++;
        } else {
            idx[*seqtotal] = READFASTA2RAFASTApos(rf);
            }
        }
    idx = realloc(idx, sizeof(long)*(*seqtotal));
    errmsg(ERROR_WARNING,
                "Marked %d members as low quality singletons", total);
    freeREADFASTA(rf);
    return idx;
    }

/* END OF SINGLETON FILTER CODE */
/* ---------------------------- */

typedef struct {
    int success;
    int total;
    } SCOREHISTORY;


/* THIS SECTION *REALLY* NEEDS TIDYING */

static void cluster_database(char *avapath, char *dbpath,
            FILE *dbfp, char *submatpath,
            int anytimecount, int historycount, int joinscore,
            BOOLEAN thinkaloud, BOOLEAN showalignments,
            BOOLEAN reportonbreak, BOOLEAN showhistogram,
            char *outputdir, BOOLEAN revcomp,
            char *avapathrc, int gapopen, int gapextend){
    int seqtotal;
    register SUBMAT *submat = newSUBMAT(submatpath);
    register long *rfi = indexFastaWithFilter(dbpath, &seqtotal,
                                              joinscore, submat);
    register ALLVSALL *ava = openALLVSALL(avapath),
                      *avarc = NULL;
    register CLUSTERSET *cs = newCLUSTERSET(seqtotal);
    register int i, tmp;
    register int apos, bpos, currscore, succcomp = 0, currcomp = 0;
    register LIST *history = newLIST();
    register SCOREHISTORY *sh, shvolume = {0,0};
    register float anytimeratio = 0.0;
    CLUSTERSETREPORTFUNCS crfo = {
                                     initreportclustersoutput,
                                singletonreportclustersoutput,
                                    startreportclustersoutput,
                                   finishreportclustersoutput,
                                   memberreportclustersoutput
                                };
    CLUSTERSETREPORTSTATUS crs;
    CLUSTERSETREPORTSTATUSOUTPUT crso;
    ALLVSALLRESULT avaresult;
    CLUSTERDATA cd;
    cd.fp = dbfp;
    cd.index = rfi;
    cd.compcount = 0;
    cd.submat = submat;
    cd.joinscore = joinscore;
    cd.thinkaloud = thinkaloud;
    cd.showalignments = showalignments;
    cd.rctype = 1;
    cd.gapopen   = gapopen;
    cd.gapextend = gapextend;
    currscore = viewTopALLVSALL(ava);
    if(historycount && anytimecount)
        anytimeratio = (float)historycount/(float)anytimecount;
    if(revcomp){
        if(!avapathrc)
            errmsg(ERROR_FATAL, "Avapath for revcomp required");
        avarc = openALLVSALL(avapathrc);
        tmp = viewTopALLVSALL(avarc);
        if(currscore < tmp)
            currscore = tmp;
        }
    if(reportonbreak)
        signal(SIGINT, &setreportonbreakflag);
    if(outputdir){
        crso.outputdir = outputdir;
        crso.crs = &crs;
        if(mkdir(outputdir, FILE_ALLOW_U_W_UGO_RX) == 0){
            /* OLD SGI LIBRARY BUG REQUIRES CHMOD */
            chmod(outputdir, FILE_ALLOW_U_W_UGO_RX);
        } else {
            errmsg(ERROR_FATAL, "Could not create output dir [%s]\n", 
                   outputdir);
            }
    } else {
        errmsg(ERROR_WARNING, "No output directory specified;"
                              " sequences will not be written");
        }
#ifdef SHOW_PERFORMANCE
    cd.lastchange = 0;
    cd.cs = cs;
#endif /* SHOW_PERFORMANCE */
    crs.cd = &cd; /* PROVIDE CLUSTERDATA FOR CLUSTERSETREPORTSTATUS */
    errmsg(ERROR_INFO, "Database contains %d sequences", seqtotal);
    if(revcomp){
        while((cd.rctype = nextALLVSALLpair(ava, avarc, &avaresult))){
            if(avaresult.score < currscore){ /* CHANGE OF SCORE */
                currcomp -= cd.compcount;
                currcomp = -currcomp;
                errmsg(ERROR_INFO, "After %d/%d comparisons,"
                                   " changing score %d->%d",
                        succcomp, currcomp,
                        currscore, avaresult.score);
                if(historycount){ /* IF ANYTIME EXIT ALLOWED */
                    /* ADD CURRENT TO HISTORY */
                    shvolume.total += currcomp;
                    shvolume.success += succcomp;
                    sh = NEW(SCOREHISTORY);
                    sh->total = currcomp;
                    sh->success = succcomp;
                    queueLIST(history, sh);
                    /* REMOVE PREV FROM HISTORY */
                    while(occupiedLIST(history)
                        &&((shvolume.total
               -(((SCOREHISTORY*)viewEndLIST(history))->total))
                         >= historycount)){
                        sh = popLIST(history);
                        shvolume.total -= sh->total;
                        shvolume.success -= sh->success;
                        free(sh);
                        }
                    errmsg(ERROR_INFO, "Currently %d/%d in history",
                           shvolume.success, shvolume.total);
                    if(shvolume.total > historycount){
                        if((shvolume.success == 0)
                        || (anytimeratio < ((float)shvolume.total
                                           /(float)shvolume.success))){
                             errmsg(ERROR_INFO,
                                    "Anytime exit: %d/%d in history"
                                    " scoring above %d",
                                    shvolume.success,
                                    shvolume.total, joinscore);
                             break;
                             }
                         }
                     }
                if(showhistogram && succcomp){
                    statusCLUSTERSET(cs, stdout);
                    histogramCLUSTERSET(cs, stdout);
                    printf("--\n\n");
                    }
                if(local_reportonbreakflag)
                    break;
                currscore = avaresult.score;
                succcomp = 0;
                currcomp = cd.compcount;
                }
#ifdef SHOW_PERFORMANCE
            cd.wordscore = avaresult.score;
#endif /* SHOW_PERFORMANCE */
            /* printf("Result %d vs %d wordscore=%d\n", */
                 /* avaresult.from, avaresult.to, avaresult.score); */
            apos = avaresult.from-1;
            bpos = avaresult.to-1;
     /* printf("Popped: %d %d %d %d %d\n", */
     /* avaresult.score, apos, bpos, rfi[apos], rfi[bpos]); */
            if((rfi[apos] >= 0) && (rfi[bpos] >= 0))
                if(offersubmitCLUSTERSET(cs, apos, bpos,
                                 &cd, compare_EST_sequences) > 0)
                    succcomp++;
            }
    } else {
        errmsg(ERROR_FATAL, "Sorry - non-revcomp code under revision");
        }
    statusCLUSTERSET(cs, stdout);
    histogramCLUSTERSET(cs, stdout);
    if(outputdir){
        for(i = 0; i < seqtotal; i++)
            if(rfi[i] < 0)
                rfi[i] = -rfi[i];
        reportCLUSTERSET(cs, &crso, &crfo);
        fclose(crso.singletonfp);
        }
    errmsg(ERROR_INFO, "Pairwise comparisons made: %d\n", cd.compcount);
    freeLIST(history);
    free(rfi);
    freeALLVSALL(ava);
    if(revcomp)
        freeALLVSALL(avarc);
    freeCLUSTERSET(cs);
    freeSUBMAT(submat);
    return;
    }

#define ARGUMENT_COUNT 15

int estateMAIN(){
    static FILE *dbfp  = (FILE*)0;
    static char *dbpath = NULL;
    static char *avapath = NULL, *avapathrc = NULL;
    static char *submat = "nucleic";
    static int anytimecount = 5;
    static int historycount = 5000;
    static BOOLEAN thinkaloud = FALSE;
    static BOOLEAN showalignments = FALSE;
    static BOOLEAN reportonbreak = TRUE;
    static BOOLEAN showhistogram = TRUE;
    static int joinscore = 750;
    static int gapopen = -12;
    static int gapextend = -8;
    static char *outputdir = NULL;
    static BOOLEAN revcomp = TRUE;
    ARGUMENT argdata[ARGUMENT_COUNT] = {
    { ARG_FILE,  'd', "database", "fasta database path",
      &dbfp, &dbpath, "ESTCLUSTER_DATABASE_PATH", "r", TRUE},
    { ARG_STRING, 'A', "allvsall", "allvsall score file",
      &avapath, NULL, "ESTCLUSTER_ALLVSALL_PATH", NULL, TRUE},
    { ARG_STRING, 'R', "allvsallrc", "allvsall revcomp score file",
      &avapathrc, NULL, "ESTCLUSTER_ALLVSALL_REVCOMP_PATH",
                                                         NULL, FALSE},
    { ARG_STRING, 's', "submat", "substitution matrix",
      &submat, NULL, "ESTCLUSTER_SUBMAT", NULL, FALSE},
    { ARG_INT, 'C',   "anytimecount", "anytime exit count",
      &anytimecount, NULL, "ESTCLUSTER_ANYTIME_COUNT", NULL, FALSE},
    { ARG_INT, 'T',   "historycount", "history size count",
      &historycount, NULL, "ESTCLUSTER_HISTORY_COUNT", NULL, FALSE},
    { ARG_BOOLEAN, 't',   "thinkaloud", "show clustering information",
      &thinkaloud, NULL, "ESTCLUSTER_THINKALOUD", NULL, FALSE},
    { ARG_BOOLEAN, 'a',   "showalignments", "show pairwise alignments",
      &showalignments, NULL, "ESTCLUSTER_SHOWALIGNMENTS", NULL, FALSE},
    { ARG_BOOLEAN, 'H', "showhistogram",
                                    "show histogram during clustering",
      &showhistogram, NULL, "ESTCLUSTER_SHOWHISTOGRAM", NULL, FALSE},
    { ARG_INT,  'O', "gapopen", "gap open penalty",
      &gapopen, NULL, "ESTCLUSTER_GAPOPEN",  NULL, FALSE},
    { ARG_INT,  'E', "gapextend", "gap extend penalty",
      &gapextend, NULL, "ESTCLUSTER_GAPEXTEND",  NULL, FALSE},

    { ARG_BOOLEAN, 'b',   "reportonbreak", "report clusters on break",
      &reportonbreak, NULL, "ESTCLUSTER_REPORTONBREAK", NULL, FALSE},
    { ARG_INT, 'j', "jointhreshold", "cluster joining score threshold",
      &joinscore, NULL, "ESTCLUSTER_JOINTHRESHOLD", NULL, FALSE},
    { ARG_STRING, 'o', "outputdir", "sequence output directory",
      &outputdir, NULL, "ESTCLUSTER_OUTPUTDIR", NULL, FALSE},
    { ARG_BOOLEAN, 'r', "revcomp", "use reverse complement",
      &revcomp, NULL, "ESTCLUSTER_REVCOMP", NULL, FALSE}
    };
    processARGUMENT(argdata, ARGUMENT_COUNT, global_emn,
                    "rapid accurate EST clustering");
    cluster_database(avapath, dbpath, dbfp, submat,
         anytimecount, historycount, joinscore, thinkaloud,
         showalignments, reportonbreak, showhistogram,
         outputdir, revcomp, avapathrc,
         gapopen, gapextend);
    return 0;
    }

/*
*/
