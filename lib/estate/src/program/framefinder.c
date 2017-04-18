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

/* framefinder : error tolerant DNA sequence framefinding.
*/

#include <ctype.h> /* for toupper() */

#include "../general/arg.h"
#include "../general/error.h"
#include "../general/common.h"
#include "../parse/readfasta.h"
#include "../sequence/ntaa.h"
#include "../parse/parser.h"
#include "../struct/vfsm.h"
#include "../sequence/sequtil.h"

ALLOW_ERROR_MESSAGES;

/* -------------------------------------- */
/* START OF WORD PROBABILITY PARSING CODE */

typedef struct {
       VFSM *vfsm;
     double *prob;
    } FRAMEFINDER_WORDPROB;

static FRAMEFINDER_WORDPROB *newFRAMEFINDER_WORDPROB(FILE *fp){
    register FRAMEFINDER_WORDPROB *ffwp = NEW(FRAMEFINDER_WORDPROB);
    register PARSER *p = newPARSER(fp);
    register int count, wordlen;
    do { count = wordPARSER(p);
         if(count == EOF)
             errmsg(ERROR_FATAL, "No words in input");
    } while(count != 2);
    wordlen = strlen((char*)p->word[1]);
    ffwp->vfsm = newVFSM("ACGTN", wordlen, FALSE);
    ffwp->prob = calloc(ffwp->vfsm->lrw, sizeof(double));
    if(!ffwp->prob)
        errmsg(ERROR_FATAL, "Insufficient memory for vfsm");
    do {
       if(count == 2)
          ffwp->prob[word2posVFSM(ffwp->vfsm, p->word[1])]
              = strtod((char*)p->word[0], NULL);
    } while((count = wordPARSER(p)) != EOF);
    freePARSER(p);
    return ffwp;
    }

static void freeFRAMEFINDER_WORDPROB(FRAMEFINDER_WORDPROB *ffwp){
    freeVFSM(ffwp->vfsm);
    free(ffwp->prob);
    free(ffwp);
    return;
    }

/* END OF WORD PROBABILITY PARSING CODE */
/* ------------------------------------ */

/* ------------------------------------ */
/* START OF FRAMEFINDER_RESULT_SET CODE */

typedef struct {
    char *aaseq;  /* TRANSLATED SEQUENCE */
     int  aalen;  /* TRANSLATION LENGTH  */
     int  start;  /* START POSITION      */
     int  finish; /* FINISH POSITION     */
  double  score;  /* RESULT SCORE        */
    } FRAMEFINDER_RESULT;

#define FRAMEFINDER_RESULT_TYPE_TOTAL           4
#define FRAMEFINDER_RESULT_REVCOMP          (1<<0)
#define FRAMEFINDER_RESULT_LOCAL_STRICT     (1<<1)
#define FRAMEFINDER_RESULT_LENGTH_CHUNKSIZE   100

#define FRAMEFINDER_RESULT_CHECK(frs, type) (frs->result[type])

typedef struct {
    FRAMEFINDER_RESULT  *result[FRAMEFINDER_RESULT_TYPE_TOTAL];
                   int   lenalloc;
                   int   mask;
    } FRAMEFINDER_RESULT_SET;

static FRAMEFINDER_RESULT_SET *newFRAMEFINDER_RESULT_SET(
                 BOOLEAN tryrevcomp, BOOLEAN trystrictlocal){
    register int i;
    register FRAMEFINDER_RESULT_SET *frs = NEW(FRAMEFINDER_RESULT_SET);
    frs->lenalloc = FRAMEFINDER_RESULT_LENGTH_CHUNKSIZE;
    frs->mask = tryrevcomp?FRAMEFINDER_RESULT_REVCOMP:0;
    if(trystrictlocal)
        frs->mask |= FRAMEFINDER_RESULT_LOCAL_STRICT;
    for(i = 0; i < FRAMEFINDER_RESULT_TYPE_TOTAL; i++){
        if((i | frs->mask) != frs->mask){
            frs->result[i] = NULL;
            continue; /* MALLOC ONLY WHEN REQUIRED */
            }
        frs->result[i] = NEW(FRAMEFINDER_RESULT);
        frs->result[i]->aaseq = malloc(sizeof(char)*frs->lenalloc);
        }
    return frs;
    }

static void freeFRAMEFINDER_RESULT_SET(FRAMEFINDER_RESULT_SET *frs){
    register int i;
    for(i = 0; i < FRAMEFINDER_RESULT_TYPE_TOTAL; i++)
        if(FRAMEFINDER_RESULT_CHECK(frs, i)){
            free(frs->result[i]->aaseq);
            free(frs->result[i]);
            }
    free(frs);
    return;
    }

static void resizeFRAMEFINDER_RESULT_SET(FRAMEFINDER_RESULT_SET *frs,
                                         int newlength){
    register int i;
    if(frs->lenalloc <= newlength){
        frs->lenalloc = newlength
                      + FRAMEFINDER_RESULT_LENGTH_CHUNKSIZE;
        for(i = 0; i < FRAMEFINDER_RESULT_TYPE_TOTAL; i++){
            if(FRAMEFINDER_RESULT_CHECK(frs, i)){
                frs->result[i]->aaseq = realloc(frs->result[i]->aaseq,
                                               frs->lenalloc);
                }
            }
        }
    return;
    }

static int printFRAMEFINDERmodel(FILE *fp, int model){
    return fprintf(fp, "{%s,%s}",
       (model & FRAMEFINDER_RESULT_REVCOMP)?"revcomp":"forward",
       (model & FRAMEFINDER_RESULT_LOCAL_STRICT)?"strict":"local");
    }

#ifdef NOTNOW
/* NOW BROKEN */
static void infoFRAMEFINDER_RESULT_SET(FRAMEFINDER_RESULT_SET *frs){
    errmsg(ERROR_INFO, "FrameFinder Result Set:\n"
                       "     result : %p\n"
                       "      total : %d\n"
                       "   lenalloc : %d\n",
                      frs->result,
                      frs->total,
                      frs->lenalloc);
    return;
    }
#endif /* NOTNOW */

/* END OF FRAMEFINDER_RESULT_SET CODE */
/* ---------------------------------- */

/* ------------------------------ */
/* START OF FRAMEFINDER NODE CODE */

#define FRAMEFINDER_NODE_CHUNKSIZE 100

typedef struct {
unsigned char base;    /* CURRENT SEQUENCE BASE                 */
         char residue; /* CURRENT SEQUENCE BASE                 */
         char trace;   /* TRACEBACK STORAGE                     */
          int cover;   /* USED FOR PROBABILITY COVERAGE SCORING */
       double prob;    /* LOG ODDS CODING WORD PROBABILITY      */
       double score;   /* FRAME CODING POTENTIAL SCORE          */
       double info;    /* LOCAL INFORMATION CONTENT             */
    } FRAMEFINDER_NODE;

typedef struct {
      FRAMEFINDER_NODE *node;       /* PROBABILITY NODES */
                   int  length;
                   int  alloc;
FRAMEFINDER_RESULT_SET *ffrs;
                double  fspen;          /* FRAMESHIFT PENALTY     */
                double  inspen;         /*  INSERTION PENALTY     */
                double  delpen;         /*   DELETION PENALTY     */
                double  openbonus;      /* BONUS FOR FRAME OPEN   */
               BOOLEAN  tryrevcomp;     /* LOOK AT REVCOMP        */
               BOOLEAN  trystrictlocal; /* TRY STRICT LOCAL FIRST */
               BOOLEAN  thinkaloud;     /* SHOW THINKING          */
    } FRAMEFINDER;

static void checkValidParamFRAMEFINDER(FRAMEFINDER *ff){
    register char *msg = " penalty should not be positive";
    if(ff->fspen > 0)
        errmsg(ERROR_FATAL, "frameshift%s", msg);
    if(ff->inspen > 0)
        errmsg(ERROR_FATAL, "insertion%s", msg);
    if(ff->delpen > 0)
        errmsg(ERROR_FATAL, "deletion%s", msg);
    if((-ff->openbonus) < (ff->fspen+Max(ff->inspen,ff->delpen)))
        errmsg(ERROR_WARNING, "Open bonus may be too high [%.2f]",
                               ff->openbonus);
    return;
    }

static FRAMEFINDER *newFRAMEFINDER(double fspen, double inspen,
       double delpen, double openbonus,
       BOOLEAN tryrevcomp, BOOLEAN trystrictlocal, BOOLEAN thinkaloud){
    register FRAMEFINDER *ff = NEW(FRAMEFINDER);
    ff->alloc = FRAMEFINDER_NODE_CHUNKSIZE;
    ff->node = malloc(sizeof(FRAMEFINDER_NODE)*ff->alloc);
    ff->ffrs = newFRAMEFINDER_RESULT_SET(tryrevcomp, trystrictlocal);
    ff->length = 0;
    ff->fspen = fspen;
    ff->inspen = inspen;
    ff->delpen = delpen;
    ff->openbonus = openbonus;
    ff->tryrevcomp = tryrevcomp;
    ff->trystrictlocal = trystrictlocal;
    ff->thinkaloud = thinkaloud;
    checkValidParamFRAMEFINDER(ff);
    return ff;
    }

static void freeFRAMEFINDER(FRAMEFINDER *ff){
    freeFRAMEFINDER_RESULT_SET(ff->ffrs);
    free(ff->node);
    free(ff);
    return;
    }

static void resizeFRAMEFINDER(FRAMEFINDER *ff, int length){
    ff->length = length;
    if(ff->alloc <= length){
        ff->alloc  = length
                   + FRAMEFINDER_NODE_CHUNKSIZE;
        ff->node = realloc(ff->node,
                           sizeof(FRAMEFINDER_NODE)*ff->alloc);
        resizeFRAMEFINDER_RESULT_SET(ff->ffrs, (length/3)+1);
        }
    return;
    }

#ifdef NOTNOW
static void infoFRAMEFINDER(FRAMEFINDER *ff){
    errmsg(ERROR_INFO, "FrameFinder :\n"
                       "       node : %p\n"
                       "     length : %d\n"
                       "      alloc : %d\n"
                       "       ffrs : %p\n"
                       "      fspen : %.2f\n"
                       "     inspen : %.2f\n"
                       "     delpen : %.2f\n"
                       "  tryrevcomp : %s\n"
                       " thinkaloud : %s\n",
                      ff->node,
                      ff->length,
                      ff->alloc,
                      ff->ffrs,
                      ff->fspen,
                      ff->inspen,
                      ff->delpen,
                      ff->tryrevcomp?"TRUE":"FALSE",
                      ff->thinkaloud?"TRUE":"FALSE");
    infoFRAMEFINDER_RESULT_SET(ff->ffrs);
    return;
    }
#endif /* NOTNOW */

/* END OF FRAMEFINDER NODE CODE */
/* ---------------------------- */

/*
NNNNNNACGTACGTNNNNNN
01234567890123456789
012345        543210
21x12          21x12
321x123      321x123
*/

/*
 FW:L:SS
 FW:L:CO
 FW:G
 RC:L:SS
 RC:L:CO
 RC:G
*/

static void ClearFrameData(FRAMEFINDER *ff){
    register int i;
    for(i = 0; i < ff->length; i++){
        ff->node[i].trace = '\0';
        ff->node[i].cover =   0;
        ff->node[i].prob  =   0.0;
        ff->node[i].score =   0.0;
        ff->node[i].info  =   0.0;
        }
    return;
    }

static void SetFrameProbabilities(FRAMEFINDER *ff,
                               FRAMEFINDER_WORDPROB *ffwp){
    register int i, j, currframe, inbase;
    register VFSMTYPE state = 0;
    register double setprob;
    for(i = 0; i < ff->length; i++){ /* SET PROBABILITIES */
        inbase = toupper(ff->node[i].base);
        if(!SEQUTILisacgt(inbase))
            inbase = 'N';
        state = MchangeVFSMstate(ffwp->vfsm, state, inbase);
        if(MisleafVFSMstate(ffwp->vfsm, state)){
            currframe = (i+1)%3;
            setprob = ffwp->prob[Mstate2posVFSM(ffwp->vfsm, state)];
            for(j = 0; j < ffwp->vfsm->depth; j+=3){
                ff->node[i-j-2].prob += setprob;
                ff->node[i-j-2].cover++;
                }
            }
        }
    for(i = 0; i < ff->length; i++){
        if((ff->node[i].prob) && (ff->node[i].cover > 1))
            ff->node[i].prob  /=  ff->node[i].cover;
            }
    return;
    }

static void SetFrameTranslation(FRAMEFINDER *ff){
    register int i;
    for(i = 0; i < ff->length; i++){
        ff->node[i].residue = NTAATranslateBase(ff->node[i].base,
                                                ff->node[i+1].base,
                                                ff->node[i+2].base);
        }
    return;
    }

/* SetFrameInfo :
   ASSUMES MINIMUM 3(N) (LEADING AND TRAILING) HAVE BEEN ADDED.
*/
static void SetFrameInfo(FRAMEFINDER *ff){
    register int wlen = 5, start = 2, stop = ff->length-2;
    register int i, ninwin = wlen-1;
    register double frac = 1.0/(double)wlen;
    /* printf("Start at [%d] stop at [%d] (len=%d,wlen=%d)\n",  */
        /* start, stop, ff->length, wlen); */
    for(i = start; i < stop; i++){
        if(ff->node[i+start].base == 'N')
            ninwin++;
        /* printf("oh yea [%c] with [%d] [%8.4f]\n",  */
            /* ff->node[i].base, ninwin, (wlen-ninwin)*frac); */
        ff->node[i].info = (wlen-ninwin)*frac;
        if(ff->node[i-start].base == 'N')
            ninwin--;
        }
    return;
    }

static void ShowFramefinderThinking(FRAMEFINDER *ff,
                                    FRAMEFINDER_RESULT *ffr,
                                    int model){
    register int i, j;
    printf("-----\nFramefinder Thinking");
    printFRAMEFINDERmodel(stdout, model);
    printf("\n");
    printf("Start at %d, finish at %d\n", ffr->start, ffr->finish);
    printf("Maxscore %.9f\n", ffr->score);
    for(i = 0; i < ff->length; i++){
        printf("%c %3d %c [%c%c%c] prob=%8.4f score=%8.4f tb [",
                ((i >= ffr->start) && (i <= ffr->finish))?'>':'#',
                i, ff->node[i].base,
                ((i+0)%3)?' ':ff->node[i].residue,
                ((i+2)%3)?' ':ff->node[i].residue,
                ((i+1)%3)?' ':ff->node[i].residue,
                ff->node[i].prob, ff->node[i].score);
        for(j = 0; j < 6; j++){
            printf("%c", GetBit(ff->node[i].trace, j)?('0'+j):'-');
            }
        printf("] info=%8.4f\n", ff->node[i].info);
        }
    printf("--\n");
    return;
    }

#define FrameChallenge(pos, score, cand, back)    \
    if(score <= cand){                            \
        if(score == cand){                        \
           SetBitOn(ff->node[pos].trace, back);   \
        } else {                                  \
           score = cand;                          \
           SetBitOnly(ff->node[pos].trace, back); \
           }                                      \
       }

static int FindTraceContrib(int trace){
    static int pref[6] = {3,4,2,5,1,0};
    register int i;
    for(i = 0; i < 6; i++)
        if(GetBit(trace, pref[i]))
            return pref[i];
    errmsg(ERROR_FATAL, "No traceback set");
    return -1;
    }

/*
 FW:L:SS  REVCOMP|LOCAL_START_STOP
 FW:L:CO  REVCOMP|LOCAL_CONTENT
 FW:G     REVCOMP|GLOBAL
*/

static double FindBestFrameLocalStrict(FRAMEFINDER *ff,
              FRAMEFINDER_WORDPROB *ffwp, FRAMEFINDER_RESULT *ffr){
    register int i, j, k, start, maxscorepos = -1;
    register double score, cand, maxscore = 0.0, penalty;
    for(i = 0; i < ff->length-2; i++){
        score = 0.0;
        switch(i){ /* SKIP SOME LOOKBACK FOR INITIAL BASES */
            default:
                penalty = (ff->fspen + (2*ff->delpen)) 
                        * ff->node[i].info;
                cand = ff->node[i-5].score + penalty;
                FrameChallenge(i, score, cand, 5);
            case 4:
                penalty = (ff->fspen + ff->delpen) * ff->node[i].info;
                cand = ff->node[i-4].score + penalty;
                FrameChallenge(i, score, cand, 4);
            case 3:
                cand = ff->node[i-3].score;
                FrameChallenge(i, score, cand, 3);
            case 2:
                penalty = (ff->fspen + ff->inspen) * ff->node[i].info;
                cand = ff->node[i-2].score + penalty;
                FrameChallenge(i, score, cand, 2);
            case 1:
                penalty = (ff->fspen + (2*ff->inspen)) 
                        * ff->node[i].info;
                cand = ff->node[i-1].score + penalty;
                FrameChallenge(i, score, cand, 1);
            case 0:
                /* FORCE START */
                /* if((score > 0.0) || (ff->node[i].residue == 'M')) */
                    /* score += ff->node[i].prob; */
                /* else */
                    /* ff->node[i].trace = 0; */
                if(score > 0.0)
                    score += ff->node[i].prob;
                else if (ff->node[i].residue == 'M') /* OPEN BONUS */
                    score += ff->node[i].prob + ff->openbonus;
                else
                    ff->node[i].trace = 0;
                FrameChallenge(i, score, 0.0, 0);
                ff->node[i].score = score;
                if(ff->node[i+3].residue == '*'){ /* FORCE STOP */
                    if(maxscore <= score){  /* MAXIMISE LENGTH */
                        maxscore = score;
                        maxscorepos = i;
                        }
                    ff->node[i].score = 0.0;
                    }
                break;
            }
        }
    for(start = maxscorepos, k = 0; start > 0; k++){
        j = FindTraceContrib(ff->node[start].trace);
        ffr->aaseq[k] = ff->node[start].residue; 
        if(j == 0)
            break;
        start-=j; /* TRACEBACK */
        }
    ffr->start  = start;
    ffr->finish = maxscorepos;
    ffr->score  = maxscore;
    ffr->aalen = k;
    ffr->aaseq[k] = '\0';
    SEQUTILreverse(ffr->aaseq, k);
    return maxscore;
    }

static double FindBestFrameLocal(FRAMEFINDER *ff,
              FRAMEFINDER_WORDPROB *ffwp, FRAMEFINDER_RESULT *ffr){
    register int i, j, k, start, maxscorepos = -1;
    register double score, cand, maxscore = 0.0, penalty;
    for(i = 0; i < ff->length-2; i++){
        score = 0.0;
        switch(i){ /* SKIP SOME LOOKBACK FOR INITIAL BASES */
            default:
                penalty = (ff->fspen + (2*ff->delpen))
                        * ff->node[i].info;
                cand = ff->node[i-5].score + penalty;
                FrameChallenge(i, score, cand, 5);
            case 4:
                penalty = (ff->fspen + ff->delpen) * ff->node[i].info;
                cand = ff->node[i-4].score + penalty;
                FrameChallenge(i, score, cand, 4);
            case 3:
                cand = ff->node[i-3].score;
                FrameChallenge(i, score, cand, 3);
            case 2:
                penalty = (ff->fspen + ff->inspen) * ff->node[i].info;
                cand = ff->node[i-2].score + penalty;
                FrameChallenge(i, score, cand, 2);
            case 1:
                penalty = (ff->fspen
                        + (2*ff->inspen)) * ff->node[i].info;
                cand = ff->node[i-1].score + penalty;
                FrameChallenge(i, score, cand, 1);
            case 0:
                score += ff->node[i].prob;
                FrameChallenge(i, score, 0, 0);
                ff->node[i].score = score;
                if(ff->node[i+3].residue == '*') /* SKIP STOP */
                    ff->node[i].score = 0.0; /* FS ROUND STOP CODONS */
                if(maxscore <= score){ /* MAXIMISE LENGTH */
                    maxscore = score;
                    maxscorepos = i;
                    }
                break;
            }
        }
    for(start = maxscorepos, k = 0; start > 0; k++){
        j = FindTraceContrib(ff->node[start].trace);
        ffr->aaseq[k] = ff->node[start].residue; 
        if(j == 0)
            break;
        start-=j; /* TRACEBACK */
        }
    ffr->start  = start;
    ffr->finish = maxscorepos;
    ffr->score  = maxscore;
    ffr->aalen = k;
    ffr->aaseq[k] = '\0';
    SEQUTILreverse(ffr->aaseq, k);
    return maxscore;
    }

static void PrepareFrameFinder(FRAMEFINDER *ff,
                               FRAMEFINDER_WORDPROB *ffwp){
    ClearFrameData(ff);
    SetFrameProbabilities(ff, ffwp);
    SetFrameTranslation(ff);
    SetFrameInfo(ff);
    return;
    }

static void PrintBestFrame(FILE *fp, FRAMEFINDER_RESULT *ffr,
                           READFASTA *rf, int model){
    register float used = ((float)(ffr->aalen*3))
                          /((float)rf->seqlen)*100.0;
    printf(">%s [framefinder (%d,%d) score=%.2f used=%.2f%% ",
          rf->def, ffr->start, ffr->finish, ffr->score, used);
    printFRAMEFINDERmodel(stdout, model);
    printf(" ]\n");
    SEQUTILwriteFASTAblock(stdout, ffr->aaseq, ffr->aalen);
    return;
    }

static void SetFrameRevComp(FRAMEFINDER *ff){
    register int a, z;
    register unsigned char swap;
    for(a = 0, z = ff->length-1; a < z; a++, z--){
        swap = SEQUTILcomplement(ff->node[a].base);
        ff->node[a].base = SEQUTILcomplement(ff->node[z].base);
        ff->node[z].base = swap;
        }
    if(ff->length&1){ /* IF ODD, COMPLEMENT THE CENTRAL BASE */
        a = ff->length>>1;
        ff->node[a].base = SEQUTILcomplement(ff->node[a].base);
        }
    return;
    }

/* CHANGE HERE AT SOME POINT,
   SO THAT ALL REVCOMP TRANSLATIONS
   ARE GENERATED TOGETHER
*/
static void MakeCandidateTranslation(FRAMEFINDER *ff, 
                                    FRAMEFINDER_WORDPROB *ffwp){
    register int i;
    register FRAMEFINDER_RESULT *ffr;
    for(i = 0; i < FRAMEFINDER_RESULT_TYPE_TOTAL; i++){
        ffr = FRAMEFINDER_RESULT_CHECK(ff->ffrs, i);
        if(ffr){ /* GENERATE RESULT */
            if(i & FRAMEFINDER_RESULT_REVCOMP)
                SetFrameRevComp(ff);
            PrepareFrameFinder(ff, ffwp); 
            if(i & FRAMEFINDER_RESULT_LOCAL_STRICT)
                FindBestFrameLocalStrict(ff, ffwp, ffr);
            else
                FindBestFrameLocal(ff, ffwp, ffr);
            if(ff->thinkaloud)
                ShowFramefinderThinking(ff, ffr, i);
            if(i & FRAMEFINDER_RESULT_REVCOMP)
                SetFrameRevComp(ff); /* RESET TO FORWARD */
            }
        }
    return;
    }

static void ChooseBestTranslation(FRAMEFINDER *ff, READFASTA *rf){
    register int i, choice = -1;
    register double choicescore = 0.0;
    register FRAMEFINDER_RESULT *ffr;
    for(i = 0; i < FRAMEFINDER_RESULT_TYPE_TOTAL; i++){
        ffr = FRAMEFINDER_RESULT_CHECK(ff->ffrs, i);
        if(ffr){
            choice = i;
            choicescore = ffr->score;
            break;
            }
        }
    while(i < FRAMEFINDER_RESULT_TYPE_TOTAL){
        ffr = FRAMEFINDER_RESULT_CHECK(ff->ffrs, i);
        if((ffr) && (ffr->score > choicescore)){
            choice = i;
            choicescore = ffr->score;
            }
        i++;
        }
    PrintBestFrame(stdout, ff->ffrs->result[choice],  rf, choice);
    return;
    }

static void FrameFind(FRAMEFINDER *ff, FILE *wordprobfp, char *dbpath){
    register FRAMEFINDER_WORDPROB *ffwp 
              = newFRAMEFINDER_WORDPROB(wordprobfp);
    register READFASTA *rf = newREADFASTA(dbpath);
    register int i, fflen;
    while(nextREADFASTA(rf)){
        fflen = rf->seqlen+(ffwp->vfsm->depth*2);
        resizeFRAMEFINDER(ff, fflen);
        /* infoFRAMEFINDER(ff); */
        for(i = 0; i < ffwp->vfsm->depth; i++) /* SET START */
            ff->node[i].base = 'N';
        for(i = 0; i < rf->seqlen; i++) /* COPY SEQUENCE */
            ff->node[i+ffwp->vfsm->depth].base = rf->seq[i];
        for(i = 0; i < ffwp->vfsm->depth; i++) /* SET END */
            ff->node[rf->seqlen+ffwp->vfsm->depth+i].base = 'N';
        MakeCandidateTranslation(ff, ffwp);
        /* ChooseBestTranslation(ff, ffwp, rf); */
        ChooseBestTranslation(ff, rf);
        }
    freeFRAMEFINDER_WORDPROB(ffwp);
    freeREADFASTA(rf);
    return;
    }

#define ARGUMENT_COUNT 9

int estateMAIN(){
    static char *dbpath = NULL;
    static FILE *wordprobfp;
    static char *wordprobpath = NULL;
    static BOOLEAN thinkaloud = FALSE;
    static BOOLEAN tryrevcomp  = TRUE;
    static BOOLEAN trystrictlocal = TRUE;
    static double  fspen = -24.0;
    static double inspen =  -2.0;
    static double delpen =  -2.0;
    static double openbonus =  14.0;
    ARGUMENT argdata[ARGUMENT_COUNT] = {
    { ARG_FILE,    'w', "wordprob", "word probabilities",
     &wordprobfp, &wordprobpath, "FRAMEFINDER_WORDPROB", "r", TRUE},
    { ARG_STRING,  'd', "database", "fasta format DNA sequences",
     &dbpath, NULL, "FRAMEFINDER_DATABASE", NULL, TRUE},
    { ARG_BOOLEAN, 't', "thinkaloud", "show calculation information",
     &thinkaloud, NULL, "FRAMEFINDER_THINKALOUD", NULL, FALSE},
    { ARG_BOOLEAN, 'r', "tryrevcomp", "try reverse complement",
     &tryrevcomp, NULL, "FRAMEFINDER_TRYREVCOMP", NULL, FALSE},
    { ARG_BOOLEAN, 's', "trystrictlocal", "try strict local model",
     &trystrictlocal, NULL, "FRAMEFINDER_TRYSTRICTLOCAL", NULL, FALSE},
    { ARG_DOUBLE,  'F', "fspen", "frameshift penalty",
     &fspen, NULL, "FRAMEFINDER_FSPEN", NULL, FALSE},
    { ARG_DOUBLE,  'I', "inspen", "insertion penalty",
     &inspen, NULL, "FRAMEFINDER_INSPEN", NULL, FALSE},
    { ARG_DOUBLE,  'D', "delpen", "deletion penalty",
     &delpen, NULL, "FRAMEFINDER_DELPEN", NULL, FALSE},
    { ARG_DOUBLE,  'O', "openbonus", "frame open bonus",
     &openbonus, NULL, "FRAMEFINDER_OPENBONUS", NULL, FALSE}
    };
    register FRAMEFINDER *ff;
    wordprobfp = stdin;
    processARGUMENT(argdata, ARGUMENT_COUNT, global_emn,
             "framefinder : error tolerant DNA sequence framefinding");
    NTAAInitialiseData();
    ff = newFRAMEFINDER(fspen, inspen, delpen, openbonus,
                        tryrevcomp, trystrictlocal, thinkaloud);
    /* infoFRAMEFINDER(ff); */
    FrameFind(ff, wordprobfp, dbpath);
    freeFRAMEFINDER(ff);
    return 0;
    }

/* plan:
   allow both score threshold and usage threshold.

ACGTACGTACGTACGT
111222333444555.
.111222333444555
..111222333444..

ACGTACGTACGTACGTACG
111111222222333333.
.111111222222333333
..111111222222.....

  o read wordusestats using vfsm
    foreach seq
        DONE:generate frequency array using vfsm.
        DONE:blend frequencies.
        DONE:convert to probability array.
        DONE:pick optimal path through probability array.
        DONE:repeat for revcomp.
        DONE:translate optimal path from forward or revcomp.

Frame transitions:
-----------------
ACGTACGTACGTACGTA
111222333444555..
.111222333444555.
..111222333444555

Frameshift possibilities:
------------------------
curr:     1     2     3
next
    1     =  2d,1i 1d,2i
    2  1d,2i    =  2d,1i
    3  2d,1i 1d,2i    =

  NO- using hexamer freqs for not only the words in the sequence,
      but the words implied by the frameshifts under consideration.
  NO- revcomp bias on basis of poly(A) or poly(T) presence.
  NO- special algorithm for contig_proc generated files:
          ie consider units between 20(N) as rearrangeable.
  NO- Initialisation approach, so in traceback, 
      if at an M and score <= 1?, stop. 
  NO- relate frameshift penalty to average expected length
      and observed lengths of reading frame. 8->26927
  NO- Look at average length of coding sequences, and average score
      with log odds ratio score statistics.

General:
-------
    [ FrameFinder : {Fw|Rc} Region=(%d,%d) Score=%.2f Usage=%.2f ]
    - use of proper training set.
        o Obtain scores and associated data for various model
          types, and graph to see where various model types are best.
          ie. observe various models in sensitivity/selectivity space.
    - options for fshift/ins/del penalties, and boolean -t --thinkaloud
    - Introduce some sort of confidence score ?
    - clip initial and tail Xs
    - proper revcomp generation and selection.
    - force frameshift around stop codons.
    - Try setting frameshift penalty as a fraction of current score ?
    - Add function for parameter learning.
    - Try addition of pseudo-score on frameopen.
    - Add tail X trimming for non-strict models.
    - Fix position reporting in output for revcomp models.
    - Try applying a 'simulated annealing' ? approach,
      whereby repeatly apply algorithm, reducing penalties,
      deciding where to stop ?  What would this show ?
    - Fix --showthinking output to work properly

Poly(N)
-------
    - experiment with varying the frameshift penalty depending
      on the consequences - ie no penalty during a poly(N).
    - Frameshift penalty *= (current word info/max word info).
    - Use window length of 5.

Model types
-----------
    - For start: only add prob if already +ve, or if(aaseq[i] == 'M')
    - For stop: only take maxscore if(aaseq[i] == '*')
    - For global: allow scores less than zero,
                  start max[0-3], stop max[last3].
    - Use progressively less stringent models
      until a satisfactory score if found;
      [A] Forward
      [B] RevComp
        [1] Local
            (a) start & stop
            (b) content only
        [2] Joined Local (???)
        [3] Global
    - Add option to select model used.
    - Add initial pass to survey sequence quality 
      and force global/nostartstop based on model type.
    - Does joined local have time complexity of O(N^2)

TODO:
    - Rewrite dp routines to facillitate intro of global model ?
    - Point of prob score addition ?
    - Frameshift implied hexamer statistics ?
*/

