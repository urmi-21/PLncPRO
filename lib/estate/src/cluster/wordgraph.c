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

#ifndef INCLUDED_WORDGRAPH_C
#define INCLUDED_WORDGRAPH_C

#include <limits.h>
#include "../general/common.h"
#include "../general/error.h"
#include "../struct/fsm.h"
#include "../parse/rafasta.h"
#include "../sequence/sequtil.h"
#include "avastore.h"
#include "../struct/vfsm.h"
#include "../struct/list.h"

ALLOW_ERROR_MESSAGES;

/* START OF  FSM-BASED WORDGRAPH CODE */

/* addValidSequenceWordFSM: ADD WORDS TO FSM f FROM seq,
                            USING ONLY [ACGT] WORDS OF wordlen.
*/
void addValidSequenceWordFSM(FSM *f, char *seq, long wordlen){
    register unsigned char *p;
    register int ctr = 0;
    for(p = (unsigned char*)seq; *p; p++){
        ctr = f->index[*p]?ctr+1:0;
        if(ctr >= wordlen) /* ONLY ADD VALID WORDS */
            submitFSM(f, p-wordlen+1, wordlen);
        }
    return;
    }

FSM *validSequenceWordFSM(char *seq, long wordlen){
    register FSM *f = BLANK(FSM);
    f->index['A'] = f->index['C'] = f->index['G'] = f->index['T'] = 1;
    initFSM(f);
    addValidSequenceWordFSM(f, seq, wordlen);
    compileFSM(f);
    return f;
    }

int fsmWORDGRAPH(char *fastapath, char *avapath,
                 int wordlen, int wordmin){
    int total;
    register int i, j, count, scoretotal = 0;
    register RAFASTAALL *rfa = readRAFASTAALL(fastapath, &total);
    register FSM *fsm;
    register ALLVSALL *ava = newALLVSALL(avapath);
    errmsg(ERROR_INFO, "Read %d sequences", total);
    for(i = 0; i < total; i++){
        if(rfa[i].len < wordlen) /* SEQUENCE TOO SHORT */
            continue;
        setqueryALLVSALL(ava, i+1);
        fsm = validSequenceWordFSM(rfa[i].seq, wordlen);
        for(j = 0; j < i; j++){
            count = countFSM(fsm, (unsigned char*)rfa[j].seq);
            if(count >= wordmin){
                scoretotal++;
                MaddresultALLVSALL(ava, j+1, count);
                }
            }
        freeFSM(fsm);
        }
    errmsg(ERROR_INFO, "All vs All generated %d / %d possible scores",
            scoretotal, ((total*total)-total)/2);
    finishALLVSALL(ava);
    infoALLVSALL(ava, stdout);
    freeRAFASTAALL(rfa);
    freeALLVSALL(ava);
    return scoretotal;
    }

int fsmWORDGRAPHrevcomp(char *fastapath, char *avapath,
                        char *avarcpath, int wordlen, int wordmin){
    int total;
    register int i, j, count, rccount, scoretotal = 0, rcscoretotal = 0;
    register RAFASTAALL *rfa = readRAFASTAALL(fastapath, &total);
    register FSM *fsm, *rcfsm;
    register ALLVSALL *ava = newALLVSALL(avapath),
                      *rcava = newALLVSALL(avarcpath);
    errmsg(ERROR_INFO, "Read %d sequences", total);
    for(i = 0; i < total; i++){
        if(rfa[i].len < wordlen) /* SEQUENCE TOO SHORT */
            continue;
        setqueryALLVSALL(ava, i+1);
        setqueryALLVSALL(rcava, i+1);
        fsm = validSequenceWordFSM(rfa[i].seq, wordlen);
        SEQUTILrevcomp((unsigned char*)rfa[i].seq, rfa[i].len);
        rcfsm = validSequenceWordFSM(rfa[i].seq, wordlen);
        SEQUTILrevcomp((unsigned char*)rfa[i].seq, rfa[i].len);
        for(j = 0; j < i; j++){
            count = countFSM(fsm, (unsigned char*)rfa[j].seq);
            if(count >= wordmin){
                scoretotal++;
                MaddresultALLVSALL(ava, j+1, count);
                }
            rccount = countFSM(rcfsm, (unsigned char*)rfa[j].seq);
            if(rccount >= wordmin){
                rcscoretotal++;
                MaddresultALLVSALL(rcava, j+1, rccount);
                }
            }
        freeFSM(fsm);
        freeFSM(rcfsm);
        }
    errmsg(ERROR_INFO, "All vs All generated %d / %d possible scores",
            scoretotal+rcscoretotal,
            ((total*total)-total));
    errmsg(ERROR_INFO, "Comprised of %d forward and %d revcomp scores",
            scoretotal, rcscoretotal);
    finishALLVSALL(ava);
    finishALLVSALL(rcava);
    infoALLVSALL(ava, stdout);
    infoALLVSALL(rcava, stdout);
    freeRAFASTAALL(rfa);
    freeALLVSALL(ava);
    freeALLVSALL(rcava);
    return scoretotal+rcscoretotal;
    }

/*   END OF  FSM-BASED WORDGRAPH CODE */


/* START OF VFSM-BASED WORDGRAPH CODE */

typedef struct VISITEDWORDLISTNODE {
struct VISITEDWORDLISTNODE *parent; /* (NULL IF LAST IN LIST)    */
                       int  seqid;  /* SEQUENCE TO EMIT          */
                      char  count;  /* MULTIPLICITY OF THIS WORD */
    } VISITEDWORDLISTNODE;

/* ---------------------------- */
/* START OF RESULT SECTION CODE */

typedef struct {
    int mailbox;  /* INDICATE SCORE STARTED      */
    int position; /* INDICATE POSITION OF SCORE  */
    int seqid;    /* IDENTITY OF THIS SEQUENCE   */
    int score;    /* SCORE FOR THIS SEQUENCE     */
    } RESULTSECTIONNODE;

typedef struct {
RESULTSECTIONNODE  *resultv; /* STRUCTURE FOR RESULTS */
              int   resulta;
            ALLVSALL  *ava;
    } RESULTSECTION;

#define RESULTSECTIONCHUNKSIZE 100

static RESULTSECTION *newRESULTSECTION(char *avapath){
    register RESULTSECTION *rs = NEW(RESULTSECTION);
    rs->resulta = RESULTSECTIONCHUNKSIZE;
    rs->resultv = calloc(rs->resulta, sizeof(RESULTSECTIONNODE));
    rs->ava = newALLVSALL(avapath);
    return rs;
    }

static void freeRESULTSECTION(RESULTSECTION *rs){
    free(rs->resultv);
    freeALLVSALL(rs->ava);
    free(rs);
    return;
    }

static void resizeRESULTSECTION(RESULTSECTION *rs, int newsize){
    register int prev;
    if(newsize >= rs->resulta){
        prev = rs->resulta;
        rs->resulta += RESULTSECTIONCHUNKSIZE;
        rs->resultv = realloc(rs->resultv,
                        sizeof(RESULTSECTIONNODE)*rs->resulta);
        memset(rs->resultv+prev, 0, sizeof(RESULTSECTIONNODE)
                                     *RESULTSECTIONCHUNKSIZE);
        }
    return;
    }

/* END OF RESULT SECTION CODE */
/* -------------------------- */

/* ------------------------- */
/* START OF WORDHISTORY CODE */

typedef struct {
    int *whv;
    int  whc;
    int  wha;
    } WORDHISTORY;

#define WORDHISTORYCHUNKSIZE 100

static WORDHISTORY *newWORDHISTORY(){
    register WORDHISTORY *wh = NEW(WORDHISTORY);
    wh->wha = WORDHISTORYCHUNKSIZE;
    wh->whc = 0;
    wh->whv = malloc(sizeof(int)*wh->wha);
    return wh;
    }

static void freeWORDHISTORY(WORDHISTORY *wh){
    free(wh->whv);
    free(wh);
    return;
    }

static void addWORDHISTORY(WORDHISTORY *wh, int state){
    if(wh->whc >= wh->wha){
        wh->wha += WORDHISTORYCHUNKSIZE;
        wh->whv = realloc(wh->whv, sizeof(int)*wh->wha);
        }
    wh->whv[wh->whc++] = state;
    return;
    }

/* END OF WORDHISTORY CODE */
/* ----------------------- */

typedef struct {
 VISITEDWORDLISTNODE **wordarray;      /* VFSM BASE ARRAY */
                 int   wordarraylen;
         WORDHISTORY  *wordhistory;
       RESULTSECTION  *result;       /* STRUCTURE FOR RESULTS */
 VISITEDWORDLISTNODE  *vwlnbuffer;   /* VISITED WORD LIST NODE BUFFER */
                 int   vwlnbufferleft;
                LIST  *vwlnbufferhistory;
                 int   wordmin;
    } VFSMCLUSTER;

#define VWLNBUFFERCHUNKSIZE 1024

static VISITEDWORDLISTNODE *newVISITEDWORDLISTNODE(VFSMCLUSTER *vc){
    if(!vc->vwlnbufferleft){
        vc->vwlnbufferleft = VWLNBUFFERCHUNKSIZE;
        vc->vwlnbuffer = malloc(sizeof(VISITEDWORDLISTNODE)
                            *vc->vwlnbufferleft);
        MqueueLIST(vc->vwlnbufferhistory, vc->vwlnbuffer);
        }
    return &vc->vwlnbuffer[--vc->vwlnbufferleft];
    }

static VFSMCLUSTER *newVFSMCLUSTER(char *avapath, int wordmin){
    register VFSMCLUSTER *vc = NEW(VFSMCLUSTER);
    vc->wordhistory = newWORDHISTORY();
    vc->result = newRESULTSECTION(avapath);
    vc->vwlnbufferleft = 0;
    vc->vwlnbufferhistory = newLIST();
    vc->wordmin = wordmin;
    return vc;
    }

static void freeVFSMCLUSTER(VFSMCLUSTER *vc){
    freeLISTptr(vc->vwlnbufferhistory, free);
    free(vc->wordarray);
    freeWORDHISTORY(vc->wordhistory);
    freeRESULTSECTION(vc->result);
    free(vc);
    return;
    }

/* ----------------------- */
/* START OF VFSM FUNCTIONS */

static void initVFSMCLUSTER(void *info, int numstates){
    register VFSMCLUSTER *vc = (VFSMCLUSTER*)info;
    vc->wordarray = calloc(numstates, sizeof(VISITEDWORDLISTNODE*));
    vc->wordarraylen = numstates;
    if(!vc->wordarray)
        errmsg(ERROR_FATAL, "Failed to make wordarray [%d]",
                              numstates);
    return;
    }

/* startVFSMCLUSTER: START OF SEQUENCE
*/
static void startVFSMCLUSTER(void *info, int seqid, int pos){
    register VFSMCLUSTER *vc = info;
    resizeRESULTSECTION(vc->result, seqid);
    vc->wordhistory->whc = 0;
    return;
    }

static void ascendVFSMCLUSTER(VISITEDWORDLISTNODE *vwln,
                                 RESULTSECTION *rs, int seqid,
                                 int homecount, int *prepresultctr){
    register int pos;
    register RESULTSECTIONNODE *rsn;
    if(homecount == 1)
        while(vwln){
            pos = vwln->seqid;
            rsn = &rs->resultv[pos];
            if(rsn->mailbox == seqid){ /* ALREADY STARTED */
                rs->resultv[rsn->position].score += vwln->count;
            } else { /* NO RESULTS YET */
                rsn->mailbox = seqid;
                rsn->position = *prepresultctr;
                rsn = &rs->resultv[(*prepresultctr)++];
                rsn->seqid = pos;
                rsn->score = vwln->count;
                }
            vwln = vwln->parent;
            }
    else
        while(vwln){
            pos = vwln->seqid;
            rsn = &rs->resultv[pos];
            if(rsn->mailbox == seqid){ /* ALREADY STARTED */
                rs->resultv[rsn->position].score
                                 += homecount*vwln->count;
            } else { /* NO RESULTS YET */
                rsn->mailbox = seqid;
                rsn->position = *prepresultctr;
                rsn = &rs->resultv[(*prepresultctr)++];
                rsn->seqid = pos;
                rsn->score = homecount*vwln->count;
                }
            vwln = vwln->parent;
            }
    return;
    }

/* finishVFSMCLUSTER: AT END OF SEQUENCE
*/
static void finishVFSMCLUSTER(void *info, int seqid){
    register VFSMCLUSTER *vc = info;
    register int i;
    int prepresultctr = 0;
    register VISITEDWORDLISTNODE *vwln;
    for(i = 0; i < vc->wordhistory->whc; i++){ /* FOR EVERY WORD SEEN */
        vwln = vc->wordarray[vc->wordhistory->whv[i]];
        ascendVFSMCLUSTER(vwln->parent, vc->result, seqid,
                          vwln->count, &prepresultctr);
        }
    setqueryALLVSALL(vc->result->ava, seqid);
    for(i = 0; i < prepresultctr; i++){ /* STORE RESULTS */
        if(vc->result->resultv[i].score >= vc->wordmin)
            MaddresultALLVSALL(vc->result->ava,
                           vc->result->resultv[i].seqid,
                           vc->result->resultv[i].score);
        }
    return;
    }

/* observeVFSMCLUSTER: FOR EVERY VALID WORD HIT
*/
static void observeVFSMCLUSTER(void *info, int seqid,
                               int seqpos, int stateid){
    register VFSMCLUSTER *vc = info;
    register VISITEDWORDLISTNODE **ptp = &vc->wordarray[stateid],
                                  *node = *ptp;
    if(node){ /* IF SEEN THIS WORD BEFORE */
        if(node->seqid == seqid){ /* ALREADY SEEN IN THIS SEQUENCE */
            node->count++;
            if(node->count == CHAR_MAX)
                  errmsg(ERROR_FATAL,
                         "Node count overflow at [%d] for [%d]",
                                      node->count, node->seqid);
        } else { /* NOT SEEN THIS IN SEQUENCE */
            node = newVISITEDWORDLISTNODE(vc);
            node->parent = *ptp;
            node->seqid = seqid;
            node->count = 1;
            *ptp = node;
            addWORDHISTORY(vc->wordhistory, stateid);
            }
    } else { /* NEVER SEEN THIS WORD BEFORE */
        *ptp = node = newVISITEDWORDLISTNODE(vc);
        node->parent = NULL;
        node->seqid  = seqid;
        node->count  = 1;
        }
    return;
    }

/* END OF VFSM FUNCTIONS */
/* --------------------- */

/* -------------------------------- */
/* START OF VFSMCLUSTERREVCOMP CODE */

typedef struct {
      VFSMCLUSTER *vc;
    RESULTSECTION *result;
              int *wordcount;
      WORDHISTORY *wordhistory;
    } VFSMCLUSTERREVCOMP;

static VFSMCLUSTERREVCOMP *newVFSMCLUSTERREVCOMP(VFSMCLUSTER *vc,
                                                 char *avapath){
    register VFSMCLUSTERREVCOMP *vcrc = NEW(VFSMCLUSTERREVCOMP);
    vcrc->vc = vc;
    vcrc->result = newRESULTSECTION(avapath);
    vcrc->wordhistory = newWORDHISTORY();
    return vcrc;
    }

static void freeVFSMCLUSTERREVCOMP(VFSMCLUSTERREVCOMP *vcrc){
    freeRESULTSECTION(vcrc->result);
    freeWORDHISTORY(vcrc->wordhistory);
    free(vcrc);
    return;
    }

/* END OF VFSMCLUSTERREVCOMP CODE */
/* ------------------------------ */

/* ------------------------------- */
/* START OF VFSM REVCOMP FUNCTIONS */

static void initVFSMCLUSTERrevcomp(void *info, int numstates){
    register VFSMCLUSTERREVCOMP *vcrc = info;
    vcrc->wordcount = calloc(numstates, sizeof(int));
    if(!vcrc->wordcount)
        errmsg(ERROR_FATAL, "Failed to make revcomp wordcount [%d]",
                              numstates);
    return;
    }

/* startVFSMCLUSTERrevcomp: START OF SEQUENCE
*/
static void startVFSMCLUSTERrevcomp(void *info, int seqid, int pos){
    register VFSMCLUSTERREVCOMP *vcrc = info;
    resizeRESULTSECTION(vcrc->result, seqid);
    vcrc->wordhistory->whc = 0;
    return;
    }

/* finishVFSMCLUSTER: AT END OF SEQUENCE
*/
static void finishVFSMCLUSTERrevcomp(void *info, int seqid){
    register VFSMCLUSTERREVCOMP *vcrc = info;
    register VISITEDWORDLISTNODE *vwln;
    register int i, homecount;
    int prepresultctr = 0;
    for(i = 0; i < vcrc->wordhistory->whc; i++){/* FOR EACH WORD SEEN */
        homecount = vcrc->wordcount[vcrc->wordhistory->whv[i]];
        vcrc->wordcount[vcrc->wordhistory->whv[i]] = 0;
        vwln = vcrc->vc->wordarray[vcrc->wordhistory->whv[i]];
        if(vwln && (vwln->seqid == seqid)) /* DON'T COMPARE WITH SELF */
            vwln = vwln->parent;
        if(vwln) /* TIDY HERE */
            ascendVFSMCLUSTER(vwln, vcrc->result,
                              seqid, homecount, &prepresultctr);
        }
    setqueryALLVSALL(vcrc->result->ava, seqid);
    for(i = 0; i < prepresultctr; i++){ /* STORE RESULTS */
        if(vcrc->result->resultv[i].score >= vcrc->vc->wordmin)
            MaddresultALLVSALL(vcrc->result->ava,
                           vcrc->result->resultv[i].seqid,
                           vcrc->result->resultv[i].score);
        }
    return;
    }

/* observeVFSMCLUSTERrevcomp: FOR EVERY VALID WORD HIT
*/
static void observeVFSMCLUSTERrevcomp(void *info, int seqid,
                               int seqpos, int stateid){
    register VFSMCLUSTERREVCOMP *vcrc = info;
    if(vcrc->wordcount[stateid]){ /* SEEN WORD WITH THIS SEQUENCE */
        vcrc->wordcount[stateid]++;
    } else {
        vcrc->wordcount[stateid] = 1;
        addWORDHISTORY(vcrc->wordhistory, stateid);
        }
    return;
    }

/* END OF VFSM REVCOMP FUNCTIONS */
/* ----------------------------- */

int vfsmWORDGRAPH(char *fastapath, char *avapath,
                   int wordlen, int wordmin){
    register FILE *fp = fopen(fastapath, "r");
    register VFSM *vfsm = newVFSM("ACGT", wordlen, FALSE);
    register VFSMCLUSTER *vc = newVFSMCLUSTER(avapath, wordmin);
    register int scoretotal;
    PARSEFASTAVFSMFUNCS pf = {initVFSMCLUSTER,
                              startVFSMCLUSTER,
                              finishVFSMCLUSTER,
                              observeVFSMCLUSTER, NULL};
    if(!fp)
        errmsg(ERROR_FATAL, "Could not open [%s]", fastapath);
/* wordmin ISSUE TO ADDRESS HERE */
    parsefastaVFSM(fp, vfsm, vc, &pf);
    finishALLVSALL(vc->result->ava);
    infoALLVSALL(vc->result->ava, stdout);
    scoretotal = vc->result->ava->scoretotal;
    freeVFSMCLUSTER(vc);
    freeVFSM(vfsm);
    fclose(fp);
    return scoretotal;
    }

int vfsmWORDGRAPHrevcomp(char *fastapath, char *avapath,
                         char *avarcpath, int wordlen, int wordmin){
    register FILE *fp = fopen(fastapath, "r");
    register VFSM *vfsm = newVFSM("ACGT", wordlen, FALSE);
    register VFSMCLUSTER *vc = newVFSMCLUSTER(avapath, wordmin);
    register VFSMCLUSTERREVCOMP *vcrc
                               = newVFSMCLUSTERREVCOMP(vc, avarcpath);
    register int scoretotal;
    PARSEFASTAVFSMFUNCS pf = {initVFSMCLUSTER,
                              startVFSMCLUSTER,
                              finishVFSMCLUSTER,
                              observeVFSMCLUSTER, NULL};
    PARSEFASTAVFSMFUNCS pfrc = {initVFSMCLUSTERrevcomp,
                                startVFSMCLUSTERrevcomp,
                                finishVFSMCLUSTERrevcomp,
                                observeVFSMCLUSTERrevcomp, NULL};
    if(!fp)
        errmsg(ERROR_FATAL, "Could not open [%s]", fastapath);

/* wordmin ISSUE TO ADDRESS HERE */
    parsefastaVFSMrevcomp(fp, vfsm, vc, &pf, vcrc, &pfrc);

    finishALLVSALL(vc->result->ava);
    finishALLVSALL(vcrc->result->ava);
    infoALLVSALL(vc->result->ava, stdout);
    infoALLVSALL(vcrc->result->ava, stdout);
    scoretotal =    vc->result->ava->scoretotal
               +  vcrc->result->ava->scoretotal;
    freeVFSMCLUSTER(vc);
    freeVFSMCLUSTERREVCOMP(vcrc);
    freeVFSM(vfsm);
    fclose(fp);
    return scoretotal;
    }

/*   END OF VFSM-BASED WORDGRAPH CODE */

#ifdef TEST_THIS_MODULE

int main(){
    /* register char *fastapath  */
    /* = "/people/gslater/archive/stratagene_pancreas.fasta"; */
    /* register char *fastapath = "/people/gslater/ten.fasta"; */
    /* register char *fastapath = "/people/gslater/archive/" */
                          /* "soares_infant_brain_inib.fasta"; */
    register char *fastapath = "/people/gslater/archive/testdata/"
                   /* "stratagene_mouse_diaphragm937303.fasta"; */
                   "stratagene_mouse_diaphragm937303_nseg.fasta";
                   /* "soares_fetal_liver_spleen_1nfls.fasta"; */
    /* register char *fastapath = "test.fasta"; */
    /* register char *avapath = "/people/gslater/scratch/test.ava"; */
    register char *avapath = "/dev/null";
    /* register char *avapath = "/people/gslater/archive/testdata/" */
                   /* "stratagene_mouse_diaphragm937303.ava"; */
                   /* "stratagene_mouse_diaphragm937303_nseg.ava"; */
                   /* "soares_fetal_liver_spleen_1nfls.ava"; */
    register int wordlen = 12;
    register int total =
        /* fsmWORDGRAPH(fastapath, avapath, wordlen); */
        vfsmWORDGRAPH(fastapath, avapath, wordlen);
    errmsg(ERROR_INFO, "Total is %d", total);
    /* remove(avapath); */
    return 0;
    }

#endif /* TEST_THIS_MODULE */

#endif /* INCLUDED_WORDGRAPH_C */

/*
    DONE:revcomp

    Memory saving options:
        [1] Use dual pass. (max req. still similar)
        [2] Write some to disk. (slow).
        [3] Convert long lists to arrays.
            (put freed nodes onto a buffer queue)
*/

