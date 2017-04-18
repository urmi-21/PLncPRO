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

#ifndef INCLUDED_VFSM_C
#define INCLUDED_VFSM_C

#include "../general/error.h"
#include "../struct/vfsm.h"
#include "../sequence/sequtil.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* FOR memset()     */
#include <strings.h> /* FOR ffs()        */
#include <ctype.h>   /* FOR tolower()    */
#include <math.h>    /* FOR pow()        */

ALLOW_ERROR_MESSAGES;

VFSM *newVFSM(char *alphabet, int depth, BOOLEAN casesensitive){
    register VFSM *vfsm = NEW(VFSM);
    register int i;
    register double maxcheck;
    register VFSMTYPE tmp;
    vfsm->alphabet  = (unsigned char*)strdup(alphabet);
    vfsm->alphasize = strlen(alphabet);
    vfsm->depth     = depth;
    if(vfsm->depth < 3)
        errmsg(ERROR_FATAL, "Minimum VFSM depth = 3");
    maxcheck = vfsm->alphasize
          * ((pow(vfsm->alphasize, depth)-1)/(vfsm->alphasize-1));
    if(maxcheck >= VFSMTYPE_MAX)
        errmsg(ERROR_FATAL, "Too many states for VFSM [%e]", maxcheck);
    for(i = 3, tmp = vfsm->prs = vfsm->alphasize; i < depth; i++){
        tmp*=vfsm->alphasize;
        vfsm->prs += tmp;
        }
    vfsm->prs++;
    vfsm->lrs   = vfsm->prs+(tmp*vfsm->alphasize);
    vfsm->prw   = vfsm->lrs-vfsm->prs;
    vfsm->lrw   = vfsm->prw*vfsm->alphasize;
    vfsm->total = vfsm->lrw+vfsm->lrs;
    if(NotSingleBit(vfsm->alphasize))
        vfsm->ispoweroftwo = vfsm->logalphasize = 0;
    else {
        vfsm->ispoweroftwo = vfsm->prw-1;
        vfsm->logalphasize = ffs(vfsm->alphasize)-1;
        } /* IF alphasize IS POWER OF 2, SO WILL BE vfsm->prw */
    memset(vfsm->index, 0, sizeof(char)*ALPHABETSIZE);
    if(casesensitive)
        for(i = 0; i < vfsm->alphasize; i++)
            vfsm->index[toupper(vfsm->alphabet[i])] =
            vfsm->index[tolower(vfsm->alphabet[i])] = i+1;
    else
        for(i = 0; i < vfsm->alphasize; i++)
            vfsm->index[vfsm->alphabet[i]] = i+1;
    return vfsm;
    }

void freeVFSM(VFSM *vfsm){
    free(vfsm->alphabet);
    free(vfsm);
    return;
    }

/* word2posVFSM : RETURNS str (WORD) AS A LEAF NODE POSITION.
   RETURNS -1 IF GIVEN A NON-VALID STRING.
*/
VFSMTYPE word2posVFSM(VFSM *vfsm, unsigned char *word){
    register int i;
    register VFSMTYPE pos;
    for(i = 0; word[i]; i++)
        if(!vfsm->index[word[i]])
            break;
    if(word[i] || (i != vfsm->depth))
        return -1;
    if(vfsm->ispoweroftwo)
        for(i = 0, pos = 0; i < vfsm->depth; i++){
            pos <<= vfsm->logalphasize;
            pos |= vfsm->index[word[i]]-1;
            }
    else
        for(i = 0, pos = 0; i < vfsm->depth; i++){
            pos *= vfsm->alphasize;
            pos += vfsm->index[word[i]]-1;
            }
    return pos;
    }

/* pos2wordVFSM : RETURNS int (LEAF NODE POSITION) AS A STRING.
   RETURNS NULL IF GIVEN A NON-VALID pos.
*/
BOOLEAN pos2wordVFSM(VFSM *vfsm, VFSMTYPE pos, char *word){
    register int i, j;
    if((pos >= vfsm->lrw) || (pos < 0))
        return FALSE;
    if(vfsm->ispoweroftwo)
        for(i = vfsm->depth-1, j = vfsm->alphasize-1; i >= 0; i--){
            word[i] = vfsm->alphabet[pos&j];
            pos >>= vfsm->logalphasize;
            }
    else
        for(i = vfsm->depth-1; i >= 0; i--){
            word[i] = vfsm->alphabet[pos%vfsm->alphasize];
            pos/=vfsm->alphasize;
            }
    word[vfsm->depth] = '\0';
    return TRUE;
    }

VFSMTYPE changeVFSMstate(VFSM *vfsm, int state, unsigned char newch){
    if(!vfsm->index[newch])
        return 0;     /* NOT A VALID CHARACTER: RESET VFSM */
    if(state >= vfsm->lrs) /* AT BASE: MOVE TO LONGEST SUFFIX STATE */
        state = vfsm->prs+((state-vfsm->lrs) % vfsm->prw);
    /* DESCEND */
    return (state*vfsm->alphasize)+vfsm->index[newch];
    }

VFSMTYPE changeVFSMstatePOW2(VFSM *vfsm, int state,
                             unsigned char newch){
    if(!vfsm->index[newch])
        return 0;     /* NOT A VALID CHARACTER: RESET VFSM */
    if(state >= vfsm->lrs) /* AT BASE: MOVE TO LONGEST SUFFIX STATE */
        state = vfsm->prs+((state-vfsm->lrs) & vfsm->ispoweroftwo);
    /* DESCEND */
    return (state<<vfsm->logalphasize)+vfsm->index[newch];
    }

/* revcompVFSMpos : RETURN A POSITION CORRESPONDING TO THE
                    REVERSE COMPLEMENT OF THE pos.
*/
VFSMTYPE revcompVFSMpos(VFSM *vfsm, VFSMTYPE pos){
    /* REVERSE ALL BITS UNLESS THEY ARE ADJACENT */
    pos = ((pos & 0xcccccccc) >>  2) | ((pos & 0x33333333) <<  2); 
    pos = ((pos & 0xf0f0f0f0) >>  4) | ((pos & 0x0f0f0f0f) <<  4); 
    pos = ((pos & 0xff00ff00) >>  8) | ((pos & 0x00ff00ff) <<  8); 
    pos = (pos >> 16) | (pos << 16);
    return (~pos) >> (32-(vfsm->depth<<1)); /* RETURN COMPLEMENT */
    }
/* ASSUMES VFSMTYPE IS A 32-BIT INTEGER
      -NEED TO ADD A 64-BIT VERSION.
   ASSUMES ALPHABET SIZE IS 4. (ie. 'ACGT')
*/

/* START OF FASTA PARSING VFSM CODE */

static void blankVFSMinit(void *info, int numstates){
    return;
    }
static void blankVFSMstart(void *info, int seqid, int pos){
    return;
    }
static void blankVFSMfinish(void *info, int seqid){
    return;
    }
static void blankVFSMobserve(void *info, int seqid, int seqpos, 
                             int stateid){
    return;
    }
static void blankVFSMvalidch(void *info, int seqid, int seqpos, 
                             char symbol){
    return;
    }

static void parsefastaVFSMprep(PARSEFASTAVFSMFUNCS *pf){
    if(!pf->init)
        pf->init = blankVFSMinit;
    if(!pf->start)
        pf->start = blankVFSMstart;
    if(!pf->finish)
        pf->finish = blankVFSMfinish;
    if(!pf->observe)
        pf->observe = blankVFSMobserve;
    if(!pf->validch)
        pf->validch= blankVFSMvalidch;
    return;
    }

int parsefastaVFSM(FILE *fp, VFSM *vfsm, void *info,
                            PARSEFASTAVFSMFUNCS *pf){
    register char *ueeof = "Unexpected EOF in fasta file";
    register VFSMTYPE state = 0;
    register int ch, prev = '\n', pos = 0, seqid = 1, seqpos = 0;
    parsefastaVFSMprep(pf);
    pf->init(info, vfsm->lrw); /* INIT */
    while(((ch = getc(fp))!='>')&&(prev != '\n')){ /* SKIP TO START */
        if(ch == EOF) /* SKIPPING TO START */
            errmsg(ERROR_FATAL, ueeof);
        prev = ch;
        }
    while((ch = getc(fp)) != '\n') /* SKIP DEFINITON */
        if(ch == EOF)
            errmsg(ERROR_FATAL, ueeof);
    pos = ftell(fp); /* SET RAFASTA POSITION */
    pf->start(info, seqid, pos); /* START */
    do {
        ch = getc(fp);
        switch(ch){
            case EOF:
                pf->finish(info, seqid); /* FINISH */
                if(seqid > 1023) /* END CURRENT LINE OF DOTS */
                    fputc('\n', stderr);
                errmsg(ERROR_INFO,  "Parsed [%d] sequences", seqid);
                return seqid; /* COMPLETE */
            case '>':
                pf->finish(info, seqid); /* FINISH */
                if(!(seqid++ & 1023)) /* SHOW DOT EVERY 1024 SEQS */
                    fputc('.', stderr);
                state = 0;
                while((ch = getc(fp)) != '\n') /* SKIP DEFINITON */
                    if(ch == EOF)
                        errmsg(ERROR_FATAL, ueeof);
                pos = ftell(fp); /* RAFASTA POSITION */
                pf->start(info, seqid, pos); /* START */
                seqpos = 0;
                break;
            case ' ': case '\t': case '\n':
                break; /* SKIP WHITE SPACE */
            default: /* OTHER CHAR */
                if(vfsm->index[ch]){
                    pf->validch(info, seqid, seqpos, (char)ch);
                    state = MchangeVFSMstatePOW2(vfsm,
                                           state, toupper(ch));
                    if(MisleafVFSMstate(vfsm, state)) /* OBSERVE */
                        pf->observe(info, seqid, seqpos,
                            Mstate2posVFSM(vfsm, state));
                } else {
                    state = 0;
                    }
                seqpos++;
                break;
            }
    } while(TRUE);
    /* WILL ONLY RETURN ONCE EOF IS REACHED */
    }

#define REVCOMPSEQCHUNKSIZE 100

static void revcompVFSMvisit(VFSM *vfsm, PARSEFASTAVFSMFUNCS *pf,
                             unsigned char *seq, int len,
                             void *info, int seqid, int pos){
    register VFSMTYPE state = 0;
    register int i;
    pf->start(info, seqid, pos);
    seq[len] = '\0';
    for(i = len-1; i >= 0; i--){
        if(seq[i] == '-'){
            state = 0;
        } else {
            state = MchangeVFSMstatePOW2(vfsm, state, seq[i]);
            if(MisleafVFSMstate(vfsm, state)){ /* OBSERVE */
                pf->observe(info, seqid, len-i+1,
                            Mstate2posVFSM(vfsm, state));
                }
            }
        }
    pf->finish(info, seqid);
    return;
    }

int parsefastaVFSMrevcomp(FILE *fp, VFSM *vfsm,
                          void *info,   PARSEFASTAVFSMFUNCS *pf,
                          void *inforc, PARSEFASTAVFSMFUNCS *pfrc){
    register char *ueeof = "Unexpected EOF in fasta file";
    register VFSMTYPE state = 0;
    register int ch, prev = '\n', pos = 0, seqid = 1, seqpos = 0;
    register int rcseqlen = 0, rcseqalloc = REVCOMPSEQCHUNKSIZE;
    register char *rcseq = malloc(sizeof(char)*rcseqalloc);
    parsefastaVFSMprep(pf);
    parsefastaVFSMprep(pfrc);
    pf->init(info, vfsm->lrw); /* INIT */
    pfrc->init(inforc, vfsm->lrw); /* INIT */
    while(((ch = getc(fp))!='>')&&(prev != '\n')){ /* SKIP TO START */
        if(ch == EOF) /* SKIPPING TO START */
            errmsg(ERROR_FATAL, ueeof);
        prev = ch;
        }
    while((ch = getc(fp)) != '\n') /* SKIP DEFINITON */
        if(ch == EOF)
            errmsg(ERROR_FATAL, ueeof);
    pos = ftell(fp); /* SET RAFASTA POSITION */
    pf->start(info, seqid, pos); /* START */
    do {
        ch = getc(fp);
        switch(ch){
            case EOF:
                pf->finish(info, seqid); /* FINISH */
                revcompVFSMvisit(vfsm, pfrc, (unsigned char*)rcseq,
                                 rcseqlen, inforc, seqid, pos);
                rcseqlen = 0;
                if(seqid > 1023) /* END CURRENT LINE OF DOTS */
                    fputc('\n', stderr);
                errmsg(ERROR_INFO,  "Parsed [%d] sequences", seqid);
                free(rcseq);
                return seqid; /* COMPLETE */
            case '>':
                pf->finish(info, seqid); /* FINISH */
                revcompVFSMvisit(vfsm, pfrc, (unsigned char*)rcseq,
                                 rcseqlen, inforc, seqid, pos);
                rcseqlen = 0;
                if(!(seqid++ & 1023)) /* SHOW DOT EVERY 1024 SEQS */
                    fputc('.', stderr);
                state = 0;
                while((ch = getc(fp)) != '\n') /* SKIP DEFINITON */
                    if(ch == EOF)
                        errmsg(ERROR_FATAL, ueeof);
                pos = ftell(fp); /* RAFASTA POSITION */
                pf->start(info, seqid, pos); /* START */
                seqpos = 0;
                break;
            case ' ': case '\t': case '\n':
                break; /* SKIP WHITE SPACE */
            default: /* OTHER CHAR */
                if(vfsm->index[ch]){
                    if(rcseqlen >= rcseqalloc)
                        rcseq = realloc(rcseq, sizeof(char)
                              *(rcseqalloc+=REVCOMPSEQCHUNKSIZE));    
                    rcseq[rcseqlen++] = SEQUTILcomplement(ch);
                    state = MchangeVFSMstatePOW2(vfsm, state, toupper(ch));
                    if(MisleafVFSMstate(vfsm, state)){  /* OBSERVE */
                        pf->observe(info, seqid, seqpos,
                            Mstate2posVFSM(vfsm, state));
                          }
                    } else {
                        if(rcseqlen >= rcseqalloc)
                            rcseq = realloc(rcseq, sizeof(char)
                                  *(rcseqalloc+=REVCOMPSEQCHUNKSIZE));
                        rcseq[rcseqlen++] = '-';
                        state = 0;
                        }
                    seqpos++;
                    break;
            }
    } while(TRUE);
    /* WILL ONLY RETURN ONCE EOF IS REACHED */
    }

/*   END OF FASTA PARSING VFSM CODE */

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */
static void infoVFSM(VFSM *vfsm){
    printf("VFSM info:\n"

           " alphabet      = [%s]\n"
           " alphabet size = %d\n"
           " depth         = %d\n"
           " prs           = %" VFSMTYPE_PRINT "\n"
           " prw           = %" VFSMTYPE_PRINT "\n"
           " lrs           = %" VFSMTYPE_PRINT "\n"
           " lrw           = %" VFSMTYPE_PRINT "\n"
           " total         = %" VFSMTYPE_PRINT "\n"
           " ispoweroftwo  = %d\n",

           vfsm->alphabet,
           vfsm->alphasize,
           vfsm->depth,
           vfsm->prs,
           vfsm->prw,
           vfsm->lrs,
           vfsm->lrw,
           vfsm->total,
           vfsm->ispoweroftwo
           );
    return;
    }

int main(){
    register int i, j;
    register VFSM *vfsm = newVFSM("ACGT", 12, TRUE);
    register unsigned char *seq = (unsigned char*)
                 "CGATCGATCTGATCGTAGNTAGCTCGATCGATGNAGCTAGC";
    register VFSMTYPE state = 0;
    char word[12];
    infoVFSM(vfsm);
/*
    for(i = 0; i < vfsm->lrw; i++){
        pos2wordVFSM(vfsm, i, word);
        printf("[%d] is [%s]\n", i, word);
        if(word2posVFSM(vfsm, word) != i)
            errmsg(ERROR_FATAL, "Bad word2posVFSM conversion");
        }
*/
    for(i = j = 0; seq[i]; i++){
        /* state = changeVFSMstate(vfsm, state, seq[i]); */
        if(!vfsm->index[seq[i]]){
            state = 0;
            continue;
            }
        state = MchangeVFSMstatePOW2(vfsm, state, seq[i]);
        if(pos2wordVFSM(vfsm, state, word))
            printf("Word number [%2d]=[%s]\n", ++j, word);
        if(pos2wordVFSM(vfsm, revcompVFSMpos(vfsm, state), word))
            printf("Revcomp     [%2d]=[%s]\n", j, word);
        }
    freeVFSM(vfsm);
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_VFSM_C */

