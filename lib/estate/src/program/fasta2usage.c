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

/* fasta2usage : Calculate word frequency scores from
                 fasta format sequences. Version 2.1
*/

#include <ctype.h>
#include "../general/arg.h"
#include "../general/error.h"
#include "../general/common.h"
#include "../struct/vfsm.h"

#include "../parse/readfasta.h"
#include "../general/dict.h"

ALLOW_ERROR_MESSAGES;

typedef struct {
   int *count;
   int  numstates;
   int jump;
   } WORDCOUNTINFO;

typedef struct {
    int jump;
    int *count;
    int composition[ALPHABETSIZE];
    } WORDUSAGEINFO;

static void observeVFSMWordUsage(void *info, int seqid,
                                 int seqpos, int stateid){
    register WORDUSAGEINFO *wui = (WORDUSAGEINFO*)info;
    if(!((seqpos+1) % wui->jump))
        wui->count[stateid]++;
    return;
    }

static void validcharVFSMWordUsage(void *info, int seqid,
                                  int seqpos, char symbol){
    register WORDUSAGEINFO *wui = (WORDUSAGEINFO*)info;
    wui->composition[toupper(symbol)]++;
    return;
    }

static void printComposition(int *comp){
    register int i;
    printf("#Composition");
    for(i = 0; i < ALPHABETSIZE; i++)
        if(comp[i])
            printf(" %c=%d", i, comp[i]);
    printf("\n");
    return;
    }

static void walkshowingcount(char *name, void *value, int count){
    printf("%s %d\n", name, (int)value);
    return;
    }

static void countWithDictionary(char *path, int wordlen,
                                int jump, char *alphabet){
    register READFASTA *rf = newREADFASTA(path);
    register DICT *d = newDICT();
    register int i, validlen;
    register char tmpch;
    int comp[ALPHABETSIZE], index[ALPHABETSIZE];
    memset(comp, 0, sizeof(int)*ALPHABETSIZE);
    memset(index, 0, sizeof(int)*ALPHABETSIZE);
    for(i = 0; alphabet[i]; i++)
        index[(unsigned char)alphabet[i]] = alphabet[i];
    while(nextREADFASTA(rf)){
        for(i = validlen = 0; i < rf->seqlen; i++)
            if(index[(unsigned char)rf->seq[i]]){
                validlen++;
                comp[index[(unsigned char)rf->seq[i]]]++;
                    if((validlen >= wordlen) && (!((i+1) % jump))){
                        tmpch = rf->seq[i+1];
                        rf->seq[i+1] = '\0';
                        countDICT(d, rf->seq+i-wordlen+1);
                        rf->seq[i+1] = tmpch;
                        }
            } else {
                validlen = 0;
                }
        }
    printComposition(comp);
    sortWalkDICT(d, walkshowingcount);
    freeREADFASTA(rf);
    freeDICT(d);
    return;
    }

static void countWordUsage(char *path, FILE *fp, int wordlen,
                           int jump, char *alphabet){
    register VFSM *vfsm = newVFSM(alphabet, wordlen, FALSE);
    register char *word = malloc(sizeof(char)*wordlen);
    register int i;
    register double maxcheck;
    WORDUSAGEINFO wui;
    PARSEFASTAVFSMFUNCS pf = {NULL, NULL, NULL,
                              observeVFSMWordUsage,
                              validcharVFSMWordUsage};
    memset(wui.composition, 0, sizeof(int)*ALPHABETSIZE);
    wui.jump = jump;
    maxcheck = (double)vfsm->lrw*sizeof(int);
    if(maxcheck < VFSMTYPE_MAX)
        wui.count = calloc(vfsm->lrw, sizeof(int));
    else
        wui.count = NULL;
    if(wui.count){
        parsefastaVFSM(fp, vfsm, &wui, &pf);
        printComposition(wui.composition);
        for(i = 0; i < vfsm->lrw; i++){
            if(wui.count[i]){
                pos2wordVFSM(vfsm, i, word);
                printf("%d %s\n", wui.count[i], word);
                }
            }
        free(wui.count);
    } else {
        errmsg(ERROR_WARNING, "Insufficient memory for VFSM");
        errmsg(ERROR_WARNING, "Switching to use DICT");
        countWithDictionary(path, wordlen, jump, alphabet);
        }
    freeVFSM(vfsm);
    free(word);
    return;
    }

static void countComposition(char *path, int jump, char *alphabet){
    register READFASTA *rf = newREADFASTA(path);
    register int i;
    int comp[ALPHABETSIZE], index[ALPHABETSIZE], single[ALPHABETSIZE];
    memset(comp, 0, sizeof(int)*ALPHABETSIZE);
    memset(index, 0, sizeof(int)*ALPHABETSIZE);
    memset(single, 0, sizeof(int)*ALPHABETSIZE);
    for(i = 0; alphabet[i]; i++)
        index[(unsigned char)alphabet[i]] = alphabet[i];
    while(nextREADFASTA(rf)){
        for(i = 0; i < rf->seqlen; i++)
            if(index[(unsigned char)rf->seq[i]]){
                comp[index[(unsigned char)rf->seq[i]]]++;
                    if(!((i+1) % jump))
                        single[index[(unsigned char)rf->seq[i]]]++;
                }
        }
    printComposition(comp);
    for(i = 0; i < ALPHABETSIZE; i++)
        if(single[i])
            printf("%c %d\n", i, single[i]);
    freeREADFASTA(rf);
    return;
    }

static void countPairUsage(char *path, int jump, char *alphabet){
    register READFASTA *rf = newREADFASTA(path);
    register int i, j;
    int comp[ALPHABETSIZE], index[ALPHABETSIZE],
        pair[ALPHABETSIZE][ALPHABETSIZE];
    register unsigned char prev;
    memset(comp, 0, sizeof(int)*ALPHABETSIZE);
    memset(index, 0, sizeof(int)*ALPHABETSIZE);
    for(i = 0; i < ALPHABETSIZE; i++)
        memset(pair[i], 0, sizeof(int)*ALPHABETSIZE);
    for(i = 0; alphabet[i]; i++)
        index[(unsigned char)alphabet[i]] = alphabet[i];
    while(nextREADFASTA(rf)){
        if(index[(unsigned char)rf->seq[0]]){
            prev = rf->seq[0];
            comp[index[prev]]++;
        } else {
            prev = 0;
            }
        for(i = 1; i < rf->seqlen; i++){
            if(index[(unsigned char)rf->seq[i]]){
                comp[index[(unsigned char)rf->seq[i]]]++;
                if(prev)
                    if(!((i+1) % jump))
                        pair[index[prev]]
                            [index[(unsigned char)rf->seq[i]]]++;
                prev = rf->seq[i];
            } else {
                prev = 0;
                }
            }
        }
    printComposition(comp);
    for(i = 0; i < ALPHABETSIZE; i++)
        for(j = 0; j < ALPHABETSIZE; j++)
        if(pair[i][j])
            printf("%c%c %d\n", i, j, pair[i][j]);
    freeREADFASTA(rf);
    return;
    }

#define ARGUMENT_COUNT 5

int estateMAIN(){
    static char *alphabet = "ACGT";
    static char *dbpath = NULL;
    static FILE *dbfp;
    static int   wordlen = 3;
    static int   jump    = 3;
    static BOOLEAN forcedict = FALSE;
    ARGUMENT argdata[ARGUMENT_COUNT] = {
    { ARG_FILE,  'd', "database", "fasta format sequences",
      &dbfp, &dbpath, "FASTA2USAGE_DATABASE", "r",  TRUE},
    { ARG_INT,   'w', "wordlen",  "word length",
     &wordlen, NULL, "FASTA2USAGE_WORDLEN", NULL, FALSE},
    { ARG_INT,   'j', "jump", "jump between words",
     &jump, NULL, "FASTA2USAGE_JUMP",   NULL, FALSE},
    { ARG_BOOLEAN,  'f', "forcedict", "force use of dictionary",
     &forcedict, NULL, "FASTA2USAGE_FORCEDICT",   NULL, FALSE},
    { ARG_STRING,  'a', "alphabet", "alphabet to score",
     &alphabet, NULL, "FASTA2USAGE_ALPHABET",  NULL, FALSE}
    };
    dbfp = stdin;
    processARGUMENT(argdata, ARGUMENT_COUNT, global_emn,
       "calculate word usage counts from fasta format sequences");
    if(!strcmp(alphabet, "protein"))
        alphabet = "ARNDCQEGHILKMFPSTWYV";
    errmsg(ERROR_INFO, "Parsing [%s]", dbpath);
    if(wordlen < 1)
        errmsg(ERROR_FATAL, "Wordlen must be more greater than 0");
    if(forcedict)
        countWithDictionary(dbpath, wordlen, jump, alphabet);
    else
        switch(wordlen){
            case 1:
                countComposition(dbpath, jump, alphabet);
                break;
            case 2:
                countPairUsage(dbpath, jump, alphabet);
                break;
            default:
                countWordUsage(dbpath, dbfp, wordlen, jump, alphabet);
                break;
            }
    return 0;
    }

