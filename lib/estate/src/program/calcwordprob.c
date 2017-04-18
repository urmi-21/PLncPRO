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

/* calcwordprob : expand DNA word counts to guess for Ns 
                  with log odds probababilities.
*/

#include "../general/arg.h"
#include "../general/error.h"
#include "../general/common.h"
#include "../parse/parser.h"
#include "../struct/vfsm.h"
#include <math.h>

ALLOW_ERROR_MESSAGES;

static void expandforNrecur(int pos, int count, int uid,
                         VFSM *vout, int *voutcount, int *voutmail){
    register int i, currpos;
    register char c;
    register char *word = malloc(sizeof(char)*vout->depth);
    pos2wordVFSM(vout, pos, word);
    voutmail[pos] = uid;
    voutcount[pos] += count;
    for(i = 0; i < vout->depth; i++){
        if(word[i] != 'N'){
            c = word[i];
            word[i] = 'N';
            currpos = word2posVFSM(vout, (unsigned char*)word);
            if(voutmail[currpos] != uid){
                expandforNrecur(currpos, count, uid,
                                vout, voutcount, voutmail);
                }
            word[i] = c;
            }
        }
    free(word);
    return;
    }

static void expandforN(VFSM *vin, VFSM *vout, 
                   int *vincount, int *voutcount, int *voutmail){
    register int i;
    register char *word = malloc(sizeof(char)*vin->depth);
    for(i = 0; i < vin->lrw; i++){ /* FOR EACH [ACGT] word */
        pos2wordVFSM(vin, i, word);
        expandforNrecur(word2posVFSM(vout, (unsigned char*)word),
                        vincount[i], i+1, vout, voutcount, voutmail);
        }
    free(word);
    return;
    }

/* CalculateChanceWordProbability : calculate the chance
                                    of word occuring by chance.
*/
static double CalculateChanceWordProbability(unsigned char *word,
                                  int len, double *compprob){
    register unsigned char *a, *z;
    register double cp = compprob[*word];
    for(a = word+1, z = word+len; a < z; a++)
        cp *= compprob[*a];
    return cp;
    }

static void printProbs(VFSM *vout, int *voutcount, int *composition){
    register int i;
    register char *word = malloc(sizeof(char)*vout->depth);
    register int counttotal = voutcount[vout->lrw-1];
    register double rawprob, chanceprob, oddsprob, logoddsprob;
    double compprob[ALPHABETSIZE] = {0.0};
    register double logbase = log(2.0);
    compprob['N'] += compprob['A'] = composition['A'];
    compprob['N'] += compprob['C'] = composition['C'];
    compprob['N'] += compprob['G'] = composition['G'];
    compprob['N'] += compprob['T'] = composition['T'];
    compprob['A'] /= compprob['N'];
    compprob['C'] /= compprob['N'];
    compprob['G'] /= compprob['N'];
    compprob['T'] /= compprob['N'];
    compprob['N'] /= compprob['N'];
    errmsg(ERROR_INFO, "Composition percentage :"
           " A=%2.2f C=%2.2f G=%2.2f T=%2.2f N=%2.2f",
                     compprob['A']*100.0,
                     compprob['C']*100.0,
                     compprob['G']*100.0,
                     compprob['T']*100.0,
                     compprob['N']*100.0);
    for(i = 0; i < vout->lrw; i++){
        pos2wordVFSM(vout, i, word); 
        if(voutcount[i])
            rawprob = ((double)voutcount[i])/((double)counttotal); 
        else /* AVOID -Inf PROBABILITIES WHERE COUNT == 0 */
            rawprob = 0.48/((double)counttotal); 
        chanceprob = CalculateChanceWordProbability((unsigned char*)
                                   word, vout->depth, compprob);
        oddsprob = rawprob/chanceprob;
        logoddsprob = log(oddsprob)/logbase;
        /* printf("%.9f %9d %s cp=%.9f op=%.9f lop=%.9f\n", */
              /* rawprob, voutcount[i], word, */
              /* chanceprob, oddsprob, logoddsprob); */
        printf("% .9f %s\n", logoddsprob, word);
        }
    free(word);
    return;
    }

static void expandWordCounts(FILE *datafp){
    register VFSM *vin, *vout;
    register PARSER *p = newPARSER(datafp);
    register int i, count, wordlen, *vincount, *voutcount, *voutmail;
    int composition[ALPHABETSIZE] = {0};
    do { 
        count = wordPARSER(p);
        if((count) && (!strcmp((char*)p->word[0], "#Composition"))){
            for(i = 1; i < count; i++) /* READ COMPOSITION INFO */
                composition[p->word[i][0]] = atoi((char*)p->word[i]+2); 
            continue;
            }
        if(count == EOF)
            errmsg(ERROR_FATAL, "No words in input");
    } while(count != 2);
    wordlen = strlen((char*)p->word[1]);
    vin  = newVFSM("ACGT",  wordlen, FALSE);
    vout = newVFSM("ACGTN", wordlen, FALSE);

    vincount = malloc(vin->lrw*sizeof(int));
    if(!vincount)
        errmsg(ERROR_FATAL, "Insufficient memory for vfsm (1)");
    voutcount = calloc(vout->lrw, sizeof(int));
    if(!voutcount)
        errmsg(ERROR_FATAL, "Insufficient memory for vfsm (2)");
    voutmail = calloc(vout->lrw, sizeof(int));
    if(!voutmail)
        errmsg(ERROR_FATAL, "Insufficient memory for vfsm (3)");

    do {
       if(count == 2)
          vincount[word2posVFSM(vin, p->word[1])]
                                 = atoi((char*)p->word[0]);
    } while((count = wordPARSER(p)) != EOF);

    expandforN(vin, vout, vincount, voutcount, voutmail);

    free(voutmail);

    freeVFSM(vin);
    free(vincount);

    printProbs(vout, voutcount, composition);

    freeVFSM(vout);
    free(voutcount);

    freePARSER(p);
    return;
    }

#define ARGUMENT_COUNT 1

int estateMAIN(){
    static FILE *datafp;
    static char *datapath = NULL; 
    ARGUMENT argdata[ARGUMENT_COUNT] = {
    { ARG_FILE,  'w', "wordcount", "word count data file",
     &datafp, &datapath, "CALCWORDPROB_WORDCOUNT", "r",  TRUE},  
    };
    datafp = stdin;
    processARGUMENT(argdata, ARGUMENT_COUNT, global_emn, 
      "generate log-odds probabilities from wordcounts (with Ns)");
    expandWordCounts(datafp);
    return 0;
    }

