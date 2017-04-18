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

/* wordsearch : find dictionary words in sequence database.
   Guy St.C. Slater..  Version 2.0  November 1997.
*/

/* This is just here as a test of the finite state machine library.
*/

#include "../general/arg.h"
#include "../general/error.h"
#include "../general/common.h"
#include "../general/stringlist.h"

ALLOW_ERROR_MESSAGES;

#include "../struct/fsm.h"

#include "../parse/readfasta.h"
#include "../parse/ioutil.h"
#include "../parse/parser.h"

#include <ctype.h>

static BOOLEAN checkword(unsigned char *word, int len, BOOLEAN *aav){
    register int i;
    for(i = 0; i < len; i++){
        word[i] = toupper(word[i]);
        if(!aav[word[i]])
            return FALSE;
        }
    return TRUE;
    }

static unsigned char *local_aas = (unsigned char*)
                                  "ARNDCQEGHILKMFPSTWYV";

static STRINGLIST *getValidWords(FILE *fp, int minlen){
    register int i, len;
    register STRINGLIST *sl = newSTRINGLIST();
    register PARSER *p = newPARSER(fp);    
    BOOLEAN aav[ALPHABETSIZE];
    for(i = 0; i < ALPHABETSIZE; i++)
        aav[i] = FALSE;
    for(i = strlen((char*)local_aas)-1; i >= 0; i--)
        aav[local_aas[i]] = TRUE;
    while((len = linePARSER(p)) != EOF){
        if(len < minlen)
            continue;
        if(!checkword(p->line, len, aav))
            continue;
        addSTRINGLIST(sl, (char*)p->line, len);
        }
    if(sl->c)
        errmsg(ERROR_INFO, "Selected [%d] valid words for FSM", 
                            sl->c);
    else
        errmsg(ERROR_FATAL, "Found no suitable words in dictionary");
    freePARSER(p);
    return sl;
    }

static void searchSequence(FSM *f, unsigned char *seq, char *def){
    register FSMNODE *n = f->root;
    register short ch;
    do {
        if(n[ch = f->index[*seq]].data.n){
           printf(">found [%s] length = %d in [%.40s]\n", 
                     (char*)n[ch].data.n,
                     strlen(n[ch].data.n), def);
           fflush(stdout);
           }
        n = n[ch].next;
    } while(*++seq);
    return;
    }
 
static void wordsearchFSMdatabase(STRINGLIST *sl, char *path){
    register FSM *fsm = makeFSMptr(sl->v, sl->c, (void**)sl->v);
    register READFASTA *rf = newREADFASTA(path);
    infoFSM(fsm);
    errmsg(ERROR_INFO, "Compiled FSM");
    while(nextREADFASTA(rf))
        searchSequence(fsm, (unsigned char*)rf->seq, rf->def);
    freeREADFASTA(rf);
    freeFSM(fsm);
    return;
    }

#define ARGUMENT_COUNT 2

int estateMAIN(){
    static FILE *sdb  = NULL;
    static FILE *wdb;
    static int minlen = 7;
    static char *sdbname;

    ARGUMENT argdata[ARGUMENT_COUNT] = {
    {ARG_FILE, 'd', "database", "fasta format database",
    &sdb, &sdbname, "WORDSEARCH_DATABASE_PATH", "r", TRUE},
    {ARG_INT,  'l', "minlen", "minimum length",
    &minlen, NULL, "WORDSEARCH_MINIMUM_LENGTH",  "[0,10]", FALSE}
    };

    register STRINGLIST *sl;
    wdb = stdin;
    processARGUMENT(argdata, ARGUMENT_COUNT, global_emn,
              "dictionary word search in sequence databases");
    sl = getValidWords(wdb, minlen);
    wordsearchFSMdatabase(sl, sdbname);
    freeSTRINGLIST(sl);
    return 0;
    }

