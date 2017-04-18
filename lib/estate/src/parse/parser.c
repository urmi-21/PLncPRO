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

#ifndef INCLUDED_PARSER_C
#define INCLUDED_PARSER_C

#include <stdlib.h>
#include <string.h> /* FOR strdup() */
#include "parser.h"
#include "../general/common.h"

PARSER *newPARSER(FILE *fp){
    PARSER *p = NEW(PARSER);
    p->fp = fp;
    p->linealloc = PARSER_CHUNK_SIZE;
    p->line = malloc(sizeof(char)*p->linealloc);
    p->wordalloc = PARSER_CHUNK_SIZE;
    p->word = malloc(sizeof(char**)*p->wordalloc);
    return p;
    }

void freePARSER(PARSER *p){
    free(p->word);
    free(p->line);
    free(p);
    return;
    }

/* linePARSER : read next line of input into line.
                returns the line length.
*/
int linePARSER(PARSER *p){
    register int ch, pos = 0;
    for(;;){
        switch(ch = getc(p->fp)){
            case '\n':
                p->line[pos] = '\0';
                return pos;
            case EOF:
                p->line[pos] = '\0';
                return pos?pos:EOF;
            default:
                if(p->linealloc <= pos){
                    p->linealloc+=PARSER_CHUNK_SIZE;
                    p->line = realloc(p->line, 
                                      sizeof(char)*p->linealloc);
                    }
                p->line[pos++] = ch;
                break;
            }
        }
    /* return; */
    }

/* wordPARSERoffset : SPLIT LAST LINE OF length INTO WORDS
                      BEGINNING AT start.
*/
int wordPARSERoffset(PARSER *p, int length, int start){
    register unsigned char *ptr;
    register int pos = 0;
    switch(length){
        case 0:
            return 0;
        case EOF:
            return EOF;
        default:
            break;
        }
    if((p->line[start] != ' ') && (p->line[start] != '\t'))
        p->word[pos++] = p->line+start;
    for(ptr = p->line+start; *ptr; ptr++){
        if((*ptr == ' ') || (*ptr == '\t')){
            *ptr = '\0';
            do { ptr++;
            } while((*ptr == ' ') || (*ptr == '\t'));
            if(!*ptr)
                break;
            if(p->wordalloc <= pos){
                p->wordalloc+=PARSER_CHUNK_SIZE;
                p->word = realloc(p->word,
                                  sizeof(char**)*p->wordalloc);
                }
            p->word[pos++] = ptr;
            }
        }
    p->word[pos] = NULL;
    return pos;
    }

/* wordPARSER : read next line and splits into words.
                returns the number of words read.
*/
int wordPARSER(PARSER *p){
    return wordPARSERoffset(p, linePARSER(p), 0);
    }

/* linePARSERfile : READ IN A WHOLE FILE, RETURNING ARRAY OF LINES,
                    RECORDING NUMBER IN total. IGNORES BLANK LINES.
*/
char **linePARSERfile(FILE *fp, int *total){
    register char **ptp;
    register PARSER *p = newPARSER(fp);
    int linesread  = 0, linesalloc = PARSER_CHUNK_SIZE, len;
    ptp = malloc(sizeof(char*)*linesalloc);
    while((len = wordPARSER(p)) != EOF){
        if(!*p->line) /* IGNORE BLANK LINES */
            continue;
        if(linesread >= linesalloc-2){
            linesalloc += PARSER_CHUNK_SIZE;
            ptp = realloc(ptp, sizeof(char*)*linesalloc); 
            }
        ptp[linesread++] = strdup((char*)p->line);
        }
    ptp[linesread] = (char*)0;
    freePARSER(p); 
    if(total)
        *total = linesread;
    return ptp;
    }

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

int main(){
    register PARSER *p = newPARSER(stdin);
    register int i, len;
    /* while((len = linePARSER(p)) != EOF) */
        /* printf("LINE [%d][%s]\n", len, p->line); */
    while((len = wordPARSER(p)) != EOF){
        if((len == 1) && (!strcmp((char*)p->word[0], "quit")))
            break;
        printf("WORD [%d] - ", len);
        for(i = 0; i < len; i++)
            printf("[%s] ", p->word[i]);
        printf("\n");
        }
    freePARSER(p);
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */

#endif /* INCLUDED_PARSER_C */

