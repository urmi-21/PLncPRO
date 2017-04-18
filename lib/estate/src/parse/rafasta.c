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

#ifndef INCLUDED_RAFASTA_C
#define INCLUDED_RAFASTA_C

#include <stdio.h>
#include <ctype.h> /* FOR isspace() */
#include <string.h> /* FOR memmove() */
#include "../general/common.h"
#include "../general/error.h"
#include "../parse/rafasta.h"
#include "../parse/ioutil.h"

ALLOW_ERROR_MESSAGES;

/* getRAFASTAindices : RETURNS AN ARRAY WITH total MEMBERS 
                       OF START POINTS FOR ALL THE SEQUENCES IN fp.
*/
long *getRAFASTAindices(FILE *fp, int *total){
    register int ch, alloc = RAFASTA_INDEX_CHUNKSIZE, count = 0;
    register long *index = malloc(sizeof(int)*alloc);
    register BOOLEAN linestart = TRUE, defline = TRUE;
    do {
        switch(ch = getc(fp)){
            case EOF:
                index = realloc(index, sizeof(int)*count);
                if(total)
                    *total = count;
                return index;
            case '\n':
                linestart = TRUE;
                if(defline == TRUE){
                    if(count >= alloc){
                        alloc <<= 1;
                        index = realloc(index, sizeof(int)*alloc);
                        if(!index)
                          errmsg(ERROR_FATAL, "Insufficient memory"); 
                        }
                    index[count] = ftell(fp);
                    count++;
                    defline   = FALSE;
                    }
                break;
            case '>':
                if(linestart == TRUE)
                    defline = TRUE;
                /* FALL THROUGH */
            default:
                linestart = FALSE;
                break;
            }
    } while(TRUE);
    /* MUST REACH EOF TO EXIT LOOP */
    }

/* printRAFASTAseq : READ SEQUENCE FROM pos IN fp TO out.
*/
int printRAFASTAseq(FILE *fp, long pos, FILE *out){
    register int ch, seqlen = 0;
    register BOOLEAN linestart = FALSE;
    fseek(fp, pos, SEEK_SET);
    do {
        switch(ch = getc(fp)){
            case '\n':
                linestart = TRUE;
                putc('\n', out);
                break;
            case '>':
                if(linestart != TRUE)
                    break; 
                /* FALL THROUGH */ 
            case EOF:
                return seqlen;
            default:
                if(isalpha(ch)){
                    putc(ch, out);
                    seqlen++;
                    }
                linestart = FALSE;
                break;
            } 
    } while(TRUE);
    /* MUST REACH EOF OR END OF SEQUENCE TO EXIT LOOP */
    }

/* getRAFASTAseq : RETURN SEQUENCE FROM pos IN fp.
                   length WILL CONTAIN THE SEQUENCE LENGTH.
*/
char *getRAFASTAseq(FILE *fp, long pos, int *length){
    register int ch;
    register int alloc = RAFASTA_SEQ_CHUNKSIZE;
    register BOOLEAN linestart = FALSE;
    register char *seq = malloc(sizeof(char)*alloc);
    fseek(fp, pos, SEEK_SET);
    *length = 0;
    do {
        switch(ch = getc(fp)){
            case '\n':
                linestart = TRUE;
                break;
            case '>':
                if(linestart != TRUE)
                    break; 
                /* FALL THROUGH */ 
            case EOF:
                seq[*length] = '\0'; 
                return seq;
            default:
                if(isalpha(ch)){
                    if((*length) >= alloc){
                        alloc <<= 1;
                        seq = realloc(seq, sizeof(char)*alloc);
                        }  
                    seq[(*length)++] = toupper(ch);
                    }
                linestart = FALSE;
                break;
            } 
    } while(TRUE);
    /* MUST REACH EOF OR END OF SEQUENCE TO EXIT LOOP */
    }

/* locateRAFASTAdefn :  RETURNS POSITION OF START OF DEFINITION
                      FOR SEQUENCE STARTING AT pos IN FILE fp. 
*/
long locateRAFASTAdefn(FILE *fp, long pos){
    register int i, ch = 0;
    register long start, end = pos;
    register long backtrack = RAFASTA_BACKTRACK_GUESS;
    register long location  = RAFASTA_NO_DEFINITION;
    register BOOLEAN linestart;

    do { 
         start = end;
         end = start-backtrack;
         backtrack <<= 1; /* DOUBLE BACKTRACK SIZE */
         if(end <= 0){
             end = 0;
             linestart = TRUE;
         } else {
             linestart = FALSE;
             }
         fseek(fp, end, SEEK_SET);
         for(i = end; i < start; i++){
             switch(ch = getc(fp)){
                 case EOF:
                     errmsg(ERROR_FATAL, "No sequence definition"); 
                 case '\n':
                     linestart = TRUE;
                     break;
                 case '>':
                     if(linestart == TRUE)
                         location = ftell(fp);
                     /* FALL THROUGH */ 
                 default:
                     linestart = FALSE;
                     break;
                 }
             }        
         if(location != RAFASTA_NO_DEFINITION)
             return location;
    } while(end++); /* INCASE JUMP BACK TO OVERLAP POINT */
    return RAFASTA_NO_DEFINITION;
    }

/* printRAFASTAdef : PRINT FASTA DEFINITION FOR pos IN fp TO out. 
*/
void printRAFASTAdef(FILE *fp, long pos, FILE *out){
    register int i;
    register long location = locateRAFASTAdefn(fp, pos);  
    if(location == RAFASTA_NO_DEFINITION){
        errmsg(ERROR_WARNING, "Could not find definition"); 
        return;
        }
    fseek(fp, location, SEEK_SET);
    for(i = pos-location-1; i; i--) 
        putc(getc(fp), out); /* OK FOLLOWING locateRAFASTAdefn CALL */
    putc('\n', out);
    return; 
    }

/* getRAFASTAdef : RETURN THE DEFINITION PRECEEDING pos IN fp.
*/
char *getRAFASTAdef(FILE *fp, long pos){
    register int len;
    register long location = locateRAFASTAdefn(fp, pos);  
    register char *def;
    if(location == RAFASTA_NO_DEFINITION){
        errmsg(ERROR_WARNING, "Could not find definition"); 
        return NULL;
        }
    len = pos-location-1; 
    fseek(fp, location, SEEK_SET);
    def = malloc(sizeof(char)*(len+1));
    fread(def, sizeof(char), len, fp);
    def[len] = '\0';
    return def; 
    }

/* writeRAFASTAindices : READS FASTA FILE FROM infp, 
                         WRITING INDEX FILE TO outfp.
                         RETURNS NUMBER OF SEQUENCES READ.
*/
int writeRAFASTAindices(FILE *infp, FILE *outfp){
    register int ch, total = 0;
    register BOOLEAN linestart = TRUE, defline = TRUE;
    long pos;
    do {
        switch(ch = getc(infp)){
            case EOF:
                return total;
            case '\n':
                linestart = TRUE;
                if(defline == TRUE){
                   pos = ftell(infp);
                   fwrite(&pos, sizeof(long), 1, outfp);
                   total++;
                   defline = FALSE;
                   }
            break;
            case '>':
                if(linestart == TRUE)
                    defline = TRUE;
                /* FALLTHROUGH */
            default:
                linestart = FALSE;
            }
    } while(TRUE);
    /* MUST REACH EOF TO RETURN FROM FUNCTION */
    }

long *readRAFASTAindices(char *path, int *total){
    long count;
    register long *l = (long*)readFILE(path, &count);
    if(!l)
        errmsg(ERROR_FATAL, "Missing index file \"%s\"", path);
    *total = count/sizeof(long);
    return l;
    }

RAFASTAALL *readRAFASTAALL(char *path, int *total){
    long length;
    register  int seqcount = 0;
    register char *filedata = readFILE(path, &length);
    register char *p = filedata, *dst = p;
    register RAFASTAALL *rfa = InitAlloc(RAFASTAALL);
    if(!filedata)
        errmsg(ERROR_FATAL, "Error reading fasta file \"%s\"", path);
    while(*p != '>') /* SKIP TO START OF FILE */
        if(!*p)
            errmsg(ERROR_FATAL, "No sequences in [%s]", path);
    p++; /* PASS '>' */
    do {
        StepRealloc(rfa, seqcount, RAFASTAALL);
        rfa[seqcount].def = p++;
        while(*p != '\n') /* SKIP TO END OF DEFINITION */
            if(!*p++)
                errmsg(ERROR_FATAL, "Unexpected EOF in [%s]", path);
        *p++ = '\0'; /* END DEFINITION */
        while(isspace(*p)) /* SKIP TO START OF SEQUENCE */
            if(!*p++)
                errmsg(ERROR_FATAL, "Unexpected EOF in [%s]", path);
        rfa[seqcount].seq = dst = p;
        rfa[seqcount].pos = p-filedata;
        do {
            switch(*p){
                case '>':
                    *dst = '\0';
                    rfa[seqcount].len = dst-rfa[seqcount].seq;
                    break;
                case '\n': /*FALLTHROUGH*/
                case '\t': /*FALLTHROUGH*/
                case ' ': 
                    break;
                case '\0':
                    *dst = '\0';
                    rfa[seqcount].len = dst-rfa[seqcount++].seq;
                    StepRealloc(rfa, seqcount, RAFASTAALL);
                    memset(&rfa[seqcount], 0, sizeof(RAFASTAALL));
                    if(total)
                        *total = seqcount;
                    memmove(filedata, rfa->def, /* MAKE FREEABLE */
                           sizeof(char)*(strlen(rfa->def)+1));
                    rfa->def = filedata;
                    return rfa;
                default:
                    *dst++ = toupper(*p);
                    break;
                }
        } while(*p++ != '>');
        seqcount++;
    } while(TRUE);
    /* NEVER REACHED */
    }

void freeRAFASTAALL(RAFASTAALL *rfa){
    free(rfa->def);
    free(rfa);
    return;
    }


#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

int main(int argc, char **argv){
    FILE *fp = (argc>1)?fopen(argv[1], "r"):stdin;
    int i, len, seqlen; 
    long *idx = getRAFASTAindices(fp, &len); 
    char *def;
    for(i = 0; i < len; i++){
        def = getRAFASTAdef(fp, idx[i]); 
        printf(">%s\n", def);
        free(def);
        /* printRAFASTAseq(fp, idx[i], stdout); */
        printf("[%s]\n", getRAFASTAseq(fp, idx[i], &seqlen) );
        }
    fclose(fp); 
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_RAFASTA_C */

