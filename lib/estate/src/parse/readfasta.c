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

#ifndef INCLUDED_READFASTA_C
#define INCLUDED_READFASTA_C

#include <string.h> /* FOR memmove() */
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "readfasta.h"
#include "../general/error.h"

ALLOW_ERROR_MESSAGES;

/* makeREADFASTA : CREATE THE BASIC COMPONENTS OF A FASTA READER.
*/
READFASTA *makeREADFASTA(){
    register READFASTA *r = BLANK(READFASTA);
    r->alloc = READFASTA_CHUNKSIZE;
    r->buf = (char*)malloc(sizeof(char)*r->alloc);
    return r;
    }

/* initREADFASTA : INITIALSE BEFORE READING
*/
void initREADFASTA(READFASTA *r, char *path){
    register char *a = r->buf, *z;
    r->full = fread(r->buf,sizeof(char),r->alloc,r->fp); /* READ */
    if(*r->buf != '>')
        for(z = a+r->full-1;
           ( (a[0] != '\n') || (a[1] != '>') ); a++)
            if(a == z)
                errmsg(ERROR_FATAL, "No sequences in [%s]", path);
    r->indent = a-r->buf;
    return;
    }

/* newREADFASTA : PREPARE FOR READING FASTA FORMAT FILE path.
*/
READFASTA *newREADFASTA(char *path){
    register READFASTA *r = makeREADFASTA();
    struct stat s;
    if(!strcmp(path, "-"))
        r->fp = stdin;
    else if(!(r->fp = fopen(path, "r")) )  /* OPEN FILE */
             errmsg(ERROR_FATAL, "Could not open database [%s]", path);
    initREADFASTA(r, path);
    fstat(fileno(r->fp), &s);
    r->stop = r->filesize = s.st_size;
    return r;
    }

void freeREADFASTA(READFASTA *r){
    fclose(r->fp);
    free(r->buf);
    free(r);
    return;
    }

/* getMoreFastaData : READ SOME MORE OF THE FASTA FILE.
   BEHAVES AS IF r->stop IS THE END OF THE FILE.
   RETURNS FALSE IF END OF FILE (OR DIVISION) HAS BEEN REACHED. 
*/

static BOOLEAN getMoreFastaData(READFASTA *r, char **src, char **dst){
    register char *before;
    register long  diff = 0, newread = 0, request, remaining;

   if(    ((float)READFASTA_OCCUPANCY/100)
       <= ((float)(r->full-r->indent)/r->alloc) ){ /* REALLOCATE */
        r->alloc+=READFASTA_CHUNKSIZE;
        before = r->buf;
        r->buf=(char*)realloc(r->buf, sizeof(char*)*r->alloc);
        diff = before-r->buf;
    } else {   /* SHIFT BACK */
        diff = r->indent;
        r->indent = 0;
        r->full-=diff; /* UP FROM BELOW NEXT LINE */
        memmove(r->buf, r->buf+diff, r->full);
        }
    (*src)-=diff;
    (*dst)-=diff;

    request = r->alloc-r->full;
    remaining = r->stop-ftell(r->fp);
    if(remaining <= 0)
        return FALSE;
    newread=fread(*src, sizeof(char), Min(request,remaining), r->fp);
    r->full+=newread;
    return newread?TRUE:FALSE;
    }

long nextREADFASTA(READFASTA *r){
    char *dst, *src = r->buf+r->indent;
    do {  /* GET DEFINITION LINE */ 
        if(src >= r->buf+r->full)
           if(!getMoreFastaData(r, &src, &dst))
               return 0; 
    } while(*src++ != '\n');
    *(src-1) = '\0';
    r->deflen = src-r->buf-r->indent-2;
    dst = src; 
    do {                       /* GET SEQUENCE DATA   */
        if(src >= r->buf+r->full)
            if(!getMoreFastaData(r, &src, &dst))
                *src = '\0';
        switch(*src = toupper(*src)){ 
            case '>':
                if(*(src-1) != '\n')
                    break; /* ( ELSE FALL THROUGH ) */
            case '\0':
                *dst = '\0';
                r->def = r->buf+r->indent+1; 
                r->pos = ftell(r->fp)-r->full+r->indent;
                r->seq = r->def+r->deflen+1;  
                r->seqlen = dst-r->seq;
                r->indent = src-r->buf;
                r->count++;
                return r->seqlen; 
            case 'A': case 'B': case 'C': case 'D': case 'E':
            case 'F': case 'G': case 'H': case 'I': case 'J':
            case 'K': case 'L': case 'M': case 'N': case 'O':
            case 'P': case 'Q': case 'R': case 'S': case 'T':
            case 'U': case 'V': case 'W': case 'X': case 'Y':
            case 'Z': 
                  *dst++ = *src;
                  break;
            }
        src++; 
    } while( TRUE );
    }

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

int main(int argc, char **argv){
    READFASTA *r = newREADFASTA( ((argc>1)?argv[1]:"-") );
    register int i = 0;

    while(nextREADFASTA(r))
         i++;
    printf("Found %d sequences\n", i);

/*  while(nextREADFASTA(r))
        printf("%7d:%.70s\n", ++i, r->def);
*/

/*   while(nextREADFASTA(r)){
        putc('>', stdout);
        fwrite(r->def, sizeof(char), r->deflen, stdout);  
        printf("[%d]\n", r->seqlen);
        writeFASTAblock(r->seq, r->seqlen);
        }
*/
    putc('\n', stdout);
    freeREADFASTA(r);
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_READFASTA_C */

