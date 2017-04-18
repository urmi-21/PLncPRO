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

#ifndef INCLUDED_IOUTIL_C
#define INCLUDED_IOUTIL_C

#include <stdio.h>
#include <stdlib.h>
#include "ioutil.h"
#include <sys/types.h>
#include <sys/stat.h>

void catFILE(FILE *in, FILE *out){
    register int ch;
    while((ch = getc(in)) != EOF)
        putc(ch, out);
    return;
    }

void bitprint(FILE *fp, int i, int n){
    while(n--) putc('0' + GetBit(i,n), fp);
    return;
    }

void reversestring(char *str, int len){
    register char *a, *z, swap;
    for(a = str, z = str+len-1; a < z; a++, z--)
        Swap(*a,*z,swap);
    return;
    }

char *readFILE(char *path, long *length){
    register FILE *fp;
    register char *p;
    struct stat s;
    if(!(fp = fopen(path, "r")))
        return NULL;
    if(fstat(fileno(fp), &s) == -1){
        fclose(fp);
        return NULL;
        }
    if(!(p = (char*)malloc(sizeof(char)*(s.st_size+1)))){
        fclose(fp);
        return NULL;
        }
    if(!(*length = fread(p, sizeof(char), s.st_size, fp))){
        free(p);
        fclose(fp);
        return NULL;
        }
    p[*length] = '\0';
    fclose(fp);
    return p;
    }

/* seekFILE : MOVES fp TO END OF THE NEXT OCCURENCE OF s */
/*            RETURNS EOF IF NOT FOUND                   */
/*            RETURNS POSITIONS MOVED OTHERWISE EOF      */
int seekFILE(FILE *fp, char *s){ 
    register char *ptr = s;   
    register int ch; 
    register int start = ftell(fp);
    ch = getc(fp);
    do {
        if(ch == EOF)
            return EOF; 
        /* else (*ptr == ch)?(ptr++):(ptr=s); */
        /* ADDED CASE FOR MISMATCH ON FIRST PATTERN CHARACTER */
        else (*ptr == ch)?(ptr++):(ptr=(*s == ch)?s+1:s);
    ch = getc(fp);
    } while(*ptr);
    return( ftell(fp) - start );
    }

/* seekFILEpipe : MOVES fp TO END OF THE NEXT OCCURENCE OF s
                  VERSION TO WORK WITH PIPES - WILL NOT CALL ftell()
*/
int seekFILEpipe(FILE *fp, char *s){ 
    register char *ptr = s;   
    register int ch, ctr = 1;
    ch = getc(fp);
    do {
        if(ch == EOF)
            return EOF;
        else (*ptr == ch)?(ptr++):(ptr=(*s == ch)?s+1:s);
    ch = getc(fp);
    ctr++;
    } while(*ptr);
    return ctr;
    }

#include <unistd.h> /* FOR gethostname() */
#include <time.h>   /* FOR time() */

#define MAXHOSTNAME_LENGTH 100

char *tmpFILEpath(char *base){
    register pid_t pid = getpid(), ppid = getppid();
    register time_t sectime = time(NULL);
    register char *tmppath;
    char hostname[MAXHOSTNAME_LENGTH];
    gethostname(hostname, MAXHOSTNAME_LENGTH);
    hostname[MAXHOSTNAME_LENGTH-1] = '\0';
    if(!base)
        base = "/tmp";
    tmppath = malloc(sizeof(char)
                   *(strlen(base)+strlen(hostname)+100));
    sprintf(tmppath, "%s/%s_%d_%d_%d.tmpfile",
                     base, hostname,
                     (int)pid, (int)ppid, (int)sectime);
    return tmppath;
    }

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

int main(){
    long i;
    register char *data = readFILE("/home/guy/.signature", &i);
    register char *tmppath = tmpFILEpath(NULL);
    printf("DATA ---[\n%s\n]---\n", data);
    free(data); 
    printf("have tmp path [%s]\n", tmppath);
    free(tmppath);
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_IOUTIL_C */

