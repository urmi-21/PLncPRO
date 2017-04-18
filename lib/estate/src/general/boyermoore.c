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

#ifndef INCLUDED_BOYERMOORE_C
#define INCLUDED_BOYERMOORE_C

#include "boyermoore.h"
#include <stdlib.h>

#ifndef ALPHABETSIZE
#define ALPHABETSIZE (1<<8) 
#endif /* ALPHABETSIZE */

/* BMindex : MAKES AN INDEX FOR BOYERMOORE FROM PATTERN pat OF 
   LENGTH len.  THE CALLER IS REQUIRED TO FREE THIS INDEX.
   FOR BOYER-MOORE WITH SINGLE BACKTRACKING INDEX.
*/
int *BMindex(unsigned char *pat, int len){
    register int i;
    register unsigned char *a, *z;
    register int *index = (int*)malloc(sizeof(int)*(ALPHABETSIZE+1));
    int skip[ALPHABETSIZE];
    for(i = 0; i < ALPHABETSIZE; i++){
        index[i] = '\0';
        skip[i] = len;
        }
    for(i = 0, a = pat, z = pat+len; a < z; a++)
        if(!index[*a])
            index[*a] = ++i; 
    for(a = pat; a < z; a++)
        skip[index[*a]] = len-(a-pat)-1;
    for(i = 0; i < ALPHABETSIZE; i++)
        index[i] = skip[index[i]];
    index[ALPHABETSIZE] = len; /* STORE LENGTH AT END */
    return index;
    }

/* BMsearch : LOOKS FOR PATTERN a IN STRING a.        
              BASED ON SEDGEWICK p288. FIXED BUG for(... j>=0 ...) 
              SIMPLIFIED IMPLEMENTATION, USING POINTERS.
              RETURNS POINTER TO MATCH IN a, OR NULL IF p IS ABSENT.
*/ 

#define BOYER_MOORE_WITH_POINTERS 
#ifdef  BOYER_MOORE_WITH_POINTERS 

unsigned char *BMsearch(int *index, unsigned char *p,
                                    unsigned char *a, int alen){ 
    register int s, t, m=index[ALPHABETSIZE];
    register unsigned char *ap, *pp, *az, *pz;
    for(ap = a+m-1, az = a+alen, pp = pz = p+m-1; pp>=p; ap--, pp--)
        while(*ap != *pp){
            t = index[*ap];
            s = pz-pp+1;
            ap += (s>t)?s:t;
            if(ap >= az)
                return NULL; /* FAIL */
            pp = pz;
            }
    return ap+1; /* SUCCEED */
    } 

#else /* BOYER_MOORE_WITH_POINTERS */

char *BMsearch(int *index, unsigned char *p, unsigned char *a,
               int alen){ 
    register int i, j, t, m=index[ALPHABETSIZE];
    for(i=m-1,j=m-1; j>=0;i--,j--)
        while(a[i] != p[j]){
            t = index[a[i]];
            i += (m-j>t)?m-j:t;
            if(i>=alen)
                return NULL; /* FAIL */
            j = m-1;
            }
    return a+i+1; /* SUCCEED */
    } 
#endif /* BOYER_MOORE_WITH_POINTERS */

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

#include <stdio.h>

int main(){
#define LONGTEST
#ifdef LONGTEST
    char *pat = "sting";
    char *txt = "a string searching example consisting of ..."; 
#else /* LONGTEST */ 
    char *pat = "and"; 
    char *txt = "ndndand."; 
#endif /* LONGTEST */
    int *index = BMindex(pat, strlen(pat)); 
    char *result = BMsearch(index, pat, txt, strlen(txt) );
    free(index);

    printf("Result = \"%s\"\n", result?result:"<failed>"); 
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_BOYERMOORE_C */

