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

#ifndef INCLUDED_STRINGLIST_C
#define INCLUDED_STRINGLIST_C

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "common.h"
#include "stringlist.h"

#define STRINGLISTCHUNKSIZE 10

STRINGLIST *newSTRINGLIST(){ /* CREATE STRINGLIST */
    register STRINGLIST *s = NEW(STRINGLIST); 
    s->v = (char**)malloc(sizeof(char*)*(s->a = STRINGLISTCHUNKSIZE));
    s->c = 0;
    return s;
    }

void freeSTRINGLIST(STRINGLIST *s){
    if(!s)return;
    while(s->c){
        free(s->v[--s->c]); 
        }
    free(s->v);
    free(s);
    return;
    }

/* addSTRINGLIST: ADD line OF len TO s.
                  NEED NOT END IN NULL.
*/
void addSTRINGLIST(STRINGLIST *s, char *line, int len){
    if(s->c >= s->a)
        s->v = (char**)realloc(s->v, 
                       sizeof(char*)*(s->a+=STRINGLISTCHUNKSIZE) );
    s->v[s->c] = (char*)malloc(sizeof(char)*(len+1));
    memcpy(s->v[s->c], line, sizeof(char)*len);
    s->v[s->c++][len] = '\0';
    return;
    }

STRINGLIST *copySTRINGLIST(STRINGLIST *s){
    register STRINGLIST *ns = NEW(STRINGLIST); 
    register int len; 
    ns->v = (char**)malloc(sizeof(char*)*(ns->a = s->a) );
    for(ns->c = 0; ns->c < s->c; ns->c++){
          len = strlen(s->v[ns->c]);
          ns->v[ns->c] = (char*)malloc(sizeof(char)*len);
          memcpy(ns->v[ns->c], s->v[ns->c], sizeof(char)*(len+1));
          }
    return ns; 
    }

STRINGLIST *joinSTRINGLIST(STRINGLIST *a, STRINGLIST *b){
    register int len = a->c+b->c;
    if(!a)return b;
    if(!b)return a;
    a->v = (char**)realloc(a->v, (a->a = sizeof(char*)*len)); 
    memcpy(a->v+a->c, b->v, sizeof(char*)*b->c);
    a->c = len;
    free(b->v);
    free(b); 
    return a;
    }

int printSTRINGLIST(FILE *fp, STRINGLIST *s){
    register int i, total = 0;
    if(!s) return 0; 
    for(i = 0; i < s->c; i++)
#ifdef DEBUG
        total += fprintf(fp, "[%3d] : [%s]\n", i, s->v[i]);
#else
        /* total += fputs(s->v[i], fp); */
        total += fprintf(fp, "%s\n", s->v[i]);
#endif /* DEBUG */
    return total;  
    }

STRINGLIST *toSTRINGLIST(char *line, int len, char delim){
    register char *a, *b, *z = line+len;
    STRINGLIST *s = newSTRINGLIST();
    for(a = b = line; a < z; a++){
        while(isspace(*a))a++;
        b = a;
        while(a < z){
            if(*a == delim)
                break; 
            a++;
            }
        addSTRINGLIST(s, b, a-b);
        }
    return s; 
    }

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

int main(){
    STRINGLIST *s = newSTRINGLIST();
    addSTRINGLIST(s, "ONE\n", 4);
    addSTRINGLIST(s, "TWO\n", 4);
    s = joinSTRINGLIST(s, copySTRINGLIST(s));
    s = joinSTRINGLIST(s, toSTRINGLIST("THREE\n FOUR\n ", 12, ' ')); 
    printSTRINGLIST(stdout, s); 
    return 0;
    }

/* SHOULD OUTPUT:
   "ONE TWO ONE TWO THREE FOUR"
*/


/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_STRINGLIST_C */


