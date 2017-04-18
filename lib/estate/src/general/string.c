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

#ifndef INCLUDED_STRING_C
#define INCLUDED_STRING_C

#include <stdlib.h>
#include <string.h>
#include "string.h"
#include "common.h"

#define STRINGCHUNKSIZE 100

STRING *newSTRING(){
    register STRING *s = NEW(STRING); 
    s->alloc = sizeof(char)*STRINGCHUNKSIZE;
    s->str  = (char*)malloc(s->alloc);
    s->len = 0;
    *s->str = '\0';
    return s;
    }

void freeSTRING(STRING *s){
    if(s->str)
        free(s->str);
    free(s);
    return;
    }

STRING *makeSTRING(char *string, int len){
    register STRING *s = NEW(STRING); 
    s->str = (char*)malloc(sizeof(char)*(len+1)); 
    memcpy(s->str, string, sizeof(char)*(len));
    s->alloc = s->len = len; 
    s->str[s->len] = '\0';
    return s;
    }

STRING *copySTRING(STRING *a){
    register STRING *s = NEW(STRING); 
    s->str = (char*)malloc(sizeof(char)*(s->alloc = a->alloc));
    memcpy(s->str, a->str, sizeof(char)*(s->len = a->len));
    return s;
    }

STRING *joinSTRING(STRING *a, STRING *b){
    a->alloc = sizeof(char)*(a->len+b->len+1); 
    a->str = (char*)realloc(a->str, a->alloc); 
    memcpy(a->str+a->len, b->str, sizeof(char)*(b->len+1));
    a->len = a->alloc;
    freeSTRING(b);
    return a;
    }

void linetoSTRING(STRING *s, char *line, int len){
    s->len+=len;
    if(s->len >= s->alloc)
        s->str = (char*)realloc(s->str, (sizeof(char)*
               (s->alloc = (s->len+STRINGCHUNKSIZE))) );
    memcpy(s->str+s->len-len, line, (sizeof(char)*len));
    s->str[s->len] = '\0';
    return;
    }

void addSTRING(STRING *s, char ch){
    if((s->len+1) >= s->alloc)
        s->str = (char*)realloc(s->str,(s->alloc+=STRINGCHUNKSIZE));
    s->str[s->len++] = ch;
    s->str[s->len] = '\0';
    return;
    }

STRING *dlinetoSTRING(STRING *s, char *line, int len, char delim){
    if(!s)
        return makeSTRING(line, len);
    addSTRING(s, delim);
    linetoSTRING(s, line, len); 
    return s;
    }

/* reverseSTRING : REVERSES s; FOR EVEN OR ODD LENGTH STRINGS */
void reverseSTRING(STRING *s){
    register char *a, *z, swap;
    for(a = s->str, z = s->str+s->len-1; a < z; a++, z--)
        Swap(*a,*z,swap);
    return;
    }

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

int main(){
    STRING *p, *s = newSTRING();
    linetoSTRING(s, "hello", 5);
    printf("s->str [%s] s->len[%d]\n", s->str, s->len);
    linetoSTRING(s, " world ", 7);
    p = joinSTRING(s, makeSTRING("the end", 7));
    printf("RESULT [%s]\n", p->str); 
    s = dlinetoSTRING(s, "...and this.", 12, '\0');
    printf("Finally [%s]\n", p->str); 
    freeSTRING(s); 
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_STRING_C */

