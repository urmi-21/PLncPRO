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

#ifndef INCLUDED_DICT_C
#define INCLUDED_DICT_C

#include <stdlib.h>
#include "dict.h"

#define DICTWORDCHUNKSIZE 100

/* newDICT : MAKE A NEW DICTIONARY */
DICT *newDICT(){  
    register DICT *n = NEW(DICT); 
    n->d = NULL;
    n->c = '\0';
    n->a = n;     
    return n;
    }

/* insertDICT : INSERT A NEW NODE TO THE DICTIONARY IN THE LINKED 
                LIST AT A POSITION IMMEDIATELY AFTER THE ONE GIVEN. */ 
static DICT *insertDICT(DICT *s, char c){
    DICT *n = NEW(DICT); 
    n->d = NULL;
    n->c = c;
    n->a = s->a;
    return (s->a = n); 
    } 

/* rotateDICT : ROTATE A LIST, AND RETURN THE POSITION HOLDING 
                CHARACTER c, INSERTING IT IF NECESSARY */ 
static DICT *rotateDICT(DICT *n, char c){
    register DICT *init = n; 
    if(!n->c) 
        if(n == n->a) 
            if(!n->d){
                n->c = c;
                return n;
                }
    do{ if(n->c == c)
            return n;
        if(n->c < n->a->c){        
            if( (c > n->c) && (c < n->a->c) )  
                return insertDICT(n, c); 
        } else if ( (c > n->c) || (c < n->a->c) )
                return insertDICT(n, c); 
    } while(init != (n = n->a));
    return (DICT*)0; /* SHOULD NEVER BE REACHED */
    }

/* spinDICT : SPIN A LIST TO LOOK FOR CHARACTER c, 
              RETURN NULL IF ABSENT */ 
static DICT *spinDICT(DICT *n, char c){
    register DICT *init = n; 
    do{ if(n->c == c)
            return n;
        if(n->c < n->a->c){ 
            if( (c > n->c) && (c < n->a->c) )  
                return NULL;
        } else if ( (c > n->c) || (c < n->a->c) )
                return NULL;
    } while(init != (n = n->a));
    return NULL;
    }

/* uniqAddDICT : ADD NAME:VALUE PAIR TO DICTIONARY, IF NAME PRESENT 
                 RETURNS VALUE, OTHERWISE RETURNS NULL. */
void *uniqAddDICT(DICT *n, char *w, void *v){
    do {
        n = rotateDICT(n, *w);
        if(!n->d)
            n->d = newDICT();
        n = n->d;
    } while(*++w); /* ASSUME STRLEN > 0 */
    n = rotateDICT(n, '\0');
    if(n->d)
        return n->d; /* FAIL, RETURN PREV VALUE */
    n->d = v;
    return NULL; /* SUCCEED, RETURN NULL */
    }

/* uniqAddDICTlen : ADD NAME:VALUE PAIR TO DICTIONARY, IF NAME PRESENT 
                    RETURNS VALUE, OTHERWISE RETURNS NULL. 
                    WILL NOT LOOK FOR STRING TERMINATING NULL
*/
void *uniqAddDICTlen(DICT *n, char *w, void *v, int len){
    register char *z = w+len;
    do {
        n = rotateDICT(n, *w);
        if(!n->d)
            n->d = newDICT();
        n = n->d;
        w++;
    } while(w != z); /* ASSUME STRLEN > 0 */
    n = rotateDICT(n, '\0');
    if(n->d)
        return n->d; /* FAIL, RETURN PREV VALUE */
    n->d = v;
    return NULL; /* SUCCEED, RETURN NULL */
    }

/* replaceAddDICT : ADD NAME:VALUE PAIR TO DICTIONARY, REPLACING
                    THE PREVIOUS VALUE, WHICH IS RETURNED */
void *replaceAddDICT(DICT *n, char *w, void *v){
    register void *p;
    do {
        n = rotateDICT(n, *w);
        if(!n->d)
            n->d = newDICT();
        n = n->d;
    } while(*++w);
    n = rotateDICT(n, '\0');
    p = n->d;
    n->d = v;
    return p;
    }

/* mergeAddDICT : ADD OR JOIN w+v TO n USING f.
                  DATA IS PASSED TO f AS NAME VALUE PAIRS,
*/
void mergeAddDICT(DICT *n, char *w, void *v, JOINFUNC f){
    do {
        n = rotateDICT(n, *w);
        if(!n->d)
            n->d = newDICT();
        n = n->d;
    } while(*++w);
    n = rotateDICT(n, '\0');
    if(!n->d){
        n->d = v; /* IF NO VALUE, THEN ADD ON */
        return;
        }
    n->d = f((void*)n->d, v);
    return;
    }

/* countDICT : INCREMENT VALUE FOR w, RETURN THE PREVIOUS VALUE 
*/
long countDICT(DICT *n, char *w){
    register long t;
    do {
        n = rotateDICT(n, *w);
        if(!n->d)
            n->d = newDICT();
        n = n->d;
    } while(*++w);
    n = rotateDICT(n, '\0');
    t = (long)n->d;
    n->d = (void*)((long)n->d+1);
    return t;
    }

/* lookupDICT : RETURNS VALUE CORRESPONDING TO w, OR RETURNS NULL
                IF w CANNOT BE FOUND IN THE DICTIONARY */
void *lookupDICT(DICT *n, char *w){
    register DICT *t;
    do {
        if((t = spinDICT(n, *w)))
            n = t->d;
        else
            return NULL;
    } while(*++w);
    return (t = spinDICT(n, '\0'))?t->d:NULL;
    }

static void walkDICTrecur(DICT *n, char **w, int p,
                         int *a, long *t, WALKFUNC f){
    register DICT *init = n;
    if(p > (*a)-2){
        (*a)+=DICTWORDCHUNKSIZE;
        *w = (char*)realloc(*w, sizeof(char)*(*a));
        }
    do {
      (*w)[p] = n->c;
      if(n->c)
          walkDICTrecur(n->d, w, p+1, a, t, f);
      else if(n->d)
               f(*w, n->d, (*t)++);
    } while(init != (n = n->a));
    return;
    }
/* walkDICT : CALL f AT EVERY NODE IN n */
long  walkDICT(DICT *n, WALKFUNC f){
    int a = DICTWORDCHUNKSIZE;
    long t = 0;
    char *w = (char*)malloc(sizeof(char*)*a);
    walkDICTrecur(n, &w, 0, &a, &t, f);
    free(w);
    return t;
    }

static void walkDICTinforecur(DICT *n, char **w, int p,
                         int *a, long *t, WALKFUNCINFO f, void *info){
    register DICT *init = n;
    if(p > (*a)-2){
        (*a)+=DICTWORDCHUNKSIZE;
        *w = (char*)realloc(*w, sizeof(char)*(*a));
        }
    do {
      (*w)[p] = n->c;
      if(n->c)
          walkDICTinforecur(n->d, w, p+1, a, t, f, info);
      else if(n->d)
               f(*w, n->d, (*t)++, info);
    } while(init != (n = n->a));
    return;
    }

/* walkDICTptr : CALL f AT EVERY NODE IN n, SUPPLYING WALKFUNC with v */
long  walkDICTinfo(DICT *n, WALKFUNCINFO f, void *info){
    int a = DICTWORDCHUNKSIZE;
    long t = 0;
    char *w = (char*)malloc(sizeof(char*)*a);
    walkDICTinforecur(n, &w, 0, &a, &t, f, info);
    free(w);
    return t;
    }

/* -------------------------------------------------------------- */
/* free functions
*/
static void freeDICTrecur(DICT *n, DICT *s){
    if(n->c)
        freeDICTrecur(n->d->a, n->d); /* DOWN */
    if(n != s)
        freeDICTrecur(n->a, s);       /* ACROSS */
    free(n);
    return;
    }

/* freeDICT : calls freeing recursive function.
*/
void freeDICT(DICT *n){
    freeDICTrecur(n->a, n);
    return;
    }

static void freeDICTptr_recur(DICT *n, DICT *s, FREEFUNC ff){
    if(n->c)
        freeDICTptr_recur(n->d->a, n->d, ff); /* DOWN */
    else 
        ff(n->d);
    if(n != s)
        freeDICTptr_recur(n->a, s, ff);       /* ACROSS */
    free(n);
    return;
    }

/* freeDICTptr : FREE DICT AND THE RECORDS CONTAINED WITHIN. 
*/
void freeDICTptr(DICT *n, FREEFUNC ff){
    freeDICTptr_recur(n->a, n, ff);
    return;
    }

/* -------------------------------------------------------------- */

static void sortWalkDICTrecur(DICT *n, char **w, int p,
                         int *a, long *t, WALKFUNC f){
    register DICT *init;
    if(p > (*a)-2){
        (*a)+=DICTWORDCHUNKSIZE;
        *w = (char*)realloc(*w, sizeof(char)*(*a));
        } 
    while(n->c < n->a->c) /* GO TO START OF LIST (FOR SORT) */
        n = n->a;
    init = (n = n->a);
    do {
      (*w)[p] = n->c;
      if(n->c)
          sortWalkDICTrecur(n->d, w, p+1, a, t, f);
      else if(n->d)
               f(*w, n->d, (*t)++);
    } while(init != (n = n->a));
    return;
    }

/* sortWalkDICT : CALL f AT EVERY NODE IN n IN ALPHABETICAL ORDER */
long sortWalkDICT(DICT *n, WALKFUNC f){
    int a = DICTWORDCHUNKSIZE;
    long t = 0;
    char *w = (char*)malloc(sizeof(char*)*a);
    sortWalkDICTrecur(n, &w, 0, &a, &t, f);
    free(w);
    return t;
    }

typedef struct {
    JOINFUNC  join;
        DICT *target;
    } DICTJOINER;

/* joindictfunc : WALK FUNC FOR JOINING DICTIONARY.
*/
static void joindictfunc(char *name, void *value, 
                          int count, void *info){
    register DICTJOINER *j = (DICTJOINER*)info; 
    mergeAddDICT(j->target, name, value, j->join);
    return;
    }

/* joinDICT : COMBINE TWO DICTIONARIES.
*/
DICT *joinDICT(DICT *a, DICT *b, JOINFUNC f){
    DICTJOINER j;
    j.target = a;
    j.join   = f;
    walkDICTinfo(b, joindictfunc, &j);
    freeDICT(b);
    return a;
    }

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */
#include <stdio.h>

static void walktest(char *name, void *value, int count){
    puts(name);
    return;
    }

int main(){
    DICT *d = newDICT();
    char line[100];
    while(fgets(line, 100, stdin))
        if(*line) /* LINE LENGTH MUST BE > 0 */
            listAddDICT(d, line, line, strlen(line));
    sortWalkDICT(d, walktest);
    freeDICTlist(d);
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */

#endif /* INCLUDED_DICT_C */

