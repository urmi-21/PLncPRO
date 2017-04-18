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

/* --------------------------------------------------------------- */
/* common.h : SHARED DEFINITIONS                                   */
/* --------------------------------------------------------------- */

#ifndef INCLUDED_COMMON_H
#define INCLUDED_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ALPHABETSIZE (sizeof(char)<<8)

#ifndef BOOLEAN
typedef unsigned char BOOLEAN;
#endif /* BOOLEAN */

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif /* TRUE */

/* ----------------------- */
/* BIT MANIPULATION MACROS */

/* BitOn : TURN THE nth BIT OF i ON */
#define BitOn(i,n) ((i)|(1<<(n)))

/* BitOff : TURN THE nth BIT OF i OFF */
#define BitOff(i,n) ((i)&~(1<<(n)))

/* BitToggle : CHANGE VALUE OF nth BIT OF i */
#define BitToggle(i,n) ((i))^(1<<(n))

/* SetBitOn : SET THE nth BIT OF i ON */
#define SetBitOn(i,n) (i)|=(1<<(n))

/* SetBitOnly : SET i TO HAVE ONLY THE nth BIT ON */
#define SetBitOnly(i,n) (i)=(1<<(n))

/* SetBitOff : SET THE nth BIT OF i OFF */
#define SetBitOff(i,n) (i)&=(~(1<<(n)))

/* SetBitToggle : CHANGE VALUE OF nth BIT OF i */
#define SetBitToggle(i,n) (i)^+(1<<(n))

/* GetBit : GET BIT n OF i */
#define GetBit(i,n) (((i)>>(n))&1)

/* SHOWS Z BITS WHICH OCCUR Y BITS FROM THE RIGHT IN X 
 * eg. GetBits(x, 0, 1) returns right most bit of x   */
#define GetBits(x,y,z) (((x)>>(y))&~(~0<<(z)))

/* SET Z BITS WHICH OCCUR Y BITS FROM THE RIGHT IN X TO W */
#define SetBits(w,x,y,z) (((x)&~(~(~0<<(z))<<(y)))|((w)<<(y)))

/* NotSingleBit : TEST IF A SINGLE BIT OF n IS ON */
#define NotSingleBit(n) (n&(n-1))

/* RemoveRightBit : REMOVE RIGHT MOST ON BIT OF n */
#define RemoveRightBit(n) (n&=(n-1))

/* isOdd : TEST IF IS AN ODD NUBMER */ 
#define isOdd(n) ((n)&1)

/* ------------------------ */
/* MEMORY ALLOCATION MACROS */

#define NEW(type)  (type*)malloc(sizeof(type))
#define BLANK(type) (type*)calloc(1, sizeof(type))

/* Set of macros for dynamic array allocation without keeping
   track of number of elements already allocated.

     InitAlloc : Start allocation.
     PrepAlloc : Start allocation with ctr elements
   StepRealloc : Increment allocation for ptr with ctr elements.
   AutoRealloc : Realloc ptr with ctr elements for inc extra elements. 

Will maintain allocation of ((ctr % ALLOCSTEP)*ALLOCSTEP) members.
Efficient, as (ALLOCSTEP+1) must be a power of two.
*/
#define ALLOCDIMENTION 8
#define ALLOCSTEP ((1<<(ALLOCDIMENTION))-1)

#define InitAlloc(type) (malloc(sizeof(type)*ALLOCSTEP))

#define PrepAlloc(type, ctr) (malloc(sizeof(type)   \
                             *(((ctr)|ALLOCSTEP)+1)))

#define StepRealloc(ptr, ctr, type)  \
    if(!((ctr) & ALLOCSTEP))         \
       (ptr) = realloc((ptr), sizeof(type)*(((ctr)|ALLOCSTEP)+1))

#define AutoRealloc(ptr, ctr, inc, type)                       \
    if(((ctr) | ALLOCSTEP) != ( ((ctr)+(inc)) | ALLOCSTEP) )   \
        (ptr) = realloc((ptr), sizeof(type)*                   \
                            ( ( ((ctr)+(inc)) | ALLOCSTEP)+1) )

/* ----------------- */
/* COMPARISON MACROS */

#define Min(x,y) (((x)<(y))?(x):(y))
#define Max(x,y) (((x)>(y))?(x):(y))
#define Challenge(a,b) if((a)<(b)) (a)=(b)

#define Swap(x,y,temp) ((temp)=(x),(x)=(y),(y)=(temp))
/* the XOR method is neater, but slower */

#define Negate(n) (n)=-(n)

/* FILE PERMISSIONS */
#define FILE_ALLOW_ALL_R   S_IRUSR|S_IRGRP|S_IROTH
#define FILE_ALLOW_ALL_W   S_IWUSR|S_IWGRP|S_IWOTH
#define FILE_ALLOW_ALL_X   S_IXUSR|S_IXGRP|S_IXOTH

#define FILE_ALLOW_ALL_RW  FILE_ALLOW_ALL_R|FILE_ALLOW_ALL_W
#define FILE_ALLOW_ALL_RX  FILE_ALLOW_ALL_R|FILE_ALLOW_ALL_X
#define FILE_ALLOW_ALL_RWX FILE_ALLOW_ALL_RW|FILE_ALLOW_ALL_X
#define FILE_ALLOW_U_W_UGO_RX FILE_ALLOW_ALL_R|FILE_ALLOW_ALL_X|S_IWUSR
#define FILE_ALLOW_U_W_UGO_R FILE_ALLOW_ALL_R|S_IWUSR

/* ------------- */
/* OFFSET MACROS */
#include <stddef.h> /* FOR offsetof() */
#ifndef offsetof
#define offsetof(type, mem) ((size_t) \
                ((char *)&((type *)0)->mem - (char *)(type *)0))
#endif /* offsetof */
#define OFFSET_ITEM(type, offset, instance) \
    (*(type*)((char*)(instance) + (offset)))

/* TYPE CONVERSION UNION STRUCTURE */

typedef struct {
    char type; /*  USE LETTER TO INDICATE TYPE */
    union {
        char  c;
        char *s;
         int  i;
        void *v;
       float  f;
      double  d;
        } u;
    } ANYTYPE;

/* GENERAL FUNCTION TYPEDEFS */
typedef void   (*FREEFUNC)     (void *v);
typedef void * (* NEWFUNC)     ();
typedef void * (*JOINFUNC)     (void *a, void *b);
typedef  int   (*PRINTFUNC)    (FILE *fp, void *v);
typedef void   (*PARSEFUNC)    (FILE *fp, ANYTYPE *a);
typedef void   (*WALKFUNC)     (char *name, void *value, int count);
typedef void   (*WALKFUNCINFO) (char *name, void *value, int count,
                                void *info);
typedef int    (*COMPFUNC)     (void *low, void *high);
typedef BOOLEAN (*PQCOMPFUNC)  (void *low, void *high);
typedef    void (*PQPOPFUNC)   (void *data, int count,
                                int total, void *info);

/* OTHER GENERAL TYPES */
typedef struct { char *name; void *value; } NAMEVALUEPAIR;

#endif /* INCLUDED_COMMON_H */


