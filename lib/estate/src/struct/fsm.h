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

/* Efficient Finite State Machine Routines.
   Guy St.C. Slater..  March 1997.
*/

#ifndef INCLUDED_FSM_H
#define INCLUDED_FSM_H

#include "../general/common.h"

typedef struct FSMNODE {
    struct FSMNODE *next;
    union { long l; void *n; } data;
    } FSMNODE;

typedef struct {
unsigned char   index[256]; /* MAP ANY CHARACTER TO THE ARRAY   */
         char   width;      /* ALPHABET SIZE FOR THIS FSM       */
      FSMNODE  *root;       /* BASE OF THE FINITE STATE MACHINE */
          int   alloc;      /* NUMBER OF NODES ALLOCATED        */ 
          int   count;      /* NUMBER OF NODES USED             */ 
    } FSM; 

 FSM *makeFSM(char **words, long total); 
 FSM *wordFSM(unsigned char *str, long wlen); 
long  countFSM(FSM *f, unsigned char *str);
void  infoFSM(FSM *f);
void  freeFSM(FSM *f);
void  freeFSMptr(FSM *f, FREEFUNC ff);
 FSM *copyFSM(FSM *f);

FSM *newFSM(char *alpha);
void initFSM(FSM *f);
void compileFSM(FSM *f);
void submitFSM(FSM *f, unsigned char *s, long length);
void submitFSMptr(FSM *f, unsigned char *s, long length, void *v);
void submitFSMptrJOIN(FSM *f, unsigned char *s, long length,
                      void *v, JOINFUNC combine);
void compileFSMptr(FSM *f);
 FSM *makeFSMptr(char **words, long total, void **ptr);
 FSM *makeFSMptrJOIN(char **words, long total, 
                     void **ptr, JOINFUNC combine);

/* HISTOGRAM ROUTINES */

typedef struct { 
    long count; /* OCCURENCES OF WORD IN QUERY.            */
    long id;    /* IDENTITY OF SEQUENCE DURING LAST VISIT. */
    long match; /* NUMBER OF MATCHES MADE.                 */
    } FSMHISTNODE;

typedef struct {
    FSM         *fsm;    /* THE BASIC FINITE STATE MACHINE. */
    FSMHISTNODE *hist;   /* THE DISTANCE CALCULATION NODES. */
    } FSMHIST;

FSMHIST *makeFSMHIST(FSM *f);
  void   freeFSMHIST(FSMHIST *h);
  long   pairsFSMHIST(FSMHIST *h, unsigned char *seq, long id);

# endif /* INCLUDED_FSM_H */

