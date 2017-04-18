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

/* Virtual Finite State Machine Routines.
   Guy St.C. Slater..  September 1998.  Version 2.0
*/

#ifndef INCLUDED_VFSM_H
#define INCLUDED_VFSM_H

#include <stdio.h>
#include <limits.h> /* FOR UINT_MAX etc */
#include "../general/common.h"

typedef unsigned int VFSMTYPE;
#define VFSMTYPE_MAX UINT_MAX
#define VFSMTYPE_PRINT "u"

typedef struct {
         char  index[ALPHABETSIZE]; /* LOOKUP INDEX               */
unsigned char *alphabet;            /* ALPHABET                   */
          int  alphasize;           /* ALPHABET SIZE (BASE)       */
          int  logalphasize;        /* LOG ALPHABET SIZE (BASE)   */
          int  depth;               /* DEPTH                      */
          int  ispoweroftwo;        /* USED FOR &(base-1) TRICKS  */
     VFSMTYPE  prs;                 /* PENULTIMATE ROW START      */
     VFSMTYPE  prw;                 /* PENULTIMATE ROW WIDTH      */
     VFSMTYPE  lrs;                 /* LAST ROW START             */
     VFSMTYPE  lrw;                 /* LAST ROW WIDTH (LEAFCOUNT) */
     VFSMTYPE  total;               /* TOTAL NUMBER OF STATES     */
    } VFSM;

VFSM *newVFSM(char *alphabet, int depth, BOOLEAN casesensitive);
void freeVFSM(VFSM *vfsm);
VFSMTYPE word2posVFSM(VFSM *vfsm, unsigned char *word);
BOOLEAN pos2wordVFSM(VFSM *vfsm, VFSMTYPE pos, char *word);
VFSMTYPE changeVFSMstate(VFSM *vfsm, int state, unsigned char newch);
VFSMTYPE changeVFSMstatePOW2(VFSM *vfsm, int state,
                             unsigned char newch);

/* NB. THIS FUNCTION ASSUMES 32 BIT VFSMTYPE */
VFSMTYPE revcompVFSMpos(VFSM *vfsm, VFSMTYPE pos);

#define MisleafVFSMstate(v, s) ((s) >= (v)->lrs)
#define Mstate2posVFSM(v, s) ((s)-(v)->lrs)

#define MchangeVFSMstate(v, s, c)                 \
    ((s) < (v)->lrs)                              \
    ? ((s)*(v)->alphasize)+(v)->index[(c)] \
    : (((v)->prs+(((s)-(v)->lrs) % (v)->prw))     \
    * (v)->alphasize)+(v)->index[(c)];

#define MchangeVFSMstatePOW2(v, s, c)                   \
    ((s) < (v)->lrs)                                    \
    ?  ((s)<<(v)->logalphasize)+(v)->index[(c)]  \
    :  (((v)->prs+(((s)-(v)->lrs) & (v)->ispoweroftwo)) \
    << (v)->logalphasize)+(v)->index[(c)];

/* ispoweroftwo AND logalphasize
   WILL BOTH BE ZERO IF ALPHABETSIZE IS NOT A POWER OF TWO.
*/


/* START OF FASTA PARSING VFSM CODE */

typedef struct {
    void (*   init)(void *info, int numstates);
    void (*  start)(void *info, int seqid, int pos);
    void (* finish)(void *info, int seqid);
    void (*observe)(void *info, int seqid, int seqpos, int stateid);
    void (*validch)(void *info, int seqid, int seqpos, char symbol);
    } PARSEFASTAVFSMFUNCS;

int parsefastaVFSM(FILE *fp, VFSM *vfsm, void *info,
                            PARSEFASTAVFSMFUNCS *pf);
int parsefastaVFSMrevcomp(FILE *fp, VFSM *vfsm,
                          void *info,   PARSEFASTAVFSMFUNCS *pf,
                          void *inforc, PARSEFASTAVFSMFUNCS *pfrc);

/*   END OF FASTA PARSING VFSM CODE */

# endif /* INCLUDED_VFSM_H */

