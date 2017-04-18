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

/* readfasta: Fast Fasta format sequence database reading code.   
   Guy St.C. Slater.  December 1996.
*/
#ifndef INCLUDED_READFASTA_H
#define INCLUDED_READFASTA_H

#include <stdio.h>
#include "../general/common.h"

#define READFASTA_CHUNKSIZE   65536 /* HOW MUCH IS READ AT A TIME  */
#define READFASTA_OCCUPANCY      75 /* %BUFFER FULL FOR REALLOC    */

/* READFASTA_TASK_SIZE SHOULD BE MORE THAN READFASTA_CHUNKSIZE */

typedef struct {
    FILE *fp;       /* FILE POINTER                   */
    long  pos;      /* POSITION IN THE FILE AT START  */
    char *buf;      /* DATA BUFFER                    */
    long  alloc;    /* LENGTH ALLOCATED FOR BUFFER    */
    long  indent;   /* LENGTH TO START OF UNREAD DATA */
    long  full;     /* LENGTH OF DATA IN BUFFER       */
    char *seq;      /* POINTS TO START OF SEQUENCE    */
    long  seqlen;   /* LENGTH OF SEQUENCE             */
    char *def;      /* POINTS TO START OF DEFINITION  */
    long  deflen;   /* LENGTH OF DEFINITION LINE      */
    long  filesize; /* LENGTH OF THE FILE             */
    long  stop;     /* POSITION TO STOP READING       */
    long  count;    /* NUMBER OF SEQUENCES READ       */
} READFASTA;

READFASTA *newREADFASTA(char *path);
READFASTA *makeREADFASTA();
void freeREADFASTA(READFASTA *r);
long nextREADFASTA(READFASTA *r);
long  posREADFASTA(READFASTA *r, long pos);
void initREADFASTA(READFASTA *r, char *path);

#define READFASTA2RAFASTApos(r) ((r)->pos + (r)->deflen + 2)

#endif /* INCLUDED_READFASTA_H */

