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

/* rafasta : RANDOM ACCESS TO FASTA FORMAT FILES ROUTINES.
   Pointer is stored to start of sequence, not definition.
   Guy St.C. Slater.. August 1997.
*/
#ifndef INCLUDED_RAFASTA_H
#define INCLUDED_RAFASTA_H

#define RAFASTA_INDEX_CHUNKSIZE 1024
#define RAFASTA_SEQ_CHUNKSIZE   1024
#define RAFASTA_BACKTRACK_GUESS  144 
#define RAFASTA_NO_DEFINITION     -1

long *   getRAFASTAindices(FILE *fp, int *total);
 int   printRAFASTAseq(FILE *fp, long pos, FILE *out);
char *   getRAFASTAseq(FILE *fp, long pos, int *length);
long  locateRAFASTAdefn(FILE *fp, long pos);
void   printRAFASTAdef(FILE *fp, long pos, FILE *out);
char *   getRAFASTAdef(FILE *fp, long pos);

 int  writeRAFASTAindices(FILE *infp, FILE *outfp);
long *readRAFASTAindices(char *path, int *total);

typedef struct {
    char *def;
    char *seq;
     int  len;
    long  pos;
    } RAFASTAALL;

RAFASTAALL *readRAFASTAALL(char *path, int *total);
      void  freeRAFASTAALL(RAFASTAALL *rfa);

#endif /* INCLUDED_RAFASTA_H */

