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

#ifndef INCLUDED_EDITDIST_C
#define INCLUDED_EDITDIST_C

#include <stdio.h>
#include <stdlib.h>
#include "editdist.h"
#include "../general/common.h"

EDITDIST *newEDITDIST(char *qyseq, int qylen){
    register size_t qmem = sizeof(char)*(qylen+1);
    register size_t vmem = sizeof(int)*(qylen+1);
    register EDITDIST *ed = NEW(EDITDIST); 
    ed->query = (char*)malloc(qmem);
    ed->vectorA = (int*)malloc(vmem);
    ed->vectorB = (int*)malloc(vmem);
    memcpy(ed->query, qyseq, qmem);
    ed->length  = qylen;
    return ed;
    }

void freeEDITDIST(EDITDIST *ed){
    free(ed->query);
    free(ed->vectorA);
    free(ed->vectorB);
    free(ed);
    return;
    }

/* calcEDITDIST: CALCULATE EDIT DISTANCE.
                 THE ALGORITHM HERE IS O(mn) COULD BE O(kn).
                 SEE Ukkonnen1995 FOR REFERENCE TO REPLACEMENT.
*/
int calcEDITDIST(EDITDIST *ed, char *dbseq, int dblen){
    register int i, j;
    register int *curr = ed->vectorA, *prev = ed->vectorB, *swap;
    for(i = 0; i <= ed->length; i++)
        prev[i] = i;
    for(i = 1; i <= dblen; i++){
        curr[0] = i;
        for(j = 1; j <= ed->length; j++)
            if((dbseq[i-1] == ed->query[j-1])
                || (curr[j-1] < prev[j-1])
                || (prev[j] < prev[j-1]) )
                curr[j] = prev[j-1];
            else
                curr[j] = prev[j-1] + 1;
        Swap(curr,prev,swap);
        }
    return prev[ed->length];
    }

/* calcEDITDIST: CALCULATE EDIT DISTANCE ONLINE.
*/
int calcEDITDISTonline(EDITDIST *ed, FILE *fp, long pos){
    register int i, j, ch;
    register int *curr = ed->vectorA, *prev = ed->vectorB, *swap;
    for(i = 0; i <= ed->length; i++)
        prev[i] = i;
    fseek(fp, pos, SEEK_SET);
    i = 1;
    do {
        switch(ch = getc(fp)){
            case 'A': case 'B': case 'C': case 'D': case 'E':
            case 'F': case 'G': case 'H': case 'I': case 'J':
            case 'K': case 'L': case 'M': case 'N': case 'O':
            case 'P': case 'Q': case 'R': case 'S': case 'T':
            case 'U': case 'V': case 'W': case 'X': case 'Y':
            case 'Z': /* VALID CHARACTER */
                curr[0] = i;
                for(j = 1; j <= ed->length; j++)
                    if((ch == ed->query[j-1])
                        || (curr[j-1] < prev[j-1])
                        || (prev[j] < prev[j-1]) )
                        curr[j] = prev[j-1];
                    else
                        curr[j] = prev[j-1] + 1;
                Swap(curr,prev,swap);
                break;
            case '>': case EOF:
                return prev[ed->length];
            }
        i++;
    } while(TRUE);
    /* WILL ONLY RETURN AFTER "\n>" OR EOF */
    }

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

int main(){
/* ADD CHECKING CODE HERE */
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_EDITDIST_C */

