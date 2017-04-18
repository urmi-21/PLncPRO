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

#ifndef INCLUDED_NTAA_C
#define INCLUDED_NTAA_C

#include <ctype.h>
#include <string.h>
#include "ntaa.h"
#include "../general/error.h"

ALLOW_ERROR_MESSAGES;

/* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
 0 - 0000 - blank          REPRESENTATION OF NUCLEOTIDES  
 1 G 0001 C                -----------------------------
 2 A 0010 T 
 3 R 0011 Y (AG)       THE STANDARD (IUB) NUCLEOTIDES ARE CODES 
 4 T 0100 A            ARE USED IN THE ARRAY "-GARTKWDCSMVYBHN" 
 5 K 0101 M (GT)       SO THAT:        
 6 W 0110 W (AT)        
 7 D 0111 H (AGT)       1: FOUR BITS PER BASE ARE USED. 
 8 C 1000 G             2: A BIT IS SET FOR EACH BASE REPRESENTED.  
 9 S 1001 S (CG)        3: ANY BASE REVERSED IS IT'S COMPLEMENT.   
10 M 1010 K (AC) 
11 V 1011 B (ACG)  
12 Y 1100 R (CT)       TRANSLATION IS CONSTANT TIME, (LOOKUP TABLE)
13 B 1101 V (CGT)      TAKING INTO ACCOUNT OF ALL OF THESE CODES. 
14 H 1110 D (ACT)
15 N 1111 N (ATGC)
 :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

#define NTAA_GENCODESIZE                 (4*4*4)
#define NTAA_PIMASETSIZE                     18

static char local_ntaa_init = 0;
static unsigned char local_ntaa_nt [NTAA_BASESETSIZE] 
                                  = "-GARTKWDCSMVYBHN";
static unsigned char local_ntaa_aa2d [NTAA_ALPHABETSIZE];
static unsigned char local_ntaa_code  [NTAA_GENCODESIZE]
 = "GGGGEEDDVVVVAAAARRSSKKNNMIIITTTTW*CC**YYLLFFSSSSRRRRQQHHLLLLPPPP";
static long local_ntaa_aamask[NTAA_AASETSIZE];

unsigned char global_ntaa_nt2d [NTAA_ALPHABETSIZE];
long global_ntaa_trans  [NTAA_TRANSLATIONNUMBER];
const unsigned char global_ntaa_aa [NTAA_AASETSIZE] 
                  = "-ARNDCQEGHILKMFPSTWYV*XXXXXXXXXXXXXXXXXX";
const unsigned char global_ntaa_aapima [NTAA_AASETSIZE] 
                  = "-ARNDCQEGHILKMFPSTWYV*ablkonihdmcepjfrxX";
const char global_ntaa_aanam  [NTAA_AASETSIZE] = 
      "\0ARNDCQEGHILKMFPSTWYV*\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";
const unsigned char global_ntaa_ntc [NTAA_BASESETSIZE]
                                   = "-CTYAMWHGSKBRVDN";
unsigned char global_ntaa_nt2cd  [NTAA_ALPHABETSIZE];

static void NTAAInitialiseNucleotideData(){
    register int i; 
    memset(global_ntaa_nt2d, 0, sizeof(global_ntaa_nt2d));
    memset(global_ntaa_nt2cd, 0, sizeof(global_ntaa_nt2cd));
    for(i = 0; i < NTAA_BASESETSIZE; i++)
        global_ntaa_nt2d[local_ntaa_nt[i]] = 
            global_ntaa_nt2d[tolower(local_ntaa_nt[i])] = i;
    for(i = 0; i < NTAA_BASESETSIZE; i++)
        global_ntaa_nt2cd[global_ntaa_ntc[i]] = 
            global_ntaa_nt2cd[tolower(global_ntaa_ntc[i])] = i;
    global_ntaa_nt2d['X'] = global_ntaa_nt2d['x'] 
                          = global_ntaa_nt2d['N'];
    global_ntaa_nt2d['U'] = global_ntaa_nt2d['u'] 
                          = global_ntaa_nt2d['T'];
    global_ntaa_nt2cd['X'] = global_ntaa_nt2cd['x'] 
                           = global_ntaa_nt2cd['N'];
    global_ntaa_nt2cd['U'] = global_ntaa_nt2cd['u'] 
                           = global_ntaa_nt2cd['T'];
    return;
    }

void NTAAInitialisePeptideData(){
    register int i, j; 
    unsigned char pimagrp[NTAA_PIMASETSIZE][6] = {
     "aIV", "bLM",  "dFWY", "lND",   "kDE",  "oEQ",  
     "nKR", "iST",  "hAG",  "cab",   "edH",  "mlk", 
     "pon", "jihP", "fCcd", "rHmpi", "xfrj", "Xx*" };
    memset(local_ntaa_aa2d, 0, sizeof(local_ntaa_aa2d));
    memset(local_ntaa_aamask, 0, sizeof(local_ntaa_aamask));
    for(i = 0; i < NTAA_AASETSIZE; i++)
        local_ntaa_aa2d[global_ntaa_aapima[i]] = i;
    for(i = 1; i < 22; i++) 
        local_ntaa_aamask[i] = (1L<<(i-1));  /* FIRST IS ZERO  */
    for(i = 0; i < NTAA_PIMASETSIZE; i++){
        local_ntaa_aamask[local_ntaa_aa2d[
                          pimagrp[i][0]]] 
      = local_ntaa_aamask[local_ntaa_aa2d[
                          pimagrp[i][1]]]; 
        for(j = 2; pimagrp[i][j]; j++)
            local_ntaa_aamask[local_ntaa_aa2d[
                              pimagrp[i][0]]] 
         |= local_ntaa_aamask[local_ntaa_aa2d[
                              pimagrp[i][j]]]; 
        }
    return;
    }

void NTAAInitialiseTranslationData(){ 
    register char a,b,c,x,y,z;
    register long i,t;
    memset(global_ntaa_trans, 0, sizeof(global_ntaa_trans));
    for(x = 0; x < 16; x++){ 
        for(y = 0; y < 16; y++){
            for(z = 0; z < 16; z++){
                t = 0;
                for(a = 0; a<4; a++){
                    if(x!=(x|1<<a))
                        continue;
                    for(b = 0; b<4; b++){
                        if(y!=(y|1<<b))
                            continue;
                        for(c = 0; c<4; c++){
                            if(z!=(z|1<<c))
                                continue;
                            t = (t|local_ntaa_aamask[
                                local_ntaa_aa2d[
                                local_ntaa_code[
                                ((a<<4)|(b<<2)|c)]]]);
                            }
                        }
                    }
                for(i=0; local_ntaa_aamask[i] 
                   != (t|local_ntaa_aamask[i]); i++);
                global_ntaa_trans[x|(y<<4)|(z<<8)] = i; 
                }
            }
       }
    return; 
    }

void NTAAInitialiseData(){
    if(local_ntaa_init)
        return;
    local_ntaa_init = 1;
    NTAAInitialiseNucleotideData();
    NTAAInitialisePeptideData();
    NTAAInitialiseTranslationData(); 
    return;
    }

void revcompNTAA(char *seq, int len){
    register unsigned char *a, *z, swap;
    for(a = (unsigned char*)seq, z = (unsigned char*)seq+len-1;
        a < z; a++, z--){
        swap = NTAAComplement(*a);
        *a   = NTAAComplement(*z);
        *z   = swap;
        }
    return;
    }

/* translateNTAAseq : RETURNS LENGTH OF aaseq GENERATED.
*/
int translateNTAAseq(unsigned char *dnaseq, int dnalen, int frame,
                     unsigned char *aaseq){
    register unsigned char *dp, *ap = aaseq, *end;
    if((frame > 0) && (frame < 4)){
        end = dnaseq+dnalen-2;
        for(dp = dnaseq+frame-1; dp < end; dp+=3)
            *ap++ = NTAATranslate(dp);
        *ap= '\0';
        return ap-aaseq;
        }       
    if((frame < 0) && (frame > -4)){ 
        for(dp = dnaseq+dnalen+frame-2; dp >= dnaseq; dp-=3)
            *ap++ = NTAATranslateRC(dp);
        *ap= '\0';
        return ap-aaseq;
        }       
    errmsg(ERROR_FATAL, "Invalid reading frame [%d]", frame); 
    return 0; /* NEVER REACHED */
    }


#include <stdio.h>
#include "ntaa.h"

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

int main(){
    unsigned char *seq = "ATG";
    unsigned char seq2[17];
    strcpy(seq2, "-GARTKWDCSMVYBHN");

    NTAAInitialiseData();
    printf("[%s] -> [%c] revcomp -> [%c]\n", 
            seq, NTAATranslate(seq), NTAATranslateRC(seq) );
    seq = "GGX"; 
    printf("[%s] -> [%c] revcomp -> [%c]\n", 
            seq, NTAATranslate(seq), NTAATranslateRC(seq) );
    seq = "NGT"; 
    printf("[%s] -> [%c] revcomp -> [%c]\n", 
            seq, NTAATranslate(seq), NTAATranslateRC(seq) );
    printf("[%s] -> [%c] revcomp -> [%c]\n", 
            seq, NTAATranslate(seq), NTAATranslate(seq) );
    printf("B4:Seq2: [%s]\n", seq2);
    revcompNTAA(seq2, strlen(seq2));
    printf("RC:Seq2: [%s]\n", seq2);
    return 0;
    }
/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_NTAA_C */

