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

/* ntaa: Nucleotide and peptide translation and ambiguity routines 
   Guy St.C. Slater.  Janurary 1996.  Version 2.1
*/

#ifndef INCLUDED_NTAA_H
#define INCLUDED_NTAA_H

#define NTAA_AASETSIZE                        40
#define NTAA_ALPHABETSIZE                 (1<<8) 
#define NTAA_TRANSLATIONNUMBER              4096
#define NTAA_BASESETSIZE                  (1<<4)

extern const unsigned char global_ntaa_aa [NTAA_AASETSIZE];
extern unsigned const char global_ntaa_aapima [NTAA_AASETSIZE];
extern unsigned char global_ntaa_nt2d [NTAA_ALPHABETSIZE];
extern long global_ntaa_trans [NTAA_TRANSLATIONNUMBER];
extern const unsigned char global_ntaa_ntc [NTAA_BASESETSIZE];
extern unsigned char global_ntaa_nt2cd [NTAA_ALPHABETSIZE];
extern const char global_ntaa_aanam  [NTAA_AASETSIZE];

#define NTAAComplement(c) \
    (global_ntaa_ntc[global_ntaa_nt2d[c]])

#define NTAACodonRC2AA(n) (global_ntaa_trans[            \
                          (global_ntaa_nt2cd[n[2]])      \
                         |(global_ntaa_nt2cd[n[1]]<<4)   \
                         |(global_ntaa_nt2cd[n[0]]<<8)])

#define NTAACodon2AA(n) (global_ntaa_trans[           \
                        (global_ntaa_nt2d[n[0]])      \
                       |(global_ntaa_nt2d[n[1]]<<4)   \
                       |(global_ntaa_nt2d[n[2]]<<8)])

#define NTAATranslate(n)     (global_ntaa_aa[NTAACodon2AA(n)])
#define NTAATranslatePima(n) (global_ntaa_aapima[NTAACodon2AA(n)])
#define NTAATranslateRC(n)   (global_ntaa_aa[NTAACodonRC2AA(n)])
#define NTAATransNonAmb(d)   (global_ntaa_aanam[NTAACodon2AA(d)])
#define NTAATransNonAmbRC(d) (global_ntaa_aanam[NTAACodonRC2AA(d)])
#define NTAANonRed(c) (global_ntaa_aanam[c])

#define NTAABase2AA(a,b,c)                                     \
               (global_ntaa_trans[(global_ntaa_nt2d[(a)])      \
                                 |(global_ntaa_nt2d[(b)]<<4)   \
                                 |(global_ntaa_nt2d[(c)]<<8)])
#define NTAATranslateBase(a,b,c) (global_ntaa_aa[NTAABase2AA(a,b,c)])

void NTAAInitialiseData();   /* MUST BE INITIALISED AT RUN TIME */

void revcompNTAA(char *seq, int len);
 int translateNTAAseq(unsigned char *dnaseq, int dnalen, int frame,
                      unsigned char *aaseq);

#endif /* INCLUDED_NTAA_H */

