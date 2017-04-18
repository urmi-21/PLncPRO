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

/* parseseqdb: A parser for EMBL, GenBank and Swissprot databases.
   Guy St.C. Slater.  November 1997.
*/

#ifndef INCLUDED_PARSESEQDB_H
#define INCLUDED_PARSESEQDB_H

#include <stdio.h>
#include <time.h>
#include "../general/common.h"
#include "../general/stringlist.h"
#include "../general/string.h"
#include "../parse/parser.h"

/* SECTIONS REQUIRED */
#define DBPARSE_ACCESSION      (1L<< 0)
#define DBPARSE_BASECOUNT      (1L<< 1)
#define DBPARSE_CLASSIFICATION (1L<< 2)
#define DBPARSE_COMMENT        (1L<< 3)
#define DBPARSE_DATE           (1L<< 4)
#define DBPARSE_DBXREF         (1L<< 5)
#define DBPARSE_DEFINITION     (1L<< 6)
#define DBPARSE_DIVISION       (1L<< 7)
#define DBPARSE_FEATURETABLE   (1L<< 8)
#define DBPARSE_GENENAME       (1L<< 9)
#define DBPARSE_IDENTIFIER     (1L<<10)
#define DBPARSE_KEYWORD        (1L<<11)
#define DBPARSE_LENGTH         (1L<<12)
#define DBPARSE_NID            (1L<<13)
#define DBPARSE_ORGANELLE      (1L<<14)
#define DBPARSE_ORGANISM       (1L<<15)
#define DBPARSE_REFERENCE      (1L<<16)
#define DBPARSE_SEGMENT        (1L<<17)
#define DBPARSE_SEQUENCE       (1L<<18)
#define DBPARSE_TYPE           (1L<<19)
#define DBPARSE_ALL            (~(0L))

typedef struct {
 STRING  *feature;
 STRING  *location;
 STRING **name;
 STRING **value;
    int   count;
    } DBFEATURE;

typedef struct {
     STRING  *author;
     STRING  *title;
     STRING  *journal;
     STRING  *position;
 STRINGLIST **dbxrefv;
        int   dbxrefc;
     STRING  *comment;
    } DBREFERENCE;

typedef struct {
  STRINGLIST  *accession;
        int   *basecount;
  STRINGLIST  *classification;
      STRING  *comment;
      time_t   date;
  STRINGLIST **dbxrefv;
         int   dbxrefc;
      STRING  *definition;
        char  *division;
   DBFEATURE **featuretablev;
         int   featuretablec;
      STRING  *genename;
        char  *identifier;
        long   index;
  STRINGLIST  *keyword;
        int   length;
        char  *nid;
        char  *organism;
      STRING  *organelle;
        long   recordlength;
 DBREFERENCE **referencev;
         int   referencec;
         int  *segment;
      STRING  *sequence;
        char   type;
   } DBENTRY;

DBENTRY *GetNextDatabaseEntry(FILE *fp, long mask);
 int printDBENTRY(FILE *fp, DBENTRY *e);
void freeDBENTRY(DBENTRY *e);

#endif /* INCLUDED_PARSESEQDB_H */

