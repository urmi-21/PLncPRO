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

/* parsedbest: CONCISE DBEST DATABASE PARSER
   Guy St.C. Slater..  Version 1.1  August 1998.
*/

#ifndef INCLUDED_PARSEDBEST_C
#define INCLUDED_PARSEDBEST_C

#include "parsedbest.h"

#include "../general/error.h"
#include "../general/string.h"
#include "../struct/list.h"

#include <stdlib.h>  /* FOR atoi()    */
#include <string.h>  /* FOR strdup()  */
#include <strings.h> /* FOR ffs()     */
#include <ctype.h>   /* FOR isspace() */

ALLOW_ERROR_MESSAGES;

/* DEBUGGING OPTIONS ETC. */
/* #define DBEST_TEST_PARSER */
/* #define DBEST_PRINT_MISSED_FIELDS */
/* #define DBEST_SHOW_FORMAT */
#define DBEST_MONITOR_USAGE

/* PARSEDBEST_MULTIPLE_CHUNKSIZE MUST BE POWER OF 2 */
#define PARSEDBEST_MULTIPLE_CHUNKSIZE 8

/* ############################################################### */
/* START OF TAG PROCESSING FUNCTIONS                               */

/* PUT GENERIC FUNCTIONS HERE THESE IN THE IOUTIL SECTION    */
/* MAYBE REMOVE THE USE OF STRING OBJECT HERE FOR EFFICIENCY */

static void parseDBESTword(FILE *fp, ANYTYPE *a){
    register int ch;
    register STRING *s = newSTRING();
    a->type = 's';
    while((ch = getc(fp)) != EOF){
        switch(ch){
            case ' ': /* SKIP WHITE SPACE */
            case '\t':
                break;
            case '\n':
                a->u.s = s->str;
                free(s); /* FREE STRING SHELL ONLY */
                return;
            default:
                addSTRING(s, ch);
                break;
            }
        }
    return;
    }

static void parseDBESTline(FILE *fp, ANYTYPE *a){
    register int ch;
    register STRING *s;
    a->type = 's';
    while(isspace(ch = getc(fp))) /* SKIP LEADING WHITE SPACE */
        if(ch == EOF){
            a->u.s = NULL;
            return;
            }
    s = newSTRING();
    do {
        switch(ch){
            case '\t': /* REPLACE TABS WITH SPACES */
                addSTRING(s, ' ');
                break;
            case '\n':
                a->u.s = s->str;
                free(s); /* FREE STRING SHELL ONLY */
                return;
            default:
                addSTRING(s, ch);
                break;
            }
    } while((ch = getc(fp)) != EOF);
    a->u.s = NULL;
    freeSTRING(s);
    return;
    }

static void parseDBESTblock(FILE *fp, ANYTYPE *a){
    register STRING *s = newSTRING();
    register int ch, prev = ' ';
    a->type = 's';
    while( (ch = getc(fp)) != EOF){
        switch(ch){
            case '\n':
                if(prev == '\n'){
                    a->u.s = s->str; 
                    free(s); /* FREE SHELL ONLY */
                    return; 
                    }
                break;
            case ' ':
            case '\t':
                break;
            default:
                addSTRING(s, ch);
                break;
            }
        prev = ch;
        }
    freeSTRING(s);
    a->u.s = NULL;
    return;
    }

static void parseDBESTpara(FILE *fp, ANYTYPE *a){
    register STRING *s = newSTRING();
    register int ch, prev = ' ';
    a->type = 's';
    while(isspace(ch = getc(fp)));
    do {
        switch(ch){
            case '\n':
                if(prev == '\n'){ /* SECTION END */
                    a->u.s = s->str; 
                    free(s); /* FREE SHELL ONLY */
                    return; 
                    }
                break;
            case ' ': /*FALLTHROUGH*/
            case '\t':
                if(prev == '\n'){ /* SKIP SPACE AT START OF LINE */
                    while(isspace(ch = getc(fp))){
                        if((ch == '\n') && (prev == '\n')){
                            a->u.s = s->str;
                            free(s); /* FREE SHELL ONLY */
                            return;
                            }
                        prev = ch;
                        }
                    addSTRING(s, ' ');
                    }
                addSTRING(s, ch);
                break;
            default:
                if(prev == '\n'){  /* TAG END */
                    ungetc(ch, fp);
                    a->u.s = s->str;  
                    free(s); /* FREE SHELL ONLY */
                    return;
                    }
                addSTRING(s, ch);
                break;
            }
        prev = ch;
    } while( (ch = getc(fp)) != EOF);
    freeSTRING(s);
    a->u.s = NULL;
    return;
    }

static void parseDBESTint(FILE *fp, ANYTYPE *a){
    register char *p;
    parseDBESTword(fp, a);
    p = a->u.s;
    a->type = 'i';
    a->u.i  = atoi(a->u.s);
    free(p);
    return;
    }

static void parseDBESTdouble(FILE *fp, ANYTYPE *a){
    register char *p;
    parseDBESTword(fp, a);
    p = a->u.s;
    a->type = 'd';
    a->u.d  = strtod(p, NULL);
    /* NEED TO USE strtod() FOR e NOTATION */
    free(p);
    return;
    }

static int printDBESTstring(FILE *fp, void *v){
    return fprintf(fp, "%s\n", (char*)v);
    }

static int printDBESTint(FILE *fp, void *v){
    return fprintf(fp, "%d\n", (int)v);
    }

static int printDBESTdouble(FILE *fp, void *v){
    register double *d = (double*)&v;
    return fprintf(fp, "%.3e\n", *d);
/* TRY return fprintf(fp, "%.3e\n", *(double *) v); */
    }

#define DBEST_BLOCK_WIDTH 60
#define DBEST_BLOCK_WRAP  10

static char *local_dbest_indent = "                ";

static int printDBESTblock(FILE *fp, void *v){
    register int total = 0, rem, prt;
    register char *a, *z; 
    for(a = (char*)v, z = a+strlen(a); a < z; a+=DBEST_BLOCK_WIDTH){
        if(total) /* TIDY HERE */
            total += fwrite(local_dbest_indent, sizeof(char), 16, fp);
        rem = z-a;
        prt = Min(DBEST_BLOCK_WIDTH, rem);
        total += fwrite(a,    sizeof(char),  prt, fp);
        total += fwrite("\n", sizeof(char),  1, fp);
        }
    return total + fwrite("\n", sizeof(char),  1, fp);
    }

static int printDBESTpara(FILE *fp, void *v){
    register char *a = (char*)v, *z = a+strlen(a), *go, *stop = a-1;
    register char *mem;
    register int total = 0;
    do {
        go = stop+1;
        if((stop+=DBEST_BLOCK_WIDTH) > z) /* LESS THAN A LINE */
            stop = z;
        else {
            mem = stop;
            while(!isspace(*stop)){ /* MOVE BACK TO LAST SPACE */
                if(stop == go){     /* NO SPACES ON LINE       */
                    stop = mem;
                    break;
                    }
                stop--;
                }
            }
        if(total) /* TIDY HERE */
            total += fputs(local_dbest_indent, fp);
        total += fwrite(go, sizeof(char), stop-go, fp);
        total += fputc('\n', fp);
    } while(stop < z);
    return total;
    }

/* ############################################################### */
/* TAG TYPE DEFINITIONS                                           */

typedef struct {     /* NEW FUNCTIONS CALLED FROM WITHIN PARSE */
    PARSEFUNC parse;
    PRINTFUNC print;
     FREEFUNC free;  /* MAY BE NULL IF NOT REQUIRED       */
       size_t size;  /* SIZE OF TYPE (OR POINTER TO TYPE) */
    } DBEST_TAG_TYPE;

#define DBEST_TAG_TYPE_TOTAL 6

#define DBEST_TAG_WORD  0
#define DBEST_TAG_LINE  1
#define DBEST_TAG_BLOCK 2
#define DBEST_TAG_PARA  3
#define DBEST_TAG_INT   4
#define DBEST_TAG_FLOAT 5

static DBEST_TAG_TYPE local_tag_type [DBEST_TAG_TYPE_TOTAL] = {
{ parseDBESTword,   printDBESTstring, free, sizeof(char*)  },
{ parseDBESTline,   printDBESTstring, free, sizeof(char*)  },
{ parseDBESTblock,  printDBESTblock,  free, sizeof(char*)  },
{ parseDBESTpara,   printDBESTpara,   free, sizeof(char*)  },
{ parseDBESTint,    printDBESTint,    NULL, sizeof(int)    },
{ parseDBESTdouble, printDBESTdouble, NULL, sizeof(double) }
    };

/* ############################################################### */
/* SECTION INSTRUCTION DEFINITIONS                                 */

typedef struct {
       char *tag;      /* TAG TO INSERT INTO FSM                   */
        int  offset;   /* STRUCTURE MEMBER OFFSET                  */
     size_t  size;     /* SIZEOF OBJECT IF TO BE MADE BLANK        */
        int  start;    /* START OF MEMBER LIST IN local_inst       */
        int  stop;     /*   END OF MEMBER LIST IN local_inst       */
        int  multiple; /* COUNTER OFFSET IF MULTIPLE, -1 IF SINGLE */
    } PARSE_DBEST_SECTION_INST;

#define DBEST_SECTIONS 13

static PARSE_DBEST_SECTION_INST local_section_inst[DBEST_SECTIONS] = {
/* 00 */ {"\nIDENTIFIERS",    offsetof(DBESTRECORD, identifiers),
          sizeof(DBEST_IDENTIFIERS),     0,  7, -1},
/* 01 */ {"\nCLONE INFO",     offsetof(DBESTRECORD, cloneinfo),
          sizeof(DBEST_CLONE_INFO),      8, 15, -1},
/* 02 */ {"\nPRIMERS",        offsetof(DBESTRECORD, primers),
          sizeof(DBEST_PRIMERS),        16, 18, -1},
/* 03 */ {"\nSEQUENCE",       offsetof(DBESTRECORD, sequence),
          sizeof(DBEST_SEQUENCE),       19, 22, 
                              offsetof(DBESTRECORD, sequence_count)},
/* 04 */ {"\nCOMMENTS",       offsetof(DBESTRECORD, comments),
          sizeof(DBEST_COMMENTS),       23, 23, -1},
/* 05 */ {"\nPUTATIVE ID",    offsetof(DBESTRECORD, putativeid),
          sizeof(DBEST_PUTATIVE_ID),    24, 24, -1},
/* 06 */ {"\nLIBRARY",        offsetof(DBESTRECORD, library),
          sizeof(DBEST_LIBRARY),        25, 40, -1},
/* 07 */ {"\nSUBMITTER",      offsetof(DBESTRECORD, submitter),
          sizeof(DBEST_SUBMITTER),      41, 47, -1},
/* 08 */ {"\nCITATIONS",      offsetof(DBESTRECORD, citations),
          sizeof(DBEST_CITATIONS),      48, 53,
                              offsetof(DBESTRECORD, citations_count)},
/* 09 */ {"\nMAP DATA",       offsetof(DBESTRECORD, mapdata),
          sizeof(DBEST_MAP_DATA),       54, 57,
                              offsetof(DBESTRECORD, mapdata_count)},
/* 10 */ {"\nMAP METHODS",    offsetof(DBESTRECORD, mapmethods),
          sizeof(DBEST_MAP_METHODS),    58, 59,
                              offsetof(DBESTRECORD, mapmethods_count)},
/* 11 */ {"\nMAP SUBMITTERS", offsetof(DBESTRECORD, mapsubmitters),
          sizeof(DBEST_MAP_SUBMITTERS), 60, 66,
                         offsetof(DBESTRECORD, mapsubmitters_count)},
/* 12 */ {"\nNEIGHBORS",      offsetof(DBESTRECORD, neighbours),
          sizeof(DBEST_NEIGHBOURS),     67, 68, 
                              offsetof(DBESTRECORD, neighbours_count)}
    };

/* ############################################################### */
/* TAG INSTRUCTION DEFINITIONS                                     */

typedef struct {
       char *tag;     /* TAG TO INSERT INTO FSM          */
        int  section; /* SECTION TO WHICH THIS BELONGS   */
        int  offset;  /* STRUCTURE MEMBER OFFSET         */
        int  type;    /* FIELD TYPE TO USE               */
    } PARSE_DBEST_INST_TAG;

/* if(tag == NULL), is the default section instruction.
*/

#define DBEST_INSTRUCTIONS 69

static PARSE_DBEST_INST_TAG local_inst [DBEST_INSTRUCTIONS] = {
/* IDENTIFIERS */
/* 00 */  {"\ndbEST Id:",        DBEST_SECTION_IDENTIFIERS,
          offsetof(DBEST_IDENTIFIERS, dbEST_Id),     DBEST_TAG_INT},
/* 01 */  {"\nEST name:",        DBEST_SECTION_IDENTIFIERS,
          offsetof(DBEST_IDENTIFIERS, EST_name),     DBEST_TAG_WORD},
/* 02 */  {"\nGenBank Acc:",     DBEST_SECTION_IDENTIFIERS,
          offsetof(DBEST_IDENTIFIERS, GenBank_Acc),  DBEST_TAG_WORD},
/* 03 */  {"\nGenBank gi:",      DBEST_SECTION_IDENTIFIERS,
          offsetof(DBEST_IDENTIFIERS, GenBank_gi),   DBEST_TAG_INT},
/* 04 */  {"\nGDB Id:",          DBEST_SECTION_IDENTIFIERS,
          offsetof(DBEST_IDENTIFIERS, GDB_Id),       DBEST_TAG_INT},
/* 05 */  {"\nGDB Dsegment:",    DBEST_SECTION_IDENTIFIERS,
          offsetof(DBEST_IDENTIFIERS, GDB_Dsegment), DBEST_TAG_WORD},
/* 06 */  {"\nDatabase:",        DBEST_SECTION_IDENTIFIERS,
          offsetof(DBEST_IDENTIFIERS, Database), DBEST_TAG_LINE},
/* 07 */  {"\nSecondary Acc:",   DBEST_SECTION_IDENTIFIERS,
          offsetof(DBEST_IDENTIFIERS, Secondary_Acc), DBEST_TAG_WORD},
          /* ONLY FOUND IN dbEST_Id:727058 */
/* CLONE_INFO */
/* 08 */  {"\nClone Id:",    DBEST_SECTION_CLONE_INFO,
          offsetof(DBEST_CLONE_INFO, Clone_Id), DBEST_TAG_LINE},
/* 09 */  {"\nSource:",      DBEST_SECTION_CLONE_INFO,
          offsetof(DBEST_CLONE_INFO, Source),   DBEST_TAG_PARA},
/* 10 */  {"\nId as DNA:",  DBEST_SECTION_CLONE_INFO,
          offsetof(DBEST_CLONE_INFO, Id_as_DNA), DBEST_TAG_WORD},
/* 11 */  {"\nId in host:",  DBEST_SECTION_CLONE_INFO,
          offsetof(DBEST_CLONE_INFO, Id_in_host), DBEST_TAG_INT},
/* 12 */  {"\nOther ESTs on clone:",  DBEST_SECTION_CLONE_INFO,
          offsetof(DBEST_CLONE_INFO, Other_ESTs_on_clone), 
                                         DBEST_TAG_PARA},
/* 13 */  {"\nInsert length:",  DBEST_SECTION_CLONE_INFO,
          offsetof(DBEST_CLONE_INFO, Insert_length), DBEST_TAG_INT},
/* 14 */  {"\nPlate:",  DBEST_SECTION_CLONE_INFO,
          offsetof(DBEST_CLONE_INFO, Plate), DBEST_TAG_LINE},
/* 15 */  {"\nDNA type:",  DBEST_SECTION_CLONE_INFO,
          offsetof(DBEST_CLONE_INFO, DNA_type), DBEST_TAG_WORD},
/* PRIMERS */
/* 16 */  {"\nSequencing:",  DBEST_SECTION_PRIMERS,
          offsetof(DBEST_PRIMERS, Sequencing), DBEST_TAG_PARA},
/* 17 */  {"\nPCR forward:",  DBEST_SECTION_PRIMERS,
          offsetof(DBEST_PRIMERS, PCR_forward), DBEST_TAG_LINE},
/* 18 */  {"\nPCR backward:",  DBEST_SECTION_PRIMERS,
          offsetof(DBEST_PRIMERS, PCR_backward), DBEST_TAG_LINE},
/* SEQUENCE */
/* 19 */  { NULL,            DBEST_SECTION_SEQUENCE,
          offsetof(DBEST_SEQUENCE, sequence),      DBEST_TAG_BLOCK},
/* 20 */  {"\nQuality:", DBEST_SECTION_SEQUENCE,
          offsetof(DBEST_SEQUENCE, Quality), DBEST_TAG_LINE},
/* 21 */  {"\nEntry Created:", DBEST_SECTION_SEQUENCE,
          offsetof(DBEST_SEQUENCE, Entry_Created), DBEST_TAG_LINE},
/* 22 */  {"\nLast Updated:",   DBEST_SECTION_SEQUENCE,
          offsetof(DBEST_SEQUENCE, Last_Updated),  DBEST_TAG_LINE},
/* COMMENTS */
/* 23 */  {NULL,   DBEST_SECTION_COMMENTS,
          offsetof(DBEST_COMMENTS, comments),  DBEST_TAG_PARA},
/* PUTATIVE ID */
/* 24 */  {NULL,   DBEST_SECTION_PUTATIVE_ID,
          offsetof(DBEST_PUTATIVE_ID, putative_id),  DBEST_TAG_PARA},
/* LIBRARY */
/* 25 */  {"\ndbEST lib id:",   DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, dbEST_lib_id),  DBEST_TAG_INT},
/* 26 */  {"\nLib Name:",       DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, Lib_Name),  DBEST_TAG_PARA},
/* 27 */  {"\nOrganism:",       DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, Organism),  DBEST_TAG_LINE},
/* 28 */  {"\nStrain:",         DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, Strain),  DBEST_TAG_LINE},
/* 29 */  {"\nCultivar:",       DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, Cultivar), DBEST_TAG_LINE},
/* 30 */  {"\nSex:",            DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, Sex), DBEST_TAG_LINE},
/* 31 */  {"\nOrgan:",          DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, Organ),  DBEST_TAG_LINE},
/* 32 */  {"\nTissue type:",    DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, Tissue_type),  DBEST_TAG_PARA},
/* 33 */  {"\nCell type:", DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, Cell_type),  DBEST_TAG_LINE},
/* 34 */  {"\nCell line:", DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, Cell_line),  DBEST_TAG_WORD},
/* 35 */  {"\nDevelop. stage:", DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, Develop_stage),  DBEST_TAG_PARA},
/* 36 */  {"\nLab host:", DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, Lab_host),  DBEST_TAG_PARA},
/* 37 */  {"\nVector:",         DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, Vector),  DBEST_TAG_LINE},
/* 38 */  {"\nR. Site 1:",      DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, R_Site_1),  DBEST_TAG_WORD},
/* 39 */  {"\nR. Site 2:",      DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, R_Site_2),  DBEST_TAG_WORD},
/* 40 */  {"\nDescription:",    DBEST_SECTION_LIBRARY,
          offsetof(DBEST_LIBRARY, Description), DBEST_TAG_PARA},
/* SUBMITTER */
/* 41 */  {"\nName:",  DBEST_SECTION_SUBMITTER,
          offsetof(DBEST_SUBMITTER, Name),  DBEST_TAG_PARA},
/* 42 */  {"\nLab:",  DBEST_SECTION_SUBMITTER,
          offsetof(DBEST_SUBMITTER, Lab),  DBEST_TAG_PARA},
/* 43 */  {"\nInstitution:",  DBEST_SECTION_SUBMITTER,
          offsetof(DBEST_SUBMITTER, Institution),  DBEST_TAG_PARA},
/* 44 */  {"\nAddress:",  DBEST_SECTION_SUBMITTER,
          offsetof(DBEST_SUBMITTER, Address),  DBEST_TAG_PARA},
/* 45 */  {"\nTel:",  DBEST_SECTION_SUBMITTER,
          offsetof(DBEST_SUBMITTER, Tel),  DBEST_TAG_LINE},
/* 46 */  {"\nFax:",  DBEST_SECTION_SUBMITTER,
          offsetof(DBEST_SUBMITTER, Fax),  DBEST_TAG_LINE},
/* 47 */  {"\nE-mail:",  DBEST_SECTION_SUBMITTER,
          offsetof(DBEST_SUBMITTER, E_mail),  DBEST_TAG_LINE},
/* CITATIONS */
/* 48 */  {"\nMedline UID:",  DBEST_SECTION_CITATIONS,
          offsetof(DBEST_CITATIONS, Medline_UID),  DBEST_TAG_INT},
/* 49 */  {"\nTitle:",        DBEST_SECTION_CITATIONS,
          offsetof(DBEST_CITATIONS, Title),        DBEST_TAG_PARA},
/* 50 */  {"\nAuthors:",      DBEST_SECTION_CITATIONS,
          offsetof(DBEST_CITATIONS, Authors),      DBEST_TAG_PARA},
/* 51 */  {"\nCitation:",     DBEST_SECTION_CITATIONS,
          offsetof(DBEST_CITATIONS, Citation),     DBEST_TAG_LINE},
/* 52 */  {"\nYear:",         DBEST_SECTION_CITATIONS,
          offsetof(DBEST_CITATIONS, Year),     DBEST_TAG_INT},
/* 53 */  {"\nStatus:",       DBEST_SECTION_CITATIONS,
          offsetof(DBEST_CITATIONS, Status),     DBEST_TAG_LINE},
/* MAP DATA */
/* 54 */  {"\nMap:",          DBEST_SECTION_MAP_DATA,
          offsetof(DBEST_MAP_DATA, Map),       DBEST_TAG_LINE},
/* 55 */  {"\nMethod:",       DBEST_SECTION_MAP_DATA,
          offsetof(DBEST_MAP_DATA, Method),    DBEST_TAG_LINE},
/* 56 */  {"\nCitation:",     DBEST_SECTION_MAP_DATA,
          offsetof(DBEST_MAP_DATA, Citation),  DBEST_TAG_LINE},
/* 57 */  {"\nSubmitter:",     DBEST_SECTION_MAP_DATA,
          offsetof(DBEST_MAP_DATA, Submitter), DBEST_TAG_LINE},
/* MAP METHODS */
/* 58 */  {"\nMethod Name:",     DBEST_SECTION_MAP_METHODS,
          offsetof(DBEST_MAP_METHODS, Method_Name), DBEST_TAG_LINE},
/* 59 */  {"\nDescrn:",     DBEST_SECTION_MAP_METHODS,
          offsetof(DBEST_MAP_METHODS, Descrn), DBEST_TAG_PARA},
/* MAP SUBMITTERS */
/* 60 */  {"\nName:",  DBEST_SECTION_MAP_SUBMITTERS,
          offsetof(DBEST_MAP_SUBMITTERS, Name),  DBEST_TAG_LINE},
/* 61 */  {"\nLab:",  DBEST_SECTION_MAP_SUBMITTERS,
          offsetof(DBEST_MAP_SUBMITTERS, Lab),  DBEST_TAG_LINE},
/* 62 */  {"\nInstitution:",  DBEST_SECTION_MAP_SUBMITTERS,
          offsetof(DBEST_MAP_SUBMITTERS, Institution),  DBEST_TAG_PARA},
/* 63 */  {"\nAddress:",  DBEST_SECTION_MAP_SUBMITTERS,
          offsetof(DBEST_MAP_SUBMITTERS, Address),  DBEST_TAG_PARA},
/* 64 */  {"\nTel:",  DBEST_SECTION_MAP_SUBMITTERS,
          offsetof(DBEST_MAP_SUBMITTERS, Tel),  DBEST_TAG_LINE},
/* 65 */  {"\nFax:",  DBEST_SECTION_MAP_SUBMITTERS,
          offsetof(DBEST_MAP_SUBMITTERS, Fax),  DBEST_TAG_LINE},
/* 66 */  {"\nE-mail:",  DBEST_SECTION_MAP_SUBMITTERS,
          offsetof(DBEST_MAP_SUBMITTERS, E_mail),  DBEST_TAG_LINE},
/* NEIGHBOURS */
/* 67 */  {"\nNeighbor:",  DBEST_SECTION_NEIGHBOURS,
          offsetof(DBEST_NEIGHBOURS, Neighbor),  DBEST_TAG_PARA},
/* 68 */  {"\nPvalue:",  DBEST_SECTION_NEIGHBOURS,
          offsetof(DBEST_NEIGHBOURS, Pvalue),  DBEST_TAG_FLOAT}
/* ## */
    };

/* ############################################################### */
/* DBESTRECORD NEW/FREE/PRINT FUNCTIONS                            */

static DBESTRECORD *newDBESTRECORD(DBESTPARSER *dp){
    register DBESTRECORD *r = BLANK(DBESTRECORD);
    register int prev, mask = dp->interest, type; 
    register void *v;
    do {
        prev = mask;
        mask &= (mask-1); /* ONE ITERATION PER SECTION */
        type = ffs(prev^mask)-2;
        v = calloc(1, local_section_inst[type].size);
        OFFSET_ITEM(void*, local_section_inst[type].offset, r) = v;
    } while(mask);
    r->position = ftell(dp->fp);
    r->path = dp->paths[dp->path_count];
    return r;
    }

static void freeDBESTRECORDsection(void *data, 
                PARSE_DBEST_SECTION_INST *inst, int count){
    register int i;
    register void *v;
    register FREEFUNC ff;
    register int mctr, mstop;
    for(mctr = 0, mstop = count*inst->size;  /* MULTIPLE COUNTING */
                        mctr < mstop; mctr+=inst->size){
        for(i = inst->start; i <= inst->stop; i++)
            if((v 
                = OFFSET_ITEM(void*, local_inst[i].offset+mctr, data)))
                if((ff = local_tag_type[local_inst[i].type].free))
                    ff(v);
        }
    free(data);
    return;
    }

void freeDBESTRECORD(DBESTRECORD *r, DBESTPARSER *dp){
    register int prev, mask = dp->interest, type; 
    register void *v;
    do { 
        prev = mask;
        mask &= (mask-1); /* ONE ITERATION PER SECTION */
        type = ffs(prev^mask)-2;
        v = OFFSET_ITEM(void*, local_section_inst[type].offset, r);
        freeDBESTRECORDsection(v, &local_section_inst[type],
              (local_section_inst[type].multiple == -1)?1:
              OFFSET_ITEM(int, local_section_inst[type].multiple, r));
    } while(mask);
    free(r);
    return; 
    }

#ifdef DBEST_PRINT_MISSED_FIELDS
static int printDBESTRECORDsection(FILE *fp, void *data, 
                 PARSE_DBEST_SECTION_INST *inst, int count){
    register int i, total = 0;
    register int mctr, mstop;
    register void *v;
    for(i = inst->start; i <= inst->stop; i++)
        if((v = OFFSET_ITEM(void*, local_inst[i].offset, data))){
            total = fprintf(fp, "%s\n", inst->tag+1);
            break; /* ONLY PRINT SECTION NAME IF CONTAINS DATA */
            }
    if(total == 0)
        return 0;
    for(mctr = 0, mstop = count*inst->size;  /* MULTIPLE COUNTING */
                        mctr < mstop; mctr+=inst->size){
        for(i = inst->start; i <= inst->stop; i++){
            if((v = 
               OFFSET_ITEM(void*, local_inst[i].offset+mctr, data))){
                total += fprintf(fp, "%-15s ", local_inst[i].tag?
                                           local_inst[i].tag+1:" ");
                total += 
                      local_tag_type[local_inst[i].type].print(fp, v);
                }
            }
        total += fprintf(fp, "\n");
        }
    return total;
    }

static int printDBESTRECORD(FILE *fp, DBESTRECORD *r, DBESTPARSER *dp){
    register int prev, mask = dp->interest, type, total = 0; 
    register void *v;
    do {
        prev = mask;
        mask &= (mask-1); /* ONE ITERATION PER SECTION */
        type = ffs(prev^mask)-2;
        v = OFFSET_ITEM(void*, local_section_inst[type].offset, r);
        total += printDBESTRECORDsection(fp, v,
              &local_section_inst[type], 
              (local_section_inst[type].multiple == -1)?1:
              OFFSET_ITEM(int, local_section_inst[type].multiple, r));
    } while(mask);
    printf("||\n");
    return total;
    }
#endif /* DBEST_PRINT_MISSED_FIELDS */

/* ############################################################### */

typedef struct {
            /* s(ection)|e(nd)|p(arse)|c(ombination)|n(eighbours) */
    char  action;  
     int  section; /* IF action=='s', SECTION BEING ENTERED       */
     int *id;      /* IF action=='p', IS ARRAY OF DBEST_SECTIONS  */
    } PARSE_DBEST_INST;

static void freePARSE_DBEST_INST(void *v){
    register PARSE_DBEST_INST *pdi = (PARSE_DBEST_INST*)v;
    if(pdi->id)
        free(pdi->id);
    free(pdi);
    return;
    }

static void *joinPARSE_DBEST_INST(void *a, void *b){
    register PARSE_DBEST_INST *ia = a, *ib = b;
    register int i;
    if((ia->action != 'p') || (ib->action != 'p'))
        errmsg(ERROR_FATAL, "Non parseable redundant tags"); 
    if((ia->action == 'n') && (ib->action == 'n'))
        errmsg(ERROR_FATAL, "Cannot combine neighbour tags"); 
    for(i = 0; i < DBEST_SECTIONS; i++){
        if((ia->id[i] != -1) && (ib->id[i] != -1))
            errmsg(ERROR_FATAL, "Redundant tags from same section");
        if(ib->id[i] != -1)
            ia->id[i] = ib->id[i];
        }
    freePARSE_DBEST_INST(ib);
    return a;
    }

/* ############################################################### */
/* DBESTPARSER NEW/FREE FUNCTIONS                                  */
    
static BOOLEAN nextDBESTPARSERfile(DBESTPARSER *dp){
    if(dp->fp)
        fclose(dp->fp);
    if(dp->path_count >= dp->path_total)
        return FALSE; /* NO MORE FILES TO PARSE */
    errmsg(ERROR_INFO, "Parsing dbEST file [%s]",
                        dp->paths[dp->path_count]);
    dp->fp = fopen(dp->paths[dp->path_count++], "r");
    if(!dp->fp)
        errmsg(ERROR_FATAL, "Could not reopen dbEST file, [%s]",
                             dp->paths[dp->path_count-1]);
    return TRUE;
    }

static void initDBESTPARSERfile(DBESTPARSER *dp, 
                                       char **paths, int total){
    register int i;
    register FILE *fp;
    if(!total)
        errmsg(ERROR_FATAL, "No dbEST files to parse");
    dp->path_total = total;
    dp->paths = malloc(sizeof(char*)*total);
    for(i = 0; i < total; i++){
        dp->paths[i] = strdup(paths[i]);
        /* TRY OPENING EACH FILE AT START TO CHECK OK */
        fp = fopen(paths[i], "r");
        if(!fp)
            errmsg(ERROR_FATAL, "Could not open dbEST file, [%s]",
                   paths[i]);
        fclose(fp);
        }
    errmsg(ERROR_INFO, "All dbEST file names valid");
    nextDBESTPARSERfile(dp);
    return;
    }
    
static void freeDBESTPARSERfile(DBESTPARSER *dp){
    register int i;
    for(i = 0; i < dp->path_total; i++)
        free(dp->paths[i]);
    free(dp->paths);
    return;
    }

DBESTPARSER *newDBESTPARSER(char **paths, int total, 
                                   int interest){
    register DBESTPARSER *dp = BLANK(DBESTPARSER);
    register LIST *tag = newLIST(), *inst = newLIST();
    register int i = 0, j = 0, inst_total = 0;
    int notused;
    register PARSE_DBEST_INST *pdi;
    register void **namelist, **datalist;
    initDBESTPARSERfile(dp, paths, total);
    if(!interest)
        errmsg(ERROR_FATAL, "No dbEST sections to parse");
    dp->count = 0;
    dp->interest = interest;

    pdi = BLANK(PARSE_DBEST_INST); /* END OF RECORD TAG */
    pdi->action = (interest & DBEST_SECTION_NEIGHBOURS)?'E':'e';
    queueLIST(tag, "\n||");
    queueLIST(inst, pdi);

    for(i = 0; i < DBEST_SECTIONS; i++){ /* PARSE SECTIONS */
        pdi = BLANK(PARSE_DBEST_INST);
        /* IF INTERESTED AND HAS SECTION TAG */
        if((local_inst[local_section_inst[i].start].section & interest) 
         && (!local_inst[local_section_inst[i].start].tag) )
            pdi->action = 'c';
        else
            pdi->action = 's';
        pdi->section = i;
        queueLIST(tag, local_section_inst[i].tag);
        queueLIST(inst, pdi);
        }

    inst_total = 1+DBEST_SECTIONS;

    for(i = 0; i < DBEST_INSTRUCTIONS; i++){ /* PARSE TAGS */
        if(local_inst[i].section & interest){ /* IF OF INTEREST */
            if(local_inst[i].tag){ /* IF NOT COMBINATION TAG */
                pdi = BLANK(PARSE_DBEST_INST);
                pdi->action = 'p';
                pdi->id = malloc(sizeof(int)*DBEST_SECTIONS);
                for(j = 0; j < DBEST_SECTIONS; j++)
                    pdi->id[j] = -1;
                pdi->id[ffs(local_inst[i].section)-2] = i;
                queueLIST(tag, local_inst[i].tag);
                queueLIST(inst, pdi);
                inst_total++;
                }
            }
        }

    /* INSTRUCTION FOR NEIGHBOURS */
    if(interest & DBEST_SECTION_NEIGHBOURS){
        pdi = BLANK(PARSE_DBEST_INST);
        pdi->action = 'n';
        pdi->section = TRUE;
        queueLIST(tag, "\nTop 15 nucleotide matches");
        queueLIST(inst, pdi);
        pdi = BLANK(PARSE_DBEST_INST);
        pdi->action = 'n';
        pdi->section = FALSE;
        queueLIST(tag, "\nTop 15 protein matches");
        queueLIST(inst, pdi);
        inst_total+=2;
        }

    namelist = toarrayLIST(tag,  &notused);
    datalist = toarrayLIST(inst, &notused);
    dp->fsm = makeFSMptrJOIN((char**)namelist,
                   inst_total, datalist, joinPARSE_DBEST_INST);
    free(namelist);
    free(datalist);

    errmsg(ERROR_INFO,
           "Built dbest parser to find %d field types", inst_total);
    return dp;
    }

void freeDBESTPARSER(DBESTPARSER *dp){
    freeDBESTPARSERfile(dp);
    freeFSMptr(dp->fsm, freePARSE_DBEST_INST);
    free(dp);
    return;
    }

/* ############################################################### */
/* CENTRAL PARSING FUNCTIONS                                       */

static FSMNODE *parseField(DBESTPARSER *dp, DBESTRECORD *r,
                           int id, int section){
    register void *v, **ptp;
    register PARSEFUNC pf = local_tag_type[local_inst[id].type].parse;
    register int *multiple, offset = 0;
    register size_t oldsize, addsize;
    register double *dptr;
    ANYTYPE a;
    v = OFFSET_ITEM(void*, local_section_inst[section].offset, r);
    if(local_section_inst[section].multiple == -1){ /* SINGLE TYPE */
        ptp = &OFFSET_ITEM(void*, local_inst[id].offset, v);
        if(*ptp)
            errmsg(ERROR_FATAL, "Unexpected 2nd entry for [%s:%s]", 
                          local_section_inst[section].tag+1,
                          local_inst[id].tag+1);
    } else { /* MULTIPLE TYPE */
        multiple = &OFFSET_ITEM(int,
                    local_section_inst[section].multiple, r);
        offset = local_inst[id].offset
               + (local_section_inst[section].size*(*multiple));
        ptp = &OFFSET_ITEM(void*, offset, v);
        if(*ptp){ /* IF ENTRY ALREADY PRESENT */
            if(!((*multiple) & (PARSEDBEST_MULTIPLE_CHUNKSIZE-1))){
                oldsize = local_section_inst[section].size
                                *((*multiple)+1);
                addsize = local_section_inst[section].size
                                *(PARSEDBEST_MULTIPLE_CHUNKSIZE);
                v = realloc(v, oldsize+addsize);
                memset((char*)v+oldsize, 0, addsize);
                OFFSET_ITEM(void*,
                   local_section_inst[section].offset, r) = v;
                }
            (*multiple)++;
            offset += local_section_inst[section].size;
            ptp = &OFFSET_ITEM(void*, offset, v);
            }
        }
    a.u.v = *ptp;   /* ALLOW MODIFIABILILTY */
    pf(dp->fp, &a); /* CALL PARSE FUNCTION  */
    if(a.type == 'd'){
        dptr = &OFFSET_ITEM(double, offset, v);
        *dptr = a.u.d;
    } else {
        *ptp = a.u.v; /* SEEMS TO BE OK FOR OTHER TYPES (CHANGE) */
        }
    return dp->fsm->root[dp->fsm->index['\n']].next; /* RESET FSM */
    }

/*
     v IS ADDRESS OF SECTION STRUCTURE IN r.
   ptp IS POINTER TO ADDRESS OF ITEM IN v.
*/

static void setDBESTRECORDmultiples(DBESTRECORD *r){
    register int i, *multiple;
    for(i = 0; i < DBEST_SECTIONS; i++)
       if(local_section_inst[i].multiple != -1){
           multiple = &OFFSET_ITEM(int, 
                    local_section_inst[i].multiple, r);
           (*multiple)++;
           }
    return;
    }

/* ############################################################### */

DBESTRECORD *parseDBESTRECORD(DBESTPARSER *dp){
    register int section = DBEST_SECTION_UNKNOWN, ch = '\n', id;
    register unsigned char c;
    register FSMNODE *n = dp->fsm->root[dp->fsm->index['\n']].next;
    register PARSE_DBEST_INST *inst;
    register DBESTRECORD *r = newDBESTRECORD(dp);
#ifdef DBEST_PRINT_MISSED_FIELDS
    register BOOLEAN prevspace = TRUE;
    register BOOLEAN shownid = FALSE;
#endif /* DBEST_PRINT_MISSED_FIELDS */
    while((ch = getc(dp->fp)) != EOF){
        c = dp->fsm->index[ch];
        inst = n[c].data.n;
        n = n[c].next;
        if(inst){
            switch(inst->action){
                case 's': /* SECTION */
                    section = inst->section;
                    break;
                case 'E': /* END OF RECORD (WITH NEIGHBOURS) */
                    if(r->neighbours_protein == -1)
                        r->neighbours_protein =
                            r->neighbours_count+1; 
 /*FALLTHROUGH*/case 'e': /* END OF RECORD */
                    section = DBEST_SECTION_UNKNOWN;
                    dp->count++;
                    setDBESTRECORDmultiples(r);
                    return r;
                case 'c': /* COMBINATION */
                    section = inst->section;
                    id = local_section_inst[section].start;
                    n = parseField(dp, r, id, section);
                    break;
                case 'p': /* PARSE FIELD */
                    id = inst->id[section];
                    if(id != -1) /* NEVER HERE ??? */
                        n = parseField(dp, r, id, section); 
                    break;
                case 'n': /* NEIGHBOURS */
                    if(inst->section){ /* ASSUME PROTEIN 1ST */
                        if(r->neighbours_protein == -1)
                            r->neighbours_protein = 
                               r->neighbours_count+1; 
                    } else { /* NUCLEOTIDE SIGN */
                        r->neighbours_protein = -1;
                        }
                    break;
                default:
                    errmsg(ERROR_FATAL, "Unknown action [%c]",
                                         inst->action);
                    break;
                }
            }
#ifdef DBEST_PRINT_MISSED_FIELDS
/* DEBUGGING CODE TO CHECK THAT PARSER UNDERSTANDS ALL OF DBEST */
              else {
                if(n == dp->fsm->root[dp->fsm->index[ch]].next)
                    if(isspace(ch)){
                        if(prevspace == FALSE)
                            printf(" ");
                        prevspace = TRUE;
                    } else {
                        if(!shownid){
                            printf("\nID[%d]\n", 
                                r->identifiers->dbEST_Id);
                            shownid = TRUE;
                            }
                        printf("%c", ch);
                        prevspace = FALSE;
                        }
                }
#endif /* DBEST_PRINT_MISSED_FIELDS */
        }
    freeDBESTRECORD(r, dp);
    if(nextDBESTPARSERfile(dp))
        return parseDBESTRECORD(dp); /* RECURSE */
    else
        return NULL;
    }

/* ############################################################### */

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

#ifdef DBEST_MONITOR_USAGE

/* monitorDBESTusage: SHOWS THE AMMOUNT EACH FIELD IS USED.
*/
static void monitorDBESTusage(char **paths, int total){
    register DBESTPARSER *dp = newDBESTPARSER(paths, total, 
                                  DBEST_SECTION_ALL);
    register DBESTRECORD *r;
    register int i, j, section;
    int count[DBEST_INSTRUCTIONS], mean[DBEST_INSTRUCTIONS];
    register void *sect, *inst;
    for(i = 0; i < DBEST_INSTRUCTIONS; i++) 
        count[i] = mean[i] = 0;
    while((r = parseDBESTRECORD(dp))){
        for(i = 0; i < DBEST_INSTRUCTIONS; i++){
            section = ffs(local_inst[i].section)-2;
            sect = OFFSET_ITEM(void*, 
                   local_section_inst[section].offset, r);
            inst = OFFSET_ITEM(void*, local_inst[i].offset, sect);
            if(inst)
                count[i]++;
/*
            if(local_inst[section].multiple != -1){
                multiple = OFFSET_ITEM(int, 
                    local_section_inst[section].multiple, r);
                for(j = 0; j < multiple; j++){
                    }
                }
*/
            }
        freeDBESTRECORD(r, dp);
        }
    printf("Results:\n");
    for(i = 0; i < DBEST_SECTIONS; i++){
        printf("%s", local_section_inst[i].tag+1);
        if(local_section_inst[i].multiple != -1)
            printf(" [multiple]");
        printf("\n");
        for(j = local_section_inst[i].start; 
            j <= local_section_inst[i].stop; j++){
            if(local_inst[j].tag)
                printf("%s ", local_inst[j].tag+1);
            else 
                printf("    ");
            printf(" %2.2f%%   (%d)\n", 
                  ((float)count[j]/(float)dp->count)*100.0, count[j]);
            }
        } 
    errmsg(ERROR_INFO, "Monitored usage in [%d] files", dp->count);
    freeDBESTPARSER(dp);
    return;
    }
#endif /* DBEST_MONITOR_USAGE */

/* ############################################################### */

#ifdef DBEST_SHOW_FORMAT
/* printDBESTformat : PRINT OUT WHOLE DATABASE FORMAT
                      AS UNDERSTOOD BY THIS PARSER
*/
static int printDBESTformat(){
    register int i, j;
    for(i = 0; i < DBEST_SECTIONS; i++){
        printf("%s", local_section_inst[i].tag+1);
        if(local_section_inst[i].multiple != -1)
            printf(" [multiple]");
        printf("\n");
        for(j = local_section_inst[i].start; 
            j <= local_section_inst[i].stop; j++){
            if(local_inst[j].tag)
                printf("%s ", local_inst[j].tag+1);
            else 
                printf("    ");
            printf("<");
            switch(local_inst[j].type){
                case DBEST_TAG_WORD:
                    printf("word");
                    break;
                case DBEST_TAG_LINE:
                    printf("line");
                    break;
                case DBEST_TAG_BLOCK:
                    printf("block");
                    break;
                case DBEST_TAG_PARA:
                    printf("paragraph");
                    break;
                case DBEST_TAG_INT:
                    printf("int");
                    break;
                case DBEST_TAG_FLOAT:
                    printf("float");
                    break;
                default:
                    errmsg(ERROR_FATAL, "Unknown tag type");
                    break;
                }
            printf(">\n");
            }
        printf("\n");
        }
    return;
    }
#endif /* DBEST_SHOW_FORMAT */

/* ############################################################### */

int main(int argc, char **argv){
#ifdef DBEST_SHOW_FORMAT
    printDBESTformat();
#endif /* DBEST_SHOW_FORMAT */

#ifdef DBEST_MONITOR_USAGE
    monitorDBESTusage(argv+1, argc-1);
#endif /* DBEST_MONITOR_USAGE */

#ifdef DBEST_TEST_PARSER
    register DBESTPARSER *dp = newDBESTPARSER(argv+1, argc-1,
                                 0 
                              |DBEST_SECTION_IDENTIFIERS
                              /* |DBEST_SECTION_CLONE_INFO */
                              /* |DBEST_SECTION_PRIMERS */
                              /* |DBEST_SECTION_SEQUENCE */
                              /* |DBEST_SECTION_COMMENTS */
                              /* |DBEST_SECTION_PUTATIVE_ID */
                              /* |DBEST_SECTION_LIBRARY */
                              /* |DBEST_SECTION_SUBMITTER */
                              /* |DBEST_SECTION_CITATIONS  */
                              /* |DBEST_SECTION_MAP_DATA */
                              /* |DBEST_SECTION_MAP_METHODS */
                              /* |DBEST_SECTION_MAP_SUBMITTERS */
                              /* |DBEST_SECTION_NEIGHBOURS */
                              |DBEST_SECTION_ALL
                              );
    register DBESTRECORD *r;
    while(r = parseDBESTRECORD(dp)){
#ifdef DBEST_PRINT_MISSED_FIELDS
        printf("\n### ... parsed this record:\n");
        printDBESTRECORD(stdout, r, dp);
#endif /* DBEST_PRINT_MISSED_FIELDS */
        freeDBESTRECORD(r, dp);
        }
    printf("Parsed [%d] record%c\n", dp->count, (dp->count>1)?'s':' ');
    freeDBESTPARSER(dp);
#endif /* DBEST_TEST_PARSER */
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */

#endif /* INCLUDED_PARSEDBEST_C */

/*
To do:
    Add "^PolyA Tail:" field for PRIMERS section.
    Calculate number of different entries (if<100) for each field.
    Write script reformat interpreter ?
*/

