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

#ifndef INCLUDED_PARSEDBEST_H
#define INCLUDED_PARSEDBEST_H

#include <stdio.h>
#include "../general/common.h"
#include "../struct/fsm.h"

/* SECTION MASK COMPONENTS */

#define DBEST_SECTION_UNKNOWN        (1<< 0)
#define DBEST_SECTION_IDENTIFIERS    (1<< 1)
#define DBEST_SECTION_CLONE_INFO     (1<< 2)
#define DBEST_SECTION_PRIMERS        (1<< 3)
#define DBEST_SECTION_SEQUENCE       (1<< 4)
#define DBEST_SECTION_COMMENTS       (1<< 5)
#define DBEST_SECTION_PUTATIVE_ID    (1<< 6)
#define DBEST_SECTION_LIBRARY        (1<< 7)
#define DBEST_SECTION_SUBMITTER      (1<< 8)
#define DBEST_SECTION_CITATIONS      (1<< 9)
#define DBEST_SECTION_MAP_DATA       (1<<10)
#define DBEST_SECTION_MAP_METHODS    (1<<11)
#define DBEST_SECTION_MAP_SUBMITTERS (1<<12)
#define DBEST_SECTION_NEIGHBOURS     (1<<13)
#define DBEST_SECTION_SPECIAL        (1<<14)

#define DBEST_SECTION_ALL            \
        DBEST_SECTION_IDENTIFIERS    \
       |DBEST_SECTION_CLONE_INFO     \
       |DBEST_SECTION_PRIMERS        \
       |DBEST_SECTION_SEQUENCE       \
       |DBEST_SECTION_COMMENTS       \
       |DBEST_SECTION_PUTATIVE_ID    \
       |DBEST_SECTION_LIBRARY        \
       |DBEST_SECTION_SUBMITTER      \
       |DBEST_SECTION_CITATIONS      \
       |DBEST_SECTION_MAP_DATA       \
       |DBEST_SECTION_MAP_METHODS    \
       |DBEST_SECTION_MAP_SUBMITTERS \
       |DBEST_SECTION_NEIGHBOURS    

/* SECTION STRUCTURES */

typedef struct {
     int  dbEST_Id; 
    char *EST_name;
    char *GenBank_Acc;
     int  GenBank_gi; 
    char *GDB_Id;
    char *GDB_Dsegment;
    char *Database;
    char *Secondary_Acc;
    } DBEST_IDENTIFIERS;

typedef struct {
    char *Clone_Id;
    char *Source; 
    char *Id_as_DNA;
     int *Id_in_host;
    char *Other_ESTs_on_clone;
     int  Insert_length;
    char *Plate;
    char *DNA_type;
    } DBEST_CLONE_INFO;

typedef struct {
    char *Sequencing;
    char *PCR_forward;
    char *PCR_backward;
    } DBEST_PRIMERS;

typedef struct {
    char *sequence;
    char *Quality;
    char *Entry_Created;
    char *Last_Updated;
    } DBEST_SEQUENCE;

typedef struct {
    int   dbEST_lib_id;
   char  *Lib_Name;
   char  *Organism;
   char  *Strain;
   char  *Cultivar; 
   char  *Sex;
   char  *Organ;
   char  *Tissue_type;
   char  *Cell_type;
   char  *Cell_line;
   char  *Develop_stage;
   char  *Lab_host;
   char  *Vector; 
   char  *R_Site_1;
   char  *R_Site_2;
   char  *Description;
    } DBEST_LIBRARY;

typedef struct {
    char *Name;
    char *Lab;
    char *Institution;
    char *Address;
    char *Tel;
    char *Fax;
    char *E_mail;
    } DBEST_SUBMITTER;

typedef DBEST_SUBMITTER DBEST_MAP_SUBMITTERS;

typedef struct {
    int  Medline_UID;
   char *Title;
   char *Authors;
   char *Citation;
    int  Year;
   char *Status;
    } DBEST_CITATIONS;

typedef struct {
    char *Map;
    char *Method;
    char *Citation;
    char *Submitter;
    } DBEST_MAP_DATA;

typedef struct {
    char *Method_Name;
    char *Descrn;
    } DBEST_MAP_METHODS;


typedef struct {
    char *comments;
    } DBEST_COMMENTS;

typedef struct {
    char *putative_id;
    } DBEST_PUTATIVE_ID;

typedef struct {
    char *Neighbor;
  double  Pvalue;
    } DBEST_NEIGHBOURS;

/* WHOLE RECORD STRUCTURE */

typedef struct {
    char *path;       /* POINTER TO CURRENT FILE PATH     */
    long position;    /* POSITION OF RECORD START IN FILE */
    DBEST_IDENTIFIERS    *identifiers;
    DBEST_CLONE_INFO     *cloneinfo;
    DBEST_PRIMERS        *primers;
    DBEST_SEQUENCE       *sequence;
               int       *sequence_count; /* (MULTIPLE QUALITY) */
    DBEST_COMMENTS       *comments;
    DBEST_PUTATIVE_ID    *putativeid;
    DBEST_LIBRARY        *library;
    DBEST_SUBMITTER      *submitter;
    DBEST_CITATIONS      *citations;
                int       citations_count;
    DBEST_MAP_DATA       *mapdata;
                int       mapdata_count;
    DBEST_MAP_METHODS    *mapmethods;
                int       mapmethods_count;
    DBEST_MAP_SUBMITTERS *mapsubmitters;
                int       mapsubmitters_count;
    DBEST_NEIGHBOURS     *neighbours;
                int       neighbours_count;
                int       neighbours_protein; /* NUMBER PROTEIN HITS */
    } DBESTRECORD;

/* PARSER STRUCTURE */

typedef struct {
        FILE  *fp;         /* CURRENT FILE BEING PARSED  */
         FSM  *fsm;        /* FSM TO PARSE FILE          */
         int   interest;   /* SECTIONS TO BE PARSED      */
         int   count;      /* NUMBER OF RECORDS PARSED   */
        char **paths;      /* LIST OF FILES TO PARSE     */
         int   path_total; /* NUMBER OF FILES TO PARSE   */
         int   path_count; /* CURRENT FILE BEING PARSED  */
    } DBESTPARSER;

DBESTPARSER *  newDBESTPARSER(char **paths, int total, int interest);
       void   freeDBESTPARSER(DBESTPARSER *dp);
DBESTRECORD *parseDBESTRECORD(DBESTPARSER *dp);
       void   freeDBESTRECORD(DBESTRECORD *r, DBESTPARSER *dp);

#endif /* INCLUDED_PARSEDBEST_H */


