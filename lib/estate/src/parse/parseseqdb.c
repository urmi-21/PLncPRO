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

#ifndef INCLUDED_PARSESEQDB_C
#define INCLUDED_PARSESEQDB_C

#include <stdio.h>
#include <ctype.h>
#include <string.h> /* FOR strdup() */
#include <time.h> /* FOR strptime() */

#ifdef linux
/* LINUX SEEMS FUSSY ABOUT THIS WITH -Wall : DON'T KNOW WHY */
char *strptime(char *buf, const char *format, const struct tm *tm);
#endif /* linux */

#include "../general/common.h"
#include "parseseqdb.h"
#include "parser.h"
#include "ioutil.h"

#define DBLINEhash(a,b,c) ((a)|((b)<<8)|((c)<<16)) 
#define MAXDATELINE 144

static void freeDBREFERENCE(DBREFERENCE *r){
    int i;
    if(r->author)
        freeSTRING(r->author);
    if(r->title)
        freeSTRING(r->title);
    if(r->journal) 
        freeSTRING(r->journal);
    if(r->position) 
        freeSTRING(r->position);
    if(r->comment) 
        freeSTRING(r->comment);
    if(r->dbxrefv){
        for(i = 0; i < r->dbxrefc; i++)
            freeSTRINGLIST(r->dbxrefv[i]);
        free(r->dbxrefv); 
        }
    free(r);
    return;
    }

static void freeDBFEATURE(DBFEATURE *f){
    int i;
    if(f->feature)
        freeSTRING(f->feature);
    if(f->location)
        freeSTRING(f->location);
    for(i = 0; i < f->count; i++){
            if(f->name[i])
                freeSTRING(f->name[i]);
            if(f->value[i])
                freeSTRING(f->value[i]);
            }
    if(f->name)  
        free(f->name);
    if(f->value)  
        free(f->value);
    free(f);
    return;
    }

void freeDBENTRY(DBENTRY *e){
    int i;
    if(e->identifier)
        free(e->identifier); 
    if(e->definition)
        freeSTRING(e->definition);
    if(e->division) 
        free(e->division); 
    if(e->accession)
        freeSTRINGLIST(e->accession);
    if(e->keyword)
        freeSTRINGLIST(e->keyword);
    if(e->organism) 
        free(e->organism);
    if(e->classification) 
        freeSTRINGLIST(e->classification);
    if(e->comment)
        freeSTRING(e->comment);
    if(e->genename)
        freeSTRING(e->genename);
    if(e->nid)
        free(e->nid);
    if(e->referencev){
        for(i = 0; i < e->referencec; i++)
            freeDBREFERENCE(e->referencev[i]);
        free(e->referencev);
        }
    if(e->dbxrefv){
        for(i = 0; i < e->dbxrefc; i++)
            freeSTRINGLIST(e->dbxrefv[i]);
        free(e->dbxrefv); 
        }
    if(e->organelle)
        freeSTRING(e->organelle);
    if(e->basecount)
        free(e->basecount);
    if(e->sequence)
        freeSTRING(e->sequence);
    if(e->segment)
        free(e->segment);
    if(e->featuretablev){
        for(i = 0; i < e->featuretablec; i++)
            freeDBFEATURE(e->featuretablev[i]);
        free(e->featuretablev);
        }
    free(e);
    return;
    }

static int printDBFEATURE(FILE *fp, DBFEATURE *f, int num){
    int i, total = 0;
    total += fprintf(fp, "##FEATURE [%d]\n", num);
    if(f->feature)
        total += fprintf(fp, "  FEATURE [%s]\n", f->feature->str);
    if(f->location)
        total += fprintf(fp, "  LOCATION [%s]\n", f->location->str);
    for(i = 0; i < f->count; i++)
        if(f->name[i])
            total += fprintf(fp, "[%d] NAME [%s] VALUE [%s]\n",
                i, f->name[i]->str, f->value[i]?f->value[i]->str:"-");
    return total;
    }

static int printDBREFERENCE(FILE *fp, DBREFERENCE *r, int num){
    int i, total = 0;
    total += fprintf(fp, ">--REFERENCE [%d]\n", num);
    if(r->author)
        total += fprintf(fp, "  AUTHOR    [%s]\n", r->author->str);
    if(r->title)
        total += fprintf(fp, "  TITLE     [%s]\n", r->title->str);
    if(r->journal)
        total += fprintf(fp, "  JOURNAL   [%s]\n", r->journal->str);
    if(r->comment)
        total += fprintf(fp, "  COMMENT   [%s]\n", r->comment->str);
    if(r->position)
        total += fprintf(fp, "  POSITION  [%s]\n", r->position->str);
    for(i = 0; i < r->dbxrefc; i++){
        total += fprintf(fp, "  REFDBXREF [%d]\n", i);
        total += printSTRINGLIST(fp, r->dbxrefv[i]);
        }
    return total;
    }

int printDBENTRY(FILE *fp, DBENTRY *e){
    int i, total = 0;
    struct tm *t;
    char line[MAXDATELINE];
    total += fprintf(fp, "{ ---\n");
    if(e->identifier)
        total += fprintf(fp, "IDENTIFIER [%s]\n", e->identifier);
    if(e->definition)
        total += fprintf(fp, "DEFINITION [%s]\n", e->definition->str);
    if(e->length)
        total += fprintf(fp, "LENGTH     [%d]\n", e->length);
    if(e->date){
        t = localtime(&e->date);
        strftime(line, MAXDATELINE, "%A %d-%B-%Y", t);
        total += fprintf(fp, "DATE [%s]\n", line);
        }
    if(e->division)
        total += fprintf(fp, "DIVISION [%s]\n", e->division);
    if(e->type)
        total += fprintf(fp, "TYPE [%s]\n",
                   (e->type=='N')?"Nucleotide":"Peptide");
    if(e->accession){
        total += fprintf(fp, "ACCESSION NUMBERS:\n");
        total += printSTRINGLIST(fp, e->accession);
        }
    if(e->keyword){
        total += fprintf(fp,"KEYWORDS: \n");
        total += printSTRINGLIST(fp, e->keyword);
        }
    if(e->organism)
        total += fprintf(fp, "ORGANISM [%s]\n", e->organism);
    if(e->genename)
        total += fprintf(fp, "GENENAME [%s]\n", e->genename->str);
    if(e->classification){
        total += fprintf(fp,"CLASSIFICATION: \n");
        total += printSTRINGLIST(fp, e->classification);
        }
    if(e->comment)
        total += fprintf(fp, "COMMENT [%s]\n", e->comment->str);
    if(e->nid)
        total += fprintf(fp, "NID [%s]\n", e->nid);
    for(i = 0; i < e->referencec; i++)
        total += printDBREFERENCE(fp, e->referencev[i], i);
    for(i = 0; i < e->dbxrefc; i++){
        total += fprintf(fp, "DBXREF [%d]\n", i);
        total += printSTRINGLIST(fp, e->dbxrefv[i]);
        }
    if(e->organelle)
        total += fprintf(fp, "ORGANELLE [%s]\n", e->organelle->str);
    if(e->basecount)
        total += fprintf(fp,
                  "BASECOUNT [%d A][%d C][[%d G][%d T] [%d other\n",
                  e->basecount[0], e->basecount[1],
                  e->basecount[2], e->basecount[3], e->basecount[4]);
    if(e->sequence)
        total += fprintf(fp, "SEQUENCE [%s]\n", e->sequence->str);
    if(e->segment)
        total += fprintf(fp, "SEGMENT [%d] of [%d]\n",
                 e->segment[0], e->segment[1]);
    if(e->featuretablev)
        for(i = 0; i < e->featuretablec; i++)
            if(e->featuretablev[i])
                total += printDBFEATURE(fp, e->featuretablev[i], i);
    total += fprintf(fp, "--- }\n");
    return total;
    }

static time_t ConvertDBDate(char *date){
    struct tm t;
    strptime(date, "%d-%b-%Y ", &t);
    return mktime(&t);
    }

static void ParseEmblIdentifier(DBENTRY *e, char *line, int len){
    register STRINGLIST *s = toSTRINGLIST(line, len, ' ');
    register int i;
    e->identifier = strdup(s->v[0]);
    e->length = atoi(s->v[s->c-2]);
    if(s->c == 6){
        e->type = 'N';
        e->division = strdup(s->v[3]);
        i = strlen(e->division)-1;
        if(e->division[i] == ';')
           e->division[i] = '\0';
    } else {
        e->type = 'P';
    }
    freeSTRINGLIST(s);
    return;
    }

static void ParseGenBankLocus(DBENTRY *e, char *line, int len){
    STRINGLIST *s = toSTRINGLIST(line, len, ' ');
    e->identifier = strdup(s->v[0]);
    e->length = atoi(s->v[1]);
    e->division = strdup(s->v[s->c-2]);
    e->date = ConvertDBDate(s->v[s->c-1]);
    freeSTRINGLIST(s);
    e->type = 'N';  /* AS GENBANK ALWAYS PROTEIN */
    return;
    }

static STRINGLIST *ProcessDBLIST(char *str, int len){
    if(str[len-1] == '.'){
       if(len < 2)
           return NULL;
       else
           str[--len] = '\0';
       }
    return toSTRINGLIST(str, len, ';');
    }

static int CleanSequenceLine(char *str){
    char *s, *d;
    s = d = str;
    while((*d = toupper(*s))){
        if(isalpha(*d))
            d++;
        s++;
        }
    return d-str;
    }

static STRING *CleanSTRING(STRING *s){
    char *p, *d;
    if( (s->str[0] == '\0')
        && ((s->str[1] == '\0') && (s->str[1] == ';')) ){
        freeSTRING(s);
        return (STRING*)0;
        }
    switch(s->len){
        case 1:
            if(*s->str != ';')
                break;
        case 0:
            freeSTRING(s);
            return (STRING*)0;
        }
    if(s->str[s->len-1] == ';')
       s->str[s->len-1] = '\0';
    if(*s->str != '\"')
        return s;
    p = s->str+1;
    d = s->str;
    while((*d++ = *p++));
    s->len = d-s->str;
    if(d[-2] == '\"'){
        d[-2] = '\0';
        s->len-=2;
        }
    return s;
    }

static int *ParseSegment(char *line, int len){
    STRINGLIST *s = toSTRINGLIST(line, len, ' ');
    int *sg = malloc(sizeof(int)*2);
    if(s->c == 3){
        sg[0] = atoi(s->v[0]);
        sg[1] = atoi(s->v[2]);
        return sg;
        }
    free(sg);
    return (int*)0;
    }

static int *ParseBaseCount(char *line, int len){
    register STRINGLIST *s = toSTRINGLIST(line, len, ' ');
    register int *bc = calloc(5, sizeof(int));
    switch(s->c){
        case 10:  /* 0,2,4,6,8 */
            bc[4] = atoi(s->v[8]);
/*FALLTHROUGH*/
        case 8:   /* 0,2,4,6 */
            bc[0] = atoi(s->v[0]);
            bc[1] = atoi(s->v[2]);
            bc[2] = atoi(s->v[4]);
            bc[3] = atoi(s->v[6]);
            break;
        case 13:  /* 3,5,7,9,11 */
            bc[0] = atoi(s->v[ 3]);
            bc[1] = atoi(s->v[ 5]);
            bc[2] = atoi(s->v[ 7]);
            bc[3] = atoi(s->v[ 9]);
            bc[4] = atoi(s->v[11]);
            break;
        default:
            free(bc);
            freeSTRINGLIST(s);
            return (int*)0;
        }
    freeSTRINGLIST(s);
    return bc;
    }

void CleanFeatureTable(DBFEATURE *f){
    int i;
    register char *p, *s;
    for(i = 0; i < f->count; i++){ /* CLEAN "" FROM VALUE ENDS */
        if(f->value[i] && (f->value[i]->str[0] == '\"')){
            p = f->value[i]->str;
            s = p+1;
            while((*p++ = *s++));
            while(p > f->value[i]->str){
                  if(*p == '\"'){
                      *p = '\0';
                      break;
                      }
                  p--;
                  }
            }
        }
    return;
    }

DBENTRY *GetNextDatabaseEntry(FILE *fp, long mask){
    DBENTRY *e = BLANK(DBENTRY);
    DBREFERENCE *r = (DBREFERENCE*)0;
    DBFEATURE *f = (DBFEATURE*)0;
    unsigned char *p;
    STRING *acc = newSTRING(), *kw = newSTRING(), *cf = newSTRING();
    STRING **cftv = (STRING**)0;
    int i, len, indent = 0;
    long token, prev = DBLINEhash('X','X',' ');
    long reftoken, refprev = 0;
    PARSER *parser = newPARSER(fp);
    e->index = ftell(fp);
    e->sequence = newSTRING();
    while((len = linePARSER(parser)) != EOF){
        if(isspace(*parser->line))
            token = prev;
        else
            token = DBLINEhash(parser->line[0],
                               parser->line[1],
                               parser->line[2]);
        switch(token){
            case DBLINEhash('L','O','C'):
                 indent = 12;
                 if((DBPARSE_IDENTIFIER|DBPARSE_DIVISION
                    |DBPARSE_LENGTH|DBPARSE_DATE) & mask)
                     ParseGenBankLocus(e, (char*)parser->line+indent,
                                          len-indent);
                 break;
            case DBLINEhash('I','D',' '):
                 indent = 5;
                 if((DBPARSE_IDENTIFIER|DBPARSE_DIVISION
                    |DBPARSE_LENGTH|DBPARSE_TYPE) & mask)
                     ParseEmblIdentifier(e, (char*)parser->line+indent,
                                            len-indent);
                 break;
            case DBLINEhash('A','C','C'):
            case DBLINEhash('A','C',' '):
                 if(DBPARSE_ACCESSION & mask)
                     linetoSTRING(acc, (char*)parser->line+indent,
                                       len-indent);
                 break;
            case DBLINEhash('D','E','F'):
            case DBLINEhash('D','E',' '):
                 if(DBPARSE_DEFINITION & mask)
                     e->definition = dlinetoSTRING(e->definition,
                         (char*)parser->line+indent, len-indent, ' ');
                 break;
            case DBLINEhash('G','N',' '):
                 if(DBPARSE_GENENAME & mask)
                    e->genename = dlinetoSTRING(e->genename,
                        (char*)parser->line+indent, len-indent, ' ');
                 break;
            case DBLINEhash('D','T',' '):
                 if(DBPARSE_DATE & mask)
                     e->date 
                         = ConvertDBDate((char*)parser->line+indent);
                 break;
            case DBLINEhash('N','I',' '):
            case DBLINEhash('N','I','D'):
                 if(DBPARSE_NID & mask)
                     e->nid = strdup((char*)parser->line+indent);
                 break;
            case DBLINEhash('C','C',' '):
            case DBLINEhash('C','O','M'):
                 if(DBPARSE_COMMENT & mask)
                     e->comment = dlinetoSTRING(e->comment,
                       (char*)parser->line+indent, len-indent, '\n'); 
                 break;
            case DBLINEhash('O','S',' '):
                 if(DBPARSE_ORGANISM & mask)
                     e->organism = strdup((char*)parser->line+indent);
                 break;
            case DBLINEhash('O','C',' '):
                 if(DBPARSE_CLASSIFICATION & mask)
                     cf = dlinetoSTRING(cf, (char*)parser->line+indent, 
                                            len-indent, ' '); 
                 break;
            case DBLINEhash('O','G',' '):
                 if(DBPARSE_ORGANELLE & mask)
                     e->organelle = dlinetoSTRING(e->organelle, 
                                (char*)parser->line+indent, len-indent, ' '); 
                 break;
            case DBLINEhash('S','O','U'):
                 if((DBPARSE_ORGANISM|DBPARSE_CLASSIFICATION) & mask){
                     if(parser->line[0] != 'S'){
                         if(e->organism){
                             linetoSTRING(cf,
                                          (char*)parser->line+indent, 
                                          len-indent);
                         } else {
                             e->organism =
                                   strdup((char*)parser->line+indent);
                             }
                         }
                     }
                 break;
            case DBLINEhash('K','E','Y'):
            case DBLINEhash('K','W',' '):
                 if(DBPARSE_KEYWORD & mask)
                     linetoSTRING(kw, (char*)parser->line+indent,
                                  len-indent);
                 break; 
            case DBLINEhash('R','E','F'): 
                 if(DBPARSE_REFERENCE & mask){
                     if(token == prev){
                         if(isspace(parser->line[2])){
                             reftoken = refprev;
                         } else {
                             reftoken 
                                = DBLINEhash(parser->line[2],
                                             parser->line[3],
                                             parser->line[4]); 
                             r = e->referencev[e->referencec-1]; 
                         }
                         switch(reftoken){
                             case DBLINEhash('A','U','T'):
                                 r->author = dlinetoSTRING(r->author, 
                                           (char*)parser->line+indent,
                                           len-indent, ' '); 
                                 break;
                             case DBLINEhash('T','I','T'):
                                 r->title = dlinetoSTRING(r->title, 
                                            (char*)parser->line+indent,
                                            len-indent, ' '); 
                                 break;
                             case DBLINEhash('J','O','U'):
                                 r->journal = dlinetoSTRING(r->journal, 
                                            (char*)parser->line+indent,
                                            len-indent, ' '); 
                                 break;
                             case DBLINEhash('R','E','M'):
                                 r->comment = dlinetoSTRING(r->comment, 
                                            (char*)parser->line+indent,
                                            len-indent, ' '); 
                                 break;
                             case DBLINEhash('M','E','D'):
                                 if(r->dbxrefv){
                                     r->dbxrefc++;
                                     StepRealloc(r->dbxrefv, r->dbxrefc,
                                                STRINGLIST**);
                                 } else {
                                     r->dbxrefv 
                                         = InitAlloc(STRINGLIST**);
                                     r->dbxrefc = 1;
                                 }
                                 r->dbxrefv[r->dbxrefc-1] 
                                                 = newSTRINGLIST(); 
                                 addSTRINGLIST(
                                            r->dbxrefv[r->dbxrefc-1], 
                                            "MEDLINE", 7); 
                                 addSTRINGLIST(
                                            r->dbxrefv[r->dbxrefc-1], 
                                            (char*)parser->line+indent,
                                            len-indent);
                                 break;
                             } 
                         refprev = reftoken;
                         break; /* PREVENT FROM FALLING THROUGH */ 
                         }
                     refprev = DBLINEhash('X','X',' '); 
                     } /* ALLOW TO FALL THROUGH */
            case DBLINEhash('R','N',' '):
                 if(DBPARSE_REFERENCE & mask){
                     if(e->referencev){
                        StepRealloc(e->referencev, 
                                    e->referencec, DBREFERENCE*);
                        e->referencec++;
                     } else {
                        e->referencev = InitAlloc(DBREFERENCE**); 
                        e->referencec = 1; 
                        }
                     e->referencev[e->referencec-1]
                                = BLANK(DBREFERENCE);
                     }
                 break; 
            case DBLINEhash('F','E','A'):
                 token = DBLINEhash('F','T',' ');
                 break; 
            case DBLINEhash('F','T',' '):
                 if(DBPARSE_FEATURETABLE & mask){
                     if(!e->featuretablev) /* NEW FEATURE TABLE */
                         e->featuretablev = InitAlloc(DBFEATURE**); 
                     if(isspace(parser->line[5])){ 
                         if(parser->line[21] == '/'){/* NEW SUBRECORD */
                             StepRealloc(f->name,  f->count, STRING**);
                             StepRealloc(f->value, f->count, STRING**);
                             for(p = parser->line+21; *p; p++)
                                 if(*p == '=')
                                     break;
                             if(*p){ /* IF NAME VALUE PAIR */
                                 *p = '\0';
                                 f->value[f->count] = makeSTRING(
                                       (char*)p+1, 
                                       len-(p-parser->line)-1);
                                 cftv = &f->value[f->count];
                             } else {
                                 f->value[f->count] = newSTRING(); 
                             }
                             f->name[f->count] = makeSTRING(
                                         (char*)parser->line+22,
                                         len-22); /* CATCH 22 ? */
                             f->count++;  
                         } else { /* ADD TO CURRENT VALUE */
                             linetoSTRING(*cftv,
                                   (char*)parser->line+21, len-21); 
                         }
                     } else { /* NEW FT ENTRY */
                         StepRealloc(e->featuretablev, 
                                     e->featuretablec, DBFEATURE**);
                         e->featuretablev[e->featuretablec]
                             = BLANK(DBFEATURE);
                         f = e->featuretablev[e->featuretablec++];
                         for(p = parser->line+5; (!isspace(*p)); p++);
                         *p = '\0'; /* GET 1ST WORD */
                         f->feature = makeSTRING((char*)parser->line+5,
                                              p-parser->line-5);
                         f->location = makeSTRING(
                                       (char*)parser->line+21, len-21);
                         cftv = &f->location;
                         f->name = InitAlloc(STRING**);
                         f->value = InitAlloc(STRING**); 
                         }
                     } /* END OF FEATURETABLE PARSER */
                 break;
            case DBLINEhash('R','A',' '): 
                 if(DBPARSE_REFERENCE & mask)
                     e->referencev[e->referencec-1]->author
                          = makeSTRING((char*)parser->line+indent,
                                      len-indent);
                 break;
            case DBLINEhash('R','T',' '):
                 if(DBPARSE_REFERENCE & mask)
                    e->referencev[e->referencec-1]->title
                       = dlinetoSTRING(
                         e->referencev[e->referencec-1]->title,
                         (char*)parser->line+indent, len-indent, ' '); 
                 break;
            case DBLINEhash('R','C',' '):
                 if(DBPARSE_REFERENCE & mask)
                    e->referencev[e->referencec-1]->comment
                       = dlinetoSTRING(
                         e->referencev[e->referencec-1]->comment,
                         (char*)parser->line+indent, len-indent, ' '); 
                 break;
            case DBLINEhash('R','L',' '):
                 if(DBPARSE_REFERENCE & mask)
                    e->referencev[e->referencec-1]->journal
                       = makeSTRING((char*)parser->line+indent,
                                    len-indent); 
                 break;
            case DBLINEhash('R','P',' '):
                 if(DBPARSE_REFERENCE & mask)
                    e->referencev[e->referencec-1]->position
                       = makeSTRING((char*)parser->line+indent,
                                    len-indent); 
                 break;
            case DBLINEhash('R','X',' '):
                 if(DBPARSE_REFERENCE & mask){
                    r = e->referencev[e->referencec-1];
                    if(r->dbxrefv){
                        r->dbxrefc++;
                        StepRealloc(r->dbxrefv, r->dbxrefc, 
                                    STRINGLIST**);
                    } else {
                        r->dbxrefv = InitAlloc(STRINGLIST**);
                        r->dbxrefc = 1;
                        }
                    r->dbxrefv[r->dbxrefc-1] 
                     = ProcessDBLIST((char*)parser->line+indent,
                                     len-indent);
                    }
                 break;
            case DBLINEhash('D','R',' '):
                 if(DBPARSE_DBXREF & mask){
                     if(e->dbxrefv){
                          e->dbxrefc++;
                          StepRealloc(e->dbxrefv, e->dbxrefc, 
                                     STRINGLIST**); 
                     } else { 
                          e->dbxrefv = InitAlloc(STRINGLIST**);
                          e->dbxrefc = 1; 
                         }
                     e->dbxrefv[e->dbxrefc-1] =
                          ProcessDBLIST((char*)parser->line+indent,
                                        len-indent); 
                     }
                 break;
            case DBLINEhash('S','E','G'):
                if(DBPARSE_SEGMENT & mask)
                    e->segment = ParseSegment(
                              (char*)parser->line+indent, len-indent); 
            case DBLINEhash('B','A','S'):
                if(DBPARSE_BASECOUNT & mask)
                    e->basecount = ParseBaseCount(
                              (char*)parser->line+indent, len-indent); 
                break; 
            case DBLINEhash('O','R','I'):
                 if((DBPARSE_SEQUENCE & mask) && (prev == token)){ 
                     i = CleanSequenceLine((char*)parser->line);
                     linetoSTRING(e->sequence, (char*)parser->line, i);
                     }
                 break;
            case DBLINEhash('S','Q',' '):
                 if(prev != token){
                     if(DBPARSE_BASECOUNT & mask)
                         e->basecount = ParseBaseCount(
                                        (char*)parser->line+indent,
                                        len-indent);
                     } else {
                         if(DBPARSE_SEQUENCE & mask){
                             i = CleanSequenceLine(
                                   (char*)parser->line+indent);
                             linetoSTRING(e->sequence, 
                                   (char*)parser->line+indent, i); 
                             }
                         }
                 break;
            case DBLINEhash('/','/','\0'): 
            case DBLINEhash('/','/',' '): 
                 if(!e->sequence->len){
                     freeSTRING(e->sequence); 
                     e->sequence = (STRING*)0; 
                     }
                 if(kw->len) 
                     e->keyword = ProcessDBLIST(kw->str, kw->len);
                 freeSTRING(kw);
                 if(cf->len)
                     e->classification 
                                = ProcessDBLIST(cf->str, cf->len);
                 freeSTRING(cf);
                 if(acc->len)
                     e->accession = ProcessDBLIST(acc->str, acc->len);
                 freeSTRING(acc);
                 e->recordlength = (ftell(fp) - e->index);
                 for(i = 0; i < e->referencec; i++){
                     if(e->referencev[i]->author)
                         e->referencev[i]->author = 
                           CleanSTRING(e->referencev[i]->author);
                     if(e->referencev[i]->title)  
                         e->referencev[i]->title = 
                           CleanSTRING(e->referencev[i]->title);
                     if(e->referencev[i]->comment)  
                         e->referencev[i]->comment =  
                           CleanSTRING(e->referencev[i]->comment); 
                     }
                 for(i = 0; i < e->featuretablec; i++)
                     CleanFeatureTable(e->featuretablev[i]); 
                 freePARSER(parser);
                 return e;
                 break;
            default: prev = DBLINEhash('X','X',' ');
            }
        prev = token;
        }
    freePARSER(parser);
    return (DBENTRY*)0;
    } 

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

int main(int argc, char **argv){
    FILE *fp = (argc>1)?fopen(argv[1], "r"):stdin;
    DBENTRY *e;
    if(!fp) exit(1);
    while((e = GetNextDatabaseEntry(fp, 
          DBPARSE_ALL))){
        /* DBPARSE_ACCESSION */
         /* DBPARSE_ORGANISM */
        /* |DBPARSE_DEFINITION|DBPARSE_SEQUENCE */
           /* DBPARSE_ACCESSION|DBPARSE_IDENTIFIER */
          /* |DBPARSE_REFERENCE */
          /* |DBPARSE_DEFINITION|DBPARSE_ORGANISM */
          /* |DBPARSE_ORGANELLE|DBPARSE_LENGTH */
          /* |DBPARSE_FEATURETABLE|DBPARSE_SEQUENCE  */
          /* ))){ */
        printDBENTRY(stdout, e);
        freeDBENTRY(e);
        }
    fclose(fp);
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_PARSESEQDB_C */ 
 
/* TODO: NEED TO ADD SWISSPROT FEATURE TABLE PARSING */

