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

/* flat2coding : extract coding sequences from fasta format file.
*/

#include <ctype.h> /* isdigit() */
#include <string.h> /* for strstr() */
#include "../general/arg.h"
#include "../general/error.h"
#include "../general/common.h"
#include "../parse/parseseqdb.h"
#include "../sequence/sequtil.h"

ALLOW_ERROR_MESSAGES;

#define CDSToken_Position     (1L<<0) 
#define CDSToken_BracketOpen  (1L<<1) 
#define CDSToken_BracketClose (1L<<2) 
#define CDSToken_Comma        (1L<<3)  
#define CDSToken_DoubleDots   (1L<<4)  
#define CDSToken_Complement   (1L<<5)  
#define CDSToken_Join         (1L<<6)  
#define CDSToken_End          (1L<<7)  

typedef struct {
    long  token;
    void *value;
    } CDSTOKEN;

BOOLEAN CDStoTOKENS(STRING *cds, CDSTOKEN *ct){
    char *p;
    int loc, len = 0;
    for(p = cds->str; *p; p++){
        switch(*p){
            case ' ':
            case '\n':
            case '\t':
                break;
            case '(':
                ct[len++].token = CDSToken_BracketOpen;
                break;
            case ')':
                ct[len++].token = CDSToken_BracketClose;
                break;
            case ',':
                ct[len++].token = CDSToken_Comma;
                break;
            case '.':
                if(*++p == '.')
                    ct[len++].token = CDSToken_DoubleDots;
                else
                    return FALSE;
                break;
            case 'c':
                if(strncmp(p, "complement", 10))
                    return FALSE;
                ct[len++].token = CDSToken_Complement;
                p+=9;
                break;
            case 'j':
                if(strncmp(p, "join", 4))
                    return FALSE;
                ct[len++].token = CDSToken_Join;
                p+=3;
                break;
            case '0': case '1': case '2': case '3': case '4':
            case '5': case '6': case '7': case '8': case '9':
                for(loc = *p++-'0'; isdigit(*p); p++)
                    loc = (loc*10)+*p-'0';
                p--;
                ct[len].token = CDSToken_Position;
                ct[len++].value = (void*)loc;
                break;
            default:
                return FALSE;
            }
        }
    ct[len++].token = CDSToken_End;
    return TRUE;
    }

BOOLEAN CDSTOKENStoSTRING(STRING *s, STRING *seq,
                          CDSTOKEN *ct, long end, char *acc){
    int i, len, start, finish, prev;
    char *from;
    BOOLEAN state;
    STRING tmp;
    for(i = 0; ct[i].token != end; i++)
        switch(ct[i].token){
            case CDSToken_Position:
                if( (ct[i+1].token != CDSToken_DoubleDots)
                  ||(ct[i+2].token != CDSToken_Position) )
                    return FALSE;
                start  = (int)ct[i].value;
                finish = (int)ct[i+2].value;
                if(start < 1){
                    errmsg(ERROR_WARNING, "Negative start in [%s]\n",
                                    acc);
                    return FALSE;
                    }
                if(finish > seq->len){
                    errmsg(ERROR_WARNING,
                                  "Finish after end in [%s]\n", acc);
                    return FALSE;
                    }
                if(start > finish){
                    errmsg(ERROR_WARNING, "Reversed region in [%s]\n",
                                    acc);
                    return FALSE;
                    }
                from = seq->str+start-1;
                len  = finish-start+1;
                linetoSTRING(s, from, len);
                i+=2;
                break;
            case CDSToken_Join:
                if(ct[i+1].token != CDSToken_BracketOpen)
                    return FALSE;
                return CDSTOKENStoSTRING(s, seq, ct+i+2,
                                  CDSToken_BracketClose, acc);
                break;
            case CDSToken_Complement:
                if(ct[i+1].token != CDSToken_BracketOpen)
                    return FALSE;
                prev = s->len;
                state = CDSTOKENStoSTRING(s, seq, ct+i+2,
                                  CDSToken_BracketClose, acc);
                tmp.str = s->str+prev;
                tmp.len = s->len-prev;
                SEQUTILrevcomp((unsigned char*)tmp.str, tmp.len);
                return state;
            case CDSToken_Comma:
                break;
            case CDSToken_BracketOpen:
            case CDSToken_DoubleDots:
            case CDSToken_BracketClose:
            default:
                return FALSE;
                break;
            }
    return TRUE;
    }

/* ValidEndsCDS : CHECK FOR START AND IN-FRAME STOP CODON AT 
                  BEGINNING AND END OF THE SEQUENCE.         */
BOOLEAN ValidEndsCDS(STRING *seq){
    char *end = seq->str+seq->len;
    if(seq->len % 3) return FALSE;        /* CHECK FRAME */
    if(seq->str[0] != 'A') return FALSE;  /* CHECK START */
    if(seq->str[1] != 'T') return FALSE;
    if(seq->str[2] != 'G') return FALSE;
    if(end[-3] != 'T') return FALSE;      /* CHECK STOP  */
    if(end[-2] == 'A'){
        if( (end[-1] != 'G') && (end[-1] != 'A') )
            return FALSE;
    } else if(end[-2] == 'G'){
        if(end[-1] != 'A')
            return FALSE;
    } else return FALSE;
    return TRUE;
    }
            
static BOOLEAN checkCDSquality(DBFEATURE *f){
    register int i;
    for(i = 0; i < f->count; i++){
        if(!strncmp(f->name[i]->str, "partial", 7)){
            return FALSE;
            }
        if(strstr(f->value[i]->str, "putative"))
            return FALSE;
        if(strstr(f->value[i]->str, "hypothetical"))
            return FALSE;
        if(strstr(f->value[i]->str, "pot.")) /* POTENTIAL */
            return FALSE;
        }
    return TRUE;
    }

int printDBSELECTION(FILE *fp, DBENTRY *e, char *organism){
    int i, num = 0;
    BOOLEAN good;
    register char *p;
    DBFEATURE *f;
    CDSTOKEN *ct;
    STRING *coding;
    if(!e->accession)return 0; 
    if(!e->definition)return 0; 
    if(!e->featuretablec)return 0;
    if(!e->organism)return 0;
    if(!e->sequence)return 0;
    if(!e->length)return 0;
    if(e->organelle)return 0; /* IGNORE IF ORGANELLE */
    if(strcasecmp(e->organism, organism))return 0; /* CHECK ORGANISM */
    for(i = 0; i < e->featuretablec; i++){
        f = e->featuretablev[i];
        if(!strcmp(f->feature->str, "CDS")){
            if(!checkCDSquality(f))
                continue;
            ct = (CDSTOKEN*)malloc(sizeof(CDSTOKEN)*f->location->len);
            if(CDStoTOKENS(f->location, ct)){
                coding = newSTRING();
                if(CDSTOKENStoSTRING(coding, e->sequence, 
                           ct, CDSToken_End, e->accession->v[0]))
                    if(ValidEndsCDS(coding)){
                        good = TRUE; 
                        for(p = coding->str; *p; p++)
                            switch(*p){
                                case 'A': case 'C': 
                                case 'T': case 'G':
                                    break;  
                                default:
                                    good = FALSE;
                                }       
                        if(good){
                            fprintf(fp, ">%s_%d Length %d bp\n", 
                                e->accession->v[0], ++num, coding->len);
                            SEQUTILwriteFASTAblock(fp, coding->str,
                                                       coding->len);
                            }       
                        }       
                freeSTRING(coding);
                }
            free(ct);
            }
        }
    return num;
    }

static void parseDatabase(FILE *fp, char *organism){
    register int total = 0, cds = 0;
    register DBENTRY *e;
    while((e = GetNextDatabaseEntry(fp, 
           DBPARSE_ACCESSION|DBPARSE_IDENTIFIER
          |DBPARSE_DEFINITION|DBPARSE_ORGANISM
          |DBPARSE_ORGANELLE|DBPARSE_LENGTH
          |DBPARSE_FEATURETABLE|DBPARSE_SEQUENCE 
          ))){     
          /* printDBENTRY(stdout, e); */
          cds += printDBSELECTION(stdout, e, organism); 
          total++;
          freeDBENTRY(e);
          }       
    errmsg(ERROR_INFO, "Selected %d CDS from %d sequences\n",
                       cds, total);
    return;
    }

#define ARGUMENT_COUNT 2 

int estateMAIN(){
    static char *dbpath;
    static char *organism = "Homo sapiens (human)";
    static FILE *dbfp;
    ARGUMENT argdata[ARGUMENT_COUNT] = {
    { ARG_FILE,  'd', "database", "EMBL/GenBank flat file",
    &dbfp, &dbpath, "FLAT2CODING_DATABASE",  "r", TRUE},
    { ARG_STRING, 'o', "organism", "exact organism name",
    &organism, NULL, "FLAT2CODING_ORGANISM",   NULL,  FALSE}
    };
    dbfp = stdin;
    processARGUMENT(argdata, ARGUMENT_COUNT, global_emn, 
                    "extract coding sequences from a database");
    parseDatabase(dbfp, organism);
    return 0;
    }

/**/

