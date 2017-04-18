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

#ifndef INCLUDED_ARG_C
#define INCLUDED_ARG_C

/* #define ARG_CHECK_FOR_FORMAT_ERRORS */

#include <stdlib.h>
#include <string.h> /* FOR strdup() */
#include <ctype.h> /* FOR isprint() */

#include "arg.h"
#include "common.h"
#include "error.h"
#include "config.h"
#include "stringlist.h"
#include "dict.h"

ALLOW_ERROR_MESSAGES;

#define ARGUMENT_YEAR      1999
#define ARGUMENT_PACKAGE   "ESTate"
#define ARGUMENT_VERSION   "0.5.0"
#define ARGUMENT_INSTITUTE "UK-HGMP-RC"
#define ARGUMENT_URL       "http://www.hgmp.mrc.ac.uk/~gslater/"
#define ARGUMENT_EMAIL     "mailto:gslater@hgmp.mrc.ac.uk"
#define ARGUMENT_NAME      "Guy St.C. Slater."

static char *local_arg_names[ARGUMENT_TOTAL] =
  { "boolean", "letter", "short",  "int",  "long",  
    "float",   "double", "string", "path", "list"};

/* RULES FOR PROCESSING ARGUMENTS.
   Arguments are sought in this order.
    [1] Command line.
    [2] Web form data.
    [3] Config file.
    [4] Environment variable.
    [5] Hard_coded default.

   Other rules:
    [1] typical usage [-n <value> | --name <value>]
    [2] boolean usage [-d|--dna]
    [3] accept any number of booleans after hyphen,
        but only a single name for a value.
    [4] universal arguments -h|--help, -c|--config, -v|--version
    [5] '\0' or (char*)0 shows arg as unavailable.
    [6] format for hint string:
         NULL || "[ [<int>,<int>] || [expr] 
         || [opta,optb,etc] || [mode]"
         [<num>,<num>] : min-max for numbers.
         [expr] : grep like char range for char argument.
         [opta,optb,optn] : for strings
         [mode] : mode to use with fopen for FILE
    [6] If type is FILE, convert next arg of '-' to stdin.
    [7] Lists are colon deliminated. (apple:cat:bannana)
    [8] Any unflagged arguments are used to fill (in reverse order)
        from the end of the arg list.
    [9] For any arg flagged as vital, data must be found,
        (other than the hard-coded default).
   [10] Show current value, not default.

prepareARGUMENT(); PREPROCESS ARGS INTO LOCAL DATA STRUCT
  checkARGUMENT(); CHECK VALID ARGS. (ONLY WITH #ifdef DEBUG)
  usageARGUMENT(); PRINT USAGE
  parseARGUMENT(); MAIN WORK
processARGUMENT(); CALLED FROM estateMAIN();

*/

/* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

typedef struct {
    union {
    STRINGLIST *argstring;
       BOOLEAN  argchar[ALPHABETSIZE];
        double  argdouble[2];
          long  arglong[2];
          char *argfile;
    } hint;
     char *valuestr;
     char *spareval;
} PROCARGINFO;

static struct {
         int    argc; /* SET BY main() */
         char **argv; /* SET BY main() */
         char  *name;
         char  *defn;
          int   snametable[ALPHABETSIZE];
         DICT  *lnamedict;
   PROCARGINFO *info;
      BOOLEAN   fullusagerequest;
          int   sparectr;
          int   vitalcount;
    } local_arg;

int main(int argc, char **argv){
    int estateMAIN();
    register int status; 
    local_arg.argc = argc;
    local_arg.argv = argv;
    status = estateMAIN(); /* CALL MAIN PROGRAM FUNCTION */
    freeCONFIG();
    return status;
    }

/* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */
/* CODE FOR ARGUMENT PREPARATION                                  */

/* argchar = [expr]
      expr = se,expr
        se = char
           | charcharchar
           | chara-charb
*/
static void parseARGHINTchar(char *hint, BOOLEAN *result){
    register unsigned char *p = (unsigned char*)hint+1, c, prev = '\0';
    register int i;
    register BOOLEAN inrange = FALSE;
    for(i = 0; i < ALPHABETSIZE; i++)
        result[i] = FALSE;
    while(*p){
        if(*p == ','){
            inrange = FALSE;
        } else if(*p == '-'){
            inrange = TRUE;
        } else if(*p == ']'){
            return;
        } else if(isalpha(*p)){
            if(inrange)
                for(c = prev+1; c <= *p; c++)
                    result[c] = TRUE; 
            else
                result[*p] = TRUE; 
            inrange = FALSE;
            prev = *p; 
            }
        p++;
        }
    return;
    }

static void parseARGHINTstring(char *hint, STRINGLIST **argstring){
    register int len = strlen(hint);
    *argstring = toSTRINGLIST(hint+1, len-2, ',');
    return;
    }

static void parseARGHINTlong(char *hint, long *arglong){
    register int len = strlen(hint);
    register STRINGLIST *sl;
    sl = toSTRINGLIST(hint+1, len-1, ',');
    arglong[0] = atol(sl->v[0]);
    arglong[1] = atol(sl->v[1]);
    freeSTRINGLIST(sl);
    return;
    }

static void parseARGHINTdouble(char *hint, double *argdouble){
    register int len = strlen(hint);
    register STRINGLIST *sl;
    sl = toSTRINGLIST(hint+1, len-1, ',');
    argdouble[0] = atof(sl->v[0]);
    argdouble[1] = atof(sl->v[1]);
    freeSTRINGLIST(sl);
    return;
    }

static void parseARGHINTfile(char *hint, char **argfile){
    *argfile = strdup(hint);
    return;
    }

static PROCARGINFO *preparePROCARGINFO(ARGUMENT *a, int count){
    register int i;
    register PROCARGINFO *pah = malloc(sizeof(PROCARGINFO)*count);
    local_arg.lnamedict  = newDICT();
    memset(local_arg.snametable, 0, sizeof(int)*ALPHABETSIZE);
    local_arg.snametable['h'] = -1;
    local_arg.snametable['c'] = -2;
    local_arg.snametable['v'] = -3;
    uniqAddDICT(local_arg.lnamedict, "help",    (void*)-1);
    uniqAddDICT(local_arg.lnamedict, "config",  (void*)-2);
    uniqAddDICT(local_arg.lnamedict, "version", (void*)-3);
    for(i = 0; i < count; i++){
#ifdef ARG_CHECK_FOR_FORMAT_ERRORS
        if(local_arg.snametable[a[i].sname])
             errmsg(ERROR_FATAL, "sname used twice, [%c]", a[i].sname); 
        if(lookupDICT(local_arg.lnamedict, a[i].lname))
             errmsg(ERROR_FATAL, "lname used twice, [%c]", a[i].sname); 
#endif /* ARG_CHECK_FOR_FORMAT_ERRORS */
        local_arg.snametable[a[i].sname] = i+1;
        uniqAddDICT(local_arg.lnamedict, a[i].lname, (void*)(i+1));
        pah[i].valuestr = NULL;
        if(a[i].hint){ /* PARSE HINT STRINGS */
            switch(a[i].type){
                case ARG_CHAR:
                    parseARGHINTchar(a[i].hint, pah[i].hint.argchar);
                    break;
                case ARG_SHORT:
                case ARG_INT:
                case ARG_LONG:
                    parseARGHINTlong(a[i].hint, 
                                     pah[i].hint.arglong);
                    break;
                case ARG_FLOAT:
                case ARG_DOUBLE:
                    parseARGHINTdouble(a[i].hint,
                                       pah[i].hint.argdouble);
                    break;
                case ARG_STRING:
                    parseARGHINTstring(a[i].hint, 
                                       &pah[i].hint.argstring);
                    break;
                case ARG_FILE:
                    parseARGHINTfile(a[i].hint, 
                                     &pah[i].hint.argfile);
                    break;
                default:
                    break;
                }
            }
        }
    return pah;
    }

static void freePROCARGINFO(PROCARGINFO *pah, 
                               ARGUMENT *a, int count){
    register int i;
    for(i = 0; i < count; i++){
        if(a[i].hint){ /* PARSE HINT STRINGS */
            switch(a[i].type){
                case ARG_STRING:
                    freeSTRINGLIST(pah[i].hint.argstring);
                    break;
                case ARG_FILE:
                    free(pah[i].hint.argfile);
                    break;
                default:
                    break;
                }
            }
        }
    freeDICT(local_arg.lnamedict);
    return;
    }

/* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

#ifdef ARG_CHECK_FOR_FORMAT_ERRORS

static BOOLEAN checkARGNUMBERhint(char *str){
    register int len = strlen(str);
    register STRINGLIST *sl;
    if( (*str != '[') && (str[len] != ']'))
        return FALSE;
    sl = toSTRINGLIST(str+1, len-1, ',');
    if(sl->c != 2)
        return FALSE;
    if( (!strlen(sl->v[0])) || (!strlen(sl->v[1])) )
        return FALSE;
    freeSTRINGLIST(sl);
    return TRUE;
    }

static BOOLEAN checkARGSTRINGhint(char *str){
    register int len = strlen(str);
    register STRINGLIST *sl;
    register int i;
    if( (*str != '[') && (str[len] != ']'))
        return FALSE;
    sl = toSTRINGLIST(str+1, len-1, ',');
    for(i = 0; i < sl->c; i++)
        if(!strlen(sl->v[i]))
            return FALSE;
    freeSTRINGLIST(sl);
    return TRUE;
    }

static BOOLEAN checkARGCHARhint(char *str){
    register char *p = str+1;
    register BOOLEAN lastpunc = TRUE;
    if(*str != '[')
        return FALSE;
    while(*p){
        if(*p == ','){
            if(lastpunc)
                return FALSE;
            lastpunc = TRUE;
        } else if(*p == '-'){
            if(isalpha(p[-2]))
                return FALSE;
            if(lastpunc)
                return FALSE;
            lastpunc = TRUE;   
        } else if(isalpha(*p)){
            lastpunc = FALSE;
        } else if(*p == ']'){
            if(lastpunc)
                return FALSE;
            else 
                return TRUE;
        } else 
              return FALSE;
        p++;
        }
    return FALSE;
    }

static BOOLEAN checkARGFILEhint(char *str){
    register int i;
    char *argfilehinttypes[12] = 
      {"r",  "w",  "r+",  "w+",  "a",  "a+",
       "rb", "wb", "r+b", "w+b", "ab", "a+b"};
    for(i = 0; i < 12; i++) 
        if(!strcmp(str, argfilehinttypes[i]))
            return TRUE;
    return FALSE;
    }

static void checkARGUMENT(ARGUMENT *a, int count){
    register int i;
    for(i = 0; i < count; i++){
        if(!a[i].desc)
            errmsg(ERROR_FATAL, "No description for arg [%d]\n", i);
        if(!a[i].address)
            errmsg(ERROR_FATAL, "No address for [%s]\n", a[i].desc);
        if(!a[i].lname)
            errmsg(ERROR_FATAL, "No lname for [%s]\n", a[i].desc);
        if(!isprint(a[i].sname))
            errmsg(ERROR_FATAL, "No printable sname for [%s]\n",
                                 a[i].desc);
        if(!a[i].vname)
            errmsg(ERROR_FATAL, "No variable name for [%s]\n", 
                                a[i].desc);
        if(a[i].hint){
            switch(a[i].type){
                case ARG_CHAR:
                    if(!checkARGCHARhint(a[i].hint))
                        errmsg(ERROR_FATAL, "Bad char hint %s", 
                                a[i].hint);
                    break;
                case ARG_STRING:
                    if(!checkARGSTRINGhint(a[i].hint))
                        errmsg(ERROR_FATAL, "Bad string hint %s",
                                a[i].hint);
                    break;
                case ARG_SHORT:
                case ARG_INT:
                case ARG_LONG:
                case ARG_FLOAT:
                case ARG_DOUBLE:
                    if(!checkARGNUMBERhint(a[i].hint))
                        errmsg(ERROR_FATAL, "Bad number hint %s",
                                a[i].hint);
                    break;
                case ARG_FILE:
                    if(!checkARGFILEhint(a[i].hint))
                        errmsg(ERROR_FATAL, "Bad file hint %s",
                                a[i].hint);
                    break;
                default:
                    errmsg(ERROR_FATAL,
                           "Hint string inappropriate for %s",
                           local_arg_names[a[i].type]);
                    break;
                }
            }
        }
    errmsg(ERROR_INFO, "Argument data is valid");
    return;
    }

#endif /* ARG_CHECK_FOR_FORMAT_ERRORS */

/* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */
/* CODE FOR DISPLAYING USAGE INFORMATION                          */

static void usageheaderARGUMENT(){
    register int i, llen = 66;
    for(i = 0;  i < llen; i++)
        errputc('-');
    errprintf("\n%s version %s.     (C) %d. %s. %s.\n%s  %s\n",
            ARGUMENT_PACKAGE, ARGUMENT_VERSION, ARGUMENT_YEAR,
            ARGUMENT_NAME, ARGUMENT_INSTITUTE, ARGUMENT_URL,
            ARGUMENT_EMAIL);
    for(i = 0;  i < llen; i++)
        errputc('-');
    errputc('\n');
    return;
    }

static void shortusageARGUMENT(ARGUMENT *a, int count){
    register int i, j, total, namelen = strlen(local_arg.name);
    usageheaderARGUMENT();
    errprintf("%s : %s\n\n  Try: `%s --help` for more information.\n",
                           local_arg.name, local_arg.defn,
                           local_arg.name);
    total = errprintf("Usage: %s [ -hvc ] [ -", local_arg.name);
    for(i = 0; i < count; i++){
        errputc(a[i].sname);
        total++;
        }
    total += errprintf(" ]");
    for(i = 0; i < count; i++){
        if(a[i].vital){
            if((total+strlen(a[i].desc)+3) > 80){
                errputc('\n');
                for(total = j = namelen+4; j > 0; j--) 
                    errputc(' ');
                }
            total += errprintf(" <%s>", a[i].desc);
            }
        }
    errprintf("\n\n");
    return;
    }

static void printARGUMENTvalue(ARGUMENT_TYPE type, void *address){
    switch(type){
        case ARG_BOOLEAN:
            errprintf((*(BOOLEAN *)address)?"True":"False");
            break;
        case ARG_CHAR:
            errprintf("%c", (*(char *)address) );
            break;
        case ARG_SHORT:
            errprintf("%d", (*(short *)address) );
            break;
        case ARG_INT:
            errprintf("%d", (*(int *)address) );
            break;
        case ARG_LONG:
            errprintf("%ld", (*(long *)address) );
            break;
        case ARG_FLOAT:
            errprintf("%.3f", (*(float *)address) );
            break;
        case ARG_DOUBLE:
            errprintf("%.3f", (*(double *)address) );
            break;
        case ARG_STRING:
            if((*(char **)address))
                errprintf("%s", (*(char **)address) );
            else
                errprintf("( NULL )");
            break;
        case ARG_FILE:
            if((*(FILE **)address) == stdin)
                errprintf("(stdin)");
            else if((*(FILE **)address) == stdout)
                errprintf("(stdout)");
#ifdef ARG_CHECK_FOR_FORMAT_ERRORS
            else 
                errmsg(ERROR_FATAL, "Unknown FILE arg default");
#endif /* ARG_CHECK_FOR_FORMAT_ERRORS */
            break;
        case ARG_LIST:
#ifdef ARG_CHECK_FOR_FORMAT_ERRORS
            errmsg(ERROR_FATAL,
                   "Hard coded default for ARG_LIST not implemented");
#endif /* ARG_CHECK_FOR_FORMAT_ERRORS */
            break;
#ifdef ARG_CHECK_FOR_FORMAT_ERRORS
        default:
            errmsg(ERROR_FATAL, "Bad Argument type [%d]", a[pos].type);
            break;
#endif /* ARG_CHECK_FOR_FORMAT_ERRORS */
        }
    return;
    }

static void fullusageARGUMENT(ARGUMENT *a, int count){
    register int i;
    register char *tu, *ta = "  ", *tb = "[]";
    shortusageARGUMENT(a, count);
    errprintf("Command line options:\n\n"
              " [ -h | --help ] : show this help\n"
              " [ -v | --version ] : show version info\n"
              " [ -c | --config <path> ] : specify config file\n\n");
    for(i = 0; i < count; i++){
        tu = a[i].vital?ta:tb;
        errprintf(" %c -%c | --%s %s<%s>%s %c : %s\n", 
                 tu[0], a[i].sname, a[i].lname, 
                 (a[i].type == ARG_BOOLEAN)?"[ ":"",
                 local_arg_names[a[i].type], 
                 (a[i].type == ARG_BOOLEAN)?" ]":"",
                 tu[1], a[i].desc);
        errprintf(" Default set with $%s\n", a[i].vname);
        if(local_arg.info[i].valuestr)
            errprintf(" Currently set to \"%s\"\n", 
                                local_arg.info[i].valuestr);
        else  {
              if(a[i].vital)
                  errprintf(" (not set : *** value required ***)\n");
              else {
                  errprintf(" (using hard coded default : \"");
                  printARGUMENTvalue(a[i].type, a[i].address);
                  errprintf("\"\n");
                  }
              }
        if(a[i].hint){
            switch(a[i].type){
                case ARG_CHAR:
                    errprintf(" Acceptable chars: %s", a[i].hint);
                    break;
                case ARG_STRING:
                    errprintf(" Acceptable strings: %s", a[i].hint);
                    break;
                case ARG_SHORT:
                case ARG_INT:
                case ARG_LONG:
                    errprintf(" Acceptable numbers from [%ld] to [%ld]", 
                        local_arg.info[i].hint.arglong[0],
                        local_arg.info[i].hint.arglong[1]);
                    break;
                case ARG_FLOAT:
                case ARG_DOUBLE:
                    errprintf(
                        " Acceptable numbers from [%2.2f] to [%2.2f]", 
                        local_arg.info[i].hint.argdouble[0],
                        local_arg.info[i].hint.argdouble[1]);
                    break;
                case ARG_FILE:
                    errprintf(" Opening file with \"%s\"", 
                        local_arg.info[i].hint.argfile);
                    break;
                default:
                    break;
                }
            errputc('\n');
            }
        errputc('\n');
        }
    return;
    }

/* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

static int handleGeneralArgument(ARGUMENT *a, int count,
                                 int argid, char *next){
    switch(argid){
        case -1: /* -h|--help */
            local_arg.fullusagerequest = TRUE;
            break;
        case -2: /* -c|--config <name> */
            if(!next)
                errmsg(ERROR_FATAL, "No config path given");
            readCONFIG(local_arg.argv[0], next);
            return 1; /* USED NEXT */
            break;
        case -3: /* -v | --version */
            usageheaderARGUMENT();
            exit(1);
            break;
#ifdef ARG_CHECK_FOR_FORMAT_ERRORS
        default:
            errmsg(ERROR_FATAL, "Unknown argid [%d]", argid);
            break;
#endif /* ARG_CHECK_FOR_FORMAT_ERRORS */
        }
    return 0;
    }

static char *parseboolARGUMENT(char *str){
    register char *truestr = "TRUE", *falsestr = "FALSE";
    if(!str)
        return NULL;
    if(!strcasecmp(str, "YES"  )) return truestr;
    if(!strcasecmp(str, "ON"   )) return truestr;
    if(!strcasecmp(str, "TRUE" )) return truestr;
    if(!strcasecmp(str, "YEAH" )) return truestr;
    if(!strcasecmp(str, "NO"   )) return falsestr;
    if(!strcasecmp(str, "OFF"  )) return falsestr;
    if(!strcasecmp(str, "FALSE")) return falsestr;
    if(!strcasecmp(str, "NOPE" )) return falsestr;
    return NULL; /* NOT A BOOLEAN VALUE */
    }

static void setARGUMENTstr(ARGUMENT *a, int count, 
                           int argid, char *value){
    if(local_arg.info[argid].valuestr){
        shortusageARGUMENT(a, count);
        errmsg(ERROR_FATAL, "Multiple arguments given for \"%s\"",
                            a[argid].desc); 
        }
    local_arg.info[argid].valuestr = value;
    return;
    }

static int handleARGUMENTid(ARGUMENT *a, int count,
                            int argid, char *next){
    register char *boolresult; 
    if(argid > 0){ /* IS PROGRAM SPECIFIC ARGUMENT */
        argid--;
        if(a[argid].type == ARG_BOOLEAN){
            boolresult = parseboolARGUMENT(next); 
            if(boolresult){
                setARGUMENTstr(a, count, argid, boolresult);
                return 1;
                }
            setARGUMENTstr(a, count, argid, "TRUE");
            return 0;
        } else {
            if(!next){
                shortusageARGUMENT(a, count);
                errmsg(ERROR_FATAL, "Expected <%s> for %s", 
                                    local_arg_names[a[argid].type],
                                    a[argid].desc);
                exit(1);
                }
            setARGUMENTstr(a, count, argid, next);
            return 1;
            }
        }
    /* OTHERWISE IS GENERAL ARGUMENT */
    return handleGeneralArgument(a, count, argid, next);
    }

static int parseARGUMENTdoublehyphen(ARGUMENT *a, int count, 
                               char *argname, char *next){
    register int argid = (int)lookupDICT(local_arg.lnamedict,
                                         argname);
    if(!argid){
        shortusageARGUMENT(a, count);
        errmsg(ERROR_FATAL, "Unknown argument [--%s]", argname);
        exit(1);
        }
    return handleARGUMENTid(a, count, argid, next);
    }

static int parseARGUMENTsinglehyphen(ARGUMENT *a, int count, 
                          unsigned char *argname, char *next){
    register int i, argid, nextmove = 0;
    for(i = 0; argname[i]; i++){
        argid = local_arg.snametable[argname[i]];
        if(!argid){
            shortusageARGUMENT(a, count);
            errmsg(ERROR_FATAL, "Unknown option \"-%c\"", argname[i]);
            }
        nextmove += handleARGUMENTid(a, count, argid, next);
        if(nextmove > 1){
            shortusageARGUMENT(a, count);
            errmsg(ERROR_FATAL, "Ambigious option \"-%s %s\"", 
                                                   argname, next);
            }
        }
    return nextmove;
    }

static int parseARGUMENTzerohyphen(ARGUMENT *a, int count, 
                                     char *argname){
    if(local_arg.sparectr >= local_arg.vitalcount){
        shortusageARGUMENT(a, count);
        errmsg(ERROR_FATAL, "Too many unflagged arguments");
        }
    local_arg.info[local_arg.sparectr++].spareval = argname;
    return 0;
    }

static void parseARGUMENT(ARGUMENT *a, int count){
    register int i;
    for(i = 1; i < local_arg.argc; i++){
        if(local_arg.argv[i][0] == '-'){
            if(local_arg.argv[i][1] == '-'){ /* DOUBLE HYPHEN */
                i+= parseARGUMENTdoublehyphen(a, count,
                          local_arg.argv[i]+2,
                ((i+1)==local_arg.argc)?NULL:local_arg.argv[i+1]);
            } else { /* SINGLE HYPHEN */
                i+= parseARGUMENTsinglehyphen(a, count,
                    (unsigned char*)local_arg.argv[i]+1,
                ((i+1)==local_arg.argc)?NULL:local_arg.argv[i+1]);
                }
        } else { /* NO HYPHEN */
            i+= parseARGUMENTzerohyphen(a, count, local_arg.argv[i]);
            }
        }
    return;
    }

/* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

static void convertARGUMENT(ARGUMENT *a, int count,
                            int pos, char *str){
    register int i;
    register long l;
    register double d;
    register FILE *fp;
    switch(a[pos].type){
        case ARG_BOOLEAN:
            (*(BOOLEAN *)a[pos].address) 
                  = strcmp(str,"TRUE")?FALSE:TRUE;
            break;
        case ARG_CHAR:
            if((!str[0]) || str[1])
                errmsg(ERROR_FATAL, "Bad char value [%s]\n", str);
            (*(char *)a[pos].address) = str[0];
            break; 
        case ARG_SHORT:
            i = atoi(str);
            (*(short *)a[pos].address) = i;
            break;
        case ARG_INT:
            i = atoi(str);
            (*(int *)a[pos].address) = i;
            break;
        case ARG_LONG:
            l = atol(str);
            (*(long *)a[pos].address) = l;
            break;
        case ARG_FLOAT:
            d = atol(str);
            (*(float *)a[pos].address) = d;
            break;
        case ARG_DOUBLE:
            d = atol(str);
            (*(double *)a[pos].address) = d;
            break;
        case ARG_STRING:
            (*(char **)a[pos].address) = str;
            break;
        case ARG_FILE:
            if((str[0] == '-') && (str[1] == '\0'))
                fp = stdin;
            else {
                 fp = fopen(str, a[pos].hint);
                 if(!fp)
                     errmsg(ERROR_FATAL, "Could not open file [%s]",
                             str);
                 }
            (*(FILE **)a[pos].address) = fp;
            break;
        case ARG_LIST:
            (*(STRINGLIST**)a[pos].address) = toSTRINGLIST(str, 
                                                  strlen(str), ':');
            break;
#ifdef ARG_CHECK_FOR_FORMAT_ERRORS
        default:
            errmsg(ERROR_FATAL, "Bad Argument type [%d]", a[pos].type);
            break;
#endif /* ARG_CHECK_FOR_FORMAT_ERRORS */
        }
    return;
    }

/* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

static void confirmARGUMENTchar(unsigned char c, PROCARGINFO *info,
                                char *name){
    if(!info->hint.argchar[c])
        errmsg(ERROR_FATAL, "Unacceptable char --%s [%c]", name, c);
    return;
    }

static void confirmARGUMENTstring(char *s, PROCARGINFO *info, 
                                  char *name){
    register int i;
    for(i = 0; i < info->hint.argstring->c; i++){
        if(!strcmp(s, info->hint.argstring->v[i]))
            return;
        }
    errmsg(ERROR_FATAL, "Unsuitable string for --%s [%s]", name, s);
    return;
    }

static void confirmARGUMENTlong(long l, PROCARGINFO *info,
                                char *name){
    if(l < info->hint.arglong[0])
        errmsg(ERROR_FATAL, "Value for --%s is too low [%d]", name, l);
    if(l > info->hint.arglong[1])
        errmsg(ERROR_FATAL, "Value for --%s is too high [%d]", name, l);
    return;
    }

static void confirmARGUMENTdouble(double d, PROCARGINFO *info,
                                char *name){
    if(d < info->hint.argdouble[0])
        errmsg(ERROR_FATAL, "Value for --%s is too low [%.2f]",
               name, d);
    if(d > info->hint.argdouble[1])
        errmsg(ERROR_FATAL, "Value for --%s is too high [%.2f]",
               name, d);
    return;
    }

static void confirmARGUMENT(ARGUMENT *a, int count, int pos){
    if(a[pos].hint){ /* PARSE HINT STRINGS */
        switch(a[pos].type){
            case ARG_CHAR:
                confirmARGUMENTchar( (*(char*)a[pos].address),
                                   &local_arg.info[pos], 
                                   a[pos].lname);
                break;
            case ARG_STRING:
                confirmARGUMENTstring( (*(char**)a[pos].address),
                                     &local_arg.info[pos], 
                                     a[pos].lname);
                break;
            case ARG_SHORT:
            case ARG_INT:
            case ARG_LONG:
                confirmARGUMENTlong( (*(long*)a[pos].address),
                                   &local_arg.info[pos],
                                   a[pos].lname);
                break;
            case ARG_FLOAT:
            case ARG_DOUBLE:
                confirmARGUMENTdouble( (*(double*)a[pos].address),
                                     &local_arg.info[pos],
                                     a[pos].lname);
                break;
            default:
                break;
            }
        }
    return;
    }

static void fillremainingARGUMENT(ARGUMENT *a, int count){
    register int i, j;
    for(i = local_arg.sparectr-1, j = count; i >= 0; i--){ /* SPARES */
        while(!a[--j].vital);
        setARGUMENTstr(a, count, j, local_arg.info[i].spareval);
        }
    for(i = 0; i < count; i++) /* TRY TO FILL WITH CONFIG VALUE */
        if(!local_arg.info[i].valuestr)
            local_arg.info[i].valuestr = getCONFIG(a[i].vname);
    return;
    }

/* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

static void checkvitalARGUMENT(ARGUMENT *a, int count){
    register int i;
    for(i = 0; i < count; i++) /* CHECK VITAL ARGS */
        if((!local_arg.info[i].valuestr) && (a[i].vital)){
            shortusageARGUMENT(a, count);
            errmsg(ERROR_FATAL, "No value for -%c | --%s <%s>", 
                            a[i].sname, a[i].lname, a[i].desc);
            }
    return;
    }

static void transferARGUMENT(ARGUMENT *a, int count){
    register int i;
    for(i = 0; i < count; i++) /* CONVERT VALUE */
        if(local_arg.info[i].valuestr){
            convertARGUMENT(a, count, i, local_arg.info[i].valuestr);
            if(a[i].hint)
                confirmARGUMENT(a, count, i);
            if(a[i].string) /* STORE ORIGINAL STRING */
                *a[i].string = local_arg.info[i].valuestr;
            }
    return;
    }

/* processARGUMENT : CALLED FROM estateMAIN() TO PARSE COMMAND LINE 
*/
void processARGUMENT(ARGUMENT *a, int count, char *name, char *defn){
    register int i, namestart = 0, namelen = strlen(name);
    local_arg.name = malloc(sizeof(char)*namelen);
    local_arg.fullusagerequest = FALSE;
    local_arg.sparectr = 0;
    local_arg.vitalcount = 0;
    for(i = 0; i < count; i++)
        if(a[i].vital)
            local_arg.vitalcount++;
    while(name[namestart]){
        if(name[namestart] == '/')
            break;
        namestart++;
        }
    namestart = name[namestart]?(namestart+1):0;
    strncpy(local_arg.name, name+namestart,
                            sizeof(char)*(namelen-namestart-2));
    local_arg.name[namelen-namestart-2] = '\0';
    local_arg.defn = defn;
    local_arg.info = preparePROCARGINFO(a, count);
#ifdef ARG_CHECK_FOR_FORMAT_ERRORS
    checkARGUMENT(a, count);
#endif /* ARG_CHECK_FOR_FORMAT_ERRORS */
    parseARGUMENT(a, count); /* MAP COMMAND ARGS TO PROG ARGS */
    readCONFIG(local_arg.argv[0], getenv("ESTATECONFIGDIR"));
    fillremainingARGUMENT(a, count);
    if(local_arg.fullusagerequest){
        fullusageARGUMENT(a, count);
        exit(1);
        }
    checkvitalARGUMENT(a, count);
    transferARGUMENT(a, count);
    freePROCARGINFO(local_arg.info, a, count);
    free(local_arg.info);
    free(local_arg.name);
    return;
    }

/* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

#ifdef NOTNOW

static void showARGUMENTwww(){
    register char *myurl = ARGUMENT_MYURL;
    register int urllen = (strlen(myurl)+strlen(local_arg.name)+10);
    register char *url = malloc(sizeof(char)*urllen);
    register int comlen = urllen+70;
    register char *command = malloc(sizeof(char)*comlen);
    errmsg(ERROR_INFO, "Starting web help for [%s]", 
                        local_arg.name);
    sprintf(url, "%s#%s.html", ARGUMENT_MYURL, local_arg.name);
    sprintf(command, "netscape -noraise -remote "
                     "\'openURL(%s,new-window)\' 2> /dev/null", url);
    if(system(command))
        errmsg(ERROR_FATAL, "Help available at:\n%s", url);
    free(command);
    free(url);
    exit(1);
    return;
    } 
#endif /* NOTNOW */

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

#define ARGUMENT_COUNT 11

int estateMAIN(){
    static BOOLEAN b1 = FALSE;
    static BOOLEAN b2 = TRUE;
    static char c    = 'x';
    static short s   = 5;
    static int i     = 7;
    static char *intstr;
    static long l    = 12345;
    static float f   = 0.0;
    static double d  = 0.0;
    static FILE *fp;
    static char *fppath = NULL;
    static char *string = "def string";
    static STRINGLIST *list = NULL;

    register int ctr;

    ARGUMENT argdata[ARGUMENT_COUNT] = {
    {ARG_BOOLEAN, 'b', "bool1var", "test boolean one",
     &b1, NULL, "ARGTEST_BOOL1", NULL, FALSE},

    {ARG_BOOLEAN, 'B', "bool2var", "test boolean two",
     &b2, NULL, "ARGTEST_BOOL2", NULL, FALSE},

    {ARG_CHAR,    'C', "charvar", "test character",
     &c, NULL, "ARGTEST_CHAR", "[a-e,quot,p,r-s]", TRUE},

    {ARG_SHORT,   's', "shortvar", "test short",
     &s, NULL, "ARGTEST_SHORT",  NULL, FALSE},

    {ARG_INT,     'i', "intvar", "test integer",
     &i, &intstr, "ARGTEST_INT", "[15,30]", FALSE},

    {ARG_LONG,    'l', "longvar", "test long",
     &l, NULL,    "ARGTEST_LONG", "[1000,2000]", TRUE},

    {ARG_FLOAT,  'f', "floatvar", "test float",
     &f, NULL, "ARGTEST_FLOAT", NULL, FALSE},

    {ARG_DOUBLE,  'd', "dblvar", "test double",
     &d, NULL, "ARGTEST_DOUBLE", "[0,100]", TRUE},

    {ARG_FILE,    'F', "filepath", "test path",
     &fp, &fppath, "ARGTEST_FILE", "w",  FALSE},

    {ARG_STRING,  'S', "strvar",  "test string",
     &string, NULL, "ARGTEST_STRING", "[this,that]", FALSE},

    {ARG_LIST,  'L', "listvar", "test list", 
     &list,   NULL, "ARGTEST_LIST", NULL,   FALSE}
    };
    fp = stdout;
    processARGUMENT(argdata, ARGUMENT_COUNT, 
                    global_emn, "test program");
    errmsg(ERROR_INFO, "Argument processing results:\n"
                       "    bool1var = [%s]\n"
                       "    bool2var = [%s]\n" 
                       "    charvar  = [%c]\n" 
                       "    shortvar = [%d]\n" 
                       "    intvar   = [%d] <%s>\n" 
                       "    longvar  = [%d]\n" 
                       "    floatvar = [%2.2f]\n" 
                       "    dblvar   = [%.3e]\n" 
                       "    filepath = [%p] <%s>\n" 
                       "    string   = [%s]",

                       b1?"TRUE":"FALSE",
                       b2?"TRUE":"FALSE",
                       c,
                       s,
                       i, intstr,
                       l,
                       f,
                       d,
                       fp, fppath,
                       string);
    printf("    list     = {");
    if(list)
        for(ctr = 0; ctr < list->c; ctr++)
            printf("%s, ", list->v[ctr]);
    printf("}\n");
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_ARG_C */

/* TODO:
    Add min/max checking for integer types from limits.
    Join values from multi arguments for list.
    Need to free memory from ARG_LIST at program end.
*/

