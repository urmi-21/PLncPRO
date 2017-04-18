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

#ifndef INCLUDED_CONFIG_C
#define INCLUDED_CONFIG_C

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/param.h> /* FOR realpath() */
#include <unistd.h>   /* FOR stat() */
#include <sys/stat.h> /* FOR stat() */
#include "common.h"
#include "string.h"
#include "config.h"
#include "error.h"

/* LOOK FOR ".ESTaterc" FILE IN : (IN THIS ORDER)
     -GIVEN PATH
     -CURRENT DIRECTORY
     -HOME DIRECTORY
     -EXECUTABLE DIRECTORY 
     -ESTate/etc DIRECTORY
*/

#define CONFIGTOKEN_END    0 
#define CONFIGTOKEN_NAME   (1L<<0)
#define CONFIGTOKEN_EQUALS (1L<<1)
#define CONFIGTOKEN_VALUE  (1L<<2)
#define CONFIGTOKEN_ERROR  (1L<<3)

#define CONFIG_FILE_NAME ".ESTaterc"
#define CONFIG_HOME_EVN  "HOME"
#define CONFIG_ETC_DIR   "../etc/"

ALLOW_ERROR_MESSAGES;

DICT *global_DICT_config = NULL;
#define MAXCONFIGLINE             140


/* getCONFIG : ATTEMPT TO GET VALUE AS FOLLOWS:
               configdata > envvariable
*/
char *getCONFIG(char *name){
      register char *value;
      value = lookupDICT(global_DICT_config, name);
      return (value)?(value):getenv(name);
      }

/* IMPLEMENTATION OF which (1) 
*/
static BOOLEAN whichPATH(char *prog, char *result, int maxresultlen){
    register char c, *s, *p, *path = getenv("PATH");
    register int proglen = strlen(prog);
    struct stat fileinfo;
    for(s = p = path; *p; p++){
        if(*p == ':'){
            c = *p;
            *p = '\0';
            if((proglen+(p-s)+2) > maxresultlen)
                errmsg(ERROR_FATAL, "Path too long for whichPATH");
            sprintf(result, "%s/%s", s, prog);
            *p = c;
            if((stat(result, &fileinfo) == 0) /* FILE EXISTS */
              && (fileinfo.st_mode & (S_IXUSR|S_IXGRP|S_IXOTH)))
                /* AND IS SOMEHOW EXECUTABLE */
                return TRUE; /* SUCCESS */
            s = p+1;
            }
        }
    return FALSE; /* FAIL */
    }

static FILE *findCONFIG(char *argvzero){
    FILE *fp;
    register char *etc = "../etc/";
    register char *home = getenv(CONFIG_HOME_EVN);
    register char *path, *p;
    register char *configfilename = CONFIG_FILE_NAME;
    register int lenc = strlen(configfilename), lenh = strlen(home),
        lena = argvzero?strlen(argvzero):0,
        lene = strlen(etc), lend;
    char realexecpath[MAXPATHLEN]; 
    if((fp = fopen(configfilename, "r"))) /* CURRENT DIRECTORY */
        return fp;
    path = (char*)malloc(sizeof(char)*lenc+Max(lenh,lena)+lene);
    sprintf(path, "%s/%s", home, configfilename);

    fp = fopen(path, "r"); /* HOME DIRECTORY */
    if(fp){
        free(path);
        return fp;
        } 
    if(!realpath(argvzero, realexecpath)){ /* IF CAN'T GET REAL PATH */
        if(!(whichPATH(argvzero, realexecpath, MAXPATHLEN))
            && (realpath(argvzero, realexecpath))){
            free(path);
            return (FILE*)0;
            }
        }
    lend = strlen(realexecpath);
    path = realloc(path, sizeof(char)*(lend+lenc+lenh)); 
    for(p = realexecpath+lend; p > realexecpath; p--)
        if(*p == '/'){
            p++;
            break;
            }
    *p = '\0';
    sprintf(path, "%s%s", realexecpath, configfilename);
    fp = fopen(path, "r");  /* EXEC DIRECTORY */
    if(!fp){
        sprintf(path, "%s%s%s", realexecpath, etc, configfilename+1);
        fp = fopen(path, "r"); /* ETC DIRECTORY */
        }
    free(path);
    return fp;
    }

static int getNextConfigToken(FILE *fp, void **v,
                              int *ln, long *ep, int pt){
    register int ch;
    register STRING *data = newSTRING();
    register BOOLEAN moved;
    *v = data;
    while((ch = getc(fp)) != EOF){
        switch(ch){
            case '=':
                if(pt != CONFIGTOKEN_NAME){
                    freeSTRING(data);
                    *v = "equals should be preceded with a name";
                    return CONFIGTOKEN_ERROR;
                    }
                return CONFIGTOKEN_EQUALS;
            case '\"':
                if(pt != CONFIGTOKEN_EQUALS){
                    freeSTRING(data);
                    *v = "value should be preceded with equals";
                    return CONFIGTOKEN_ERROR;
                    }
                while((ch = getc(fp)) != '\"'){
                    switch(ch){
                        case '\n': 
                            (*ln)++; 
                            *ep = ftell(fp);
                            break;
                        case '\\': 
                            switch(ch = getc(fp)){
                                case '\n': /* CONTINUE LINE */
                                    (*ln)++; 
                                    *ep = ftell(fp);
                                    moved = FALSE;
                                    while(isspace(ch = getc(fp)) ){
                                        if(ch == EOF)
                                            return CONFIGTOKEN_END;
                                        moved = TRUE;
                                        } 
                                    if(moved)
                                       addSTRING(data, ' ');
                                    ungetc(ch, fp);
                                    break;
                                case 'n': 
                                    addSTRING(data, '\n');
                                    break;
                                case 't': 
                                    addSTRING(data, '\t');
                                    break;
                                default:
                                    addSTRING(data, ch);
                                    break;
                                }
                            break;
                        case EOF:
                            return CONFIGTOKEN_END;
                        default:
                            addSTRING(data, ch);
                            break;
                        }
                    }
                return CONFIGTOKEN_VALUE;
            case '#': 
                while((ch = getc(fp)) != '\n')
                    if(ch == EOF)
                        return CONFIGTOKEN_END;
                    /* ELSE FALL THROUGH AS IF END OF LINE */
/*FALLTHROUGH*/
            case '\n': 
                (*ln)++;
                *ep = ftell(fp);
                break; 
            case ' ': case '\t': 
                break;
            default:
                if(pt != CONFIGTOKEN_VALUE){
                    freeSTRING(data);
                    *v = "Unexpected location for name";
                    return CONFIGTOKEN_ERROR;
                    }
                addSTRING(data, ch);
                while((ch = getc(fp)) != EOF){
                    if( (!isalnum(ch)) && (ch != '_') )
                        break;
                    addSTRING(data, ch);
                    }
                if(ch != '='){
                    freeSTRING(data);
                    *v = "equals sign missing";
                    return CONFIGTOKEN_ERROR;
                    }
                ungetc(ch, fp);
                return CONFIGTOKEN_NAME;
            }
        }
    return CONFIGTOKEN_END;
    }

static void expandValueString(STRING *s){
    register STRING *ns = newSTRING();
    register char *p = s->str, *end = s->str+s->len, *go;
    register char *expvalue;
    while(p < end){
        if( (*p == '$') && (*(p+1) == '{') ){
            go = p+=2;
            while(*p != '}'){
                if(p >= end)
                    errmsg(ERROR_FATAL, "Bad Expression [%s]", 
                                        s->str);
                p++;
                } 
            *p++ = '\0'; 
            if(!(expvalue = getCONFIG(go)))
                errmsg(ERROR_FATAL, "Unknown variable [%s]", go);
            linetoSTRING(ns, expvalue, strlen(expvalue)); 
            }
        addSTRING(ns, *p++);
        }
    s->str = ns->str; s->len = ns->len; s->alloc = ns->alloc;
    free(ns);
    return; 
    }

void readCONFIG(char *argvzero, char *path){
    register FILE *fp = (FILE*)0;
    register int i, ch, tok = CONFIGTOKEN_VALUE;
    register STRING *name = (STRING*)0, *value;
    register char *p;
    register long linepos;
    void *v;
    int linenum = 1;
    long errpos = 0;
    if(global_DICT_config) 
        return; /* IF ALREADY CALLED */
    if(path){
        fp = fopen(path, "r");
        if(!fp)
            errmsg(ERROR_INFO, "Could not find config file [%s]", path);
        fp = findCONFIG(argvzero);     
    } else {
        fp = findCONFIG(argvzero);     
        }
    global_DICT_config = newDICT();
    if(!fp){ /* IF NOT FOUND CONFIG FILE */
        errmsg(ERROR_INFO, "No configuration file found");
        return;
        }
    while((tok = getNextConfigToken(fp, &v, &linenum, &errpos, tok))){
        p = (char*)v;
        switch(tok){
            case CONFIGTOKEN_NAME:
                name = (STRING*)v;
                break;
            case CONFIGTOKEN_EQUALS:
                break;
            case CONFIGTOKEN_VALUE:
                value = (STRING*)v;
                expandValueString(value);
                if(!name)
                    errmsg(ERROR_FATAL, "Value without name [%s]",
                                         value->str); 
                if(uniqAddDICT(global_DICT_config, 
                                name->str, value->str))
                    errmsg(ERROR_FATAL, "Repeated Name [%s]", 
                                         name->str);
                freeSTRING(name);
                name = (STRING*)0;
                free(value);
                break;
            case CONFIGTOKEN_ERROR:
                linepos = ftell(fp)-errpos;
                fseek(fp, errpos, SEEK_SET);
                while((ch = getc(fp)) != EOF){
                    errprintf("%c", ch);
                    if(ch == '\n')
                        break;      
                    }
                for(i = 1; i < linepos; i++)
                    errputc(' ');
                errprintf("^--- %s\n", p); 
                errmsg(ERROR_FATAL, "Config error on line %d.", 
                                    linenum); 
            default:
                errmsg(ERROR_FATAL, "Syntax error in config file");
                exit(1);
            }
        }
    fclose(fp);
    return;
    }

void freeCONFIG(){
    if(global_DICT_config)
        freeDICT(global_DICT_config); 
    return;
    }

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

int main(int argc, char **argv){
    readCONFIG(*argv, NULL);
    printf("RESULT [%s]\n", getCONFIG("SWISSPROTFASTA") );
    freeCONFIG();
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_CONFIG_C */

