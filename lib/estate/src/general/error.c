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

/* Error message writing code.
   Guy St.C. Slater.  April 1997.
*/

#ifndef INCLUDED_ERROR_C
#define INCLUDED_ERROR_C

#include <stdio.h>
#include <string.h>
#include "error.h"

int global_error_state = ERROR_STATE_STDERR;
FILE *local_error_fp = NULL;

void set_local_error_fp(){
    switch(global_error_state){
        case ERROR_STATE_STDERR:
            local_error_fp = stdout;
            break;
        case ERROR_STATE_STDOUT:
            local_error_fp = stderr;
            break;
        /* case ERROR_STATE_FILE: */
            /* local_error_fp = ???; */
            /* break;  */
        }
    return;
    }

void errmsg(int type, int line, char *file, char *format, ...){
    char *errtype[ERROR_TYPES] = {"Fatal", "Warning", "Info"};
    va_list ap;
    va_start(ap, format);
    if(!local_error_fp)
        set_local_error_fp();
    fprintf(local_error_fp, "#%s: (%.*s) ", errtype[type], 
                            (int)(strlen(file)-2), file); 
    vfprintf(local_error_fp, format, ap);
    va_end(ap);
#ifdef DEBUG
    if(!type) 
        fprintf(local_error_fp, "(Line %d)\n", line);
#endif
    putc('\n', local_error_fp);
    if(type)
       return;    
    exit(type);
    }

int errprintf(char *format, ...){
    register int total;
    va_list ap;  
    va_start(ap, format);
    if(!local_error_fp)
        set_local_error_fp();
    total = vfprintf(local_error_fp, format, ap);
    va_end(ap);
    return total;
    }

void errputc(char c){ 
    if(!local_error_fp)
        set_local_error_fp();
    fputc(c, local_error_fp);
    return;
    }

#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */
ALLOW_ERROR_MESSAGES;

int main(){
    errmsg(ERROR_INFO, "testing error library");
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */ 

#endif /* INCLUDED_ERROR_C */

