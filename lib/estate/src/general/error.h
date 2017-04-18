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

/* Error Message handling. 
   Guy St.C. Slater..  April 1997.
*/

#ifndef INCLUDED_ERROR_H
#define INCLUDED_ERROR_H

#include <stdio.h>
#include <stdarg.h>

extern int global_error_state;

#define ALLOW_ERROR_MESSAGES static char *global_emn = __FILE__

#define ERROR_ARGS     __LINE__, global_emn
#define ERROR_TYPES    3
#define ERROR_FATAL    0, ERROR_ARGS
#define ERROR_WARNING  1, ERROR_ARGS
#define ERROR_INFO     2, ERROR_ARGS

#define ERROR_STATE_STDERR 0
#define ERROR_STATE_STDOUT 1
#define ERROR_STATE_FILE   2

void errmsg(int type, int line, char *file, char *format, ...);
 int errprintf(char *format, ...);
void errputc(char c);

#endif /* INCLUDED_ERROR_H */

