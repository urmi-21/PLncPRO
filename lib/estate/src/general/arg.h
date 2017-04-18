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

/* arg.h: Flexible command line interpretation code.
   Guy St.C. Slater..  November 1997.  Version 2.1.
*/

#ifndef INCLUDED_ARG_H
#define INCLUDED_ARG_H

#include "common.h"
#include <stdio.h>

typedef enum {
    ARG_BOOLEAN, ARG_CHAR,   ARG_SHORT,  ARG_INT,  ARG_LONG,
    ARG_FLOAT,   ARG_DOUBLE, ARG_STRING, ARG_FILE, ARG_LIST
} ARGUMENT_TYPE;

#define ARGUMENT_TOTAL 10

typedef struct {
    ARGUMENT_TYPE   type;    /* THE TYPE OF ARGUMENT TO FIND   */
    unsigned char   sname;   /* NAME OF COMMAND LINE CHAR      */
             char  *lname;   /* NAME OF COMMAND LINE STRING    */
             char  *desc;    /* DESCRIPTION OF THE OPTION      */
             void  *address; /* ADDRESS TO HOLD RESULT         */
             char **string;  /* ADDRESS TO POINT ORIG STRING   */
             char  *vname;   /* env/.rc/cgi VARIABLE NAME      */
             char  *hint;    /* STRING WITH OTHER INFO         */
          BOOLEAN   vital;   /* MUST A VALUE BE FOUND FOR THIS */
    } ARGUMENT;

void processARGUMENT(ARGUMENT *a, int count, char *name, char *defn);

#endif /* INCLUDED_ARGS_H */

