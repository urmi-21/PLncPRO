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

/* Config: Location and parsing of runtime configuration data. 
   Guy St.C. Slater..  February 1997.
*/

#ifndef INCLUDED_CONFIG_H
#define INCLUDED_CONFIG_H

#include "dict.h"

void readCONFIG(char *argvzero, char *path);

extern DICT *global_DICT_config;  
char *getCONFIG(char *name);
void  freeCONFIG();

#endif /* INCLUDED_CONFIG_H */

