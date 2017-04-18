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

/* Clusterset: code to store and build a set of clusters.
   Guy St.C. Slater.. November 1998.

   See Sedgewick, Chapter 30 for more info on the
   path-compressed, weight balanced union-find algorithm used.
*/
#ifndef INCLUDED_CLUSTER_H
#define INCLUDED_CLUSTER_H

typedef struct {
    int *tree;           /* UNION-FIND TREE                    */
    int *output;         /* LIST STORAGE FOR OUTPUT            */
    int  total;          /* TOTAL NUMBER OF MEMBERS ALLOWED    */
    int  alloced;        /* MEMBERS FOR WHICH MEMORY ALLOCATED */
    int  singletontotal; /* NUMBER OF SINGLETONS REMAINING     */
    int  clustertotal;   /* CURRENT NUMBER OF CLUSTERS FORMED  */
    } CLUSTERSET;

/* CLUSTERSETCOMPFUNC: RETURNS ZERO IS COMPARISON AT OR ABOVE THRESHOLD
*/
typedef int (*CLUSTERSETCOMPFUNC) (void  *info, int a, int b);

/* CLUSTERSETREPORTFUNC: REPORT ON MEMBER OF A CLUSTERSET.
*/
/* typedef void (*CLUSTERSETREPORTFUNC) (void  *info, int memberid,  */
                                /* int clusterid, int count, int total); */

CLUSTERSET *       newCLUSTERSET(int total);
      void        freeCLUSTERSET(CLUSTERSET *cs);
      void    increaseCLUSTERSET(CLUSTERSET *cs, int newtotal);
       int      lookupCLUSTERSET(CLUSTERSET *cs, int member);
       int      relateCLUSTERSET(CLUSTERSET *cs, int a, int b);
       int      submitCLUSTERSET(CLUSTERSET *cs, int a, int b);
       int offersubmitCLUSTERSET(CLUSTERSET *cs, int a, int b,
                                 void *info, CLUSTERSETCOMPFUNC cf);
      void      statusCLUSTERSET(CLUSTERSET *cs, FILE *fp);
      void   histogramCLUSTERSET(CLUSTERSET *cs, FILE *fp);

      /* void      reportCLUSTERSET(CLUSTERSET *cs, void *info,  */
                    /* CLUSTERSETREPORTFUNC rf, BOOLEAN reportsingletons);
                     * */

typedef struct {
    void (*     init)(void *info, int membertotal, int clustertotal,
                                  int singletontotal);
    void (*singleton)(void *info, int singletonid);
    void (*    start)(void *info, int clusterid);
    void (*   finish)(void *info);
    void (*   member)(void *info, int memberid);
    } CLUSTERSETREPORTFUNCS;
/* NULL FUNCTIONS WILL NOT BE CALLED */

void reportCLUSTERSET(CLUSTERSET *cs, void *info,
                      CLUSTERSETREPORTFUNCS *rf);

#endif /* INCLUDED_CLUSTER_H */

