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

#ifndef INCLUDED_CLUSTER_C
#define INCLUDED_CLUSTER_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* FOR memset() */

#include "../general/common.h"
#include "../parse/ioutil.h"
#include "clusterset.h"

#define CLUSTERSETCHUNKSIZE 100

/* newCLUSTERSET: MAKE A NEW CLUSTERSET STRUCTURE WITH TOTAL of members.
*/
CLUSTERSET *newCLUSTERSET(int total){
    register CLUSTERSET *cs = NEW(CLUSTERSET);
    cs->clustertotal = 0;
    cs->total = cs->alloced = cs->singletontotal = total;
    cs->tree = calloc(total+1, sizeof(int));
    cs->output = malloc((total+1)*sizeof(int));
    return cs;
    }

/* freeCLUSTERSET: DESTROY CLUSTERSET STRUCTURE cs.
*/
void freeCLUSTERSET(CLUSTERSET *cs){
    free(cs->tree);
    free(cs->output);
    free(cs);
    return;
    }

/* increaseCLUSTERSET: GROW CLUSTERSET TO ALLOW total MEMBERS.
*/
void increaseCLUSTERSET(CLUSTERSET *cs, int newtotal){
    register int diff = newtotal-cs->total, actualtotal;
    register size_t datasize;
    if(diff < 1)
        return; /* NOTHING TO BE DONE */
    cs->total += diff;
    cs->singletontotal += diff;
    if(newtotal >= cs->alloced){ /* IF NEED TO REALLOC */
        actualtotal = diff+CLUSTERSETCHUNKSIZE; /* BUFFER */
        datasize = sizeof(int)*actualtotal;
        cs->tree = realloc(cs->tree, datasize);
        cs->output = realloc(cs->output, datasize);
        memset(cs->tree+cs->alloced, 0,
               sizeof(int)*(actualtotal-cs->alloced));
        cs->alloced = actualtotal;
        }
    return;
    }

/*
   cs->tree[x] will contain -1 if x is a singleton.
                            y if x is a member of y.
*/

/* idlookupCLUSTERSET: RETURNS RAW CLUSTERSET ID FOR member.
*/
static int idlookupCLUSTERSET(CLUSTERSET *cs, int member){
    register int temp, comp, id = comp = member+1;
    while(cs->tree[id] > 0) /* LOOKUP ID */
        id = cs->tree[id];
    while(cs->tree[comp] > 0){ /* COMPRESS PATHS */
        temp = comp;
        comp = cs->tree[comp];
        cs->tree[temp] = id;
        }
    return id;
    }

/* lookupCLUSTERSET: CHECK TO SEE IF member IS IN A SET IN cs.
         RETURNS ZERO IF member IS A SINGLETON,
         OTHERWISE RETURNS THE CURRENT SET ID OF member.
*/
int lookupCLUSTERSET(CLUSTERSET *cs, int member){
    register int id = idlookupCLUSTERSET(cs, member);
    return cs->tree[id]?id:0;
    }

/* relateCLUSTERSET: TEST TO SEE IF a AND b ARE MEMBERS OF THE SAME SET.
         RETURNS ZERO IF EITHER OR BOTH A OR B ARE SINGLETONS.
                 -1 IF A AND B ARE MEMBERS OF DIFFERENT SETS.
                 CURRENT SET ID IF A AND B ARE MEMBERS OF THE SAME SET.
*/
int relateCLUSTERSET(CLUSTERSET *cs, int a, int b){
    register int ida, idb;
    ida = lookupCLUSTERSET(cs, a); 
    if(!ida) /* A IS SINGLETON */
        return 0;
    idb = lookupCLUSTERSET(cs, b);
    if(!idb) /* B IS A SINGLETON */
        return 0;
    return (ida == idb)?ida:-1; /* RETURN ID IF SAME, ELSE -1 */
    }

/* submitCLUSTERSET: JOIN A AND B IF EITHER OR BOTH ARE SINGLETONS.
          RETURNS ZERO IF ALREADY MEMBERS OF SAME SET.
          RETURNS ID OF NEW SET IF JOINED.
*/
int submitCLUSTERSET(CLUSTERSET *cs, int a, int b){
    register int ida = idlookupCLUSTERSET(cs, a),
                 idb = idlookupCLUSTERSET(cs, b);
    if(cs->tree[ida]){ /* (x,?) */
        if(cs->tree[idb]){ /* ([xy],[xy]) */
            if(ida == idb){ /* (x,x) */
                return 0; /* BOTH ALREADY IN SAME CLUSTER */
            } else { /* (x,y) */
                cs->clustertotal--; /* JOINING TWO CLUSTERS */
                }
        } else { /* (x,S) */
            cs->singletontotal--; /* B IS A SINGLETON */
            }
    } else { /* (S,?) */
        if(cs->tree[idb]){ /* (S,x) */
            cs->singletontotal--; /* A IS A SINGLETON */
        } else {  /* (S,S) */
            cs->singletontotal-=2; /* BOTH ARE SINGLETONS */
            cs->clustertotal++;    /* START A NEW CLUSTER */
            }
       }
    if(cs->tree[ida] < cs->tree[idb]){ /* ADD TO B */
        cs->tree[ida] += cs->tree[idb] - 1;
        cs->tree[idb] = ida;
        return idb;
        } /* ELSE ADD TO A */
    cs->tree[idb] += cs->tree[ida] - 1;
    cs->tree[ida] = idb;
    return ida;
    }

/* offersubmitCLUSTERSET: WILL JOIN A AND B IF EITHER OR BOTH
       ARE SINGLETONS AND IF CLUSTERSETCOMPFUNC RETURNS ZERO.
       RETURNS ZERO IF A AND B ARE ALREADY MEMBERS OF THE SAME SET.
                 (NO CALL TO cf WILL HAVE BEEN MADE)
       RETURNS -1 IF CLUSTERSETCOMPFUNC RETURNS NON-ZERO.
       RETURN CURRENT SET ID IN WHICH A AND B HAVE BEEN ADDED
*/
int offersubmitCLUSTERSET(CLUSTERSET *cs, int a, int b,
                void *info, CLUSTERSETCOMPFUNC cf){
    register int ida = idlookupCLUSTERSET(cs, a),
                 idb = idlookupCLUSTERSET(cs, b);
    /* fputc('?', stderr); */
    if(cs->tree[ida]){ /* (x,?) */
        if(cs->tree[idb]){ /* ([xy],[xy]) */
            if(ida == idb){ /* (x,x) */
                return 0; /* BOTH ALREADY IN SAME CLUSTER */
            } else { /* (x,y) */
                /* fputc('~', stderr); */
                if(cf(info, a, b)) /* ATTEMPT TO JOIN */
                    return -1;
                cs->clustertotal--; /* JOINING TWO CLUSTERS */
                }
        } else { /* (x,S) */
            /* fputc('<', stderr); */
            if(cf(info, a, b))
                return -1;
            cs->singletontotal--; /* B IS A SINGLETON */
            }
    } else { /* (S,?) */
        /* fputc('>', stderr); */
        if(cf(info, a, b))
            return -1;
        if(cs->tree[idb]){ /* (S,x) */
            cs->singletontotal--; /* A IS A SINGLETON */
        } else {  /* (S,S) */
            cs->singletontotal-=2; /* BOTH ARE SINGLETONS */
            cs->clustertotal++;    /* START A NEW CLUSTER */
            }
        }
    if(cs->tree[ida] < cs->tree[idb]){ /* ADD TO B */
        cs->tree[ida] += cs->tree[idb] - 1;
        cs->tree[idb] = ida;
        return idb;
        } /* ELSE ADD TO A */
    cs->tree[idb] += cs->tree[ida] - 1;
    cs->tree[ida] = idb;
    return ida;
    }

/* ----------------------------------------- */
/* START OF BLANK reportCLUSTERSET FUNCTIONS */

static void blankREPORTinit(void *info, int membertotal,
                            int clustertotal, int singletontotal){
    return;
    }

static void blankREPORTstart(void *info, int clusterid){
    return;
    }

static void blankREPORTfinish(void *info){
    return;
    }

static void blankREPORTmember(void *info, int memberid){
    return;
    }

/* END OF BLANK reportCLUSTERSET FUNCTIONS */
/* --------------------------------------- */

/* reportCLUSTERSET: REPORTS ON ALL CLUSTERS IN cs USING FUNCTIONS
                     IN rf WHEN AVAILABLE WITH info.
                     (AS A RESULT, WILL COMPRESS ALL CLUSTERSET PATHS)
*/
void reportCLUSTERSET(CLUSTERSET *cs, void *info, 
                      CLUSTERSETREPORTFUNCS *rf){
    register int i, currid, id, singleton = -1, cluster = -1;
    memset(cs->output, 0, sizeof(int)*(cs->total+1)); /* CLEAR OUTPUT */
    if(!rf->init)      rf->init      = blankREPORTinit;
    if(!rf->start)     rf->start     = blankREPORTstart;
    if(!rf->finish)    rf->finish    = blankREPORTfinish;
    if(!rf->member)    rf->member    = blankREPORTmember;
    for(i = 0; i < cs->total; i++){ /* FOR EACH ELEMENT */
        id = lookupCLUSTERSET(cs, i); /* LOOK UP ID */
        if(id){ /* MEMBER OF A CLUSTER */
            if((i+1) == id){ /* IS THE ROOT OF THE CLUSTER */
                if(!(cs->output[id] || (cluster == id))){ /* UNSEEN */
                    cs->output[id] = cluster; /* ADD ROOT */ 
                    cluster = id;
                    }
            } else { /* IS A MEMBER OF THE CLUSTER */
                /* ADD JUST AFTER CLUSTERLIST ROOT */
                if((!cs->output[id]) && (cluster != id)){ /* UNSEEN */
                    cs->output[id] = cluster; /* ADD ROOT EARLY */
                    cluster = id;
                    }
                /* ADD AFTER ROOT */
                cs->output[i+1] = cs->output[id];
                cs->output[id] = i+1;
                }
        } else { /* IS A SINGLETON */
            cs->output[i+1] = singleton; /* BUILD SINGLETON LIST */
            singleton = i+1;
            }
        }
    rf->init(info, cs->total, cs->clustertotal, cs->singletontotal); 
    if(rf->singleton) /* SHOW SINGLETONS */
        for(i = singleton; i >= 0; i=cs->output[i])
            rf->singleton(info, i);
    currid = -1;
    for(i = cluster; i >= 0; i=cs->output[i]){ /* SHOW CLUSTERS */
        if((id = cs->tree[i]) < 0)
            id = i;
        if(id != currid){
            if(currid != -1)
                rf->finish(info);
            rf->start(info, id);
            }
        rf->member(info, i);
        currid = id;
        }
    if(currid != -1)
        rf->finish(info);
    return;
    }

typedef struct {
    int currcount;
    int *hist;
    } CLUSTERSETHISTINFO;

static void histREPORTstart(void *info, int clusterid){
    CLUSTERSETHISTINFO *chi = (CLUSTERSETHISTINFO*)info;
    chi->currcount = 0;
    return;
    }

static void histREPORTfinish(void *info){
    CLUSTERSETHISTINFO *chi = (CLUSTERSETHISTINFO*)info;
    chi->hist[chi->currcount]++;
    return;
    }

static void histREPORTmember(void *info, int memberid){
    CLUSTERSETHISTINFO *chi = (CLUSTERSETHISTINFO*)info;
    chi->currcount++;
    return;
    }

static void printHistBar(FILE *fp, int num, BOOLEAN reverse){
    static char bar[32];
    register int i, width = num;
    if(width > 32)
        width = 32;
    for(i = 0; i < width; i++)
        bar[i] = '#';
    while(i < 32){
        bar[i] = '.';
        i++;
        }
    if(num > 32){
        if(reverse){
            i = sprintf(bar, "+%d", num>>5);
            bar[i] = '<';
            reversestring(bar+1, i-1);
        } else {
            i = sprintf(bar, "+%d", num>>5);
            bar[i] = '>';
            }
        }
    if(reverse){
        for(i = 31; i >= 0; i--)
           fputc(bar[i], fp);
    } else {
        fwrite(bar, sizeof(char), 32, fp);
        }
    return;
    }

void histogramCLUSTERSET(CLUSTERSET *cs, FILE *fp){
    register int i;
    CLUSTERSETHISTINFO chi;
    CLUSTERSETREPORTFUNCS crf = {NULL, NULL,
                                 histREPORTstart,
                                 histREPORTfinish,
                                 histREPORTmember};
    chi.hist = calloc(cs->total, sizeof(int));
    reportCLUSTERSET(cs, &chi, &crf);
    fprintf(fp, "Histogram of distribution of cluster sizes\n"
                "                                  size : number\n");
    printHistBar(fp, 1, TRUE);
    fprintf(fp, "%7d %-7d", 1, cs->singletontotal);
    printHistBar(fp, cs->singletontotal, FALSE);
    fputc('\n', fp);
    for(i = 0; i < cs->total; i++){
        if(chi.hist[i]){
            printHistBar(fp, i, TRUE);
            fprintf(fp, "%7d %-7d", i, chi.hist[i]);
            printHistBar(fp, chi.hist[i], FALSE);
            fputc('\n', fp);
            }
        }
    free(chi.hist);
    fputc('\n', fp);
    return;
    }

/* statusCLUSTERSET: WRITES INFO ABOUT CURRENT STATE OF cs TO fp.
*/
void statusCLUSTERSET(CLUSTERSET *cs, FILE *fp){
    fprintf(fp, "Clusterset has %d members\n"
                 "  Current number of clusters: %d\n"
                 "  Number of remaining singletons %d.\n\n",
                 cs->total,
                 cs->clustertotal, cs->singletontotal);
    return;
    }

#ifdef TEST_THIS_MODULE

#define VERTEXTOTAL 13
#define EDGETOTAL   17

typedef struct {
    int currid;      /* CLUSTERID BEING REPORTED */
    int membercount; /* MEMBERS SEEN IN CURRENT CLUSTER */
    int clustercount; /* CLUSTERS SEEN */
    } CLUSTERSETREPORTSTATUS;

static void testREPORTinit(void *info, int membertotal,
                            int clustertotal, int singletontotal){
    printf(">init membertotal=%d clustertotal=%d singletontotal=%d\n",
                  membertotal,   clustertotal,   singletontotal);
    return;
    }

#ifdef NOTNOW
static void testREPORTsingleton(void *info, int singletonid){
    printf(">singleton singletonid=%d\n", singletonid);
    return;
    }
#endif /* NOTNOW */

static void testREPORTstart(void *info, int clusterid){
    printf(">start clusterid=%d\n", clusterid);
    return;
    }

static void testREPORTfinish(void *info){
    printf(">finish\n");
    return;
    }

static void testREPORTmember(void *info, int memberid){
    printf(">member memberid=%d\n", memberid);
    return;
    }

int main(){
    register CLUSTERSET *cs = newCLUSTERSET(5);
    char *edge[EDGETOTAL] = 
        {"AG","AB","AC","LM","JM","JL","JK","ED","FD","HI",
         "FE","AF","GE","GC","GH","JG","LG"};
    register int i;
    CLUSTERSETREPORTFUNCS crf = {testREPORTinit, 
                                 /* testREPORTsingleton, */ NULL,
                                 testREPORTstart,
                                 testREPORTfinish,
                                 testREPORTmember};
    CLUSTERSETREPORTSTATUS crs;
    statusCLUSTERSET(cs, stdout);
    increaseCLUSTERSET(cs, VERTEXTOTAL);
    statusCLUSTERSET(cs, stdout);
    reportCLUSTERSET(cs, &crs, &crf);

    for(i = 0; i < EDGETOTAL; i++){
        submitCLUSTERSET(cs, edge[i][0]-'A', edge[i][1]-'A');
        statusCLUSTERSET(cs, stdout);
        reportCLUSTERSET(cs, &crs, &crf);
        }
    freeCLUSTERSET(cs);

    return 0;
    }

#endif /* TEST_THIS_MODULE */

#endif /* INCLUDED_CLUSTER_C */

/*
*/

