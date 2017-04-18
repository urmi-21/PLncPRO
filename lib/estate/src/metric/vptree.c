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

#ifndef INCLUDED_VPTREE_C
#define INCLUDED_VPTREE_C

#include <stdio.h>
#include <stdlib.h>

#include <sys/types.h> /* FOR MKDIR   */
#include <sys/stat.h>  /* FOR MKDIR   */
#include <unistd.h>    /* FOR SYMLINK */
#include <signal.h>    /* FOR signal FOR oursearchVPTREE */

#include "../general/common.h"
#include "../general/error.h"
#include "vptree.h"
#include "metric.h"
#include "../parse/rafasta.h"
#include "../parse/ioutil.h"
#include "../struct/pqueue.h"
#include "../struct/list.h"

ALLOW_ERROR_MESSAGES;

/*
   NOTES:
   -----
   DISTANCE IS REUSED; IT IS USED FOR STORING DISTANCES
   FROM VANTAGE POINT PRIOR TO NODE SORTING.

   PARTITION WILL BE SET TO VPT_NO_PARTITION IF A SEQUENTIAL SEARCH
   IS REQUIRED OF ALL CHILD NODES. (PARTITIONING NOT PERFORMED).

   IMPLEMENTED SO THAT THE INNER SUBSPACE
   IS ALWAYS ON THE LEFT OF THE PARTITION.

   LOCATION IS REUSED; WHEN vpt->cache != NULL, LOCATION POINTS
   TO A SEQUENCE IN THE CACHE, NOT THE DATABASE POSITION.
*/

#define VPT_NO_PARTITION -1

typedef struct {
   long location;  /* LOCATION OF THE SEQUENCE IN DATABASE    */
    int partition; /* POSITION OF THE PARTITION NODE          */
  short distance;  /* DISTANCE TO THE PARTITION NODE (RADIUS) */
    } VPTREENODE;

typedef struct {
    VPTREENODE *root;    /* ROOT OF THE TREE */
           int  total;   /* NUMBER OF SEQUENCES IN THE TREE */
          char *dbpath;  /* PATH TO DATABASE */
          FILE *dbdata;  /* DATABASE */
        METRIC *metric;  /* TREE COMPARISON METRIC */
    } VPTREE;

typedef struct {
    char *seq;         /* DATABASE SEQUENCE */
     int  length;      /* SEQUENCE LENGTH */
    long  location;    /* CORRESPONDING FILE LOCATION */
    } VPTREECACHE;

typedef struct {
         VPTREE *vpt;     /* THE TREE */
    VPTREECACHE *cache;   /* SEQUENCE CACHE */
            int  comps;   /* COMPARISON COUNT */
            int  skips;   /* SKIP COUNT (REDUNDANCY) */
    } VPTREEBUILD;

/* newVPTREE : CREATE A NEW VP TREE.
*/
static VPTREE *newVPTREE(char *dbpath, METRIC *metric){
    register VPTREE *vpt = NEW(VPTREE);
    register size_t datasize;
    register int i;
    register long *idx;
    vpt->dbpath = strdup(dbpath);
    vpt->metric = metric;
    if(!(vpt->dbdata = fopen(dbpath, "r")))
        errmsg(ERROR_FATAL, "Could not open [%s]", dbpath);
    if(*dbpath != '/')
        errmsg(ERROR_FATAL, "Absolute path for [%s] required", dbpath);
    errmsg(ERROR_INFO, "Counting sequences in [%s]", dbpath);
    idx = (long*)getRAFASTAindices(vpt->dbdata, &vpt->total);
    errmsg(ERROR_INFO, "Making vptree from %d sequences", vpt->total);
    datasize = vpt->total * sizeof(VPTREENODE);
    if(!(vpt->root = (VPTREENODE*)malloc(datasize)))
        errmsg(ERROR_FATAL, "Insufficient memory for VPTREE size=%d",
               vpt->total);
    memset(vpt->root, 0, datasize);
    for(i = 0; i < vpt->total; i++)
        vpt->root[i].location = idx[i];
    free(idx);
    return vpt;
    }

/* freeVPTREE : FREE AND DESTROY A VPTREE.
*/
static void freeVPTREE(VPTREE *vpt){
    freeMETRIC(vpt->metric);
    fclose(vpt->dbdata);
    free(vpt->dbpath);
    free(vpt->root);
    free(vpt);
    return;
    }

/* _______________________________________________________________
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DATABASE SEQUENCE CACHE ROUTINES .

   HOW IT WORKS:
       o ANY SECTION WITH LESS THAN VPTREE_CACHE_THRESHOLD
         SEQUENCES IN IT WILL BE CACHED BEFORE CLUSTERING.
       o (VPTREENODE*)->location WILL REPLACED WITH THE LOCATION
         OF THE SEQUENCE ON THE ARRAY.
*/

#define VPTREE_CACHE_THRESHOLD 10000 /* ZERO TO TURN CACHE OFF */

void newVPTREECACHE(VPTREEBUILD *b, int start, int length){
    register int i;
    b->cache = malloc(sizeof(VPTREECACHE)*length);
    errmsg(ERROR_INFO, "Caching [%d] seqs", length);
    for(i = 0; i < length; i++){
        b->cache[i].location = b->vpt->root[start+i].location;
        b->vpt->root[start+i].location = i;
        b->cache[i].seq = getRAFASTAseq(b->vpt->dbdata,
             b->cache[i].location, &b->cache[i].length);
        }
    return;
    }

void freeVPTREECACHE(VPTREEBUILD *b, int start, int length){
    register int i, j;
    errmsg(ERROR_INFO, "Emptying cache of [%d] seqs", length);
    for(i = 0; i < length; i++){
        j = b->vpt->root[start+i].location;
        b->vpt->root[start+i].location = b->cache[j].location;
        /* b->vpt->root[start+i].location = b->cache[i].location;  */
        free(b->cache[j].seq);
        }
    free(b->cache);
    b->cache = NULL;
    return;
    }

/* _______________________________________________________________
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* compVPTREE : COMPARISON OF DISTANCES USED BY SORTING/MEDIAN ROUTINE.
*/
static int compVPTREE(VPTREENODE *a, VPTREENODE *b){
    return (a->distance - b->distance);
    }

/* partitionVPTREE : FIND THE OPTIMAL PARTITION POINT FOR THE
                     DISTRIBUTION OF SCORES GIVEN.  THIS WILL
                     BE THE POINT CLOSEST TO THE CENTRE, AT WHICH
                     AT WHICH THE DISTANCE AT THE PARTITION IS
                     DIFFERENT TO THE DISTANCE AT THE POINT
                     IMMEDIATELY TO THE LEFT OF THE PARTITION.
   SECTION LENGTH MUST BE >=2 TO USE THIS FUNCTION.
*/
static int partitionVPTREE(VPTREE *vpt, int start, int length){
    register int centre = (length+1)>>1;
    register int sw = ((length)>>1)-1;
    register int i, split = -1;
    register VPTREENODE *n = vpt->root+start;
    for(i = 1; i <= sw; i++){
        if(n[centre].distance != n[centre-i].distance){
            split = centre-i+1; /* PARTITION NEARER */
            break;
            }
        if(n[centre].distance != n[centre+i].distance){
            split = centre+i;   /* PARTITION FORWARD */
            break;
            }
        }
    if(split == -1){ /* FAILED TO FIND PARTITION SO FAR */
        if( (length == (length|1)) /* ODD */
         && (n[centre].distance != n[1].distance) )
            split = 2;    /* CHECK ODD COLUMN */
        else if(n[centre].distance != 0)
            split = 1;    /* LINEAR PARTITIONING */
        else {
            /* errmsg(ERROR_WARNING, "Could not partition");  */
            split = VPT_NO_PARTITION;
            }
        }
    return split;
    }

/* swapVPTREE : SWAP NODES a AND b.
*/
static void swapVPTREE(VPTREENODE *a, VPTREENODE *b){
    register int temp;
    Swap(a->location,  b->location,  temp);
    Swap(a->partition, b->partition, temp);
    Swap(a->distance,  b->distance,  temp);
    return;
    }

/* selectVPTREE : SELECT A NODE TO USE AS THE VANTAGE POINT
                  FOR THE SUBSPACE OF start...length
                  ADD MORE INTELLIGENT SELECTION ALGORITHM LATER.
*/
static void selectVPTREE(VPTREE *vpt, int start, int length){
    register int centre = (length+1)>>1;
    swapVPTREE(vpt->root+start, vpt->root+start+centre);
    /* swapVPTREE(vpt->root+start, vpt->root+start+length-1);  */
    return;
    }

/* scoreVPTREE : ASSIGN DISTANCE SCORE FOR MEMBER UPTO length
                 SEQUENCES AWAY FROM start IN VPTREE vpt.
*/
static void scoreVPTREE(VPTREEBUILD *b, int start, int length){
    int len;
    register char *seq;
    register void *sc;
    register int i, z;
    if(b->cache){
        seq = b->cache[b->vpt->root[start].location].seq;
        sc = b->vpt->metric->prep(seq, b->vpt->metric->param);
        for(i = start+1, z = start+length; i < z; i++)
            b->vpt->root[i].distance = b->vpt->metric->cache(sc,
                   b->cache[b->vpt->root[i].location].seq,
                   b->cache[b->vpt->root[i].location].length);
        b->vpt->metric->free(sc);
    } else {
        seq = getRAFASTAseq(b->vpt->dbdata,
                            b->vpt->root[start].location, &len);
        sc = b->vpt->metric->prep(seq, b->vpt->metric->param);
        for(i = start+1, z = start+length; i < z; i++)
            b->vpt->root[i].distance = b->vpt->metric->dist(sc,
                        b->vpt->root[i].location, b->vpt->dbdata);
        b->vpt->metric->free(sc);
        free(seq);
        }
    return;
    }


/* buildVPTREErecur : RECURSIVE FUNCTION FOR VPTREE BUILDING
*/
static void buildVPTREErecur(VPTREEBUILD *b, int start, int length){
    register VPTREENODE *n = b->vpt->root+start; /* SUBTASK ROOT */
    register int p;
    register BOOLEAN madecache = FALSE;
    switch(length){
        case 0:  break; /* SHOULD NEVER HAPPEN HERE */
        case 1:
            n->partition = 0;
            n->distance  = 0;
            break;
        case 2: /* AVOID SELECTION IF WIDTH IS 2 */
            selectVPTREE(b->vpt, start, 2);
            scoreVPTREE(b, start, 2);
            b->comps++; /* NO NEED FOR SORT */
            if(n->distance == 0){
                n->partition = VPT_NO_PARTITION;
                errmsg(ERROR_WARNING, "Found identical PAIR");
                b->skips++;
            } else {
                n->partition = 1;
                buildVPTREErecur(b, start+1, 1);
                }
            break;
        default:
            if(!b->cache)
                if(length < VPTREE_CACHE_THRESHOLD){
                    newVPTREECACHE(b, start, length);
                    madecache = TRUE; /* CACHE IN USE */
                    }
            selectVPTREE(b->vpt, start, length);
            scoreVPTREE(b, start, length);
            b->comps+=length;
            qsort(n+1, length-1, sizeof(VPTREENODE),
                 (int(*)(const void *, const void *))compVPTREE);
            /* NEED TO REPLACE WITH A MERGESORT ROUTINE HERE */
            /* CANNOT USE MEDIAN ALGORITHM AS NON-PERFECT BALANCING */
            n->partition = partitionVPTREE(b->vpt, start, length);
            if(n->partition == VPT_NO_PARTITION){
                errmsg(ERROR_WARNING, "Found identical GROUP of %d",
                       length);
                b->skips+=length;
            } else {
                p = n->partition;
                n->distance = n[p].distance;
                if(p > 1) /* FIXES PARTITION=1 BUG */
                    buildVPTREErecur(b, start+1, p-1);
                buildVPTREErecur(b, start+p, length-p);
                }
            break;
        }
    if(madecache) /* CACHE OFF */
        freeVPTREECACHE(b, start, length);
    return;
    }

/* newVPTREEBUILD : MAKE NEW BUILD DATA STRUCTURE.
*/
VPTREEBUILD *newVPTREEBUILD(VPTREE *vpt){
    VPTREEBUILD *b = NEW(VPTREEBUILD);
    b->comps = 0;
    b->skips = 0;
    b->vpt = vpt;
    b->cache = NULL;
    return b;
    }

/* buildVPTREE : BUILD A NEW VPTREE.
*/
static void buildVPTREE(VPTREE *vpt){
    VPTREEBUILD *b = newVPTREEBUILD(vpt);
    buildVPTREErecur(b, 0, vpt->total);
    errmsg(ERROR_INFO, "Built vptree of %d elements "
                       "with %d comparisons and %d skips",
                       vpt->total, b->comps, b->skips);
    free(b);
    return;
    }

#define VPTREE_FILE_DATABASE 0
#define VPTREE_FILE_INDEX    1
#define VPTREE_FILE_METRIC   2

/* getVPTREEpath : MAKE A PATH FOR A VPTREE COMPONENT.
                   USED BY readVPTREE AND writeVPTREE.
*/
static char *getVPTREEpath(char *idxpath, char component){
    register int idxlen = strlen(idxpath);
    char *files[3] = {"database", "index", "metric"};
    register char *fullpath = malloc(sizeof(char)
                         *(idxlen+strlen
                               (files[(unsigned char)component])+2));
    sprintf(fullpath, "%s/%s", idxpath,
                                files[(unsigned char)component]);
    return fullpath;
    }

/* writeVPTREE : WRITE A VPTREE TO DISK.
                 STRUCTURE: DIRECTORY <path>
                 CONTAINING THE FOLLOWING FILES:
                 index    : THE VPTREE INDEX
                 database : LOGICAL LINK TO THE DATABASE
                 metric   : DESCRIPTION OF THE METRIC USED
*/
/* NEED TO SHARE CODE BETWEEN read AND write FUNCTIONS */
static void writeVPTREE(VPTREE *vpt, char *idxpath){
    struct stat info;
    register BOOLEAN dirmade = TRUE;
    register char *cpath;
    register FILE *fp;
    if(mkdir(idxpath, FILE_ALLOW_U_W_UGO_RX) == 0) /* MAKE DIR */
        chmod(idxpath, FILE_ALLOW_U_W_UGO_RX); /* SG LIBRARY BUG FIX */
    else {
        if(stat(idxpath, &info))
           dirmade = FALSE;
        else if(S_ISDIR(info.st_mode)){
                 cpath = malloc(sizeof(char)*(strlen(idxpath)+8));
                 sprintf(cpath, "%s_update", idxpath);
                 errmsg(ERROR_WARNING,
                      "Writing to update directory [%s]", cpath);
                 writeVPTREE(vpt, cpath);
                 free(cpath);
                 return;
             } else
                 dirmade = FALSE;
        if(dirmade == FALSE)
            errmsg(ERROR_FATAL, "Could not make directory [%s]",
                                idxpath);
        }
    cpath = getVPTREEpath(idxpath, VPTREE_FILE_DATABASE);
    if(symlink(vpt->dbpath, cpath)) /* MAKE THE DB LINK */
        errmsg(ERROR_FATAL, "Could not make link to [%s]", vpt->dbpath);
    chmod(cpath, FILE_ALLOW_U_W_UGO_RX);
    free(cpath);

    cpath = getVPTREEpath(idxpath, VPTREE_FILE_INDEX);
    if(!(fp = fopen(cpath, "wb")))
        errmsg(ERROR_FATAL,
               "Could not open [%s] to write vptree", cpath);
    fwrite(vpt->root, sizeof(VPTREENODE), vpt->total, fp);
    fchmod(fileno(fp), FILE_ALLOW_U_W_UGO_R);
    fclose(fp);
    free(cpath);

    cpath = getVPTREEpath(idxpath, VPTREE_FILE_METRIC);
    if(!(fp = fopen(cpath, "w")))
        errmsg(ERROR_FATAL,
               "Could not open [%s] to write metric", cpath);
    writeMETRIC(fp, vpt->metric);
    fchmod(fileno(fp), FILE_ALLOW_U_W_UGO_R);
    fclose(fp);
    free(cpath);

    return;
    }

/* readVPTREE : READ A VPTREE FROM DISK.
*/
static VPTREE *readVPTREE(char *idxpath){
    register VPTREE *vpt = NEW(VPTREE);
    register char *cpath;
    register FILE *fp;
    struct stat s;

    cpath = getVPTREEpath(idxpath, VPTREE_FILE_INDEX);
    fp = fopen(cpath, "rb");
    if(!fp)
        errmsg(ERROR_FATAL,
               "Could not open [%s] to read vptree", cpath);
    if(fstat(fileno(fp), &s) == -1)
        errmsg(ERROR_FATAL,
               "Could not read file info for [%s]\n", cpath);
    vpt->total = s.st_size/sizeof(VPTREENODE);
    if(!(vpt->root = (VPTREENODE*)malloc(s.st_size)))
        errmsg(ERROR_FATAL,
            "Could not allocate %d bytes for [%s]\n",
                                             s.st_size, cpath);
    fread(vpt->root, sizeof(VPTREENODE), vpt->total, fp);
    fclose(fp);
    free(cpath);

    cpath = getVPTREEpath(idxpath, VPTREE_FILE_METRIC);
    if(!(fp = fopen(cpath, "r")))
        errmsg(ERROR_FATAL,
               "Could not open [%s] to read metric", cpath);
    vpt->metric = readMETRIC(fp);
    fclose(fp);
    free(cpath);

    cpath = getVPTREEpath(idxpath, VPTREE_FILE_DATABASE);
    vpt->dbpath = cpath;  /* THUS DO NOT FREE cpath HERE */
    if(!(vpt->dbdata = fopen(cpath, "r")))
        errmsg(ERROR_FATAL, "Could not open database, [%s]", cpath);
    return vpt;
    }

/* makeVPTREE : BUILD A VPTREE FROM dbpath USING metric
                AND WRITE IT OUT TO idxpath.
*/
void makeVPTREE(char *dbpath, char *idxpath, METRIC *metric){
    VPTREE *vpt = newVPTREE(dbpath, metric);
    buildVPTREE(vpt);
    writeVPTREE(vpt, idxpath);
    freeVPTREE(vpt);
    return;
    }
/* ___________________________________________________________ */

typedef struct {
     int start;     /* START OF THIS REGION         */
     int length;    /* LENGTH OF THIS REGION        */
   short depth;     /* DEPTH IN THE TREE            */
 BOOLEAN infront;   /* IS A FORWARD SUBSPACE        */
     int potential; /* NEAREST POSSIBLE DESCENDANT  */
     int distance;  /* DISTANCE FROM THE QUERY      */
#define VPTREE_PERFORMANCE_DATA
#ifdef VPTREE_PERFORMANCE_DATA
     int obsvcomps;  /* COMPARISONS WHEN ENCOUNTERED */
     int usedcomps;  /* COMPARISONS WHEN DIST MEASD  */
#endif /* VPTREE_PERFORMANCE_DATA */
     } VPTREEOUR;

static struct {
    void (*sig)(int);
    int comps;  /* NUMBER OF COMPARISONS   */
    int skips;  /* NUMBER OF SKIPS MADE    */
    int hits;   /* NUMBER OF HITS REPORTED */
    int elim;   /* RANGE ELIMINATED        */
    int last;   /* DISTANCE OF LAST HIT    */
    } local_oursearchdata;

/* newVPTREEOUR : CREATE NEW TASK FOR PRIORITY QUEUE STORAGE.
*/
static VPTREEOUR *newVPTREEOUR(VPTREE *vpt, int start, int length,
                             short depth, int infront, int potential){
    register VPTREEOUR *our = NEW(VPTREEOUR);
    our->start = start;
    our->length = length;
    our->depth = depth;
    our->infront = infront;
    our->potential  = potential;
    our->distance = 0; /* NECESSARY HERE ??? */
#ifdef VPTREE_PERFORMANCE_DATA
    our->obsvcomps = local_oursearchdata.comps;
#endif /* VPTREE_PERFORMANCE_DATA */
    return our;
    }

static void oursearchinterrupt(int id){
    printf("Search stopped after %d comps and %d skips\n"
           "Found %d sequences within distance %d\n",
           local_oursearchdata.comps, local_oursearchdata.skips,
           local_oursearchdata.hits,  local_oursearchdata.elim);
    if(local_oursearchdata.sig)
        local_oursearchdata.sig(id);
    else
        exit(0);
    return;
    }

/* ourcompVPTREEobserved : COMPARISON FUNCTION FOR observed PQUEUE.
                           PRIORITY: POTENTIAL >> INFRONT >> DEPTH
*/
static BOOLEAN ourcompVPTREEobserved(void *low, void *high){
    register VPTREEOUR *a = low, *b = high;
    if(a->potential == b->potential)
        return (a->infront == b->infront)
              ?(a->depth > b->depth)
              :(a->infront > b->infront);
    return (a->potential < b->potential);
    }

/*
TODO:
POTENTIAL >> DEPTH >> INFRONT
INFRONT IS ARBITARY.
LOWER RADIUS (node distance) SHOULD HAVE HIGHER PRIORITY.

*/

/* ourcompVPTREEused : COMPARISON FUNCTION FOR used PQUEUE.
                       PRIORITY SIMPLY ACCORDING LOWEST DISTANCE.
*/
static BOOLEAN ourcompVPTREEused(void *low, void *high){
    register VPTREEOUR *a = low, *b = high;
    return (a->distance < b->distance);
    }

/* showVPTREEhit : SHOW A SEQUENCE WHICH HAS BEEN HIT.
*/
static void showVPTREEhit(VPTREE *vpt, VPTREEOUR *r){
    register char *defn = getRAFASTAdef(vpt->dbdata,
                                   vpt->root[r->start].location);
    register int i;
    local_oursearchdata.hits++;
    if(local_oursearchdata.last == r->distance)
        i = printf("   = ");
    else {
        i = printf("%4d ", local_oursearchdata.hits);
        local_oursearchdata.last = r->distance;
        }
    i += printf("| %3d | %7d | %3d | ",
        r->distance, local_oursearchdata.comps,
        local_oursearchdata.skips);
#ifdef VPTREE_PERFORMANCE_DATA
     i+= printf("%4d | %4d | ", r->obsvcomps, r->usedcomps);
#endif /* VPTREE_PERFORMANCE_DATA */
    if(i < 79){
        defn[79-i] = '\0';
        puts(defn);
    } else putchar('\n');
    fflush(stdout);
    free(defn);
    return;
    }


/* oursearchVPTREE : MY OBSERVED-USED-REPORTED SEARCH STRATEGY.
                  ATTEMPTS TO REPORT HITS FOUND DURING SEARCH,
                  AS EARLY AS POSSIBLE.
*/
static void oursearchVPTREE(VPTREE *vpt, void *sc){
    register PQUEUE *observed = newPQUEUE(0, ourcompVPTREEobserved);
    register PQUEUE *used = newPQUEUE(0, ourcompVPTREEused);
    register VPTREEOUR *o, *no, *r;
    register VPTREENODE *n;
    register int i, p, infront, pp = 0;
    local_oursearchdata.sig  = signal(SIGINT, &oursearchinterrupt);
    local_oursearchdata.hits = 0;
    local_oursearchdata.comps = 0;
    local_oursearchdata.skips = 0;
    local_oursearchdata.elim = -1;
    local_oursearchdata.last = -1;
    o = newVPTREEOUR(vpt, 0, vpt->total, 0, TRUE, 0); /*TOTAL SPACE*/

    do {
        do {
            if(used->root) {
                r = (VPTREEOUR*)used->root->data;
            } else
                break;
            if(r->distance <= o->potential){
                showVPTREEhit(vpt, r);
                local_oursearchdata.elim = r->distance;
                r = popPQUEUE(used);
                free(r);
            } else
                  break;
            } while(TRUE);

        o->distance = vpt->metric->dist(sc,
                      vpt->root[o->start].location, vpt->dbdata);
        local_oursearchdata.comps++;
#ifdef VPTREE_PERFORMANCE_DATA
        o->usedcomps = local_oursearchdata.comps;
#endif /* VPTREE_PERFORMANCE_DATA */
        pushPQUEUE(used, o);
        n = vpt->root+o->start;
        p = n->partition;

        if(p == VPT_NO_PARTITION){
   /* printf("Leaf Bucket Encountered: size=%d\n", o->length); */
            local_oursearchdata.skips+=(o->length-1);
            for(i = 1; i < o->length; i++){ /* BUG HERE */
                /* errmsg(ERROR_WARNING, "Expect error"); */
                no = newVPTREEOUR(vpt, o->start+i, VPT_NO_PARTITION,
                                       o->depth+1, TRUE,
                                       o->distance);
                no->distance = o->distance;
#ifdef VPTREE_PERFORMANCE_DATA
                o->usedcomps = local_oursearchdata.comps;
#endif /* VPTREE_PERFORMANCE_DATA */
                pushPQUEUE(used, no);
                }
        } else {
            switch(o->length){
                case 0: /* IGNORE */
                    break;
                case 1: /* LEAF NODE */
                    break;
                case 2:
                    no = newVPTREEOUR(vpt, o->start+1, 1,
                                     o->depth+1, TRUE, o->potential);
                    pushPQUEUE(observed, no);
                    break;
                default:
                    if(o->distance < n->distance){
                        infront = TRUE;
                        pp = n->distance - o->distance;
                        if(pp < o->potential+1)  /* INHERIT */
                            pp = o->potential+1; /* OPT */
                    } else {
                        infront = FALSE;
                        pp = o->distance - n->distance+1; /* OPT */
                        if(pp < o->potential) /* INHERIT */
                            pp = o->potential;
                        }
                    if(p > 1){
                        no = newVPTREEOUR(vpt, o->start+1, p-1,
                                           o->depth+1, infront,
                                           infront?o->potential:pp);
                        pushPQUEUE(observed, no);
                        }
                    no = newVPTREEOUR(vpt, o->start+p, o->length-p,
                                           o->depth+1, !infront,
                                           infront?pp:o->potential);
                    pushPQUEUE(observed, no);
                    break;
                }
            }
    } while((o = popPQUEUE(observed)));
    printf("Finished after [%d] comps, [%d] skips\n",
              local_oursearchdata.comps, local_oursearchdata.skips);
    freePQUEUE(observed);
    freePQUEUE(used);
    signal(SIGINT, local_oursearchdata.sig); /* RETURN TO PREV VALUE */
    return;
    }

/* ___________________________________________________________ */

void searchVPTREE(char *idxpath, char *query){
    register VPTREE *vpt = readVPTREE(idxpath);
    register void *sc;
    printf("RANK | DIST|  COMPS  |SKIPS"
#ifdef VPTREE_PERFORMANCE_DATA
           "| C@OB | C@US "
#endif /* VPTREE_PERFORMANCE_DATA */
           "| DESCRIPTION\n"
           "-----|-----|---------|-----"
#ifdef VPTREE_PERFORMANCE_DATA
           "|------|------"
#endif /* VPTREE_PERFORMANCE_DATA */
           "|------------\n");
    sc = vpt->metric->prep(query, vpt->metric->param);
    oursearchVPTREE(vpt, sc);
    vpt->metric->free(sc);
    freeVPTREE(vpt);
    return;
    }

void displayVPTREE(char *idxpath){
    register VPTREE *vpt = readVPTREE(idxpath);
    register char *defn;
    register int i;
    for(i = 0; i < vpt->total; i++){
        defn = getRAFASTAdef(vpt->dbdata,
                             vpt->root[i].location);
        printf("%5d %5d %5d %.55s\n",
           i+1, vpt->root[i].distance, vpt->root[i].partition, defn);
        free(defn);
        }
    freeVPTREE(vpt);
    return;
    }

/* |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| */
/* PARALLEL CODE FOR VPTREE BUILDING                              */

/* #define ALLOW_PARALLEL_IMPLENTATION */
#ifdef ALLOW_PARALLEL_IMPLENTATION
#include "parallel.h"

/* STRATEGY:
       o Use the SGI m_fork() library.
       o buildVPTREErecur() can be used for complete
         clustering of a subtask. And can be used in parallel.
       o Prep should always be called individually by each processor,
         as some metrics (eg FSMHIST mailboxes) do not parallelise.
       o define an alternative subtask processor,
         scoreVPTREEparallel().
       o Use three phases to the parallelisation.
            1: Serial clustering, but parallelised score
               generation.  When length below SCORETASK_LIMIT,
               place sub task in QUEUE:initial, else process.
               This task should proceed until there are atleast
               numprocs tasks in QUEUE:initial, but end
               soon as the qsort etc, is being done in serial.
            2: Parallel:For each element in QUEUE:initial,
               score, if subtask length exceeds SPLITTASK limit,
               put into QUEUE:initial, else into QUEUE:final
            3: Call buildVPTREErecur() with each small
               remaining subtask in the QUEUE:final
               for completion with out further parallelisation.
       o Not parallelised sequences counting, as is I/O limited.
*/

/* SCORETASK_LIMIT : DON'T USE PARALLEL SCORING IF TASK SIZE SMALLER
*/
#define VPTREE_PARALLEL_SCORETASK_LIMIT 1000


/* SCORETASK_SIZE : SIZE FOR A SUBTASK IN PARALLEL PHASE 1.
*/
#define VPTREE_PARALLEL_SCORETASK_SIZE 1000

/* SPLITTASK_LIMIT: IF SUBTASK IS BELOW THIS SIZE,
                      USE buildVPTREErecur (SERIAL).
*/
#define VPTREE_PARALLEL_SPLITTASK_LIMIT 1000

/* scoreVPTREE_PP1_section : PARALLEL SCORING FUNCTION.
*/
static void scoreVPTREE_PP1_section(VPTREE *vpt,
                   int start, int length, int numprocs, char *seq){
    register void *sc = vpt->metric->prep(seq, vpt->metric->param);
    register int myid = m_get_myid();
    register int i, z, count, from, to;
    register VPTREENODE *root = vpt->root;
    register FILE *fp = fopen(vpt->dbpath, "r");
    register METRIC m; /* TIDY ALL THIS HERE WITH A copyVPT FUNC ?? */
    m.dist = vpt->metric->dist;
    m.free = vpt->metric->free;
    do {
       count = nextPARALLELcounter();
       from = count*VPTREE_PARALLEL_SCORETASK_SIZE+1;
       if(length < from)
           break;
       to =   from+VPTREE_PARALLEL_SCORETASK_SIZE-1;
       if(length < to)
           to = length;
       fprintf(stderr, "Proc [%d] Count[%d] from[%d]-to[%d]\n",
              m_get_myid(), count, from, to);
       for(i = start+from, z = start+to; i < z; i++)
           root[i].distance = m.dist(sc, root[i].location, fp);
    } while(TRUE);
    fprintf(stderr, "Proc [%d] finished\n", myid);
    m.free(sc);
    fclose(fp);
    return;
    }

/* scoreVPTREE_PP1 : ASSIGN DISTANCE SCORE FOR MEMBER UPTO length
                 SEQUENCES AWAY FROM start IN VPTREE vpt.
                 FOR PARALLEL PHASE 1.
*/
static void scoreVPTREE_PP1(VPTREE *vpt, int start, int length,
                            int numprocs){
    int len;
    register char *seq = getRAFASTAseq(vpt->dbdata,
                                  vpt->root[start].location, &len);
    m_fork(scoreVPTREE_PP1_section, vpt, start, length, numprocs, seq);
    free(seq);
    errmsg(ERROR_WARNING, "Finished parallel score of [%d]", length);
    return;
    }

typedef struct {
    int start;
    int len;
    int id;
    } VPTREETASK;

/* buildVPTREErecur : RECURSIVE FUNCTION FOR VPTREE BUILDING.
*/
static LIST *buildVPTREErecurPP1(VPTREE *vpt, int start, int length,
             int *comps, int *skips, int *id, int numprocs){
    register VPTREENODE *n = vpt->root+start; /* SUBTASK ROOT */
    register int p;
    register LIST *leftover = newLIST();
    register VPTREETASK *task;
    switch(length){
        case 0:  break; /* SHOULD NEVER HAPPEN HERE */
        case 1:
            n->partition = 0;
            n->distance  = 0;
            break;
        case 2: /* AVOID SELECTION IF WIDTH IS 2 */
            selectVPTREE(vpt, start, 2);
            task = NEW(VPTREETASK); /* STORE FOR LATER */
            task->start = start;
            task->len   = 2;
            task->id    = (*id)++;
            queueLIST(leftover, task);
            (*comps)++; /* NO NEED FOR SORT */
            if(n->distance == 0){
                n->partition = VPT_NO_PARTITION;
                errmsg(ERROR_WARNING, "Found identical PAIR");
                (*skips)++;
            } else {
                n->partition = 1;
                buildVPTREErecurPP1(vpt, start+1, 1,
                                    comps, skips, id, numprocs);
                }
            break;
        default:
            selectVPTREE(vpt, start, length);
            if(length > VPTREE_PARALLEL_SCORETASK_LIMIT)
                scoreVPTREE_PP1(vpt, start, length, numprocs);
            else {
                task = NEW(VPTREETASK);
                task->start = start;
                task->len   = length;
                task->id    = (*id)++;
                queueLIST(leftover, task);
                }
            (*comps)+=length;
            qsort(n+1, length-1, sizeof(VPTREENODE),
                 (int(*)(const void *, const void *))compVPTREE);
            /* NEED TO REPLACE WITH A MERGESORT ROUTINE HERE */
            /* CANNOT USE MEDIAN ALGORITHM AS NON-PERFECT BALANCING */
            n->partition = partitionVPTREE(vpt, start, length);
            if(n->partition == VPT_NO_PARTITION){
                errmsg(ERROR_WARNING, "Found identical GROUP of %d",
                       length);
                (*skips)+=length;
            } else {
                p = n->partition;
                n->distance = n[p].distance;
                if(p > 1) /* FIX PARTITION=1 BUG */
                    buildVPTREErecurPP1(vpt, start+1, p-1,
                                        comps, skips, id, numprocs);
                buildVPTREErecurPP1(vpt, start+p, length-p,
                                    comps, skips, id, numprocs);
                }
            break;
        }
    return leftover;
    }

static void buildVPTREE_PP2(VPTREE *vpt, LIST *todo,
                            int numprocs, int total){
    register int id;
    register VPTREETASK *vt;
    while((id = nextPARALLELcounter()) < total){
        printf("proc %d, id %d\n", m_get_myid(), id);
        /* do { vt = todo->a->n->v; */
        /* } while(vt->id != id);  */
        /* vt = popLIST(todo);  */
        /* printf("proc %d got %d (%d) [from %d for %d]\n",  */
                /* m_get_myid(), vt->id,  id, vt->start, vt->len); */
        /* free(vt); */
        }
    return;
    }


/* buildVPTREEparallel : BUILD A NEW VPTREE IN PARALLEL.
*/
static void buildVPTREEparallel(VPTREE *vpt, int numprocs){
    int comps = 0, skips = 0, id = 0;
    register LIST *todo;
    errmsg(ERROR_INFO, "Starting Phase 1 parallel build");
    todo = buildVPTREErecurPP1(vpt, 0, vpt->total,
                        &comps, &skips, &id, numprocs);
    errmsg(ERROR_WARNING, "Phase 2 has %d subtasks", id);
    m_fork(buildVPTREE_PP2, vpt, todo, numprocs, id);

    errmsg(ERROR_INFO, "Built vptree of %d elements "
                       "with %d comparisons and %d skips",
                       vpt->total, comps, skips);
    return;
    }

/* makeVPTREEparallel : BUILD A VPTREE FROM dbpath
                IN PARALLEL WITH numprocs
                USING metric AND WRITE IT OUT TO idxpath.
*/
void makeVPTREEparallel(char *dbpath, char *idxpath,
                        METRIC *metric, int numprocs){
    VPTREE *vpt = newVPTREE(dbpath, metric);
    if(numprocs == 1)
        buildVPTREE(vpt); /* SERIAL IMPLEMENTATION */
    else
        buildVPTREEparallel(vpt, numprocs);
    writeVPTREE(vpt, idxpath);
    freeVPTREE(vpt);
    return;
    }

#endif /* ALLOW_PARALLEL_IMPLENTATION */

/* ___________________________________________________________ */

/* |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| */


#ifdef TEST_THIS_MODULE
/* ### TEST CODE ##################################### */

int main(){
    /* makeVPTREE("small.fasta", "a.vpt",  */
              /* newMETRIC(METRIC_EDITDIST, newMETRICPARAMwd(4)) ); */
    makeVPTREE("/home/guy/checkout/estate/example/estsmall.fasta",
              "test.vpt", newMETRIC(METRIC_EDITDIST, NULL ) );
    return 0;
    }

/* ### TEST CODE ##################################### */
#endif /* TEST_THIS_MODULE */

#endif /* INCLUDED_VPTREE_C */

/*  TODO:
    Collaborative subspace elimination.
    Remove passing of the numprocs variable in parallel routines.

    Add a build progress indicatior.
    Try the near-metric, D(a,b) = S(a,a)+S(b,b)-2S(a,b).
    Add a decent vantage point picking algorithm.
*/

