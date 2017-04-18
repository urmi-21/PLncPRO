# ------------------------------------------------------------------
# Makefile for ESTate source code.
# ------------------------------------------------------------------

CC        = cc
# CC        = gcc

LDFLAGS   = -O
# LDFLAGS   = -g -Wall -Werror -O
# LDFLAGS   = -O3 -fomit-frame-pointer

# IT IS BEST TO TURN ON COMPILER OPTIMIZATION, IT IS REQUIRED
# TO MAKE THE DYNAMIC PROGRAMMING ROUTINES WORK EFFICIENTLY.

# THE LOCATION WHERE BINARIES WILL BE PLACED AFTER A make install
BINDIR    = ../bin/

# YOU SHOULD NOT NEED TO EDIT ANYTHING BELOW THIS LINE
# ----------------------------------------------------

TESTFLAGS = -DTEST_THIS_MODULE -DDEBUG
LOADLIBS  =

CHECK     = $(GENERAL)common.h
TESTCC    = $(CC) $(LDFLAGS) $(TESTFLAGS) -o $(@F)
BUILDCC   = $(CC) $(LDFLAGS) -o $(@F)

.SUFFIXES: .o .c .test

# ------------------------------------------------------------------

GENERAL   = general/
PROGRAM   = program/
STRUCT    = struct/
PARSE     = parse/
DYNAMIC   = dynamic/
SEQUENCE  = sequence/
CLUSTER   = cluster/
METRIC    = metric/

GENERAL_OBJ = $(GENERAL)arg.o     $(GENERAL)config.o     \
              $(GENERAL)dict.o    $(GENERAL)error.o      \
              $(GENERAL)string.o  $(GENERAL)stringlist.o

# ------------------------------------------------------------------

all:	 program
everything: all test

general: $(GENERAL_OBJ)

general_test: arg.test config.test dict.test \
              error.test string.test stringlist.test
dynamic_test: affine.test submat.test editdist.test
struct_test: fsm.test list.test pqueue.test vfsm.test
parse_test: ioutil.test parsedbest.test parser.test \
            rafasta.test readfasta.test parseseqdb.test
sequence_test:	sequtil.test ntaa.test
cluster_test:	clusterset.test avastore.test wordgraph.test
metric_test:	metric.test vptree.test

test:		general_test dynamic_test struct_test parse_test \
                sequence_test cluster_test metric_test

PROGRAM_LIST=	alignsearch calcwordprob estcluster fasta2usage \
                flat2coding framefinder precluster revcomp      \
                translate vpbuild vpsearch wordsearch
program:	$(PROGRAM_LIST)

install:	all
		mkdir -p $(BINDIR); cp $(PROGRAM_LIST) $(BINDIR)

clean:;	touch a.o $(GENERAL)a.o a.test $(PROGRAM_LIST); \
        rm *.o */*.o *.test $(PROGRAM_LIST)

# ------------------------------------------------------------------
# DEFINE SUFFIX RULES

.c.o:		$(CHECK)
		$(CC) $(LDFLAGS) -c $< -o $@

.c.test:	$(CHECK)
		$(TESTCC) $<

.o:		$(CHECK)
		$(BUILDCC) $<

# ------------------------------------------------------------------
# RULES FOR BUILDING PROGRAMS


CALCWORDPROB_OBJ = $(PROGRAM)calcwordprob.c $(GENERAL_OBJ) \
                  $(STRUCT)vfsm.o $(PARSE)parser.o \
                  $(SEQUENCE)sequtil.o
calcwordprob: 	$(CALCWORDPROB_OBJ)
		$(BUILDCC) $(CALCWORDPROB_OBJ) -lm

ALIGNSEARCH_OBJ = $(PROGRAM)alignsearch.c $(PARSE)readfasta.o \
                 $(PARSE)rafasta.o \
                 $(DYNAMIC)submat.o \
                 $(DYNAMIC)affine.o $(SEQUENCE)sequtil.o \
                 $(PARSE)ioutil.o $(STRUCT)pqueue.o \
                 $(GENERAL_OBJ)
alignsearch: 	$(ALIGNSEARCH_OBJ)
		$(BUILDCC) $(ALIGNSEARCH_OBJ)

ESTCLUSTER_OBJ = $(PROGRAM)estcluster.c $(GENERAL_OBJ) \
              $(CLUSTER)avastore.o $(STRUCT)pqueue.o \
              $(PARSE)rafasta.o $(PARSE)ioutil.o \
              $(DYNAMIC)affine.o $(DYNAMIC)submat.o \
              $(CLUSTER)clusterset.o $(SEQUENCE)sequtil.o \
              $(PARSE)readfasta.o $(STRUCT)list.o
estcluster: 	$(ESTCLUSTER_OBJ)
		$(BUILDCC) $(ESTCLUSTER_OBJ)

FRAMEFINDER_OBJ = $(PROGRAM)framefinder.c $(GENERAL_OBJ) \
                  $(PARSE)readfasta.o $(SEQUENCE)ntaa.o \
                  $(PARSE)parser.o $(STRUCT)vfsm.o \
                  $(SEQUENCE)sequtil.o
framefinder: 	$(FRAMEFINDER_OBJ)
		$(BUILDCC) $(FRAMEFINDER_OBJ) -lm

FASTA2USAGE_OBJ = $(PROGRAM)fasta2usage.c $(GENERAL_OBJ) \
                  $(STRUCT)vfsm.o $(SEQUENCE)sequtil.o \
                  $(PARSE)readfasta.o
fasta2usage: 	$(FASTA2USAGE_OBJ)
		$(BUILDCC) $(FASTA2USAGE_OBJ) -lm

FLAT2CODING_OBJ= $(PROGRAM)flat2coding.c $(GENERAL_OBJ) \
                 $(PARSE)parseseqdb.o $(PARSE)parser.o \
                 $(SEQUENCE)sequtil.o
flat2coding: 	$(CHECK) $(FLAT2CODING_OBJ)
		$(BUILDCC) $(FLAT2CODING_OBJ)

REVCOMP_OBJ = $(PROGRAM)revcomp.c $(GENERAL_OBJ) \
              $(PARSE)readfasta.c $(SEQUENCE)sequtil.o
revcomp:	$(REVCOMP_OBJ)
	 	$(BUILDCC) $(REVCOMP_OBJ)

TRANSLATE_OBJ = $(PROGRAM)translate.c $(GENERAL_OBJ) \
                $(PARSE)readfasta.c $(SEQUENCE)ntaa.o  \
                $(SEQUENCE)sequtil.o
translate:	$(TRANSLATE_OBJ)
	 	$(BUILDCC) $(TRANSLATE_OBJ)

PRECLUSTER_OBJ = $(PROGRAM)precluster.c $(GENERAL_OBJ) \
                $(CLUSTER)wordgraph.c $(CLUSTER)avastore.o  \
                $(PARSE)ioutil.o $(PARSE)rafasta.o \
                $(STRUCT)fsm.o $(STRUCT)vfsm.o     \
                $(STRUCT)list.o $(STRUCT)pqueue.o  \
                $(SEQUENCE)sequtil.o
precluster:	$(PRECLUSTER_OBJ)
	 	$(BUILDCC) $(PRECLUSTER_OBJ) -lm

VPBUILD_OBJ= $(PROGRAM)vpbuild.c $(GENERAL_OBJ) \
             $(METRIC)vptree.c $(METRIC)metric.o \
             $(DYNAMIC)editdist.o $(STRUCT)fsm.o \
             $(PARSE)rafasta.o $(STRUCT)pqueue.o \
             $(STRUCT)list.o $(PARSE)ioutil.o
vpbuild:	$(VPBUILD_OBJ)
		$(BUILDCC) $(VPBUILD_OBJ)

VPSEARCH_OBJ = $(PROGRAM)vpsearch.c $(GENERAL_OBJ) \
             $(METRIC)vptree.c $(METRIC)metric.o \
             $(DYNAMIC)editdist.o $(STRUCT)fsm.o \
             $(PARSE)rafasta.o $(STRUCT)pqueue.o \
             $(STRUCT)list.o $(PARSE)ioutil.o    \
             $(SEQUENCE)sequtil.o
vpsearch:	$(VPSEARCH_OBJ)
		$(BUILDCC) $(VPSEARCH_OBJ)

WORDSEARCH_OBJ = $(PROGRAM)wordsearch.c $(GENERAL_OBJ) \
		 $(STRUCT)fsm.o $(PARSE)ioutil.o       \
                 $(PARSE)parser.o $(PARSE)readfasta.o
wordsearch: 	$(WORDSEARCH_OBJ)	
		$(BUILDCC) $(WORDSEARCH_OBJ)

# ------------------------------------------------------------------
# RULES FOR GENERAL TEST OBJECTS

error.test:		$(GENERAL)error.test
string.test:		$(GENERAL)string.test
stringlist.test:	$(GENERAL)stringlist.test
dict.test:		$(GENERAL)dict.test

ARG_TEST_OBJ =  $(GENERAL)arg.c \
                $(GENERAL)config.o     \
                $(GENERAL)dict.o   $(GENERAL)error.o   \
                $(GENERAL)string.o $(GENERAL)stringlist.o
arg.test:	$(ARG_TEST_OBJ)
		$(TESTCC) $(ARG_TEST_OBJ)

CONFIG_TEST_OBJ = $(GENERAL)config.c \
                  $(GENERAL)dict.o \
                  $(GENERAL)error.o $(GENERAL)string.o     \
                  $(GENERAL)stringlist.o
config.test:	$(CONFIG_TEST_OBJ)
		$(TESTCC) $(CONFIG_TEST_OBJ)

# ------------------------------------------------------------------
# RULES FOR STRUCT TEST OBJECTS

fsm.test:	$(STRUCT)fsm.test
list.test:	$(STRUCT)list.test

PQUEUE_TEST_OBJ = $(STRUCT)pqueue.c $(GENERAL)error.o
pqueue.test:	$(PQUEUE_TEST_OBJ)	
	 	$(TESTCC) $(PQUEUE_TEST_OBJ)	

VFSM_TEST_OBJ = $(STRUCT)vfsm.c $(GENERAL)error.o \
                $(SEQUENCE)sequtil.o
vfsm.test:	$(VFSM_TEST_OBJ)	
	 	$(TESTCC) $(VFSM_TEST_OBJ) -lm

# ------------------------------------------------------------------
# RULES FOR PARSE TEST OBJECTS

ioutil.test:	$(PARSE)ioutil.test
parser.test:	$(PARSE)parser.test

PARSEDBEST_TEST_OBJ = $(PARSE)parsedbest.c $(GENERAL)error.o \
                      $(STRUCT)fsm.o $(STRUCT)list.o \
                      $(GENERAL)string.o
parsedbest.test:	$(PARSEDBEST_TEST_OBJ)
			$(TESTCC) $(PARSEDBEST_TEST_OBJ)

READFASTA_TEST_OBJ = 	$(PARSE)readfasta.c $(GENERAL)error.o
readfasta.test:	$(READFASTA_TEST_OBJ)
		$(TESTCC) $(READFASTA_TEST_OBJ)

RAFASTA_TEST_OBJ =	$(PARSE)rafasta.c $(GENERAL)error.o \
                        $(PARSE)ioutil.o
rafasta.test:	$(RAFASTA_TEST_OBJ)
		$(TESTCC) $(RAFASTA_TEST_OBJ)

PARSESEQDB_TEST_OBJ =	$(PARSE)parseseqdb.c $(GENERAL)string.o \
                        $(GENERAL)stringlist.o $(PARSE)ioutil.o \
                        $(PARSE)parser.o
parseseqdb.test:	$(PARSESEQDB_TEST_OBJ)
		$(TESTCC) $(PARSESEQDB_TEST_OBJ)

# ------------------------------------------------------------------
# RULES FOR DYNAMIC TEST OBJECTS

SUBMAT_TEST_OBJ = $(DYNAMIC)submat.c $(GENERAL)error.o
submat.test:	$(SUBMAT_TEST_OBJ)	
	 	$(TESTCC) $(SUBMAT_TEST_OBJ)	

AFFINE_TEST_OBJ = $(DYNAMIC)affine.c \
                    $(GENERAL)error.o $(DYNAMIC)submat.o
affine.test:	$(AFFINE_TEST_OBJ)	
	 	$(TESTCC) $(AFFINE_TEST_OBJ)	

editdist.test:	$(DYNAMIC)editdist.test

# ------------------------------------------------------------------
# RULES FOR SEQUENCE TEST OBJECTS

SEQUTIL_TEST_OBJ = $(SEQUENCE)sequtil.c
sequtil.test:	$(SEQUTIL_TEST_OBJ)
	 	$(TESTCC) $(SEQUTIL_TEST_OBJ)

NTAA_TEST_OBJ = $(SEQUENCE)ntaa.c $(GENERAL)error.o
ntaa.test:	$(NTAA_TEST_OBJ)
	 	$(TESTCC) $(NTAA_TEST_OBJ)

# ------------------------------------------------------------------
# RULES FOR CLUSTER TEST OBJECTS

AVASTORE_TEST_OBJ = $(CLUSTER)avastore.c $(STRUCT)pqueue.o \
                    $(GENERAL)error.o
avastore.test:	$(AVASTORE_TEST_OBJ)
	 	$(TESTCC) $(AVASTORE_TEST_OBJ)

CLUSTERSET_TEST_OBJ = $(CLUSTER)clusterset.c $(PARSE)ioutil.o
clusterset.test:	$(CLUSTERSET_TEST_OBJ)
	 	$(TESTCC) $(CLUSTERSET_TEST_OBJ)

WORDGRAPH_TEST_OBJ = $(CLUSTER)wordgraph.c $(GENERAL)error.o \
                $(PARSE)ioutil.o $(PARSE)rafasta.o \
                $(STRUCT)fsm.o $(STRUCT)vfsm.o \
                $(STRUCT)list.o $(STRUCT)pqueue.o \
                $(CLUSTER)avastore.o $(SEQUENCE)sequtil.o
wordgraph.test:	$(WORDGRAPH_TEST_OBJ)
	 	$(TESTCC) $(WORDGRAPH_TEST_OBJ) -lm

# ------------------------------------------------------------------
# RULES FOR METRIC TEST OBJECTS

METRIC_TEST_OBJ = $(METRIC)metric.c $(GENERAL)error.o \
                  $(DYNAMIC)editdist.o $(STRUCT)fsm.o
metric.test:	$(METRIC_TEST_OBJ)
		$(TESTCC) $(METRIC_TEST_OBJ)

VPTREE_TEST_OBJ = $(METRIC)vptree.c $(METRIC)metric.o \
                  $(GENERAL)error.o \
                  $(DYNAMIC)editdist.o $(STRUCT)fsm.o \
                  $(PARSE)rafasta.o $(STRUCT)pqueue.o \
                  $(STRUCT)list.o $(PARSE)ioutil.o
vptree.test:	$(VPTREE_TEST_OBJ)
		$(TESTCC) $(VPTREE_TEST_OBJ)

# ------------------------------------------------------------------

