#!/bin/bash

#Steps
#1 Download data from: 
#2 Setup a blastdb and replace ~/uniprotdb/prot.db with your blastdb
#3 execute the required commands to build prediction models
#NOTE: make sure all relative paths to files and data are correct

plncpro build -p plncpro_data/plant_new_fasta/monocot/train/monocot_pct_train.fa -n plncpro_data/plant_new_fasta/monocot/train/monocot_lnct_train.fa -o monocot_model -m monocot.model -d ~/uniprotdb/prot.db -t 6

plncpro build -p plncpro_data/plant_new_fasta/amt/train/amt_pct_train.fa -n plncpro_data/plant_new_fasta/amt/train/amt_lnct_train.fa -o amt_model -m amt.model -d ~/uniprotdb/prot.db -t 6

plncpro build -p plncpro_data/plant_new_fasta/at/train/at_pct_train.fa -n plncpro_data/plant_new_fasta/at/train/at_lnct_train.fa -o at_model -m at.model -d ~/uniprotdb/prot.db -t 6

plncpro build -p plncpro_data/plant_new_fasta/cr/train/cr_pct_train.fa -n plncpro_data/plant_new_fasta/cr/train/cr_lnct_train.fa -o cr_model -m cr.model -d ~/uniprotdb/prot.db -t 6


plncpro build -p plncpro_data/plant_new_fasta/gm/train/gm_pct_train.fa -n plncpro_data/plant_new_fasta/gm/train/gm_lnct_train.fa -o gm_model -m gm.model -d ~/uniprotdb/prot.db -t 6


plncpro build -p plncpro_data/plant_new_fasta/os/train/os_pct_train.fa -n plncpro_data/plant_new_fasta/os/train/os_lnct_train.fa -o os_model -m os.model -d ~/uniprotdb/prot.db -t 6


plncpro build -p plncpro_data/plant_new_fasta/pp/train/pp_pct_train.fa -n plncpro_data/plant_new_fasta/pp/train/pp_lnct_train.fa -o pp_model -m pp.model -d ~/uniprotdb/prot.db -t 6


plncpro build -p plncpro_data/plant_new_fasta/sm/train/sm_pct_train.fa -n plncpro_data/plant_new_fasta/sm/train/sm_lnct_train.fa -o sm_model -m sm.model -d ~/uniprotdb/prot.db -t 6


plncpro build -p plncpro_data/plant_new_fasta/st/train/st_pct_train.fa -n plncpro_data/plant_new_fasta/st/train/st_lnct_train.fa -o st_model -m st.model -d ~/uniprotdb/prot.db -t 6



plncpro build -p plncpro_data/plant_new_fasta/vv/train/vv_pct_train.fa -n plncpro_data/plant_new_fasta/vv/train/vv_lnct_train.fa -o vv_model -m vv.model -d ~/uniprotdb/prot.db -t 6


plncpro build -p plncpro_data/plant_new_fasta/zm/train/zm_pct_train.fa -n plncpro_data/plant_new_fasta/zm/train/zm_lnct_train.fa -o zm_model -m zm.model -d ~/uniprotdb/prot.db -t 6


plncpro build -p plncpro_data/hg24/train/hg24_pct_train_5000.fa -n plncpro_data/hg24/train/hg24_lnct_train_5000.fa -o hg_model -m hg.model -d ~/uniprotdb/prot.db -t 6


plncpro build -p plncpro_data/mm8/train/m8_pct_train_2500.fa -n plncpro_data/mm8/train/m8_lnct_train_2500.fa -o mm_model -m mm.model -d ~/uniprotdb/prot.db -t 6
