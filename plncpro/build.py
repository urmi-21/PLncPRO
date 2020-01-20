'''
PLNCPRO
This python scripts build predictive model from training set.
#################################################################################
#  This file automatically extracts features for given labeled input sequences    #
#  The input is given as command line in the following order            #
#       positiveexample.fa and negativeexamples.fa                #
#  The final output of this script is the feature file containg both pos and -ve#
#  features. We label lncrna as 0 and pcrna as 1                #
#  give lncrna.fa as first input and pcrna.fa as second                #
################################################################################# 

Author : Urminder Singh
email: urmind13_sit@jnu.ac.in
UrMi 21/4/16
'''

import sys, getopt
import os
from Bio import SeqIO

def printhelp():
                        
     



    print("*********************************************Help*********************************************")
    print("                   DESCRIPTION")
    print("This script generates classification model from codin and non coding transcripts")
    print("Arguments:")
    print("-h     print this message")
    print("-p,--pos     path to file containing protein coding examples")
    print("-n,--neg     path to file containing non coding examples")
    print("-m,--model     output model name")
    print("-o,--outdir     output directory name to store all results")    
    print("-d     path to blast database")
    print("                   OPTIONAL")
    print("-t     number of threads[default: 4]")
    print("-k     number of trees[default: 1000]")
    print("-r     clean up intermediate files")
    print("-v     show more messages")    
    print("--min_len     specifiy min_length to filter input files")
    print("--noblast     Don't use blast features")
    print("--no_ff     Don't use framefinder features")
    print("--qcov_hsp     specify qcov parameter for blast[default:30]")
    print("--pos_blastres     path to blast output for positive input file")
    print("--neg_blastres     path to blast output for negative input file")
    



def main(args = sys.argv,home=None):
    
    ######################################

    ############Define variables##############
    pos_flag=False
    neg_flag=False
    db_flag=False
    model_flag=False
    removefiles_flag=False
    noblast_flag=False
    noblast_flag_val="false"
    pos_blastres_flag=False
    neg_blastres_flag=False
    no_ff_flag=False
    no_ff_flagval="false"
    min_len_flag=False
    ##set default values
    min_length=1
    num_threads=4
    num_trees=1000
    max_qcov_hsp=30
    max_targets=10
    outdir=os.getcwd()
    out_dir_flag=False
    pos_file=""
    neg_file=""
    model_name=""
    blastdb=""
    pos_blastres_file=""
    neg_blastres_file=""
    vflag=False
    ###################################################################
    path_sep=os.pathsep
    ############################Set input options######################
    try:
        opts, args = getopt.getopt(sys.argv[2:],"ht:n:rm:p:o:d:k:v",["pos=","neg=","noblast","no_ff","qcov_hsp=","threads=","db=","outdir=","model=","num_trees=","remove_temp","pos_blastres=","neg_blastres=","min_len="])
    except getopt.GetoptError:
        printhelp()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            printhelp()
            sys.exit()
            
        elif opt in ("-p", "--pos"):
            #print 'pos found'
            pos_file=os.path.abspath(arg)
            #print (pos_file)
            #print arg
        elif opt in ("-n", "--neg"):
            #print 'neg found'
            neg_file=os.path.abspath(arg)
            #print neg_file
        elif opt in ("-m", "--model"):
            #print 'model found'
            model_name=arg
            #print model_name
        elif opt in ("-d", "--db"):
            #print 'db found'
            blastdb=(arg)
            db_flag=True
            #print blastdb
        elif opt in ("-t", "--threads"):
            #print 'threads found'
            num_threads=int(arg)
            #print num_threads
        elif opt in ("-k", "--num_trees"):
            #print 'ntree found'
            num_trees=int(arg)
            #print num_trees
        elif opt in ("-r", "--remove_temp"):
            #print 'remove found'
            removefiles_flag=True
            #print removefiles_flag
        elif opt in ("--noblast"):
            #print 'noblast found'
            noblast_flag=True
            noblast_flag_val="true"
            #print noblast_flag
        elif opt in ("--no_ff"):
            #print 'no_ff found'
            no_ff_flag=True
            no_ff_flagval="true"
            #print arg
        elif opt in ("--qcov_hsp"):
            #print 'qcov_hsp found'
            max_qcov_hsp=arg
            #print arg
        elif opt in ("--pos_blastres"):
            #print 'pos_blastres found'
            pos_blastres_flag=True
            pos_blastres_file=arg
            #print arg
        elif opt in ("--neg_blastres"):
            #print 'neg_blastres found'
            neg_blastres_flag=True
            neg_blastres_file=arg
            #print arg
        elif opt in ("-o","--outdir"):
            #print 'outdir found'
            outdir=os.getcwd()+'/'+arg
            outdir=os.path.abspath(arg)
            out_dir_flag=True
            #print outdir
            #sys.exit(0)
        elif opt in ("-v","--verbose"):
            #print 'v found'
            vflag=True
        elif opt in ("--min_len"):
            #print 'ml found'
            min_length=int(arg)
            min_len_flag=True

    ###Check all necessary inputs
    if out_dir_flag==False:
        print("please give output directory name -o out_dir")
        sys.exit(0)
    if model_name=="":
        print("please specify output model name -m model_name")
        sys.exit(0)

    #check if model exists
    model_file=outdir+'/'+model_name

    if (os.path.isfile(model_file) ):
        print(('Error... model file already exists: '+model_file+'\nExiting...'))
        sys.exit(0)
        
    #check pos,neg exists
    if not (os.path.isfile(pos_file) ):
        print(('Please check pos file...Error file:'+pos_file+ ' doesn\'t exist'))
        sys.exit(0)
        
    if not (os.path.isfile(neg_file) ):
        print(('Please check neg file...Error file:'+neg_file+' doesn\'t exist'))
        sys.exit(0)

##check blast database
    if noblast_flag==False:
        if (neg_blastres_flag==False or pos_blastres_flag==False):
            if db_flag==False:
                print('Please specify blast database...Error')
                sys.exit(0)
                
    ##check for blastres files
    if pos_blastres_flag==True:
        if not (os.path.isfile(pos_blastres_file) ):
            print(('Please check pos_blastres file...Error file: '+pos_blastres_file+ ' doesn\'t exist'))
            sys.exit(0)
            
            
    if neg_blastres_flag==True:
        if not (os.path.isfile(neg_blastres_file) ):
            print(('Please check neg_blastres file...Error file: '+neg_blastres_file+ ' doesn\'t exist'))
            sys.exit(0)

    if vflag:
        print(outdir)
        print(neg_file)
        print((os.path.dirname(os.path.realpath(neg_file))))


####################################################################################################################
##########################################START Reading FILES#######################################################

    ##remove sequences with min length
    if min_len_flag==True:
        print('Removing short sequences........')
        output_handle = open(neg_file+'_temp'+str(min_length), "w")
        ctr=0
        short_ctr=0
        for record in SeqIO.parse(neg_file, "fasta"):
            length=len(record)
            if length>=min_length:
                ctr=ctr+1
                SeqIO.write(record, output_handle, "fasta")
            else:
                short_ctr=short_ctr+1

        output_handle = open(pos_file+'_temp'+str(min_length), "w")
        for record in SeqIO.parse(pos_file, "fasta"):
            length=len(record)
            if length>=min_length:
                ctr=ctr+1
                SeqIO.write(record, output_handle, "fasta")
            else:
                short_ctr=short_ctr+1
        print(('Short Sequences <',str(min_length),str(short_ctr),' removed,',str(ctr),' retained'))
        neg_file=neg_file+'_temp'+str(min_length)
        pos_file=pos_file+'_temp'+str(min_length)
        print('New files with filtered sequences:')
        print(neg_file)
        print(pos_file)
    
############################################Read neg file############################################################
    if vflag:
        print('Reading Negative File...\nExtracting Features...')
    os.system("python "+home+"/bin/extractfeatures.py "+neg_file+" 0")

    ########Run framefinder
    ff_featurefile_neg=neg_file+"_ffout_framefinderfeatures"
    if no_ff_flag==True:
        print('Skipping framefinder...')
        os.system("echo '' > "+neg_file+"_ffout")

    else:
        os.system(home+"/lib/framefinder/framefinder -r False -w "+home+"/lib/framefinder/framefinder.model "+neg_file+" > "+neg_file+"_ffout")
        #parse framefinder results and write framefinder feature files
        if vflag:
            print('Extracting Framefinder Features...')
        os.system("python "+home+"/bin/ffparse.py "+neg_file+"_ffout")

    ######Run BLASTX########
    if noblast_flag==False:
        if neg_blastres_flag==True:
            if vflag:
                print('Parsing Negative Blast results...\nExtracting Features...')
            os.system("python "+home+"/bin/blastparse_mt3.py "+neg_blastres_file)
            if vflag:
                print('Merging all Features...')
            os.system("python "+home+"/bin/mergefeatures.py "+neg_file+"_features "+no_ff_flagval+" "+ff_featurefile_neg+" "+noblast_flag_val+" "+neg_blastres_file+"_blastfeatures "+" "+neg_file+"_all_features")
        else:
            #filename for blastres
            blastres_neg=neg_file+"_blastres"
            if vflag:
                print('Running BLASTX...This might take some time depending on your input.')
            bcommand="blastx -query "+neg_file+" -db "+blastdb+" -outfmt '6 qseqid sseqid pident evalue qcovs qcovhsp score bitscore qframe sframe' -out "+blastres_neg+" -qcov_hsp_perc "+ str(max_qcov_hsp)+" -num_threads "+str(num_threads)
            print((str(bcommand)))
            os.system(str(bcommand))
            if vflag:
                print('Parsing Blast Results...')
                print("python "+home+"/bin/blastparse_mt3.py "+blastres_neg)

            os.system("python "+home+"/bin/blastparse_mt3.py "+blastres_neg)
            print('Merging all Features...')
            print("python "+home+"/bin/mergefeatures.py "+neg_file+"_features "+no_ff_flagval+" "+ff_featurefile_neg+" "+noblast_flag_val+" "+blastres_neg+"_blastfeatures "+" "+neg_file+"_all_features")
            os.system("python "+home+"/bin/mergefeatures.py "+neg_file+"_features "+no_ff_flagval+" "+ff_featurefile_neg+" "+noblast_flag_val+" "+blastres_neg+"_blastfeatures "+" "+neg_file+"_all_features")

    else:
        if vflag:
            print('Skipping Blast...')
        blastres_neg=neg_file+"_blastres"
        os.system("echo 'X    X    0    0    0    0    0    0    0    0' > "+blastres_neg)
        os.system("python "+home+"/bin/blastparse_mt3.py "+blastres_neg)
        print('Merging all Features...')
        os.system("python "+home+"/bin/mergefeatures.py "+neg_file+"_features "+no_ff_flagval+" "+ff_featurefile_neg+" "+noblast_flag_val+" "+blastres_neg+"_blastfeatures "+" "+neg_file+"_all_features")    

########################################################################################################################
############################################Read pos file###############################################################
    if vflag:
        print('Reading Positive File...\nExtracting Features...')
    os.system("python "+home+"/bin/extractfeatures.py "+pos_file+" 1")

    ########Run framefinder
    ff_featurefile_pos=pos_file+"_ffout_framefinderfeatures"
    if no_ff_flag==True:
        print('Skipping framefinder...')
        os.system("echo '' > "+pos_file+"_ffout")
    else:
        os.system(home+"/lib/framefinder/framefinder -r False -w "+home+"/lib/framefinder/framefinder.model "+pos_file+" > "+pos_file+"_ffout")
        #parse framefinder results and write framefinder feature files
        if vflag:
            print('Extracting Framefinder Features...')
        os.system("python "+home+"/bin/ffparse.py "+pos_file+"_ffout")

    ######Run BLASTX########
    if noblast_flag==False:
        if pos_blastres_flag==True:
            if vflag:
                print('Parsing Blast results...\nExtracting Features...')
            os.system("python "+home+"/bin/blastparse_mt3.py "+pos_blastres_file)
            if vflag:
                print('Merging all Features...')
            os.system("python "+home+"/bin/mergefeatures.py "+pos_file+"_features "+no_ff_flagval+" "+ff_featurefile_pos+" "+noblast_flag_val+" "+pos_blastres_file+"_blastfeatures "+" "+pos_file+"_all_features")
        else:
            #filename for blastres
            blastres_pos=pos_file+"_blastres"
            if vflag:
                print('Running BLASTX...This might take some time depending on your input.')
            bcommand="blastx -query "+pos_file+" -db "+blastdb+" -outfmt '6 qseqid sseqid pident evalue qcovs qcovhsp score bitscore qframe sframe' -out "+blastres_pos+" -qcov_hsp_perc "+ str(max_qcov_hsp)+" -num_threads "+str(num_threads)
            print((str(bcommand)))
            os.system(str(bcommand))
            if vflag:
                print('Parsing Blast Results...')

            os.system("python "+home+"/bin/blastparse_mt3.py "+blastres_pos)
            print('Merging all Features...')
            os.system("python "+home+"/bin/mergefeatures.py "+pos_file+"_features "+no_ff_flagval+" "+ff_featurefile_pos+" "+noblast_flag_val+" "+blastres_pos+"_blastfeatures "+" "+pos_file+"_all_features")

    else:
        if vflag:
            print('Skipping Blast...')
        blastres_pos=pos_file+"_blastres"
        os.system("echo 'X    X    0    0    0    0    0    0    0    0' > "+blastres_pos)
        os.system("python "+home+"/bin/blastparse_mt3.py "+blastres_pos)
        print('Merging all Features...')
        os.system("python "+home+"/bin/mergefeatures.py "+pos_file+"_features "+no_ff_flagval+" "+ff_featurefile_pos+" "+noblast_flag_val+" "+blastres_pos+"_blastfeatures "+" "+pos_file+"_all_features")    

    #########################################################################################################################
    #########################################Merge all features into single file#############################################
    ##mergeboth files in one
    if vflag:
        print('Merging Pos and Neg features...')
    os.system("cat "+neg_file+"_all_features > "+neg_file+"_final_features")
    os.system("cat "+pos_file+"_all_features | tail -n+2 >> "+neg_file+"_final_features")
    #########################################################################################################################
    ##########################################Build Model####################################################################

    if vflag:
        print('Building Model...')
    os.system("python "+home+"/bin/rf/buildmodel.py "+neg_file+"_final_features "+model_name+" "+str(num_trees)+" "+str(num_threads))

    ################Remove Temp Files##################
    if removefiles_flag==True:
        print('Removing temp files...')
        #print "rm -f "+pos_file+"_ffout_framefinderfeatures"
        os.system("rm -f "+pos_file+"_ffout_framefinderfeatures")
        os.system("rm -f "+pos_file+"_blastres_blastfeatures")
        os.system("rm -f "+pos_file+"_features")
        os.system("rm -f "+pos_file+"_ffout")
    
        os.system("rm -f "+neg_file+"_ffout_framefinderfeatures")
        os.system("rm -f "+neg_file+"_blastres_blastfeatures")
        os.system("rm -f "+neg_file+"_features")
        os.system("rm -f "+neg_file+"_ffout")
        if min_len_flag==True:    
            os.system("rm -f "+neg_file)
            os.system("rm -f "+pos_file)    

    #########################################Move files to out_dir##########################################
    neg_files_dir=os.path.dirname(os.path.realpath(neg_file))
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    os.system("mv "+neg_file+"_all_features"+" "+outdir+"/")
    os.system("mv "+neg_file+"_blastres"+" "+outdir)
    os.system("mv "+neg_file+"_final_features"+" "+outdir)
    os.system("mv "+pos_file+"_all_features"+" "+outdir+"/")
    os.system("mv "+pos_file+"_blastres"+" "+outdir)

    if removefiles_flag==False:
        os.system("mv "+pos_file+"_blastres_blastfeatures"+" "+outdir)
        os.system("mv "+pos_file+"_features"+" "+outdir)
        os.system("mv "+pos_file+"_ffout"+" "+outdir)
        os.system("mv "+pos_file+"_ffout_framefinderfeatures"+" "+outdir)
        os.system("mv "+neg_file+"_blastres_blastfeatures"+" "+outdir)
        os.system("mv "+neg_file+"_features"+" "+outdir)
        os.system("mv "+neg_file+"_ffout"+" "+outdir)
        os.system("mv "+neg_file+"_ffout_framefinderfeatures"+" "+outdir)
    ##move model file
    os.system("mv "+neg_files_dir+"/"+model_name+" "+outdir)

    print(('All outputs saved to: '+ outdir))
    print('END')
    
if __name__ == "__main__":
    main()
