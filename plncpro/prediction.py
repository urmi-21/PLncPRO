'''
PLNCPRO
#################################################################################
#  This file automatically extracts features for given fasta file of sequences    #
#  It also take a model name as argument and then predicts the classes of sequences
#  in the input file.                                #
#  We label lncrna as 0 and pcrna as 1                        #
#  Please give same parameters for blast; those given at model building time     #
#################################################################################
Author : Urminder Singh
email: urmind13_sit@jnu.ac.in
UrMi 21/4/16
''' 
import sys, getopt
import os
from Bio import SeqIO

###########Define Funcs###############
def printhelp():
    print ("*********************************************Help*********************************************")
    print ("DESCRIPTION")
    print ("This script classifies transcripts as coding or non coding transcripts")
    print ("Arguments:")
    print ("-h     print this message")
    print ("-p     output file name to store prediction results")
    print ("-i     path to file containing input sequences")
    print ("-m     path to the model file")
    print ("-o     output directory name to store all results")
    print ("-d     path to blast database")
    print ("                   OPTIONAL")
    print ("-t     number of threads[default: 4]")
    print ("-l     path to the files containg labels(this outputs performance of the classifier)")
    print ("-r     clean up intermediate files")
    print ("-v     show more messages")    
    print ("--min_len     specifiy min_length to filter input files")
    print ("--noblast     Don't use blast features")
    print ("--no_ff     Don't use framefinder features")
    print ("--qcov_hsp     specify qcov parameter for blast[default:30]")
    print ("--blastres     path to blast output for input file")
    
    
def main(args = sys.argv,home=None):
    
    
    ######################################
    ############Define variables##############
    in_flag=False
    out_flag=False
    db_flag=False
    model_flag=False
    removefiles_flag=False
    noblast_flag=False
    noblast_flag_val="false"
    blastres_flag=False
    no_ff_flag=False
    no_ff_flagval="false"
    min_len_flag=False
    out_file_flag=False
    ##set default values
    min_length=1
    num_threads=4
    max_qcov_hsp=30
    max_targets=10
    outdir=os.getcwd()
    out_dir_flag=False
    in_file=""
    out_file=""
    model_file=""
    blastdb=""
    blastres_file=""
    vflag=False
    lflag=False
    lflag_val="false"
    label_file=""
    #####################################
    try:
        opts, args = getopt.getopt(sys.argv[2:],"ht:i:rm:p:o:d:vl:",["prediction_out=","infile=","noblast","no_ff","qcov_hsp=","threads=","db=","outdir=","model=","remove_temp","blastres=","labels=","min_len="])
    except getopt.GetoptError:
        printhelp()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            printhelp()
            sys.exit()
        elif opt in ("-p", "--prediction_out"):
            #print 'pos found'
            out_file=os.path.abspath(arg)
            outfile=arg
            out_file_flag=True
            #print (out_file)
            #print arg
        elif opt in ("-i", "--infile"):
            #print 'neg found'
            in_file=os.path.abspath(arg)
            #print in_file
        elif opt in ("-m", "--model"):
            #print 'model found'
            model_file=arg
            #print model_file
        elif opt in ("-d", "--db"):
            #print 'db found'
            blastdb=(arg)
            db_flag=True
            #print blastdb
        elif opt in ("-t", "--threads"):
            #print 'threads found'
            num_threads=int(arg)
            #print num_threads
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
        elif opt in ("--blastres"):
            #print 'blastres found'
            blastres_flag=True
            blastres_file=arg
            #print arg
        elif opt in ("-o","--outdir"):
            #print 'outdir found'
            outdir=os.getcwd()+'/'+arg
            outdir=os.path.abspath(arg)
            out_dir_flag=True
        elif opt in ("-v","--verbose"):
            #print 'v found'
            vflag=True
        elif opt in ("-l","--labels"):
            #print 'l found'
            lflag=True
            lflag_val="true"
            label_file=arg
        elif opt in ("--min_len"):
            #print 'ml found'
            min_length=int(arg)
            min_len_flag=True
    

    ###Check all necessary inputs

    if out_file_flag==False:
        print ("please give output file name for predictions -p out_file")
        sys.exit(0)
    if out_dir_flag==False:
        print ("please give output directory name -o out_dir")
        sys.exit(0)

    if model_file=="":
        print ("please specify model file -m model")
        sys.exit(0)
        #check if model exists

    if not (os.path.isfile(model_file) ):
        print(('Please check model file...Error file:'+model_file+ ' doesn\'t exist'))
        sys.exit(0)

    #check pos,neg exists
    if not (os.path.isfile(in_file) ):
        print(('Please check input file...Error file:'+in_file+ ' doesn\'t exist'))
        sys.exit(0)

    if lflag==True:
        if not (os.path.isfile(label_file) ):
            print(('Please check label file...Error file:'+label_file+ ' doesn\'t exist'))
            sys.exit(0)

    ##check blast database
    if noblast_flag==False:
        if (blastres_flag==False):
            if db_flag==False:
                print ('Please specify blast database...Error')
                sys.exit(0)

##check for blastres files
    if blastres_flag==True:
        if not (os.path.isfile(blastres_file) ):
            print(('Please check blastres file...Error file: '+blastres_file+ ' doesn\'t exist'))
            sys.exit(0)

####################################################################################################################
##########################################START Reading FILES#######################################################

    ##remove sequences with min length
    if min_len_flag==True:
        print ('Removing short sequences........')
        print(('minlen:',min_length))
        output_handle = open(in_file+'_temp_'+str(min_length), "w")
        ctr=0
        short_ctr=0
        temp_rec=[]
        for record in SeqIO.parse(in_file, "fasta"):
            length=len(record.seq)
            #print length
            if length>=min_length:
                ctr=ctr+1
                temp_rec.append(record)
                #SeqIO.write(record, output_handle, "fasta")
                #print 'writing',record.id
            else:
                short_ctr=short_ctr+1
    
        SeqIO.write(temp_rec, output_handle, "fasta")
        output_handle.close()        #important else gives errors
    
        print(('Short Sequences <',str(min_length),str(short_ctr),' removed,',str(ctr),' retained'))
        in_file=in_file+'_temp_'+str(min_length)
        print ('New file with filtered sequences:')    
        print (in_file)


########################################################################################################################
############################################Read in file###############################################################
    if vflag:
        print ('Reading Input File...\nExtracting Features...')
    
    os.system("python "+home+"/bin/extractfeatures.py "+in_file+" 1")

    ########Run framefinder
    ff_featurefile_pos=in_file+"_ffout_framefinderfeatures"
    if no_ff_flag==True:
        print ('Skipping framefinder...')
        os.system("echo '' > "+in_file+"_ffout")

    else:
        #print "lib/framefinder/framefinder -r False -w lib/framefinder/framefinder.model "+in_file+" > "+in_file+"_ffout"
        os.system(home+"/lib/framefinder/framefinder -r False -w "+home+"/lib/framefinder/framefinder.model "+in_file+" > "+in_file+"_ffout")
        #parse framefinder results and write framefinder feature files
        if vflag:
            print ('Extracting Framefinder Features...')
        os.system("python "+home+"/bin/ffparse.py "+in_file+"_ffout")

    ######Run BLASTX########
    if noblast_flag==False:
        if blastres_flag==True:
            if vflag:
                print ('Parsing Blast results...\nExtracting Features...')
                
            os.system("python "+home+"/bin/blastparse_mt3.py "+blastres_file)
            
            if vflag:
                print ('Merging all Features...')
                
            os.system("python "+home+"/bin/mergefeatures.py "+in_file+"_features "+no_ff_flagval+" "+ff_featurefile_pos+" "+noblast_flag_val+" "+blastres_file+"_blastfeatures "+" "+in_file+"_all_features")
            
        else:
            #filename for blastres
            blastres_pos=in_file+"_blastres"
            if vflag:
                print ('Running BLASTX...This might take some time depending on your input.')
            bcommand="blastx -query "+in_file+" -db "+blastdb+" -outfmt '6 qseqid sseqid pident evalue qcovs qcovhsp score bitscore qframe sframe' -out "+blastres_pos+" -qcov_hsp_perc "+ str(max_qcov_hsp)+" -num_threads "+str(num_threads)
            print((str(bcommand)))
            os.system(str(bcommand))
            if vflag:
                print ('Parsing Blast Results...')

            os.system("python "+home+"/bin/blastparse_mt3.py "+blastres_pos)
            print ('Merging all Features...')
            os.system("python "+home+"/bin/mergefeatures.py "+in_file+"_features "+no_ff_flagval+" "+ff_featurefile_pos+" "+noblast_flag_val+" "+blastres_pos+"_blastfeatures "+" "+in_file+"_all_features")

    else:
        if vflag:
            print ('Skipping Blast...')
        blastres_pos=in_file+"_blastres"
        os.system("echo 'X    X    0    0    0    0    0    0    0    0' > "+blastres_pos)
        os.system("python "+home+"/bin/blastparse_mt3.py "+blastres_pos)
        print ('Merging all Features...')
        os.system("python "+home+"/bin/mergefeatures.py "+in_file+"_features "+no_ff_flagval+" "+ff_featurefile_pos+" "+noblast_flag_val+" "+blastres_pos+"_blastfeatures "+" "+in_file+"_all_features")    

#########################################################################################################################
##############################################Start prediction###########################################################
##mergeboth files in one
    if vflag:
        print ('Predicting...')
        
    print((str("python "+home+"/bin/rf/predict.py "+in_file+"_all_features "+model_file+" "+out_file+" "+lflag_val+" "+label_file)))
    os.system("python "+home+"/bin/rf/predict.py "+in_file+"_all_features "+model_file+" "+out_file+" "+lflag_val+" "+label_file)

################Remove Temp Files##################
    if removefiles_flag==True:
        print ('Removing temp files...')
        #print "rm -f "+in_file+"_ffout_framefinderfeatures"
        os.system("rm -f "+in_file+"_ffout_framefinderfeatures")
        os.system("rm -f "+in_file+"_blastres_blastfeatures")
        os.system("rm -f "+in_file+"_features")
        os.system("rm -f "+in_file+"_ffout")
        if min_len_flag==True:    
            os.system("rm -f "+in_file)

#########################################Move files to out_dir##########################################
    files_dir=os.path.dirname(os.path.realpath(in_file))
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    os.system("mv "+in_file+"_all_features"+" "+outdir+"/")
    os.system("mv "+in_file+"_blastres"+" "+outdir)
    os.system("mv "+out_file+" "+outdir)

    print(('All outputs saved to: '+ outdir))
    print ('END')
    
    
    
    
if __name__ == "__main__":
    main()
    
    