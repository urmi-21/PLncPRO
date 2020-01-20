'''
PLNCPRO
This file takes as input
and it extracts dasta seq by label after reading prediction file
Author : Urminder Singh
email: urmind13_sit@jnu.ac.in
UrMi 3/5/16
'''

import sys
import getopt
#import math
#import re
from Bio import SeqIO
from Bio.Seq import Seq

def main(args = sys.argv,home=None):
	#set defaults
	cutoff=0
	min_len=0
	max_len=float('Inf')
	label='0'

	try:
		opts, args = getopt.getopt(sys.argv[2:],"hf:o:p:l:s:r:",["ifile=","ofile=","min=","max="])
		#print opts
	except getopt.GetoptError:
		print('predstoseq.py -f <input fastafile> -o <outputfile> -p <predictionfile> -l <required label default:0> -s <class_prob_cutoff, default:0> -m <min_length, default:0> <max_length, default:inf>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('Use this to extract lncRNA or mRNA sequences, as predicted by PLNCPRO, from the input fasta file')
			print('Usage:')
			print('predstoseq.py -f <input fastafile> -o <outputfile> -p <predictionfile> -l <required label default:0> -s <class_prob_cutoff[range 0-1], default:0> -m <min_length, default:0> <max_length, default:inf>')
			sys.exit()
		elif opt in ("-f", "--ifile"):
			#print 'infile found'
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
			#print outputfile
		elif opt in ("-p"):
			predsfile= arg
			#print predsfile
		elif opt in ("-l"):
			label=arg
			#print label
		elif opt in ("-s"):
			#print 's: '+arg
			cutoff=float(arg)
			#print cutoff
		elif opt in ("-r"):
			min_len=int(arg)
		elif opt in ("--min"):
			min_len=int(arg)
			#print 'max found'
		elif opt in ("--max"):
			max_len=int(arg)
			#print 'max found'
	
	if cutoff>1 or cutoff<0:
		print('please enter probability value in range [0,1]')
		sys.exit()
	if max_len<=min_len:
		print('Error: check min and max len')
		sys.exit()
	if label != '0' and label != '1':
		print('Please check lablel is 0 or 1')
		sys.exit()

	print('**********Extracting Sequences***************')
	print(('class prob cutoff='+ str(cutoff)))
	print(('min length cutoff='+ str(min_len)))
	print(('max length cutoff='+ str(max_len)))

	#label=sys.argv[4]
	#open preds file
	idlist=[]
	with open(predsfile) as f:
		content=f.readlines()

	for l in content:
	
		if label=='1':
			tocheck=float(l.split('\t')[2])
		elif label=='0':
			tocheck=float(l.split('\t')[3])
		else:
			print('check label\nError')
			sys.exit(0)
		
		if l.split('\t')[1]==label:
			if cutoff <= tocheck :
				idlist.append(l.split('\t')[0])
	print(('Total sequences in prediction file with label '+label+' and class prob >= '+str(cutoff)+', were: '+str(len(idlist))))
	ctr=0
	min_len_filter=0
	max_len_filter=0
	#extract seq from fasta
	output_handle = open(outputfile, "w")
	for record in SeqIO.parse(inputfile, "fasta"):
		found=0
		seqid=record.id
		#print record.seq
		if seqid in idlist :
			#write to file
			#print 'found'
			if len(record.seq)>max_len:
				max_len_filter=max_len_filter+1
			elif len(record.seq)<min_len:
				min_len_filter=min_len_filter+1
			else:
				#print len(record.seq)
				ctr=ctr+1
				SeqIO.write(record, output_handle, "fasta")
				found =1

			
		
	
	
	print(('Total filtered due to length < '+str(min_len)+', were: '+str(min_len_filter)))
	print(('Total filtered due to length > '+str(max_len)+', were: '+str(max_len_filter)))
	print(('Sequences not found in the fasta file were: '+str(len(idlist)-ctr-min_len_filter-max_len_filter)))
	if len(idlist)-ctr-min_len_filter-max_len_filter>0:
		print('WARNING:')
		print(('Please check input fasta as '+str(len(idlist)-ctr-min_len_filter-max_len_filter)+ ' sequences in the prediction file did not match to any sequences in fasta file'))
	print(('Total sequences written: '+str(ctr)))
	print(('File '+outputfile+' saved!'))

if __name__ == "__main__":
	main()
