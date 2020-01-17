'''
PLNCPRO
This file merges all the features created by extractfeatures.py, blastparse.py and ffparse.py
The program takes as input the three features files in the following order
extractfeatures.py blastparse.py ffparse.py
the program writes a final mergedfeature file
arguments 
1-->feature file frome extractfeatures.py
2-->framefinder flag
3-->framefinder feature file
4-->no blast flag
5-->blastres feature file
6-->outfile name
Author : Urminder Singh
email: urmind13_sit@jnu.ac.in
UrMi 14/1/2016
'''

import sys
#import math
#import re


#read exfeature file
featuredata=[]
with open(sys.argv[1]) as f:
	featuredata=f.readlines()
no_framefinder_flag=sys.argv[2]
if no_framefinder_flag=='false':
	#read framefinder file feature file
	with open(sys.argv[3]) as f:
		content=f.readlines()
	#check that featuredata and content are of same length

	if len(featuredata)!=len(content):
		print('error in framefinder file\nPlease check framefinder output file...')
		sys.exit(0)
	elif len(featuredata)==len(content):
		#print 'data_good'
		#merge content in featuredata
		for i in range(len(featuredata)):
			#print content[i]
			#merge at the end
			featuredata[i]=featuredata[i].split('\n')[0]+'\t'+content[i].split('\t')[1]+'\t'+content[i].split('\t')[2] #1-->score and 2-->cov
			#featuredata[i]=featuredata[i].split('\n')[0]+'\t'+content[i].split('\t')[2]

#print featuredata[0]
#for x in featuredata:
#	print x
		
#open blastfeature file and merge
#here some ids may not be present in blast results so match both featuredata and blastfeaturefile
#by query id and then merge coressponding rows
#read framefinder file feature file

no_blast_flag=sys.argv[4]
if no_blast_flag=='false':
	with open(sys.argv[5]) as f:
		blastfeature=f.readlines()

	found_flag=0
	##write firts line
	#featuredata[0]=featuredata[0].split('\n')[0]+'\t'+blastfeature[0].split('\t')[1]+'\t'+blastfeature[0].split('\t')[2]+'\t'+blastfeature[0].split('\t')[3].split('\n')[0]
	for i in range(0,len(featuredata)):
		found_flag=0
		#qid is second col, first col is label
		qid=featuredata[i].split('\t')[1]
		#if no_ff_flag
		qid=qid.split('\n')[0]
		#print qid
		for x in blastfeature:
			#print 'pair:',qid,x
			if qid in x:
				featuredata[i]=featuredata[i].split('\n')[0]+'\t'+x.split('\t')[1]+'\t'+x.split('\t')[2]+'\t'+x.split('\t')[3]+'\t'+x.split('\t')[4].split('\n')[0]
				#featuredata[i]=featuredata[i].split('\n')[0]+'\t'+x.split('\t')[4].split('\n')[0] #only numhits
				found_flag=1
				#print 'here found***********************************************************************'
				break;
			#if no entry found
			found_flag=0
		if found_flag==0:
			featuredata[i]=featuredata[i].split('\n')[0]+'\t'+str(-1)+'\t'+str(0)+'\t'+str(2)+'\t'+str(0)
			#featuredata[i]=featuredata[i].split('\n')[0]+'\t'+str(0)#for numhits
	

#finally write the merged feature file
fname=sys.argv[6]
f=open(fname,'w')
for x in featuredata:
	f.write(x)
	f.write('\n')
#print fname+', File written'
