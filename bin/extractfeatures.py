'''
PLNCPRO
This file reads a fasta or multiple fasta file and extracts
various features. It takes the input fasta file as input and creates an output file
Author : Urminder Singh
email: urmind13_sit@jnu.ac.in
UrMi 13/1/16
'''
##Note to self:This file is a mess please refractor

import sys
#import math
#import re
from Bio import SeqIO
from Bio.Seq import Seq
import regex as re
import time
start_time = time.time()
########################################
##Define 64 Trimers the values are not used by PlncPRO I think these values were DNA bending propensities.
trimervalues = {}
trimervalues["AAT"]=-0.280
trimervalues["ATT"]=-0.280
trimervalues["AAA"]=-0.274
trimervalues["TTT"]=-0.274
trimervalues["CCA"]=-0.246
trimervalues["TGG"]=-0.246
trimervalues["AAC"]=-0.205
trimervalues["GTT"]=-0.205
trimervalues["ACT"]=-0.183
trimervalues["AGT"]=-0.183
trimervalues["CCG"]=-0.136
trimervalues["CGG"]=-0.136
trimervalues["ATC"]=-0.110
trimervalues["GAT"]=-0.110
trimervalues["AAG"]=-0.081
trimervalues["CTT"]=-0.081
trimervalues["CGC"]=-0.077
trimervalues["GCG"]=-0.077
trimervalues["AGG"]=-0.057
trimervalues["CCT"]=-0.057
trimervalues["GAA"]=-0.037
trimervalues["TTC"]=-0.037
trimervalues["ACG"]=-0.033
trimervalues["CGT"]=-0.033  
trimervalues["ACC"]=-0.032  
trimervalues["GGT"]=-0.032  
trimervalues["GAC"]=-0.013  
trimervalues["GTC"]=-0.013  
trimervalues["CCC"]=0.012  
trimervalues["GGG"]=0.012  
trimervalues["ACA"]=-0.006  
trimervalues["TGT"]=-0.006  
trimervalues["CGA"]=-0.003  
trimervalues["TCG"]=-0.003  
trimervalues["GGA"]=0.013  
trimervalues["TCC"]=0.013  
trimervalues["CAA"]=0.015  
trimervalues["TTG"]=0.015  
trimervalues["AGC"]=0.017  
trimervalues["GCT"]=0.017  
trimervalues["GTA"]=0.025  
trimervalues["TAC"]=0.025  
trimervalues["AGA"]=0.027  
trimervalues["TCT"]=0.027  
trimervalues["CTC"]=0.031  
trimervalues["GAG"]=0.031  
trimervalues["CAC"]=0.040  
trimervalues["GTG"]=0.040  
trimervalues["TAA"]=0.068  
trimervalues["TTA"]=0.068  
trimervalues["GCA"]=0.076  
trimervalues["TGC"]=0.076  
trimervalues["CTA"]=0.090  
trimervalues["TAG"]=0.090 
trimervalues["GCC"]=0.107  
trimervalues["GGC"]=0.107  
trimervalues["ATG"]=0.134  
trimervalues["CAT"]=0.134  
trimervalues["CAG"]=0.175  
trimervalues["CTG"]=0.175  
trimervalues["ATA"]=0.182  
trimervalues["TAT"]=0.182  
trimervalues["TCA"]=0.194  
trimervalues["TGA"]=0.194
###############################
###only five trimers
'''trimervalues = {}
trimervalues["CGA"]=0
trimervalues["GCG"]=0
trimervalues["CGC"]=0
trimervalues["CCG"]=0
trimervalues["TCG"]=0
'''

################################
#function to calculate entopy of the DNA sequence, give input counts of A G C and T
def getlog2(val):
	if val == 0:
		return 0
	else:
		return math.log(val,2)
def returnentropy(a,g,c,t):
	pa=0.0
	pt=0.0
	pc=0.0
	pg=0.0
	size=a+g+c+t
	pa=a/size
	pt=t/size
	pc=c/size
	pg=g/size
	entrpy=(-1)*((pa*getlog2(pa))+(pt*getlog2(pt))+(pc*getlog2(pc))+(pg*getlog2(pg)))
	return entrpy

def returncount(data,char):
	ctr=0
	for x in data:
		if x == char:
			ctr=ctr+1
	return ctr

def returntrimercount(data):
	
	ctr=0	
	length=len(data)
	#calculate GC1 GC2 GC3
	GC1=0
	GC2=0
	GC3=0
	#print length
	data_str=''.join(data)
	for x in trimervalues:
		trimer_str=''.join(x)
		matches = re.findall(r''+trimer_str, data_str, overlapped=True)
		trimervalues[trimer_str]=len(matches)
	
	#normalise trimer count by len-2
	for x in trimervalues:
		trimervalues[x]=trimervalues[x]/(length-2)*100
		#print x,trimervalues[x]
		
	return ctr

def gethexmatrix(seq):
	a1=0
	a2=0
	a3=0
	a4=0
	a5=0
	a6=0
	g1=0
	g2=0
	g3=0
	g4=0
	g5=0
	g6=0
	t1=0
	t2=0
	t3=0
	t4=0
	t5=0
	t6=0
	c1=0
	c2=0
	c3=0
	c4=0
	c5=0
	c6=0

	hexmer=[]
	for i in range(len(seq)-5):
		hexstr=''
		for j in range(i,i+6):
			#print seq[j],
			hexstr=hexstr+seq[j]
			
		hexmer.append(hexstr)
		#print ''
	#print len(hexmer)
	for x in hexmer:
		if x[0]=='A':
			a1=a1+1
		elif x[0]=='G':
			g1=g1+1
		elif x[0]=='C':
			c1=c1+1
		elif x[0]=='T':
			t1=t1+1

		if x[1]=='A':
			a2=a2+1
		elif x[1]=='G':
			g2=g2+1
		elif x[1]=='C':
			c2=c2+1
		elif x[1]=='T':
			t2=t2+1

		if x[2]=='A':
			a3=a3+1
		elif x[2]=='G':
			g3=g3+1
		elif x[2]=='C':
			c3=c3+1
		elif x[2]=='T':
			t3=t3+1

		if x[3]=='A':
			a4=a4+1
		elif x[3]=='G':
			g4=g4+1
		elif x[3]=='C':
			c4=c4+1
		elif x[3]=='T':
			t4=t4+1

		if x[4]=='A':
			a5=a5+1
		elif x[4]=='G':
			g5=g5+1
		elif x[4]=='C':
			c5=c5+1
		elif x[4]=='T':
			t5=t5+1


		if x[5]=='A':
			a6=a6+1
		elif x[5]=='G':
			g6=g6+1
		elif x[5]=='C':
			c6=c6+1
		elif x[5]=='T':
			t6=t6+1


	#print a1/len(hexmer),a2/len(hexmer),a3/len(hexmer),a4/len(hexmer),a5/len(hexmer),a6/len(hexmer)
	#print g1/len(hexmer),g2/len(hexmer),g3/len(hexmer),g4/len(hexmer),g5/len(hexmer),g6/len(hexmer)
	#print c1/len(hexmer),c2/len(hexmer),c3/len(hexmer),c4/len(hexmer),c5/len(hexmer),c6/len(hexmer)
	#print t1/len(hexmer),t2/len(hexmer),t3/len(hexmer),t4/len(hexmer),t5/len(hexmer),t6/len(hexmer)
	#print '*********'
	res=str(a1/len(hexmer))+'\t'+str(a2/len(hexmer))+'\t'+str(a3/len(hexmer)) +'\t'+str(a4/len(hexmer)) +'\t'+ str(a5/len(hexmer))+'\t'+str(a6/len(hexmer))
	res=res+'\t'+str(g1/len(hexmer))+'\t'+str(g2/len(hexmer))+'\t'+str(g3/len(hexmer)) +'\t'+str(g4/len(hexmer)) +'\t'+ str(g5/len(hexmer))+'\t'+str(g6/len(hexmer))
	res=res+'\t'+str(c1/len(hexmer))+'\t'+str(c2/len(hexmer))+'\t'+str(c3/len(hexmer)) +'\t'+str(c4/len(hexmer)) +'\t'+ str(c5/len(hexmer))+'\t'+str(c6/len(hexmer))
	res=res+'\t'+str(t1/len(hexmer))+'\t'+str(t2/len(hexmer))+'\t'+str(t3/len(hexmer)) +'\t'+str(t4/len(hexmer)) +'\t'+ str(t5/len(hexmer))+'\t'+str(t6/len(hexmer))
	return res


###############################################################################
###########################START MAIN##########################################
###########################END OF FUNCTIONS####################################
###############################################################################
#do for each sequence
#label should be either 1 or 0
label=int(sys.argv[2])

if label != 1 and label!=0:
	print('please check label value; 0 or 1')
	sys.exit(0)
fname=sys.argv[1]+'_features'
f=open(fname,'w')
firstline_flag=0
#print 'extracting features...'
#print 'Total Sequences: '+ str(len(list(SeqIO.parse(sys.argv[1], "fasta"))))
for record in SeqIO.parse(sys.argv[1], "fasta"):
	seqid=record.id
	#calculate composition
	length=len(record)
	
		
	'''countA=returncount(record.seq,'A')
	countG=returncount(record.seq,'G')
	countC=returncount(record.seq,'C')
	countT=returncount(record.seq,'T')
	perA=countA/length*100
	perG=countG/length*100
	perC=countC/length*100
	perT=countT/length*100
	perGC=perG+perC
	perAT=perA+perT
	'''
	#print perA,perG,perC,perT,perA+perG+perC+perT,perGC,perAT
	#calculate entropy
	#entropy=returnentropy(countA,countG,countC,countT)
	#print'entropy= '+str(returnentropy(countA,countG,countC,countT))
	#calculate trimers	
	returntrimercount(record)
	
	#create string of trimer values
	tstr=''
	tstrval=''
	for x in trimervalues:
		tstr=tstr+'\t'+str(trimervalues[x])
		tstrval=tstrval+'\t'+str(x)
	
	#remove extra space
	tstr=tstr[1:]
	tstrval=tstrval[1:]
	
	#hexvals='a1\ta2\ta3\ta4\ta5\ta6\tg1\tg2\tg3\tg4\tg5\tg6\tc1\tc2\tc3\tc4\tc5\tc6\tt1\tt2\tt3\tt4\tt5\tt6'
	#hexres=gethexmatrix(record.seq)
	#print hexres
	#print str(seqid),str(length),str(perA),str(perG),str(perC),str(perT),str(perGC),str(perAT),str(entropy),tstr
	#print feature names first
	if firstline_flag==0:
		#f.write(str('Label')+'\t'+str('seqid')+'\t'+str('length')+'\t'+str('perA')+'\t'+str('perG')+'\t'+str('perC')+'\t'+str('perT')+'\t'+str('perGC')+'\t'+str('perAT')+'\t'+str('entropy')+'\t'+tstrval+'\t'+hexvals+'\n')
		#f.write(str('Label')+'\t'+str('seqid')+'\t'+str('length')+'\t'+hexvals+'\n')
		#f.write(str('Label')+'\t'+str('seqid')+'\t'+str('length')+'\t'+str('perA')+'\t'+str('perG')+'\t'+str('perC')+'\t'+str('perT')+'\t'+str('perGC')+'\t'+str('perAT')+'\t'+str('entropy')+'\t'+tstrval+'\n')
		f.write(str('Label')+'\t'+str('seqid')+'\t'+tstrval+'\t'+str('length')+'\n')
		#f.write(str('Label')+'\t'+str('seqid')+'\t'+str('length')+'\n')
		
		firstline_flag=1
		
	
	#f.write(str(label)+'\t'+str(seqid)+'\t'+str(length)+'\t'+str(perA)+'\t'+str(perG)+'\t'+str(perC)+'\t'+str(perT)+'\t'+str(perGC)+'\t'+str(perAT)+'\t'+str(entropy)+'\t'+tstr+'\t'+hexres+'\n')
	#f.write(str(label)+'\t'+str(seqid)+'\t'+str(length)+'\t'+hexres+'\n')
	#f.write(str(label)+'\t'+str(seqid)+'\t'+str(length)+'\t'+str(perA)+'\t'+str(perG)+'\t'+str(perC)+'\t'+str(perT)+'\t'+str(perGC)+'\t'+str(perAT)+'\t'+str(entropy)+'\t'+tstr+'\n')
	f.write(str(label)+'\t'+str(seqid)+'\t'+tstr+'\t'+str(length)+'\n')
	#f.write(str(label)+'\t'+str(seqid)+'\t'+str(length)+'\n')
	
#print 'Done! '+str(fname)+' saved!'

#print("--- %s seconds ---" % (time.time() - start_time))
