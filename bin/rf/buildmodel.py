'''
This file reads two feature files and builds rf model classifier.
The label is the first column of each file and first line contains the feature names.
Files are give as arguments and a third argument id fo model name used to save model to file
UrMi 23/01/2016
'''

from sklearn.ensemble import RandomForestClassifier
from numpy import genfromtxt, savetxt
import numpy as np
import math
import sys
import pickle
import os

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
def main():
	#print os.path.dirname(sys.argv[1])
	#iterations to build model
	num_iterations=4
	#model parameters as arguments
	num_trees=int(sys.argv[3])
	num_threads=int(sys.argv[4])
	#read files and save in two diff datasets
	#read first line from both files and check for consistency
	features=genfromtxt(open(sys.argv[1],'r'), delimiter='\t', dtype=None,encoding=None)[0]
	features=features[2:]
	#read file
	dataset = genfromtxt(open(sys.argv[1],'r'), delimiter='\t', dtype='f8')[1:]
	target = [x[0] for x in dataset]
	train = [x[2:] for x in dataset]
	
	
	#create and train the random forest
	#multi-core CPUs can use: rf = RandomForestClassifier(n_estimators=100, n_jobs=2)
	rf = RandomForestClassifier(n_estimators=num_trees,oob_score=True,n_jobs=num_threads)
	rf.fit(train, target)
	oobscore=rf.oob_score_
	#print 'oob=',str(oobscore)
	model_fname=os.path.dirname(sys.argv[1])+'/'+sys.argv[2]
	with open(model_fname, 'wb') as f:
		pickle.dump(rf, f)
	#build model num_iterations times and save with highest oob score
	for i in range(num_iterations):
		rf = RandomForestClassifier(n_estimators=num_trees,oob_score=True,n_jobs=num_threads)
		rf.fit(train, target)
		thisoobscore=rf.oob_score_
		#print 'oob=',str(thisoobscore)
		if thisoobscore>oobscore:
			oobscore=thisoobscore
			#print 'writing to file...'
			with open(model_fname, 'wb') as f:
				pickle.dump(rf, f)
	print(bcolors.HEADER+'****************************************************************Model Details****************************************************************'		+ bcolors.ENDC)
	#print '\t\t\t\t\tFinal oob score=:',str(oobscore)
	print(bcolors.OKBLUE + "\t\t\t\t\tFinal oob score=: "+str(oobscore) + bcolors.ENDC)
	print(bcolors.OKBLUE + "\t\t\t\t\tUse this model with the predict program\n\t\t\t\t\tRelative feature importance:\n"+ bcolors.ENDC)
	for f in sorted(zip([round(x, 4) for x in rf.feature_importances_], features),reverse=True):
		print(f, end=' ')
'''in prediction file
with open('path/to/file', 'rb') as f:
    rf = cPickle.load(f)
preds = rf.predict(new_X)
'''
	
if __name__=="__main__":
	main()
