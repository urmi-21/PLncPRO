'''
This file reads feature file and predicts classes. this file reads a rf classifier model and the does prediction.
The label is the first column of each file and first line contains the feature names.
Files are give as arguments and a second argument id fo model name used to read already built model from file
argv[1]-->input fa file
argv[2]-->rf model file
argv[3]-->outputfile name
argv[4]-->if known labels are present
argv[5]--> known label filename
UrMi 23/01/2016
'''

from sklearn.ensemble import RandomForestClassifier
from numpy import genfromtxt, savetxt
import numpy as np
import math
import sys
import pickle
import os
from sklearn.metrics import classification_report
from sklearn.metrics import roc_curve, auc
#import matplotlib.pyplot as plt

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
    
    #read first line from both files and check for consistency
    features=genfromtxt(open(sys.argv[1],'r'), delimiter='\t', dtype=None,encoding=None)[0]
    features=features[2:]
    text=genfromtxt(open(sys.argv[1],'r'), delimiter='\t', dtype=None,encoding=None)[1:]
    seqids=[x[1] for x in text]    
    #read file
    dataset = genfromtxt(open(sys.argv[1],'r'), delimiter='\t', dtype='f8')[1:]
    target = [x[0] for x in dataset] # not required for prediction
    X = [x[2:] for x in dataset]
    with open(sys.argv[2], 'rb') as f:
        rf = pickle.load(f)
        preds = rf.predict(X)
        preds_prob= rf.predict_proba(X)

    #print preds
    #print preds_prob    
    #print 'lens:'
    #print len(preds),len(preds_prob)
    num_pos=0
    num_neg=0
    for x in preds:
        if x==0:
            num_neg=num_neg+1        
        else:
            num_pos=num_pos+1

    print('\t\tPredicted as neg: '+str(num_neg)+': '+str(num_neg/len(preds)))
    print('\t\tPredicted as pos: '+str(num_pos)+': '+str(num_pos/len(preds)))

    #write predictions to file
    #fname=os.path.dirname(sys.argv[1])+'/'+sys.argv[3]
    fname=sys.argv[3]
    f=open(fname,'w')
    for i in range(len(seqids)):
        f.write(str(seqids[i])+'\t'+str(int(preds[i]))+'\t'+str(float(preds_prob[i][1]))+'\t'+str(float(preds_prob[i][0]))+'\n')
    #print 'Predictions written to file: '+fname
    if('true' in sys.argv[4]):
        print('\t\tCalculating Performance of Classifier...')
        #read file with known labels
        label=[]
        with open(sys.argv[5]) as f:
            content=f.readlines()
        for x in content:
            label.append(int(x.split('\n')[0]))

        #check if num labels == num examples
        #print 'len preds=',str(len(preds))
        #print 'len label=',str(len(label))
        if not(len(preds) == len(label)):
            print('Error, Num of labels does not match Num of examples')
            sys.exit(1)
        
        #start calculation
        tp=0
        tn=0
        fp=0
        fn=0
        p=0
        n=0
        for i in range(0,len(preds)):
            if label[i]==1 and preds[i]==1:
                tp=tp+1
            elif label[i]==1 and preds[i]==0:
                fn=fn+1
            elif label[i]==0 and preds[i]==0:
                tn=tn+1
            elif label[i]==0 and preds[i]==1:
                fp=fp+1
            else:
                print('Error in Predictions!!!')
                sys.exit(1)
        p=tp+fn
        n=tn+fp
        print(bcolors.FAIL)

        print('\t\t\tPrediction')
        print('\t'+str('_______________________________________________________________'))
        print('\t\t'+str('Pos')+'\t|\t'+str('Neg')+'\t|\t'+str('Total'))
        print('\t\t'+str(tp)+'\t|\t'+str(fn)+'\t|\t'+str(p))
        print('\t\t'+str(fp)+'\t|\t'+str(tn)+'\t|\t'+str(n))
        print('\t'+str('_______________________________________________________________'))

        print('\t\tp:',str(p))
        print('\t\tn:',str(n))
        print('\t\ttp:',str(tp))
        print('\t\ttn:',str(tn))
        print('\t\tfp:',str(fp))
        print('\t\tfn:',str(fn))
        
        try:        
            sens=tp/p
        except ZeroDivisionError:
            #print "division by zero...can't calculate Sens"
            sens='NaN'
        try:
            spc=tn/n
        except ZeroDivisionError:
            #print "division by zero...can't calculate Spc"
            spc='NaN'
        acc=((tp+tn)/(p+n))*100
        mcc=(tp*tn)-(fp*fn)
        den=math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
        try:
            mcc=mcc/den
        except ZeroDivisionError:
            #print "division by zero...can't calculate MCC"
            mcc='NaN'
        print('\t\tAcc=',str(acc))
        print('\t\tSens=',str(sens))
        print('\t\tSpc=',str(spc))
        print('\t\tMCC=',str(mcc))        
        
        #for i in range(len(seqids)):
        #    print seqids[i],preds[i],label[i]
        # Determine the false positive and true positive rates
        #print rf.predict_proba(X)[:,1]
        fpr, tpr, _ = roc_curve(label, rf.predict_proba(X)[:,1])
        #print 'len prob fpr','tpr'        
        #print len(fpr),len(tpr)
        # Calculate the AUC
        roc_auc = auc(fpr, tpr)
        print('\t\tROC AUC: %0.4f' % roc_auc)
        print(bcolors.ENDC)

#main function does everything. Keep this
if __name__=="__main__":
    main()
