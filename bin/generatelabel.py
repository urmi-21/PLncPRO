'''
This file just generates labels i.e 0 or 1 and preints them on terminal
It may be helpful to generate labels for known data and save them using '>'
[1]-->number of 1's to print
[2]-->number of 0's to print
e.g. python generatelabel.py 200 200 > labels.txt
'''
import sys
for i in range(int(sys.argv[1])):
	print(1)
for i in range(int(sys.argv[2])):
	print(0)
