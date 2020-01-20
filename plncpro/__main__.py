#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 13:31:41 2020

@author: usingh

plncpro main script
"""



import sys
import plncpro.prediction
import plncpro.build
import plncpro.predstoseq



######################################
class bcolors:
	FAIL = '\033[91m'
	ENDC = '\033[0m'
    
######################################
print((bcolors.FAIL))
print ('\t\t\t\t  _____  _                  _____    _____     ____  ')
print ('\t\t\t\t |  __ \| |                |  __ \  |  __ \   / __ \ ')
print ('\t\t\t\t | |__) | |  _ __     ___  | |__) | | |__) | | |  | |')
print ('\t\t\t\t |  ___/| | |  _ \   / __| |  ___/  |  _  /  | |  | |')
print ('\t\t\t\t | |    | | | | | | | (__  | |      | | \ \  | |__| |')
print ('\t\t\t\t |_|    |_| |_| |_|  \___| |_|      |_|  \_\  \____/ ')
print((bcolors.ENDC))



def main():
    """The main routine."""
    validCommands=['predict','build','predtoseq']
    if len(sys.argv) < 2:
        print("Valid commands are:\t"+"\t".join(validCommands))
        sys.exit(1)
        
    thisCommand=sys.argv[1]
        
    if thisCommand not in validCommands:
        print("Valid commands are:\t"+"\t".join(validCommands))
    elif thisCommand==validCommands[0]:
        print ("calling predict")
        plncpro.prediction.main()
        
    elif thisCommand==validCommands[1]:
        plncpro.build.main()
        
    elif thisCommand==validCommands[2]:
        print ("calling predtos")
    

if __name__ == "__main__":
    main()