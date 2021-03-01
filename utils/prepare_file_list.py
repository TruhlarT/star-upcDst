#!/usr/bin/python

#--------------------------------------------------------------
# macro to unite upcDst file into form, where there is exactly 1 file for 1 runNumber
# The input list must be sorted via runNumbers
# For that purpose use: sort -t "/" -k 7,7 "tmpDefault.list" 
# usage:
# ./upcUnite.py  inputSource
# ./upcUnite.py 
#
#--------------------------------------------------------------

import os
from subprocess import Popen
import re
import sys
import time


#_____________________________________________________________________________
if __name__ == "__main__":

    # Setting the input list of upcDst files
    args = sys.argv
    args.pop(0) # cut first input from terminal = ./SubmitPlugin
    if len(args) == 1:
        inputSource = args.pop(0) # read second input from terminal = inputSource
    else:
        print "missing input"
        exit()


    with open("doneNames.list") as f:
        doneList = f.read().splitlines()

    fileSet = set()

    processed_name = "muDstPart.list"
    file_prN = open(processed_name, "w")

    for line in open(inputSource, "r"):
        file_name = line[:line.find(".root")]
        file_name = file_name[line.rfind("/") + 1:]
        if file_name not in doneList:
            if file_name not in fileSet:
                file_prN.write(line)
                fileSet.add( file_name)
            else:
                print "Already in set: " + file_name

    file_prN.close()
    processed_name = "processedNames.list"
    file_prN = open(processed_name, "w")

    for name in fileSet:
        file_prN.write(name+"\n");


    file_prN.close()












