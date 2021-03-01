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


    fileSet = set()

    for line in open(inputSource, "r"):
        file_name = line[:line.find("root") + 4]
        #print file_name
        fileSet.add(file_name)


    processed_name = "part1Extra.list"
    file_prN = open(processed_name, "w")

    for line in open("part1.list","r"):
        file_name = line[line.rfind("/") + 1:line.find(".root") + 5]
        if file_name in fileSet:
            #print file_name
            file_prN.write(line);


    file_prN.close()












