#!/usr/bin/python


import os
from subprocess import Popen, PIPE
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



    processed_name = "badMuDst2.list"
    file_prN = open(processed_name, "w")

    badSet = set()

    for line in open(inputSource, "r"):
        badSet.add(line[:len(line)-1])


    for line in open("part1.list"):
        file_name = line[line.find("/st_rp")-8:line.find("/st_rp")]
        if file_name in badSet:
            file_prN.write(line)

    file_prN.close()












