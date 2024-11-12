#!/usr/bin/python


import os
from subprocess import Popen
from os.path import exists
import re
import sys
import time


#_____________________________________________________________________________
if __name__ == "__main__":

    # Setting the input list of upcDst files
    args = sys.argv
    args.pop(0) # cut first input from terminal
    if len(args) == 1:
        inputSource = args.pop(0) # read second input from terminal = inputSource
    else:
        print "missing input"
        exit()

    for line in open(inputSource, "r"):
        file_name = line[:line.find(".root")]
        file_name = file_name[line.rfind("/") - 8:line.rfind("/")]
        file_path = os.path.join("lists", file_name + ".list")
        with open(file_path, "a") as file_prN:
            file_prN.write(line)
            file_prN.close()












