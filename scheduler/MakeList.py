#!/usr/bin/python

from glob import glob
from subprocess import Popen, PIPE
import os
import sys

from argument_parser import argument_parser

#_____________________________________________________________________________
if __name__ == "__main__":

    #merge output files by chunks of a given size

    #config file from command line argumet
    args = sys.argv
    args.pop(0)
    config = args.pop(0)

    parser = argument_parser()
    parser.add_parameter("top")
    parser.add_parameter("add_input")

    parser.parse(config)

    top = parser.get("top")
    qlist = parser.get("add_input")

    basedir = parser.get("top") + parser.get("add_input")
    print basedir
    cmd = "ls " + basedir + "/*.root"
    out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()[0].split("\n")

    file_list = open(basedir+"/picoDst.list", "w")

    for fline in out:
        if len(fline) == 0: continue
        file_list.write(fline+"\n")


    print "List created: " + basedir+ "/picoDst.list"
























