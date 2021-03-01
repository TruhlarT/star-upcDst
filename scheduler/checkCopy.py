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
    parser.add_parameter("outlist")
    parser.add_parameter("add_input", list)

    parser.parse(config)

    top = parser.get("top")
    qlist = parser.get("add_input")
    outlist = parser.get("outlist")


    for q in qlist:
        log_path = top + q[0] + "/logs"
        cmd = "ls " + log_path + "/*.out"
        logs = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()[0].split("\n")
        #remove empty elements (last was empty string)
        logs = [i for i in logs if i]
        
        nFiles = 0
        for log in logs:
            for line in open(log, "r"):
                if line.find("RunFilterMaker, nFiles:") != -1:
                    nFiles = nFiles + int(line[len("RunFilterMaker, nFiles: "):])
                    break

        count = len(open(q[1]).readlines(  ))
        if nFiles != count :
            print "There is a problem! Check logs for " + q[0]
            print "Logs: " + str(nFiles) + " List: " + str(count)
























