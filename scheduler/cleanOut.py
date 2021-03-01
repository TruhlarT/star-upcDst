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

        for log in logs:
            firstLine=[]
            with open(log,'r') as myFile:
                for indx, line in enumerate(myFile, 1):
                    if line.find("We END on") != -1 or line.find("We are now leaving with error on") != -1:
                        firstLine.append(indx)

            myFile.close()
            if len(firstLine) > 1:
                fout = open(log_path + "/final.out", "w")
                with open(log,'r') as myFile:
                    for indx, line in enumerate(myFile, 1):
                        if indx > firstLine[len(firstLine)-2]:
                            fout.write(line)
                fout.close()
                log_newname = log[:log.find(".out")] + "_old.out"
                os.rename(log,log_newname)
                os.rename(log_path + "/final.out",log) 
                print "Changed log: " + log    






















