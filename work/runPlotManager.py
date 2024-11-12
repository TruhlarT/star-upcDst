#!/usr/bin/python

from glob import glob
from subprocess import Popen, PIPE
import os
import sys
import re
import shutil
import filecmp

#_____________________________________________________________________________
if __name__ == "__main__":

    #get the name of root file to process
    args = sys.argv
    args.pop(0) # cut first input from terminal = ./runPlotManager.py
            
    currDir = os.getcwd()
    
    basedir = "/gpfs01/star/pwg/truhlar/Run17_P20ic/"
    if len(args) == 1:
        tag = str(args.pop(0))
        if tag == "l":
            infile = "AnalysisOutput.root"
            location = currDir + "/"
        else:
            infile = "StRP_production_0000.root"
            location = basedir + tag + "/merged/"
    else: 
        #get outputs directory from submit used for production
        for line in open("submit.xml").xreadlines(  ): 
            if "Location" in line:
                basedir = line.lstrip().split(">")[1].split("/sched")[0]
                tag = basedir.rsplit("/",1)[1]
        infile = "StRP_production_0000.root"
        location = basedir + "/merged/"

    #Copy libstar-upc.so
    shutil.copyfile(currDir + "/../build/libstar-upc.so", "libstar-upc.so")

    print "  Input file:     ", infile
    print "  File location:     ", location
    #print "  Tag :      ", tag
    src_file = currDir + "/include/RunDef.h"
    backup_file = currDir + "/include/RunDef_backup.h"
    new_file = currDir + "/config/" + tag + ".h"

    # Check if the source file and new file are different
    if tag != "l" and not filecmp.cmp(src_file, new_file, shallow=False):
            # Copy RunDef.h to backup
            shutil.copyfile(src_file, backup_file)
            # Replace RunDef.h with the new file
            shutil.copyfile(new_file, src_file)
            print "  RunDef.h updated by ", new_file

    make = Popen("make", stdout=PIPE, stderr=PIPE).communicate()
    if len(make[1]) != 0:
        print make[1]
    #exit()


    #embedFile = embedDir + "/merged/StRP_production_0000.root"
    root_cmd = "./PlotManager " + location + infile #"\",\"" + embedFile + "\")"
    root = Popen(root_cmd.split(), stdout=PIPE, stderr=PIPE).communicate()

    print root[0], root[1]
    if tag != "l" and not filecmp.cmp(src_file, new_file, shallow=False):
        shutil.copyfile(currDir + "/include/RunDef_backup.h", currDir + "/include/RunDef.h")
    print "All done."
