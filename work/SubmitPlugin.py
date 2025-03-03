#!/usr/bin/python

#--------------------------------------------------------------
# macro to submit rpDst analysis 
#
# usage:
# ./SubmitPlugin.py tag inputSource
# ./SubmitPlugin.py tag
#
#--------------------------------------------------------------

import os
from subprocess import Popen
import re
import sys
import time
import shutil


#_____________________________________________________________________________
if __name__ == "__main__":

    #name of config file from command line argument
    args = sys.argv
    args.pop(0) # cut first input from terminal = ./SubmitPlugin
    if len(args) == 2:
        tag = args.pop(0) # read third input from terminal = tag (outputDir)
        inputSource = args.pop(0) # read second input from terminal = inputSource
    elif len(args) == 1:
        tag = args.pop(0) # read third input from terminal = tag (outputDir)
        #inputSource = "/star/u/truhlar/star-upcDst/work/lists/forwardRPTest.list"
        inputSource = "/star/u/truhlar/star-upcDst/work/lists/extendedUpc.list"
        #inputSource = "/gpfs01/star/pwg/truhlar/Embedding/PiplusPiminus/output.list"   
        #inputSource = "/gpfs01/star/pwg/truhlar/Embedding/KplusKminus/output.list"   
        #inputSource = "/gpfs01/star/pwg/truhlar/Embedding/ProtonProtonbar/output.list"   
    else:
        print "Wrong input arguments"
        exit()

    exe = "AnalysisManager"

    #working directory is assumed the same as for this macro
    workDir = os.getcwd()

    #directory with source codes, one level up
    srcDir = workDir[:workDir.rfind("/")]

    top = "/gpfs01/star/pwg/truhlar/Run17_P20ic"
    basedir = top + "/" + tag
    inputSource = "<input URL=\"filelist:" + inputSource + "\"/>"
    print("Input source: "+inputSource)
    print("Output dir: "+basedir)
    print("Source dir: "+srcDir)
    shutil.copyfile(srcDir + "/work/include/RunDef.h", srcDir + "/work/config/" + tag + ".h")
    #create jobs output folders
    Popen("mkdir -p "+basedir, shell=True).communicate()
    #proceed with scheduler submit
    Popen("mkdir "+basedir+"/logs", shell=True).communicate()
    Popen("mkdir "+basedir+"/sched", shell=True).communicate()
    #create xml
    xmlname = "submit.xml"
    xml = open(xmlname, "w")
    for line in open("scheduler_template.xml", "r"):
        line = re.sub(r"__BASEDIR__", basedir, line)
        line = re.sub(r"__INPUT__", inputSource, line)
        line = re.sub(r"__SRCDIR__", srcDir, line)
        line = re.sub(r"__EXE__", exe, line)
        xml.write(line)
    xml.close()
    #submit for a given catalog query
    os.chdir("sched")
    Popen("star-submit ../"+xmlname, shell=True).communicate()
    print "All done."




















