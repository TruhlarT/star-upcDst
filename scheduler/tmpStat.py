#!/usr/bin/python

from glob import glob
from subprocess import Popen, PIPE
import sys
import os
from math import ceil

from argument_parser import argument_parser

#_____________________________________________________________________________
if __name__ == "__main__":

    #command line argumets
    resubmit = False
    args = sys.argv
    while args != []:
        arg = args.pop(0)
        if arg == "-r": resubmit = True
        elif arg == "-c": config = args.pop(0)

    #get outputs directory from config used for production
    parser = argument_parser()
    parser.add_parameter("top")
    parser.add_parameter("add_input", list)
    parser.parse(config)
    qlist = parser.get("add_input")
    basedirList = []
    for runNumber in  qlist:
        basedirList.append(parser.get("top") + runNumber[0]    )

    #submitted jobs
    joblist = []
    jobdist = []
    for basedir in basedirList:
        for job in glob(basedir + "/sched/*_*.csh"):
            joblist.append( job.split("sched/sched")[1].split(".csh")[0] )
            jobdist.append(basedir)

    print "Submitted:", len(joblist)

    #running jobs
    running = []
    cmd = "condor_q -global | "+os.getlogin()
    out = Popen(cmd.split(), stdout=PIPE).communicate()[0].split("\n")
    for i in out:
        i1 = i.find("sched/sched")+len("sched/sched")
        i2 = i.find(".csh")
        if i2 < 0: continue
        jobid = i[i1:i2]
        #select jobs belonging to this production
        if jobid not in joblist: continue
        running.append(jobid)

    print "Running:", len(running)
    #done jobs
    donelist = []
    totsiz = 0
    #list all root files with size
    for basedir in basedirList:
        cmd = "ls -s " + basedir + "/*.root"
        out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()[0].split("\n")
        for fline in out:
            if len(fline) == 0: continue
            #get name and size for each file
            siznam = fline.lstrip().split(" ")
            size = int(siznam[0])
            name = siznam[1]
            #test for zero-sized outputs
            if size == 0:
                print "Size: " + fline
                continue
            #test for segmentation violation in logs
            cmd = "grep \"segmentation violation\" " + basedir + "/logs/" + name.split("/")[-1][:-len(".root")] + ".out"
            grepTest = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()[0].split("\n")
            if len(grepTest) != 1:
                print "Grep: " + fline
                continue
            #job id from file name
            totsiz += size
            donelist.append(name.split("/")[-1][:-len(".root")])

    #missing jobs
    missing = []
    for job in joblist:
        if job in running or job in donelist:
            continue
        missing.append(job)
        print "Missing: " + job

    print "Errors:", len(missing)

    print "Done:", len(donelist)
    print "Output total size:".ljust(20), totsiz, "K"
    print "".ljust(20), ceil(float(totsiz)/1024), "M"
    print "".ljust(20), ceil(float(totsiz)/1024**2), "G"

 
    exit()
    #resubmit missing jobs if requested
    if resubmit is True and len(missing) > 0:
        print "Resubmitting the missing jobs"
        ijob = 0
        for job in missing:
            print "Clearing the log files"
            log_path = jobdist[joblist.index(job)] + "/logs/"
            cmd = "rm " + log_path + job + ".out"
            out = Popen(cmd.split(), stdout=PIPE, stderr=PIPE).communicate()
            cmd = "rm " + log_path + job + ".err"
            out = Popen(cmd.split(), stdout=PIPE, stderr=PIPE).communicate()
            print "Job to resubmit:", job, "idx in list:", ijob
            session_num = job.split("_")
            cmd = "star-submit -r " + session_num[1] + " " + session_num[0] + ".session.xml"
            out = Popen(cmd.split(), stdout=PIPE, stderr=PIPE).communicate()
            print out[0]
            print out[1]
            ijob += 1
        print "Resubmittion done."



















