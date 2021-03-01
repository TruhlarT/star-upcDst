#!/usr/bin/bash

#here is catalog command to get a list of muDst for P18ib PR data from rhicf

#catalog query
query="trgsetupname=pp500_production_2017,filetype=daq_reco_MuDst,filename~st_rp,storage=nfs"
nfiles=0

#target directory for the files
targetdir="/gpfs01/star/pwg/truhlar/"

#output for retrieval request
files_request="newMuDst.list"

#-- end of config --


cat /dev/null > $files_request

hpss_list=$(get_file_list.pl -keys path,filename -cond $query -delim "/" -limit $nfiles)

#loop over files
for ifile in $hpss_list
do
    echo $ifile $targetdir$ifile >> $files_request
done

#submit the request
#hpss_user.pl -f $files_request

