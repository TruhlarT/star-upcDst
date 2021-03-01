#!/bin/bash
echo "Getting files using: get_file_list"

get_file_list.pl -delim "/" -keys node,path,filename -cond production=P18ih,trgsetupname=pp500_production_2017,filetype=daq_reco_MuDst,filename~st_rp,storage!=hpss -limit 0 > newMuDst.list

cat newMuDst.list | wc -l

sed -e 's/^/root:\/\//' newMuDst.list > editednewMuDst.list

python prepare_file_list.py editednewMuDst.list

cat doneNames.list > tmp1.list
cat processedNames.list >> tmp1.list

sort tmp1.list | uniq -d | wc -l
echo "0 Good, else bad...."
echo "if 0, add new part to proccessedMuDst.list and also names!"

rm newMuDst.list
rm editednewMuDst.list
rm tmp1.list
