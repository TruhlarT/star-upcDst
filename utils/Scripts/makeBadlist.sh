#!/bin/bash
# A sample Bash script, by Trhlar

while IFS='' read -r line || [[ -n "$line" ]]; do
    cmd="/star/data01/pwg_tasks/upc02/productionRun17/Part1/18*/sched/sched${line}.list"
	cat $cmd >> badMuDst.list
done < "$1"


  
