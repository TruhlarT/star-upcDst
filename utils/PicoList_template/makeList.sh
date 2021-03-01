#!/bin/bash
# A sample Bash script, by Trhlar


while IFS='' read -r line || [[ -n "$line" ]]; do
   runNumber=$(echo "$line" | cut -d '/' -f 10 )
	#echo $runNumber
	echo $line >> $runNumber.list	
done < "$1"
