#!/bin/bash


cat $1 | cut -d '/' -f 2- > edited$1
#sed -i 's|home/starlib|root://xrdstar.rcf.bnl.gov:1095//home/starlib|g' edited$1
sed -i 's/::/\//g' edited$1
sed -i -e 's/^/\//' edited$1

cat edited$1 > tmp.list
cat proccessedMuDst.list >> tmp.list
sort tmp.list | uniq -d > duplicate.list
cat edited$1 >> duplicate.list
sort duplicate.list | uniq -u > muDstPart.list
rm tmp.list
rm duplicate.list