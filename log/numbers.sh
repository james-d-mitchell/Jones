#!/bin/bash 

# Simple script to extract the actual numbers produced that are in log files, 
# e.g. 
#
# log/numbers.sh log/lovelace/jones 
#
# will create file log/laptop/jones.txt with the number in the jones* log 
# files in the directory log/laptop.

set -e

if [ -f $1.txt ]; then
  rm -f $1.txt
fi
touch $1.txt
for file in $1*
do
   tail -1 $file >> $1.txt
done
