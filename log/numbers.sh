#!/bin/bash 

# Simple script to extract the actual numbers produced that are in log files, 
# e.g. 
#
# log/numbers.sh log/laptop/jones 
#
# will create file log/laptop/jones.txt with the number in the jones* log 
# files in the directory log/laptop.

set -e
for file in $1*
do
   tail -1 $file >> $1.txt
done
