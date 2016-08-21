#!/bin/bash 
#
# Simple script to make log files, the args should be:
#
# 1) the program (i.e. jones, motzkin, kauffman)
# 2) the directory to write the log files into
# 3) the maximum value to run the code on.
#
# Should be called from inside the log directory 
set -e
if [ "$#" -ne 3 ]; then
  echo "You must enter exactly 3 arguments (program, directory, maximum)"
  exit 1
fi

for i in $(eval echo {1..$3});
do
  ../$1 $i -v >> $2/$1`printf %02d $i`.log
done
