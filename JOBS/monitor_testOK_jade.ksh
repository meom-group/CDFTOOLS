#!/bin/ksh
set -x
n=0
while (( $n < $1 )) ; do
 sleep 3
 n=$( ls -l ????/OK_MONITOR 2> /dev/null | wc -l )
done
