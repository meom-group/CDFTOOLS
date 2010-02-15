#!/bin/ksh
set -x
n=0
while (( $n < 6 )) ; do
 sleep 3
 n=$( ls -l OK? 2> /dev/null | wc -l )
done
