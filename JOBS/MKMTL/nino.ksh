#!/bin/ksh
dir=$( basename `pwd` )
CONFIG=${dir%-DIAGS}

# if no NINO.txt files in the dir, skip
 ls *NINO* 1> /dev/null 2>&1
 if [ $? != 0 ] ; then echo no nino to deal with ... ; exit ; fi

cat *NINO* > ../${CONFIG}-MONITOR/${CONFIG}_nino.mtl
