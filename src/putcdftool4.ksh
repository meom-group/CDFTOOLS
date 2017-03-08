#!/bin/ksh


mkdir -p MODIFIED

for file in *.?90 ; do 
cat $file | sed -e '/History/ a\
\ \ !!         :  4.0  : 03/2017  : J.M. Molines  ' -e '/CDFTOOLS_3.0/ c\
\ \ !! CDFTOOLS_4.0 , MEOM 2017 ' -e '/Copyright/ c\
\ \ !! Copyright (c) 2017, J.-M. Molines ' > MODIFIED/$file

done
