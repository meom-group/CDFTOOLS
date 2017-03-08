#!/bin/ksh


mkdir -p MODIFIED

for key in $( cat noclass_lst.txt | awk '{print $1","$NF}' ) ; do
   echo $key
   file=${key%,*}
   class=${key#*,}
   echo $file  '--->' $class

cat $file | sed -e '/History/ a\
\ \ !!         :  4.0  : 03/2017  : J.M. Molines  ' -e "/!! Software governed/ a\
\ \ !! @class $class" -e '/CDFTOOLS_3.0/ c\
\ \ !! CDFTOOLS_4.0 , MEOM 2017 ' -e '/Copyright/ c\
\ \ !! Copyright (c) 2017, J.-M. Molines ' > MODIFIED/$file

done
