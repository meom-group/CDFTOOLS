#!/bin/ksh

#   $Rev$
#   $Date$
# suppose that we run this script in $CDFTOOLS/DOC
CDFTOOLS=../

grep subsection cdftools_user.tex | grep -v addcontent | grep underline | sed -e 's@\\subsection\*{\\underline{@@'  -e 's/:}}//' \
        -e 's/\\//g' | sort > list_man
here=$(pwd)
cd $CDFTOOLS
ls -1 *90 | sed -e 's/.f90//' | sort > $here/list_prog
cd $here

n=01
for f in $( cat list_prog ); do
  grep -q $f list_man
  if [ $? == 1 ] ; then
    printf "\n %02d  %s \t %s \n " $n $f 'missing in manual'
    n=$(( n + 1 ))
  fi
done
    printf "\n"

for f in $( cat list_man ); do
  grep -q $f list_prog
  if [ $? == 1 ] ; then
    printf "%s \t %s \n \n" $f 'missing in CDFTOOLS ??'
  fi
done

\rm -f list_prog list_man

