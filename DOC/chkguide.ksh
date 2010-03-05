#!/bin/ksh

#   $Rev$
#   $Date$
CDFTOOLS=../

grep subsection cdftools_prog.tex | grep -v addcontent | grep underline | sed -e 's@\\subsection\*{\\underline{@@'  -e 's/}}//' \
        -e 's/\\//g' | sed -e 's/^ //' |  sed -e 's/^.*F/F/' |  sed -e 's/^.*S/S/' | sort > list_man

here=$(pwd)
cd $CDFTOOLS
grep FUNCTION cdfio.f90 | grep -v END | grep -v -e '!' | sed -e 's/^.*F/F/' | sort > tttmp
grep FUNCTION eos.f90 | grep -v END | grep -v -e '!' | sed -e 's/^.*F/F/' | sort >> tttmp
grep SUBROUTINE cdfio.f90 | grep -v END | grep -v -e '!' | sed -e 's/^.*S/S/' | sort >> tttmp
grep SUBROUTINE eos.f90 | grep -v END | grep -v -e '!' | sed -e 's/^.*S/S/' | sort >> tttmp
cat tttmp | sort > $here/list_prog  ; \rm tttmp
cd $here

 cat list_prog
n=01
for f in $( cat list_prog | awk '{ print $0 }' ); do
echo $f
#  g=$( echo $f | awk '{print $2}' )
#echo $g
#  grep -q $g list_man
#  if [ $? == 1 ] ; then
#    printf "\n %02d  %s \t %s \n " $n $f 'missing in manual'
#    n=$(( n + 1 ))
#  fi
done
    printf "\n"
exit

for f in $( cat list_man ); do
  grep -q $f list_prog
  if [ $? == 1 ] ; then
    printf "%s \t %s \n \n" $f 'missing in CDFTOOLS ??'
  fi
done

\rm -f list_prog list_man

