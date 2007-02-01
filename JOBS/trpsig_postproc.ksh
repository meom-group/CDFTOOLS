#!/bin/ksh

#set -x

# This script is used to perform the computation of Density Class Transport after
# their computation on zahir with monthly mean
# (1) It first  reorganize the results on a section/year basis.
# (2) Then it computes the mean of 12 months for trpsig 
# (3) finally it translates the mean results in a txt file for each section

### (1) : reorganization of the results.
here=$( basename `pwd` )

if [ ${here##*-} != TRPSIG ] ; then 
   echo this script mus be run in a -TRPSIG directory.
   exit 1
fi

CONFIG=${here%-TRPSIG}

# Fonction reorganize : put the files in the right position
reorganize() {

for year in $( ls -d ???? ) ; do
 cd $year

 for f in *.bimg ; do
   tmp=${f#${CONFIG}_yy????m??_}
   dir=${tmp%_*.bimg}
   if [ ! -d ../$dir ] ; then
      mkdir ../$dir
   fi

   if [ ! -d ../$dir/$year ] ; then
      mkdir ../$dir/$year
   fi
  mv $f ../$dir/$year
 done

 if [ ! -d ../LST ] ; then  mkdir ../LST ; fi
 if [ ! -d ../TRPSIG ] ; then  mkdir ../TRPSIG ; fi

 mv *.lst ../LST
 mv *.txt ../TRPSIG
 cd ../
 rmdir $year

done

}

### (2) :  now compute the mean foreach  year and sections
#        we suppose that section name starts with 2 digits 01_ 02_ 03_ etc ...
mean() {
# section name are codes as 2 digit_Capitalized_Name (eg: 01_Denmark_Strait or 07_Bab_el_Mandeb )
for stnam in [0-9][0-9]_[A-Z]* ; do
  cd $stnam
  printf "%s"  "Working for station $stnam  "

  # note that bimgmoy4 and bimgcaltrans exec are in cdftools-2.0 (extension ...)
  for d in ???? ; do
    printf "%4d " $(( d ))
    cd $d
    bimgmoy4 ${CONFIG}_y*trpsig.bimg > /dev/null
    mv moy.bimg ${CONFIG}_y${d}_${stnam}_trpsig.bimg
   # (2.2) : translate results into txt file foreach section
    bimgcaltrans ${CONFIG}_y${d}_${stnam}_trpsig.bimg > $TMPDIR/${CONFIG}_y${d}_${stnam}_trpsig.txt
    
    \rm [mv]*bimg
   cd ../
  done
  printf "\n"
  cd ../
done

}

##################### Main script: call functions

reorganize


mean



