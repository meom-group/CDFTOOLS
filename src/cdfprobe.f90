PROGRAM cdfprobe
  !!======================================================================
  !!                     ***  PROGRAM  cdfprobe  ***
  !!=====================================================================
  !!  ** Purpose : Display time series of a variable at a given point
  !!
  !!
  !! History : 2.1  : 12/2006  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class file_informations
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)    :: narg, iargc, ijarg     ! browse line
  INTEGER(KIND=4)    :: iilook, ijlook, iklook ! point to look at
  CHARACTER(LEN=256) :: cf_in, cldum , cv_in   ! file name  variable name
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfprobe -f IN-file -v IN-var -i ilook -j jlook [-k klook]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'      Displays a series of pair of values (time, value) corresponding to the'
     PRINT *,'      IN-var variable in IN-file, at location (ilook, jlook,[klook]). The'
     PRINT *,'      standard output can be piped to a graphical tool such as ''graph'' to'
     PRINT *,'      easily plot the time evolution of IN-var.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file : input file to look for' 
     PRINT *,'       -i ilook   : i position of the probe.'
     PRINT *,'       -j jlook   : j position of the probe.'
     PRINT *,'       -v IN-var  : name of the cdf variabled to be displayed'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-k klook] : Use the probe at level klook, instead of the first level.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       2 columns ( time , value ) ASCII output on display'
     PRINT *,'       time are given in days since the begining of the run.'
     STOP 
  ENDIF

  ! Browse command line
  iklook=1
  ijarg=1
  DO WHILE ( ijarg <= narg)
     CALL getarg(ijarg,cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f' ) ; CALL getarg(ijarg,cf_in) ; ijarg=ijarg+1
     CASE ( '-v' ) ; CALL getarg(ijarg,cv_in) ; ijarg=ijarg+1
     CASE ( '-i' ) ; CALL getarg(ijarg,cldum) ; ijarg=ijarg+1 ; READ(cldum,*) iilook
     CASE ( '-j' ) ; CALL getarg(ijarg,cldum) ; ijarg=ijarg+1 ; READ(cldum,*) ijlook
     CASE ( '-k' ) ; CALL getarg(ijarg,cldum) ; ijarg=ijarg+1 ; READ(cldum,*) iklook
     CASE DEFAULT  ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( chkfile(cf_in) ) STOP 99 ! missing file

     CALL gettimeseries(cf_in, cv_in, iilook, ijlook, klev=iklook)

END PROGRAM cdfprobe
