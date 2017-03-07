PROGRAM cdfprobe
  !!======================================================================
  !!                     ***  PROGRAM  cdfprobe  ***
  !!=====================================================================
  !!  ** Purpose : Display time series of a variable at a given point
  !!
  !!
  !! History : 2.1  : 12/2006  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class file_informations
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)    :: narg, iargc            ! browse line
  INTEGER(KIND=4)    :: iilook, ijlook, ilevel ! point to look at
  CHARACTER(LEN=256) :: cf_in, cldum , cv_in   ! file name  variable name
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfprobe IN-file ilook jlook cdfvar [level]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'      Display a 2 columns output time (in days), value.'  
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       IN-file : input file to look for' 
     PRINT *,'       ilook jlook : i,j position of the probe.'
     PRINT *,'       cdfvar : name of the cdf variabled to be displayed'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [level] : This optional last argument is used' 
     PRINT *,'               to specify a model level, instead of first.'
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
  CALL getarg(1, cf_in )
  CALL getarg(2, cldum ) ; READ(cldum,*) iilook
  CALL getarg(3, cldum ) ; READ(cldum,*) ijlook
  CALL getarg(4, cv_in ) 

  IF ( chkfile(cf_in) ) STOP ! missing file

  IF ( narg == 5 ) THEN
     CALL getarg(5, cldum) ;  READ(cldum,*) ilevel
     CALL gettimeseries(cf_in, cv_in, iilook, ijlook, klev=ilevel)
  ELSE
     CALL gettimeseries(cf_in, cv_in, iilook, ijlook             )
  ENDIF

END PROGRAM cdfprobe
