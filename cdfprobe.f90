PROGRAM cdfprobe
  !!----------------------------------------------------------------
  !!                 ***   Program cdfprobe   ***
  !!  
  !!    Purpose : display time series of a variable at a given point
  !!
  !!    Usage :  cdfprobe  ncfile  i j cdfvar 
  !!
  !!  history:
  !!      Original: J.M. Molines Dec. 2006 
  !!----------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  USE cdfio
  IMPLICIT NONE
  INTEGER :: narg, iargc
  INTEGER :: ilook, jlook, ilevel
  CHARACTER(LEN=80) :: cfile, cdum , cvar

  narg=iargc()
  IF ( narg /= 4  ) THEN
     PRINT *,' USAGE: cdfprobe cdf_file i j cdfvar [level]'
     PRINT *,'      Display a 2 columns output time(d) value '
     STOP 
  ENDIF

  ! Browse command line
  CALL getarg(1, cfile)
  CALL getarg(2, cdum) ; READ(cdum,*) ilook
  CALL getarg(3, cdum) ; READ(cdum,*) jlook
  CALL getarg(4, cvar) 
  IF ( narg == 5 ) THEN
     CALL getarg(5, cdum) ;  READ(cdum,*) ilevel
     CALL gettimeseries(cfile,cvar,ilook,jlook,klev=ilevel)
  ELSE
     CALL gettimeseries(cfile,cvar,ilook,jlook)
  ENDIF

END PROGRAM cdfprobe
