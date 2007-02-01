PROGRAM cdfprofile
  !!---------------------------------------------------------------------
  !!               ***  PROGRAM cdfprofile  ***
  !!
  !!  **  Purpose: extract a verticcal profile from a CDFfile
  !!  
  !!  **  Method:  read (i,j) position of point to extract
  !!               read varname
  !!               print profile
  !!
  !!
  !! history :
  !!   Original :  J.M. Molines June 2005
  !!---------------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER :: narg, iargc, istatus
  INTEGER :: jk
  INTEGER :: ilook, jlook
  INTEGER :: npiglo, npjglo, npk

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: v2d
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: depth, profile

  CHARACTER(LEN=80) :: cdum, cfile, cvar, cdep

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg /= 4  ) THEN
     PRINT *,' Usage : cdfprofile  I J file varname '
     PRINT *,' Output on standard output'
     STOP
  ENDIF


  CALL getarg (1, cdum)
  READ(cdum,*) ilook
  CALL getarg (2, cdum)
  READ(cdum,*) jlook
  CALL getarg(3, cfile)
  CALL getarg(4, cvar)

  npiglo= getdim (cfile,'x')
  npjglo= getdim (cfile,'y')
  npk   = getdim (cfile,'depth',cdep)

  ! Allocate arrays
  ALLOCATE( v2d (npiglo,npjglo), depth(npk) ,profile(npk) )

  depth(:) = getvar1d(cfile,cdep,npk,istatus)

  DO jk=1,npk
     v2d (:,:)= getvar(cfile, cvar,  jk ,npiglo,npjglo)
     profile(jk) = v2d(ilook,jlook)
  END DO
  PRINT *, "FILE : ", TRIM(cfile)
  PRINT *, "    ", TRIM(cdep),"         ", TRIM(cvar),"(",ilook,",",jlook,")"
  DO jk=1, npk
     PRINT *, depth(jk), profile(jk)
  END DO

END PROGRAM cdfprofile
