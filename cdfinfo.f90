PROGRAM cdfinfo
  !!======================================================================
  !!                     ***  PROGRAM  cdfinfo  ***
  !!=====================================================================
  !!  ** Purpose : Give very basic informations for Netcdf File
  !!
  !!  ** Method  : to be improved
  !!
  !! History : 2.1  : 09/2010  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jvar                     ! dummy loop index
  INTEGER(KIND=4)                               :: ierr                     ! working integer
  INTEGER(KIND=4)                               :: narg, iargc              ! 
  INTEGER(KIND=4)                               :: npiglo, npjglo, npk ,npt ! size of the domain
  INTEGER(KIND=4)                               :: nvars                    ! Number of variables in a file

  CHARACTER(LEN=256)                            :: cf_in                    ! file name
  CHARACTER(LEN=256)                            :: cv_dep                   ! depth name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names                 ! array of var name
  
  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvar                  ! variable attributes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfinfo ''model cdf file'' '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Gives very basic information about the file given in arguments.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        model output file in netcdf.' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'        On standard ouput, gives the size of the domain, the depth '
     PRINT *,'        dimension name, the number of variables.'
     PRINT *,'      '
     STOP
  ENDIF

  CALL getarg (1, cf_in)
  IF ( chkfile(cf_in) ) STOP ! missing file

  npiglo = getdim (cf_in,cn_x)
  npjglo = getdim (cf_in,cn_y)
  npk    = getdim (cf_in,cn_z, cdtrue=cv_dep, kstatus=ierr)

  IF (ierr /= 0 ) THEN
     npk   = getdim (cf_in,'z',cdtrue=cv_dep,kstatus=ierr)
     IF (ierr /= 0 ) THEN
       npk   = getdim (cf_in,'sigma',cdtrue=cv_dep,kstatus=ierr)
        IF ( ierr /= 0 ) THEN 
          npk = getdim (cf_in,'nav_lev',cdtrue=cv_dep,kstatus=ierr)
            IF ( ierr /= 0 ) THEN 
              npk = getdim (cf_in,'levels',cdtrue=cv_dep,kstatus=ierr)
              IF ( ierr /= 0 ) THEN 
                PRINT *,' assume file with no depth'
                npk=0
              ENDIF
            ENDIF
        ENDIF
     ENDIF
  ENDIF

  npt    = getdim (cf_in,cn_t)

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  PRINT *,' Depth dimension name is ', TRIM(cv_dep)

  nvars = getnvar(cf_in)
  PRINT *,' nvars =', nvars

  ALLOCATE (cv_names(nvars)  )
  ALLOCATE (stypvar(nvars)  )

  ! get list of variable names 
  cv_names(:)=getvarname(cf_in, nvars, stypvar)

  DO jvar = 1, nvars
   PRINT *, 'variable# ',jvar,' is : ',TRIM(cv_names(jvar))
  END DO

END PROGRAM cdfinfo
