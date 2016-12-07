PROGRAM cdfstatcoord
  !!======================================================================
  !!                     ***  PROGRAM  cdfstatcoord  ***
  !!=====================================================================
  !!  ** Purpose : Compute statistics about the grid metric versus latitude
  !!
  !!  ** Method  : bins e1 and e2 by latitudes and takes the mean value 
  !!               of each bin
  !!
  !! History : 2.1  : 07/2007  : J.M. Molines : Original code (T. Penduff idea)
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

  INTEGER(KIND=4)                           :: narg, iargc          ! browse lines
  INTEGER(KIND=4)                           :: npiglo, npjglo       ! size of the domain
  INTEGER(KIND=4)                           :: ngood                ! point counter

  REAL(KIND=4), PARAMETER                   :: pp_binsize=2.        ! bin size
  REAL(KIND=4), PARAMETER                   :: pp_latmin=-80.       ! minimum latitude
  REAL(KIND=4), PARAMETER                   :: pp_latmax=90.        ! maximum latitude
  REAL(KIND=4)                              :: rlat, rlat1, rlat2   ! working variables
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: e1, e2, gphi, zmask  ! metrics and mask

  REAL(KIND=8)                              :: de1mean, de2mean     ! mean value of horiz metrics

  CHARACTER(LEN=256)                        :: cf_coo, cf_msk       ! file names
  CHARACTER(LEN=256)                        :: cv_msk='tmask'       ! mask variable name

  LOGICAL, DIMENSION(:,:), ALLOCATABLE      :: lgood                ! flag for point selection
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg < 2 ) THEN
     PRINT *,' usage : cdfstatcoord COOR-file MSK-file [ MSK-var ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Computes and displays statistics about grid metrics vs latitude.'
     PRINT *,'       Bins e1 and e2 by latitude bins, and compute the mean of each bin.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       COOR-file : coordinates file with e1 e2 metrics' 
     PRINT *,'       MSK-file  : mask file '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [MSK-var] : mask variable name. Default is ', TRIM(cv_msk) 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none apart those requested on command line.' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       Standard output'
     STOP
  ENDIF

  CALL getarg (1, cf_coo)
  CALL getarg (2, cf_msk)
  IF ( narg == 3 ) CALL getarg(3, cv_msk)
  
  IF ( chkfile(cf_coo) .OR. chkfile(cf_msk) ) STOP ! missing files

  npiglo= getdim (cf_coo, cn_x)
  npjglo= getdim (cf_coo, cn_y)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo

  ALLOCATE ( e1(npiglo,npjglo)  , e2(npiglo,npjglo) )
  ALLOCATE ( gphi(npiglo,npjglo), zmask(npiglo,npjglo), lgood(npiglo,npjglo)  )

  ! read grid metrics and latitude
  e1   = getvar(cf_coo, cn_ve1t,  1, npiglo, npjglo)
  e2   = getvar(cf_coo, cn_ve2t,  1, npiglo, npjglo)
  gphi = getvar(cf_coo, cn_gphit, 1, npiglo, npjglo)
  ! read zmask (1)
  zmask = getvar(cf_msk, cv_msk,  1, npiglo, npjglo)

  rlat = pp_latmin + pp_binsize/2.
  DO WHILE ( rlat <= pp_latmax ) 
     rlat1 = rlat - pp_binsize/2. ; rlat2 = rlat + pp_binsize/2.
     lgood = .FALSE.
     WHERE ( rlat1 <= gphi .AND. gphi < rlat2  .AND. zmask /= 0 )  lgood=.TRUE.
     ngood = COUNT(lgood)
     IF ( ngood /= 0 ) THEN
        de1mean = SUM( e1, mask=lgood) / ngood  
        de2mean = SUM( e2, mask=lgood) / ngood
     ELSE
        de1mean = -999.
        de2mean = -999.
     ENDIF
     PRINT '(f8.3, 3f15.3,i8)', rlat, de1mean, de2mean ,de1mean/de2mean, ngood
     rlat = rlat + pp_binsize
  ENDDO

END PROGRAM cdfstatcoord
