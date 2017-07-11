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
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class statistics
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: narg, iargc,ijarg    ! browse lines
  INTEGER(KIND=4)                           :: npiglo, npjglo       ! size of the domain
  INTEGER(KIND=4)                           :: ngood                ! point counter

  REAL(KIND=4), PARAMETER                   :: pp_binsize=2.        ! bin size
  REAL(KIND=4), PARAMETER                   :: pp_latmin=-80.       ! minimum latitude
  REAL(KIND=4), PARAMETER                   :: pp_latmax=90.        ! maximum latitude
  REAL(KIND=4)                              :: rlat, rlat1, rlat2   ! working variables
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: e1, e2, gphi, zmask  ! metrics and mask

  REAL(KIND=8)                              :: de1mean, de2mean     ! mean value of horiz metrics

  CHARACTER(LEN=256)                        :: cf_coo, cf_msk       ! file names
  CHARACTER(LEN=256)                        :: cv_msk               ! mask variable name
  CHARACTER(LEN=256)                        :: cldum                ! mask variable name

  LOGICAL, DIMENSION(:,:), ALLOCATABLE      :: lgood                ! flag for point selection
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  cv_msk = cn_tmask

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfstatcoord -c COOR-file -m MSK-file [-v MSK-var ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute and displays statistics about grid metrics vs latitude.'
     PRINT *,'       Bins e1 and e2 by latitude bins, and compute the mean of each bin.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -c COOR-file : coordinates file with e1 e2 metrics' 
     PRINT *,'       -m MSK-file  : mask file '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-v MSK-var] : mask variable name. Default is ', TRIM(cv_msk) 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none apart those requested on command line.' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       Standard output'
     STOP 
  ENDIF

  ijarg=1
  DO WHILE (ijarg <= narg )
     CALL getarg(ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum ) 
     CASE ( '-c' ) ; CALL getarg (ijarg, cf_coo) ; ijarg=ijarg+1
     CASE ( '-m' ) ; CALL getarg (ijarg, cf_msk) ; ijarg=ijarg+1
     CASE ( '-v' ) ; CALL getarg (ijarg, cv_msk) ; ijarg=ijarg+1
     CASE DEFAULT  ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( chkfile(cf_coo) .OR. chkfile(cf_msk) ) STOP 99 ! missing files

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
     PRINT '(a)', 'Latitude       e1mean         e2mean           e1/e2    npoints'
     PRINT '(a)', '-----------------------------------------------------------------'
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
