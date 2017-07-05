PROGRAM cdfdifmask
  !!======================================================================
  !!                     ***  PROGRAM  cdfdifmask  ***
  !!=====================================================================
  !!  ** Purpose : Build the difference between 2 mask files
  !!
  !!
  !! History : 2.1  : ??????   : ???          : Original code
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: jk, jvar               ! dummy loop index
  INTEGER(KIND=4)                            :: ierr                   ! working integer
  INTEGER(KIND=4)                            :: narg, iargc            ! browsing command line
  INTEGER(KIND=4)                            :: npiglo, npjglo, npk    ! size of the domain
  INTEGER(KIND=4)                            :: ncout                  ! ncid of output file
  INTEGER(KIND=4), DIMENSION(4)              :: ipk                    ! outptut variables : levels,
  INTEGER(KIND=4), DIMENSION(4)              :: id_varout              ! ncdf varid's

  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: zmask, zmask2          ! 2D mask at current level
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: tim                    ! dummy time variable

  CHARACTER(LEN=256)                         :: cf_out='mask_diff.nc'  ! Output file name
  CHARACTER(LEN=256)                         :: cf_msk1, cf_msk2       ! name of input files
  CHARACTER(LEN=256)                         :: cv_in                  ! variable name

  TYPE(variable), DIMENSION(4)               :: stypvar                ! data structure

  LOGICAL                                    :: lchk                   ! checking file existence
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfdifmask  mask1 mask2'
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the difference between 2 mask files.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       mask1, mask2 : model files to be compared.' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'       variables : tmask, umask, vmask, fmask'
     STOP
  ENDIF
  CALL getarg (1, cf_msk1)
  CALL getarg (2, cf_msk2)

  lchk =           chkfile ( cf_msk1 )
  lchk = lchk .OR. chkfile ( cf_msk2 )
  IF ( lchk ) STOP 99 ! missing file

  npiglo = getdim (cf_msk1, cn_x)
  npjglo = getdim (cf_msk1, cn_y)
  npk    = getdim (cf_msk1, 'z' )  ! mask file have a z depth dim instead of depth ...

  ipk(:) = npk
  stypvar(:)%cunits            = '1/0'
  stypvar(:)%rmissing_value    = 9999.
  stypvar(:)%valid_min         = 0.
  stypvar(:)%valid_max         = 1.
  stypvar(:)%conline_operation = 'N/A'
  stypvar(:)%caxis             = 'TZYX'
  stypvar(:)%cprecision        = 'by'

  stypvar(1)%cname='tmask' ;  stypvar(1)%clong_name='tmask' ;  stypvar(1)%cshort_name='tmask'
  stypvar(2)%cname='umask' ;  stypvar(2)%clong_name='umask' ;  stypvar(2)%cshort_name='umask'
  stypvar(3)%cname='vmask' ;  stypvar(3)%clong_name='vmask' ;  stypvar(3)%cshort_name='vmask'
  stypvar(4)%cname='fmask' ;  stypvar(4)%clong_name='fmask' ;  stypvar(4)%cshort_name='fmask'

  ncout = create      (cf_out, cf_msk1, npiglo, npjglo, npk,      cdep='z', cdepvar='nav_lev')
  ierr  = createvar   (ncout,  stypvar, 4,      ipk,    id_varout                            )
  ierr  = putheadervar(ncout,  cf_msk1, npiglo, npjglo, npk,      cdep='nav_lev'             )

  ALLOCATE (zmask(npiglo,npjglo), zmask2(npiglo,npjglo))

  DO jvar=1,4
     cv_in = stypvar(jvar)%cname
     PRINT *, ' making difference for ', TRIM(cv_in)
     DO jk=1, npk
        PRINT * ,'jk = ', jk
        zmask(:,:)  = getvar(cf_msk1, cv_in, jk, npiglo, npjglo)
        zmask2(:,:) = getvar(cf_msk2, cv_in, jk, npiglo, npjglo)
        zmask(:,:)  = zmask2(:,:) - zmask(:,:)
        ierr        = putvar(ncout, id_varout(jvar), zmask, jk, npiglo, npjglo)
     END DO  ! loop to next level
  END DO

  ALLOCATE ( tim(1) )
  tim(:) = 0.
  ierr   = putvar1d(ncout, tim, 1, 'T')
  ierr   = closeout(ncout)

END PROGRAM cdfdifmask
