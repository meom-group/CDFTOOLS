PROGRAM cdfdifmask
  !!======================================================================
  !!                     ***  PROGRAM  cdfdifmask  ***
  !!=====================================================================
  !!  ** Purpose : Build the difference between 2 mask files
  !!
  !!
  !! History : 2.1  : ??????   : ???          : Original code
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class file_operations
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: jk, jvar               ! dummy loop index
  INTEGER(KIND=4)                            :: ierr                   ! working integer
  INTEGER(KIND=4)                            :: narg, iargc, ijarg     ! browsing command line
  INTEGER(KIND=4)                            :: npiglo, npjglo, npk    ! size of the domain
  INTEGER(KIND=4)                            :: ncout                  ! ncid of output file
  INTEGER(KIND=4), DIMENSION(4)              :: ipk                    ! outptut variables : levels,
  INTEGER(KIND=4), DIMENSION(4)              :: id_varout              ! ncdf varid's

  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: zmask, zmask2          ! 2D mask at current level
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: tim                    ! dummy time variable

  CHARACTER(LEN=256)                         :: cf_out='mask_diff.nc'  ! Output file name
  CHARACTER(LEN=256)                         :: cf_msk1, cf_msk2       ! name of input files
  CHARACTER(LEN=256)                         :: cv_in                  ! variable name
  CHARACTER(LEN=256)                         :: cldum                  ! working char variable

  TYPE(variable), DIMENSION(4)               :: stypvar                ! data structure

  LOGICAL                                    :: lchk                   ! checking file existence
  LOGICAL                                    :: lnc4 = .FALSE.         ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfdifmask -m  mask1 mask2 [-o OUT-file] [-nc4]'
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the difference between 2 NEMO mask files.' 
     PRINT *,'       This difference is not easy to perform with nco.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -m mask1 mask2 : model mask files to be compared.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file ]: specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ]     : Use netcdf4 output with chunking and deflation level 1'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option is used.'
     PRINT *,'       variables : ',TRIM(cn_tmask),' ', TRIM(cn_umask),' ',TRIM(cn_vmask),' ',TRIM(cn_fmask)
     STOP
  ENDIF

  ijarg=1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-m'   ) ; CALL getarg(ijarg, cf_msk1) ; ijarg=ijarg+1
        ; CALL getarg(ijarg, cf_msk2) ; ijarg=ijarg+1
        ! options
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out ) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP
     END SELECT
  END DO

  lchk =           chkfile ( cf_msk1 )
  lchk = lchk .OR. chkfile ( cf_msk2 )
  IF ( lchk ) STOP ! missing file

  npiglo = getdim (cf_msk1, cn_x)
  npjglo = getdim (cf_msk1, cn_y)
  npk    = getdim (cf_msk1, 'z' )  ! mask file have a z depth dim instead of depth ...

  ALLOCATE (zmask(npiglo,npjglo), zmask2(npiglo,npjglo))
  ALLOCATE ( tim(1) )

  CALL CreateOutput

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

  ierr   = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    DO jvar = 1, 4
       stypvar(jvar)%ichunk       = (/npiglo,MAX(1,npjglo/30),1,1 /)
    ENDDO
    ipk(:) = npk
    stypvar(:)%cunits            = '1/0'
    stypvar(:)%rmissing_value    = 9999.
    stypvar(:)%valid_min         = 0.
    stypvar(:)%valid_max         = 1.
    stypvar(:)%conline_operation = 'N/A'
    stypvar(:)%caxis             = 'TZYX'
    stypvar(:)%cprecision        = 'by'

    stypvar(1)%cname=cn_tmask ;  stypvar(1)%clong_name=cn_tmask ;  stypvar(1)%cshort_name=cn_tmask
    stypvar(2)%cname=cn_umask ;  stypvar(2)%clong_name=cn_umask ;  stypvar(2)%cshort_name=cn_umask
    stypvar(3)%cname=cn_vmask ;  stypvar(3)%clong_name=cn_vmask ;  stypvar(3)%cshort_name=cn_vmask
    stypvar(4)%cname=cn_fmask ;  stypvar(4)%clong_name=cn_fmask ;  stypvar(4)%cshort_name=cn_fmask

    ncout = create      (cf_out, cf_msk1, npiglo, npjglo, npk,      cdep='z', cdepvar='nav_lev', ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar, 4,      ipk,    id_varout                            , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_msk1, npiglo, npjglo, npk,      cdep='nav_lev'                           )

    tim(:) = 0.
    ierr   = putvar1d(ncout, tim, 1, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfdifmask
