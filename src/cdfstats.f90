PROGRAM cdfstats
  !!======================================================================
  !!                     ***  PROGRAM  cdfstats  ***
  !!=====================================================================
  !!  ** Purpose : Compute RMS/CORREL between 2 files. 
  !!               Seasonal cycle removed
  !!
  !!  ** Method  :
  !!
  !! History : 2.1  : 2009     : M.A. Balmaseda : original code  from cdfrmsssh.f90
  !!         :  4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils,  ONLY : SetGlobalAtt
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class statistics
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: jt, jm                 ! dummy loop index
  INTEGER(KIND=4)                           :: narg, iargc, ijarg, ij ! browse command line
  INTEGER(KIND=4)                           :: npiglo, npjglo         ! size of the domain
  INTEGER(KIND=4)                           :: npk, nt                ! size of the domain
  INTEGER(KIND=4), PARAMETER                :: jpvar=4                ! number of output variables
  INTEGER(KIND=4), DIMENSION(jpvar)         :: ipk, id_varout

  TYPE(variable), DIMENSION(jpvar)          :: stypvar                ! structure for attribute

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: u, v                   ! input variables from file 1 and 2
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask, e1t, e2t        ! mask and metrics

  REAL(KIND=8), DIMENSION(1)                :: dtim                   ! time for output (dummy)
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_er, dl_uv           ! rms, correlation
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_sn, dl_sg           ! signal/noise signal ratios
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_du, dl_dv           !  variable anomaly
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_u2, dl_v2           ! quadratic sum
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_um, dl_vm           ! linear sum
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_area                ! cell areas
  REAL(KIND=8)                              :: dl_spma                ! total area of the ocean ( info only)
  REAL(KIND=8)                              :: dl_spmu, dl_spmv       ! global mean (info only)
  REAL(KIND=8)                              :: dl_fct, dl_fcts        ! scaling coefficients

  CHARACTER(LEN=256)                        :: cf_in, cf_ref          ! input and reference file names
  CHARACTER(LEN=256)                        :: cf_msk, cf_hgr         ! current mask and hgr file
  CHARACTER(LEN=256)                        :: cf_out = 'stats.nc'    ! output file
  CHARACTER(LEN=256)                        :: cglobal                ! Global attribute
  CHARACTER(LEN=256)                        :: cldum                  ! dummy string for arguments
  CHARACTER(LEN=20)                         :: cv_name1, cv_name2     ! variable name
  CHARACTER(LEN=2)                          :: cy                     ! (1 or 12 ) 

  INTEGER(KIND=4)                           :: ncy                    ! 1/12 for annual/seasonal statistics
  INTEGER(KIND=4)                           :: ncout                  ! ID of netcdf output file
  INTEGER(KIND=4)                           :: ierr                   ! error status for ncdf

  LOGICAL                                   :: lchk                   ! flag for checking missing files
  LOGICAL                                   :: lnc4 = .FALSE.         ! Use nc4 with chunking and deflation
  !!--------------------------------------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfstats -f IN-file -r REF-file -ncy ncy [-v1 VAR-name1] ...'
     PRINT *,'             ... [-v2 VAR-name2] [-m MSH-MSK-file ] [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'         This tool computes some statistics (rms, correlation, signal/noise'
     PRINT *,'         ratio and signal ratio [ratio of std deviation]) between two files.'
     PRINT *,'         In this tool, the files are supposed to hold monthly averages values,'
     PRINT *,'         for many years. Specifying ncy=12, allows to remove the seasonal'
     PRINT *,'         cycle of the data.'
     PRINT *,'      '
     PRINT *,'         This program was initially written for SSH statistics between model'
     PRINT *,'         output and AVISO files (default variable names are ',TRIM(cn_sossheig)
     PRINT *,'         for this reason ). It can now be used with any variables.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        -f IN-file  : Model data file.'
     PRINT *,'        -r REF-file : Reference data file (usually observation file) '
     PRINT *,'        -ncy ncy    : ncy can only be 1 or 12. If set to 12, annual cycle is '
     PRINT *,'                 removed  from the data.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        [-v1 VAR-name1] : Variable name in the model file. Default is ',TRIM(cn_sossheig),'.'
     PRINT *,'        [-v2 VAR-name2] : Variable name in the reference file. If not specified'
     PRINT *,'              assumes that it is the same than in the model file.'
     PRINT *,'        [-m MSH-MSK-file] : specify a mesh_mask file holding the tmaskutil'
     PRINT *,'              and the horizontal metrics. If this option is not used,'
     PRINT *,'              mask are taken in ',TRIM(cn_fmsk), ' and horizontal metric'
     PRINT *,'              is taken in ',TRIM(cn_fhgr)
     PRINT *,'        [-o OUT-file] : specify the output file name instead of ',TRIM(cf_out)
     PRINT *,'        [-nc4 ]  : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'              This option is effective only if cdftools are compiled with'
     PRINT *,'              a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ' , TRIM(cn_fmsk),' and ', TRIM(cn_fhgr)
     PRINT *,'           or mesh_mask file specified in -m option'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'        netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables are : ' 
     PRINT *,'               rms    : RMS between the input files'
     PRINT *,'               correl : CORREL between the input files'
     PRINT *,'               rrat   : Signal to noise ratio '
     PRINT *,'               srat   : Signal ratio (stdev ratio) '
     PRINT *,'      '
     STOP 
  ENDIF

  ! default values
  cf_msk   = cn_fmsk
  cf_hgr   = cn_fhgr
  cv_name1 = cn_sossheig  
  cv_name2 = cn_sossheig  
  CALL SetGlobalAtt( cglobal )  ! global attribute for history : command line

  ! Browse command line
  ijarg = 1 
  DO WHILE (ijarg <= narg )
     CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-f ' ) ; CALL getarg (ijarg, cf_in   ) ; ijarg=ijarg+1
     CASE ( '-r ' ) ; CALL getarg (ijarg, cf_ref  ) ; ijarg=ijarg+1
     CASE ( '-ncy') ; CALL getarg (ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) ncy
        ! option
     CASE ( '-m ' ) ; CALL getarg (ijarg, cf_msk  ) ; ijarg=ijarg+1
        ;             cf_hgr = cf_msk
     CASE ( '-v1' ) ; CALL getarg (ijarg, cv_name1) ; ijarg=ijarg+1 ; cv_name2=cv_name1
     CASE ( '-v2' ) ; CALL getarg (ijarg, cv_name2) ; ijarg=ijarg+1 
     CASE ( '-o ' ) ; CALL getarg (ijarg, cf_out  ) ; ijarg=ijarg+1
     CASE ( '-nc4') ; lnc4 = .TRUE.
     CASE DEFAULT   ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  ! Security check for files
  lchk = chkfile ( cf_in  )
  lchk = chkfile ( cf_ref ) .OR. lchk
  lchk = chkfile ( cf_msk ) .OR. lchk
  lchk = chkfile ( cf_hgr ) .OR. lchk
  IF (lchk ) STOP 99 ! missing files

  ! log arguments do far
  PRINT *,'IN-file   : ', TRIM(cf_in )
  PRINT *,'REF-file  : ', TRIM(cf_ref)
  PRINT *,'NCY       : ', ncy
  PRINT *,'VAR-name1 : ', TRIM(cv_name1)
  PRINT *,'VAR_name2 : ', TRIM(cv_name2)
  PRINT *,'MASK file : ', TRIM(cf_msk  )
  PRINT *,'HGR file  : ', TRIM(cf_hgr  )

  ! define domain size from IN-file
  npiglo = getdim (cf_in, cn_x )
  npjglo = getdim (cf_in, cn_y )
  nt     = getdim (cf_in, cn_t )  ! read time dimension
  npk    = 1

  PRINT *, 'NPIGLO =', npiglo
  PRINT *, 'NPJGLO =', npjglo
  PRINT *, 'NPK    =', npk

  ! Allocate arrays from domain size
  ALLOCATE(     u(npiglo,npjglo),   v(npiglo,npjglo)                     )
  ALLOCATE( tmask(npiglo,npjglo), e1t(npiglo,npjglo), e2t(npiglo,npjglo) )

  ALLOCATE( dl_er(npiglo,npjglo), dl_uv(npiglo,npjglo) )
  ALLOCATE( dl_sn(npiglo,npjglo), dl_sg(npiglo,npjglo) )
  ALLOCATE( dl_u2(npiglo,npjglo), dl_v2(npiglo,npjglo) )
  ALLOCATE( dl_du(npiglo,npjglo), dl_dv(npiglo,npjglo) )
  ALLOCATE( dl_um(npiglo,npjglo), dl_vm(npiglo,npjglo) )
  ALLOCATE( dl_area(npiglo,npjglo)                     )

  CALL CreateOutput

  ! Read mask and metrics  
  tmask= getvar(cf_msk, cn_tmaskutil, 1, npiglo, npjglo)
  e1t  = getvar(cf_hgr, cn_ve1t,     1, npiglo, npjglo)
  e2t  = getvar(cf_hgr, cn_ve2t,     1, npiglo, npjglo)

  dl_area(:,:) = tmask(:,:)*e1t(:,:)*e2t(:,:)*1.d0   ! masked cell area
  dl_spma      = SUM(dl_area)                        ! model ocean area
  dl_fct       = 1.d0/float(nt)
  dl_fcts      = 1.d0*float(ncy)*dl_fct

  PRINT *, 'dl_spma = ',dl_spma, SUM(tmask)


  PRINT *,' output file created'
  dl_er(:,:) = 0.d0   ! rms
  dl_uv(:,:) = 0.d0   ! correlation
  dl_u2(:,:) = 0.d0   ! variance var1
  dl_v2(:,:) = 0.d0   ! variance var2
  dl_sn(:,:) = 0.d0   ! signal to noise ratio (rms/sdv)
  dl_sg(:,:) = 0.d0   ! signal ratio  (sdv(1)/sdv(2))

  DO jm = 1, ncy      ! loop on month (ncy=12)  or no loop if annual file (ncy=1)
     dl_um(:,:) = 0.d0
     dl_vm(:,:) = 0.d0
     PRINT *,' computing mean for month ',jm
     DO jt = jm, nt, ncy
        u(:,:) = getvar(cf_in , cv_name1, 1, npiglo, npjglo, ktime=jt)
        v(:,:) = getvar(cf_ref, cv_name2, 1, npiglo, npjglo, ktime=jt)

        dl_um(:,:) = dl_um(:,:) + u(:,:)*tmask(:,:)*1.d0
        dl_vm(:,:) = dl_vm(:,:) + v(:,:)*tmask(:,:)*1.d0
     ENDDO

     dl_um(:,:) = dl_um(:,:)*dl_fcts
     dl_vm(:,:) = dl_vm(:,:)*dl_fcts

     PRINT *,'MIN MAX UM ',MINVAL(dl_um), MAXVAL(dl_um)
     PRINT *,'MIN MAX VM ',MINVAL(dl_vm), MAXVAL(dl_vm)
     PRINT *,'computing 2nd order statistics'
     DO jt = jm, nt, ncy
        u(:,:) = getvar(cf_in , cv_name1, 1, npiglo, npjglo, ktime = jt)
        v(:,:) = getvar(cf_ref, cv_name2, 1, npiglo, npjglo, ktime = jt)

        ! anomaly
        dl_du(:,:) = (u(:,:) - dl_um(:,:))*tmask(:,:)
        dl_dv(:,:) = (v(:,:) - dl_vm(:,:))*tmask(:,:)

        ! REM no used if print below commented out (jmm ?)
        !         dl_spmu = SUM(u(:,:)*dl_area(:,:))/dl_spma
        !         dl_spmv = SUM(v(:,:)*dl_area(:,:))/dl_spma
        !         PRINT *,' jt, dl_spmu, dl_spmv ', jt, dl_spmu,dl_spmv

        dl_u2(:,:)=dl_u2(:,:) +  dl_du(:,:)             *  dl_du(:,:)
        dl_v2(:,:)=dl_v2(:,:) +  dl_dv(:,:)             *  dl_dv(:,:)
        dl_er(:,:)=dl_er(:,:) + (dl_du(:,:)-dl_dv(:,:)) * (dl_du(:,:)-dl_dv(:,:))
        dl_uv(:,:)=dl_uv(:,:) +  dl_du(:,:)             *  dl_dv(:,:)
     ENDDO
  ENDDO    ! loop on month

  dl_u2(:,:) = dl_u2(:,:)*dl_fct
  dl_v2(:,:) = dl_v2(:,:)*dl_fct 
  dl_uv(:,:) = dl_uv(:,:)*dl_fct
  dl_er(:,:) = SQRT(dl_er(:,:)*dl_fct)

  WHERE (tmask(:,:) > 0 )  dl_uv(:,:) = dl_uv(:,:)/SQRT(dl_u2(:,:)*dl_v2(:,:))
  WHERE (tmask(:,:) > 0 )  dl_sn(:,:) = dl_er(:,:)/SQRT(dl_v2(:,:))
  WHERE (tmask(:,:) > 0 )  dl_sg(:,:) = SQRT(dl_u2(:,:)/dl_v2(:,:))

  ! some print on standard output
  PRINT *,'MIN MAX RMS          ', MINVAL(dl_er), MAXVAL(dl_er)
  PRINT *,'MIN MAX CORREL       ', MINVAL(dl_uv), MAXVAL(dl_uv)
  PRINT *,'MIN MAX SIGNAL/NOISE ', MINVAL(dl_sn), MAXVAL(dl_sn)
  PRINT *,'MIN MAX SIGNAL RATIO ', MINVAL(dl_sg), MAXVAL(dl_sg)

  ! output on NC file
  ierr = putvar(ncout, id_varout(1), REAL(dl_er), 1, npiglo, npjglo)
  ierr = putvar(ncout, id_varout(2), REAL(dl_uv), 1, npiglo, npjglo)
  ierr = putvar(ncout, id_varout(3), REAL(dl_sn), 1, npiglo, npjglo)
  ierr = putvar(ncout, id_varout(4), REAL(dl_sg), 1, npiglo, npjglo)

  ierr      = closeout(ncout             )

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
    ! prepare output file 
    !   common features to all variables
    ipk    (:)                   = 1
    stypvar(:)%conline_operation = 'N/A'
    stypvar(:)%caxis             = 'TYX'

    !   specific features
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = 'rms'
    stypvar(1)%cunits            = 'm'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = 0.
    stypvar(1)%valid_max         = 100.
    stypvar(1)%clong_name        = 'RMS_'//TRIM(cv_name1)//'_'//TRIM(cv_name2)//'_'//cy
    stypvar(1)%cshort_name       = 'rms'

    stypvar(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname             = 'correl'
    stypvar(2)%cunits            = 'ndim'
    stypvar(2)%rmissing_value    = 0.
    stypvar(2)%valid_min         = -1.
    stypvar(2)%valid_max         = 1.
    stypvar(2)%clong_name        = 'CORREL_'//TRIM(cv_name1)//'_'//TRIM(cv_name2)//'_'//cy
    stypvar(2)%cshort_name       = 'correl'

    stypvar(3)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(3)%cname             = 'rrat'
    stypvar(3)%cunits            = 'N/A'
    stypvar(3)%rmissing_value    = 0.
    stypvar(3)%valid_min         = 0.
    stypvar(3)%valid_max         = 100.
    stypvar(3)%clong_name        = 'RMS/signal_'//TRIM(cv_name1)//'_'//TRIM(cv_name2)//'_'//cy
    stypvar(3)%cshort_name       = 'rrat'

    stypvar(4)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(4)%cname             = 'srat'
    stypvar(4)%cunits            = 'N/A'
    stypvar(4)%rmissing_value    = 0.
    stypvar(4)%valid_min         = 0.
    stypvar(4)%valid_max         = 100.
    stypvar(4)%clong_name        = 'sdvm/sdvo_'//TRIM(cv_name1)//'_'//TRIM(cv_name2)//'_'//cy
    stypvar(4)%cshort_name       = 'srat'

    PRINT *,' creating output file'
    ncout = create   (cf_out, cf_in,   npiglo, npjglo, npk                              , ld_nc4=lnc4)
    ierr  = createvar(ncout,  stypvar, jpvar,  ipk,    id_varout, cdglobal=TRIM(cglobal), ld_nc4=lnc4)
    ierr  = putheadervar(ncout,  cf_in, npiglo, npjglo, npk, ld_xycoo=.TRUE. )

    dtim = 1.d0
    ierr = putvar1d(ncout,dtim,1,'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfstats
