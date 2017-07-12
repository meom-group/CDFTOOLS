PROGRAM cdfvtrp
  !!======================================================================
  !!                     ***  PROGRAM  cdfvtrp  ***
  !!=====================================================================
  !!  ** Purpose : Compute verticaly integrated transport.
  !!
  !!  ** Method  : Read the velocity components, and computed the verticaly
  !!               averaged transport at each grid cell ( velocity location).
  !!
  !! History : 2.1  : 01/2005  : J.M. Molines : Original code
  !!                : 01/2008  : P. Mathiot for -lbathy option
  !!           3.0  : 05/2011  : J.M. Molines : Doctor norm + Lic., merge
  !!                             with cdftrp_bathy
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class transport
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: ji, jj, jk, jt       ! dummy loop index
  INTEGER(KIND=4)                            :: it                   ! time index
  INTEGER(KIND=4)                            :: ierr                 ! working integer
  INTEGER(KIND=4)                            :: narg, iargc, ijarg   ! command line 
  INTEGER(KIND=4)                            :: npiglo, npjglo       ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt             ! size of the domain
  INTEGER(KIND=4)                            :: ncout                ! ncid of output file
  INTEGER(KIND=4)                            :: nvarout = 2          ! number of output variables
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout       ! for variable output

  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: e31d                 ! e3t metrics (full step)
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: e1u, e1v             ! horizontal metrics
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: e2u, e2v             !  "            "
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: e3u, e3v             ! vertical metrics
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: tmask, hdepw         ! tmask and bathymetry
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zdhdx, zdhdy         ! bottom slope
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zalpha               ! angle of rotation
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zu, zv               ! velocity components

  REAL(KIND=8), DIMENSION(:),    ALLOCATABLE :: dtim                 ! time counter
  REAL(KIND=8), DIMENSION(:,:),  ALLOCATABLE :: dwku , dwkv          ! working arrays
  REAL(KIND=8), DIMENSION(:,:),  ALLOCATABLE :: dtrpu, dtrpv         ! barotropic transport 

  TYPE (variable), DIMENSION(:), ALLOCATABLE :: stypvar              ! structure for attribute

  CHARACTER(LEN=256)                         :: cf_ufil              ! input U- file
  CHARACTER(LEN=256)                         :: cf_vfil              ! input V- file
  CHARACTER(LEN=256)                         :: cf_out='trp.nc'      ! output file
  CHARACTER(LEN=256)                         :: cv_soastrp='soastrp' ! Along Slope TRansPort
  CHARACTER(LEN=256)                         :: cv_socstrp='socstrp' ! Cross Slope TRansPort
  CHARACTER(LEN=256)                         :: cldum                ! dummy character variable

  LOGICAL                                    :: lfull  = .FALSE.     ! flag for full step
  LOGICAL                                    :: lbathy = .FALSE.     ! flag for slope current
  LOGICAL                                    :: lchk   = .FALSE.     ! flag for missing files
  LOGICAL                                    :: lnc4   = .FALSE.     ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfvtrp -u U-file -v V-file [-full] [-bathy] [-vvl] ...'
     PRINT *,'               ... [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the vertically integrated transports at each grid cell.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -u U-file : netcdf gridU file' 
     PRINT *,'       -v V-file : netcdf gridV file' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-full ]  : To be used in case of full step configuration.'
     PRINT *,'                Default is partial steps.'
     PRINT *,'       [-bathy ] : When used, cdfvtrp also compute the along slope and cross'
     PRINT *,'                slope transport components.'
     PRINT *,'                Bathymetry is read from ',TRIM(cn_fzgr),' file.'
     PRINT *,'       [-vvl  ] : Use time-varying vertical metrics'
     PRINT *,'       [-o OUT-file  ] : specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ]     : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'                This option is effective only if cdftools are compiled with'
     PRINT *,'                a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cn_fhgr),' and ',TRIM(cn_fzgr)
     PRINT *,'        ',TRIM(cn_fmsk),' is required only with -bathy option.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'       variables : ' 
     PRINT *,'           ', TRIM(cn_sozoutrp),' : zonal transport.'
     PRINT *,'           ', TRIM(cn_somevtrp),' : meridional transport.'
     PRINT *,'          If option -bathy is used :'
     PRINT *,'           ', TRIM(cv_soastrp),' : along slope transport'
     PRINT *,'           ', TRIM(cv_socstrp),' : cross slope transport'
     PRINT *,'      '
     STOP 
  ENDIF

  ! scan command line and set flags
  ijarg = 1 
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum ) 
     CASE ('-u'     ) ; CALL getarg(ijarg, cf_ufil) ; ijarg=ijarg+1
     CASE ('-v'     ) ; CALL getarg(ijarg, cf_vfil) ; ijarg=ijarg+1
        ! options
     CASE ('-full'  ) ; lfull  = .TRUE.
     CASE ('-vvl'   ) ; lg_vvl = .TRUE.
     CASE ('-o'     ) ; CALL getarg(ijarg, cf_out) ; ijarg=ijarg+1
     CASE ('-nc4'   ) ; lnc4   = .TRUE.
     CASE ('-bathy' ) ; lbathy = .TRUE. ; nvarout = 4
     CASE DEFAULT     ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO
  
  ! file existence check
  lchk = lchk .OR. chkfile ( cn_fzgr )
  lchk = lchk .OR. chkfile ( cn_fhgr )
  lchk = lchk .OR. chkfile ( cf_ufil )
  lchk = lchk .OR. chkfile ( cf_vfil )
  IF ( lbathy ) lchk = lchk .OR. chkfile ( cn_fmsk )
  IF ( lchk ) STOP 99   ! missing files

  IF ( lg_vvl) THEN
     cn_fe3u = cf_ufil
     cn_fe3v = cf_vfil
     cn_ve3u = cn_ve3uvvl
     cn_ve3v = cn_ve3vvvl
  ENDIF

  ALLOCATE ( ipk(nvarout), id_varout(nvarout), stypvar(nvarout) )

  npiglo = getdim (cf_ufil, cn_x)
  npjglo = getdim (cf_ufil, cn_y)
  npk    = getdim (cf_ufil, cn_z)
  npt    = getdim (cf_ufil, cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ! Allocate arrays
  ALLOCATE ( e1v(npiglo,npjglo), e3v(npiglo,npjglo)    )
  ALLOCATE ( e2u(npiglo,npjglo), e3u(npiglo,npjglo)    )
  ALLOCATE ( zu(npiglo,npjglo), zv(npiglo,npjglo)      )
  ALLOCATE ( dwku(npiglo,npjglo), dwkv(npiglo,npjglo)  )
  ALLOCATE ( dtrpu(npiglo,npjglo), dtrpv(npiglo,npjglo))
  ALLOCATE ( e31d(npk), dtim(npt)                      )

  IF ( lbathy ) THEN ! allocate extra arrays
     ALLOCATE ( e1u(npiglo, npjglo), e2v(npiglo, npjglo))
     ALLOCATE ( tmask(npiglo,npjglo), hdepw(npiglo, npjglo) )
     ALLOCATE ( zdhdx(npiglo,npjglo), zdhdy(npiglo, npjglo) )
     ALLOCATE ( zalpha(npiglo,npjglo) )
  ENDIF

  CALL CreateOutput

  e1v(:,:) = getvar(cn_fhgr, cn_ve1v, 1, npiglo, npjglo)
  e2u(:,:) = getvar(cn_fhgr, cn_ve2u, 1, npiglo, npjglo)
  IF ( lfull )  e31d(:)  = getvare3(cn_fzgr, cn_ve3t1d, npk )

  IF ( lbathy ) THEN  ! read extra metrics
    e1u(:,:)   = getvar(cn_fhgr, cn_ve1u,  1, npiglo, npjglo)
    e2v(:,:)   = getvar(cn_fhgr, cn_ve2v,  1, npiglo, npjglo)
    tmask(:,:) = getvar(cn_fmsk, cn_tmask, 1, npiglo, npjglo)
    hdepw(:,:) = getvar(cn_fzgr, cn_hdepw, 1, npiglo, npjglo)
  ENDIF

  DO jt = 1, npt
     dtrpu(:,:)= 0.d0
     dtrpv(:,:)= 0.d0
     IF ( lg_vvl ) THEN ; it = jt
     ELSE ;               it = 1
     ENDIF

     DO jk = 1,npk
        PRINT *,'level ',jk
        ! Get velocities at jk
        zu(:,:)= getvar(cf_ufil, cn_vozocrtx, jk ,npiglo, npjglo, ktime=jt)
        zv(:,:)= getvar(cf_vfil, cn_vomecrty, jk ,npiglo, npjglo, ktime=jt)

        ! get e3v at level jk
        IF ( lfull ) THEN
           e3v(:,:) = e31d(jk)
           e3u(:,:) = e31d(jk)
        ELSE
           e3v(:,:) = getvar(cn_fe3v, cn_ve3v, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl)
           e3u(:,:) = getvar(cn_fe3u, cn_ve3u, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl)
        ENDIF
        dwku(:,:) = zu(:,:)*e2u(:,:)*e3u(:,:)*1.d0
        dwkv(:,:) = zv(:,:)*e1v(:,:)*e3v(:,:)*1.d0
        ! integrates vertically 
        dtrpu(:,:) = dtrpu(:,:) + dwku(:,:)
        dtrpv(:,:) = dtrpv(:,:) + dwkv(:,:)

     END DO  ! loop to next level

     ierr = putvar(ncout, id_varout(1) ,REAL(dtrpu(:,:)), 1, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(2) ,REAL(dtrpv(:,:)), 1, npiglo, npjglo, ktime=jt)

     IF ( lbathy ) THEN
       ! compute transport component at T point 
       dwku(:,:) = 0.d0        ! U direction
       DO jj=1, npjglo
          DO ji= 2,npiglo
            dwku(ji,jj) = 0.5 * ( dtrpu(ji,jj) + dtrpu(ji-1,jj) )
          ENDDO
          ! E-W periodicity :
          dwku(1,jj) = dwku(npiglo-1, jj)
       ENDDO
       dwkv(:,:) = 0.d0        ! V direction
       DO jj=2, npjglo
          DO ji= 1,npiglo
            dwkv(ji,jj) = 0.5 * ( dtrpv(ji,jj) + dtrpv(ji,jj-1) )
          ENDDO
       ENDDO

       ! compute bathymetric slope at T point (centered scheme)
       zdhdx = 0.e0            ! U direction
       DO jj=1,npjglo          
          DO ji=2, npiglo-1
            zdhdx(ji,jj) = ( hdepw(ji+1,jj) - hdepw(ji-1,jj)) / ( e1u(ji,jj) + e1u(ji-1,jj) ) * tmask(ji,jj)
          END DO
       END DO

       zdhdy = 0.e0            ! V direction
       DO jj=2,npjglo-1        
          DO ji=1, npiglo
            zdhdy(ji,jj) = ( hdepw(ji,jj+1) - hdepw(ji,jj-1)) / ( e2v(ji,jj) + e2v(ji,jj-1) ) * tmask(ji,jj)
          END DO
       END DO
      
       ! compute the angle between the bathymetric slope and model coordinates
       zalpha(:,:) = ATAN2( zdhdx, zdhdy ) * tmask(:,:)
 
       ! apply the rotation on the transport
       dtrpu(:,:) = (  dwku(:,:) * COS(zalpha) + dwkv(:,:)* SIN(zalpha) ) * tmask(:,:)
       dtrpv(:,:) = ( -dwku(:,:) * SIN(zalpha) + dwkv(:,:)* COS(zalpha) ) * tmask(:,:)
      
       ierr = putvar(ncout, id_varout(3) ,REAL(dtrpu(:,:)), 1, npiglo, npjglo, ktime=jt)
       ierr = putvar(ncout, id_varout(4) ,REAL(dtrpv(:,:)), 1, npiglo, npjglo, ktime=jt)
     ENDIF
  END DO

  ierr = closeout (ncout)
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
    ! define  variables for output 
    ipk(:) = 1   ! all 2D variables
    stypvar%rmissing_value    = 0.
    stypvar%valid_min         = -100.
    stypvar%valid_max         = 100.
    stypvar%cunits            = 'm3/s'
    stypvar%conline_operation = 'N/A'
    stypvar%caxis             = 'TYX'

    stypvar(1)%ichunk         = (/npiglo,MAX(1,npjglo/30),1,1 /) ; stypvar(2)%ichunk      = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname       = cn_sozoutrp                         ; stypvar(2)%cname       = cn_somevtrp
    stypvar(1)%clong_name  = 'Zonal_barotropic_transport'        ; stypvar(2)%clong_name  = 'Meridional_barotropic_transport'
    stypvar(1)%cshort_name = cn_sozoutrp                         ; stypvar(2)%cshort_name = cn_somevtrp

    IF ( lbathy ) THEN
       stypvar(3)%ichunk      = (/npiglo,MAX(1,npjglo/30),1,1 /)
       stypvar(3)%cname       = cv_soastrp                       ; stypvar(4)%cname       = cv_socstrp
       stypvar(3)%clong_name  = 'Along_Slope_Barotropic_Transp'  ; stypvar(4)%clong_name  = 'Cross_Slope_Barotropic_Transp'
       stypvar(3)%cshort_name = cv_soastrp                       ; stypvar(4)%cshort_name = cv_socstrp
    ENDIF

    ! create output fileset
    ncout = create      (cf_out, cf_ufil, npiglo,  npjglo, 1         , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar, nvarout, ipk,    id_varout , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_ufil, npiglo,  npjglo, 1                       )
  
    dtim  = getvar1d(cf_ufil, cn_vtimec, npt     )
    ierr  = putvar1d(ncout,   dtim,      npt, 'T')
  END SUBROUTINE CreateOutput

END PROGRAM cdfvtrp
