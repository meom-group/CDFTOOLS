PROGRAM cdfgeo_uv
  !!======================================================================
  !!                     ***  PROGRAM  cdfgeo_uv  ***
  !!=====================================================================
  !!  ** Purpose : Compute the ug and vg component of the geostrophic 
  !!               velocity from the SSH field
  !!
  !!  ** Method  : ug = -g/f * d(ssh)/dy
  !!               vg =  g/f * d(ssh)/dx
  !!
  !!  **  Note : ug is located on a V grid point
  !!             vg                 U grid point
  !!
  !!
  !! History : 2.1  : 02/2008  : J. Juanno    : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic. and bug fix
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class derived_fields
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jt     ! dummy loop index
  INTEGER(KIND=4)                           :: npiglo, npjglo ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt       ! size of the domain
  INTEGER(KIND=4)                           :: narg, iargc    ! browse line
  INTEGER(KIND=4)                           :: ijarg          ! browse line
  INTEGER(KIND=4)                           :: ncoutu         ! ncid for ugeo file
  INTEGER(KIND=4)                           :: ncoutv         ! ncid for vgeo file
  INTEGER(KIND=4)                           :: ierr           ! error status
  INTEGER(KIND=4)                           :: ioption=0      ! Option for C-grid interpolation
  INTEGER(KIND=4), DIMENSION(1)             :: ipk            ! levels of output vars
  INTEGER(KIND=4), DIMENSION(1)             :: id_varoutu     ! varid for ugeo
  INTEGER(KIND=4), DIMENSION(1)             :: id_varoutv     ! varid for vgeo

  REAL(KIND=4)                              :: grav           ! gravity
  REAL(KIND=4)                              :: ffu, ffv       ! coriolis param f at U and V point
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1u, e2v, ff   ! horiz metrics, coriolis (f-point)
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2u, e1v       ! horiz metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: glamu, gphiu   ! longitude latitude u-point
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: glamv, gphiv   ! longitude latitude v-point
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: un, vn         ! velocity components
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zun, zvn       ! velocity components
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsshn          ! ssh
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zwrk           ! working array for interpolation
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: umask, vmask   ! mask at u and v points

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim           ! time counter

  CHARACTER(LEN=256)                        :: cf_tfil        ! input file name
  CHARACTER(LEN=256)                        :: cf_uout='ugeo.nc' 
  CHARACTER(LEN=256)                        :: cf_vout='vgeo.nc'
  CHARACTER(LEN=256)                        :: cldum          ! dummy character variable
  CHARACTER(LEN=256)                        :: cl_global      ! global attribute

  TYPE(variable), DIMENSION(1)              :: stypvaru       ! attributes for ugeo
  TYPE(variable), DIMENSION(1)              :: stypvarv       ! attributes for vgeo

  LOGICAL                                   :: lchk           ! file existence flag
  LOGICAL                                   :: lnc4 = .FALSE. ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  grav = 9.81  ! gravity

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfgeo-uv -f T-file [-o UOUT-file VOUT-file ] [-ssh SSH-var] [-nc4] [-C option]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'         Compute the geostrophic velocity components from the gradient of the'
     PRINT *,'         SSH read in the input file. '
     PRINT *,'      '
     PRINT *,'         Without any -C option, the zonal component is located on a C-grid '
     PRINT *,'         V point, the meridional one is located on a C-grid U point. See the'
     PRINT *,'         use of the -C option in order to have (Ugeo, Vgeo) at (U,V) points on'
     PRINT *,'         the C-grid.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'      -f T-file : netcdf file with SSH (input).' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'      [-o UOUT-file VOUT-file]: specify the names of the output files.'
     PRINT *,'              Default are: ',TRIM(cf_uout),' ',TRIM(cf_vout),'.'
     PRINT *,'      [-nc4]: Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'              This option is effective only if cdftools are compiled with'
     PRINT *,'              a netcdf library supporting chunking and deflation.'
     PRINT *,'      [-C option]: Using this option, the output velocity component are at the'
     PRINT *,'              correct (U,V) points on the C-grid. Two options are available :'
     PRINT *,'           option = 1 : SSH is interpolated on the F point prior derivation.'
     PRINT *,'           option = 2 : Ugeo and Vgeo are interpolated on the C-grid after'
     PRINT *,'              derivation.'
     PRINT *,'              Both options should give very similar results...'
     PRINT *,'      [-ssh SSH-var ]: give the name of the SSH variable if not ',TRIM(cn_sossheig)
     PRINT *,'  '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cn_fhgr),' and ',TRIM(cn_fzgr)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       - netcdf file : ', TRIM(cf_uout) ,' (default)'
     PRINT *,'           variables : ', TRIM(cn_vozocrtx)
     PRINT *,'           Unless -C option is used : '
     PRINT *,'             *** CAUTION:  this variable is located on V-point ***'
     PRINT *,'       - netcdf file : ', TRIM(cf_vout) ,' (default)'
     PRINT *,'           variables : ', TRIM(cn_vomecrty)
     PRINT *,'           Unless -C option is used : '
     PRINT *,'             *** CAUTION:  this variable is located on U-point ***'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ('-f'   ) ; CALL getarg(ijarg, cf_tfil     ) ; ijarg = ijarg + 1 
        ! options
     CASE ('-o'   ) ; CALL getarg(ijarg, cf_uout     ) ; ijarg = ijarg + 1
        ;             CALL getarg(ijarg, cf_vout     ) ; ijarg = ijarg + 1
     CASE ('-nc4' ) ; lnc4 = .TRUE.
     CASE ('-C'   ) ; CALL getarg(ijarg, cldum       ) ; ijarg = ijarg + 1
        ;             READ(cldum, * ) ioption
     CASE ('-ssh' ) ; CALL getarg(ijarg, cn_sossheig ) ; ijarg = ijarg + 1 
     CASE DEFAULT   ; PRINT *, ' ERROR : ',TRIM(cldum),' :  unknown option.' ; STOP 99
     END SELECT
  ENDDO

  SELECT CASE ( ioption )
  CASE ( 0 )   ; PRINT *,' *** UGEO on V-point, VGEO on U-point ***' ; cl_global=' (Ugeo, Vgeo) are on (V,U) points of the C-grid'
  CASE ( 1 )   ; PRINT *,' *** Use SSH interpolation ***'            ; cl_global=' (Ugeo, Vgeo) are on (U,V) points of the C-grid (SSH interp)'
  CASE ( 2 )   ; PRINT *,' *** Use Ugeo Vgeo interpolation ***'      ; cl_global=' (Ugeo, Vgeo) are on (U,V) points of the C-grid (velocity interp)'
  CASE DEFAULT ;  PRINT *, ' +++ ERROR: -C can use only option 1 or 2 +++' ; STOP 99
  END SELECT
  
  lchk = chkfile(cn_fhgr)
  lchk = chkfile(cn_fzgr) .OR. lchk
  lchk = chkfile(cf_tfil) .OR. lchk
  IF ( lchk ) STOP 99 ! missing file

  npiglo = getdim(cf_tfil, cn_x)
  npjglo = getdim(cf_tfil, cn_y)
  npk    = getdim(cf_tfil, cn_z) 
  npt    = getdim(cf_tfil, cn_t) 

  PRINT *, ' NPIGLO= ', npiglo
  PRINT *, ' NPJGLO= ', npjglo
  PRINT *, ' NPK   = ', npk
  PRINT *, ' NPT   = ', npt

  ! Allocate the memory
  ALLOCATE ( e1u(npiglo,npjglo), e2v(npiglo,npjglo) )
  IF( ioption == 1 ) ALLOCATE ( e1v(npiglo,npjglo), e2u(npiglo,npjglo) )
  ALLOCATE ( ff(npiglo,npjglo), dtim(npt)  )
  ALLOCATE ( glamu(npiglo,npjglo), gphiu(npiglo,npjglo)  )
  ALLOCATE ( glamv(npiglo,npjglo), gphiv(npiglo,npjglo)  )
  ALLOCATE ( un(npiglo,npjglo), vn(npiglo,npjglo)  )
  IF ( ioption == 2 ) ALLOCATE ( zun(npiglo,npjglo), zvn(npiglo,npjglo)  )
  ALLOCATE ( zsshn(npiglo,npjglo) )
  ALLOCATE ( umask(npiglo,npjglo), vmask(npiglo,npjglo) )
  ALLOCATE ( zwrk (npiglo,npjglo) )

  ! Read the metrics from the mesh_hgr file
  e1u   = getvar(cn_fhgr, cn_ve1u,  1, npiglo, npjglo)
  e1v   = getvar(cn_fhgr, cn_ve1v,  1, npiglo, npjglo)
  e2u   = getvar(cn_fhgr, cn_ve2u,  1, npiglo, npjglo)
  e2v   = getvar(cn_fhgr, cn_ve2v,  1, npiglo, npjglo)
  ff    = getvar(cn_fhgr, cn_vff,   1, npiglo, npjglo) 

  glamu = getvar(cn_fhgr, cn_glamu, 1, npiglo, npjglo)
  gphiu = getvar(cn_fhgr, cn_gphiu, 1, npiglo, npjglo)
  glamv = getvar(cn_fhgr, cn_glamv, 1, npiglo, npjglo)
  gphiv = getvar(cn_fhgr, cn_gphiv, 1, npiglo, npjglo)

  CALL CreateOutputUV

  ! Read ssh
  DO jt=1,npt
     zwrk(:,:) = getvar(cf_tfil, cn_sossheig, 1, npiglo, npjglo, ktime=jt)
     IF ( ioption == 1 ) THEN
        PRINT *, ' *** interpolation of SSH ...'
        DO jj=1, npjglo -1
           DO ji=1, npiglo -1
              zsshn(ji,jj) = 0.25*( zwrk (ji,jj  ) +  zwrk (ji+1,jj  ) + &
                   &                zwrk (ji,jj+1) +  zwrk (ji+1,jj+1) )
           ENDDO
        ENDDO
     ELSE
        zsshn(:,:) = zwrk(:,:) 
     ENDIF

     IF ( jt == 1 ) THEN
        ! compute the masks (do not depend on time
        umask=0. ; vmask = 0
        DO jj = 1, npjglo 
           DO ji = 1, npiglo - 1
              umask(ji,jj) = zwrk(ji,jj)*zwrk(ji+1,jj)
              IF (umask(ji,jj) /= 0.) umask(ji,jj) = 1.
           END DO
        END DO

        DO jj = 1, npjglo - 1
           DO ji = 1, npiglo
              vmask(ji,jj) = zwrk(ji,jj)*zwrk(ji,jj+1)
              IF (vmask(ji,jj) /= 0.) vmask(ji,jj) = 1.
           END DO
        END DO
        ! e1u and e1v are modified to simplify the computation below
        ! note that geostrophy is not available near the equator ( f=0)
        IF ( ioption == 0 .OR. ioption == 2 ) THEN  ! SSH at T point
           DO jj=2, npjglo - 1
              DO ji=2, npiglo - 1
                 ffu = ff(ji,jj) + ff(ji,  jj-1)
                 IF ( ffu /= 0. ) THEN  ; e1u(ji,jj)= 2.* grav * umask(ji,jj) / ( ffu ) / e1u(ji,jj)
                 ELSE                   ; e1u(ji,jj)= 0.  ! spvalue
                 ENDIF

                 ffv = ff(ji,jj) + ff(ji-1,jj  )
                 IF ( ffv /= 0. ) THEN  ; e2v(ji,jj)= 2.* grav * vmask(ji,jj) / ( ffv ) / e2v(ji,jj)
                 ELSE                   ; e2v(ji,jj)= 0.  ! spvalue
                 ENDIF
              END DO
           END DO
        ELSE    ! SSH at F point
           DO jj=2, npjglo - 1
              DO ji=2, npiglo - 1
                 ffu = ff(ji,jj) + ff(ji,  jj-1)
                 IF ( ffu /= 0. ) THEN ; e2u(ji,jj)= 2.* grav * umask(ji,jj) / ( ffu ) / e2u(ji,jj)
                 ELSE                  ; e2u(ji,jj)= 0.  ! spvalue
                 ENDIF

                 ffv = ff(ji,jj) + ff(ji-1,jj  )
                 IF ( ffv /= 0. ) THEN ; e1v(ji,jj)= 2.* grav * vmask(ji,jj) / ( ffv ) / e1v(ji,jj)
                 ELSE                  ; e1v(ji,jj)= 0.  ! spvalue
                 ENDIF
              END DO
           ENDDO
        ENDIF
     END IF

     ! Calculation of geostrophic velocity :
     un(:,:) = 0.
     vn(:,:) = 0.

     IF ( ioption == 0 .OR. ioption == 2 ) THEN
        DO jj = 1,npjglo - 1
           DO ji = 1,npiglo -1
              vn(ji,jj) =   e1u(ji,jj) * ( zsshn(ji+1,jj  ) - zsshn(ji,jj) ) 
              un(ji,jj) = - e2v(ji,jj) * ( zsshn(ji  ,jj+1) - zsshn(ji,jj) ) 
           END DO
        END DO
     ELSE    ! SSH at F point
        DO jj = 2,npjglo 
           DO ji = 2,npiglo 
              vn(ji,jj) =   e1v(ji,jj) * ( zsshn(ji,jj) - zsshn(ji-1,jj) ) 
              un(ji,jj) = - e2u(ji,jj) * ( zsshn(ji,jj) - zsshn(ji,jj-1) ) 
           END DO
        END DO
     ENDIF

     IF ( ioption == 2 ) THEN ! interpolate ugeo, vgeo on (U,V) point
        PRINT *, ' *** interpolation of velocities ...'
        DO jj=2,npjglo -1
           DO ji=2, npiglo -1
              zun(ji,jj) = 0.25*( un(ji,jj) + un(ji, jj-1) + un(ji+1,jj) + un (ji+1, jj-1) )
              zvn(ji,jj) = 0.25*( vn(ji,jj) + vn(ji-1, jj) + vn(ji,jj+1) + vn (ji-1, jj+1) )
           ENDDO
        ENDDO
        un(:,:) = zun(:,:)
        vn(:,:) = zvn(:,:)
     ENDIF

     ! write un and vn  ...
     ierr = putvar(ncoutu, id_varoutu(1), un(:,:), 1, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncoutv, id_varoutv(1), vn(:,:), 1, npiglo, npjglo, ktime=jt)

  END DO  ! time loop

  ierr = closeout(ncoutu)
  ierr = closeout(ncoutv)

CONTAINS

  SUBROUTINE CreateOutputUV
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutputUV  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ipk(1)                        = 1
    stypvaru(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvaru(1)%cname             = TRIM(cn_vozocrtx)
    stypvaru(1)%cunits            = 'm/s'
    stypvaru(1)%rmissing_value    = 0.
    stypvaru(1)%valid_min         = 0.
    stypvaru(1)%valid_max         = 20.
    stypvaru(1)%clong_name        = 'Zonal_Geostrophic_Velocity'
    stypvaru(1)%cshort_name       = TRIM(cn_vozocrtx)
    stypvaru(1)%conline_operation = 'N/A'
    stypvaru(1)%caxis             = 'TYX'

    stypvarv(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvarv(1)%cname             = TRIM(cn_vomecrty)
    stypvarv(1)%cunits            = 'm/s'
    stypvarv(1)%rmissing_value    = 0.
    stypvarv(1)%valid_min         = 0.
    stypvarv(1)%valid_max         = 20.
    stypvarv(1)%clong_name        = 'Meridional_Geostrophic_Velocity'
    stypvarv(1)%cshort_name       = TRIM(cn_vomecrty)
    stypvarv(1)%conline_operation = 'N/A'
    stypvarv(1)%caxis             = 'TYX'

    ! create output filesets
    ncoutu = create      (cf_uout, cf_tfil,  npiglo, npjglo, 0         , ld_nc4=lnc4                     )
    ierr   = createvar   (ncoutu,  stypvaru, 1,      ipk,    id_varoutu, ld_nc4=lnc4, cdglobal=cl_global )
    IF ( ioption == 0 ) THEN
       ! U geo  ! @ V-point !
       ierr   = putheadervar(ncoutu,  cf_tfil,  npiglo, npjglo, 0, pnavlon=glamv, pnavlat=gphiv)
    ELSE
       ierr   = putheadervar(ncoutu,  cf_tfil,  npiglo, npjglo, 0, pnavlon=glamu, pnavlat=gphiu)
    ENDIF

    dtim = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr = putvar1d(ncoutu,  dtim,      npt, 'T')

    ! V geo  ! @ U-point !
    ncoutv = create      (cf_vout, cf_tfil,  npiglo, npjglo, 0         , ld_nc4=lnc4                     )
    ierr   = createvar   (ncoutv,  stypvarv, 1,      ipk,    id_varoutv, ld_nc4=lnc4, cdglobal=cl_global )
    IF ( ioption == 0 ) THEN
       ! V geo  ! @ U-point !
       ierr   = putheadervar(ncoutv,  cf_tfil,  npiglo, npjglo, 0, pnavlon=glamu, pnavlat=gphiu)
    ELSE
       ierr   = putheadervar(ncoutv,  cf_tfil,  npiglo, npjglo, 0, pnavlon=glamv, pnavlat=gphiv)
    ENDIF

    dtim = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr = putvar1d(ncoutv,  dtim,      npt, 'T')

  END SUBROUTINE CreateOutputUV

END PROGRAM cdfgeo_uv

