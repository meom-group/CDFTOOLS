PROGRAM cdfmhst
  !!======================================================================
  !!                     ***  PROGRAM  cdfmhst  ***
  !!=====================================================================
  !!  ** Purpose : Compute Meridional Heat Salt  Transport.
  !!
  !!  ** Method  : Starts from the mean VT, VS fields computed by cdfvT.
  !!               Zonal and vertical integration are performed for these
  !!               quantities. If a sub-basin mask is provided, then a
  !!               meridional H/S transoport is computed for each sub basin.
  !!
  !! History : 2.1  : 01/2005  : J.M. Molines  : Original code
  !!                : 04/2005  : A.M. Treguier : adaptation to regional config
  !!                : 04/2007  : J.M. Molines  : add netcdf output
  !!           3.0  : 05/2011  : J.M. Molines  : Doctor norm + Lic.
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

  INTEGER(KIND=4)                              :: jj, jk, jt       ! dummy loop index
  INTEGER(KIND=4)                              :: jbasins, jvar    ! dummy loop index
  INTEGER(KIND=4)                              :: narg, iargc      ! command line 
  INTEGER(KIND=4)                              :: ijarg            ! argument counter
  INTEGER(KIND=4)                              :: npiglo, npjglo   ! size of the domain
  INTEGER(KIND=4)                              :: npk, npt         ! size of the domain
  INTEGER(KIND=4)                              :: numouth = 10     ! logical unit for heat 
  INTEGER(KIND=4)                              :: numouts = 11     ! logical unit for salt
  INTEGER(KIND=4)                              :: npvar=1          ! number of variables type
  INTEGER(KIND=4)                              :: nbasins          ! number of basins
  INTEGER(KIND=4)                              :: ierr             ! error status
  INTEGER(KIND=4)                              :: ncout            ! ncid of output file
  INTEGER(KIND=4)                              :: ivar             ! variable index
  INTEGER(KIND=4), DIMENSION(2)                :: iloc             ! working array
  INTEGER(KIND=4), DIMENSION(:),  ALLOCATABLE  :: ipk, id_varout   ! for output variables

  REAL(KIND=4), PARAMETER                      :: pprau0 = 1000.   ! reference density
  REAL(KIND=4), PARAMETER                      :: pprcp  = 4000.   ! specific heat
  REAL(KIND=4), PARAMETER                      :: ppspval= 9999.99 ! missing value

  REAL(KIND=4), DIMENSION(1)                   :: gdep             ! dummy depth array
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: tim              ! time counter
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: e31d             ! 1D e3t for full step
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: zmask            ! mask
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e1v, e3v, gphiv  ! metrics and latitude
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: zvt, zvs         ! transport components
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rdumlon          ! dummy longitude = 0.
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rdumlat          ! latitude for i = north pole

  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dzonal_heat_glo  ! zonal integral
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dzonal_heat_atl  ! zonal integral 
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dzonal_heat_pac
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dzonal_heat_ind
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dzonal_heat_aus
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dzonal_heat_med
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dzonal_salt_glo
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dzonal_salt_atl
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dzonal_salt_pac
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dzonal_salt_ind
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dzonal_salt_aus
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dzonal_salt_med
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dmtrp            ! transport in PW ir kT/s
  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: dwkh, dtrph      ! working variables
  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: dtrps, dwks      ! working variables

  TYPE(variable), DIMENSION(:),    ALLOCATABLE :: stypvar          ! structure for attributes

  CHARACTER(LEN=256)                           :: cf_vtfil         ! input VT file name
  CHARACTER(LEN=256)                           :: cf_outh='zonal_heat_trp.dat'
  CHARACTER(LEN=256)                           :: cf_outs='zonal_salt_trp.dat'
  CHARACTER(LEN=256)                           :: cf_outnc='mhst.nc'
  CHARACTER(LEN=256)                           :: cv_zomht='zomht' ! MHT variable name
  CHARACTER(LEN=256)                           :: cv_zomst='zomst' ! MST variable name
  CHARACTER(LEN=256)                           :: cldum            ! dummy character variable
  CHARACTER(LEN=4),  DIMENSION(5)              :: cbasin=(/'_glo','_atl','_inp','_ind','_pac'/)
  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: cvarname         ! varname arrays

  LOGICAL                                      :: llglo = .FALSE.  ! flag for sub basin file
  LOGICAL                                      :: lchk  = .FALSE.  ! flag for missing files
  LOGICAL                                      :: lfull = .FALSE.  ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmhst  VT-file  [MST] [-full]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the meridional heat/salt transport as a function of '
     PRINT *,'       latitude. If the file ',TRIM(cn_fbasins),' is provided, the meridional '
     PRINT *,'       heat/salt transport for each sub-basin is also computed.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       VT-file  : netcdf file containing the mean value of the products' 
     PRINT *,'                  U.S, U.T, V.S and V.T (obtained with cdfvT).'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [MST ] : output flag for meridional salt transport on netcdf files.'
     PRINT *,'                If not specified, only the MHT is output.' 
     PRINT *,'       [-full ] : to be set for full step case.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ', TRIM(cn_fhgr),', ',TRIM(cn_fzgr),' and ',TRIM(cn_fmsk)
     PRINT *,'        If ',TRIM(cn_fbasins),' is also available, sub-basin meridional transports'
     PRINT *,'        are also computed.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       ASCII files : ', TRIM(cf_outh),' : Meridional Heat Transport'
     PRINT *,'                     ', TRIM(cf_outs),' : Meridional Salt Transport'
     PRINT *,'       netcdf file : ', TRIM(cf_outnc)
     PRINT *,'           variables : ( [... ] : MST option ) '
     PRINT *,'                       ', TRIM(cv_zomht),cbasin(1),' : Meridional Heat Transport (global)'
     PRINT *,'                     [ ', TRIM(cv_zomst),cbasin(1),' : Meridional Salt Transport (global) ] '
     PRINT *,'       If ',TRIM(cn_fbasins),' is available, per basin meridional transport '
     PRINT *,'       are also available:' 
              DO jbasins=2, 5
     PRINT *,'                       ', TRIM(cv_zomht),cbasin(jbasins),' : Meridional Heat Transport'
     PRINT *,'                     [ ', TRIM(cv_zomst),cbasin(jbasins),' : Meridional Salt Transport ]'
              END DO
     STOP
  ENDIF

  npvar = 1    ! default value ( no MST output)
  ijarg = 1
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum) ; ijarg = ijarg+1
     SELECT CASE ( cldum)
     CASE ( 'MST' )   ; npvar    =2
     CASE ( '-full' ) ; lfull    = .TRUE.
     CASE DEFAULT     ; cf_vtfil = cldum
     END SELECT
  END DO

  ! check for missing files
  lchk = lchk .OR. chkfile( cn_fhgr )
  lchk = lchk .OR. chkfile( cn_fzgr )
  lchk = lchk .OR. chkfile( cn_fmsk )
  lchk = lchk .OR. chkfile( cf_vtfil)
  IF ( lchk ) STOP ! missing files

  ! check for sub basin file and set appropriate variables
  nbasins = 1 
  IF ( .NOT. chkfile(cn_fbasins) ) THEN
     llglo   = .TRUE.
     nbasins = 5 
  ENDIF

  npiglo = getdim (cf_vtfil, cn_x)
  npjglo = getdim (cf_vtfil, cn_y)
  npk    = getdim (cf_vtfil, cn_z)
  npt    = getdim (cf_vtfil, cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ! Allocate arrays
  ALLOCATE ( tim(npt) )
  ALLOCATE ( dwkh(npiglo,npjglo), zmask(npiglo,npjglo), zvt(npiglo,npjglo) )
  ALLOCATE ( dwks(npiglo,npjglo), zvs(npiglo,npjglo) )
  ALLOCATE ( e1v(npiglo,npjglo), e3v(npiglo,npjglo), gphiv(npiglo,npjglo))
  ALLOCATE ( dtrph(npiglo,npjglo))
  ALLOCATE ( dtrps(npiglo,npjglo))
  ALLOCATE ( dzonal_heat_glo(npjglo), dzonal_heat_atl(npjglo), dzonal_heat_pac(npjglo) )
  ALLOCATE ( dzonal_heat_ind(npjglo), dzonal_heat_aus(npjglo), dzonal_heat_med(npjglo) )
  ALLOCATE ( dzonal_salt_glo(npjglo), dzonal_salt_atl(npjglo), dzonal_salt_pac(npjglo) )
  ALLOCATE ( dzonal_salt_ind(npjglo), dzonal_salt_aus(npjglo), dzonal_salt_med(npjglo) )
  ALLOCATE ( dmtrp(npjglo) )
  ALLOCATE ( rdumlon(1,npjglo), rdumlat(1,npjglo))

  IF ( lfull ) ALLOCATE ( e31d(npk) )

  e1v(:,:)   = getvar(cn_fhgr, cn_ve1v,  1, npiglo, npjglo)
  gphiv(:,:) = getvar(cn_fhgr, cn_gphiv, 1, npiglo, npjglo)
  gdep(:)    = 0.  ! dummy depth for netcdf output

  IF ( lfull ) e31d = getvare3(cn_fzgr, cn_ve3t, npk )

  iloc         = MAXLOC( gphiv )
  rdumlat(1,:) = gphiv(iloc(1),:)
  rdumlon(:,:) = 0.   ! set the dummy longitude to 0

  ! prepare output netcdf output file
  ! Allocate output variables
  ALLOCATE(stypvar(nbasins*npvar),  cvarname(nbasins*npvar) )
  ALLOCATE(    ipk(nbasins*npvar), id_varout(nbasins*npvar) )

  ipk(:)=1               ! all output variables have only 1 level !
  DO jbasins = 1,nbasins

     cvarname(jbasins)                  = TRIM(cv_zomht)//TRIM(cbasin(jbasins))
     stypvar(jbasins)%cname             = cvarname(jbasins)
     stypvar(jbasins)%cunits            = 'PW'
     stypvar(jbasins)%rmissing_value    = ppspval
     stypvar(jbasins)%valid_min         = -10.
     stypvar(jbasins)%valid_max         = 20
     stypvar(jbasins)%clong_name        = 'Meridional Heat Transport '//TRIM(cbasin(jbasins))
     stypvar(jbasins)%cshort_name       = cvarname(jbasins)
     stypvar(jbasins)%conline_operation = 'N/A'
     stypvar(jbasins)%caxis             = 'TY'

     IF ( npvar == 2 ) THEN
        ! MST
        ivar = nbasins+jbasins
        cvarname(ivar)                   = TRIM(cv_zomst)//TRIM(cbasin(jbasins))
        stypvar(ivar )%cname             = cvarname(ivar)
        stypvar(ivar )%cunits            = 'T/sec'
        stypvar(ivar )%rmissing_value    = ppspval
        stypvar(ivar )%valid_min         = -10.e9
        stypvar(ivar )%valid_max         = 20.e9
        stypvar(ivar )%clong_name        = 'Meridional Salt Transport '//TRIM(cbasin(jbasins))
        stypvar(ivar )%cshort_name       = cvarname(ivar)
        stypvar(ivar )%conline_operation = 'N/A'
        stypvar(ivar )%caxis             = 'TY'
     ENDIF
  END DO

  ! create output fileset
  ncout = create      (cf_outnc, cf_vtfil, 1,             npjglo, 1, cdep='depthv'                              )
  ierr  = createvar   (ncout,    stypvar,  nbasins*npvar, ipk,    id_varout                                     )
  ierr  = putheadervar(ncout,    cf_vtfil, 1,             npjglo, 1, pnavlon=rdumlon, pnavlat=rdumlat, pdep=gdep)

  tim   = getvar1d (cf_vtfil, cn_vtimec, npt     )
  ierr  = putvar1d (ncout,    tim,       npt, 'T')

  OPEN(numouth,FILE=cf_outh,FORM='FORMATTED', RECL=256)  ! to avoid wrapped line with ifort
  OPEN(numouts,FILE=cf_outs,FORM='FORMATTED', RECL=256)  ! to avoid wrapped line with ifort

  DO jt=1, npt
     dtrph(:,:) = 0.d0
     dtrps(:,:) = 0.d0
     DO jk = 1,npk
        PRINT *,'level ',jk
        ! Get temperature and salinity at jk
        zvt(:,:)= getvar(cf_vtfil, cn_vomevt, jk, npiglo, npjglo, ktime=jt)
        zvs(:,:)= getvar(cf_vtfil, cn_vomevs, jk, npiglo, npjglo, ktime=jt)

        ! get e3v at level jk
        IF ( lfull ) THEN
           e3v(:,:) = e31d(jk)
        ELSE
           e3v(:,:)  = getvar(cn_fzgr, 'e3v_ps', jk, npiglo, npjglo, ldiom=.TRUE.)
        ENDIF
        dwkh(:,:) = zvt(:,:)*e1v(:,:)*e3v(:,:)*1.d0
        dwks(:,:) = zvs(:,:)*e1v(:,:)*e3v(:,:)*1.d0

        ! integrates vertically 
        dtrph(:,:) = dtrph(:,:) + dwkh(:,:) * pprau0 * pprcp
        dtrps(:,:) = dtrps(:,:) + dwks(:,:)  

     END DO  ! loop to next level

     ! global 
     zmask(:,:) = getvar(cn_fmsk, 'vmask', 1, npiglo, npjglo)
     DO jj=1,npjglo
        dzonal_heat_glo(jj) = SUM( dtrph(2:npiglo-1,jj)*zmask(2:npiglo-1,jj) )
        dzonal_salt_glo(jj) = SUM( dtrps(2:npiglo-1,jj)*zmask(2:npiglo-1,jj) )
     END DO

     IF ( llglo ) THEN
        ! Zonal mean with mask
        ! Atlantic 
        zmask(:,:) = getvar(cn_fbasins, 'tmaskatl', 1, npiglo, npjglo)
        DO jj=1,npjglo
           dzonal_heat_atl(jj) = SUM( dtrph(:,jj)*zmask(:,jj) )
           dzonal_salt_atl(jj) = SUM( dtrps(:,jj)*zmask(:,jj) )
        END DO

        ! Pacific
        zmask(:,:) = getvar(cn_fbasins, 'tmaskpac', 1, npiglo, npjglo)
        DO jj=1,npjglo
           dzonal_heat_pac(jj) = SUM( dtrph(:,jj)*zmask(:,jj) )
           dzonal_salt_pac(jj) = SUM( dtrps(:,jj)*zmask(:,jj) )
        END DO

        ! Indian
        zmask(:,:) = getvar(cn_fbasins, 'tmaskind', 1, npiglo, npjglo)
        DO jj=1,npjglo
           dzonal_heat_ind(jj) = SUM( dtrph(:,jj)*zmask(:,jj) )
           dzonal_salt_ind(jj) = SUM( dtrps(:,jj)*zmask(:,jj) )
        END DO

        ! Austral
        dzonal_heat_aus = 0.d0
        dzonal_salt_aus = 0.d0
        !    zmask(:,:)=getvar(cn_fbasins,'tmaskant',1,npiglo,npjglo)
        !    DO jj=1,npjglo
        !       dzonal_heat_aus(jj)= SUM( dtrph(:,jj)*zmask(:,jj))
        !       dzonal_salt_aus(jj)= SUM( dtrps(:,jj)*zmask(:,jj))
        !    END DO

        !   ! Med
        dzonal_heat_med = 0.d0
        dzonal_salt_med = 0.d0

        !    zmask(:,:)=getvar(cn_fbasins,'tmaskmed',1,npiglo,npjglo)
        !    DO jj=1,npjglo
        !       dzonal_heat_med(jj)= SUM( dtrph(:,jj)*zmask(:,jj))
        !       dzonal_salt_med(jj)= SUM( dtrps(:,jj)*zmask(:,jj))
        !    END DO
     ENDIF


     DO jvar=1,npvar   !  MHT [ and MST ]  (1 or 2 )
        IF ( jvar == 1 ) THEN
           ! MHT
           ivar=1
           dmtrp(:) = dzonal_heat_glo(:)/1.d15                        ! GLO
           WHERE ( dmtrp == 0 ) dmtrp = ppspval
           ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), 1, 1, npjglo, ktime=jt)
           ivar=ivar+1
           IF ( nbasins == 5 ) THEN
              dmtrp(:) = dzonal_heat_atl(:)/1.d15                      ! ATL
              WHERE ( dmtrp == 0 ) dmtrp = ppspval         
              ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), 1, 1, npjglo, ktime=jt)
              ivar=ivar+1
              dmtrp(:) = (dzonal_heat_ind(:) + dzonal_heat_pac(:))/1.d15  ! INP
              WHERE ( dmtrp == 0 ) dmtrp = ppspval
              ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), 1, 1, npjglo, ktime=jt)
              ivar=ivar+1
              dmtrp(:) = dzonal_heat_ind(:)/1.d15                      ! IND
              WHERE ( dmtrp == 0 ) dmtrp = ppspval
              ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), 1, 1, npjglo, ktime=jt)
              ivar=ivar+1
              dmtrp(:) = dzonal_heat_pac(:)/1.d15                      ! PAC
              WHERE ( dmtrp == 0 ) dmtrp = ppspval
              ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), 1, 1, npjglo, ktime=jt)
              ivar=ivar+1
           ENDIF
        ELSE
           ! MST
           dmtrp(:) = dzonal_salt_glo(:)/1.d6                        ! GLO
           WHERE ( dmtrp == 0 ) dmtrp = ppspval
           ierr=putvar(ncout, id_varout(ivar), REAL(dmtrp), 1, 1, npjglo, ktime=jt)
           ivar=ivar+1
           IF ( nbasins == 5 ) THEN
              dmtrp(:) = dzonal_salt_atl(:)/1.d6                      ! ATL
              WHERE ( dmtrp == 0 ) dmtrp = ppspval
              ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), 1, 1, npjglo, ktime=jt)
              ivar=ivar+1
              dmtrp(:) = (dzonal_salt_ind(:) + dzonal_salt_pac(:))/1.d6  ! INP
              WHERE ( dmtrp == 0 ) dmtrp = ppspval
              ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), 1, 1, npjglo, ktime=jt)
              ivar=ivar+1
              dmtrp(:) = dzonal_salt_ind(:)/1.d6                      ! IND
              WHERE ( dmtrp == 0 ) dmtrp = ppspval
              ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), 1, 1, npjglo, ktime=jt)
              ivar=ivar+1
              dmtrp(:) = dzonal_salt_pac(:)/1.d6                      ! PAC
              WHERE ( dmtrp == 0 ) dmtrp = ppspval
              ierr=putvar(ncout, id_varout(ivar), REAL(dmtrp), 1, 1, npjglo, ktime=jt)
              ivar=ivar+1
           ENDIF
        ENDIF
     END DO

     WRITE(numouth,*)'! Zonal heat transport (integrated alon I-model coordinate) (in Pw)'
     IF ( llglo ) THEN
        WRITE(numouth,*)'! J        Global          Atlantic         Pacific          Indian           Mediteranean     Austral  '
        WRITE(numouth,*)' ! time : ', jt
        DO jj=npjglo, 1, -1
           WRITE(numouth,9000) jj, &
                rdumlat(1,jj),  dzonal_heat_glo(jj)/1d15 , &
                dzonal_heat_atl(jj)/1d15, &
                dzonal_heat_pac(jj)/1d15, &
                dzonal_heat_ind(jj)/1d15, &
                dzonal_heat_med(jj)/1d15, &
                dzonal_heat_aus(jj)/1d15
        END DO
     ELSE
        WRITE(numouth,*)'! J        Global        '
        WRITE(numouth,*)' ! time : ', jt
        DO jj=npjglo, 1, -1
           WRITE(numouth,9000) jj, &
                rdumlat(1,jj),  dzonal_heat_glo(jj)/1d15  
        END DO
     ENDIF
     !               
     WRITE(numouts,*)' ! Zonal salt transport (integrated alon I-model coordinate) (in 10^6 kg/s)'
     IF ( llglo ) THEN
        WRITE(numouts,*)' ! J        Global          Atlantic         Pacific          Indian           Mediteranean     Austral  '
        WRITE(numouts,*)' ! time : ', jt
        !               
        DO jj=npjglo, 1, -1
           WRITE(numouts,9001) jj, &
                rdumlat(1,jj),  dzonal_salt_glo(jj)/1d6 , &
                dzonal_salt_atl(jj)/1d6, &
                dzonal_salt_pac(jj)/1d6, &
                dzonal_salt_ind(jj)/1d6, &
                dzonal_salt_med(jj)/1d6, &
                dzonal_salt_aus(jj)/1d6
        END DO
     ELSE
        WRITE(numouts,*)' J        Global  '
        WRITE(numouts,*)' ! time : ', jt
        DO jj=npjglo, 1, -1
           WRITE(numouts,9001) jj, &
                rdumlat(1,jj),  dzonal_salt_glo(jj)/1d6  
        ENDDO
     ENDIF

  ENDDO  ! time loop
  ierr = closeout(ncout)
  CLOSE(numouth)
  CLOSE(numouts)

9000 FORMAT(I4,6(1x,f9.3,1x,f8.4))
9001 FORMAT(I4,6(1x,f9.2,1x,f9.3))

END PROGRAM cdfmhst
