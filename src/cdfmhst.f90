PROGRAM cdfmhst
  !!======================================================================
  !!                     ***  PROGRAM  cdfmhst  ***
  !!=====================================================================
  !!  ** Purpose : Compute Meridional Heat Salt  Transport.
  !!
  !!  ** Method  : Starts from the mean VT, VS fields computed by cdfvT.
  !!               Optionally, it can read V, T  and S files in separate files
  !!               Zonal and vertical integration are performed for these
  !!               quantities. If a sub-basin mask is provided, then a
  !!               meridional H/S transoport is computed for each sub basin.
  !!               Meridional H/S transports at different depths
  !!
  !! History : 2.1  : 01/2005  : J.M. Molines  : Original code
  !!                : 04/2005  : A.M. Treguier : adaptation to regional config
  !!                : 04/2007  : J.M. Molines  : add netcdf output
  !!           3.0  : 05/2011  : J.M. Molines  : Doctor norm + Lic.
  !!                  10/2012  : M.A. Balmaseda: opt separate V,T,S infiles
  !!                                             opt z dimension 
  !!                                             add inp0
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

  INTEGER(KIND=4)                              :: ji, jj, jk, jt   ! dummy loop index
  INTEGER(KIND=4)                              :: jbasins, jvar    ! dummy loop index
  INTEGER(KIND=4)                              :: it               ! time index for vvl
  INTEGER(KIND=4)                              :: narg, iargc      ! command line 
  INTEGER(KIND=4)                              :: ik0              ! working integer
  INTEGER(KIND=4)                              :: ijarg            ! argument counter
  INTEGER(KIND=4)                              :: npiglo, npjglo   ! size of the domain
  INTEGER(KIND=4)                              :: npk, npt, npko   ! size of the domain
  INTEGER(KIND=4)                              :: numouth = 10     ! logical unit for heat 
  INTEGER(KIND=4)                              :: numouts = 11     ! logical unit for salt
  INTEGER(KIND=4)                              :: npvar=1          ! number of variables type
  INTEGER(KIND=4)                              :: nbasins          ! number of basins
  INTEGER(KIND=4)                              :: nbasinso         ! basins in output
  INTEGER(KIND=4)                              :: ierr             ! error status
  INTEGER(KIND=4)                              :: ncout            ! ncid of output file
  INTEGER(KIND=4)                              :: ivar             ! variable index
  INTEGER(KIND=4), DIMENSION(2)                :: iloc             ! working array
  INTEGER(KIND=4), DIMENSION(:),  ALLOCATABLE  :: ipk, id_varout   ! for output variables

  REAL(KIND=4), PARAMETER                      :: pprau0 = 1000.   ! reference density
  REAL(KIND=4), PARAMETER                      :: pprcp  = 4000.   ! specific heat
  REAL(KIND=4), PARAMETER                      :: ppspval= 9999.99 ! missing value

  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: gdep             ! depth array
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: tim              ! time counter
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: e31d             ! 1D e3t for full step
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: zmask            ! mask
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e1v, e3v, gphiv  ! metrics and latitude
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: zvt, zvs         ! transport components
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: zv, zt, zs       ! v,t,s
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

  CHARACTER(LEN=256)                           :: cf_vtfil         ! input VT or V file name
  CHARACTER(LEN=256)                           :: cf_vfil          ! separate V files
  CHARACTER(LEN=256)                           :: cf_tfil          ! separate T,S files
  CHARACTER(LEN=256)                           :: cf_sfil          ! separate T,S files
  CHARACTER(LEN=256)                           :: cf_outh='zonal_heat_trp.dat'
  CHARACTER(LEN=256)                           :: cf_outs='zonal_salt_trp.dat'
  CHARACTER(LEN=256)                           :: cf_outnc='mhst.nc'
  CHARACTER(LEN=256)                           :: cv_zomht='zomht' ! MHT variable name
  CHARACTER(LEN=256)                           :: cv_zomst='zomst' ! MST variable name
  CHARACTER(LEN=256)                           :: cldum            ! dummy character variable
  CHARACTER(LEN=5),  DIMENSION(6)              :: cbasin=(/'_glo ','_atl ','_inp ','_ind ','_pac ','_inp0'/)
  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: cvarname         ! varname arrays
  CHARACTER(LEN=256)                           :: cldimension      ! varname arrays

  LOGICAL                                      :: llglo = .FALSE.  ! flag for sub basin file
  LOGICAL                                      :: lchk  = .FALSE.  ! flag for missing files
  LOGICAL                                      :: lfull = .FALSE.  ! flag for full step
  LOGICAL                                      :: lsepf = .FALSE.  ! flag for separate files
  LOGICAL                                      :: lzdim = .FALSE.  ! flag for saving cumulated MHT
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmhst  -vt VT-file | (-v V-file -t T-file [-s S-file]) [-MST] ...'
     PRINT *,'              ... [-b BASIN-mask] [-full] [-Zdim] [-o OUT-file] [-vvl] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the meridional heat/salt transport as a function of latitude.'
     PRINT *,'       If the basin mask-file ',TRIM(cn_fbasins),' is provided, the meridional'
     PRINT *,'       heat/salt transport for each sub-basin is also computed.'
     PRINT *,'      '
     PRINT *,'       This program is designed to read the time-average products V.T into'
     PRINT *,'       the corresponding VT-file (produced by cdfvT). If VT-files are not'
     PRINT *,'       available, or if working with model snap-shot, V-file and T-file can'
     PRINT *,'       be specified as an alternative.'
     PRINT *,'      '
     PRINT *,'       In the evaluation of the ''meridional'' component, a ''zonal'' integration'
     PRINT *,'       is performed. Note that in CDFTOOLS, as far as ''zonal'' integrals'
     PRINT *,'       are concerned, they in-fact correspond to ''along-model I'' integrals.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -vt VT-file  : netcdf file containing the mean value of the products' 
     PRINT *,'                  U.S, U.T, V.S and V.T (obtained with cdfvT).'
     PRINT *,'         or   '     
     PRINT *,'       -v V-file  : specify V-file '
     PRINT *,'       -t T-file  : specify T-file '
     PRINT *,'       If not using ''-vt'' option, both V-file and T-file MUST be specified.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-s S-file] : specify S-file (salinity) if the salinity field is not in'
     PRINT *,'              T-file.'
     PRINT *,'       [-MST ]   : Indicates the the meridional Salt transport will also be '
     PRINT *,'              computed.  If not specified, only the MHT is output.' 
     PRINT *,'       [-full ] : to be set for full step case.'
     PRINT *,'       [-Zdim ] : to be set to output vertical structure of Heat/salt transport'
     PRINT *,'       [-b BASIN-mask ] : use BASIN-mask instead of default', TRIM(cf_outnc)
     PRINT *,'       [-o OUT-file ] : change name of the output file. Default:', TRIM(cf_outnc)
     PRINT *,'       [-vvl ] : use time-varying vertical metrics.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ', TRIM(cn_fhgr),', ',TRIM(cn_fzgr),' and ',TRIM(cn_fmsk)
     PRINT *,'        If ',TRIM(cn_fbasins),' is also available, sub-basin meridional transports'
     PRINT *,'        are also computed.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       ASCII files : ', TRIM(cf_outh),' : Meridional Heat Transport'
     PRINT *,'                     ', TRIM(cf_outs),' : Meridional Salt Transport'
     PRINT *,'       netcdf file : ', TRIM(cf_outnc),' unless -o option is used.'
     PRINT *,'           variables : ( [... ] : MST option ) '
     PRINT *,'                       ', TRIM(cv_zomht),cbasin(1),' : Meridional Heat Transport (global)'
     PRINT *,'                     [ ', TRIM(cv_zomst),cbasin(1),' : Meridional Salt Transport (global) ] '
     PRINT *,'       If ',TRIM(cn_fbasins),' is available, per basin meridional transport '
     PRINT *,'       are also available:' 
     DO jbasins=2, 6
        PRINT *,'                       ', TRIM(cv_zomht),cbasin(jbasins),' : Meridional Heat Transport'
        PRINT *,'                     [ ', TRIM(cv_zomst),cbasin(jbasins),' : Meridional Salt Transport ]'
     PRINT *,'      '
     END DO
     STOP
  ENDIF

  cf_vtfil='none'
  cf_vfil='none'
  cf_tfil='none'
  cf_sfil='none'
  lsepf  = .TRUE.
  npvar   = 1    ! default value ( no MST output)
  ijarg   = 1
  ! browse command line and detect the file name as argument different from any option
  ! count the number of files. Assume VT if 1 only, V TS if 2 and V T S if 3
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum) ; ijarg = ijarg+1
     SELECT CASE ( cldum)
     CASE ( '-vt'   ) ; CALL getarg(ijarg, cf_vtfil  ) ; ijarg = ijarg+1
     CASE ( '-v'    ) ; CALL getarg(ijarg, cf_vfil   ) ; ijarg = ijarg+1
     CASE ( '-t'    ) ; CALL getarg(ijarg, cf_tfil   ) ; ijarg = ijarg+1
     CASE ( '-s'    ) ; CALL getarg(ijarg, cf_sfil   ) ; ijarg = ijarg+1
     CASE ( 'MST' )   ; npvar    = 2
     CASE ( '-full' ) ; lfull    = .TRUE.
     CASE ( '-Zdim' ) ; lzdim    = .TRUE.
     CASE ( '-b'    ) ; CALL getarg(ijarg, cn_fbasins) ; ijarg = ijarg+1
     CASE ( '-o'    ) ; CALL getarg(ijarg, cf_outnc  ) ; ijarg = ijarg+1
     CASE ( '-vvl'  ) ; lg_vvl   = .TRUE.
     CASE DEFAULT     ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP
     END SELECT
  END DO

  ! security check
  IF ( TRIM(cf_vtfil) /= 'none' ) lsepf = .FALSE.

  IF ( lsepf ) THEN  ! need to have cf_vfil and cf_tfil at least
     lchk = chkfile ( cf_vfil )
     lchk = lchk .OR. chkfile( cf_tfil )
     IF ( lchk ) STOP 'Missing V-file of T-file '
     IF ( TRIM(cf_sfil) == 'none' ) cf_sfil = cf_tfil
  ENDIF

  ! check for missing files
  lchk = lchk .OR. chkfile( cn_fhgr )
  lchk = lchk .OR. chkfile( cn_fzgr )
  lchk = lchk .OR. chkfile( cn_fmsk )
  lchk = lchk .OR. chkfile( cf_vtfil)
  IF ( lsepf ) THEN
     lchk = lchk .OR. chkfile( cf_tfil)
     lchk = lchk .OR. chkfile( cf_sfil)
  ENDIF

  IF ( lchk ) STOP ! missing files
  IF ( lg_vvl ) cn_fe3v = cf_vtfil  ! REM: in case of vvl e3v is either in VT file or in V file 

  ! check for sub basin file and set appropriate variables
  IF ( .NOT. chkfile(cn_fbasins ) ) THEN
     llglo    = .TRUE.
     nbasins  = 5
     nbasinso = 6
  ELSE
     PRINT *,' Only compute for GLOBAL basin'
     llglo    = .FALSE.
     nbasins  = 1
     nbasinso = 1
  ENDIF

  npiglo = getdim (cf_vtfil, cn_x)
  npjglo = getdim (cf_vtfil, cn_y)
  npk    = getdim (cf_vtfil, cn_z)
  npt    = getdim (cf_vtfil, cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  IF ( lzdim ) THEN  ! output file hold heat/salt transport by levels
     npko = npk
     cldimension = 'TZY'
  ELSE
     npko = 1
     cldimension = 'TY'
  ENDIF

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
  ALLOCATE ( rdumlon(1,npjglo), rdumlat(1,npjglo))
  ALLOCATE ( dmtrp(npjglo) )
  ALLOCATE ( gdep (npko) )

  IF ( lsepf ) THEN
     ALLOCATE ( zv(npiglo,npjglo), zt(npiglo,npjglo), zs(npiglo,npjglo) )
  ENDIF

  IF ( lfull ) ALLOCATE ( e31d(npk) )

  e1v(:,:)   = getvar(cn_fhgr, cn_ve1v,  1, npiglo, npjglo)
  gphiv(:,:) = getvar(cn_fhgr, cn_gphiv, 1, npiglo, npjglo)

  IF ( lzdim ) THEN
     gdep(:)   = getvare3(cn_fzgr, cn_gdept, npk )
  ELSE
     gdep(:)    = 0.  ! dummy depth for netcdf output
  ENDIF

  IF ( lfull ) e31d = getvare3(cn_fzgr, cn_ve3t, npk )

  iloc         = MAXLOC( gphiv )
  rdumlat(1,:) = gphiv(iloc(1),:)
  rdumlon(:,:) = 0.   ! set the dummy longitude to 0

  ! prepare output netcdf output file
  ! Allocate output variables
  ALLOCATE(stypvar(nbasinso*npvar),  cvarname(nbasinso*npvar) )
  ALLOCATE(    ipk(nbasinso*npvar), id_varout(nbasinso*npvar) )

  CALL CreateOutput

  OPEN(numouth,FILE=cf_outh,FORM='FORMATTED', RECL=256)  ! to avoid wrapped line with ifort
  OPEN(numouts,FILE=cf_outs,FORM='FORMATTED', RECL=256)  ! to avoid wrapped line with ifort


  DO jt=1, npt
     IF ( lg_vvl ) THEN ; it=jt
     ELSE               ; it=1
     ENDIF
     dtrph(:,:) = 0.d0
     dtrps(:,:) = 0.d0
     DO jk = 1,npk
        PRINT *,'level ',jk
        ! Get temperature and salinity at jk
        IF ( lsepf ) THEN
           zv(:,:)= getvar(cf_vfil, cn_vomecrty, jk, npiglo, npjglo, ktime=jt)
           zt(:,:)= getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
           zs(:,:)= getvar(cf_sfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)
           zvt(:,:)=0.
           zvs(:,:)=0.
           DO ji=1, npiglo
              DO jj = 1, npjglo -1
                 zvt(ji,jj) = 0.5 * ( zt(ji,jj) + zt(ji,jj+1) )*zv(ji,jj)  ! temper at Vpoint
                 zvs(ji,jj) = 0.5 * ( zs(ji,jj) + zs(ji,jj+1) )*zv(ji,jj)  ! salinity at Vpoint
              END DO
           END DO

        ELSE
           zvt(:,:)= getvar(cf_vtfil, cn_vomevt, jk, npiglo, npjglo, ktime=jt)
           zvs(:,:)= getvar(cf_vtfil, cn_vomevs, jk, npiglo, npjglo, ktime=jt)
        ENDIF
        ! get e3v at level jk
        IF ( lfull ) THEN ; e3v(:,:) = e31d(jk)
        ELSE              ; e3v(:,:) = getvar(cn_fe3v, cn_ve3v, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl )
        ENDIF
        dwkh(:,:) = zvt(:,:)*e1v(:,:)*e3v(:,:)*1.d0
        dwks(:,:) = zvs(:,:)*e1v(:,:)*e3v(:,:)*1.d0

        ! integrates vertically 
        dtrph(:,:) = dtrph(:,:) + dwkh(:,:) * pprau0 * pprcp
        dtrps(:,:) = dtrps(:,:) + dwks(:,:)  

        !global
        zmask(:,:) = getvar(cn_fmsk, cn_vmask, 1, npiglo, npjglo)
        DO jj=1,npjglo
           dzonal_heat_glo(jj) = SUM( dtrph(2:npiglo-1,jj)*zmask(2:npiglo-1,jj) )
           dzonal_salt_glo(jj) = SUM( dtrps(2:npiglo-1,jj)*zmask(2:npiglo-1,jj) )
        END DO

        IF ( llglo ) THEN
           ! Zonal mean with mask
           ! Atlantic 
           zmask(:,:) = getvar(cn_fbasins, cn_tmaskatl, 1, npiglo, npjglo)
           DO jj=1,npjglo
              dzonal_heat_atl(jj) = SUM( dtrph(:,jj)*zmask(:,jj) )
              dzonal_salt_atl(jj) = SUM( dtrps(:,jj)*zmask(:,jj) )
           END DO

           ! Pacific
           zmask(:,:) = getvar(cn_fbasins, cn_tmaskpac, 1, npiglo, npjglo)
           DO jj=1,npjglo
              dzonal_heat_pac(jj) = SUM( dtrph(:,jj)*zmask(:,jj) )
              dzonal_salt_pac(jj) = SUM( dtrps(:,jj)*zmask(:,jj) )
           END DO

           ! Indian
           zmask(:,:) = getvar(cn_fbasins, cn_tmaskind, 1, npiglo, npjglo)
           DO jj=1,npjglo
              dzonal_heat_ind(jj) = SUM( dtrph(:,jj)*zmask(:,jj) )
              dzonal_salt_ind(jj) = SUM( dtrps(:,jj)*zmask(:,jj) )
           END DO

           ! Austral
           dzonal_heat_aus = 0.d0
           dzonal_salt_aus = 0.d0
           !    zmask(:,:)=getvar(cn_fbasins,cn_tmaskant,1,npiglo,npjglo)
           !    DO jj=1,npjglo
           !       dzonal_heat_aus(jj)= SUM( dtrph(:,jj)*zmask(:,jj))
           !       dzonal_salt_aus(jj)= SUM( dtrps(:,jj)*zmask(:,jj))
           !    END DO

           !   ! Med
           dzonal_heat_med = 0.d0
           dzonal_salt_med = 0.d0

           !    zmask(:,:)=getvar(cn_fbasins,cn_tmaskmed,1,npiglo,npjglo)
           !    DO jj=1,npjglo
           !       dzonal_heat_med(jj)= SUM( dtrph(:,jj)*zmask(:,jj))
           !       dzonal_salt_med(jj)= SUM( dtrps(:,jj)*zmask(:,jj))
           !    END DO
        ENDIF

        IF (  lzdim  .OR. ( jk == npk ) ) THEN  !output this level
           IF ( lzdim ) THEN ; ik0 = jk   
           ELSE              ; ik0 = 1 
           ENDIF

           DO jvar=1,npvar   !  MHT [ and MST ]  (1 or 2 )
              IF ( jvar == 1 ) THEN
                 ! MHT
                 ivar=1
                 dmtrp(:) = dzonal_heat_glo(:)/1.d15                            ! GLO
                 WHERE ( dmtrp == 0 ) dmtrp = ppspval
                 ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), ik0, 1, npjglo, ktime=jt) ; ivar=ivar+1
                 IF ( nbasins == 5 ) THEN
                    dmtrp(:) = dzonal_heat_atl(:)/1.d15                         ! ATL
                    WHERE ( dmtrp == 0 ) dmtrp = ppspval         
                    ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), ik0, 1, npjglo, ktime=jt); ivar=ivar+1
                    dmtrp(:) = (dzonal_heat_ind(:) + dzonal_heat_pac(:))/1.d15  ! INP
                    WHERE ( dmtrp == 0 ) dmtrp = ppspval
                    ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), ik0, 1, npjglo, ktime=jt); ivar=ivar+1
                    dmtrp(:) = dzonal_heat_ind(:)/1.d15                         ! IND
                    WHERE ( dmtrp == 0 ) dmtrp = ppspval
                    ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), ik0, 1, npjglo, ktime=jt); ivar=ivar+1
                    dmtrp(:) = dzonal_heat_pac(:)/1.d15                         ! PAC
                    WHERE ( dmtrp == 0 ) dmtrp = ppspval
                    ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), ik0, 1, npjglo, ktime=jt); ivar=ivar+1
                    ! now inp0
                    dmtrp(:) = ( dzonal_heat_glo(:) - dzonal_heat_atl(:) )/1.d15 ! INP0
                    WHERE ( dmtrp == 0 ) dmtrp = ppspval
                    ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), ik0, 1, npjglo, ktime=jt); ivar=ivar+1
                 ENDIF
              ELSE
                 ! MST
                 dmtrp(:) = dzonal_salt_glo(:)/1.d6                              ! GLO
                 WHERE ( dmtrp == 0 ) dmtrp = ppspval
                 ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), ik0, 1, npjglo, ktime=jt); ivar=ivar+1
                 IF ( nbasins == 5 ) THEN
                    dmtrp(:) = dzonal_salt_atl(:)/1.d6                           ! ATL
                    WHERE ( dmtrp == 0 ) dmtrp = ppspval
                    ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), ik0, 1, npjglo, ktime=jt); ivar=ivar+1
                    dmtrp(:) = (dzonal_salt_ind(:) + dzonal_salt_pac(:))/1.d6    ! INP
                    WHERE ( dmtrp == 0 ) dmtrp = ppspval
                    ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), ik0, 1, npjglo, ktime=jt); ivar=ivar+1
                    dmtrp(:) = dzonal_salt_ind(:)/1.d6                           ! IND
                    WHERE ( dmtrp == 0 ) dmtrp = ppspval
                    ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), ik0, 1, npjglo, ktime=jt); ivar=ivar+1
                    dmtrp(:) = dzonal_salt_pac(:)/1.d6                           ! PAC
                    WHERE ( dmtrp == 0 ) dmtrp = ppspval
                    ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), ik0, 1, npjglo, ktime=jt); ivar=ivar+1
                    ! now inp0
                    dmtrp(:) = ( dzonal_salt_glo(:) - dzonal_salt_atl(:) )/1.d6  ! INP0  
                    WHERE ( dmtrp == 0 ) dmtrp = ppspval
                    ierr = putvar(ncout, id_varout(ivar), REAL(dmtrp), ik0, 1, npjglo, ktime=jt); ivar=ivar+1
                 ENDIF
              ENDIF
           END DO

        ENDIF   ! end loop on check for ouptut
     END DO  ! loop to next level


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
    ipk(:)=npko               ! all output variables either 1 or npko levels
    DO jbasins = 1,nbasinso
       cvarname(jbasins)                  = TRIM(cv_zomht)//TRIM(cbasin(jbasins))
       stypvar(jbasins)%cname             = cvarname(jbasins)
       stypvar(jbasins)%cunits            = 'PW'
       stypvar(jbasins)%rmissing_value    = ppspval
       stypvar(jbasins)%valid_min         = -10.
       stypvar(jbasins)%valid_max         = 20
       stypvar(jbasins)%clong_name        = 'Meridional Heat Transport '//TRIM(cbasin(jbasins))
       stypvar(jbasins)%cshort_name       = cvarname(jbasins)
       stypvar(jbasins)%conline_operation = 'N/A'
       stypvar(jbasins)%caxis             = cldimension

       IF ( npvar == 2 ) THEN
          ! MST
          ivar = nbasinso+jbasins
          cvarname(ivar)                   = TRIM(cv_zomst)//TRIM(cbasin(jbasins))
          stypvar(ivar )%cname             = cvarname(ivar)
          stypvar(ivar )%cunits            = 'T/sec'
          stypvar(ivar )%rmissing_value    = ppspval
          stypvar(ivar )%valid_min         = -10.e9
          stypvar(ivar )%valid_max         = 20.e9
          stypvar(ivar )%clong_name        = 'Meridional Salt Transport '//TRIM(cbasin(jbasins))
          stypvar(ivar )%cshort_name       = cvarname(ivar)
          stypvar(ivar )%conline_operation = 'N/A'
          stypvar(ivar )%caxis             = cldimension
       ENDIF
    END DO

    ! create output fileset
    ncout = create      (cf_outnc, cf_vtfil, 1,             npjglo, npko, cdep='depthv'                              )
    ierr  = createvar   (ncout,    stypvar,  nbasinso*npvar, ipk,    id_varout                                       )
    ierr  = putheadervar(ncout,    cf_vtfil, 1,             npjglo, npko, pnavlon=rdumlon, pnavlat=rdumlat, pdep=gdep)

    tim   = getvar1d (cf_vtfil, cn_vtimec, npt     ) 
    ierr  = putvar1d (ncout,    tim,       npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfmhst
