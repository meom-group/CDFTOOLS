PROGRAM cdfmocsig
  !!======================================================================
  !!                     ***  PROGRAM  cdfmocsig  ***
  !!=====================================================================
  !!  ** Purpose : Compute the Meridional Overturning Cell (MOC)
  !!               using density bins. 
  !!
  !!  ** Method  : The MOC is computed from the V velocity field, collected in density bins,
  !!               (reference depth is given as the 3rd argument) and integrated
  !!               throughout the density bins, then zonally averaged with
  !!               eventual masking for oceanic basins.
  !!               In the present version the masking corresponds to the global
  !!               configuration. MOC for Global, Atlantic, Indo-Pacific, Indian,Pacific ocean
  !!               Results are saved on mocsig.nc file with variables name respectively
  !!               zomsfglo, zomsfatl, zomsfinp, zomsfind, zomsfpac.
  !!               If no new_maskglo.nc file found, then the mask.nc file is used and
  !!               only zomsfglo is computed.

  !!
  !! History : 2.1  : 11/2005  : A.M. Treguier : Original code from cdfmoc
  !!                : 03/2010  : C. Dufour     : Choice of depth reference
  !!                                             improvements 
  !!           3.0  : 04/2011  : J.M. Molines  : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER(KIND=2), DIMENSION (:,:,:), ALLOCATABLE ::  ibmask              ! nbasins x npiglo x npjglo
  INTEGER(KIND=2), DIMENSION (:,:),   ALLOCATABLE ::  itmask              ! tmask from salinity field

  INTEGER(KIND=4)                                 :: jbasin, jj, jk       ! dummy loop index
  INTEGER(KIND=4)                                 :: ji, jt, jbin         ! dummy loop index
  INTEGER(KIND=4)                                 :: nbins                ! number of  density  bins
  INTEGER(KIND=4)                                 :: npglo, npatl, npinp  ! basins index (mnemonics)
  INTEGER(KIND=4)                                 :: npind, nppac         !  "      "
  INTEGER(KIND=4)                                 :: nbasins              ! number of basins
  INTEGER(KIND=4)                                 :: ierr                 ! working integer
  INTEGER(KIND=4)                                 :: narg, iargc, iarg    ! command line  browsing 
  INTEGER(KIND=4)                                 :: ijarg, ii            !  "             "
  INTEGER(KIND=4)                                 :: npiglo,npjglo        ! size of the domain
  INTEGER(KIND=4)                                 :: npk, npt             ! size of the domain
  INTEGER(KIND=4)                                 :: ncout                ! ncid of output file
  INTEGER(KIND=4), DIMENSION(2)                   :: iloc                 ! working array
  INTEGER(KIND=4), DIMENSION(:),      ALLOCATABLE :: ipk, id_varout       ! output variable levels and id
  INTEGER(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: ibin                 ! remaping density in bin number

  REAL(KIND=4), DIMENSION (:,:),      ALLOCATABLE :: e1v, gphiv           ! horizontal metrics, latitude
  REAL(KIND=4), DIMENSION (:,:),      ALLOCATABLE :: zt, zs               ! temperature, salinity
  REAL(KIND=4), DIMENSION (:,:),      ALLOCATABLE :: zv, zveiv            ! velocity and bolus velocity
  REAL(KIND=4), DIMENSION (:,:),      ALLOCATABLE :: e3v                  ! vertical metrics
  REAL(KIND=4), DIMENSION (:,:),      ALLOCATABLE :: rdumlon              ! dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:),      ALLOCATABLE :: rdumlat              ! latitude for i = north pole
  REAL(KIND=4), DIMENSION (:,:),      ALLOCATABLE :: zttmp                ! arrays to call sigmai and mask it
  REAL(KIND=4), DIMENSION (:),        ALLOCATABLE :: sigma                ! density coordinate 
  REAL(KIND=4), DIMENSION (:),        ALLOCATABLE :: e31d                 ! vertical level (full step)
  REAL(KIND=4), DIMENSION (:),        ALLOCATABLE :: tim                  ! time counter
  REAL(KIND=4)                                    :: pref=0.              ! depth reference for pot. density 
  REAL(KIND=4)                                    :: sigmin               ! minimum density for bining
  REAL(KIND=4)                                    :: sigstp               ! density step for bining

  REAL(KIND=8), DIMENSION(:,:,:),     ALLOCATABLE :: dmoc                 ! nbasins x npjglo x npk
  REAL(KIND=8), DIMENSION(:,:),       ALLOCATABLE :: dens                 ! density
  REAL(KIND=8), DIMENSION(:,:),       ALLOCATABLE :: dmoc_tmp             ! temporary transport array

  CHARACTER(LEN=256)                              :: cf_vfil              ! meridional velocity file
  CHARACTER(LEN=256)                              :: cf_tfil              ! temperature/salinity file
  CHARACTER(LEN=256)                              :: cf_moc='mocsig.nc'   ! output file
  CHARACTER(LEN=255)                              :: cglobal              ! Global attribute
  CHARACTER(LEN=256)                              :: cldum                ! dummy char variable

  TYPE(variable), DIMENSION(:), ALLOCATABLE       :: stypvar              ! output var properties

  LOGICAL, DIMENSION(3)                           :: lbin                 ! flag for bin specifications
  LOGICAL                                         :: lbas   = .FALSE.     ! flag for basins file
  LOGICAL                                         :: lprint = .FALSE.     ! flag for extra print
  LOGICAL                                         :: leiv   = .FALSE.     ! flag for Eddy Induced Velocity (GM)
  LOGICAL                                         :: lfull  = .FALSE.     ! flag for full step
  LOGICAL                                         :: lchk   = .FALSE.     ! flag for missing file
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmocsig  V_file T_file depth_ref [-eiv] [-full]  ... '
     PRINT *,'         ...  [-sigmin sigmin] [-sigstp sigstp] [-nbins nbins] [-v] '
     PRINT *,'     PURPOSE : '
     PRINT *,'       Computes the MOC in density-latitude coordinates. The global value'
     PRINT *,'       is always computed. Values for oceanic sub-basins are calculated'
     PRINT *,'       if the file ', TRIM(cn_fbasins), ' is provided.'
     PRINT *,'       Last arguments is the reference depth for potential density, in m.'
     PRINT *,'       Actually only 0 1000 or 2000 are available with standard values for'
     PRINT *,'       density bins. If you specify another reference depth, you must also'
     PRINT *,'       specify the minimum density, the bin size and the number of bins,'
     PRINT *,'       with the options -sigmin, -sigstp, -nbins'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        V_file  : Netcdf gridV file'
     PRINT *,'        T_file  : Netcdf gridT file'
     PRINT *,'        depth_ref : reference depth for density '
     PRINT *,'               for depth values of 0 1000 or 2000, pre-defined limits for'
     PRINT *,'               minimum density, number of density bins and width of density'
     PRINT *,'               bins are provided. For other reference depth, you must use'
     PRINT *,'               -sigmin, -sigstp and -nbins options (see below).'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-eiv ] : takes into account VEIV Meridional eddy induced velocity'
     PRINT *,'                 -> To be used only if Gent and McWilliams parameterization '
     PRINT *,'                    has been used '
     PRINT *,'       [ -full ] : Works with full step instead of standard partial steps'
     PRINT *,'       [ -sigmin ] : Specify minimum of density for bining'
     PRINT *,'       [ -sigstp ] : Specify density step for bining'
     PRINT *,'       [ -nbins ]  : Specify the number of density bins you want'
     PRINT *,'       [ -v  ]     : Verbose option for more info during execution'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        Files ', TRIM(cn_fzgr),', ',TRIM(cn_fhgr),', ', TRIM(cn_fmsk)
     PRINT *,'        File ', TRIM(cn_fbasins),' is optional [sub basins masks]'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_moc) 
     PRINT *,'       variables ',TRIM( cn_zomsfglo),' : Global ocean '
     PRINT *,'       variables ',TRIM( cn_zomsfatl),' : Atlantic Ocean '
     PRINT *,'       variables ',TRIM( cn_zomsfinp),' : Indo Pacific '
     PRINT *,'       variables ',TRIM( cn_zomsfind),' : Indian Ocean alone'
     PRINT *,'       variables ',TRIM( cn_zomsfpac),' : Pacific Ocean alone'
     PRINT *,'       If file ',TRIM(cn_fbasins),' is not present, ',TRIM(cn_fmsk),' file'
     PRINT *,'       is used and only ',TRIM( cn_zomsfglo),' is produced.'
     STOP
  ENDIF

  cglobal = 'Partial step computation'
  lbin=(/.TRUE.,.TRUE.,.TRUE./)
  ijarg = 1 ; ii = 0
  DO WHILE ( ijarg <= narg )
     CALL getarg (ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ('-full')
        lfull   = .TRUE.
        cglobal = 'Full step computation'
     CASE ('-eiv')
        leiv    = .TRUE.
     CASE ('-sigmin')
        CALL getarg (ijarg, cldum) ; ijarg=ijarg+1 ; READ(cldum,*) sigmin
        lbin(1) = .FALSE.
     CASE ('-nbins')
        CALL getarg (ijarg, cldum) ; ijarg=ijarg+1 ; READ(cldum,*) nbins
        lbin(2) = .FALSE.
     CASE ('-sigstp')
        CALL getarg (ijarg, cldum) ; ijarg=ijarg+1 ; READ(cldum,*) sigstp
        lbin(3) = .FALSE.
     CASE ('-v')
        lprint = .TRUE.
     CASE DEFAULT
        ii=ii+1
        SELECT CASE (ii)
        CASE ( 1 ) ; cf_vfil = cldum
        CASE ( 2 ) ; cf_tfil = cldum
        CASE ( 3 ) ; READ(cldum,*) pref
        CASE DEFAULT
           STOP 'ERROR : Too many arguments ...'
        END SELECT
     END SELECT
  END DO

  ! check file existence
  lchk = lchk .OR. chkfile ( cn_fhgr )
  lchk = lchk .OR. chkfile ( cn_fzgr )
  lchk = lchk .OR. chkfile ( cn_fmsk )
  lchk = lchk .OR. chkfile ( cf_vfil )
  lchk = lchk .OR. chkfile ( cf_tfil )
  IF ( lchk ) STOP  ! missing file(s)

  ! re-use lchk for binning control : TRUE if no particular binning specified
  lchk = lbin(1) .OR. lbin(2) .OR. lbin(3) 

  npiglo = getdim (cf_vfil,cn_x)
  npjglo = getdim (cf_vfil,cn_y)
  npk    = getdim (cf_vfil,cn_z)
  npt    = getdim (cf_vfil,cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  !setting up the building command in global attribute
  CALL SetGlobalAtt(cglobal, 'A')  ! append command name to global attribute

  !  Detects newmaskglo file
  lbas = .NOT. chkfile (cn_fbasins )

  IF (lbas) THEN
     nbasins = 5
  ELSE
     nbasins = 1
  ENDIF

  ALLOCATE ( stypvar(nbasins), ipk(nbasins), id_varout(nbasins) )

  IF ( lchk )  THEN  ! use default bins definition according to pref 
     ! Define density parameters
     SELECT CASE ( INT(pref) )
     CASE ( 0 )
        nbins  = 52
        sigmin  = 23.
        sigstp = 0.1
     CASE ( 1000 )
        nbins  = 88
        sigmin  = 24.
        sigstp = 0.1
     CASE ( 2000)
        nbins  = 158
        sigmin  = 30.
        sigstp = 0.05
     CASE DEFAULT
        PRINT *,' This value of depth_ref (',pref,') is not implemented as standard'
        PRINT *,' You must use the -sigmin, -sigstp and -nbins options to precise'
        PRINT *,' the density bining you want to use.'
        STOP
     END SELECT
  ENDIF
  PRINT '(a,f6.1,a)',           '  For reference depth ', pref, ' m, '
  PRINT '(a,f5.2,a,f5.2,a,i3)', '  You are using -sigmin ', sigmin,' -sigstp ', sigstp,' -nbins ', nbins

  ALLOCATE ( sigma(nbins) )  

  ! define densities at middle of bins
  DO ji=1,nbins
     sigma(ji)  = sigmin +(ji-0.5)*sigstp
  ENDDO
  IF (lprint) PRINT *, ' min density:',sigma(1), ' max density:', sigma(nbins)

  !global   ; Atlantic  ; Indo-Pacif ; Indian  ; Pacif
  npglo= 1  ; npatl=2   ;  npinp=3   ; npind=4 ; nppac=5

  ! Common to all variables :
  stypvar%cunits            = 'Sverdrup'
  stypvar%rmissing_value    = 99999.
  stypvar%valid_min         = -1000.
  stypvar%valid_max         =  1000.
  stypvar%scale_factor      = 1.
  stypvar%add_offset        = 0.
  stypvar%savelog10         = 0.
  stypvar%conline_operation = 'N/A'
  stypvar%caxis             = 'TZY'

  ipk(:) = npk  

  ! Global basin
  stypvar(npglo)%cname       = cn_zomsfglo
  stypvar(npglo)%clong_name  = 'Meridional_Overt.Cell_Global'
  stypvar(npglo)%cshort_name = cn_zomsfglo

  IF (lbas) THEN
     stypvar(npatl)%cname       = cn_zomsfatl
     stypvar(npatl)%clong_name  = 'Meridional_Overt.Cell_Atlantic'
     stypvar(npatl)%cshort_name = cn_zomsfatl

     stypvar(npinp)%cname       = cn_zomsfinp
     stypvar(npinp)%clong_name  = 'Meridional_Overt.Cell_IndoPacif'
     stypvar(npinp)%cshort_name = cn_zomsfinp

     stypvar(npind)%cname       = cn_zomsfind
     stypvar(npind)%clong_name  = 'Meridional_Overt.Cell_Indian'
     stypvar(npind)%cshort_name = cn_zomsfind

     stypvar(nppac)%cname       = cn_zomsfpac
     stypvar(nppac)%clong_name  = 'Meridional_Overt.Cell_pacif'
     stypvar(nppac)%cshort_name = cn_zomsfpac
  ENDIF

  ! Allocate arrays
  ALLOCATE ( ibmask(nbasins,npiglo,npjglo) )
  ALLOCATE ( zv (npiglo,npjglo), zt(npiglo,npjglo), zs(npiglo,npjglo))
  ALLOCATE ( e3v(npiglo,npjglo) )
  ALLOCATE ( ibin(npiglo, npjglo) )
  ALLOCATE ( e1v(npiglo,npjglo), gphiv(npiglo,npjglo) )
  ALLOCATE ( dmoc(nbasins, npjglo, nbins) )
  ALLOCATE ( dmoc_tmp(nbins,npiglo) )
  ALLOCATE ( rdumlon(1,npjglo) , rdumlat(1,npjglo))
  ALLOCATE ( dens(npiglo,npjglo))
  ALLOCATE ( itmask(npiglo,npjglo), zttmp(npiglo,npjglo))
  ALLOCATE ( tim(npt), e31d(npk)  )

  IF ( leiv ) THEN
     ALLOCATE ( zveiv (npiglo,npjglo))
  END IF

  e1v(:,:)   = getvar(cn_fhgr,   cn_ve1v,  1, npiglo, npjglo) 
  gphiv(:,:) = getvar(cn_fhgr,   cn_gphiv, 1, npiglo, npjglo)

  IF ( lfull ) e31d(:)  = getvare3(cn_fzgr, cn_ve3t, npk)

  iloc         = MAXLOC(gphiv)
  rdumlat(1,:) = gphiv(iloc(1),:)
  rdumlon(:,:) = 0.               ! set the dummy longitude to 0

  ! create output fileset
! ncout = create      (cf_moc, cf_vfil, 1,       npjglo, nbins,  cdep='sigma')
  ncout = create      (cf_moc, 'none', 1,       npjglo, nbins,  cdep='sigma')
  ierr  = createvar   (ncout,  stypvar, nbasins, ipk ,id_varout, cdglobal=cglobal)
  ierr  = putheadervar(ncout,  cf_vfil, 1,       npjglo, nbins,  pnavlon=rdumlon, pnavlat=rdumlat, pdep=sigma)

  tim  = getvar1d(cf_vfil, cn_vtimec, npt     )
  ierr = putvar1d(ncout,   tim,       npt, 'T')

  ! reading the masks
  ibmask(npglo,:,:) = getvar(cn_fmsk, 'vmask', 1, npiglo, npjglo)

  IF ( lbas ) THEN
     ibmask(npatl,:,:) = getvar(cn_fbasins, 'tmaskatl', 1, npiglo, npjglo)
     ibmask(npind,:,:) = getvar(cn_fbasins, 'tmaskind', 1, npiglo, npjglo)
     ibmask(nppac,:,:) = getvar(cn_fbasins, 'tmaskpac', 1, npiglo, npjglo)
     ibmask(npinp,:,:) = ibmask(nppac,:,:) + ibmask(npind,:,:)
     ! ensure that there are no overlapping on the masks
     WHERE(ibmask(npinp,:,:) > 0 ) ibmask(npinp,:,:) = 1
     ! change global mask for GLOBAL periodic condition
     ibmask(1,1,     :) = 0.
     ibmask(1,npiglo,:) = 0.
  ENDIF

  DO jt=1, npt
     ! initialize moc to 0
     dmoc(:,:,:) = 0.d0

     DO jk=1,npk-1
        !               for testing purposes only loop from 2 to 400
        IF (lprint) PRINT *,' working at depth ',jk
        ! Get velocities v at jj
        zv(:,:) = getvar(cf_vfil, cn_vomecrty, jk, npiglo, npjglo)
        IF ( leiv ) THEN
           zveiv(:,:) = getvar(cf_vfil, cn_vomeeivv, jk, npiglo,npjglo)
           zv(:,:)    = zv(:,:) + zveiv(:,:)
        END IF
        zt(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo)
        zs(:,:) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo)

        ! get e3v at latitude jj
        IF ( lfull ) THEN
          e3v(:,:) = e31d(jk)
        ELSE
          e3v(:,:) = getvar(cn_fzgr, 'e3v_ps', jk, npiglo, npjglo, ldiom=.TRUE.)
        ENDIF
        !
        !  finds density 
        itmask =  1
        WHERE ( zs == 0 ) itmask = 0
        dens  = sigmai(zt, zs, pref, npiglo, npjglo)
        zttmp = dens* itmask ! convert to single precision 
        ! find bin numbers
        ibin(:,:) = INT( (zttmp-sigmin)/sigstp )
        ibin(:,:) = MAX( ibin(:,:), 1    )
        ibin(:,:) = MIN( ibin(:,:), nbins)

        DO jj=2,npjglo-1
           dmoc_tmp = 0
           !  converts transport in "k" to transport in "sigma"
           !  indirect adresssing - do it once and not for each basin!
           DO ji=2,npiglo-1
              dmoc_tmp(ibin(ji,jj),ji)=dmoc_tmp(ibin(ji,jj),ji) - e1v(ji,jj)*e3v(ji,jj)*zv(ji,jj)
           END DO
           ! integrates 'zonally' (along i-coordinate) 
           ! add to dmoc the contributions from level jk  at all densities jbin
           DO jbin =1,nbins  
              DO ji=2,npiglo-1
                 DO jbasin= 1, nbasins
                    ! For all basins 
                    dmoc(jbasin,jj,jbin)=dmoc(jbasin,jj,jbin ) + dmoc_tmp(jbin,ji) * ibmask(jbasin,ji,jj)
                 ENDDO
              END DO
           END DO
           !               end of loop on latitude for filling dmoc
        END DO
        !  end of loop on depths for calculating transports     
     END DO

     ! integrates across bins from highest to lowest density
     dmoc(:,:,nbins) = dmoc(:,:,nbins)/1.e6
     DO jk=nbins-1, 1, -1
        dmoc(:,:,jk) = dmoc(:,:,jk+1) + dmoc(:,:,jk)/1.e6
     END DO  ! loop to next bin

     ! netcdf output  
     DO jbasin = 1, nbasins
        DO jk = 1, nbins
           ierr = putvar (ncout, id_varout(jbasin), REAL(dmoc(jbasin,:,jk)), jk, 1, npjglo)
        END DO
     END DO

  ENDDO  ! time loop

  ierr = closeout(ncout)

END PROGRAM cdfmocsig
   
