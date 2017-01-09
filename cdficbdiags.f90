PROGRAM cdficbdiag
  !!======================================================================
  !!                     ***  PROGRAM  cdficediag  ***
  !!=====================================================================
  !!  ** Purpose : Compute the Ice volume, area and extend for each 
  !!               hemisphere
  !!
  !!  ** Method  : Use the icemod files for input and determine the
  !!               hemisphere with sign of the coriolis parameter.
  !!
  !! History : 2.1  : 01/2006  : J.M. Molines : Original code
  !!         : 2.1  : 07/2009  : R. Dussin    : Add Ncdf output
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !! Modified: 3.0  : 08/2011  : P.   Mathiot : Add LIM3 option
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id: cdficediags.f90 759 2014-07-21 22:01:28Z molines $
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: jk, jj, jt           ! dummy loop index
  INTEGER(KIND=4)                            :: ierr                 ! working integer
  INTEGER(KIND=4)                            :: narg, iargc          ! command line 
  INTEGER(KIND=4)                            :: npiglo, npjglo, npt  ! size of the domain
  INTEGER(KIND=4)                            :: nvpk                 ! vertical levels in working variable
  INTEGER(KIND=4)                            :: nperio = 4           ! boundary condition ( periodic, north fold)
  INTEGER(KIND=4)                            :: ikx=1, iky=1, ikz=0  ! dims of netcdf output file
  INTEGER(KIND=4)                            :: nboutput=4           ! number of values to write in cdf output
  INTEGER(KIND=4)                            :: ncout                ! for netcdf output
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout

  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: e1, e2               ! metrics
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: tmask, ff            ! npiglo x npjglo
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: ricbmass, ricbmelt ! thickness, leadfrac (concentration)
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: rdumlon, rdumlat     ! dummy lon lat for output
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: tim                  ! time counter

  REAL(KIND=8)                               :: dmasss, dmelts        ! volume, area extend South hemisphere
  REAL(KIND=8)                               :: dextends, dextends2  ! volume, area extend South hemisphere
  REAL(KIND=8)                               :: dmassn, dmeltn        ! volume, area extend North hemisphere
  REAL(KIND=8)                               :: dextendn, dextendn2  ! volume, area extend North hemisphere

  TYPE(variable), DIMENSION(:),  ALLOCATABLE :: stypvar              ! structure of output
  !
  CHARACTER(LEN=256)                         :: cf_ifil              ! input ice file
  CHARACTER(LEN=256)                         :: cf_out='icbdiags.nc' ! output file
  CHARACTER(LEN=256)                         :: cldum                ! dummy string
  !
  LOGICAL                                    :: lchk  = .false.      ! missing file flag
  LOGICAL                                    :: llim3 = .false.      ! LIM3 flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdficediag ICE-file [-lim3] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute the spatially integrated icb mass and melt flux.'
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       ICB-file : a single netcdf icb file' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cn_fhgr),' and ',TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : [NS]Mass  (Kg )'
     PRINT *,'                     [NS]Melt    (Kg/s )'
     PRINT *,'               N = northern hemisphere'
     PRINT *,'               S = southern hemisphere'
     PRINT *,'       standard output'
     STOP
  ENDIF

  CALL getarg (1, cf_ifil)

  lchk = lchk .OR. chkfile(cn_fhgr) 
  lchk = lchk .OR. chkfile(cn_fmsk) 
  lchk = lchk .OR. chkfile(cf_ifil)

  IF ( lchk ) STOP ! missing file

  npiglo = getdim (cf_ifil,cn_x)
  npjglo = getdim (cf_ifil,cn_y)
  npt    = getdim (cf_ifil,cn_t)

  ALLOCATE ( tmask(npiglo,npjglo) ,ff(npiglo,npjglo) )
  ALLOCATE ( ricbmass(npiglo,npjglo) )
  ALLOCATE ( ricbmelt(npiglo,npjglo) )
  ALLOCATE ( e1(npiglo,npjglo),e2(npiglo,npjglo) )
  ALLOCATE ( tim(npt) )

  ALLOCATE ( stypvar(nboutput), ipk(nboutput), id_varout(nboutput) )
  ALLOCATE ( rdumlon(1,1), rdumlat(1,1) )

  rdumlon(:,:) = 0.
  rdumlat(:,:) = 0.

  ipk(:) = 1

  ! define new variables for output 
  stypvar%scale_factor      = 1.
  stypvar%add_offset        = 0.
  stypvar%savelog10         = 0.
  stypvar%conline_operation = 'N/A'
  stypvar%caxis             = 'T'

  stypvar(1)%cname          = 'NMass'
  stypvar(1)%cunits         = 'Kg'
  stypvar(1)%clong_name     = 'Icb_Mass_in_Northern_Hemisphere'
  stypvar(1)%cshort_name    = 'NMass'

  stypvar(2)%cname          = 'NMelt'
  stypvar(2)%cunits         = 'Kg/s'
  stypvar(2)%clong_name     = 'Icb_melt_in_Northern_Hemisphere'
  stypvar(2)%cshort_name    = 'NMelt'

  stypvar(3)%cname          = 'SVMass'
  stypvar(3)%cunits         = 'Kg'
  stypvar(3)%clong_name     = 'Icb_Mass_in_Southern_Hemisphere'
  stypvar(3)%cshort_name    = 'SMass'

  stypvar(4)%cname          = 'SMelt'
  stypvar(4)%cunits         = 'Kg/s'
  stypvar(4)%clong_name     = 'Icb_Melt_in_Southern_Hemisphere'
  stypvar(4)%cshort_name    = 'SMelt'


  e1(:,:) = getvar(cn_fhgr, cn_ve1t,  1, npiglo, npjglo)
  e2(:,:) = getvar(cn_fhgr, cn_ve2t,  1, npiglo, npjglo)
  ff(:,:) = getvar(cn_fhgr, cn_gphit, 1, npiglo, npjglo) ! only the sign of ff is important

  ! modify the mask for periodic and north fold condition (T pivot, F Pivot ...)
  ! in fact should be nice to use jperio as in the code ...
  tmask(:,:)=getvar(cn_fmsk,'tmask',1,npiglo,npjglo)
  SELECT CASE (nperio)
  CASE (0) ! closed boundaries
     ! nothing to do
  CASE (4) ! ORCA025 type boundary
     tmask(1:2,:)=0.
     tmask(:,npjglo)=0.
     tmask(npiglo/2+1:npiglo,npjglo-1)= 0.
  CASE (6)
     tmask(1:2,:)=0.
     tmask(:,npjglo)=0.
  CASE DEFAULT
     PRINT *,' Nperio=', nperio,' not yet coded'
     STOP
  END SELECT

  ricbmass(:,:)=0.
  ricbmelt(:,:)=0.


  IF (chkvar(cf_ifil, cn_iicbmass)) STOP

  IF (chkvar(cf_ifil, cn_iicbmelt)) STOP
  !
  DO jt = 1, npt
     IF (TRIM(cn_iicbmass) .NE. 'missing') ricbmass(:,:) = getvar(cf_ifil, cn_iicbmass, 1, npiglo, npjglo, ktime=jt)
     ricbmelt(:,:) = getvar(cf_ifil, cn_iicbmelt, 1, npiglo, npjglo, ktime=jt)

     ! North : ff > 0 
     dmassn     = SUM( ricbmass (:,:)* e1(:,:) * e2(:,:) * tmask (:,:), (ff > 0 ) )
     dmeltn    = SUM( ricbmelt (:,:)* e1(:,:) * e2(:,:) * tmask (:,:), (ff > 0 ) )

     ! South : ff < 0
     dmasss     = SUM( ricbmass (:,:)* e1(:,:) * e2(:,:) * tmask (:,:), (ff < 0) )
     dmelts    = SUM( ricbmelt (:,:)* e1(:,:) * e2(:,:) * tmask (:,:), (ff < 0 )) 

     PRINT *,' TIME = ', jt,' ( ',tim(jt),' )'
     PRINT *,' Northern Hemisphere ' 
     PRINT *,'          NMass (Kg)  ', dmassn
     PRINT *,'          NMelt (Kg/s)    ', dmeltn
     PRINT *
     PRINT *,' Southern Hemisphere ' 
     PRINT *,'          Mass (Kg)  ', dmasss
     PRINT *,'          Melt (Kg/s)    ', dmelts

     IF ( jt == 1 ) THEN
        ! create output fileset
        ncout = create      (cf_out, 'none',  ikx,      iky, ikz,     cdep='depthw'                   )
        ierr  = createvar   (ncout,  stypvar, nboutput, ipk, id_varout                                )
        ierr  = putheadervar(ncout,  cf_ifil, ikx,      iky, ikz,     pnavlon=rdumlon, pnavlat=rdumlat)

        tim   = getvar1d(cf_ifil, cn_vtimec, npt     )
        ierr  = putvar1d(ncout,   tim,       npt, 'T')
     ENDIF

     ! netcdf output 
     ierr = putvar0d(ncout,id_varout(1), REAL(dmassn), ktime=jt)
     ierr = putvar0d(ncout,id_varout(2), REAL(dmeltn), ktime=jt)
     ierr = putvar0d(ncout,id_varout(3), REAL(dmasss), ktime=jt)
     ierr = putvar0d(ncout,id_varout(4), REAL(dmelts), ktime=jt)

  END DO ! time loop
  ierr = closeout(ncout)

END PROGRAM cdficbdiag
