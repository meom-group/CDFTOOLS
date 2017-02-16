MODULE icediag

CONTAINS

    SUBROUTINE icediags(e1, e2, tmask, ff, ricethick, riceldfra, dvoln, darean, dextendn, &
                        dextendn2, dvols, dareas, dextends, dextends2)
        REAL(KIND=4), DIMENSION(:,:), INTENT(IN) :: e1, e2               ! metrics
        REAL(KIND=4), DIMENSION(:,:), INTENT(IN) :: tmask, ff            ! npiglo x npjglo
        REAL(KIND=4), DIMENSION(:,:), INTENT(IN) :: ricethick, riceldfra ! thickness, leadfrac (concentration)
        REAL(KIND=8), INTENT(OUT)        :: dvols, dareas        ! volume, area extend South hemisphere
        REAL(KIND=8), INTENT(OUT)        :: dextends, dextends2  ! volume, area extend South hemisphere
        REAL(KIND=8), INTENT(OUT)        :: dvoln, darean        ! volume, area extend North hemisphere
        REAL(KIND=8), INTENT(OUT)        :: dextendn, dextendn2


        ! North : ff > 0
        dvoln     = SUM( ricethick (:,:)* e1(:,:) * e2(:,:) * riceldfra (:,:) * tmask (:,:), (ff > 0 ) )
        darean    = SUM(                  e1(:,:) * e2(:,:) * riceldfra (:,:) * tmask (:,:), (ff > 0 ) )
        dextendn  = SUM(                  e1(:,:) * e2(:,:) * riceldfra (:,:) * tmask (:,:), (riceldfra > 0.15 .AND. ff > 0 ) )
        ! JMM added 22/01/2007 : to compute same extent than the NSIDC
        dextendn2 = SUM(                  e1(:,:) * e2(:,:)                   * tmask (:,:), (riceldfra > 0.15 .AND. ff > 0 ) )

        ! South : ff < 0
        dvols     = SUM( ricethick (:,:)* e1(:,:) * e2(:,:) * riceldfra (:,:) * tmask (:,:), (ff < 0 ) )
        dareas    = SUM(                  e1(:,:) * e2(:,:) * riceldfra (:,:) * tmask (:,:), (ff < 0 ) )
        dextends  = SUM(                  e1(:,:) * e2(:,:) * riceldfra (:,:) * tmask (:,:), (riceldfra > 0.15 .AND. ff < 0  ) )
        dextends2 = SUM(                  e1(:,:) * e2(:,:)                   * tmask (:,:), (riceldfra > 0.15 .AND. ff < 0  ) )

        dvoln = dvoln / 1d9
        darean = darean / 1d9
        dextendn = dextendn / 1d9
        dextendn2 = dextendn2 / 1d9

        dvols = dvols / 1d9
        dareas = dareas / 1d9
        dextends = dextends / 1d9
        dextends2 = dextends2 / 1d9
    END SUBROUTINE icediags

END MODULE icediag


PROGRAM cdficediag
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
  USE icediag
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: jt           ! dummy loop index
  INTEGER(KIND=4)                            :: ierr                 ! working integer
  INTEGER(KIND=4)                            :: narg, iargc          ! command line
  INTEGER(KIND=4)                            :: ijarg, ireq         ! command line
  INTEGER(KIND=4)                            :: npiglo, npjglo, npt  ! size of the domain
  INTEGER(KIND=4)                            :: nperio = 4           ! boundary condition ( periodic, north fold)
  INTEGER(KIND=4)                            :: ikx=1, iky=1, ikz=0  ! dims of netcdf output file
  INTEGER(KIND=4)                            :: nboutput=8           ! number of values to write in cdf output
  INTEGER(KIND=4)                            :: ncout                ! for netcdf output
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout

  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: e1, e2               ! metrics
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: tmask, ff            ! npiglo x npjglo
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: ricethick, riceldfra ! thickness, leadfrac (concentration)
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: rdumlon, rdumlat     ! dummy lon lat for output
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: tim                  ! time counter

  REAL(KIND=8)                               :: dvols, dareas        ! volume, area extend South hemisphere
  REAL(KIND=8)                               :: dextends, dextends2  ! volume, area extend South hemisphere
  REAL(KIND=8)                               :: dvoln, darean        ! volume, area extend North hemisphere
  REAL(KIND=8)                               :: dextendn, dextendn2  ! volume, area extend North hemisphere

  TYPE(variable), DIMENSION(:),  ALLOCATABLE :: stypvar              ! structure of output
  !
  CHARACTER(LEN=2048)                         :: cf_ifil              ! input ice file
  CHARACTER(LEN=2048)                         :: cf_out='icediags.nc' ! output file
  CHARACTER(LEN=2048)                         :: cldum                ! dummy string
  CHARACTER(LEN=2048)                         :: cn_mask='tmask'      ! mask variable name
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
     PRINT *,'        Compute the ice volume, area and extent for each hemisphere.'
     PRINT *,'        The extent is computed in a similar way to NSIDC for easy '
     PRINT *,'        comparison : the extent is the surface of the grid cells covered'
     PRINT *,'        by ice when the ice concentration is above 0.15'
     PRINT *,'      '
     PRINT *,'        For compatibility with previous version, another estimate of '
     PRINT *,'        the extend is computed using grid cell surfaces weighted by the'
     PRINT *,'        ice concentration, but it will be deprecated soon.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       ICE-file : netcdf icemod file (LIM2 by default)' 
     PRINT *,'      '
     PRINT *,'     OPTION :'
     PRINT *,'       [-lim3 ] : LIM3 variable name convention is used'
     PRINT *,'       [-o OUT-file ] : specify output file instead of ',TRIM(cf_out)
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cn_fhgr),' and ',TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : [NS]Volume  (10^9 m3 )'
     PRINT *,'                     [NS]Area    (10^9 m2 )'
     PRINT *,'                     [NS]Extent  (10^9 m2 ) -- obsolete --'
     PRINT *,'                     [NS]Exnsidc (10^9 m2 )'
     PRINT *,'               N = northern hemisphere'
     PRINT *,'               S = southern hemisphere'
     PRINT *,'       standard output'
     STOP
  ENDIF

  CALL getarg (1, cf_ifil)
  
  ijarg = 1 ; ireq = 0
  
  DO WHILE ( ijarg <= narg )
    CALL getarg (ijarg, cldum   ) ; ijarg = ijarg + 1 
    SELECT CASE ( cldum )
    CASE ( '-lim3'    ) ; llim3 = .true.
    CASE ( '-lim2' ) ; llim3 = .false.
    CASE ( '-mask' ) ;
        CALL getarg (ijarg, cn_mask) ; ijarg=ijarg+1
    CASE ( '-maskfile' ) ;
        CALL getarg (ijarg, cn_fmsk) ; ijarg=ijarg+1
    CASE ('-o') 
        CALL getarg (ijarg, cf_out) ; ijarg=ijarg+1
    CASE DEFAULT 
      ireq=ireq+1
      SELECT CASE ( ireq )
      CASE ( 1 )    ; cf_ifil=cldum
      CASE DEFAULT  ; PRINT *,' Too many arguments'
      END SELECT
    END SELECT
  END DO

  lchk = lchk .OR. chkfile(cn_fhgr)
  lchk = lchk .OR. chkfile(cn_fmsk)
  lchk = lchk .OR. chkfile(cf_ifil)

  IF ( lchk ) STOP ! missing file

  npiglo = getdim (cf_ifil,cn_x)
  npjglo = getdim (cf_ifil,cn_y)
  npt    = getdim (cf_ifil,cn_t)

  ALLOCATE ( tmask(npiglo,npjglo) ,ff(npiglo,npjglo) )
  ALLOCATE ( ricethick(npiglo,npjglo) )
  ALLOCATE ( riceldfra(npiglo,npjglo) )
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

  stypvar(1)%cname          = 'NVolume'
  stypvar(1)%cunits         = '10^9 m3'
  stypvar(1)%clong_name     = 'Ice_volume_in_Northern_Hemisphere'
  stypvar(1)%cshort_name    = 'NVolume'

  stypvar(2)%cname          = 'NArea'
  stypvar(2)%cunits         = '10^9 m2'
  stypvar(2)%clong_name     = 'Ice_area_in_Northern_Hemisphere'
  stypvar(2)%cshort_name    = 'NArea'

  stypvar(3)%cname          = 'NExtent'
  stypvar(3)%cunits         = '10^9 m2'
  stypvar(3)%clong_name     = 'Ice_extent_in_Northern_Hemisphere'
  stypvar(3)%cshort_name    = 'NExtent'

  stypvar(4)%cname          = 'NExnsidc'
  stypvar(4)%cunits         = '10^9 m2'
  stypvar(4)%clong_name     = 'Ice_extent_similar_to_NSIDC_in_Northern_Hemisphere'
  stypvar(4)%cshort_name    = 'NExnsidc'

  stypvar(5)%cname          = 'SVolume'
  stypvar(5)%cunits         = '10^9 m3'
  stypvar(5)%clong_name     = 'Ice_volume_in_Southern_Hemisphere'
  stypvar(5)%cshort_name    = 'SVolume'

  stypvar(6)%cname          = 'SArea'
  stypvar(6)%cunits         = '10^9 m2'
  stypvar(6)%clong_name     = 'Ice_area_in_Southern_Hemisphere'
  stypvar(6)%cshort_name    = 'SArea'

  stypvar(7)%cname          = 'SExtent'
  stypvar(7)%cunits         = '10^9 m2'
  stypvar(7)%clong_name     = 'Ice_extent_in_Southern_Hemisphere'
  stypvar(7)%cshort_name    = ''

  stypvar(8)%cname          = 'SExnsidc'
  stypvar(8)%cunits         = '10^9 m2'
  stypvar(8)%clong_name     = 'Ice_extent_similar_to_NSIDC_in_Southern_Hemisphere'
  stypvar(8)%cshort_name    = 'SExnsidc'

  e1(:,:) = getvar(cn_fhgr, cn_ve1t,  1, npiglo, npjglo)
  e2(:,:) = getvar(cn_fhgr, cn_ve2t,  1, npiglo, npjglo)
  ff(:,:) = getvar(cn_fhgr, cn_gphit, 1, npiglo, npjglo) ! only the sign of ff is important

  ! modify the mask for periodic and north fold condition (T pivot, F Pivot ...)
  ! in fact should be nice to use jperio as in the code ...
  tmask(:,:)=getvar(cn_fmsk,cn_mask,1,npiglo,npjglo)
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

  ricethick(:,:)=0.
  riceldfra(:,:)=0.

  IF (llim3) THEN
     cn_iicethic = cn_iicethic3
     cn_ileadfra = cn_ileadfra3
  END IF

  ! Check variable
  IF (chkvar(cf_ifil, cn_iicethic)) THEN
     cn_iicethic='missing'
     PRINT *,'' 
     PRINT *,' WARNING, ICE THICKNESS IS SET TO 0. '
     PRINT *,' '
  END IF

  IF (chkvar(cf_ifil, cn_ileadfra)) STOP
  !
  DO jt = 1, npt
     IF (TRIM(cn_iicethic) .NE. 'missing') ricethick(:,:) = getvar(cf_ifil, cn_iicethic, 1, npiglo, npjglo, ktime=jt)
     riceldfra(:,:) = getvar(cf_ifil, cn_ileadfra, 1, npiglo, npjglo, ktime=jt)

     CALL icediags(e1, e2, tmask, ff, ricethick, riceldfra, dvoln, darean, dextendn, dextendn2, dvols, dareas, dextends, dextends2)

     PRINT *,' TIME = ', jt,' ( ',tim(jt),' )'
     PRINT *,' Northern Hemisphere ' 
     PRINT *,'          NVolume (10^9 m3)  ', dvoln
     PRINT *,'          NArea (10^9 m2)    ', darean
     PRINT *,'          NExtend (10^9 m2)  ', dextendn
     PRINT *,'          NExnsidc (10^9 m2) ', dextendn2
     PRINT *
     PRINT *,' Southern Hemisphere ' 
     PRINT *,'          SVolume (10^9 m3)  ', dvols
     PRINT *,'          SArea (10^9 m2)    ', dareas
     PRINT *,'          SExtend (10^9 m2)  ', dextends
     PRINT *,'          SExnsidc (10^9 m2) ', dextends2

     IF ( jt == 1 ) THEN
        ! create output fileset
        ncout = create      (cf_out, 'none',  ikx,      iky, ikz,     cdep='depthw'                   )
        ierr  = createvar   (ncout,  stypvar, nboutput, ipk, id_varout                                )
        ierr  = putheadervar(ncout,  cf_ifil, ikx,      iky, ikz,     pnavlon=rdumlon, pnavlat=rdumlat)

        tim   = getvar1d(cf_ifil, cn_vtimec, npt     )
        ierr  = putvar1d(ncout,   tim,       npt, 'T')
     ENDIF

     ! netcdf output 
     ierr = putvar0d(ncout,id_varout(1), REAL(dvoln     ), ktime=jt)
     ierr = putvar0d(ncout,id_varout(2), REAL(darean    ), ktime=jt)
     ierr = putvar0d(ncout,id_varout(3), REAL(dextendn  ), ktime=jt)
     ierr = putvar0d(ncout,id_varout(4), REAL(dextendn2 ), ktime=jt)
     ierr = putvar0d(ncout,id_varout(5), REAL(dvols     ), ktime=jt)
     ierr = putvar0d(ncout,id_varout(6), REAL(dareas    ), ktime=jt)
     ierr = putvar0d(ncout,id_varout(7), REAL(dextends  ), ktime=jt)
     ierr = putvar0d(ncout,id_varout(8), REAL(dextends2 ), ktime=jt)
  END DO ! time loop
  ierr = closeout(ncout)

END PROGRAM cdficediag


