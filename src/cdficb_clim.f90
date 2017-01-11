PROGRAM cdficb_clim
  !!======================================================================
  !!                     ***  PROGRAM  cdficb_clim  ***
  !!=====================================================================
  !!  ** Purpose : Compute the iceberg mass and melt
  !!
  !!  ** Method  : Use the icb files for input and determine the
  !!
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: jk, jj, jt, ji       ! dummy loop index
  INTEGER(KIND=4)                            :: ierr                 ! working integer
  INTEGER(KIND=4)                            :: narg, iargc          ! command line 
  INTEGER(KIND=4)                            :: npiglo, npjglo, npt, numFiles  ! size of the domain
  INTEGER(KIND=4)                            :: nvpk                 ! vertical levels in working variable
  INTEGER(KIND=4)                            :: nperio = 4           ! boundary condition ( periodic, north fold)
  INTEGER(KIND=4)                            :: ikx, iky, ikz=0      ! dims of netcdf output file
  INTEGER(KIND=4)                            :: nboutput=2           ! number of values to write in cdf output
  INTEGER(KIND=4)                            :: ncout                ! for netcdf output
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout, itimeVar

  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: e1, e2               ! metrics
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: tmask, ff            ! npiglo x npjglo
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: ricbmass, ricbmelt   ! icbmass icbmelt
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: rdumlon, rdumlat     ! dummy lon lat for output
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: tim                  ! time counter

  TYPE(variable), DIMENSION(:),  ALLOCATABLE :: stypvar              ! structure of output
  !
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cf_icb            ! input icb file
  CHARACTER(LEN=256)                         :: cf_out='icbdiags.nc' ! output file
  CHARACTER(LEN=256)                         :: cldum                ! dummy string
  !
  LOGICAL                                    :: lchk  = .false.      ! missing file flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdficb_clim 12-ICB-monthly-means-files'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute the 2D field of icb mass and icb melt.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       ICE-file : netcdf icb file' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cn_fhgr),' and ',TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : Mass  (Kg/m2 )'
     PRINT *,'                     Melt  (Kg/m2/s )'
     STOP
  ENDIF

  CALL getarg(1,cldum)
  READ(cldum,*) numFiles
  
  IF (numFiles < 12) STOP

  ALLOCATE(cf_icb(numFiles))

  DO ji= 1, numFiles
        CALL getarg(ji+1,cf_icb(ji))
        lchk = lchk .OR. chkfile(cf_icb(ji))
  END DO


  lchk = lchk .OR. chkfile(cn_fhgr) 
  lchk = lchk .OR. chkfile(cn_fmsk) 

  IF ( lchk ) STOP ! missing file

  npiglo = getdim (cf_icb(1),cn_x)
  npjglo = getdim (cf_icb(1),cn_y)
  npt    = 12
  ikx = npiglo
  iky = npjglo

  ALLOCATE ( tmask(npiglo,npjglo) ,ff(npiglo,npjglo) )
  ALLOCATE ( ricbmass(npiglo,npjglo) )
  ALLOCATE ( ricbmelt(npiglo,npjglo) )
  ALLOCATE ( e1(npiglo,npjglo),e2(npiglo,npjglo) )
  ALLOCATE ( tim(npt),itimeVar(npt) )
 
  itimeVar = (/(ji,ji=1,12)/)  

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

  stypvar(1)%cname          = 'Mass'
  stypvar(1)%cunits         = 'Kg/m2'
  stypvar(1)%clong_name     = 'Icb mass per unit of area'
  stypvar(1)%cshort_name    = 'Mass'

  stypvar(2)%cname          = 'Melt'
  stypvar(2)%cunits         = 'Kg/m2/s'
  stypvar(2)%clong_name     = 'Icb melt flux'
  stypvar(2)%cshort_name    = 'Melt'


  e1(:,:) = getvar(cn_fhgr, cn_ve1t,  1, npiglo, npjglo)
  e2(:,:) = getvar(cn_fhgr, cn_ve2t,  1, npiglo, npjglo)
  ff(:,:) = getvar(cn_fhgr, cn_gphit, 1, npiglo, npjglo) ! only the sign of ff is important

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

  itimeVar = (/(ji,ji=1,12)/)
  ! Check variable
  IF (chkvar(cf_icb(1), cn_iicbmass)) THEN
     cn_iicbmass='missing'
     PRINT *,'' 
     PRINT *,' WARNING, ICEBERG MASS IS SET TO 0. '
     PRINT *,' '
  END IF

  !
  DO jt = 1, npt
     IF (TRIM(cn_iicbmass) /= 'missing') ricbmass(:,:) = getvar(cf_icb(jt), cn_iicbmass, 1, npiglo, npjglo, ktime=1)
     ricbmelt(:,:) = getvar(cf_icb(jt), cn_iicbmelt, 1, npiglo, npjglo, ktime=1)

    IF ( jt == 1 ) THEN
         ! create output fileset
        ncout = create      (cf_out, 'none',  ikx,      iky, ikz,     cdep='depthw'                   )
        ierr  = createvar   (ncout,  stypvar, nboutput, ipk, id_varout                                )
        ierr  = putheadervar(ncout,  cf_icb(1), ikx,      iky, ikz)

        tim   = getvar1d(cf_icb(1), cn_vtimec, npt     )
        tim = (/(ji,ji=1,12)/)
        ierr  = putvar1d(ncout,   tim,       npt, 'T')
     ENDIF

     ! netcdf output
     !ierr = putvar0d(ncout,cn_vtimec,jt,ktime=jt) 
     ierr = putvar(ncout,id_varout(1),REAL(ricbmass(:,:)),1,npiglo,npjglo, ktime=jt)
     ierr = putvar(ncout,id_varout(2),REAL(ricbmelt(:,:)),1,npiglo,npjglo, ktime=jt)

  END DO ! time loop
  ierr = closeout(ncout)

END PROGRAM cdficb_clim
