PROGRAM cdfvertmean
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfvertmean  ***
  !!
  !!  **  Purpose  :  Compute the vertical average of a scalar quantity
  !!                  between z layers
  !!                  PARTIAL STEPS
  !!  
  !!  **  Method   :  compute the sum ( V  * e1 *e2 * e3 *mask )
  !!                  for the mixed layer stored into gridT file
  !!
  !!
  !! history ;
  !!  Original :  J.M. Molines ( 2008) January
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk, jvar
  INTEGER   :: k1,k2
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: nvars, ivar

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  e3,  zs   !:  metrics, salinity
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  hdep              !:  mxl depth
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask             !:  npiglo x npjglo
  REAL(KIND=8), DIMENSION (:,:),   ALLOCATABLE ::  zvol2d             !:  npiglo x npjglo
  REAL(KIND=4),DIMENSION(:), ALLOCATABLE       ::  gdep             !:  

  REAL(KIND=8)      :: zvol,  dep_up, dep_down
  REAL(KIND=8), DIMENSION (:,:),   ALLOCATABLE ::  zvertmean         !:  mxl salt content

  CHARACTER(LEN=80) :: cfilet 
  CHARACTER(LEN=80) :: coordzgr='mesh_zgr.nc',cmask='mask.nc'
  CHARACTER(LEN=80)               :: ctype='T'
  CHARACTER(LEN=80)               :: cdum
  CHARACTER(LEN=80)               :: cvarnam, cdep, ce3, cvmask
  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: cvarname    !: name of input variables
  TYPE(variable), DIMENSION(:),ALLOCATABLE    :: typvarin     !: stucture for attributes

  ! Output stuff
  INTEGER                         :: ncout, ierr
  INTEGER,           DIMENSION(1) :: ipk, id_varout  !: only one output variable
  REAL(KIND=4),      DIMENSION(1) :: tim,dep       !: time output
  CHARACTER(LEN=80)               :: cfileout='vertmean.nc'

  TYPE(variable), DIMENSION(1)    :: typvar         !: stucture for attributes
  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfvertmean  datafile varname T|U|V|W z1 z2  '
     PRINT *,' Computes the vertical mean value of variable between z1 and z2'
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_zgr.nc ,mask.nc '
     PRINT *,'  must be in the current directory'
     PRINT *,' Output ncdf file vertmean.nc, variable 2D sovertmean'
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  CALL getarg (2, cvarnam)
  CALL getarg (3, ctype)
  CALL getarg (4, cdum) ; READ(cdum,*) dep_up
  CALL getarg (5, cdum) ; READ(cdum,*) dep_down

  IF (dep_down < dep_up ) THEN
    PRINT *,'Give depth limits in increasing order !'
    STOP
  ENDIF

  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')

  nvars = getnvar(cfilet)
  ALLOCATE( cvarname(nvars), typvarin(nvars) )
  cvarname(:)=getvarname(cfilet,nvars,typvarin)
  ivar=1
  DO jvar=1,nvars
    IF ( TRIM(cvarname(jvar)) == TRIM(cvarnam) ) THEN
      EXIT
    ENDIF
   ivar=ivar+1
  ENDDO
  IF ( ivar == nvars+1 ) THEN
     PRINT *,' Variable ',TRIM(cvarnam),' not found in ', TRIM(cfilet)
     STOP
  ENDIF
  
  
  dep(1) = 0.
  ipk(:) = 1
  typvar(1)%name= 'sovertmean'
  typvar(1)%units=typvarin(ivar)%units
  typvar(1)%missing_value=typvarin(ivar)%missing_value
  typvar(1)%valid_min= typvarin(ivar)%valid_min
  typvar(1)%valid_max= typvarin(ivar)%valid_max
  typvar(1)%long_name='vertical average of '//TRIM(typvarin(ivar)%long_name)
  typvar(1)%short_name='sovertmean'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo) , zvertmean(npiglo, npjglo) )
  ALLOCATE ( zs(npiglo,npjglo) ,hdep(npiglo,npjglo)  )
  ALLOCATE ( e3(npiglo,npjglo) ,zvol2d(npiglo,npjglo) )
  ALLOCATE ( gdep(npk) )

  ! Initialize output file
  ncout = create(cfileout, cfilet, npiglo,npjglo,1)
  ierr=createvar(ncout ,typvar,1, ipk,id_varout )
  ierr=putheadervar(ncout, cfilet,npiglo, npjglo,1,pdep=dep)
  tim=getvar1d(cfilet,'time_counter',1)
  ierr=putvar1d(ncout,tim,1,'T')

  ! Read vertical depth at w point
  SELECT CASE ( ctype)
  CASE( 'T','U','V','t','u','v');  cdep='gdepw' ; ce3='e3t_ps' 
  CASE( 'W' ,'w')               ;  cdep='gdept' ; ce3='e3w_ps' 
  CASE DEFAULT ; PRINT *,'Point type ', TRIM(ctype),' not known! ' ; STOP
  END SELECT
     gdep(:) = getvare3(coordzgr,cdep,npk)

  ! set mask variable name
  SELECT CASE (ctype )
  CASE ('T','t','W','w') ; cvmask='tmask'
  CASE ('U','u')         ; cvmask='umask'
  CASE ('V','v')         ; cvmask='vmask'
  END SELECT

  ! Look for k1 and k2 as nearest level of dep_up and dep_down
  k1=1; k2=npk
  DO jk=1,npk
   IF ( gdep(jk) <= dep_up ) k1=jk
   IF ( gdep(jk) <= dep_down ) k2=jk
  ENDDO
    
  PRINT *, dep_up, dep_down, k1, k2 , gdep(k1), gdep(k2)
     

  zvol=0.d0
  zvertmean(:,:)=0.d0

  DO jk = k1, k2
     ! Get temperatures at jk
     zs(:,:)= getvar(cfilet, cvarnam,  jk ,npiglo,npjglo)
     zmask(:,:)=getvar(cmask,cvmask,jk,npiglo,npjglo)

     ! get e3 at level jk ( ps...)
     e3(:,:) = getvar(coordzgr, ce3, jk,npiglo,npjglo, ldiom=.true.)
     IF (jk == k1 ) THEN
       hdep(:,:)=gdep(jk)+e3(:,:)
       e3(:,:)=MIN(e3,hdep-REAL(dep_up))
     ENDIF
     IF ( jk == k2 ) THEN
      e3(:,:)=MIN(e3,REAL(dep_down)-gdep(jk))
     ENDIF

     zvol=SUM( e3 * zmask)
     zvol2d=e3(:,:)*zmask+ zvol2d(:,:)
     zvertmean(:,:)=zvertmean(:,:)+ zs*e3*zmask

     IF (zvol /= 0 )THEN
        !   go on !
     ELSE
        !   no more layer below !    
        EXIT   ! get out of the jk loop
     ENDIF

  END DO

  ! Output to netcdf file : kg/m2
  WHERE ( zvol2d /= 0 )  zvertmean=zvertmean/zvol2d
  ierr = putvar(ncout, id_varout(1) ,REAL(zvertmean), 1,npiglo, npjglo)
  ierr=closeout(ncout)

END PROGRAM cdfvertmean
