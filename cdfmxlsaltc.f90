PROGRAM cdfmxlsaltc
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfmxlsaltc  ***
  !!
  !!  **  Purpose  :  Compute the heat content in the mixed layer
  !!                  PARTIAL STEPS
  !!  
  !!  **  Method   :  compute the sum ( S  * e1 *e2 * e3 *mask )
  !!                  for the mixed layer stored into gridT file
  !!
  !!
  !! history ;
  !!  Original :  J.M. Molines ( 2006) April 
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  e3,  zs   !:  metrics, salinity
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmxl              !:  mxl depth
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask             !:  npiglo x npjglo
  REAL(KIND=4),DIMENSION(:), ALLOCATABLE       ::  gdepw             !:  

  REAL(KIND=8), PARAMETER :: rprho0=1020., rpcp=4000.
  REAL(KIND=8)      :: zvol, zsum, zvol2d, zsum2d, zsurf
  REAL(KIND=8), DIMENSION (:,:),   ALLOCATABLE ::  zmxlsaltc         !:  mxl salt content

  CHARACTER(LEN=80) :: cfilet 
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc',cmask='mask.nc'

  ! Output stuff
  INTEGER                         :: ncout, ierr
  INTEGER,           DIMENSION(1) :: ipk, id_varout  !: only one output variable
  REAL(KIND=4),      DIMENSION(1) :: tim,dep       !: time output
  CHARACTER(LEN=80)               :: cfileout='mxlsaltc.nc'

  TYPE(variable), DIMENSION(1)    :: typvar         !: stucture for attributes
  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmxlsaltc  gridTfile   '
     PRINT *,' Computes the heat content in the mixed layer (Joules)'
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc '
     PRINT *,'  must be in the current directory'
     PRINT *,' Output ncdf file mxlsaltc.nc, variable 2D somxlsaltc'
     STOP
  ENDIF

  CALL getarg (1, cfilet)

  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')

  dep(1) = 0.
  ipk(:) = 1
  typvar(1)%name= 'somxlsaltc'
  typvar(1)%units='kg/m2'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= 0.
  typvar(1)%valid_max= 1.e9
  typvar(1)%long_name='Mixed_Layer_Salt_Content'
  typvar(1)%short_name='somxlsaltc'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo) , zmxlsaltc(npiglo, npjglo) )
  ALLOCATE ( zs(npiglo,npjglo) ,zmxl(npiglo,npjglo)  )
  ALLOCATE ( e3(npiglo,npjglo) )
  ALLOCATE ( gdepw(npk) )

  ! Initialize output file
  ncout = create(cfileout, cfilet, npiglo,npjglo,1)
  ierr=createvar(ncout ,typvar,1, ipk,id_varout )
  ierr=putheadervar(ncout, cfilet,npiglo, npjglo,1,pdep=dep)
  tim=getvar1d(cfilet,'time_counter',1)
  ierr=putvar1d(ncout,tim,1,'T')

  ! Read vertical depth at w point
  gdepw(:) = getvare3(coordzgr,'gdepw',npk)

  ! Read Mixed Layer Depth in the gridT file
  ! Note that it is usually a mean value (5-day mean for instance). Therefore,
  ! in general it is not a particular value of gdepw.
  zmxl(:,:)=getvar(cfilet,'somxl010',1,npiglo,npjglo)

  zvol=0.d0
  zsum=0.d0
  zmxlsaltc(:,:)=0.d0
  DO jk = 1,npk
     ! Get temperatures at jk
     zs(:,:)= getvar(cfilet, 'vosaline',  jk ,npiglo,npjglo)
     zmask(:,:)=getvar(cmask,'tmask',jk,npiglo,npjglo)

     ! get e3 at level jk ( ps...)
     e3(:,:) = getvar(coordzgr, 'e3t_ps', jk,npiglo,npjglo, ldiom=.true.)

     !  e3 is used as a flag for the mixed layer; It is 0 outside the mixed layer
     e3(:,:)=MAX ( 0., MIN(e3,zmxl-gdepw(jk) ) )
     WHERE ( e3 == 0 ) zmask = 0.

     zvol2d=SUM( e3 * zmask)
     zmxlsaltc(:,:)=zmxlsaltc(:,:)+ zs*e3*zmask

     IF (zvol2d /= 0 )THEN
        !   go on !
     ELSE
        !   no more layer below !    
        EXIT   ! get out of the jk loop
     ENDIF

  END DO

  ! Output to netcdf file : kg/m2
  zmxlsaltc=zmxlsaltc*rprho0
  ierr = putvar(ncout, id_varout(1) ,REAL(zmxlsaltc), 1,npiglo, npjglo)
  ierr=closeout(ncout)

END PROGRAM cdfmxlsaltc
