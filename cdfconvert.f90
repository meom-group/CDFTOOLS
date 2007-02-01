PROGRAM cdfconvert
  !!-------------------------------------------------------------------
  !!              PROGRAM CDFCONVERT
  !!              ******************
  !!
  !!  **  Purpose: Convert a set of dimgfile (Clipper like)
  !!               to a set of CDF files (Drakkar like )
  !!  
  !!  **  Method: Read tag then open the respective T S 2D U V files to create 
  !!              gridT, gridU and gridV files.
  !!              Requires  mesh_hgr.nc and mesh_zgr.nc  files
  !!
  !! history:
  !!    Original:  J.M. Molines (Jan. 2007 )
  !!-------------------------------------------------------------------
  !!
  !! * Modules used
  USE cdfio 

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk, jvar
  INTEGER   :: narg, iargc, nvar
  INTEGER   :: npiglo,npjglo, npk                           !: size of the domain

  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: v2d, glam, gphi
  REAL(KIND=4) , DIMENSION (:), ALLOCATABLE :: dep
  REAL(KIND=4) ,DIMENSION(1)                  :: timean

  CHARACTER(LEN=80) :: ctag, confcase

  ! Dimg stuff
  INTEGER   :: irecl, ii, nt, ndim
  INTEGER   :: numu=10, numv=11, numt=12,  nums=14, num2d=15, numssh=16, numuu=17, numvv=18
  CHARACTER(LEN=80) :: cdimgu, cdimgv,cdimgt, cdimgs, cdimg2d !: file name dimg
  CHARACTER(LEN=80) :: cdimguu, cdimgvv, cdimgssh             !: file name dimg (optional)
  CHARACTER(LEN=80) :: cheader
  CHARACTER(LEN=4) :: cver
  REAL(KIND=4) :: x1,y1, dx,dy, spval
  LOGICAL :: lexist

  ! Netcdf Stuff
  CHARACTER(LEN=80) :: cfilu ,cfilv ,cfilt, cfilbsf                   !: file name nc
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc', coordzgr='mesh_zgr.nc'
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvar
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipk, id_varout
  INTEGER    :: ncout
  INTEGER    :: istatus

  !!  Read command line
  narg= iargc()
  IF ( narg /= 2 ) THEN
     PRINT *,' Usage : cdfconvert ''Clipper tag '' ''CLIPPER confcase'' '
     PRINT *,'    Output on gridT.nc, gridU.nc and gridV.nc '
     PRINT *,'    mesh_hgr and mesh_zgr must be in the current directory'
     STOP
  ENDIF
  !!
  CALL getarg (1, ctag)
  CALL getarg (2, confcase)

  !! Build dimg file names
  cdimgu=TRIM(confcase)//'_U_'//TRIM(ctag)//'.dimg'
  cdimgv=TRIM(confcase)//'_V_'//TRIM(ctag)//'.dimg'
  cdimgt=TRIM(confcase)//'_T_'//TRIM(ctag)//'.dimg'
  cdimgs=TRIM(confcase)//'_S_'//TRIM(ctag)//'.dimg'
  cdimg2d=TRIM(confcase)//'_2D_'//TRIM(ctag)//'.dimg'

  cdimgssh=TRIM(confcase)//'_SSH_'//TRIM(ctag)//'.dimg'
  cdimguu=TRIM(confcase)//'_UU_'//TRIM(ctag)//'.dimg'
  cdimgvv=TRIM(confcase)//'_VV_'//TRIM(ctag)//'.dimg'

  cfilu=TRIM(confcase)//'_'//TRIM(ctag)//'_gridU.nc'
  cfilv=TRIM(confcase)//'_'//TRIM(ctag)//'_gridV.nc'
  cfilt=TRIM(confcase)//'_'//TRIM(ctag)//'_gridT.nc'
  cfilbsf=TRIM(confcase)//'_'//TRIM(ctag)//'_PSI.nc'

  ! open (and check ?? if they exists )
  irecl=isdirect(cdimgu)  ; OPEN( numu,FILE=cdimgu, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
  irecl=isdirect(cdimgv)  ; OPEN( numv,FILE=cdimgv, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
  irecl=isdirect(cdimgt)  ; OPEN( numt,FILE=cdimgt, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
  irecl=isdirect(cdimgs)  ; OPEN( nums,FILE=cdimgs, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
  irecl=isdirect(cdimg2d) ; OPEN( num2d,FILE=cdimg2d, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
  
  READ(numt,REC=1) cver, cheader, ii, npiglo, npjglo, npk

  ALLOCATE (v2d(npiglo, npjglo), glam(npiglo,npjglo), gphi(npiglo,npjglo), dep(npk) )
  READ(numt,REC=1) cver, cheader, ii, npiglo, npjglo, npk, nt, ndim, &
                        x1,y1,dx,dy,spval,   &
                        (dep(jk),jk=1,npk), &
                        timean(1)
  ! transform Clipper days to drakkar seconds ...
  timean(1)=timean(1)*86400.

  ! Build gridT file with votemper, vosaline, sossheig, ... fluxes ...
  INQUIRE(FILE=cdimgssh, EXIST=lexist)
  IF ( lexist ) THEN
    irecl=isdirect(cdimgssh); OPEN( numssh,FILE=cdimgssh, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
    nvar=10 
  ELSE
    nvar=9
  ENDIF

  ALLOCATE ( typvar(nvar), ipk(nvar), id_varout(nvar) )
  jvar=1
  ipk(jvar)      = npk
  typvar(jvar)%name='votemper'
  typvar(jvar)%units='C'
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= -2.
  typvar(jvar)%valid_max= 40.
  typvar(jvar)%long_name='Potential Temperature'
  typvar(jvar)%short_name='votemper'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TZYX'
  jvar=jvar+1

  ipk(jvar)      = npk
  typvar(jvar)%name='vosaline'
  typvar(jvar)%units='PSU'
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= 0.
  typvar(jvar)%valid_max= 45.
  typvar(jvar)%long_name='Salinity'
  typvar(jvar)%short_name='vosaline'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TZYX'
  jvar=jvar+1

  IF ( lexist ) THEN
  ipk(jvar)      = 1
  typvar(jvar)%name='sossheig'
  typvar(jvar)%units='m'
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= -10.
  typvar(jvar)%valid_max= 10.
  typvar(jvar)%long_name='Sea_Surface_height'
  typvar(jvar)%short_name='sossheig'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TYX'
  jvar=jvar+1
  ENDIF

  ipk(jvar)      = 1
  typvar(jvar)%name='somxl010'         ! rec 12 of dimg file 2D
  typvar(jvar)%units='m'
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= 0.
  typvar(jvar)%valid_max= 7000.
  typvar(jvar)%long_name='Mixed_Layer_Depth_on_0.01_rho_crit'
  typvar(jvar)%short_name='somxl010'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TYX'
  jvar=jvar+1

  ipk(jvar)      = 1
  typvar(jvar)%name='sohefldo'         ! rec 4 of dimg file 2D
  typvar(jvar)%units='W/m2'
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= -1000.
  typvar(jvar)%valid_max= 1000.
  typvar(jvar)%long_name='Net_Downward_Heat_Flux'
  typvar(jvar)%short_name='sohefldo'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TYX'
  jvar=jvar+1

  ipk(jvar)      = 1
  typvar(jvar)%name='soshfldo'         ! rec 8 of dimg file 2D (qsr)
  typvar(jvar)%units='W/m2'
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= -1000.
  typvar(jvar)%valid_max= 1000.
  typvar(jvar)%long_name='Short_Wave_Radiation'
  typvar(jvar)%short_name='soshfldo'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TYX'
  jvar=jvar+1

  ipk(jvar)      = 1
  typvar(jvar)%name='sowaflup'         ! rec 5 of dimg file 2D (emp)
  typvar(jvar)%units='kg/m2/s'         ! conversion required from CLIPPER /86400.
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= -1000.
  typvar(jvar)%valid_max= 1000.
  typvar(jvar)%long_name='Net_Upward_Water_Flux'
  typvar(jvar)%short_name='sowaflup'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TYX'
  jvar=jvar+1

  ipk(jvar)      = 1
  typvar(jvar)%name='sowafldp'         ! rec 10 of dimg file 2D (erp)
  typvar(jvar)%units='kg/m2/s'         ! conversion required from CLIPPER /jvar.
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= -1000.
  typvar(jvar)%valid_max= 1000.
  typvar(jvar)%long_name='Surface_Water_Flux:Damping'
  typvar(jvar)%short_name='sowafldp'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TYX'
  jvar=jvar+1

  ipk(jvar)      = 1
  typvar(jvar)%name='soicecov'         ! rec 13 of dimg file 2D (erp)
  typvar(jvar)%units='%'         
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= 0.
  typvar(jvar)%valid_max= 1.
  typvar(jvar)%long_name='Ice Cover'
  typvar(jvar)%short_name='soicecov'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TYX'
  jvar=jvar+1

  ipk(jvar)      = 1
  typvar(jvar)%name='sohefldp'         ! rec 9 of dimg file 2D (erp)
  typvar(jvar)%units='W/m2'         
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= -10.
  typvar(jvar)%valid_max= 10.
  typvar(jvar)%long_name='Surface Heat Flux: Damping'
  typvar(jvar)%short_name='sohefldp'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TYX'

  glam=getvar(coordhgr,'glamt',1,npiglo,npjglo)
  gphi=getvar(coordhgr,'gphit',1,npiglo,npjglo)
  dep=getvare3(coordzgr,'gdept',npk)

  ncout =create(cfilt, 'none',npiglo,npjglo,npk,cdep='deptht' )
  istatus= createvar(ncout ,typvar,nvar, ipk,id_varout )
  istatus= putheadervar(ncout, 'none', npiglo, npjglo,npk,&
                     pnavlon=glam,pnavlat=gphi,pdep=dep )

  jvar=1
  ! T
  DO jk=1, npk
   READ(numt,REC=jk+1) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
   istatus=putvar(ncout, id_varout(jvar),v2d, jk, npiglo, npjglo)
  END DO
  jvar=jvar+1

  print *, 'Done for T'

  ! S
  DO jk=1, npk
   READ(nums,REC=jk+1) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
   istatus=putvar(ncout, id_varout(jvar),v2d, jk, npiglo, npjglo)
  END DO
  jvar=jvar+1
  print *, 'Done for S'

  IF ( lexist ) THEN
  ! SSH
  READ(numssh,REC=2) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  istatus=putvar(ncout, id_varout(jvar),v2d, 1, npiglo, npjglo)
  jvar=jvar+1
  print *, 'Done for SSH'
  ENDIF

  ! MXL
  READ(num2d,REC=12) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  istatus=putvar(ncout, id_varout(jvar),v2d, 1, npiglo, npjglo)
  jvar=jvar+1
  print *, 'Done for MXL'

  ! QNET
  READ(num2d,REC=4) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  istatus=putvar(ncout, id_varout(jvar),v2d, 1, npiglo, npjglo)
  jvar=jvar+1
  print *, 'Done for QNET'

  ! QSR
  READ(num2d,REC=8) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  istatus=putvar(ncout, id_varout(jvar),v2d, 1, npiglo, npjglo)
  jvar=jvar+1
  print *, 'Done for QSR'

  ! EMP
  READ(num2d,REC=5) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  v2d=v2d/86400. ! to change units
  istatus=putvar(ncout, id_varout(jvar),v2d, 1, npiglo, npjglo)
  jvar=jvar+1
  print *, 'Done for EMP'

  ! ERP
  READ(num2d,REC=10) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  v2d=v2d/86400. ! to change units
  istatus=putvar(ncout, id_varout(jvar),v2d, 1, npiglo, npjglo)
  jvar=jvar+1
  print *, 'Done for ERP'

  ! FREEZE
  READ(num2d,REC=13) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  istatus=putvar(ncout, id_varout(jvar),v2d, 1, npiglo, npjglo)
  jvar=jvar+1
  print *, 'Done for FREEZE'

  ! QRP
  READ(num2d,REC=9) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  istatus=putvar(ncout, id_varout(jvar),v2d, 1, npiglo, npjglo)
  jvar=jvar+1
  print *, 'Done for QRP'

  istatus=putvar1d(ncout,timean,1,'T')
  istatus=CLOSEOUT(ncout)
  DEALLOCATE ( typvar, ipk, id_varout )


!!!!! GRID U !!!!!
 ! Build gridU file with vozocrtx, sozotaux
  INQUIRE(FILE=cdimguu, EXIST=lexist)
  IF ( lexist ) THEN
    irecl=isdirect(cdimguu); OPEN( numuu,FILE=cdimguu, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
    nvar=3 
  ELSE
    nvar=2
  ENDIF
  ALLOCATE ( typvar(nvar), ipk(nvar), id_varout(nvar) )
  
  jvar=1
  ipk(jvar)      = npk
  typvar(jvar)%name='vozocrtx'
  typvar(jvar)%units='m/s'
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= 0.
  typvar(jvar)%valid_max= 20.
  typvar(jvar)%long_name='Zonal Velocity '
  typvar(jvar)%short_name='vozocrtx'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TZYX'
  jvar=jvar+1

  ipk(jvar)      = 1
  typvar(jvar)%name='sozotaux'
  typvar(jvar)%units='N/m2'
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= 0.
  typvar(jvar)%valid_max= 20.
  typvar(jvar)%long_name='Zonal Wind Stress'
  typvar(jvar)%short_name='sozotaux'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TYX'
  jvar=jvar+1

  IF ( lexist ) THEN
  ipk(jvar)      = npk
  typvar(jvar)%name='vozocrtx_sqd'
  typvar(jvar)%units='m2/s2'
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= 0.
  typvar(jvar)%valid_max= 100.
  typvar(jvar)%long_name='MS_Zonal_Velocity'
  typvar(jvar)%short_name='vozocrtx_sqd'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TZYX'
  ENDIF

  glam=getvar(coordhgr,'glamu',1,npiglo,npjglo)
  gphi=getvar(coordhgr,'gphiu',1,npiglo,npjglo)
  dep=getvare3(coordzgr,'gdept',npk)

  ncout =create(cfilu, 'none',npiglo,npjglo,npk,cdep='depthu' )
  istatus= createvar(ncout ,typvar,nvar, ipk,id_varout )
  istatus= putheadervar(ncout, 'none', npiglo, npjglo,npk,&
                     pnavlon=glam,pnavlat=gphi,pdep=dep )

  jvar=1
  DO jk=1, npk
   READ(numu,REC=jk+1) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
   istatus=putvar(ncout, id_varout(jvar),v2d, jk, npiglo, npjglo)
  END DO
  jvar=jvar+1
  print *, 'Done for U'

  READ(num2d, REC=2 ) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  istatus=putvar(ncout, id_varout(jvar),v2d, 1, npiglo, npjglo)
  jvar=jvar+1
  print *, 'Done for TAUX'

  IF ( lexist ) THEN
  DO jk=1, npk
   READ(numuu,REC=jk+1) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
   istatus=putvar(ncout, id_varout(jvar),v2d, jk, npiglo, npjglo)
  END DO
  print *, 'Done for UU'
  ENDIF


  istatus=putvar1d(ncout,timean,1,'T')
  istatus=CLOSEOUT(ncout)
  DEALLOCATE ( typvar, ipk, id_varout )

!!!!! GRID V !!!!!
  ! Build gridV file with vomecrty, sometauy
  INQUIRE(FILE=cdimgvv, EXIST=lexist)
  IF ( lexist ) THEN
    irecl=isdirect(cdimgvv); OPEN( numvv,FILE=cdimgvv, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
    nvar=3 
  ELSE
    nvar=2
  ENDIF
  ALLOCATE ( typvar(nvar), ipk(nvar), id_varout(nvar) )

  jvar=1
  ipk(jvar)      = npk
  typvar(jvar)%name='vomecrty'
  typvar(jvar)%units='m/s'
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= 0.
  typvar(jvar)%valid_max= 20.
  typvar(jvar)%long_name='Meridinal  Velocity '
  typvar(jvar)%short_name='vomecrty'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TZYX'
  jvar=jvar+1

  ipk(jvar)      = 1
  typvar(jvar)%name='sometauy'
  typvar(jvar)%units='N/m2'
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= 0.
  typvar(jvar)%valid_max= 20.
  typvar(jvar)%long_name='Meridional Wind Stress'
  typvar(jvar)%short_name='sometauy'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TYX'
  jvar=jvar+1

  IF ( lexist ) THEN
  ipk(jvar)      = npk
  typvar(jvar)%name='vomecrty_sqd'
  typvar(jvar)%units='m2/s2'
  typvar(jvar)%missing_value=0.
  typvar(jvar)%valid_min= 0.
  typvar(jvar)%valid_max= 100.
  typvar(jvar)%long_name='MS_Meridional_Velocity'
  typvar(jvar)%short_name='vomecrty_sqd'
  typvar(jvar)%online_operation='N/A'
  typvar(jvar)%axis='TZYX'
  ENDIF


  glam=getvar(coordhgr,'glamv',1,npiglo,npjglo)
  gphi=getvar(coordhgr,'gphiv',1,npiglo,npjglo)
  dep=getvare3(coordzgr,'gdept',npk)
  
  ncout =create(cfilv, 'none',npiglo,npjglo,npk,cdep='depthv' )
  istatus= createvar(ncout ,typvar,nvar, ipk,id_varout )
  istatus= putheadervar(ncout, 'none', npiglo, npjglo,npk,&
                     pnavlon=glam,pnavlat=gphi,pdep=dep )

  DO jk=1, npk
   READ(numv,REC=jk+1) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
   istatus=putvar(ncout, id_varout(1),v2d, jk, npiglo, npjglo)
  END DO
  print *, 'Done for V'

  READ(num2d, REC=3 ) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  istatus=putvar(ncout, id_varout(2),v2d, 1, npiglo, npjglo)
  print *, 'Done for TAUY'

  IF ( lexist ) THEN
  DO jk=1, npk
   READ(numvv,REC=jk+1) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
   istatus=putvar(ncout, id_varout(jvar),v2d, jk, npiglo, npjglo)
  END DO
  print *, 'Done for VV'
  ENDIF

  istatus=putvar1d(ncout,timean,1,'T')
  istatus=CLOSEOUT(ncout)

  DEALLOCATE ( typvar, ipk, id_varout )

!!!!! PSI !!!!!
  ! Build PSI file with sobarstf
  nvar=1  
  ALLOCATE ( typvar(nvar), ipk(nvar), id_varout(nvar) )
  ipk(1)      = 1
  typvar(1)%name='sobarstf'
  typvar(1)%units='m3/s'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= -3.e8
  typvar(1)%valid_max= 3.e8
  typvar(1)%long_name='Barotropic_Stream_Function'
  typvar(1)%short_name='sobarstf'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TYX'

  glam=getvar(coordhgr,'glamf',1,npiglo,npjglo)
  gphi=getvar(coordhgr,'gphif',1,npiglo,npjglo)
  dep=getvare3(coordzgr,'gdept',1)
  
  ncout =create(cfilbsf, 'none',npiglo,npjglo,1,cdep='depthu' )
  istatus= createvar(ncout ,typvar,nvar, ipk,id_varout )
  istatus= putheadervar(ncout, 'none', npiglo, npjglo,1,&
                     pnavlon=glam,pnavlat=gphi,pdep=dep )

   READ(num2d,REC=7) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
   istatus=putvar(ncout, id_varout(1),v2d, 1, npiglo, npjglo)
  print *, 'Done for PSI'

  istatus=putvar1d(ncout,timean,1,'T')
  istatus=CLOSEOUT(ncout)

  DEALLOCATE ( typvar, ipk, id_varout )

CONTAINS
        INTEGER FUNCTION isdirect(clname)
!!!                     FUNCTION ISDIRECT
!!!                     *****************
!!!
!!!    PURPOSE : This integer function returns the record length if clname
!!!              is a valid dimg file, it returns 0 either.
!!!
!!!    METHOD : Open the file and look for the key characters (@!01) for
!!!             identification.
!!!
!!!    AUTHOR : Jean-Marc Molines (Apr. 1998)
!!! -------------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) ::  clname
      CHARACTER(LEN=4)  ::  cver
      CHARACTER(LEN=80) ::  clheader
!
      INTEGER :: irecl

!
      OPEN(100,FILE=clname, FORM   ='UNFORMATTED', ACCESS ='DIRECT', RECL   =88)
       READ(100,REC=1) cver ,clheader,irecl
      CLOSE(100)
!
      IF (cver ==  '@!01' ) THEN
       isdirect=irecl
      ELSE
       isdirect=0
      END IF
!
      END FUNCTION isdirect
END PROGRAM cdfconvert
