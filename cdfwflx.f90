PROGRAM cdfwflx
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfwflx  ***
  !!
  !!  **  Purpose  :  Produce a file with the water flux separated into 4 components:
  !!                   E (soevap), P (soprecip), R (sorunoff), dmp (sowafldp).
  !!                   The total water flux is E -P -R + dmp. Units in this program
  !!                  are mm/days. 
  !!  
  !!  **  Method   :  Evap is computed from the latent heat flux : evap=-qla/Lv 
  !!                  Runoff is read from the climatological input file
  !!                  dmp is read from the file (sowafldp)
  !!                  Precip is then computed as the difference between the
  !!                  total water flux (sowaflup) and the E-R+dmp. In the high latitudes
  !!                  this precip includes the effect of snow (storage/melting). Therefore
  !!                  it may differ slightly from the input precip file.
  !!
  !! history ;
  !!  Original :  J.M. Molines (January 2008 )
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jj, jk ,ji                  !: dummy loop index
  INTEGER   :: narg, iargc                 !: command line 
  INTEGER   :: npiglo,npjglo               !: size of the domain

  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  zmask, zwk                 !:  work array
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  evap, precip, runoff, wdmp !: water flux components
  REAL(KIND=4)                                   ::  Lv=2.5e6                    !: latent HF <--> evap conversion

  CHARACTER(LEN=256) :: cfilet , cfiler

  INTEGER    :: istatus
  ! output stuff
  INTEGER, PARAMETER              :: jpvarout=5
  INTEGER                         :: ncout, ierr
  INTEGER,    DIMENSION(jpvarout) :: ipk, id_varout  !: only one output variable
  REAL(KIND=4),      DIMENSION(1) :: tim,dep       !: time output
  CHARACTER(LEN=256)               :: cfileout='wflx.nc'

  TYPE(variable), DIMENSION(jpvarout) :: typvar        !: structure for attributes


  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfwflx  Tfile Runoff file'
     PRINT *,' Computes the water fluxes components'
     PRINT *,' Output on wflx.nc, soevap,soprecip,sorunoff,sowadmp,sowaflup'
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  CALL getarg (2, cfiler)
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')

  ! prepare output variables
  dep(1) = 0.
  ipk(:)= 1  ! all variables ( output are 2D)

  typvar(1)%name= 'soevap'
  typvar(2)%name= 'soprecip'
  typvar(3)%name= 'sorunoff'
  typvar(4)%name= 'sowadmp'
  typvar(5)%name= 'sowaflux'
  typvar%units='mm/day'
  typvar%missing_value=0.
  typvar%valid_min= -100.
  typvar%valid_max= 100.
  typvar(1)%long_name='Evaporation'
  typvar(2)%long_name='Precipitation'
  typvar(3)%long_name='Runoff'
  typvar(4)%long_name='SSS damping'
  typvar(5)%long_name='Total water flux'
  typvar(1)%short_name='soevap'
  typvar(2)%short_name='soprecip'
  typvar(3)%short_name='sorunoff'
  typvar(4)%short_name='sowadmp'
  typvar(5)%short_name='sowaflux'
  typvar%online_operation='N/A'
  typvar%axis='TYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo


  ALLOCATE ( zmask(npiglo,npjglo), zwk(npiglo,npjglo))
  ALLOCATE ( evap(npiglo,npjglo), precip(npiglo,npjglo), runoff(npiglo,npjglo), wdmp(npiglo,npjglo) )

 ! read vosaline for masking purpose
     zwk(:,:) =  getvar(cfilet, 'vosaline',  1 ,npiglo,npjglo)
     zmask=1. ; WHERE ( zwk == 0 ) zmask=0.

 ! Evap : 
     evap(:,:)= -1.* getvar(cfilet, 'solhflup',  1 ,npiglo,npjglo)/Lv*86400. *zmask(:,:)  ! mm/days
      print *,'Evap done'
 ! Wdmp
     wdmp(:,:)= getvar(cfilet, 'sowafldp',  1 ,npiglo,npjglo)*86400.*zmask(:,:)           ! mm/days
      print *,'Damping done'
 ! Runoff
     runoff(:,:)= getvar(cfiler, 'sorunoff',  1 ,npiglo,npjglo)*86400.*zmask(:,:)         ! mm/days
      print *,'Runoff done'
 ! total water flux
     zwk(:,:) =  getvar(cfilet, 'sowaflup',  1 ,npiglo,npjglo)*86400.*zmask(:,:)          ! mm/days
      print *,'Total water flux done'
 ! Precip:
     precip(:,:)= evap(:,:)-runoff(:,:)+wdmp(:,:)-zwk(:,:)                     ! mm/day
      print *,'Precip done'

 ! Write output file
 !
 ncout = create(cfileout, cfilet, npiglo,npjglo,1)
 ierr = createvar(ncout ,typvar ,jpvarout, ipk,id_varout )
 ierr= putheadervar(ncout, cfilet,npiglo, npjglo,1,pdep=dep)
 tim=getvar1d(cfilet,'time_counter',1)
 ierr = putvar(ncout, id_varout(1) ,evap,   1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(2) ,precip, 1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(3) ,runoff, 1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(4) ,wdmp,   1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(5) ,zwk,    1,npiglo, npjglo)
 ierr=putvar1d(ncout,tim,1,'T')

 ierr=closeout(ncout)

   END PROGRAM cdfwflx
