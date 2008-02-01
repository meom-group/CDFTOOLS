PROGRAM cdfbuoyflx
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfbuoyflx  ***
  !!
  !!  **  Purpose  :  Produce a file with the water flux separated into 4 components:
  !!                   E (evap), P (precip), R (runoff), dmp (sssdmp).
  !!                   The total water flux is E -P -R + dmp. Units in this program
  !!                  are mm/days.  (Up to that it is the same than cdfwflx)
  !!
  !!                   It also produces un the same file the component of the heat flux 
  !!                   Latent Heat FLux, Sensible Heat flux, Long Wave HF, Short Wave HF, Net HF 
  !!                  
  !!                  Buoyancy fluxes are also computed, as a net value but also with the 
  !!                  contribution of each term.
  !!  
  !!  **  Method   :  Evap is computed from the latent heat flux : evap=-qla/Lv 
  !!                  Runoff is read from the climatological input file
  !!                  dmp is read from the file (sowafldp)
  !!                  Precip is then computed as the difference between the
  !!                  total water flux (sowaflup) and the E-R+dmp. In the high latitudes
  !!                  this precip includes the effect of snow (storage/melting). Therefore
  !!                  it may differ slightly from the input precip file.
  !!                  
  !!                  Heat fluxes are directly copied from the gridT files, same name, same units
  !!                  We also add sst and SSS for convenience.
  !!                  
  !!                  Buoyancy fluxes are also computed as :
  !!                     BF = -1/rho (alpha x TF  - beta SF )
  !!                       (TF = thermal part, SF = haline part )
  !!                     TF = 1/(rho x Cp)* Q
  !!                     SF = 1/(1-SSS) x (E-P) x SSS
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
  USE eos

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jj, jk ,ji                  !: dummy loop index
  INTEGER   :: narg, iargc                 !: command line 
  INTEGER   :: npiglo,npjglo               !: size of the domain

  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  zmask,zcoefq,zcoefw              !:  work array
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  zalbet, zbeta                    !:  work array
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  evap, precip, runoff, wdmp, wnet !: water flux components
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  wice, precip_runoff               !: water flux components
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  qlat, qsb,    qlw,    qsw,  qnet !: heat flux components
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  b_evap, b_precip, b_runoff, b_wdmp, bw_net  !: BF water flux components
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  b_qlat, b_qsb,    b_qlw,    b_qsw , bh_net  !: BF heat flux components
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  sst, sss, buoyancy_fl                      !: Total buoyancy flux

  ! Physical constants
  REAL(KIND=4)                                   ::  Lv=2.5e6                   !: latent HF <--> evap conversion
  REAL(KIND=4)                                   ::  Cp = 4000.                 !: specific heat of water 

  CHARACTER(LEN=80) :: cfilet , cfiler

  INTEGER    :: istatus
  ! output stuff
  INTEGER, PARAMETER              :: jpvarout=25
  INTEGER                         :: ncout, ierr
  INTEGER,    DIMENSION(jpvarout) :: ipk, id_varout  !: only one output variable
  REAL(KIND=4),      DIMENSION(1) :: tim,dep       !: time output
  CHARACTER(LEN=80)               :: cfileout='buoyflx.nc'

  TYPE(variable), DIMENSION(jpvarout) :: typvar        !: structure for attributes


  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfbuoyflx  Tfile Runoff file'
     PRINT *,' produces the water fluxes components'
     PRINT *,' produces the heat fluxes components'
     PRINT *,' produces the net fluxes'
     PRINT *,' produces the buoyancy water fluxes components'
     PRINT *,' produces the buoyancy heat fluxes components'
     PRINT *,' produces the buoyancy net fluxes'
     PRINT *,' produces the sss and sst '
     PRINT *,' Output on buoyflx.nc , 25 variables (2D) '
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  CALL getarg (2, cfiler)
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')

  ! prepare output variables
  dep(1) = 0.
  ipk(:)= 1  ! all variables ( output are 2D)
  typvar%online_operation='N/A'
  typvar%axis='TYX'

  !  
  ! 1--> 7 water fluxes                     ;   ! 8 --> 12    heat fluxes
  typvar(1)%name= 'evap'                  ;  typvar(8)%name= 'latent'      
  typvar(2)%name= 'precip'                ;  typvar(9)%name= 'sensible'      
  typvar(3)%name= 'runoff'                ;  typvar(10)%name= 'longwave'      
  typvar(4)%name= 'sssdmp'                ;  typvar(11)%name= 'solar'      
  typvar(5)%name= 'watnet'                ;  typvar(12)%name= 'heatnet'      
  typvar(6)%name= 'wice'            
  typvar(7)%name= 'precip_runoff'  
  typvar(1:7)%units='mm/day'              ;  typvar(8:12)%units='W/m2'      
  typvar(1:7)%missing_value=0.            ;  typvar(8:12)%missing_value=0.      
  typvar(1:7)%valid_min= -100.            ;  typvar(8:12)%valid_min= -500.      
  typvar(1:7)%valid_max= 100.             ;  typvar(8:12)%valid_max= 500.      
  typvar(1)%long_name='Evaporation'       ;  typvar(8)%long_name='Latent Heat flux'      
  typvar(2)%long_name='Precipitation'     ;  typvar(9)%long_name='Sensible Heat flux'       
  typvar(3)%long_name='Runoff'            ;  typvar(10)%long_name='Long Wave Heat flux'      
  typvar(4)%long_name='SSS damping'       ;  typvar(11)%long_name='Short Wave Heat flux'
  typvar(5)%long_name='Total water flux'  ;  typvar(12)%long_name='Net Heat Flux'      
  typvar(6)%long_name='Ice congelation and melting'  
  typvar(7)%long_name='Precip and runoff together' 
  typvar(1)%short_name='evap'             ;  typvar(8)%short_name='latent'      
  typvar(2)%short_name='precip'           ;  typvar(9)%short_name='sensible'      
  typvar(3)%short_name='runoff'           ;  typvar(10)%short_name='longwave'      
  typvar(4)%short_name='sssdmp'           ;  typvar(11)%short_name='solar'      
  typvar(5)%short_name='watnet'           ;  typvar(12)%short_name='heatnet'       
  typvar(6)%short_name='wice'  
  typvar(7)%short_name='precip_runoff' 

  ! 13--> 17  buoy water fluxes             ;   ! 18 --> 22    buoy heat fluxes
  typvar(13)%name= 'evap_b'               ;  typvar(18)%name= 'latent_b'
  typvar(14)%name= 'precip_b'             ;  typvar(19)%name= 'sensible_b'
  typvar(15)%name= 'runoff_b'             ;  typvar(20)%name= 'longwave_b'
  typvar(16)%name= 'sssdmp_b'             ;  typvar(21)%name= 'solar_b'
  typvar(17)%name= 'watnet_b'             ;  typvar(22)%name= 'heatnet_b'
  typvar(13:17)%units='1e-6 kg/m2/s'      ;  typvar(18:22)%units='1e-6 kg/m2/s'
  typvar(13:17)%missing_value=0.          ;  typvar(18:22)%missing_value=0.
  typvar(13:17)%valid_min= -100.          ;  typvar(18:22)%valid_min= -500.
  typvar(13:17)%valid_max= 100.           ;  typvar(18:22)%valid_max= 500.
  typvar(13)%long_name='buoy flx evap'    ;  typvar(18)%long_name='buoy Latent Heat flux'
  typvar(14)%long_name='buoy flx precip'  ;  typvar(19)%long_name='buoy Sensible Heat flux'
  typvar(15)%long_name='buoy flx runoff'  ;  typvar(20)%long_name='buoy Long Wave Heat flux'
  typvar(16)%long_name='buoy flx damping' ;  typvar(21)%long_name='buoy Short Wave Heat flux'
  typvar(17)%long_name='buoy haline flx'  ;  typvar(22)%long_name='buoy thermo Flux'
  typvar(13)%short_name='evap_b'          ;  typvar(18)%short_name='latent_b'
  typvar(14)%short_name='precip_b'        ;  typvar(19)%short_name='sensible_b'
  typvar(15)%short_name='runoff_b'        ;  typvar(20)%short_name='longwave_b'
  typvar(16)%short_name='sssdmp_b'        ;  typvar(21)%short_name='solar_b'
  typvar(17)%short_name='watnet_b'        ;  typvar(22)%short_name='heatnet_b'

  ! total buoyancy flux
  typvar(23)%name= 'buoyancy_fl'
  typvar(23)%units='1e-6 kg/m2/s'
  typvar(23)%missing_value=0.
  typvar(23)%valid_min= -100.
  typvar(23)%valid_max= 100.
  typvar(23)%long_name='buoyancy flux'
  typvar(23)%short_name='buoyancy_fl'

  ! SSS                                       ; SST
  typvar(24)%name= 'sss'                      ;   typvar(25)%name= 'sst'
  typvar(24)%units='PSU'                      ;   typvar(25)%units='Celsius'
  typvar(24)%missing_value=0.                 ;   typvar(25)%missing_value=0.
  typvar(24)%valid_min= 0.                    ;   typvar(25)%valid_min= -2.
  typvar(24)%valid_max= 45                    ;   typvar(25)%valid_max= 45
  typvar(24)%long_name='Sea Surface Salinity' ;   typvar(25)%long_name='Sea Surface Temperature'
  typvar(24)%short_name='sss  '               ;   typvar(25)%short_name='sst'

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo


  ALLOCATE ( zmask(npiglo,npjglo), wnet(npiglo,npjglo), zalbet(npiglo,npjglo), zbeta(npiglo, npjglo) )
  ALLOCATE ( zcoefq(npiglo,npjglo), zcoefw(npiglo,npjglo) )
  ALLOCATE ( evap(npiglo,npjglo), precip(npiglo,npjglo), runoff(npiglo,npjglo), wdmp(npiglo,npjglo) )
  ALLOCATE ( wice(npiglo,npjglo), precip_runoff(npiglo,npjglo) )
  ALLOCATE ( qlat(npiglo,npjglo), qsb(npiglo,npjglo), qlw(npiglo,npjglo), qsw(npiglo,npjglo), qnet(npiglo,npjglo) )
  ALLOCATE ( b_evap(npiglo,npjglo), b_precip(npiglo,npjglo), b_runoff(npiglo,npjglo), b_wdmp(npiglo,npjglo),bw_net(npiglo,npjglo) ) 
  ALLOCATE ( b_qlat(npiglo,npjglo), b_qsb(npiglo,npjglo),    b_qlw(npiglo,npjglo),    b_qsw(npiglo,npjglo), bh_net(npiglo,npjglo))
  ALLOCATE ( buoyancy_fl(npiglo,npjglo), sst(npiglo,npjglo), sss(npiglo,npjglo) )

 ! read sss for masking purpose and sst
     sss(:,:) =  getvar(cfilet, 'vosaline',  1 ,npiglo,npjglo)
     zmask=1. ; WHERE ( sss == 0 ) zmask=0.
     sst(:,:) =  getvar(cfilet, 'votemper',  1 ,npiglo,npjglo)

 ! Evap : 
     qlat(:,:)= getvar(cfilet, 'solhflup',  1 ,npiglo,npjglo)*zmask(:,:)                  ! W/m2 
     evap(:,:)= -1.* qlat(:,:) /Lv*86400. *zmask(:,:)                                     ! mm/days
      print *,'Evap done'
 ! Wdmp
     wdmp(:,:)= getvar(cfilet, 'sowafldp',  1 ,npiglo,npjglo)*86400.*zmask(:,:)           ! mm/days
      print *,'Damping done'
 ! Runoff
     runoff(:,:)= getvar(cfiler, 'sorunoff',  1 ,npiglo,npjglo)*86400.*zmask(:,:)         ! mm/days
      print *,'Runoff done'
 ! total water flux (emps)
     wnet(:,:) =  getvar(cfilet, 'sowaflcd',  1 ,npiglo,npjglo)*86400.*zmask(:,:)          ! mm/days
      print *,'Total water flux done'
 ! fsalt = contribution of ice freezing and melting to salinity ( + = freezing, - = melting )Q
     wice(:,:) = getvar(cfilet, 'iowaflup',  1 ,npiglo,npjglo)*86400.*zmask(:,:)          ! mm/days
      print *,'ice contribution done'
 ! Precip:
     precip(:,:)= evap(:,:)-runoff(:,:)+wdmp(:,:)-wnet(:,:)+wice(:,:)                     ! mm/day
      print *,'Precip done'
 ! Precip+runoff : (as a whole ) (interpolated on line)
     precip_runoff(:,:)= evap(:,:)+wdmp(:,:)-wnet(:,:)+wice(:,:)                          ! mm/day
      print *,'Precip done'
  ! other heat fluxes
     qsb(:,:)= getvar(cfilet, 'sosbhfup',  1 ,npiglo,npjglo)*zmask(:,:)                  ! W/m2 
      print *,'qsb done'
     qlw(:,:)= getvar(cfilet, 'solwfldo',  1 ,npiglo,npjglo)*zmask(:,:)                  ! W/m2 
      print *,'qlw done'
     qsw(:,:)= getvar(cfilet, 'soshfldo',  1 ,npiglo,npjglo)*zmask(:,:)                  ! W/m2 
      print *,'qsw done'
     qnet(:,:)= getvar(cfilet,'sohefldo',  1 ,npiglo,npjglo)*zmask(:,:)                  ! W/m2 
      print *,'qnet done'

  ! buoyancy flux
    zalbet(:,:)= albet ( sst, sss, 0., npiglo,npjglo)
    zbeta (:,:)= beta  ( sst, sss, 0., npiglo,npjglo)
    zcoefq(:,:)= -zbeta * zalbet /Cp * 1.e6
    zcoefw(:,:)=  zbeta*sss/(1-sss/1000.)/86400. *1.e6   ! division by 86400 to get back water fluxes in kg/m2/s
    buoyancy_fl=0. ; bh_net=0. ; b_qlat=0. ; b_qlw=0. ; b_qsw=0. ; b_qsb=0.
                     bw_net=0. ; b_evap=0. ; b_precip=0.; b_wdmp=0. ; b_runoff=0.
    WHERE (sss /= 0 ) 
     bh_net(:,:)= zcoefq * qnet
     b_qlat(:,:)= zcoefq * qlat
     b_qlw (:,:)= zcoefq * qlw
     b_qsw (:,:)= zcoefq * qsw
     b_qsb (:,:)= zcoefq * qsb

     bw_net(:,:)= zcoefw * wnet
     b_evap(:,:)= zcoefw * evap
     b_precip(:,:)= -zcoefw * precip
     b_runoff(:,:)= -zcoefw * runoff
     b_wdmp(:,:)= zcoefw * wdmp

!    buoyancy_fl(:,:) = zcoefq * qnet +zcoefw * wnet
     buoyancy_fl(:,:) = bh_net + bw_net
    END WHERE 

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
 ierr = putvar(ncout, id_varout(5) ,wnet,   1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(6) ,wice,   1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(7) ,precip_runoff,   1,npiglo, npjglo)

 ierr = putvar(ncout, id_varout(8)  ,qlat,  1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(9)  ,qsb,   1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(10)  ,qlw,   1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(11)  ,qsw,   1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(12) ,qnet,  1,npiglo, npjglo)

 ierr = putvar(ncout, id_varout(13) ,b_evap,   1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(14) ,b_precip, 1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(15) ,b_runoff, 1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(16) ,b_wdmp,   1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(17) ,bw_net,   1,npiglo, npjglo)

 ierr = putvar(ncout, id_varout(18) ,b_qlat,  1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(19) ,b_qsb,   1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(20) ,b_qlw,   1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(21) ,b_qsw,   1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(22) ,bh_net,  1,npiglo, npjglo)

 ierr = putvar(ncout, id_varout(23) ,buoyancy_fl,  1,npiglo, npjglo)

 ierr = putvar(ncout, id_varout(24) ,sss,  1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(25) ,sst,  1,npiglo, npjglo)

 ierr=putvar1d(ncout,tim,1,'T')

 ierr=closeout(ncout)

   END PROGRAM cdfbuoyflx
