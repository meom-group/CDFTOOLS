PROGRAM cdfbuoyflx
  !!======================================================================
  !!                     ***  PROGRAM  cdfbuoyflx  ***
  !!=====================================================================
  !!  ** Purpose :  Produce a file with the water flux separated into 4 components:
  !!                E (evap), P (precip), R (runoff), dmp (sssdmp).
  !!                The total water flux is E -P -R + dmp. Units in this program
  !!                are mm/days.  (Up to that it is the same than cdfwflx)
  !!
  !!                It also produces un the same file the component of the heat flux
  !!                Latent Heat FLux, Sensible Heat flux, Long Wave HF, Short Wave HF,
  !!                Net HF
  !!
  !!                Buoyancy fluxes are also computed, as a net value but also with the
  !!                contribution of each term.
  !!
  !!  ** Method  : Evap is computed from the latent heat flux : evap=-qla/Lv
  !!               Runoff is read from the climatological input file
  !!               dmp is read from the file (sowafldp)
  !!               Precip is then computed as the difference between the
  !!               total water flux (sowaflup) and the E-R+dmp. In the high latitudes
  !!               this precip includes the effect of snow (storage/melting). Therefore
  !!               it may differ slightly from the input precip file.
  !!
  !!               Heat fluxes are directly copied from the gridT files, same name, same units
  !!               We also add sst and SSS for convenience.
  !!
  !!               Buoyancy fluxes are also computed as :
  !!                  BF = g/rho ( alpha x TF  - beta x SF ) 
  !!                       (TF = thermal part, SF = haline part )
  !!                  TF = Qnet/cp
  !!                  SF = rho x (E-P) x SSS
  !!  ** Reference :
  !!              Atmosphere, Ocean and Climate Dynamics: An Introductory Text. By John Marshall,
  !!              R. Alan Plumb ( Academic Press, 2008 ) Eq. 11.4 p 225. 
  !!
  !! History : 2.1  : 01/2008  : J.M. Molines : Original code
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!           3.0  : 09/2015  : J.M. Molines : add nc4 capabilities, optional output file
  !!                                            short output,
  !!                                            different management of read fluxes (XIOS ...)
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class forcing
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: np_varout=25
  INTEGER(KIND=4)                           :: ncout, ierr
  INTEGER(KIND=4)                           :: jt                                ! dummy loop index
  INTEGER(KIND=4)                           :: narg, iargc, ijarg                ! command line 
  INTEGER(KIND=4)                           :: npiglo, npjglo, npt               ! size of the domain
  INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:):: ipk, id_varout  

  ! Physical constants
  REAL(KIND=4)                              :: Lv = 2.5e6                        ! latent HF <--> evap conversion
  REAL(KIND=4)                              :: Cp = 4000.                        ! specific heat of water 
  REAL(KIND=4)                              :: Rho = 1026.                       ! reference density
  REAL(KIND=4)                              :: Grav = 9.81                       ! Gravity

  REAL(KIND=4)                              :: zsps                              ! Missing value for salinity
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: zdep                         ! deptht
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask, zcoefq, zcoefw             ! work array
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zalbet, zbeta                     ! work array
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: evap, precip, runoff, wdmp, wnet  ! water flux components
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: wice, precip_runoff               ! water flux components
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: qlat, qsb, qlw, qsw, qnet         ! heat flux components
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: b_evap, b_precip, b_runoff        ! BF water flux components
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: b_wdmp, bw_net                    ! BF water flux components
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: b_qlat, b_qsb, b_qlw              ! BF heat flux components
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: b_qsw , bh_net                    ! BF heat flux components
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsst, zsss, buoyancy_fl           ! Total buoyancy flux

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim                              ! time counter

  CHARACTER(LEN=256)                        :: cf_tfil ,cf_flxfil, cf_rnfil      ! input file gridT, flx and runoff
  CHARACTER(LEN=256)                        :: cf_out='buoyflx.nc'               ! output file
  CHARACTER(LEN=256)                        :: cldum                             ! dummy character variable
  CHARACTER(LEN=256)                        :: cv_sss                            ! Actual name for SSS
  CHARACTER(LEN=256)                        :: cv_sst                            ! Actual name for SST

  TYPE(variable), ALLOCATABLE, DIMENSION(:) :: stypvar                           ! structure for attributes

  LOGICAL                                   :: lchk =.FALSE.                     ! flag for missing files
  LOGICAL                                   :: lnc4 =.FALSE.                     ! flag for netcdf4 output
  LOGICAL                                   :: lsho =.FALSE.                     ! flag for short output
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfbuoyflx  -t T-file [-r RNF-file] [-f FLX-file ] [-sss SSS-name]'
     PRINT *,'     ... [-sst SST-name] [-nc4] [-o OUT-file] [-short ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute (or read) the heat and water fluxes components.'
     PRINT *,'       Compute (or read) the net heat and water fluxes.'
     PRINT *,'       Compute the buoyancy heat and water fluxes components.'
     PRINT *,'       Compute the net buoyancy fluxes.'
     PRINT *,'       Save sss and sst. '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file   : netcdf file with temperature and salinity '
     PRINT *,'      '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-r RNF-file ] : Specify a run-off file if runoff not in T-file '
     PRINT *,'                         nor in FLX-file'
     PRINT *,'       [-f FLX-file ] : Use this option if fluxes are not saved in gridT files'
     PRINT *,'       [-sss SSS-name ] : Use this option if SSS variable name in T-file '
     PRINT *,'                          differ from ',TRIM(cn_vosaline)
     PRINT *,'       [-sst SST-name ] : Use this option if SST variable name in T-file '
     PRINT *,'                          differ from ',TRIM(cn_votemper)
     PRINT *,'       [-nc4 ] Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'               This option is effective only if cdftools are compiled with'
     PRINT *,'               a netcdf library supporting chunking and deflation.'
     PRINT *,'       [-o OUT-file ] Default is ', TRIM(cf_out)
     PRINT *,'       [-short ] With this option only save the buoyancy flux without '
     PRINT *,'                  all the components of the flux.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : 25 variables (2D) or 1 variable in case of -short option'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      '
     STOP 
  ENDIF
  ijarg   = 1
  cf_flxfil='none'
  cf_rnfil='none'
  cv_sss=cn_vosaline
  cv_sst=cn_votemper

  DO   WHILE ( ijarg <= narg ) 
     CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg (ijarg, cf_flxfil) ; ijarg = ijarg + 1
        ;              lchk = lchk .OR. chkfile (cf_flxfil)
     CASE ( '-t'   ) ; CALL getarg (ijarg, cf_tfil) ; ijarg = ijarg + 1
        ;              lchk = lchk .OR. chkfile (cf_tfil  )
     CASE ( '-r'   ) ; CALL getarg (ijarg, cf_rnfil) ; ijarg = ijarg + 1
        ;              lchk = lchk .OR. chkfile (cf_rnfil )
     CASE ( '-sss' ) ; CALL getarg (ijarg, cv_sss) ; ijarg = ijarg + 1
     CASE ( '-sst' ) ; CALL getarg (ijarg, cv_sst) ; ijarg = ijarg + 1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE ( '-o'   ) ; CALL getarg (ijarg, cf_out) ; ijarg = ijarg + 1
     CASE ('-short') ; lsho = .TRUE. ; np_varout = 1
     CASE DEFAULT    ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO
  IF (lchk ) STOP 99 ! missing files

  IF ( cf_flxfil == 'none' ) THEN
     cf_flxfil = cf_tfil
  ENDIF
  ! If no runoff file specified, assume that run off are in flx file [ which must be read by the way ... ]
  IF ( cf_rnfil == 'none' ) THEN
     cf_rnfil = cf_flxfil
  ENDIF

  ! Look for Missing value for salinity
  zsps = getspval(cf_tfil, cn_vosaline)

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npt    = getdim (cf_tfil,cn_t)

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npt    =', npt

  CALL CreateOutput

  ! always allocated
  ALLOCATE ( zmask(npiglo,npjglo), wnet(npiglo,npjglo), zalbet(npiglo,npjglo), zbeta(npiglo, npjglo) )
  ALLOCATE ( zcoefq(npiglo,npjglo), zcoefw(npiglo,npjglo), qnet(npiglo,npjglo) )
  ALLOCATE ( bw_net(npiglo,npjglo) ,bh_net(npiglo,npjglo)) 
  ALLOCATE ( buoyancy_fl(npiglo,npjglo), zsst(npiglo,npjglo), zsss(npiglo,npjglo) )

  ! allocated only for full output
  IF ( .NOT. lsho ) THEN
     ALLOCATE ( evap(npiglo,npjglo), precip(npiglo,npjglo), runoff(npiglo,npjglo), wdmp(npiglo,npjglo) )
     ALLOCATE ( wice(npiglo,npjglo), precip_runoff(npiglo,npjglo) )
     ALLOCATE ( qlat(npiglo,npjglo), qsb(npiglo,npjglo), qlw(npiglo,npjglo), qsw(npiglo,npjglo) )
     ALLOCATE ( b_evap(npiglo,npjglo), b_precip(npiglo,npjglo), b_runoff(npiglo,npjglo), b_wdmp(npiglo,npjglo) ) 
     ALLOCATE ( b_qlat(npiglo,npjglo), b_qsb(npiglo,npjglo),    b_qlw(npiglo,npjglo),    b_qsw(npiglo,npjglo)  )
  ENDIF

  DO jt = 1, npt
     ! read sss for masking purpose and sst
     zsss(:,:) = getvar(cf_tfil, cv_sss, 1, npiglo, npjglo, ktime=jt)
     zmask=1. ; WHERE ( zsss == zsps ) zmask=0.
     zsst(:,:) = getvar(cf_tfil, cv_sst, 1, npiglo, npjglo, ktime=jt)

     ! total water flux (emps)
     wnet(:,:) = getvar(cf_flxfil, cn_sowaflup, 1, npiglo, npjglo, ktime=jt )*86400.*zmask(:,:)          ! mm/days
     qnet(:,:)=  getvar(cf_flxfil, cn_sohefldo, 1, npiglo, npjglo, ktime=jt )*zmask(:,:)    ! W/m2 

     ! buoyancy flux
     zalbet(:,:)= albet ( zsst, zsss, 0., npiglo, npjglo)
     zbeta (:,:)= beta  ( zsst, zsss, 0., npiglo, npjglo)
     zcoefq(:,:)= Grav/Rho *( zbeta * zalbet /Cp ) * 1.e6
     zcoefw(:,:)= Grav* zbeta * zsss / 86400. /1000 * 1.e6   ! division by 86400 and 1000 to get back water fluxes in m/s

     buoyancy_fl=0. ; bh_net=0. ; bw_net=0.
     WHERE ( zmask == 1 ) 
        bh_net(:,:)= zcoefq * qnet
        bw_net(:,:)= zcoefw * wnet
        buoyancy_fl(:,:) = ( bh_net - bw_net ) 
     END WHERE

     IF ( .NOT. lsho ) THEN
        ! Evap : 
        qlat(:,:)= getvar(cf_flxfil, cn_solhflup, 1, npiglo, npjglo, ktime=jt) *zmask(:,:)    ! W/m2 
        evap(:,:)= -1.* qlat(:,:) /Lv*86400. *zmask(:,:)                                    ! mm/days

        ! Wdmp
        wdmp(:,:)= getvar(cf_flxfil, cn_sowafldp, 1, npiglo, npjglo, ktime=jt)*86400.*zmask(:,:) ! mm/days

        ! Runoff  ! take care : not a model output (time_counter may disagree ... jmm
        runoff(:,:)= getvar(cf_rnfil, 'sorunoff', 1, npiglo, npjglo)*86400.*zmask(:,:)         ! mm/days

        ! fsalt = contribution of ice freezing and melting to salinity ( + = freezing, - = melting )Q
        wice(:,:) = getvar(cf_flxfil, cn_iowaflup, 1, npiglo, npjglo, ktime=jt )*86400.*zmask(:,:)          ! mm/days

        ! Precip:
        precip(:,:)= evap(:,:)-runoff(:,:)+wdmp(:,:)-wnet(:,:)+wice(:,:)                     ! mm/day

        ! Precip+runoff : (as a whole ) (interpolated on line)
        precip_runoff(:,:)= evap(:,:)+wdmp(:,:)-wnet(:,:)+wice(:,:)                          ! mm/day

        ! other heat fluxes
        qsb(:,:)= getvar(cf_flxfil, cn_sosbhfup,  1, npiglo, npjglo, ktime = jt )*zmask(:,:)    ! W/m2 
        qlw(:,:)= getvar(cf_flxfil, cn_solwfldo,  1, npiglo, npjglo, ktime = jt )*zmask(:,:)    ! W/m2 
        qsw(:,:)= getvar(cf_flxfil, cn_soshfldo,  1, npiglo, npjglo, ktime = jt )*zmask(:,:)    ! W/m2 


        b_qlat=0.  ; b_qlw=0.    ; b_qsw=0.   ; b_qsb=0.
        b_evap=0.  ; b_precip=0. ; b_wdmp=0.  ; b_runoff=0.

        WHERE (zsss /= 0 ) 
           b_qlat(:,:)= zcoefq * qlat
           b_qlw (:,:)= zcoefq * qlw
           b_qsw (:,:)= zcoefq * qsw
           b_qsb (:,:)= zcoefq * qsb

           b_evap(:,:)= zcoefw * evap
           b_precip(:,:)= -zcoefw * precip
           b_runoff(:,:)= -zcoefw * runoff
           b_wdmp(:,:)= zcoefw * wdmp

           !    buoyancy_fl(:,:) = zcoefq * qnet +zcoefw * wnet
        END WHERE
     ENDIF

     ! Write output file
     IF ( lsho ) THEN
        ierr = putvar(ncout, id_varout(1),buoyancy_fl, 1,npiglo, npjglo, ktime=jt )
     ELSE
        ierr = putvar(ncout, id_varout(1), evap,   1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(2), precip, 1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(3), runoff, 1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(4), wdmp,   1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(5), wnet,   1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(6), wice,   1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(7), precip_runoff,   1,npiglo, npjglo, ktime=jt )

        ierr = putvar(ncout, id_varout(8), qlat,   1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(9), qsb,    1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(10),qlw,    1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(11),qsw,    1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(12),qnet,   1, npiglo, npjglo, ktime=jt )

        ierr = putvar(ncout, id_varout(13),b_evap,  1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(14),b_precip,1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(15),b_runoff,1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(16),b_wdmp,  1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(17),bw_net,  1, npiglo, npjglo, ktime=jt )

        ierr = putvar(ncout, id_varout(18),b_qlat,  1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(19),b_qsb,   1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(20),b_qlw,   1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(21),b_qsw,   1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(22),bh_net,  1, npiglo, npjglo, ktime=jt )

        ierr = putvar(ncout, id_varout(23),buoyancy_fl, 1,npiglo, npjglo, ktime=jt )

        ierr = putvar(ncout, id_varout(24), zsss,  1, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(25), zsst,  1, npiglo, npjglo, ktime=jt )
     ENDIF
  END DO  ! time loop

  dtim = getvar1d(cf_tfil, cn_vtimec, npt     )
  ierr = putvar1d(ncout,   dtim,      npt, 'T')

  ierr=closeout(ncout)

CONTAINS
  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Set up all things required for the output file, create
    !!               the file and write the header part.
    !!
    !! ** Method  :  Use global module variables
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4) :: jv   ! dummy loop index
    !!----------------------------------------------------------------------
    ! prepare output variables
    ALLOCATE (zdep(1), dtim(npt) )
    zdep(1) = 0.
    ALLOCATE (ipk(np_varout), id_varout(np_varout), stypvar(np_varout) )
    ipk(:)  = 1  ! all variables ( output are 2D)
    stypvar%conline_operation = 'N/A'
    stypvar%caxis             = 'TYX'
    DO jv = 1, np_varout
       stypvar(jv)%ichunk = (/npiglo,MAX(1,npjglo/30), 1, 1 /)
    ENDDO

    IF ( lsho ) THEN
       ! total buoyancy flux
       stypvar(1)%cname= 'buoyancy_fl'
       stypvar(1)%cunits='1e-6 m2/s3'
       stypvar(1)%rmissing_value=0.
       stypvar(1)%valid_min= -100.
       stypvar(1)%valid_max= 100.
       stypvar(1)%clong_name='buoyancy flux'
       stypvar(1)%cshort_name='buoyancy_fl'
    ELSE
       ! 1--> 7 water fluxes                     ;   ! 8 --> 12    heat fluxes
       stypvar(1)%cname= 'evap'                  ;  stypvar(8)%cname= 'latent'      
       stypvar(2)%cname= 'precip'                ;  stypvar(9)%cname= 'sensible'      
       stypvar(3)%cname= 'runoff'                ;  stypvar(10)%cname= 'longwave'      
       stypvar(4)%cname= 'sssdmp'                ;  stypvar(11)%cname= 'solar'      
       stypvar(5)%cname= 'watnet'                ;  stypvar(12)%cname= 'heatnet'      
       stypvar(6)%cname= 'wice'            
       stypvar(7)%cname= 'precip_runoff'  

       stypvar(1:7)%cunits='mm/day'              ;  stypvar(8:12)%cunits='W/m2'      
       stypvar(1:7)%rmissing_value=0.            ;  stypvar(8:12)%rmissing_value=0.      
       stypvar(1:7)%valid_min= -100.             ;  stypvar(8:12)%valid_min= -500.      
       stypvar(1:7)%valid_max= 100.              ;  stypvar(8:12)%valid_max= 500.      
       stypvar(1)%clong_name='Evaporation'       ;  stypvar(8)%clong_name='Latent Heat flux'      
       stypvar(2)%clong_name='Precipitation'     ;  stypvar(9)%clong_name='Sensible Heat flux'       
       stypvar(3)%clong_name='Runoff'            ;  stypvar(10)%clong_name='Long Wave Heat flux'      
       stypvar(4)%clong_name='SSS damping'       ;  stypvar(11)%clong_name='Short Wave Heat flux'
       stypvar(5)%clong_name='Total water flux'  ;  stypvar(12)%clong_name='Net Heat Flux'      
       stypvar(6)%clong_name='Ice congelation and melting'  
       stypvar(7)%clong_name='Precip and runoff together' 

       stypvar(1)%cshort_name='evap'             ;  stypvar(8)%cshort_name='latent'      
       stypvar(2)%cshort_name='precip'           ;  stypvar(9)%cshort_name='sensible'      
       stypvar(3)%cshort_name='runoff'           ;  stypvar(10)%cshort_name='longwave'      
       stypvar(4)%cshort_name='sssdmp'           ;  stypvar(11)%cshort_name='solar'      
       stypvar(5)%cshort_name='watnet'           ;  stypvar(12)%cshort_name='heatnet'       
       stypvar(6)%cshort_name='wice'  
       stypvar(7)%cshort_name='precip_runoff' 

       ! 13--> 17  buoy water fluxes             ;   ! 18 --> 22    buoy heat fluxes
       stypvar(13)%cname= 'evap_b'               ;  stypvar(18)%cname= 'latent_b'
       stypvar(14)%cname= 'precip_b'             ;  stypvar(19)%cname= 'sensible_b'
       stypvar(15)%cname= 'runoff_b'             ;  stypvar(20)%cname= 'longwave_b'
       stypvar(16)%cname= 'sssdmp_b'             ;  stypvar(21)%cname= 'solar_b'
       stypvar(17)%cname= 'watnet_b'             ;  stypvar(22)%cname= 'heatnet_b'

       stypvar(13:17)%cunits='1e-6 m2/s3'      ;  stypvar(18:22)%cunits='1e-6 m2/s3'
       stypvar(13:17)%rmissing_value=0.          ;  stypvar(18:22)%rmissing_value=0.
       stypvar(13:17)%valid_min= -100.           ;  stypvar(18:22)%valid_min= -500.
       stypvar(13:17)%valid_max= 100.            ;  stypvar(18:22)%valid_max= 500.

       stypvar(13)%clong_name='buoy flx evap'    ;  stypvar(18)%clong_name='buoy Latent Heat flux'
       stypvar(14)%clong_name='buoy flx precip'  ;  stypvar(19)%clong_name='buoy Sensible Heat flux'
       stypvar(15)%clong_name='buoy flx runoff'  ;  stypvar(20)%clong_name='buoy Long Wave Heat flux'
       stypvar(16)%clong_name='buoy flx damping' ;  stypvar(21)%clong_name='buoy Short Wave Heat flux'
       stypvar(17)%clong_name='buoy haline flx'  ;  stypvar(22)%clong_name='buoy thermo Flux'

       stypvar(13)%cshort_name='evap_b'          ;  stypvar(18)%cshort_name='latent_b'
       stypvar(14)%cshort_name='precip_b'        ;  stypvar(19)%cshort_name='sensible_b'
       stypvar(15)%cshort_name='runoff_b'        ;  stypvar(20)%cshort_name='longwave_b'
       stypvar(16)%cshort_name='sssdmp_b'        ;  stypvar(21)%cshort_name='solar_b'
       stypvar(17)%cshort_name='watnet_b'        ;  stypvar(22)%cshort_name='heatnet_b'

       ! total buoyancy flux
       stypvar(23)%cname= 'buoyancy_fl'
       stypvar(23)%cunits='1e-6 m2/s3'
       stypvar(23)%rmissing_value=0.
       stypvar(23)%valid_min= -100.
       stypvar(23)%valid_max= 100.
       stypvar(23)%clong_name='buoyancy flux'
       stypvar(23)%cshort_name='buoyancy_fl'

       ! SSS                                         ; SST
       stypvar(24)%cname= 'sss'                      ;   stypvar(25)%cname= 'sst'
       stypvar(24)%cunits='PSU'                      ;   stypvar(25)%cunits='Celsius'
       stypvar(24)%rmissing_value=0.                 ;   stypvar(25)%rmissing_value=0.
       stypvar(24)%valid_min= 0.                     ;   stypvar(25)%valid_min= -2.
       stypvar(24)%valid_max= 45                     ;   stypvar(25)%valid_max= 45
       stypvar(24)%clong_name='Sea Surface Salinity' ;   stypvar(25)%clong_name='Sea Surface Temperature'
       stypvar(24)%cshort_name='sss  '               ;   stypvar(25)%cshort_name='sst'
    ENDIF

    ncout = create      (cf_out, cf_tfil, npiglo,    npjglo, 1,          ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar, np_varout, ipk,    id_varout , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_tfil, npiglo,    npjglo, 1,   pdep=zdep )
  END SUBROUTINE CreateOutput


END PROGRAM cdfbuoyflx
