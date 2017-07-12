PROGRAM cdf16bit
  !!======================================================================
  !!                     ***  PROGRAM  cdf16bit  ***
  !!======================================================================
  !!  ** Purpose : Transform the 32bit precision input file into a 16bit prec.
  !!               Uses constant scale_factor and add_offset.
  !!
  !!  ** Method  : Store the results on a 'cdf16bit.nc' file similar to the input file.
  !!               Scale factor and offset are pre-defined for authorized cdf varname
  !!
  !! History : 2.1 !  11/2006    J.M. Molines   : Original code
  !!           3.0 !  12/2010    J.M. Molines   : Full Doctor form + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   sf_ao         : Scale Factor Add Offset
  !!   check_scaling : verify that the sf_ao do not produce saturation
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class data_transformation
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jt, jvar, jv   ! dummy loop index
  INTEGER(KIND=4)                               :: ierr               ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg ! 
  INTEGER(KIND=4)                               :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                               :: nvars              ! Number of variables in a file
  INTEGER(KIND=4)                               :: ncout              ! ncid of output file
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var             ! arrays of var id's
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                ! arrays of var id's
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout          ! arrays of var id's
  INTEGER(KIND=2), DIMENSION(:,:),  ALLOCATABLE :: i2d                ! 16 bit 2D array fro conversion

  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d                ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: zmax, zmin         ! min and max of the field at level(jk)
  REAL(KIND=4)                                  :: sf, ao             ! scale_factor, add_offset
  REAL(KIND=4)                                  :: zchkmax, zchkmin   ! scale_factor, add_offset checking values
  REAL(KIND=4)                                  :: zzmax, zzmin       ! min and max of the full 3D field
  REAL(KIND=4)                                  :: spval              ! missing value, fill_value, spval ...

  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dtim               ! time of file

  CHARACTER(LEN=256)                            :: cf_in              ! input file
  CHARACTER(LEN=256)                            :: cf_out='cdf16bit.nc' ! outputfile
  CHARACTER(LEN=256)                            :: cldum              ! dummy string
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names           ! array of var name

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar            ! Type variable is defined in cdfio.

  LOGICAL                                       :: l_chk    = .FALSE. ! logical flags to save line options
  LOGICAL                                       :: l_verbose= .FALSE. ! logical flags to save line options
  LOGICAL                                       :: lnc4     = .FALSE. ! Use nc4 with chunking and deflation
  LOGICAL, DIMENSION(:,:),          ALLOCATABLE :: lmask              ! 2D logical land/sea mask (true on ocean)
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdf16bit -f 32BIT-file [-check] [-verbose] [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Convert input 32 bit precision file into 16 bit precision file using' 
     PRINT *,'       add_offset and scale_factor. '
     PRINT *,'       Note that predifined values for these two parameters are defined '
     PRINT *,'       according to the variable name. If variable name is not supported,'
     PRINT *,'       no conversion is performed.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f 32BIT-file : input 32 bit file to be converted' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-check ]   : control than the scale factors are adequate.'
     PRINT *,'       [-verbose ] : give information level by level.' 
     PRINT *,'       [-o OUT-file] : Specify output file name instead of ', TRIM(cf_out)
     PRINT *,'       [-nc4 ] : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : same names than in input file'
     STOP 
  ENDIF
  !!
  ijarg = 1
  DO WHILE ( ijarg <= narg)  
     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE (  cldum  )  
     CASE ( '-f'      ) ; CALL getarg (ijarg, cf_in ) ; ijarg = ijarg + 1
        ! options
     CASE ( '-check'  ) ; l_chk     =.TRUE.
     CASE ( '-verbose') ; l_chk     =.TRUE. 
        ;                 l_verbose =.TRUE.
     CASE ( '-o'      ) ; CALL getarg (ijarg, cf_out) ; ijarg = ijarg + 1
     CASE ( '-nc4'    ) ; lnc4      =.TRUE.
     CASE DEFAULT       ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  IF ( chkfile(cf_in)  ) STOP 99 ! missing file

  ! get domain dimension from input file
  npiglo = getdim (cf_in, cn_x)
  npjglo = getdim (cf_in, cn_y)
  npk    = getdim (cf_in, cn_z, kstatus=ierr)

  IF (ierr /= 0 ) THEN
     npk   = getdim (cf_in,'z',kstatus=ierr)
     IF (ierr /= 0 ) THEN ;  PRINT *,' ERROR : depth dimension name not suported' ; STOP 99
     ENDIF
  ENDIF
  npt    = getdim (cf_in, cn_t)

  ! Allocate memory 
  ALLOCATE( v2d(npiglo,npjglo), i2d(npiglo,npjglo), lmask(npiglo, npjglo) )
  ALLOCATE( zmin(npk) , zmax(npk) , dtim(npt))

  ! Get the number of variables held in the file, allocate arrays
  nvars = getnvar(cf_in)
  PRINT *,' nvars =', nvars

  ALLOCATE (cv_names(nvars) )
  ALLOCATE (stypvar(nvars) )
  ALLOCATE (id_var(nvars),ipk(nvars),id_varout(nvars) )

  ! get list of variable names and collect attributes in stypvar (optional)
  cv_names(:)=getvarname(cf_in,nvars,stypvar) 

  CALL CreateOutput

  ! Loop on all variables of the file
  DO jvar = 1,nvars
     IF (cv_names(jvar) == 'none' ) THEN
        ! skip these variable  they are copied in ncout by putheader above
     ELSE
        sf=stypvar(jvar)%scale_factor
        ao=stypvar(jvar)%add_offset
        PRINT *,' Working with ', TRIM(cv_names(jvar)), ipk(jvar), sf, ao
        spval=stypvar(jvar)%rmissing_value 
        DO jt = 1, npt
           DO jk = 1, ipk(jvar)
              v2d(:,:)= getvar(cf_in, cv_names(jvar), jk ,npiglo, npjglo, ktime=jt )
              IF ( sf == 1. .AND. ao == 0 ) THEN
                 ! write FLOATS
                 IF ( stypvar(jvar)%savelog10 == 1 ) THEN
                    WHERE ( v2d /= spval )
                       v2d(:,:)= LOG10(v2d)
                    ELSEWHERE
                       v2d = 0.
                    END WHERE
                 ENDIF
                 ierr = putvar(ncout, id_varout(jvar) ,v2d, jk, npiglo, npjglo, ktime=jt)
                 ! skip remaining of the do-loop, treat next level
                 CYCLE
              ENDIF
              IF ( stypvar(jvar)%savelog10 == 0 ) THEN
                 ! take care of not converting 'special values'
                 WHERE( v2d /= spval ) 
                    i2d(:,:)=NINT((v2d(:,:)-ao)/sf)
                 ELSEWHERE
                    i2d(:,:)=0
                 END WHERE
              ELSE  ! store log10  ao and sf refer to the log10 of the variable
                 WHERE( v2d /= spval ) 
                    i2d(:,:)=NINT((LOG10(v2d(:,:))-ao)/sf)
                 ELSEWHERE
                    i2d(:,:)=0
                 END WHERE
              ENDIF
              CALL checkscaling
              ! write SHORT to the file
              ierr = putvar(ncout, id_varout(jvar) ,i2d, jk, npiglo, npjglo, ktime=jt)
           END DO  ! loop to next level
        END DO  ! next time loop
     END IF
  END DO ! loop to next var in file

  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE sf_ao (kvar)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE sf_ao  ***
    !!
    !! ** Purpose :  Set the scale_factor and add_offset for the variable kvar
    !!               Also set the flag savelog10 when the log10 of the variable
    !!               is stored, instead of the proper variable. 
    !!
    !! ** Method  : Recognize the variable name and set pre-defined values.
    !!              Give the min and max value for a given variable, and remap
    !!              it on -32000 +32000 (Integer*2. The max I2 is 32767 (2^15 -1)
    !!              Taking 32000 leaves  allows a slight overshoot ( 1%) 
    !!
    !! ** Comments : With select case (which gives a much more readable code,
    !!               the CASE statement requires a constant matching pattern,
    !!               thus avoiding the use of dynamically adjusted names as
    !!               defined in modcdfnames ...
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) :: kvar  ! variable number

    CHARACTER(LEN=256) :: clvarname
    REAL(KIND=4)       :: zvmin, zvmax
    !!----------------------------------------------------------------------

    clvarname=cv_names(kvar)
    SELECT CASE (clvarname)
       ! gridT
    CASE ('votemper')             ! Potential temperature (Deg C)
       zvmin= -3.   ;  zvmax = 42.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('vosaline')             ! Salinity (PSU)
       zvmin= 0.   ;  zvmax = 42.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('sossheig')             ! Sea Surface Heigh (m)
       zvmin= -2.5   ;  zvmax = 2.5
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('somxl010')             ! Mixed layer depth (m)
       zvmin= 0.   ;  zvmax = 5000.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('sohefldo')             ! Total Heat flux Down (W/m2)
       zvmin= -1500.   ;  zvmax = 500.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('soshfldo')             ! Solar Heat flux Down (W/m2)
       zvmin= -0.1   ;  zvmax = 500.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('sowaflup')             ! Evaporation - Precipitation Up ( kg/m2/s)
       zvmin= -0.1   ;  zvmax = 0.1
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('sowafldp')             ! SSS damping term Up (kg/m2/s )
       zvmin= -10.  ;  zvmax = 15.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('iowaflup')             ! ???
       zvmin= -1.  ;  zvmax = 0.1
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('sowaflcd')             ! Concentration Dilution water flux (kg/m2/s)
       zvmin=-1.   ;  zvmax = 15.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('solhflup')             ! Latent Heat Flux Up (W/m2)
       zvmin=-800.   ;  zvmax = 150.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('solwfldo')             ! Long Wave radiation Heat flux Down (W/m2)
       zvmin=-200.   ;  zvmax = 50.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('sosbhfup')             ! Sensible Heat Flux Up (W/m2)
       zvmin=-800.   ;  zvmax = 100.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

       ! gridU
    CASE ('vozocrtx')             ! Zonal Velocity U (m/s)
       zvmin= -3.0   ;  zvmax = 3.0
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('sozotaux')             !  Zonal Wind Stress (N/m2)
       zvmin= -1.5   ;  zvmax = 1.5
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

       ! gridV
    CASE ('vomecrty')             ! Meridional Velocity V (m/s)
       zvmin= -3.0   ;  zvmax = 3.0
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('sometauy')             !  Meridional  Wind Stress (N/m2)
       zvmin= -1.5   ;  zvmax = 1.5
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

       ! gridW
    CASE ('vovecrtz')             ! Vertical Velocity W (m/s)
       zvmin= -1.e-2   ;  zvmax = 1.e-2
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('votkeavt')             ! Vertical mixing coef log(avt) log(m2/s)
       zvmin= -8.   ;  zvmax = 2.0
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=1.

       !icemod 
    CASE ('isnowthi')             ! Snow Thickness (m)
       zvmin=0.   ;  zvmax = 5.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('iicethic')             ! Ice Thickness (m)
       zvmin=0.   ;  zvmax = 15.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('iiceprod')             ! Ice Production (m/kt) (step ice)
       zvmin=-0.05   ;  zvmax = 0.05
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('ileadfra')             ! Ice Lead Fraction (%) (In fact, ice concentration)
       zvmin= 0  ;  zvmax = 1.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('iicetemp')             ! Ice Temperature (Deg C )
       zvmin= -50.  ;  zvmax = 0.1
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('ioceflxb')             !Ocean Ice flux  (W/m2) 
       zvmin= -100.  ;  zvmax = 2500.0
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('iicevelu')             ! Zonal Ice Velocity (m/s) (at U point)
       zvmin= -2.  ;  zvmax = 2.0
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('iicevelv')             ! Meridional Ice Velocity (m/s) (at V point)
       zvmin= -2.  ;  zvmax = 2.0
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('isstempe')             ! Sea Surface Temperature (Deg C)
       zvmin= -3.  ;  zvmax = 42.0
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('isssalin')             ! Sea Surface Salinity (PSU)
       zvmin= 0.  ;  zvmax = 42.0
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('iocetflx')             ! Total Flux at Ocean Surface (W/m2)
       zvmin= -1500.  ;  zvmax = 500.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('iocesflx')             ! Solar Flux at Ocean Surface (W/m2)
       zvmin= 0.  ;  zvmax = 500.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('iocwnsfl')             ! Non Solar Flux at Ocean surface (W/m2)
       zvmin= -1500.  ;  zvmax = 200.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('iocesafl')             ! Salt Flux at Ocean Surface (kg/m2/kt)
       zvmin= -300.  ;  zvmax = 300.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('iocestru')             ! Zonal Ice Ocean Stress (N/m2)
       zvmin= -1.5  ;  zvmax = 1.5
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('iocestrv')             ! Meridional Ice Ocean Stress (N/m2)
       zvmin= -1.5  ;  zvmax = 1.5
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('iicesflx')             ! Solar FLux at ice/ocean Surface (W/m2)
       zvmin= -1.0  ;  zvmax = 500.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('iicenflx')             ! Non Solar FLux at ice/ocean Surface (W/m2)
       zvmin= -1500.  ;  zvmax = 300.
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

    CASE ('isnowpre')             ! Snow Precipitation (kg/day)
       zvmin= 0.  ;  zvmax = 0.0001
       stypvar(kvar)%add_offset=(zvmin + zvmax) /2.
       stypvar(kvar)%scale_factor= (zvmax-stypvar(kvar)%add_offset)/32000.
       stypvar(kvar)%savelog10=0.

       ! TRC
    CASE ('cfc11')               ! Concentration tracer 1
       zvmin= 0.  ;  zvmax = 0.0001
       stypvar(kvar)%add_offset=0.
       stypvar(kvar)%scale_factor= 1.
       stypvar(kvar)%savelog10=1.

    CASE ('bombc14')               ! Concentration tracer 1
       zvmin= 0.  ;  zvmax = 0.0001
       stypvar(kvar)%add_offset=0.
       stypvar(kvar)%scale_factor= 1.
       stypvar(kvar)%savelog10=1.

    CASE DEFAULT
       PRINT *, TRIM(clvarname),'  is not recognized !'
       PRINT *, 'No conversion will be performed'
       stypvar(kvar)%scale_factor=1.0
       stypvar(kvar)%add_offset=0.
       stypvar(kvar)%savelog10=0.
    END SELECT

  END SUBROUTINE sf_ao

  SUBROUTINE checkscaling()
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE checkscaling  ***
    !!
    !! ** Purpose :  Check if the scale_factor and add_offset are ok for 
    !!               the current v2d field
    !!
    !! ** Method  :  - Needs -check and/or -verbose line option to be activated.
    !!               - Find the min and max of 3D field (called every level, and determine min/max)
    !!               - if -verbose option set, give details at every levels
    !!               - When last level is done, give the diagnostics in case of conflict  
    !!
    !!----------------------------------------------------------------------
    IF ( l_chk ) THEN  ! with this option, check if the max value of the field  can be
       !mapped on I2 with actual values of Scale_factor and Add_offset
       lmask=.TRUE. ; WHERE (v2d == spval ) lmask=.FALSE.
       ! Works with log10 of v2d in case of savelog10=1
       IF (stypvar(jvar)%savelog10 == 1 ) THEN
          WHERE( v2d /= 0. ) v2d=LOG10(v2d)
       ENDIF
       zmax(jk)=MAXVAL(v2d,lmask) ; zmin(jk)=MINVAL(v2d,lmask)

       ! Additional output if verbose mode
       IF ( l_verbose ) THEN
          zchkmax=(zmax(jk) - ao )/sf ; zchkmin = (zmin(jk) -ao ) /sf
          IF ( zchkmax >= 2**15 ) THEN
             PRINT *,TRIM(cv_names(jvar)), ' LEVEL ', jk ,' MIN = ',zmin(jk),' MAX = ', zmax(jk)
             PRINT *,' W A R N I N G ! : maximum too high for (sf,ao) pair.', TRIM(cf_in)
             PRINT *,' Optimal value for this level AO = ', (zmin(jk) + zmax(jk) )/2.,' [ ', &
                  &  stypvar(jvar)%add_offset,']'
             PRINT *,' Optimal value for this level SF = ', (zmax(jk) - (zmin(jk) + zmax(jk) )/2. )/32000., &
                  & ' [ ',stypvar(jvar)%scale_factor,' ] '
          END IF

          IF ( zchkmin < -2**15 ) THEN
             PRINT *,TRIM(cv_names(jvar)), ' LEVEL ', jk ,' MIN = ',zmin(jk),' MAX = ', zmax(jk)
             PRINT *,' W A R N I N G ! : minimum too low for (sf,ao) pair.', TRIM(cf_in)
             PRINT *,' Optimal value for this level AO = ', (zmin(jk) + zmax(jk) )/2.,' [ ', &
                  &   stypvar(jvar)%add_offset,']'
             PRINT *,' Optimal value for this level SF = ', (zmax(jk) - (zmin(jk) + zmax(jk) )/2. )/32000., &
                  & ' [ ',stypvar(jvar)%scale_factor,' ] '
          END IF
       END IF   ! verbose mode

       ! Print a warning if necessary after the last level of var has been processed
       IF ( jk == ipk(jvar) ) THEN
          zzmax=MAXVAL(zmax(1:ipk(jvar))) ; zzmin=MINVAL(zmin(1:ipk(jvar)))
          zchkmax=(zzmax - ao )/sf ; zchkmin = (zzmin -ao ) /sf
          IF ( zchkmax >= 2**15 ) THEN
             PRINT *,' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
             PRINT *,TRIM(cv_names(jvar)), ' MIN = ',zzmin,' MAX = ',zzmax,TRIM(cf_in)
             PRINT *,' WARNING ! : maximum too high for (sf,ao) pair.'
             PRINT *,' Optimal value for this level AO = ', (zzmin + zzmax )/2.,' [ ', &
                  &   stypvar(jvar)%add_offset,']'
             PRINT *,' Optimal value for this level SF = ', (zzmax - (zzmin + zzmax )/2. )/32000., &
                  & ' [ ',stypvar(jvar)%scale_factor,' ] '
             PRINT *,' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
          END IF

          IF ( zchkmin < -2**15 ) THEN
             PRINT *,' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
             PRINT *,TRIM(cv_names(jvar)), ' MIN = ',zzmin,' MAX = ', zzmax,TRIM(cf_in)
             PRINT *,' WARNING ! : minimum too low for (sf,ao) pair.'
             PRINT *,' Optimal value for  AO = ', (zzmin + zzmax )/2.,' [ ', & 
                  &   stypvar(jvar)%add_offset,']'
             PRINT *,' Optimal value for  SF = ', (zzmax - (zzmin + zzmax )/2. )/32000., &
                  & ' [ ',stypvar(jvar)%scale_factor,' ] '
             PRINT *,' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
          ENDIF
       END IF   ! last level
    END IF   ! check mode
  END SUBROUTINE checkscaling

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    id_var(:)  = (/(jv, jv=1,nvars)/)
    ! ipk gives the number of level or 0 if not a T[Z]YX  variable
    ipk(:)     = getipk (cf_in,nvars)

    ! flags variable not to be treated by changing their name to none
    WHERE( ipk == 0 ) cv_names='none'
    stypvar(:)%cname=cv_names

    ! create output fileset
    ! fills the scale_factor and add_offset attribute according to variable name
    ! if the variable is not documented, then, sf=1, ao=0. and no conversion
    ! is performed for this variable (It stays in REAL*4 )
    DO jvar=1,nvars
       IF (cv_names(jvar) /= 'none' ) CALL sf_ao(jvar)
       stypvar(jvar)%ichunk = (/npiglo,MAX(1,npjglo/30),1,1 /)
    END DO

    ! create output file taking the sizes in cf_in
    ncout =create(cf_out, cf_in,npiglo,npjglo,npk           , ld_nc4=lnc4 )
    ! The variables are created as FLOAT or SHORT depending on the scale_factor AND add_offset attribute
    ierr= createvar(ncout , stypvar,  nvars, ipk, id_varout , ld_nc4=lnc4 )
    ierr= putheadervar(ncout , cf_in, npiglo, npjglo, npk)

    ! Get time and write time
    dtim=getvar1d(cf_in,cn_vtimec,npt) 
    ierr=putvar1d(ncout,dtim,npt,'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdf16bit
