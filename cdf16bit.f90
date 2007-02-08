PROGRAM cdf16bit
  !!-----------------------------------------------------------------------
  !!                 ***  PROGRAM cdf16bit  ***
  !!
  !!  **  Purpose: Transform the 32bit precision input file into a 16bit prec.
  !!               Uses constant scale_factor and add_offset
  !!               Store the results on a 'cdf16bit.nc' file similar to the input file.
  !!               Scale factor and offset are pre-defined for authorized cdf varname
  !!
  !!  **  Method:  read, transform and write 
  !!               Optional checks can be performed.
  !!
  !! history :
  !!     Original code :   J.M. Molines (Nov 2006 )
  !!-----------------------------------------------------------------------
  !!
  !! * Module used
  USE cdfio 

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk,jt,jvar, jv, jarg                         !: dummy loop index
  INTEGER   :: ierr                                         !: working integer
  INTEGER   :: narg, iargc                                  !: 
  INTEGER   :: npiglo,npjglo, npk                           !: size of the domain
  INTEGER   ::  nvars                                       !: Number of variables in a file
  INTEGER , DIMENSION(:), ALLOCATABLE :: id_var , &         !: arrays of var id's
       &                             ipk    , &             !: arrays of vertical level for each var
       &                             id_varout              !: cdf id for varout
  INTEGER(KIND=2),DIMENSION(:,:), ALLOCATABLE :: i2d       !: 16 bit 2D array fro conversion

  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE  :: v2d        !: Array to read a layer of data
  REAL(KIND=4), DIMENSION(1)                  :: tim        !: time of file
  REAL(KIND=4)                                :: sf, ao     !: scale_factor, add_offset
  REAL(KIND=4)                      :: zchkmax, zchkmin     !: scale_factor, add_offset checking values
  REAL(KIND=4),DIMENSION(:), ALLOCATABLE ::zmax, zmin       !: min and max of the field at level(jk)
  REAL(KIND=4)                      :: zzmax, zzmin         !: min and max of the full 3D field
  REAL(KIND=4)                      :: spval                !: missing value, fill_value, spval ...

  CHARACTER(LEN=80) :: cfile ,cfileout, cdum                !: file name
  CHARACTER(LEN=80) ,DIMENSION(:), ALLOCATABLE:: cvarname   !: array of var name
  
  TYPE (variable), DIMENSION(:), ALLOCATABLE :: typvar      !: Type variable is defined in cdfio.
                                                            !: It is used for attributes
  INTEGER    :: ncout
  INTEGER    :: istatus

  LOGICAL    :: l_chk=.false., l_verbose=.false.           !: logical flags to save line options
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: lmask            !: 2D logical land/sea mask (true on ocean)
  !!

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdf16bit 32bit_file [ -check ] [ -verbose]'
     PRINT *,'   If -check is used, control than the scale factors are adequate'
     PRINT *,'   If -verbose is used, give information level by level.'
     STOP
  ENDIF
  !!
  CALL getarg (1, cfile)

  ! Check for options and reflect options on logical flags
  IF (narg > 1 ) THEN 
    DO jarg=2,narg
     CALL getarg(jarg,cdum)
     IF ( cdum == '-check' ) THEN
        l_chk=.true.
     ELSE IF ( cdum == '-verbose' ) THEN
        l_chk=.true. ; l_verbose=.true.
     ELSE
        PRINT *,' OPTION ',TRIM(cdum),' not supported.' ; STOP
     ENDIF
    END DO
  ENDIF

  ! get domain dimension from input file
  npiglo= getdim (cfile,'x')
  npjglo= getdim (cfile,'y')
  npk   = getdim (cfile,'depth',kstatus=istatus)

  IF (istatus /= 0 ) THEN
     npk   = getdim (cfile,'z',kstatus=istatus)
     IF (istatus /= 0 ) STOP 'depth dimension name not suported'
  ENDIF

  ! Allocate memory 
  ALLOCATE( v2d(npiglo,npjglo), i2d(npiglo,npjglo), lmask(npiglo, npjglo) )
  ALLOCATE( zmin(npk) , zmax(npk) )

  ! Get the number of variables held in the file, allocate arrays
  nvars = getnvar(cfile)
  PRINT *,' nvars =', nvars

  ALLOCATE (cvarname(nvars) )
  ALLOCATE (typvar(nvars) )
  ALLOCATE (id_var(nvars),ipk(nvars),id_varout(nvars) )

  ! get list of variable names and collect attributes in typvar (optional)
  cvarname(:)=getvarname(cfile,nvars,typvar) 

  id_var(:)  = (/(jv, jv=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cfile,nvars)

  ! flags variable not to be treated by changing their name to none
  WHERE( ipk == 0 ) cvarname='none'
  typvar(:)%name=cvarname

  ! create output fileset
  cfileout='cdf16bit.nc'

  ! fills the scale_factor and add_offset attribute according to variable name
  ! if the variable is not documented, then, sf=1, ao=0. and no conversion
  ! is performed for this variable (It stays in REAL*4 )
  DO jvar=1,nvars
     IF (cvarname(jvar) /= 'none' ) CALL sf_ao(jvar)
  END DO

  ! create output file taking the sizes in cfile
  ncout =create(cfileout, cfile,npiglo,npjglo,npk)
  ! The variables are created as FLOAT or SHORT depending on the scale_factor AND add_offset attribute
  ierr= createvar(ncout , typvar,  nvars, ipk, id_varout )
  ierr= putheadervar(ncout , cfile, npiglo, npjglo, npk)

 ! Get time and write time
  tim=getvar1d(cfile,'time_counter',1) ; ierr=putvar1d(ncout,tim,1,'T')
  ! Loop on all variables of the file
  DO jvar = 1,nvars
     IF (cvarname(jvar) == 'none' ) THEN
        ! skip these variable  they are copied in ncout by putheader above
     ELSE
        sf=typvar(jvar)%scale_factor
        ao=typvar(jvar)%add_offset
        PRINT *,' Working with ', TRIM(cvarname(jvar)), ipk(jvar), sf, ao
        spval=typvar(jvar)%missing_value 
        DO jk = 1, ipk(jvar)
           v2d(:,:)= getvar(cfile, cvarname(jvar), jk ,npiglo, npjglo )
           IF ( sf == 1. .AND. ao == 0 ) THEN
            ! write FLOATS
            ierr = putvar(ncout, id_varout(jvar) ,v2d, jk, npiglo, npjglo)
            ! skip remaining of the do-loop, treat next level
            CYCLE
           ENDIF
           IF ( typvar(jvar)%savelog10 == 0 ) THEN
             ! take care of not converting 'special values'
             WHERE( v2d /= spval ) 
               i2d(:,:)=NINT((v2d(:,:)-ao)/sf)
             ELSEWHERE
               i2d(:,:)=0
             END WHERE
           ELSE  ! store log10  ao and sf refer to the log10 of the variable
             WHERE( v2d /= spval ) 
               i2d(:,:)=NINT((log10(v2d(:,:))-ao)/sf)
             ELSEWHERE
               i2d(:,:)=0
             END WHERE
           ENDIF
           CALL checkscaling
           ! write SHORT to the file
           ierr = putvar(ncout, id_varout(jvar) ,i2d, jk, npiglo, npjglo)
        END DO  ! loop to next level
     END IF
  END DO ! loop to next var in file

  istatus = closeout(ncout)

  CONTAINS
  SUBROUTINE sf_ao (kvar)
  !! --------------------------------------------------------------------------------
  !!                  ***   Subroutine sfao  ***
  !!
  !!   **  Purpose : set the scale_factor and add_offset for the variable kvar
  !!                 Also set the flag savelog10 when the log10 of the variable
  !!                 is stored, instead of the proper variable.
  !!
  !!   **  Method : recognize the variable name and set pre-defined values.
  !!                Give the min and max value for a given variable, and remap
  !!                it on -32000 +32000 (Integer*2. The max I2 is 32767 (2^15 -1)
  !!                Taking 32000 leaves  allows a slight overshoot ( 1%).
  !!                
  !!
  !!    history:
  !!         Original : J.M. MOLINES (Nov. 2006)
  !! --------------------------------------------------------------------------------
  !* Arguments
  INTEGER :: kvar  !: variable number

  ! * Local variables
  CHARACTER(LEN=80) :: clvarname
  REAL(KIND=4) :: zvmin, zvmax

  clvarname=cvarname(kvar)
  SELECT CASE (clvarname)
  ! gridT
    CASE ('votemper')             ! Potential temperature (Deg C)
      zvmin= -3.   ;  zvmax = 42.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('vosaline')             ! Salinity (PSU)
      zvmin= 0.   ;  zvmax = 42.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('sossheig')             ! Sea Surface Heigh (m)
      zvmin= -2.5   ;  zvmax = 2.5
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('somxl010')             ! Mixed layer depth (m)
      zvmin= 0.   ;  zvmax = 5000.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('sohefldo')             ! Total Heat flux Down (W/m2)
      zvmin= -1500.   ;  zvmax = 500.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('soshfldo')             ! Solar Heat flux Down (W/m2)
      zvmin= -0.1   ;  zvmax = 500.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('sowaflup')             ! Evaporation - Precipitation Up ( kg/m2/s)
      zvmin= -0.1   ;  zvmax = 0.1
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('sowafldp')             ! SSS damping term Up (kg/m2/s )
      zvmin= -10.  ;  zvmax = 15.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('iowaflup')             ! ???
      zvmin= -1.  ;  zvmax = 0.1
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('sowaflcd')             ! Concentration Dilution water flux (kg/m2/s)
      zvmin=-1.   ;  zvmax = 15.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('solhflup')             ! Latent Heat Flux Up (W/m2)
      zvmin=-800.   ;  zvmax = 150.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('solwfldo')             ! Long Wave radiation Heat flux Down (W/m2)
      zvmin=-200.   ;  zvmax = 50.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('sosbhfup')             ! Sensible Heat Flux Up (W/m2)
      zvmin=-800.   ;  zvmax = 100.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    ! gridU
    CASE ('vozocrtx')             ! Zonal Velocity U (m/s)
      zvmin= -3.0   ;  zvmax = 3.0
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('sozotaux')             !  Zonal Wind Stress (N/m2)
      zvmin= -1.5   ;  zvmax = 1.5
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    ! gridV
    CASE ('vomecrty')             ! Meridional Velocity V (m/s)
      zvmin= -3.0   ;  zvmax = 3.0
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('sometauy')             !  Meridional  Wind Stress (N/m2)
      zvmin= -1.5   ;  zvmax = 1.5
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    ! gridW
    CASE ('vovecrtz')             ! Vertical Velocity W (m/s)
      zvmin= -1.e-2   ;  zvmax = 1.e-2
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('votkeavt')             ! Vertical mixing coef log(avt) log(m2/s)
      zvmin= -8.   ;  zvmax = 2.0
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=1.

    !icemod 
    CASE ('isnowthi')             ! Snow Thickness (m)
      zvmin=0.   ;  zvmax = 5.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('iicethic')             ! Ice Thickness (m)
      zvmin=0.   ;  zvmax = 15.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('iiceprod')             ! Ice Production (m/kt) (step ice)
      zvmin=-0.05   ;  zvmax = 0.05
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('ileadfra')             ! Ice Lead Fraction (%) (In fact, ice concentration)
      zvmin= 0  ;  zvmax = 1.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('iicetemp')             ! Ice Temperature (Deg C )
      zvmin= -50.  ;  zvmax = 0.1
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('ioceflxb')             !Ocean Ice flux  (W/m2) 
      zvmin= -100.  ;  zvmax = 2500.0
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('iicevelu')             ! Zonal Ice Velocity (m/s) (at U point)
      zvmin= -2.  ;  zvmax = 2.0
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('iicevelv')             ! Meridional Ice Velocity (m/s) (at V point)
      zvmin= -2.  ;  zvmax = 2.0
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('isstempe')             ! Sea Surface Temperature (Deg C)
      zvmin= -3.  ;  zvmax = 42.0
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('isssalin')             ! Sea Surface Salinity (PSU)
      zvmin= 0.  ;  zvmax = 42.0
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('iocetflx')             ! Total Flux at Ocean Surface (W/m2)
      zvmin= -1500.  ;  zvmax = 500.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('iocesflx')             ! Solar Flux at Ocean Surface (W/m2)
      zvmin= 0.  ;  zvmax = 500.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('iocwnsfl')             ! Non Solar Flux at Ocean surface (W/m2)
      zvmin= -1500.  ;  zvmax = 200.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('iocesafl')             ! Salt Flux at Ocean Surface (kg/m2/kt)
      zvmin= -300.  ;  zvmax = 300.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('iocestru')             ! Zonal Ice Ocean Stress (N/m2)
      zvmin= -1.5  ;  zvmax = 1.5
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('iocestrv')             ! Meridional Ice Ocean Stress (N/m2)
      zvmin= -1.5  ;  zvmax = 1.5
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('iicesflx')             ! Solar FLux at ice/ocean Surface (W/m2)
      zvmin= -1.0  ;  zvmax = 500.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('iicenflx')             ! Non Solar FLux at ice/ocean Surface (W/m2)
      zvmin= -1500.  ;  zvmax = 300.
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE ('isnowpre')             ! Snow Precipitation (kg/day)
      zvmin= 0.  ;  zvmax = 0.0001
      typvar(kvar)%add_offset=(zvmin + zvmax) /2.
      typvar(kvar)%scale_factor= (zvmax-typvar(kvar)%add_offset)/32000.
      typvar(kvar)%savelog10=0.

    CASE DEFAULT
      PRINT *, TRIM(clvarname),'  is not recognized !'
      PRINT *, 'No conversion will be performed'
      typvar(kvar)%scale_factor=1.0
      typvar(kvar)%add_offset=0.
      typvar(kvar)%savelog10=0.
  END SELECT
  END SUBROUTINE sf_ao

  SUBROUTINE checkscaling
  !!---------------------------------------------------------------------------
  !!           ***  Subroutine checkscaling   ***
  !!
  !!    * Purpose : Check if the scale_factor and add_offset are ok for the current v2d field
  !!
  !!    * Method : - Needs -check and/or -verbose line option to be activated.
  !!               - Find the min and max of 3D field (called every level, and determine min/max)
  !!               - if -verbose option set, give details at every levels
  !!               - When last level is done, give the diagnostics in case of conflict
  !!
  !!  history:
  !!         Original : J.M. Molines ( Nov. 2006)
  !!---------------------------------------------------------------------------
  !! * All variables are global from the main program
  !!
     IF ( l_chk ) THEN  ! with this option, check if the max value of the field  can be
                        !mapped on I2 with actual values of Scale_factor and Add_offset
       lmask=.true. ; WHERE (v2d == spval ) lmask=.false.
       ! Works with log10 of v2d in case of savelog10=1
       IF (typvar(jvar)%savelog10 == 1 ) THEN
         WHERE( v2d /= 0. ) v2d=LOG10(v2d)
       ENDIF
       zmax(jk)=MAXVAL(v2d,lmask) ; zmin(jk)=MINVAL(v2d,lmask)
 
       ! Additional output if verbose mode
       IF ( l_verbose ) THEN
        zchkmax=(zmax(jk) - ao )/sf ; zchkmin = (zmin(jk) -ao ) /sf
        IF ( zchkmax >= 2**15 ) THEN
           PRINT *,TRIM(cvarname(jvar)), ' LEVEL ', jk ,' MIN = ',zmin(jk),' MAX = ', zmax(jk)
           PRINT *,' W A R N I N G ! : maximum too high for (sf,ao) pair.', TRIM(cfile)
           PRINT *,' Optimal value for this level AO = ', (zmin(jk) + zmax(jk) )/2.,' [ ', &
                 &  typvar(jvar)%add_offset,']'
           PRINT *,' Optimal value for this level SF = ', (zmax(jk) - (zmin(jk) + zmax(jk) )/2. )/32000., &
                 & ' [ ',typvar(jvar)%scale_factor,' ] '
         END IF

        IF ( zchkmin < -2**15 ) THEN
           PRINT *,TRIM(cvarname(jvar)), ' LEVEL ', jk ,' MIN = ',zmin(jk),' MAX = ', zmax(jk)
           PRINT *,' W A R N I N G ! : minimum too low for (sf,ao) pair.', TRIM(cfile)
           PRINT *,' Optimal value for this level AO = ', (zmin(jk) + zmax(jk) )/2.,' [ ', &
                 &   typvar(jvar)%add_offset,']'
           PRINT *,' Optimal value for this level SF = ', (zmax(jk) - (zmin(jk) + zmax(jk) )/2. )/32000., &
                 & ' [ ',typvar(jvar)%scale_factor,' ] '
        END IF
       END IF   ! verbose mode

       ! Print a warning if necessary after the last level of var has been processed
       IF ( jk == ipk(jvar) ) THEN
         zzmax=MAXVAL(zmax(1:ipk(jvar))) ; zzmin=MINVAL(zmin(1:ipk(jvar)))
         zchkmax=(zzmax - ao )/sf ; zchkmin = (zzmin -ao ) /sf
         IF ( zchkmax >= 2**15 ) THEN
           PRINT *,' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
           PRINT *,TRIM(cvarname(jvar)), ' MIN = ',zzmin,' MAX = ',zzmax,TRIM(cfile)
           PRINT *,' WARNING ! : maximum too high for (sf,ao) pair.'
           PRINT *,' Optimal value for this level AO = ', (zzmin + zzmax )/2.,' [ ', &
                 &   typvar(jvar)%add_offset,']'
           PRINT *,' Optimal value for this level SF = ', (zzmax - (zzmin + zzmax )/2. )/32000., &
                 & ' [ ',typvar(jvar)%scale_factor,' ] '
           PRINT *,' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
         END IF

         IF ( zchkmin < -2**15 ) THEN
           PRINT *,' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
           PRINT *,TRIM(cvarname(jvar)), ' MIN = ',zzmin,' MAX = ', zzmax,TRIM(cfile)
           PRINT *,' WARNING ! : minimum too low for (sf,ao) pair.'
           PRINT *,' Optimal value for  AO = ', (zzmin + zzmax )/2.,' [ ', & 
                 &   typvar(jvar)%add_offset,']'
           PRINT *,' Optimal value for  SF = ', (zzmax - (zzmin + zzmax )/2. )/32000., &
                 & ' [ ',typvar(jvar)%scale_factor,' ] '
           PRINT *,' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
         ENDIF
       END IF   ! last level
     END IF   ! check mode
  END SUBROUTINE checkscaling

END PROGRAM cdf16bit
