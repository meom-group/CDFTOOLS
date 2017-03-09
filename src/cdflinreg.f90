PROGRAM cdflinreg
  !!======================================================================
  !!                     ***  PROGRAM  cdflinreg  ***
  !!=====================================================================
  !!  ** Purpose : Compute linear regression coef from a bunch of input 
  !!               cdf files given as argument.
  !!               Store the results on a 'similar' cdf file.
  !!
  !!  ** Method  : compute a and b such as yr = a . t + b
  !!              yr is the estimation of the field value, t is the time (in days ).
  !!              a= cov(y,t) / var(t)
  !!              b= moy(y) - a . moy(t)
  !!              R2 pearson value [0,1], giving the quality of the adjustment is also given
  !!              R2= a*a*var(t)/var(y)
  !!              cov(y,t)= moy(y*t) - moy(y)*moy(t)
  !!              var(t)  = moy(t*t) - moy(t)*moy(t)
  !!              var(y)  = moy(y*y) - moy(y)*moy(y)
  !!
  !! History : 2.1  : 01/2008  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!--------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class statistics
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER                    :: jptmax=1000               ! maximum number of time frame
  INTEGER(KIND=4)                               :: jk, jfil, jvar, jv, jt    ! dummy loop index
  INTEGER(KIND=4)                               :: ierr, ijvar               ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg        ! 
  INTEGER(KIND=4)                               :: npiglo, npjglo, npk, npt  ! size of the domain
  INTEGER(KIND=4)                               :: idep, idep_max     ! possible depth index, maximum
  INTEGER(KIND=4)                               :: nvars                     ! Number of variables in a file
  INTEGER(KIND=4)                               :: nfiles                    ! Number of files on the command line
  INTEGER(KIND=4)                               :: ntframe                   ! Cumul of time frame
  INTEGER(KIND=4)                               :: ncout                     ! ncid of output file
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE    :: ipki                      ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE    :: ipko                      ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE    :: id_varout                 ! id of output variable

  REAL(KIND=4)                                  :: zspval = -99999.          ! special value/ missing value
  REAL(KIND=4), DIMENSION(2)                    :: timean                    ! trick : timean(1) hold moy(t) (days)
  !                                                                          !         timean(2) hold moy(t2) (days)**2
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE       :: tim                       ! time counter

  REAL(KIND=8)                                  :: dt, dt2                   ! variables for cumulated time values
  REAL(KIND=8)                                  :: dtotal_time               ! cumulated time
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE     :: dy, dyy, dyt              ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE     :: dv2d                      ! Array to read a layer of data
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE     :: dmean, dmean2, dmean3     ! 
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE     :: dareg, dbreg, dpear       ! slope, origin ordinate, pearson coef

  CHARACTER(LEN=256)                            :: cf_in                     ! file names
  CHARACTER(LEN=256)                            :: cf_out='linreg.nc'        ! file names
  CHARACTER(LEN=256)                            :: cv_dep                    ! depth variable name
  CHARACTER(LEN=256)                            :: cldum                     ! depth variable name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_namesi                 ! array of var name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_nameso                 ! array of var22 name for output
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cf_lst                    ! list of input file names
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: clv_dep                   ! array of possible depth name (or 3rd dimension)

  TYPE (variable), DIMENSION(:), ALLOCATABLE    :: stypvari, stypvaro        ! data structure

  LOGICAL                                       :: lcaltmean                 ! flag for timemean computation
  LOGICAL                                       :: lnc4   = .FALSE.          ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdflinreg -l LIST-files  [-o OUT-file] [-nc4] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute the linear regression coefficients for a bunch of'
     PRINT *,'        input files. '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -l  LIST-files: A list of netcdf model file of same kind' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file ] : specify name of the output file instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ] : Use netcdf4 output with chunking and deflation level 1'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless option -o is used.'
     PRINT *,'         variables : for each input variables, there are 3 computed field'
     PRINT *,'            - slope coefficient'
     PRINT *,'            - barycenter '
     PRINT *,'            - Pearson Coefficient'
     STOP
  ENDIF

  ijarg = 1
  DO WHILE (ijarg <= narg )
     CALL getarg(ijarg, cldum  ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-l'   ) ; CALL GetFileList
        ! options
     CASE ( '-o'   ) ; CALL getarg (ijarg, cf_out ) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP
     END SELECT
  ENDDO

  !! Initialisation from 1st file (all file are assume to have the same geometry)
  cf_in=cf_lst(1)  
  IF ( chkfile(cf_in) ) STOP ! missing file

  npiglo = getdim (cf_in,cn_x                             )
  npjglo = getdim (cf_in,cn_y                             )

  ! looking for npk among various possible name
  idep_max=8
  ALLOCATE ( clv_dep(idep_max) )
  clv_dep(:) = (/cn_z,'z','sigma','nav_lev','levels','ncatice','icbcla','icbsect'/)
  idep=1  ; ierr=1000
  DO WHILE ( ierr /= 0 .AND. idep <= idep_max )
     npk  = getdim (cf_in, clv_dep(idep), cdtrue=cv_dep, kstatus=ierr)
     idep = idep + 1
  ENDDO

  IF ( ierr /= 0 ) THEN  ! none of the dim name was found
     PRINT *,' assume file with no depth'
     npk=0
  ENDIF

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk

  ALLOCATE( dy   (npiglo,npjglo), dyt   (npiglo,npjglo), dyy   (npiglo,npjglo), dv2d(npiglo,npjglo) )
  ALLOCATE( dmean(npiglo,npjglo), dmean2(npiglo,npjglo), dmean3(npiglo,npjglo)                      )
  ALLOCATE( dareg(npiglo,npjglo), dbreg (npiglo,npjglo) ,dpear (npiglo,npjglo)                      )
  ALLOCATE( tim  (jptmax) )

  nvars = getnvar(cf_in)
  PRINT *,' nvars =', nvars

  ALLOCATE (cv_namesi(nvars), cv_nameso(3*nvars) )
  ALLOCATE (stypvari (nvars), stypvaro (3*nvars) )
  ALLOCATE (ipki     (nvars), ipko     (3*nvars) )
  ALLOCATE (                  id_varout(3*nvars) )

  ! get list of variable names and collect attributes in stypvari (optional)
  cv_namesi(:) = getvarname(cf_in, nvars, stypvari)

  CALL CreateOutput

  lcaltmean=.TRUE. ; dt=0.d0 ; dt2=0.d0
  DO jvar = 1,nvars
     ijvar=(jvar-1)*3 + 1
     IF (cv_namesi(jvar) == cn_vlon2d .OR. &
          cv_namesi(jvar) == cn_vlat2d ) THEN
        ! skip these variable
     ELSE
        PRINT *,' Working with ', TRIM(cv_namesi(jvar)), ipki(jvar), jvar
        DO jk = 1, ipki(jvar)
           dy(:,:) = 0.d0 ; dyt(:,:) = 0.d0 ; dyy(:,:) =0.d0 ; dtotal_time = 0.;  ntframe=0
           DO jfil = 1, nfiles
              cf_in=cf_lst(jfil)
              IF ( jvar == 1 ) THEN
                 IF ( chkfile(cf_in) ) STOP ! missing file
              ENDIF
              npt = getdim (cf_in,cn_t)
              ntframe=ntframe+npt
              IF ( lcaltmean )  THEN
                 ! read time and convert seconds to years
                 tim(ntframe-npt+1:ntframe)=getvar1d(cf_in,cn_vtimec,npt)/86400.d0/365.
              END IF

              DO jt=1,npt
                 ! If forcing fields is without depth dimension
                 dv2d(:,:) = getvar(cf_in, cv_namesi(jvar), jk ,npiglo, npjglo, ktime=jt )
                 dy(:,:)   = dy(:,:)  + dv2d(:,:)
                 dyy(:,:)  = dyy(:,:) + dv2d(:,:)*dv2d(:,:)
                 dyt(:,:)  = dyt(:,:) + dv2d(:,:)*tim(ntframe-npt+jt)
              ENDDO
           END DO
           ! finish with level jk ; compute mean (assume zspval is 0 )
           dt          = SUM(tim(1:ntframe)                )
           dt2         = SUM(tim(1:ntframe)*tim(1:ntframe) )
           dmean(:,:)  = dy (:,:) / ntframe
           dmean2(:,:) = dyt(:,:) / ntframe
           dmean3(:,:) = dyy(:,:) / ntframe

           IF (lcaltmean )  THEN
              timean(1)= dt/ntframe
              timean(2)= dt2/ntframe
              ierr=putvar1d(ncout,timean,1,'T')
           END IF

           !compute dareg, dbreg, dpear
           WHERE (dmean /= 0 ) 
              dareg(:,:) = ( dmean2(:,:) - dmean(:,:) *timean(1) ) / ( timean(2) -timean(1)*timean(1) )
              dbreg(:,:) = dmean(:,:) - dareg(:,:)*timean(1) 
              dpear(:,:) = dareg(:,:)*dareg(:,:)*( timean(2) -timean(1)*timean(1))/( dmean3(:,:) -dmean(:,:)*dmean(:,:) )
              WHERE (dpear < 0 ) dpear=0 ; WHERE (dpear > 1 ) dpear=1
           ELSEWHERE
              dareg=zspval ; dbreg=zspval ; dpear=zspval
           ENDWHERE

           ierr = putvar(ncout, id_varout(ijvar  ), REAL(dareg), jk, npiglo, npjglo)
           ierr = putvar(ncout, id_varout(ijvar+1), REAL(dbreg), jk, npiglo, npjglo)
           ierr = putvar(ncout, id_varout(ijvar+2), REAL(dpear), jk, npiglo, npjglo)
           lcaltmean = .FALSE. ! tmean already computed
        END DO  ! loop to next level
     END IF
  END DO ! loop to next var in file

  ierr = closeout(ncout)

CONTAINS 

  SUBROUTINE GetFileList
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetFileList  ***
    !!
    !! ** Purpose :  Set up a file list given on the command line as 
    !!               blank separated list
    !!
    !! ** Method  :  Scan the command line until a '-' is found
    !!----------------------------------------------------------------------
    INTEGER (KIND=4)  :: ji
    INTEGER (KIND=4)  :: icur
    !!----------------------------------------------------------------------
    !!
    nfiles=0
    ! need to read a list of file ( number unknow ) 
    ! loop on argument till a '-' is found as first char
    icur=ijarg                          ! save current position of argument number
    DO ji = icur, narg                  ! scan arguments till - found
       CALL getarg ( ji, cldum )
       IF ( cldum(1:1) /= '-' ) THEN ; nfiles = nfiles+1
       ELSE                          ; EXIT
       ENDIF
    ENDDO
    ALLOCATE (cf_lst(nfiles) )
    DO ji = icur, icur + nfiles -1
       CALL getarg(ji, cf_lst( ji -icur +1 ) ) ; ijarg=ijarg+1
    END DO
  END SUBROUTINE GetFileList

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ! ipki gives the number of level or 0 if not a T[Z]YX  variable
    ipki(:) = getipk (cf_in, nvars, cdep=cv_dep)

    ! for each input variable, create 3 output variables
    DO jvar = 1, nvars
       ijvar=(jvar -1)*3 +1
       ! AREG
       ipko     (ijvar)                    = ipki(jvar)
       cv_nameso(ijvar)                    = TRIM(cv_namesi(jvar))//'_areg'
       stypvaro(ijvar)%ichunk              = (/npiglo,MAX(1,npjglo/30),1,1 /)
       stypvaro(ijvar)%cname               = TRIM(stypvari(jvar)%cname)//'_areg'   ! name
       stypvaro(ijvar)%cunits              = TRIM(stypvari(jvar)%cunits)//'/year'  ! unit
       stypvaro(ijvar)%rmissing_value      = zspval                                ! missing_value
       stypvaro(ijvar)%valid_min           = -100.                                 ! valid_min = zero
       stypvaro(ijvar)%valid_max           =  100.                                 ! valid_max *valid_max
       stypvaro(ijvar)%scale_factor        = 1.
       stypvaro(ijvar)%add_offset          = 0.
       stypvaro(ijvar)%savelog10           = 0.
       stypvaro(ijvar)%clong_name          = TRIM(stypvari(jvar)%clong_name)//'_linear_slope'
       stypvaro(ijvar)%cshort_name         = TRIM(stypvari(jvar)%cshort_name)//'_areg'     
       stypvaro(ijvar)%conline_operation   = TRIM(stypvari(jvar)%conline_operation) 
       stypvaro(ijvar)%caxis               = TRIM(stypvari(jvar)%caxis) 
       ! BREG
       ipko     (ijvar+1)                  = ipki(jvar)
       cv_nameso(ijvar+1)                  = TRIM(cv_namesi(jvar))//'_breg'
       stypvaro(ijvar+1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
       stypvaro(ijvar+1)%cname             = TRIM(stypvari(jvar)%cname)//'_breg'   ! name
       stypvaro(ijvar+1)%cunits            = TRIM(stypvari(jvar)%cunits)           ! unit
       stypvaro(ijvar+1)%rmissing_value    = zspval                                ! missing_value
       stypvaro(ijvar+1)%valid_min         = -100.                                 ! valid_min = zero
       stypvaro(ijvar+1)%valid_max         =  100.                                 ! valid_max *valid_max
       stypvaro(ijvar+1)%scale_factor      = 1.
       stypvaro(ijvar+1)%add_offset        = 0.
       stypvaro(ijvar+1)%savelog10         = 0.
       stypvaro(ijvar+1)%clong_name        = TRIM(stypvari(jvar)%clong_name)//'_b' 
       stypvaro(ijvar+1)%cshort_name       = TRIM(stypvari(jvar)%cshort_name)//'_breg' 
       stypvaro(ijvar+1)%conline_operation = TRIM(stypvari(jvar)%conline_operation)
       stypvaro(ijvar+1)%caxis             = TRIM(stypvari(jvar)%caxis)
       ! R2 pearson
       ipko     (ijvar+2)                  = ipki(jvar)
       cv_nameso(ijvar+2)                  = TRIM(cv_namesi(jvar))//'_r2'
       stypvaro(ijvar+2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
       stypvaro(ijvar+2)%cname             = TRIM(stypvari(jvar)%cname)//'_r2'    ! name
       stypvaro(ijvar+2)%cunits            = 'no unit'                            ! unit
       stypvaro(ijvar+2)%rmissing_value    = zspval                               ! missing_value
       stypvaro(ijvar+2)%valid_min         = 0.                                   ! valid_min = zero
       stypvaro(ijvar+2)%valid_max         = 1.                                   ! valid_max *valid_max
       stypvaro(ijvar+2)%scale_factor      = 1.
       stypvaro(ijvar+2)%add_offset        = 0.
       stypvaro(ijvar+2)%savelog10         = 0.
       stypvaro(ijvar+2)%clong_name        = TRIM(stypvari(jvar)%clong_name)//'_r2_Pearson'
       stypvaro(ijvar+2)%cshort_name       = TRIM(stypvari(jvar)%cshort_name)//'_r2'
       stypvaro(ijvar+2)%conline_operation = TRIM(stypvari(jvar)%conline_operation)
       stypvaro(ijvar+2)%caxis             = TRIM(stypvari(jvar)%caxis)
    END DO

    ! discart useless variables assigning 'none' name
    WHERE( ipki == 0 ) cv_namesi = 'none'
    WHERE( ipko == 0 ) cv_nameso = 'none'
    ! synchronize variables names
    stypvari(:)%cname = cv_namesi
    stypvaro(:)%cname = cv_nameso

    ! create output file taking the sizes in cf_in
    ncout  = create      (cf_out, cf_in,    npiglo,  npjglo, npk, cdep=cv_dep , ld_nc4=lnc4 )
    ierr   = createvar   (ncout,  stypvaro, 3*nvars, ipko,   id_varout        , ld_nc4=lnc4 )
    ierr   = putheadervar(ncout,  cf_in,    npiglo,  npjglo, npk, cdep=cv_dep )

  END SUBROUTINE CreateOutput

END PROGRAM cdflinreg
