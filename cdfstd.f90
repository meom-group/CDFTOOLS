PROGRAM cdfstd
  !!======================================================================
  !!                     ***  PROGRAM  cdfstd  ***
  !!=====================================================================
  !!  ** Purpose : Compute Standard deviation values for all the 
  !!               variables in a bunch of cdf files given as argument
  !!               Store the results on a 'similar' cdf file.
  !!
  !!  ** Method  : Compute mean, mean squared, then the variance and
  !!               the standard deviation
  !!
  !! History : 2.1  : 04/2006  : F. Castruccio : Original code (from cdfmoy)
  !!           3.0  : 01/2011  : J.M. Molines  : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jfil, jt        ! dummy loop index
  INTEGER(KIND=4)                               :: jvar, jv            ! dummy loop index
  INTEGER(KIND=4)                               :: narg, iargc, ijarg  ! browse line
  INTEGER(KIND=4)                               :: npiglo, npjglo      ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt            ! size of the domain
  INTEGER(KIND=4)                               :: nvars               ! number of variables in a file
  INTEGER(KIND=4)                               :: ntframe             ! cumul of time frame
  INTEGER(KIND=4)                               :: ncout               ! ncid of stdev file output
  INTEGER(KIND=4)                               :: ncou2               ! ncid of mean file output (optional)
  INTEGER(KIND=4)                               :: ierr                ! error status
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var              ! varid's of input variables
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk, id_varout      ! levels and varid's of output vars
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varoutm          ! varid's of mean var output (optional)

  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d                 ! 2d data array
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: tim                 ! tim counter
  REAL(KIND=4), DIMENSION(1)                    :: timean              ! mean time

  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dtab, dtab2         ! cumulated values and squared values
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dstd                ! standard deviation
  REAL(KIND=8)                                  :: dtotal_time         ! cumulated time

  CHARACTER(LEN=256)                            :: cf_in               ! input file
  CHARACTER(LEN=256)                            :: cf_out='cdfstd.nc'  ! std dev output file
  CHARACTER(LEN=256)                            :: cf_moy='cdfmoy.nc'  ! mean output file (optional)
  CHARACTER(LEN=256)                            :: cv_dep              ! depth variable name
  CHARACTER(LEN=256)                            :: cldum               ! dummy string
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_namesi           ! array of var name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_nameso           ! array of var name for output

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvari            ! attributes of input variables
  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvaro            ! attributes of output variables

  LOGICAL                                       :: lcaltmean           ! time mean computation flag
  LOGICAL                                       :: lsave=.false.       ! mean value save flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfstd list_of files [-save]' 
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the standard deviation of the variables belonging to a set of' 
     PRINT *,'       files given as arguments.  This computation is direct and does not '
     PRINT *,'       required a pre-processing with any of the cdfmoy tools.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       List on netcdf files of the same type, forming a time-series' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [ -save ] : Save the mean value of the field, in addition to the '
     PRINT *,'                   std deviation' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       - netcdf file : ', TRIM(cf_out) 
     PRINT *,'           variables :  IN-var_std, same units than input variables.'
     PRINT *,'       - netcdf file : ', TRIM(cf_moy),' in case of -save option.' 
     PRINT *,'           variables :  IN-var, same units than input variables.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'        cdfmoy, cdfrmsssh, cdfstdevw'
     STOP
  ENDIF

  ! look for -save option and one of the file name 
  ijarg = 1
  DO WHILE ( ijarg <= narg ) 
    CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
    SELECT CASE ( cldum )
    CASE ( '-save' ) 
      lsave = .true.
    CASE DEFAULT 
       CALL getarg (ijarg, cf_in) ; ijarg = ijarg + 1
       EXIT  ! got the first file
    END SELECT
  END DO

  IF ( chkfile(cf_in) ) STOP ! missing file

  npiglo = getdim (cf_in, cn_x)
  npjglo = getdim (cf_in, cn_y)
  npk    = getdim (cf_in, cn_z, cdtrue=cv_dep, kstatus=ierr)

  IF (ierr /= 0 ) THEN
     npk   = getdim (cf_in,'z',kstatus=ierr)
     IF (ierr /= 0 ) THEN
        npk   = getdim (cf_in,'sigma',cdtrue=cv_dep,kstatus=ierr)
       IF (ierr /= 0 ) THEN
         PRINT *,' assume file with no depth'
         npk=0
       ENDIF
     ENDIF
  ENDIF

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk

  ALLOCATE( dtab(npiglo,npjglo), dtab2(npiglo,npjglo), v2d(npiglo,npjglo) )
  ALLOCATE( dstd(npiglo,npjglo)                                           )

  nvars = getnvar(cf_in)
  PRINT *,' nvars =', nvars

  ALLOCATE (cv_namesi(nvars), cv_nameso(nvars)          )
  ALLOCATE (stypvari(nvars), stypvaro(nvars)            )
  ALLOCATE (id_var(nvars), ipk(nvars), id_varout(nvars) )
  IF ( lsave ) ALLOCATE (id_varoutm(nvars)              )

  cv_namesi(:) = getvarname(cf_in, nvars, stypvari)

  id_var(:)  = (/(jv, jv=1,nvars)/)
  ipk(:)     = getipk(cf_in, nvars, cdep=cv_dep)
  DO jvar = 1, nvars 
     cv_nameso(jvar) = TRIM(cv_namesi(jvar))//'_std' 
  ENDDO

  WHERE( ipk == 0 ) cv_nameso='none'

  DO jvar = 1, nvars
     stypvaro(jvar)             = stypvari(jvar)
     stypvaro(jvar)%cname       = cv_nameso(jvar)
     stypvaro(jvar)%clong_name  = 'Std Deviation of '//TRIM(cv_namesi(jvar))
     stypvaro(jvar)%cshort_name = cv_nameso(jvar)
  END DO

  ! create output fileset
  ncout = create      (cf_out, cf_in,    npiglo, npjglo, npk, cdep=cv_dep )
  ierr  = createvar   (ncout,  stypvaro, nvars,  ipk,    id_varout        )
  ierr  = putheadervar(ncout,  cf_in,    npiglo, npjglo, npk, cdep=cv_dep )

  IF ( lsave )  THEN
    ! create output fileset for mean values
    ncou2 = create      (cf_moy, cf_in,    npiglo, npjglo, npk, cdep=cv_dep )
    ierr  = createvar   (ncou2,  stypvari, nvars,  ipk,    id_varoutm       )
    ierr  = putheadervar(ncou2,  cf_in,    npiglo, npjglo, npk, cdep=cv_dep )
  ENDIF

  lcaltmean=.TRUE.
  DO jvar = 1,nvars
     IF ( cv_namesi(jvar) == cn_vlon2d .OR. &
          cv_namesi(jvar) == cn_vlat2d .OR. &
          cv_nameso(jvar) == 'none'  ) THEN
        ! skip these variable
     ELSE
        PRINT *,' Working with ', TRIM(cv_namesi(jvar)), ipk(jvar)

        DO jk = 1, ipk(jvar)
           PRINT *,'level ',jk

           dtab(:,:) = 0.d0; dtab2(:,:) = 0.d0; dtotal_time = 0.d0 
           ntframe = 0
           DO jfil = 1, narg
              CALL getarg (jfil, cf_in)
              IF ( chkfile(cf_in) ) STOP ! missing file

              IF ( lcaltmean )  THEN
                 npt = getdim (cf_in, cn_t)
                 ALLOCATE (tim(npt) )
                 tim(:)      = getvar1d(cf_in, cn_vtimec, npt)
                 dtotal_time = dtotal_time + SUM(DBLE(tim))
                 DEALLOCATE ( tim ) 
              END IF

              DO jt=1,npt
                 ntframe = ntframe + 1
                 v2d(  :,:) = getvar(cf_in, cv_namesi(jvar), jk, npiglo, npjglo, ktime=jt)
                 dtab( :,:) = dtab( :,:) + v2d(:,:)*1.d0
                 dtab2(:,:) = dtab2(:,:) + v2d(:,:)*v2d(:,:)*1.d0
              END DO
           END DO

           ! finish with level jk ; compute mean (assume spval is 0 )
           dtab( :,:) = dtab( :,:) / ntframe
           dtab2(:,:) = dtab2(:,:) / ntframe

           WHERE ( dtab2 - dtab*dtab >= 0 ) 
             dstd = SQRT(dtab2 - dtab*dtab)
           ELSE WHERE
             dstd = 0.d0
           END WHERE

           ! store variable on output file
                        ierr = putvar(ncout, id_varout(jvar),  REAL(dstd), jk, npiglo, npjglo, kwght=ntframe)
           IF ( lsave ) ierr = putvar(ncou2, id_varoutm(jvar), REAL(dtab), jk, npiglo, npjglo, kwght=ntframe)

           IF ( lcaltmean )  THEN
              timean(1) = dtotal_time / ntframe
                           ierr = putvar1d(ncout, timean, 1, 'T')
              IF ( lsave ) ierr = putvar1d(ncou2, timean, 1, 'T')
              lcaltmean = .FALSE. ! tmean already computed 
           END IF
        END DO  ! loop to next level
     END IF

  END DO ! loop to next var in file

  ierr = closeout(ncout)
  ierr = closeout(ncou2)

END PROGRAM cdfstd
