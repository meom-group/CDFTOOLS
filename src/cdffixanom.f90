PROGRAM cdffixanom
  !!======================================================================
  !!                     ***  PROGRAM  cdffixanom  ***
  !!=====================================================================
  !!  ** Purpose : Tool dedicated to the local correction of initial 
  !!               conditions files, where deepest part are some time 
  !!               Spurious. Below a given level, the values of the variable
  !!               are replaced by the upper level.
  !!
  !! History : 4.0  : 04/2020  : J.M. Molines (from cdfmax/cdfmltmask)
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2020
  !! $Id$
  !! Copyright (c) 2020, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class file_informations
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk,  jt
  INTEGER(KIND=4)                               :: narg, iargc, ijarg
  INTEGER(KIND=4)                               :: ni, nj, nk, nt       ! size of the global domain
  INTEGER(KIND=4)                               :: niz, njz             ! size of the domain
  INTEGER(KIND=4)                               :: ikref                ! revel of reference
  INTEGER(KIND=4)                               :: iimin=1, iimax=0     ! i-limit of the domain
  INTEGER(KIND=4)                               :: ijmin=1, ijmax=0     ! j-limit of the domain
  INTEGER(KIND=4)                               :: itmin=1, itmax=0     ! t-limit of the domain
  INTEGER(KIND=4)                               :: istatus,ierr         ! working integer

  REAL(KIND=4)                                  :: zspval               ! missing value or spval
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE     :: v2d, vref            ! running 2d var, and reference var.

  CHARACTER(LEN=256)                            :: cf_in                ! input file name
  CHARACTER(LEN=256)                            :: cf_out               ! input file name
  CHARACTER(LEN=256)                            :: cv_in='none'         ! current variable name
  CHARACTER(LEN=256)                            :: cldum                ! dummy char variable

  LOGICAL                                       :: lout = .FALSE.       ! Flag set when name of output file given

  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdffixanom -f IN-file -v IN-var -reflev kref  ...'
     PRINT *,'      ... [-w imin imax jmin jmax] [-time tmin tmax ] ...'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        In the subdomain, replace the variable value by values at kref,'
     PRINT *,'        from kref+1 to the last wet point. Action is performed on a copy'
     PRINT *,'        of the input file.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file : input file '
     PRINT *,'       -v IN-var  : input variable'
     PRINT *,'       -reflev kref : Give reference level'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-w imin imax jmin jmax] : restrict to sub area specified by the '
     PRINT *,'                  given limits. '
     PRINT *,'       [-time tmin tmax ] : restrict to the indicated time windows.'
     PRINT *,'       [-o OUT-file] : specify the output file name instead of default.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       output is done on <IN-file>_fixanom, or OUT-file if -o option used.'
     STOP 
  ENDIF

  ijarg=1
  ikref=0
  cv_in = 'none'
  DO  WHILE (ijarg <=  narg)
     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE (cldum )
     CASE ( '-f'      ) ; CALL getarg(ijarg, cf_in) ; ijarg = ijarg + 1
     CASE ( '-v'      ) ; CALL getarg(ijarg, cv_in) ; ijarg = ijarg + 1 
     CASE ( '-reflev' ) ; CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ikref
     CASE ( '-w'      ) ; CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) iimin
        ;                 CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) iimax
        ;                 CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmin
        ;                 CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmax
     CASE ( '-time'   ) ; CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) itmin
        ;                 CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) itmax
     CASE ( '-o'      ) ; CALL getarg(ijarg, cf_out); ijarg = ijarg + 1 ; lout=.TRUE.

     CASE DEFAULT     ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  ! Check if reference level given
  IF ( ikref == 0 ) THEN
    PRINT *,' *** ERROR : You must specify a reference level with -reflev '
    STOP 99
  ENDIF
  
  ! Check for variable name
  IF ( TRIM(cv_in) == 'none' ) THEN
    PRINT *,' *** ERROR : You must specify a variable name with -v option.'
    STOP 99
  ENDIF

  ! Check for missing files
  IF ( chkfile(cf_in) ) STOP 99 ! missing file

  IF ( .NOT. lout ) THEN
     cf_out   = TRIM(cf_in)//'_fixanom'
  ENDIF
  CALL system(' dd bs=10000000 if='//TRIM(cf_in)//' of='//TRIM(cf_out) )
  cf_in = cf_out

  PRINT *,' Working on copy : ', TRIM(cf_out)

  ni=0 ; nj=0; nk=0; nt=0 

  ni = getdim(cf_in, cn_x, cldum, istatus)
  IF ( istatus == 1 ) THEN 
     ni = getdim(cf_in, 'lon', cldum, istatus)
     IF ( istatus == 1 ) THEN
        PRINT *,' No X or lon dim found ' ; STOP 99
     ENDIF
  ENDIF
  IF ( iimax == 0 ) iimax = ni

  nj = getdim(cf_in, cn_y, cldum, istatus)
  IF ( istatus == 1 ) THEN 
     nj = getdim(cf_in, 'lat', cldum, istatus)
     IF ( istatus == 1 ) THEN
        PRINT *,' No y or lat dim found ' ; STOP 99
     ENDIF
  ENDIF
  IF ( ijmax == 0 ) ijmax = nj

  nk=getdim(cf_in, cn_z, cldum, istatus)
  IF ( istatus == 1 ) THEN 
     nk = getdim(cf_in, 'z', cldum, istatus)
     IF ( istatus == 1 ) THEN 
        nk = getdim(cf_in, 'lev', cldum, istatus)
        IF ( istatus == 1 ) THEN
           PRINT *,' No dep or z or lev  dim found ' 
           nk = 1
        ENDIF
     ENDIF
  ENDIF
  IF ( nk <= ikref ) THEN
    PRINT *,' *** ERROR : reference level >= total levels in files'
    STOP 99
  ENDIF

  nt = getdim(cf_in, cn_t, cldum, istatus)

  IF ( istatus == 1 ) THEN 
     nt = getdim(cf_in, 'step', cldum, istatus)
     IF ( istatus == 1 ) THEN
        PRINT *,' No time or step dim found ' 
     ENDIF
  ENDIF

  ! fix the size of the zoomed area, or the whole domain if no zoom
  niz = iimax - iimin + 1
  njz = ijmax - ijmin + 1

  IF (nt == 0 ) nt = 1  ! assume a 1 time frame file
  IF ( itmax == 0 ) itmax = nt

  zspval = getatt(cf_in, cv_in, cn_missing_value)

  ALLOCATE (vref(niz,njz), v2d(niz,njz) )

  ! Loop on time
  DO jt=itmin, itmax
   ! read ref value 
   vref(:,:) = getvar(cf_in, cv_in, ikref, niz, njz, kimin=iimin, kjmin=ijmin, ktime=jt)
   ! Loop on level
   DO jk=ikref+1, nk
     v2d(:,:)  = getvar(cf_in, cv_in, jk, niz, njz, kimin=iimin, kjmin=ijmin, ktime=jt)
     WHERE (v2d /= zspval ) v2d=vref
     ierr = putvar(cf_out, cv_in, jk, niz, njz, kimin=iimin, kjmin=ijmin, ptab=v2d,ktime=jt)
   ENDDO
  ENDDO

END PROGRAM cdffixanom
