PROGRAM cdfcmp
  !!======================================================================
  !!                     ***  PROGRAM  cdfcmp  ***
  !!======================================================================
  !!  ** Purpose : Find the differences between one same variable in two different files
  !!               Indicate where are located these differences
  !!               Indicate the relative difference
  !!
  !!  ** Method  : Compare var1 and var2
  !!               If it differs, print in standard output where are located diff
  !!               Spatial sub-area restriction can be defined
  !!
  !! History : 3.0 !  08/2012    A. Lecointre   : Original code + Full Doctor form + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class file_operations
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk,jj,ji, jvar, jjvar ! dummy loop index
  INTEGER(KIND=4)                               :: ierr                  ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg    ! argument on line 
  INTEGER(KIND=4)                               :: npiglo, npjglo        ! size fo the domain
  INTEGER(KIND=4)                               :: iimin=1, iimax=0      ! i-limit of the domain
  INTEGER(KIND=4)                               :: ijmin=1, ijmax=0      ! j-limit of the domain
  INTEGER(KIND=4)                               :: ikmin=1, ikmax=0      ! k-limit of the domain
  INTEGER(KIND=4)                               :: nvars                 ! Number of variables in a file
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                   ! arrays of var id's
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: var1, var2            ! variables to compare
  REAL(KIND=4)                                  :: dvar                  ! relative difference
  CHARACTER(LEN=256)                            :: cf1_in,cf2_in         ! input file name
  CHARACTER(LEN=256)                            :: cv_in                 ! variable name
  CHARACTER(LEN=256)                            :: cldum                 ! working string
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names              ! array of var name
  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar               ! Type variable is defined in cdfio.

  !!--------------------------------------------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfcmp -f1 IN-file1 -f2 IN-file2 -v IN-var ...'
     PRINT *,'     ... [-lev kmin kmax ] [-zoom imin imax jmin jmax] ...'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Find where IN-var is different between IN-file1 and IN-file2 '
     PRINT *,'        Options allow to restrict the finding to a sub area in space'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f1 IN-file1 : input file1'
     PRINT *,'       -f2 IN-file2 : input file2'
     PRINT *,'       -v  IN-var   : input variable'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-lev kmin kmax ] : restrict to level between kmin and kmax. '
     PRINT *,'       [-zoom imin imax jmin jmax] : restrict to sub area specified'
     PRINT *,'                                     by the given limits. '
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       output is done on standard output.'
     STOP 
  ENDIF
  !!
  ijarg  = 1
  DO  WHILE (ijarg <=  narg)
     CALL getarg(ijarg,cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-f1' ) ; CALL getarg(ijarg, cf1_in) ; ijarg = ijarg + 1
     CASE ( '-f2' ) ; CALL getarg(ijarg, cf2_in) ; ijarg = ijarg + 1
     CASE ( '-v'  ) ; CALL getarg(ijarg, cv_in ) ; ijarg = ijarg + 1
     CASE ( '-lev') ; CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) ikmin
        ;             CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) ikmax
     CASE ('-zoom') ; CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) iimin
        ;             CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) iimax
        ;             CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmin
        ;             CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmax
     CASE DEFAULT   ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  IF ( chkfile(cf1_in) .OR. chkfile(cf2_in) ) STOP 99 ! missing file
  IF ( chkvar(cf1_in, cv_in) .OR. chkvar(cf2_in, cv_in) ) STOP 99  ! missing var

  npiglo = getdim (cf1_in, cn_x)
  npjglo = getdim (cf1_in, cn_y)
  IF ( iimax == 0 ) iimax = npiglo
  IF ( ijmax == 0 ) ijmax = npjglo

  ! get the number of vertical levels of cv_in variable
  nvars = getnvar(cf1_in)
  ALLOCATE (ipk(nvars),cv_names(nvars),stypvar(nvars))
  cv_names(:)=getvarname(cf1_in,nvars,stypvar)
  ipk(:) = getipk (cf1_in,nvars)
  DO jvar=1,nvars
     IF ( cv_names(jvar) == cv_in ) jjvar=jvar
  ENDDO
  IF ( ikmax == 0 ) ikmax = ipk(jjvar)

  ! Allocate memory.
  ALLOCATE(var1(npiglo, npjglo))
  ALLOCATE(var2(npiglo, npjglo))

  PRINT *,' Working with ', TRIM(cv_in),' defined on ', ipk(jjvar),' level(s)'
  DO jk = ikmin, ikmax
     PRINT *,'# -------------------------------------------'
     PRINT '(A19,I3)','# Checking level: ',jk
     PRINT *,'# i     j    k      var1      var2  %reldiff'
     var1(:,:)=9999.0
     var2(:,:)=9999.0
     var1(:,:) = getvar(cf1_in, cv_in,  jk, npiglo, npjglo)
     var2(:,:) = getvar(cf2_in, cv_in,  jk, npiglo, npjglo)
     DO jj=ijmin, ijmax
        DO ji=iimin, iimax
           IF ( var1(ji,jj) /= var2(ji,jj) ) THEN
              dvar = 100.0*(var1(ji,jj)-var2(ji,jj))/var1(ji,jj)
              PRINT '(I4,2X,I4,2X,I3,2X,F8.3,2X,F8.3,2X,F8.3)',ji,jj,jk,var1(ji,jj),var2(ji,jj),dvar
           ENDIF
        ENDDO
     ENDDO
  ENDDO

END PROGRAM cdfcmp
