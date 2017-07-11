PROGRAM cdfnan
  !!======================================================================
  !!                     ***  PROGRAM  cdfnan  ***
  !!=====================================================================
  !!  ** Purpose : Replace the nan values by spval or another value 
  !!               given in argument
  !!
  !! History : 2.1  : 05/2010  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
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

  INTEGER(KIND=4)                               :: jk, jvar, jfil, jt     ! dummy loop index
  INTEGER(KIND=4)                               :: ierr                   ! working integer
  INTEGER(KIND=4)                               :: idep, idep_max         ! possible depth index, maximum
  INTEGER(KIND=4)                               :: narg, iargc, ijarg     ! browse line
  INTEGER(KIND=4)                               :: nfiles                 ! number of files to process
  INTEGER(KIND=4)                               :: ncid                   ! ncid of input file for rewrite
  INTEGER(KIND=4)                               :: npiglo, npjglo         ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt               ! size of the domain
  INTEGER(KIND=4)                               :: nvars                  ! Number of variables in a file
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                    ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var                 ! arrays of var id

  REAL(KIND=4)                                  :: zspval, replace        ! spval, replace value
  REAL(KIND=4)                                  :: rabsmax                ! spval, replace value
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: tab                    ! Arrays for data

  CHARACTER(LEN=256)                            :: cldum                  ! dummy string for getarg
  CHARACTER(LEN=256)                            :: cf_inout               ! file name
  CHARACTER(LEN=256)                            :: cunits, clname, csname ! attributes
  CHARACTER(LEN=256)                            :: cv_dep                 ! name of dep dimension
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names               ! array of var name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cf_lst                 ! list of file to process
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: clv_dep                ! array of possible depth name (or 3rd dimension)

  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvar                ! type for attributes

  LOGICAL                                       :: l_replace = .FALSE.    ! flag for replace value
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfnan -l LST-files [-r value] [-absmax max] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Detect NaN values in the input files, and change them to either spval'
     PRINT *,'       (missing_value) or the value given with -r option.'
     PRINT *,'      '
     PRINT *,'       When absolute value is larger than huge or that the value given with'
     PRINT *,'       the -absmax option, it is also replaced as NaN are replaced.'
     PRINT *,'      '
     PRINT *,'     CAUTION :'
     PRINT *,'      ################################'
     PRINT *,'      # INPUT FILES ARE OVER-WRITTEN #'
     PRINT *,'      ################################'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -l LST-files : A blank-separated list of the name of the files to '
     PRINT *,'              process. All files in the list must have the same geometry and'
     PRINT *,'              must contain the same variables.'
     PRINT *,'              CAUTION : input files are over-written!'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-r value] : use value instead of missing_value for replacing NaN.'
     PRINT *,'       [-absmax rabsmax ] : replace values whose absolute value is greater'
     PRINT *,'                    than max.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : input file is rewritten without NaN.' 
     PRINT *,'         variables : same name as input.' 
     PRINT *,'      '
     STOP 
  ENDIF

  rabsmax=HUGE(0.0)

  ijarg=1 ; nfiles=0
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum) ;   ijarg = ijarg + 1
     SELECT CASE (cldum)
     CASE ('-l'      ) ; CALL GetFileList
        ! options
     CASE ('-value'  ) ; CALL getarg( ijarg, cldum) ; ijarg = ijarg+1 ; READ(cldum,*) replace ; l_replace=.TRUE.
     CASE ('-absmax' ) ; CALL getarg( ijarg, cldum) ; ijarg = ijarg+1 ; READ(cldum,*) rabsmax 
     CASE DEFAULT      ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  cf_inout = cf_lst(1)
  IF ( chkfile (cf_inout) )  STOP 99 ! missing file

  npiglo = getdim (cf_inout, cn_x              )
  npjglo = getdim (cf_inout, cn_y              )

  ! looking for npk among various possible name
  idep_max=8
  ALLOCATE ( clv_dep(idep_max) )
  clv_dep(:) = (/cn_z,'z','sigma','nav_lev','levels','ncatice','icbcla','icbsect'/)
  idep=1  ; ierr=1000
  DO WHILE ( ierr /= 0 .AND. idep <= idep_max )
     npk  = getdim (cf_inout, clv_dep(idep), cdtrue=cv_dep, kstatus=ierr)
     idep = idep + 1
  ENDDO

  IF ( ierr /= 0 ) THEN  ! none of the dim name was found
     PRINT *,' assume file with no depth'
     npk=0
  ENDIF

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( tab(npiglo,npjglo) )

  nvars = getnvar(cf_inout)

  ALLOCATE (cv_names(nvars), id_var(nvars),ipk(nvars), stypvar(nvars))

  cv_names(:) = getvarname(cf_inout,nvars,stypvar)
  ipk(:)      = getipk(cf_inout,nvars)
  id_var(:)   = getvarid(cf_inout,nvars)

  DO jfil = 1, nfiles
     cf_inout = cf_lst(jfil)
     PRINT *, 'Change NaN on file ', cf_inout
     ncid = ncopen(cf_inout)
     npt  = getdim (cf_inout,cn_t)

     DO jvar = 1,nvars
        IF (    cv_names(jvar) == cn_vlon2d   .OR. &
             &  cv_names(jvar) == cn_vlat2d   .OR. &
             &  cv_names(jvar) == cn_vtimec   .OR. &
             &  cv_names(jvar) == cn_vdeptht  .OR. &
             &  cv_names(jvar) == cn_vdepthu  .OR. &
             &  cv_names(jvar) == cn_vdepthv         )  THEN
           ! skip these variable
        ELSE
           IF ( l_replace ) THEN
              zspval=replace
           ELSE
              ierr = getvaratt (cf_inout, cv_names(jvar), cunits, zspval, clname, csname)
           ENDIF

           DO jt=1,npt
              DO jk = 1, ipk(jvar) 
                 tab(:,:) = getvar(cf_inout, cv_names(jvar), jk, npiglo, npjglo, ktime=jt )
                 !                   WHERE( isnan(tab(:,:)) ) tab(:,:) = zspval
                 ! isnan function is not available on xlf90 compiler
                 ! we replace it by the following test that gives the same results
                 ! reference : http://www.unixguide.net/ibm/faq/faq3.03.shtml
                 WHERE( tab(:,:) /=  tab(:,:) ) tab(:,:) = zspval
                 WHERE( tab(:,:) < -rabsmax )   tab(:,:) = zspval
                 WHERE( tab(:,:) >  rabsmax )   tab(:,:) = zspval
                 ierr = putvar(ncid, id_var(jvar), tab, jk, npiglo, npjglo, ktime=jt)
              ENDDO
           END DO
        ENDIF
     ENDDO
  ENDDO

  ierr = closeout(ncid)

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


END PROGRAM cdfnan
