PROGRAM cdfnan
  !!======================================================================
  !!                     ***  PROGRAM  cdfnan  ***
  !!=====================================================================
  !!  ** Purpose : Replace the nan values by spval or another value 
  !!               given in argument
  !!
  !! History : 2.1  : 05/2010  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
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

  INTEGER(KIND=4)                               :: jk, jvar, jt           ! dummy loop index
  INTEGER(KIND=4)                               :: ierr                   ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg     ! browse line
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
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names               ! array of var name

  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvar                ! type for attributes

  LOGICAL                                       :: l_replace = .false.    ! flag for replace value
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfnan list_of_model_output_files [-value replace] [-absmax rabsmax ] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Detect NaN values in the input files, and change them to '
     PRINT *,'       either spval (missing_value) or the value given as option.'
     PRINT *,'       Does the same for absolute values > huge(0.0)'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       list of model output files. They must be of same type and have'
     PRINT *,'       similar sizes. CAUTION : input files are rewritten !'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-value replace ] : use replace instead of missing_value for'
     PRINT *,'                           changing NaN.'
     PRINT *,'       [-absmax rabsmax ] : replace values whose absolute value is greater '
     PRINT *,'                           than rabsmax.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : input file is rewritten without NaN.' 
     PRINT *,'         variables : same name as input.' 
     STOP
  ENDIF

  rabsmax=huge(0.0)
  ijarg=1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum) ;   ijarg = ijarg + 1
     SELECT CASE (cldum)
     CASE ('-value' )
        CALL getarg( ijarg, cldum) ; ijarg = ijarg+1 ; 
        READ(cldum,*) replace ; l_replace=.true.
     CASE ('-absmax' )
        CALL getarg( ijarg, cldum) ; ijarg = ijarg+1 ; 
        READ(cldum,*) rabsmax 
     CASE DEFAULT
        cf_inout=TRIM(cldum)
     END SELECT
  END DO
  IF ( chkfile (cf_inout) )  STOP ! missing file

  npiglo = getdim (cf_inout, cn_x              )
  npjglo = getdim (cf_inout, cn_y              )
  npk    = getdim (cf_inout, cn_z, kstatus=ierr)

  IF (ierr /= 0 ) THEN
     npk   = getdim (cf_inout,'z',kstatus=ierr)
     IF (ierr /= 0 ) THEN
        PRINT *, 'ASSUME NO VERTICAL DIMENSIONS !'
        npk=0
     ENDIF
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

  !re scan argument list 
  ijarg = 1
  DO WHILE (ijarg <=  narg ) 
     CALL getarg (ijarg, cf_inout) ; ijarg = ijarg + 1 

     SELECT CASE ( cf_inout)
     CASE ('-value' )
        ! replace already read, just skip
        ijarg = ijarg + 1
     CASE ('-absmax' )
        !  already read, just skip
        ijarg = ijarg + 1
     CASE DEFAULT  ! reading files
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
                    WHERE( tab(:,:) < -rabsmax ) tab(:,:) = zspval
                    WHERE( tab(:,:) >  rabsmax ) tab(:,:) = zspval
                    ierr = putvar(ncid, id_var(jvar), tab, jk, npiglo, npjglo, ktime=jt)
                 ENDDO
              END DO
           ENDIF
        ENDDO
     END SELECT
  ENDDO

  ierr = closeout(ncid)

END PROGRAM cdfnan
