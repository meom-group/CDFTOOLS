PROGRAM cdfbottom
  !!-------------------------------------------------------------------
  !!             ***  PROGRAM cdfbottom  ***
  !!
  !!  **  Purpose: Extract the bottom value for the 3D variables
  !!               which are in the input file. Store the results
  !!               on a similar file, with the same variable name.
  !!  
  !!  **  Method: Try to avoid 3 d arrays 
  !!      Uses the corresponding mask file to determine the bottom.
  !!      If no mask found it assumes that 0.0000 values corresponds
  !!      to masked values.
  !!
  !! history: 
  !!     Original :   J.M. Molines (Nov. 2005)
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk , jv, jvar                       !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: nvars                               !: number of variables in the input file
  INTEGER, DIMENSION(:), ALLOCATABLE :: ipk,ipko,& !: outptut variables : number of levels,
       &               id_var,     id_varout       !: ncdf varid's
  real(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: zfield ,&   !: Array to read a layer of data
       &                                         zbot , &    !: array to store the bottom value
       &                                         zmask       !: 2D mask at current level
  REAL(KIND=4),DIMENSION(1)                   ::  tim

  CHARACTER(LEN=80) :: cfile, cdum ,cmask='mask.nc',cfileout='bottom.nc' !:
  CHARACTER(LEN=1) :: ctype=' '
  CHARACTER(LEN=5) :: cvmask=' '                                   !: name of the mask variable
  CHARACTER(LEN=80) ,DIMENSION(:), ALLOCATABLE   :: cvarname       !: array of var name
  
  TYPE (variable), DIMENSION(:), ALLOCATABLE :: typvar             !: structure for variable attribute

  INTEGER    :: ncout
  INTEGER    :: istatus

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfbottom  ncfile [ T | U | V | F]'
     PRINT *,'     grid point type is optional: if not given'
     PRINT *,'  it does''nt require the mask.nc file and'
     PRINT *,'  assumes that data points with 0 are land points '
     PRINT *,' Output on bottom.nc, variables as in the input file'
     STOP
  ENDIF

  CALL getarg (1, cfile)
  npiglo= getdim (cfile,'x')
  npjglo= getdim (cfile,'y')
  npk   = getdim (cfile,'depth')

  ALLOCATE (zfield(npiglo,npjglo), zbot(npiglo,npjglo),zmask(npiglo,npjglo))

  IF (narg == 2 ) THEN
     CALL getarg (2, ctype )

     SELECT CASE ( ctype )
     CASE ( 'T', 't', 'S', 's' )
        cvmask='tmask'
     CASE ( 'U', 'u' )
        cvmask='umask'
     CASE ( 'V', 'v' )
        cvmask='vmask'
     CASE ( 'F', 'f' )
        cvmask='fmask'
        PRINT *, 'Be carefull with fmask ... !!!'
     CASE DEFAULT
        PRINT *, ' ERROR : This type of point ', ctype,' is not known !'
        STOP
     END SELECT

  ENDIF

  ! look for the number of variables in the input file
  nvars = getnvar(cfile)
  PRINT *,' nvars =', nvars

  ALLOCATE (cvarname(nvars) ,typvar(nvars))
  ALLOCATE (id_var(nvars),ipk(nvars),id_varout(nvars) ,ipko(nvars) )

  cvarname(:)=getvarname(cfile,nvars,typvar)

  id_var(:)  = (/(jv, jv=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cfile,nvars)
  ipko(:)= 1  ! all variables output are 2D

  WHERE( ipk <= 1 ) cvarname='none'
! typvar%name=cvarname
! typvar%axis='TYX'
  DO jvar=1,nvars
    typvar(jvar)%name=cvarname(jvar)
    typvar(jvar)%axis='TYX'
    cdum=typvar(jvar)%long_name
    typvar(jvar)%long_name='Bottom '//TRIM(cdum)
  END DO
  ! create output fileset
  ! create output file taking the sizes in cfile

  ncout =create(cfileout, cfile,npiglo,npjglo,1)

  ierr= createvar(ncout , typvar,  nvars, ipko, id_varout )

  ierr= putheadervar(ncout , cfile, npiglo, npjglo, 1)

  DO jvar = 1,nvars
     zfield = 0.
     zbot   = 0.

     IF (cvarname(jvar) == 'none' ) THEN
        ! skip these variable
     ELSE
        PRINT *, ' WORKING with ', TRIM( cvarname(jvar) ), ipk(jvar)
        DO jk = 1, ipk(jvar)
           zmask = 1.
           zfield(:,:) = getvar(cfile, cvarname(jvar),  jk ,npiglo, npjglo)
           IF (jk == 1 .AND. jvar == 1. )  THEN
              tim=getvar1d(cfile,'time_counter',1)
           ENDIF
           IF ( cvmask == ' ' ) THEN
              WHERE ( zfield /= 0 )
                 zbot = zfield
              END WHERE
           ELSE
              zmask(:,:) = getvar(cmask, cvmask, jk, npiglo, npjglo)
              WHERE ( zmask /= 0 )
                 zbot = zfield
              END WHERE
           ENDIF
        END DO
        ierr = putvar(ncout, id_varout(jvar) ,zbot, 1,npiglo, npjglo)
     ENDIF

  END DO
  ierr=putvar1d(ncout,tim,1,'T')

  istatus = closeout(ncout)
END PROGRAM cdfbottom
