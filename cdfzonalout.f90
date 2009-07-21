PROGRAM cdfzonalout
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfzonalout  ***
  !!
  !!  **  Purpose  :  Output zonal mean/integral as ascii files
  !!  
  !!  **  Method   :  
  !!     Read zonalmean or zonalsum file, determine 1D variable and dump them on an ASCII file
  !!
  !!
  !! history ;
  !!  Original :  J.M. Molines (Feb. 2006) 
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jbasin, jj, jk ,ji ,jvar ,jjvar     !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: nvars , mvar                        !: number of variables in the file
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipk, ijvar, ipko, id_varout    !: jpbasin x nvar

  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlon              !: dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlat              !: latitude for i = north pole
  REAL(KIND=4), DIMENSION (1)                    ::  tim

  REAL(KIND=8), DIMENSION (:,:,:),   ALLOCATABLE ::  zv

  CHARACTER(LEN=256) :: cfilev 
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE   :: cvarname             !: array of var name for input
  TYPE(variable), DIMENSION(:),ALLOCATABLE :: typvar

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfzonalout  file  '
     STOP
  ENDIF

  CALL getarg (1, cfilev)

  nvars  = getnvar(cfilev)
  ALLOCATE ( cvarname(nvars)  ,ipk(nvars), ijvar(nvars), typvar(nvars)  )
  cvarname(1:nvars) = getvarname(cfilev,nvars,typvar)
  ipk(1:nvars) = getipk(cfilev,nvars)

  ! Open standard output with reclen 2048 for avoid wrapping with ifort
  OPEN(6,FORM='FORMATTED',RECL=2048)
  ! look for 1D var ( f(lat) )
  mvar = 0
  DO jvar = 1,nvars
     ! skip variables such as nav_lon, nav_lat, time_counter deptht ...
     IF (ipk(jvar) == 0 .OR. ipk(jvar) > 1 ) THEN
        cvarname(jvar)='none'
     ELSE
        mvar = mvar + 1       ! count for valid input variables
        ijvar(mvar) = jvar    ! use indirect adressing for those variables
     ENDIF
  END DO
  WRITE(6,*) 'Number of 1D variables :', mvar
  DO jjvar=1,mvar
    jvar=ijvar(jjvar)
    WRITE(6,*) '     ',TRIM(cvarname(jvar))
  ENDDO
   

  npiglo= getdim (cfilev,'x')
  npjglo= getdim (cfilev,'y')
  npk   = getdim (cfilev,'depth')


  WRITE(6,*)  'npiglo=', npiglo
  WRITE(6,*)  'npjglo=', npjglo
  WRITE(6,*)  'npk   =', npk

  ! Allocate arrays
  ALLOCATE ( zv(npiglo,npjglo,mvar) )
  ALLOCATE ( dumlon(1,npjglo) , dumlat(1,npjglo))


  dumlat(:,:) = getvar(cfilev,'nav_lat',1,1,npjglo)

  ! main computing loop
  DO jjvar = 1, mvar
     jvar = ijvar(jjvar)
        ! Get variables and mask at level jk
        zv(:,:,jjvar)       = getvar(cfilev,   cvarname(jvar),  1 ,1,npjglo)

  END DO ! next variable

   WRITE(6,*) ' J  LAT ', (TRIM(cvarname(ijvar(jjvar))),' ',jjvar=1,mvar)
  DO jj=npjglo,1,-1
    WRITE(6,*)  jj, dumlat(1,jj), zv(1,jj,1:mvar)
  ENDDO


END PROGRAM cdfzonalout
