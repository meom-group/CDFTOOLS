PROGRAM cdfsigi
  !!-------------------------------------------------------------------
  !!             ***  PROGRAM cdfsigi  ***
  !!
  !!  **  Purpose: Compute sigmai 3D field from gridT file
  !!                Store the results on a 'similar' cdf file.
  !!  
  !!  **  Method: read temp and salinity, compute sigma-i
  !!              using depth given in argument (meters or dbar)
  !!
  !! history: 
  !!     Original :   J.M. Molines (Nov 2004 ) for ORCA025
  !!                  J.M. Molines Apr 2005 : use modules
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio
  USE eos

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk , jt                             !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: 
  INTEGER   :: npiglo,npjglo, npk, npt             !: size of the domain
  INTEGER, DIMENSION(1) ::  ipk, &                 !: outptut variables : number of levels,
       &                    id_varout              !: ncdf varid's
  real(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: ztemp, zsal ,&   !: Array to read a layer of data
       &                                         zsigi , &        !: potential density (sig-i)
       &                                         zmask            !: 2D mask at current level
  REAL(KIND=4),DIMENSION(:),ALLOCATABLE       :: tim
  REAL(KIND=4)                                :: prof=0.! in meters
  REAL(KIND=4)                                :: spval  !: missing value

  CHARACTER(LEN=80) :: cfilet ,cfileout='sigi.nc' !:
  CHARACTER(LEN=80) :: cdum

  TYPE(variable) , DIMENSION(1) :: typvar         !: structure for attributes

  INTEGER    :: ncout
  INTEGER    :: istatus
  INTEGER, DIMENSION (2) :: ismin, ismax
  REAL(KIND=4)   :: sigmin, sigmax

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfsigi  gridT Ref_dep(m)'
     PRINT *,' Output on sigi.nc, variable vosigmai'
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  CALL getarg (2, cdum)
  READ(cdum,*) prof
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')
  npt   = getdim (cfilet,'time')

  ipk(:)= npk  ! all variables (input and output are 3D)
  typvar(1)%name= 'vosigmai'
  typvar(1)%units='kg/m3'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= 0.001
  typvar(1)%valid_max= 45.
  typvar(1)%long_name='Potential_density:refered to '//TRIM(cdum)//' m'
  typvar(1)%short_name='vosigmai'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'npt   =', npt

  ALLOCATE (ztemp(npiglo,npjglo), zsal(npiglo,npjglo), zsigi(npiglo,npjglo) ,zmask(npiglo,npjglo))
  ALLOCATE (tim(npt) )

  ! create output fileset

  ncout =create(cfileout, cfilet, npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfilet,npiglo, npjglo,npk)

  spval=getatt(cfilet,'vosaline','missing_value')

  tim=getvar1d(cfilet,'time_counter',npt)
  ierr=putvar1d(ncout,tim,npt,'T')
  DO jt = 1, npt
     PRINT *,'time: ',jt
  DO jk = 1, npk
     zmask(:,:)=1.

     ztemp(:,:)= getvar(cfilet, 'votemper',  jk ,npiglo, npjglo,ktime=jt)
     zsal(:,:) = getvar(cfilet, 'vosaline',  jk ,npiglo, npjglo,ktime=jt)

     WHERE(zsal == spval ) zmask = 0

     zsigi(:,:) = sigmai ( ztemp,zsal,prof,npiglo,npjglo )* zmask(:,:)
     IF ( npiglo /= 1 .AND. npjglo /= 1 ) THEN
       sigmin=minval(zsigi(2:npiglo-1,2:npjglo-1) ,zmask(2:npiglo-1,2:npjglo-1)==1)
       sigmax=maxval(zsigi(2:npiglo-1,2:npjglo-1) ,zmask(2:npiglo-1,2:npjglo-1)==1)
       ismin= minloc(zsigi(2:npiglo-1,2:npjglo-1) ,zmask(2:npiglo-1,2:npjglo-1)==1)
       ismax= maxloc(zsigi(2:npiglo-1,2:npjglo-1) ,zmask(2:npiglo-1,2:npjglo-1)==1)
       PRINT *,'Level ',jk,': min = ', sigmin,' at ', ismin(1), ismin(2)
       PRINT *,'               : max = ', sigmax,' at ', ismax(1), ismax(2)
     ENDIF

     ierr = putvar(ncout, id_varout(1) ,zsigi, jk,npiglo, npjglo,ktime=jt)

  END DO  ! loop to next level
  END DO  ! loop on time

  istatus = closeout(ncout)
END PROGRAM cdfsigi
