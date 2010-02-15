PROGRAM cdfsig0
  !!-------------------------------------------------------------------
  !!             ***  PROGRAM cdfsig0  ***
  !!
  !!  **  Purpose: Compute sigma0 3D field from gridT file
  !!                Store the results on a 'similar' cdf file.
  !!  
  !!  **  Method: Try to avoid 3 d arrays 
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
  INTEGER   :: jk, jt                              !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: 
  INTEGER   :: npiglo,npjglo, npk, npt             !: size of the domain
  INTEGER, DIMENSION(1) ::  ipk, &                 !: outptut variables : number of levels,
       &                    id_varout              !: ncdf varid's
  real(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: ztemp, zsal ,&   !: Array to read a layer of data
       &                                         zsig0 , &        !: potential density (sig-0)
       &                                         zmask            !: 2D mask at current level
  REAL(KIND=4),DIMENSION(:),ALLOCATABLE   ::  tim

  CHARACTER(LEN=256) :: cfilet ,cfileout='sig0.nc' !:

  TYPE(variable) , DIMENSION(1) :: typvar         !: structure for attributes

  INTEGER    :: ncout
  INTEGER    :: istatus
  INTEGER, DIMENSION (2) :: ismin, ismax
  REAL(KIND=4)   :: sigmin, sigmax

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfsig0  gridT '
     PRINT *,' Output on sig0.nc, variable vosigma0'
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')
  npt   = getdim (cfilet,'time')

  ipk(:)= npk  ! all variables (input and output are 3D)
  typvar(1)%name= 'vosigma0'
  typvar(1)%units='kg/m3'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= 0.001
  typvar(1)%valid_max= 40.
  typvar(1)%long_name='Potential_density:sigma-0'
  typvar(1)%short_name='vosigma0'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'npt   =', npt

  ALLOCATE (ztemp(npiglo,npjglo), zsal(npiglo,npjglo), zsig0(npiglo,npjglo) ,zmask(npiglo,npjglo))
  ALLOCATE (tim(npt))

  ! create output fileset

  ncout =create(cfileout, cfilet, npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfilet,npiglo, npjglo,npk)
  tim=getvar1d(cfilet,'time_counter',npt)
  ierr=putvar1d(ncout,tim,npt,'T')

  DO jt=1,npt
    PRINT *,' TIME = ', jt, tim(jt)/86400.,' days'
  DO jk = 1, npk
     zmask(:,:)=1.

     ztemp(:,:)= getvar(cfilet, 'votemper',  jk ,npiglo, npjglo,ktime=jt)
     zsal(:,:) = getvar(cfilet, 'vosaline',  jk ,npiglo, npjglo,ktime=jt)

     WHERE(zsal == 0 ) zmask = 0

     zsig0(:,:) = sigma0 ( ztemp,zsal,npiglo,npjglo )* zmask(:,:)
     
!     sigmin=minval(zsig0(2:npiglo-1,2:npjglo-1) ,zmask(2:npiglo-1,2:npjglo-1)==1)
!     sigmax=maxval(zsig0(2:npiglo-1,2:npjglo-1) ,zmask(2:npiglo-1,2:npjglo-1)==1)
!     ismin= minloc(zsig0(2:npiglo-1,2:npjglo-1) ,zmask(2:npiglo-1,2:npjglo-1)==1)
!     ismax= maxloc(zsig0(2:npiglo-1,2:npjglo-1) ,zmask(2:npiglo-1,2:npjglo-1)==1)
!     PRINT *,'Level ',jk,': min = ', sigmin,' at ', ismin(1), ismin(2)
!     PRINT *,'               : max = ', sigmax,' at ', ismax(1), ismax(2)

     ierr = putvar(ncout, id_varout(1) ,zsig0, jk,npiglo, npjglo,ktime=jt)

  END DO  ! loop to next level
  END DO  ! next time frame

  istatus = closeout(ncout)
END PROGRAM cdfsig0
