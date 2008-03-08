PROGRAM cdfsiginsitu
  !!-------------------------------------------------------------------
  !!             ***  PROGRAM cdfsiginsitu  ***
  !!
  !!  **  Purpose: Compute sigmainsitu  3D field from gridT file
  !!                Store the results on a 'similar' cdf file.
  !!  
  !!  **  Method: read temp and salinity, compute sigmainsitu
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
  INTEGER   :: jk,jt                               !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: 
  INTEGER   :: npiglo,npjglo, npk ,npt             !: size of the domain
  INTEGER, DIMENSION(1) ::  ipk, &                 !: outptut variables : number of levels,
       &                    id_varout              !: ncdf varid's
  real(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: ztemp, zsal ,&   !: Array to read a layer of data
       &                                         zsigi , &        !: potential density (sig-i)
       &                                         zmask            !: 2D mask at current level
  REAL(KIND=4),DIMENSION(:), ALLOCATABLE      :: prof ,tim    !: prof (m) and time (sec)
  REAL(KIND=4)                                :: spval  !: missing value

  CHARACTER(LEN=80) :: cfilet ,cfileout='siginsitu.nc' !:
  CHARACTER(LEN=80) :: cdum

  TYPE(variable) , DIMENSION(1) :: typvar         !: structure for attributes

  INTEGER    :: ncout
  INTEGER    :: istatus
  INTEGER, DIMENSION (2) :: ismin, ismax
  REAL(KIND=4)   :: sigmin, sigmax

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfsiginsitu gridT '
     PRINT *,' Output on siginsitu.nc, variable vosigmainsitu'
     PRINT *,' Depths are taken from  input file '
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')
  npt   = getdim (cfilet,'time')

  ipk(:)= npk  ! all variables (input and output are 3D)
  typvar(1)%name= 'vosigmainsitu'
  typvar(1)%units='kg/m3'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= 0.001
  typvar(1)%valid_max= 45.
  typvar(1)%long_name='in situ density '
  typvar(1)%short_name='vosigmainsitu'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'npt   =', npt

  ALLOCATE (ztemp(npiglo,npjglo), zsal(npiglo,npjglo), zsigi(npiglo,npjglo) ,zmask(npiglo,npjglo))
  ALLOCATE (prof(npk) , tim(npt) )

  ! create output fileset

  ncout =create(cfileout, cfilet, npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfilet,npiglo, npjglo,npk)
  prof(:)=getvar1d(cfilet,'deptht',npk)
  tim=getvar1d(cfilet,'time_counter',npt)
  ierr=putvar1d(ncout,tim,npt,'T')

  spval=getatt(cfilet,'vosaline','missing_value')
  DO jt=1,npt
     PRINT *,'time ',jt, tim(jt)/86400.,' days'
  DO jk = 1, npk
     zmask(:,:)=1.

     ztemp(:,:)= getvar(cfilet, 'votemper',  jk ,npiglo, npjglo)
     zsal(:,:) = getvar(cfilet, 'vosaline',  jk ,npiglo, npjglo)

     WHERE(zsal == spval ) zmask = 0

     zsigi(:,:) = sigmai ( ztemp,zsal,prof(jk),npiglo,npjglo )* zmask(:,:)
     ierr = putvar(ncout, id_varout(1) ,zsigi, jk,npiglo, npjglo)

  END DO  ! loop to next level
  END DO  ! loop to next time

  istatus = closeout(ncout)
END PROGRAM cdfsiginsitu
