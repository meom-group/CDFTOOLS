PROGRAM cdfdifmask
  !!-------------------------------------------------------------------
  !!             ***  PROGRAM cdfdifmask  ***
  !!
  !!  **  Purpose: Build mask file from a salinity output
  !!  
  !!  **  Method:  Read vosaline and set tmask to 1 where sal is not 0
  !!               then umask, vmask and fmask are deduced from tmask
  !!               REM: the result may be locally different for fmask than
  !!                   fmask produced online as there are computed on line
  !!
  !! history: 
  !!     Original :   J.M. Molines November 2005
  !!-------------------------------------------------------------------
  !!  $Rev: 255 $
  !!  $Date: 2009-07-21 17:49:27 +0200 (Tue, 21 Jul 2009) $
  !!  $Id: cdfdifmask.f90 255 2009-07-21 15:49:27Z molines $
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk,jt, jvar                   !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc , ntags                 !: 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER, DIMENSION(4) ::  ipk, &                 !: outptut variables : number of levels,
       &                    id_varout              !: ncdf varid's
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE ::     zmask,zmask2            !: 2D mask at current level

  CHARACTER(LEN=256)                 :: cvar       !: array of var name
  CHARACTER(LEN=256)                 :: cfile1, cfile2, cline,cfileout='mask_diff.nc'
  TYPE(variable), DIMENSION(4) :: typvar
  REAL(KIND=4) ,DIMENSION(1)                  :: timean

  INTEGER    :: ncout, npt
  INTEGER    :: istatus
  REAL(4) :: ss

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfdifmask  mask1 mask2'
     STOP
  ENDIF

  CALL getarg (1, cfile1)
  CALL getarg (2, cfile2)
  npiglo= getdim (cfile1,'x')
  npjglo= getdim (cfile1,'y')
  npk   = getdim (cfile1,'depth')

   print *, npiglo, npjglo, npk

  ipk(1:4)      = npk
  typvar(1)%name='tmask'
  typvar(2)%name='umask'
  typvar(3)%name='vmask'
  typvar(4)%name='fmask'
  typvar(1:4)%units='1/0'
  typvar(1:4)%missing_value=9999.
  typvar(1:4)%valid_min= 0.
  typvar(1:4)%valid_max= 1.
  typvar(1)%long_name='tmask'
  typvar(2)%long_name='umask'
  typvar(3)%long_name='vmask'
  typvar(4)%long_name='fmask'
  typvar(1)%short_name='tmask'
  typvar(2)%short_name='umask'
  typvar(3)%short_name='vmask'
  typvar(4)%short_name='fmask'
  typvar(1:4)%online_operation='N/A'
  typvar(1:4)%axis='TZYX'
  typvar(1:4)%precision='i2'

  ncout =create(cfileout, cfile1,npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,4, ipk,id_varout )
  ierr= putheadervar(ncout, cfile1, npiglo, npjglo,npk)


  ALLOCATE (zmask(npiglo,npjglo),zmask2(npiglo,npjglo))

  npt= 0
  DO jvar=1,4
     cvar=typvar(jvar)%name
  DO jk=1, npk
     zmask(:,:)= getvar(cfile1, cvar,  jk ,npiglo, npjglo)
     zmask2(:,:)= getvar(cfile2, cvar,  jk ,npiglo, npjglo)
     zmask(:,:)= zmask2(:,:) - zmask(:,:)
     ierr=putvar(ncout,id_varout(jvar), zmask, jk ,npiglo, npjglo)
  END DO  ! loop to next level
  END DO
  timean(:)=0.
  ierr=putvar1d(ncout,timean,1,'T')
  istatus = closeout(ncout)


END PROGRAM cdfdifmask
