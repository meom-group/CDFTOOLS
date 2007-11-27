PROGRAM cdfmkmask
  !!-------------------------------------------------------------------
  !!             ***  PROGRAM cdfmkmask  ***
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
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk,jt                         !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc , ntags                 !: 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER, DIMENSION(4) ::  ipk, &                 !: outptut variables : number of levels,
       &                    id_varout              !: ncdf varid's
  real(KIND=4) , DIMENSION (:,:), ALLOCATABLE ::     zmask,zmask2            !: 2D mask at current level

  CHARACTER(LEN=80) ,DIMENSION(4)   :: cvarname   !: array of var name
  CHARACTER(LEN=80)                 :: cfilet, cline,cfileout='mask_sal.nc'
  TYPE(variable), DIMENSION(4) :: typvar
  REAL(KIND=4) ,DIMENSION(1)                  :: timean

  INTEGER    :: ncout, npt
  INTEGER    :: istatus
  REAL(4) :: ss

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmkmask  gridT '
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')

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

  ncout =create(cfileout, cfilet,npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,4, ipk,id_varout )
  ierr= putheadervar(ncout, cfilet, npiglo, npjglo,npk)


  ALLOCATE (zmask(npiglo,npjglo),zmask2(npiglo,npjglo))

  npt= 0
  DO jk=1, npk
     zmask(:,:)= getvar(cfilet, 'vosaline',  jk ,npiglo, npjglo)
     WHERE (zmask > 0 ) zmask = 1
     ierr=putvar(ncout,id_varout(1), zmask, jk ,npiglo, npjglo)
     ! now umask
     zmask2=0.
     DO ji=1,npiglo-1
       DO jj=1,npjglo
        zmask2(ji,jj)=zmask(ji,jj)*zmask(ji+1,jj)
       END DO
     END DO
     ierr=putvar(ncout,id_varout(2), zmask2, jk ,npiglo, npjglo)

    ! now vmask
     zmask2=0.
     DO ji=1,npiglo
       DO jj=1,npjglo-1
        zmask2(ji,jj)=zmask(ji,jj)*zmask(ji,jj+1)
       END DO
     END DO
     ierr=putvar(ncout,id_varout(3), zmask2, jk ,npiglo, npjglo)

     !now fmask
     zmask2=0.
     DO ji=1,npiglo-1
       DO jj=1,npjglo-1
        zmask2(ji,jj)=zmask(ji,jj)*zmask(ji,jj+1)*zmask(ji+1,jj)*zmask(ji+1,jj+1)
       END DO
     END DO
     ierr=putvar(ncout,id_varout(4), zmask2, jk ,npiglo, npjglo)
  END DO  ! loop to next level
  timean(:)=0.
  ierr=putvar1d(ncout,timean,1,'T')
  istatus = closeout(ncout)


END PROGRAM cdfmkmask
