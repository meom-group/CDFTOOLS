PROGRAM cdfmsk
  !!-------------------------------------------------------------------
  !!             ***  PROGRAM cdfmsk  ***
  !!
  !!  **  Purpose: Computes the number of land points from the mask
  !!  
  !!  **  Method: Try to avoid 3 d arrays 
  !!
  !! history: 
  !!     Original :   J.M. Molines May 2005
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
  INTEGER, DIMENSION(1) ::  ipk, &                 !: outptut variables : number of levels,
       &                    id_varout              !: ncdf varid's
  real(KIND=4) , DIMENSION (:,:), ALLOCATABLE ::     zmask            !: 2D mask at current level

  CHARACTER(LEN=80) ,DIMENSION(1)   :: cvarname   !: array of var name
  CHARACTER(LEN=80)                 :: cfilet

  INTEGER    :: ncout, npt
  INTEGER    :: istatus
  REAL(4) :: ss

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmsk  maskfile '
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'z')

  ALLOCATE (zmask(npiglo,npjglo))

  npt= 0
  DO jk=1, npk
     zmask(:,:)= getvar(cfilet, 'tmask',  jk ,npiglo, npjglo)
     ss=sum(zmask)
     print *, jk, ss, ss/npiglo/npjglo*100,' % H'
     npt = npt + ss
  END DO  ! loop to next level
  PRINT *, ' Number of Ocean points :', npt,'  ',(1.*npt)/npiglo/npjglo/npk*100,' %'
  PRINT *, ' Number of Land points :', npiglo*npjglo*npk - npt,'  ',(npiglo*npjglo*npk -1.*npt)/npiglo/npjglo/npk*100,' %'

END PROGRAM cdfmsk
