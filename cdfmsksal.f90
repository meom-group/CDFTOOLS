PROGRAM cdfmsksal
  !!-------------------------------------------------------------------
  !!             ***  PROGRAM cdfmsksal  ***
  !!
  !!  **  Purpose: Computes the number of land points from the gridT file
  !!  
  !!  **  Method: Try to avoid 3 d arrays gridT file Work with vosaline
  !!
  !! history: 
  !!     Original :   J.M. Molines November 2005
  !!-------------------------------------------------------------------
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
  CHARACTER(LEN=80)                 :: cfilet, cline

  INTEGER    :: ncout, npt
  INTEGER    :: istatus
  REAL(4) :: ss

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmsksal  gridT '
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'z')

  ALLOCATE (zmask(npiglo,npjglo))

  npt= 0
!  DO jk=1, npk
   DO jk=1, 1
     zmask(:,:)= getvar(cfilet, 'vosaline',  jk ,npiglo, npjglo)
     WHERE (zmask > 0 ) zmask = 1
     ss=sum(zmask)
     print *, jk, ss, ss/npiglo/npjglo*100,' % H'
     npt = npt + ss
  END DO  ! loop to next level
  OPEN (10, FILE='tmask.bimg',FORM='UNFORMATTED')
  cline='tmask(1) from '//trim(cfilet)//' (cdfmsksal)'
  WRITE(10) cline
  WRITE(10) cline
  WRITE(10) cline
  WRITE(10) cline
  WRITE(10) npiglo, npjglo,1,1,1,0
  WRITE(10) 1.,1.,1.,1., 9999.
  WRITE(10) 0.
  WRITE(10) 0.
  WRITE(10)(( zmask(ji,jj),ji=1,npiglo),jj=1,npjglo)
  CLOSE(10)

  PRINT *, ' Number of Ocean points :', npt,'  ',(1.*npt)/npiglo/npjglo/npk*100,' %'
  PRINT *, ' Number of Land points :', npiglo*npjglo*npk - npt,'  ',(npiglo*npjglo*npk -1.*npt)/npiglo/npjglo/npk*100,' %'

END PROGRAM cdfmsksal
