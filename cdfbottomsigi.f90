PROGRAM cdfbottomsigi
  !!-------------------------------------------------------------------
  !!             ***  PROGRAM cdfbottomsigi  ***
  !!
  !!  **  Purpose: Compute the bottom sigi from gridT file
  !!                Store the results on a 'similar' cdf file.
  !!  
  !!  **  Method: Try to avoid 3 d arrays 
  !!      uses vosaline do determine the bottom points
  !!
  !! history: 
  !!     Original :   J.M. Molines (Nov. 2005)
  !!              :   P. Mathiot (2008) from cdfbottomsig0 to cdfbottomdigi
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
  INTEGER   :: jk                                  !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER, DIMENSION(1) ::  ipk, &                 !: outptut variables : number of levels,
       &                    id_varout              !: ncdf varid's
  real(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: ztemp, zsal ,&   !: Array to read a layer of data
       &                                         ztemp0 , &       !: temporary array to read temp
       &                                         zsal0  , &       !: temporary array to read sal
       &                                         zsigi , &        !: potential density (sig-0)
       &                                         zmask            !: 2D mask at surface
  REAL(KIND=4),DIMENSION(1)                   ::  tim
  REAL(KIND=4)  :: zref
  CHARACTER(LEN=80) :: cdum, cfilet ,cfileout='botsigi.nc' !:
  TYPE (variable), DIMENSION(1) :: typvar         !: structure for attributes

  INTEGER    :: ncout
  INTEGER    :: istatus
  INTEGER, DIMENSION (2) :: ismin, ismax
  REAL(KIND=4)   :: sigmin, sigmax

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfbottomsigi  gridT zref'
     PRINT *,' Output on botsigi.nc, variable sobotsigi'
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  CALL getarg (2, cdum); READ(cdum,*) zref

  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')

  ipk(:)= 1  ! all variables (input and output are 3D)
  typvar(1)%name= 'sobotsigi'
  typvar(1)%units='kg/m3'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= 0.001
  typvar(1)%valid_max= 50.
  typvar(1)%long_name='Bottom_Potential_density'
  typvar(1)%short_name='sobotsigi'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TYX'

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE (ztemp(npiglo,npjglo), zsal(npiglo,npjglo), zsigi(npiglo,npjglo) ,zmask(npiglo,npjglo))
  ALLOCATE (ztemp0(npiglo,npjglo), zsal0(npiglo,npjglo) )

  ! create output fileset

  ncout =create(cfileout, cfilet, npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfilet,npiglo, npjglo,npk)

  zsal = 0.
  ztemp = 0.
  zmask = 1.
  DO jk = 1, npk
     PRINT *,'level ',jk
     zsal0(:,:) = getvar(cfilet, 'vosaline',  jk ,npiglo, npjglo)
     ztemp0(:,:)= getvar(cfilet, 'votemper',  jk ,npiglo, npjglo)
     IF (jk == 1  )  THEN
        tim=getvar1d(cfilet,'time_counter',1)
        WHERE( zsal0 == 0. ) zmask=0.
     END IF
     WHERE ( zsal0 /= 0 )
       zsal=zsal0 ; ztemp=ztemp0
     END WHERE

  ENDDO
     
     zsigi(:,:) = sigmai ( ztemp,zsal,zref,npiglo,npjglo )* zmask(:,:)

     sigmin=minval(zsigi(2:npiglo-1,2:npjglo-1) ,zmask(2:npiglo-1,2:npjglo-1)==1)
     sigmax=maxval(zsigi(2:npiglo-1,2:npjglo-1) ,zmask(2:npiglo-1,2:npjglo-1)==1)
     ismin= minloc(zsigi(2:npiglo-1,2:npjglo-1) ,zmask(2:npiglo-1,2:npjglo-1)==1)
     ismax= maxloc(zsigi(2:npiglo-1,2:npjglo-1) ,zmask(2:npiglo-1,2:npjglo-1)==1)
     PRINT *,'Bottom density : min = ', sigmin,' at ', ismin(1), ismin(2)
     PRINT *,'               : max = ', sigmax,' at ', ismax(1), ismax(2)

     ierr = putvar(ncout, id_varout(1) ,zsigi, 1,npiglo, npjglo)

        ierr=putvar1d(ncout,tim,1,'T')

  istatus = closeout(ncout)
END PROGRAM cdfbottomsigi
