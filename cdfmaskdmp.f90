PROGRAM cdfmaskdmp
  !!-------------------------------------------------------------------
  !!             ***  PROGRAM cdfmaskdmp  ***
  !!
  !!  **  Purpose: Compute 3D mask for AABW relaxation from T and S climatologies
  !!                Store the results on a cdf file.
  !!  
  !!  **  Method: read temp and salinity, compute sigma-2
  !!              compute coefs, create mask
  !!
  !! history: 
  !!     Original :   R. Dussin (sept 2010) for ORCA025
  !!
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio
  USE eos

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk , jt , jj , ji                   !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: 
  INTEGER   :: npiglo,npjglo, npk, npt             !: size of the domain
  INTEGER, DIMENSION(1) ::  ipk, &                 !: outptut variables : number of levels,
       &                    id_varout              !: ncdf varid's
  real(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: ztemp, zsal ,&   !: Array to read a layer of data
       &                                         zsigi , &        !: potential density (sig-i)
       &                                         zmask , &        !: 2D mask at current level
       &                                         zwdmp , &        !: damping mask at current level
       &                                         zlat             !: latitude
  REAL(KIND=4),DIMENSION(:),ALLOCATABLE       :: tim , zdep
  REAL(KIND=4)                                :: spval  !: missing value

  CHARACTER(LEN=256) :: cfilet , cfiles, cfilemask='mask.nc', cfileout='mask_dmp.nc' !:
  CHARACTER(LEN=256) :: cdum

  TYPE(variable) , DIMENSION(1) :: typvar         !: structure for attributes

  INTEGER    :: ncout
  INTEGER    :: istatus
  ! default parameters
  REAL(KIND=4)   :: prof=2000.
  REAL(KIND=4)   :: snmax=37.16 , swidth=0.025
  REAL(KIND=4)   :: hmin=1000. , hwidth=100.
  REAL(KIND=4)   :: latmax=-20. , latwidth=2.
  REAL(KIND=4)   :: riri, fifi, loulou

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmaskdmp fileT fileS [prof snmax swidth hmin hwidth latmax latwidth]'
     PRINT *,' default is : cdfmaskdmp fileT fileS 2000. 37.16 0.025 1000. 100. -20. 2. '
     PRINT *,' mask.nc must be in your directory '
     PRINT *,' Output on mask_dmp.nc, variable wdmp'
     STOP
  ENDIF

  IF ( narg > 2 .AND. narg < 9 ) THEN
     PRINT *,'wrong number of arguments'
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  CALL getarg (2, cfiles)
  IF ( narg == 9 ) THEN

    CALL getarg (3, cdum)
    READ(cdum,*) prof
    CALL getarg (4, cdum)
    READ(cdum,*) snmax
    CALL getarg (5, cdum)
    READ(cdum,*) swidth
    CALL getarg (6, cdum)
    READ(cdum,*) hmin
    CALL getarg (7, cdum)
    READ(cdum,*) hwidth
    CALL getarg (8, cdum)
    READ(cdum,*) latmax
    CALL getarg (9, cdum)
    READ(cdum,*) latwidth
 
  ENDIF

  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')
  npt   = getdim (cfilet,'time')

  ipk(:)= npk  ! all variables (input and output are 3D)
  typvar(1)%name='wdmp'
  typvar(1)%missing_value=1.e+20
  typvar(1)%axis='TZYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'npt   =', npt

  ALLOCATE (ztemp(npiglo,npjglo), zsal(npiglo,npjglo), zsigi(npiglo,npjglo) ,zmask(npiglo,npjglo) , zlat(npiglo,npjglo))
  ALLOCATE (zwdmp(npiglo,npjglo))
  ALLOCATE (tim(npt) , zdep(npk))

  ! create output fileset

  ncout =create(cfileout, cfilet, npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfilet,npiglo, npjglo,npk)

  tim=getvar1d(cfilet,'time_counter',npt)
  zdep=getvar1d(cfilet,'deptht',npk)
  zlat(:,:) = getvar(cfilet, 'nav_lat',  1 ,npiglo, npjglo)

  ierr=putvar1d(ncout,tim,npt,'T')

  DO jt = 1, npt
     PRINT *,'time: ',jt
  DO jk = 1, npk
     PRINT *, 'jk = ', jk

     ztemp(:,:)= getvar(cfilet, 'votemper',  jk ,npiglo, npjglo,ktime=jt)
     zsal(:,:) = getvar(cfiles, 'vosaline',  jk ,npiglo, npjglo,ktime=jt)

     zmask(:,:) = getvar(cfilemask, 'tmask',  jk ,npiglo, npjglo)

     zsigi(:,:) = sigmai ( ztemp,zsal,prof,npiglo,npjglo )* zmask(:,:)

     DO jj=1,npjglo
        DO ji=1,npiglo

           riri=tanh((zdep(jk)-hmin)/hwidth)/2. + 0.5
           fifi=tanh((zsigi(ji,jj)-snmax)/swidth)/2. + 0.5
           loulou=tanh(-(zlat(ji,jj)-latmax)/latwidth)/2. + 0.5

           zwdmp(ji,jj)=riri * fifi * loulou

        ENDDO
     ENDDO

     zwdmp(:,:) = zwdmp(:,:) * zmask(:,:)

     ierr = putvar(ncout, id_varout(1) ,zwdmp, jk,npiglo, npjglo,ktime=jt)

  END DO  ! loop to next level
  END DO  ! loop on time

  istatus = closeout(ncout)
END PROGRAM cdfmaskdmp
