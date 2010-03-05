PROGRAM cdfspice
  !!---------------------------------------------------------------------------------
  !!             ***  PROGRAM cdfspice  ***
  !!
  !!  **  Purpose: Compute spiciness 3D field from gridT file
  !!               Store the results on a 'similar' cdf file.
  !!  
  !!  **  Method: Try to avoid 3 d arrays 
  !!              Following Flament (2002) "A state variable for characterizing water
  !!              masses and their diffusive stability: spiciness."
  !!              Progress in Oceanography Volume 54, 2002, Pages 493-501.   
  !!
  !!  **  Definition:  spiciness = sum(i=0,5)[sum(j=0,4)[b(i,j)*theta^i*(s-35)^j]] 
  !!                   with:  b     -> coefficients
  !!                          theta -> potential temperature
  !!                          s     -> salinity 
  !!
  !!  **  Example: spice(15,33) = 0.54458641375             
  !!
  !! history: 
  !!     Original :   C.O. Dufour (Mar 2010) 
  !!----------------------------------------------------------------------------------
  !!----------------------------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk, jt, ji, jj                      !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: 
  INTEGER   :: npiglo,npjglo, npk, npt             !: size of the domain
  INTEGER, DIMENSION(1) ::  ipk, &                 !: outptut variables : number of levels,
       &                    id_varout              !: ncdf varid's
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: ztemp, zsal ,&             !: Array to read a layer of data
       &                                         ztempt, zsalt, zsalref ,&  !: temporary arrays
       &                                         zspi , &                   !: potential density (sig-0)
       &                                         zmask                      !: 2D mask at current level

  REAL(KIND=8) , DIMENSION (6,5) ::              beta             !: coefficients of spiciness formula
  REAL(KIND=4),DIMENSION(:),ALLOCATABLE   ::  tim

  CHARACTER(LEN=256) :: cfilet ,cfileout='spice.nc' !:

  TYPE(variable) , DIMENSION(1) :: typvar         !: structure for attributes

  INTEGER    :: ncout
  INTEGER    :: istatus

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfspice  gridT '
     PRINT *,' Output on spice.nc, variable vospice'
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')
  npt   = getdim (cfilet,'time')

  ipk(:)= npk  ! all variables (input and output are 3D)
  typvar(1)%name= 'vospice'
  typvar(1)%units='kg/m3'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= -300.
  typvar(1)%valid_max= 300.
  typvar(1)%long_name='spiciness'
  typvar(1)%short_name='vospice'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'npt   =', npt

  ALLOCATE (ztemp(npiglo,npjglo), zsal(npiglo,npjglo), zspi(npiglo,npjglo) ,zmask(npiglo,npjglo))
  ALLOCATE (ztempt(npiglo,npjglo), zsalt(npiglo,npjglo), zsalref(npiglo,npjglo))
  ALLOCATE (tim(npt))

  ! create output fileset

  ncout =create(cfileout, cfilet, npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfilet,npiglo, npjglo,npk)
  tim=getvar1d(cfilet,'time_counter',npt)
  ierr=putvar1d(ncout,tim,npt,'T')

  ! Define coefficients to compute spiciness
  beta(1,1) = 0           ; beta(1,2) = 7.7442e-01  ; beta(1,3) = -5.85e-03   ; beta(1,4) = -9.84e-04   ; beta(1,5) = -2.06e-04
  beta(2,1) = 5.1655e-02  ; beta(2,2) = 2.034e-03   ; beta(2,3) = -2.742e-04  ; beta(2,4) = -8.5e-06    ; beta(2,5) = 1.36e-05
  beta(3,1) = 6.64783e-03 ; beta(3,2) = -2.4681e-04 ; beta(3,3) = -1.428e-05  ; beta(3,4) = 3.337e-05   ; beta(3,5) = 7.894e-06
  beta(4,1) = -5.4023e-05 ; beta(4,2) = 7.326e-06   ; beta(4,3) = 7.0036e-06  ; beta(4,4) = -3.0412e-06 ; beta(4,5) = -1.0853e-06
  beta(5,1) = 3.949e-07   ; beta(5,2) = -3.029e-08  ; beta(5,3) = -3.8209e-07 ; beta(5,4) = 1.0012e-07  ; beta(5,5) = 4.7133e-08
  beta(6,1) = -6.36e-10   ; beta(6,2) = -1.309e-09  ; beta(6,3) = 6.048e-09   ; beta(6,4) = -1.1409e-09 ; beta(6,5) = -6.676e-10

  ! Compute spiciness
  DO jt=1,npt
    PRINT *,' TIME = ', jt, tim(jt)/86400.,' days'
  DO jk = 1, npk
     zmask(:,:)=1.

     ztemp(:,:)= getvar(cfilet, 'votemper',  jk ,npiglo, npjglo,ktime=jt)
     zsal(:,:) = getvar(cfilet, 'vosaline',  jk ,npiglo, npjglo,ktime=jt)

     WHERE(zsal == 0 ) zmask = 0

     ! spiciness at time jt, at level jk  
     zspi(:,:) = 0
     zsalref(:,:) = zsal(:,:) - 35.
     ztempt(:,:) = 1.
     DO ji=1,6
       zsalt(:,:) = 1.
       DO jj=1,5
         zspi(:,:) = zspi(:,:) + beta(ji,jj)*ztempt(:,:)*zsalt(:,:)
         zsalt(:,:) = zsalt(:,:)*zsalref(:,:)
       END DO
       ztempt(:,:) = ztempt(:,:)*ztemp(:,:)     
     END DO

     ierr = putvar(ncout, id_varout(1) ,REAL(zspi), jk,npiglo, npjglo,ktime=jt)

  END DO  ! loop to next level
  END DO  ! next time frame

  istatus = closeout(ncout)
END PROGRAM cdfspice
