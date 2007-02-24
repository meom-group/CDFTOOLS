PROGRAM cdfeke
  !!-------------------------------------------------------------------
  !!              PROGRAM CDFEKE
  !!              **************
  !!
  !!  **  Purpose: Compute EKE from mean files :
  !!                mean gridU , MS gridU mean gridV  MS gridV
  !!  
  !!  **  Method: Try to avoid 3 d arrays 
  !!
  !! history:
  !!    Original:  J.M. Molines (Nov 2004 ) for ORCA025
  !!               J.M. Molines (Apr 2005) : use of modules
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !!
  !! * Modules used
  USE cdfio 

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk
  INTEGER   :: narg, iargc                                  !: 
  INTEGER   :: npiglo,npjglo, npk                                !: size of the domain
  INTEGER, DIMENSION(1) ::  ipk, id_varout
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: u, v, u2, v2,  eke
  REAL(KIND=4)                                :: ua, va
  REAL(KIND=4) ,DIMENSION(1)                  :: timean

  CHARACTER(LEN=80) :: cfileu ,cfileu2,cfilev, cfilev2, cfilet, cfileout='eke.nc'            !: file name
  TYPE(variable), DIMENSION(1) :: typvar

  INTEGER    :: ncout
  INTEGER    :: istatus, ierr

  !!  Read command line
  narg= iargc()
  IF ( narg /= 5 ) THEN
     PRINT *,' Usage : cdfeke ''gridU gridU2 gridV gridV2 gridT2'' '
     PRINT *,'   Grid T2 is only required for the Tgrid of output field'
     PRINT *,'   We suggest to give a gridT2 file, which is smaller '
     PRINT *,'   Output on eke.nc ,variable voeke'
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cfileu)
  CALL getarg (2, cfileu2)
  CALL getarg (3, cfilev)
  CALL getarg (4, cfilev2)
  CALL getarg (5, cfilet)

  npiglo = getdim (cfileu,'x')
  npjglo = getdim (cfileu,'y')
  npk    = getdim (cfileu,'depth')

  ipk(1)      = npk
  typvar(1)%name='voeke'
  typvar(1)%units='m2/s2'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= 0.
  typvar(1)%valid_max= 10000.
  typvar(1)%long_name='Eddy_Kinetic_Energy'
  typvar(1)%short_name='voeke'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( u(npiglo,npjglo), u2(npiglo,npjglo), v(npiglo,npjglo) ,v2(npiglo,npjglo) )
  ALLOCATE( eke(npiglo,npjglo) )

  ncout =create(cfileout, cfilet,npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfilet, npiglo, npjglo,npk)

  DO jk = 1, npk
    u(:,:) = getvar(cfileu,'vozocrtx',jk ,npiglo, npjglo)
    v(:,:) = getvar(cfilev,'vomecrty',jk ,npiglo, npjglo)
    u2(:,:) = getvar(cfileu2,'vozocrtx_sqd',jk ,npiglo, npjglo)
    v2(:,:) = getvar(cfilev2,'vomecrty_sqd',jk ,npiglo, npjglo)

    ua = 0. ; va = 0. ; eke(:,:) = 0.
    DO ji=2, npiglo
      DO jj=2,npjglo
        ua = 0.5* ((u2(ji,jj)-u(ji,jj)*u(ji,jj))+ (u2(ji-1,jj)-u(ji-1,jj)*u(ji-1,jj)))
        va = 0.5* ((v2(ji,jj)-v(ji,jj)*v(ji,jj))+ (v2(ji,jj-1)-v(ji,jj-1)*v(ji,jj-1)))
        eke(ji,jj) = 0.5 * ( ua + va )
      END DO
    END DO
    ierr=putvar(ncout,id_varout(1), eke, jk ,npiglo, npjglo)
  END DO
    timean=getvar1d(cfileu,'time_counter',1)
    ierr=putvar1d(ncout,timean,1,'T')
    istatus = closeout(ncout)

END PROGRAM cdfeke
