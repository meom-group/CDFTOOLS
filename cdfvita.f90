PROGRAM cdfvita
  !!-------------------------------------------------------------------
  !!              PROGRAM CDFVITA
  !!              **************
  !!
  !!  **  Purpose: Compute surface velocity on t grid
  !!                 gridU ,  gridV   gridT (reference)
  !!  
  !!  **  Method: Try to avoid 3 d arrays 
  !!
  !! history:
  !!    Original:  J.M. Molines (Nov 2006 ) for ORCA025
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
  INTEGER, DIMENSION(3) ::  ipk, id_varout
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: u, v, ua, va, vmod
  REAL(KIND=4) ,DIMENSION(1)                  :: timean

  CHARACTER(LEN=80) :: cfileu ,cfilev,  cfilet, cfileout='vita.nc'            !: file name
  TYPE(variable), DIMENSION(3) :: typvar

  INTEGER    :: ncout
  INTEGER    :: istatus, ierr

  !!  Read command line
  narg= iargc()
  IF ( narg /= 3 ) THEN
     PRINT *,' Usage : cdfvita ''gridU  gridV  gridT2'' '
     PRINT *,'   Grid T2 is only required for the Tgrid of output field'
     PRINT *,'   We suggest to give a gridT2 file, which is smaller '
     PRINT *,'   Output on vita.nc ,variables sovitua sovitva sovitmod'
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cfileu)
  CALL getarg (2, cfilev)
  CALL getarg (3, cfilet)

  npiglo = getdim (cfileu,'x')
  npjglo = getdim (cfileu,'y')
  npk    = getdim (cfileu,'depth')

  ipk(1)      = npk
  typvar(1)%name='sovitua'
  typvar(1)%units='m/s'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= 0.
  typvar(1)%valid_max= 10000.
  typvar(1)%long_name='Surface Zonal Velocity T point'
  typvar(1)%short_name='sovitua'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TYX'

  ipk(2)      = npk
  typvar(2)%name='sovitva'
  typvar(2)%units='m/s'
  typvar(2)%missing_value=0.
  typvar(2)%valid_min= 0.
  typvar(2)%valid_max= 10000.
  typvar(2)%long_name='Surface Meridional Velocity T point'
  typvar(2)%short_name='sovitva'
  typvar(2)%online_operation='N/A'
  typvar(2)%axis='TYX'

  ipk(3)      = npk
  typvar(3)%name='sovitmod'
  typvar(3)%units='m/s'
  typvar(3)%missing_value=0.
  typvar(3)%valid_min= 0.
  typvar(3)%valid_max= 10000.
  typvar(3)%long_name='Surface  Velocity module T point'
  typvar(3)%short_name='sovitmod'
  typvar(3)%online_operation='N/A'
  typvar(3)%axis='TYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( u(npiglo,npjglo),  v(npiglo,npjglo)  )
  ALLOCATE( ua(npiglo,npjglo), va(npiglo,npjglo), vmod(npiglo,npjglo) )

  ncout =create(cfileout, cfilet,npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,3, ipk,id_varout )
  ierr= putheadervar(ncout, cfilet, npiglo, npjglo,npk)

  DO jk = 1, npk
    u(:,:) = getvar(cfileu,'vozocrtx',jk ,npiglo, npjglo)
    v(:,:) = getvar(cfilev,'vomecrty',jk ,npiglo, npjglo)

    ua = 0. ; va = 0. ; ua(:,:) = 0. ; va(:,:)=0. ; vmod(:,:)=0.
    DO ji=2, npiglo
      DO jj=2,npjglo
        ua(ji,jj) = 0.5* (u(ji,jj)+ u(ji-1,jj))
        va(ji,jj) = 0.5* (v(ji,jj)+ v(ji,jj-1))
        vmod(ji,jj) = SQRT( ua(ji,jj)*ua(ji,jj) + va(ji,jj)*va(ji,jj) )
      END DO
    END DO
    ierr=putvar(ncout,id_varout(1), ua, jk ,npiglo, npjglo)
    ierr=putvar(ncout,id_varout(2), va, jk ,npiglo, npjglo)
    ierr=putvar(ncout,id_varout(3), vmod, jk ,npiglo, npjglo)
  END DO
    timean=getvar1d(cfileu,'time_counter',1)
    ierr=putvar1d(ncout,timean,1,'T')
    istatus = closeout(ncout)

END PROGRAM cdfvita
