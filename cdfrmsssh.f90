PROGRAM cdfrmsssh
  !!-------------------------------------------------------------------
  !!             ***  PROGRAM cdfrmsssh  ***
  !!
  !!  **  Purpose: Compute RMS SSH
  !!  
  !!  **  Method: Try to avoid 3 d arrays 
  !!
  !! history :
  !!  Original :  J.M. Molines (Nov 2004 ) for ORCA025
  !!              J.M. Molines Apr 2005 : use modules
  !!-------------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk
  INTEGER   :: narg, iargc                          !: 
  INTEGER   :: npiglo,npjglo, npk                   !: size of the domain
  INTEGER, DIMENSION(1) ::  ipk, id_varout
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: u, u2,  rms
  REAL(KIND=4) ,DIMENSION(1)                  :: timean

  CHARACTER(LEN=80) :: cfile ,cfile2 ,cfileout='rms.nc'            !: file name

  TYPE(variable), DIMENSION(1) :: typvar          !: structure for attribute

  INTEGER    :: ncout
  INTEGER    :: istatus, ierr

  !!  Read command line
  narg= iargc()
  IF ( narg /= 2 ) THEN
     PRINT *,' Usage : cdfrmsssh ''gridX gridX2'' '
     PRINT *,'   Output on rms.nc , variable sossheig_rms '
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cfile)
  CALL getarg (2, cfile2)

  npiglo= getdim (cfile,'x')
  npjglo= getdim (cfile,'y')
  npk   = getdim (cfile,'depth')

  ipk(1) = 1
  typvar(1)%name= 'sossheig_rms'
  typvar(1)%units='m'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= 0.
  typvar(1)%valid_max= 100.
  typvar(1)%long_name='RMS_Sea_Surface_height'
  typvar(1)%short_name='sossheig_rms'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( u(npiglo,npjglo), u2(npiglo,npjglo) )
  ALLOCATE( rms(npiglo,npjglo) )

  ncout =create(cfileout, cfile,npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfile, npiglo, npjglo, npk)

  DO jk = 1, ipk(1)
     u(:,:) = getvar(cfile,'sossheig',jk, npiglo, npjglo)
     u2(:,:) = getvar(cfile2,'sossheig_sqd',jk, npiglo, npjglo)

     rms(:,:) = 0.
     DO ji=2, npiglo
        DO jj=2,npjglo
           rms(ji,jj)  =  SQRT((u2(ji,jj)-u(ji,jj)*u(ji,jj)))
        END DO
     END DO
     ierr=putvar(ncout,id_varout(1), rms, jk, npiglo, npjglo)
  END DO
  timean=getvar1d(cfile,'time_counter',1)
  ierr=putvar1d(ncout,timean,1,'T')
  istatus = closeout(ncout)

END PROGRAM cdfrmsssh
