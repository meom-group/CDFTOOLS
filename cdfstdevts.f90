PROGRAM cdfstdevts
  !!--------------------------------------------------------------------
  !!         ***  PROGRAM cdfstdevts  ***
  !!
  !!  **  Purpose : Compute standard deviation of TS fields
  !!  
  !!  **  Method  : Start from T2 files computed with cdfmoy_sal2_temp2
  !!
  !! history :
  !!     Original : J.M. Molines (nov 2004) for ORCA025
  !!                J.M. Molines (Apr 2005) : use of modules
  !!--------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk, jvar
  INTEGER   :: narg, iargc                                  !: 
  INTEGER   :: npiglo,npjglo, npk                                !: size of the domain
  INTEGER, DIMENSION(2) ::  ipk, id_varout
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: u, u2,  stdev
  REAL(KIND=4) ,DIMENSION(1)                  :: timean

  CHARACTER(LEN=256) :: cfile ,cfile2 ,cfileout='stdevts.nc'            !: file name
  CHARACTER(LEN=256), DIMENSION(2) :: cvar, cvar2

  TYPE(variable), DIMENSION(2)    :: typvar          !: structure for attributes

  INTEGER    :: ncout
  INTEGER    :: istatus, ierr

  !!  Read command line
  narg= iargc()
  IF ( narg /= 2 ) THEN
     PRINT *,' Usage : cdfstdevts ''gridX gridX2'' '
     PRINT *,'   Output on stdevts.nc variable votemper_stdev vosaline_stdev'
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cfile)
  CALL getarg (2, cfile2)

  npiglo= getdim (cfile,'x')
  npjglo= getdim (cfile,'y')
  npk   = getdim (cfile,'depth')

  cvar(1)='votemper'  ; cvar2(1)='votemper_sqd'
  cvar(2)='vosaline'  ; cvar2(2)='vosaline_sqd'

  ipk(1) = npk
  typvar(1)%name= 'votemper_stdev'
  typvar(1)%units='DegC'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= 0.
  typvar(1)%valid_max= 20.
  typvar(1)%long_name='stdev_temperature'
  typvar(1)%short_name='votemper_stdev'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZYX'

  ipk(2) = npk
  typvar(2)%name= 'vosaline_stdev'
  typvar(2)%units='PSU'
  typvar(2)%missing_value=0.
  typvar(2)%valid_min= 0.
  typvar(2)%valid_max= 10.
  typvar(2)%long_name='STDEV_salinity'
  typvar(2)%short_name='vosaline_stdev'
  typvar(2)%online_operation='N/A'
  typvar(2)%axis='TZYX'



  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( u(npiglo,npjglo), u2(npiglo,npjglo) )
  ALLOCATE( stdev(npiglo,npjglo) )

  ncout =create(cfileout, cfile,npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,2, ipk,id_varout )
  ierr= putheadervar(ncout, cfile, npiglo, npjglo, npk)

  DO jvar=1,2
  DO jk = 1, ipk(jvar)
     u(:,:) = getvar(cfile,cvar(jvar),jk, npiglo, npjglo)
     u2(:,:) = getvar(cfile2,cvar2(jvar),jk, npiglo, npjglo)

     stdev(:,:) = 0.
     DO ji=2, npiglo
        DO jj=2,npjglo
           stdev(ji,jj)  =  ((u2(ji,jj)-u(ji,jj)*u(ji,jj)))
        END DO
     END DO
     ierr=putvar(ncout,id_varout(jvar), sqrt(real(stdev)), jk, npiglo, npjglo)
  END DO
  timean=getvar1d(cfile,'time_counter',1)
  END DO
  ierr=putvar1d(ncout,timean,1,'T')
  istatus = closeout(ncout)

END PROGRAM cdfstdevts
