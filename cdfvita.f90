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
  INTEGER   :: ji,jj,jk, jlev
  INTEGER   :: narg, iargc, ijarg                                  !: 
  INTEGER   :: npiglo,npjglo, npk                                !: size of the domain
  INTEGER   :: nlev, nvar, ik
  INTEGER, DIMENSION(:),ALLOCATABLE ::  ipk, id_varout, nklev
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvar
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: u, v, ua, va, vmod
  REAL(KIND=4) ,DIMENSION(1)                  :: timean
  REAL(KIND=4) ,DIMENSION(:),     ALLOCATABLE :: gdept, gdeptall

  CHARACTER(LEN=256) :: cfileu ,cfilev, cfilew,  cfilet, cfileout='vita.nc'            !: file name
  CHARACTER(LEN=256) :: cdum

  INTEGER  :: ncout
  INTEGER  :: istatus, ierr
  LOGICAL  :: lvertical = .false.

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfvita gridU  gridV  gridT2 [-w gridW ] [-lev level_list]'
     PRINT *,'   Grid T2 is only required for the Tgrid of output field'
     PRINT *,'   if optionnal -w gridW file is given, then the W component '
     PRINT *,'   is also interpolated'
     PRINT *,'   We suggest to give a gridT2 file, which is smaller '
     PRINT *,'  [-lev level_list ] : specify a list of level to be used '
     PRINT *,'       (default option is to use all input levels).'
     PRINT *,'       This option MUST be the last on the command line !!'
     PRINT *,'   Output on vita.nc ,variables sovitua sovitva sovitmod [ sovitwa ]'
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  nlev = 0
  ijarg=1
  DO WHILE ( ijarg <= narg )
     CALL getarg( ijarg, cdum ) ; ijarg=ijarg+1
     SELECT CASE ( cdum )
     CASE ( '-lev' )
        nlev= narg - ijarg + 1
        ALLOCATE (nklev(nlev) )
        DO jlev = 1, nlev
           CALL getarg( ijarg, cdum ) ; ijarg=ijarg+1 ; READ(cdum,* ) nklev(jlev)
        ENDDO
     CASE ( '-w' )
        CALL getarg( ijarg, cfilew ) ; ijarg=ijarg+1
        lvertical=.true.
     CASE DEFAULT
        cfileu=cdum
        CALL getarg( ijarg, cfilev ) ; ijarg=ijarg+1
        CALL getarg( ijarg, cfilet ) ; ijarg=ijarg+1
     END SELECT
  ENDDO

  ! adjust number of variable according to -w option
  IF ( lvertical ) THEN
     nvar = 4 
  ELSE
     nvar = 3
  ENDIF

  ALLOCATE ( ipk(nvar), id_varout(nvar), typvar(nvar) )

  npiglo = getdim (cfileu,'x')
  npjglo = getdim (cfileu,'y')
  npk    = getdim (cfileu,'depth')

  IF ( nlev == 0 ) THEN ! take all levels
     nlev = npk
     ALLOCATE (nklev(nlev) )
     DO jlev = 1, nlev
        nklev(jlev) = jlev
     ENDDO
  ENDIF

  ALLOCATE ( gdept(nlev) )

  ipk(1)      = nlev
  typvar(1)%name='sovitua'
  typvar(1)%units='m/s'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= 0.
  typvar(1)%valid_max= 10000.
  typvar(1)%long_name='Zonal Velocity T point'
  typvar(1)%short_name='sovitua'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TYX'

  ipk(2)      = nlev
  typvar(2)%name='sovitva'
  typvar(2)%units='m/s'
  typvar(2)%missing_value=0.
  typvar(2)%valid_min= 0.
  typvar(2)%valid_max= 10000.
  typvar(2)%long_name='Meridional Velocity T point'
  typvar(2)%short_name='sovitva'
  typvar(2)%online_operation='N/A'
  typvar(2)%axis='TYX'

  ipk(3)      = nlev
  typvar(3)%name='sovitmod'
  typvar(3)%units='m/s'
  typvar(3)%missing_value=0.
  typvar(3)%valid_min= 0.
  typvar(3)%valid_max= 10000.
  typvar(3)%long_name='Velocity module T point'
  typvar(3)%short_name='sovitmod'
  typvar(3)%online_operation='N/A'
  typvar(3)%axis='TYX'

  IF ( lvertical ) THEN
     ipk(4)      = nlev
     typvar(4)%name='sovitwa'
     typvar(4)%units='mm/s'
     typvar(4)%missing_value=0.
     typvar(4)%valid_min= 0.
     typvar(4)%valid_max= 10000.
     typvar(4)%long_name='Vertical Velocity at T point'
     typvar(4)%short_name='sovitwa'
     typvar(4)%online_operation='N/A'
     typvar(4)%axis='TYX'
  ENDIF


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'nlev  =', nlev

  ALLOCATE( u(npiglo,npjglo),  v(npiglo,npjglo) , gdeptall(npk) )
  ALLOCATE( ua(npiglo,npjglo), va(npiglo,npjglo), vmod(npiglo,npjglo) )

  gdeptall(:) = getvar1d(cfilet,'deptht',npk)

  DO jlev = 1, nlev
     ik = nklev(jlev)
     gdept(jlev) = gdeptall(ik)
  ENDDO

  ncout =create(cfileout, cfilet, npiglo, npjglo, nlev)
  ierr= createvar(ncout, typvar, nvar, ipk, id_varout )
  ierr= putheadervar(ncout, cfilet, npiglo, npjglo, nlev, pdep=gdept)

  DO jlev = 1, nlev
     ik = nklev(jlev)
     u(:,:) = getvar(cfileu,'vozocrtx',ik ,npiglo, npjglo)
     v(:,:) = getvar(cfilev,'vomecrty',ik ,npiglo, npjglo)

     ua = 0. ; va = 0. ; ua(:,:) = 0. ; va(:,:)=0. ; vmod(:,:)=0.
     DO ji=2, npiglo
        DO jj=2,npjglo
           ua(ji,jj) = 0.5* (u(ji,jj)+ u(ji-1,jj))
           va(ji,jj) = 0.5* (v(ji,jj)+ v(ji,jj-1))
           vmod(ji,jj) = SQRT( ua(ji,jj)*ua(ji,jj) + va(ji,jj)*va(ji,jj) )
        END DO
     END DO
     ierr=putvar(ncout,id_varout(1), ua, jlev ,npiglo, npjglo)
     ierr=putvar(ncout,id_varout(2), va, jlev ,npiglo, npjglo)
     ierr=putvar(ncout,id_varout(3), vmod, jlev ,npiglo, npjglo)
  END DO

  IF ( lvertical ) THEN
     ! reuse u an v arrays to store Wk and Wk+1
     DO jlev=1, nlev-1
        u(:,:) = getvar(cfilew,'vovecrtz',nklev(jlev)   ,npiglo, npjglo)
        v(:,:) = getvar(cfilew,'vovecrtz',nklev(jlev)+1 ,npiglo, npjglo)
        ua(:,:)=0.5*(u(:,:) + v(:,:))*1000.  ! mm/sec
        ierr=putvar(ncout,id_varout(4), ua, jlev ,npiglo, npjglo)
     END DO

     IF (nlev == npk ) THEN
        ua(:,:)=0.e0
     ELSE
        u(:,:) = getvar(cfilew,'vovecrtz',nklev(nlev)   ,npiglo, npjglo)
        v(:,:) = getvar(cfilew,'vovecrtz',nklev(nlev)+1 ,npiglo, npjglo)
        ua(:,:)=0.5*(u(:,:) + v(:,:))*1000.  ! mm/sec
     ENDIF
     ierr=putvar(ncout,id_varout(4), ua, nlev ,npiglo, npjglo)
  ENDIF

  timean=getvar1d(cfileu,'time_counter',1)
  ierr=putvar1d(ncout,timean,1,'T')
  istatus = closeout(ncout)

END PROGRAM cdfvita
