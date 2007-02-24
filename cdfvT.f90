PROGRAM cdfvT
  !!-------------------------------------------------------------------
  !!                ***  PROGRAM cdfvT  ***
  !!
  !!  **  Purpose: 
  !!  
  !!  **  Method: Try to avoid 3 d arrays 
  !!
  !! history :
  !!   Original : J.M. Molines (Nov 2004 ) for ORCA025
  !!              J.M. Molines  (apr 2005 ) : use of modules
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
  INTEGER, DIMENSION(4) ::  ipk, id_varout
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: zcumulut, zcumulus  !: Arrays for cumulated values
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: zcumulvt, zcumulvs  !: Arrays for cumulated values
  REAL(KIND=8)                                :: total_time
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: ztemp, zsal ,&       !: Array to read a layer of data
       &                                         zvitu, zvitv,   &
       &                                         zworku, zworkv,   &
       &                                         rmean
  REAL(KIND=4),DIMENSION(1)                   :: timean, tim

  CHARACTER(LEN=80) :: cfilet,cfileu,cfilev ,cfileout='vt.nc', config , ctag !:
  TYPE (variable), DIMENSION(4)     :: typvar     !: structure for attributes
  LOGICAL :: lexist                               !: to inquire existence of files

  INTEGER    :: ncout
  INTEGER    :: istatus

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfvT CONFIG ''list_of_tags'' '
     PRINT *,'        CONFIG is the CONFIG name (eg: ORCA025-G32 ) '
     PRINT *,'        list_of_tags is the list of the time tags (y....m.. d..)'
     PRINT *,'        on which the mean values of UT, US, VT, VS are computes'
     PRINT *,'   Output on vt.nc variables vozout, vozous, vomevt, vomevs  '
     STOP
  ENDIF

  ntags = narg -1   ! first argument is the config name
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, config)
  CALL getarg (2, ctag)
  WRITE(cfilet,'(a,"_",a,"_gridT.nc")') TRIM(config),TRIM(ctag)
  INQUIRE(FILE=cfilet,EXIST=lexist)
  IF ( .NOT. lexist ) THEN
    WRITE(cfilet,'(a,"_",a,"_grid_T.nc")') TRIM(config),TRIM(ctag)
    INQUIRE(FILE=cfilet,EXIST=lexist)
    IF ( .NOT. lexist ) THEN
      PRINT *,' ERROR : missing gridT or even grid_T file '
      STOP
    ENDIF
  ENDIF

  PRINT *,TRIM(cfilet)
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')

  ipk(:)= npk  ! all variables (input and output are 3D)
  ! define output variables
  typvar(1)%name= 'vomevt'
  typvar(2)%name= 'vomevs'
  typvar(3)%name= 'vozout'
  typvar(4)%name= 'vozous'

  typvar(1)%units='m.DegC.s-1'
  typvar(2)%units='m.PSU.s-1'
  typvar(3)%units='m.DegC.s-1'
  typvar(4)%units='m.PSU.s-1'

  typvar%missing_value=0.
  typvar%valid_min= -100.
  typvar%valid_max= 100.

  typvar(1)%long_name='Meridional_VT'
  typvar(2)%long_name='Meridional_VS'
  typvar(3)%long_name='Zonal_UT'
  typvar(4)%long_name='Zonal_US'

  typvar(1)%short_name='vomevt'
  typvar(2)%short_name='vomevs'
  typvar(3)%short_name='vozout'
  typvar(4)%short_name='vozous'

  typvar%online_operation='N/A'
  typvar%axis='TZYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( zcumulut(npiglo,npjglo), zcumulus(npiglo,npjglo) )
  ALLOCATE( zcumulvt(npiglo,npjglo), zcumulvs(npiglo,npjglo) )
  ALLOCATE( zvitu(npiglo,npjglo),zvitv(npiglo,npjglo) )
  ALLOCATE( zworku(npiglo,npjglo),zworkv(npiglo,npjglo) )
  ALLOCATE( ztemp(npiglo,npjglo) ,zsal(npiglo,npjglo) )
  ALLOCATE( rmean(npiglo,npjglo))


  ! create output fileset

  ncout =create(cfileout, cfilet, npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,4, ipk,id_varout )
  ierr= putheadervar(ncout, cfilet,npiglo, npjglo, npk )

  DO jk = 1, npk
     PRINT *,'level ',jk
     zcumulut(:,:) = 0.d0 ;  zcumulvt(:,:) = 0.d0 ; total_time = 0.
     zcumulus(:,:) = 0.d0 ;  zcumulvs(:,:) = 0.d0 

     DO jt = 2, narg
        CALL getarg (jt, ctag)

        WRITE(cfilet,'(a,"_",a,"_gridT.nc")') TRIM(config),TRIM(ctag)
        INQUIRE(FILE=cfilet,EXIST=lexist)
        IF ( .NOT. lexist ) THEN
          WRITE(cfilet,'(a,"_",a,"_grid_T.nc")') TRIM(config),TRIM(ctag)
          INQUIRE(FILE=cfilet,EXIST=lexist)
          IF ( .NOT. lexist ) THEN
             PRINT *,' ERROR : missing gridT or even grid_T file '
             STOP
          ENDIF
        ENDIF

        WRITE(cfileu,'(a,"_",a,"_gridU.nc")') TRIM(config),TRIM(ctag)
        INQUIRE(FILE=cfileu,EXIST=lexist)
        IF ( .NOT. lexist ) THEN
          WRITE(cfileu,'(a,"_",a,"_grid_U.nc")') TRIM(config),TRIM(ctag)
          INQUIRE(FILE=cfileu,EXIST=lexist)
          IF ( .NOT. lexist ) THEN
             PRINT *,' ERROR : missing gridU or even grid_U file '
             STOP
          ENDIF
        ENDIF

        WRITE(cfilev,'(a,"_",a,"_gridV.nc")') TRIM(config),TRIM(ctag)
        INQUIRE(FILE=cfilev,EXIST=lexist)
        IF ( .NOT. lexist ) THEN
          WRITE(cfileu,'(a,"_",a,"_grid_V.nc")') TRIM(config),TRIM(ctag)
          INQUIRE(FILE=cfileu,EXIST=lexist)
          IF ( .NOT. lexist ) THEN
             PRINT *,' ERROR : missing gridV or even grid_V file '
             STOP
          ENDIF
        ENDIF

        IF (jk == 1  )  THEN
           tim=getvar1d(cfilet,'time_counter',1)
           total_time = total_time + tim(1)
        END IF

        zvitu(:,:)= getvar(cfileu, 'vozocrtx' , jk ,npiglo, npjglo )
        zvitv(:,:)= getvar(cfilev, 'vomecrty' , jk ,npiglo, npjglo )
        ztemp(:,:)= getvar(cfilet, 'votemper',  jk ,npiglo, npjglo )
        zsal(:,:) = getvar(cfilet, 'vosaline',  jk ,npiglo, npjglo )

        ! temperature
        zworku(:,:) = 0. ; zworkv(:,:) = 0.
        DO ji=1, npiglo-1
           DO jj = 1, npjglo -1
              zworku(ji,jj) = 0.5 * ( ztemp(ji,jj) + ztemp(ji+1,jj) )  ! temper at Upoint
              zworkv(ji,jj) = 0.5 * ( ztemp(ji,jj) + ztemp(ji,jj+1) )  ! temper at Vpoint
           END DO
        END DO

        zcumulut(:,:) = zcumulut(:,:) + zworku(:,:) * zvitu(:,:)
        zcumulvt(:,:) = zcumulvt(:,:) + zworkv(:,:) * zvitv(:,:)

        ! salinity
        zworku(:,:) = 0. ; zworkv(:,:) = 0.
        DO ji=1, npiglo-1
           DO jj = 1, npjglo -1
              zworku(ji,jj) = 0.5 * ( zsal(ji,jj) + zsal(ji+1,jj) )  ! salinity  at Upoint
              zworkv(ji,jj) = 0.5 * ( zsal(ji,jj) + zsal(ji,jj+1) )  ! salinity  at Vpoint
           END DO
        END DO

        zcumulus(:,:) = zcumulus(:,:) + zworku(:,:) * zvitu(:,:)
        zcumulvs(:,:) = zcumulvs(:,:) + zworkv(:,:) * zvitv(:,:)

     END DO

     ! finish with level jk ; compute mean (assume spval is 0 )
     rmean(:,:) = zcumulvt(:,:)/ntags
     ierr = putvar(ncout, id_varout(1) ,rmean, jk,npiglo, npjglo )

     rmean(:,:) = zcumulvs(:,:)/ntags
     ierr = putvar(ncout, id_varout(2) ,rmean, jk,npiglo, npjglo )

     rmean(:,:) = zcumulut(:,:)/ntags
     ierr = putvar(ncout, id_varout(3) ,rmean, jk,npiglo, npjglo )

     rmean(:,:) = zcumulus(:,:)/ntags
     ierr = putvar(ncout, id_varout(4) ,rmean, jk,npiglo, npjglo )

     IF (jk == 1 )  THEN
        timean(1)= total_time/ntags
        ierr=putvar1d(ncout,timean,1,'T')
     END IF

  END DO  ! loop to next level

  istatus = closeout(ncout)

END PROGRAM cdfvT
