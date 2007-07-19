PROGRAM cdfvsig
  !!-------------------------------------------------------------------
  !!                ***  PROGRAM cdfvsig  ***
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
  USE eos

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk,jt                         !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc , ntags                 !: 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: zcumulusig, zcumulsigu, zcumulu  !: Arrays for cumulated values
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: zcumulvsig, zcumulsigv, zcumulv  !: Arrays for cumulated values
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: zcumulwsig, zcumulsigw, zcumulw  !: Arrays for cumulated values
  REAL(KIND=8)                                :: total_time
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: umask, vmask, wmask !: mask of the velocity points
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: ztemp, zsal         !: Array to read a layer of data
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: ztempu, zsalu       !: Array to read a layer of data
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: ztempup, zsalup     !: Array to read a layer of data
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: ztempv, zsalv       !: Array to read a layer of data
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: ztempw, zsalw  ,&   !: Array to read a layer of data
       &                                         zvitu, zvitv, zvitw , &
       &                                         zworku, zworkv, zworkw , &
       &                                         rmean
  REAL(KIND=4),DIMENSION(1)                   :: timean, tim

  CHARACTER(LEN=80) :: config , ctag  !:
  CHARACTER(LEN=80) :: cfilet,cfileu,cfilev, cfilew 
  CHARACTER(LEN=80) :: cfilmask='mask.nc'
  CHARACTER(LEN=80) :: cfilusig='usig.nc',  cfilvsig='vsig.nc',  cfilwsig='wsig.nc' !:
  INTEGER, DIMENSION(3) ::  ipkusig, id_varoutusig,&
                            ipkvsig, id_varoutvsig,&
                            ipkwsig, id_varoutwsig

  TYPE (variable), DIMENSION(3)     :: typvarusig ,&     !: structure for attributes
                                       typvarvsig ,&     !: structure for attributes
                                       typvarwsig        !: structure for attributes
  LOGICAL :: lexist                               !: to inquire existence of files

  INTEGER    :: ncoutusig, ncoutvsig, ncoutwsig
  INTEGER    :: istatus

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfvsig CONFIG ''list_of_tags'' '
     PRINT *,'        CONFIG is the CONFIG name (eg: ORCA025-G32 ) '
     PRINT *,'        list_of_tags is the list of the time tags (y....m.. d..)'
     PRINT *,'        on which the mean values of Usigma, Vsigma, Wsigma are computes'
     PRINT *,'   Output on usig.nc variables vousig vosigu vozocrtx '
     PRINT *,'   Output on vsig.nc variables vovsig vosigv vomecrty '
     PRINT *,'   Output on wsig.nc variables vowsig vosigw vovecrtz '
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

  ipkusig(:)= npk  ! all variables (input and output are 3D)
  ipkvsig(:)= npk  ! all variables (input and output are 3D)
  ipkwsig(:)= npk  ! all variables (input and output are 3D)

  ! define output variables
  typvarusig(1)%name= 'vousig'    ; typvarvsig(1)%name= 'vovsig'   ;  typvarwsig(1)%name= 'vowsig'
  typvarusig(2)%name= 'vosigu'    ; typvarvsig(2)%name= 'vosigv'   ;  typvarwsig(2)%name= 'vosigw'
  typvarusig(3)%name= 'vozocrtx'  ; typvarvsig(3)%name= 'vomecrty' ;  typvarwsig(3)%name= 'vovecrtz'

  typvarusig(1)%units='kg.m-2.s-1' ; typvarvsig(1)%units='kg.m-2.s-1' ; typvarwsig(1)%units='kg.m-2.s-1'
  typvarusig(2)%units='kg.m-3'     ; typvarvsig(2)%units='kg.m-3'     ; typvarwsig(2)%units='kg.m-3'
  typvarusig(3)%units='m.s-1'      ; typvarvsig(3)%units='m.s-1'      ; typvarwsig(3)%units='m.s-1'

  typvarusig%missing_value=0.      ;  typvarvsig%missing_value=0.     ; typvarwsig%missing_value=0.
  typvarusig%valid_min= -100.      ;  typvarvsig%valid_min= -100.     ; typvarwsig%valid_min= -100.
  typvarusig%valid_max= 100.       ;  typvarvsig%valid_max= 100.      ; typvarwsig%valid_max= 100.

  typvarusig(1)%long_name='Mean U x sigma0' ; typvarvsig(1)%long_name='Mean V x sigma0'  
                                                                               typvarwsig(1)%long_name='Mean W x sigma0'
  typvarusig(2)%long_name='Mean sigma0 at U point' ; typvarvsig(2)%long_name='Mean sigma0 at V point' 
                                                                               typvarwsig(2)%long_name='Mean sigma0 at W point'
  typvarusig(3)%long_name='Mean zonal velocity' ;  typvarvsig(3)%long_name='Mean meridional velocity'  
                                                                               typvarwsig(3)%long_name='Mean vertical velocity'

  typvarusig(1)%short_name= 'vousig' ; typvarvsig(1)%short_name= 'vovsig' ; typvarwsig(1)%short_name= 'vowsig'
  typvarusig(2)%short_name= 'vosigu' ; typvarvsig(2)%short_name= 'vosigv' ; typvarwsig(2)%short_name= 'vosigw'
  typvarusig(3)%short_name= 'vozocrtx' ; typvarvsig(3)%short_name= 'vomecrty' ; typvarwsig(3)%short_name= 'vovecrtz'

  typvarusig%online_operation='N/A' ; typvarvsig%online_operation='N/A';  typvarwsig%online_operation='N/A'
  typvarusig%axis='TZYX'            ; typvarvsig%axis='TZYX'           ;   typvarwsig%axis='TZYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( zcumulusig(npiglo,npjglo), zcumulsigu(npiglo,npjglo), zcumulu(npiglo,npjglo) )
  ALLOCATE( zcumulvsig(npiglo,npjglo), zcumulsigv(npiglo,npjglo), zcumulv(npiglo,npjglo) )
  ALLOCATE( zcumulwsig(npiglo,npjglo), zcumulsigw(npiglo,npjglo), zcumulw(npiglo,npjglo) )
  ALLOCATE( zvitu(npiglo,npjglo),zvitv(npiglo,npjglo),zvitw(npiglo,npjglo) )
  ALLOCATE( zworku(npiglo,npjglo),zworkv(npiglo,npjglo),zworkw(npiglo,npjglo) )
  ALLOCATE( ztemp(npiglo,npjglo) ,zsal(npiglo,npjglo) )
  ALLOCATE( ztempup(npiglo,npjglo) ,zsalup(npiglo,npjglo) )
  ALLOCATE( ztempu(npiglo,npjglo) ,zsalu(npiglo,npjglo) )
  ALLOCATE( ztempw(npiglo,npjglo) ,zsalw(npiglo,npjglo) )
  ALLOCATE( rmean(npiglo,npjglo))
  ALLOCATE( umask(npiglo,npjglo), vmask(npiglo,npjglo), wmask(npiglo,npjglo) )


  ! create output fileset s ! JMM BUG : cfileu v w not defined here !

  ncoutusig =create(cfilusig, cfileu, npiglo,npjglo,npk)
  ncoutvsig =create(cfilvsig, cfilev, npiglo,npjglo,npk)
  ncoutwsig =create(cfilwsig, cfilew, npiglo,npjglo,npk)

  ierr= createvar(ncoutusig ,typvarusig,3, ipkusig,id_varoutusig )
  ierr= createvar(ncoutvsig ,typvarvsig,3, ipkvsig,id_varoutvsig )
  ierr= createvar(ncoutwsig ,typvarwsig,3, ipkwsig,id_varoutwsig )

  ierr= putheadervar(ncoutusig, cfileu,npiglo, npjglo, npk )
  ierr= putheadervar(ncoutvsig, cfilev,npiglo, npjglo, npk )
  ierr= putheadervar(ncoutwsig, cfilew,npiglo, npjglo, npk )

  DO jk = 1, npk
     PRINT *,'level ',jk
     total_time = 0.
     zcumulusig(:,:) = 0.d0 ;  zcumulvsig(:,:) = 0.d0 ; zcumulwsig(:,:) = 0.d0 
     zcumulsigu(:,:) = 0.d0 ;  zcumulsigv(:,:) = 0.d0 ; zcumulsigw(:,:) = 0.d0
     zcumulu(:,:) = 0.d0    ;  zcumulv(:,:) = 0.d0    ; zcumulw(:,:) = 0.d0

     umask(:,:)= getvar(cfilmask, 'umask' , jk ,npiglo, npjglo )
     vmask(:,:)= getvar(cfilmask, 'vmask' , jk ,npiglo, npjglo )
     wmask(:,:)= getvar(cfilmask, 'tmask' , jk ,npiglo, npjglo )

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
          WRITE(cfilev,'(a,"_",a,"_grid_V.nc")') TRIM(config),TRIM(ctag)
          INQUIRE(FILE=cfilev,EXIST=lexist)
          IF ( .NOT. lexist ) THEN
             PRINT *,' ERROR : missing gridV or even grid_V file '
             STOP
          ENDIF
        ENDIF

        WRITE(cfilew,'(a,"_",a,"_gridW.nc")') TRIM(config),TRIM(ctag)
        INQUIRE(FILE=cfilew,EXIST=lexist)
        IF ( .NOT. lexist ) THEN
          WRITE(cfilew,'(a,"_",a,"_grid_W.nc")') TRIM(config),TRIM(ctag)
          INQUIRE(FILE=cfilew,EXIST=lexist)
          IF ( .NOT. lexist ) THEN
             PRINT *,' ERROR : missing gridW or even grid_W file '
             STOP
          ENDIF
        ENDIF

        IF (jk == 1  )  THEN
           tim=getvar1d(cfilet,'time_counter',1)
           total_time = total_time + tim(1)
        END IF

        zvitu(:,:)= getvar(cfileu, 'vozocrtx' , jk ,npiglo, npjglo )
        zvitv(:,:)= getvar(cfilev, 'vomecrty' , jk ,npiglo, npjglo )
        zvitw(:,:)= getvar(cfilew, 'vovecrtz' , jk ,npiglo, npjglo )
        ztemp(:,:)= getvar(cfilet, 'votemper',  jk ,npiglo, npjglo )
        zsal(:,:) = getvar(cfilet, 'vosaline',  jk ,npiglo, npjglo )

        ! density horizontal flux
        zworku(:,:) = 0. ; zworkv(:,:) = 0.
        DO ji=1, npiglo-1
           DO jj = 1, npjglo -1
              ztempu(ji,jj) = 0.5 * ( ztemp(ji,jj) + ztemp(ji+1,jj) )  ! temper at Upoint
              zsalu(ji,jj)  = 0.5 * ( zsal(ji,jj)  +  zsal(ji+1,jj) )  ! temper at Upoint
              ztempv(ji,jj) = 0.5 * ( ztemp(ji,jj) + ztemp(ji,jj+1) )  ! temper at Vpoint
              zsalv(ji,jj)  = 0.5 * ( zsal(ji,jj)  +  zsal(ji,jj+1) )  ! temper at Upoint
           END DO
        END DO
        zworku(:,:)=sigma0(ztempu, zsalu,npiglo, npjglo) * umask(:,:)
        zworkv(:,:)=sigma0(ztempv, zsalv,npiglo, npjglo) * vmask(:,:)

        zcumulusig(:,:) = zcumulusig(:,:) + zworku(:,:) * zvitu(:,:)
        zcumulvsig(:,:) = zcumulvsig(:,:) + zworkv(:,:) * zvitv(:,:)
        zcumulsigu(:,:) = zcumulsigu(:,:) + zworku(:,:) 
        zcumulsigv(:,:) = zcumulsigv(:,:) + zworkv(:,:)
        zcumulu(:,:)    = zcumulu(:,:)    +  zvitu(:,:)
        zcumulv(:,:)    = zcumulv(:,:)    +  zvitv(:,:)
        
        IF (jk > 1 ) THEN !  compute wsig
          ztempw=0.5*(ztempup + ztemp)
          zsalw=0.5*(zsalup + zsal)
          zworkw=sigma0(ztempw, zsalw, npiglo,npjglo) * wmask (:,:) ! yes w mask is ok from up to down
          zcumulwsig(:,:)=zcumulwsig(:,:) + zworkw(:,:) * zvitw(:,:)
          zcumulsigw(:,:)=zcumulsigw(:,:) + zworkw(:,:)
          zcumulw(:,:)    = zcumulw(:,:)  + zvitw(:,:)
        ENDIF
        ! save upper T and S for next jk vertical interp at w point
        ztempup=ztemp
        zsalup=zsal

     END DO  ! time loop

     ! finish with level jk ; compute mean (assume spval is 0 )
     rmean(:,:) = zcumulusig(:,:)/ntags ;  ierr = putvar(ncoutusig, id_varoutusig(1) ,rmean, jk,npiglo, npjglo )
     rmean(:,:) = zcumulsigu(:,:)/ntags ;  ierr = putvar(ncoutusig, id_varoutusig(2) ,rmean, jk,npiglo, npjglo )
     rmean(:,:) = zcumulu(:,:)/ntags    ;  ierr = putvar(ncoutusig, id_varoutusig(3) ,rmean, jk,npiglo, npjglo )

     rmean(:,:) = zcumulvsig(:,:)/ntags ;  ierr = putvar(ncoutvsig, id_varoutvsig(1) ,rmean, jk,npiglo, npjglo )
     rmean(:,:) = zcumulsigv(:,:)/ntags ;  ierr = putvar(ncoutvsig, id_varoutvsig(2) ,rmean, jk,npiglo, npjglo )
     rmean(:,:) = zcumulv(:,:)/ntags    ;  ierr = putvar(ncoutvsig, id_varoutvsig(3) ,rmean, jk,npiglo, npjglo )

     rmean(:,:) = zcumulwsig(:,:)/ntags ;  ierr = putvar(ncoutwsig, id_varoutwsig(1) ,rmean, jk,npiglo, npjglo )
     rmean(:,:) = zcumulsigw(:,:)/ntags ;  ierr = putvar(ncoutwsig, id_varoutwsig(2) ,rmean, jk,npiglo, npjglo )
     rmean(:,:) = zcumulw(:,:)/ntags    ;  ierr = putvar(ncoutwsig, id_varoutwsig(3) ,rmean, jk,npiglo, npjglo )

     IF (jk == 1 )  THEN
        timean(1)= total_time/ntags
        ierr=putvar1d(ncoutusig,timean,1,'T')
        ierr=putvar1d(ncoutvsig,timean,1,'T')
        ierr=putvar1d(ncoutwsig,timean,1,'T')
     END IF

  END DO  ! loop to next level

  istatus = closeout(ncoutusig)
  istatus = closeout(ncoutvsig)
  istatus = closeout(ncoutwsig)

END PROGRAM cdfvsig
