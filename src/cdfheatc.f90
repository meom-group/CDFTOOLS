PROGRAM cdfheatc
  !!======================================================================
  !!                     ***  PROGRAM  cdfheatc  ***
  !!=====================================================================
  !!  ** Purpose : Compute the heat content of the ocean : 1 single value
  !!
  !!  ** Method  : compute the sum ( rho cp T  * e1t *e2t * e3t * tmask )
  !!
  !! History : 2.1  : 03/2006  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: jk, jt              ! dummy loop index
  INTEGER(KIND=4)                           :: ik                  ! working integer
  INTEGER(KIND=4)                           :: ierr                ! working integer
  INTEGER(KIND=4)                           :: iimin=0, iimax=0    ! domain limitation for computation
  INTEGER(KIND=4)                           :: ijmin=0, ijmax=0    ! domain limitation for computation
  INTEGER(KIND=4)                           :: ikmin=0, ikmax=0    ! domain limitation for computation
  INTEGER(KIND=4)                           :: mxloption=0         ! mixed layer option    
  INTEGER(KIND=4)                           :: narg, iargc, ijarg  ! command line 
  INTEGER(KIND=4)                           :: npiglo, npjglo      ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt            ! size of the domain
  INTEGER(KIND=4)                           :: nvpk                ! vertical levels in working variable

  REAL(KIND=4), PARAMETER                   :: pprho0=1020.        ! water density (kg/m3)
  REAL(KIND=4), PARAMETER                   :: ppcp=4000.          ! calorific capacity (J/kg/m3)

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2t                 ! horizontal metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: weight              ! horizontal metrics
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: e3t               ! vertical metric
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: temp                ! temperature
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: tmask             ! tmask
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rmxldep             ! mixed layer depth
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: gdepw               ! depth
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                 ! time counter
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: e31d                ! vertical metrics in case of full step

  REAL(KIND=8)                              :: dvol                ! 3D volume of the ocean
  REAL(KIND=8)                              :: dsum                ! weighted sum 3D
  REAL(KIND=8)                              :: dvol2d              ! volume of a layer
  REAL(KIND=8)                              :: dsum2d              ! weigthed sum per layer
  REAL(KIND=8)                              :: dsurf               ! surface of a layer

  CHARACTER(LEN=256)                        :: cf_tfil             ! input gridT file
  CHARACTER(LEN=256)                        :: cldum               ! dummy character variable

  LOGICAL                                   :: lfull=.FALSE.       ! flag for full step computation
  LOGICAL                                   :: lchk                ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfheatc  T-file ...'
     PRINT *,'    ... [imin imax jmin jmax kmin kmax] [-full] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Computes the heat content in the specified area (Joules)'
     PRINT *,'        A sub-domain can be specified in option.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file : a file with temperature and salinity' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [imin imax jmin jmax kmin kmax] : limit of a sub domain where'
     PRINT *,'                      the heat content will be calculated.'
     PRINT *,'                   - if imin = 0 then ALL i are taken'
     PRINT *,'                   - if jmin = 0 then ALL j are taken'
     PRINT *,'                   - if kmin = 0 then ALL k are taken'
     PRINT *,'       [-full ] : assume full step model output instead of default'
     PRINT *,'                  partial steps.'
     PRINT *,'       [-mxloption ] : pass 1 to compute only in the mixed layer, -1 to exclude'
     PRINT *,'                       it from the calculations '
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       Files ',TRIM(cn_fhgr),', ',TRIM(cn_fzgr),' and ',TRIM(cn_fmsk) 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : to be done ....'
     PRINT *,'       Standard output'
     STOP
  ENDIF

  PRINT *,'I am debugging the correct one'
  ijarg = 1 
  CALL getarg (ijarg, cf_tfil) ; ijarg = ijarg + 1

  lchk = chkfile(cn_fhgr)
  lchk = chkfile(cn_fzgr) .OR. lchk
  lchk = chkfile(cn_fmsk) .OR. lchk
  lchk = chkfile(cf_tfil) .OR. lchk
  IF ( lchk ) STOP ! missing files

  DO WHILE ( ijarg <= narg ) 
     CALL getarg ( ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-full' ) ; lfull = .true.
     CASE ( '-mxloption' ) ;
        CALL getarg ( ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) mxloption
     CASE DEFAULT
        PRINT *,' Reading 6 values : imin imax jmin jmax kmin kmax'
                                                          READ(cldum,*) iimin
        CALL getarg ( ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) iimax
        CALL getarg ( ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmin
        CALL getarg ( ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmax
        CALL getarg ( ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ikmin
        CALL getarg ( ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ikmax
     END SELECT
  END DO

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  IF (iimin /= 0 ) THEN ; npiglo = iimax - iimin + 1;  ELSE ; iimin=1 ; ENDIF
  IF (ijmin /= 0 ) THEN ; npjglo = ijmax - ijmin + 1;  ELSE ; ijmin=1 ; ENDIF
  IF (ikmin /= 0 ) THEN ; npk    = ikmax - ikmin + 1;  ELSE ; ikmin=1 ; ENDIF

  nvpk   = getvdim(cf_tfil,cn_votemper)
  IF (nvpk == 2 ) nvpk = 1
  IF (nvpk == 3 ) nvpk = npk
  IF (ikmax == 0) ikmax = nvpk

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt
  PRINT *, 'nvpk   = ', nvpk

  ! Allocate arrays
  PRINT *, 'Allocate TMASK'
  ALLOCATE ( tmask(npiglo,npjglo,nvpk))
  PRINT *, 'Allocate temp'
  ALLOCATE ( temp (npiglo,npjglo))
  PRINT *, 'Allocate e2t and e3t'
  ALLOCATE ( e2t(npiglo,npjglo), e3t(npiglo,npjglo,nvpk))
  PRINT *, 'Allocate weight'
  ALLOCATE ( weight(npiglo,npjglo))

  IF (mxloption /= 0) THEN
      PRINT *, 'Allocate rmxldep'
      ALLOCATE ( rmxldep(npiglo,npjglo))
  ENDIF

  PRINT *, 'Allocate gdepw'
  ALLOCATE ( gdepw(npk), tim(npt))
  IF ( lfull ) ALLOCATE ( e31d(npk))

  e2t(:,:) = getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin) * getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin)
  gdepw(:) = getvare3(cn_fzgr, cn_gdepw,  npk)
  tim  (:) = getvare3(cf_tfil, cn_vtimec, npt)
  
  IF ( lfull ) e31d(:) = getvare3(cn_fzgr, cn_ve3t, npk)

  DO jt=1,npt
     dvol = 0.d0
     dsum = 0.d0
     PRINT * ,'TIME : ', tim(jt)
     IF (mxloption /= 0) rmxldep(:,:) = getvar(cf_tfil, cn_somxl010, 1, npiglo, npjglo, ktime=jt)

     DO jk = ikmin,ikmax
        ik = jk + ikmin -1
        ! Get velocities v at ik
        temp( :,:)   = getvar(cf_tfil, cn_votemper, ik, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=jt)
        IF ( jt == 1 ) THEN
            PRINT *, 'Read mask'
            tmask(:,:,jk)   = getvar(cn_fmsk, 'tmask', ik, npiglo, npjglo, kimin=iimin, kjmin=ijmin)

            ! get e3t at level ik ( ps...)
            PRINT *, 'Load e3t'
            IF ( lfull ) THEN
               e3t(:,:, jk) = e31d(jk)
            ELSE
               e3t(:,:, jk) = getvar(cn_fzgr, 'e3t_0', ik, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ldiom=.TRUE.)
            ENDIF
        ENDIF
        PRINT *, 'Compute weight'
        weight = e2t * tmask(:,:,jk)

        dsurf  = sum(weight)
        
        SELECT CASE ( mxloption ) 
         CASE ( 1 ) 
            weight(:,:) = MAX ( 0., MIN( e3t(:,:, jk),rmxldep-gdepw(ik) ) ) * weight
         CASE ( -1 )
            weight(:,:) = MIN ( e3t(:,:, jk), MAX( 0.,gdepw(ik)+e3t(:,:,jk)-rmxldep ) ) * weight
         CASE ( 0 )
            weight(:,:)= e3t(:,:, jk) * weight
         END SELECT


        dvol2d = sum(weight)
        dvol   = dvol + dvol2d

        dsum2d = sum(weight * temp)
        dsum   = dsum + dsum2d

        IF (dvol2d /= 0 )THEN
           PRINT *, ' Heat Content  at level ',ik,'(',gdepw(ik),' m) ',pprho0*ppcp*dsum2d, 'surface = ',dsurf/1.e6,' km^2'
        ELSE
           PRINT *, ' No points in the water at level ',ik,'(',gdepw(ik),' m) '
        ENDIF

     END DO
     
     PRINT * ,' Total Heat content        : ', pprho0*ppcp*dsum ,' Joules'
     PRINT * ,' Total Heat content/volume : ', pprho0*ppcp*dsum/dvol ,' Joules/m3 '
  END DO

END PROGRAM cdfheatc
