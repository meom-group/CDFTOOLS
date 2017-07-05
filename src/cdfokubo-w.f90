  PROGRAM cdfokubow
  !!---------------------------------------------------------------------------
  !!         ***  PROGRAM  cdfokubow ***
  !!
  !!  **  Purpose: Compute the okubow weiss parameter on F-points for given gridU gridV files and variables (like cdfcurl routine)
  !!
  !! history :
  !!   Original :  B. Djath (August 2012)
  !!---------------------------------------------------------------------
  !!  $Rev: 256 $
  !!  $Date: 2012-08-31 19:49:27 +0200 (ven. 31 aout 2012) $
  !! 
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! Copyright (c) 2012, B. Djath
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jt                            ! dummy loop index
  INTEGER(KIND=4)                           :: ilev                                  ! level to be processed
  INTEGER(KIND=4)                           :: npiglo, npjglo                        ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt                              ! size of the domain
  INTEGER(KIND=4)                           :: narg, iargc                           ! browse command line
  INTEGER(KIND=4)                           :: ncout, ierr                           ! browse command line
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout                        ! output variable properties

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2v, e1u, e1f, e2f, e1t, e2t          ! horizontql metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: un, vn                                ! velocity field
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zun, zvn                              ! working arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: okubow, fmask, tmask                  ! curl and fmask
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rotn, cisah1, cisah2t, cisah2         ! curl and fmask
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                                   ! time counter

  CHARACTER(LEN=256)                        :: cf_ufil, cf_vfil                      ! file names
  CHARACTER(LEN=256)                        :: cf_out = 'okubow.nc'                  ! output file name
  CHARACTER(LEN=256)                        :: cv_u, cv_v                            ! variable names
  CHARACTER(LEN=256)                        :: cldum                                 ! dummy string

  TYPE (variable), DIMENSION(1)             :: stypvar                               ! structure for attibutes

  LOGICAL                                   :: lforcing = .FALSE.                    ! forcing flag
  LOGICAL                                   :: lchk     = .FALSE.                    ! flag for missing files
  LOGICAL                                   :: lperio   = .FALSE.                    ! flag for E-W periodicity
  !!----------------------------------------------------------------------------------
  CALL ReadCdfNames() 

  narg = iargc()
  IF ( narg /= 5 ) THEN
     PRINT *,' usage : cdfokubow U-file V-file U-var V-var lev'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute Okubo-Weiss parameter of a vector field, at a specified level.'  
     PRINT *,'       If level is specified as 0, assume that the input files are'
     PRINT *,'       forcing files, presumably on A-grid. In this latter case, the'
     PRINT *,'       vector field is interpolated on the C-grid. In any case, the'
     PRINT *,'       curl is computed on the F-point.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       U-file : zonal component of the vector field.'
     PRINT *,'       V-file : meridional component of the vector field.'
     PRINT *,'       U-var  : zonal component variable name'
     PRINT *,'       V-var  : meridional component variable name.'
     PRINT *,'       lev    : level to be processed. If set to 0, assume forcing file '
     PRINT *,'                in input.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ', TRIM(cn_fhgr),' and ', TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : sokubow (s^-2)'
     STOP
  ENDIF

  CALL getarg(1, cf_ufil)
  CALL getarg(2, cf_vfil)
  CALL getarg(3, cv_u   )
  CALL getarg(4, cv_v   )
  CALL getarg(5, cldum  ) ;  READ(cldum,*) ilev

  lchk = chkfile(cn_fhgr ) .OR. lchk
  lchk = chkfile(cn_fmsk ) .OR. lchk
  lchk = chkfile(cf_ufil ) .OR. lchk
  lchk = chkfile(cf_vfil ) .OR. lchk
  IF ( lchk ) STOP 99 ! missing files

  ! define new variables for output
  stypvar(1)%cname             = 'sokubow'
  stypvar(1)%cunits            = 's-2'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = -1000.
  stypvar(1)%valid_max         =  1000.
  stypvar(1)%clong_name        = 'Okubo_Weiss_param (okubow)'
  stypvar(1)%cshort_name       = 'sokubow'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TYX'

  ipk(1) = 1  !  2D

  npiglo = getdim(cf_ufil,cn_x)
  npjglo = getdim(cf_ufil,cn_y)
  npk    = getdim(cf_ufil,cn_z)
  npt    = getdim(cf_ufil,cn_t) 

  PRINT *, 'npiglo = ',npiglo
  PRINT *, 'npjglo = ',npjglo
  PRINT *, 'npk    = ',npk
  PRINT *, 'npt    = ',npt
  PRINT *, 'ilev   = ',ilev

  !test if lev exists
  IF ( (npk==0) .AND. (ilev > 0) ) THEN
     PRINT *, 'Problem : npk = 0 and lev > 0 STOP'
     STOP 99
  END IF

  ! if forcing field 
  IF ( ilev==0 .AND. npk==0 ) THEN
     lforcing=.true.
     npk = 1 ; ilev=1
     PRINT *, 'npk =0, assume 1'
  END IF

  IF ( npt==0 ) THEN
     PRINT *, 'npt=0, assume 1'
     npt=1
  END IF
  ! check files and determines if the curl will be 2D of 3D
  ! ????????????

  ! Allocate the memory
  ALLOCATE ( e1u(npiglo,npjglo) , e1f(npiglo,npjglo) )
  ALLOCATE ( e2v(npiglo,npjglo) , e2f(npiglo,npjglo) )
  ALLOCATE ( e1t(npiglo,npjglo) , e2t(npiglo,npjglo) )
  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo)  )
  ALLOCATE ( zun(npiglo,npjglo) , zvn(npiglo,npjglo) )
  ALLOCATE ( cisah1(npiglo,npjglo) , cisah2(npiglo,npjglo) )
  ALLOCATE ( cisah2t(npiglo,npjglo) , tmask(npiglo,npjglo) )
  ALLOCATE ( okubow(npiglo,npjglo) , fmask(npiglo,npjglo) )
  ALLOCATE ( rotn(npiglo,npjglo) , tim(npt) )

  e1u =  getvar(cn_fhgr, cn_ve1u, 1, npiglo, npjglo)
  e1f =  getvar(cn_fhgr, cn_ve1f, 1, npiglo, npjglo)
  e2v =  getvar(cn_fhgr, cn_ve2v, 1, npiglo, npjglo)
  e2f =  getvar(cn_fhgr, cn_ve2f, 1, npiglo, npjglo)
  e1t =  getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo)
  e2t =  getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo)

  ! use zun and zvn to store f latitude and longitude for output
  zun = getvar(cn_fhgr, cn_glamf, 1, npiglo, npjglo)
  zvn = getvar(cn_fhgr, cn_gphif, 1, npiglo, npjglo)

  ! look for  E-W periodicity
  IF ( zun(1,1) == zun(npiglo-1,1) ) lperio = .TRUE.

  ! create output fileset
  ncout = create      (cf_out, cf_ufil, npiglo, npjglo, 0                           )
  ierr  = createvar   (ncout , stypvar, 1,      ipk,    id_varout                   )
  ierr  = putheadervar(ncout,  cf_ufil, npiglo, npjglo, 0, pnavlon=zun, pnavlat=zvn )

  tim  = getvar1d(cf_ufil, cn_vtimec, npt      )
  ierr = putvar1d(ncout,   tim,       npt,  'T')
  
  DO jt=1,npt
     IF (MOD(jt,100)==0 ) PRINT *, jt,'/',npt
        ! if files are forcing fields
        zun(:,:) =  getvar(cf_ufil, cv_u, ilev ,npiglo,npjglo, ktime=jt)
        zvn(:,:) =  getvar(cf_vfil, cv_v, ilev ,npiglo,npjglo, ktime=jt)
        tmask(:,:) = getvar(cn_fmsk, 'tmask', ilev , npiglo, npjglo)

     IF ( lforcing ) THEN ! for forcing file u and v are on the A grid
        DO ji=1, npiglo-1
           un(ji,:) = 0.5*(zun(ji,:) + zun(ji+1,:))
        END DO
        !
        DO jj=1, npjglo-1
           vn(:,jj) = 0.5*(zvn(:,jj) + zvn(:,jj+1))
        END DO
        ! end compute u and v on U and V point
     ELSE
       un(:,:) = zun(:,:)
       vn(:,:) = zvn(:,:)
     END IF

     ! compute the mask
     IF ( jt==1 ) THEN
        DO jj = 1, npjglo - 1
           DO ji = 1, npiglo - 1
              fmask(ji,jj)=0.
              fmask(ji,jj)= un(ji,jj)*un(ji,jj+1) * vn(ji,jj)*vn(ji+1,jj)
              IF (fmask(ji,jj) /= 0.) fmask(ji,jj)=1.
           ENDDO
        ENDDO
     END IF

      rotn(:,:) = 0. ; cisah1(:,:) = 0. ; cisah2t(:,:) = 0. ; cisah2(:,:) = 0. ;okubow(:,:) = 0.
     DO jj = 1, npjglo -1 
        DO ji = 1, npiglo -1   ! vector opt.
           rotn(ji,jj) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ) - e2v(ji,jj) * vn(ji,jj)    &
                &         - e1u(ji  ,jj+1) * un(ji  ,jj+1) + e1u(ji,jj) * un(ji,jj)  ) &
                &         * fmask(ji,jj) / ( e1f(ji,jj) * e2f(ji,jj) )       ! quantity on f grid

           cisah1(ji,jj) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ) - e2v(ji,jj) * vn(ji,jj)    &
                &         + e1u(ji  ,jj+1) * un(ji  ,jj+1) - e1u(ji,jj) * un(ji,jj)  )   &
                &         * fmask(ji,jj) / ( e1f(ji,jj) * e2f(ji,jj) )       ! quantity on f grid

           cisah2t(ji,jj) = (  e1u(ji+1,jj  ) * un(ji+1,jj  ) - e1u(ji,jj) * un(ji,jj)    &
                &         - e2v(ji  ,jj+1) * vn(ji  ,jj+1) + e2v(ji,jj) * vn(ji,jj)  )    &
                &         * tmask(ji,jj) / ( e1t(ji,jj) * e2t(ji,jj) )      ! quantity on T grid
 
           cisah2(ji,jj)  = 0.25 * fmask(ji,jj) * ( cisah2t(ji,jj) * cisah2t(ji,jj)        &
                &         + cisah2t(ji+1,jj) * cisah2t(ji+1,jj) +  cisah2t(ji,jj+1)        &
                &         * cisah2t(ji,jj+1)  + cisah2t(ji+1,jj+1) * cisah2t(ji+1,jj+1) )       ! quantity computed on f grid

               okubow(ji,jj) = cisah1(ji,jj) * cisah1(ji,jj) + cisah2(ji,jj) - rotn(ji,jj)*rotn(ji,jj)

        END DO
     END DO

     IF ( lperio ) okubow(npiglo,:) = okubow(2, :)
     ! write rotn on file at level k and at time jt
     ierr = putvar(ncout, id_varout(1), okubow, 1, npiglo, npjglo, ktime=jt)
  END DO
  ierr = closeout(ncout)

END PROGRAM cdfokubow

