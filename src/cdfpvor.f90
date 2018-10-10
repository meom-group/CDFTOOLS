PROGRAM cdfpvor
  !!======================================================================
  !!                     ***  PROGRAM  cdfpvor  ***
  !!=====================================================================
  !!  ** Purpose : Compute the Ertel Potential vorticity
  !!
  !!  ** Method  : Formula :
  !!      Qpot = drho/dz * ( f + xsi ) = Qstr + Qrel
  !!      * f is the Coriolis factor, computed from the latitudes of the T-grid :
  !!        f(i,j) = 2 * omega * sin ( phit(i,j) * pi / 180 )
  !!
  !!      * xsi is the relative vorticity (vertical component of the velocity curl),
  !!        computed from the relative vorticity of the F-points interpolated at
  !!        the T-points :
  !!        xsif(i,j) = ( ue(i,j) - ue(i,j+1) - ve(i,j) + ve(i+1,j) ) / areaf(i,j)
  !!         with : ue(i,j) = U(i,j) * e1u(i,j)
  !!                ve(i,j) = V(i,j) * e2v(i,j)
  !!             areaf(i,j) = e1f(i,j) * e2f(i,j)
  !!        xsi(i,j) = ( xsif(i-1,j-1) + xsif(i-1,j) + xsif(i,j-1) + xsif(i,j) ) / 4
  !!                 = (  ue(i-1,j-1) + ue(i,j-1) - ue(i-1,j+1) - ue(i,j+1)
  !!                    - ve(i-1,j-1) - ve(i-1,j) + ve(i+1,j-1) + ve(i+1,j) )
  !!                    / 4 / areat(i,j)
  !!        with : areat(i,j) = e1t(i,j) * e2t(i,j)
  !!        units : U, V in m.s-1
  !!           e1u, e2v, e1f, e2f in m
  !!           f, xsi in s-1
  !!           Qpot, Qrel, Qstr in 1.e-7 kg.m-4.s-1
  !!
  !! History : 2.1  : 12/2005  : A.M. Treguier : Original code
  !!           3.0  : 05/2011  : J.M. Molines  : Doctor norm + Lic., merge with cdfpv
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!-------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class derived_fields
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                             :: jk, jj, ji, jt       ! dummy loop index
  INTEGER(KIND=4)                             :: it                   ! time index for vvl
  INTEGER(KIND=4)                             :: ierr                 ! working integer
  INTEGER(KIND=4)                             :: narg, iargc          ! command line
  INTEGER(KIND=4)                             :: ijarg                ! command line
  INTEGER(KIND=4)                             :: npiglo, npjglo       ! size of the domain
  INTEGER(KIND=4)                             :: npk, npt             ! size of the domain
  INTEGER(KIND=4)                             :: iup=1, idown=2, itmp ! working interger
  INTEGER(KIND=4)                             :: ncout                ! ncid for output file
  INTEGER(KIND=4)                             :: nvar=3               ! number of output variable
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE  :: ipk, id_varout       ! output variable id's

  REAL(KIND=4)                                :: zpi, zomega, rau0sg  ! physical constant
  REAL(KIND=4)                                :: zsps                 ! Missing value for salinity
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: gdepw                ! deptht
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: e31d                 ! metric for full step
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zn2                  ! Brunt Vaissala frequency
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: tmask                ! tmask from salinity
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e3w                  ! vertical metric
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1u, e2v             ! horizontal metric
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1t, e2t             ! horizontal metric at T point
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: gphit                ! latitude of t point
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: ztemp, zsal, zwk     ! array to ead 2 layer of data

  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: dtim                 ! time counter
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dun, dvn             ! velocity component and flx
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: drotn                ! curl of the velocity
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: d2fcor               ! coriolis term at T point
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dstretch             ! stretching vorticity
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dareat               ! area of T cells

  CHARACTER(LEN=256)                          :: cf_tfil              ! input T file
  CHARACTER(LEN=256)                          :: cf_sfil              ! salinity file (option)
  CHARACTER(LEN=256)                          :: cf_ufil              ! input U file
  CHARACTER(LEN=256)                          :: cf_vfil              ! input V file
  CHARACTER(LEN=256)                          :: cf_out='pvor.nc'     ! output file
  CHARACTER(LEN=256)                          :: cf_e3w               ! file holding e3w in the case of vvl
  CHARACTER(LEN=256)                          :: cldum                ! dummy character variable

  TYPE(variable), DIMENSION(:), ALLOCATABLE   :: stypvar              ! structure for attribute

  LOGICAL                                     :: lfull   = .FALSE.    ! flag for full step
  LOGICAL                                     :: lertel  = .TRUE.     ! flag for large scale pv
  LOGICAL                                     :: lchk    = .FALSE.    ! flag for missing files
  LOGICAL                                     :: lnc4    = .FALSE.    ! flag for netcdf4 chunking and deflation
  LOGICAL                                     :: l_metric= .TRUE.     ! flag for netcdf4 chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfpvor -t T-file  -u U-file -v V-file [-s S-file] [-full] [-lspv] ...'
     PRINT *,'           ... [-nometric] [-o OUT-file] [-nc4] [-vvl W-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the Ertel potential vorticity and save the relative  ' 
     PRINT *,'       vorticity, the stretching and the total potential vorticity. '
     PRINT *,'       Qtot = ( f + xsi ) . D(rho)/D(z)  = Qstrech + Qrel           '
     PRINT *,'       With -lspv option, compute only Qstretch or Large Scale P V '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file : netcdf file for temperature and salinity.           ' 
     PRINT *,'            If salinity not in T-file use -s option.'
     PRINT *,'       -u U-file : netcdf file for zonal component of the velocity.    '
     PRINT *,'       -v V-file : netcdf file for meridional component of the velocity.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-s S-file] : specify salinity file if not T-file.'
     PRINT *,'       [-full ] : indicate a full step configuration.                ' 
     PRINT *,'       [-lspv ] : calculate only the large scale potential vorticity.'
     PRINT *,'                  ( replace the old cdflspv tool).'
     PRINT *,'                  If used only T-file is required, no need for velocities.'
     PRINT *,'       [-nc4 ] :  use netcdf4 with chunking and deflation '
     PRINT *,'       [-o OUT-file] : use output file instead of default ',TRIM(cf_out)
     PRINT *,'       [-vvl W-file] : use time-varying vertical metrics. W-file holds the'
     PRINT *,'                  time-varying e3w vertical metrics.'
     PRINT *,'       [-nometric] : Do not use metrics for computing the vorticity.'
     PRINT *,'                 Assume all metrics to be 1, arbitrary.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fhgr),' and ',TRIM(cn_fzgr)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : vorelvor (1.e-7 kg.m-4.s-1 ) relative vorticity'
     PRINT *,'                     vostrvor (1.e-7 kg.m-4.s-1 ) stretching vorticity'
     PRINT *,'                     vototvor (1.e-7 kg.m-4.s-1 ) total potential vorticity'
     PRINT *,'                  Ertel PV are located at T points.'
     PRINT *,'           '
     PRINT *,'       With option -lspv :'
     PRINT *,'       netcdf file : lspv.nc'
     PRINT *,'         variables :  volspv  (1.e-7 kg.m-4.s-1 ) large scale potential '
     PRINT *,'              vorticity. LSPV is located at W points.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfcurl ( compute only the curl on 1 level)'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1 
  cf_sfil = 'none'
  DO WHILE ( ijarg <= narg )
     CALL getarg( ijarg, cldum ) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-t'      ) ; CALL getarg( ijarg, cf_tfil) ; ijarg = ijarg + 1
     CASE ( '-u'      ) ; CALL getarg( ijarg, cf_ufil) ; ijarg = ijarg + 1
     CASE ( '-v'      ) ; CALL getarg( ijarg, cf_vfil) ; ijarg = ijarg + 1
      ! options
     CASE ( '-s'      ) ; CALL getarg( ijarg, cf_sfil) ; ijarg = ijarg + 1
     CASE ( '-full'   ) ; lfull  = .TRUE.
     CASE ( '-lspv'   ) ; lertel = .FALSE. ; nvar = 1  ; cf_out = 'lspv.nc'
     CASE ( '-nc4'    ) ; lnc4   = .TRUE.
     CASE ( '-o'      ) ; CALL getarg( ijarg, cf_out ) ; ijarg = ijarg + 1
     CASE ( '-vvl'    ) ; lg_vvl = .TRUE.
        ;                 CALL getarg( ijarg, cf_e3w ) ; ijarg = ijarg + 1
     CASE ('-nometric') ; l_metric   = .FALSE.
        ;               ; cf_out = 'pvor_grid.nc'
     CASE DEFAULT       ;  PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  IF ( cf_sfil == 'none' ) cf_sfil=cf_tfil

  IF ( l_metric ) THEN
    lchk = lchk .OR. chkfile( cn_fzgr)
    lchk = lchk .OR. chkfile( cn_fhgr)
  ENDIF

  lchk = lchk .OR. chkfile( cf_tfil)
  lchk = lchk .OR. chkfile( cf_sfil)
  IF ( lertel ) THEN
     lchk = lchk .OR. chkfile( cf_ufil)
     lchk = lchk .OR. chkfile( cf_vfil)
  ENDIF
  IF ( lg_vvl ) THEN
     lchk = lchk .OR. chkfile( cf_e3w )
  ENDIF
  IF ( lchk   ) STOP 99 ! missing file

  ! Look for MissingValue for salinity
  zsps = getspval(cf_sfil, cn_vosaline)

  IF ( lg_vvl ) THEN
     cn_fe3w = cf_e3w
     cn_ve3w = cn_ve3wvvl
  ENDIF

  npiglo = getdim (cf_tfil, cn_x)
  npjglo = getdim (cf_tfil, cn_y)
  npk    = getdim (cf_tfil, cn_z)
  npt    = getdim (cf_tfil, cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE ( e1u(npiglo,npjglo),   e1t(npiglo,npjglo)     )
  ALLOCATE ( e2v(npiglo,npjglo),   e2t(npiglo,npjglo)     )
  ALLOCATE ( gphit(npiglo,npjglo), d2fcor(npiglo,npjglo)  )
  ALLOCATE ( tmask(npiglo,npjglo)                         )
  ALLOCATE ( gdepw(npk), dtim(npt)                        )

  ALLOCATE ( dstretch(npiglo,npjglo)                      )
  IF ( lertel ) THEN
     ALLOCATE ( dareat(npiglo,npjglo)                     )
     ALLOCATE ( dun(npiglo,npjglo),    dvn(npiglo,npjglo) )
     ALLOCATE ( drotn(npiglo,npjglo)                      )
  ENDIF
  IF ( lfull ) ALLOCATE ( e31d(npk)                       )

  IF ( l_metric ) THEN
    e1u   = getvar(cn_fhgr, cn_ve1u,  1, npiglo, npjglo)
    e1t   = getvar(cn_fhgr, cn_ve1t,  1, npiglo, npjglo)
    e2v   = getvar(cn_fhgr, cn_ve2v,  1, npiglo, npjglo)
    e2t   = getvar(cn_fhgr, cn_ve2t,  1, npiglo, npjglo)
    gphit = getvar(cn_fhgr, cn_gphit, 1, npiglo, npjglo)
  ELSE
    e1u   = 1.
    e1t   = 1.
    e2v   = 1.
    e2t   = 1.
    gphit = 10.  ! will give a constant value for f (10 deg N, arbitrary)
  ENDIF

  rau0sg      = 1020./9.81
  zpi         = ACOS(-1.)
  zomega      = 2.0  * zpi /(3600*24)
  d2fcor(:,:) = 2.d0 * zomega * SIN(gphit(:,:)*zpi/180.0)

  IF ( lertel ) THEN
     dareat(:,:) = 4.d0 * e1t(:,:) * e2t(:,:)  ! factor of 4 to normalize relative vorticity
  ENDIF

  gdepw(:) = getvare3(cn_fzgr, cn_gdepw, npk)

  ALLOCATE (ztemp(npiglo,npjglo,2), zsal(npiglo,npjglo,2)) 
  ALLOCATE (zwk(npiglo,npjglo,2)                         )
  ALLOCATE (zn2(npiglo,npjglo) ,    e3w(npiglo,npjglo)   )
  ALLOCATE ( stypvar(nvar), ipk(nvar), id_varout(nvar)   )

  ! create output fileset
  CALL CreateOutput

  IF ( lfull  ) THEN
    IF ( l_metric ) THEN
      e31d = getvare3( cn_fzgr, cn_ve3w1d, npk )
    ELSE
      e31d = 1.
    ENDIF
  ENDIF

  DO jt=1,npt
     IF ( lg_vvl ) THEN ; it=jt
     ELSE               ; it=1
     ENDIF
     !  2 levels of T and S are required : iup,idown (with respect to W level)
     !  Compute from bottom to top (for vertical integration)
     PRINT *,'time=',jt,'(days:',dtim(jt)/86400.,')'
     ztemp(:,:,idown) = getvar(cf_tfil, cn_votemper, npk-1, npiglo, npjglo, ktime=jt)
     zsal( :,:,idown) = getvar(cf_sfil, cn_vosaline, npk-1, npiglo, npjglo, ktime=jt)

     ! -------------------------------- LOOP OVER LEVELS
     DO jk = npk-1, 1, -1 
        PRINT *,'            level ',jk
        IF ( lertel ) THEN 
           ! ------------------------------------RELATIVE VORTICITY FIRST
           dun(:,:) = getvar(cf_ufil, cn_vozocrtx, jk ,npiglo, npjglo, ktime=jt)
           dvn(:,:) = getvar(cf_vfil, cn_vomecrty, jk ,npiglo, npjglo, ktime=jt)
           dun(:,:) = dun(:,:)*e1u(:,:)
           dvn(:,:) = dvn(:,:)*e2v(:,:)
           !     relative vorticity at T point
           drotn(:,:) = 0.d0
           DO jj = 2, npjglo -1 
              DO ji = 2, npiglo -1    
                 drotn(ji,jj) = ( dun(ji-1,jj-1)  + dun(ji,jj-1)  &
                        &        -dun(ji-1,jj+1)  - dun(ji,jj+1)  &
                        &        -dvn(ji-1,jj-1)  - dvn(ji-1,jj)  &
                        &        +dvn(ji+1,jj-1)  + dvn(ji+1,jj)) &
                        / dareat(ji,jj) 
              END DO
           END DO
        ENDIF

        !  now  tmask and Vaisala Frequency bn2
        IF ( jk > 1) THEN 
           tmask(:,:)=1.
           ztemp(:,:,iup) = getvar(cf_tfil, cn_votemper, jk-1 ,npiglo, npjglo, ktime=jt)
           zsal(:,:,iup)  = getvar(cf_sfil, cn_vosaline, jk-1 ,npiglo, npjglo, ktime=jt)
           WHERE(zsal(:,:,idown) == zsps ) tmask = 0

           IF ( l_metric ) THEN
              IF ( lfull ) THEN ; e3w(:,:) = e31d(jk)
              ELSE              ; e3w(:,:) = getvar(cn_fe3w, cn_ve3w, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl )
              ENDIF
           ELSE                 ; e3w(:,:) = 1.
           ENDIF
              

           WHERE (e3w == 0 ) e3w = 1.

           zwk(:,:,iup) = eosbn2 ( ztemp, zsal, gdepw(jk), e3w, npiglo, npjglo ,iup, idown)* tmask(:,:)
           !
           IF ( lertel ) THEN ! put zn2 at T level (k )
              WHERE ( zwk(:,:,idown) == 0 ) 
                 zn2(:,:) =  zwk(:,:,iup)
              ELSEWHERE                     
                 zn2(:,:) = 0.5 * ( zwk(:,:,iup) + zwk(:,:,idown) ) * tmask(:,:)
              END WHERE
           ELSE               ! keep bn2 at w points
              zn2(:,:) = zwk(:,:,iup) * tmask(:,:)
           ENDIF
        ENDIF
        !
        !   now rotn will be converted to relative vorticity and zn2 to stretching
        dstretch(:,:) = d2fcor(:,:)* rau0sg * zn2(:,:)

        IF ( lertel ) THEN
           drotn(:,:)    = drotn(:,:) * rau0sg * zn2(:,:)
           ! write the three variables on file at level k
           ierr = putvar(ncout, id_varout(1), REAL( drotn          )*1.e7, jk, npiglo, npjglo, ktime=jt)
           ierr = putvar(ncout, id_varout(2), REAL( dstretch       )*1.e7, jk, npiglo, npjglo, ktime=jt)
           ierr = putvar(ncout, id_varout(3), REAL((drotn+dstretch))*1.e7, jk, npiglo, npjglo, ktime=jt)
        ELSE
           ! save absolute value of dstretch, as in olf cdflspv
           ierr = putvar(ncout, id_varout(1), REAL( ABS(dstretch)  )*1.e7, jk, npiglo, npjglo, ktime=jt)
        ENDIF

        itmp = idown ; idown = iup ; iup = itmp

     END DO  ! loop to next level

     !  set zero at bottom and surface
     zwk(:,:,1) = 0.e0
     ierr = putvar(ncout, id_varout(1), zwk(:,:,1), 1,   npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(1), zwk(:,:,1), npk, npiglo, npjglo, ktime=jt)

     IF (lertel ) THEN
        ierr = putvar(ncout, id_varout(2), zwk(:,:,1), 1,   npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(3), zwk(:,:,1), 1,   npiglo, npjglo, ktime=jt)

        ierr = putvar(ncout, id_varout(2), zwk(:,:,1), npk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(3), zwk(:,:,1), npk, npiglo, npjglo, ktime=jt)
     ENDIF
  END DO   ! loop on time

  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ipk(:)= npk                   ! Those three variables are  3D
    stypvar%cunits            = '1.e-7 kg.m-4.s-1'
    stypvar%rmissing_value    = 0.
    stypvar%valid_min         = -10000.
    stypvar%valid_max         = 10000.
    stypvar%conline_operation = 'N/A'
    stypvar%caxis             = 'TZYX'

    DO ji = 1, nvar
       stypvar(ji)%ichunk = (/npiglo,MAX(1,npjglo/30), 1, 1 /)
    ENDDO

    IF (lertel ) THEN
       ! define variable name and attribute
       stypvar(1)%cname = 'vorelvor' ; stypvar(1)%clong_name = 'Relative_component_of_Ertel_PV'
       stypvar(2)%cname = 'vostrvor' ; stypvar(2)%clong_name = 'Stretching_component_of_Ertel_PV'
       stypvar(3)%cname = 'vototvor' ; stypvar(3)%clong_name = 'Ertel_potential_vorticity'

       stypvar(1)%cshort_name = 'vorelvor'
       stypvar(2)%cshort_name = 'vostrvor'
       stypvar(3)%cshort_name = 'vototvor'

       ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk       , ld_nc4=lnc4 )
       ierr  = createvar   (ncout,  stypvar, nvar,   ipk,    id_varout , ld_nc4=lnc4 )
       ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk       )
    ELSE
       stypvar(1)%cname = 'volspv' ; stypvar(1)%clong_name = 'Large Scale Potential_vorticity'
       stypvar(1)%cshort_name = 'volspv'
       ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk, cdep=TRIM(cn_vdepthw) , ld_nc4=lnc4 )
       ierr  = createvar   (ncout,  stypvar, nvar,   ipk,    id_varout                  , ld_nc4=lnc4 )
       ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk, pdep=gdepw            )
    ENDIF

    dtim = getvar1d(cf_ufil, cn_vtimec, npt      )
    ierr = putvar1d(ncout,   dtim,      npt, 'T' )

  END SUBROUTINE CreateOutput

END PROGRAM cdfpvor
