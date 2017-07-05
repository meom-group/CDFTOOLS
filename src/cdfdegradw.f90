PROGRAM cdfdegradw
  !!======================================================================
  !!                   ***  PROGRAM  cdfdegradw  ***
  !!=====================================================================
  !!  ** Purpose : Degrade horizontal resolution of NEMO gridW output
  !!
  !!  ** Method  : Degradation procedure ensuring the conservation
  !!               of water fluxes.
  !!
  !! History : 1.0 : 11/2011  : X. Meunier : Original code
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: nflsdc           ! array for fluid sub-domain counter
  INTEGER(KIND=4), DIMENSION(:),     ALLOCATABLE :: ipk, id_varout, id_varout2
  INTEGER(KIND=4)                           :: ji, jj, jk, jt, jri, jrj   ! dummy loop index
  INTEGER(KIND=4)                           :: ik, ii, ij         !
  INTEGER(KIND=4)                           :: ix, iy             !
  INTEGER(KIND=4)                           :: ierr               ! error status
  INTEGER(KIND=4)                           :: nri=1, nrj=1       ! scale degrading factor along i and j  
  INTEGER(KIND=4)                           :: iimin=2            ! indice from where the filtering begins
  INTEGER(KIND=4)                           :: ijmin=2            ! indice from where the filtering begins
  INTEGER(KIND=4)                           :: narg, iargc, ijarg ! browse command line
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                           :: ncout, ncout2      ! ncid of output file
  INTEGER(KIND=4)                           :: msks               ! mask sum on each coarse cells
  INTEGER(KIND=4)                           :: nvpk               ! vertical levels in working variable
  INTEGER(KIND=4)                           :: ikx=1, iky=1       ! dims of netcdf output file

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: glamf, gphif       ! lon, lat
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1, e2, zv         ! metrics, variable
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rdumlon, rdumlat   ! dummy lon/lat for output file
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rdummymean         ! array for mean value on output file
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask              ! 2D mask at current level
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                ! time counter
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: gdep               ! depth 
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: zdep               ! depth of the whole vertical levels
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: e31d               ! 1d vertical spacing
  REAL(KIND=4)                              :: zspval             ! missing value

  REAL(KIND=8)                              :: dvol, dsum, dsurf  ! cumulated values
  REAL(KIND=8)                              :: dvol2d, dsum2d     !

  CHARACTER(LEN=256)                        :: cf_wfil            ! input filename
  CHARACTER(LEN=256)                        :: cf_ncout,cf_ncout2 ! output file name
  !
  CHARACTER(LEN=256)                        :: cv_nam             ! current variable name
  CHARACTER(LEN=256)                        :: cv_dep             ! deptht name
  CHARACTER(LEN=20)                         :: cv_e1, cv_e2       ! horizontal metrics names
  CHARACTER(LEN=20)                         :: cv_e3, cv_e31d     ! vertical metrics names
  CHARACTER(LEN=20)                         :: cv_msk             ! mask variable name
  !
  CHARACTER(LEN=256)                        :: clunits            ! attribute of output file : units
  CHARACTER(LEN=256)                        :: cllong_name        !     "      long name
  CHARACTER(LEN=256)                        :: cglobal            !     "      global 
  CHARACTER(LEN=256)                        :: clshort_name       !     "      short name
  CHARACTER(LEN=256)                        :: cldum              ! dummy char variable

  TYPE(variable), DIMENSION(:), ALLOCATABLE :: stypvar            ! structure of output
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: stypvar2           ! structure of sub-domain counter output

  LOGICAL                                   :: lfull     = .false.! full step  flag
  LOGICAL                                   :: lchk               ! flag for missing files
 
 !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfdegradw IN-Wfile IN-var ri rj [i0 j0]'
     PRINT *,'       ... [-full]'
     PRINT *,'       '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Degrad the horizontal resolution of NEMO W-grid ouput,       ' 
     PRINT *,'       for each z-level and time step, with a ratio of ri along     '
     PRINT *,'       x direction and rj along y direction. If specified, the input'
     PRINT *,'       grid is considered starting from the indices i0 and j0.      ' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       IN-Wfile  : netcdf W-file.' 
     PRINT *,'       IN-var    : name of netcdf variable to work with' 
     PRINT *,'       ri        : degradation ratio for x-direction   ' 
     PRINT *,'       rj        : degradation ratio for y-direction   ' 
     PRINT *,'      '
     PRINT *,'     OPTIONS : '
     PRINT *,'       [i0 j0] : spatial indices from where starting the procedure   '
     PRINT *,'                 of degradation.                                    ' 
     PRINT *,'       [-full] : flag for full steps grid, instead of default partial'
     PRINT *,'                 steps.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       Files ', TRIM(cn_fhgr),', ', TRIM(cn_fzgr),', ', TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : degraded_cdfvar.nc '
     PRINT *,'       netcdf file : flsdc.nc'
     PRINT *,'      '
    STOP
  ENDIF

  cglobal = 'Partial step computation'
  ijarg = 1 ; ii = 0
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 
     SELECT CASE (cldum) 
     CASE ('-full' )
        lfull = .true. 
        cglobal = 'full step computation'
     CASE DEFAULT 
         ii=ii+1
        SELECT CASE (ii) 
        CASE ( 1 ) ; cf_wfil  = cldum 
        CASE ( 2 ) ; cv_nam = cldum 
        CASE ( 3 ) ; READ(cldum,*) nri  
        CASE ( 4 ) ; READ(cldum,*) nrj  
        CASE ( 5 ) ; READ(cldum,*) iimin
        CASE ( 6 ) ; READ(cldum,*) ijmin
        CASE DEFAULT
          PRINT *, ' ERROR : Too many arguments ...'
          STOP 99
        END SELECT
     END SELECT
  END DO

  lchk = chkfile(cn_fhgr)
  lchk = chkfile(cn_fzgr) .OR. lchk
  lchk = chkfile(cn_fmsk) .OR. lchk
  lchk = chkfile(cf_wfil) .OR. lchk
  IF ( lchk ) STOP 99 ! missing file

  cv_dep = 'none'
  npiglo = getdim (cf_wfil, cn_x)
  npjglo = getdim (cf_wfil, cn_y)
  npk    = getdim (cf_wfil, cn_z, cdtrue=cv_dep, kstatus=ierr)

  IF (ierr /= 0 ) THEN
     npk   = getdim (cf_wfil, 'z', cdtrue=cv_dep, kstatus=ierr)
     IF (ierr /= 0 ) THEN
        npk   = getdim (cf_wfil, 'sigma', cdtrue=cv_dep, kstatus=ierr)
        IF ( ierr /= 0 ) THEN
           npk = getdim (cf_wfil, 'nav_lev', cdtrue=cv_dep, kstatus=ierr)
           IF ( ierr /= 0 ) THEN
              PRINT *,' assume file with no depth'
              npk=0
           ENDIF
        ENDIF
     ENDIF
  ENDIF

  npt    = getdim (cf_wfil, cn_t)
  nvpk   = getvdim(cf_wfil, cv_nam)

  IF (npk   == 0 ) THEN ; npk = 1;  ENDIF ! no depth dimension ==> 1 level
  IF (iimin < 2 )  THEN
     PRINT *,'iimin value is too low'
     STOP 99
  END IF
  IF (ijmin  < 2 ) THEN
     PRINT *,'ijmin value is too low'
     STOP 99
  END IF
  npiglo = ( (npiglo - iimin ) / nri )*nri
  npjglo = ( (npjglo - ijmin ) / nrj )*nrj
  ikx = npiglo / nri
  iky = npjglo / nrj

  PRINT *,'ikx    : ',ikx
  PRINT *,'iky    : ',iky

  IF (nvpk == 2 ) npk = 1
  ! IF (nvpk == 3 ) nvpk = npk

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt
  PRINT *, 'nvpk   = ', nvpk
  PRINT *, 'depth dim name is ',TRIM(cv_dep)

  ! Allocate arrays
  ALLOCATE (zmask(npiglo,npjglo) )
  ALLOCATE (zv(npiglo,npjglo) )
  ALLOCATE (e1(npiglo,npjglo), e2(npiglo,npjglo) )
  ALLOCATE (gdep (npk), e31d(npk), tim(npt) )
  ALLOCATE (zdep(npk) )

  cv_e1    = cn_ve1t
  cv_e2    = cn_ve2t
  cv_e3    = 'e3w_ps'
  cv_e31d  = cn_ve3w
  cv_msk = 'tmask'
  cv_dep   = cn_gdepw
 
  e1(:,:) = getvar  (cn_fhgr, cv_e1,  1,  npiglo, npjglo, kimin=iimin, kjmin=ijmin)
  e2(:,:) = getvar  (cn_fhgr, cv_e2,  1,  npiglo, npjglo, kimin=iimin, kjmin=ijmin)
  IF ( lfull )  e31d(:) = getvare3(cn_fzgr, cv_e31d, npk)
  zdep(:) = getvare3(cn_fzgr, cv_dep, npk)
  gdep(:) = zdep(  1  :npk)

  ALLOCATE ( stypvar( 1 ), ipk( 1 ), id_varout( 1 ) )
  ALLOCATE ( stypvar2( 1 ), id_varout2( 1 ) )
  ALLOCATE ( rdumlon(ikx,iky), rdumlat(ikx,iky), rdummymean(ikx,iky) ) 
  ALLOCATE ( nflsdc(ikx,iky) )
  ALLOCATE ( glamf(npiglo+1,npjglo+1), gphif(npiglo+1,npjglo+1))

  glamf = getvar(cn_fhgr, cn_glamf, 1, npiglo+1 , npjglo+1 , kimin=iimin-1 , kjmin=ijmin-1 )
  gphif = getvar(cn_fhgr, cn_gphif, 1, npiglo+1 , npjglo+1 , kimin=iimin-1 , kjmin=ijmin-1 )
  DO jj=1,iky
     DO ji=1,ikx
        rdumlat(ji,jj) = 0.5 * ( gphif((ji-1)*nri+1,(jj-1)*nrj+1) + gphif( ji   *nri+1,jj*nrj+1) )
        rdumlon(ji,jj) = 0.5 * ( glamf( ji   *nri+1,(jj-1)*nrj+1) + glamf((ji-1)*nri+1,jj*nrj+1) )
     END DO
  END DO

  ipk(:) = npk ! mean for each level

  ierr=getvaratt (cf_wfil, cv_nam, clunits, zspval, cllong_name, clshort_name)

  ! define new variables for output 
  stypvar%rmissing_value    = 99999.
  stypvar%scale_factor      = 1.
  stypvar%add_offset        = 0.
  stypvar%savelog10         = 0.
  stypvar%conline_operation = 'N/A'

  stypvar%cunits         = TRIM(clunits)
  stypvar%valid_min      = -100000.
  stypvar%valid_max      = 100000.
  stypvar%cname          = 'degraded_'//TRIM(cv_nam)
  stypvar%clong_name     = 'degraded '//TRIM(cllong_name)
  stypvar%cshort_name    = 'degraded_'//TRIM(clshort_name)

  stypvar2%rmissing_value    = 99999.
  stypvar2%scale_factor      = 1.
  stypvar2%add_offset        = 0.
  stypvar2%savelog10         = 0.
  stypvar2%conline_operation = 'N/A'

  stypvar2%valid_min      = 0.
  stypvar2%valid_max      = 100000.
  stypvar2%cname          = 'flsdc'
  stypvar2%clong_name     = 'fluid sub-domain counter'
  stypvar2%cshort_name    = 'fdsdc'

  IF ( nvpk == 2 ) THEN
     stypvar%caxis       = 'TYX'
     stypvar2%caxis      = 'YX'
  ELSE
     stypvar%caxis       = 'TZYX'
     stypvar2%caxis      = 'ZYX'
  END IF

  ! create output fileset
  cf_ncout  = 'degraded_'//TRIM(cv_nam)//'.nc'
  cf_ncout2 = 'flsdc.nc'
  ncout  = create      (cf_ncout ,   'none',  ikx,   iky,   npk, cdep=cv_dep)
  ncout2 = create      (cf_ncout2,   'none',  ikx,   iky,   npk, cdep=cv_dep)
  ierr  = createvar   (ncout ,      stypvar , 1 , ipk,   id_varout  , cdglobal=TRIM(cglobal) ) 
  ierr  = createvar   (ncout2,      stypvar2, 1 , ipk,   id_varout2 , cdglobal=TRIM(cglobal) ) 
  ierr  = putheadervar(ncout ,      cf_wfil,  ikx, iky, npk, pnavlon=rdumlon, pnavlat=rdumlat, pdep=gdep, cdep=cv_dep)
  ierr  = putheadervar(ncout2,      cn_fmsk,  ikx, iky, npk, pnavlon=rdumlon, pnavlat=rdumlat, pdep=gdep, cdep=cv_dep)
  tim   = getvar1d(cf_wfil, cn_vtimec, npt)
  ierr  = putvar1d(ncout,  tim,       npt, 'T')

  DO jk = 1, npk
     ik = jk
     ! Get variable v at ik
     zmask(:,:) = getvar(cn_fmsk, cv_msk, ik, npiglo, npjglo, kimin=iimin, kjmin=ijmin          )
     DO jt=1,npt
        ! 
        zv   (:,:) = getvar(cf_wfil,  cv_nam,   ik, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=jt)
        DO jj=1,iky
           DO ji=1,ikx
              dvol = 0.d0
              dsum = 0.d0
              msks = 0
              DO jrj=1,nrj
                 DO jri=1,nri
                    ix = (ji-1)*nri+jri
                    iy = (jj-1)*nrj+jrj
                    dvol2d = e1(ix,iy) * e2(ix,iy) * zmask(ix,iy)
                    msks = msks + zmask(ix,iy)
                    dvol = dvol + dvol2d
                    dsum = dsum + zv(ix,iy) * dvol2d
                 END DO
              END DO
              IF ( dvol > 0.d0 ) THEN
                 rdummymean(ji,jj) = dsum/dvol
              ELSE
                 rdummymean(ji,jj) = 99999.
              END IF
              nflsdc(ji,jj) = msks
           END DO
        END DO
        ierr = putvar(ncout , id_varout(1), rdummymean , jk, ikx, iky, ktime=jt )
     END DO ! time loop
     ierr = putvar(ncout2, id_varout2(1), DBLE(nflsdc), jk, ikx, iky )
  END DO  

  ierr = closeout(ncout )
  ierr = closeout(ncout2)

9004 FORMAT(f9.2,' ',i2,' ',f9.2)

END PROGRAM cdfdegradw
