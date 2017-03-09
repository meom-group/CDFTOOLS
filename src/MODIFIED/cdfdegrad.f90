PROGRAM cdfdegrad
  !!======================================================================
  !!                   ***  PROGRAM  cdfdegrad  ***
  !!=====================================================================
  !!  ** Purpose : Degrade horizontal resolution of NEMO output
  !!
  !!  ** Method  : Degradation procedure ensuring the conservation
  !!               of tracers fluxes, preserving the tracer budget (T)
  !!               preserve the mass fluxes ( U,V,W)
  !!
  !! History : 1.0 : 11/2011  : X. Meunier : Original code
  !!           3.0 : 02/2017  : J.-M. Molines : merge of all cdfdegradx tools
  !!                            debuging but still not completed
  !!         : 4.0 : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class file_operations
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: nflsdc        ! array for fluid sub-domain counter
  INTEGER(KIND=4), DIMENSION(:),     ALLOCATABLE :: ipk, id_varout ! netcdf output
  INTEGER(KIND=4)                           :: ji, jj, jk, jt     ! dummy loop index
  INTEGER(KIND=4)                           :: jri, jrj           ! dummy loop index
  INTEGER(KIND=4)                           :: ik, ii, ij         !
  INTEGER(KIND=4)                           :: it                 ! time index for vvl
  INTEGER(KIND=4)                           :: ix, iy             !
  INTEGER(KIND=4)                           :: idep, idep_max     ! index for finding vertical dimension
  INTEGER(KIND=4)                           :: ierr               ! error status
  INTEGER(KIND=4)                           :: nri=1, nrj=1       ! scale degrading factor along i and j  
  INTEGER(KIND=4)                           :: iimin=2            ! indice from where the filtering begins
  INTEGER(KIND=4)                           :: ijmin=2            ! indice from where the filtering begins
  INTEGER(KIND=4)                           :: narg, iargc, ijarg ! browse command line
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                           :: ncout              ! ncid of output file
  INTEGER(KIND=4)                           :: msks               ! mask sum on each coarse cells
  INTEGER(KIND=4)                           :: nvpk               ! vertical levels in working variable
  INTEGER(KIND=4)                           :: ikx, iky           ! dims of netcdf output file

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask              ! 2D mask at current level
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: glamf, gphif       ! lon, lat
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1, e2, e3, zv     ! metrics, variable
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rdumlon, rdumlat   ! dummy lon/lat for output file
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rdummymean         ! array for mean value on output file
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                ! time counter
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: gdep               ! depth 
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: zdep               ! depth of the whole vertical levels
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: e31d               ! 1d vertical spacing
  REAL(KIND=4)                              :: zspval             ! missing value

  REAL(KIND=8)                              :: dvol, dsum, dsurf  ! cumulated values
  REAL(KIND=8)                              :: dvol2d, dsum2d     !

  CHARACTER(LEN=256)                        :: cf_in              ! input filename
  CHARACTER(LEN=256)                        :: cf_e3              ! vertical metrics filename
  CHARACTER(LEN=256)                        :: cf_out=''          ! output file name
  !
  CHARACTER(LEN=256)                        :: cv_nam             ! current variable name
  CHARACTER(LEN=20)                         :: cv_e1, cv_e2       ! horizontal metrics names
  CHARACTER(LEN=20)                         :: cv_e3, cv_e31d     ! vertical metrics names
  CHARACTER(LEN=20)                         :: cv_msk             ! mask variable name
  CHARACTER(LEN=256)                        :: cv_dep             ! deptht name
  CHARACTER(LEN=256), DIMENSION(:),ALLOCATABLE :: clv_dep             ! deptht name
  !
  CHARACTER(LEN=256)                        :: clunits            ! attribute of output file : units
  CHARACTER(LEN=256)                        :: cldum              ! dummy char variable
  CHARACTER(LEN=256)                        :: cllong_name        !     "      long name
  CHARACTER(LEN=256)                        :: cglobal            !     "      global 
  CHARACTER(LEN=256)                        :: clshort_name       !     "      short name
  CHARACTER(LEN=256)                        :: ctyp=''            ! Point on C grid

  TYPE(variable), DIMENSION(:), ALLOCATABLE :: stypvar            ! structure of output

  LOGICAL                                   :: lfull  = .FALSE.   ! full step  flag
  LOGICAL                                   :: lstart = .FALSE.   ! change staring point flag
  LOGICAL                                   :: lchk               ! flag for missing files
  LOGICAL                                   :: ll_pt  = .FALSE.   ! flag for T point
  LOGICAL                                   :: ll_pu  = .FALSE.   ! flag for U point
  LOGICAL                                   :: ll_pv  = .FALSE.   ! flag for V point
  LOGICAL                                   :: ll_pw  = .FALSE.   ! flag for W point

  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfdegrad -f IN-file -v IN-var -ratio ri rj -p T|U|V|W  [-start i0 j0]'
     PRINT *,'       ... [-full] [-vvl] [-o OUT-file]'
     PRINT *,'       '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Degrad the horizontal resolution of a NEMO ouput file,       ' 
     PRINT *,'       for each z-level and time step, with a ratio of ri along     '
     PRINT *,'       x direction and rj along y direction. If specified, the input'
     PRINT *,'       grid is considered starting from the indices i0 and j0.      ' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file  : netcdf input file on grid point corresponding to -p option.' 
     PRINT *,'       -v IN-var    : name of netcdf variable to work with' 
     PRINT *,'       -ratio ri rj : degradation ratio for x-direction and y-direction.' 
     PRINT *,'       -p T|U|V|W   : position of variable on C-grid.'
     PRINT *,'      '
     PRINT *,'     OPTIONS : '
     PRINT *,'       [-start i0 j0] : spatial indices from where the procedure of   '
     PRINT *,'                        degradation starts. ' 
     PRINT *,'       [-o OUT-file ]:  output filename instead of ''degraded_<IN-var>.nc'' '
     PRINT *,'       [-vvl ] : use time-varying vertical metrics.'
     PRINT *,'       [-full] : flag for full steps grid, instead of default partial'
     PRINT *,'                 steps.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       Files ', TRIM(cn_fhgr),', ', TRIM(cn_fzgr),', ', TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : degraded_<IN-var>.nc unless -o option is used'
     PRINT *,'          variables : degraded_<IN-VAR> '
     PRINT *,'                      flsdc : fluid subdomain counter '
     PRINT *,'      '
     STOP
  ENDIF

  cglobal = 'Partial step computation'
  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 
     SELECT CASE (cldum) 
     CASE ('-f'    ) ; CALL getarg(ijarg, cf_in   ) ; ijarg = ijarg + 1
     CASE ('-v'    ) ; CALL getarg(ijarg, cv_nam  ) ; ijarg = ijarg + 1
     CASE ('-ratio') ; CALL getarg(ijarg, cldum   ) ; ijarg = ijarg + 1 ; READ(cldum,*) nri
        CALL getarg(ijarg, cldum   ) ; ijarg = ijarg + 1 ; READ(cldum,*) nrj
     CASE ('-p'    ) ; CALL getarg(ijarg, ctyp    ) ; ijarg = ijarg + 1
        ! options
     CASE ('-o'    ) ; CALL getarg(ijarg, cf_out  ) ; ijarg = ijarg + 1
     CASE ('-full' ) ; lfull  = .TRUE. ; cglobal = 'full step computation'
     CASE ('-vvl'  ) ; lg_vvl = .TRUE. 
     CASE ('-start') ; lstart = .TRUE.
        CALL getarg(ijarg, cldum   ) ; ijarg = ijarg + 1 ; READ(cldum,*) iimin
        CALL getarg(ijarg, cldum   ) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmin
     CASE DEFAULT    ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option !' ; STOP
     END SELECT
  END DO

  lchk = chkfile(cn_fhgr)
  lchk = chkfile(cn_fzgr) .OR. lchk
  lchk = chkfile(cn_fmsk) .OR. lchk
  lchk = chkfile(cf_in  ) .OR. lchk
  IF ( lchk   ) STOP ! missing file

  SELECT CASE ( ctyp )
  CASE ( 'T' )
     ll_pt = .TRUE.
     IF ( .NOT. lstart ) THEN ; iimin=2 ; ijmin=2 ; ENDIF
     cf_e3    = cn_fe3t
     cv_e1    = cn_ve1t
     cv_e2    = cn_ve2t
     cv_e3    = cn_ve3t
     cv_e31d  = cn_ve3t
     cv_msk   = cn_tmask
  CASE ( 'U' )
     ll_pu = .TRUE.
     IF ( .NOT. lstart ) THEN ; iimin=1 ; ijmin=2 ; ENDIF
     cf_e3    = cn_fe3u
     cv_e1    = cn_ve1u
     cv_e2    = cn_ve2u
     cv_e3    = cn_ve3u
     cv_e31d  = cn_ve3u
     cv_msk   = cn_umask
  CASE ( 'V' )
     ll_pv = .TRUE.
     IF ( .NOT. lstart ) THEN ; iimin=2 ; ijmin=1 ; ENDIF
     cf_e3    = cn_fe3v
     cv_e1    = cn_ve1v
     cv_e2    = cn_ve2v
     cv_e3    = cn_ve3v
     cv_e31d  = cn_ve3v
     cv_msk   = cn_vmask
  CASE ( 'W' )
     ll_pw = .TRUE.
     IF ( .NOT. lstart ) THEN ; iimin=2 ; ijmin=2 ; ENDIF
     cf_e3    = cn_fe3w
     cv_e1    = cn_ve1t
     cv_e2    = cn_ve2t
     cv_e3    = cn_ve3w
     cv_e31d  = cn_ve3w
     cv_msk   = cn_tmask
  CASE DEFAULT ; PRINT *,' ERROR : C-grid point ', TRIM(ctyp),' not recognized' ; STOP
  END SELECT

  IF ( lg_vvl ) cf_e3 = cf_in  

  cv_dep = 'none'
  npiglo = getdim (cf_in  , cn_x)
  npjglo = getdim (cf_in  , cn_y)

  ! looking for npk among various possible name
  idep_max=8
  ALLOCATE ( clv_dep(idep_max) )
  clv_dep(:) = (/cn_z,'z','sigma','nav_lev','levels','ncatice','icbcla','icbsect'/)
  idep=1  ; ierr=1000

  DO WHILE ( ierr /= 0 .AND. idep <= idep_max )
     npk  = getdim (cf_in  , clv_dep(idep), cdtrue=cv_dep, kstatus=ierr)
     idep = idep + 1
  ENDDO

  IF ( ierr /= 0 ) THEN  ! none of the dim name was found
     PRINT *,' assume file with no depth'
     npk=0
  ENDIF

  npt    = getdim (cf_in  , cn_t   )
  nvpk   = getvdim(cf_in  , cv_nam )

  IF (npk   == 0 ) THEN ; npk = 1; ENDIF ! no depth dimension ==> 1 level
  ! retail horizontal fine resolution sizes to be multiple of the coarsening ratio

  SELECT CASE ( ctyp )
  CASE ( 'T' , 'W' ) 
     IF (iimin < 2  ) THEN ; PRINT *,'iimin value is too low' ; STOP ; ENDIF
     IF (ijmin < 2  ) THEN ; PRINT *,'ijmin value is too low' ; STOP ; ENDIF
     npiglo = ( (npiglo - iimin ) / nri )*nri
     npjglo = ( (npjglo - ijmin ) / nrj )*nrj
     ikx = npiglo / nri
     iky = npjglo / nrj
     ALLOCATE (e1(npiglo,npjglo), e2(npiglo,npjglo), e3(npiglo,npjglo) )
  CASE ( 'U' ) 
     IF (iimin < 1  ) THEN ; PRINT *,'iimin value is too low' ; STOP ; ENDIF
     IF (ijmin < 2  ) THEN ; PRINT *,'ijmin value is too low' ; STOP ; ENDIF
     npiglo = ( (npiglo - iimin - 1 ) / nri )*nri + 1
     npjglo = ( (npjglo - ijmin ) / nrj )*nrj
     ikx = ( npiglo - 1 ) / nri + 1
     iky = ( npjglo     ) / nrj
     ALLOCATE (e2(npiglo,npjglo), e3(npiglo,npjglo) )
  CASE ( 'V' ) 
     IF (iimin < 2  ) THEN ; PRINT *,'iimin value is too low' ; STOP ; ENDIF
     IF (ijmin < 1  ) THEN ; PRINT *,'ijmin value is too low' ; STOP ; ENDIF
     npiglo = ( (npiglo - iimin ) / nri )*nri
     npjglo = ( (npjglo - ijmin - 1 ) / nrj )*nrj + 1
     ikx = ( npiglo     ) / nri
     iky = ( npjglo - 1 ) / nrj + 1
     ALLOCATE (e1(npiglo,npjglo), e3(npiglo,npjglo) )
  END SELECT

  PRINT *,'ikx    : ',ikx
  PRINT *,'iky    : ',iky

  IF (nvpk == 2 ) npk = 1

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt
  PRINT *, 'nvpk   = ', nvpk
  PRINT *, 'depth dim name is ',TRIM(cv_dep)

  ! Allocate arrays
  ALLOCATE (zmask(npiglo,npjglo) )
  ALLOCATE (zv(npiglo,npjglo) )
  ALLOCATE (gdep (npk), e31d(npk), tim(npt) )
  ALLOCATE (zdep(npk) )

  IF ( ll_pt .OR. ll_pw .OR. ll_pv )  e1(:,:) = getvar  (cn_fhgr, cv_e1,  1,  npiglo, npjglo, kimin=iimin, kjmin=ijmin)
  IF ( ll_pt .OR. ll_pw .OR. ll_pu )  e2(:,:) = getvar  (cn_fhgr, cv_e2,  1,  npiglo, npjglo, kimin=iimin, kjmin=ijmin)

  IF ( lfull )  e31d(:) = getvare3(cn_fzgr, cv_e31d, npk)
  zdep(:) = getvare3(cn_fzgr, cv_dep, npk)
  gdep(:) = zdep(  1  :npk)

  ALLOCATE ( stypvar( 2 ), ipk( 2 ), id_varout( 2 ) )
  ALLOCATE ( rdumlon(ikx,iky), rdumlat(ikx,iky), rdummymean(ikx,iky) ) 
  ALLOCATE ( nflsdc(ikx,iky) )

  ! JMM: follow original X. Meunier code for u/v/t
  !      Can be done simpler
  SELECT CASE ( ctyp )
  CASE ( 'T' , 'W' )
     ALLOCATE ( glamf(npiglo+1,npjglo+1), gphif(npiglo+1,npjglo+1))
     glamf = getvar(cn_fhgr, cn_glamf, 1, npiglo+1 , npjglo+1 , kimin=iimin-1 , kjmin=ijmin-1 )
     gphif = getvar(cn_fhgr, cn_gphif, 1, npiglo+1 , npjglo+1 , kimin=iimin-1 , kjmin=ijmin-1 )

     DO jj=1,iky
        DO ji=1,ikx
           rdumlat(ji,jj) = 0.5 * ( gphif((ji-1)*nri+1,(jj-1)*nrj+1) + gphif( ji   *nri+1,jj*nrj+1) )
           rdumlon(ji,jj) = 0.5 * ( glamf( ji   *nri+1,(jj-1)*nrj+1) + glamf((ji-1)*nri+1,jj*nrj+1) )
        END DO
     END DO
  CASE ( 'U' )
     ALLOCATE ( glamf(npiglo,npjglo+1), gphif(npiglo,npjglo+1))
     glamf = getvar(cn_fhgr, cn_glamf, 1, npiglo , npjglo+1 , kimin=iimin , kjmin=ijmin-1 )
     gphif = getvar(cn_fhgr, cn_gphif, 1, npiglo , npjglo+1 , kimin=iimin , kjmin=ijmin-1 )

     DO ji=1,ikx
        DO jj=1,iky
           rdumlat(ji,jj) = 0.5 * ( gphif((ji-1)*nri+1,(jj-1)*nrj+1) + gphif((ji-1)*nri+1,jj*nrj+1) )
           rdumlon(ji,jj) = 0.5 * ( glamf((ji-1)*nri+1,(jj-1)*nrj+1) + glamf((ji-1)*nri+1,jj*nrj+1) )
        END DO
     END DO
  CASE ( 'V' )
     ALLOCATE ( glamf(npiglo+1,npjglo), gphif(npiglo+1,npjglo))
     glamf = getvar(cn_fhgr, cn_glamf, 1, npiglo+1 , npjglo , kimin=iimin-1 , kjmin=ijmin )
     gphif = getvar(cn_fhgr, cn_gphif, 1, npiglo+1 , npjglo , kimin=iimin-1 , kjmin=ijmin )
     DO jj=1,iky
        DO ji=1,ikx
           rdumlat(ji,jj) = 0.5 * ( gphif((ji-1)*nri+1,(jj-1)*nrj+1) + gphif(ji*nri+1,(jj-1)*nrj+1) )
           rdumlon(ji,jj) = 0.5 * ( glamf((ji-1)*nri+1,(jj-1)*nrj+1) + glamf(ji*nri+1,(jj-1)*nrj+1) )
        END DO
     END DO
  END SELECT
  DEALLOCATE ( glamf, gphif )

  CALL CreateOutput

  DO jt=1,npt
     IF ( lg_vvl ) THEN ; it=jt
     ELSE               ; it=1
     ENDIF
     DO jk = 1, npk
        ik = jk
        ! Get variable v at ik
        zmask(:,:) = getvar(cn_fmsk, cv_msk, ik, npiglo, npjglo, kimin=iimin, kjmin=ijmin          )

        IF (ll_pw ) THEN
           ! do nothing
        ELSE  ! read vertical metrics
           IF ( lfull ) THEN ; e3(:,:) = e31d(jk)
           ELSE              ; e3(:,:) = getvar(cf_e3, cv_e3,  ik, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=it, ldiom=.NOT.lg_vvl )
           END IF
        ENDIF
        ! 
        zv   (:,:) = getvar(cf_in  ,  cv_nam,   ik, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=jt)
        SELECT CASE ( ctyp )
        CASE ( 'T' )  ! weighted volume averaging
           DO jj=1,iky
              DO ji=1,ikx
                 dvol = 0.d0
                 dsum = 0.d0
                 msks = 0
                 ! volume averaging
                 DO jrj=1,nrj
                    DO jri=1,nri
                       ix = (ji-1)*nri+jri
                       iy = (jj-1)*nrj+jrj
                       dvol2d = e1(ix,iy) * e2(ix,iy) * e3(ix,iy) * zmask(ix,iy)
                       msks = msks + zmask(ix,iy)
                       dvol = dvol + dvol2d
                       dsum = dsum + zv(ix,iy) * dvol2d
                    END DO
                 END DO
                 IF ( dvol > 0.d0 ) THEN ; rdummymean(ji,jj) = dsum/dvol
                 ELSE                    ; rdummymean(ji,jj) = 99999.
                 END IF
                 nflsdc(ji,jj) = msks
              END DO
           END DO
        CASE ( 'U' )  ! lateral volume flux averaging
           DO ji=1,ikx
              ix = (ji-1)*nri+1
              DO jj=1,iky
                 dvol = 0.d0
                 dsum = 0.d0
                 msks = 0
                 DO jrj=1,nrj
                    iy = (jj-1)*nrj+jrj
                    ! JMM :note that X. Meunier code uses e1 instead of e2 ( bug ?)
                    ! JMM : for conservation of fluxes I think that dvol2d must not be masked. (bug ?)
                    dvol2d = e2(ix,iy) * e3(ix,iy) * zmask(ix,iy)
                    msks = msks + zmask(ix,iy)
                    dvol = dvol + dvol2d
                    dsum = dsum + zv(ix,iy) * dvol2d
                 END DO
                 IF ( dvol > 0.d0 ) THEN
                    rdummymean(ji,jj) = dsum/dvol
                 ELSE
                    rdummymean(ji,jj) = 99999.
                 END IF
                 nflsdc(ji,jj) = msks
              END DO
           END DO
        CASE ( 'V' ) ! lateral volume flux averaging
           DO jj=1,iky
              iy = (jj-1)*nrj+1
              DO ji=1,ikx
                 dvol = 0.d0
                 dsum = 0.d0
                 msks = 0
                 DO jri=1,nri
                    ix = (ji-1)*nri+jri
                    ! JMM :note that X. Meunier code uses e2 instead of e1 ( bug ?)
                    ! JMM : for conservation of fluxes I think that dvol2d must not be masked. (bug ?)
                    dvol2d = e1(ix,iy) * e3(ix,iy) * zmask(ix,iy)
                    msks = msks + zmask(ix,iy)
                    dvol = dvol + dvol2d
                    dsum = dsum + zv(ix,iy) * dvol2d
                 END DO
                 IF ( dvol > 0.d0 ) THEN
                    rdummymean(ji,jj) = dsum/dvol
                 ELSE
                    rdummymean(ji,jj) = 99999.
                 END IF
                 nflsdc(ji,jj) = msks
              END DO
           END DO
        CASE ( 'W' )  ! surface cell averaging 
           DO jj=1,iky
              DO ji=1,ikx
                 dvol = 0.d0
                 dsum = 0.d0
                 msks = 0
                 DO jrj=1,nrj
                    DO jri=1,nri
                       ix = (ji-1)*nri+jri
                       iy = (jj-1)*nrj+jrj
                       ! JMM : for conservation of fluxes I think that dvol2d must not be masked. (bug ?)
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
        END SELECT
        ierr = putvar(ncout , id_varout(1), rdummymean ,  jk, ikx, iky, ktime=jt )
        ierr = putvar(ncout , id_varout(2), DBLE(nflsdc), jk, ikx, iky, ktime=jt )  
     END DO ! k-loop
  END DO    ! time-loop

  ierr = closeout(ncout )

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
    ipk(:) = npk ! mean for each level

    ierr=getvaratt (cf_in  , cv_nam, clunits, zspval, cllong_name, clshort_name)

    ! define new variables for output 
    stypvar(1)%rmissing_value    = 99999.
    stypvar(1)%scale_factor      = 1.
    stypvar(1)%add_offset        = 0.
    stypvar(1)%savelog10         = 0.
    stypvar(1)%conline_operation = 'N/A'

    stypvar(1)%cunits         = TRIM(clunits)
    stypvar(1)%valid_min      = -100000.
    stypvar(1)%valid_max      =  100000.
    stypvar(1)%cname          = 'degraded_'//TRIM(cv_nam)
    stypvar(1)%clong_name     = 'degraded '//TRIM(cllong_name)
    stypvar(1)%cshort_name    = 'degraded_'//TRIM(clshort_name)

    stypvar(2)%rmissing_value    = 99999.
    stypvar(2)%scale_factor      = 1.
    stypvar(2)%add_offset        = 0.
    stypvar(2)%savelog10         = 0.
    stypvar(2)%conline_operation = 'N/A'

    stypvar(2)%cunits         = '[1]'
    stypvar(2)%valid_min      = 0.
    stypvar(2)%valid_max      = 100000.
    stypvar(2)%cname          = 'flsdc'
    stypvar(2)%clong_name     = 'fluid sub-domain counter'
    stypvar(2)%cshort_name    = 'fdsdc'

    IF ( nvpk == 2 ) THEN
       stypvar(:)%caxis       = 'TYX'
    ELSE
       stypvar(:)%caxis       = 'TZYX'
    END IF

    ! create output fileset
    IF ( cf_out == '' ) cf_out= 'degraded_'//TRIM(cv_nam)//'.nc'  ! set default name if not specified on command line
    ncout = create      (cf_out ,     'none',  ikx,   iky,   npk, cdep=cv_dep)
    ierr  = createvar   (ncout ,      stypvar , 2 , ipk,   id_varout  , cdglobal=TRIM(cglobal) ) 
    ierr  = putheadervar(ncout ,      cf_in  ,  ikx, iky, npk, pnavlon=rdumlon, pnavlat=rdumlat, pdep=gdep, cdep=cv_dep)
    tim   = getvar1d(cf_in  , cn_vtimec, npt)
    ierr  = putvar1d(ncout,  tim,       npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfdegrad
