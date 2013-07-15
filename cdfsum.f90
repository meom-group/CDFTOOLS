PROGRAM cdfsum
  !!======================================================================
  !!                     ***  PROGRAM  cdfsum  ***
  !!=====================================================================
  !!  ** Purpose : Compute the sum of a variable over the ocean, or
  !!               part of the ocean
  !!
  !!  ** Method  : this code is for partial steps configuration
  !!              sum = sum ( V * e1 *e2 * e3 *mask )
  !!              CAUTION : this version is still tricky, as it does not
  !!              compute the same thing in case of forcing field or
  !!              model field. Need clarification ( JMM)
  !!
  !! History : 2.1  : 11/2008  : P. Mathiot   : Original code (from cdfmean)
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
  INTEGER(KIND=4)                           :: ik                  ! dummy loop index
  INTEGER(KIND=4)                           :: iimin=0, iimax=0    ! domain limitation for computation
  INTEGER(KIND=4)                           :: ijmin=0, ijmax=0    ! domain limitation for computation
  INTEGER(KIND=4)                           :: ikmin=0, ikmax=0    ! domain limitation for computation
  INTEGER(KIND=4)                           :: ierr                ! working integer
  INTEGER(KIND=4)                           :: narg, iargc         ! command line 
  INTEGER(KIND=4)                           :: npiglo, npjglo      ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt            ! size of the domain
  INTEGER(KIND=4)                           :: nvpk                ! vertical levels in working variable
  INTEGER(KIND=4)                           :: numout=10           ! logical unit
  INTEGER(KIND=4)                           :: ncout               ! for netcdf output
  INTEGER(KIND=4), DIMENSION(2)             :: ipk, id_varout

  REAL(KIND=4)                              :: zspval             ! missing value
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1, e2, e3,  zv     ! metrics, velocity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask               ! npiglo x npjglo
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: gdep                ! depth 
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                 ! time
  REAL(KIND=4), DIMENSION(1,1)              :: rdumlon, rdumlat    ! dummy latitude and longitude
  REAL(KIND=4), DIMENSION(1,1)              :: rdummy              ! dummy 2d variable for result

  REAL(KIND=8)                              :: dvol, dvol2d        ! volume of the ocean/ layer
  REAL(KIND=8)                              :: dsurf               ! surface of the ocean
  REAL(KIND=8)                              :: dsum, dsum2d        ! global sum /layer sum
  REAL(KIND=8)                              :: dsumt               ! global sum over time

  CHARACTER(LEN=256)                        :: cldum               ! dummy string
  CHARACTER(LEN=256)                        :: cf_in               ! file name 
  CHARACTER(LEN=256)                        :: cf_out='cdfsum.nc'  ! output file name 
  CHARACTER(LEN=256)                        :: cv_dep              ! depth name
  CHARACTER(LEN=256)                        :: cv_in               ! variable name
  CHARACTER(LEN=20)                         :: cv_e1, cv_e2, cv_e3 ! name of the horiz/vert metrics
  CHARACTER(LEN=20)                         :: cv_msk              ! name of mask variable
  CHARACTER(LEN=20)                         :: cvartype            ! variable type
  CHARACTER(LEN=256)                        :: clunits            ! attribute of output file : units
  CHARACTER(LEN=256)                        :: cllong_name        !     "      long name
  CHARACTER(LEN=256)                        :: clshort_name       !     "      short name
  CHARACTER(LEN=256)                        :: cglobal            !     "      global 

  TYPE(variable), DIMENSION(2)              :: stypvar             ! structure of output

  LOGICAL                                   :: lforcing            ! forcing flag
  LOGICAL                                   :: lchk                ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfsum IN-file IN-var T| U | V | F | W  ... '
     PRINT *,'             ... [imin imax jmin jmax kmin kmax] [-full ] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Computes the sum value of the field (3D, weighted)' 
     PRINT *,'       This sum can be optionally limited to a sub-area.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       IN-file : netcdf input file.' 
     PRINT *,'       IN-var  : netcdf variable to work with.'
     PRINT *,'       T| U | V | F | W : C-grid point where IN-var is located.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [imin imax jmin jmax kmin kmax] : limit of the sub area to work with.' 
     PRINT *,'              if imin=0 all i are taken'
     PRINT *,'              if jmin=0 all j are taken'
     PRINT *,'              if kmin=0 all k are taken'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'      ', TRIM(cn_fhgr),', ',TRIM(cn_fzgr),' and ',TRIM(cn_fmsk) 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       Standard output.'
     PRINT *,'       netcdf file : ',TRIM(cf_out),' with 2 variables : vertical profile of sum'
     PRINT *,'                     and 3D sum.'
     STOP
  ENDIF

  CALL getarg (1, cf_in)
  CALL getarg (2, cv_in)
  CALL getarg (3, cvartype)

  lchk = chkfile(cn_fhgr)
  lchk = chkfile(cn_fzgr) .OR. lchk
  lchk = chkfile(cn_fmsk) .OR. lchk
  lchk = chkfile(cf_in  ) .OR. lchk
  IF ( lchk ) STOP ! missing file

  IF (narg > 3 ) THEN
     IF ( narg /= 9 ) THEN
        PRINT *, ' ERROR : You must give 6 optional values (imin imax jmin jmax kmin kmax)'
        STOP
     ELSE
        ! input optional iimin iimax ijmin ijmax
        CALL getarg ( 4,cldum) ; READ(cldum,*) iimin
        CALL getarg ( 5,cldum) ; READ(cldum,*) iimax
        CALL getarg ( 6,cldum) ; READ(cldum,*) ijmin
        CALL getarg ( 7,cldum) ; READ(cldum,*) ijmax
        CALL getarg ( 8,cldum) ; READ(cldum,*) ikmin
        CALL getarg ( 9,cldum) ; READ(cldum,*) ikmax
     ENDIF
  ENDIF

  npiglo = getdim (cf_in,cn_x)
  npjglo = getdim (cf_in,cn_y)
  npk    = getdim (cf_in,cn_z)
  nvpk   = getvdim(cf_in,cv_in)
  npt    = getdim (cf_in,cn_t)

  IF (iimin /= 0 ) THEN ; npiglo = iimax - iimin + 1;  ELSE ; iimin = 1 ;  ENDIF
  IF (ijmin /= 0 ) THEN ; npjglo = ijmax - ijmin + 1;  ELSE ; ijmin = 1 ;  ENDIF
  IF (ikmin /= 0 ) THEN ; npk    = ikmax - ikmin + 1;  ELSE ; ikmin = 1 ;  ENDIF

  IF (nvpk == 2 ) nvpk = 1
  IF (nvpk == 3 ) nvpk = npk

  PRINT *, 'Size of the extracted area :'
  PRINT *, '  npiglo = ', npiglo
  PRINT *, '  npjglo = ', npjglo
  PRINT *, '  npk    = ', npk
  PRINT *, '  nvpk   = ', nvpk
  PRINT *, '  npt    = ', npt

  lforcing=.FALSE.
  IF ( (npk == 0) ) THEN
     lforcing = .TRUE.
     npk      = 1
     PRINT *, 'W A R N I N G : you used a forcing field'
  END IF

  IF (lforcing)  OPEN(unit=numout, file='cdfsum.txt' , form='formatted', status='new', iostat=ierr)

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo) )
  ALLOCATE ( zv   (npiglo,npjglo) )
  ALLOCATE ( e1   (npiglo,npjglo), e2(npiglo,npjglo), e3(npiglo,npjglo) )
  ALLOCATE ( gdep (npk), tim(npt) )

  SELECT CASE (TRIM(cvartype))
  CASE ( 'T' )
     cv_e1  = cn_ve1t
     cv_e2  = cn_ve2t
     cv_e3  = 'e3t_ps'
     cv_msk = 'tmask'
     cv_dep = cn_gdept
  CASE ( 'U' )
     cv_e1  = cn_ve1u
     cv_e2  = cn_ve2u
     cv_e3  = 'e3t_ps'
     cv_msk = 'umask'
     cv_dep = cn_gdept
  CASE ( 'V' )
     cv_e1  = cn_ve1v
     cv_e2  = cn_ve2v
     cv_e3  = 'e3t_ps'
     cv_msk = 'vmask'
     cv_dep = cn_gdept
  CASE ( 'F' )
     cv_e1  = cn_ve1f
     cv_e2  = cn_ve2f
     cv_e3  = 'e3t_ps'
     cv_msk = 'fmask'
     cv_dep = cn_gdept
  CASE ( 'W' )
     cv_e1  = cn_ve1t
     cv_e2  = cn_ve2t
     cv_e3  = 'e3w_ps'
     cv_msk = 'tmask'
     cv_dep = cn_gdepw
  CASE DEFAULT
     PRINT *, 'this type of variable is not known :', TRIM(cvartype)
     STOP
  END SELECT

  e1(:,:) = getvar  (cn_fhgr, cv_e1, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin)
  e2(:,:) = getvar  (cn_fhgr, cv_e2, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin)
  gdep(:) = getvare3(cn_fzgr, cv_dep,   npk                                   )

 rdumlon = 0. ; rdumlat = 0.
 ipk(1) = nvpk  ! vertical profile
 ipk(2) = 1     ! 3D sum
  ierr=getvaratt (cf_in, cv_in, clunits, zspval, cllong_name, clshort_name)

  ! define new variables for output 
  stypvar%rmissing_value    = 99999.
  stypvar%valid_min         = -1000.
  stypvar%valid_max         = 1000.
  stypvar%scale_factor      = 1.
  stypvar%add_offset        = 0.
  stypvar%savelog10         = 0.
  stypvar%conline_operation = 'N/A'

  stypvar(1)%cname          = 'sum_'//TRIM(cv_in)
  stypvar(1)%cunits         = TRIM(clunits)//'.m2'
  stypvar(1)%clong_name     = 'sum'//TRIM(cllong_name)
  stypvar(1)%cshort_name    = 'sum'//TRIM(clshort_name)
  stypvar(1)%caxis          = 'ZT'

  stypvar(2)%cname          = 'sum_3D'//TRIM(cv_in)
  stypvar(2)%cunits         = TRIM(clunits)//'.m3'
  stypvar(2)%clong_name     = 'sum_3D'//TRIM(cllong_name)
  stypvar(2)%cshort_name    = 'sum_3D'//TRIM(clshort_name)
  stypvar(2)%caxis          = 'T'

  ncout = create      (cf_out,     'none',  1,     1  ,   nvpk, cdep=cv_dep)
  ierr  = createvar   (ncout,      stypvar, 2    , ipk,   id_varout        )
  ierr  = putheadervar(ncout,      cf_in,   1,     1, npk, pnavlon=rdumlon, pnavlat=rdumlat, pdep=gdep(1:nvpk), cdep=cv_dep)
  tim   = getvar1d(cf_in, cn_vtimec, npt)
  ierr  = putvar1d(ncout,  tim,      npt, 'T')



  dsumt = 0.d0
  DO jt = 1,npt
     dvol = 0.d0
     dsum = 0.d0
     zv   = 0.
     DO jk = 1,nvpk
        ik = jk + ikmin -1
        ! Get velocities v at ik
        zv   (:,:) = getvar(cf_in,   cv_in,  ik, npiglo, npjglo, ktime=jt,   kimin=iimin, kjmin=ijmin)
        zmask(:,:) = getvar(cn_fmsk, cv_msk, ik, npiglo, npjglo,             kimin=iimin, kjmin=ijmin)
        !    zmask(:,npjglo)=0.

        ! get e3 at level ik ( ps...)
        e3(:,:) = getvar(cn_fzgr, cv_e3, ik, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ldiom=.TRUE.)      
        ! 
        IF (.NOT. lforcing) THEN
           dsurf  = SUM(DBLE(e1 * e2      * zmask))
           dvol2d = SUM(DBLE(e1 * e2 * e3 * zmask))
           dvol   = dvol + dvol2d
           dsum2d = SUM(DBLE(zv))
           dsum   = dsum + dsum2d
           IF (dvol2d /= 0 )THEN
              PRINT *, ' Sum value at level ', ik, '(',gdep(ik),' m) ', dsum2d
              rdummy(1,1)= REAL(dsum2d)
           ELSE
              PRINT *, ' No points in the water at level ', ik, '(',gdep(ik),' m) '
              rdummy(1,1)= 99999.
           ENDIF
        ELSE
           dsurf  = SUM(DBLE(     e1 * e2 * zmask))
           dsum2d = SUM(DBLE(zv * e1 * e2 * zmask))
           dsum   = dsum + dsum2d
           PRINT *, ' Sum value at time ',jt,' = ', dsum2d
           PRINT *, '          Surface  = ', dsurf/1.d6,' km^2'
           PRINT *, '       mean value  = ', dsum2d/dsurf
           WRITE (numout,'(i4," ",1e12.6)') jt, dsum2d
           rdummy(1,1) = REAL(dsum2d)
        END IF
        ierr = putvar( ncout, id_varout(1), rdummy, jk, 1,1, ktime=jt)
     END DO
     dsumt = dsumt + dsum
     IF (.NOT. lforcing) PRINT * ,' Sum value over the ocean: ', dsumt
     rdummy(1,1) = REAL(dsumt)
     ierr = putvar( ncout, id_varout(2), rdummy, 1, 1, 1, ktime=jt)
  END DO  ! time loop
  
  PRINT *, ' mean Sum over time ', dsumt/npt

  CLOSE(numout)
  ierr=closeout(ncout)

END PROGRAM cdfsum
