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
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class integration
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: jk, jt              ! dummy loop index
  INTEGER(KIND=4)                           :: ik, it              ! working integer
  INTEGER(KIND=4)                           :: iimin=0, iimax=0    ! domain limitation for computation
  INTEGER(KIND=4)                           :: ijmin=0, ijmax=0    ! domain limitation for computation
  INTEGER(KIND=4)                           :: ikmin=0, ikmax=0    ! domain limitation for computation
  INTEGER(KIND=4)                           :: ierr                ! working integer
  INTEGER(KIND=4)                           :: narg, iargc,ijarg   ! command line 
  INTEGER(KIND=4)                           :: npiglo, npjglo      ! size of the domain
  INTEGER(KIND=4)                           :: npk, npkk, npt      ! size of the domain
  INTEGER(KIND=4)                           :: nvpk                ! vertical levels in working variable
  INTEGER(KIND=4)                           :: numout=10           ! logical unit
  INTEGER(KIND=4)                           :: ncout               ! for netcdf output
  INTEGER(KIND=4), DIMENSION(2)             :: ipk, id_varout

  REAL(KIND=4)                              :: zspval             ! missing value
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: gdep                ! depth 
  REAL(KIND=4), DIMENSION(1,1)              :: rdummy              ! dummy 2d variable for result
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1, e2, e3,  zv     ! metrics, velocity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask               ! npiglo x npjglo

  REAL(KIND=8)                              :: dvol, dvol2d        ! volume of the ocean/ layer
  REAL(KIND=8)                              :: dsurf               ! surface of the ocean
  REAL(KIND=8)                              :: dsum, dsum2d        ! global sum /layer sum
  REAL(KIND=8)                              :: dsumt               ! global sum over time
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim                ! time

  CHARACTER(LEN=256)                        :: cldum               ! dummy string
  CHARACTER(LEN=256)                        :: cf_in               ! file name 
  CHARACTER(LEN=256)                        :: cf_e3               ! name of vertical metrics file
  CHARACTER(LEN=256)                        :: cf_out='cdfsum.nc'  ! output file name 
  CHARACTER(LEN=256)                        :: cv_dep              ! depth name in mesh_zgr
  CHARACTER(LEN=256)                        :: cdep                ! depth name in output file
  CHARACTER(LEN=256)                        :: cv_in               ! variable name
  CHARACTER(LEN=20)                         :: cv_e1, cv_e2, cv_e3 ! name of the horiz/vert metrics
  CHARACTER(LEN=20)                         :: cv_msk              ! name of mask variable
  CHARACTER(LEN=20)                         :: cl_vmsk             ! name of external mask variable (-M option)
  CHARACTER(LEN=20)                         :: cvartype            ! variable type
  CHARACTER(LEN=256)                        :: clunits            ! attribute of output file : units
  CHARACTER(LEN=256)                        :: cllong_name        !     "      long name
  CHARACTER(LEN=256)                        :: clshort_name       !     "      short name
  CHARACTER(LEN=256)                        :: cglobal            !     "      global 

  TYPE(variable), DIMENSION(2)              :: stypvar             ! structure of output

  LOGICAL                                   :: lforcing            ! forcing flag
  LOGICAL                                   :: lchk  = .FALSE.     ! flag for missing files
  LOGICAL                                   :: lnc4  = .FALSE.     ! flag for netcdf4 
  LOGICAL                                   :: lerror= .FALSE.     ! flag for missing arguments
  LOGICAL                                   :: lfmsk = .FALSE.     ! flag for using non standard mask file
  LOGICAL                                   :: lnwgh = .FALSE.     ! flag for no weight
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfsum -f IN-file -v IN-var -p C-type [-nowght]... '
     PRINT *,'          ... [-w imin imax jmin jmax kmin kmax] [-full ] ...'
     PRINT *,'          ... [-o OUT-file] [-nc4] [-M MSK-file VAR-mask ] [-vvl] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the sum value of the field (3D weighted, unless -nowght used).'
     PRINT *,'       This sum can be optionally limited to a 3D sub-area.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file : netcdf input file.' 
     PRINT *,'       -v IN-var  : netcdf variable to work with.'
     PRINT *,'       -p C-type  : one of T|U|V|F|W indicating rhe  C-grid point where'
     PRINT *,'               IN-var is located.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-w imin imax jmin jmax kmin kmax]: set the 3D window limiting sub area.'
     PRINT *,'              if imin=0 all i are taken'
     PRINT *,'              if jmin=0 all j are taken'
     PRINT *,'              if kmin=0 all k are taken'
     PRINT *,'       [-nowght] : compute the aritmetic sum of the variable over the grid' 
     PRINT *,'               points, (not the area/volume weighted sum ). ' 
     PRINT *,'       [-full ]:  Use full steps instead of default partial steps'
     PRINT *,'       [-vvl  ]:  use time -varying  vertical metrics'
     PRINT *,'       [-o OUT-file ] : name of the output file instead of', TRIM(cf_out)
     PRINT *,'       [-nc4 ] : use netcdf4 chunking and deflation.'
     PRINT *,'       [-M MSK-file VAR-mask] : Allow the use of a non standard mask file '
     PRINT *,'              with VAR-mask, instead of ',TRIM(cn_fmsk),' and the variable'
     PRINT *,'              associated with the grid point set by -p argument.'
     PRINT *,'              This option is a usefull alternative to the -w option, when the '
     PRINT *,'              area of interest is not ''box-like''. However, for vertical '
     PRINT *,'              selection, both -w and -M can be used together.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'      ', TRIM(cn_fhgr),', ',TRIM(cn_fzgr),' and ',TRIM(cn_fmsk),'. If'
     PRINT *,'         -M option is used, the specified mask file is required instead '
     PRINT *,'         ', TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       Standard output.'
     PRINT *,'       netcdf file : ',TRIM(cf_out),' unless modified with -o option. '
     PRINT *,'           - 2 variables : vertical profile of sum and 3D sum.'
     PRINT *,'                         names are sum_<varname> and sum3D_<varname>.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO: '
     PRINT *,'       cdfmean '
     STOP 
  ENDIF

  ijarg=1
  ! mandatory arguments are set to none by default for further check
  cf_in='none' ; cv_in='none'; cvartype='none'
  DO WHILE ( ijarg <= narg )
     CALL getarg( ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE (cldum)
     CASE ( '-f '   ) ; CALL getarg(ijarg, cf_in   ) ; ijarg=ijarg+1
     CASE ( '-v '   ) ; CALL getarg(ijarg, cv_in   ) ; ijarg=ijarg+1
     CASE ( '-p '   ) ; CALL getarg(ijarg, cvartype) ; ijarg=ijarg+1
        ! options
     CASE ( '-o '   ) ; CALL getarg(ijarg, cf_out  ) ; ijarg=ijarg+1
     CASE ( '-vvl ' ) ; lg_vvl=.TRUE.
     CASE ( '-nc4 ' ) ; lnc4  =.TRUE.
     CASE ( '-nowght ' ) ; lnwgh  =.TRUE.
     CASE ( '-w'    ) ; 
        ;               CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*)  iimin
        ;               CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*)  iimax
        ;               CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*)  ijmin
        ;               CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*)  ijmax
        ;               CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*)  ikmin
        ;               CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*)  ikmax
     CASE ( '-M'    ) ;    lfmsk =.TRUE.
        ;               CALL getarg(ijarg, cn_fmsk ) ; ijarg=ijarg+1 
        ;               CALL getarg(ijarg, cl_vmsk ) ; ijarg=ijarg+1 
     CASE DEFAULT     ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO
  
  IF ( cf_in == 'none' )  THEN
     PRINT *,' You must specify an input file with -f option'
     lerror= lerror .OR. .TRUE. 
  ENDIF
  IF ( cv_in == 'none' )  THEN
     PRINT *,' You must specify an input variable with -v option'
     lerror= lerror .OR. .TRUE. 
  ENDIF
  IF ( cvartype == 'none' )  THEN
     PRINT *,' You must specify a point type with -p option'
     lerror= lerror .OR. .TRUE. 
  ENDIF
  IF (lerror ) STOP 99

  IF ( .NOT. lnwgh )  lchk = chkfile(cn_fhgr)
  IF ( .NOT. lnwgh )  lchk = chkfile(cn_fzgr) .OR. lchk
  lchk = chkfile(cn_fmsk) .OR. lchk
  lchk = chkfile(cf_in  ) .OR. lchk
  IF ( lchk ) STOP 99 ! missing file

  npiglo = getdim (cf_in,cn_x )
  npjglo = getdim (cf_in,cn_y )
  npk    = getdim (cf_in,cn_z )
  nvpk   = getvdim(cf_in,cv_in)
  npt    = getdim (cf_in,cn_t )
  npkk = npk

  IF (iimin /= 0 ) THEN ; npiglo = iimax - iimin + 1;  ELSE ; iimin = 1 ;  ENDIF
  IF (ijmin /= 0 ) THEN ; npjglo = ijmax - ijmin + 1;  ELSE ; ijmin = 1 ;  ENDIF
  IF (ikmin /= 0 ) THEN ; npkk   = ikmax - ikmin + 1;  ELSE ; ikmin = 1 ; ikmax = npk ;  ENDIF

  IF (nvpk == 2 ) nvpk = 1
  IF (nvpk == 3 ) nvpk = npkk

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
  ALLOCATE ( gdep (npk), dtim(npt) )

  SELECT CASE (TRIM(cvartype))
  CASE ( 'T' )
     cv_e1  = cn_ve1t
     cv_e2  = cn_ve2t
     cf_e3  = cn_fe3t
     cv_e3  = cn_ve3t
     cv_msk = cn_tmask
     cv_dep = cn_gdept
     cdep   = cn_vdeptht
     IF (lg_vvl) cv_e3 = cn_ve3tvvl
  CASE ( 'U' )
     cv_e1  = cn_ve1u
     cv_e2  = cn_ve2u
     cf_e3  = cn_fe3u
     cv_e3  = cn_ve3u
     cv_msk = cn_umask
     cv_dep = cn_gdept
     cdep   = cn_vdepthu
     IF (lg_vvl) cv_e3 = cn_ve3uvvl
  CASE ( 'V' )
     cv_e1  = cn_ve1v
     cv_e2  = cn_ve2v
     cf_e3  = cn_fe3v
     cv_e3  = cn_ve3v
     cv_msk = cn_vmask
     cv_dep = cn_gdept
     cdep   = cn_vdepthv
     IF (lg_vvl) cv_e3 = cn_ve3vvvl
  CASE ( 'F' )
     cv_e1  = cn_ve1f
     cv_e2  = cn_ve2f
     cf_e3  = cn_fe3t
     cv_e3  = cn_ve3t  !!!! CAUTION : e3f not defined
     cv_msk = cn_fmask
     cv_dep = cn_gdept
     cdep   = cn_vdeptht
     IF (lg_vvl) cv_e3 = cn_ve3tvvl
  CASE ( 'W' )
     cv_e1  = cn_ve1t
     cv_e2  = cn_ve2t
     cf_e3  = cn_fe3w
     cv_e3  = cn_ve3w
     cv_msk = cn_tmask
     cv_dep = cn_gdepw
     cdep   = cn_vdepthw
     IF (lg_vvl) cv_e3 = cn_ve3wvvl
  CASE DEFAULT
     PRINT *, 'this type of variable is not known :', TRIM(cvartype)
     STOP 99
  END SELECT

  IF ( lg_vvl ) cf_e3 = cf_in

  ! set cv_mask to on-line specified name if -M option used
  IF ( lfmsk ) cv_msk = cl_vmsk

  IF ( lnwgh)  THEN
     e1(:,:) = 1
     e2(:,:) = 1
  ELSE
     e1(:,:) = getvar  (cn_fhgr, cv_e1, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin)
     e2(:,:) = getvar  (cn_fhgr, cv_e2, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin)
  ENDIF
  gdep(:) = getvare3(cn_fzgr, cv_dep,   npk                                   )

  CALL CreateOutput

  dsumt = 0.d0
  DO jt = 1,npt
     IF ( lg_vvl ) THEN ; it = jt
     ELSE               ; it = 1
     ENDIF
     dvol = 0.d0
     dsum = 0.d0
     zv   = 0.
     DO jk = 1,nvpk
        ik = jk + ikmin -1
        ! Get field  at ik
        zv   (:,:) = getvar(cf_in,   cv_in,  ik, npiglo, npjglo, ktime=jt,   kimin=iimin, kjmin=ijmin)
        zmask(:,:) = getvar(cn_fmsk, cv_msk, ik, npiglo, npjglo,             kimin=iimin, kjmin=ijmin)
        !    zmask(:,npjglo)=0.
        zv = zv * zmask

        ! get e3 at level ik ( ps...)
        IF ( lnwgh ) THEN
           e3(:,:) = 1.
        ELSE
           e3(:,:) = getvar(cf_e3, cv_e3, ik, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=it, ldiom=.NOT.lg_vvl)      
        ENDIF
        ! 
        IF (.NOT. lforcing) THEN
           dsurf  = SUM(DBLE(e1 * e2      * zmask))
           dvol2d = SUM(DBLE(e1 * e2 * e3 * zmask))
           dvol   = dvol + dvol2d
           dsum2d = SUM(DBLE(zv * e1 * e2 ))
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
           WRITE (numout,'(i4," ",1e13.6)') jt, dsum2d
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

  CONTAINS
  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose : Create netcdf output file  
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(1,1)              :: zdumlon, zdumlat    ! dummy latitude and longitude


    !!
    ! define new variables for output 
    zdumlon = 0. ; zdumlat = 0.
    ierr=getvaratt (cf_in, cv_in, clunits, zspval, cllong_name, clshort_name)
    stypvar%rmissing_value    = 99999.
    stypvar%valid_min         = -1000.
    stypvar%valid_max         = 1000.
    stypvar%scale_factor      = 1.
    stypvar%add_offset        = 0.
    stypvar%savelog10         = 0.
    stypvar%conline_operation = 'N/A'

    ipk(1) = nvpk  ! vertical profile
    stypvar(1)%ichunk         = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname          = 'sum_'//TRIM(cv_in)
    IF ( lnwgh ) THEN
      stypvar(1)%cunits         = TRIM(clunits)
    ELSE
      stypvar(1)%cunits         = TRIM(clunits)//'.m2'
    ENDIF
    stypvar(1)%clong_name     = 'sum'//TRIM(cllong_name)
    stypvar(1)%cshort_name    = 'sum'//TRIM(clshort_name)
    stypvar(1)%caxis          = 'ZT'

    ipk(2) = 1     ! 3D sum
    stypvar(2)%ichunk         = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname          = 'sum_3D'//TRIM(cv_in)
    stypvar(2)%cunits         = TRIM(clunits)//'.m3'
    stypvar(2)%clong_name     = 'sum_3D'//TRIM(cllong_name)
    stypvar(2)%cshort_name    = 'sum_3D'//TRIM(clshort_name)
    stypvar(2)%caxis          = 'T'

    ncout = create      (cf_out,     'none',  1,     1  ,   nvpk, cdep=cdep , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,      stypvar, 2    , ipk,   id_varout       , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,      cf_in,   1,     1, npkk,              &
                    &  pnavlon=zdumlon, pnavlat=zdumlat,                   &
                    &  pdep=gdep(ikmin:ikmax),                             &
                    &  cdep=cdep                                           )
    dtim  = getvar1d(cf_in, cn_vtimec, npt     )
    ierr  = putvar1d(ncout,  dtim,     npt, 'T')
  
  END SUBROUTINE CreateOutput

END PROGRAM cdfsum
