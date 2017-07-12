PROGRAM cdflap
  !!======================================================================
  !!                     ***  PROGRAM  cdflap  ***
  !!=====================================================================
  !!  ** Purpose : Compute the laplacian of a field given in argument
  !!
  !!  ** Method  : Use a classical centered stencil
  !!               Gradient along the coasline are set to 0
  !!
  !! History :  3.0  : 09/2014  : J.M. Molines 
  !!         :  4.0  : 03/2017  : J.M. Molines  
  !!          
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class derived_fields
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                              :: ji, jj, jk, jt, jv ! dummy loop index
  INTEGER(KIND=4)                              :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                              :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                              :: npkv               ! vertical size of the variable
  INTEGER(KIND=4)                              :: narg, iargc, ijarg ! browse line
  INTEGER(KIND=4)                              :: ncout              ! ncid of output file
  INTEGER(KIND=4)                              :: nvars              ! nmber of variables in input file
  INTEGER(KIND=4)                              :: ierr               ! error status
  INTEGER(KIND=4)                              :: iioff1, iioff2     ! i-offset for mask
  INTEGER(KIND=4)                              :: ijoff1, ijoff2     ! j-offset for mask
  INTEGER(KIND=4)                              :: ii1, ii2, ij1, ij2 ! working index
  INTEGER(KIND=4), DIMENSION(1)                :: ipk, id_varout     ! levels and varid's of output vars
  INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: ipkin              ! ipk of all variables in input files

  REAL(KIND=4)                                 :: rmissing           ! missing value or_FillValue of input variable
  REAL(KIND=4)                                 :: grav=9.81          ! Gravity acceleration (m2/s)
  REAL(KIND=4)                                 :: omega              ! earth rotation rate
  REAL(KIND=4)                                 :: rpi                ! 3.14159...
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e1_i1, e1_i2       ! along i horizontal metric
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e2_j1, e2_j2       ! along j horizontal metric
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rmski, rmskj       ! relevant mask
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: v2d                ! input variable
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: ff                 ! Coriolis
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: gphi               ! Coriolis

  REAL(KIND=8)                                 :: dspval=99999.d0    ! output laplacian
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dtim               ! time counter
  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: dlap               ! output laplacian

  CHARACTER(LEN=256)                           :: cf_in              ! input file name
  CHARACTER(LEN=256)                           :: cv_in              ! input variable name
  CHARACTER(LEN=256)                           :: cv_units           ! units of input variable name
  CHARACTER(LEN=256)                           :: cv_lat             ! name of latitude variable in hgr
  CHARACTER(LEN=3)                             :: ct_in              ! input variable type [ T U V F ] on C-grid
  CHARACTER(LEN=256)                           :: cf_out='lap.nc'    !output file name
  CHARACTER(LEN=256)                           :: cln_in             ! Long name of input variable
  CHARACTER(LEN=256)                           :: csn_in             ! Short name of input variable
  CHARACTER(LEN=256)                           :: cldum              ! dummy string
  CHARACTER(LEN=10)                            :: ce1_i1, ce1_i2     ! name of relevant horizontal i-metric
  CHARACTER(LEN=10)                            :: ce2_j1, ce2_j2     ! name of relevant horizontal i-metric
  CHARACTER(LEN=10)                            :: cmask_i, cmask_j   ! name of relevant mask variable
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_nam             ! array of var name

  TYPE(variable), DIMENSION(1)                 :: stypvar            ! output attributes
  TYPE(variable), DIMENSION(:), ALLOCATABLE    :: sdum               ! input attributes

  LOGICAL                                      :: lchk               ! missing files flag
  LOGICAL                                      :: l_overf2=.FALSE.   ! overf flag
  LOGICAL                                      :: l_metric=.TRUE.    ! use metric flag
  LOGICAL                                      :: lnc4    = .FALSE.  ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdflap -f IN-file -v IN-var -t IN-type [-overf2] [-nometric] ...'
     PRINT *,'               ...[-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the Laplacian of the variable IN-var in file IN-file. Assume'
     PRINT *,'       that the data are on a C-grid model (as NEMO).'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file : netcdf file in input'
     PRINT *,'       -v IN-var  : name of the variable to process '
     PRINT *,'       -t IN-TYPE : Position of the variable on the C-grid [ T U V F ]'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-overf2] : save laplacien/f/f*g (where f is the local coriolis '
     PRINT *,'            parameter, and g is the accelaration due to gravity --9.81 m/s2-- )'
     PRINT *,'            For the SSH field, this is a proxy for geostrophic vorticity'
     PRINT *,'       [-nometric] : compute laplacian without considering metrics '
     PRINT *,'       [-o OUT-file] : specify output file name instead of ',TRIM(cf_out)
     PRINT *,'            This option must be used after the -overf2 or -nometric option, as'
     PRINT *,'            output file name is redefined when using these options.'
     PRINT *,'       [-nc4 ] : Use netcdf4 output with chunking and deflation level 1..'
     PRINT *,'            This option is effective only if cdftools are compiled with'
     PRINT *,'            a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fhgr)," ",TRIM(cn_fzgr) ,' and ', TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : lap<var> (unit/m2)'
     PRINT *,'       if option -overf2 is used, netcdf file is lapoverf2.nc and '
     PRINT *,'       variable is lap<var>overf2'
     STOP 
  ENDIF

  ijarg = 1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1
     SELECT CASE (cldum )
     CASE ( '-f'       ) ; CALL getarg(ijarg,cf_in ) ; ijarg=ijarg+1
     CASE ( '-v'       ) ; CALL getarg(ijarg,cv_in ) ; ijarg=ijarg+1
     CASE ( '-t'       ) ; CALL getarg(ijarg,ct_in ) ; ijarg=ijarg+1
        !options
     CASE ( '-overf2'  ) ; l_overf2 = .TRUE.
        ;                  cf_out   = 'lapoverf2.nc'
     CASE ( '-nometric') ; l_metric = .FALSE.
        ;                  cf_out   = 'lapgrid.nc'
     CASE ( '-o'       ) ; CALL getarg(ijarg,cf_out) ; ijarg=ijarg+1
     CASE ( '-nc4'     ) ; lnc4     = .TRUE.
     CASE DEFAULT        ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.'; STOP 99
     END SELECT
  ENDDO
  PRINT *, ' TYP ', ct_in

  ! check if files exists
  lchk=.FALSE.
  IF ( l_metric ) THEN
     lchk = chkfile (cn_fhgr)
     lchk = chkfile (cn_fzgr) .OR. lchk
     lchk = chkfile (cn_fmsk) .OR. lchk
  ENDIF
  lchk = chkfile (cf_in  ) .OR. lchk
  IF ( lchk ) STOP 99 ! missing files

  npiglo = getdim(cf_in,cn_x)
  npjglo = getdim(cf_in,cn_y)
  npk    = getdim(cf_in,cn_z)
  npt    = getdim(cf_in,cn_t)

  !IF ( npk == 0 ) THEN 
  !  PRINT *, 'Input file with no vertical dimension'
  !  PRINT *, 'Set npk to 1 '
  !  npk = 1
  !ENDIF

  PRINT *, 'NPIGLO : ', npiglo
  PRINT *, 'NPJGLO : ', npjglo
  PRINT *, 'NPK    : ', npk
  PRINT *, 'NPT    : ', npT

  nvars = getnvar(cf_in)   ! get the number of variables in files
  ALLOCATE(cv_nam(nvars), ipkin(nvars) ,sdum(nvars) )

  ! get list of variable names and collect attributes in stypvar
  cv_nam(:) = getvarname(cf_in,nvars, sdum)

  ! determine the vertical size of the working variable
  ipkin(:)  = getipk (cf_in,nvars)
  DO jv =1, nvars
     IF ( cv_nam(jv) == cv_in ) THEN
        npkv=ipkin(jv)
        EXIT
     ENDIF
  ENDDO

  ierr = getvaratt (cf_in, cv_in, cv_units, rmissing, cln_in, csn_in)

  ! fix mesh mask variables according to variable type
  ! HOT ! need to write down stencil to fully understand those settings ...
  SELECT CASE (ct_in )
  CASE ( 'T' )
     ! needs umask, vmask, e1u, e1t, e2v, e2t 
     cmask_i = cn_umask  ; cmask_j = cn_vmask
     ce1_i1  = cn_ve1u   ; ce1_i2  = cn_ve1t
     ce2_j1  = cn_ve2v   ; ce2_j2  = cn_ve2t
     iioff1  = 0         ; iioff2  = 1
     ijoff1  = 0         ; ijoff2  = 1
     cv_lat  =cn_gphit
  CASE ( 'U' )
     ! needs tmask, fmask, e1t, e1u, e2f, e2u 
     cmask_i = cn_tmask  ; cmask_j = cn_fmask
     ce1_i1  = cn_ve1t   ; ce1_i2  = cn_ve1u
     ce2_j1  = cn_ve2f   ; ce2_j2  = cn_ve2u
     iioff1  = 1         ; iioff2  = 0
     ijoff1  = 0         ; ijoff2  = 1
     cv_lat  =cn_gphiu
  CASE ( 'V' )
     ! needs fmask, tmask, e1f, e1v, e2t, e2v 
     cmask_i = cn_fmask  ; cmask_j = cn_tmask
     ce1_i1  = cn_ve1f   ; ce1_i2  = cn_ve1v
     ce2_j1  = cn_ve2t   ; ce2_j2  = cn_ve2v
     iioff1  = 0         ; iioff2  = 1
     ijoff1  = 1         ; ijoff2  = 0
     cv_lat  =cn_gphiv
  CASE ( 'F' )
     ! needs vmask, umask, e1v, e1f, e2u, e2f 
     cmask_i = cn_vmask  ; cmask_j = cn_umask
     ce1_i1  = cn_ve1v   ; ce1_i2  = cn_ve1f
     ce2_j1  = cn_ve2u   ; ce2_j2  = cn_ve2f
     iioff1  = 1         ; iioff2  = 0
     ijoff1  = 1         ; ijoff2  = 0
     cv_lat  =cn_gphif
  CASE DEFAULT
     PRINT *, ' TYPE ', TRIM(ct_in),' unknown on C-grid'
     STOP 99
  END SELECT

  ! Allocate the memory
  ALLOCATE ( e1_i1(npiglo,npjglo), e2_j1(npiglo,npjglo) )
  ALLOCATE ( e1_i2(npiglo,npjglo), e2_j2(npiglo,npjglo) )
  ALLOCATE ( rmski(npiglo,npjglo), rmskj(npiglo,npjglo) )
  IF ( l_overf2 ) THEN
     ALLOCATE (  ff(npiglo,npjglo), gphi(npiglo,npjglo) )
  ENDIF

  ALLOCATE ( dlap(npiglo,npjglo), v2d(npiglo, npjglo)  )
  ALLOCATE ( dtim(npt) )

  ! Read the metrics from the mesh_hgr file
  IF ( l_metric ) THEN
     e1_i1 = getvar(cn_fhgr, ce1_i1, 1, npiglo, npjglo)
     e1_i2 = getvar(cn_fhgr, ce1_i2, 1, npiglo, npjglo)
     e2_j1 = getvar(cn_fhgr, ce2_j1, 1, npiglo, npjglo)
     e2_j2 = getvar(cn_fhgr, ce2_j2, 1, npiglo, npjglo)
  ELSE
     e1_i1 = 1.
     e1_i2 = 1.
     e2_j1 = 1.
     e2_j2 = 1.
  ENDIF
  !
  IF ( l_overf2 ) THEN
     ! to prepare computation of g/f*LAP(SSH) 
     gphi    = getvar(cn_fhgr, cv_lat, 1, npiglo, npjglo)
     rpi     = ACOS(-1.)
     omega   = 2 * rpi / 86400.
     ff     = 2.*omega* SIN (gphi*rpi/180.)
     WHERE (ff == 0 ) ff = EPSILON(1.0)
  ENDIF

  ! define new variables for output
  CALL CreateOutput

  ! Main time loop
  DO jt = 1, npt
     ! Main level loop from top to bottom
     DO jk = 1, npkv
        PRINT *,'jt = ', jt,' jk = ', jk

        ! variable at level jk time jt
        v2d(:,:) =  getvar(cf_in, cv_in, jk, npiglo, npjglo, ktime=jt)
        ! relevant mask
        IF ( l_metric ) THEN
           rmski = getvar(cn_fmsk, cmask_i, jk, npiglo, npjglo)
           rmskj = getvar(cn_fmsk, cmask_j, jk, npiglo, npjglo)
        ELSE
           ! something can be done with regard to field's  missing value
           rmski = 1.
           rmskj = 1.
        ENDIF
        ! Compute laplacian
        dlap(:,:) = 0.d0
        DO jj = 2, npjglo -1
           ij1 = jj + ijoff1
           ij2 = jj - ijoff2
           DO ji = 2, npiglo -1
              ii1 = ji + iioff1
              ii2 = ji - iioff2
              dlap(ji,jj) = ((v2d(ji+1,jj  )-v2d(ji  ,jj  ))*1.d0/e1_i1(ii1,jj)*rmski(ii1,jj)                &
                   &       - (v2d(ji  ,jj  )-v2d(ji-1,jj  ))*1.d0/e1_i1(ii2,jj)*rmski(ii2,jj))/e1_i2(ji,jj)  &
                   &       +((v2d(ji  ,jj+1)-v2d(ji  ,jj  ))*1.d0/e2_j1(ji,ij1)*rmskj(ji,ij1)                &
                   &       - (v2d(ji  ,jj  )-v2d(ji  ,jj-1))*1.d0/e2_j1(ji,ij2)*rmskj(ji,ij2))/e2_j2(ji,jj)
           END DO
        END DO

        IF ( l_overf2 ) THEN
           WHERE (dlap == 0.d0 )
              dlap = dspval
           ELSEWHERE
              dlap = grav*dlap /ff/ff
           ENDWHERE
        ENDIF

        ! write level jk 
        ierr = putvar(ncout, id_varout(1), REAL(dlap), jk, npiglo, npjglo, ktime=jt)

     END DO  ! loop to next level
  END DO ! loop on time

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
    ipk(1) = npkv
    IF ( l_overf2) THEN ; stypvar(1)%cname      = 'lap'//TRIM(cv_in)//'overf2'
       ;                  stypvar(1)%cunits     = TRIM(cv_units)//'/m'
       ;                  stypvar(1)%clong_name = 'Laplacian of '//TRIM(cln_in)
    ELSE                ; stypvar(1)%cname      = 'lap'//TRIM(cv_in)
       ;                  stypvar(1)%cunits     = TRIM(cv_units)//'/m2'
       ;                  stypvar(1)%clong_name = 'g /f/f* Laplacian of '//TRIM(cln_in)
    ENDIF

    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%rmissing_value    = dspval
    stypvar(1)%valid_min         = -10.
    stypvar(1)%valid_max         = 10.
    stypvar(1)%cshort_name       = stypvar(1)%cname
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TZYX'

    ! create output fileset
    ncout = create      (cf_out, cf_in, npiglo, npjglo, npk         , ld_nc4=lnc4 ) 
    ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_in, npiglo, npjglo, npk         )

    dtim = getvar1d(cf_in , cn_vtimec, npt     )
    ierr = putvar1d(ncout , dtim     , npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdflap

