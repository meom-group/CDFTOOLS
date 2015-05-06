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
  !!          
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

  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e1_i1, e1_i2       ! along i horizontal metric
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e2_j1, e2_j2       ! along j horizontal metric
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rmski, rmskj       ! relevant mask
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: v2d                ! input variable
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: ff                 ! Coriolis
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: gphi               ! Coriolis
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: tim                ! time counter
  REAL(KIND=4)                                 :: rmissing           ! missing value or_FillValue of input variable
  REAL(KIND=4)                                 :: grav=9.81          ! Gravity acceleration (m2/s)
  REAL(KIND=4)                                 :: omega              ! earth rotation rate
  REAL(KIND=4)                                 :: rpi                ! 3.14159...

  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: dlap               ! output laplacian
  REAL(KIND=8)                                 :: dspval=99999.d0    ! output laplacian

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
  LOGICAL                                      :: l_overf2=.FALSE.    ! overf flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg < 2 ) THEN
     PRINT *,' usage : cdflap IN-file IN-var  IN-type [-overf2]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the Laplacian of the variable IN-var in file IN-file'
     PRINT *,'       Assumes that the data are on a C-grid model (as NEMO) '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       IN-file : netcdf file in input'
     PRINT *,'       IN-var  : name of the variable to process '
     PRINT *,'       IN-TYPE : Position of the variable on the C-grid [ T U V F ]'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -overf2 : save laplacien/f/f*g (where f is the local coriolis '
     PRINT *,'            parameter, and g is the accelaration due to gravity --9.81 m/s2-- )'
     PRINT *,'            For the SSH field, this is a proxy for geostrophic vorticity'
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
  CALL getarg(ijarg, cf_in) ; ijarg = ijarg + 1
  CALL getarg(ijarg, cv_in) ; ijarg = ijarg + 1
  CALL getarg(ijarg, ct_in) ; ijarg = ijarg + 1
  IF ( narg == 4 ) THEN
  ! assume -overf2 option !
     l_overf2=.true.
     cf_out='lapoverf2.nc'
  ENDIF
  PRINT *, ' TYP ', ct_in

  ! check if files exists
  lchk = chkfile (cn_fhgr)
  lchk = chkfile (cn_fzgr) .OR. lchk
  lchk = chkfile (cn_fmsk) .OR. lchk
  lchk = chkfile (cf_in  ) .OR. lchk
  IF ( lchk ) STOP ! missing files

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
  ! determine the vertical size of the working variable
  nvars = getnvar(cf_in)   ! get the number of variables in files
  ALLOCATE(cv_nam(nvars), ipkin(nvars) ,sdum(nvars) )

  ! get list of variable names and collect attributes in stypvar
  cv_nam(:) = getvarname(cf_in,nvars, sdum)
  ipkin(:)  = getipk (cf_in,nvars)
  DO jv =1, nvars
   IF ( cv_nam(jv) == cv_in ) THEN
     npkv=ipkin(jv)
     EXIT
   ENDIF
  ENDDO

  ierr = getvaratt (cf_in, cv_in, cv_units, rmissing, cln_in, csn_in)

  ! define new variables for output
  ipk(1) = npkv
  IF ( l_overf2) THEN
     stypvar(1)%cname             = 'lap'//TRIM(cv_in)//'overf2'
  ELSE
     stypvar(1)%cname             = 'lap'//TRIM(cv_in)
  ENDIF
  stypvar(1)%cunits            = TRIM(cv_units)//'/m2'
  IF ( l_overf2) stypvar(1)%cunits  = TRIM(cv_units)//'/m'
  stypvar(1)%rmissing_value    = dspval
  stypvar(1)%valid_min         = -10.
  stypvar(1)%valid_max         = 10.
  stypvar(1)%clong_name        = 'Laplacian of '//TRIM(cln_in)
  IF ( l_overf2) stypvar(1)%clong_name        = 'g /f/f* Laplacian of '//TRIM(cln_in)
  stypvar(1)%cshort_name       = stypvar(1)%cname
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TZYX'

  ! fix mesh mask variables according to variable type
  ! HOT ! need to write down stencil to fully understand those settings ...
  SELECT CASE (ct_in )
  CASE ( 'T' )
      ! needs umask, vmask, e1u, e1t, e2v, e2t 
      cmask_i = 'umask'  ; cmask_j = 'vmask'
      ce1_i1  = cn_ve1u  ; ce1_i2  = cn_ve1t
      ce2_j1  = cn_ve2v  ; ce2_j2  = cn_ve2t
      iioff1  = 0        ; iioff2  = 1
      ijoff1  = 0        ; ijoff2  = 1
      cv_lat  =cn_gphit
  CASE ( 'U' )
      ! needs tmask, fmask, e1t, e1u, e2f, e2u 
      cmask_i = 'tmask'  ; cmask_j = 'fmask'
      ce1_i1  = cn_ve1t  ; ce1_i2  = cn_ve1u
      ce2_j1  = cn_ve2f  ; ce2_j2  = cn_ve2u
      iioff1  = 1        ; iioff2  = 0
      ijoff1  = 0        ; ijoff2  = 1
      cv_lat  =cn_gphiu
  CASE ( 'V' )
      ! needs fmask, tmask, e1f, e1v, e2t, e2v 
      cmask_i = 'fmask'  ; cmask_j = 'tmask'
      ce1_i1  = cn_ve1f  ; ce1_i2  = cn_ve1v
      ce2_j1  = cn_ve2t  ; ce2_j2  = cn_ve2v
      iioff1  = 0        ; iioff2  = 1
      ijoff1  = 1        ; ijoff2  = 0
      cv_lat  =cn_gphiv
  CASE ( 'F' )
      ! needs vmask, umask, e1v, e1f, e2u, e2f 
      cmask_i = 'vmask'  ; cmask_j = 'umask'
      ce1_i1  = cn_ve1v  ; ce1_i2  = cn_ve1f
      ce2_j1  = cn_ve2u  ; ce2_j2  = cn_ve2f
      iioff1  = 1        ; iioff2  = 0
      ijoff1  = 1        ; ijoff2  = 0
      cv_lat  =cn_gphif
  CASE DEFAULT
      PRINT *, ' TYPE ', TRIM(ct_in),' unknown on C-grid'
      STOP
  END SELECT
      

  ! Allocate the memory
  ALLOCATE ( e1_i1(npiglo,npjglo), e2_j1(npiglo,npjglo) )
  ALLOCATE ( e1_i2(npiglo,npjglo), e2_j2(npiglo,npjglo) )
  ALLOCATE ( rmski(npiglo,npjglo), rmskj(npiglo,npjglo) )
  IF ( l_overf2 ) THEN
     ALLOCATE (  ff(npiglo,npjglo), gphi(npiglo,npjglo) )
  ENDIF

  ALLOCATE ( dlap(npiglo,npjglo), v2d(npiglo, npjglo)  )
  ALLOCATE ( tim(npt) )

  ! Read the metrics from the mesh_hgr file
  e1_i1 = getvar(cn_fhgr, ce1_i1, 1, npiglo, npjglo)
  e1_i2 = getvar(cn_fhgr, ce1_i2, 1, npiglo, npjglo)
  e2_j1 = getvar(cn_fhgr, ce2_j1, 1, npiglo, npjglo)
  e2_j2 = getvar(cn_fhgr, ce2_j2, 1, npiglo, npjglo)
!
  IF ( l_overf2 ) THEN
    ! to prepare computation of g/f*LAP(SSH) 
     gphi    = getvar(cn_fhgr, cv_lat, 1, npiglo, npjglo)
     rpi     = acos(-1.)
     omega   = 2 * rpi / 86400.
     ff     = 2.*omega* sin (gphi*rpi/180.)
     WHERE (ff == 0 ) ff = epsilon(1.0)
  ENDIF

  ! create output fileset
  ncout = create      (cf_out, cf_in, npiglo, npjglo, npk         ) 
  ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout )
  ierr  = putheadervar(ncout,  cf_in, npiglo, npjglo, npk         )

  tim  = getvar1d(cf_in , cn_vtimec, npt     )
  ierr = putvar1d(ncout , tim      , npt, 'T')

  ! Main time loop
  DO jt = 1, npt
     ! Main level loop from top to bottom
     DO jk = 1, npkv
        PRINT *,'jt = ', jt,' jk = ', jk

        ! variable at level jk time jt
        v2d(:,:) =  getvar(cf_in, cv_in, jk, npiglo, npjglo, ktime=jt)
        ! relevant mask
        rmski = getvar(cn_fmsk, cmask_i, jk, npiglo, npjglo)
        rmskj = getvar(cn_fmsk, cmask_j, jk, npiglo, npjglo)
      

        ! Compute laplacian
        dlap(:,:) = 0.d0
        DO jj = 2, npjglo -1
           ij1 = jj + ijoff1
           ij2 = jj - ijoff2
           DO ji = 2, npiglo -1
              ii1 = ji + iioff1
              ii2 = ji - iioff2
              dlap(ji,jj) =  ( ( v2d(ji+1,jj) - v2d(ji,jj) )*1.d0 /e1_i1(ii1,jj)* rmski(ii1,jj)    &
           &  - ( v2d(ji,jj) - v2d(ji-1,jj) )*1.d0/e1_i1(ii2,jj)*  rmski(ii2,jj) ) / e1_i2(ji,jj)  &
           &  + ( ( v2d(ji,jj+1) - v2d(ji,jj) )*1.d0 /e2_j1(ji,ij1)* rmskj(ji,ij1)                 &
           &  - ( v2d(ji,jj) - v2d(ji,jj-1) )*1.d0/e2_j1(ji,ij2)*  rmskj(ji,ij2) ) / e2_j2(ji,jj)
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

END PROGRAM cdflap

