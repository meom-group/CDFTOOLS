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

  INTEGER(KIND=4)                              :: ji, jj, jk, jt     ! dummy loop index
  INTEGER(KIND=4)                              :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                              :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                              :: narg, iargc, ijarg ! browse line
  INTEGER(KIND=4)                              :: ncout              ! ncid of output file
  INTEGER(KIND=4)                              :: ierr               ! error status
  INTEGER(KIND=4)                              :: iioff1, iioff2     ! i-offset for mask
  INTEGER(KIND=4)                              :: ijoff1, ijoff2     ! j-offset for mask
  INTEGER(KIND=4)                              :: ii1, ii2, ij1, ij2 ! working index
  INTEGER(KIND=4), DIMENSION(1)                :: ipk, id_varout     ! levels and varid's of output vars

  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e1_i1, e1_i2       ! along i horizontal metric
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e2_j1, e2_j2       ! along j horizontal metric
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rmski, rmskj       ! relevant mask
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: v2d                ! input variable
  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: rlap               ! output laplacian
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: tim                ! time counter
  REAL(KIND=4)                                 :: rmissing           ! missing value or_FillValue of input variable

  CHARACTER(LEN=256)                           :: cf_in              ! input file name
  CHARACTER(LEN=256)                           :: cv_in              ! input variable name
  CHARACTER(LEN=256)                           :: cv_units           ! units of input variable name
  CHARACTER(LEN=256)                           :: ct_in              ! input variable type [ T U V F ] on C-grid
  CHARACTER(LEN=256)                           :: cf_out='lap.nc'    !output file name
  CHARACTER(LEN=256)                           :: cln_in             ! Long name of input variable
  CHARACTER(LEN=256)                           :: csn_in             ! Short name of input variable
  CHARACTER(LEN=256)                           :: cldum              ! dummy string
  CHARACTER(LEN=10)                            :: ce1_i1, ce1_i2     ! name of relevant horizontal i-metric
  CHARACTER(LEN=10)                            :: ce2_j1, ce2_j2     ! name of relevant horizontal i-metric
  CHARACTER(LEN=10)                            :: cmask_i, cmask_j   ! name of relevant mask variable

  TYPE(variable), DIMENSION(1)                 :: stypvar            ! output attributes

  LOGICAL                                      :: lchk               ! missing files flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg < 2 ) THEN
     PRINT *,' usage : cdflap IN-file IN-var  IN-type '
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
     PRINT *,'       none so far.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fhgr)," ",TRIM(cn_fzgr) ,' and ', TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : lap (unit/m2)'
     STOP
  ENDIF

  ijarg = 1
  CALL getarg(ijarg, cf_in) ; ijarg = ijarg + 1
  CALL getarg(ijarg, cv_in) ; ijarg = ijarg + 1
  CALL getarg(ijarg, ct_in) ; ijarg = ijarg + 1

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

  ierr = getvaratt (cf_in, cv_in, cv_units, rmissing, cln_in, csn_in)

  ! define new variables for output
  ipk(1)                       = npk 
  stypvar(1)%cname             = 'lap'
  stypvar(1)%cunits            = TRIM(cv_units)//'/m2'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = -10.
  stypvar(1)%valid_max         = 10.
  stypvar(1)%clong_name        = 'Laplacian of '//TRIM(cln_in)
  stypvar(1)%cshort_name       = 'lap' 
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TZYX'

  ! fix mesh mask variables according to variable type
  ! HOT ! need to write down stencil to fully understand those settings ...
  SELECT CASE (ct_in )
  CASE ( ' T ' )
      ! needs umask, vmask, e1u, e1t, e2v, e2t 
      cmask_i = 'umask'  ; cmask_j = 'vmask'
      ce1_i1  = cn_ve1u  ; ce1_i2  = cn_ve1t
      ce2_j1  = cn_ve2v  ; ce2_j2  = cn_ve2t
      iioff1  = 0        ; iioff2  = 1
      ijoff1  = 0        ; ijoff2  = 1
  CASE ( ' U ' )
      ! needs tmask, fmask, e1t, e1u, e2f, e2u 
      cmask_i = 'tmask'  ; cmask_j = 'fmask'
      ce1_i1  = cn_ve1t  ; ce1_i2  = cn_ve1u
      ce2_j1  = cn_ve2f  ; ce2_j2  = cn_ve2u
      iioff1  = 1        ; iioff2  = 0
      ijoff1  = 0        ; ijoff2  = 1
  CASE ( ' V ' )
      ! needs fmask, tmask, e1f, e1v, e2t, e2v 
      cmask_i = 'fmask'  ; cmask_j = 'tmask'
      ce1_i1  = cn_ve1f  ; ce1_i2  = cn_ve1v
      ce2_j1  = cn_ve2t  ; ce2_j2  = cn_ve2v
      iioff1  = 0        ; iioff2  = 1
      ijoff1  = 1        ; ijoff2  = 0
  CASE ( ' F ' )
      ! needs vmask, umask, e1v, e1f, e2u, e2f 
      cmask_i = 'vmask'  ; cmask_j = 'umask'
      ce1_i1  = cn_ve1v  ; ce1_i2  = cn_ve1f
      ce2_j1  = cn_ve2u  ; ce2_j2  = cn_ve2f
      iioff1  = 1        ; iioff2  = 0
      ijoff1  = 1        ; ijoff2  = 0
  END SELECT
      

  ! Allocate the memory
  ALLOCATE ( e1_i1(npiglo,npjglo), e2_j1(npiglo,npjglo) )
  ALLOCATE ( e1_i2(npiglo,npjglo), e2_j2(npiglo,npjglo) )
  ALLOCATE ( rmski(npiglo,npjglo), rmskj(npiglo,npjglo) )

  ALLOCATE ( rlap(npiglo,npjglo), v2d(npiglo, npjglo)  )
  ALLOCATE ( tim(npt) )

  ! Read the metrics from the mesh_hgr file
  e1_i1 = getvar(cn_fhgr, ce1_i1, 1, npiglo, npjglo)
  e1_i2 = getvar(cn_fhgr, ce1_i2, 1, npiglo, npjglo)
  e2_j1 = getvar(cn_fhgr, ce2_j1, 1, npiglo, npjglo)
  e2_j2 = getvar(cn_fhgr, ce2_j2, 1, npiglo, npjglo)

  ! create output fileset
  ncout = create      (cf_out, cf_in, npiglo, npjglo, npk         ) 
  ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout )
  ierr  = putheadervar(ncout,  cf_in, npiglo, npjglo, npk         )

  tim  = getvar1d(cf_in , cn_vtimec, npt     )
  ierr = putvar1d(ncout , tim      , npt, 'T')

  ! Main time loop
  DO jt = 1, npt
     ! Main level loop from top to bottom
     DO jk = 1, npk
        PRINT *,'jt = ', jt,' jk = ', jk

        ! variable at level jk time jt
        v2d(:,:) =  getvar(cf_in, cv_in, jk, npiglo, npjglo, ktime=jt)
        ! relevant mask
        rmski = getvar(cn_fmsk, cmask_i, jk, npiglo, npjglo)
        rmskj = getvar(cn_fmsk, cmask_j, jk, npiglo, npjglo)

        ! Compute laplacian
        rlap(:,:) = 0.
        DO jj = 2, npjglo -1
           ij1 = jj + ijoff1
           ij2 = jj - ijoff2
           DO ji = 2, npiglo -1
              ii1 = ji + iioff1
              ii2 = ji - iioff2
              rlap(ji,jj) =  ( ( v2d(ji+1,jj) - v2d(ji,jj) ) /e1_i1(ii1,jj)* rmski(ii1,jj) - ( v2d(ji,jj) - v2d(ji-1,jj) )/e1_i1(ii2,jj)*  rmski(ii2,jj) ) / e1_i2(ji,jj)  &
                &          + ( ( v2d(ji,jj+1) - v2d(ji,jj) ) /e2_j1(ij1,jj)* rmskj(ji,ij1) - ( v2d(ji,jj) - v2d(ji,jj-1) )/e2_j1(ji,ij2)*  rmskj(ji,ij2) ) / e2_j2(ji,jj)
           END DO
        END DO

        ! write level jk 
        ierr = putvar(ncout, id_varout(1), rlap, jk, npiglo, npjglo, ktime=jt)

     END DO  ! loop to next level
  END DO ! loop on time

  ierr = closeout(ncout)

END PROGRAM cdflap
