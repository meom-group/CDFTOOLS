PROGRAM cdfdiv
  !!======================================================================
  !!                     ***  PROGRAM  cdfdiv  ***
  !!=====================================================================
  !!  ** Purpose : Compute the divergence for given gridU gridV files 
  !!               and variables
  !!
  !!  ** Method  : Use the same stencil than in NEMO code for computing
  !!               vertical velocities
  !!
  !! History :  3.0  : 10/2011  : P. Mathiot : first version, based on cdfw.f90
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
  INTEGER(KIND=4)                              :: itmp               ! working integer for level swap
  INTEGER(KIND=4), DIMENSION(1)                :: ipk, id_varout     ! levels and varid's of output vars

  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e1t, e2t           ! horizontal T metrics
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e1v, e2u           ! horizontal V and U metrics
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e3v, e3u, e3t      ! vertical metrics
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: glamt, gphit       ! T longitude latitude
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: un, vn             ! horizontal velocity component
  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: hdivn              ! horizontal divergence
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: gdept              ! depth of T points
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: tim                ! time counter
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: e31d               ! vertical metrics (full step)

  CHARACTER(LEN=256)                           :: cf_ufil            ! U file name
  CHARACTER(LEN=256)                           :: cf_vfil            ! V file name
  CHARACTER(LEN=256)                           :: cf_out='div.nc'      ! W file name ( output)
  CHARACTER(LEN=256)                           :: cldum              ! dummy string

  TYPE(variable), DIMENSION(1)                 :: stypvar            ! output attributes

  LOGICAL                                      :: lchk               ! missing files flag
  LOGICAL                                      :: lfull=.FALSE.      ! full step flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg < 2 ) THEN
     PRINT *,' usage : cdfdiv U-file V-file [ U-var V-var ] [ -full ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the divergence of the flow from the U and V velocity components'
     PRINT *,'       Limitation: coded only for C grid (be carefful with forcing field)' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       U-file : netcdf file with the zonal velocity component.' 
     PRINT *,'       V-file : netcdf file with the meridional velocity component.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [ U-var V-var ] : names of the zonal and meridional velocity '
     PRINT *,'                         components. Default are ', TRIM(cn_vozocrtx),' and ', TRIM(cn_vomecrty)
     PRINT *,'       [ -full ] : in case of full step configuration. Default is partial step.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fhgr),' and ',TRIM(cn_fzgr) 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : div (s-1)'
     STOP
  ENDIF

  ijarg = 1
  CALL getarg(ijarg, cf_ufil) ; ijarg = ijarg + 1
  CALL getarg(ijarg, cf_vfil) ; ijarg = ijarg + 1

  DO WHILE (ijarg <= narg )
     CALL getarg(ijarg, cldum) ;
     SELECT CASE ( cldum )
     CASE ( '-full' ) 
        lfull = .TRUE.
        ijarg = ijarg + 1
     CASE DEFAULT
        CALL getarg(ijarg, cn_vozocrtx) ; ijarg = ijarg + 1
        CALL getarg(ijarg, cn_vomecrty) ; ijarg = ijarg + 1
     END SELECT
  END DO

  PRINT *, cn_vozocrtx, cn_vomecrty

  lchk = chkfile (cn_fhgr)
  lchk = chkfile (cn_fzgr) .OR. lchk
  lchk = chkfile (cf_ufil) .OR. lchk
  lchk = chkfile (cf_vfil) .OR. lchk
  IF ( lchk ) STOP ! missing files

  npiglo = getdim(cf_ufil,cn_x)
  npjglo = getdim(cf_ufil,cn_y)
  npk    = getdim(cf_ufil,cn_z)
  npt    = getdim(cf_ufil,cn_t)

  ! define new variables for output
  ipk(1)                       = npk 
  stypvar(1)%cname             = 'div'
  stypvar(1)%cunits            = 's-1'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = -10.
  stypvar(1)%valid_max         = 10.
  stypvar(1)%clong_name        = 'Divergence field'
  stypvar(1)%cshort_name       = 'div' 
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TZYX'

  ! Allocate the memory
  ALLOCATE ( e1v(npiglo,npjglo), e2u(npiglo,npjglo) )
  ALLOCATE ( e1t(npiglo,npjglo), e2t(npiglo,npjglo) )
  ALLOCATE ( e3u(npiglo,npjglo), e3v(npiglo,npjglo), e3t(npiglo,npjglo) )
  ALLOCATE ( glamt(npiglo,npjglo), gphit(npiglo,npjglo)  )
  ALLOCATE ( un(npiglo,npjglo), vn(npiglo,npjglo), hdivn(npiglo,npjglo) )
  ALLOCATE ( gdept(npk), tim(npt) )
  IF ( lfull ) ALLOCATE ( e31d (npk) )

  ! Read the metrics from the mesh_hgr file
  e2u = getvar(cn_fhgr, cn_ve2u, 1, npiglo, npjglo)
  e1v = getvar(cn_fhgr, cn_ve1v, 1, npiglo, npjglo)
  e1t = getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo)
  e2t = getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo)

  ! and the coordinates   from the mesh_hgr file
  glamt = getvar(cn_fhgr, cn_glamt, 1, npiglo, npjglo)
  gphit = getvar(cn_fhgr, cn_gphit, 1, npiglo, npjglo)

  ! Read the depth of the w points (in the file, it is not a vector but a 1x1xnpk array)
  gdept(:) = getvare3(cn_fzgr, cn_gdept, npk)
  IF ( lfull ) e31d(:) = getvare3(cn_fzgr, cn_ve3t, npk)

  ! create output fileset
  ncout = create      (cf_out, cf_ufil, npiglo, npjglo, npk, cdep=cn_vdepthw     )
  ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout                )
  ierr  = putheadervar(ncout,  'dummy', npiglo, npjglo, npk, glamt, gphit, gdept )

  tim  = getvar1d(cf_ufil, cn_vtimec, npt     )
  ierr = putvar1d(ncout  , tim      , npt, 'T')

  ! Main time loop
  DO jt = 1, npt
     ! Main level loop from top to bottom
     DO jk = 1, npk
        PRINT *,'jt = ', jt,' jk = ', jk

        ! velocities at level jk
        un(:,:) =  getvar(cf_ufil, cn_vozocrtx, jk, npiglo, npjglo, ktime=jt)
        vn(:,:) =  getvar(cf_vfil, cn_vomecrty, jk, npiglo, npjglo, ktime=jt)

        IF ( lfull ) THEN
           e3u(:,:) = e31d(jk)
           e3v(:,:) = e31d(jk)
           e3t(:,:) = e31d(jk)
        ELSE
           ! e3 metrics at level jk ( Partial steps)
           e3u(:,:) = getvar(cn_fzgr, 'e3u_ps', jk, npiglo, npjglo, ldiom=.TRUE.) 
           e3v(:,:) = getvar(cn_fzgr, 'e3v_ps', jk, npiglo, npjglo, ldiom=.TRUE.) 
           e3t(:,:) = getvar(cn_fzgr, 'e3t_ps', jk, npiglo, npjglo, ldiom=.TRUE.) 
        ENDIF

        ! Compute divergence :
        DO jj = 2, npjglo -1
           DO ji = 2, npiglo -1
              hdivn(ji,jj) =   &
                &  (  e2u(ji,jj)*e3u(ji,jj) * un(ji,jj) - e2u(ji-1,jj  )*e3u(ji-1,jj  )  * un(ji-1,jj )     &       
                &   + e1v(ji,jj)*e3v(ji,jj) * vn(ji,jj) - e1v(ji  ,jj-1)*e3v(ji  ,jj-1)  * vn(ji  ,jj-1)  ) &
                & / ( e1t(ji,jj)*e2t(ji,jj) * e3t(ji,jj) )
           END DO
        END DO

        ! write level jk 
        ierr = putvar(ncout, id_varout(1), hdivn, jk, npiglo, npjglo, ktime=jt)

     END DO  ! loop to next level
  END DO ! loop on time

  ierr = closeout(ncout)

END PROGRAM cdfdiv

