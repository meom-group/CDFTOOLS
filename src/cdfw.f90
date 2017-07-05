PROGRAM cdfw
  !!======================================================================
  !!                     ***  PROGRAM  cdfw  ***
  !!=====================================================================
  !!  ** Purpose : Compute the 3D w for given gridU gridV files 
  !!               and variables
  !!
  !!  ** Method  : Use the equation on continuity: Integrate the 
  !!               horizontal divergence from bottom to the top.
  !!               ( Use the same routines than in the NEMO code )
  !!
  !! History : 2.1  : 06/2005  : J.M. Molines : Original code
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

  INTEGER(KIND=4)                              :: ji, jj, jk, jt     ! dummy loop index
  INTEGER(KIND=4)                              :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                              :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                              :: narg, iargc, ijarg ! browse line
  INTEGER(KIND=4)                              :: ncout              ! ncid of output file
  INTEGER(KIND=4)                              :: ierr               ! error status
  INTEGER(KIND=4)                              :: itop = 1           ! top array index
  INTEGER(KIND=4)                              :: ibot = 2           ! bottom array index
  INTEGER(KIND=4)                              :: itmp               ! working integer for level swap
  INTEGER(KIND=4), DIMENSION(1)                :: ipk, id_varout     ! levels and varid's of output vars

  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE  :: wn                 ! vertical velocity on the top
  !                                                                  ! and bottom of a cell.
  !                                                                  ! wn(top) is computed
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e1t, e2t           ! horizontal T metrics
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e1v, e2u           ! horizontal V and U metrics
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e3u, e3v, e3t      ! vertical metrics (partial steps)
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: glamt, gphit       ! T longitude latitude
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: un, vn             ! horizontal velocity component
  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: hdivn              ! horizontal divergence
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: gdepw              ! depth of W points
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: tim                ! time counter
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: e31d               ! vertical metrics ( full step)

  CHARACTER(LEN=256)                           :: cf_ufil            ! U file name
  CHARACTER(LEN=256)                           :: cf_vfil            ! V file name
  CHARACTER(LEN=256)                           :: cf_out='w.nc'      ! W file name ( output)
  CHARACTER(LEN=256)                           :: cldum              ! dummy string

  TYPE(variable), DIMENSION(1)                 :: stypvar            ! output attributes

  LOGICAL                                      :: lchk               ! missing files flag
  LOGICAL                                      :: lfull=.FALSE.      ! full step flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg < 2 ) THEN
     PRINT *,' usage : cdfw U-file V-file [ U-var V-var ] [ -full]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the vertical velocity from the vertical integration of'
     PRINT *,'       of the horizontal divergence of the velocity.'
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
     PRINT *,'         variables : ', TRIM(cn_vovecrtz),' (m/s)'
     STOP
  ENDIF

  ijarg = 1
  CALL getarg(ijarg, cf_ufil) ; ijarg = ijarg + 1
  CALL getarg(ijarg, cf_vfil) ; ijarg = ijarg + 1

  DO WHILE (ijarg <= narg )
     CALL getarg(ijarg, cldum) 
     SELECT CASE ( cldum )
     CASE ( '-full' ) 
        lfull = .TRUE.
        ijarg = ijarg + 1
     CASE DEFAULT
        CALL getarg(ijarg, cn_vozocrtx) ; ijarg = ijarg + 1
        CALL getarg(ijarg, cn_vomecrty) ; ijarg = ijarg + 1
     END SELECT
  END DO

  lchk = chkfile (cn_fhgr)
  lchk = chkfile (cn_fzgr) .OR. lchk
  lchk = chkfile (cf_ufil) .OR. lchk
  lchk = chkfile (cf_vfil) .OR. lchk
  IF ( lchk ) STOP 99 ! missing files

  npiglo = getdim(cf_ufil,cn_x)
  npjglo = getdim(cf_ufil,cn_y)
  npk    = getdim(cf_ufil,cn_z)
  npt    = getdim(cf_ufil,cn_t)

  ! define new variables for output
  ipk(1)                       = npk 
  stypvar(1)%cname             = TRIM(cn_vovecrtz)
  stypvar(1)%cunits            = 'm/s'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = -1.
  stypvar(1)%valid_max         = 1.
  stypvar(1)%clong_name        = 'Vertical_Velocity'
  stypvar(1)%cshort_name       = TRIM(cn_vovecrtz)
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TZYX'

  ! Allocate the memory
  ALLOCATE ( e1v(npiglo,npjglo), e2u(npiglo,npjglo) )
  ALLOCATE ( e1t(npiglo,npjglo), e2t(npiglo,npjglo) )
  ALLOCATE ( e3u(npiglo,npjglo), e3v(npiglo,npjglo), e3t(npiglo,npjglo) )
  ALLOCATE ( glamt(npiglo,npjglo), gphit(npiglo,npjglo)  )
  ALLOCATE ( un(npiglo,npjglo), vn(npiglo,npjglo), hdivn(npiglo,npjglo) )
  ALLOCATE ( wn(npiglo,npjglo,2) )
  ALLOCATE ( gdepw(npk), tim(npt) )
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
  gdepw(:) = getvare3(cn_fzgr, cn_gdepw, npk)
  IF ( lfull ) e31d(:) = getvare3(cn_fzgr, cn_ve3t, npk)

  ! create output fileset
  ncout = create      (cf_out, cf_ufil, npiglo, npjglo, npk, cdep=cn_vdepthw     )
  ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout                )
  ierr  = putheadervar(ncout,  'dummy', npiglo, npjglo, npk, glamt, gphit, gdepw )

  tim  = getvar1d(cf_ufil, cn_vtimec, npt     )
  ierr = putvar1d(ncout, tim,       npt, 'T')

  DO jt = 1, npt
     wn(:,:,:) = 0.
     ! Main level loop from bottom to top
     DO jk = npk-1, 1, -1
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

        ! Computation from the bottom
        wn(:,:,itop) = wn(:,:,ibot) - e3t(:,:) * hdivn(:,:)

        ! write wn  on file at level jk (This coculd be epensive at it writes from the bottom ...
        ierr = putvar(ncout, id_varout(1), wn(:,:,itop), jk, npiglo, npjglo, ktime=jt)

        ! swap top and bottom index
        itmp=itop ; itop=ibot ; ibot=itmp 

     END DO  ! loop to next level
  END DO ! loop on time

  ierr = closeout(ncout)


END PROGRAM cdfw

