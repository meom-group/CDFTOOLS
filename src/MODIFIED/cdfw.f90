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
  !!         :  4.0  : 03/2017  : J.M. Molines  
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
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

  INTEGER(KIND=4)                              :: ji, jj, jk, jt     ! dummy loop index
  INTEGER(KIND=4)                              :: it                 ! time index for vvl
  INTEGER(KIND=4)                              :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                              :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                              :: narg, iargc        ! browse line
  INTEGER(KIND=4)                              :: ijarg, ireq        ! browse line
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
  CHARACTER(LEN=256)                           :: cf_tfil            ! Used for VVL e3t
  CHARACTER(LEN=256)                           :: cf_out='w.nc'      ! W file name ( output)
  CHARACTER(LEN=256)                           :: cldum              ! dummy string

  TYPE(variable), DIMENSION(1)                 :: stypvar            ! output attributes

  LOGICAL                                      :: lchk               ! missing files flag
  LOGICAL                                      :: lfull = .FALSE.    ! full step flag
  LOGICAL                                      :: lnc4  = .FALSE.    ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg < 2 ) THEN
     PRINT *,' usage : cdfw U-file V-file [U-var V-var] [-full] [-o OUT-file] [-nc4] ...'
     PRINT *,'           ...[-vvl T-file ]'
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
     PRINT *,'       [ -o OUT-file ] : specify the output file name instead of ', TRIM(cf_out)
     PRINT *,'       [ -nc4 ]     : Use netcdf4 output with chunking and deflation level 1'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
     PRINT *,'       [ -vvl T-file ] : Use time varying vertical metrics (e3t), provided '
     PRINT *,'                in T-file'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fhgr),' and ',TRIM(cn_fzgr) 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out),' unless -o option is used.'
     PRINT *,'         variables : ', TRIM(cn_vovecrtz),' (m/s)'
     STOP
  ENDIF

  ijarg = 1 ; ireq=0
  DO WHILE (ijarg <= narg )
     CALL getarg(ijarg, cldum) ;  ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-full' ) ; lfull  = .TRUE.
     CASE ( '-vvl'  ) ; CALL getarg(ijarg, cf_tfil) ;  ijarg = ijarg + 1
                        lg_vvl = .TRUE. 
     CASE ( '-o'    ) ; CALL getarg(ijarg, cf_out ) ;  ijarg = ijarg + 1
     CASE ( '-nc4'  ) ; lnc4   = .TRUE.
     CASE DEFAULT
        ireq=ireq+1
        SELECT CASE ( ireq )
        CASE ( 1 )  ; cf_ufil     = cldum
        CASE ( 2 )  ; cf_vfil     = cldum
        CASE ( 3 )  ; cn_vozocrtx = cldum
        CASE ( 4 )  ; cn_vomecrty = cldum
        CASE DEFAULT 
           PRINT *, ' ERROR: Too many ''free'' arguments !'
           STOP 1
        END SELECT
     END SELECT
  END DO

  lchk = chkfile (cn_fhgr)
  lchk = chkfile (cn_fzgr) .OR. lchk
  lchk = chkfile (cf_ufil) .OR. lchk
  lchk = chkfile (cf_vfil) .OR. lchk
  IF ( lchk ) STOP ! missing files

  IF ( lg_vvl ) THEN 
     cn_fe3u = cf_ufil
     cn_fe3v = cf_vfil
     cn_fe3t = cf_tfil
  ENDIF

  npiglo = getdim(cf_ufil,cn_x)
  npjglo = getdim(cf_ufil,cn_y)
  npk    = getdim(cf_ufil,cn_z)
  npt    = getdim(cf_ufil,cn_t)

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

  CALL CreateOutput

  DO jt = 1, npt
     IF ( lg_vvl ) THEN ;  it=jt
     ELSE               ;  it=1
     ENDIF
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
           e3u(:,:) = getvar(cn_fe3u, cn_ve3u, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl ) 
           e3v(:,:) = getvar(cn_fe3v, cn_ve3v, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl ) 
           e3t(:,:) = getvar(cn_fe3t, cn_ve3t, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl ) 
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
    ! define new variables for output
    ipk(1)                       = npk 
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = TRIM(cn_vovecrtz)
    stypvar(1)%cunits            = 'm/s'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = -1.
    stypvar(1)%valid_max         = 1.
    stypvar(1)%clong_name        = 'Vertical_Velocity'
    stypvar(1)%cshort_name       = TRIM(cn_vovecrtz)
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TZYX'

    ! create output fileset
    ncout = create      (cf_out, cf_ufil,  npiglo, npjglo, npk, cdep=cn_vdepthw  , ld_nc4=lnc4 )
    ierr  = createvar   (ncout , stypvar,  1,      ipk,    id_varout,              ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  'dummy',  npiglo, npjglo, npk, glamt, gphit, gdepw            )

    tim  = getvar1d(cf_ufil, cn_vtimec, npt     )
    ierr = putvar1d(ncout, tim,       npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfw

