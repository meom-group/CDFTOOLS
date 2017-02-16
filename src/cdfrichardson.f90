PROGRAM cdfrichardson
  !!======================================================================
  !!                     ***  PROGRAM  cdfrichardson  ***
  !!=====================================================================
  !!  ** Purpose : Compute the Richardson NUmber
  !!               using same algoritm than NEMO
  !!
  !!  ** Method  : Try to avoid 3 d arrays : work with 2 levels at a time
  !!              The Richardson number is computed as 
  !!              Ri = N^2/ dz(U)**2 
  !!              and dz(U)** [ squared vertical velocity derivative] is :
  !!              dz(ub)*dz(ub) + dz(vb)*dz(vb)
  !!
  !! History : 2.0  : 11/2004  : J.M. Molines : Original code
  !!           2.1  : 04/2005  : J.M. Molines : use cdfio
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames   ! for cdf variable names
  USE eos
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                              :: ji, jj, jk, jt           ! dummy loop index
  INTEGER(KIND=4)                              :: ierr                     ! working integer
  INTEGER(KIND=4)                              :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                              :: npiglo, npjglo, npk, npt ! size of the domain
  INTEGER(KIND=4)                              :: iup = 1, idown = 2, itmp ! for swapping the levels
  INTEGER(KIND=4)                              :: ncout                    ! ncid of output file
  INTEGER(KIND=4), DIMENSION(2)                :: ipk, id_varout           ! level and id of output variables

  REAL(KIND=4)                                 :: zpi                      ! 3.14...
  REAL(KIND=4)                                 :: rspval=0.                ! missing_value
  REAL(KIND=4)                                 :: zcoef, zdku, zdkv, zzri  ! working real
  REAL(KIND=4), DIMENSION (:,:,:), ALLOCATABLE :: ztemp, zsal, zwk         ! Array to read 2 layer of data
  REAL(KIND=4), DIMENSION (:,:,:), ALLOCATABLE :: zu, zv                   ! Array to read 2 layer of velocities
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zri                      ! Richardson number
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zmask, e3w               ! mask and metric
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE :: gdep, tim, e3w1d         ! depth and time

  CHARACTER(LEN=2048)                           :: cldum                    ! dummy char variable
  CHARACTER(LEN=2048)                           :: cf_tfil                  ! input T file name
  CHARACTER(LEN=2048)                           :: cf_ufil                  ! input U file name
  CHARACTER(LEN=2048)                           :: cf_vfil                  ! input V file name
  CHARACTER(LEN=2048)                           :: cf_out = 'richardson.nc' ! output file name
  CHARACTER(LEN=2048)                           :: cglobal                  ! global attribute
  CHARACTER(LEN=80)                            :: cv_e3w  = 'e3w_ps'       ! e3w variable name (partial step)
  CHARACTER(LEN=80)                            :: cv_ric  = 'voric'        ! cdf variable name for N2
  CHARACTER(LEN=80)                            :: cv_dep                   ! cdf variable name for depth

  TYPE(variable), DIMENSION(1)                 :: stypvar                  ! variable attribute

  LOGICAL                                      :: l_w=.FALSE.              ! flag for vertical location of ric
  LOGICAL                                      :: lchk=.TRUE.              ! check missing files
  LOGICAL                                      :: lfull=.FALSE.            ! full step flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfrichardson  gridT gridU gridV [ W ] [-full]'
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the Richardson Number (Ri) according to' 
     PRINT *,'       temperature, salinity and velocity components'
     PRINT *,'       given in the input files.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       gridT : input gridT file for temperature and salinity' 
     PRINT *,'       gridU : input gridU file for zonal velocity component'
     PRINT *,'       gridV : input gridV file for meridional velocity component'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [ W ] : keep N2 at W points. Default is to interpolate N2' 
     PRINT *,'             at T point on the vertical'
     PRINT *,'       [ -full ] : indicate a full step configuration instead of'
     PRINT *,'                the default partial steps.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fzgr),' is needed for this program.' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'       variables : ', TRIM(cv_ric)
     STOP
  ENDIF

  cglobal = 'Partial step computation'

  ijarg = 1
  CALL getarg (ijarg, cf_tfil) ; ijarg = ijarg + 1
  CALL getarg (ijarg, cf_ufil) ; ijarg = ijarg + 1
  CALL getarg (ijarg, cf_vfil) ; ijarg = ijarg + 1

  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE (cldum)
     CASE ('W','w') ; l_w   = .TRUE.
     CASE ('-full') ; lfull = .TRUE. ; cglobal = 'full step computation'
     CASE DEFAULT   ; PRINT *,' Option not understood :', TRIM(cldum) ; STOP
     END SELECT
  END DO

  lchk = chkfile (cn_fzgr )
  lchk = lchk .OR. chkfile (cf_tfil  )
  lchk = lchk .OR. chkfile (cf_ufil  )
  lchk = lchk .OR. chkfile (cf_vfil  )
  IF ( lchk  ) STOP  ! missing files 

  npiglo = getdim (cf_tfil, cn_x)
  npjglo = getdim (cf_tfil, cn_y)
  npk    = getdim (cf_tfil, cn_z)
  npt    = getdim (cf_tfil, cn_t)

  ipk(1)                       = npk  !  3D
  stypvar(1)%cname             = cv_ric
  stypvar(1)%cunits            = 'no'
  stypvar(1)%rmissing_value    = rspval
  stypvar(1)%valid_min         = 0.
  stypvar(1)%valid_max         = 50000.
  stypvar(1)%clong_name        = 'Richardson Number'
  stypvar(1)%cshort_name       = cv_ric
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TZYX'

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  ALLOCATE (ztemp(npiglo,npjglo,2), zsal(npiglo,npjglo,2) )
  ALLOCATE (zu(npiglo,npjglo,2),      zv(npiglo,npjglo,2) )
  ALLOCATE (zwk(npiglo,npjglo,2), zmask(npiglo,npjglo)    )
  ALLOCATE (zri(npiglo,npjglo), e3w(npiglo,npjglo)        )
  ALLOCATE (gdep(npk), tim(npt)               )
  zwk(:,:,:) = rspval
  zri(:,:)   = rspval

  IF ( lfull ) ALLOCATE (e3w1d(npk) )

  cv_dep=cn_gdept
  IF (l_w) cv_dep=cn_gdepw

  gdep(:) = getvare3(cn_fzgr, cv_dep, npk) 

  ! create output fileset
  ncout = create      (cf_out,   cf_tfil,  npiglo, npjglo, npk)
  ierr  = createvar   (ncout ,   stypvar, 1,      ipk,    id_varout, cdglobal=TRIM(cglobal))
  ierr  = putheadervar(ncout,    cf_tfil,  npiglo, npjglo, npk, pdep=gdep)

  zpi=ACOS(-1.)

  tim  = getvar1d(cf_tfil, cn_vtimec, npt   )
  ierr = putvar1d(ncout,  tim,       npt,'T')

  IF ( lfull )  e3w1d(:) = getvare3(cn_fzgr, cn_ve3w, npk)

  gdep(:) = getvare3(cn_fzgr, cn_gdepw, npk)

  DO jt=1,npt
     !  2 levels of T and S are required : iup,idown (with respect to W level)
     !  Compute from bottom to top (for vertical integration)
     ztemp(:,:,idown) = getvar(cf_tfil, cn_votemper,  npk-1, npiglo, npjglo, ktime=jt)
     zsal( :,:,idown) = getvar(cf_tfil, cn_vosaline,  npk-1, npiglo, npjglo, ktime=jt)
     zu( :,:,idown)   = getvar(cf_ufil, cn_vozocrtx,  npk-1, npiglo, npjglo, ktime=jt)
     zv( :,:,idown)   = getvar(cf_vfil, cn_vomecrty,  npk-1, npiglo, npjglo, ktime=jt)

     DO jk = npk-1, 2, -1 
        PRINT *,'level ',jk
        zmask(:,:)=1.
        ztemp(:,:,iup)= getvar(cf_tfil, cn_votemper, jk-1, npiglo, npjglo, ktime=jt)
        WHERE(ztemp(:,:,idown) == 0 ) zmask = 0
        zsal(:,:,iup) = getvar(cf_tfil, cn_vosaline, jk-1, npiglo, npjglo, ktime=jt)
        zu( :,:,iup)  = getvar(cf_ufil, cn_vozocrtx, jk-1, npiglo, npjglo, ktime=jt)
        zv( :,:,iup)  = getvar(cf_vfil, cn_vomecrty, jk-1, npiglo, npjglo, ktime=jt)


        IF ( lfull ) THEN
           e3w(:,:) = e3w1d(jk)
        ELSE
           e3w(:,:) = getvar(cn_fzgr, cv_e3w  , jk, npiglo, npjglo, ldiom=.TRUE.)
        ENDIF

        zwk(:,:,iup) = eosbn2(ztemp, zsal, gdep(jk), e3w, npiglo, npjglo ,iup, idown)* zmask(:,:)

        DO jj = 2, npjglo - 1
           DO ji = 2, npiglo - 1
              zcoef = 0.5 / e3w(ji,jj)
              !                                            ! shear of horizontal velocity
              zdku = zcoef * (  zu(ji-1,jj,iup  ) + zu(ji,jj,iup    )   &
                   &           -zu(ji-1,jj,idown) - zu(ji,jj,idown  )  )
              zdkv = zcoef * (  zv(ji,jj-1,iup  ) + zv(ji,jj,iup    )   &
                   &           -zv(ji,jj-1,idown) - zv(ji,jj,idown  )  )
              !                                            ! richardson number (minimum value set to zero)
              zzri = zwk(ji,jj,iup) / ( zdku*zdku + zdkv*zdkv + 1.e-20 )
              zwk(ji,jj,iup) = MAX( zzri, 0.e0 )
           ENDDO
        ENDDO

        IF ( .NOT. l_w ) THEN
           ! now put zri at T level (k )
           WHERE ( zwk(:,:,idown) == 0 ) 
              zri(:,:) =  zwk(:,:,iup)
           ELSEWHERE
              zri(:,:) = 0.5 * ( zwk(:,:,iup) + zwk(:,:,idown) ) * zmask(:,:)
           END WHERE
        ELSE
           zri(:,:) = zwk(:,:,iup)
        ENDIF

        WHERE ( zri < 0  .AND. zri /= rspval )  zri = rspval
        ierr = putvar(ncout, id_varout(1), zri, jk, npiglo, npjglo, ktime=jt )
        itmp = idown ; idown = iup ; iup = itmp

     END DO  ! loop to next level
  END DO

  ierr = closeout(ncout)

END PROGRAM cdfrichardson
