PROGRAM cdfbn2
  !!======================================================================
  !!                     ***  PROGRAM  cdfbn2  ***
  !!=====================================================================
  !!  ** Purpose : Compute the Brunt Vaissala frequency
  !!               using same algoritm than NEMO
  !!
  !!  ** Method  : Try to avoid 3 d arrays : work with 2 levels a a time
  !!              The brunt-vaisala frequency is computed using the
  !!              polynomial expression of McDougall (1987):
  !!              N^2 = grav * beta * ( alpha/beta*dk[ t ] - dk[ s ] )/e3w
  !!              N2 is then insterpolated at T levels
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

  INTEGER(KIND=4)                              :: jk, jt                   ! dummy loop index
  INTEGER(KIND=4)                              :: ierr                     ! working integer
  INTEGER(KIND=4)                              :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                              :: npiglo, npjglo, npk, npt ! size of the domain
  INTEGER(KIND=4)                              :: iup = 1, idown = 2, itmp ! for swapping the levels
  INTEGER(KIND=4)                              :: ncout                    ! ncid of output file
  INTEGER(KIND=4), DIMENSION(2)                :: ipk, id_varout           ! level and id of output variables

  REAL(KIND=4)                                 :: zpi                      ! 3.14...
  REAL(KIND=4), DIMENSION (:,:,:), ALLOCATABLE :: ztemp, zsal, zwk         ! Array to read 2 layer of data
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zn2                      ! Brunt Vaissala Frequency (N2)
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zmask, e3w               ! mask and metric
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE :: gdep, tim, e3w1d         ! depth and time

  CHARACTER(LEN=256)                           :: cf_tfil, cldum, cv_dep   ! input file name, ...
  CHARACTER(LEN=256)                           :: cf_out = 'bn2.nc'        ! output file name
  CHARACTER(LEN=256)                           :: cglobal                  ! global attribute
  CHARACTER(LEN=80)                            :: cv_e3w  = 'e3w_ps'       ! e3w variable name (partial step)
  CHARACTER(LEN=80)                            :: cv_bn2  = 'vobn2'        ! cdf variable name for N2

  TYPE(variable), DIMENSION(1)                 :: stypvar                  ! variable attribute

  LOGICAL                                      :: l_w=.false.              ! flag for vertical location of bn2
  LOGICAL                                      :: lchk=.true.              ! check missing files
  LOGICAL                                      :: lfull=.false.            ! full step flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfbn2  T-file [W] [-full]'
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the Brunt-Vaissala frequency (N2) according to' 
     PRINT *,'       temperature and salinity given in the input file.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file : netcdf input gridT file for temperature and salinity.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [ W ] : keep N2 at W points. Default is to interpolate N2' 
     PRINT *,'             at T point on the vertical.'
     PRINT *,'       [ -full ] : indicate a full step configuration instead of'
     PRINT *,'                the default partial steps.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fzgr),' is needed for this program.' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ', TRIM(cv_bn2)
     STOP
  ENDIF

  cglobal = 'Partial step computation'

  ijarg = 1
  CALL getarg (ijarg, cf_tfil) ; ijarg = ijarg + 1

  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE (cldum)
     CASE ('W','w') ; l_w   = .true.
     CASE ('-full') ; lfull = .true. ; cglobal = 'full step computation'
     CASE DEFAULT   ; PRINT *,' Option not understood :', TRIM(cldum) ; STOP
     END SELECT
  END DO

  lchk = chkfile (cn_fzgr )
  lchk = lchk .OR. chkfile (cf_tfil  )
  IF ( lchk  ) STOP  ! missing files 

  npiglo = getdim (cf_tfil, cn_x)
  npjglo = getdim (cf_tfil, cn_y)
  npk    = getdim (cf_tfil, cn_z)
  npt    = getdim (cf_tfil, cn_t)

  ipk(1)                       = npk  !  3D
  stypvar(1)%cname             = cv_bn2
  stypvar(1)%cunits            = 's-1'
  stypvar(1)%rmissing_value    = -1000.
  stypvar(1)%valid_min         = 0.
  stypvar(1)%valid_max         = 50000.
  stypvar(1)%clong_name        = 'Brunt_Vaissala_Frequency'
  stypvar(1)%cshort_name       = cv_bn2
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TZYX'

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  ALLOCATE (ztemp(npiglo,npjglo,2), zsal(npiglo,npjglo,2) )
  ALLOCATE (zwk(npiglo,npjglo,2), zmask(npiglo,npjglo)    )
  ALLOCATE (zn2(npiglo,npjglo), e3w(npiglo,npjglo)        )
  ALLOCATE (gdep(npk), tim(npt)                           )
  IF ( lfull ) ALLOCATE (e3w1d(npk) )

  cv_dep=cn_gdept
  IF (l_w) cv_dep=cn_gdepw

  gdep(:) = getvare3(cn_fzgr, cv_dep, npk)

  ! create output fileset
  ncout = create      (cf_out,   cf_tfil,  npiglo, npjglo, npk)
  ierr  = createvar   (ncout ,   stypvar,  1,      ipk,    id_varout, cdglobal=TRIM(cglobal))
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

     DO jk = npk-1, 2, -1 
        PRINT *,'level ',jk
        zmask(:,:)=1.
        ztemp(:,:,iup)= getvar(cf_tfil, cn_votemper, jk-1, npiglo, npjglo, ktime=jt)
        WHERE(ztemp(:,:,idown) == 0 ) zmask = 0
        zsal(:,:,iup) = getvar(cf_tfil, cn_vosaline, jk-1, npiglo, npjglo, ktime=jt)

        IF ( lfull ) THEN
           e3w(:,:) = e3w1d(jk)
        ELSE
           e3w(:,:) = getvar(cn_fzgr, cv_e3w  , jk, npiglo, npjglo, ldiom=.true.)
        ENDIF

        zwk(:,:,iup) = eosbn2(ztemp, zsal, gdep(jk), e3w, npiglo, npjglo ,iup, idown)* zmask(:,:)

        IF ( .NOT. l_w ) THEN
           ! now put zn2 at T level (k )
           WHERE ( zwk(:,:,idown) == 0 ) 
              zn2(:,:) =  zwk(:,:,iup)
           ELSEWHERE
              zn2(:,:) = 0.5 * ( zwk(:,:,iup) + zwk(:,:,idown) ) * zmask(:,:)
           END WHERE
        ELSE
           zn2(:,:) = zwk(:,:,iup)
        ENDIF

        WHERE ( zn2 == 0 ) zn2 = -1000.
        ierr = putvar(ncout, id_varout(1), zn2, jk, npiglo, npjglo, ktime=jt )
        itmp = idown ; idown = iup ; iup = itmp

     END DO  ! loop to next level
  END DO

  ierr = closeout(ncout)

END PROGRAM cdfbn2
