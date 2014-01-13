PROGRAM cdfeke
  !!======================================================================
  !!                     ***  PROGRAM cdfeke   ***
  !!=====================================================================
  !!  ** Purpose : Compute Eddy Kinetic Energy 
  !!
  !!  ** Method  : Use gridU gridU2, gridV gridV2 files produced by
  !!               cdfmoy. Velocities are interpolated both on T points
  !!               and the variance is computed. If -mke option is used
  !!               the program also outputs MKE field
  !!
  !! History : pre  : 11/2004  : J.M. Molines : Original code
  !!           2.1  : 04/2005  : J.M. Molines : use modules
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: ji, jj, jk, jt     ! dummy loop index
  INTEGER(KIND=4)                            :: narg, iargc        ! command line browsing
  INTEGER(KIND=4)                            :: npiglo, npjglo     ! size of the domain (horiz)
  INTEGER(KIND=4)                            :: npk, npt           ! size of the domain vert and time
  INTEGER(KIND=4)                            :: ncout              ! ncid of output file
  INTEGER(KIND=4)                            :: ierr               ! Error status
  INTEGER(KIND=4)                            :: ivar               ! variable counter
  INTEGER(KIND=4)                            :: ip_eke, ip_mke     ! variable index
  INTEGER(KIND=4), DIMENSION(2)              :: ipk, id_varout     ! 

  REAL(KIND=4)                               :: ua, va             ! working arrays
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: tim                ! time variable
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: uc, vc, u2, v2     ! velocities etc...
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: eke                ! velocities etc...
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: rmke               ! Mean Kinetic energy

  CHARACTER(LEN=256)                         :: cf_out='eke.nc'    ! file name
  CHARACTER(LEN=256)                         :: cf_ufil, cf_u2fil  ! file name
  CHARACTER(LEN=256)                         :: cf_vfil, cf_v2fil  !   "
  CHARACTER(LEN=256)                         :: cf_tfil            !   "
  CHARACTER(LEN=256)                         :: cdum               ! dummy character variable

  TYPE(variable), DIMENSION(2)               :: stypvar            !

  LOGICAL                                    :: lchk               ! checking files existence
  LOGICAL                                    :: lperio=.FALSE.     ! checking E-W periodicity
  LOGICAL                                    :: leke=.TRUE.        ! compute EKE
  LOGICAL                                    :: lmke=.FALSE.       ! compute MKE
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfeke U-file [U2-file]  V-file [V2-file] T-file [-mke ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute the Eddy Kinetic Energy from previously computed'
     PRINT *,'        mean values and mean squared values of velocity components.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS : both ''General Use'' or ''Reduced Use'' are acceptable'
     PRINT *,'      * General Use: 5 files are given in argument, and EKE is computed'
     PRINT *,'       U-file  : gridU type file with mean U component.' 
     PRINT *,'       U2-file : gridU2 type file with mean U2 component.' 
     PRINT *,'       V-file  : gridV type file with mean V component.' 
     PRINT *,'       V2-file : gridV2 type file with mean V2 component.' 
     PRINT *,'       T-file  : any gridT or gridT2 (smaller) file, used for EKE header.'
     PRINT *,'       '
     PRINT *,'      * Reduced Use: no U2/V2 file, only MKE is computed from U and V file.'
     PRINT *,'       U-file  : gridU type file with mean U component.' 
     PRINT *,'       V-file  : gridV type file with mean V component.' 
     PRINT *,'       T-file  : any gridT or gridT2 (smaller) file, used for MKE header.'
     PRINT *,'             '
     PRINT *,'     OPTION :'
     PRINT *,'       -mke  : output MKE field together with EKE. If used, must be the last option.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : voeke (m2/s)'
     PRINT *,'         variables : vomke (m2/s) if required'
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  cf_u2fil='none'
  cf_v2fil='none'
  ip_eke=0
  ip_mke=0
  IF ( narg == 3 ) THEN
    ! suppose only mean U and Mean V files as input, and only compute MKE
    CALL getarg (1, cf_ufil )
    CALL getarg (2, cf_vfil )
    CALL getarg (3, cf_tfil )
    lmke = .TRUE.
    leke = .FALSE.
  ELSE 
    CALL getarg (1, cf_ufil )
    CALL getarg (2, cf_u2fil)
    CALL getarg (3, cf_vfil )
    CALL getarg (4, cf_v2fil)
    CALL getarg (5, cf_tfil )
  ENDIF

  IF ( narg == 6 ) THEN
    CALL getarg (narg, cdum )
    IF ( cdum == '-mke'  ) THEN 
      lmke=.TRUE.
    ELSE
      PRINT *, ' Option ', TRIM(cdum),' unknown. We stop.'
      STOP
    ENDIF
  ENDIF

  lchk =           chkfile (cf_ufil )
  lchk = lchk .OR. chkfile (cf_vfil )
  lchk = lchk .OR. chkfile (cf_tfil )
  IF (  leke ) THEN
    lchk = lchk .OR. chkfile (cf_u2fil)
    lchk = lchk .OR. chkfile (cf_v2fil)
  ENDIF

  IF ( lchk ) STOP ! missing files

  npiglo = getdim (cf_ufil,cn_x)
  npjglo = getdim (cf_ufil,cn_y)
  npk    = getdim (cf_ufil,cn_z)
  npt    = getdim (cf_ufil,cn_t)
  
  ivar = 1
  IF ( leke ) THEN 
    ip_eke=ivar
    ipk(ip_eke)                       = npk
    stypvar(ip_eke)%cname             = 'voeke'
    stypvar(ip_eke)%cunits            = 'm2/s2'
    stypvar(ip_eke)%rmissing_value    = 0.
    stypvar(ip_eke)%valid_min         = 0.
    stypvar(ip_eke)%valid_max         = 10000.
    stypvar(ip_eke)%clong_name        = 'Eddy_Kinetic_Energy'
    stypvar(ip_eke)%cshort_name       = 'voeke'
    stypvar(ip_eke)%conline_operation = 'N/A'
    stypvar(ip_eke)%caxis             = 'TZYX'
    ivar = ivar + 1
  ENDIF

  IF ( lmke ) THEN
    ip_mke=ivar
    ipk(ip_mke)                       = npk
    stypvar(ip_mke)%cname             = 'vomke'
    stypvar(ip_mke)%cunits            = 'm2/s2'
    stypvar(ip_mke)%rmissing_value    = 0.
    stypvar(ip_mke)%valid_min         = 0.
    stypvar(ip_mke)%valid_max         = 10000.
    stypvar(ip_mke)%clong_name        = 'Mean_Kinetic_Energy'
    stypvar(ip_mke)%cshort_name       = 'vomke'
    stypvar(ip_mke)%conline_operation = 'N/A'
    stypvar(ip_mke)%caxis             = 'TZYX'
  ENDIF

  ivar = MAX( ip_mke, ip_eke) ! set ivar to the effective number of variables to be output

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE( uc(npiglo,npjglo), vc(npiglo,npjglo) )
  ALLOCATE( tim(npt) )

  IF ( leke ) THEN
    ALLOCATE( eke(npiglo,npjglo)  )
    ALLOCATE( u2(npiglo,npjglo), v2(npiglo,npjglo) )
    eke(:,:) = 0.e0
  ENDIF
  IF (lmke ) THEN
    ALLOCATE( rmke(npiglo,npjglo)  )
    rmke(:,:) = 0.e0
  ENDIF

  ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk       )
  ierr  = createvar   (ncout,  stypvar, ivar,   ipk,    id_varout )
  ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk       )

  ! check for E_W periodicity
  uc(:,:) = getvar(cf_tfil, cn_vlon2d, 1, npiglo, npjglo )
  IF ( uc(1,1) ==  uc(npiglo-1,1) ) THEN 
     lperio = .TRUE. 
     PRINT *,' E-W periodicity detected '
  ENDIF

  DO jt = 1, npt  ! input file is likely to contain only one time frame but who knows ...
    DO jk = 1, npk
      uc(:,:) = getvar(cf_ufil,  cn_vozocrtx,               jk, npiglo, npjglo, ktime=jt )
      vc(:,:) = getvar(cf_vfil,  cn_vomecrty,               jk, npiglo, npjglo, ktime=jt )
      IF ( leke ) THEN
        u2(:,:) = getvar(cf_u2fil, TRIM(cn_vozocrtx)//'_sqd', jk ,npiglo, npjglo, ktime=jt )
        v2(:,:) = getvar(cf_v2fil, TRIM(cn_vomecrty)//'_sqd', jk ,npiglo, npjglo, ktime=jt )
      ENDIF

      ua = 0. ; va = 0. ; eke(:,:) = 0.
      IF ( leke ) THEN
        DO ji=2, npiglo
          DO jj=2,npjglo
              ua = 0.5* ((u2(ji,jj)-uc(ji,jj)*uc(ji,jj))+ (u2(ji-1,jj)-uc(ji-1,jj)*uc(ji-1,jj)))
              va = 0.5* ((v2(ji,jj)-vc(ji,jj)*vc(ji,jj))+ (v2(ji,jj-1)-vc(ji,jj-1)*vc(ji,jj-1)))
              eke(ji,jj) = 0.5 * ( ua + va )
          END DO
        END DO
      ENDIF

      IF ( lmke ) THEN
        DO ji=2, npiglo
          DO jj=2,npjglo
              rmke(ji,jj)=  0.5* (0.5*( uc(ji,jj)*uc(ji,jj) + uc(ji-1,jj)*uc(ji-1,jj)) + &
              &                   0.5*( vc(ji,jj)*vc(ji,jj) + vc(ji,jj-1)*vc(ji,jj-1)) )
          END DO
        END DO
      ENDIF

      IF ( lperio ) eke(1,:) = eke(npiglo-1,:)
      IF ( leke ) THEN 
         IF ( lperio ) eke(1,:) = eke(npiglo-1,:)
         ierr=putvar(ncout,id_varout(ip_eke), eke,  jk ,npiglo, npjglo, ktime=jt )
      ENDIF
      IF ( lmke ) THEN 
         IF ( lperio ) rmke(1,:) = rmke(npiglo-1,:)
         ierr=putvar(ncout,id_varout(ip_mke), rmke, jk ,npiglo, npjglo, ktime=jt )
      ENDIF
    END DO
  END DO ! time loop

  tim  = getvar1d(cf_ufil, cn_vtimec, npt     )
  ierr = putvar1d(ncout,   tim,       npt, 'T')

  ierr = closeout(ncout)

END PROGRAM cdfeke
