PROGRAM cdfeke
  !!======================================================================
  !!                     ***  PROGRAM cdfeke   ***
  !!=====================================================================
  !!  ** Purpose : Compute Eddy Kinetic Energy 
  !!
  !!  ** Method  : Use gridU gridU2, gridV gridV2 files produced by
  !!               cdfmoy. Velocities are interpolated both on T points
  !!               and the variance is computed
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
  INTEGER(KIND=4), DIMENSION(1)              :: ipk, id_varout     ! 

  REAL(KIND=4)                               :: ua, va             ! working arrays
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: tim                ! time variable
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: uc, vc, u2, v2     ! velocities etc...
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: eke                ! velocities etc...

  CHARACTER(LEN=256)                         :: cf_out='eke.nc'    ! file name
  CHARACTER(LEN=256)                         :: cf_ufil, cf_u2fil  ! file name
  CHARACTER(LEN=256)                         :: cf_vfil, cf_v2fil  !
  CHARACTER(LEN=256)                         :: cf_tfil            !

  TYPE(variable), DIMENSION(1)               :: stypvar            !

  LOGICAL                                    :: lchk               ! checking files existence
  LOGICAL                                    :: lperio=.FALSE.     ! checking E-W periodicity
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  !!  Read command line
  narg= iargc()
  IF ( narg /= 5 ) THEN
     PRINT *,' usage : cdfeke U-file U2-file V-file V2-file T2-file'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute the Eddy Kinetic Energy from previously computed'
     PRINT *,'        mean values and mean squared values of velocity components.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       U-file  : gridU type file with mean U component.' 
     PRINT *,'       U2-file : gridU2 type file with mean U2 component.' 
     PRINT *,'       V-file  : gridV type file with mean V component.' 
     PRINT *,'       V2-file : gridV2 type file with mean V2 component.' 
     PRINT *,'       T2-file : any gridT or gridT2 (smaller) file, used for EKE header.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : voeke (m2/s)'
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cf_ufil )
  CALL getarg (2, cf_u2fil)
  CALL getarg (3, cf_vfil )
  CALL getarg (4, cf_v2fil)
  CALL getarg (5, cf_tfil )

  lchk =           chkfile (cf_ufil )
  lchk = lchk .OR. chkfile (cf_u2fil)
  lchk = lchk .OR. chkfile (cf_vfil )
  lchk = lchk .OR. chkfile (cf_v2fil)
  lchk = lchk .OR. chkfile (cf_tfil )
  IF ( lchk ) STOP ! missing files

  npiglo = getdim (cf_ufil,cn_x)
  npjglo = getdim (cf_ufil,cn_y)
  npk    = getdim (cf_ufil,cn_z)
  npt    = getdim (cf_ufil,cn_t)

  ipk(1)                       = npk
  stypvar(1)%cname             = 'voeke'
  stypvar(1)%cunits            = 'm2/s2'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = 0.
  stypvar(1)%valid_max         = 10000.
  stypvar(1)%clong_name        = 'Eddy_Kinetic_Energy'
  stypvar(1)%cshort_name       = 'voeke'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TZYX'

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE( uc(npiglo,npjglo), u2(npiglo,npjglo), vc(npiglo,npjglo), v2(npiglo,npjglo) )
  ALLOCATE( eke(npiglo,npjglo) , tim(npt) )

  ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk       )
  ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout )
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
      u2(:,:) = getvar(cf_u2fil, TRIM(cn_vozocrtx)//'_sqd', jk ,npiglo, npjglo, ktime=jt )
      v2(:,:) = getvar(cf_v2fil, TRIM(cn_vomecrty)//'_sqd', jk ,npiglo, npjglo, ktime=jt )

      ua = 0. ; va = 0. ; eke(:,:) = 0.
      DO ji=2, npiglo
        DO jj=2,npjglo
          ua = 0.5* ((u2(ji,jj)-uc(ji,jj)*uc(ji,jj))+ (u2(ji-1,jj)-uc(ji-1,jj)*uc(ji-1,jj)))
          va = 0.5* ((v2(ji,jj)-vc(ji,jj)*vc(ji,jj))+ (v2(ji,jj-1)-vc(ji,jj-1)*vc(ji,jj-1)))
          eke(ji,jj) = 0.5 * ( ua + va )
        END DO
      END DO
      IF ( lperio ) eke(1,:) = eke(npiglo-1,:)
      ierr=putvar(ncout,id_varout(1), eke, jk ,npiglo, npjglo, ktime=jt )
    END DO
  END DO ! time loop

  tim  = getvar1d(cf_ufil, cn_vtimec, npt     )
  ierr = putvar1d(ncout,   tim,       npt, 'T')

  ierr = closeout(ncout)

END PROGRAM cdfeke
