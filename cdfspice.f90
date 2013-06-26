PROGRAM cdfspice
  !!======================================================================
  !!                     ***  PROGRAM  cdfspice  ***
  !!=====================================================================
  !!  ** Purpose : Compute spiciness 3D field from gridT file
  !!               Store the results on a 'similar' cdf file.
  !!
  !!  ** Method  :  spiciness = sum(i=0,5)[sum(j=0,4)[b(i,j)*theta^i*(s-35)^j]]
  !!                   with:  b     -> coefficients
  !!                          theta -> potential temperature
  !!                          s     -> salinity
  !!
  !!  **  Example:
  !!       spice(15,33)=   0.5445863      0.544586321373410  calcul en double
  !!       spice(15,33)=   0.5445864      (calcul en simple precision)
  !!
  !!  ** References : Flament (2002) "A state variable for characterizing 
  !!              water masses and their diffusive stability: spiciness."
  !!              Progress in Oceanography Volume 54, 2002, Pages 493-501.
  !!
  !! History : 2.1  : 03/2010  : C.O. Dufour  : Original code
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

  INTEGER(KIND=4)                           :: ji, jj, jk, jt     ! dummy loop index
  INTEGER(KIND=4)                           :: ierr               ! error status
  INTEGER(KIND=4)                           :: narg, iargc        ! browse command line
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                           :: ncout              ! ncid of output file
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout     ! level and  varid's

  REAL(KIND=4)                              :: zspval             ! missing value
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                ! time counter

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dtemp              ! temperature
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dtempt             ! temperature
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dsal               ! salinity
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dsalt              ! salinity
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dsalref            ! reference salinity
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dspi               ! spiceness
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dmask              ! 2D mask at current level
  REAL(KIND=8), DIMENSION(6,5)              :: dbet               ! coefficients of spiciness formula

  CHARACTER(LEN=256)                        :: cf_tfil            ! input filename
  CHARACTER(LEN=256)                        :: cf_out='spice.nc'  ! output file name

  TYPE (variable), DIMENSION(1)             :: stypvar            ! structure for attributes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfspice T-file '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the spiceness corresponding to temperatures and salinities'
     PRINT *,'       given in the input file.' 
     PRINT *,'      '
     PRINT *,'       spiciness = sum(i=0,5)[sum(j=0,4)[b(i,j)*theta^i*(s-35)^j]]'
     PRINT *,'                 with:  b     -> coefficients'
     PRINT *,'                        theta -> potential temperature'
     PRINT *,'                        s     -> salinity'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file : netcdf file with temperature and salinity (gridT)' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : vospice'
     PRINT *,'      '
     PRINT *,'     REFERENCE :'
     PRINT *,'       Flament (2002) "A state variable for characterizing '
     PRINT *,'             water masses and their diffusive stability: spiciness."'
     PRINT *,'             Progress in Oceanography Volume 54, 2002, Pages 493-501.'
     STOP
  ENDIF
  IF ( narg == 0 ) THEN
     PRINT *,'usage : cdfspice  gridT '
     PRINT *,'    Output on spice.nc, variable vospice'
     STOP
  ENDIF

  CALL getarg (1, cf_tfil)

  IF ( chkfile(cf_tfil) ) STOP ! missing files

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  ipk(:)                       = npk 
  stypvar(1)%cname             = 'vospice'
  stypvar(1)%cunits            = 'kg/m3'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = -300.
  stypvar(1)%valid_max         = 300.
  stypvar(1)%clong_name        = 'spiciness'
  stypvar(1)%cshort_name       = 'vospice'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TZYX'

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE (dtemp(npiglo,npjglo), dsal (npiglo,npjglo) )
  ALLOCATE (dspi( npiglo,npjglo), dmask(npiglo,npjglo) )
  ALLOCATE (dtempt(npiglo,npjglo), dsalt(npiglo,npjglo))
  ALLOCATE (dsalref(npiglo,npjglo))
  ALLOCATE (tim(npt))

  ! create output fileset
  ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk       )
  ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout )
  ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk       )

  tim  = getvar1d(cf_tfil, cn_vtimec, npt     )
  ierr = putvar1d(ncout,   tim,       npt, 'T')

  zspval = getatt(cf_tfil, cn_vosaline, 'missing_value')

  ! Define coefficients to compute spiciness (R*8)
  dbet(1,1) = 0.d0        ; dbet(1,2) = 7.7442d-01  ; dbet(1,3) = -5.85d-03   ; dbet(1,4) = -9.84d-04   ; dbet(1,5) = -2.06d-04
  dbet(2,1) = 5.1655d-02  ; dbet(2,2) = 2.034d-03   ; dbet(2,3) = -2.742d-04  ; dbet(2,4) = -8.5d-06    ; dbet(2,5) = 1.36d-05
  dbet(3,1) = 6.64783d-03 ; dbet(3,2) = -2.4681d-04 ; dbet(3,3) = -1.428d-05  ; dbet(3,4) = 3.337d-05   ; dbet(3,5) = 7.894d-06
  dbet(4,1) = -5.4023d-05 ; dbet(4,2) = 7.326d-06   ; dbet(4,3) = 7.0036d-06  ; dbet(4,4) = -3.0412d-06 ; dbet(4,5) = -1.0853d-06
  dbet(5,1) = 3.949d-07   ; dbet(5,2) = -3.029d-08  ; dbet(5,3) = -3.8209d-07 ; dbet(5,4) = 1.0012d-07  ; dbet(5,5) = 4.7133d-08
  dbet(6,1) = -6.36d-10   ; dbet(6,2) = -1.309d-09  ; dbet(6,3) = 6.048d-09   ; dbet(6,4) = -1.1409d-09 ; dbet(6,5) = -6.676d-10

  ! Compute spiciness
  DO jt=1,npt
     PRINT *,' TIME = ', jt, tim(jt)/86400.,' days'
     DO jk = 1, npk
        dmask(:,:) = 1.

        dtemp(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
        dsal( :,:) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)

        WHERE(dsal == zspval ) dmask = 0

        ! spiciness at time jt, at level jk  
        dspi(:,:)    = 0.d0
        dsalref(:,:) = dsal(:,:) - 35.d0
        dtempt(:,:)  = 1.d0
        DO ji=1,6
           dsalt(:,:) = 1.d0
           DO jj=1,5
              dspi( :,:) = dspi (:,:) + dbet   (ji,jj) * dtempt(:,:) * dsalt(:,:)
              dsalt(:,:) = dsalt(:,:) * dsalref( :,: )
           END DO
           dtempt(:,:) = dtempt(:,:) * dtemp(:,:)     
        END DO

        ierr = putvar(ncout, id_varout(1), REAL(dspi*dmask), jk, npiglo, npjglo, ktime=jt)

     END DO  ! loop to next level
  END DO  ! next time frame

  ierr = closeout(ncout)

END PROGRAM cdfspice
