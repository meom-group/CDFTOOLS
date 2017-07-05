PROGRAM cdf2matlab
  !!======================================================================
  !!                     ***  PROGRAM  cdf2matlab  ***
  !!=====================================================================
  !!  ** Purpose : Reshapes ORCA grids to be matlab-friendly
  !!
  !!  ** Method  : transform input file with monotonically increasing
  !!               longitudes.
  !!
  !! History : 2.1  : 01/2011  : R. Dussin    : Original code
  !!           3.0  : 03/2011  : J.M. Molines : Doctor norm + Lic.
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

  INTEGER(KIND=4)                           :: ji, jj              ! dummy loop index
  INTEGER(KIND=4)                           :: narg, iargc         ! 
  INTEGER(KIND=4)                           :: npiglo, npjglo, npk ! size of the domain
  INTEGER(KIND=4)                           :: npiglox2            ! new model size in x
  INTEGER(KIND=4)                           :: ilev, iindex, itmp
  INTEGER(KIND=4)                           :: ncout
  INTEGER(KIND=4)                           :: ierr
  INTEGER(KIND=4), DIMENSION(3)             :: ipk, id_varout

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zlon, zlat, zvar    ! input arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zlonout, zlatout    ! output arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zlonwork, zlatwork  ! working arrays arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zvarout, zvarwork   ! working arrays arrays
  REAL(KIND=4), DIMENSION(1)                :: tim

  CHARACTER(LEN=256)                        :: cf_in               ! input file name
  CHARACTER(LEN=256)                        :: cf_out='output.nc'  ! output file name
  CHARACTER(LEN=256)                        :: cv_in               ! input variable  name
  CHARACTER(LEN=256)                        :: cldum               !  dummy character variable

  TYPE(variable), DIMENSION(3)              :: stypvar             ! structure for attribute
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg /= 3 ) THEN
     PRINT *,' usage : cdf2matlab IN-file IN-var level '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Convert global nemo input file (ORCA configurations) into' 
     PRINT *,'       a file with monotonically increasing longitudes.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       IN-file : input model file.' 
     PRINT *,'       IN-var  : netcdf variable name to process.'
     PRINT *,'       level   : level to process.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : same name than in input file.'
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cf_in)
  CALL getarg (2, cv_in)
  CALL getarg (3, cldum) ; READ(cldum,*) ilev

  IF ( chkfile (cf_in) ) STOP 99 ! missing file

  npiglo = getdim (cf_in,cn_x)
  npjglo = getdim (cf_in,cn_y)
  npk    = getdim (cf_in,cn_z)

  ipk(:)                       = 1
  stypvar(1)%cname             = 'lon'
  stypvar(1)%cunits            = 'degrees'
  stypvar(1)%valid_min         = -180.
  stypvar(1)%valid_max         = 540.
  stypvar(1)%clong_name        = 'longitude'
  stypvar(1)%cshort_name       = 'lon'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'YX'

  stypvar(2)%cname             = 'lat'
  stypvar(2)%cunits            = 'degrees'
  stypvar(2)%rmissing_value    = 0.
  stypvar(2)%valid_min         = -90.
  stypvar(2)%valid_max         = 90.
  stypvar(2)%clong_name        = 'latitude'
  stypvar(2)%cshort_name       = 'lat'
  stypvar(2)%conline_operation = 'N/A'
  stypvar(2)%caxis             = 'YX'

  stypvar(3)%cname             = cv_in
  stypvar(3)%cunits            = ''
  stypvar(3)%rmissing_value    = 0.
  stypvar(3)%clong_name        = ''
  stypvar(3)%cshort_name       = cv_in
  stypvar(3)%conline_operation = 'N/A'
  stypvar(3)%caxis             = 'TYX'

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk

  npiglox2 = 2 * npiglo

  ALLOCATE( zvar(npiglo,npjglo), zlon(npiglo,npjglo), zlat(npiglo,npjglo) )
  ALLOCATE( zvarout(npiglox2,npjglo), zlonout(npiglox2,npjglo), zlatout(npiglox2,npjglo) )

  ncout = create    (cf_out,   cf_in,   npiglox2, npjglo, 1         )
  ierr  = createvar (ncout,    stypvar, 3,        ipk,    id_varout )

  zlon(:,:) = getvar(cf_in, cn_vlon2d, 1,    npiglo, npjglo)
  zlat(:,:) = getvar(cf_in, cn_vlat2d, 1,    npiglo, npjglo)
  zvar(:,:) = getvar(cf_in, cv_in,     ilev, npiglo, npjglo)

  DO jj=1,npjglo
     iindex = MINLOC( ABS(zlon(:,jj) + 180 ),1 ) ! find the discontinuity in lon array
     itmp   = npiglo - iindex + 1

     zlonout(1:itmp,jj) = zlon(iindex:npiglo,jj) ; zlonout(itmp+1:npiglo,jj) = zlon(1:iindex-1,jj)
     zlonout(npiglo+1:npiglo+itmp  ,jj) = zlon(iindex:npiglo,jj) + 360. 
     zlonout(npiglo+itmp+1:npiglox2,jj) = zlon(1:iindex-1,   jj) + 360.

     zlatout(1:itmp,jj) = zlat(iindex:npiglo,jj) ; zlatout(itmp+1:npiglo,jj) = zlat(1:iindex-1,jj)
     zlatout(npiglo+1:npiglo+itmp,  jj) = zlat(iindex:npiglo,jj)  
     zlatout(npiglo+itmp+1:npiglox2,jj) = zlat(1:iindex-1,   jj) 

     zvarout(1:itmp,jj) = zvar(iindex:npiglo,jj) ; zvarout(itmp+1:npiglo,jj) = zvar(1:iindex-1,jj)
     zvarout(npiglo+1:npiglo+itmp,  jj) = zvar(iindex:npiglo,jj)  
     zvarout(npiglo+itmp+1:npiglox2,jj) = zvar(1:iindex-1,   jj) 
  END DO

  ! Special treatement for ORCA2

  IF ( ( npiglo == 182 ) .AND. ( npjglo == 149 ) ) THEN
     PRINT *, 'Assuming that this config is ORCA2 !'

     ALLOCATE( zvarwork(npiglox2,npjglo), zlonwork(npiglox2,npjglo), zlatwork(npiglox2,npjglo) )

     !! init the arryas
     zlonwork(:,:) = zlonout(:,:)
     zlatwork(:,:) = zlatout(:,:)
     zvarwork(:,:) = zvarout(:,:)

     !! swap values to keep lon increasing 
     zlonwork(131,:) = zlonout(130,:) ; zlonwork(npiglo+131,:) = zlonout(npiglo+130,:)
     zlatwork(131,:) = zlatout(130,:) ; zlatwork(npiglo+131,:) = zlatout(npiglo+130,:)
     zvarwork(131,:) = zvarout(130,:) ; zvarwork(npiglo+131,:) = zvarout(npiglo+130,:)

     zlonwork(130,:) = zlonout(131,:) ; zlonwork(npiglo+130,:) = zlonout(npiglo+131,:)
     zlatwork(130,:) = zlatout(131,:) ; zlatwork(npiglo+130,:) = zlatout(npiglo+131,:)
     zvarwork(130,:) = zvarout(131,:) ; zvarwork(npiglo+130,:) = zvarout(npiglo+131,:)

     !! swapping the arrays
     zlonout(:,:) = zlonwork(:,:)
     zlatout(:,:) = zlatwork(:,:)
     zvarout(:,:) = zvarwork(:,:)

  ENDIF
 
  ierr = putvar(ncout,id_varout(1), zlonout, 1, npiglox2, npjglo)
  ierr = putvar(ncout,id_varout(2), zlatout, 1, npiglox2, npjglo)
  ierr = putvar(ncout,id_varout(3), zvarout, 1, npiglox2, npjglo)

  tim  = getvar1d(cf_in, cn_vtimec, 1     )
  ierr = putvar1d(ncout, tim,       1, 'T')
  ierr = closeout(ncout)

  PRINT *, 'Tip : in matlab, do not plot the last line (e.g. maximum northern latitude) '

END PROGRAM cdf2matlab
