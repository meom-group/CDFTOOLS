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
  !!         : 4.0  : 03/2017  : J.M. Molines  

  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class data_transformation
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jt          ! dummy loop index
  INTEGER(KIND=4)                           :: narg, iargc,ijarg   ! 
  INTEGER(KIND=4)                           :: npiglo, npjglo, npk ! size of the domain
  INTEGER(KIND=4)                           :: npt                 ! number of time frame
  INTEGER(KIND=4)                           :: npiglox2            ! new model size in x
  INTEGER(KIND=4)                           :: ilev, iindex, itmp
  INTEGER(KIND=4)                           :: ncout
  INTEGER(KIND=4)                           :: ierr
  INTEGER(KIND=4), DIMENSION(3)             :: ipk, id_varout

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zlon, zlat, zvar    ! input arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zlonout, zlatout    ! output arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zlonwork, zlatwork  ! working arrays arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zvarout, zvarwork   ! working arrays arrays

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim

  CHARACTER(LEN=256)                        :: cf_in               ! input file name
  CHARACTER(LEN=256)                        :: cf_out='output.nc'  ! output file name
  CHARACTER(LEN=256)                        :: cv_in               ! input variable  name
  CHARACTER(LEN=256)                        :: cldum               !  dummy character variable

  TYPE(variable), DIMENSION(3)              :: stypvar             ! structure for attribute

  LOGICAL                                   :: lnc4 = .FALSE.      ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdf2matlab -f IN-file -v IN-var -k level [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Convert global nemo input file (ORCA configurations) into' 
     PRINT *,'       a file with monotonically increasing longitudes.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file : input model file.' 
     PRINT *,'       -v IN-var  : netcdf variable name to process.'
     PRINT *,'       -k level   : level to process.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        [-o OUT-file] : specify output file name instead of ',TRIM(cf_out)
     PRINT *,'        [-nc4 ]: Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
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
  ijarg=1
  DO WHILE (ijarg <= narg )
     CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'  ) ; CALL getarg (ijarg, cf_in ) ; ijarg=ijarg+1
     CASE ( '-v'  ) ; CALL getarg (ijarg, cv_in ) ; ijarg=ijarg+1
     CASE ( '-k'  ) ; CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) ilev
        ! options
     CASE ( '-o'  ) ; CALL getarg (ijarg, cf_out) ; ijarg=ijarg+1
     CASE ( '-nc4') ; lnc4 = .TRUE.
     CASE DEFAULT   ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO
  IF ( chkfile (cf_in) ) STOP 99  ! missing file

  npiglo = getdim (cf_in,cn_x)
  npjglo = getdim (cf_in,cn_y)
  npk    = getdim (cf_in,cn_z)
  npt    = getdim (cf_in,cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt
  npiglox2 = 2 * npiglo

  ALLOCATE( zvar(npiglo,npjglo), zlon(npiglo,npjglo), zlat(npiglo,npjglo), dtim(npt)  )
  ALLOCATE( zvarout(npiglox2,npjglo), zlonout(npiglox2,npjglo), zlatout(npiglox2,npjglo) )

  CALL CreateOutput
  zlon(:,:) = getvar(cf_in, cn_vlon2d, 1,    npiglo, npjglo)
  zlat(:,:) = getvar(cf_in, cn_vlat2d, 1,    npiglo, npjglo)
  DO jt=1, npt

  zvar(:,:) = getvar(cf_in, cv_in,     ilev, npiglo, npjglo, ktime=jt )

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

  ierr = putvar(ncout,id_varout(1), zlonout, 1, npiglox2, npjglo, ktime=jt)
  ierr = putvar(ncout,id_varout(2), zlatout, 1, npiglox2, npjglo, ktime=jt)
  ierr = putvar(ncout,id_varout(3), zvarout, 1, npiglox2, npjglo, ktime=jt)
  ENDDO

  ierr = closeout(ncout)

  PRINT *, 'Tip : in matlab, do not plot the last line (e.g. maximum northern latitude) '

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
    ipk(:)                       = 1
    stypvar(1)%ichunk            = (/npiglox2,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = 'lon'
    stypvar(1)%cunits            = 'degrees'
    stypvar(1)%valid_min         = -180.
    stypvar(1)%valid_max         = 540.
    stypvar(1)%clong_name        = 'longitude'
    stypvar(1)%cshort_name       = 'lon'
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'YX'

    stypvar(2)%ichunk            = (/npiglox2,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname             = 'lat'
    stypvar(2)%cunits            = 'degrees'
    stypvar(2)%rmissing_value    = 0.
    stypvar(2)%valid_min         = -90.
    stypvar(2)%valid_max         = 90.
    stypvar(2)%clong_name        = 'latitude'
    stypvar(2)%cshort_name       = 'lat'
    stypvar(2)%conline_operation = 'N/A'
    stypvar(2)%caxis             = 'YX'

    stypvar(3)%ichunk            = (/npiglox2,MAX(1,npjglo/30),1,1 /)
    stypvar(3)%cname             = cv_in
    stypvar(3)%cunits            = ''
    stypvar(3)%rmissing_value    = 0.
    stypvar(3)%clong_name        = ''
    stypvar(3)%cshort_name       = cv_in
    stypvar(3)%conline_operation = 'N/A'
    stypvar(3)%caxis             = 'TYX'

    ncout = create    (cf_out,   cf_in,   npiglox2, npjglo, 1         , ld_nc4=lnc4 )
    ierr  = createvar (ncout,    stypvar, 3,        ipk,    id_varout , ld_nc4=lnc4 )

    dtim = getvar1d(cf_in, cn_vtimec, npt     )
    ierr = putvar1d(ncout, dtim,      npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdf2matlab
