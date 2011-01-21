PROGRAM cdf2matlab
  !!-------------------------------------------------------------------
  !!             ***  PROGRAM cdf2matlab  ***
  !!
  !!  **  Purpose: Reshapes ORCA grids to be matlab-friendly
  !!  
  !!  **  Method: Try to avoid 3 d arrays 
  !!
  !! history :
  !!  Original :  R. Dussin (Jan 2011 )
  !!             
  !!-------------------------------------------------------------------
  !!  $Rev: 256 $
  !!  $Date: 2009-07-21 17:49:27 +0200 (mar. 21 juil. 2009) $
  !!  $Id: cdf2matlab.f90 256 2009-07-21 15:49:27Z molines $
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji, jj
  INTEGER   :: narg, iargc                           !: 
  INTEGER   :: npiglo, npjglo, npk, npiglox2         !: size of the domain
  INTEGER   :: zlev , zindex, ztmp
  INTEGER, DIMENSION(3) :: ipk, id_varout
  REAL(KIND=4) , DIMENSION(:,:), ALLOCATABLE :: zlon, zlat, zvar          ! input arrays
  REAL(KIND=4) , DIMENSION(:,:), ALLOCATABLE :: zlonout, zlatout, zvarout ! output arrays
  REAL(KIND=4) , DIMENSION(:,:), ALLOCATABLE :: zlonwork, zlatwork, zvarwork ! working arrays arrays
  REAL(KIND=4) , DIMENSION(1)                :: timean

  CHARACTER(LEN=256) :: cfile, cvarin, cdum, cfileout='output.nc'        !: file name

  TYPE(variable), DIMENSION(3) :: typvar          !: structure for attribute

  INTEGER    :: ncout
  INTEGER    :: istatus, ierr

  !!  Read command line
  narg= iargc()
  IF ( narg /= 3 ) THEN
     PRINT *,' Usage : cdf2matlab file variable level '
     PRINT *,'   Output on output.nc '
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cfile)
  CALL getarg (2, cvarin)
  CALL getarg (3, cdum) ; READ(cdum,*) zlev

  npiglo= getdim (cfile,'x')
  npjglo= getdim (cfile,'y')
  npk   = getdim (cfile,'depth')

  ipk(:) = 1
  typvar(1)%name= 'lon'
  typvar(1)%units='degrees'
  typvar(1)%valid_min= -180.
  typvar(1)%valid_max= 540.
  typvar(1)%long_name='longitude'
  typvar(1)%short_name='lon'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='YX'

  typvar(2)%name= 'lat'
  typvar(2)%units='degrees'
  typvar(2)%missing_value=0.
  typvar(2)%valid_min= -90.
  typvar(2)%valid_max= 90.
  typvar(2)%long_name='latitude'
  typvar(2)%short_name='lat'
  typvar(2)%online_operation='N/A'
  typvar(2)%axis='YX'

  typvar(3)%name= cvarin
  typvar(3)%units=''
  typvar(3)%missing_value=0.
  typvar(3)%long_name=''
  typvar(3)%short_name=cvarin
  typvar(3)%online_operation='N/A'
  typvar(3)%axis='TYX'

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  npiglox2 = 2 * npiglo

  ALLOCATE( zvar(npiglo,npjglo), zlon(npiglo,npjglo), zlat(npiglo,npjglo) )
  ALLOCATE( zvarout(npiglox2,npjglo),zlonout(npiglox2,npjglo),zlatout(npiglox2,npjglo) )

  ncout =create(cfileout, cfile,npiglox2,npjglo,1)
  ierr= createvar(ncout ,typvar,3, ipk,id_varout )

  zlon(:,:) = getvar(cfile,'nav_lon',1, npiglo, npjglo)
  zlat(:,:) = getvar(cfile,'nav_lat',1, npiglo, npjglo)
  zvar(:,:) = getvar(cfile,cvarin,zlev, npiglo, npjglo)

  DO jj=1,npjglo

     zindex = MINLOC( ABS(zlon(:,jj) + 180 ),1 ) ! find the discontinuity in lon array
     ztmp   = npiglo - zindex + 1

     zlonout(1:ztmp,jj) = zlon(zindex:npiglo,jj) ; zlonout(ztmp+1:npiglo,jj) = zlon(1:zindex-1,jj)
     zlonout(npiglo+1:npiglo+ztmp,jj) = zlon(zindex:npiglo,jj) + 360. 
     zlonout(npiglo+ztmp+1:npiglox2,jj) = zlon(1:zindex-1,jj) + 360.

     zlatout(1:ztmp,jj) = zlat(zindex:npiglo,jj) ; zlatout(ztmp+1:npiglo,jj) = zlat(1:zindex-1,jj)
     zlatout(npiglo+1:npiglo+ztmp,jj) = zlat(zindex:npiglo,jj)  
     zlatout(npiglo+ztmp+1:npiglox2,jj) = zlat(1:zindex-1,jj) 

     zvarout(1:ztmp,jj) = zvar(zindex:npiglo,jj) ; zvarout(ztmp+1:npiglo,jj) = zvar(1:zindex-1,jj)
     zvarout(npiglo+1:npiglo+ztmp,jj) = zvar(zindex:npiglo,jj)  
     zvarout(npiglo+ztmp+1:npiglox2,jj) = zvar(1:zindex-1,jj) 

  END DO

  ! Special treatement for ORCA2

  IF ( ( npiglo .EQ. 182 ) .AND. ( npjglo .EQ. 149 ) ) THEN
     PRINT *, 'Assuming that config is ORCA2'

     ALLOCATE( zvarwork(npiglox2,npjglo),zlonwork(npiglox2,npjglo),zlatwork(npiglox2,npjglo) )

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

  ierr=putvar(ncout,id_varout(1), zlonout, 1, npiglox2, npjglo)
  ierr=putvar(ncout,id_varout(2), zlatout, 1, npiglox2, npjglo)
  ierr=putvar(ncout,id_varout(3), zvarout, 1, npiglox2, npjglo)

  timean=getvar1d(cfile,'time_counter',1)
  ierr=putvar1d(ncout,timean,1,'T')
  istatus = closeout(ncout)

END PROGRAM cdf2matlab
