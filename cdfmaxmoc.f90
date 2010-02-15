PROGRAM cdfmaxmoc
  !!---------------------------------------------------------------------------------------------------
  !!              ***  PROGRAM cdfmaxmoc  ***
  !! 
  !!   ** Purpose : Compute the maximum of the overturning fonction from a file calculated by cdfmoc
  !!
  !!   ** Method : maxovt 'ovtfile' latmin latmax depmin depmax
  !!               return ovtmaximum and ovt minimum in the defined range.
  !!               Also give location of those extrema
  !!               works for Atlantic and Global MOC
  !!
  !!  * history:
  !!              July 2005 : original : J.M. Molines
  !!              November :  modified and adapted to cdf output R. Dussin.
  !!---------------------------------------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !!
  USE cdfio
  IMPLICIT NONE
  !
  INTEGER :: jj, jk                    ! dummy loop index
  INTEGER :: npjglo, npk               ! size of the overturning
  INTEGER :: narg, iargc               ! line command stuff
  INTEGER :: jmin, jmax, kmin, kmax    ! (latitude, depth) window where to look at extrema
  INTEGER :: jlatmin, jlatmax, kdmin, kdmax ! index of found extrema
  INTEGER :: iminloc(3), imaxloc(3)    ! temporary array to use with minloc/maxloc
  ! added to write in netcdf
  INTEGER :: kx=1, ky=1, kz=1          ! dims of netcdf output file
  INTEGER :: nboutput=6                ! number of values to write in cdf output
  INTEGER :: ncout, ierr               ! for netcdf output
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipk, id_varout
  !
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE ::  zomoc   ! zonal MOC (1,npjglo,jpk)
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::    rlat    ! latitude (1, npjglo)
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE   ::   gdepw    ! depth read in the header
  REAL(KIND=4)                              ::   ovtmax, ovtmin  ! 
  REAL(KIND=4)                              ::   rlatmin, rlatmax, depmin , depmax
  ! added to write in netcdf
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::  dumlon, dumlat
  REAL(KIND=4), DIMENSION (1)               ::  tim ! time counter
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvar  ! structure of output
  !
  CHARACTER(LEN=256) :: cdum, cfile, comment, cbasin, cvar
  ! added to write in netcdf
  CHARACTER(LEN=256) :: cfileoutnc='maxmoc.nc' , cflagcdf
  ! added to write in netcdf
  LOGICAL :: lwrtcdf=.FALSE.

  ! * main program
  narg=iargc()
  IF (narg >= 6 .AND. narg <= 7 ) THEN
     CALL getarg(1,cfile)
     CALL getarg(2,cbasin)
     CALL getarg(3,cdum)
     READ(cdum,*) rlatmin
     CALL getarg(4,cdum)
     READ(cdum,*) rlatmax
     CALL getarg(5,cdum)
     READ(cdum,*) depmin
     CALL getarg(6,cdum)
     READ(cdum,*) depmax
     IF (narg==7) THEN
        CALL getarg(7,cdum)
        READ(cdum,*) cflagcdf
     ENDIF
  ELSE
     PRINT *,' USAGE: cdfmaxmoc ''ovt_file.nc'' cbasin latmin latmax depmin depmax [cdfout]'
     PRINT *,'        cbasin is one of atl glo inp ind or pac '
     PRINT *,' Output on standard output by default'
     PRINT *,' Output on netcdf is available adding cdfout as last argument'
     STOP
  ENDIF

  IF(cflagcdf=='cdfout') THEN
     lwrtcdf=.TRUE.
  ENDIF

  npjglo=getdim(cfile,'y')
  npk=getdim(cfile,'depth')

  ALLOCATE ( zomoc (1,npjglo,npk) ,gdepw(npk), rlat(1,npjglo))
  gdepw(:)  = -getvar1d(cfile,'depthw',npk)
  rlat(:,:) = getvar(cfile,'nav_lat',1,1,npjglo)

  SELECT CASE (cbasin)
  CASE ('atl')
     cvar='zomsfatl'
  CASE ('glo')
     cvar='zomsfglo'
  CASE ('pac')
     cvar='zomsfpac'
  CASE ('inp')
     cvar='zomsfinp'
  CASE ('ind')
     cvar='zomsfind'
  CASE DEFAULT
     STOP 'basin not found'
  END SELECT

  IF(lwrtcdf) THEN

     ALLOCATE ( typvar(nboutput), ipk(nboutput), id_varout(nboutput) )
     ALLOCATE (dumlon(1,1) , dumlat(1,1) )

     dumlon(:,:)=0.
     dumlat(:,:)=0.

     DO jj=1,nboutput
        ipk(jj)=1
     ENDDO

     ! define new variables for output 
     typvar(1)%name='maxmoc'
     typvar(1)%units='Sverdrup'
     typvar%missing_value=99999.
     typvar(1)%valid_min= -1000.
     typvar(1)%valid_max= 1000.
     typvar%scale_factor= 1.
     typvar%add_offset= 0.
     typvar%savelog10= 0.
     typvar(1)%long_name='Maximum_Overturing'
     typvar(1)%short_name='maxmoc'
     typvar%online_operation='N/A'
     typvar%axis='T'

     typvar(2)%name='minmoc'
     typvar(2)%units='Sverdrup'
     typvar(2)%valid_min= -1000.
     typvar(2)%valid_max= 1000.
     typvar(2)%long_name='Minimum_Overtuning'
     typvar(2)%short_name='minmoc'

     typvar(3)%name='latmaxmoc'
     typvar(3)%units='Degrees'
     typvar(3)%valid_min= -90.
     typvar(3)%valid_max= 90.
     typvar(3)%long_name='Latitude_of_Maximum_Overturing'
     typvar(3)%short_name='latmaxmoc'

     typvar(4)%name='latminmoc'
     typvar(4)%units='Degrees'
     typvar(4)%valid_min= -1000.
     typvar(4)%valid_max= 1000.
     typvar(4)%long_name='Latitude_of_Minimum_Overtuning'
     typvar(4)%short_name='latminmoc'

     typvar(5)%name='depthmaxmoc'
     typvar(5)%units='Meters'
     typvar(5)%valid_min= -10000.
     typvar(5)%valid_max= 0.
     typvar(5)%long_name='Depth_of_Maximum_Overturing'
     typvar(5)%short_name='depthmaxmoc'

     typvar(6)%name='depthminmoc'
     typvar(6)%units='Meters'
     typvar(6)%valid_min= -10000.
     typvar(6)%valid_max= 0.
     typvar(6)%long_name='Depth_of_Minimum_Overtuning'
     typvar(6)%short_name='depthminmoc'

  ENDIF

  DO jk=1,npk
     zomoc (:,:,jk) = getvar(cfile,cvar,jk,1,npjglo)
  END DO

  ! look for jmin-jmax :
  DO jj=1, npjglo
     IF ( rlat(1,jj) <= rlatmin )  jmin = jj
     IF ( rlat(1,jj) <= rlatmax )  jmax = jj
  END DO

  ! look for kmin kmax
  DO jk=1,npk
     IF ( gdepw(jk) <= depmin ) kmin = jk
     IF ( gdepw(jk) <= depmax ) kmax = jk
  END DO

  ! look for max/min overturning
  ovtmax = MAXVAL(zomoc(1,jmin:jmax,kmin:kmax))
  ovtmin = MINVAL(zomoc(1,jmin:jmax,kmin:kmax))

  ! find location of min/max
  iminloc =MINLOC(zomoc(:,jmin:jmax,kmin:kmax))
  imaxloc =MAXLOC(zomoc(:,jmin:jmax,kmin:kmax))

  ! results from minloc/maxloc is relative to the sub -array given as arguments
  jlatmin= iminloc(2)+jmin -1 ; jlatmax = imaxloc(2)+jmin -1
  kdmin  = iminloc(3)+kmin -1 ; kdmax   = imaxloc(3)+kmin -1

  ! PRINT * , 'latmin  = ', rlat(1,jmin), 'latmax= ', rlat(1,jmax)
  ! PRINT *,  'Dep min = ', gdepw(kmin), 'Dep max = ',gdepw(kmax)
  PRINT *,' Maximum ', ovtmax ,' Sv latitude ', rlat(1,jlatmax),' depth = ', gdepw(kdmax)
  PRINT *,' Minimum ', ovtmin ,' Sv latitude ', rlat(1,jlatmin),' depth = ', gdepw(kdmin)

  IF(lwrtcdf) THEN

     ! create output fileset
     ncout =create(cfileoutnc,'none',kx,ky,kz,cdep='depthw')
     ierr= createvar(ncout,typvar,nboutput,ipk,id_varout )
     ierr= putheadervar(ncout, cfile,kx, &
          ky,kz,pnavlon=dumlon,pnavlat=dumlat,pdep=gdepw)
     tim=getvar1d(cfile,'time_counter',1)
     ierr=putvar1d(ncout,tim,1,'T')

     ! netcdf output 
     ierr = putvar0d(ncout,id_varout(1), REAL(ovtmax) )
     ierr = putvar0d(ncout,id_varout(2), REAL(ovtmin) )
     ierr = putvar0d(ncout,id_varout(3), REAL(rlat(1,jlatmax)) )
     ierr = putvar0d(ncout,id_varout(4), REAL(rlat(1,jlatmin)) )
     ierr = putvar0d(ncout,id_varout(5), REAL(gdepw(kdmax)) )
     ierr = putvar0d(ncout,id_varout(6), REAL(gdepw(kdmin)) )

     ierr = closeout(ncout)

  ENDIF

END PROGRAM cdfmaxmoc
