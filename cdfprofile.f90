PROGRAM cdfprofile
  !!---------------------------------------------------------------------
  !!               ***  PROGRAM cdfprofile  ***
  !!
  !!  **  Purpose: extract a verticcal profile from a CDFfile
  !!  
  !!  **  Method:  read (i,j) position of point to extract
  !!               read varname
  !!               print profile
  !!
  !!
  !! history :
  !!   Original :  J.M. Molines June 2005
  !!---------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER :: narg, iargc, istatus
  INTEGER :: jk
  INTEGER :: ilook, jlook
  INTEGER :: npiglo, npjglo, npk
  ! added to write in netcdf
  INTEGER :: kx=1, ky=1, kz            ! dims of netcdf output file
  INTEGER :: jj, nboutput=1                ! number of values to write in cdf output
  INTEGER :: ncout, ierr               ! for netcdf output
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipk, id_varout

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: v2d, lon, lat
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: depth, profile
  ! added to write in netcdf
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::  dumlon, dumlat
  REAL(KIND=4), DIMENSION (1)               ::  tim ! time counter
  REAL(KIND=4), DIMENSION (1,1)             ::  dummymean 
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvar  ! structure of output


  CHARACTER(LEN=256) :: cdum, cfile, cvar, cdep
  ! added to write in netcdf
  CHARACTER(LEN=256) :: cfileoutnc='profile.nc'


  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg /= 4  ) THEN
     PRINT *,' Usage : cdfprofile  I J file varname '
     PRINT *,' Output on standard output and netcdf'
     STOP
  ENDIF


  CALL getarg (1, cdum)
  READ(cdum,*) ilook
  CALL getarg (2, cdum)
  READ(cdum,*) jlook
  CALL getarg(3, cfile)
  CALL getarg(4, cvar)

  npiglo= getdim (cfile,'x')
  npjglo= getdim (cfile,'y')
  npk   = getdim (cfile,'depth',cdep)

  ! Allocate arrays
  ALLOCATE( v2d (npiglo,npjglo), depth(npk) ,profile(npk) )
  ALLOCATE ( typvar(nboutput), ipk(nboutput), id_varout(nboutput) )
  ALLOCATE (dumlon(1,1) , dumlat(1,1) ,lon(npiglo,npjglo), lat(npiglo,npjglo))

  lon(:,:)= getvar(cfile, 'nav_lon',  1 ,npiglo,npjglo)
  lat(:,:)= getvar(cfile, 'nav_lat',  1 ,npiglo,npjglo)

  dumlon(:,:)=lon(ilook,jlook)
  dumlat(:,:)=lat(ilook,jlook)

  DO jj=1,nboutput
     ipk(jj)=npk
  ENDDO

  ! define new variables for output 
  typvar(1)%name=TRIM(cvar)
  typvar(1)%units='Sverdrup'
  typvar%missing_value=99999.
  typvar(1)%valid_min= -1000.
  typvar(1)%valid_max= 1000.
  typvar%scale_factor= 1.
  typvar%add_offset= 0.
  typvar%savelog10= 0.
  !typvar(1)%long_name=
  typvar(1)%short_name=TRIM(cvar)
  typvar%online_operation='N/A'
  typvar%axis='TZ'

  depth(:) = getvar1d(cfile,cdep,npk,istatus)
  kz=npk

  ! create output fileset
  ncout =create(cfileoutnc,'none',kx,ky,npk,cdep='depth')
  ierr= createvar(ncout,typvar,nboutput,ipk,id_varout )
  ierr= putheadervar(ncout, cfile,kx, &
       ky,kz,pnavlon=dumlon,pnavlat=dumlat,pdep=depth)
  tim=getvar1d(cfile,'time_counter',1)
  ierr=putvar1d(ncout,tim,1,'T')

  DO jk=1,npk
     v2d (:,:)= getvar(cfile, cvar,  jk ,npiglo,npjglo)
     profile(jk) = v2d(ilook,jlook)
     ! netcdf output 
     dummymean(1,1)=profile(jk)
     ierr = putvar(ncout, id_varout(1), dummymean, jk, kx, ky )

  END DO
  PRINT *, "FILE : ", TRIM(cfile)
  PRINT *, "    ", TRIM(cdep),"         ", TRIM(cvar),"(",ilook,",",jlook,")"
  DO jk=1, npk
     PRINT *, depth(jk), profile(jk)
  END DO

  ierr = closeout(ncout)


END PROGRAM cdfprofile
