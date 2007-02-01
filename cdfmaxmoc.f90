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
  !!              November :  modified and adapted to cdf output
  !!---------------------------------------------------------------------------------------------------
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
  !
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE ::  zomoc   ! zonal MOC (1,npjglo,jpk)
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::    rlat    ! latitude (1, npjglo)
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE   ::   gdepw    ! depth read in the header
  REAL(KIND=4)                              ::   ovtmax, ovtmin  ! 
  REAL(KIND=4)                              ::   rlatmin, rlatmax, depmin , depmax
  !
  CHARACTER(LEN=80) :: cdum, cfile, comment, cbasin, cvar

  ! * main program
  narg=iargc()
  IF (narg /= 6 ) THEN
     PRINT *,' USAGE: cdfmaxmoc ''ovt_file,nc'' cbasin latmin latmax depmin depmax'
     PRINT *,'       cbasin is one of atl glo inp ind or pac '
     PRINT *,' Output on standard output'
     STOP
  ENDIF

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
  ovtmax = maxval(zomoc(1,jmin:jmax,kmin:kmax))
  ovtmin = minval(zomoc(1,jmin:jmax,kmin:kmax))

  ! find location of min/max
  iminloc =minloc(zomoc(:,jmin:jmax,kmin:kmax))
  imaxloc =maxloc(zomoc(:,jmin:jmax,kmin:kmax))

  ! results from minloc/maxloc is relative to the sub -array given as arguments
  jlatmin= iminloc(2)+jmin -1 ; jlatmax = imaxloc(2)+jmin -1
  kdmin  = iminloc(3)+kmin -1 ; kdmax   = imaxloc(3)+kmin -1

  ! PRINT * , 'latmin  = ', rlat(1,jmin), 'latmax= ', rlat(1,jmax)
  ! PRINT *,  'Dep min = ', gdepw(kmin), 'Dep max = ',gdepw(kmax)
  PRINT *,' Maximum ', ovtmax ,' Sv latitude ', rlat(1,jlatmax),' depth = ', gdepw(kdmax)
  PRINT *,' Minimum ', ovtmin ,' Sv latitude ', rlat(1,jlatmin),' depth = ', gdepw(kdmin)

END PROGRAM cdfmaxmoc
