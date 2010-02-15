PROGRAM cdfsigtrp
  !!---------------------------------------------------------------------
  !!               ***  PROGRAM cdfsigtrp  ***
  !!
  !!  **  Purpose: Compute density class Mass Transports  across a section
  !!               PARTIAL STEPS version
  !!  
  !!  **  Method:
  !!        -The begining and end point of the section are given in term of f-points index.
  !!        -The program works for zonal or meridional sections.
  !!        -The section definitions are given in an ASCII FILE dens_section.dat
  !!            foreach sections, 2 lines : (i) : section name (String, no blank)
  !!                                       (ii) : imin imax jmin jmax for the section
  !!        -Only vertical slices corrsponding to the sections are read in the files.
  !!            read metrics, depth, etc
  !!            read normal velocity (either vozocrtx oy vomecrty )
  !!            read 2 rows of T and S ( i i+1  or j j+1 )
  !!                compute the mean value at velocity point
  !!                compute sigma0 (can be easily modified for sigmai )
  !!            compute the depths of isopyncal surfaces
  !!            compute the transport from surface to the isopycn
  !!            compute the transport in each class of density
  !!            compute the total transport (for information)
  !!
  !! history :
  !!   Original :  J.M. Molines March 2006
  !!---------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio
  USE eos

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: nbins                                  !: number of density classes
  INTEGER   :: ji, jk, jclass, jsec,jiso , jbin,jarg  !: dummy loop index
  INTEGER   :: ipos                                   !: working variable
  INTEGER   :: narg, iargc, nxtarg                    !: command line 
  INTEGER   :: npk, nk                                !: vertical size, number of wet layers in the section
  INTEGER   :: numbimg=10                             !: optional bimg logical unit
  INTEGER   :: numout=11                              !: ascii output

  INTEGER                            :: nsection                    !: number of sections (overall)
  INTEGER ,DIMENSION(100)            :: imina, imaxa, jmina, jmaxa  !: sections limits
  INTEGER                            :: imin, imax, jmin, jmax      !: working section limits
  INTEGER                            :: npts                        !: working section number of h-points

  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE :: gdept, gdepw !: depth of T and W points 
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zs, zt       !: salinity and temperature from file 
  REAL(KIND=4), DIMENSION (:,:,:), ALLOCATABLE :: tmpm, tmpz   !: temporary arrays

  ! double precision for cumulative variables and densities
  REAL(KIND=8), DIMENSION (:),     ALLOCATABLE ::  eu                 !: either e1v or e2u
  REAL(KIND=8), DIMENSION (:,:),   ALLOCATABLE ::  zu, e3 , zmask     !: velocities e3 and umask
  REAL(KIND=8), DIMENSION (:,:),   ALLOCATABLE ::  zsig ,gdepu        !: density, depth of vel points
  REAL(KIND=8)                                 :: sigma_min, sigma_max,dsigma   !: Min and Max for sigma bining
  REAL(KIND=8)                                 :: sigma,zalfa                   !: current working sigma
  REAL(KIND=8), DIMENSION (:),   ALLOCATABLE   :: sigma_lev                     !: built array with sigma levels
  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE   :: hiso                          !: depth of isopycns

  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE   :: zwtrp,zwtrpp, zwtrpm,zwtrpbin, trpbin!: transport arrays
  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE   :: zwtrpsum, zwtrpsump, zwtrpsumm,heightvein 
  INTEGER, DIMENSION (:), ALLOCATABLE   :: zxxx 

  CHARACTER(LEN=256), DIMENSION (100)           :: csection                     !: section name
  CHARACTER(LEN=256) :: cfilet, cfileu, cfilev, cfilesec='dens_section.dat'      !: files name
  CHARACTER(LEN=256) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc'          !: coordinates files
  CHARACTER(LEN=256) :: cfilout='trpsig.txt'                                     !: output file
  CHARACTER(LEN=256) :: cdum                                                     !: dummy string

  LOGICAL    :: l_merid                     !: flag is true for meridional working section
  LOGICAL    :: l_print=.FALSE.             !: flag  for printing additional results
  LOGICAL    :: l_bimg=.FALSE.              !: flag  for bimg output

  !!  * Initialisations

  !  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg < 3  ) THEN
     PRINT *,' Usage : cdfsigtrp  gridTfile  gridUfile gridVfile   sigma_min sigma_max nbins [options]'
     PRINT *,'            sigma_min, sigma_max : limit for density bining '
     PRINT *,'                           nbins : number of bins to use '
     PRINT *,'     Possible options :'
     PRINT *,'         -print :additional output is send to std output'
     PRINT *,'         -bimg : 2D (x=lat/lon, y=sigma) output on bimg file for hiso, cumul trp, trp'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc must be in the current directory'
     PRINT *,' File  section.dat must also be in the current directory '
     PRINT *,' Output on trpsig.txt'
     STOP
  ENDIF

  !! Read arguments
  CALL getarg (1, cfilet)
  CALL getarg (2, cfileu)
  CALL getarg (3, cfilev)

	
  DO jarg=7, narg
     CALL getarg(jarg,cdum)
     SELECT CASE (cdum)
     CASE ('-print' )
        l_print = .TRUE.
     CASE ('-bimg')
        l_bimg = .TRUE.
     CASE DEFAULT
        PRINT *,' Unknown option ', TRIM(cdum),' ... ignored'
     END SELECT
  END DO

	
  ! Initialise sections from file 
  CALL section_init(cfilesec, csection,imina,imaxa,jmina,jmaxa, nsection)

  npk = getdim (cfilet,'depth')

  ALLOCATE ( gdept(npk), gdepw(npk) )
  ALLOCATE ( zwtrpsum  (nsection, 1), zwtrpsump (nsection, 1), zwtrpsumm (nsection, 1) )
  ALLOCATE ( heightvein(nsection, 1)) 
  ! read gdept, gdepw : it is OK even in partial cells, as we never use the bottom gdep
  gdept(:) = getvare3(coordzgr,'gdept',npk)
  gdepw(:) = getvare3(coordzgr,'gdepw',npk)
    print *, 'gdept', gdept

	
  !! *  Main loop on sections
PRINT*,'ok'

  DO jsec=1,nsection
     l_merid=.FALSE.
     imin=imina(jsec) ; imax=imaxa(jsec) ; jmin=jmina(jsec) ; jmax=jmaxa(jsec)
     IF (imin == imax ) THEN        ! meridional section
        l_merid=.TRUE.
        npts=jmax-jmin

     ELSE IF ( jmin == jmax ) THEN  ! zonal section
        npts=imax-imin 
	PRINT *,' Section ',TRIM(csection(jsec)),' is zonal',npts,imin,imax,jmin,jmax
     ELSE
        PRINT *,' Section ',TRIM(csection(jsec)),' is neither zonal nor meridional :('
        PRINT *,' We skip this section .'
        CYCLE
     ENDIF

	

!gh
	nbins=0

     print *,' allocate deb'
     ALLOCATE ( zu(npts, npk), zt(npts,npk), zs(npts,npk) ,zsig(npts,0:npk) )
     ALLOCATE ( eu(npts), e3(npts,npk), gdepu(npts, npk), zmask(npts,npk) )
     ALLOCATE ( tmpm(1,npts,2), tmpz(npts,1,2) )
     ALLOCATE ( zwtrp(npts, nbins+1) , zwtrpp(npts,nbins+1), zwtrpm(npts,nbins+1) )
!     ALLOCATE ( zxxx(npts+npk+1)  )
     ALLOCATE ( zxxx(2)  )
     print *,' allocate fin'

     zt = 0. ; zs = 0. ; zu = 0. ; gdepu= 0. ; zmask = 0.  ; zsig=0.d0
     print *,' allocate raz'

     IF (l_merid ) THEN   ! meridional section at i=imin=imax
        tmpm(:,:,1)=getvar(coordhgr, 'e2u', 1,1,npts, kimin=imin, kjmin=jmin+1)
        eu(:)=tmpm(1,:,1)  ! metrics varies only horizontally
        DO jk=1,npk
           ! initiliaze gdepu to gdept()
           gdepu(:,jk) = gdept(jk)

           ! vertical metrics (PS case)
           tmpm(:,:,1)=getvar(coordzgr,'e3u_ps',jk,1,npts, kimin=imin, kjmin=jmin+1, ldiom=.true.)
           e3(:,jk)=tmpm(1,:,1)
           tmpm(:,:,1)=getvar(coordzgr,'e3w_ps',jk,1,npts, kimin=imin, kjmin=jmin+1, ldiom=.true.)
           tmpm(:,:,2)=getvar(coordzgr,'e3w_ps',jk,1,npts, kimin=imin+1, kjmin=jmin+1, ldiom=.true.)
           IF (jk >= 2 ) THEN
              DO ji=1,npts
                 gdepu(ji,jk)= gdepu(ji,jk-1) + MIN(tmpm(1,ji,1), tmpm(1,ji,2))
              END DO
           ENDIF

           ! Normal velocity
           tmpm(:,:,1)=getvar(cfileu,'vozocrtx',jk,1,npts, kimin=imin, kjmin=jmin+1)
           zu(:,jk)=tmpm(1,:,1)

           ! salinity and deduce umask for the section
           tmpm(:,:,1)=getvar(cfilet,'vosaline',jk,1,npts, kimin=imin  , kjmin=jmin+1)
           tmpm(:,:,2)=getvar(cfilet,'vosaline',jk,1,npts, kimin=imin+1, kjmin=jmin+1)
           zmask(:,jk)=tmpm(1,:,1)*tmpm(1,:,2)
           WHERE ( zmask(:,jk) /= 0 ) zmask(:,jk)=1
           ! do not take special care for land value, as the corresponding velocity point is masked
           zs(:,jk) = 0.5 * ( tmpm(1,:,1) + tmpm(1,:,2) )

           ! limitation to 'wet' points
           IF ( SUM(zs(:,jk))  == 0 ) THEN
              nk=jk ! first vertical point of the section full on land
              EXIT  ! as soon as all the points are on land
           ENDIF

           END DO

     ELSE                   ! zonal section at j=jmin=jmax
        tmpz(:,:,1)=getvar(coordhgr, 'e1v', 1,npts,1,kimin=imin, kjmin=jmin)
        eu=tmpz(:,1,1)
        DO jk=1,npk
           ! initiliaze gdepu to gdept()
           gdepu(:,jk) = gdept(jk)

           ! vertical metrics (PS case)
           tmpz(:,:,1)=getvar(coordzgr,'e3v_ps',jk, npts, 1, kimin=imin+1, kjmin=jmin, ldiom=.true.)
           e3(:,jk)=tmpz(:,1,1)
           tmpz(:,:,1)=getvar(coordzgr,'e3w_ps',jk,npts,1, kimin=imin+1, kjmin=jmin, ldiom=.true.)
           tmpz(:,:,2)=getvar(coordzgr,'e3w_ps',jk,npts,1, kimin=imin+1, kjmin=jmin+1, ldiom=.true.)
           IF (jk >= 2 ) THEN
              DO ji=1,npts
                 gdepu(ji,jk)= gdepu(ji,jk-1) + MIN(tmpz(ji,1,1), tmpz(ji,1,2))
              END DO
           ENDIF

           ! Normal velocity
           tmpz(:,:,1)=getvar(cfilev,'vomecrty',jk,npts,1, kimin=imin+1, kjmin=jmin)
           zu(:,jk)=tmpz(:,1,1)

           ! salinity and mask
           tmpz(:,:,1)=getvar(cfilet,'vosaline',jk, npts, 1, kimin=imin+1, kjmin=jmin)
           tmpz(:,:,2)=getvar(cfilet,'vosaline',jk, npts, 1, kimin=imin+1, kjmin=jmin+1)
           zmask(:,jk)=tmpz(:,1,1)*tmpz(:,1,2)
           WHERE ( zmask(:,jk) /= 0 ) zmask(:,jk)=1
           ! do not take special care for land value, as the corresponding velocity point is masked
           zs(:,jk) = 0.5 * ( tmpz(:,1,1) + tmpz(:,1,2) )

           ! limitation to 'wet' points
           IF ( SUM(zs(:,jk))  == 0 ) THEN
              nk=jk ! first vertical point of the section full on land
              EXIT  ! as soon as all the points are on land
           ENDIF

           END DO

     ENDIF
     print *,' lecture done '

	zxxx=MINLOC(zu)
        print *,'zxxx', zxxx, jsec
	heightvein(jsec,1)=gdept(zxxx(2))

	     ! compute transport between surface and isopycn 
     IF (l_print) PRINT *,' TRP SURF -->  ISO (SV)'
     DO jiso = 1, nbins + 1
        DO ji=1,npts
        print *,'debut de raz', ji, nbins, jiso, npts
           zwtrp (ji,jiso) = 0.d0
           zwtrpp(ji,jiso) = 0.d0
           zwtrpm(ji,jiso) = 0.d0
        print *,'debut de bcle jk'
           DO jk=1, nk
                 zwtrp(ji,jiso)= zwtrp(ji,jiso) + eu(ji)*e3(ji,jk)*zu(ji,jk)
              IF ( zu(ji,jk) >= 0 ) THEN
                 zwtrpp(ji,jiso)= zwtrpp(ji,jiso) + eu(ji)*e3(ji,jk)*zu(ji,jk)
              ELSE  
	         zwtrpm(ji,jiso)= zwtrpm(ji,jiso) + eu(ji)*e3(ji,jk)*zu(ji,jk)
              ENDIF

           END DO
        print *,'fin de bcle jk'
        END DO

	zwtrpsum (jsec,1)= SUM( zwtrp (:,jiso) )/1.e6
	zwtrpsump(jsec,1)= SUM( zwtrpp(:,jiso) )/1.e6
	zwtrpsumm(jsec,1)= SUM( zwtrpm(:,jiso) )/1.e6
     END DO
     ! free memory for the next section
     print *,' dealloc '
     DEALLOCATE ( zu,zt, zs ,zsig ,gdepu)
     DEALLOCATE ( eu, e3 ,tmpm, tmpz,zmask )
     DEALLOCATE ( zwtrp ,zwtrpp,zwtrpm ,zxxx)
     print *,' dealloc done '

  END DO   ! next section

	PRINT*,'ok end'
  !! Global Output
PRINT*,numout,cfilout
     WRITE(1,9006)  TRIM(cfilet(1:ipos-1))
	PRINT*,'ok'
       WRITE(1,9008) ' Section : ', ' transport     ',' transport v<0 :',' transport v<0 :',' vein height   '
	PRINT*,'ok'
     DO jsec=1,nsection
       WRITE(1,9007) csection(jsec),zwtrpsum  (jsec,1),zwtrpsump (jsec,1),zwtrpsumm (jsec,1),heightvein(jsec,1)
     ENDDO
    
 	OPEN( numout, FILE=cfilout )
	PRINT*,'ok'
       ipos=INDEX(cfilet,'_gridT.nc')
	PRINT*,'ok'
       WRITE(numout,9006)  TRIM(cfilet(1:ipos-1))
	PRINT*,'ok'
       WRITE(numout,9008) ' Section : ', ' transport     ',' transport v<0 :',' transport v<0 :',' vein height   '
	PRINT*,'ok'
     DO jsec=1,nsection
       WRITE(numout,9007) csection(jsec),zwtrpsum  (jsec,1),zwtrpsump (jsec,1),zwtrpsumm (jsec,1),heightvein(jsec,1)
      ! WRITE(numout,9004)' transport     ',  (zwtrpsum  (jsec,1),jsec=1,nsection)
      ! WRITE(numout,9004)' transport v<0 ',  (zwtrpsump (jsec,1),jsec=1,nsection)
      ! WRITE(numout,9004)' transport v>0 ',  (zwtrpsumm (jsec,1),jsec=1,nsection)
      ! WRITE(numout,9004)' vein height   ',  (heightvein(jsec,1),jsec=1,nsection)
     ENDDO
     CLOSE( numout )

	PRINT*,'ok'

9000 FORMAT(i7,25f8.3)
9001 FORMAT(i7,25f8.0)
9002 FORMAT(f7.3,25f8.0)
9003 FORMAT(f7.3,25f8.3)
!9004 FORMAT(a15,e16.7)
9004 FORMAT('#',a15, 20(2x,e16.7,2x))
9005 FORMAT('#',a15, 20(2x,a12,2x) )
9006 FORMAT('# ',a)
9007 FORMAT('#',a8, 4(2x,e16.7,2x))
9008 FORMAT('#',a8, 4(2x,a15,2x) )

CONTAINS
  SUBROUTINE section_init(cdfile,cdsection,kimin,kimax,kjmin,kjmax,knumber)
    IMPLICIT NONE
    ! Arguments
    INTEGER, DIMENSION(100) :: kimin,kimax, kjmin,kjmax
    INTEGER, INTENT(OUT) :: knumber
    CHARACTER(LEN=256), DIMENSION(100) :: cdsection
    CHARACTER(LEN=*), INTENT(IN) :: cdfile

    ! Local variables
    INTEGER :: ii, numit=10, jsec
    CHARACTER(LEN=256) :: cline

    OPEN(numit, FILE=cdfile)
    ii=0

    DO
       READ(numit,'(a)') cline
       IF (INDEX(cline,'EOF') == 0 ) THEN
          READ(numit,*)    ! skip one line
          ii = ii + 1
       ELSE
          EXIT
       ENDIF
    END DO

    knumber=ii
    IF ( knumber > 100 ) THEN
      PRINT *,' ERROR : no more than 100 sections are allowed'
      STOP
    ENDIF
    REWIND(numit)
    DO jsec=1,knumber
       READ(numit,'(a)') cdsection(jsec)
       READ(numit,*) kimin(jsec), kimax(jsec), kjmin(jsec), kjmax(jsec)
    END DO

    CLOSE(numit)

  END SUBROUTINE section_init


END PROGRAM cdfsigtrp
