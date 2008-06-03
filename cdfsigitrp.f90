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
  !!  $Rev: 131 $
  !!  $Date: 2007-12-14 09:21:24 +0100 (Fri, 14 Dec 2007) $
  !!  $Id: cdfsigtrp.f90 131 2007-12-14 08:21:24Z molines $
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
  INTEGER ,DIMENSION(:), ALLOCATABLE :: imina, imaxa, jmina, jmaxa  !: sections limits
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
  REAL(KIND=4)                                 :: rdep                   !: dpeth ref
  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE   :: zwtrp, zwtrpbin, trpbin       !: transport arrays

  CHARACTER(LEN=80), DIMENSION (:), ALLOCATABLE :: csection                     !: section name
  CHARACTER(LEN=80) :: cfilet, cfileu, cfilev, cfilesec='dens_section.dat'      !: files name
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc'          !: coordinates files
  CHARACTER(LEN=80) :: cfilout='trpsig.txt'                                     !: output file
  CHARACTER(LEN=80) :: cdum                                                     !: dummy string

  LOGICAL    :: l_merid                     !: flag is true for meridional working section
  LOGICAL    :: l_print=.FALSE.             !: flag  for printing additional results
  LOGICAL    :: l_bimg=.FALSE.              !: flag  for bimg output

  !!  * Initialisations

  !  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg < 6  ) THEN
     PRINT *,' Usage : cdfsigtrp  gridTfile  gridUfile gridVfile   sigma_min sigma_max nbins zref [options]'
     PRINT *,'            sigma_min, sigma_max : limit for density bining '
     PRINT *,'                           nbins : number of bins to use '
     PRINT *, '           zref : depth of sigmai'
     PRINT *,'     Possible options :'
     PRINT *,'         -print :additional output is send to std output'
     PRINT *,'         -bimg : 2D (x=lat/lon, y=sigma) output on bimg file for hiso, cumul trp, trp'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc must be in the current directory'
     PRINT *,' File  section.dat must also be in the current directory '
     PRINT *,' Output on trpsig.txt and on standard output '
     STOP
  ENDIF

  !! Read arguments
  CALL getarg (1, cfilet)
  CALL getarg (2, cfileu)
  CALL getarg (3, cfilev)
  CALL getarg (4,cdum)           ;        READ(cdum,*) sigma_min
  CALL getarg (5,cdum)           ;        READ(cdum,*) sigma_max
  CALL getarg (6,cdum)           ;        READ(cdum,*) nbins
  CALL getarg (7,cdum)           ;        READ(cdum,*) rdep
  DO jarg=8, narg
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

  ! Allocate and build sigma levels and section array
  ALLOCATE ( sigma_lev (nbins+1) , trpbin(nsection,nbins)  )

  sigma_lev(1)=sigma_min
  dsigma=( sigma_max - sigma_min) / nbins
  DO jclass =2, nbins+1
     sigma_lev(jclass)= sigma_lev(1) + (jclass-1) * dsigma
  END DO

  ! Look for vertical size of the domain
  npk = getdim (cfilet,'depth')
  ALLOCATE ( gdept(npk), gdepw(npk) )

  ! read gdept, gdepw : it is OK even in partial cells, as we never use the bottom gdep
  gdept(:) = getvare3(coordzgr,'gdept', npk)
  gdepw(:) = getvare3(coordzgr,'gdepw', npk)

  !! *  Main loop on sections

  DO jsec=1,nsection
     l_merid=.FALSE.
     imin=imina(jsec) ; imax=imaxa(jsec) ; jmin=jmina(jsec) ; jmax=jmaxa(jsec)
     IF (imin == imax ) THEN        ! meridional section
        l_merid=.TRUE.
        npts=jmax-jmin

     ELSE IF ( jmin == jmax ) THEN  ! zonal section
        npts=imax-imin 

     ELSE
        PRINT *,' Section ',TRIM(csection(jsec)),' is neither zonal nor meridional :('
        PRINT *,' We skip this section .'
        CYCLE
     ENDIF

     ALLOCATE ( zu(npts, npk), zt(npts,npk), zs(npts,npk) ,zsig(npts,0:npk) )
     ALLOCATE ( eu(npts), e3(npts,npk), gdepu(npts, npk), zmask(npts,npk) )
     ALLOCATE ( tmpm(1,npts,2), tmpz(npts,1,2) )
     ALLOCATE ( zwtrp(npts, nbins+1) , hiso(npts,nbins+1), zwtrpbin(npts,nbins) )

     zt = 0. ; zs = 0. ; zu = 0. ; gdepu= 0. ; zmask = 0.  ; zsig=0.d0

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

           ! temperature
           tmpm(:,:,1)=getvar(cfilet,'votemper',jk,1,npts, kimin=imin, kjmin=jmin+1)
           tmpm(:,:,2)=getvar(cfilet,'votemper',jk,1,npts, kimin=imin+1, kjmin=jmin+1)
           zt(:,jk) = 0.5 * ( tmpm(1,:,1) + tmpm(1,:,2) )

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

           ! temperature
           tmpz(:,:,1)=getvar(cfilet,'votemper',jk, npts, 1, kimin=imin+1, kjmin=jmin)
           tmpz(:,:,2)=getvar(cfilet,'votemper',jk, npts, 1, kimin=imin+1, kjmin=jmin+1)
           zt(:,jk) = 0.5 * ( tmpz(:,1,1) + tmpz(:,1,2) )
        END DO

     ENDIF

     ! compute density only for wet points
     zsig(:,1:nk)=sigmai( zt, zs, rdep, npts, nk)*zmask(:,:)
     zsig(:,0)=zsig(:,1)-1.e-4   ! dummy layer for easy interpolation

     ! Some control print 
     IF ( l_print ) THEN
        PRINT *,' T (deg C)' 
        DO jk=1,nk
           PRINT 9000, jk,  (zt(ji,jk),ji=1,npts)
        END DO

        PRINT *,' S (PSU)'
        DO jk=1,nk
           PRINT 9000, jk,  (zs(ji,jk),ji=1,npts)
        END DO

        PRINT *,' SIG (kg/m3 - 1000 )'
        DO jk=1,nk
           PRINT 9000, jk,  (zsig(ji,jk),ji=1,npts)
        END DO

        PRINT *,' VELOCITY (cm/s ) '
        DO jk=1,nk
           PRINT 9000, jk,  (zu(ji,jk)*100,ji=1,npts)
        END DO

        PRINT *,' GDEPU (m) '
        DO jk=1,nk
           PRINT 9001,jk,  (gdepu(ji,jk)*zmask(ji,jk),ji=1,npts)
        END DO

        PRINT *, 'E3 (m)'
        DO jk=1,nk
           PRINT 9001,jk,  (e3(ji,jk)*zmask(ji,jk),ji=1,npts)
        END DO
     END IF

     ! compute depth of isopynals (nbins+1 )
     IF (l_print )  PRINT *,' DEP ISO ( m )'
     DO  jiso =1, nbins+1
        sigma=sigma_lev(jiso)
!!!  REM : I and K loop can be inverted if necessary
        DO ji=1,npts
           hiso(ji,jiso) = gdept(npk)
           DO jk=1,nk 
              IF ( zsig(ji,jk) < sigma ) THEN
              ELSE
                 ! interpolate between jk-1 and jk
                 zalfa=(sigma - zsig(ji,jk-1)) / ( zsig(ji,jk) -zsig(ji,jk-1) )
                 IF (ABS(zalfa) > 1.1 ) THEN   ! case zsig(0) = zsig(1)-1.e-4
                    hiso(ji,jiso)= 0.
                 ELSE
                    hiso(ji,jiso)= gdepu(ji,jk)*zalfa + (1.-zalfa)* gdepu(ji,jk-1)
                 ENDIF
                 EXIT
              ENDIF
           END DO
        END DO
        IF (l_print) PRINT 9002, sigma,(hiso(ji,jiso),ji=1,npts)
     END DO

     ! compute transport between surface and isopycn 
     IF (l_print) PRINT *,' TRP SURF -->  ISO (SV)'
     DO jiso = 1, nbins + 1
        sigma=sigma_lev(jiso)
        DO ji=1,npts
           zwtrp(ji,jiso) = 0.d0
           DO jk=1, nk-1
              IF ( gdepw(jk+1) < hiso(ji,jiso) ) THEN
                 zwtrp(ji,jiso)= zwtrp(ji,jiso) + eu(ji)*e3(ji,jk)*zu(ji,jk)
              ELSE  ! last box ( fraction)
                 zwtrp(ji,jiso)= zwtrp(ji,jiso) + eu(ji)*(hiso(ji,jiso)-gdepw(jk))*zu(ji,jk)
                 EXIT  ! jk loop
              ENDIF
           END DO
        END DO
        IF (l_print) PRINT  9003, sigma,(zwtrp(ji,jiso)/1.e6,ji=1,npts)
     END DO

     ! binned transport : difference between 2 isopycns
     IF (l_print) PRINT *,' TRP bins (SV)'
     DO jbin=1, nbins
        sigma=sigma_lev(jbin)
        DO ji=1, npts
           zwtrpbin(ji,jbin) = zwtrp(ji,jbin+1) -  zwtrp(ji,jbin) 
        END DO
        trpbin(jsec,jbin)=SUM(zwtrpbin(:,jbin) )
        IF (l_print) PRINT  9003, sigma,(zwtrpbin(ji,jbin)/1.e6,ji=1,npts), trpbin(jsec,jbin)/1.e6
     END DO
     PRINT *,' Total transport in all bins :',TRIM(csection(jsec)),' ',SUM(trpbin(jsec,:) )/1.e6
     

     ! output of the code for 1 section
     IF (l_bimg) THEN
       ! (along section, depth ) 2D variables
        cdum=TRIM(csection(jsec))//'_trpdep.bimg'
        OPEN(numbimg,FILE=cdum,FORM='UNFORMATTED')
        cdum=' 3 dimensions in this isopycnal file '
        WRITE(numbimg) cdum
        cdum=' 1: T ;  2: S ; 3: sigma ; 4: Velocity '
        WRITE(numbimg) cdum
        WRITE(cdum,'(a,4i5.4)') ' from '//TRIM(csection(jsec)), imin,imax,jmin,jmax
        WRITE(numbimg) cdum
        cdum=' file '//TRIM(cfilet)
        WRITE(numbimg) cdum
        WRITE(numbimg) npts,nk,1,1,4,0
        WRITE(numbimg) 1.,-float(nk),1.,1., 0.
        WRITE(numbimg) 0.
        WRITE(numbimg) 0.
        ! temperature
        WRITE(numbimg) (( REAL(zt(ji,jk)), ji=1,npts) , jk=nk,1,-1 )
        ! salinity
        WRITE(numbimg) (( REAL(zs(ji,jk)), ji=1,npts) , jk=nk,1,-1 )
        ! sigma
        WRITE(numbimg) (( REAL(zsig(ji,jk)), ji=1,npts) , jk=nk,1,-1 )
        ! Velocity
        WRITE(numbimg) (( REAL(zu(ji,jk)), ji=1,npts) , jk=nk,1,-1 )
        CLOSE(numbimg)

        ! (along section, sigma ) 2D variables
        cdum=TRIM(csection(jsec))//'_trpsig.bimg'
        OPEN(numbimg,FILE=cdum,FORM='UNFORMATTED')
        cdum=' 3 dimensions in this isopycnal file '
        WRITE(numbimg) cdum
        cdum=' 1: hiso ;  2: bin trp ; 3: cumulated  trp '
        WRITE(numbimg) cdum
        WRITE(cdum,'(a,4i5.4)') ' from '//TRIM(csection(jsec)), imin,imax,jmin,jmax
        WRITE(numbimg) cdum
        cdum=' file '//TRIM(cfilet)
        WRITE(numbimg) cdum
        WRITE(numbimg) npts,nbins,1,1,3,0
        WRITE(numbimg) 1.,-REAL(sigma_lev(nbins)),1.,REAL(dsigma), 0.
        WRITE(numbimg) 0.
        WRITE(numbimg) 0.
        ! hiso
        WRITE(numbimg) (( REAL(hiso(ji,jiso)), ji=1,npts) , jiso=nbins,1,-1)
        ! binned transport
        WRITE(numbimg) (( REAL(zwtrpbin(ji,jiso))/1.e6, ji=1,npts) , jiso=nbins,1,-1)
        ! cumulated transport
        WRITE(numbimg) (( REAL(zwtrp(ji,jiso))/1.e6, ji=1,npts) , jiso=nbins,1,-1)
        CLOSE(numbimg)
     ENDIF

     ! free memory for the next section
     DEALLOCATE ( zu,zt, zs ,zsig ,gdepu, hiso, zwtrp, zwtrpbin )
     DEALLOCATE ( eu, e3 ,tmpm, tmpz,zmask )

  END DO   ! next section

  !! Global Output
     OPEN( numout, FILE=cfilout)
       ipos=INDEX(cfilet,'_gridT.nc')
       WRITE(numout,9006)  TRIM(cfilet(1:ipos-1))
       WRITE(numout,9005) ' sigma  ', (csection(jsec),jsec=1,nsection)
     DO jiso=1,nbins
       WRITE(numout,9004) sigma_lev(jiso), (trpbin(jsec,jiso),jsec=1,nsection)
     ENDDO
     CLOSE(numout)

9000 FORMAT(i7,25f8.3)
9001 FORMAT(i7,25f8.0)
9002 FORMAT(f7.3,25f8.0)
9003 FORMAT(f7.3,25f8.3)
9004 FORMAT(f9.4, 20e16.7)
9005 FORMAT('#',a9, 20(2x,a12,2x) )
9006 FORMAT('# ',a)

CONTAINS
  SUBROUTINE section_init(cdfile,cdsection,kimin,kimax,kjmin,kjmax,knumber)
    IMPLICIT NONE
    ! Arguments
    INTEGER, DIMENSION(:),ALLOCATABLE :: kimin,kimax, kjmin,kjmax
    INTEGER, INTENT(OUT) :: knumber
    CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: cdsection
    CHARACTER(LEN=*), INTENT(IN) :: cdfile

    ! Local variables
    INTEGER :: ii, numit=10, jsec
    CHARACTER(LEN=80) :: cline

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
    ALLOCATE( cdsection(knumber) )
    ALLOCATE( kimin(knumber), kimax(knumber), kjmin(knumber), kjmax(knumber) )
    REWIND(numit)
    DO jsec=1,knumber
       READ(numit,'(a)') cdsection(jsec)
       READ(numit,*) kimin(jsec), kimax(jsec), kjmin(jsec), kjmax(jsec)
    END DO

    CLOSE(numit)

  END SUBROUTINE section_init


END PROGRAM cdfsigtrp
