PROGRAM cdftemptrp_full
  !!---------------------------------------------------------------------
  !!               ***  PROGRAM cdftemptrp_full  ***
  !!
  !!  **  Purpose: Compute temperature class Mass Transports  across a section
  !!                FULL STEPS version
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
  !!            read 2 rows of T ( i i+1  or j j+1 )
  !!                compute the mean value at velocity point
  !!            compute the depths of isothermal surfaces
  !!            compute the transport from surface to the isotherm
  !!            compute the transport in each class of temperature
  !!            compute the total transport (for information)
  !!
  !! history :
  !!   Original :  F. Castruccio ( Fall 2006)
  !!
  !!---------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: nbins                                  !: number of density classes
  INTEGER   :: ji, jk, jclass, jsec, jiso, jbin, jarg !: dummy loop index
  INTEGER   :: ipos                                   !: working variable
  INTEGER   :: narg, iargc, nxtarg                    !: command line 
  INTEGER   :: npk, nk                                !: vertical size, number of wet layers in the section
  INTEGER   :: numbimg=10                             !: optional bimg logical unit
  INTEGER   :: numout=11                              !: ascii output

  INTEGER                            :: nsection                    !: number of sections (overall)
  INTEGER ,DIMENSION(:), ALLOCATABLE :: imina, imaxa, jmina, jmaxa  !: sections limits
  INTEGER                            :: imin, imax, jmin, jmax      !: working section limits
  INTEGER                            :: npts                        !: working section number of h-points

  REAL(KIND=4), DIMENSION (1,1)                :: ztmp         !: used to read gdepx(1,1,k,1)
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE :: gdept, gdepw !: depth of T and W points 
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE :: e3t          !: depth of T and W points 
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zt           !: temperature from file 
  REAL(KIND=4), DIMENSION (:,:,:), ALLOCATABLE :: tmpm, tmpz   !: temporary arrays

  ! double precision for cumulative variables
  REAL(KIND=8), DIMENSION (:),     ALLOCATABLE ::  eu                 !: either e1v or e2u
  REAL(KIND=8), DIMENSION (:,:),   ALLOCATABLE ::  zu, e3 , zmask     !: velocities e3 and umask
  REAL(KIND=8), DIMENSION (:,:),   ALLOCATABLE ::  ztemp, gdepu       !: temp., depth of vel points
  REAL(KIND=8)                                 :: temp_min, temp_max,dtemp   !: Min and Max for temp. bining
  REAL(KIND=8)                                 :: temp,zalfa                   !: current working temp.
  REAL(KIND=8), DIMENSION (:),   ALLOCATABLE   :: temp_lev                     !: built array with temp. levels
  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE   :: hiso                          !: depth of isotherms

  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE   :: zwtrp, zwtrpbin, trpbin       !: transport arrays

  CHARACTER(LEN=80), DIMENSION (:), ALLOCATABLE :: csection                     !: section name
  CHARACTER(LEN=80) :: cfilet, cfileu, cfilev, cfilesec='temp_section.dat'      !: files name
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc'          !: coordinates files
  CHARACTER(LEN=80) :: cfilout='trptemp.txt'                                    !: output file
  CHARACTER(LEN=80) :: cdum                                                     !: dummy string

  LOGICAL    :: l_merid                     !: flag is true for meridional working section
  LOGICAL    :: l_print=.FALSE.             !: flag  for printing additional results
  LOGICAL    :: l_bimg=.FALSE.              !: flag  for bimg output

  !!  * Initialisations

  !  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg < 6  ) THEN
     PRINT '(255a)',' Usage : cdftemptrp-full  gridTfile  gridUfile gridVfile   temp_max temp_min nbins [options]'
     PRINT '(255a)','            temp_max, temp_min : limit for temperature bining '
     PRINT '(255a)','                           nbins : number of bins to use '
     PRINT '(255a)','     Possible options :'
     PRINT '(255a)','         -print :additional output is send to std output'
     PRINT '(255a)','         -bimg : 2D (x=lat/lon, y=temp) output on bimg file for hiso, cumul trp, trp'
     PRINT '(255a)',' Files mesh_hgr.nc, mesh_zgr.nc must be in the current directory'
     PRINT '(255a)',' File  temp_section.dat must also be in the current directory '
     PRINT '(255a)',' Output on trptemp.txt'
     STOP
  ENDIF

  !! Read arguments
  CALL getarg (1, cfilet)
  CALL getarg (2, cfileu)
  CALL getarg (3, cfilev)
  CALL getarg (4,cdum)           ;        READ(cdum,*) temp_max
  CALL getarg (5,cdum)           ;        READ(cdum,*) temp_min
  CALL getarg (6,cdum)           ;        READ(cdum,*) nbins

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

  ! Allocate and build temp. levels and section array
  ALLOCATE ( temp_lev (nbins+1) , trpbin(nsection,nbins)  )

  temp_lev(1)=temp_max
  dtemp=( temp_max - temp_min) / nbins
  DO jclass =2, nbins+1
     temp_lev(jclass)= temp_lev(1) - (jclass-1) * dtemp
  END DO

  ! Look for vertical size of the domain
  npk = getdim (cfilet,'depth')
  ALLOCATE ( gdept(npk), gdepw(npk), e3t(npk) )
  
  ! read gdept, gdepw 
     gdept(:) = getvare3(coordzgr, 'gdept',npk)
     gdepw(:) = getvare3(coordzgr, 'gdepw',npk)
     e3t(:)   = getvare3(coordzgr, 'e3t',npk)

  !! *  Main loop on sections
  
  write(*,*) 'nsection',nsection
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

     ALLOCATE ( zu(npts, npk), zt(npts,npk), ztemp(npts,0:npk))
     ALLOCATE ( eu(npts), e3(npts,npk), gdepu(npts, npk), zmask(npts,npk) )
     ALLOCATE ( tmpm(1,npts,2), tmpz(npts,1,2) )
     ALLOCATE ( zwtrp(npts, nbins+1) , hiso(npts,nbins+1), zwtrpbin(npts,nbins) )

     zt = 0. ; zu = 0. ; gdepu= 0. ; zmask = 0. ; ztemp=0.d0

     IF (l_merid ) THEN   ! meridional section at i=imin=imax
        tmpm(:,:,1)=getvar(coordhgr, 'e2u', 1,1,npts, kimin=imin, kjmin=jmin+1)
        eu(:)=tmpm(1,:,1)  ! metrics varies only horizontally
        DO jk=1,npk
           ! initiliaze gdepu to gdept()
           gdepu(:,jk) = gdept(jk)

           ! vertical metrics (Full step )
           e3(:,jk)=e3t(jk)

           ! Normal velocity
           tmpm(:,:,1)=getvar(cfileu,'vozocrtx',jk,1,npts, kimin=imin, kjmin=jmin+1)
           zu(:,jk)=tmpm(1,:,1)

           ! temperature
           tmpm(:,:,1)=getvar(cfilet,'votemper',jk,1,npts, kimin=imin, kjmin=jmin+1)
           tmpm(:,:,2)=getvar(cfilet,'votemper',jk,1,npts, kimin=imin+1, kjmin=jmin+1)
           zmask(:,jk)=tmpm(1,:,1)*tmpm(1,:,2)
           WHERE ( zmask(:,jk) /= 0 ) zmask(:,jk)=1
           ! do not take special care for land value, as the corresponding velocity point is masked
           zt(:,jk) = 0.5 * ( tmpm(1,:,1) + tmpm(1,:,2) )

           ! limitation to 'wet' points
           IF ( SUM(zt(:,jk))  == 0 ) THEN
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

           ! vertical metrics (Full step case)
           e3(:,jk)=e3t(jk)

           ! Normal velocity
           tmpz(:,:,1)=getvar(cfilev,'vomecrty',jk,npts,1, kimin=imin+1, kjmin=jmin)
           zu(:,jk)=tmpz(:,1,1)

           ! temperature
           tmpz(:,:,1)=getvar(cfilet,'votemper',jk, npts, 1, kimin=imin+1, kjmin=jmin)
           tmpz(:,:,2)=getvar(cfilet,'votemper',jk, npts, 1, kimin=imin+1, kjmin=jmin+1)
           zmask(:,jk)=tmpz(:,1,1)*tmpz(:,1,2)
           WHERE ( zmask(:,jk) /= 0 ) zmask(:,jk)=1
           ! do not take special care for land value, as the corresponding velocity point is masked
           zt(:,jk) = 0.5 * ( tmpz(:,1,1) + tmpz(:,1,2) )
       
           ! limitation to 'wet' points
           IF ( SUM(zt(:,jk))  == 0 ) THEN
              nk=jk ! first vertical point of the section full on land
              EXIT  ! as soon as all the points are on land
           ENDIF
           
        END DO

     ENDIF

     ! temp. only for wet points
     ztemp(:,1:nk)=zt(:,:)
     ztemp(:,0)=ztemp(:,1)-1.e-4   ! dummy layer for easy interpolation

     ! Some control print 
     IF ( l_print ) THEN
        PRINT *,' T (deg C)' 
        DO jk=1,nk
           PRINT 9000, jk,  (ztemp(ji,jk),ji=1,npts)
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

     ! compute depth of isotherms (nbins+1 )
     IF (l_print )  PRINT *,' DEP ISO ( m )'
     DO  jiso =1, nbins+1
        temp=temp_lev(jiso)
!!!  REM : I and K loop can be inverted if necessary
        DO ji=1,npts
           hiso(ji,jiso) = gdept(npk)
           DO jk=1,nk 
              IF ( ztemp(ji,jk) > temp ) THEN
              ELSE
                 ! interpolate between jk-1 and jk
                 zalfa=(temp - ztemp(ji,jk-1)) / ( ztemp(ji,jk) -ztemp(ji,jk-1) )
                 IF (ABS(zalfa) > 1.1 ) THEN   ! case ztemp(0) = ztemp(1)-1.e-4
                    hiso(ji,jiso)= 0.
                 ELSE
                    hiso(ji,jiso)= gdepu(ji,jk)*zalfa + (1.-zalfa)* gdepu(ji,jk-1)
                 ENDIF
                 EXIT
              ENDIF
           END DO
        END DO
        IF (l_print) PRINT 9002, temp,(hiso(ji,jiso),ji=1,npts)
     END DO

     ! compute transport between surface and isotherm 
     IF (l_print) PRINT *,' TRP SURF -->  ISO (SV)'
     DO jiso = 1, nbins + 1
        temp=temp_lev(jiso)
        DO ji=1,npts
           zwtrp(ji,jiso) = 0.d0
           DO jk=1, nk
              IF ( gdepw(jk+1) < hiso(ji,jiso) ) THEN
                 zwtrp(ji,jiso)= zwtrp(ji,jiso) + eu(ji)*e3(ji,jk)*zu(ji,jk)
              ELSE  ! last box ( fraction)
                 zwtrp(ji,jiso)= zwtrp(ji,jiso) + eu(ji)*(hiso(ji,jiso)-gdepw(jk))*zu(ji,jk)
                 EXIT  ! jk loop
              ENDIF
           END DO
        END DO
        IF (l_print) PRINT  9003, temp,(zwtrp(ji,jiso)/1.e6,ji=1,npts)
     END DO

     ! binned transport : difference between 2 isotherms
     IF (l_print) PRINT *,' TRP bins (SV)'
     DO jbin=1, nbins
        temp=temp_lev(jbin)
        DO ji=1, npts
           zwtrpbin(ji,jbin) = zwtrp(ji,jbin+1) -  zwtrp(ji,jbin) 
        END DO
        trpbin(jsec,jbin)=SUM(zwtrpbin(:,jbin) )
        IF (l_print) PRINT  9003, temp,(zwtrpbin(ji,jbin)/1.e6,ji=1,npts), trpbin(jsec,jbin)/1.e6
     END DO
     PRINT *,' Total transport in all bins :',TRIM(csection(jsec)),' ',SUM(trpbin(jsec,:) )/1.e6
     

     ! output of the code for 1 section
     IF (l_bimg) THEN
       ! (along section, depth ) 2D variables
        cdum=TRIM(csection(jsec))//'_trpdep.bimg'
        OPEN(numbimg,FILE=cdum,FORM='UNFORMATTED')
        cdum=' 3 dimensions in this file '
        WRITE(numbimg) cdum
        cdum=' 1: T ; 2: Velocity '
        WRITE(numbimg) cdum
        WRITE(cdum,'(a,4i5.4)') ' from '//TRIM(csection(jsec)), imin,imax,jmin,jmax
        WRITE(numbimg) cdum
        cdum=' file '//TRIM(cfilet)
        WRITE(numbimg) cdum
        WRITE(numbimg) npts,nk,1,1,2,0
        WRITE(numbimg) 1.,-float(nk),1.,1., 0.
        WRITE(numbimg) 0.
        WRITE(numbimg) 0.
        ! temperature
        WRITE(numbimg) (( REAL(ztemp(ji,jk)), ji=1,npts) , jk=nk,1,-1 )
        ! Velocity
        WRITE(numbimg) (( REAL(zu(ji,jk)), ji=1,npts) , jk=nk,1,-1 )
        CLOSE(numbimg)

        ! (along section, temp ) 2D variables
        cdum=TRIM(csection(jsec))//'_trptemp.bimg'
        OPEN(numbimg,FILE=cdum,FORM='UNFORMATTED')
        cdum=' 3 dimensions in this file '
        WRITE(numbimg) cdum
        cdum=' 1: hiso ;  2: bin trp '
        WRITE(numbimg) cdum
        WRITE(cdum,'(a,4i5.4)') ' from '//TRIM(csection(jsec)), imin,imax,jmin,jmax
        WRITE(numbimg) cdum
        cdum=' file '//TRIM(cfilet)
        WRITE(numbimg) cdum
        WRITE(numbimg) npts,nbins,1,1,2,0
        WRITE(numbimg) 1.,-REAL(temp_lev(nbins)),1.,REAL(dtemp), 0.
        WRITE(numbimg) 0.
        WRITE(numbimg) 0.
        ! hiso
        WRITE(numbimg) (( REAL(hiso(ji,jiso)), ji=1,npts) , jiso=nbins,1,-1)
        ! binned transport
        WRITE(numbimg) (( REAL(zwtrpbin(ji,jiso))/1.e6, ji=1,npts) , jiso=nbins,1,-1)
        CLOSE(numbimg)
     ENDIF

     ! free memory for the next section
     DEALLOCATE ( zu, zt, ztemp, gdepu, hiso, zwtrp, zwtrpbin )
     DEALLOCATE ( eu, e3 ,tmpm, tmpz, zmask )

  END DO   ! next section

  !! Global Output
     OPEN( numout, FILE=cfilout)
       ipos=INDEX(cfilet,'_gridT.nc')
       WRITE(numout,9006)  TRIM(cfilet(1:ipos-1))
       WRITE(numout,9005) ' temp.  ', (csection(jsec),jsec=1,nsection)
     DO jiso=1,nbins
       WRITE(numout,9004) temp_lev(jiso), (trpbin(jsec,jiso),jsec=1,nsection)
     ENDDO
     CLOSE(numout)

9000 FORMAT(i7,40f8.3)
9001 FORMAT(i7,40f8.0)
9002 FORMAT(f7.3,40f8.0)
9003 FORMAT(f7.3,40f8.3)
9004 FORMAT(f9.4, 40e16.7)
9005 FORMAT('#',a9, 40(2x,a12,2x) )
9006 FORMAT('# ',a)

CONTAINS
  SUBROUTINE section_init(cdfile,cdsection,kimin,kimax,kjmin,kjmax,knumber)
    IMPLICIT NONE
    ! Arguments
    INTEGER, DIMENSION(:), ALLOCATABLE :: kimin,kimax, kjmin,kjmax
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


END PROGRAM cdftemptrp_full
