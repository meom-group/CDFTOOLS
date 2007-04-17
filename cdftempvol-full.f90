PROGRAM cdftempvol_full
  !!---------------------------------------------------------------------
  !!               ***  PROGRAM cdftempvol_full  ***
  !!
  !!  **  Purpose: Compute water volume in a given domain between isotherms 
  !!                FULL STEPS version
  !! 
  !!  **  Method   :  compute the sum ( e1 * e2 * e3 * mask )
  !!        -The box boundary are given by imin, imax, jmin, jmax
  !!            read metrics, depth, etc
  !!            read T and SSH
  !!            compute the depths of isothermal surfaces
  !!            compute the volume from surface to the isotherm
  !!            compute the volume in each class of temperature
  !!            compute the total volume
  !!
  !! history :
  !!   Original :  F. Castruccio (Fall 2006)
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
  INTEGER   :: ji, jj, jk, jclass, jiso, jbin, jarg   !: dummy loop index
  INTEGER   :: ipos                                   !: working variable
  INTEGER   :: narg, iargc                            !: command line 
  INTEGER   :: npiglo,npjglo                          !: size of the domain
  INTEGER   :: npk, nk                                !: vertical size, number of wet layers in the section
  INTEGER   :: numbimg=10                             !: optional bimg logical unit
  INTEGER   :: numout=11                              !: ascii output

  INTEGER                            :: imin, imax, jmin, jmax      !: working box limits

  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE :: gdept, gdepw !: depth of T and W points 
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE :: e3t          !: depth of T and W points 
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: e1t, e2t     !: lon, lat of T from file 
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zt           !: temperature from file
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zssh         !: SSH from file 
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: tmp          !: temporary array

  ! double precision for cumulative variables
  REAL(KIND=8), DIMENSION (:,:),   ALLOCATABLE ::  e1, e2             !: either e1t or e2t
  REAL(KIND=8), DIMENSION (:,:),   ALLOCATABLE ::  e3 , zmask         !: e3 and zmask
  REAL(KIND=8), DIMENSION (:,:),   ALLOCATABLE ::  ztemp, gdep        !: temp., depth of temp. points
  REAL(KIND=8)                                 :: temp_min, temp_max, dtemp    !: Min and Max for temp. bining
  REAL(KIND=8)                                 :: temp,zalfa                   !: current working temp.
  REAL(KIND=8), DIMENSION (:),   ALLOCATABLE   :: temp_lev                     !: built array with temp. levels
  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE   :: hiso                         !: depth of isotherms

  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE   :: zwvol, zwvolbin, volbin2      !: volume arrays
  REAL(KIND=8), DIMENSION (:), ALLOCATABLE     :: volbin                        !: volume arrays

  CHARACTER(LEN=80) :: cfilet                                                   !: files name
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc'          !: coordinates files
  CHARACTER(LEN=80) :: cfilout='voltemp.txt'                                    !: output file
  CHARACTER(LEN=80) :: cdum                                                     !: dummy string

  LOGICAL    :: l_print=.FALSE.             !: flag  for printing additional results
  LOGICAL    :: l_print2=.FALSE.             !: flag  for printing additional results
  LOGICAL    :: l_bimg=.FALSE.              !: flag  for bimg output

  !!  * Initialisations

  !  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg < 6  ) THEN
     PRINT '(255a)',' Usage : cdftempvol-full  gridTfile  imin, imax, jmin, jmax temp_max temp_min nbins [options]'
     PRINT '(255a)','            imin, imax, jmin, jmax : horizontal limit of the box'
     PRINT '(255a)','            temp_max, temp_min : limit for temperature bining '
     PRINT '(255a)','                           nbins : number of bins to use '
     PRINT '(255a)','     Possible options :'
     PRINT '(255a)','         -print :additional output is send to std output'
     PRINT '(255a)','         -bimg : 2D (x=lat/lon, y=temp) output on bimg file for hiso, cumul trp, trp'
     PRINT '(255a)',' Files mesh_hgr.nc, mesh_zgr.nc must be in the current directory'
     PRINT '(255a)',' Output on voltemp.txt'
     STOP
  ENDIF

  !! Read arguments
  CALL getarg (1, cfilet)
  CALL getarg (2,cdum)           ;        READ(cdum,*) imin
  CALL getarg (3,cdum)           ;        READ(cdum,*) imax
  CALL getarg (4,cdum)           ;        READ(cdum,*) jmin
  CALL getarg (5,cdum)           ;        READ(cdum,*) jmax
  CALL getarg (6,cdum)           ;        READ(cdum,*) temp_max
  CALL getarg (7,cdum)           ;        READ(cdum,*) temp_min
  CALL getarg (8,cdum)           ;        READ(cdum,*) nbins

  DO jarg=9, narg
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

  ! Allocate and build temp. levels and section array
  ALLOCATE ( temp_lev (nbins+1) )

  temp_lev(1)=temp_max
  dtemp=( temp_max - temp_min) / nbins
  DO jclass =2, nbins+1
     temp_lev(jclass)= temp_lev(1) - (jclass-1) * dtemp
  END DO

  ! Look for size of the domain
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk = getdim (cfilet,'depth')
  IF (imin /= 0 ) THEN ; npiglo=imax -imin + 1;  ELSE ; imin=1 ; ENDIF
  IF (jmin /= 0 ) THEN ; npjglo=jmax -jmin + 1;  ELSE ; jmin=1 ; ENDIF
  ALLOCATE ( gdept(npk), gdepw(npk), e1t(npiglo,npjglo), e2t(npiglo,npjglo), e3t(npk) )
  ALLOCATE ( volbin(nbins), volbin2(npjglo,nbins) )
  volbin=0.d0 ; volbin2=0.d0 

  ! read dimensions 
  gdept(:) = getvare3(coordzgr, 'gdept',npk)
  gdepw(:) = getvare3(coordzgr, 'gdepw',npk)
  e1t(:,:) = getvar(coordhgr, 'e1t', 1,npiglo,npjglo,kimin=imin,kjmin=jmin)
  e2t(:,:) = getvar(coordhgr, 'e2t', 1,npiglo,npjglo,kimin=imin,kjmin=jmin)
  e3t(:)   = getvare3(coordzgr, 'e3t',npk)


  !! *  Main loop
  
  DO jj=jmin,jmax

     ALLOCATE ( zt(npiglo,npk), tmp(npiglo,1), ztemp(npiglo,0:npk), zssh(npiglo,1) ) 
     ALLOCATE ( e1(npiglo,npk), e2(npiglo,npk), e3(npiglo,npk), gdep(npiglo, npk), zmask(npiglo,npk) )
     ALLOCATE ( zwvol(npiglo, nbins+1) , hiso(npiglo,nbins+1), zwvolbin(npiglo,nbins) )

     zssh= 0. ; gdep= 0. ; zmask = 0. ; ztemp=0.d0 ; e1=0.d0 ; e2=0.d0 ; e3=0.d0
     zwvol=0.d0 ; zwvolbin=0.d0
 
     zssh(:,:)=getvar(cfilet,'sossheig',1, npiglo, 1 , kimin=imin+1 , kjmin=jj)

     DO jk=1,npk
        ! initiliaze gdep to gdept()
        gdep(:,jk) = gdept(jk)

        ! metrics (Full step case)
        e1(:,jk)=e1t(:,jj-jmin+1)
        e2(:,jk)=e2t(:,jj-jmin+1)
        e3(:,jk)=e3t(jk)
     

        ! temperature
        tmp(:,:)=getvar(cfilet,'votemper',jk, npiglo, 1, kimin=imin+1, kjmin=jj)
        zmask(:,jk)=tmp(:,1)
        WHERE ( zmask(:,jk) /= 0 ) zmask(:,jk)=1

        zt(:,jk) = tmp(:,1)

        ! limitation to 'wet' points
        IF ( SUM(zt(:,jk))  == 0 ) THEN
           nk=jk ! first vertical point of the section full on land
           EXIT  ! as soon as all the points are on land
        ENDIF

     END DO

     ! temp. only for wet points
     ztemp(:,1:nk)=zt(:,:)
     ztemp(:,0)=ztemp(:,1)-1.e-4   ! dummy layer for easy interpolation


     ! Some control print 
     IF ( l_print2 ) THEN
        PRINT *,' T (deg C)' 
        DO jk=1,nk
           PRINT 9000, jk,  (ztemp(ji,jk),ji=1,npiglo)
        END DO

        PRINT *,' SSH (m)'  
        PRINT 9000, 1, (zssh(ji,1),ji=1,npiglo) 

        PRINT *,' GDEP (m) '
        DO jk=1,nk
           PRINT 9001,jk,  (gdep(ji,jk)*zmask(ji,jk),ji=1,npiglo)
        END DO

        PRINT *, 'E1 (m)'
        DO jk=1,nk
           PRINT 9001,jk,  (e1(ji,jk)*zmask(ji,jk),ji=1,npiglo)
        END DO

        PRINT *, 'E2 (m)'
        DO jk=1,nk
           PRINT 9001,jk,  (e2(ji,jk)*zmask(ji,jk),ji=1,npiglo)
        END DO

        PRINT *, 'E3 (m)'
        DO jk=1,nk
           PRINT 9001,jk,  (e3(ji,jk)*zmask(ji,jk),ji=1,npiglo)
        END DO
     END IF

     ! compute depth of isotherms (nbins+1 )
     IF (l_print )  PRINT *,' DEP ISO ( m )'
     DO  jiso =1, nbins+1
        temp=temp_lev(jiso)
!!!  REM : I and K loop can be inverted if necessary
        DO ji=1,npiglo
           hiso(ji,jiso) = gdept(npk)
           DO jk=1,nk 
              IF ( ztemp(ji,jk) > temp ) THEN
              ELSE
                 ! interpolate between jk-1 and jk
                 zalfa=(temp - ztemp(ji,jk-1)) / ( ztemp(ji,jk) -ztemp(ji,jk-1) )
                 IF (ABS(zalfa) > 1.1 ) THEN   ! case ztemp(0) = ztemp(1)-1.e-4
                    hiso(ji,jiso)= 0.
                 ELSE
                    hiso(ji,jiso)= gdep(ji,jk)*zalfa + (1.-zalfa)* gdep(ji,jk-1)
                 ENDIF
                 EXIT
              ENDIF
           END DO
        END DO
        IF (l_print) PRINT 9002, temp,(hiso(ji,jiso),ji=1,npiglo)
     END DO

     ! compute volume between surface and isotherm 
     IF (l_print) PRINT *,' VOL SURF -->  ISO (1.e12 M3)'
     DO jiso = 1, nbins + 1
        temp=temp_lev(jiso)
        DO ji=1,npiglo
           !zwvol(ji,jiso) = e1(ji,1)*e2(ji,1)*zssh(ji,1)
           DO jk=1, nk
              IF ( gdepw(jk+1) < hiso(ji,jiso) ) THEN
                 zwvol(ji,jiso)= zwvol(ji,jiso) + e1(ji,jk)*e2(ji,jk)*e3(ji,jk)
              ELSE  ! last box ( fraction)
                 zwvol(ji,jiso)= zwvol(ji,jiso) + e1(ji,jk)*e2(ji,jk)*(hiso(ji,jiso)-gdepw(jk))
                 EXIT  ! jk loop
              ENDIF
           END DO
        END DO
        IF (l_print) PRINT  9003, temp,(zwvol(ji,jiso)/1.e12,ji=1,npiglo)
     END DO

     ! binned volume : difference between 2 isotherms
     IF (l_print) PRINT *,' VOL bins (SV)'
     DO jbin=1, nbins
        temp=temp_lev(jbin)
        DO ji=1, npiglo
           zwvolbin(ji,jbin) = zwvol(ji,jbin+1) -  zwvol(ji,jbin) 
        END DO
        volbin2(jj-jmin+1,jbin)=SUM(zwvolbin(:,jbin) )
        IF (l_print) PRINT  9003, temp,(zwvolbin(ji,jbin)/1.e12,ji=1,npiglo), volbin2(jj-jmin+1,jbin)/1.e12
        volbin(jbin)=volbin(jbin)+volbin2(jj-jmin+1,jbin)
     END DO
     PRINT *,' Total volume in all bins (1e.15 M3):',SUM(volbin2(jj-jmin+1,:) )/1.e15
     

!    ! output of the code for 1 section
!    IF (l_bimg) THEN
!      ! (along section, depth ) 2D variables
!       cdum='Tdep.bimg'
!       OPEN(numbimg,FILE=cdum,FORM='UNFORMATTED')
!       cdum=' 3 dimensions in this file '
!       WRITE(numbimg) cdum
!       cdum=' 1: T  '
!       WRITE(numbimg) cdum
!       WRITE(cdum,'(a,4i5.4)') 'in box ', imin,imax,jmin,jmax
!       WRITE(numbimg) cdum
!       cdum=' file '//TRIM(cfilet)
!       WRITE(numbimg) cdum
!       WRITE(numbimg) npiglo,nk,1,1,2,0
!       WRITE(numbimg) 1.,-float(nk),1.,1., 0.
!       WRITE(numbimg) 0.
!       WRITE(numbimg) 0.
!       ! temperature
!       WRITE(numbimg) (( REAL(ztemp(ji,jk)), ji=1,npiglo) , jk=nk,1,-1 )
!       CLOSE(numbimg)

!       ! (along section, temp ) 2D variables
!       cdum='Volume_water_Tdep.bimg'
!       OPEN(numbimg,FILE=cdum,FORM='UNFORMATTED')
!       cdum=' 3 dimensions in this file '
!       WRITE(numbimg) cdum
!       cdum=' 1: hiso ;  2: bin vol '
!       WRITE(numbimg) cdum
!       WRITE(cdum,'(a,4i5.4)') ' in box ', imin,imax,jmin,jmax
!       WRITE(numbimg) cdum
!       cdum=' file '//TRIM(cfilet)
!       WRITE(numbimg) cdum
!       WRITE(numbimg) npiglo,nbins,1,1,2,0
!       WRITE(numbimg) 1.,-REAL(temp_lev(nbins)),1.,REAL(dtemp), 0.
!       WRITE(numbimg) 0.
!       WRITE(numbimg) 0.
!       ! hiso
!       WRITE(numbimg) (( REAL(hiso(ji,jiso)), ji=1,npiglo) , jiso=nbins,1,-1)
!       ! binned transport
!       WRITE(numbimg) (( REAL(zwvolbin(ji,jiso))/1.e15, ji=1,npiglo) , jiso=nbins,1,-1)
!       CLOSE(numbimg)
!    ENDIF

     ! free memory for the next section
     DEALLOCATE ( zt, tmp, ztemp, zssh )
     DEALLOCATE ( e1, e2, e3, gdep, zmask )
     DEALLOCATE ( zwvol, hiso, zwvolbin )

     PRINT *,' Total volume in all bins (1e.15 M3):',SUM(volbin(:)/1.e15 )

  END DO   ! next section

  !! Global Output
     OPEN( numout, FILE=cfilout)
       ipos=INDEX(cfilet,'_gridT.nc')
       WRITE(numout,9006)  TRIM(cfilet(1:ipos-1))
       WRITE(numout,*) ' temp.  '
     DO jiso=1,nbins
       WRITE(numout,9004) temp_lev(jiso), volbin(jiso)
     ENDDO
     CLOSE(numout)

9000 FORMAT(i7,60f8.3)
9001 FORMAT(i7,60f8.0)
9002 FORMAT(f7.3,60f8.0)
9003 FORMAT(f7.3,60f8.3)
9004 FORMAT(f9.4, 60e16.7)
9005 FORMAT('#',a9, 60(2x,a12,2x) )
9006 FORMAT('# ',a)


END PROGRAM cdftempvol_full
