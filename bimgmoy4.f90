PROGRAM bimgmoy4
!!
!!! ----------------------------------------------------------
!!!                  PROGRAM BIMGMOY4
!!!                  ****************
!!!
!!!     PURPOSE:
!!!     --------
!!!     This program read a list of bimg files (either sequential
!!!     or direct files) and compute the mean value which is dumped
!!!     to the file named moy.bimg (or moy.dimg). It also computes the
!!!     second order moment (ie the mean square of the files), and the
!!!     variance.
!!!
!!      METHOD:
!!      -------
!!      This program assume that all the files have the same geometry.
!!      It also assumes that if the first file is a direct file all
!!      files are direct. The variance is computed as 
!!        var = (mean_squared) - (mean) squared
!!       File can contain time frames (as many as you want) BUT
!!      they must be 3D files, ie you cannot have both nk AND ndim non zero.
!!
!!      USAGE:
!!      ------
!!        bimgmoy2 'file list'
!!     
!!      OUTPUT:
!!      -------
!!      On moy.bimg/moy.dimg , moy2.bimg/moy2.dimg AND var.bimg/var.dimg
!!      
!!      EXTERNAL:
!!      --------
!!        isdirect : normally in libbimg.a
!!        iargc, getarg  : libU77
!!     
!!      AUTHOR:
!!      ------
!!           J.M. Molines, May 1998
!!!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
!! 0.0 Declarations:
!! -----------------
        IMPLICIT NONE
!
        INTEGER narg, iargc,  nrecl
        INTEGER ni, nj, nk, ndim, icod, nt, irec, nk2, nt2, nframe
        INTEGER ji, jj, jk, jt, jfich, irecl
!
	CHARACTER*80 cline1, cline2, cline3, cline4
	CHARACTER*80 clfil1
        CHARACTER*4 VER
!
	REAL(KIND=4),DIMENSION(:,:), ALLOCATABLE :: v2d
	REAL(KIND=4),DIMENSION(:), ALLOCATABLE ::  h1d, time_tag
        REAL x1, y1, dx, dy, spval, time_tag1,time_mean
!
! ... REAL*8 are necessary on workstations to avoid truncation error
!     in the variance computation.
        REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: v3d
!
        LOGICAL lbimgflag
!!
!! 1.0 Initialisations:
!! --------------------
	narg=iargc()
	IF (narg .LT. 1) THEN
	 print *,'USAGE: bimgmoy4  file* [-bimg]'
	 print *,' sortie sur moy.bimg | moy.dimg'
	 print *,' sortie sur moy2.bimg|moy2.dimg'
         print *,' If -bimg is specified (optional) the output '
         print *,'  is forced to be bimg regardless of the input format'
	STOP
	END IF
        lbimgflag = .false.
! ... check for -bimg option 
	CALL getarg(narg,clfil1)
       IF (clfil1 .EQ. '-bimg' ) THEN
         narg      = narg - 1
         lbimgflag = .true.
       END IF
        CALL getarg(1,clfil1)

        nrecl=ISDIRECT(clfil1)
       IF ( nrecl .EQ. 0 ) THEN

       OPEN(10,FILE=clfil1,FORM='UNFORMATTED')
       READ(10) cline1
       READ(10) cline2
       READ(10) cline3
       READ(10) cline4
!
       READ(10) ni,nj,nk,nt,ndim,icod
       PRINT *, ni,nj,nk,nt,ndim,icod
       CLOSE(10)
       ELSE
       OPEN(10,FILE=clfil1,FORM='UNFORMATTED',  ACCESS='DIRECT', RECL=nrecl)
       READ(10,REC=1) VER,cline1,irecl, ni,nj,nk,nt,ndim
       CLOSE(10)
       ENDIF
       nk=max(nk,ndim)


! ALLOCATE ARRAYS
       ALLOCATE( v2d(ni,nj), v3d(ni,nj,nk), h1d(nk), time_tag(nt) )
 
!      v3d=0.d0
       time_mean = 0.
       nframe = 0
!
!! 2.0 Loop on files, accumulate values and squared values
!! -------------------------------------------------------
!!
       DO jfich = 1, narg
       CALL getarg(jfich,clfil1)
       PRINT *, ni,nj,nk,nt,ndim,icod
       print *,trim(clfil1)
! ... check if the file is direct or not
        nrecl=ISDIRECT(clfil1)
       IF ( nrecl .EQ. 0 ) THEN
!
! ... The file is a bimg file ...
! 
       OPEN(10,FILE=clfil1,FORM='UNFORMATTED')
       rewind(10)
       READ(10) cline1
       READ(10) cline2
       READ(10) cline3
       READ(10) cline4
       print *, trim(cline1)
       print *, trim(cline2)
       print *, trim(cline3)
       print *, trim(cline4)
!
!
       READ(10) ni,nj,nk,nt,ndim,icod
! ... Stop if the file hold more than 1 time frame or more than 1 dim
       IF ( (ndim .NE. 1) .AND.  (nk .NE. 1)) THEN
        print *,' This program only works with files'
        print *,' having both nk and ndim not zero'
        print *,' Sorry .... :( '
        STOP 
       END IF
       IF (ndim .NE. 1) THEN
        nk = ndim
        nk2=1
       ELSE
        nk2=nk
       ENDIF
!
       READ(10) x1,y1,dx,dy,spval
       READ(10) (h1d(jk),jk=1,nk2)
       DO jt =1,nt
       READ(10) time_tag(jt)
       IF (jt .EQ. 1 .AND. jfich .EQ. 1 ) time_tag1=time_tag(1)
       time_mean = time_mean + time_tag(jt)/nt
!
        DO jk=1,nk
          READ(10)((v2d(ji,jj),ji=1,ni),jj=1,nj)
           DO ji=1,ni
            DO jj=1,nj
              IF (v2d(ji,jj) .NE. spval) THEN
             v3d  (ji,jj,jk) = v3d  (ji,jj,jk) + dble(v2d(ji,jj))
              ELSE
             v3d  (ji,jj,jk) = dble(spval)
              END IF
            END DO
           END DO
        END DO
	 nframe = nframe + 1
       END DO
       CLOSE(10)
       ELSE
!
! ... The file is a dimg file
!
       OPEN(10,FILE=clfil1,FORM='UNFORMATTED', ACCESS='DIRECT', RECL=nrecl)
       READ(10,REC=1) VER,cline1,irecl,  &
     &         ni,nj,nk,nt,ndim,  &
     &         x1,y1,dx,dy,spval,  &
     &         (h1d(jk),jk=1,nk),  &
     &         (time_tag(jt),jt=1,nt)
       IF ( (ndim .NE. 1) .AND.  (nk .NE. 1)) THEN
        print *,' This program only works with files'
        print *,' having both nk and ndim not zero'
        print *,' Sorry .... :( '
        STOP 
       END IF
       IF (ndim .NE. 1)  nk = ndim

       DO jt = 1, nt
       IF (jt .EQ. 1 .AND. jfich .EQ. 1 ) time_tag1=time_tag(1)
       time_mean = time_mean + time_tag(jt)/nt
        DO jk=1,nk
          irec = 2 + (jt -1)*nk + (jk -1 )
         READ(10,REC=irec)((v2d(ji,jj),ji=1,ni),jj=1,nj)
           DO ji=1,ni
            DO jj=1,nj
              IF (v2d(ji,jj) .NE. spval) THEN
             v3d  (ji,jj,jk) = v3d  (ji,jj,jk) + dble(v2d(ji,jj))
              ELSE
             v3d  (ji,jj,jk) = dble(spval)
              END IF
            END DO
           END DO
        END DO
	 nframe = nframe + 1
       END DO
       CLOSE(10)
       END IF
! ... Loop on files
       END DO
!
!!
!! 3.0 Compute mean value and mean squared value, and Variance
!! ---------------------------------------------
       DO jk=1,nk
       DO ji=1,ni
       DO jj=1,nj
        IF(v3d(ji,jj,jk) .NE. spval) THEN
           v3d  (ji,jj,jk)=v3d  (ji,jj,jk) / float(nframe)
        END IF
       END DO
       END DO
       END DO
!!
!! 4.0 Output to bimg or dimg file (depending of the input files)
!!     The bimg format can be forced by the option -bimg
!! --------------------------------------------------------------
       IF (ndim .NE. 1 ) THEN 
         nk = 1
         nk2 = ndim 
       ELSE 
         nk2 = nk 
       END IF
! ... There is only one time frame in the output
	 nt2 = nt
         nt  = 1
       time_mean=time_mean/narg
       cline3=cline3(2:)
	IF (nrecl .EQ. 0 .OR. lbimgflag ) THEN
        OPEN(10,FILE='moy.bimg',FORM='UNFORMATTED')
!
! ... Mean file
!
       WRITE(cline1,100) nframe
100    FORMAT('MEAN  values from ',i3.3,' dumps')
       WRITE(cline2,101)  time_tag1, time_tag(nt2)
101       FORMAT('computed between day ',f8.0,' and day ',f8.0)
       WRITE(cline4,'(16hMade by bimgmoy2)')
       WRITE(10) cline1
       WRITE(10) cline2
       WRITE(10) cline3
       WRITE(10) cline4
       WRITE(10) ni,nj,nk,nt,ndim,icod
       WRITE(10) x1,y1,dx,dy,spval
       WRITE(10)(h1d(jk),jk=1,nk)
       WRITE(10) time_mean
       DO jk=1,nk2
       WRITE(10)((REAL(v3d(ji,jj,jk)),ji=1,ni),jj=1,nj)
       END DO
       CLOSE(10)

       ELSE
! ... DIMG file
      OPEN(10,FILE='moy.dimg',FORM='UNFORMATTED', ACCESS='DIRECT', RECL = nrecl)
!
! ... Mean file
!
       WRITE(cline1,103) nframe,time_tag1, time_tag(nt2)
103    FORMAT('MEAN  values from ',i3.3,' dumps (',f8.0,'->',f8.0,')')
       WRITE(10,REC=1) VER,cline1,nrecl,   &
     &         ni,nj,nk,nt,ndim,   &
     &         x1,y1,dx,dy,spval,   &
     &         (h1d(jk),jk=1,nk),   &
     &         (time_mean,jt=1,nt)
       DO jk=1,nk2
        irec = jk+1
       WRITE(10,REC=irec)((REAL(v3d(ji,jj,jk)),ji=1,ni),jj=1,nj)
       END DO
       CLOSE(10)
       ENDIF

       print *,' bimgmoy2 successfull!'

    CONTAINS
        FUNCTION isdirect(clname)
!!! -------------------------------------------------------------------------
!!!                     FUNCTION ISDIRECT
!!!                     *****************
!!!
!!!    PURPOSE : This integer function returns the record length if clname
!!!              is a valid dimg file, it returns 0 either.
!!!
!!!    METHOD : Open the file and look for the key characters (@!01) for
!!!             identification.
!!!
!!!    AUTHOR : Jean-Marc Molines (Apr. 1998)
!!! -------------------------------------------------------------------------
!! 1.0 Declarations:
!! -----------------
      IMPLICIT NONE
      INTEGER isdirect
      CHARACTER*(*) clname
      CHARACTER*4 VER
      CHARACTER*80 clheader
!
      INTEGER irecl
!!
!! 2.0 Look for VER:
!! ----------------
!!
      OPEN(100,FILE=clname,   &
     &         FORM   ='UNFORMATTED',   &
     &         ACCESS ='DIRECT',   &
     &         RECL   =88)
       READ(100,REC=1) VER,clheader,irecl
	print *,'VER',VER
      CLOSE(100)
!
      IF (VER .EQ. '@!01' ) THEN
       isdirect=irecl
      ELSE
       isdirect=0
      END IF
!
      END FUNCTION isdirect

 END PROGRAM bimgmoy4


