PROGRAM bimgcaltrans
  !!--------------------------------------------------------------
  !!             *** PROGRAM  bimgcaltrans  ***
  !! 
  !!   ** Purpose: Compute density class transport from bimg files
  !!               produced by cdfsigtrp
  !!
  !!   History :
  !!       Original : J.M. Molines (22/03/2006)
  !!--------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: npi,npj,npk,npt, npdim
  INTEGER :: ji,jj,jdim
  INTEGER :: narg, iargc
  INTEGER :: numin=10

  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: v2d
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: sig
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: trp
  REAL(KIND=4) :: x1, sigmin, dx,dsig,spval, bidon

  CHARACTER(LEN=80) :: cfile, comm

  !! 
  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' USAGE : bimgcaltrans  bimg-file'
     STOP
  ENDIF

  CALL getarg(1,cfile)
  OPEN(numin,FILE=cfile,FORM='UNFORMATTED')
  READ(numin) comm
  READ(numin) comm
  READ(numin) comm
  READ(numin) comm
  READ(numin) npi,npj,npk,npt,npdim

  ALLOCATE ( v2d(npi,npj), sig(npj), trp(npj) )
  READ(numin) x1,sigmin, dx,dsig,spval
  READ(numin)   ! skip h
  READ(numin)   ! skip t
  READ(numin)   ! skip 1rst dim = hiso
  READ(numin) ((v2d(ji,jj),ji=1,npi), jj=1,npj)
  
  !! Build sig
  sig(1)=-sigmin
  DO jj=2,npj
    sig(jj)=sig(1)-(jj-1)*dsig
  END DO

  trp(:)=0.d0
  DO jj=1,npj
    trp(jj)=SUM(v2d(:,jj))
  END DO
  
  DO jj=npj,1,-1
    PRINT 9004, sig(jj),trp(jj)
  END DO

9004 FORMAT(f9.4, 20e16.7)

END PROGRAM  bimgcaltrans

