!
! 
PROGRAM coordinates2nc
  !!-------------------------------------------------------------------------
  !!                        PROGRAM coordinates2nc
  !!                        **********************
  !! ** Purpose:
  !!         transform a "clipper" coordinates file into an
  !!         ioipsl coordinate file
  !!
  !! History :
  !!   February 2003 : Anne de Miranda
  !!--------------------------------------------------------------------------
  USE ioipsl
  !  USE histcom

  IMPLICIT NONE


  INTEGER(4)           :: jpi , jpj
  INTEGER(4),PARAMETER ::  wp = 8
  INTEGER(4)           ::  numcoo, nummsh
  INTEGER(4)           ::  nrecl8
  INTEGER(4)           ::  ngrid, ncid
  INTEGER(4)           ::  ni,nj, nk
  INTEGER(4)           ::  ilev, itime, narg, iargc

  REAL(wp)             ::  zjpiglo,zjpjglo,znbsel,zrecl8
  REAL(wp)             :: znt,zdim,xx1,yy1,ddx,ddy,sspval

  REAL(wp), DIMENSION(:,:),ALLOCATABLE ::   &
       zlamt, zphit, zdta
  REAL(wp), DIMENSION(1):: zdep
  REAL(wp) zdt,zdate0 
  REAL(wp), DIMENSION(:,:),  ALLOCATABLE ::   &
       glamt, glamu, glamv, glamf, gphit, gphiu, gphiv, gphif, &
       e1t, e1u, e1v, e1f, e2t, e2u, e2v, e2f

  CHARACTER (len=21) ::   clname
  CHARACTER(80)        :: cltextco, cfilcoo, cfilout

  LOGICAL clog

  numcoo=10
  nummsh=11
  narg=iargc()
  IF ( narg /= 2 ) THEN
    PRINT *, ' USAGE: coordinates2nc ''clipper coordinate'' '' nc coordinates'' '
    STOP
  END IF

  CALL getarg(1, cfilcoo)
  CALL getarg(2,cfilout)
  !
  ! ... Read coordinates (only the used metric is read)
  !
  nrecl8=200
  OPEN(numcoo,FILE=cfilcoo,status='old' & 
       ,form='unformatted',      access='direct',recl=nrecl8)
  READ(numcoo,rec=1)      cltextco, zrecl8, zjpiglo, zjpjglo
  CLOSE(numcoo)
  !
  PRINT *, cltextco, zrecl8, zjpiglo, zjpjglo
  nrecl8=zrecl8
  PRINT*, nrecl8
  OPEN(numcoo,FILE=cfilcoo,status='old',form='unformatted', &
       access='direct',recl=nrecl8)
  READ (numcoo, rec = 1) cltextco, zrecl8, zjpiglo, zjpjglo, znbsel, znt, &
       zdim, xx1, yy1, ddx, ddy, sspval
  ni=zjpiglo
  nj=zjpjglo
  PRINT *, 'ni=',ni,'nj=',nj
  ALLOCATE(glamt(ni,nj))
  ALLOCATE(glamu(ni,nj))
  ALLOCATE(glamv(ni,nj))
  ALLOCATE(glamf(ni,nj))
  ALLOCATE(gphit(ni,nj))
  ALLOCATE(gphiu(ni,nj))
  ALLOCATE(gphiv(ni,nj))
  ALLOCATE(gphif(ni,nj))
  ALLOCATE(e1t(ni,nj))
  ALLOCATE(e1u(ni,nj))
  ALLOCATE(e1v(ni,nj))
  ALLOCATE(e1f(ni,nj))
  ALLOCATE(e2t(ni,nj))
  ALLOCATE(e2u(ni,nj))
  ALLOCATE(e2v(ni,nj))
  ALLOCATE(e2f(ni,nj))

  READ(numcoo,rec=2) glamt
  READ(numcoo,rec=3) glamu
  READ(numcoo,rec=4) glamv
  READ(numcoo,rec=5) glamf
  READ(numcoo,rec=6) gphit
  READ(numcoo,rec=7) gphiu
  READ(numcoo,rec=8) gphiv
  READ(numcoo,rec=9) gphif
  READ(numcoo,rec=10) e1t
  READ(numcoo,rec=11) e1u
  READ(numcoo,rec=12) e1v
  READ(numcoo,rec=13) e1f
  READ(numcoo,rec=14) e2t
  READ(numcoo,rec=15) e2u
  READ(numcoo,rec=16) e2v
  READ(numcoo,rec=17) e2f

  PRINT *, 'Reading ',TRIM(cfilcoo),' OK.'

  clname = cfilout
  ilev = 1
  itime=0
  zdate0=0.
  zdt=0.
  clog=.FALSE.

  ALLOCATE(zlamt(ni,nj))
  ALLOCATE(zphit(ni,nj))

  zlamt(:,:) = 0.
  zphit(:,:) = 0.
  zdep(1) = 0.

  jpi = zjpiglo
  jpj = zjpjglo
  PRINT* , jpi, jpj

  !  CALL histbeg ( clname, jpi, zlamt, jpj, zphit, &
  !       1, ni, 1, nj, 0, zdate0, zdt, ngrid, ncid)

  PRINT *,ngrid, ncid

  !  CALL restini(clname,jpi,jpj,zlamt,zphit,ilev,zdep,clname   &
  !       ,itime,zdate0,zdt,nummsh)
  CALL restini('NONE',jpi,jpj,glamt,gphit,1,zdep,TRIM(clname),itime,zdate0,zdt,nummsh)

  ! Horizontal grid-point position
  !-----------------------------
  CALL restput(nummsh,'glamt',jpi,jpj,1,0,glamt)
  CALL restput(nummsh,'glamu',jpi,jpj,1,0,glamu)
  CALL restput(nummsh,'glamv',jpi,jpj,1,0,glamv)
  CALL restput(nummsh,'glamf',jpi,jpj,1,0,glamf)
  !
  CALL restput(nummsh,'gphit',jpi,jpj,1,0,gphit)
  CALL restput(nummsh,'gphiu',jpi,jpj,1,0,gphiu)
  CALL restput(nummsh,'gphiv',jpi,jpj,1,0,gphiv)
  CALL restput(nummsh,'gphif',jpi,jpj,1,0,gphif)

  ! Horizontal scale factors
  ! ---------------------------
  CALL restput(nummsh,'e1t',jpi,jpj,1,0,e1t)
  CALL restput(nummsh,'e1u',jpi,jpj,1,0,e1u)
  CALL restput(nummsh,'e1v',jpi,jpj,1,0,e1v)
  CALL restput(nummsh,'e1f',jpi,jpj,1,0,e1f)
  !
  CALL restput(nummsh,'e2t',jpi,jpj,1,0,e2t)
  CALL restput(nummsh,'e2u',jpi,jpj,1,0,e2u)
  CALL restput(nummsh,'e2v',jpi,jpj,1,0,e2v)
  CALL restput(nummsh,'e2f',jpi,jpj,1,0,e2f)


END PROGRAM coordinates2nc
