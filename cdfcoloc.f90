PROGRAM cdfcoloc
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfcoloc  ***
  !!
  !!  **  Purpose  :  colocalisation for Greg Holloway
  !!  
  !!  **  Method   :  Use the weight files computed with cdfweight
  !!
  !! history ;
  !!  Original :  J.M. Molines (16 may 2007 )
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER, PARAMETER :: jptyp=5
  INTEGER :: narg, iargc, nid=0
  INTEGER :: ji,jj, jk, jid, jtyp
  INTEGER :: i1, j1, i2, j2, i3, j3, i4, j4, k1, k2
  INTEGER :: imin,  jmin
  INTEGER :: iloc, jloc, kloc
  INTEGER :: npiglo, npjglo, iquadran, npk, npkv
  INTEGER :: numgreg=10, ios, id, idep, numbin=20, numout=30

  REAL(KIND=8)  :: alpha, beta, gamma, hN, scale
  REAL(KIND=8)  :: xmin, ymin

  REAL(KIND=4)                      :: vup, vdo, wup, wdo
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: v3d
  REAL(KIND=8), DIMENSION(:,:)  , ALLOCATABLE :: vinterp, v2d,  e

  INTEGER(KIND=2), DIMENSION(:,:,:), ALLOCATABLE :: mask

  CHARACTER(LEN=80) :: cdum
  ! file name
  CHARACTER(LEN=80) :: coord='coordinates.nc', czgr='mesh_zgr.nc', cmask='mask.nc'
  CHARACTER(LEN=80) :: cweight,  cgridt, cgridu, cgridv, cctyp, cfil, cvar, cvmask
  CHARACTER(LEN=3), DIMENSION(jptyp) :: ctype

  !!
  ctype=(/'U','V','Sx','Sy','H'/)
  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg /= 4  ) THEN
     PRINT *,' Usage : cdfcoloc  root_weight gridT gridU gridV '
     PRINT *,' return an ascii file iyxzUVSxSyH.txt'
     PRINT *,' Example : cdfcoloc weight gridT gridU gridV '
     STOP
  ENDIF

  CALL getarg (1, cweight )
  CALL getarg (2, cgridt )
  CALL getarg (3, cgridu )
  CALL getarg (4, cgridv )

  npiglo= getdim (cgridt,'x')
  npjglo= getdim (cgridt,'y')
  npk=    getdim (cgridt,'depth')

  ALLOCATE (v3d(npiglo, npjglo, npk), mask(npiglo, npjglo, npk) )
  ALLOCATE (v2d(npiglo, npjglo), e(npiglo,npjglo) )

  WRITE(cweight,'("weight_",a,".bin")') 'T'
  OPEN(numbin, FILE=cweight,FORM='unformatted')
  ! Determine the number of records in the weight file
  DO WHILE (1 /= 2 )
     READ(numbin, END=100)
     nid=nid+1
  END DO
100 CONTINUE
  PRINT *, nid
  CLOSE (numbin)
  ! ALLOCATE ( T(nid), S(nid), U(nid), V(nid), SSH(nid), H(nid), Sx(nid),Sy(nid) )
  ALLOCATE ( vinterp(nid,jptyp) )

  DO jtyp=1,jptyp
     cctyp=ctype(jtyp)
     SELECT CASE ( cctyp)
     CASE ('T')
        WRITE(cweight,'("weight_",a,".bin")') 'T'
        cvar='votemper' ; cvmask='tmask'
        cfil=cgridt
        npkv=npk
        scale=1.
     CASE ('S')
        WRITE(cweight,'("weight_",a,".bin")') 'T'
        cvar='vosaline' ; cvmask='tmask'
        cfil=cgridt
        npkv=npk
        scale=1.
     CASE ('U')
        WRITE(cweight,'("weight_",a,".bin")') 'U'
        cvar='vozocrtx' ; cvmask='umask'
        cfil=cgridu
        npkv=npk
        scale=100.
     CASE ('V')
        WRITE(cweight,'("weight_",a,".bin")') 'V'
        cvar='vomecrty' ; cvmask='vmask'
        cfil=cgridv
        npkv=npk
        scale=100.
     CASE ('SSH')
        WRITE(cweight,'("weight_",a,".bin")') 'T'
        cvar='sossheig' ; cvmask='tmask'
        cfil=cgridt
        npkv=1
        scale=100.
     CASE ('H')
        WRITE(cweight,'("weight_",a,".bin")') 'T'
        cvar='hdepw' ; cvmask='tmask'
        cfil=czgr
        npkv=1
        scale=1.
     CASE ('Sx') 
        WRITE(cweight,'("weight_",a,".bin")') 'U'
        cfil='none' ; cvmask='umask'
        npkv=1
        e(:,:)=getvar(coord,'e1u',1, npiglo,npjglo)     
        v2d(:,:)=getvar(czgr,'hdepw',1, npiglo,npjglo)     
        DO ji=2, npiglo-1
            v3d(ji,:,1)=(v2d(ji+1,:) -v2d(ji,:))/ e(ji,:)
        END DO
        scale=100.
     CASE ('Sy')
        WRITE(cweight,'("weight_",a,".bin")') 'V'
        cfil='none' ; cvmask='vmask'
        npkv=1
        e(:,:)=getvar(coord,'e2v',1, npiglo,npjglo)
        v2d(:,:)=getvar(czgr,'hdepw',1, npiglo,npjglo)
        DO jj=2, npjglo-1
            v3d(:,jj,1)=(v2d(:,jj+1) -v2d(:,jj))/ e(:,jj)
        END DO
        scale=100.
     END SELECT

     OPEN(numbin, FILE=cweight,FORM='unformatted')


     IF (cfil /= 'none' ) THEN
      DO jk=1, npkv
        v3d(:,:,jk)=getvar(cfil,cvar,jk, npiglo,npjglo)
      END DO
     ENDIF

     DO jk=1, npkv
        mask(:,:,jk)=getvar(cmask,cvmask,jk, npiglo,npjglo)
     END DO

     REWIND(numbin)
     PRINT *,'START coloc ...'
     DO jid=1,nid
        READ(numbin) id, ymin, xmin, idep ,imin, jmin, kloc, iquadran, hN, alpha, beta, gamma
        vinterp(jid,jtyp)=interp()
        IF ( vinterp (jid,jtyp) > -990 )  vinterp (jid,jtyp) = vinterp (jid,jtyp) * scale
     END DO

     CLOSE(numbin)
  END DO
    OPEN(numout,FILE='iyxzUVSxSyH.txt')
!       PRINT '(I5, I10, 2f9.4, I5, f10.3,1x,a,1x,I3)',jid, id, ymin,xmin,idep, vinterp(jid,jtyp),TRIM(cctyp),kloc
     WRITE(cweight,'("weight_",a,".bin")') 'T'
     OPEN(numbin, FILE=cweight,FORM='unformatted')
     DO jid=1, nid
        READ(numbin) id, ymin, xmin, idep, imin, jmin, kloc, iquadran, hN
        IF ( xmin > 180.) xmin=xmin -360
        CALL rotation( vinterp(jid,:), hN)
        WRITE(numout, '(I5, 2f10.4,    I6,       5f10.4)')  &
    &                   id, ymin,xmin,idep, (vinterp(jid,jtyp),jtyp=1,jptyp)
     END DO
     CLOSE(numbin)
     CLOSE(numout)

CONTAINS

  FUNCTION interp
    REAL(KIND=8) :: interp

    IF (imin == -1000 ) THEN
      interp=-999.
      RETURN
    ENDIF

    SELECT CASE (iquadran)
    CASE (1)
       i1=imin    ; j1 = jmin
       i2=imin +1 ; j2 = jmin
       i3=imin +1 ; j3 = jmin + 1
       i4=imin    ; j4 = jmin + 1
    CASE (2)
       i1=imin    ; j1 = jmin
       i2=imin    ; j2 = jmin - 1
       i3=imin +1 ; j3 = jmin - 1
       i4=imin +1 ; j4 = jmin
    CASE (3)
       i1=imin    ; j1 = jmin
       i2=imin -1 ; j2 = jmin
       i3=imin -1 ; j3 = jmin - 1
       i4=imin    ; j4 = jmin - 1
    CASE (4)
       i1=imin    ; j1 = jmin
       i2=imin    ; j2 = jmin + 1
       i3=imin -1 ; j3 = jmin + 1
       i4=imin -1 ; j4 = jmin
    END SELECT
    k1=kloc    ; k2 = kloc + 1

    IF (npkv == 1 ) THEN
       k1 = 1  ; k2 = 0 ; wdo = 0.
    ENDIF

    wup = mask(i1,j1,k1)*(1-alpha)*(1-beta) + mask(i2,j2,k1) * alpha *(1-beta) + &
         &  mask(i3,j3,k1)*alpha*beta         + mask(i4,j4,k1) * (1-alpha)*beta

    vup = v3d(i1,j1,k1)*(1-alpha)*(1-beta) + v3d(i2,j2,k1) * alpha *(1-beta) + &
         &  v3d(i3,j3,k1)*alpha*beta         + v3d(i4,j4,k1) * (1-alpha)*beta

    IF (k2 /= 0 ) THEN
       wdo = mask(i1,j1,k2)*(1-alpha)*(1-beta) + mask(i2,j2,k2) * alpha *(1-beta) + &
            &  mask(i3,j3,k2)*alpha*beta         + mask(i4,j4,k2) * (1-alpha)*beta

       vdo = v3d(i1,j1,k2)*(1-alpha)*(1-beta) + v3d(i2,j2,k2) * alpha *(1-beta) + &
            &  v3d(i3,j3,k2)*alpha*beta         + v3d(i4,j4,k2) * (1-alpha)*beta
    ENDIF

    IF ( wup == 0 ) THEN
       interp=-999.
    ELSE IF ( wdo == 0 ) THEN
       interp= vup/wup
    ELSE
       interp= (1 - gamma) * vup/wup + gamma * vdo/wdo
    ENDIF

  END FUNCTION interp

  SUBROUTINE rotation (pval, pcourse)
    REAL(KIND=8), DIMENSION(jptyp), INTENT(inout) :: pval
    REAL(KIND=8), INTENT(in) :: pcourse

   ! local variable
   REAL(KIND=8), DIMENSION(2) :: zcomp 
   REAL(KIND=8) :: zconv, zcourse, zpi
  
   
   zpi=ACOS(-1.d0) ; zconv=zpi/180.d0
   zcourse=pcourse*zconv  ! heading in radians
   ! Velocity : u=1, v=2
   zcomp(1:2)=pval(1:2)
   IF ( zcomp(1) > -990. .AND. zcomp(2) > -990. ) THEN
     pval(1) =  zcomp(1)*cos(zcourse) +zcomp(2)*sin(zcourse)
     pval(2) = -zcomp(1)*sin(zcourse) +zcomp(2)*cos(zcourse)
   ENDIF

   ! Bottom slope : Sx=3, Sy=4
   zcomp(1:2)=pval(3:4)
   IF ( zcomp(1) > -990. .AND. zcomp(2) >  -990. ) THEN
     pval(3) =  zcomp(1)*cos(zcourse) +zcomp(2)*sin(zcourse)
     pval(4) = -zcomp(1)*sin(zcourse) +zcomp(2)*cos(zcourse)
   ENDIF

  END SUBROUTINE rotation

END PROGRAM cdfcoloc
