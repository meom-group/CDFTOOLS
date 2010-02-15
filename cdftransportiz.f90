PROGRAM cdftransportiz
  !!---------------------------------------------------------------------
  !!               ***  PROGRAM cdftransportiz  ***
  !!
  !!  **  Purpose: Compute Transports across a section
  !!               PARTIAL STEPS version
  !!  
  !!  **  Method: Try to avoid 3 d arrays.
  !!             The begining and end point of the section are given in term of f-points index.
  !!             This program computes the transport across this section for
  !!               (1) Mass transport ( Sv)
  !!               (2) Heat Transport (PW)
  !!               (3) Salt Transport (kT/sec)
  !!             The transport is > 0 left handside of the line
  !!             This program use a zig-zag line going through U and V-points.
  !!             It takes as input : VT files, gridU, gridV files.
  !!             The mesh_hgr.nc, mesh_hzr.nc are required.
  !!             It is conveniebt to use an ASCII file as the standard input to give
  !!             the name and the imin imax jmin jmax for eaxh section required
  !!             The last name of this ASCII file must be EOF
  !!
  !!
  !! history :
  !!   Original :  J.M. Molines (jan. 2005)
  !!               J.M. Molines Apr 2005 : use modules
  !!               J.M. Molines Apr 2007 : merge with Julien Jouanno version (std + file output)
  !!               R. Dussin (Jul. 2009) : add cdf output
  !!---------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER :: nclass   !: number of depth class
  INTEGER ,DIMENSION (:),ALLOCATABLE ::  imeter  !: limit beetween depth level, in m (nclass -1)
  INTEGER ,DIMENSION (:),ALLOCATABLE :: ilev0,ilev1 !: limit in levels  ! nclass
  INTEGER   :: jk, jclass, jj                      !: dummy loop index
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: imin, imax, jmin, jmax, ik 
  INTEGER   :: numout = 10, numvtrp=11, numhtrp=12, numstrp=14
  ! added to write in netcdf
  INTEGER :: kx=1, ky=1, kz=1          ! dims of netcdf output file
  INTEGER :: nboutput=9                ! number of values to write in cdf output
  INTEGER :: ncout, ierr               ! for netcdf output
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipk, id_varout
  ! broken line stuff
  INTEGER, PARAMETER :: jpseg=10000
  INTEGER :: i0,j0,i1,j1, i, j
  INTEGER :: n,nn,k, jseg
  INTEGER :: norm_u, norm_v, ist, jst

  REAL(KIND=4) ::  rxi0,ryj0, rxi1, ryj1
  REAL(KIND=4) ::   ai,bi, aj,bj,d
  REAL(KIND=4) ::    rxx(jpseg),ryy(jpseg)
  REAL(KIND=4), DIMENSION(jpseg) :: gla, gphi

  REAL(KIND=8), DIMENSION(jpseg) :: voltrp, heatrp, saltrp
  REAL(KIND=8)                   :: voltrpsum, heatrpsum, saltrpsum
  COMPLEX yypt(jpseg), yypti

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         e1v, e3v ,gphiv, zv, zvt, zvs !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         e2u, e3u ,gphiu, zu, zut, zus !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         glamu, glamv
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE ::         gdepw
  REAL(KIND=4)                                 ::   rd1, rd2
  REAL(KIND=4)                                 ::  udum, vdum

  REAL(KIND=8),   DIMENSION (:,:), ALLOCATABLE :: zwku,zwkv,    zwkut,zwkvt,   zwkus,zwkvs
  REAL(KIND=8),   DIMENSION (:,:,:), ALLOCATABLE :: ztrpu, ztrpv, ztrput,ztrpvt, ztrpus,ztrpvs
  ! added to write in netcdf
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::  dumlon, dumlat
  REAL(KIND=4), DIMENSION (1)               ::  tim ! time counter
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvar  ! structure of output
  !
  CHARACTER(LEN=256) :: cfilet ,cfileout='section_trp.dat', &
       &                       cfileu, cfilev, csection , &
       &                       cfilvtrp='vtrp.txt', cfilhtrp='htrp.txt', cfilstrp='strp.txt'
  CHARACTER(LEN=256) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc', cdum
  CHARACTER(LEN=256) ,DIMENSION(4)   :: cvarname   !: array of var name for output

  INTEGER    ::  nxtarg
  LOGICAL    :: ltest=.FALSE.
  ! added to write in netcdf
  CHARACTER(LEN=256) :: cfileoutnc 
  ! added to write in netcdf
  LOGICAL :: lwrtcdf=.TRUE.

  ! constants
  REAL(KIND=4)   ::  rau0=1000.,  rcp=4000.

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg < 3  ) THEN
     PRINT *,' Usage : cdftransportiz [-test  u v ]  VTfile gridUfile gridVfile   ''limit of level'' '
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc must be in te current directory'
     PRINT *,' Option -test vt u v is used for testing purposes, with constant flow field'
     PRINT *,' Output on standard output and on an ascii file called section_trp.dat'
     STOP
  ENDIF


  CALL getarg (1, cfilet)
  IF ( cfilet == '-test')  THEN
     ltest = .TRUE.
     CALL getarg (2, cdum)
     READ(cdum,*) udum
     CALL getarg (3, cdum)
     READ(cdum,*) vdum
     CALL getarg (4, cfilet)
     CALL getarg (5, cfileu)
     CALL getarg (6, cfilev)
     nxtarg=6
  ELSE
     CALL getarg (2, cfileu)
     CALL getarg (3, cfilev)
     nxtarg=3
  ENDIF
  nclass = narg -nxtarg + 1

  ALLOCATE (  imeter(nclass -1), ilev0(nclass), ilev1(nclass) )

  DO jk=1, nclass -1
     CALL getarg(nxtarg+jk,cdum)
     READ(cdum,*) imeter(jk)
  END DO

  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  IF(lwrtcdf) THEN

     ALLOCATE ( typvar(nboutput), ipk(nboutput), id_varout(nboutput) )
     ALLOCATE (dumlon(1,1) , dumlat(1,1) )

     dumlon(:,:)=0.
     dumlat(:,:)=0.

     DO jj=1,nboutput
        ipk(jj)=1
     ENDDO

     ! define new variables for output 
     typvar(1)%name='vtrp'
     typvar(1)%units='Sverdrup'
     typvar%missing_value=99999.
     typvar(1)%valid_min= -1000.
     typvar(1)%valid_max= 1000.
     typvar%scale_factor= 1.
     typvar%add_offset= 0.
     typvar%savelog10= 0.
     typvar(1)%long_name='Mass_Transport'
     typvar(1)%short_name='vtrp'
     typvar%online_operation='N/A'
     typvar%axis='T'

     typvar(2)%name='htrp'
     typvar(2)%units='PW'
     typvar(2)%valid_min= -1000.
     typvar(2)%valid_max= 1000.
     typvar(2)%long_name='Heat_Transport'
     typvar(2)%short_name='htrp'

     typvar(3)%name='strp'
     typvar(3)%units='kt/s'
     typvar(3)%valid_min= -1000.
     typvar(3)%valid_max= 1000.
     typvar(3)%long_name='Salt_Transport'
     typvar(3)%short_name='strp'

     typvar(4)%name='lonmin'
     typvar(4)%units='deg'
     typvar(4)%valid_min= -180.
     typvar(4)%valid_max= 180.
     typvar(4)%long_name='minimum_longitude_of_section'
     typvar(4)%short_name='lonmin'

     typvar(5)%name='lonmax'
     typvar(5)%units='deg'
     typvar(5)%valid_min= -180.
     typvar(5)%valid_max= 180.
     typvar(5)%long_name='maximum_longitude_of_section'
     typvar(5)%short_name='lonmax'

     typvar(6)%name='latmin'
     typvar(6)%units='deg'
     typvar(6)%valid_min= -90.
     typvar(6)%valid_max= 90.
     typvar(6)%long_name='minimum_latitude_of_section'
     typvar(6)%short_name='latmin'

     typvar(7)%name='latmax'
     typvar(7)%units='deg'
     typvar(7)%valid_min= -90.
     typvar(7)%valid_max= 90.
     typvar(7)%long_name='maximum_latitude_of_section'
     typvar(7)%short_name='latmax'

     typvar(8)%name='top'
     typvar(8)%units='meters'
     typvar(8)%valid_min= 0.
     typvar(8)%valid_max= 10000.
     typvar(8)%long_name='min_depth_of_the_section'
     typvar(8)%short_name='top'

     typvar(9)%name='bottom'
     typvar(9)%units='meters'
     typvar(9)%valid_min= 0.
     typvar(9)%valid_max= 10000.
     typvar(9)%long_name='max_depth_of_the_section'
     typvar(9)%short_name='bottom'

  ENDIF

  ! Allocate arrays
  ALLOCATE( zu (npiglo,npjglo), zut(npiglo,npjglo), zus(npiglo,npjglo) )
  ALLOCATE( zv (npiglo,npjglo), zvt(npiglo,npjglo), zvs(npiglo,npjglo) )
  !
  ALLOCATE ( zwku (npiglo,npjglo), zwkut(npiglo,npjglo), zwkus(npiglo,npjglo) )
  ALLOCATE ( zwkv (npiglo,npjglo), zwkvt(npiglo,npjglo), zwkvs(npiglo,npjglo) )
  !
  ALLOCATE ( ztrpu (npiglo,npjglo,nclass), ztrpv (npiglo,npjglo,nclass))
  ALLOCATE ( ztrput(npiglo,npjglo,nclass), ztrpvt(npiglo,npjglo,nclass))
  ALLOCATE ( ztrpus(npiglo,npjglo,nclass), ztrpvs(npiglo,npjglo,nclass))
  !
  ALLOCATE ( e1v(npiglo,npjglo),e3v(npiglo,npjglo))
  ALLOCATE ( e2u(npiglo,npjglo),e3u(npiglo,npjglo))
  !
  ALLOCATE ( gphiu(npiglo,npjglo),  gphiv(npiglo,npjglo) )
  ALLOCATE ( glamu(npiglo,npjglo),  glamv(npiglo,npjglo) )
  ALLOCATE ( gdepw(npk) )
  !

  e1v(:,:) = getvar(coordhgr, 'e1v', 1,npiglo,npjglo)
  e2u(:,:) = getvar(coordhgr, 'e2u', 1,npiglo,npjglo)

  glamv(:,:) =  getvar(coordhgr, 'glamv', 1,npiglo,npjglo)
  glamu(:,:) =  getvar(coordhgr, 'glamu', 1,npiglo,npjglo)

  gphiv(:,:) = getvar(coordhgr, 'gphiv', 1,npiglo,npjglo)
  gphiu(:,:) = getvar(coordhgr, 'gphiu', 1,npiglo,npjglo)

  gdepw(:) = getvare3(coordzgr, 'gdepw',npk)

  ! look for nearest level to imeter
  ik = 1

  ilev0(1)      = 1
  ilev1(nclass) = npk-1

  DO jk = 1, nclass -1
     DO WHILE ( gdepw(ik)  < imeter(jk) )
        ik = ik +1
     END DO

     rd1= ABS(gdepw(ik-1) - imeter(jk) )
     rd2= ABS(gdepw(ik) - imeter(jk) )
     IF ( rd2 < rd1 ) THEN
        ilev1(jk) = ik -1  ! t-levels
        ilev0(jk+1) = ik
     ELSE 
        ilev1(jk) = ik -2  ! t-levels
        ilev0(jk+1) = ik -1
     END IF
  END DO
  PRINT *, 'Limits :  '
  DO jk = 1, nclass
     PRINT *, ilev0(jk),ilev1(jk), gdepw(ilev0(jk)), gdepw(ilev1(jk)+1)
  END DO

  !! compute the transport
  ztrpu (:,:,:)= 0
  ztrpv (:,:,:)= 0

  ztrput(:,:,:)= 0
  ztrpvt(:,:,:)= 0

  ztrpus(:,:,:)= 0
  ztrpvs(:,:,:)= 0
  DO jclass = 1, nclass
     DO jk = ilev0(jclass),ilev1(jclass)
        PRINT *,'level ',jk
        ! Get velocities, temperature and salinity fluxes at jk
        IF ( ltest ) THEN
           zu (:,:)= udum
           zv (:,:)= vdum
           zut(:,:)= udum
           zvt(:,:)= vdum
           zus(:,:)= udum
           zvs(:,:)= vdum
        ELSE
           zu (:,:)= getvar(cfileu, 'vozocrtx',  jk ,npiglo,npjglo)
           zv (:,:)= getvar(cfilev, 'vomecrty',  jk ,npiglo,npjglo)
           zut(:,:)= getvar(cfilet, 'vozout',  jk ,npiglo,npjglo)
           zvt(:,:)= getvar(cfilet, 'vomevt',  jk ,npiglo,npjglo)
           zus(:,:)= getvar(cfilet, 'vozous',  jk ,npiglo,npjglo)
           zvs(:,:)= getvar(cfilet, 'vomevs',  jk ,npiglo,npjglo)
        ENDIF

        ! get e3u, e3v  at level jk
        e3v(:,:) = getvar(coordzgr, 'e3v_ps', jk,npiglo,npjglo, ldiom=.TRUE.)
        e3u(:,:) = getvar(coordzgr, 'e3u_ps', jk,npiglo,npjglo, ldiom=.TRUE.)

        zwku (:,:) = zu (:,:)*e2u(:,:)*e3u(:,:)
        zwkv (:,:) = zv (:,:)*e1v(:,:)*e3v(:,:)
        zwkut(:,:) = zut(:,:)*e2u(:,:)*e3u(:,:)
        zwkvt(:,:) = zvt(:,:)*e1v(:,:)*e3v(:,:)
        zwkus(:,:) = zus(:,:)*e2u(:,:)*e3u(:,:)
        zwkvs(:,:) = zvs(:,:)*e1v(:,:)*e3v(:,:)

        ! integrates vertically 
        ztrpu (:,:,jclass) = ztrpu (:,:,jclass) + zwku (:,:)
        ztrpv (:,:,jclass) = ztrpv (:,:,jclass) + zwkv (:,:)
        ztrput(:,:,jclass) = ztrput(:,:,jclass) + zwkut(:,:) * rau0*rcp
        ztrpvt(:,:,jclass) = ztrpvt(:,:,jclass) + zwkvt(:,:) * rau0*rcp
        ztrpus(:,:,jclass) = ztrpus(:,:,jclass) + zwkus(:,:)  
        ztrpvs(:,:,jclass) = ztrpvs(:,:,jclass) + zwkvs(:,:)  

     END DO  ! loop to next level
  END DO    ! next class

  OPEN(numout,FILE=cfileout)
  ! also dump the results on txt files without any comments, some users  like it !
  OPEN(numvtrp,FILE=cfilvtrp)
  OPEN(numhtrp,FILE=cfilhtrp)
  OPEN(numstrp,FILE=cfilstrp)
  DO 
     PRINT *, ' Give name of section '
     READ(*,'(a)') csection
     IF (TRIM(csection) == 'EOF' ) THEN ; CLOSE(numout) ; CLOSE(numvtrp) ; CLOSE(numhtrp) ; CLOSE(numstrp) ; ENDIF
        IF (TRIM(csection) == 'EOF' ) EXIT
        PRINT *, ' Give imin, imax, jmin, jmax '
        READ(*,*) imin, imax, jmin, jmax
        !! Find the broken line between P1 (imin,jmin) and P2 (imax, jmax)
        !! ---------------------------------------------------------------
        ! ... Initialization
        i0=imin; j0=jmin; i1=imax;  j1=jmax
        rxi1=i1;  ryj1=j1; rxi0=i0; ryj0=j0

        ! .. Compute equation:  ryj = aj rxi + bj
        IF ( (rxi1 -rxi0) /=  0 ) THEN
           aj = (ryj1 - ryj0 ) / (rxi1 -rxi0)
           bj = ryj0 - aj * rxi0
        ELSE
           aj=10000.
           bj=0.
        END IF

        ! .. Compute equation:  rxi = ai ryj + bi
        IF ( (ryj1 -ryj0) /=  0 ) THEN
           ai = (rxi1 - rxi0 ) / ( ryj1 -ryj0 )
           bi = rxi0 - ai * ryj0
        ELSE
           ai=10000.
           bi=0.
        END IF

        ! ..  Compute the integer pathway:
        n=0
        ! .. Chose the strait line with the smallest slope
        IF (ABS(aj) <=  1 ) THEN
           ! ... Here, the best line is y(x)
           ! ... If i1 < i0 swap points and remember it has been swapped
           IF (i1 <  i0 ) THEN
              i  = i0 ; j  = j0
              i0 = i1 ; j0 = j1
              i1 = i  ; j1 = j
           END IF

           IF ( j1 >= j0 ) THEN
              ist = 1     ; jst = 1
              norm_u =  1 ;  norm_v = -1
           ELSE
              ist = 1     ; jst = 0
              norm_u = -1 ; norm_v = -1
           END IF

           ! ... compute the nearest j point on the line crossing at i
           DO i=i0,i1
              n=n+1
              IF (n > jpseg) STOP 'n > jpseg !'
              j=NINT(aj*i + bj )
              yypt(n) = CMPLX(i,j)
           END DO
        ELSE
           ! ... Here, the best line is x(y)
           ! ... If j1 < j0 swap points and remember it has been swapped
           IF (j1 <  j0 ) THEN
              i  = i0 ; j  = j0
              i0 = i1 ; j0 = j1
              i1 = i  ; j1 = j
           END IF
           IF ( i1 >=  i0 ) THEN
              ist = 1    ;  jst = 1
              norm_u = 1 ;  norm_v = -1
           ELSE
              ist = 0
              jst = 1
              norm_u = 1
              norm_v = 1
           END IF

           ! ... compute the nearest i point on the line crossing at j
           DO j=j0,j1
              n=n+1
              IF (n > jpseg) STOP 'n>jpseg !'
              i=NINT(ai*j + bi)
              yypt(n) = CMPLX(i,j)
           END DO
        END IF

        !!
        !! Look for intermediate points to be added.
        !  ..  The final positions are saved in rxx,ryy
        rxx(1)=REAL(yypt(1))
        ryy(1)=IMAG(yypt(1))
        nn=1

        DO k=2,n
           ! .. distance between 2 neighbour points
           d=ABS(yypt(k)-yypt(k-1))
           ! .. intermediate points required if d > 1
           IF ( d > 1 ) THEN
              CALL interm_pt(yypt,k,ai,bi,aj,bj,yypti)
              nn=nn+1
              IF (nn > jpseg) STOP 'nn>jpseg !'
              rxx(nn)=REAL(yypti)
              ryy(nn)=IMAG(yypti)
           END IF
           nn=nn+1
           IF (nn > jpseg) STOP 'nn>jpseg !'
           rxx(nn)=REAL(yypt(k))
           ryy(nn)=IMAG(yypt(k))
        END DO

        ! Now extract the transport through a section 
        ! ... Check whether we need a u velocity or a v velocity
        !   Think that the points are f-points and delimit either a U segment
        !   or a V segment (ist and jst are set in order to look for the correct
        !   velocity point on the C-grid
        PRINT *, TRIM(csection)
        PRINT *, 'IMIN IMAX JMIN JMAX', imin, imax, jmin, jmax
        WRITE(numout,*)'% Transport along a section by levels' ,TRIM(csection)
        WRITE(numout,*) '% nada IMIN IMAX JMIN JMAX'
        DO jclass=1,nclass
           voltrpsum = 0.
           heatrpsum = 0.
           saltrpsum = 0.

           DO jseg = 1, nn-1
              i0=rxx(jseg)
              j0=ryy(jseg)
              IF ( rxx(jseg) ==  rxx(jseg+1) ) THEN
                 gla(jseg)=glamu(i0,j0+jst)   ; gphi(jseg)=gphiu(i0,j0+jst)
                 voltrp(jseg)= ztrpu (i0,j0+jst,jclass)*norm_u
                 heatrp(jseg)= ztrput(i0,j0+jst,jclass)*norm_u
                 saltrp(jseg)= ztrpus(i0,j0+jst,jclass)*norm_u
              ELSE IF ( ryy(jseg) == ryy(jseg+1) ) THEN
                 gla(jseg)=glamv(i0+ist,j0)  ;  gphi(jseg)=gphiv(i0+ist,j0)
                 voltrp(jseg)=ztrpv (i0+ist,j0,jclass)*norm_v
                 heatrp(jseg)=ztrpvt(i0+ist,j0,jclass)*norm_v
                 saltrp(jseg)=ztrpvs(i0+ist,j0,jclass)*norm_v
              ELSE
                 PRINT *,' ERROR :',  rxx(jseg),ryy(jseg),rxx(jseg+1),ryy(jseg+1)
              END IF
              voltrpsum = voltrpsum+voltrp(jseg)
              heatrpsum = heatrpsum+heatrp(jseg)
              saltrpsum = saltrpsum+saltrp(jseg)
           END DO   ! next segment
           IF (jclass == 1 ) PRINT *, 'FROM (LON LAT): ', gla(1),gphi(1),' TO (LON LAT)', gla(nn-1), gphi(nn-1)
           PRINT *, gdepw(ilev0(jclass)), gdepw(ilev1(jclass)+1)
           PRINT *, ' Mass transport : ', voltrpsum/1.e6,' SV'
           PRINT *, ' Heat transport : ', heatrpsum/1.e15,' PW'
           PRINT *, ' Salt transport : ', saltrpsum/1.e6,' kT/s'
           IF (jclass == 1 ) THEN 
              WRITE(numout,*)  '% nada LONmin LATmin LONmax LATmax'
              WRITE(numout,*)  '% Top(m)  Bottom(m)  MassTrans(Sv) HeatTrans(PW) SaltTrans(kt/s)'
              WRITE(numout,*) 0 ,imin, imax, jmin, jmax
              WRITE(numout,9003) 0. ,gla(1),gphi(1), gla(nn-1), gphi(nn-1)
           ENDIF
           WRITE(numout,9002) gdepw(ilev0(jclass)), gdepw(ilev1(jclass)+1), voltrpsum/1.e6, heatrpsum/1.e15, saltrpsum/1.e6
           WRITE(numvtrp,'(e12.6)') voltrpsum
           WRITE(numhtrp,'(e12.6)') heatrpsum
           WRITE(numstrp,'(e12.6)') saltrpsum

           IF(lwrtcdf) THEN

              ! create output fileset
              cfileoutnc=TRIM(csection)//'_transports.nc'
              ncout =create(cfileoutnc,'none',kx,ky,kz,cdep='depthw')
              ierr= createvar(ncout,typvar,nboutput,ipk,id_varout )
              ierr= putheadervar(ncout, cfilet,kx, &
                   ky,kz,pnavlon=dumlon,pnavlat=dumlat,pdep=gdepw)
              tim=getvar1d(cfilet,'time_counter',1)
              ierr=putvar1d(ncout,tim,1,'T')

              ! netcdf output 
              ierr = putvar0d(ncout,id_varout(1), REAL(voltrpsum/1.e6) )
              ierr = putvar0d(ncout,id_varout(2), REAL(heatrpsum/1.e15) )
              ierr = putvar0d(ncout,id_varout(3), REAL(saltrpsum/1.e6) )
              ierr = putvar0d(ncout,id_varout(4), REAL(gla(1)) )
              ierr = putvar0d(ncout,id_varout(5), REAL(gla(nn-1)) )
              ierr = putvar0d(ncout,id_varout(6), REAL(gphi(1)) )
              ierr = putvar0d(ncout,id_varout(7), REAL(gphi(nn-1)) )
              ierr = putvar0d(ncout,id_varout(8), REAL(gdepw(ilev0(jclass))) )
              ierr = putvar0d(ncout,id_varout(9), REAL(gdepw(ilev1(jclass)+1)) )

              ierr = closeout(ncout)

           ENDIF


        END DO ! next class

     END DO ! infinite loop : gets out when input is EOF 

9000 FORMAT(I4,6(f9.3,f8.4))
9001 FORMAT(I4,6(f9.2,f9.3))
9002 FORMAT(f9.0,f9.0,f9.2,f9.2,f9.2)
9003 FORMAT(f9.2,f9.2,f9.2,f9.2,f9.2)

   CONTAINS
     SUBROUTINE interm_pt (ydpt,k,pai,pbi,paj,pbj,ydpti)
       !! -----------------------------------------------------
       !!           SUBROUTINE INTERM_PT
       !!           ********************
       !!
       !!   PURPOSE:
       !!   --------
       !!     Find the best intermediate points on a pathway.
       !!
       !!    ARGUMENTS:
       !!    ----------
       !!      ydpt : complex vector of the positions of the nearest points
       !!         k : current working index
       !!       pai ,pbi : slope and original ordinate of x(y)
       !!       paj ,pbj : slope and original ordinate of y(x)
       !!      ydpti : Complex holding the position of intermediate point
       !!
       !!    AUTHOR:
       !!    -------
       !!      19/07/1999 : Jean-Marc MOLINES
       !!      14/01/2005 : J M M in F90
       !!
       !!--------------------------------------------------------------
       !!
       !! 0. Declarations:
       !! ----------------
       IMPLICIT NONE
       COMPLEX, INTENT(in) :: ydpt(*)
       COMPLEX, INTENT(out) :: ydpti
       REAL(KIND=4), INTENT(IN) ::  pai,pbi,paj,pbj
       INTEGER ,INTENT(in) :: k
       ! ... local
       COMPLEX :: ylptmp1, ylptmp2
       REAL(KIND=4) ::  za0,zb0,za1,zb1,zd1,zd2
       REAL(KIND=4) ::  zxm,zym
       REAL(KIND=4) ::  zxp,zyp
       !!
       !! 1. Compute intermediate points
       !! ------------------------------
       !
       ! ... Determines whether we use y(x) or x(y):
       IF (ABS(paj) <=  1) THEN
          ! ..... y(x)
          ! ... possible intermediate points:
          ylptmp1=ydpt(k-1)+(1.,0.)
          ylptmp2=ydpt(k-1)+CMPLX(0.,SIGN(1.,paj))
          !
          ! ... M is the candidate point:
          zxm=REAL(ylptmp1)
          zym=IMAG(ylptmp1)
          za0=paj
          zb0=pbj
          !
          za1=-1./za0
          zb1=zym - za1*zxm
          ! ... P is the projection of M on the strait line
          zxp=-(zb1-zb0)/(za1-za0)
          zyp=za0*zxp + zb0
          ! ... zd1 is the distance MP
          zd1=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
          !
          ! ... M is the candidate point:
          zxm=REAL(ylptmp2)
          zym=IMAG(ylptmp2)
          za1=-1./za0
          zb1=zym - za1*zxm
          ! ... P is the projection of M on the strait line
          zxp=-(zb1-zb0)/(za1-za0)
          zyp=za0*zxp + zb0
          ! ... zd2 is the distance MP
          zd2=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
          ! ... chose the smallest (zd1,zd2)
          IF (zd2 <=  zd1) THEN
             ydpti=ylptmp2
          ELSE
             ydpti=ylptmp1
          END IF
          !
       ELSE
          !
          ! ... x(y)
          ylptmp1=ydpt(k-1)+CMPLX(SIGN(1.,pai),0.)
          ylptmp2=ydpt(k-1)+(0.,1.)
          zxm=REAL(ylptmp1)
          zym=IMAG(ylptmp1)
          za0=pai
          zb0=pbi
          !
          za1=-1./za0
          zb1=zxm - za1*zym
          zyp=-(zb1-zb0)/(za1-za0)
          zxp=za0*zyp + zb0
          zd1=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
          !
          zxm=REAL(ylptmp2)
          zym=IMAG(ylptmp2)
          za1=-1./za0
          zb1=zxm - za1*zym
          zyp=-(zb1-zb0)/(za1-za0)
          zxp=za0*zyp + zb0
          zd2=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
          IF (zd2 <=  zd1) THEN
             ydpti=ylptmp2
          ELSE
             ydpti=ylptmp1
          END IF
       END IF
     END SUBROUTINE interm_pt

   END PROGRAM cdftransportiz
