PROGRAM cdfrhoproj
  !! --------------------------------------------------------------
  !!               ***   PROGRAM RHO_VERT_INT ***
  !!  ** Purpose:  This program is used to project any scalar on the A grid
  !!               onto given isopycnic surfaces.
  !!
  !!  ** Method:   Linear interpolation is used on the vertical to define
  !!             the depth of the given isopycn and linear interpolation
  !!             is also performed on the scalar to determine its value at
  !!             this depth.
  !!
  !!  ** Usage :
  !!         cdfrhoproj [-s0 sig0] 'rho file' 'scalar file (*)'
  !!
  !!  * history:
  !!        Original : J.M. Molines for SPEM in Dynamo (1996)
  !!        Modif    : J-O. Beismann for OPA (1999)
  !!        Modif    : J.M. Molines for normalization Clipper (March 2000)
  !!                 : J.M. Molines in cdftools, f90 dor DRAKKAR (Nov. 2005)
  !! ---------------------------------------------------------------------
  !! * Used modules
  USE cdfio

  !! * Local declaration
  IMPLICIT NONE

  INTEGER :: npiglo, npjglo, npk, npkk ,npt  ,nvars
  INTEGER :: narg, iargc
  INTEGER :: ji,jj,jk,jkk,jfich,k0, jvar
  INTEGER :: istartarg = 1
  INTEGER :: ncout, ierr
  INTEGER, DIMENSION(2) :: ipk, id_varout !: for output variables
  !
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE ::  v3d, alpha
  REAL(KIND=4), DIMENSION(:,:)  , ALLOCATABLE ::  v2dint, zint, v2d
  REAL(KIND=4), DIMENSION(:)    , ALLOCATABLE ::  zi, time_tag, h1d
  REAL(KIND=4)                                ::  x1z,y1z,dxz,dyz,P1,P2
  REAL(KIND=4)                                ::  spval=999999.
  REAL(KIND=4)                                ::  spvalz=0.

  CHARACTER(LEN=80) ::  cline, cfilZI, cfildata, cfilRHOMOD, cvar, cfilout, ctype='T'
  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: czvar     !: temporary arry for variable name in file
 
  TYPE(variable), DIMENSION(2)  :: typvar      !: structure for attributes
  TYPE(variable), DIMENSION(:), ALLOCATABLE  :: typzvar      !: structure for attributes
  !
  LOGICAL :: lsingle=.false.

  !! * Read command line
  narg=iargc()
  IF (narg < 3 ) THEN
         PRINT *, &
     &' >>>> usage: cdfrhoproj [-s0 sig0 ] <cvar>  <rhofile> <file*> [ T | U | V | W ]'
         PRINT *,'   Interpolated files will be file.nc.interp'
         PRINT *,'   Isopycnal value are read on a text file ''rho_lev'' '
         PRINT *,'   unless the option  -s0 is specified with one particular value'
         PRINT *,'   cvar specify the name of the cdf variable to interpolate '
         PRINT *,'   Model density are taken on file ''rhofile'' '
         PRINT *,'   File is the netcdf file holding cvar '
         PRINT *,'   Last argument is optional (T by default) and indicate the '
         PRINT *,'   C-Grid point corresponding to file.'
         PRINT *,'   cvar will be interpolated on T point previous projection on isopycnals'
         STOP
  ENDIF
  !   seek a -s0 option
  CALL getarg(1,cline)
  IF (cline == '-s0' ) THEN
     npkk = 1
     lsingle=.true.
     istartarg = 3
     CALL getarg(2,cline)
     ALLOCATE (zi(npkk) )
     READ(cline,*) zi(1)
  END IF
  ! read ZI if not single 
  IF ( .NOT.  lsingle ) THEN
     cfilZI='rho_lev'
     OPEN(10,file=cfilZI)
     READ(10,*) npkk
     ALLOCATE (zi(npkk) )
     DO jkk=1,npkk
        READ(10,*) zi(jkk)
        PRINT *,zi(jkk)
     END DO
     CLOSE(10)
  ENDIF
  !   Seek for ctype (last argument either T U V or W )
  CALL getarg(narg,ctype)
  SELECT CASE ( ctype)
  CASE ( 'T','t','U','u','V','v','W','w'  )
     narg = narg -1  ! last argument is not a file name
  CASE DEFAULT
     ctype='T'
  END SELECT
   
  ! Read variable name
  CALL getarg(istartarg,cvar)
  ! Read Rho file
  CALL getarg(istartarg+1,cfilRHOMOD)
  npiglo=getdim(cfilRHOMOD,'x')
  npjglo=getdim(cfilRHOMOD,'y')
  npk   =getdim(cfilRHOMOD,'depth')
  npt   =getdim(cfilRHOMOD,'time')

  CALL getarg(istartarg+2, cfildata)
  nvars=getnvar(cfildata)
  ALLOCATE(czvar(nvars), typzvar(nvars))

  czvar(:)=getvarname(cfildata,nvars,typzvar)

  ALLOCATE( v3d(npiglo,npjglo,npk), alpha(npiglo, npjglo, npkk) )
  ALLOCATE( v2dint(npiglo, npjglo), v2d(npiglo,npjglo), zint(npiglo,npjglo)  )
  ALLOCATE( time_tag(npt), h1d(npk) )

  time_tag(:)=getvar1d(cfilRHOMOD,'time_counter', npt)
  h1d(:)=getvar1d(cfilRHOMOD,'deptht',npk)
   
  DO jk=1,npk
     v3d(:,:,jk) = getvar(cfilRHOMOD,'vosigma0',jk,npiglo,npjglo)
  END DO

  !! ** Compute interpolation coefficients as well as the level used
  !!    to interpolate between
  DO ji=1,npiglo
     DO jj = 1, npjglo
        jk = 1
        DO jkk=1,npkk
        !  Assume that rho (z) is increasing downward (no inversion)
        !     Caution with sigma0 at great depth !
           DO WHILE (zi(jkk) >=  v3d(ji,jj,jk) .AND. jk <= npk &
     &                .AND. v3d(ji,jj,jk) /=  spvalz )
              jk=jk+1
           END DO
           jk=jk-1
           k0=jk
           IF (jk .EQ. 0) THEN
              jk=1
              alpha(ji,jj,jkk) = 0.
           ELSE IF (v3d(ji,jj,jk+1) .EQ. spvalz ) THEN
              k0=0
              alpha(ji,jj,jkk) = 0.
           ELSE 
           ! ... alpha is always in [0,1]. Adding k0 ( >=1 ) for saving space for k0
              alpha(ji,jj,jkk)= &
     &               (zi(jkk)-v3d(ji,jj,jk))/(v3d(ji,jj,jk+1)-v3d(ji,jj,jk)) +k0
           ENDIF
        END DO
     END DO
  END DO

  !! ** Loop on the scalar files to project on choosen isopycnics surfaces
  DO jfich=istartarg+2,narg

     CALL getarg(jfich,cfildata)
     PRINT *,'working with ', TRIM(cfildata)

     IF (npt /= 1 ) THEN
        PRINT *,' This program has to be modified for multiple'
        PRINT *,' time frames.'
        STOP ' Error : npt # 1'
     ENDIF
 
     DO jk=1,npk
        v2d(:,:) = getvar(cfildata,cvar,jk,npiglo,npjglo)
        SELECT CASE ( ctype )
        CASE ('T', 't' )
           v3d(:,:,jk) = v2d(:,:)
        CASE ('U','u' )
           DO ji=2,npiglo
              DO jj=1, npjglo
                 v3d(ji,jj,jk)=0.5*( v2d(ji,jj) + v2d(ji-1,jj) )  ! put variable on T point
              END DO
           END DO
        CASE ('V','v' )
           DO jj=2,npjglo
              DO ji=1, npiglo
                 v3d(ji,jj,jk)=0.5*( v2d(ji,jj) + v2d(ji,jj-1) )  ! put variable on T point
              END DO
           END DO
         CASE('W','w' )
          STOP 'Case W not done yet :( '
         END SELECT
     END DO

     ! ... open output file and write header
     ipk(:)=npkk
     DO jvar=1,nvars
       IF ( cvar == typzvar(jvar)%name ) THEN 
          typvar(1)=typzvar(jvar)
          EXIT
       ENDIF
     END DO
     typvar(1)%long_name=TRIM(typvar(1)%long_name)//' on iso sigma'
     typvar(1)%axis='TRYX'

     typvar(2)%name= 'vodepiso'
     typvar(2)%units='m'
     typvar(2)%missing_value=999999.
     typvar(2)%valid_min= 0.
     typvar(2)%valid_max= 7000.
     typvar(2)%long_name='Depth_of_Isopycnals'
     typvar(2)%short_name='vodepiso'
     typvar(2)%online_operation='N/A'
     typvar(2)%axis='TRYX'


     cfilout=TRIM(cfildata)//'.interp'
 
     ncout = create(cfilout,cfilRHOMOD ,npiglo,npjglo,npkk)
     ierr = createvar(ncout, typvar,2,ipk, id_varout )
     ierr = putheadervar(ncout , cfilRHOMOD, npiglo, npjglo, npkk,pdep=zi)
 
     DO jkk=1,npkk
        DO ji=1,npiglo
           DO jj=1,npjglo
             ! k0 is retrieved from alpha, taking the integer part.
             ! The remnant is alpha. 
             k0=INT(alpha(ji,jj,jkk))
             alpha(ji,jj,jkk) =  alpha(ji,jj,jkk) - k0
             IF (k0 /= 0) THEN
              P1=v3d(ji,jj,k0)
              P2=v3d(ji,jj,k0+1)
                IF (P1 /= spvalz .AND. P2 /= spvalz) THEN
                   v2dint(ji,jj)=alpha(ji,jj,jkk)*P2  &
   &                         +(1-alpha(ji,jj,jkk))*P1
                    zint (ji,jj)=alpha(ji,jj,jkk)*h1d(k0+1) &
   &                         +(1-alpha(ji,jj,jkk))*h1d(k0)
                ELSE 
                   v2dint(ji,jj)=spval
                   zint  (ji,jj)=spval
               ENDIF
             ELSE 
               v2dint(ji,jj)=spval
               zint  (ji,jj)=spval
             ENDIF
           END DO
        END DO
        ierr = putvar(ncout,id_varout(1), v2dint,jkk,npiglo,npjglo)
        ierr = putvar(ncout,id_varout(2), zint  ,jkk,npiglo,npjglo)
     END DO
     ierr = putvar1d(ncout,time_tag,1,'T')
     ierr = closeout(ncout)
  END DO  ! loop on scalar files
        PRINT *,'Projection on isopycns completed successfully'
END  PROGRAM cdfrhoproj
