PROGRAM cdfsigintegr
  !! --------------------------------------------------------------
  !!               ***   PROGRAM cdfsigintegr  ***
  !!  ** Purpose:  This program is used to integrate quantities between isopycnals
  !!
  !!  ** Method:   Linear interpolation is used on the vertical to define
  !!             the depth of the given isopycn.
  !!             Then, the integral is performed from the top of the ocean down to the given
  !!             isopycnal. Finaly, by making the difference between 2 isopycnals we obtain
  !!             the required quantity.
  !!
  !!  ** Usage :
  !!         cdfsigintegr 'rho file' 'scalar file (*)'
  !!
  !!  * history:
  !!        Original : J.M. Molines December 2007 (From cdfrhoproj )
  !! ---------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------

  !! * Used modules
  USE cdfio

  !! * Local declaration
  IMPLICIT NONE

  INTEGER :: npiglo, npjglo, npk, npkk ,npt  ,nvars
  INTEGER :: narg, iargc
  INTEGER :: ji,jj,jk,jkk,jfich,k0, jvar
  INTEGER :: ncout, ierr
  INTEGER, DIMENSION(3) :: ipk, id_varout !: for output variables
  !
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE ::  v3d, alpha
  REAL(KIND=4), DIMENSION(:,:)  , ALLOCATABLE ::  v2d, e3, tmask, dum
  REAL(KIND=4), DIMENSION(:)    , ALLOCATABLE ::  zi, time_tag, h1d, gdepw
  REAL(KIND=4)                                ::  spval=999999.
  REAL(KIND=4)                                ::  spvalz=0.
  REAL(KIND=4), DIMENSION(:,:,:)  , ALLOCATABLE ::  zint   
  REAL(KIND=8), DIMENSION(:,:,:)  , ALLOCATABLE ::  v2dint   !: double precision for integration

  CHARACTER(LEN=256) ::  cfilZI, cfildata, cfilRHOMOD, cvar, cfilout, ctype='T'
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: czvar     !: temporary arry for variable name in file
  CHARACTER(LEN=256) :: coordzgr='mesh_zgr.nc'          !: coordinates files

 
  TYPE(variable), DIMENSION(3)  :: typvar      !: structure for attributes
  TYPE(variable), DIMENSION(:), ALLOCATABLE  :: typzvar      !: structure for attributes
  !

  !! * Read command line
  narg=iargc()
  IF (narg < 3 ) THEN
         PRINT *, &
     &' >>>> usage: cdfsigintegr  <cvar>  <rhofile> <file*> [ T | U | V | W ]'
         PRINT *,'   Interpolated files will be file.nc.interp'
         PRINT *,'   Isopycnal value are read on a text file ''rho_lev'' (minimum 2)'
         PRINT *,'   cvar specify the name of the cdf variable to interpolate '
         PRINT *,'   Model density are taken on file ''rhofile'' '
         PRINT *,'   File is the netcdf file holding cvar '
         PRINT *,'   Last argument is optional (T by default) and indicate the '
         PRINT *,'   C-Grid point corresponding to file.'
         PRINT *,'   cvar will be interpolated on T point previous projection on isopycnals'
         STOP
  ENDIF
     cfilZI='rho_lev'
     OPEN(10,file=cfilZI)
     READ(10,*) npkk
     ALLOCATE (zi(npkk) )
     DO jkk=1,npkk
        READ(10,*) zi(jkk)
        PRINT *,zi(jkk)
     END DO
     CLOSE(10)

  !   Seek for ctype (last argument either T U V or W )
  CALL getarg(narg,ctype)
  SELECT CASE ( ctype)
  CASE ( 'T','t','U','u','V','v','W','w'  )
     narg = narg -1  ! last argument is not a file name
  CASE DEFAULT
     ctype='T'
  END SELECT
   
  ! Read variable name
  CALL getarg(1,cvar)
  ! Read Rho file
  CALL getarg(2,cfilRHOMOD)
  npiglo=getdim(cfilRHOMOD,'x')
  npjglo=getdim(cfilRHOMOD,'y')
  npk   =getdim(cfilRHOMOD,'depth')
  npt   =getdim(cfilRHOMOD,'time')

  spvalz=getspval(cfilRHOMOD,'vosigma0')

  CALL getarg(3, cfildata)
  nvars=getnvar(cfildata)
  ALLOCATE(czvar(nvars), typzvar(nvars))

  czvar(:)=getvarname(cfildata,nvars,typzvar)

  ALLOCATE( v3d(npiglo,npjglo,npk), alpha(npiglo, npjglo, npkk) ,e3(npiglo,npjglo) )
  ALLOCATE( v2dint(npiglo, npjglo,2), v2d(npiglo,npjglo), zint(npiglo,npjglo,2)  )
  ALLOCATE( time_tag(npt), h1d(npk) ,gdepw(npk) ,tmask(npiglo,npjglo), dum(npiglo,npjglo) )

  gdepw(:) = getvare3(coordzgr,'gdepw', npk)

  time_tag(:)=getvar1d(cfilRHOMOD,'time_counter', npt)
  h1d(:)=getvar1d(cfilRHOMOD,'deptht',npk)
   
  ! Note, if working with vertical slabs, one may avoid 3D array, but may be slow ...
  tmask=1.
  DO jk=1,npk
     v3d(:,:,jk) = getvar(cfilRHOMOD,'vosigma0',jk,npiglo,npjglo)
     IF ( jk == 1 ) THEN
       WHERE (v3d(:,:,jk) == spvalz ) tmask=0.
     ENDIF
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
  DO jfich=3,narg

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
     ipk(1)=npkk-1 
     ipk(2)=npkk-1 
     ipk(3)=npkk 
     DO jvar=1,nvars
       IF ( cvar == typzvar(jvar)%name ) THEN 
          typvar(1)=typzvar(jvar)
          EXIT
       ENDIF
     END DO
     typvar(1)%name='inv'//TRIM(typvar(1)%name)
     typvar(1)%long_name=TRIM(typvar(1)%long_name)//' integrated on sigma bin'
     typvar(1)%units=TRIM(typvar(1)%units)//'.m'
     typvar(1)%missing_value=spval
     typvar(1)%axis='TRYX'

     typvar(2)%name= 'isothick'
     typvar(2)%units='m'
     typvar(2)%missing_value=spval
     typvar(2)%valid_min= 0.
     typvar(2)%valid_max= 7000.
     typvar(2)%long_name='Thickness_of_Isopycnals'
     typvar(2)%short_name='isothick'
     typvar(2)%online_operation='N/A'
     typvar(2)%axis='TRYX'

     typvar(3)%name= 'vodepiso'
     typvar(3)%units='m'
     typvar(3)%missing_value=spval
     typvar(3)%valid_min= 0.
     typvar(3)%valid_max= 7000.
     typvar(3)%long_name='Depth_of_Isopycnals'
     typvar(3)%short_name='vodepiso'
     typvar(3)%online_operation='N/A'
     typvar(3)%axis='TRYX'


     cfilout=TRIM(cfildata)//'.integr'
 
     ncout = create(cfilout,cfilRHOMOD ,npiglo,npjglo,npkk)
     ierr = createvar(ncout, typvar,3,ipk, id_varout )
     ierr = putheadervar(ncout , cfilRHOMOD, npiglo, npjglo, npkk,pdep=zi)

     ! Compute integral from surface to isopycnal
     DO jkk=1,npkk
        ! determine isopycnal surface
        DO ji=1,npiglo
           DO jj=1,npjglo
             ! k0 is retrieved from alpha, taking the integer part.
             k0=INT(alpha(ji,jj,jkk)) ; alpha(ji,jj,jkk) =  alpha(ji,jj,jkk) - k0
             IF (k0 /= 0) THEN
                    zint (ji,jj,1)=alpha(ji,jj,jkk)*h1d(k0+1) &
   &                         +(1-alpha(ji,jj,jkk))*h1d(k0)
             ELSE 
               zint  (ji,jj,1)=0.  !spval ! 
             ENDIF
           END DO
        END DO
        ! integrate from jk=1 to zint
        v2dint(:,:,1) = 0.d0

        DO jk=1,npk-1
          ! get metrixs at level jk
          e3(:,:)=getvar(coordzgr,'e3t_ps',jk,npiglo,npjglo,ldiom=.true.)
          DO ji=1,npiglo
            DO jj=1,npjglo
              IF ( gdepw(jk)+e3(ji,jj) < zint(ji,jj,1) ) THEN  ! full cell
               v2dint(ji,jj,1)=v2dint(ji,jj,1) + e3(ji,jj)* v3d(ji,jj,jk)
              ELSE IF (( zint(ji,jj,1) <= gdepw(jk)+e3(ji,jj) ) .AND. (zint(ji,jj,1) > gdepw(jk)) ) THEN
               v2dint(ji,jj,1)=v2dint(ji,jj,1)+ (zint(ji,jj,1) - gdepw(jk) )* v3d(ji,jj,jk)
              ELSE   ! below the isopycnal 
               ! do nothing for this i j point
              ENDIF
            END DO
          END DO
        END DO   ! end on vertical integral for isopynal jkk

        dum=zint(:,:,1)

        WHERE (tmask == 0. ) dum=spval
        ierr = putvar(ncout,id_varout(3), dum  ,jkk,npiglo,npjglo)

        IF (jkk > 1  ) THEN  ! compute the difference ie the inventory in the layer between 2 isopycnals
           dum=v2dint(:,:,1) - v2dint(:,:,2) ; WHERE ((tmask == 0.)  .OR. (dum < 0 ) ) dum = spval
           ierr = putvar(ncout,id_varout(1), dum,jkk-1,npiglo,npjglo)

           dum=zint  (:,:,1) - zint  (:,:,2) ; WHERE ((tmask == 0.)  .OR. (dum < 0 ) ) dum = spval
           ierr = putvar(ncout,id_varout(2), dum,jkk-1,npiglo,npjglo)
        ENDIF
         v2dint(:,:,2)=v2dint(:,:,1)
         zint  (:,:,2)=zint  (:,:,1)
          
     END DO
     ierr = putvar1d(ncout,time_tag,1,'T')
     ierr = closeout(ncout)
  END DO  ! loop on scalar files
        PRINT *,' integral between isopycnals completed successfully'
END  PROGRAM cdfsigintegr
