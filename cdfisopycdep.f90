PROGRAM cdfisopycdep
  !! --------------------------------------------------------------
  !!               ***   PROGRAM cdfisopycdep  ***
  !!  ** Purpose:  This program is used to determine the depth of isopycnal
  !!
  !!  ** Method:   Linear interpolation is used on the vertical to define
  !!             the depth of the given isopycn.
  !!
  !!  ** Usage :
  !!         cdfisopycdep [-s sigma] 'rho file' cdfsigmavar
  !!
  !!  * history:
  !!        Original : J.M. Molines for SPEM in Dynamo (1996)
  !!        Modif    : J-O. Beismann for OPA (1999)
  !!        Modif    : J.M. Molines for normalization Clipper (March 2000)
  !!                 : J.M. Molines in cdftools, f90 dor DRAKKAR (Nov. 2005)
  !! ---------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------

  !! * Used modules
  USE cdfio

  !! * Local declaration
  IMPLICIT NONE

  INTEGER :: npiglo, npjglo, npk, npkk ,npt  
  INTEGER :: narg, iargc
  INTEGER :: ji,jj,jk,jkk,k0
  INTEGER :: istartarg = 1
  INTEGER :: ncout, ierr
  INTEGER, DIMENSION(1) :: ipk, id_varout !: for output variables
  !
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE ::  v3d, alpha
  REAL(KIND=4), DIMENSION(:,:)  , ALLOCATABLE ::  zint
  REAL(KIND=4), DIMENSION(:)    , ALLOCATABLE ::  zi, time_tag, h1d
  REAL(KIND=4)                                ::  P1,P2
  REAL(KIND=4)                                ::  spval=999999.
  REAL(KIND=4)                                ::  spvalz=0.

  CHARACTER(LEN=256) ::  cline, cfilZI,  cfilsigma, cvar, cfilout
 
  TYPE(variable), DIMENSION(1)  :: typvar      !: structure for attributes
  !
  LOGICAL :: lsingle=.false.

  !! * Read command line
  narg=iargc()
  IF (narg < 2 ) THEN
         PRINT *, &
     &' >>>> usage: cdfisopycdep [-s sigma ] <rhofile> <cdfsigmavar> '
         PRINT *,'   Deptht of isopycnal surfaces will be in isopycdep.nc'
         PRINT *,'   Isopycnal value are read on a text file ''rho_lev'' '
         PRINT *,'   unless the option  -s is specified with one particular value.'
         PRINT *,'   Model density are taken on file ''rhofile'' with name cdfsigmavar'
         PRINT *,'   Output done on isopycdep.nc, var vodepiso'
         STOP
  ENDIF
  !   seek a -s option
  CALL getarg(1,cline)
  IF (cline == '-s' ) THEN
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
   
  ! Read Rho file
  CALL getarg(istartarg,cfilsigma)
  npiglo=getdim(cfilsigma,'x')
  npjglo=getdim(cfilsigma,'y')
  npk   =getdim(cfilsigma,'depth')
  npt   =getdim(cfilsigma,'time')

  ! Read variable name
  CALL getarg(istartarg+1,cvar)


  ALLOCATE( v3d(npiglo,npjglo,npk), alpha(npiglo, npjglo, npkk) )
  ALLOCATE( zint(npiglo,npjglo)  )
  ALLOCATE( time_tag(npt), h1d(npk) )

  time_tag(:)=getvar1d(cfilsigma,'time_counter', npt)
  h1d(:)=getvar1d(cfilsigma,'deptht',npk)
   
  DO jk=1,npk
     v3d(:,:,jk) = getvar(cfilsigma,cvar,jk,npiglo,npjglo)
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


     ! ... open output file and write header
     ipk(:)=npkk

     typvar(1)%name= 'vodepiso'
     typvar(1)%units='m'
     typvar(1)%missing_value=999999.
     typvar(1)%valid_min= 0.
     typvar(1)%valid_max= 7000.
     typvar(1)%long_name='Depth_of_Isopycnals'
     typvar(1)%short_name='vodepiso'
     typvar(1)%online_operation='N/A'
     typvar(1)%axis='TRYX'


     cfilout='isopycdep.nc'
 
     ncout = create(cfilout,cfilsigma ,npiglo,npjglo,npkk)
     ierr = createvar(ncout, typvar,1,ipk, id_varout )
     ierr = putheadervar(ncout , cfilsigma, npiglo, npjglo, npkk,pdep=zi)
 
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
                    zint (ji,jj)=alpha(ji,jj,jkk)*h1d(k0+1) &
   &                         +(1-alpha(ji,jj,jkk))*h1d(k0)
                ELSE 
                   zint  (ji,jj)=spval
               ENDIF
             ELSE 
               zint  (ji,jj)=spval
             ENDIF
           END DO
        END DO
        ierr = putvar(ncout,id_varout(1), zint  ,jkk,npiglo,npjglo)
     END DO
     ierr = putvar1d(ncout,time_tag,1,'T')
     ierr = closeout(ncout)
        PRINT *,'Projection on isopycns completed successfully'
END  PROGRAM cdfisopycdep
