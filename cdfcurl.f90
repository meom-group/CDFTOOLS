PROGRAM cdfcurl
  !!---------------------------------------------------------------------------
  !!         ***  PROGRAM  cdfcurl  ***
  !!
  !!  **  Purpose: Compute the curl on F-points for given gridU gridV files and variables
  !!
  !! history :
  !!   Original :  J.M. Molines (May 2005)
  !!   Modified :  P.   Mathiot (June 2007) update for forcing fields
  !!---------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER :: ji,jj,jk, jt, ilev
  INTEGER :: npiglo, npjglo, npk, nt
  INTEGER :: narg, iargc, ncout, ierr
  INTEGER, DIMENSION(1) ::  ipk, id_varout         ! 

  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE  :: e2v, e1u, e1f, e2f 
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: un, vn, rotn, fmask, zun, zvn
  REAL(KIND=4) ,DIMENSION(1)                 ::  tim

  LOGICAL :: lforcing    = .FALSE.

  CHARACTER(LEN=80) :: cfilu, cfilv, cvaru, cvarv, cdum
  CHARACTER(LEN=80) :: coord='mesh_hgr.nc', cfileout='curl.nc'
  TYPE (variable), DIMENSION(1) :: typvar         !: structure for attibutes

  !!
  narg = iargc()
  IF ( narg /= 5 ) THEN
     PRINT *,' USAGE : cdfcurl fileU fileV varU varV lev'
     PRINT *,'        lev is the level where the curl will be computed'
     PRINT *,'        Produce a cdf file curl.nc with socurl variable'
     PRINT *,'        Need mesh_hgr.nc'
     PRINT *,'                         '
     PRINT *,'        For forcing fields : '
     PRINT *,'        if no z dimension put lev=0 (PM)' 
     PRINT *,'        if file is in grid B or C, check the code (PM)' 
     STOP
  ENDIF

  CALL getarg(1, cfilu)
  CALL getarg(2, cfilv)
  CALL getarg(3, cvaru)
  CALL getarg(4, cvarv)
  CALL getarg(5, cdum)
  READ(cdum,*) ilev

  ! define new variables for output ( must update att.txt)
  typvar(1)%name='socurl'
  typvar(1)%units='s-1'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= -1000.
  typvar(1)%valid_max= 1000.
  typvar(1)%long_name='Relative_Vorticity (curl)'
  typvar(1)%short_name='socurl'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TYX'

  ipk(1) = 1  !  2D

  npiglo = getdim(cfilu,'x')
  npjglo = getdim(cfilu,'y')
  npk    = getdim(cfilu,'depth')
  nt     = getdim(cfilu,'time_counter') !PM

  PRINT *, 'npiglo =',npiglo
  PRINT *, 'npjglo =',npjglo
  PRINT *, 'npk    =',npk
  PRINT *, 'nt     =',nt  !PM
  PRINT *, 'ilev   =',ilev

  !test if lev exists
  IF ((npk==0) .AND. (ilev .GT. 0) ) THEN
     PRINT *, 'Problem : npk = 0 and lev > 0 STOP'
     STOP
  END IF

  ! if forcing field (PM)
  IF (ilev==0 .AND. npk==0) THEN
     lforcing = .TRUE.
     npk = 1
     PRINT *, 'npk =0, assume 1'
  END IF

  IF (nt==0) THEN
     PRINT *, 'nt=0, assume 1'
     nt=1
  END IF
  !end (PM)

  ! check files and determines if the curl will be 2D of 3D

  ! create output fileset
  ncout =create(cfileout, cfilu, npiglo,npjglo,0)
  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfilu, npiglo, npjglo, 0)


  ! Allocate the memory
  ALLOCATE ( e1u(npiglo,npjglo) , e1f(npiglo,npjglo) )
  ALLOCATE ( e2v(npiglo,npjglo) , e2f(npiglo,npjglo) )
  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo)  )
  ALLOCATE ( zun(npiglo,npjglo) , zvn(npiglo,npjglo) )
  ALLOCATE ( rotn(npiglo,npjglo) , fmask(npiglo,npjglo) )

  e1u=  getvar(coord, 'e1u', 1,npiglo,npjglo)
  e1f=  getvar(coord, 'e1f', 1,npiglo,npjglo)
  e2v=  getvar(coord, 'e2v', 1,npiglo,npjglo)
  e2f=  getvar(coord, 'e2f', 1,npiglo,npjglo)

  tim=getvar1d(cfilu,'time_counter',1)
  ierr=putvar1d(ncout,tim,1,'T')
  
  DO jt=1,nt

     ! begin PM
     IF (MOD(jt,100)==0) PRINT *, jt,'/',nt
     ! if files are forcing fields
     IF ( lforcing == .TRUE. ) THEN
        jk = jt
        zun(:,:) =  getvar(cfilu, cvaru, jk ,npiglo,npjglo, ktime=jt)
        zvn(:,:) =  getvar(cfilv, cvarv, jk ,npiglo,npjglo, ktime=jt)
        ! compute u and v on U and V point
        DO ji=1, npiglo-1
           un(ji,:) = 0.5*(zun(ji,:) + zun(ji+1,:))
        END DO
        !
        DO jj=1, npjglo-1
           vn(:,jj) = 0.5*(zvn(:,jj) + zvn(:,jj+1))
        END DO
        ! end compute u and v on U and V point
     ELSE IF (lforcing == .FALSE. ) THEN
        jk = ilev
        un(:,:) =  getvar(cfilu, cvaru, jk ,npiglo,npjglo, ktime=jt)
        vn(:,:) =  getvar(cfilv, cvarv, jk ,npiglo,npjglo, ktime=jt)
     END IF
     !end PM

     ! compute the mask
     IF (jt==1) THEN
        DO jj = 1, npjglo - 1
           DO ji = 1, npiglo - 1
              fmask(ji,jj)=0.
              fmask(ji,jj)= un(ji,jj)*un(ji,jj+1) * vn(ji,jj)*vn(ji+1,jj)
              IF (fmask(ji,jj) /= 0.) fmask(ji,jj)=1.
           ENDDO
        ENDDO
     END IF

     rotn(:,:) = 0.
     DO jj = 1, npjglo -1 
        DO ji = 1, npiglo -1   ! vector opt.
           rotn(ji,jj) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ) - e2v(ji,jj) * vn(ji,jj)    &
                &              - e1u(ji  ,jj+1) * un(ji  ,jj+1) + e1u(ji,jj) * un(ji,jj)  ) &
                &           * fmask(ji,jj) / ( e1f(ji,jj) * e2f(ji,jj) )
        END DO
     END DO
     ! 
     ! write rotn on file at level k and at time jt
     IF (lforcing== .TRUE. ) THEN
        ierr = putvar(ncout, id_varout(1) ,rotn, jt ,npiglo, npjglo, jt)
     ELSE IF (lforcing== .FALSE. ) THEN
        ierr = putvar(ncout, id_varout(1) ,rotn, 1 ,npiglo, npjglo, jt)
     END IF
  END DO
  ierr = closeout(ncout)

END PROGRAM cdfcurl

