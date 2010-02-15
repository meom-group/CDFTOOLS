PROGRAM cdfmoyuv
  !!---------------------------------------------------------------------------
  !!         ***  PROGRAM  cdfbti  ***
  !!
  !!  **  Purpose: Compute the temporal mean of u,v,u2,v2 and uv on the T point 
  !!               Useful for the programm cdfbti which compute the transfert of
  !!               barotropic instability energy
  !!
  !! history :
  !!   Original :  A. Melet (Feb 2008)
  !!---------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER :: ji,jj,jk, jt, ntframe, total_time, ilev
  INTEGER :: npiglo, npjglo, npk, nt, ntags
  INTEGER :: imin, imax, jmin, jmax, npil, npjl
  INTEGER :: narg, iargc, ncout, ierr
  INTEGER, DIMENSION(5) ::  ipk, id_varout         ! 

  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: u2d, v2d, uvn
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: un, vn, u2n, v2n
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: umask, vmask
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE  :: tabu, tabv, tabu2, tabv2, tabuv
  REAL(KIND=4) ,DIMENSION(1)                 ::  tim

  CHARACTER(LEN=256) :: cfileu, cfilev, cvaru, cvarv, config , ctag
  CHARACTER(LEN=256) :: cfileout='moyuv.nc', cdum

  TYPE (variable), DIMENSION(5) :: typvar         !: structure for attibutes

  !!
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' USAGE : cdfmoyuv config imin imax jmin jmax ''list_of_tags'' '
     PRINT *,'        Produce a cdf file moyuv.nc with uvbar variable on T point'
     PRINT *,'        and ubar, u2bar on U point, and vbar, v2bar on V point '
     PRINT *,'        for the region defined by imin imax jmin jmax'
     PRINT *,' '
     STOP
  ENDIF

  ntags = narg - 5   ! first five arguments are not tags
  !! Initialisation from 1st file (all file are assume to have the same
  !geometry)
  CALL getarg (1, config)
  CALL getarg (2, cdum) ; READ(cdum,*) imin
  CALL getarg (3, cdum) ; READ(cdum,*) imax
  CALL getarg (4, cdum) ; READ(cdum,*) jmin
  CALL getarg (5, cdum) ; READ(cdum,*) jmax
  CALL getarg (6, ctag)
  WRITE(cfileu,'(a,"_",a,"_gridU.nc")') TRIM(config),TRIM(ctag)

  PRINT *,TRIM(cfileu)        
  npiglo = getdim (cfileu,'x')
  npjglo = getdim (cfileu,'y')
  npk    = getdim (cfileu,'depth')
  nt     = getdim(cfileu,'time_counter')

  IF (imin /= 0 ) THEN ; npiglo=imax -imin + 1;  ELSE ; imin=1 ; ENDIF
  IF (jmin /= 0 ) THEN ; npjglo=jmax -jmin + 1;  ELSE ; jmin=1 ; ENDIF

  ! define new variables for output ( must update att.txt)

  typvar(1)%name='ubar'
  typvar(1)%long_name='temporal mean of u on U point'
  typvar(1)%short_name='ubar'

  typvar(2)%name='vbar'
  typvar(2)%long_name='temporal mean of v on V point'
  typvar(2)%short_name='vbar'

  typvar(3)%name='u2bar'
  typvar(3)%long_name='temporal mean of u * u on U point'
  typvar(3)%short_name='u2bar'

  typvar(4)%name='v2bar'
  typvar(4)%long_name='temporal mean of v * v on V point'
  typvar(4)%short_name='v2bar'

  typvar(5)%name='uvbar'
  typvar(5)%long_name='temporal mean of u * v on T point'
  typvar(5)%short_name='uvbar'

  typvar%units='m2.s-2'
  typvar%missing_value=0.
  typvar%valid_min= -1000.
  typvar%valid_max= 1000.
  typvar%online_operation='N/A'
  typvar%axis='TYX'

  ipk(1) = npk  
  ipk(2) = npk
  ipk(3) = npk
  ipk(4) = npk
  ipk(5) = npk

  PRINT *, 'npiglo =',npiglo
  PRINT *, 'npjglo =',npjglo
  PRINT *, 'npk    =',npk
  PRINT *, 'nt     =',nt 

  !test if lev exists
  IF ((npk==0) .AND. (ilev .GT. 0) ) THEN
     PRINT *, 'Problem : npk = 0 and lev > 0 STOP'
     STOP
  END IF

  ! create output fileset
  ncout =create(cfileout, cfileu, npiglo,npjglo,npk)
  ierr= createvar(ncout ,typvar,5, ipk,id_varout )
  ierr= putheadervar(ncout, cfileu, npiglo, npjglo,npk)
        
  ! Allocate the memory
  ALLOCATE ( u2d(npiglo,npjglo)  , v2d(npiglo,npjglo)  )
  ALLOCATE ( un(npiglo,npjglo) , tabu(npiglo,npjglo) )
  ALLOCATE ( vn(npiglo,npjglo) , tabv(npiglo,npjglo) )
  ALLOCATE ( u2n(npiglo,npjglo) , tabu2(npiglo,npjglo) )
  ALLOCATE ( v2n(npiglo,npjglo) , tabv2(npiglo,npjglo) )
  ALLOCATE ( uvn(npiglo,npjglo) , tabuv(npiglo,npjglo) )
  ALLOCATE ( umask(npiglo,npjglo) , vmask(npiglo,npjglo) )

  DO jk=1, npk
     PRINT *,'            level ',jk
     total_time = 0.;  ntframe=0 
     tabu(:,:) = 0.d0 ; tabv(:,:) = 0.d0 ; tabuv(:,:) = 0.d0 
     tabu2(:,:) = 0.d0 ; tabv2(:,:) = 0.d0 
     tim=getvar1d(cfileu,'time_counter',nt)
     total_time = total_time + SUM(tim(1:nt) )
     un(:,:)  =  0.d0
     vn(:,:)  =  0.d0
     u2n(:,:) =  0.d0
     v2n(:,:) =  0.d0
     uvn(:,:) =  0.d0

     DO jt= 6, narg
        ntframe=ntframe+1
        CALL getarg (jt, ctag)
        WRITE(cfileu,'(a,"_",a,"_gridU.nc")') TRIM(config),TRIM(ctag)
        WRITE(cfilev,'(a,"_",a,"_gridV.nc")') TRIM(config),TRIM(ctag)
        
        u2d(:,:)= getvar(cfileu, 'vozocrtx', jk , &
                  & npiglo, npjglo, kimin=imin, kjmin=jmin, ktime=1 )
        v2d(:,:)= getvar(cfilev, 'vomecrty', jk , &
                  & npiglo, npjglo, kimin=imin, kjmin=jmin, ktime=1 )
        
        tabu(:,:)  = tabu(:,:)  + u2d(:,:)
        tabu2(:,:) = tabu2(:,:) + u2d(:,:) * u2d(:,:)
        tabv(:,:)  = tabv(:,:)  + v2d(:,:)
        tabv2(:,:) = tabv2(:,:) + v2d(:,:) * v2d(:,:)
        
        DO jj = jmin+1, npjglo
           DO ji = imin+1, npiglo
              umask(ji,jj)=0.
              umask(ji,jj)=u2d(ji,jj)*u2d(ji-1,jj)
              vmask(ji,jj)=0.
              vmask(ji,jj)=v2d(ji,jj)*v2d(ji,jj-1)
              IF (umask(ji,jj) /= 0.) umask(ji,jj)=1.
              IF (vmask(ji,jj) /= 0.) vmask(ji,jj)=1.   

              tabuv(ji-imin,jj-jmin) = tabuv(ji-imin,jj-jmin) &
              &   + 0.5 * umask(ji,jj) * (u2d(ji,jj)+u2d(ji-1,jj)) &
              &   * 0.5 * vmask(ji,jj) * (v2d(ji,jj)+v2d(ji,jj-1)) 
           END DO
        END DO
     END DO
     un(:,:)  = tabu(:,:)  / ntframe
     vn(:,:)  = tabv(:,:)  / ntframe
     u2n(:,:) = tabu2(:,:) / ntframe
     v2n(:,:) = tabv2(:,:) / ntframe
     uvn(:,:) = tabuv(:,:) / ntframe
     ! sauvegarde
     ierr = putvar(ncout, id_varout(1) ,un, jk, npiglo, npjglo, &
                   &      ktime=1)
     ierr = putvar(ncout, id_varout(2) ,vn, jk, npiglo, npjglo, &
                   &      ktime=1)
     ierr = putvar(ncout, id_varout(3) ,u2n, jk, npiglo, npjglo, &
                   &      ktime=1)
     ierr = putvar(ncout, id_varout(4) ,v2n, jk, npiglo, npjglo, &
                   &      ktime=1)
     ierr = putvar(ncout, id_varout(5) ,uvn, jk, npiglo, npjglo, &
                   &      ktime=1)
  END DO
  ierr = closeout(ncout)

END PROGRAM cdfmoyuv

