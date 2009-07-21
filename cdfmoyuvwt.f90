PROGRAM cdfmoyuvwt
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
  INTEGER :: ji,jj,jk, jt, ntframe, ilev
  INTEGER :: npiglo, npjglo, npk, nt, ntags
  INTEGER :: imin, imax, jmin, jmax, npil, npjl
  INTEGER :: narg, iargc, ncout, ierr
  INTEGER, DIMENSION(11) ::  ipk, id_varout         ! 

  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: u2d, v2d
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: w2d, t2d
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: wxz, txz
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: un, vn, u2n, v2n, uvn
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: wn, tn, utn, vtn, t2n
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: wtn, zzz
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: umask, vmask, tmask, wmask
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: t1mask, w1mask
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE  :: tabu, tabv, tabu2, tabv2, tabuv
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE  :: tabw, tabt, tabut, tabvt, tabt2
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE  :: tabwt
  REAL(KIND=4) ,DIMENSION(1)                 ::  tim
  REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE  :: wtab
  REAL(KIND=8)                               :: total_time

  CHARACTER(LEN=256) :: cfileu, cfilev, cvaru, cvarv, config , ctag
  CHARACTER(LEN=256) :: cfilew, cfilet, cavarw, cvart
  CHARACTER(LEN=256) :: cfileout='moyuvwt.nc', cdum
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE  :: ctabtag

  TYPE (variable), DIMENSION(11) :: typvar         !: structure for attibutes

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
  ALLOCATE (ctabtag(ntags) )
  !! Initialisation from 1st file (all file are assume to have the same
  !geometry)
  CALL getarg (1, config)
  CALL getarg (2, cdum) ; READ(cdum,*) imin
  CALL getarg (3, cdum) ; READ(cdum,*) imax
  CALL getarg (4, cdum) ; READ(cdum,*) jmin
  CALL getarg (5, cdum) ; READ(cdum,*) jmax
  CALL getarg (6, ctag)
  WRITE(cfileu,'(a,"_",a,"_gridU.nc")') TRIM(config),TRIM(ctag)
  ctabtag(1)=ctag

  DO jt=7,narg
    CALL getarg(jt,ctag)
    ctabtag(jt-5)=ctag
  END DO

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

  typvar(6)%name='wbar'
  typvar(6)%long_name='temporal mean of w on W point'
  typvar(6)%short_name='wbar'
   
  typvar(7)%name='tbar'
  typvar(7)%long_name='temporal mean of T on T point in K'
  typvar(7)%short_name='tbar'

  typvar(8)%name='utbar'
  typvar(8)%long_name='temporal mean of u * T (in K) on T point'
  typvar(8)%short_name='utbar'

  typvar(9)%name='vtbar'
  typvar(9)%long_name='temporal mean of v * T (in K) on T point'
  typvar(9)%short_name='vtbar'
 
  typvar(10)%name='t2bar'
  typvar(10)%long_name='temporal mean of T * T on T point in K^2'
  typvar(10)%short_name='t2bar'
      
  typvar(11)%name='wtbar'
  typvar(11)%long_name='temporal mean of w * T (in K) on T point'
  typvar(11)%short_name='wtbar'

  typvar%units=''
  typvar%missing_value=0.
  typvar%valid_min= -1000.
  typvar%valid_max= 1000.
  typvar%online_operation='N/A'
  typvar%axis='TYX'

  ipk(:) = npk  

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
  ierr= createvar(ncout ,typvar,11, ipk,id_varout )
  ierr= putheadervar(ncout, cfileu, npiglo, npjglo,npk)
        
  ! Allocate the memory
  ALLOCATE ( u2d(npiglo,npjglo)  , v2d(npiglo,npjglo)  )
  ALLOCATE ( t2d(npiglo,npjglo)  , w2d(npiglo,npjglo)  )
  ALLOCATE ( un(npiglo,npjglo) , tabu(npiglo,npjglo) )
  ALLOCATE ( vn(npiglo,npjglo) , tabv(npiglo,npjglo) )
  ALLOCATE ( u2n(npiglo,npjglo) , tabu2(npiglo,npjglo) )
  ALLOCATE ( v2n(npiglo,npjglo) , tabv2(npiglo,npjglo) )
  ALLOCATE ( t2n(npiglo,npjglo) , tabt2(npiglo,npjglo) )
  ALLOCATE ( uvn(npiglo,npjglo) , tabuv(npiglo,npjglo) )
  ALLOCATE ( tn(npiglo,npjglo) , tabt(npiglo,npjglo) )
  ALLOCATE ( wn(npiglo,npjglo) , tabw(npiglo,npjglo) )
  ALLOCATE ( utn(npiglo,npjglo) , tabut(npiglo,npjglo) )
  ALLOCATE ( vtn(npiglo,npjglo) , tabvt(npiglo,npjglo) )
  ALLOCATE ( wtn(npiglo,npk) , tabwt(npiglo,npk), zzz(npiglo,1) )
   
  ALLOCATE ( umask(npiglo,npjglo) , vmask(npiglo,npjglo) )
  ALLOCATE ( tmask(npiglo,npjglo) , wmask(npiglo,npjglo) )
  ALLOCATE ( t1mask(npiglo,npk) , w1mask(npiglo,npk) )
  ALLOCATE ( txz(npiglo,npk) , wxz(npiglo,npk) )
  ALLOCATE ( wtab(npiglo, npjglo, npk) )
  
  DO jk=1, npk
     PRINT *,'            level ',jk
     total_time = 0.d0;  ntframe=0 
     tabu(:,:) = 0.d0 ; tabv(:,:) = 0.d0 ; tabuv(:,:) = 0.d0 
     tabu2(:,:) = 0.d0 ; tabv2(:,:) = 0.d0 ; tabt(:,:) = 0.d0 
     tabw(:,:) = 0.d0 ; tabut(:,:) = 0.d0 ; tabvt(:,:) = 0.d0 
     tabt2(:,:) = 0.d0
     un(:,:)  =  0.d0
     vn(:,:)  =  0.d0
     u2n(:,:) =  0.d0
     v2n(:,:) =  0.d0
     uvn(:,:) =  0.d0
     tn(:,:)  =  0.d0
     wn(:,:)  =  0.d0
     utn(:,:)  =  0.d0
     vtn(:,:)  =  0.d0
     t2n(:,:)  =  0.d0

     DO jt= 6, narg
        ntframe=ntframe+1
        ctag=ctabtag(jt-5)
        WRITE(cfileu,'(a,"_",a,"_gridU.nc")') TRIM(config),TRIM(ctag)
        WRITE(cfilev,'(a,"_",a,"_gridV.nc")') TRIM(config),TRIM(ctag)
        WRITE(cfilew,'(a,"_",a,"_gridW.nc")') TRIM(config),TRIM(ctag)
        WRITE(cfilet,'(a,"_",a,"_gridT.nc")') TRIM(config),TRIM(ctag)
        IF ( jk == 1 ) THEN
         tim=getvar1d(cfileu,'time_counter',nt)
         total_time = total_time + SUM(tim(1:nt) )
        ENDIF
        
        u2d(:,:)= getvar(cfileu, 'vozocrtx', jk , &
                  & npiglo, npjglo, kimin=imin, kjmin=jmin, ktime=1 )
        v2d(:,:)= getvar(cfilev, 'vomecrty', jk , &
                  & npiglo, npjglo, kimin=imin, kjmin=jmin, ktime=1 )
        w2d(:,:)= getvar(cfilew, 'vovecrtz', jk , &
                  & npiglo, npjglo, kimin=imin, kjmin=jmin, ktime=1 )
        t2d(:,:)= getvar(cfilet, 'votemper', jk , &
                  & npiglo, npjglo, kimin=imin, kjmin=jmin, ktime=1 )
        
        tabu(:,:)  = tabu(:,:)  + u2d(:,:)
        tabu2(:,:) = tabu2(:,:) + u2d(:,:) * u2d(:,:)
        tabv(:,:)  = tabv(:,:)  + v2d(:,:)
        tabv2(:,:) = tabv2(:,:) + v2d(:,:) * v2d(:,:)
        tabw(:,:)  = tabw(:,:)  + w2d(:,:)
        tabt(:,:)  = tabt(:,:)  + (t2d(:,:)+273.15)
        tabt2(:,:) = tabt2(:,:) + (t2d(:,:)+273.15)*(t2d(:,:)+273.15)
        
        DO jj = jmin+1, npjglo
           DO ji = imin+1, npiglo
              umask(ji,jj)=0.
              umask(ji,jj)=u2d(ji,jj)*u2d(ji-1,jj)
              vmask(ji,jj)=0.
              vmask(ji,jj)=v2d(ji,jj)*v2d(ji,jj-1)
              wmask(ji,jj)=0.
              wmask(ji,jj)=w2d(ji,jj)
              tmask(ji,jj)=0.
              tmask(ji,jj)=t2d(ji,jj)
              IF (umask(ji,jj) /= 0.) umask(ji,jj)=1.
              IF (vmask(ji,jj) /= 0.) vmask(ji,jj)=1.   
              IF (tmask(ji,jj) /= 0.) tmask(ji,jj)=1.
              IF (wmask(ji,jj) /= 0.) wmask(ji,jj)=1.

              tabuv(ji-imin,jj-jmin) = tabuv(ji-imin,jj-jmin) &
              &   + 0.5 * umask(ji,jj) * (u2d(ji,jj)+u2d(ji-1,jj)) &
              &   * 0.5 * vmask(ji,jj) * (v2d(ji,jj)+v2d(ji,jj-1)) 
              tabut(ji-imin,jj-jmin) = tabut(ji-imin,jj-jmin) &
              &   + 0.5 * umask(ji,jj) * (u2d(ji,jj)+u2d(ji-1,jj)) &
              &   * tmask(ji,jj) * (t2d(ji,jj)+273.15)
              tabvt(ji-imin,jj-jmin) = tabvt(ji-imin,jj-jmin) &
              &   + 0.5 * vmask(ji,jj) * (v2d(ji,jj)+v2d(ji,jj-1)) &
              &   * tmask(ji,jj) * (t2d(ji,jj)+273.15)
              
           END DO
        END DO
     END DO
     
     un(:,:)  = tabu(:,:)  / ntframe
     vn(:,:)  = tabv(:,:)  / ntframe
     u2n(:,:) = tabu2(:,:) / ntframe
     v2n(:,:) = tabv2(:,:) / ntframe
     uvn(:,:) = tabuv(:,:) / ntframe
     tn(:,:)  = tabt(:,:)  / ntframe
     wn(:,:)  = tabw(:,:)  / ntframe
     utn(:,:) = tabut(:,:) / ntframe
     vtn(:,:) = tabvt(:,:) / ntframe
     t2n(:,:) = tabt2(:,:) / ntframe
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
     ierr = putvar(ncout, id_varout(6) ,wn, jk, npiglo, npjglo, &
                   &      ktime=1)
     ierr = putvar(ncout, id_varout(7) ,tn, jk, npiglo, npjglo, &
                   &      ktime=1)
     ierr = putvar(ncout, id_varout(8) ,utn, jk, npiglo, npjglo, &
                   &      ktime=1)
     ierr = putvar(ncout, id_varout(9) ,vtn, jk, npiglo, npjglo, &
                   &      ktime=1)
     ierr = putvar(ncout, id_varout(10) ,t2n, jk, npiglo, npjglo, &
                   &      ktime=1)


  END DO
     ! this is done for everybody, dont bother the time for next variables
     tim(1)= total_time/ntframe
     ierr=putvar1d(ncout,tim,1,'T')

  DO jj = jmin+1, jmax   ! jj global
     print *, 'JJ=',jj
     ntframe=0
     tabwt(:,:) = 0.d0
     wtn(:,:)   =  0.d0
     wxz(:,:)   =  0.d0
     txz(:,:)   =  0.d0
    DO jt= 6, narg
       ntframe=ntframe+1
       ctag=ctabtag(jt-5)
       WRITE(cfilet,'(a,"_",a,"_gridT.nc")') TRIM(config),TRIM(ctag)
       WRITE(cfilew,'(a,"_",a,"_gridW.nc")') TRIM(config),TRIM(ctag)

       wxz(:,:)=getvarxz(cfilew,'vovecrtz',jj,npiglo,npk, kimin=imin,kkmin=1,ktime=1)
       txz(:,:)=getvarxz(cfilet,'votemper',jj,npiglo,npk, kimin=imin,kkmin=1,ktime=1)

         DO jk=1, npk-1
            w1mask(:,jk) = wxz(:,jk) * wxz(:,jk+1)
            t1mask(:,jk) = txz(:,jk)
            WHERE ( w1mask(:,jk) /= 0.) w1mask(:,jk)=1.
            WHERE ( t1mask(:,jk) /= 0.) t1mask(:,jk)=1.
            tabwt(:,jk) = tabwt(:,jk) + t1mask(:,jk)*(txz(:,jk)+273.15) &
                          & *0.5* w1mask(:,jk)* ( wxz(:,jk) + wxz(:,jk+1))
        END DO
     END DO
     wtn(:,:) = tabwt(:,:) / ntframe
     wtab(:,jj-jmin+1,:)= wtn(:,:)
  END DO       
  DO jk=1,npk
     ierr = putvar(ncout, id_varout(11) ,wtab(:,:,jk), jk, npiglo, npjglo )
  END DO

  ierr = closeout(ncout)

END PROGRAM cdfmoyuvwt

