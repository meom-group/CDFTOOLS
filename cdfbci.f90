PROGRAM cdfbci
  !!---------------------------------------------------------------------------
  !!         ***  PROGRAM  cdfbci  ***
  !!
  !!  **  Purpose: Compute the term of energetic transfert BCI 
  !!      for the baroclinic instability for given gridU gridV gridU2 gridV2 files and variables
  !!      The intput file is the result of a pre-processing by cdfmoyuvwt
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
  INTEGER :: ji,jj,jk, jt, ilev
  INTEGER :: npiglo, npjglo, npk, nt
  INTEGER :: narg, iargc, ncout, ierr
  INTEGER, DIMENSION(5) ::  ipk, id_varout         ! 

  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE  :: e2t, e1t 
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: anout, anovt, un, vn, tn
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: utn,vtn
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: tmask, umask, vmask
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: bci
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: dtdx, dtdy
  REAL(KIND=4) ,DIMENSION(1)                 :: tim

  CHARACTER(LEN=256) :: cfile
  CHARACTER(LEN=256) :: coord='mesh_hgr.nc', cfileout='bci.nc'
  TYPE (variable), DIMENSION(5) :: typvar         !: structure for attibutes

  !!
  narg = iargc()
  IF ( narg /= 1 ) THEN
     PRINT *,' USAGE : cdfbci file'
     PRINT *,'        Produce a cdf file bci.nc with bci variable'
     PRINT *,'        file is from cdfmoyuvwt'
     PRINT *,'        the mean must have been computed on a period long enough'
     PRINT *,'        for the statistics to be meaningful'
     PRINT *,'        Need mesh_hgr.nc'
     PRINT *,'                         '
     PRINT *,'        if file is in grid B or C, check the code (PM)' 
     STOP
  ENDIF

  CALL getarg(1, cfile)
  npiglo = getdim(cfile,'x')
  npjglo = getdim(cfile,'y')
  npk    = getdim(cfile,'depth')
  nt     = getdim(cfile,'time_counter')

  PRINT *, 'npiglo =',npiglo
  PRINT *, 'npjglo =',npjglo
  PRINT *, 'npk    =',npk
  PRINT *, 'nt     =',nt

  ! define new variables for output ( must update att.txt)
  typvar(1)%name='dTdx'
  typvar(1)%long_name='zonal derivate of Tbar on T point (*1000)'
  typvar(1)%short_name='dTdx'

  typvar(2)%name='dTdy'
  typvar(2)%long_name='meridional derivate of Tbar on T point (*1000)'
  typvar(2)%short_name='dTdy'
 
  typvar(3)%name='uT'
  typvar(3)%long_name='anomaly of u times anomaly of T on T point'
  typvar(3)%short_name='uT'
  
  typvar(4)%name='vT'
  typvar(4)%long_name='anomaly of v times anomaly of T on T point'
  typvar(4)%short_name='vT'

  typvar(5)%name='bci'
  typvar(5)%long_name='transfert of energy for the baroclinic instability (*1000)'
  typvar(5)%short_name='bci'

  typvar%units='1000 (u"T" dT/dx + v"T" dT/dy)'
  typvar%missing_value=0.
  typvar%valid_min= -1000.
  typvar%valid_max= 1000.
  typvar%online_operation='N/A'
  typvar%axis='TYX'

  ipk(:) = npk  
  
  !test if lev exists
  IF ((npk==0) .AND. (ilev .GT. 0) ) THEN
     PRINT *, 'Problem : npk = 0 and lev > 0 STOP'
     STOP
  END IF

  ! create output fileset
  ncout =create(cfileout, cfile, npiglo,npjglo,npk)
  ierr= createvar(ncout ,typvar,5, ipk,id_varout )
  ierr= putheadervar(ncout, cfile, npiglo, npjglo,npk)

  ! Allocate the memory
  ALLOCATE ( e1t(npiglo,npjglo) , e2t(npiglo,npjglo) )
  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo)  )
  ALLOCATE ( tn(npiglo,npjglo) )
  ALLOCATE ( utn(npiglo,npjglo) , vtn(npiglo,npjglo) )
  ALLOCATE ( umask(npiglo,npjglo) , vmask(npiglo,npjglo) )
  ALLOCATE ( tmask(npiglo,npjglo) )
  ALLOCATE ( dtdx(npiglo,npjglo)  , dtdy(npiglo,npjglo)  )
  ALLOCATE ( anout(npiglo,npjglo)  , anovt(npiglo,npjglo)  )
  ALLOCATE ( bci(npiglo,npjglo) )

  e1t=  getvar(coord, 'e1t', 1,npiglo,npjglo)
  e2t=  getvar(coord, 'e2t', 1,npiglo,npjglo)

  tim=getvar1d(cfile,'time_counter',nt)
  ierr=putvar1d(ncout,tim,1,'T')
     
  DO jk=1, npk
     PRINT *,'            level ',jk
           
     dtdx(:,:) = 0.d0
     dtdy(:,:) = 0.d0
        
     anovt(:,:) = 0.d0      
     anout(:,:) = 0.d0
     un(:,:)  = 0.d0
     vn(:,:)  = 0.d0
     tn(:,:)  = 0.d0
     utn(:,:) = 0.d0
     vtn(:,:) = 0.d0

     un(:,:)  =  getvar(cfile, 'ubar', jk ,npiglo,npjglo, ktime=1)
     vn(:,:)  =  getvar(cfile, 'vbar', jk ,npiglo,npjglo, ktime=1)
     tn(:,:)  =  getvar(cfile, 'tbar', jk ,npiglo,npjglo, ktime=1)
     utn(:,:) =  getvar(cfile, 'utbar', jk ,npiglo,npjglo, ktime=1)
     vtn(:,:) =  getvar(cfile, 'vtbar', jk ,npiglo,npjglo, ktime=1)

     ! compute the mask
     DO jj = 2, npjglo
        DO ji = 2, npiglo
           umask(ji,jj)=0.
           vmask(ji,jj)=0.
           tmask(ji,jj)=0.
           umask(ji,jj)= un(ji,jj)*un(ji-1,jj) 
           vmask(ji,jj)= vn(ji,jj)*vn(ji,jj-1)
           tmask(ji,jj)= tn(ji,jj)
           IF (umask(ji,jj) /= 0.) umask(ji,jj)=1.
           IF (vmask(ji,jj) /= 0.) vmask(ji,jj)=1.
           IF (tmask(ji,jj) /= 0.) tmask(ji,jj)=1.
        END DO
     END DO

     DO jj = 2, npjglo  
        DO ji = 2, npiglo    ! vector opt.
           ! calcul des dérivées au point T
           dtdx(ji,jj) = 1000/2 *( ( tn(ji,jj ) - tn(ji-1,jj) )        &
                &           * tmask(ji,jj)*tmask(ji-1,jj)           &
                &           / ( 0.5* ( e1t(ji,jj) + e1t(ji-1,jj) )) &
                &           +( tn(ji+1,jj ) - tn(ji,jj) )           &
                &           * tmask(ji+1,jj)*tmask(ji,jj)           &
                &           / ( 0.5* ( e1t(ji+1,jj) + e1t(ji,jj) )))

           dtdy(ji,jj) = 1000/2 *( ( tn(ji,jj) - tn(ji,jj-1) )         &
                &           * tmask(ji,jj)*tmask(ji,jj-1)           &
                &           / ( 0.5* ( e2t(ji,jj) + e2t(ji,jj-1) )) &
                &           +( tn(ji,jj+1 ) - tn(ji,jj) )           &
                &           * tmask(ji,jj+1)*tmask(ji,jj)           &
                &           / ( 0.5* ( e2t(ji,jj+1) + e2t(ji,jj) )) )

           anout(ji,jj)    = ( utn(ji,jj) &
                &                 -   1/2 * umask(ji,jj)*( un(ji,jj) + un(ji-1,jj) ) &
                &                   * tmask(ji,jj) * tn(ji,jj) )

           anovt(ji,jj)    = ( vtn(ji,jj) &
                &                 -   1/2 * vmask(ji,jj)*( vn(ji,jj) + vn(ji,jj-1) ) &
                &                   * tmask(ji,jj) * tn(ji,jj) )

           ! calcul du terme bti
           bci(ji,jj) = ( anout(ji,jj) * dtdx(ji,jj) &
                &         + anovt(ji,jj) * dtdy(ji,jj) )

        END DO
     END DO
     ! 
     ierr = putvar(ncout, id_varout(1) ,dtdx, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(2) ,dtdy, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(3) ,anout, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(4) ,anovt, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(5) ,bci, jk, npiglo, npjglo, ktime=1)
  END DO
  ierr = closeout(ncout)

END PROGRAM cdfbci

