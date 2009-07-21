PROGRAM cdfbti
  !!---------------------------------------------------------------------------
  !!         ***  PROGRAM  cdfbti  ***
  !!
  !!  **  Purpose: Compute the term of energetic transfert BTI 
  !!      for the barotropic instability for given gridU gridV gridU2 gridV2 files and variables
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
  INTEGER, DIMENSION(8) ::  ipk, id_varout         ! 

  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE  :: e2t, e1t, e1f, e2f 
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: un, vn, u2n, v2n, uvn
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: fmask, umask, vmask
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: anousqrt, anovsqrt, anouv, bti
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: dudx, dudy, dvdx, dvdy
  REAL(KIND=4) ,DIMENSION(1)                 :: tim

  CHARACTER(LEN=256) :: cfile
  CHARACTER(LEN=256) :: coord='mesh_hgr.nc', cfileout='bti.nc'
  TYPE (variable), DIMENSION(8) :: typvar         !: structure for attibutes

  !!
  narg = iargc()
  IF ( narg /= 1 ) THEN
     PRINT *,' USAGE : cdfbti file'
     PRINT *,'        Produce a cdf file bti.nc with bti variable'
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
  typvar(1)%name='dudx'
  typvar(1)%long_name='zonal derivate of u on T point'
  typvar(1)%short_name='dudx'

  typvar(2)%name='dvdx'
  typvar(2)%long_name='zonal derivate of v on T point'
  typvar(2)%short_name='dvdx'
 
  typvar(3)%name='dudy'
  typvar(3)%long_name='meridional derivate of u on T point'
  typvar(3)%short_name='dudy'
  
  typvar(4)%name='dvdy'
  typvar(4)%long_name='meridional derivate of v on T point'
  typvar(4)%short_name='dvdy'

  typvar(5)%name='anousqrt'
  typvar(5)%long_name='temporal mean of the square of the zonal speed anomaly'
  typvar(5)%short_name='anousqrt'

  typvar(6)%name='anovsqrt'
  typvar(6)%long_name='temporal mean of the square of the meridional speed anomaly'
  typvar(6)%short_name='anovsqrt'

  typvar(7)%name='anouv'
  typvar(7)%long_name='temporal mean of the Reynolds term'
  typvar(7)%short_name='anouanov'

  typvar(8)%name='bti'
  typvar(8)%long_name='transfert of energy for the barotropic instability'
  typvar(8)%short_name='bti'

  typvar%units='100000 s-1'
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
  ipk(6) = npk  
  ipk(7) = npk  
  ipk(8) = npk
  
  !test if lev exists
  IF ((npk==0) .AND. (ilev .GT. 0) ) THEN
     PRINT *, 'Problem : npk = 0 and lev > 0 STOP'
     STOP
  END IF

  ! create output fileset
  ncout =create(cfileout, cfile, npiglo,npjglo,npk)
  ierr= createvar(ncout ,typvar,8, ipk,id_varout )
  ierr= putheadervar(ncout, cfile, npiglo, npjglo,npk)

  ! Allocate the memory
  ALLOCATE ( e1t(npiglo,npjglo) , e1f(npiglo,npjglo) )
  ALLOCATE ( e2t(npiglo,npjglo) , e2f(npiglo,npjglo) )
  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo)  )
  ALLOCATE ( fmask(npiglo,npjglo) )
  ALLOCATE ( umask(npiglo,npjglo) , vmask(npiglo,npjglo) )
  ALLOCATE ( dudx(npiglo,npjglo)  , dudy(npiglo,npjglo)  )
  ALLOCATE ( dvdx(npiglo,npjglo)  , dvdy(npiglo,npjglo)  )
  ALLOCATE ( u2n(npiglo,npjglo)  , v2n(npiglo,npjglo)  )
  ALLOCATE ( uvn(npiglo,npjglo) )
  ALLOCATE ( anousqrt(npiglo,npjglo) , anovsqrt(npiglo,npjglo)  )
  ALLOCATE ( anouv(npiglo,npjglo), bti(npiglo,npjglo) )

  e1t=  getvar(coord, 'e1t', 1,npiglo,npjglo)
  e1f=  getvar(coord, 'e1f', 1,npiglo,npjglo)
  e2t=  getvar(coord, 'e2t', 1,npiglo,npjglo)
  e2f=  getvar(coord, 'e2f', 1,npiglo,npjglo)

  tim=getvar1d(cfile,'time_counter',nt)
  ierr=putvar1d(ncout,tim,1,'T')
     
  DO jk=1, npk
     PRINT *,'            level ',jk
           
     dudx(:,:) = 0.d0
     dvdx(:,:) = 0.d0
     dudy(:,:) = 0.d0
     dvdy(:,:) = 0.d0
        
     anousqrt(:,:) = 0.d0
     anovsqrt(:,:) = 0.d0      
     anouv(:,:)    = 0.d0
     
     un(:,:)  =  getvar(cfile, 'ubar', jk ,npiglo,npjglo, ktime=1)
     vn(:,:)  =  getvar(cfile, 'vbar', jk ,npiglo,npjglo, ktime=1)
     u2n(:,:) =  getvar(cfile, 'u2bar', jk ,npiglo,npjglo, ktime=1)
     v2n(:,:) =  getvar(cfile, 'v2bar', jk ,npiglo,npjglo, ktime=1)
     uvn(:,:) =  getvar(cfile, 'uvbar', jk ,npiglo,npjglo, ktime=1)

     ! compute the mask
     DO jj = 2, npjglo
        DO ji = 2, npiglo
           umask(ji,jj)=0.
           vmask(ji,jj)=0.
           umask(ji,jj)= un(ji,jj)*un(ji-1,jj) 
           vmask(ji,jj)= vn(ji,jj)*vn(ji,jj-1)
           IF (umask(ji,jj) /= 0.) umask(ji,jj)=1.
           IF (vmask(ji,jj) /= 0.) vmask(ji,jj)=1.
        ENDDO
     ENDDO
     DO jj = 1, npjglo-1
        DO ji = 1, npiglo-1
           fmask(ji,jj)=0.
           fmask(ji,jj)= un(ji,jj)*un(ji,jj+1) * vn(ji,jj)*vn(ji+1,jj)
           IF (fmask(ji,jj) /= 0.) fmask(ji,jj)=1.
        ENDDO
     ENDDO



     DO jj = 2, npjglo  
        DO ji = 2, npiglo    ! vector opt.
           ! calcul des dérivées au point T
           dudx(ji,jj) = 100000 * ( un(ji,jj ) -  un(ji-1,jj) )   &
                &           * umask(ji,jj) / e1t(ji,jj) 

           dvdy(ji,jj) = 100000 * ( vn(ji,jj) - vn(ji,jj-1) )   &
                &           * vmask(ji,jj) / e2t(ji,jj)             

           dudy(ji,jj) = 100000/4 *( ( un(ji,jj+1 ) - un(ji,jj) )   &
                &           * fmask(ji,jj) / e2f(ji,jj) &
                        + (un(ji,jj ) - un(ji,jj-1) )   &
                &           * fmask(ji,jj-1) / e2f(ji,jj-1) &
                        + (un(ji-1,jj+1 ) - un(ji-1,jj) )   &
                &           * fmask(ji-1,jj) / e2f(ji-1,jj) &         
                        + (un(ji-1,jj ) - un(ji-1,jj-1) )   &    
                &           * fmask(ji-1,jj-1) / e2f(ji-1,jj-1) )            

           dvdx(ji,jj) = 100000/4 *( ( vn(ji,jj ) - vn(ji-1,jj) )   &
                &           * fmask(ji-1,jj) / e1f(ji-1,jj) &
                        + (vn(ji+1,jj ) - vn(ji,jj) )   &
                &           * fmask(ji,jj) / e1f(ji,jj) &
                        + (vn(ji-1,jj-1 ) - vn(ji,jj-1) )   &
                &           * fmask(ji-1,jj-1) / e1f(ji-1,jj-1) &
                        + (vn(ji+1,jj-1 ) - vn(ji,jj-1) )   &
                &           * fmask(ji,jj-1) / e1f(ji,jj-1) )         

           ! calcul des termes de Reynolds
           anousqrt(ji,jj) = 1000/2 * umask(ji,jj)*( ( u2n(ji,jj) - un(ji,jj)*un(ji,jj) ) &
                &                     + ( u2n(ji-1,jj) - un(ji-1,jj)*un(ji-1,jj) ) )       
                
           anovsqrt(ji,jj) = 1000/2 * vmask(ji,jj)*( ( v2n(ji,jj) - vn(ji,jj)*vn(ji,jj) ) &
                &                     + ( v2n(ji,jj-1) - vn(ji,jj)*vn(ji,jj-1) ) ) 

           anouv(ji,jj)    = 1000 * ( uvn(ji,jj) &
                &                 -   1/2 * umask(ji,jj)*( un(ji,jj) + un(ji-1,jj) ) &
                &                   * 1/2 * vmask(ji,jj)*( vn(ji,jj) + vn(ji,jj-1) ) )
           ! calcul du terme bti
           bti(ji,jj) = -1 * ( anousqrt(ji,jj) * dudx(ji,jj) &
                &           + anovsqrt(ji,jj) * dvdy(ji,jj) &
                &           + anouv(ji,jj) * ( dvdx(ji,jj) + dudy(ji,jj) ))

        END DO
     END DO
     ! 
     ierr = putvar(ncout, id_varout(1) ,dudx, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(2) ,dvdx, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(3) ,dudy, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(4) ,dvdy, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(5) ,anousqrt, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(6) ,anovsqrt, jk, npiglo, npjglo, ktime=1) 
     ierr = putvar(ncout, id_varout(7) ,anouv, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(8) ,bti, jk, npiglo, npjglo, ktime=1)
  END DO
  ierr = closeout(ncout)

END PROGRAM cdfbti

