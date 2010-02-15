PROGRAM cdfnrjcomp
  !!---------------------------------------------------------------------------
  !!         ***  PROGRAM  cdfnrjcomp  ***
  !!
  !!  **  Purpose: Compute the terms for energy components 
  !!               (Mean Kinetic Energy, Eddy Kinetic Energy,
  !!                Mean Potential Energy, Eddy Potential Energy )
  !!               compute : tbar,ubar,vbar,anotsqrt,anousqrt,anovsqrt
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
  INTEGER, DIMENSION(6) ::  ipk, id_varout         ! 

  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: un, vn, u2n, v2n
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: tn, t2n, anotsqrt
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: umask, vmask
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: anousqrt, anovsqrt 
  REAL(KIND=4) ,DIMENSION(1)                 :: tim

  CHARACTER(LEN=256) :: cfile
  CHARACTER(LEN=256) :: cfileout='nrjcomp.nc'
  TYPE (variable), DIMENSION(6) :: typvar         !: structure for attibutes

  !!
  narg = iargc()
  IF ( narg /= 1 ) THEN
     PRINT *,' USAGE : cdfnrjcomp file'
     PRINT *,'        Produce a cdf file nrjcomp.nc with variables'
     PRINT *,'        tbar,ubar,vbar,anotsqrt,anousqrt,anovsqrt on T point'
     PRINT *,'        file is from cdfmoyuvwt'
     PRINT *,'        the mean must have been computed on a period long enough'
     PRINT *,'        for the statistics to be meaningful'
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
  typvar(1)%name='tbar'
  typvar(1)%long_name='temporal mean of the temperature on T point'
  typvar(1)%short_name='tbar'

  typvar(2)%name='ubar'
  typvar(2)%long_name='temporal mean of the zonal velocity on T point'
  typvar(2)%short_name='ubar'
 
  typvar(3)%name='vbar'
  typvar(3)%long_name='temporal mean of the meridional velocity on T point'
  typvar(3)%short_name='vbar'
  
  typvar(4)%name='anotsqrt'
  typvar(4)%long_name='temporal mean of the square of the temperature anomaly on T point (*1000)'
  typvar(4)%short_name='anotsqrt'

  typvar(5)%name='anousqrt'
  typvar(5)%long_name='temporal mean of the square of the zonal speed anomaly on T point (*1000)'
  typvar(5)%short_name='anousqrt'

  typvar(6)%name='anovsqrt'
  typvar(6)%long_name='temporal mean of the square of the meridional speed anomaly on T point (*1000)'
  typvar(6)%short_name='anovsqrt'


  typvar%units=' '
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
  ierr= createvar(ncout ,typvar,6, ipk,id_varout )
  ierr= putheadervar(ncout, cfile, npiglo, npjglo,npk)

  ! Allocate the memory
  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo)  )
  ALLOCATE ( umask(npiglo,npjglo) , vmask(npiglo,npjglo) )
  ALLOCATE ( u2n(npiglo,npjglo)  , v2n(npiglo,npjglo)  )
  ALLOCATE ( anousqrt(npiglo,npjglo) , anovsqrt(npiglo,npjglo)  )
  ALLOCATE ( tn(npiglo,npjglo)  , t2n(npiglo,npjglo)  )
  ALLOCATE ( anotsqrt(npiglo,npjglo) )
 

  tim=getvar1d(cfile,'time_counter',nt)
  ierr=putvar1d(ncout,tim,1,'T')
     
  DO jk=1, npk
     PRINT *,'            level ',jk
           
     anousqrt(:,:) = 0.d0
     anovsqrt(:,:) = 0.d0      
     anotsqrt(:,:) = 0.d0
     
     un(:,:)  =  getvar(cfile, 'ubar', jk ,npiglo,npjglo, ktime=1)
     vn(:,:)  =  getvar(cfile, 'vbar', jk ,npiglo,npjglo, ktime=1)
     u2n(:,:) =  getvar(cfile, 'u2bar', jk ,npiglo,npjglo, ktime=1)
     v2n(:,:) =  getvar(cfile, 'v2bar', jk ,npiglo,npjglo, ktime=1)
     tn(:,:)  =  getvar(cfile, 'tbar', jk ,npiglo,npjglo, ktime=1)
     t2n(:,:) =  getvar(cfile, 't2bar', jk ,npiglo,npjglo, ktime=1)

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

     DO jj = 2, npjglo  
        DO ji = 2, npiglo    ! vector opt.

           anotsqrt(ji,jj) = 1000 * ( t2n(ji,jj) - tn(ji,jj) * tn(ji,jj) )
           anousqrt(ji,jj) = 1000/2 * umask(ji,jj)*( ( u2n(ji,jj) - un(ji,jj)*un(ji,jj) ) &
                &                     + ( u2n(ji-1,jj) - un(ji-1,jj)*un(ji-1,jj) ) )       
                
           anovsqrt(ji,jj) = 1000/2 * vmask(ji,jj)*( ( v2n(ji,jj) - vn(ji,jj)*vn(ji,jj) ) &
                &                     + ( v2n(ji,jj-1) - vn(ji,jj)*vn(ji,jj-1) ) ) 

        END DO
     END DO
     ! 
     ierr = putvar(ncout, id_varout(1) ,tn, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(2) ,un, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(3) ,vn, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(4) ,anotsqrt, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(5) ,anousqrt, jk, npiglo, npjglo, ktime=1) 
     ierr = putvar(ncout, id_varout(6) ,anovsqrt, jk, npiglo, npjglo, ktime=1)
  END DO
  ierr = closeout(ncout)

END PROGRAM cdfnrjcomp

