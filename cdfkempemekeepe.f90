PROGRAM cdfkempemekeepe
  !!---------------------------------------------------------------------------
  !!         ***  PROGRAM  cdfkempemekeepe  ***
  !!
  !!  **  Purpose: Compute the term of energetic transfert  
  !!     from mean kinetic energy to mean potential energy (T1) 
  !!     and from eddy potential energy to eddy kinetic energy (T3)
  !!
  !! history :
  !!   Original :  A. Melet (Mar 2008)
  !!---------------------------------------------------------------------
  !!--------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used

  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER :: ji,jj,jk, jt, ilev, jmin
  INTEGER :: npiglo, npjglo, npk, nt
  INTEGER :: kimin,imin,kkmin,ktime
  INTEGER :: narg, iargc, ncout, ierr
  INTEGER, DIMENSION(2) ::  ipk, id_varout         ! 

  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: t1mask,w1mask,txz,wxz 
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: wtxz,wbartbarxz,anowtxz
  REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE  :: wbartbar,anowt 
  REAL(KIND=4) ,DIMENSION(1)                 :: tim

  CHARACTER(LEN=80) :: cfile
  CHARACTER(LEN=80) :: cfileout='transfertst1t3.nc'
  TYPE (variable), DIMENSION(2) :: typvar         !: structure for attibutes

  !!
  narg = iargc()
  IF ( narg /= 1 ) THEN
     PRINT *,' USAGE : cdfkempemekeepe file'
     PRINT *,'        Produce a cdf file transfertst1t3.nc with wT and anowT variables'
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
  typvar(1)%name='wT'
  typvar(1)%long_name='temporal mean of w times temporal mean of T on T point (*1000)'
  typvar(1)%short_name='wT'

  typvar(2)%name='anowT'
  typvar(2)%long_name='temporal mean of anomaly of w times ano of T on T point (*1000)'
  typvar(2)%short_name='dTdy'
 
  typvar%units='1000 m.K'
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
  ierr= createvar(ncout ,typvar,2, ipk,id_varout )
  ierr= putheadervar(ncout, cfile, npiglo, npjglo,npk)

  ! Allocate the memory
  ALLOCATE ( t1mask(npiglo,npk) , w1mask(npiglo,npk) )
  ALLOCATE ( txz(npiglo,npk) , wxz(npiglo,npk) )
  ALLOCATE ( wtxz(npiglo,npk) )
  ALLOCATE ( wbartbarxz(npiglo,npk), anowtxz(npiglo,npk) )
  ALLOCATE ( wbartbar(npiglo, npjglo, npk) )
  ALLOCATE ( anowt(npiglo, npjglo, npk) )
  
  tim=getvar1d(cfile,'time_counter',nt)
  ierr=putvar1d(ncout,tim,1,'T')

  DO jj = 1, npjglo
     print*, 'jj : ',jj
       wbartbarxz(:,:) = 0.d0
       anowtxz(:,:)    = 0.d0
       wtxz(:,:)     = 0.d0
       wxz(:,:)      = 0.d0
       txz(:,:)      = 0.d0

       txz(:,:) = getvarxz(cfile,'tbar',jj,npiglo,npk,kimin=1,kkmin=1,ktime=1)
       wxz(:,:) = getvarxz(cfile,'wbar',jj,npiglo,npk,kimin=1,kkmin=1,ktime=1)
       wtxz(:,:)= getvarxz(cfile,'wtbar',jj,npiglo,npk,kimin=1,kkmin=1,ktime=1)

       DO jk=1, npk-1
          w1mask(:,jk) = wxz(:,jk) * wxz(:,jk+1)
          t1mask(:,jk) = txz(:,jk)
          WHERE ( w1mask(:,jk) /= 0.) w1mask(:,jk)=1.
          WHERE ( t1mask(:,jk) /= 0.) t1mask(:,jk)=1. 
          wbartbarxz(:,jk) = 1000 * t1mask(:,jk) * txz(:,jk) &
                &     * 0.5 * w1mask(:,jk) * ( wxz(:,jk) + wxz(:,jk+1) )
          anowtxz(:,jk) = 1000 * ( wtxz(:,jk) - wbartbarxz(:,jk)*0.001 )      
       END DO
       wbartbar(:,jj,:) = wbartbarxz(:,:)
       anowt(:,jj,:)    = anowtxz(:,:)
  END DO
  DO jk=1,npk
     ierr = putvar(ncout, id_varout(1) ,wbartbar(:,:,jk), jk, npiglo, npjglo )
     ierr = putvar(ncout, id_varout(2) ,anowt(:,:,jk), jk, npiglo, npjglo )
  END DO
  ierr = closeout(ncout)  

END PROGRAM cdfkempemekeepe

