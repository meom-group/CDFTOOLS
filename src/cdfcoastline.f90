! This program is not validated nor finished in the standards of CDFTOOLS
PROGRAM cdfcofpoint
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfmean  ***
  !!
  !!  **  Purpose  :  Compute distance of first coast in grid point
  !!                  Determine the edge of a mask (for further use
  !!                  with iceshelf parametrization)
  !!  
  !!  **  Method   :  long iterative method (check furtehr time all mask point)
  !!
  !!
  !! history ;
  !!  Original :  P. Mathiot (June, 2009)
  !!-------------------------------------------------------------------
  !! * Modules used
  USE cdfio
  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: npiglo, npjglo, npk, nt, npi, npj      !: size of the domain
  INTEGER   :: ji, jj, jk, jt, i
  INTEGER   :: imin=0, imax=0, jmin=0, jmax=0      !: domain limitation for computation
  INTEGER   :: kmin=0, kmax=0                      !: domain limitation for computation
  INTEGER   :: narg, iargc                         !: command line 


  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask, mask, mask_out  !: npiglo x npjglo

  CHARACTER(LEN=256) :: cfile, cldum
  CHARACTER(LEN=256) :: cmask='mask.nc'

  ! output stuff
  INTEGER, DIMENSION(1) :: ipk, id_varout
  TYPE(variable), DIMENSION(1) :: stypvar
  REAL(KIND=4) ,DIMENSION(1)                  :: timean
  CHARACTER(LEN=256) :: cfileout='pointcoast.nc'
  INTEGER :: ncout, ierr

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmean  filemask [imin imax jmin jmax kmin kmax] '
     PRINT *,' imin imax jmin jmax kmin kmax can be given in option '
     PRINT *,'    if imin = 0 then ALL i are taken'
     PRINT *,'    if jmin = 0 then ALL j are taken'
     PRINT *,'    if kmin = 0 then ALL k are taken'
     PRINT *,'    output file is pointcoast.nc    '
     STOP
  ENDIF

  CALL getarg (1, cfile)

  npi= getdim (cfile,'x')
  npj= getdim (cfile,'y')
  npk   = getdim (cfile,'depth')
  nt    = getdim (cfile,'time')
  
  imin=1; jmin=1; jmax=npj; imax=npi

  IF (narg > 3 ) THEN
    IF ( narg /= 5 ) THEN
       PRINT *, ' ERROR : You must give 6 optional values (imin imax jmin jmax kmin kmax)'
       STOP
    ELSE
    ! input optional imin imax jmin jmax
      CALL getarg ( 2,cldum) ; READ(cldum,*) imin
      CALL getarg ( 3,cldum) ; READ(cldum,*) imax
      CALL getarg ( 4,cldum) ; READ(cldum,*) jmin
      CALL getarg ( 5,cldum) ; READ(cldum,*) jmax
    ENDIF
  ENDIF

  IF (npk == 0  ) THEN ; npk = 1              ; ENDIF  ! no depth dimension ==> 1 level
  IF (imin /= 0 ) THEN ; npiglo=imax -imin + 1;  ELSE ; imin=1 ; ENDIF
  IF (jmin /= 0 ) THEN ; npjglo=jmax -jmin + 1;  ELSE ; jmin=1 ; ENDIF
  IF (kmin /= 0 ) THEN ; npk   =kmax -kmin + 1;  ELSE ; kmin=1 ; ENDIF

  WRITE(6, *) 'npiglo=', npiglo
  WRITE(6, *) 'npjglo=', npjglo
  WRITE(6, *) 'npi   =', npi
  WRITE(6, *) 'npj   =', npj
  WRITE (6,*) 'npk   =', npk
  WRITE (6,*) 'nt    =', nt
  
  ! Allocate arrays
  ALLOCATE ( zmask(npi,npj), mask_out(npi,npj), mask(npi,npj) )
  mask_out(:,:)=0
  mask(:,:)= getvar(cfile, 'tmask',  1 ,npi,npj)
  i=0
  DO WHILE ( SUM(mask(imin+1:imax-1,jmin+1:jmax-1)) .NE. 0 )
     i=i+1
     zmask=0
     IF (MOD(i,10)==0) PRINT *, 'i = ',i, ' SUM(mask) = ',SUM(mask(imin+1:imax-1,jmin+1:jmax-1)) 
     DO ji=imin+1,imax-1
        DO jj=jmin+1,jmax-1
           IF ((mask(ji,jj)==0) .AND. (mask_out(ji,jj)==i-1)) THEN
              IF (mask(ji+1,jj  )==1) zmask(ji+1,jj  )=1
              IF (mask(ji  ,jj-1)==1) zmask(ji  ,jj-1)=1
              IF (mask(ji-1,jj  )==1) zmask(ji-1,jj  )=1
              IF (mask(ji  ,jj+1)==1) zmask(ji  ,jj+1)=1
           END IF
        END DO
     END DO
     WHERE (zmask==1) mask=0
     WHERE (zmask==1) mask_out=i
  END DO
        
  ! prepare file output
  ipk(1) = 1
  stypvar(1)%cname='pointcoast'
  stypvar(1)%cunits='px'
  stypvar(1)%rmissing_value=0
  stypvar(1)%valid_min= 1.
  stypvar(1)%valid_max= i
  stypvar(1)%clong_name='pointcoast'
  stypvar(1)%cshort_name='pointcoast'
  stypvar(1)%conline_operation='N/A'
  stypvar(1)%caxis='TZYX'
  stypvar(1)%cprecision='r4'
  PRINT *,' CREATE ...'
  ncout=create(cfileout, cfile,npi,npj,1)

  PRINT *,' CREATEVAR ...'
  ierr= createvar(ncout ,stypvar,1, ipk,id_varout )
  PRINT *,' PUTHEADERVAR ...'
  ierr= putheadervar(ncout, cfile, npi,npj,npk)
  ierr=putvar(ncout,id_varout(1),mask_out,1,npi,npj)
  timean(:)=0.
  ierr=putvar1d(ncout,timean,1,'T')
  ierr = closeout(ncout)

END PROGRAM cdfcofpoint
