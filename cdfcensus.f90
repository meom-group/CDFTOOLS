PROGRAM cdfcensus
  !!-----------------------------------------------------------------------------
  !!                      ***  PROGRAM cdfcensus  ***
  !!
  !!  ** Purpose:   Build an array giving the volume of water in a TS cell.
  !! 
  !!  ** Method:
  !!      T-file and S-file are scanned for a given region (eventually limited
  !!        in depth) and the volume of water in a (T,S) cell such that 
  !!          T < Tmodele < T+dt
  !!      and S < Smodele < S+ds. 
  !!    If Smodel or T model are out of the bound they are cumulated in the
  !!    nearest (T,S) cell. 
  !!      The output is done on a bimg file where S is given as
  !!    the x-direction and T the y-direction, the field value being the volume
  !!    of water. Due to a very large range in the water volume over the TS field
  !!    the field is indeed the LOG (1 + VOLUME), and even, the scale can be made
  !!    more non-linear by repeating the LOG operation, ie, for example,
  !!    field=LOG(1 + LOG (1 + VOLUME)). The parameter nlog, passed as command
  !!    argument can be used to fix the number of LOG. If nlog = 0, the true
  !!    volume is saved. 
  !!       Depending on the user purpose, limiting values tmin,tmax, and smin,smax
  !!    as well as the increments dt, ds can be adjusted. 

  !!    The ouput file is census.bimg and is always a bimg file. ---> to be changed
  !!
  !! history :
  !!     Jean-Marc MOLINES, 01/02/97 in the dynamo project for SPEM
  !!     Modifie a partir de water_mass_census_z par Anne de Miranda (27/09/99)
  !!     Rewritten in Dr. Form by Jean-Marc Molines, 11/01/02
  !!     Clothilde Langlais 01/06 CDF version and PS
  !!     J.M. Molines 03/06 : integration in CDFTOOLS-2.0
  !!     J.M. Molines 12/06 : output in netcdf, add sigma2 and sigma4 in output
  !!                          bimg output as option.
  !!-----------------------------------------------------------------------------
  !!
  ! * Module used
  USE cdfio
  USE eos

  ! * Local Variables
  IMPLICIT NONE

  INTEGER :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER :: ilog, nlog, idum
  INTEGER :: narg, iargc
  INTEGER :: it, is
  INTEGER :: ji, jj,jk
  INTEGER :: i1, i2, j1, j2, k1, k2
  INTEGER :: iarg, nt, ns

  REAL(kind=4),DIMENSION(:,:), ALLOCATABLE :: t, s, rsigma0, rsigma2, rsigma4
  REAL(kind=8),DIMENSION(:,:), ALLOCATABLE :: e1t, e2t
  REAL(kind=8),DIMENSION(:),   ALLOCATABLE :: e3t
  REAL(kind=8),DIMENSION(:,:), ALLOCATABLE :: e3t_ps
  REAL(kind=4),DIMENSION(:,:), ALLOCATABLE :: rcensus, dump, sx, ty
  REAL(kind=4),DIMENSION(:),   ALLOCATABLE :: depdum, tim

  REAL(kind=4)  :: tmin, tmax, dt, tm, xt,xs
  REAL(kind=4)  :: smin, smax, ds, sm
  REAL(kind=4)  :: tpoint, spoint, volpoint, rcmax, rcmax1

  REAL(kind=8)  ::  voltotal

  CHARACTER(LEN=80) :: cline1, cline2, cline3, cline4
  CHARACTER(LEN=80) :: cfilTS, cfildum, config
  CHARACTER(LEN=80) :: chgr='mesh_hgr.nc' , czgr='mesh_zgr.nc',cfileout='census.nc'

  INTEGER ::           ncout, ierr
  TYPE(variable), DIMENSION(4) :: typvar
  INTEGER, DIMENSION(4) ::  ipk, id_varout

  LOGICAL :: lcdf=.true. , lbimg=.false.

  ! Initialisations
  DATA tmin, tmax, dt /-2.0, 38.0, 0.05/
  DATA smin, smax, ds /25.0, 40.0, 0.02/


  voltotal=0.d0

  ! Browse command line
  narg=iargc()
!  IF ( (narg .NE. 2) .AND. (narg .NE. 7) .AND. (narg .NE. 10) ) THEN
   IF ( narg ==  0 ) THEN
     PRINT *,'>>>> usage: cdfcensus  ''TSfile'''
     PRINT *,' ''nlog'' [-zoom imin imax jmin jmax] [-klim kmin kmax] [-bimg]'
     PRINT *,' Output file is census.nc, variable volcensus,sigma0,sigma2,sigma4'
     PRINT *,' If -bimg is specified, censusopa.bimg is created instead of cdf file'
     PRINT *,' mesh_hgr and mesh_zgr.nc  must exist here ./ ' 
     STOP
  END IF

  CALL getarg(1,cfilTS)
  CALL getarg(2,cline1)
  READ(cline1,*) nlog
  PRINT *,' TS_FILE = ',TRIM(cfilTS)
  PRINT *,' NLOG  = ', nlog

  ! set domain size from TS file
  npiglo= getdim (cfilTS,'x')
  npjglo= getdim (cfilTS,'y')
  npk   = getdim (cfilTS,'deptht')
  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ! Allocate memory
  ALLOCATE (t(npiglo,npjglo),s(npiglo,npjglo))
  ALLOCATE (e1t(npiglo,npjglo),e2t(npiglo,npjglo),e3t(npk),e3t_ps(npiglo,npjglo))

  ! Read metrics
  e1t(:,:) = getvar(chgr,'e1t',1,npiglo,npjglo)
  e2t(:,:) = getvar(chgr,'e2t',1,npiglo,npjglo)
  e3t(:) = getvare3(czgr,'e3t',npk)  ! Not necessary for PS

  ! default is full domain, full depth
  i1 = 1
  i2 = npiglo
  j1 = 1
  j2 = npjglo
  k1 = 1
  k2 = npk

  ! Read additional optional argument (zoom)
  IF (narg.GE.3) THEN
     iarg = 3
     DO WHILE ( iarg .LE. narg )
        CALL getarg(iarg,cline1)
        iarg = iarg+1
        IF (cline1 .EQ. '-zoom') THEN
           CALL getarg(iarg,cline1)
           READ(cline1,*) i1
           iarg = iarg+1
           CALL getarg(iarg,cline1)
           READ(cline1,*) i2
           iarg = iarg+1
           CALL getarg(iarg,cline1)
           READ(cline1,*) j1
           iarg = iarg+1
           CALL getarg(iarg,cline1)
           READ(cline1,*) j2
           iarg = iarg+1
        ELSE IF (cline1 .EQ. '-klim' ) THEN
           CALL getarg(iarg,cline1)
           READ(cline1,*) k1
           iarg = iarg+1
           CALL getarg(iarg,cline1)
           READ(cline1,*) k2
           iarg = iarg+1
        ELSE IF (cline1 .EQ. '-bimg' ) THEN
           lbimg=.true.
           lcdf=.false.
        ELSE
           PRINT *,' Unknown option :',TRIM(cline1)
           STOP
        END IF
     END DO
  ENDIF

  ! Extra checking for over bound
  IF (i1.LT.0) i1=1
  IF (i2.GT.npiglo) i2=npiglo
  IF (j1.LT.0) j1=1
  IF (j2.GT.npjglo) j2=npjglo
  PRINT '(a,6i5)','indices:',i1,i2,j1,j2,k1,k2

  ! Compute the census on the requested domain
  PRINT *,' Water mass census on the file '
  PRINT *, TRIM(cfilTS)
  PRINT *, ' running .........'
  xt = (tmax - tmin )/dt + 1
  xs = (smax - smin )/ds + 1
  nt = NINT(xt)
  ns = NINT(xs)

  ! Allocate arrays
  ALLOCATE ( rcensus (ns,nt), dump(ns,nt) )
  ALLOCATE ( rsigma0(ns,nt), rsigma2(ns,nt), rsigma4(ns,nt) )
  ALLOCATE ( sx (ns,nt), ty(ns,nt), depdum(1) ,tim(1))
  rcensus(:,:)=0.

  ! fill up rsigma0 array with theoretical density
  DO ji=1,ns
     DO jj=1,nt
        spoint=smin+(ji-1)*ds
        tpoint=tmin+(jj-1)*dt
        sx(ji,jj)=spoint ; ty(ji,jj)=tpoint
     END DO
  END DO
  rsigma0=sigma0(ty,sx,ns,nt)
  rsigma2=sigmai(ty,sx,2000.,ns,nt)
  rsigma4=sigmai(ty,sx,4000.,ns,nt)
  depdum(1)=0.

  IF ( lcdf ) THEN
  ! create output fileset

  ipk(:)= 1                   ! Those three variables are  3D
  ! define variable name and attribute
  typvar(1)%name= 'volcensus'
  typvar(2)%name= 'sigma0'
  typvar(3)%name= 'sigma2'
  typvar(4)%name= 'sigma4'
  typvar(1)%units='m3'
  typvar(2:4)%units='kg/m3'
  typvar%missing_value=-100.
  typvar%valid_min= 0.
  typvar%valid_max= 1.e20
  typvar(1)%long_name='Volume_Census_TS'
  typvar(2)%long_name='Sigma0_TS'
  typvar(3)%long_name='Sigma2_TS'
  typvar(4)%long_name='Sigma4_TS'
  typvar(1)%short_name='volcensus'
  typvar(2)%short_name='sigma0'
  typvar(3)%short_name='sigma2'
  typvar(4)%short_name='sigma4'
  typvar%online_operation='N/A'
  typvar%axis='TYX'


  ncout =create(cfileout, cfilTS, ns,nt,1)
  ierr= createvar   (ncout ,typvar,4, ipk,id_varout )
  ierr= putheadervar(ncout, cfilTS,ns,nt,1,pnavlon=sx,pnavlat=ty,pdep=depdum)

  ENDIF

  ! Enter main loop
  DO jk=k1,k2
     t(:,:)=getvar(cfilTS, 'votemper',  jk ,npiglo, npjglo)
     s(:,:)=getvar(cfilTS, 'vosaline',  jk ,npiglo, npjglo)
     e3t_ps(:,:) = getvar('mesh_zgr.nc','e3t_ps',jk,npiglo,npjglo)

     DO ji=i1,i2
        DO jj=j1,j2
           tpoint=t(ji,jj)
           spoint=s(ji,jj)
           volpoint=e1t(ji,jj)*e2t(ji,jj)*e3t_ps(ji,jj)

           ! salinity = 0 on masked points ( OPA !!! )
           IF (spoint .NE. 0) THEN
              it=NINT( (tpoint-tmin)/dt) + 1
              is=NINT( (spoint-smin)/ds) + 1
              IF (it .LT. 1) it=1
              IF (is .LT. 1) is=1
              IF (it .GT. nt) it=nt
              IF (is .GT. ns) is=ns

              rcensus(is,it) = rcensus(is,it) + volpoint*1.e-1
              voltotal =voltotal + volpoint*1e-15
           END IF
        END DO
     END DO

  END DO  ! Main loop

  ! Computes some statistics
  rcmax=-100000.
  DO ji=1,ns
     DO jj=1,nt
        rcmax1=amax1(rcmax,rcensus(ji,jj))
        IF (rcmax1.NE.rcmax) THEN
           sm= smin+(ji-1)*ds
           tm= tmin+(jj-1)*dt
        END IF
        rcmax=rcmax1
     END DO
  END DO

  PRINT *,'  Total Volume of the domain in  10^15 m3:', REAL(voltotal)
  PRINT *,'   Volume of the most represented water mass :',rcmax
  PRINT *,'     Salinity=',sm
  PRINT *,'     Temperature=', tm

  ! use a distorsion function ( n x log ) to reduce extrema in the output file.
  DO ji=1,ns
     DO jj=1,nt
        dump(ji,jj)=rcensus(ji,jj)
        DO ilog=1,nlog
           dump(ji,jj)=ALOG10(1+dump(ji,jj))
        END DO
     END DO
  END DO

  IF ( lcdf ) THEN
  ! Output on census.nc file
  tim=getvar1d(cfilTS,'time_counter',1)
  ierr = putvar(ncout, id_varout(1) ,dump, 1,ns, nt)
  ierr = putvar(ncout, id_varout(2) ,rsigma0, 1,ns, nt)
  ierr = putvar(ncout, id_varout(3) ,rsigma2, 1,ns, nt)
  ierr = putvar(ncout, id_varout(4) ,rsigma4, 1,ns, nt)
  ierr=putvar1d(ncout,tim,1,'T')

  ierr=closeout(ncout)
  ENDIF

  IF (lbimg ) THEN
  ! Output on bimg file
  cfildum='censusopa.bimg'
  OPEN (10,file=cfildum,form='UNFORMATTED')

  WRITE(cline1,942)' Water Masses Census [10-15 m3] on',i1,i2,j1,j2
942 FORMAT(a,4i5)
  cline2='   computed from the following T-S files:'
  cline3=cfilTS
  cline4=''
  !
  WRITE(10) cline1
  WRITE(10) cline2
  WRITE(10) cline3
  WRITE(10) cline4
  WRITE(10) ns,nt,1,1,4,nlog
  WRITE(10) smin,tmin,ds,dt,0.
  WRITE(10) 0.
  WRITE(10) 0.
  WRITE(10) ((dump(ji,jj)  ,ji=1,ns),jj=1,nt)
  WRITE(10) ((rsigma0(ji,jj),ji=1,ns),jj=1,nt)
  WRITE(10) ((rsigma2(ji,jj),ji=1,ns),jj=1,nt)
  WRITE(10) ((rsigma4(ji,jj),ji=1,ns),jj=1,nt)
  CLOSE(10)
  ENDIF
  
  PRINT *,' Done.'
END PROGRAM cdfcensus
