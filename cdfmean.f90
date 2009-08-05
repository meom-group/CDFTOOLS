PROGRAM cdfmean
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfmean  ***
  !!
  !!  **  Purpose  :  Compute the Mean Value over the ocean
  !!                  PARTIAL STEPS
  !!  
  !!  **  Method   :  compute the sum ( V * e1 *e2 * e3 *mask )/ sum( e1 * e2 * e3 *mask )
  !!
  !!
  !! history ;
  !!  Original :  J.M. Molines (Oct. 2005) 
  !!              R. Dussin (Jul 2009) : add cdf output
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk, ik, jt, jj
  INTEGER   :: imin=0, imax=0, jmin=0, jmax=0      !: domain limitation for computation
  INTEGER   :: kmin=0, kmax=0                      !: domain limitation for computation
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk, nt              !: size of the domain
  INTEGER   :: nvpk                                !: vertical levels in working variable
  INTEGER   :: numout=10                           !: logical unit for output file
  ! added to write in netcdf
  INTEGER :: kx=1, ky=1                ! dims of netcdf output file
  INTEGER :: nvars=2                ! number of values to write in cdf output
  INTEGER :: ncout, ierr               ! for netcdf output
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipk, id_varout
  !
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  e1, e2, e3,  zv   !:  metrics, velocity
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask             !:   npiglo x npjglo
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE ::  gdep              !:  depth 
  ! added to write in netcdf
  REAL(KIND=4) :: threedmeanout, pmissing_value
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::  dumlon, dumlat, dummymean
  REAL(KIND=4), DIMENSION (1)               ::  tim ! time counter
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: meanout
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvar  ! structure of output
  !
  REAL(KIND=8)      :: zvol, zsum, zvol2d, zsum2d, zsurf
  CHARACTER(LEN=256) :: cfilev , cdum
  CHARACTER(LEN=256) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc',cmask='mask.nc'
  CHARACTER(LEN=256) :: cvar, cvartype
  CHARACTER(LEN=20) :: ce1, ce2, ce3, cvmask, cvtype, cdep
  CHARACTER(LEN=256) :: cfilout='out.txt'
  ! added to write in netcdf
  CHARACTER(LEN=256) :: cfileoutnc='cdfmean.nc' , cflagcdf
  CHARACTER(LEN=256) :: cdunits, cdlong_name, cdshort_name 
  ! added to write in netcdf
  LOGICAL :: lwrtcdf=.FALSE.


  INTEGER    :: istatus

  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmean  ncfile cdfvar T| U | V | F | W [imin imax jmin jmax kmin kmax] '
     PRINT *,' Computes the mean value of the field (3D, weighted) '
     PRINT *,' imin imax jmin jmax kmin kmax can be given in option '
     PRINT *,'    if imin = 0 then ALL i are taken'
     PRINT *,'    if jmin = 0 then ALL j are taken'
     PRINT *,'    if kmin = 0 then ALL k are taken'
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc '
     PRINT *,'  must be in the current directory'
     PRINT *,' Output on standard output'
     STOP
  ENDIF
  ! Open standard output with recl=256 to avoid wrapping of long lines (ifort)
  OPEN(6,FORM='FORMATTED',RECL=256)

  CALL getarg (1, cfilev)
  CALL getarg (2, cvar)
  CALL getarg (3, cvartype)

  IF (narg > 3 ) THEN
     IF ( narg < 9 .OR. narg > 10 ) THEN
        PRINT *, ' ERROR : You must give 6 optional values (imin imax jmin jmax kmin kmax)'
        STOP
     ELSE
        ! input optional imin imax jmin jmax
        CALL getarg ( 4,cdum) ; READ(cdum,*) imin
        CALL getarg ( 5,cdum) ; READ(cdum,*) imax
        CALL getarg ( 6,cdum) ; READ(cdum,*) jmin
        CALL getarg ( 7,cdum) ; READ(cdum,*) jmax
        CALL getarg ( 8,cdum) ; READ(cdum,*) kmin
        CALL getarg ( 9,cdum) ; READ(cdum,*) kmax
        IF ( narg==10 ) THEN
           CALL getarg (10,cdum) ; READ(cdum,*) cflagcdf
        ENDIF
     ENDIF
  ENDIF

  IF(cflagcdf=='cdfout') THEN
     lwrtcdf=.TRUE.
  ENDIF

  npiglo= getdim (cfilev,'x')
  npjglo= getdim (cfilev,'y')
  npk   = getdim (cfilev,'depth')
  nt    = getdim (cfilev,'time')
  nvpk  = getvdim(cfilev,cvar)
  IF (npk == 0  ) THEN ; npk = 1              ; ENDIF  ! no depth dimension ==> 1 level
  IF (imin /= 0 ) THEN ; npiglo=imax -imin + 1;  ELSE ; imin=1 ; ENDIF
  IF (jmin /= 0 ) THEN ; npjglo=jmax -jmin + 1;  ELSE ; jmin=1 ; ENDIF
  IF (kmin /= 0 ) THEN ; npk   =kmax -kmin + 1;  ELSE ; kmin=1 ; ENDIF

  IF (nvpk == 2 ) nvpk = 1
  IF (nvpk == 3 ) nvpk = npk

  WRITE(6, *) 'npiglo=', npiglo
  WRITE(6, *) 'npjglo=', npjglo
  WRITE (6,*) 'npk   =', npk
  WRITE (6,*) 'nt    =', nt
  WRITE (6,*) 'nvpk  =', nvpk

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo) )
  ALLOCATE ( zv(npiglo,npjglo) )
  ALLOCATE ( e1(npiglo,npjglo),e2(npiglo,npjglo), e3(npiglo,npjglo) )
  ALLOCATE ( gdep (npk) )
  SELECT CASE (TRIM(cvartype))
     CASE ( 'T' )
        ce1='e1t'
        ce2='e2t'
        ce3='e3t_ps'
        cvmask='tmask'
        cdep='gdept'
     CASE ( 'U' )
        ce1='e1u'
        ce2='e2u'
        ce3='e3t_ps'
        cvmask='umask'
        cdep='gdept'
     CASE ( 'V' )
        ce1='e1v'
        ce2='e2v'
        ce3='e3t_ps'
        cvmask='vmask'
        cdep='gdept'
     CASE ( 'F' )
        ce1='e1f'
        ce2='e2f'
        ce3='e3t_ps'
        cvmask='fmask'
        cdep='gdept'
     CASE ( 'W' )
        ce1='e1t'
        ce2='e2t'
        ce3='e3w_ps'
        cvmask='tmask'
        cdep='gdepw'
     CASE DEFAULT
        PRINT *, 'this type of variable is not known :', TRIM(cvartype)
        STOP
  END SELECT

  e1(:,:) = getvar(coordhgr, ce1, 1,npiglo,npjglo,kimin=imin,kjmin=jmin)
  e2(:,:) = getvar(coordhgr, ce2, 1,npiglo,npjglo,kimin=imin,kjmin=jmin)
  gdep(:) = getvare3(coordzgr,cdep,npk)

  IF(lwrtcdf) THEN
      ALLOCATE ( typvar(nvars), ipk(nvars), id_varout(nvars) )
      ALLOCATE (dumlon(kx,ky) , dumlat(kx,ky), dummymean(kx,ky) )
      ALLOCATE ( meanout(npk) )

      dumlon(:,:)=0.
      dumlat(:,:)=0.

      ipk(1)=npk ! mean for each level
      ipk(2)=1   ! 3D mean

      ierr=getvaratt (cfilev,cvar,cdunits, &
                pmissing_value, cdlong_name, cdshort_name)

      ! define new variables for output 
      typvar(1)%name='mean_'//TRIM(cvar)
      typvar%units=TRIM(cdunits)
      typvar%missing_value=99999.
      typvar%valid_min= -1000.
      typvar%valid_max= 1000.
      typvar%scale_factor= 1.
      typvar%add_offset= 0.
      typvar%savelog10= 0.
      typvar(1)%long_name='mean_'//TRIM(cdlong_name)
      typvar(1)%short_name='mean_'//TRIM(cdshort_name)
      typvar%online_operation='N/A'
      typvar%axis='ZT'

      typvar(2)%name='mean_3D'//TRIM(cvar)
      typvar(2)%long_name='mean_3D'//TRIM(cdlong_name)
      typvar(2)%short_name='mean_3D'//TRIM(cdshort_name)
      typvar%online_operation='N/A'
      typvar%axis='T'
   ENDIF


   OPEN(numout,FILE=cfilout)

   DO jt=1,nt
      zvol=0.d0
      zsum=0.d0
      DO jk = 1,nvpk
         ik = jk+kmin-1
         ! Get velocities v at ik
         zv(:,:)= getvar(cfilev, cvar,  ik ,npiglo,npjglo,kimin=imin,kjmin=jmin)
         !     zv(:,:)= getvar(cfilev, cvar,  jt ,npiglo,npjglo,kimin=imin,kjmin=jmin,ktime=jt)
         zmask(:,:)=getvar(cmask,cvmask,ik,npiglo,npjglo,kimin=imin,kjmin=jmin)
         !    zmask(:,npjglo)=0.  
         ! get e3 at level ik ( ps...)
         e3(:,:) = getvar(coordzgr, ce3, ik,npiglo,npjglo,kimin=imin,kjmin=jmin, ldiom=.TRUE.)

         ! 
         zsurf=SUM(e1 * e2 * zmask)
         zvol2d=SUM(e1 * e2 * e3 * zmask)
         zvol=zvol+zvol2d
         zsum2d=SUM(zv*e1*e2*e3*zmask)
         zsum=zsum+zsum2d
         IF (zvol2d /= 0 )THEN
            WRITE(6,*)' Mean value at level ',ik,'(',gdep(ik),' m) ',zsum2d/zvol2d, 'surface = ',zsurf/1.e6,' km^2'
            WRITE(numout,9004) gdep(ik),ik,zsum2d/zvol2d
            IF (lwrtcdf) meanout(jk)=zsum2d/zvol2d
         ELSE
            WRITE(6,*) ' No points in the water at level ',ik,'(',gdep(ik),' m) '
            IF (lwrtcdf) meanout(jk)=99999.
         ENDIF
      END DO
      WRITE(6,*) ' Mean value over the ocean: ', zsum/zvol, jt
      threedmeanout=zsum/zvol
   END DO
   CLOSE(1)
9004          FORMAT(f9.2,' ',i2,' ',f9.2)

   IF(lwrtcdf) THEN

      ! create output fileset
      ncout =create(cfileoutnc,'none',kx,ky,npk,cdep=cdep)
      ierr= createvar(ncout,typvar,nvars,ipk,id_varout )
      ierr= putheadervar(ncout, cfilev ,kx, &
           ky,npk,pnavlon=dumlon,pnavlat=dumlat,pdep=gdep,cdep=cdep)
      tim=getvar1d(cfilev,'time_counter',1)
      ierr=putvar1d(ncout,tim,1,'T')

      ! netcdf output 
      DO jk=1, nvpk
         dummymean(1,1)=meanout(jk)
         ierr = putvar(ncout, id_varout(1), dummymean, jk, kx, ky )
      END DO

      ierr=putvar0d(ncout,id_varout(2), threedmeanout )

      ierr = closeout(ncout)

   ENDIF

 END PROGRAM cdfmean
