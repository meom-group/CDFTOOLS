PROGRAM cdfzonalmean
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfzonalmean  ***
  !!
  !!  **  Purpose  :  Compute the zonal mean 
  !!  
  !!  **  Method   :  
  !!                  Results are saved on zonalmean.nc file with 
  !!                  variables name respectively as follow:
  !!                  same as input except that the 2 first char are
  !!                  changed to zo. Then a suffix is append to the
  !!                  name of the variable : glo atl inp ind and pac
  !!                  if a subbasin mask is given on input., else
  !!                  the suffix glo is used. Example :
  !!                  sosaline_glo sosaline_atl etc ...
  !!
  !!
  !! history ;
  !!  Original :  J.M. Molines (nov. 2005) 
  !!  Modified :  P.   Mathiot (June 2007) Update for forcing fields format
  !!                                       + for many time steps
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: npbasins=1, ivar = 0                !: number of subbasin, number of output var
  INTEGER   :: jbasin, jj, jk ,ji ,jvar ,jjvar,jkk !: dummy loop index
  INTEGER   :: jt                                  !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk, nt              !: size of the domain
  INTEGER   :: ncout
  INTEGER   :: nvars , mvar                        !: number of variables in the file
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipk, ijvar, ipko, id_varout    !: jpbasin x nvar
  INTEGER, DIMENSION(2)              ::  iloc

  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  e1, e2, gphi, zv !:  metrics, velocity
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlon              !: dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlat              !: latitude for i = north pole
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  zmaskvar
  REAL(KIND=4), DIMENSION (:),       ALLOCATABLE ::  gdep                !: gdept or gdepw
  REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  zmask               !:  jpbasins x npiglo x npjglo
  REAL(KIND=4), DIMENSION (1)                    ::  tim
  REAL(KIND=4)                                   :: spval=99999.

  REAL(KIND=8), DIMENSION (:,:),   ALLOCATABLE ::  zomsf , area        !: jpbasins x npjglo x npk

  CHARACTER(LEN=80) :: cfilev , cfileoutnc='zonalmean.nc', cdum
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc',cmaskfil='mask.nc',cbasinmask='none'
  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE   :: cvarname             !: array of var name for input
  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE   :: cvarnameo             !: array of var name for output
  TYPE(variable), DIMENSION(:), ALLOCATABLE   :: typvar                  !: structure for attributes
  TYPE(variable), DIMENSION(:), ALLOCATABLE   :: typvaro                 !: structure for attributes
  CHARACTER(LEN=10) :: ce1, ce2, cphi, cdep,cmask, cdepo
  CHARACTER(LEN=4),DIMENSION(5) :: cbasin=(/'_glo','_atl','_inp','_ind','_pac'/)

  LOGICAL :: lrevert_dep = .TRUE.         !: flag to revert the order of depth in the output file (plotting facility)
  LOGICAL :: lforcing    = .FALSE.

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfzonalmean  file  T | U | V | F | W [new_maskglo.nc]'
     PRINT *,' Computes the zonal mean '
     PRINT *,' If no new_maskglo specified, assume global '
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc '
     PRINT *,'  must be in the current directory'
     PRINT *,' Output on zonalmean.nc: '
     PRINT *,'      variables zoxxxx_glo  : Global ocean '
     PRINT *,'      variables zoxxxx_atl  : Atlantic Ocean '
     PRINT *,'      variables zoxxxx_inp  : Indo Pacific '
     PRINT *,'      variables zoxxxx_ind  : Indian Ocean alone'
     PRINT *,'      variables zoxxxx_pac  : Pacific Ocean alone'
     PRINT *,'  Depth variable output is negative (standard) unless '
     PRINT *,'  you recompile the tool with lrevert_dep=.false.'
     STOP
  ENDIF

  CALL getarg (1, cfilev)
  CALL getarg (2, cdum )

  ! set the metrics according to C grid point
  SELECT CASE (cdum)
  CASE ('T', 't', 'S', 's')
     ce1='e1t'
     ce2='e2t'
     cdep='gdept'
     cdepo='deptht'
     cphi='gphit'
     cmask='tmask'
  CASE ('U', 'u')
     ce1='e1u'
     ce2='e2u'
     cdep='gdepu'
     cdepo='depthu'
     cphi='gphiu'
     cmask='umask'
  CASE ('V', 'v')
     ce1='e1v'
     ce2='e2v'
     cdep='gdepv'
     cdepo='depthv'
     cphi='gphiv'
     cmask='vmask'
  CASE ('F', 'f')
     ce1='e1f'
     ce2='e2f'
     cdep='gdepf'
     cdepo='deptht'
     cphi='gphif'
     cmask='fmask'
  CASE ('W', 'w')
     ce1='e1t'
     ce2='e2t'
     cdep='gdepw'
     cdepo='depthw'
     cphi='gphit'
     cmask='tmask'
  CASE DEFAULT
     PRINT *, ' C grid:', TRIM(cdum),' point not known!'
     STOP
  END SELECT

  ! Read sub_basin file name (optional)
  IF (narg == 3 ) THEN
     CALL getarg(3, cbasinmask)
     npbasins=5
  ENDIF

  nvars  = getnvar(cfilev)
  ALLOCATE ( cvarname(nvars)          ,ipk(nvars), ijvar(nvars), typvar(nvars)  )
  ALLOCATE ( cvarnameo(npbasins*nvars),ipko(npbasins*nvars),id_varout(npbasins*nvars) )
  ALLOCATE ( typvaro(npbasins*nvars))
  cvarname(1:nvars) = getvarname(cfilev,nvars,typvar)
  ipk(1:nvars) = getipk(cfilev,nvars)

  ! buildt output filename
  ivar = 0
  mvar = 0
  DO jvar = 1,nvars
     ! skip variables such as nav_lon, nav_lat, time_counter deptht ...
     IF (ipk(jvar) == 0 ) THEN
        cvarname(jvar)='none'
     ELSE
        mvar = mvar + 1       ! count for valid input variables
        ijvar(mvar) = jvar    ! use indirect adressing for those variables
        DO jbasin=1,npbasins
           ivar=ivar + 1      ! count for output variables
           cvarnameo(ivar)='zo'//TRIM(cvarname(jvar)(3:))//TRIM(cbasin(jbasin) )
           ! intercept case of duplicate zonal name
           IF (cvarname(jvar) == 'iowaflup' ) cvarnameo(ivar)='zowaflio'//TRIM(cbasin(jbasin) )
           IF (cvarname(jvar) == 'cfc11' ) cvarnameo(ivar)='zocfc11'//TRIM(cbasin(jbasin) )
           IF (cvarname(jvar) == 'bombc14' ) cvarnameo(ivar)='zobc14'//TRIM(cbasin(jbasin) )
           IF (cvarname(jvar) == 'invcfc' ) cvarnameo(ivar)='zoinvcfc'//TRIM(cbasin(jbasin) )
           IF (cvarname(jvar) == 'invc14' ) cvarnameo(ivar)='zoinvc14'//TRIM(cbasin(jbasin) )
           IF (cvarname(jvar) == 'qtrcfc' ) cvarnameo(ivar)='zoqtrcfc'//TRIM(cbasin(jbasin) )
           IF (cvarname(jvar) == 'qtrc14' ) cvarnameo(ivar)='zoqtrc14'//TRIM(cbasin(jbasin) )
           IF (cvarname(jvar) == 'qintcfc' ) cvarnameo(ivar)='zoqintcfc'//TRIM(cbasin(jbasin) )
           IF (cvarname(jvar) == 'qintc14' ) cvarnameo(ivar)='zoqintc14'//TRIM(cbasin(jbasin) )
           typvaro(ivar)%name=cvarnameo(ivar)
           ! units can be build automatically ( the same as original variable)
           typvaro(ivar)%units=typvar(jvar)%units
           ! missing value, valid min and valid max : idem original field
           typvaro(ivar)%missing_value=spval
           typvaro(ivar)%valid_min=typvar(jvar)%valid_min
           typvaro(ivar)%valid_max=typvar(jvar)%valid_max
           ! longname : prefix=Zonal_Mean_   suffix=TRIM(cbasin(jbasin)
           typvaro(ivar)%long_name='Zonal_Mean_'//TRIM(typvar(jvar)%long_name)//TRIM(cbasin(jbasin) )
           ! shortname=name
           typvaro(ivar)%short_name=typvaro(ivar)%name
           ! online operation : N/A (as usual ...)
           typvaro(ivar)%online_operation='/N/A'
           ! axis : either TY( original 2D)  or TZY (original 3D)
           IF (ipk(jvar) == 1 ) THEN
             typvaro(ivar)%axis='TY'
           ELSE
             typvaro(ivar)%axis='TZY'
           ENDIF
     
           ipko(ivar)=ipk(jvar)
        END DO
     ENDIF
  END DO

  npiglo= getdim (cfilev,'x')
  npjglo= getdim (cfilev,'y')
  npk   = getdim (cfilev,'depth')
  nt    = getdim (cfilev,'time_counter')

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'nt    =', nt
  ! if forcing fields, npk=0, assume 1
  IF (npk==0) THEN
     npk = 1 
     lforcing= .TRUE.
     PRINT *,' It is a forcing field, assume npk=1 and gdep=0'
  END IF
  ! Allocate arrays
  ALLOCATE ( zmask(npbasins,npiglo,npjglo) )
  ALLOCATE ( zv(npiglo,npjglo) )
  ALLOCATE ( zmaskvar(npiglo,npjglo) )
  ALLOCATE ( e1(npiglo,npjglo),e2(npiglo,npjglo), gphi(npiglo,npjglo) ,gdep(npk) )
  ALLOCATE ( zomsf( npjglo, npk) ,area( npjglo, npk) )
  ALLOCATE ( dumlon(1,npjglo) , dumlat(1,npjglo))


  ! get the metrics
  e1(:,:)   = getvar(coordhgr, ce1, 1,npiglo,npjglo) 
  e2(:,:)   = getvar(coordhgr, ce2, 1,npiglo,npjglo) 
  gphi(:,:) = getvar(coordhgr, cphi, 1,npiglo,npjglo)
  IF (.NOT. lforcing) gdep(:)   = getvare3(coordzgr, cdep ,npk)
  IF (lforcing)  gdep(:)   = 0
  IF ( lrevert_dep ) gdep(:)   = -1.*  gdep(:)     ! helps for plotting the results

  ! Look for the i-index that go through the North Pole
  iloc        = MAXLOC(gphi)
  dumlat(1,:) = gphi(iloc(1),:)
  dumlon(:,:) = 0.   ! set the dummy longitude to 0

  ! create output fileset
  ncout = create(cfileoutnc, cfilev, 1,npjglo,npk,cdep=cdepo)
  ierr  = createvar(ncout ,typvaro,ivar, ipko,id_varout )
  ierr  = putheadervar(ncout, cfilev,1,npjglo,npk,pnavlon=dumlon,pnavlat=dumlat,pdep=gdep)
  tim   = getvar1d(cfilev,'time_counter',1)
  ierr  = putvar1d(ncout,tim,1,'T')

  ! reading the surface mask masks
  ! 1 : global ; 2 : Atlantic ; 3 : Indo-Pacif ; 4 : Indian ; 5 : Pacif
  zmask(1,:,:) = getvar(cmaskfil,cmask,1,npiglo,npjglo)
  IF ( cbasinmask /= 'none' ) THEN
     zmask(2,:,:) = getvar(cbasinmask,'tmaskatl',1,npiglo,npjglo)
     zmask(4,:,:) = getvar(cbasinmask,'tmaskind',1,npiglo,npjglo)
     zmask(5,:,:) = getvar(cbasinmask,'tmaskpac',1,npiglo,npjglo)
     zmask(3,:,:) = zmask(5,:,:)+zmask(4,:,:)
     ! ensure that there are no overlapping on the masks
     WHERE(zmask(3,:,:) > 0 ) zmask(3,:,:) = 1
  ENDIF

  ! main computing loop
  ivar = 0
  DO jjvar = 1, mvar
     jvar = ijvar(jjvar)
     DO jt = 1,nt
        IF (MOD(jt,100)==0) PRINT *, jt,'/',nt
        DO jkk = 1, ipk(jvar)
           PRINT *,TRIM(cvarname(jvar)), ' level ',jkk
           ! Get variables and mask at level jk
           IF (.NOT. lforcing) jk=jkk
           IF ( lforcing     ) jk=jt
           zv(:,:)       = getvar(cfilev,   cvarname(jvar),jk ,npiglo,npjglo)
           zmaskvar(:,:) = getvar(cmaskfil, cmask,           jkk ,npiglo,npjglo)
           
           ! For all basins 
           DO jbasin = 1, npbasins
              zomsf(:,:) = 0.d0
              area(:,:) = 0.d0
              ! integrates 'zonally' (along i-coordinate)
              DO ji=1,npiglo
                 DO jj=1,npjglo
                    zomsf(jj,jkk) = zomsf(jj,jkk) + e1(ji,jj)*e2(ji,jj)* zmask(jbasin,ji,jj)*zmaskvar(ji,jj)*zv(ji,jj)
                    area(jj,jkk)  =  area(jj,jkk) + e1(ji,jj)*e2(ji,jj)* zmask(jbasin,ji,jj)*zmaskvar(ji,jj)
                 END DO
              END DO

              ! compute the mean value if the area is not 0, else assign spval
              WHERE (area /= 0 ) 
                 zomsf=zomsf/area
              ELSEWHERE
                 zomsf=spval
              ENDWHERE
              ivar=  (jjvar-1)*npbasins + jbasin
              ierr = putvar (ncout, id_varout(ivar),REAL(zomsf(:,jkk)), jk,1,npjglo, ktime=jt)
           END DO  !next basin
        END DO  ! next k 
     END DO ! next time
  END DO ! next variable

  ierr = closeout(ncout)

END PROGRAM cdfzonalmean
