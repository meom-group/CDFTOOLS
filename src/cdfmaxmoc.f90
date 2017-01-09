PROGRAM cdfmaxmoc
  !!======================================================================
  !!                     ***  PROGRAM  cdfmaxmoc  ***
  !!=====================================================================
  !!  ** Purpose : Compute the maximum of the overturning fonction from 
  !!               a file calculated by cdfmoc
  !!
  !!  ** Method  : A spatial window, limited by latmin latmax depmin depmax
  !!               given on the command line, is used to determnine the 
  !!               maximum and minimum of the MOC as well as their 
  !!               respective depth and latitude. 
  !!
  !! History : 2.1  : 07/2005  : J.M. Molines : Original code
  !!                : 11/2009  : R. Dussin    : Netcdf output
  !!           3.0  : 03/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  !
  INTEGER(KIND=4)                             :: jj, jk            ! dummy loop index
  INTEGER(KIND=4)                             :: npjglo, npk       ! size of the overturning
  INTEGER(KIND=4)                             :: narg, iargc       ! line command stuff
  INTEGER(KIND=4)                             :: iarg              ! line command stuff
  INTEGER(KIND=4)                             :: ijmin, ijmax      ! latitude window where to look at extrema
  INTEGER(KIND=4)                             :: ikmin, ikmax      ! depth window where to look at extrema
  INTEGER(KIND=4)                             :: ilatmin, ilatmax  ! index of found extrema (latitude)
  INTEGER(KIND=4)                             :: idepmin, idepmax  ! index of found extrema (depth )
  INTEGER(KIND=4)                             :: nx=1, ny=1, nz=1  ! dims of netcdf output file
  INTEGER(KIND=4)                             :: nvarout=6         ! number of values to write in cdf output
  INTEGER(KIND=4)                             :: ncout, ierr       ! for netcdf output
  INTEGER(KIND=4), DIMENSION(3)               :: iminloc, imaxloc  ! work arrays for minloc and maxloc
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE  :: ipk, id_varout    ! netcdf output
  !
  REAL(KIND=4)                                :: ovtmax, ovtmin    ! max/ min of MOC ( Sv)
  REAL(KIND=4)                                :: rlatmin, rlatmax  ! latitude limits for searching
  REAL(KIND=4)                                :: rdepmin, rdepmax  ! depth limits for searching
  REAL(KIND=4), DIMENSION(1)                  :: tim               ! time counter
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: gdepw             ! depth read in the header
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rdumlon, rdumlat  ! dummy array for output
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rlat              ! latitude (1, npjglo)
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: rmoc              ! MOC (1, npjglo, jpk)
  !
  TYPE(variable), DIMENSION(:),   ALLOCATABLE :: stypvar           ! structure of output
  !
  CHARACTER(LEN=256)                          :: cf_moc            ! input file 
  CHARACTER(LEN=256)                          :: cf_ncout='maxmoc.nc' ! output file
  CHARACTER(LEN=256)                          :: cldum             ! dummy string for I/O
  CHARACTER(LEN=256)                          :: cbasin, cv_in     ! basin name and cdf variable name
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  narg=iargc()
  
  IF ( narg /= 6 ) THEN
     PRINT *,' usage : cdfmaxmoc OVT-file basin_name latmin latmax depmin depmax'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute the maximum and minimum of the overturning, from file OVT-file,' 
     PRINT *,'        for oceanic basin specified by cbasin, and in the geographical frame '
     PRINT *,'        defined by latmin latmax, depmin, depmax.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       OVT-file   : overturning file from cdfmoc, with or w/o sub basins.' 
     PRINT *,'       basin_name : name of oceanic subbasin as defined in ',TRIM(cn_fbasins)
     PRINT *,'                usually it can be one of atl, glo, inp, ind or pac'
     PRINT *,'                glo means no subbasins.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_ncout) 
     PRINT *,'         6 variables : '
     PRINT *,'            maxmoc, minmoc ( sv )      : max and min of overturning'
     PRINT *,'            latmaxmoc latminmoc ( deg) : latitudes of max and min.'
     PRINT *,'            depmaxmoc depminmoc ( m)   : depth of max amd min .'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfmoc '
     PRINT *,'      '
     STOP
  ENDIF

  CALL getarg(1, cf_moc)                          ! input moc file
  CALL getarg(2, cbasin)                          ! basin name
  CALL getarg(3, cldum ) ; READ(cldum,*) rlatmin  ! searching window : latmin
  CALL getarg(4, cldum ) ; READ(cldum,*) rlatmax  ! searching window : latmax
  CALL getarg(5, cldum ) ; READ(cldum,*) rdepmin  ! searching window : depth min
  CALL getarg(6, cldum ) ; READ(cldum,*) rdepmax  ! searching window : depth max
  
  IF ( chkfile(cf_moc) ) STOP ! missing file

  npjglo = getdim(cf_moc, cn_y)
  npk    = getdim(cf_moc, cn_z)

  ALLOCATE ( rmoc(1,npjglo,npk), gdepw(npk), rlat(1,npjglo))
  gdepw(:)  = -getvar1d(cf_moc, cn_vdepthw,  npk         )
  rlat(:,:) =  getvar  (cf_moc, cn_vlat2d,   1, 1, npjglo)

  SELECT CASE (cbasin)
  CASE ('atl') ; cv_in=cn_zomsfatl
  CASE ('glo') ; cv_in=cn_zomsfglo
  CASE ('pac') ; cv_in=cn_zomsfpac
  CASE ('inp') ; cv_in=cn_zomsfinp
  CASE ('ind') ; cv_in=cn_zomsfind
  CASE DEFAULT ; STOP 'basin not found'
  END SELECT

  ALLOCATE ( stypvar(nvarout), ipk(nvarout), id_varout(nvarout) )
  ALLOCATE ( rdumlon(1,1) , rdumlat(1,1) )

  rdumlon(:,:)=0.
  rdumlat(:,:)=0.

  DO jj=1,nvarout
     ipk(jj)=1
  ENDDO

  ! define new variables for output 
  ! all variables :
  stypvar%rmissing_value    = 99999.
  stypvar%scale_factor      = 1.
  stypvar%add_offset        = 0.
  stypvar%savelog10         = 0.
  stypvar%conline_operation = 'N/A'
  stypvar%caxis             = 'T'

  ! each pair of variables
  stypvar(1)%cname          = 'maxmoc'                         ; stypvar(2)%cname          = 'minmoc'
  stypvar(1)%clong_name     = 'Maximum_Overturing'             ; stypvar(2)%clong_name     = 'Minimum_Overtuning'
  stypvar(1)%cshort_name    = 'maxmoc'                         ; stypvar(2)%cshort_name    = 'minmoc'
  stypvar(1:2)%cunits       = 'Sverdrup'
  stypvar(1:2)%valid_min    = -1000. 
  stypvar(1:2)%valid_max    =  1000.

  stypvar(3)%cname          = 'latmaxmoc'                      ; stypvar(4)%cname          = 'latminmoc'
  stypvar(3)%clong_name     = 'Latitude_of_Maximum_Overturing' ; stypvar(4)%clong_name     = 'Latitude_of_Minimum_Overtuning'
  stypvar(3)%cshort_name    = 'latmaxmoc'                      ; stypvar(4)%cshort_name    = 'latminmoc'
  stypvar(3:4)%cunits       = 'Degrees'
  stypvar(3:4)%valid_min    = -90.
  stypvar(3:4)%valid_max    = 90.

  stypvar(5)%cname          = 'depthmaxmoc'                    ; stypvar(6)%cname          = 'depthminmoc'
  stypvar(5)%clong_name     = 'Depth_of_Maximum_Overturing'    ; stypvar(6)%clong_name     = 'Depth_of_Minimum_Overtuning'
  stypvar(5)%cshort_name    = 'depthmaxmoc'                    ; stypvar(6)%cshort_name    = 'depthminmoc'
  stypvar(5:6)%cunits       = 'Meters'
  stypvar(5:6)%valid_min    = -10000.
  stypvar(5:6)%valid_max    = 0.

  DO jk=1,npk
     rmoc(:,:,jk) = getvar(cf_moc, cv_in, jk, 1, npjglo)
  END DO

  ! define window in index limit
  ! look for ijmin-ijmax :
  DO jj=1, npjglo
     IF ( rlat(1,jj) <= rlatmin )  ijmin = jj
     IF ( rlat(1,jj) <= rlatmax )  ijmax = jj
  END DO

  ! look for ikmin ikmax
  DO jk=1,npk
     IF ( gdepw(jk) <= rdepmin ) ikmin = jk
     IF ( gdepw(jk) <= rdepmax ) ikmax = jk
  END DO

  ! look for max/min overturning
  ovtmax = MAXVAL(rmoc(1,ijmin:ijmax,ikmin:ikmax))
  ovtmin = MINVAL(rmoc(1,ijmin:ijmax,ikmin:ikmax))

  ! find location of min/max
  iminloc =MINLOC(rmoc(:,ijmin:ijmax,ikmin:ikmax))
  imaxloc =MAXLOC(rmoc(:,ijmin:ijmax,ikmin:ikmax))

  ! results from minloc/maxloc is relative to the sub -array given as arguments
  ilatmin = iminloc(2) + ijmin -1 ; ilatmax = imaxloc(2) + ijmin -1
  idepmin = iminloc(3) + ikmin -1 ; idepmax = imaxloc(3) + ikmin -1

  PRINT *,' Maximum ', ovtmax ,' Sv latitude ', rlat(1,ilatmax),' depth = ', gdepw(idepmax)
  PRINT *,' Minimum ', ovtmin ,' Sv latitude ', rlat(1,ilatmin),' depth = ', gdepw(idepmin)

  ! create output fileset
  ncout = create      (cf_ncout, 'none',  nx,      ny,  nz,      cdep=cn_vdepthw )
  ierr  = createvar   (ncout,    stypvar, nvarout, ipk, id_varout                )

  ierr  = putheadervar(ncout,    cf_moc,  nx,      ny,  nz,      &
                        pnavlon=rdumlon, pnavlat=rdumlat,   pdep=gdepw           )
               
  tim  = getvar1d(cf_moc,cn_vtimec, 1     )
  ierr = putvar1d(ncout, tim,       1, 'T')

  ! netcdf output 
  ierr = putvar0d(ncout,id_varout(1), REAL(ovtmax) )
  ierr = putvar0d(ncout,id_varout(2), REAL(ovtmin) )
  ierr = putvar0d(ncout,id_varout(3), REAL(rlat(1,ilatmax)) )
  ierr = putvar0d(ncout,id_varout(4), REAL(rlat(1,ilatmin)) )
  ierr = putvar0d(ncout,id_varout(5), REAL(gdepw(idepmax)) )
  ierr = putvar0d(ncout,id_varout(6), REAL(gdepw(idepmin)) )

  ierr = closeout(ncout)

END PROGRAM cdfmaxmoc
