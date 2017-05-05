PROGRAM cdftransport
  !!======================================================================
  !!                     ***  PROGRAM  cdftransport  ***
  !!=====================================================================
  !!  ** Purpose : Compute Transports across a section. 
  !!               By default, mass (Sv) and  heat(PW)/salt(kT/s) transports
  !!               are computed unless -noheat option is used (mass 
  !!               transport only).
  !!
  !!  ** Method  : The begining and end point of the section are given in 
  !!               term of F-points index. A broken line joining successive
  !!               F-points is defined between the begining and end point
  !!               of the section. Therefore each segment between F-points
  !!               is either a zonal or meridional segment corresponding to
  !!               V or U velocity component. Doing so, the volume conservation
  !!               is ensured as velocities are not interpolated, and stay 
  !!               on the native model grid. 
  !!                 The section name and the begin/end point of a section are
  !!               read from standard input, till 'EOF' is given as section
  !!               name. This make possible to give a bunch of sections in 
  !!               an ASCII files and use the < redirection.
  !!            SIGN CONVENTION : The transport is positive when the flow cross
  !!               the section to the right, negative otherwise. This depends
  !!               on the sense the section is described.  With this convention
  !!               The algebric sum of transports accross sections forming a 
  !!               closed area is 0. 
  !!            OPTIONS :
  !!               -full   : full step case
  !!               -noheat : only mass transport is computed.
  !!               -time   : specify the time frame to be used
  !!               -zlimit : transports can be computed in different depth layers
  !!                         defined by their depth limit
  !!            REQUIREMENT :
  !!               mesh-mask file are required in the current directory.
  !!            
  !!
  !! History : 2.1  : 01/2005  : J.M. Molines : Original code
  !!           2.1  : 07/2009  : R. Dussin : add cdf output
  !!           2.1  : 01/2010  : M.A. Balmaseda : Change integration signs 
  !!                             so that the transport across a segment is 
  !!                             independent of the chosen trajectory.
  !!           3.0  : 04/2011  : J.M. Molines : Doctor norm + Lic.
  !!         :  4.0  : 03/2017  : J.M. Molines
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!  interm_pt  : choose intermediate points on a broken line.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils       ! for global attribute
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)  
  !! @class transport
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                             :: jclass, jseg   ! dummy loop index
  INTEGER(KIND=4)                             :: ji, jj, jk     ! dummy loop index
  INTEGER(KIND=4)                             :: it             ! time index
  INTEGER(KIND=4), DIMENSION(:),  ALLOCATABLE :: ilev0, ilev1   ! limit in levels  (nclass)
  INTEGER(KIND=4), DIMENSION(:),  ALLOCATABLE :: ipk, id_varout ! Netcdf output
  INTEGER(KIND=4)                             :: ipos           ! working integer (position of ' ' in strings)
  INTEGER(KIND=4)                             :: ncout, ierr    ! for netcdf output
  INTEGER(KIND=4)                             :: nvarout=12     ! number of values to write in cdf output
  INTEGER(KIND=4)                             :: ivtrp          ! var index of volume transport (barotrope)
  INTEGER(KIND=4)                             :: iptrp          ! var index of volume transport (barotrope)
  INTEGER(KIND=4)                             :: imtrp          ! var index of volume transport (barotrope)
  INTEGER(KIND=4)                             :: ihtrp          ! var index of heat transport (barotrope)
  INTEGER(KIND=4)                             :: istrp          ! var index of sal transport (barotrope)
  INTEGER(KIND=4)                             :: ivtrpcl        ! var index of volume transport (p. class)
  INTEGER(KIND=4)                             :: iptrpcl        ! var index of volume transport (p. class)
  INTEGER(KIND=4)                             :: imtrpcl        ! var index of volume transport (p. class)
  INTEGER(KIND=4)                             :: ihtrpcl        ! var index of heat transport (p. class)
  INTEGER(KIND=4)                             :: istrpcl        ! var index of sal transport (p. class)
  INTEGER(KIND=4)                             :: ilonmin        ! var index of starting section longitude
  INTEGER(KIND=4)                             :: ilonmax        ! var index of ending section longitude
  INTEGER(KIND=4)                             :: ilatmin        ! var index of starting section latitude
  INTEGER(KIND=4)                             :: ilatmax        ! var index of ending section latitude
  INTEGER(KIND=4)                             :: itop           ! var index of top depth class
  INTEGER(KIND=4)                             :: ibot           ! var index of bottom depth class
  INTEGER(KIND=4)                             :: ikx=1, iky=1   ! dims of netcdf output file
  INTEGER(KIND=4)                             :: numout  = 10   ! logical unit for output file (overall)
  INTEGER(KIND=4)                             :: numvtrp = 11   ! logical unit for volume transport file
  INTEGER(KIND=4)                             :: numhtrp = 12   ! logical unit for heat transport file
  INTEGER(KIND=4)                             :: numstrp = 14   ! logical unit for salt trp file 
  INTEGER(KIND=4)                             :: nclass,ndep    ! number of depth class
  INTEGER(KIND=4)                             :: narg, iargc    ! command line 
  INTEGER(KIND=4)                             :: ijarg, nxtarg  !  "       "
  INTEGER(KIND=4)                             :: npiglo, npjglo ! size of the domain
  INTEGER(KIND=4)                             :: npk, npt       ! size of the domain
  INTEGER(KIND=4)                             :: iimin, iimax   ! i-limit of the section
  INTEGER(KIND=4)                             :: ijmin, ijmax   ! j-limit of the section
  INTEGER(KIND=4)                             :: ivar, itime    ! working integer
  INTEGER(KIND=4)                             :: ii, ij, ik     ! working integer
  INTEGER(KIND=4), PARAMETER                  :: jpseg=10000    ! used for broken line algorithm
  INTEGER(KIND=4)                             :: ii0, ij0       !  "        "             "
  INTEGER(KIND=4)                             :: ii1, ij1       !  "        "             "
  INTEGER(KIND=4)                             :: iitmp, ijtmp   !  "        "             "
  INTEGER(KIND=4)                             :: np, nn         ! segment counters, 
  INTEGER(KIND=4)                             :: iist, ijst     ! local point offset for velocity
  INTEGER(KIND=4)                             :: norm_u, norm_v ! normalization factor (sign of normal transport)
  INTEGER(KIND=4)                             :: idirx, idiry   ! sense of description of the section

  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1v, e2u       ! horizontal metric
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e3u, e3v       ! vertical metric
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: glamf          ! longitudes of F points
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: gphif          ! latitudes of F points
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zu, zut, zus   ! Zonal velocities and uT uS
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zv, zvt, zvs   ! Meridional velocities and uT uS
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zt, zs         ! temperature and salinity
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rdum           ! dummy (1x1) array for ncdf output
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zuobc, zvobc   ! arrays for OBC files (vertical slice)
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: tim            ! time counter
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: gdepw          ! depth at layer interface
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: e31d           ! vertical metric in case of full step
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: rclass         ! vertical metric in case of full step
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: rz_lst         ! limit beetween depth level, in m (nclass -1)
  REAL(KIND=4), DIMENSION(2)                  :: gla, gphi      ! lon/lat of the begining/end of section (f point)
  REAL(KIND=4), DIMENSION(jpseg)              :: rxx, ryy       ! working variables
  REAL(KIND=4)                                :: rxi0, ryj0     ! working variables
  REAL(KIND=4)                                :: rxi1, ryj1     ! working variables
  REAL(KIND=4)                                :: ai, bi         ! equation of line (y=ai.x +bi)
  REAL(KIND=4)                                :: aj, bj         ! equation of line (x=aj.y +bj
  REAL(KIND=4)                                :: rd, rd1, rd2   ! distance between point, between vertical layers
  REAL(KIND=4)                                :: udum, vdum     ! dummy velocity components for tests
  REAL(KIND=4)                                :: rau0=1000      ! density of pure water (kg/m3)
  REAL(KIND=4)                                :: rcp=4000.      ! heat capacity (J/kg/K)

  ! at every model point
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dwku,  dwkv    ! volume transport at each cell boundary
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dwkut, dwkvt   ! heat   transport at each cell boundary
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dwkus, dwkvs   ! salt   transport at each cell boundary
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dwkup, dwkvp   ! volume transport in the positive direction
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dwkum, dwkvm   !  volume transport in the negatibe direction
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dtrpu,  dtrpv  ! volume transport integrated in depth class
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dtrput, dtrpvt ! heat transport integrated in depth class
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dtrpus, dtrpvs ! salt transport integrated in depth class
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dtrpup, dtrpvp ! volume transport integrated in depth class (positive)
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dtrpum, dtrpvm ! volume transport integrated in depth class (negative)
  ! for a given section
  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: dvoltrpsum     ! volume transport by depth class across section
  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: dvoltrpsump    ! volume transport by depth class across section
  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: dvoltrpsumm    ! volume transport by depth class across section
  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: dheatrpsum     ! heat transport by depth class across section
  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: dsaltrpsum     ! salt transport by depth class across section
  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: dvolallegcl    ! over all leg volume transport by depth class
  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: dvolallegclp   ! over all leg volume transport by depth class +
  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: dvolallegclm   ! over all leg volume transport by depth class -
  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: dheatallegcl   ! over all leg heat transport by depth class 
  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: dsaltallegcl   ! over all leg salt transport by depth class 
  REAL(KIND=8), DIMENSION(jpseg)              :: dvoltrp        ! volume transport across each segment of a section
  REAL(KIND=8), DIMENSION(jpseg)              :: dvoltrpp       ! volume transport across each segment of a section
  REAL(KIND=8), DIMENSION(jpseg)              :: dvoltrpm       ! volume transport across each segment of a section
  REAL(KIND=8), DIMENSION(jpseg)              :: dheatrp        ! heat transport across each segment of a section
  REAL(KIND=8), DIMENSION(jpseg)              :: dsaltrp        ! salt transport across each segment of a section
  REAL(KIND=8)                                :: dvoltrpbrtp    ! volume transport integrated over the whole depth
  REAL(KIND=8)                                :: dvoltrpbrtpp   ! volume transport integrated over the whole depth
  REAL(KIND=8)                                :: dvoltrpbrtpm   ! volume transport integrated over the whole depth
  REAL(KIND=8)                                :: dheatrpbrtp    ! heat transport integrated over the whole depth
  REAL(KIND=8)                                :: dsaltrpbrtp    ! salt transport integrated over the whole depth
  REAL(KIND=8)                                :: dvolalleg      ! over all leg sum of volume transport
  REAL(KIND=8)                                :: dvolallegp     ! over all leg sum of volume transport +
  REAL(KIND=8)                                :: dvolallegm     ! over all leg sum of volume transport -
  REAL(KIND=8)                                :: dheatalleg     ! over all leg sum of heat transport 
  REAL(KIND=8)                                :: dsaltalleg     ! over all leg sum of salt transport 

  COMPLEX, DIMENSION(jpseg)                   :: yypt           ! array of points coordinates in a section
  COMPLEX                                     :: yypti          ! working point

  TYPE(variable), DIMENSION(:),   ALLOCATABLE :: stypvar        ! structure of output

  CHARACTER(LEN=256)                          :: cf_vtfil       ! VT file  (in)
  CHARACTER(LEN=256)                          :: cf_tfil        ! T file  (in)
  CHARACTER(LEN=256)                          :: cf_ufil        ! U file   (in)
  CHARACTER(LEN=256)                          :: cf_vfil        ! V file   (in)
  CHARACTER(LEN=256)                          :: cf_out='section_trp.dat'  ! output file name (ASCII)
  CHARACTER(LEN=256)                          :: cf_outnc            ! output netcdf file
  CHARACTER(LEN=256)                          :: cf_vtrp='vtrp.txt'  ! output volume transport file
  CHARACTER(LEN=256)                          :: cf_htrp='htrp.txt'  ! output heat transport file
  CHARACTER(LEN=256)                          :: cf_strp='strp.txt'  ! output salt transport file
  CHARACTER(LEN=256)                          :: csection            ! section names
  CHARACTER(LEN=256)                          :: cvarname            ! variable names (root)
  CHARACTER(LEN=256)                          :: clongname           ! variable longname (root)
  CHARACTER(LEN=512)                          :: cglobal             ! global attribute
  CHARACTER(LEN=256)                          :: cldum               ! dummy char variable
  CHARACTER(LEN=256)                          :: cline               ! dummy char variable
  CHARACTER(LEN=256)                          :: csfx='transports'   ! suffix for the netcdf outputfile
  CHARACTER(LEN=256), DIMENSION(3)            :: cldumt              ! dummy char variable

  LOGICAL                                     :: ltest   = .FALSE.   ! flag for test case
  LOGICAL                                     :: lfull   = .FALSE.   ! flag for full step case
  LOGICAL                                     :: lheat   = .TRUE.    ! flag for skipping heat/salt transport computation
  LOGICAL                                     :: lchk    = .FALSE.   ! flag for missing files
  LOGICAL                                     :: lpm     = .FALSE.   ! flag for plus/minus transport
  LOGICAL                                     :: lobc    = .FALSE.   ! flag for obc input files
  LOGICAL                                     :: l_merid = .FALSE.   ! flag for meridional obc
  LOGICAL                                     :: l_zonal = .FALSE.   ! flag for zonal obc
  LOGICAL                                     :: l_tsfil = .FALSE.   ! flag for using T file instead of VT file
  LOGICAL                                     :: l_self  = .FALSE.   ! flag for self mesh/mask files in the input
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  ! Print usage if no argument
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdftransport -u U-file -v V-file [-t T-file] [-vt VT-file] ...'
     PRINT *,'                  ... [-test  u v ] [-noheat ] [-pm ] [-obc] [-TS] ... '
     PRINT *,'                  ... [-full] [-time jt] [-vvl] [-self] ...'
     PRINT *,'                  ... [-zlimit dep_list] [-sfx suffix]'
     PRINT *,'      '
     PRINT *,'    PURPOSE :'
     PRINT *,'      Compute the transports (volume, heat and salt) accross a section.'
     PRINT *,'      The name of the section and the imin, imax, jmin, jmax for the section '
     PRINT *,'      is read from the standard input. To finish the program use the key name'
     PRINT *,'      ''EOF'' for the section name. It may be usefull to use the syntax :'
     PRINT *,'      ''cdftransport [..options.. ] < section_file.txt''. Section corresponds'
     PRINT *,'      to a line between 2 F-points A(imin, jmin) and B(imax,jmax). The order'
     PRINT *,'      of the points does matter: when travelling from A to B, transports to '
     PRINT *,'      the right are positive, transports to the left are negative.   '
     PRINT *,'      OBC U,V files can be used if -obc option is specified.'
     PRINT *,'      Extracted broken lines files can also be used with option -self.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'      -u U-file  : netcdf file with the zonal velocity component.'
     PRINT *,'      -v V-file  : netcdf file with the meridional velocity component.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'      [-vt VT-file]: netcdf file with mean values of vt, vs, ut, us for heat'
     PRINT *,'                    and salt transport. This option is mandatory unless '
     PRINT *,'                    -noheat option is used.'
     PRINT *,'      [-t T-file]: Temperature and Salinity file used with -TS option.'
     PRINT *,'      [-test u v]: use constant the u and v velocity components for sign '
     PRINT *,'                    test purpose.'
     PRINT *,'      [-noheat ] : use when heat and salt transport are not requested.'
     PRINT *,'      [-pm ]     : separate positive and negative contribution to'
     PRINT *,'                    the volume transport. This option implicitly set -noheat,'
     PRINT *,'                    and must be used before the file names.'
     PRINT *,'      [-obc ]    : indicates that input files are obc files (vertical slices)'
     PRINT *,'                   Take care that for this case, mesh files must be adapted.'
     PRINT *,'                   This option implicitly set -noheat, and must be used before'
     PRINT *,'                   the file names.'
     PRINT *,'      [-TS ]     : Indicate that UT VT US VS will be recomputed from T U V'
     PRINT *,'                   files. In this case use -t T-file option.'
     PRINT *,'      [-full ]   : use for full step configurations.'
     PRINT *,'      [-time jt ]: compute transports for time index jt. Default is 1.'
     PRINT *,'      [-vvl  ]   : use time varying  vertical metrics e3 read in the data file'
     PRINT *,'      [-zlimit dep_list] : Specify a list of depth (meters) defining the '
     PRINT *,'                   limits of classes for which transports will be computed.'
     PRINT *,'                   If not used, the transports are computed for the whole water'
     PRINT *,'                   column. Example : -zlimit 500 1000 creates 3 classes :'
     PRINT *,'                           0-500 '
     PRINT *,'                         500-1000 '
     PRINT *,'                        1000-bottom '
     PRINT *,'      [-self ] : This option indicates that input files corresponds to a '
     PRINT *,'                  broken line, hence  data files hold the metrics. In this '
     PRINT *,'                  case, the input file is considered as a V-file, and must'
     PRINT *,'                  specified with -v option.'
     PRINT *,'      [-sfx suffix] : use suffix instead of ',TRIM(csfx), ' in the netcdf'
     PRINT *,'                  output file. '
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'      Files ',TRIM(cn_fhgr),', ',TRIM(cn_fzgr),' must be in the current directory.'
     PRINT *,'            unless -self option is used.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'      - Standard output '
     PRINT *,'      - ASCII file reflecting the standard output: section_<SUFFIX>.dat'
     PRINT *,'      - ASCII files for volume, heat and salt transport: v<SUFFIX>.txt,'
     PRINT *,'          h<SUFFIX>.txt and s<SUFFIX>.txt.'
     PRINT *,'      - Netcdf files for each section. <SECTION>_<SUFFIX>.nc. '
     PRINT *,'        Default <SUFFIX> is ',TRIM(csfx),' and can be changed with -sfx option.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfsigtrp cdf_xtrac_brokenline'
     PRINT *,'      '
     STOP
  ENDIF

  itime  = 1
  nclass = 1
  ijarg  = 1
  CALL SetGlobalAtt(cglobal)

  ! Browse command line for arguments and/or options
  DO WHILE (ijarg <= narg )
     CALL getarg (ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ('-vt'    ) ; CALL getarg (ijarg, cf_vtfil) ; ijarg=ijarg+1 
     CASE ('-u'     ) ; CALL getarg (ijarg, cf_ufil ) ; ijarg=ijarg+1 
     CASE ('-v'     ) ; CALL getarg (ijarg, cf_vfil ) ; ijarg=ijarg+1 
     CASE ('-t'     ) ; CALL getarg (ijarg, cf_tfil ) ; ijarg=ijarg+1 
     CASE ('-test ' ) ; CALL getarg (ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) udum
        ;               CALL getarg (ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) vdum
        ;               ltest = .TRUE. 
     CASE ('-full'  ) ; lfull = .TRUE.
     CASE ('-noheat') ; lheat = .FALSE.
     CASE ('-time'  ) ; CALL getarg (ijarg, cldum   ) ; ijarg = ijarg + 1 ; READ(cldum,*) itime
     CASE ('-pm'    ) ; lpm   = .TRUE.
        ;               lheat = .FALSE.
     CASE ('-obc'   ) ; lobc  = .TRUE.
        ;               lheat = .FALSE.
     CASE ( '-TS'   ) ; l_tsfil = .TRUE.
     CASE ( '-vvl'  ) ; lg_vvl = .TRUE.
     CASE ( '-self' ) ; l_self = .TRUE.
     CASE ('-zlimit') ; CALL GetZlimit  ! Parse arguments of -zlimit to build rz_lst array
     CASE ('-sfx'   ) ; CALL getarg (ijarg, csfx   ) ; ijarg = ijarg + 1
     CASE DEFAULT     ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP
     END SELECT
  END DO

  it = 1
  IF ( lg_vvl ) THEN
     ! when vvl vertical metrics is stored in the data file with name cn_ve3u 
     cn_fe3u = cf_ufil
     cn_fe3v = cf_vfil
     it      = itime
  ENDIF

  IF ( l_self ) THEN
     cn_fzgr  = cf_vfil
     cn_fhgr  = cf_vfil
     cn_fe3u  = cf_vfil
     cn_fe3v  = cf_vfil
     cf_vtfil = cf_vfil
     cf_tfil  = cf_vfil
  ENDIF

  ! checking if all required files are available
  lchk = lchk .OR. chkfile(cn_fzgr)
  lchk = lchk .OR. chkfile(cn_fhgr)
  IF ( ltest ) THEN
     ! OK
  ELSE
     lchk = lchk .OR. chkfile(cf_ufil)
     lchk = lchk .OR. chkfile(cf_vfil)
     IF (l_tsfil)  THEN 
        lchk = lchk .OR. chkfile(cf_tfil )
     ELSE IF (lheat  ) THEN 
        lchk = lchk .OR. chkfile(cf_vtfil)
     ENDIF
  ENDIF
  IF ( lchk ) STOP ! missing files

  ! adjust the number of output variables according to options
  IF ( nclass > 1 ) THEN
     IF ( lheat ) THEN ; nvarout = 12
     ELSE              ; nvarout =  8
     ENDIF
     IF ( lpm   ) nvarout=nvarout+4
  ELSE
     IF ( lheat ) THEN ; nvarout =  9
     ELSE              ; nvarout =  7
     ENDIF
     IF ( lpm   ) nvarout=nvarout+2
  ENDIF

  ALLOCATE ( ilev0(nclass), ilev1(nclass), rclass(nclass) )
  rclass=(/(jclass, jclass=1,nclass)/)

  npiglo = getdim (cf_vfil,cn_x)
  npjglo = getdim (cf_vfil,cn_y)
  npk    = getdim (cf_vfil,cn_z)
  npt    = getdim (cf_vfil,cn_t)

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  IF ( lobc ) THEN  ! if lobc false, l_merid and l_zonal are false (default)
     IF ( npiglo == 1 ) THEN
        l_merid=.TRUE.
        ALLOCATE (zuobc(npjglo,npk), zvobc(npjglo,npk) )
        PRINT *,' Meridional OBC'
     ENDIF

     IF ( npjglo == 1 ) THEN
        l_zonal=.TRUE.
        ALLOCATE (zuobc(npiglo,npk), zvobc(npiglo,npk) )
        PRINT *,' Zonal OBC'
     ENDIF
  ENDIF

  ALLOCATE ( e31d(npk) )

  ! define new variables for output 
  ALLOCATE ( stypvar(nvarout), ipk(nvarout), id_varout(nvarout) )
  ALLOCATE ( rdum(1,1) )

  rdum(:,:)=0.e0

  ! Allocate arrays
  ALLOCATE (   zu(npiglo,npjglo),   zv(npiglo,npjglo) )
  ALLOCATE ( dwku(npiglo,npjglo), dwkv(npiglo,npjglo) )
  ALLOCATE ( dtrpu(npiglo,npjglo,nclass), dtrpv(npiglo,npjglo,nclass))
  ALLOCATE ( dvoltrpsum(nclass), dvolallegcl(nclass) )

  IF ( lpm ) THEN 
     ALLOCATE ( dwkup(npiglo,npjglo), dwkvp(npiglo,npjglo) )
     ALLOCATE ( dwkum(npiglo,npjglo), dwkvm(npiglo,npjglo) )
     ALLOCATE ( dtrpup(npiglo,npjglo,nclass), dtrpvp(npiglo,npjglo,nclass))
     ALLOCATE ( dtrpum(npiglo,npjglo,nclass), dtrpvm(npiglo,npjglo,nclass))
     ALLOCATE ( dvoltrpsump(nclass),  dvoltrpsumm(nclass)  )
     ALLOCATE ( dvolallegclp(nclass), dvolallegclm(nclass) )
  ENDIF

  IF ( lheat ) THEN
     ALLOCATE (   zut(npiglo,npjglo),   zus(npiglo,npjglo) )
     ALLOCATE (   zvt(npiglo,npjglo),   zvs(npiglo,npjglo) )
     ALLOCATE ( dwkut(npiglo,npjglo), dwkus(npiglo,npjglo) )
     ALLOCATE ( dwkvt(npiglo,npjglo), dwkvs(npiglo,npjglo) )
     ALLOCATE ( dtrput(npiglo,npjglo,nclass), dtrpvt(npiglo,npjglo,nclass))
     ALLOCATE ( dtrpus(npiglo,npjglo,nclass), dtrpvs(npiglo,npjglo,nclass))
     ALLOCATE ( dheatrpsum(nclass), dsaltrpsum(nclass)     )
     ALLOCATE ( dheatallegcl(nclass), dsaltallegcl(nclass) )
     IF ( l_tsfil ) THEN
        ALLOCATE (zt(npiglo+1,npjglo+1), zs(npiglo+1, npjglo+1) )
     ENDIF
  ENDIF
  !
  ALLOCATE ( e1v(npiglo,npjglo),e3v(npiglo,npjglo)       )
  ALLOCATE ( e2u(npiglo,npjglo),e3u(npiglo,npjglo)       )
  !
  ALLOCATE ( gphif(npiglo,npjglo) )
  ALLOCATE ( glamf(npiglo,npjglo) )
  ALLOCATE ( gdepw(npk) , tim(npt)                       )
  !
  ! read metrics and grid position
  e1v(:,:)   = getvar(cn_fhgr, cn_ve1v, 1, npiglo, npjglo)
  IF ( l_self ) THEN
     e2u(:,:)   = 1.  ! dummy value, not used
     glamf(:,:) = 1.
     gphif(:,:) = 1.
     ! use e31d for temporary calculation
     e31d(:)    =  getvar1d(cf_vtfil, cn_vdeptht, npk)
     gdepw(1)   = 0.
     gdepw(2:npk) = 0.5 * (e31d(1:npk-1) + e31d(2:npk)) ! This is just a proxy for gdepw
     ! max error ~ 1m for 46 lev
     e31d(:)    = 1.  ! set dummy value for e31d 
  ELSE
     e2u(:,:)   = getvar(cn_fhgr, cn_ve2u, 1, npiglo, npjglo)

     glamf(:,:) = getvar(cn_fhgr, cn_glamf, 1,npiglo, npjglo)
     gphif(:,:) = getvar(cn_fhgr, cn_gphif, 1,npiglo, npjglo)

     gdepw(:) = getvare3(cn_fzgr, cn_gdepw, npk)
     e31d(:)  = getvare3(cn_fzgr, cn_ve3t,  npk) ! used only for full step
  ENDIF


  IF ( lobc ) THEN
     ! read u, v on OBC
     IF ( l_zonal ) THEN   ! (jpiglo,jpk)
        zuobc(:,:)= getvarxz(cf_ufil, cn_vozocrtx, 1, npiglo, npk)
        zvobc(:,:)= getvarxz(cf_vfil, cn_vomecrty, 1, npiglo, npk)
     ENDIF
     IF ( l_merid ) THEN   ! (jpjglo,jpk)
        zuobc(:,:)= getvaryz(cf_ufil, cn_vozocrtx, 1, npjglo, npk)
        zvobc(:,:)= getvaryz(cf_vfil, cn_vomecrty, 1, npjglo, npk)
     ENDIF
  ENDIF

  ! look for nearest level to rz_lst and setup ilev0 and ilev1 (t-index of class limit)
  ik = 1
  ilev0(1) = 1 ; ilev1(nclass) = npk-1  ! default value if nclass=1
  DO jclass = 1, nclass -1
     DO WHILE ( gdepw(ik)  < rz_lst(jclass) )
        ik = ik +1
     END DO

     rd1 = ABS(gdepw(ik-1) - rz_lst(jclass) )
     rd2 = ABS(gdepw(ik  ) - rz_lst(jclass) )
     IF ( rd2 < rd1 ) THEN
        ilev1(jclass  ) = ik - 1  ! t-levels index
        ilev0(jclass+1) = ik
     ELSE 
        ilev1(jclass  ) = ik - 2  ! t-levels index
        ilev0(jclass+1) = ik - 1
     END IF
  END DO

  PRINT *, 'Limits :  '
  DO jclass = 1, nclass
     PRINT *, ilev0(jclass),ilev1(jclass), gdepw(ilev0(jclass)), gdepw(ilev1(jclass)+1)
  END DO

  ! compute the transports at each grid cell
  dtrpu (:,:,:)= 0.d0 ; dtrpv (:,:,:)= 0.d0    ! initialization to 0

  IF ( lpm   ) THEN
     dtrpup(:,:,:)= 0.d0 ; dtrpvp(:,:,:)= 0.d0  
     dtrpum(:,:,:)= 0.d0 ; dtrpvm(:,:,:)= 0.d0
  ENDIF
  IF ( lheat ) THEN
     dtrput(:,:,:)= 0.d0 ; dtrpvt(:,:,:)= 0.d0  
     dtrpus(:,:,:)= 0.d0 ; dtrpvs(:,:,:)= 0.d0
  ENDIF

  DO jclass = 1, nclass
     DO jk = ilev0(jclass),ilev1(jclass)
        PRINT *,'level ',jk
        ! Get velocities, temperature and salinity fluxes at jk
        IF ( ltest ) THEN
           zu (:,:) = udum ; zv (:,:) = vdum
           IF (lheat) THEN
              zut(:,:) = udum ; zvt(:,:) = vdum
              zus(:,:) = udum ; zvs(:,:) = vdum
           ENDIF
        ELSEIF ( lobc ) THEN
           IF      ( l_zonal ) THEN ; zu(:,1)=zuobc(:,jk) ; zv(:,1)=zvobc(:,jk) 
           ELSE IF ( l_merid ) THEN ; zu(1,:)=zuobc(:,jk) ; zv(1,:)=zvobc(:,jk) 
           ENDIF
        ELSE
           IF ( l_self ) THEN
              zu(:,:) = 0.
           ELSE
              zu (:,:) = getvar(cf_ufil, cn_vozocrtx, jk, npiglo, npjglo, ktime=itime)
           ENDIF
           zv (:,:) = getvar(cf_vfil, cn_vomecrty, jk, npiglo, npjglo, ktime=itime)
           IF (lheat) THEN
              IF ( l_tsfil ) THEN
                 zt(:,:) = 0. ; zs(:,:) = 0.
                 zt(1:npiglo,1:npjglo) =  getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=itime)
                 zs(1:npiglo,1:npjglo) =  getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=itime)
                 IF  (l_self ) THEN
                    zut(:,:) = 0.
                    zus(:,:) = 0.
                    zvt(1:npiglo,1:npjglo) = zv(1:npiglo,1:npjglo) *  zt(1:npiglo,1:npjglo)
                    zvs(1:npiglo,1:npjglo) = zv(1:npiglo,1:npjglo) *  zs(1:npiglo,1:npjglo)
                 ELSE
                    zut(1:npiglo,1:npjglo) = zu(1:npiglo,1:npjglo) * 0.5* ( zt(1:npiglo,1:npjglo) + zt(2:npiglo+1,1:npjglo  ))
                    zus(1:npiglo,1:npjglo) = zu(1:npiglo,1:npjglo) * 0.5* ( zs(1:npiglo,1:npjglo) + zs(2:npiglo+1,1:npjglo  ))
                    zvt(1:npiglo,1:npjglo) = zv(1:npiglo,1:npjglo) * 0.5* ( zt(1:npiglo,1:npjglo) + zt(1:npiglo,  2:npjglo+1))
                    zvs(1:npiglo,1:npjglo) = zv(1:npiglo,1:npjglo) * 0.5* ( zs(1:npiglo,1:npjglo) + zs(1:npiglo,  2:npjglo+1))
                 ENDIF
              ELSE
                 IF ( l_self ) THEN
                    zut(:,:) = 0.
                    zus(:,:) = 0.
                 ELSE
                    zut(:,:) = getvar(cf_vtfil, cn_vozout,   jk, npiglo, npjglo, ktime=itime)
                    zus(:,:) = getvar(cf_vtfil, cn_vozous,   jk, npiglo, npjglo, ktime=itime)
                 ENDIF
                 zvt(:,:) = getvar(cf_vtfil, cn_vomevt,   jk, npiglo, npjglo, ktime=itime)
                 zvs(:,:) = getvar(cf_vtfil, cn_vomevs,   jk, npiglo, npjglo, ktime=itime)
              ENDIF
           ENDIF
        ENDIF

        ! get e3u, e3v  at level jk
        IF ( lfull ) THEN 
           e3v(:,:) = e31d(jk)
           e3u(:,:) = e31d(jk)
        ELSE
           IF ( l_self) THEN
              e3u(:,:) = 1. !dummy value
              e3v(:,:) = getvar(cn_fe3v, cn_ve3v, jk, npiglo, npjglo, ktime=it, ldiom=.FALSE.)  ! In broken line name is cn_ve3v
           ELSE
              e3u(:,:) = getvar(cn_fe3u, cn_ve3u, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl)
              e3v(:,:) = getvar(cn_fe3v, cn_ve3v, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl)
           ENDIF
        ENDIF

        dwku (:,:) = zu (:,:)*e2u(:,:)*e3u(:,:)*1.d0
        dwkv (:,:) = zv (:,:)*e1v(:,:)*e3v(:,:)*1.d0

        IF ( lpm ) THEN 
           dwkup = 0.d0 ; dwkum = 0.d0
           WHERE ( zu >= 0. ) 
              dwkup(:,:) = zu (:,:)*e2u(:,:)*e3u(:,:)*1.d0
           ELSEWHERE         
              dwkum(:,:) = zu (:,:)*e2u(:,:)*e3u(:,:)*1.d0
           END WHERE

           dwkvp = 0.d0 ; dwkvm = 0.d0
           WHERE ( zv >= 0. ) 
              dwkvp(:,:) = zv (:,:)*e1v(:,:)*e3v(:,:)*1.d0
           ELSEWHERE         ;
              dwkvm(:,:) = zv (:,:)*e1v(:,:)*e3v(:,:)*1.d0
           END WHERE
        ENDIF

        IF ( lheat ) THEN
           dwkut(:,:) = zut(:,:)*e2u(:,:)*e3u(:,:)*1.d0
           dwkvt(:,:) = zvt(:,:)*e1v(:,:)*e3v(:,:)*1.d0
           dwkus(:,:) = zus(:,:)*e2u(:,:)*e3u(:,:)*1.d0
           dwkvs(:,:) = zvs(:,:)*e1v(:,:)*e3v(:,:)*1.d0
        ENDIF

        ! integrates vertically 
        dtrpu (:,:,jclass) = dtrpu (:,:,jclass) + dwku (:,:)
        dtrpv (:,:,jclass) = dtrpv (:,:,jclass) + dwkv (:,:)

        IF ( lpm   ) THEN
           dtrpup(:,:,jclass) = dtrpup(:,:,jclass) + dwkup(:,:)
           dtrpvp(:,:,jclass) = dtrpvp(:,:,jclass) + dwkvp(:,:)
           dtrpum(:,:,jclass) = dtrpum(:,:,jclass) + dwkum(:,:)
           dtrpvm(:,:,jclass) = dtrpvm(:,:,jclass) + dwkvm(:,:)
        ENDIF

        IF ( lheat ) THEN
           dtrput(:,:,jclass) = dtrput(:,:,jclass) + dwkut(:,:) * rau0*rcp
           dtrpvt(:,:,jclass) = dtrpvt(:,:,jclass) + dwkvt(:,:) * rau0*rcp
           dtrpus(:,:,jclass) = dtrpus(:,:,jclass) + dwkus(:,:)  
           dtrpvs(:,:,jclass) = dtrpvs(:,:,jclass) + dwkvs(:,:)  
        ENDIF

     END DO  ! loop to next level
  END DO    ! next class

  cf_out ='section_'//TRIM(csfx)//'.dat'
  cf_vtrp='v'//TRIM(csfx)//'.txt'
  cf_htrp='h'//TRIM(csfx)//'.txt'
  cf_strp='s'//TRIM(csfx)//'.txt'
  OPEN(numout,FILE=cf_out)
  ! also dump the results on txt files without any comments, some users  like it !
  OPEN(numvtrp,FILE=cf_vtrp)
  IF ( lheat ) THEN
     OPEN(numhtrp,FILE=cf_htrp) ; OPEN(numstrp,FILE=cf_strp)
  ENDIF

  !################################################################################
  ! enter interactive part 
  !################################################################################
  ! initialize all legs arrays and variable to 0 
  dvolalleg = 0.d0 ; dvolallegcl(:) = 0.d0
  IF ( lpm   ) THEN
     dvolallegp = 0.d0 ; dvolallegclp(:) = 0.d0
     dvolallegm = 0.d0 ; dvolallegclm(:) = 0.d0
  ENDIF
  IF ( lheat ) THEN
     dheatalleg = 0.d0 ; dheatallegcl(:) = 0.d0
     dsaltalleg = 0.d0 ; dsaltallegcl(:) = 0.d0
  ENDIF
  DO 
     PRINT *, ' Give name of section (EOF to finish)'
     READ(*,'(a)') cline
     ii = 0
     cldumt(:) = 'none'
     ipos = INDEX(cline,' ')
     DO WHILE ( ipos > 1 )
        ii = ii + 1
        cldumt(ii) = cline(1:ipos - 1 )
        cline = TRIM ( cline(ipos+1:) )
        ipos  = INDEX( cline,' ' )
        IF ( ii >= 3 ) EXIT
     END DO
     csection = TRIM(cldumt(1) )
     cvarname = TRIM(cldumt(2) )
     clongname = TRIM(cldumt(3) )

     IF (TRIM(csection) == 'EOF' ) THEN 
        CLOSE(numout) ; CLOSE(numvtrp) 
        IF ( lheat ) THEN 
           CLOSE(numhtrp) ; CLOSE(numstrp) 
        ENDIF
        EXIT  ! infinite DO loop
     ENDIF

     ! create output fileset
     CALL set_typvar( stypvar, csection, cvarname, clongname )
     cf_outnc = TRIM(csection)//'_'//TRIM(csfx)//'.nc'
     ncout    = create      (cf_outnc, 'none',    ikx,      iky, nclass, cdep='depth_class')
     ierr     = createvar   (ncout,    stypvar,   nvarout,  ipk, id_varout, cdglobal=TRIM(cglobal) )
     ierr     = putheadervar(ncout,    cf_ufil,   ikx, iky, nclass, pnavlon=rdum, pnavlat=rdum, pdep=rclass )
     tim      = getvar1d    (cf_ufil,  cn_vtimec, npt     )
     ierr     = putvar1d    (ncout,    tim(itime:itime),       1, 'T')

     PRINT *, ' Give iimin, iimax, ijmin, ijmax '
     READ(*,*) iimin, iimax, ijmin, ijmax
     !! Find the broken line between P1 (iimin,ijmin) and P2 (iimax, ijmax)
     ! ... Initialization
     ii0  = iimin ; ij0  = ijmin ; ii1  = iimax ;  ij1 = ijmax
     rxi0 = ii0   ; ryj0 = ij0   ; rxi1 = ii1   ; ryj1 = ij1

     ! compute direction of integrations and signs
     !The transport across the section is the dot product of
     !integral(line){(Mx,My)*dS} 
     !Mx=integral(u*dz)  My=integral(v*dz)) and dS=(dy,-dx)}

     !By defining the direction of the integration as 
     idirx = SIGN(1,ii1-ii0) !positive to the east or if ii1=ii0
     idiry = SIGN(1,ij1-ij0) !positive to the north or if ij1=ij0

     !Then dS=(e2u*idiry,-e1v*idirx)
     !This will produce the following sign convention:
     !    West-to-est line (dx>0, dy=0)=> -My*dx (-ve for a northward flow)
     !    South-to-north   (dy>0, dx=0)=>  Mx*dy (+ve for an eastward flow)
     norm_u =  idiry
     norm_v = -idirx

     ! .. Compute equation:  ryj = aj rxi + bj [valid in the (i,j) plane]
     IF ( (rxi1 -rxi0) /=  0 ) THEN
        aj = (ryj1 - ryj0 ) / (rxi1 -rxi0)
        bj = ryj0 - aj * rxi0
     ELSE
        aj = 10000.  ! flag value
        bj = 0.
     END IF

     ! .. Compute equation:  rxi = ai ryj + bi [valid in the (i,j) plane]
     IF ( (ryj1 -ryj0) /=  0 ) THEN
        ai = (rxi1 - rxi0 ) / ( ryj1 -ryj0 )
        bi = rxi0 - ai * ryj0
     ELSE
        ai = 10000. ! flag value
        bi = 0.
     END IF

     ! ..  Compute the integer pathway: a succession of F points
     np=0
     ! .. Chose the strait line with the smallest slope
     IF (ABS(aj) <=  1 ) THEN
        ! ... Here, the best line is y(x)
        ! ... If ii1 < ii0 swap points [ always describe section from left to right ]
        IF (ii1 <  ii0 ) THEN
           iitmp = ii0   ; ijtmp = ij0
           ii0   = ii1   ; ij0   = ij1
           ii1   = iitmp ; ij1   = ijtmp
        END IF

        ! iist,ijst is the grid offset to pass from F point to either U/V point
        IF ( ij1 >= ij0 ) THEN     ! line heading NE
           iist = 1 ; ijst = 1
        ELSE                       ! line heading SE
           iist = 1 ; ijst = 0
        END IF

        ! ... compute the nearest ji point on the line crossing at ji
        DO ji=ii0, ii1
           np=np+1
           IF (np > jpseg) STOP 'np > jpseg !'
           ij=NINT(aj*ji + bj )
           yypt(np) = CMPLX(ji,ij)
        END DO
     ELSE
        ! ... Here, the best line is x(y)
        ! ... If ij1 < ij0 swap points [ always describe section from bottom to top ]
        IF (ij1 <  ij0 ) THEN
           iitmp = ii0   ; ijtmp = ij0
           ii0   = ii1   ; ij0   = ij1
           ii1   = iitmp ; ij1   = ijtmp
        END IF

        ! iist,ijst is the grid offset to pass from F point to either U/V point
        IF ( ii1 >=  ii0 ) THEN
           iist = 1 ; ijst = 1
        ELSE
           iist = 0 ; ijst = 1
        END IF

        ! ... compute the nearest ji point on the line crossing at jj
        DO jj=ij0,ij1
           np=np+1
           IF (np > jpseg) STOP 'np > jpseg !'
           ii=NINT(ai*jj + bi)
           yypt(np) = CMPLX(ii,jj)
        END DO
     END IF

     !!
     !! Look for intermediate points to be added.
     !  ..  The final positions are saved in rxx,ryy
     rxx(1) = REAL(yypt(1))
     ryy(1) = IMAG(yypt(1))
     nn     = 1

     DO jk=2,np
        ! .. distance between 2 neighbour points
        rd=ABS(yypt(jk)-yypt(jk-1))
        ! .. intermediate points required if rd > 1
        IF ( rd > 1 ) THEN
           CALL interm_pt(yypt, jk, ai, bi, aj, bj, yypti)
           nn=nn+1
           IF (nn > jpseg) STOP 'nn>jpseg !'
           rxx(nn) = REAL(yypti)
           ryy(nn) = IMAG(yypti)
        END IF
        nn=nn+1
        IF (nn > jpseg) STOP 'nn>jpseg !'
        rxx(nn) = REAL(yypt(jk))
        ryy(nn) = IMAG(yypt(jk))
     END DO
     ! record longitude and latitude of initial en endind point of the section
     gla (1) = glamf( INT(rxx(1)),  INT(ryy(1))  ) 
     gphi(1) = gphif( INT(rxx(1)),  INT(ryy(1))  ) 
     gla (2) = glamf( INT(rxx(nn)), INT(ryy(nn)) ) 
     gphi(2) = gphif( INT(rxx(nn)), INT(ryy(nn)) ) 

     ! Now extract the transport through a section 
     ! ... Check whether we need a u velocity or a v velocity
     !   Think that the points are f-points and delimit either a U segment
     !   or a V segment (iist and ijst are set in order to look for the correct
     !   velocity point on the C-grid
     PRINT *, TRIM(csection)
     PRINT *, 'IMIN IMAX JMIN JMAX', iimin, iimax, ijmin, ijmax
     WRITE(numout,*) '% Transport along a section by levels' ,TRIM(csection)
     WRITE(numout,*) '% ---- IMIN IMAX JMIN JMAX'

     dvoltrpbrtp  = 0.d0
     dvoltrpbrtpp = 0.d0
     dvoltrpbrtpm = 0.d0
     dheatrpbrtp  = 0.d0
     dsaltrpbrtp  = 0.d0
     DO jclass=1,nclass
        dvoltrpsum(jclass) = 0.d0
        IF ( lpm   ) THEN
           dvoltrpsump(jclass) = 0.d0
           dvoltrpsumm(jclass) = 0.d0
        ENDIF
        IF ( lheat ) THEN
           dheatrpsum(jclass) = 0.d0
           dsaltrpsum(jclass) = 0.d0
        ENDIF

        ! segment jseg is a line between (rxx(jseg),ryy(jseg))  and (rxx(jseg+1),ryy(jseg+1))
        DO jseg = 1, nn-1
           ii0=rxx(jseg)
           ij0=ryy(jseg)
           IF ( rxx(jseg) ==  rxx(jseg+1) ) THEN    ! meridional segment, use U velocity
              dvoltrp(jseg)= dtrpu (ii0,ij0+ijst,jclass)*norm_u

              IF ( lpm   ) THEN 
                 IF (norm_u > 0 ) THEN
                    dvoltrpp(jseg)= dtrpup(ii0,ij0+ijst,jclass)*norm_u
                    dvoltrpm(jseg)= dtrpum(ii0,ij0+ijst,jclass)*norm_u
                 ELSE
                    dvoltrpp(jseg)= dtrpum(ii0,ij0+ijst,jclass)*norm_u
                    dvoltrpm(jseg)= dtrpup(ii0,ij0+ijst,jclass)*norm_u
                 ENDIF
              ENDIF

              IF ( lheat ) THEN
                 dheatrp(jseg)= dtrput(ii0,ij0+ijst,jclass)*norm_u
                 dsaltrp(jseg)= dtrpus(ii0,ij0+ijst,jclass)*norm_u
              ENDIF
           ELSE IF ( ryy(jseg) == ryy(jseg+1) ) THEN ! zonal segment, use V velocity
              dvoltrp(jseg)=dtrpv (ii0+iist,ij0,jclass)*norm_v

              IF ( lpm   ) THEN 
                 IF (norm_v > 0 ) THEN
                    dvoltrpp(jseg)=dtrpvp(ii0+iist,ij0,jclass)*norm_v
                    dvoltrpm(jseg)=dtrpvm(ii0+iist,ij0,jclass)*norm_v
                 ELSE
                    dvoltrpp(jseg)=dtrpvm(ii0+iist,ij0,jclass)*norm_v
                    dvoltrpm(jseg)=dtrpvp(ii0+iist,ij0,jclass)*norm_v
                 ENDIF
              ENDIF

              IF ( lheat ) THEN
                 dheatrp(jseg)=dtrpvt(ii0+iist,ij0,jclass)*norm_v
                 dsaltrp(jseg)=dtrpvs(ii0+iist,ij0,jclass)*norm_v
              ENDIF
           ELSE
              PRINT *,' ERROR :',  rxx(jseg),ryy(jseg),rxx(jseg+1),ryy(jseg+1) ! likely to never happen !
           END IF
           dvoltrpsum(jclass) = dvoltrpsum(jclass) + dvoltrp(jseg)
           IF ( lpm   ) THEN 
              dvoltrpsump(jclass) = dvoltrpsump(jclass) + dvoltrpp(jseg)
              dvoltrpsumm(jclass) = dvoltrpsumm(jclass) + dvoltrpm(jseg)
           ENDIF
           IF ( lheat ) THEN
              dheatrpsum(jclass) = dheatrpsum(jclass) + dheatrp(jseg)
              dsaltrpsum(jclass) = dsaltrpsum(jclass) + dsaltrp(jseg)
           ENDIF
        END DO   ! next segment

        ! Ascii outputs :      
        IF (jclass == 1 ) THEN   ! print header when it is the first class
           PRINT '(a,2f8.2,a,2f8.2)', 'FROM (LON LAT): ', gla(1),gphi(1),' TO (LON LAT) ', gla(2), gphi(2)
           WRITE(numout,*)  '% ---- LONmin LATmin LONmax LATmax'
           WRITE(numout,*)  '% Top(m)  Bottom(m)  MassTrans(Sv) HeatTrans(PW) SaltTrans(kt/s)'
           WRITE(numout,*) 0 ,iimin, iimax, ijmin, ijmax
           WRITE(numout,9003) 0. ,gla(1),gphi(1), gla(2), gphi(2)
        ENDIF

        PRINT *, gdepw(ilev0(jclass)), gdepw(ilev1(jclass)+1)
        PRINT *, ' Mass transport : ', dvoltrpsum(jclass)/1.e6,' SV'
        WRITE(numvtrp,'(e14.6)') dvoltrpsum(jclass) 
        IF ( lpm   ) THEN
           PRINT *, ' Positive Mass transport : ', dvoltrpsump(jclass)/1.e6,' SV'
           PRINT *, ' Negative Mass transport : ', dvoltrpsumm(jclass)/1.e6,' SV'
           WRITE(numout,9002) gdepw(ilev0(jclass)), gdepw(ilev1(jclass)+1), &
                &   dvoltrpsum(jclass)/1.e6, dvoltrpsump(jclass)/1.e6, dvoltrpsumm(jclass)/1.e6
           WRITE(numvtrp,'(e14.6)') dvoltrpsump(jclass) 
           WRITE(numvtrp,'(e14.6)') dvoltrpsumm(jclass) 
        ENDIF

        IF ( lheat ) THEN
           PRINT *, ' Heat transport : ', dheatrpsum(jclass)/1.e15,' PW'
           PRINT *, ' Salt transport : ', dsaltrpsum(jclass)/1.e6,' kT/s'
           WRITE(numout,9002) gdepw(ilev0(jclass)), gdepw(ilev1(jclass)+1), &
                &   dvoltrpsum(jclass)/1.e6, dheatrpsum(jclass)/1.e15, dsaltrpsum(jclass)/1.e6
           WRITE(numhtrp,'(e14.6)') dheatrpsum(jclass)
           WRITE(numstrp,'(e14.6)') dsaltrpsum(jclass)
        ELSE
           IF ( .NOT. lpm ) WRITE(numout,9002) gdepw(ilev0(jclass)), gdepw(ilev1(jclass)+1), dvoltrpsum(jclass)/1.e6
        ENDIF

        ! netcdf output 
        IF ( nclass > 1 ) THEN
           rdum(1,1) = REAL(dvoltrpsum(jclass)/1.e6)
           ierr = putvar(ncout,id_varout(ivtrpcl), rdum, jclass, 1, 1, 1 ) 
           IF ( lpm   ) THEN
              rdum(1,1) =  REAL(dvoltrpsump(jclass)/1.e6)
              ierr = putvar(ncout,id_varout(iptrpcl), rdum, jclass, 1, 1, 1 ) 
              rdum(1,1) =  REAL(dvoltrpsumm(jclass)/1.e6)
              ierr = putvar(ncout,id_varout(imtrpcl), rdum, jclass, 1, 1, 1 ) 
           ENDIF
           IF ( lheat ) THEN
              rdum(1,1) =  REAL(dheatrpsum(jclass)/1.e15)
              ierr = putvar(ncout,id_varout(ihtrpcl), rdum, jclass, 1, 1, 1 ) 
              rdum(1,1) =  REAL(dsaltrpsum(jclass)/1.e6)
              ierr = putvar(ncout,id_varout(istrpcl), rdum, jclass, 1, 1, 1 )
           ENDIF
        ENDIF
        rdum(1,1) = REAL(gdepw(ilev0(jclass)))
        ierr = putvar(ncout,id_varout(itop), rdum, jclass, 1, 1, 1 )
        rdum(1,1) = REAL(gdepw(ilev1(jclass)+1))
        ierr = putvar(ncout,id_varout(ibot), rdum, jclass, 1, 1, 1 )

        dvoltrpbrtp = dvoltrpbrtp +  dvoltrpsum(jclass)
        IF ( lpm  ) THEN
           dvoltrpbrtpp = dvoltrpbrtpp + dvoltrpsump(jclass)
           dvoltrpbrtpm = dvoltrpbrtpm + dvoltrpsumm(jclass)
        ENDIF
        IF ( lheat) THEN
           dheatrpbrtp = dheatrpbrtp +  dheatrpsum(jclass)
           dsaltrpbrtp = dsaltrpbrtp +  dsaltrpsum(jclass)
        ENDIF
        ! save sum over legs
        dvolallegcl(jclass) = dvolallegcl(jclass) + dvoltrpsum(jclass)
        IF ( lpm   ) THEN
           dvolallegclp(jclass) = dvolallegclp(jclass) + dvoltrpsump(jclass)
           dvolallegclm(jclass) = dvolallegclm(jclass) + dvoltrpsumm(jclass)
        ENDIF
        IF ( lheat ) THEN
           dheatallegcl(jclass) = dheatallegcl(jclass) + dheatrpsum(jclass)
           dsaltallegcl(jclass) = dsaltallegcl(jclass) + dsaltrpsum(jclass)
        ENDIF
     END DO ! next class
     ! save sum over legs 
     dvolalleg = dvolalleg + dvoltrpbrtp
     IF ( lpm   ) THEN
        dvolallegp = dvolallegp + dvoltrpbrtpp
        dvolallegm = dvolallegm + dvoltrpbrtpm
     ENDIF
     IF ( lheat ) THEN
        dheatalleg = dheatalleg + dheatrpbrtp
        dsaltalleg = dsaltalleg + dsaltrpbrtp
     ENDIF

     IF ( nclass > 1 ) THEN 
        PRINT *, ' ====================================================='
        PRINT *, ' total Mass transport : ', dvoltrpbrtp/1.e6,' SV'
        IF ( lpm   ) THEN
           PRINT *, ' total positive transport : ', dvoltrpbrtpp/1.e6,' SV'
           PRINT *, ' total negative transport : ', dvoltrpbrtpm/1.e6,' SV'
        ENDIF
        IF ( lheat ) THEN
           PRINT *, ' total Heat transport : ', dheatrpbrtp/1.e15,' PW'
           PRINT *, ' total Salt transport : ', dsaltrpbrtp/1.e6,' kT/s'
        ENDIF
     ENDIF
     ierr = putvar0d(ncout,id_varout(ivtrp), REAL(dvoltrpbrtp/1.e6)        )
     IF ( lpm   ) THEN
        ierr = putvar0d(ncout,id_varout(iptrp), REAL(dvoltrpbrtpp/1.e6)    )
        ierr = putvar0d(ncout,id_varout(imtrp), REAL(dvoltrpbrtpm/1.e6)    )
     ENDIF
     IF ( lheat ) THEN
        ierr = putvar0d(ncout,id_varout(ihtrp), REAL(dheatrpbrtp/1.e15)    )
        ierr = putvar0d(ncout,id_varout(istrp), REAL(dsaltrpbrtp/1.e6 )    )
     ENDIF
     ierr = putvar0d(ncout,id_varout(ilonmin), REAL(gla(1))  )
     ierr = putvar0d(ncout,id_varout(ilonmax), REAL(gla(2))  )
     ierr = putvar0d(ncout,id_varout(ilatmin), REAL(gphi(1)) )
     ierr = putvar0d(ncout,id_varout(ilatmax), REAL(gphi(2)) )
     ierr = closeout(ncout)
  END DO ! infinite loop : gets out when input is EOF 


  PRINT *,'   '
  PRINT *,' Overall transports (sum of all legs done so far)'
  DO jclass = 1, nclass
     PRINT *, gdepw(ilev0(jclass)), gdepw(ilev1(jclass)+1)
     PRINT *, ' Mass transport : ', dvolallegcl(jclass)/1.e6,' SV'
     IF ( lpm   ) THEN
        PRINT *, ' Positive Mass transport : ', dvolallegclp(jclass)/1.e6,' SV'
        PRINT *, ' Negative Mass transport : ', dvolallegclm(jclass)/1.e6,' SV'
     ENDIF

     IF ( lheat ) THEN
        PRINT *, ' Heat transport : ', dheatallegcl(jclass)/1.e15,' PW'
        PRINT *, ' Salt transport : ', dsaltallegcl(jclass)/1.e6,' kT/s'
     ENDIF
  ENDDO

  IF ( nclass > 1 ) THEN
     PRINT *, ' ====================================================='
     PRINT *, '    Mass transport      : ', dvolalleg/1.e6,' SV'
     IF ( lpm ) THEN
        PRINT *, '     positive transport : ', dvolallegp/1.e6,' SV'
        PRINT *, '     negative transport : ', dvolallegm/1.e6,' SV'
     ENDIF
     IF ( lheat ) THEN
        PRINT *, '     heat transport     : ', dheatalleg/1.e15,' PW'
        PRINT *, '     salt transport     : ', dsaltalleg/1.e6,' kT/s'
     ENDIF
  ENDIF

9000 FORMAT(I4,6(f9.3,f8.4))
9001 FORMAT(I4,6(f9.2,f9.3))
9002 FORMAT(f9.0,f9.0,f9.2,f9.2,f9.2)
9003 FORMAT(f9.2,f9.2,f9.2,f9.2,f9.2)

CONTAINS

  SUBROUTINE set_typvar ( sd_typvar, cdsection, cdvarname, cdlongname ) 
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE set_typvar  ***
    !!
    !! ** Purpose : Initialize typvar structure for netcdfoutput at a given section
    !!
    !! ** Method  : use varname longname to suffix variable name and attributes
    !!              If varname and/or logname are not given (ie 'none') take
    !!              standard default names
    !!              Netcdf id for variables are passed as global variables
    !!----------------------------------------------------------------------
    TYPE(variable), DIMENSION(:), INTENT(out) :: sd_typvar        ! structure of output
    CHARACTER(LEN=*),             INTENT(in ) :: cdsection
    CHARACTER(LEN=*),             INTENT(in ) :: cdvarname
    CHARACTER(LEN=*),             INTENT(in ) :: cdlongname
    !!
    INTEGER(KIND=4)                           :: ivar
    CHARACTER(LEN=255)                        :: csuffixvarnam
    CHARACTER(LEN=255)                        :: cprefixlongnam
    !!----------------------------------------------------------------------
    ! set suffixes according to variable/longname 
    IF ( cdvarname /= 'none' ) THEN
       csuffixvarnam = '_'//TRIM(cdvarname)
    ELSE
       csuffixvarnam = ''
    ENDIF

    IF ( cdlongname /= 'none' ) THEN
       cprefixlongnam = TRIM(cdlongname)//'_'
    ELSE
       cprefixlongnam = ''
    ENDIF

    ! set common values
    sd_typvar%rmissing_value=99999.
    sd_typvar%scale_factor= 1.
    sd_typvar%add_offset= 0.
    sd_typvar%savelog10= 0.
    sd_typvar%conline_operation='N/A'
    sd_typvar%caxis='T'

    ! set particular values for individual variables
    ivar = 1  ; ivtrp = ivar
    ipk(ivar) = 1
    sd_typvar(ivar)%cname       = 'vtrp'//TRIM(csuffixvarnam)
    sd_typvar(ivar)%cunits      = 'Sverdrup'
    sd_typvar(ivar)%valid_min   = -500.
    sd_typvar(ivar)%valid_max   =  500.
    sd_typvar(ivar)%clong_name  = TRIM(cprefixlongnam)//'Volume_Transport'
    sd_typvar(ivar)%cshort_name = 'vtrp'

    IF ( lpm ) THEN 
       ivar = ivar + 1 ; iptrp = ivar                                                  ;  imtrp = ivar+1
       ipk(ivar) = 1                                                                   ;  ipk(ivar+1) = 1
       sd_typvar(ivar)%cname       = 'ptrp'//TRIM(csuffixvarnam)                       ;  sd_typvar(ivar+1)%cname       = 'mtrp'//TRIM(csuffixvarnam)
       sd_typvar(ivar)%cunits      = 'Sverdrup'                                        ;  sd_typvar(ivar+1)%cunits      = 'Sverdrup'
       sd_typvar(ivar)%valid_min   = -500.                                             ;  sd_typvar(ivar+1)%valid_min   = -500.
       sd_typvar(ivar)%valid_max   =  500.                                             ;  sd_typvar(ivar+1)%valid_max   =  500.
       sd_typvar(ivar)%clong_name  = TRIM(cprefixlongnam)//'Positive_volume_transport' ;  sd_typvar(ivar+1)%clong_name  = TRIM(cprefixlongnam)//'Negative_volume_transport'
       sd_typvar(ivar)%cshort_name = 'ptrp'                                            ;  sd_typvar(ivar+1)%cshort_name = 'mtrp'
       ivar = ivar + 1
    ENDIF

    IF ( lheat ) THEN
       ivar = ivar + 1 ; ihtrp = ivar                                                  ;  istrp = ivar+1
       ipk(ivar) = 1                                                                   ;  ipk(ivar+1) = 1
       sd_typvar(ivar)%cname       = 'htrp'//TRIM(csuffixvarnam)                       ;  sd_typvar(ivar+1)%cname       = 'strp'//TRIM(csuffixvarnam)
       sd_typvar(ivar)%cunits      = 'PW'                                              ;  sd_typvar(ivar+1)%cunits      = 'kt/s'
       sd_typvar(ivar)%valid_min   = -1000.                                            ;  sd_typvar(ivar+1)%valid_min   = -1000.
       sd_typvar(ivar)%valid_max   =  1000.                                            ;  sd_typvar(ivar+1)%valid_max   =  1000.
       sd_typvar(ivar)%clong_name  = TRIM(cprefixlongnam)//'Heat_Transport'            ;  sd_typvar(ivar+1)%clong_name  = TRIM(cprefixlongnam)//'Salt_Transport'
       sd_typvar(ivar)%cshort_name = 'htrp'                                            ;  sd_typvar(ivar+1)%cshort_name = 'strp'
       ivar = ivar + 1
    ENDIF

    ivar = ivar + 1 ; ilonmin = ivar                                                   ;  ilonmax = ivar+1
    ipk(ivar) = 1                                                                      ;  ipk(ivar+1) = 1
    sd_typvar(ivar)%cname       = 'lonmin'//TRIM(csuffixvarnam)                        ;  sd_typvar(ivar+1)%cname       = 'lonmax'//TRIM(csuffixvarnam)
    sd_typvar(ivar)%cunits      = 'deg'                                                ;  sd_typvar(ivar+1)%cunits      = 'deg'
    sd_typvar(ivar)%valid_min   = -180.                                                ;  sd_typvar(ivar+1)%valid_min   = -180.
    sd_typvar(ivar)%valid_max   =  180.                                                ;  sd_typvar(ivar+1)%valid_max   =  180.
    sd_typvar(ivar)%clong_name  = TRIM(cprefixlongnam)//'begin_longitude'              ;  sd_typvar(ivar+1)%clong_name  = TRIM(cprefixlongnam)//'end_longitude'
    sd_typvar(ivar)%cshort_name = 'lonmin'                                             ;  sd_typvar(ivar+1)%cshort_name = 'lonmax'
    ivar = ivar + 1

    ivar = ivar + 1  ; ilatmin = ivar                                                  ;  ilatmax = ivar+1
    ipk(ivar) = 1                                                                      ;  ipk(ivar+1) = 1
    sd_typvar(ivar)%cname       = 'latmin'//TRIM(csuffixvarnam)                        ;  sd_typvar(ivar+1)%cname       = 'latmax'//TRIM(csuffixvarnam)
    sd_typvar(ivar)%cunits      = 'deg'                                                ;  sd_typvar(ivar+1)%cunits      = 'deg'
    sd_typvar(ivar)%valid_min   = -90.                                                 ;  sd_typvar(ivar+1)%valid_min   = -90.
    sd_typvar(ivar)%valid_max   =  90.                                                 ;  sd_typvar(ivar+1)%valid_max   =  90.
    sd_typvar(ivar)%clong_name  = TRIM(cprefixlongnam)//'begin_latitude'               ;  sd_typvar(ivar+1)%clong_name  = TRIM(cprefixlongnam)//'end_latitude'
    sd_typvar(ivar)%cshort_name = 'latmin'                                             ;  sd_typvar(ivar+1)%cshort_name = 'latmax'
    ivar = ivar + 1

    ivar = ivar + 1  ; itop = ivar                                                     ;  ibot = ivar+1
    ipk(ivar) = nclass                                                                 ;  ipk(ivar+1) = nclass
    sd_typvar(ivar)%cname       = 'top'                                                ;  sd_typvar(ivar+1)%cname       = 'bottom'
    sd_typvar(ivar)%cunits      = 'meters'                                             ;  sd_typvar(ivar+1)%cunits      = 'meters'
    sd_typvar(ivar)%valid_min   = 0.                                                   ;  sd_typvar(ivar+1)%valid_min   = 0.
    sd_typvar(ivar)%valid_max   = 10000.                                               ;  sd_typvar(ivar+1)%valid_max   = 10000.
    sd_typvar(ivar)%clong_name  = 'class_min_depth'                                    ;  sd_typvar(ivar+1)%clong_name  = 'class_max_depth'
    sd_typvar(ivar)%cshort_name = 'top'                                                ;  sd_typvar(ivar+1)%cshort_name = 'bottom'
    ivar = ivar + 1

    ivtrpcl = -1  ; ihtrpcl = -1 ; istrpcl = -1
    IF ( nclass > 1 ) THEN  ! define additional variable for vertical profile of transport (per class)
       ivar = ivar + 1  ; ivtrpcl = ivar
       ipk(ivar) = nclass                      
       sd_typvar(ivar)%cname       = 'vtrp_dep'//TRIM(csuffixvarnam) 
       sd_typvar(ivar)%cunits      = 'SV'  
       sd_typvar(ivar)%valid_min   = 0.       
       sd_typvar(ivar)%valid_max   = 10000.  
       sd_typvar(ivar)%clong_name  = TRIM(cprefixlongnam)//'Volume_Transport_per_class' 
       sd_typvar(ivar)%cshort_name = 'vtrp_dep'            

       IF ( lpm )  THEN
          ivar = ivar + 1  ; iptrpcl = ivar                                             ;  imtrpcl = ivar+1
          ipk(ivar) = nclass                                                            ;  ipk(ivar+1) = nclass
          sd_typvar(ivar)%cname       = 'ptrp_dep'//TRIM(csuffixvarnam)                 ;  sd_typvar(ivar+1)%cname       = 'mtrp_dep'//TRIM(csuffixvarnam)
          sd_typvar(ivar)%cunits      = 'SV'                                            ;  sd_typvar(ivar+1)%cunits      = 'SV'
          sd_typvar(ivar)%valid_min   = -500.                                           ;  sd_typvar(ivar+1)%valid_min   = -500.
          sd_typvar(ivar)%valid_max   =  500.                                           ;  sd_typvar(ivar+1)%valid_max   =  500.
          sd_typvar(ivar)%clong_name  = TRIM(cprefixlongnam)//'Positive_trp_per_class'  ;  sd_typvar(ivar+1)%clong_name  = TRIM(cprefixlongnam)//'Negative_trp_per_class'
          sd_typvar(ivar)%cshort_name = 'ptrp_dep'                                      ;  sd_typvar(ivar+1)%cshort_name = 'mtrp_dep'
          ivar = ivar + 1
       ENDIF

       IF ( lheat ) THEN
          ivar = ivar + 1  ; ihtrpcl = ivar                                              ;  istrpcl = ivar+1
          ipk(ivar) = nclass                                                             ;  ipk(ivar+1) = nclass
          sd_typvar(ivar)%cname       = 'htrp_dep'//TRIM(csuffixvarnam)                  ;  sd_typvar(ivar+1)%cname       = 'strp_dep'//TRIM(csuffixvarnam)
          sd_typvar(ivar)%cunits      = 'PW'                                             ;  sd_typvar(ivar+1)%cunits      = 'kt/s'
          sd_typvar(ivar)%valid_min   = -1000.                                           ;  sd_typvar(ivar+1)%valid_min   = -1000.
          sd_typvar(ivar)%valid_max   =  1000.                                           ;  sd_typvar(ivar+1)%valid_max   =  1000.
          sd_typvar(ivar)%clong_name  = TRIM(cprefixlongnam)//'Heat_Transport_per_class' ;  sd_typvar(ivar+1)%clong_name  = TRIM(cprefixlongnam)//'Salt_Transport_per_class'
          sd_typvar(ivar)%cshort_name = 'htrp_dep'                                       ;  sd_typvar(ivar+1)%cshort_name = 'strp_dep'
          ivar = ivar + 1
       ENDIF
    ENDIF

  END SUBROUTINE set_typvar

  SUBROUTINE interm_pt (ydpt, kk, pai, pbi, paj, pbj, ydpti)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE nterm_pt  ***
    !!
    !! ** Purpose : Find the best intermediate points on a pathway.
    !!
    !! ** Method  : ydpt : complex vector of the positions of the nearest points
    !!               kk  : current working index
    !!          pai, pbi : slope and original ordinate of x(y)
    !!          paj, pbj : slope and original ordinate of y(x)
    !!             ydpti : Complex holding the position of intermediate point 
    !!
    !! ** Reference : 19/07/1999 : J.M. Molines in Clipper
    !!----------------------------------------------------------------------
    COMPLEX, DIMENSION(:), INTENT(in ) :: ydpt
    COMPLEX,               INTENT(out) :: ydpti
    REAL(KIND=4),          INTENT(in ) :: pai, pbi, paj, pbj
    INTEGER(KIND=4),       INTENT(in ) :: kk
    ! ... local
    COMPLEX                            :: ylptmp1, ylptmp2
    REAL(KIND=4)                       :: za0, zb0
    REAL(KIND=4)                       :: za1, zb1
    REAL(KIND=4)                       :: zd1, zd2
    REAL(KIND=4)                       :: zxm, zym
    REAL(KIND=4)                       :: zxp, zyp
    !!----------------------------------------------------------------------
    ! ... Determines whether we use y(x) or x(y):
    IF (ABS(paj) <=  1) THEN
       ! .....  use y(x)
       ! ... possible intermediate points:
       ylptmp1=ydpt(kk-1)+(1.,0.)                 ! M1 
       ylptmp2=ydpt(kk-1)+CMPLX(0.,SIGN(1.,paj))  ! M2
       !
       ! ... M1 is the candidate point:
       zxm=REAL(ylptmp1)
       zym=IMAG(ylptmp1)
       za0=paj
       zb0=pbj
       !
       za1=-1./za0
       zb1=zym - za1*zxm
       ! ... P1 is the projection of M1 on the strait line
       zxp=-(zb1-zb0)/(za1-za0)
       zyp=za0*zxp + zb0
       ! ... zd1 is the distance M1P1
       zd1=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
       !
       ! ... M2 is the candidate point:
       zxm=REAL(ylptmp2)
       zym=IMAG(ylptmp2)
       za1=-1./za0
       zb1=zym - za1*zxm
       ! ... P2 is the projection of M2 on the strait line
       zxp=-(zb1-zb0)/(za1-za0)
       zyp=za0*zxp + zb0
       ! ... zd2 is the distance M2P2
       zd2=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
       ! ... chose the smallest (zd1,zd2)
       IF (zd2 <=  zd1) THEN
          ydpti=ylptmp2   ! use M2
       ELSE
          ydpti=ylptmp1   ! use M1
       END IF
       !
    ELSE   
       ! ...  use x(y)
       ! ... possible intermediate points:
       ylptmp1=ydpt(kk-1)+CMPLX(SIGN(1.,pai),0.)  ! M1
       ylptmp2=ydpt(kk-1)+(0.,1.)                 ! M2
       ! 
       ! ... M1 is the candidate point:
       zxm=REAL(ylptmp1)
       zym=IMAG(ylptmp1)
       za0=pai
       zb0=pbi
       !
       za1=-1./za0
       zb1=zxm - za1*zym
       zyp=-(zb1-zb0)/(za1-za0)
       zxp=za0*zyp + zb0
       zd1=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
       !
       zxm=REAL(ylptmp2)
       zym=IMAG(ylptmp2)
       za1=-1./za0
       zb1=zxm - za1*zym
       zyp=-(zb1-zb0)/(za1-za0)
       zxp=za0*zyp + zb0
       zd2=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
       IF (zd2 <=  zd1) THEN
          ydpti=ylptmp2
       ELSE
          ydpti=ylptmp1
       END IF
    END IF
  END SUBROUTINE interm_pt

  SUBROUTINE GetZlimit
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetZlimit  ***
    !!
    !! ** Purpose :  Set up a limit list given on the command line as 
    !!               blank separated list
    !!
    !! ** Method  :  Scan the command line until a '-' is found
    !!----------------------------------------------------------------------
    INTEGER (KIND=4)  :: icur
    !!----------------------------------------------------------------------
    !!
    ndep=0
    ! need to read a list of dep ( number unknow ) 
    ! loop on argument till a '-' is found as first char
    icur=ijarg                          ! save current position of argument number
    DO ji = icur, narg                  ! scan arguments till - found
       CALL getarg ( ji, cldum )
       IF ( cldum(1:1) /= '-' ) THEN ; ndep = ndep+1
       ELSE                          ; EXIT
       ENDIF
    ENDDO
    nclass = ndep+1
    ALLOCATE (rz_lst(ndep) )
    DO ji = icur, icur + ndep -1
       CALL getarg(ji, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) rz_lst( ji -icur +1 )
    END DO
  END SUBROUTINE GetZlimit

END PROGRAM cdftransport
