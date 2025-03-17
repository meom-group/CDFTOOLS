PROGRAM cdfsigtrp
  !!======================================================================
  !!                     ***  PROGRAM  cdfsigtrp  ***
  !!======================================================================
  !!  ** Purpose : Compute density class Mass transport across a section.
  !!
  !!  ** Method  :- The begining and end point of the section are given in 
  !!                term of f-points index.
  !!              - The program works for zonal or meridional sections.
  !!              - The section definitions are given in an ASCII FILE 
  !!                dens_section.dat:
  !!                 foreach sections, 2 lines :
  !!                       (i) : section name (String, no blank)
  !!                      (ii) : imin imax jmin jmax for the section
  !!              - Only vertical slices corrsponding to the sections are
  !!                read in the files.
  !!              - read metrics, depth, etc
  !!              - read normal velocity (either vozocrtx oy vomecrty )
  !!              - read 2 rows of T and S ( i i+1  or j j+1 )
  !!              - compute the mean value at velocity point
  !!              - compute sigma0 (can be easily modified for sigmai )
  !!              - compute the depths of isopyncal surfaces
  !!              - compute the transport from surface to the isopycn
  !!              - compute the transport in each class of density
  !!              - compute the total transport (for information)
  !!
  !! History : 2.1  : 03/2006  : J.M. Molines : Original code
  !!                : 07/2009  : R. Dussin    : add cdf output
  !!           3.0  : 06/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!  section_init   : initialize section names and positions
  !!  print_out      : routine which performs standard output if required
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos          ! for sigma0, sigmai
  USE modcdfnames  ! for ReadCdfNames
  USE modutils     ! for SetGlobalAtt
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class transport
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: ji, jk, jclass, jsec ! dummy loop index
  INTEGER(KIND=4)                               :: jiso, jbin           ! dummy loop index
  INTEGER(KIND=4)                               :: nbins                ! number of density classes
  INTEGER(KIND=4)                               :: ipos                 ! working variable
  INTEGER(KIND=4)                               :: narg, iargc          ! command line 
  INTEGER(KIND=4)                               :: ijarg, ireq, nreq    ! command line
  INTEGER(KIND=4)                               :: npk, nk              ! vertical size, number of wet layers
  INTEGER(KIND=4)                               :: numout=11            ! ascii output
  INTEGER(KIND=4)                               :: nsection             ! number of sections (overall)
  INTEGER(KIND=4)                               :: npiglo               ! length of broken line section
  INTEGER(KIND=4)                               :: iimin, iimax         ! working section limits
  INTEGER(KIND=4)                               :: ijmin, ijmax         ! working section limits
  INTEGER(KIND=4)                               :: npts                 ! number of points in section
  INTEGER(KIND=4)                               :: ikx=1, iky=1         ! dims of netcdf output file
  INTEGER(KIND=4)                               :: nboutput=2           ! number of values to write in cdf output
  INTEGER(KIND=4)                               :: ncout, ierr          ! for netcdf output
  INTEGER(KIND=4)                               :: iweight              ! weight of input file for further averaging
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: iimina, iimaxa       ! sections limits
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ijmina, ijmaxa       ! sections limits
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk, id_varout       ! variable levels and id

  REAL(KIND=4)                                  :: refdep =0.e0         ! reference depth (m)
  REAL(KIND=4)                                  :: zsps, zspu, zspv     ! Missing value for salinity, U and V
  REAL(KIND=4), DIMENSION(1)                    :: rdummy1, rdummy2     ! working variable
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: gdept, gdepw         ! depth of T and W points 
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: eu                   ! either e1v or e2u
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: e3t1d, e3w1d         ! vertical metrics in case of full step
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: rlonlat              ! longitudes/latitudes if the section
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zs, zt, zz           ! salinity and temperature from file 
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: rdumlon, rdumlat     ! dummy longitude and latitude for output
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zu                   ! velocity
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zmask                ! mask
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: tmpm, tmpz           ! temporary arrays

  ! double precision for cumulative variables and densities
  REAL(KIND=8)                                  :: dsigma_min           ! minimum density for bining
  REAL(KIND=8)                                  :: dsigma_max, dltsig   ! maximum density for bining, step
  REAL(KIND=8)                                  :: dsigma, dalfa        ! working sigma, interpolation coeff.
  REAL(KIND=8), DIMENSION(1)                    :: dtim                 ! time counter
  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dsigma_lev           ! built array with sigma levels
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: de3                  ! vertical metrics
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: ddepu                ! depth of vel points
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: ddepw                ! depth of W  points
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dsig                 ! density
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dhiso                ! depth of isopycns
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dwtrp, dwtrpbin      ! transport arrays
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dtrpbin              ! transport arrays

  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvar              ! structure of output

  CHARACTER(LEN=256)                            :: cf_tfil              ! temperature salinity file
  CHARACTER(LEN=256)                            :: cf_sfil              ! salinity file (option)
  CHARACTER(LEN=256)                            :: cf_ufil              ! zonal velocity file
  CHARACTER(LEN=256)                            :: cf_vfil              ! meridional velocity file
  CHARACTER(LEN=256)                            :: cf_wfil              ! W file for vvl e3w
  CHARACTER(LEN=256)                            :: cf_brk               ! broken-line file
  CHARACTER(LEN=256)                            :: cf_section='dens_section.dat'  ! input section file
  CHARACTER(LEN=256)                            :: cf_out='trpsig.txt'  ! output  ascii file
  CHARACTER(LEN=256)                            :: cf_nc                ! output netcdf file (2d)
  CHARACTER(LEN=256)                            :: cf_outnc             ! output netcdf file (1d, 0d))
  CHARACTER(LEN=256)                            :: cv_dep               ! depth variable
  CHARACTER(LEN=256)                            :: cldum                ! dummy string
  CHARACTER(LEN=256)                            :: cglobal              ! global attribute
  CHARACTER(LEN=80 )                            :: cfmt_9000            ! format string 
  CHARACTER(LEN=80 )                            :: cfmt_9001            ! format string
  CHARACTER(LEN=80 )                            :: cfmt_9002            ! format string
  CHARACTER(LEN=80 )                            :: cfmt_9003            ! format string
  CHARACTER(LEN=256)                            :: cl_vnam, cl_lname    ! working variables
  CHARACTER(LEN=256)                            :: csuffixvarname       !
  CHARACTER(LEN=256)                            :: cprefixlongname      !
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names             ! names of input variables
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: csection             ! section name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cvarname             ! output variable name (root)
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: clongname            ! output long name (root)

  LOGICAL                                       :: l_merid              ! flag for meridional section
  LOGICAL                                       :: ltemp  =.FALSE.      ! flag for use of temperature
  LOGICAL                                       :: lprint =.FALSE.      ! flag for extra print
  LOGICAL                                       :: lxtra  =.FALSE.      ! flag for extra netcdf output
  LOGICAL                                       :: lfull  =.FALSE.      ! flag for full step 
  LOGICAL                                       :: lntr   =.FALSE.      ! flag for neutral density
  LOGICAL                                       :: lchk   =.FALSE.      ! flag for missing files
  LOGICAL                                       :: lbrk   =.FALSE.      ! flag for broken lines
  LOGICAL                                       :: ll_teos10  = .FALSE. ! teos10 flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfsigtrp -t T-file -u U-file -v V-file [-s S-file] [-brk BRK-file] .'
     PRINT *,'              ... -smin sigma_min -smax sigma_max -nbins nbins [-print] ...'
     PRINT *,'              ... [-xtra] [-full ] [-vvl W-file] [-refdep ref_depth] ...'
     PRINT *,'              ... [-neutral ] [-section file ] [-temp] [-help] [-teos10]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute density class transports, according to the density class' 
     PRINT *,'       definition (minimum, maximum and number of bins) given in arguments.'
     PRINT *,'      '
     PRINT *,'       Sections information are read from the file ',TRIM(cf_section)
     PRINT *,'       which is a text file built  with pairs of lines giving: (1) section name'
     PRINT *,'       and (2) section location.'
     PRINT *,'       First line with section name may also have 2 additional strings holding'
     PRINT *,'       a prefix for variable output, and a long name to be used as attribute in'
     PRINT *,'       the output file. ' 
     PRINT *,'       Second line gives the location of the section with specification of four'
     PRINT *,'       integer values (imin imax jmin jmax), relative to the model grid.'
     PRINT *,'       Only  zonal or meridional sections are allowed.'
     PRINT *,'      '
     PRINT *,'       This program also offer the possibilty to read ''broken-line'' files,'
     PRINT *,'       holding already extracted data along a pseudo zonal section. In this '
     PRINT *,'       particular case, (''-brk'' switch), no additional information about the'
     PRINT *,'       section is required, nor metric files, as it is already available in the'
     PRINT *,'       input file.'
     PRINT *,'      '
     PRINT *,'       This program can also be used to compute transport by temperatures'
     PRINT *,'       classes, provided the temperatures decrease monotonically downward.'
     PRINT *,'       In this case, use -temp option and of course specify sigma_min, '
     PRINT *,'       sigma_max as temperatures.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       Input data file can be specified either using the 3 following switches:'
     PRINT *,'       -t T-file : netcdf file with temperature and salinity' 
     PRINT *,'          If salinity not in T-file use -s option.'
     PRINT *,'       -u U-file : netcdf file with zonal velocity component'
     PRINT *,'       -v V-file : netcdf file with meridional velocity component'
     PRINT *,'       Or, using the ''-brk'' switch.'
     PRINT *,'        -brk BRK-file : specify a ''broken-line'' file produced by the tool'
     PRINT *,'              cdf_xtrac_brokenline, which is considered already as a pseudo'
     PRINT *,'              zonal section, holding all the relevant metrics for the section.'
     PRINT *,'            '
     PRINT *,'       -smin sigma_min : minimum density for binning'
     PRINT *,'       -smax sigma_max : maximum density for binning'
     PRINT *,'       -nbins nbins : number of bins. This will fix the bin ''width'' '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-s S-file ]: specify salinity file if not T-file.'
     PRINT *,'       [-full] : for full step configuration' 
     PRINT *,'       [-vvl W-file]: use time varying vertical metrics. Need a W-file for e3w.'
     PRINT *,'       [-xtra] : produce extra netcdf output file which shows the details'
     PRINT *,'               of the sections (normal velocity, density, temperature, '
     PRINT *,'               salinity, transports, isopycnal depths. '
     PRINT *,'       [-print]: write the binned transports on standard output, for each'
     PRINT *,'               sections.'
     PRINT *,'       [-refdep ref_depth ]: give a reference depths for the computation of'
     PRINT *,'               potential density. Sigma_min, sigma_max must be adapted '
     PRINT *,'               accordingly.'
     PRINT *,'       [-neutral]: use neutral density instead of potential density '
     PRINT *,'       [-section file] : give the name of section file.'
     PRINT *,'               Default is ', TRIM(cf_section)
     PRINT *,'       [-temp] : use temperature instead of density for binning'
     PRINT *,'       [-help] : give commented example for the section file.'
     PRINT *,'       [-teos10] : use TEOS10 equation of state instead of default EOS80'
     PRINT *,'                 Temperature should be conservative temperature (CT) in deg C.'
     PRINT *,'                 Salinity should be absolute salinity (SA) in g/kg.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fhgr),', ', TRIM(cn_fzgr),' and ', TRIM(cf_section)
     PRINT *,'         If option ''-brk'' is used, there is no need for these files, the '
     PRINT *,'         metrics being already embarked into this section file.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       Netcdf file : There is 1 netcdf file per section. File name is build'
     PRINT *,'         from section name : Section_name_trpsig.nc'
     PRINT *,'         variables : sigma_class (upper limit of the bin)'
     PRINT *,'                     sigtrp : transport (Sv per bin)'
     PRINT *,'      '
     PRINT *,'       ascii file  : ', TRIM(cf_out) 
     PRINT *,'      '
     PRINT *,'      Standard output : the results are written on standard output only if '
     PRINT *,'         the -print option is used.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      cdfrhoproj, cdftransport, cdfsigintegr '
     PRINT *,'      '
     STOP 
  ENDIF

  ! browse command line
  ijarg = 1 ; ireq = 0 ; nreq=6
  cf_sfil='none'
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-t'      ) ; CALL getarg(ijarg, cf_tfil ) ; ijarg=ijarg+1 ; ireq=ireq+1
     CASE ( '-u'      ) ; CALL getarg(ijarg, cf_ufil ) ; ijarg=ijarg+1 ; ireq=ireq+1
     CASE ( '-v'      ) ; CALL getarg(ijarg, cf_vfil ) ; ijarg=ijarg+1 ; ireq=ireq+1
     CASE ( '-brk'    ) ; CALL getarg(ijarg, cf_brk  ) ; ijarg=ijarg+1 ; ireq=ireq+1
        ;                 lbrk   = .TRUE. ; nreq= 4
     CASE ( '-smin'   ) ; CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) dsigma_min ; ireq=ireq+1
     CASE ( '-smax'   ) ; CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) dsigma_max ; ireq=ireq+1
     CASE ( '-nbins'  ) ; CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) nbins      ; ireq=ireq+1
        ! options
     CASE ( '-s'      ) ; CALL getarg(ijarg, cf_sfil ) ; ijarg=ijarg+1 ; ireq=ireq+1
                          nreq = 7
     CASE ( '-full'   ) ; lfull  = .TRUE.
     CASE ( '-vvl'    ) ; lg_vvl = .TRUE.
        ;                 CALL getarg(ijarg, cf_wfil ) ; ijarg=ijarg+1 
     CASE ( '-xtra'   ) ; lxtra  = .TRUE.
     CASE ( '-print'  ) ; lprint = .TRUE.
     CASE ( '-temp'   ) ; ltemp  = .TRUE. 
     CASE ( '-help'   ) ; CALL file_example ; STOP 99
     CASE ( '-refdep' ) ; CALL getarg(ijarg, cldum      ) ; ijarg=ijarg+1 ; READ(cldum,*) refdep
     CASE ( '-section') ; CALL getarg(ijarg, cf_section ) ; ijarg=ijarg+1 
     CASE ( '-neutral') ; lntr   = .TRUE.
     CASE ( '-teos10' ) ; ll_teos10 = .TRUE. 
     CASE DEFAULT       ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  CALL eos_init ( ll_teos10 )

  IF ( cf_sfil == 'none' ) cf_sfil=cf_tfil

  IF ( ireq /= nreq ) THEN
     PRINT *,' ERROR :  not enough input arguments. See the usage message and correct.' 
     STOP 99
  ENDIF

  IF ( lbrk ) THEN
     cn_fzgr = cf_brk
     cn_fhgr = cf_brk
     cf_tfil = cf_brk
     cf_sfil = cf_brk
     cf_ufil = cf_brk
     cf_vfil = cf_brk 
     cf_section = cf_brk  ! never used in this case, just for check
  ENDIF

  ! check for file existence
  lchk = lchk .OR. chkfile( cn_fzgr    )
  lchk = lchk .OR. chkfile( cn_fhgr    )
  lchk = lchk .OR. chkfile( cf_section )
  lchk = lchk .OR. chkfile( cf_tfil    )
  lchk = lchk .OR. chkfile( cf_sfil    )
  lchk = lchk .OR. chkfile( cf_ufil    )
  lchk = lchk .OR. chkfile( cf_vfil    )
  IF ( lchk ) STOP 99 ! missing file

  ! Look for missing value for salinity, U and V
  zsps = getspval(cf_sfil, cn_vosaline )
  zspu = getspval(cf_ufil, cn_vozocrtx )
  zspv = getspval(cf_vfil, cn_vomecrty )

  IF ( lg_vvl )  THEN
     cn_fe3u = cf_ufil
     cn_fe3v = cf_vfil
     cn_fe3w = cf_wfil
     cn_ve3u = cn_ve3uvvl
     cn_ve3v = cn_ve3vvvl
     cn_ve3w = cn_ve3wvvl
  ENDIF

  IF ( ltemp)  THEN  ! temperature decrease downward. Change sign and swap min/max
     refdep = -10.   ! flag value
     dltsig     = dsigma_max  ! use dltsig as dummy variable for swapping
     dsigma_max = -dsigma_min
     dsigma_min = -dltsig
  ENDIF

  ! define global attribute with command line
  CALL SetGlobalAtt( cglobal)

  ! get the attribute iweight from vomecrty
  iweight = getatt(cf_vfil, cn_vomecrty, 'iweight')
  IF ( iweight == 0 ) iweight = 1  ! if 0 means that it is not defined.

  ALLOCATE ( stypvar(nboutput), ipk(nboutput), id_varout(nboutput) )
  ALLOCATE ( rdumlon(ikx,iky),  rdumlat(ikx,iky)                   )

  rdumlon(:,:)=0.
  rdumlat(:,:)=0.

  ipk(1)=nbins ! sigma for each level
  ipk(2)=nbins ! transport for each level
  ! initialisation of variable names etc... is done according to section name

  ! Initialise sections from file 
  ! first call to get nsection and allocate arrays 
  IF ( lbrk ) THEN 
     !npiglo = getdim (cf_brk, cn_x)
     !nsection = 1 ; iimina=1 ; iimaxa=npiglo ; ijmina=1 ; ijmaxa=1
     nsection = 1
  ELSE             
     nsection = 0 
     CALL section_init(cf_section, csection,cvarname,clongname,iimina, iimaxa, ijmina, ijmaxa, nsection)
  ENDIF
  ALLOCATE ( csection(nsection), cvarname(nsection), clongname(nsection) )
  ALLOCATE ( iimina(nsection), iimaxa(nsection), ijmina(nsection),ijmaxa(nsection) )
  IF ( lbrk )  THEN ! initialise section for true section
     npiglo = getdim (cf_brk, cn_x)
     iimina = 1 ; iimaxa = npiglo ; ijmina = 1 ; ijmaxa = 1
     ! file name for BRK file is <prefix>_<section>.nc
     ipos         = INDEX(cf_brk,'_',.TRUE.)
     csection(1)  = cf_brk(ipos+1:)
     ipos         = INDEX(csection(1),'.',.TRUE.)
     csection(1)  = csection(1)(1:ipos-1)
     cvarname(1)  = csection(1)
     clongname(1) = csection(1)
  ELSE
     CALL section_init(cf_section, csection,cvarname,clongname, iimina,iimaxa,ijmina,ijmaxa, nsection)
  ENDIF

  ! Allocate and build sigma levels and section array
  ALLOCATE ( dsigma_lev (nbins+1) , dtrpbin(nsection,nbins)  )

  dsigma_lev(1)=dsigma_min
  dltsig=( dsigma_max - dsigma_min) / nbins
  DO jclass =2, nbins+1
     dsigma_lev(jclass)= dsigma_lev(1) + (jclass-1) * dltsig
  END DO

  ! Look for vertical size of the domain
  npk = getdim (cf_tfil,cn_z)
  ALLOCATE ( gdept(npk), gdepw(npk) )
  IF ( lfull ) ALLOCATE ( e3t1d(npk), e3w1d(npk))

  ! read gdept, gdepw : it is OK even in partial cells, as we never use the bottom gdep
  IF ( .NOT. lbrk ) THEN
     gdept(:) = getvare3(cn_fzgr, cn_gdept, npk)
     gdepw(:) = getvare3(cn_fzgr, cn_gdepw, npk)
  ENDIF

  IF ( lfull )  THEN 
     e3t1d(:) = getvare3(cn_fzgr, cn_ve3t1d, npk)
     e3w1d(:) = getvare3(cn_fzgr, cn_ve3w1d, npk)
  ENDIF

  !! *  Main loop on sections
  DO jsec=1,nsection
     iimin=iimina(jsec) ; iimax=iimaxa(jsec)
     ijmin=ijmina(jsec) ; ijmax=ijmaxa(jsec)

     IF (iimin == iimax ) THEN        ! meridional section
        npts    = ijmax - ijmin       ! number of segments
        l_merid = .TRUE.

     ELSE IF ( ijmin == ijmax ) THEN  ! zonal section
        npts    = iimax - iimin       ! number of segments
        l_merid = .FALSE.

     ELSE
        PRINT *,' Section ',TRIM(csection(jsec)),' is neither zonal nor meridional :('
        PRINT *,' We skip this section .'
        CYCLE
     ENDIF

     ALLOCATE ( zu(npts,npk), zt(npts,npk), zs(npts,npk), zz(npts,npk), dsig(npts,0:npk)  )
     ALLOCATE ( eu(npts), de3(npts,npk), ddepu(npts, 0:npk), ddepw(npts,0:npk),zmask(npts,npk) )
     ALLOCATE ( tmpm(1,npts), tmpz(npts,1)                                   )
     ALLOCATE ( dwtrp(npts, nbins+1), dhiso(npts,nbins+1), dwtrpbin(npts,nbins) )
     ALLOCATE ( rlonlat(npts,1) )

     zt = 0. ; zs = 0. ; zu = 0. ; ddepu= 0.d0 ; zmask = 0.  ; dsig=0.d0

     IF (l_merid ) THEN   ! meridional section at i=iimin=iimax  ! use getvaryz
        tmpm(:,:)    = getvar(cn_fhgr, cn_ve2u,   1, 1, npts, kimin=iimin, kjmin=ijmin+1)
        eu(:)        = tmpm(1,:)  ! metrics varies only horizontally
        tmpm(:,:)    = getvar(cn_fhgr, cn_vlat2d, 1, 1, npts, kimin=iimin, kjmin=ijmin+1)
        rlonlat(:,1) = tmpm(1,:)  ! latitude in this case

        ! use zt and zs as temporaty variable for e3w
        IF ( lfull ) THEN
           DO ji=1, npts
              de3(ji,:) = e3t1d(:)
              zt( ji,:) = e3w1d(:)
              zs( ji,:) = e3w1d(:)
           ENDDO
        ELSE
            SELECT CASE ( cg_zgr_ver )
            CASE ( 'v2.0' ) ; cn_ve3u = 'e3u_ps' ; cn_ve3w='e3w_ps'
            CASE ( 'v3.0' ) ; cn_ve3u = 'e3u'    ; cn_ve3w='e3w'
            CASE ( 'v3.6' ) ; cn_ve3u = 'e3u_0'  ; cn_ve3w='e3w_0'
            END SELECT
           de3(:,:) = getvaryz(cn_fe3u, cn_ve3u, iimin,   npts, npk, kjmin=ijmin+1 )
           zt( :,:) = getvaryz(cn_fe3w, cn_ve3w, iimin,   npts, npk, kjmin=ijmin+1 )
           zs( :,:) = getvaryz(cn_fe3w, cn_ve3w, iimin+1, npts, npk, kjmin=ijmin+1 )
        ENDIF

        DO ji=1, npts
           ddepu(ji,1) = gdept(1) ! assume that first depth is set by gdept(1)
        ENDDO                    ! for vvl case, it might not be exact ...

        DO jk=2, npk
           DO ji=1,npts
              ddepu(ji,jk) = ddepu(ji,jk-1) +MIN (zt(ji,jk), zs(ji,jk) )
           ENDDO
        ENDDO
        ! normal velocity
        zu( :,:) = getvaryz(cf_ufil, cn_vozocrtx, iimin,   npts, npk, kjmin=ijmin+1 )
        WHERE( zu == zspu ) zu = 0.0

        ! salinity and deduce umask for the section
        zs( :,:) = getvaryz(cf_sfil, cn_vosaline, iimin,   npts, npk, kjmin=ijmin+1 )
        zt( :,:) = getvaryz(cf_sfil, cn_vosaline, iimin+1, npts, npk, kjmin=ijmin+1 )
        zmask = 1.
        WHERE ( zs == zsps .OR. zt == zsps ) zmask = 0.
        zs (:,:) = 0.5 * ( zs(:,:) + zt(:,:) )*zmask

        ! limitation to 'wet' points
        DO jk = 1, npk
           IF ( SUM(zs(:,jk)) == 0 ) THEN
              nk=jk
              EXIT
           ENDIF
        ENDDO
        PRINT *, ' NK = ',nk

        ! temperature
        zt(:,:) = getvaryz(cf_tfil, cn_votemper, iimin  , npts, npk, kjmin=ijmin+1 )
        zz(:,:) = getvaryz(cf_tfil, cn_votemper, iimin+1, npts, npk, kjmin=ijmin+1 )
        zt(:,:) = 0.5 * ( zt(:,:) + zz(:,:) ) *zmask

     ELSE                   ! zonal section at j=ijmin=ijmax ( include BRK-line sections)


        IF (lbrk ) THEN
           ! need to fix de3, ddepu, zu, zs, zt, nk , zmask
           IF ( lfull) THEN 
              DO ji=1, npts
                 de3(ji,:) = e3t1d(:)
              ENDDO
           ELSE
              tmpz(:,:)    = getvar(cn_fhgr, cn_ve1v,   1, npts, 1, kimin=iimin, kjmin=ijmin)
              eu(:)        = tmpz(:,1)
              tmpz(:,:)    = getvar(cn_fhgr, cn_vlon2d, 1, npts, 1, kimin=iimin, kjmin=ijmin)
              rlonlat(:,1) = tmpz(:,1)  ! longitude in this case
              de3(  :,:    ) = getvarxz(cf_brk,  cn_ve3v,     ijmin,   npts, npk, kimin=iimin )
              ddepu(:,1:npk) = getvarxz(cf_brk,  cn_depu3d,   ijmin,   npts, npk, kimin=iimin )
              ddepw(:,1:npk) = getvarxz(cf_brk,  cn_depw3d,   ijmin,   npts, npk, kimin=iimin )
           ENDIF
           zu(   :,:) = getvarxz(cf_vfil, cn_vomecrty, ijmin,   npts, npk, kimin=iimin )
           zt(   :,:) = getvarxz(cf_vfil, cn_votemper, ijmin,   npts, npk, kimin=iimin )
           zs(   :,:) = getvarxz(cf_vfil, cn_vosaline, ijmin,   npts, npk, kimin=iimin )
           zmask(:,:) = getvarxz(cf_vfil, cn_vmask,    ijmin,   npts, npk, kimin=iimin )
           ! limitation to 'wet' points
           ! JMM Broken line mask have a 9999 fil value (to be improved)
           WHERE(zmask == 9999.) zmask=0.

           nk = npk
           DO jk = 1, npk
              IF ( SUM(zmask(:,jk)) == 0 ) THEN
                 nk=jk
                 EXIT
              ENDIF
           ENDDO
           
        ELSE
        tmpz(:,:)    = getvar(cn_fhgr, cn_ve1v,   1, npts, 1, kimin=iimin, kjmin=ijmin)
        eu(:)        = tmpz(:,1)
        tmpz(:,:)    = getvar(cn_fhgr, cn_vlon2d, 1, npts, 1, kimin=iimin, kjmin=ijmin)
        rlonlat(:,1) = tmpz(:,1)  ! longitude in this case

           ! use zt and zs as temporary variable for e3w
           IF ( lfull ) THEN
              DO ji=1, npts
                 de3(ji,:) = e3t1d(:)
                 zt( ji,:) = e3w1d(:)
                 zs( ji,:) = e3w1d(:)
              ENDDO
           ELSE
              SELECT CASE ( cg_zgr_ver )
              CASE ( 'v2.0' ) ; cn_ve3v = 'e3v_ps' ; cn_ve3w='e3w_ps'
              CASE ( 'v3.0' ) ; cn_ve3v = 'e3v'    ; cn_ve3w='e3w'
              CASE ( 'v3.6' ) ; cn_ve3v = 'e3v_0'  ; cn_ve3w='e3w_0'
              END SELECT

              de3(:,:) = getvarxz(cn_fe3v, cn_ve3v, ijmin,   npts, npk, kimin=iimin+1 )
              zt( :,:) = getvarxz(cn_fe3w, cn_ve3w, ijmin,   npts, npk, kimin=iimin+1 )
              zs( :,:) = getvarxz(cn_fe3w, cn_ve3w, ijmin+1, npts, npk, kimin=iimin+1 )
           ENDIF

           DO ji=1, npts
              ddepu(ji,1) = gdept(1)
           ENDDO

           DO jk=2, npk
              DO ji=1,npts
                 ddepu(ji,jk) = ddepu(ji,jk-1) +MIN (zt(ji,jk), zs(ji,jk) )
              ENDDO
           ENDDO

           ! normal velocity
           zu( :,:) = getvarxz(cf_vfil, cn_vomecrty, ijmin,   npts, npk, kimin=iimin+1 )
           WHERE( zu == zspv ) zu = 0.0

           ! salinity and deduce umask for the section
           zs( :,:) = getvarxz(cf_sfil, cn_vosaline, ijmin,   npts, npk, kimin=iimin+1 )
           zt( :,:) = getvarxz(cf_sfil, cn_vosaline, ijmin+1, npts, npk, kimin=iimin+1 )

           zmask = 1.
           WHERE ( zs == zsps .OR. zt == zsps ) zmask = 0.
           zs (:,:) = 0.5 * ( zs(:,:) + zt(:,:) )* zmask

           ! limitation to 'wet' points
           DO jk = 1, npk
              IF ( SUM(zs(:,jk)) == 0 ) THEN
                 nk=jk
                 EXIT
              ENDIF
           ENDDO
          PRINT *,'JMM nk  ', nk

           ! temperature
           zt(:,:) = getvarxz(cf_tfil, cn_votemper, ijmin  , npts, npk, kimin=iimin+1 )
           zz(:,:) = getvarxz(cf_tfil, cn_votemper, ijmin+1, npts, npk, kimin=iimin+1 )
           zt(:,:) = 0.5 * ( zt(:,:) + zz(:,:) )
        ENDIF
     ENDIF

     ! compute density only for wet points
     IF ( lntr ) THEN 
        dsig(:,1:nk)=sigmantr( zt, zs,         npts, nk)*zmask(:,:)
     ELSE
        IF     ( refdep == -10.) THEN  ; dsig(:,1:nk)= -zt(:,:)*zmask(:,:)  ! change sign 
        ELSEIF ( refdep ==   0.) THEN  ; dsig(:,1:nk)=sigma0( zt, zs,         npts, nk)*zmask(:,:)
        ELSE                           ; dsig(:,1:nk)=sigmai( zt, zs, refdep, npts, nk)*zmask(:,:)
        ENDIF
     ENDIF

     dsig(:,0)=dsig(:,1)-1.e-4   ! dummy layer for easy interpolation

!! on each vertical column of water set density of land point to be the bottom most wet density
     DO ji = 1, npts
        DO jk = 1, nk
            IF ( zmask (ji,jk) == 0 ) THEN
               dsig(ji,jk)= dsig(ji,jk-1) +1.e-5  ! in order to have an increase
            ENDIF
        ENDDO
     ENDDO

     ! compute depth of isopynals (nbins+1 )
     DO  jiso =1, nbins+1
        dsigma=dsigma_lev(jiso)
!!!  REM : I and K loop can be inverted if necessary
        DO ji=1,npts
           !           dhiso(ji,jiso) = gdept(npk)  ! for broken line it is easier to use ddepu. Impact ?
           dhiso(ji,jiso) = ddepu(ji,npk)
           DO jk=1,nk 
              IF ( dsig(ji,jk) < dsigma ) THEN
              ELSE
                 ! interpolate between jk-1 and jk
                 dalfa=(dsigma - dsig(ji,jk-1)) / ( dsig(ji,jk) -dsig(ji,jk-1) )
                 IF (ABS(dalfa) > 1.1d0 .OR. dalfa < 0.d0 ) THEN   ! case dsig(0) = dsig(1)-1.e-4
                    dhiso(ji,jiso)= 0.d0
                 ELSE
                    dhiso(ji,jiso)= ddepu(ji,jk)*dalfa + (1.d0-dalfa)* ddepu(ji,jk-1)
                 ENDIF
                 EXIT
              ENDIF
           END DO
        END DO
     END DO

     ! compute transport between surface and isopycn 
     DO jiso = 1, nbins + 1
        dsigma=dsigma_lev(jiso)
        DO ji=1,npts
           dwtrp(ji,jiso) = 0.d0
           IF ( lbrk ) gdepw(:) = ddepw(ji,1:npk)   ! in broken line gdepw is not available as 1d array
           DO jk=1, nk-1
              IF ( gdepw(jk+1) < dhiso(ji,jiso) ) THEN
                 dwtrp(ji,jiso)= dwtrp(ji,jiso) + eu(ji)*de3(ji,jk)*zu(ji,jk)*1.d0
              ELSE  ! last box ( fraction)
                 dwtrp(ji,jiso)= dwtrp(ji,jiso) + eu(ji)*(dhiso(ji,jiso)-gdepw(jk))*zu(ji,jk)*1.d0
                 EXIT  ! jk loop
              ENDIF
           END DO
        END DO
     END DO

     ! binned transport : difference between 2 isopycns
     DO jbin=1, nbins
        dsigma=dsigma_lev(jbin)
        DO ji=1, npts
           dwtrpbin(ji,jbin) = dwtrp(ji,jbin+1) -  dwtrp(ji,jbin) 
        END DO
        dtrpbin(jsec,jbin)=SUM(dwtrpbin(:,jbin) )
     END DO

     ! output of the code for 1 section
     IF (lprint) CALL print_out(jsec)
     IF (lxtra ) CALL cdf_writ(jsec)
     PRINT *,' Total transport in all bins :',TRIM(csection(jsec)),' ',SUM(dtrpbin(jsec,:) )/1.d6

     ! free memory for the next section
     DEALLOCATE ( zu, zt, zs, zz, dsig, ddepu, ddepw, dhiso, dwtrp, dwtrpbin )
     DEALLOCATE ( eu, de3, tmpm, tmpz, zmask, rlonlat             )

  END DO   ! next section

  !! Global Output
  OPEN( numout, FILE=cf_out)
  ipos=INDEX(cf_tfil,'_gridT.nc')
  WRITE(numout,9006)  TRIM(cf_tfil(1:ipos-1))
  WRITE(numout,9005) ' sigma  ', (csection(jsec),jsec=1,nsection)
  DO jiso=1,nbins
     WRITE(numout,9004) dsigma_lev(jiso), (dtrpbin(jsec,jiso),jsec=1,nsection)
  ENDDO
  CLOSE(numout)

  cv_dep='levels'
  ! need to call section_init again in order to reset cvarname, clongname if they where modified 
  ! previously in cdf_writ(  case lxtra=true )
  IF (lxtra .AND. .NOT. lbrk) THEN
     CALL section_init(cf_section, csection,cvarname,clongname, iimina,iimaxa,ijmina,ijmaxa, nsection)
  ENDIF

  DO jsec=1,nsection
     ! setup output variables (section dependant for adaptative variable name (if possible)
     ! define new variables for output 

     CALL CreateOutput (jsec) 

     DO jiso=1,nbins
        rdummy1 = dsigma_lev(jiso)
        rdummy2 = dtrpbin(jsec,jiso)/1.d6  ! Sv
        ierr    = putvar(ncout, id_varout(1), rdummy1, jiso, ikx, iky )
        ierr    = putvar(ncout, id_varout(2), rdummy2, jiso, ikx, iky )
     END DO

     ierr = closeout(ncout)

  END DO

9004 FORMAT(f9.4, 20e16.7)
9005 FORMAT('#',a9, 20(2x,a12,2x) )
9006 FORMAT('# ',a)

CONTAINS

  SUBROUTINE section_init(cdfile, cdsection, cdvarname, cdlongname, kimin, kimax, kjmin, kjmax, knumber)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE section_init  ***
    !!
    !! ** Purpose : Read input ASCII file that defines section names and limit of
    !!              sections.
    !!
    !! ** Method  : At fisrt call only return the number of sections for further
    !!              allocation.  
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),                       INTENT(in   ) :: cdfile
    CHARACTER(LEN=256), DIMENSION(knumber), INTENT(out  ) :: cdsection
    CHARACTER(LEN=256), DIMENSION(knumber), INTENT(out  ) :: cdvarname
    CHARACTER(LEN=256), DIMENSION(knumber), INTENT(out  ) :: cdlongname
    INTEGER(KIND=4),                        INTENT(inout) :: knumber
    INTEGER(KIND=4), DIMENSION(knumber),    INTENT(out  ) :: kimin, kimax, kjmin, kjmax

    ! Local variables
    INTEGER(KIND=4)                                       :: jsec
    INTEGER(KIND=4)                                       :: ii, inum=10
    INTEGER(KIND=4)                                       :: ipos  
    CHARACTER(LEN=256)                                    :: cline
    CHARACTER(LEN=80), DIMENSION(3)                       :: cldum
    LOGICAL                                               :: llfirst
    !!----------------------------------------------------------------------
    llfirst=.FALSE.
    IF ( knumber == 0 ) llfirst=.TRUE.

    OPEN(inum, FILE=cdfile)
    REWIND(inum)
    ii = 0

    ! read the file just to count the number of sections
    DO
       READ(inum,'(a)') cline
       IF (INDEX(cline,'EOF') == 0 ) THEN
          READ(inum,*)    ! skip one line
          ii = ii + 1
       ELSE
          EXIT
       ENDIF
    END DO

    knumber=ii
    IF ( llfirst ) RETURN

    REWIND(inum)
    DO jsec=1,knumber
       READ(inum,'(a)') cline
       ii = 0
       cldum(:) = 'none'
       ipos = INDEX(cline,' ')
       DO WHILE ( ipos > 1 ) 
          ii = ii + 1
          cldum(ii) = cline(1:ipos - 1 )
          cline = TRIM ( cline(ipos+1:) )
          ipos  = INDEX( cline,' ' ) 
          IF ( ii >= 3 ) EXIT
       END DO
       cdsection(jsec) = TRIM(cldum(1) )
       cdvarname(jsec) = TRIM(cldum(2) )
       cdlongname(jsec) = TRIM(cldum(3) )
       READ(inum,*    ) kimin(jsec), kimax(jsec), kjmin(jsec), kjmax(jsec)
    END DO

    CLOSE(inum)

  END SUBROUTINE section_init

  SUBROUTINE cdf_writ( ksec)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE cdf_writ  ***
    !!
    !! ** Purpose :  Write output cdf files if required 
    !!
    !! ** Method  :  Most of the variables are global 
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),   INTENT(in) :: ksec  ! number of the section

    INTEGER(KIND=4)               :: ji, jk
    INTEGER(KIND=4)               :: ivar
    INTEGER(KIND=4)               :: icout
    INTEGER(KIND=4), DIMENSION(4) :: ipk, id_varout

    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zdum
    TYPE(variable),  DIMENSION(4) :: sl_typvar
    CHARACTER(LEN=255)            :: csuffixvarnam
    CHARACTER(LEN=255)            :: cprefixlongnam
    !!----------------------------------------------------------------------
    IF ( cvarname(ksec) /= 'none' ) THEN
       csuffixvarnam = '_'//TRIM(cvarname(ksec))
    ELSE
       csuffixvarnam = ''
    ENDIF
    IF ( clongname(ksec) /= 'none' ) THEN
       cprefixlongnam = TRIM(clongname(ksec))//'_'
    ELSE
       cprefixlongnam = ''
    ENDIF

    ALLOCATE ( zdum(npts,1))
    ! (along section, depth ) 2D variables
    cf_nc=TRIM(csection(ksec))//'_secdep.nc'
    ! define variables
    ipk(:)=nk
    sl_typvar%rmissing_value    = 0.
    sl_typvar%rmissing_value    = 0.
    sl_typvar%scale_factor      = 1.
    sl_typvar%add_offset        = 0.
    sl_typvar%savelog10         = 0.
    sl_typvar%iwght             = iweight
    sl_typvar%conline_operation = 'N/A'
    sl_typvar%caxis             = 'XZT'

    ivar=1
    sl_typvar(ivar)%cname          = 'temperature'//TRIM(csuffixvarnam)
    sl_typvar(ivar)%cunits         = 'Celsius'
    sl_typvar(ivar)%valid_min      = -2.
    sl_typvar(ivar)%valid_max      = 45.
    sl_typvar(ivar)%clong_name     = TRIM(cprefixlongnam)//'Potential_temperature'
    sl_typvar(ivar)%cshort_name    = 'temperature'

    ivar=ivar+1
    sl_typvar(ivar)%cname          = 'salinity'//TRIM(csuffixvarnam)
    sl_typvar(ivar)%cunits         = 'PSU'
    sl_typvar(ivar)%valid_min      = 0.
    sl_typvar(ivar)%valid_max      = 45.
    sl_typvar(ivar)%clong_name     = TRIM(cprefixlongnam)//'Salinity'
    sl_typvar(ivar)%cshort_name    = 'salinity'

    ivar=ivar+1
    sl_typvar(ivar)%cname          = 'density'//TRIM(csuffixvarnam)
    sl_typvar(ivar)%cunits         = 'kg/m3 -1000'
    sl_typvar(ivar)%valid_min      = 0.
    sl_typvar(ivar)%valid_max      = 45.
    sl_typvar(ivar)%clong_name     = TRIM(cprefixlongnam)//'potential_density'
    sl_typvar(ivar)%cshort_name    = 'density'

    ivar=ivar+1
    sl_typvar(ivar)%cname          = 'velocity'//TRIM(csuffixvarnam)
    sl_typvar(ivar)%cunits         = 'm/s'
    sl_typvar(ivar)%valid_min      = -3.
    sl_typvar(ivar)%valid_max      = 3.
    sl_typvar(ivar)%clong_name     = TRIM(cprefixlongnam)//'Normal_velocity'
    sl_typvar(ivar)%cshort_name    = 'velocity'

    icout = create      (cf_nc, 'none',    npts, 1, nk, cdep=cn_vdeptht                     )
    ierr  = createvar   (icout, sl_typvar, ivar, ipk, id_varout, cdglobal=TRIM(cglobal)     )
    ierr  = putheadervar(icout, cf_tfil,   npts, 1, nk, &
         &   pnavlon=rlonlat, pnavlat=rlonlat, cdep=cn_vdeptht                              )

    !    dtim = getvar1d(cf_tfil, cn_vtimec, 1     )
    !    ierr = putvar1d(icout,   dtim,     1 , 'T')

    DO jk = 1, nk
       zdum(:,1)=zt(:,jk)   ; ierr = putvar ( icout, id_varout(1), zdum, jk, npts, 1 )
       zdum(:,1)=zs(:,jk)   ; ierr = putvar ( icout, id_varout(2), zdum, jk, npts, 1 )
       zdum(:,1)=dsig(:,jk) ; ierr = putvar ( icout, id_varout(3), zdum, jk, npts, 1 )
       zdum(:,1)=zu(:,jk)   ; ierr = putvar ( icout, id_varout(4), zdum, jk, npts, 1 )
    END DO

    ierr = closeout(icout)

    ! (along section, sigma ) 2D variables
    cf_nc=TRIM(csection(ksec))//'_secsig.nc'
    ! define variables
    ipk(:)=nbins
    sl_typvar%rmissing_value    = 99999.
    sl_typvar%rmissing_value    = 99999.
    sl_typvar%scale_factor      = 1.
    sl_typvar%add_offset        = 0.
    sl_typvar%savelog10         = 0.
    sl_typvar%iwght             = iweight
    sl_typvar%conline_operation = 'N/A'
    sl_typvar%caxis             = 'XST'

    ivar=1
    ipk(ivar)=nbins-1
    sl_typvar(ivar)%cname          = 'isodep'//TRIM(csuffixvarnam)
    sl_typvar(ivar)%cunits         = 'm'
    sl_typvar(ivar)%valid_min      = 0.
    sl_typvar(ivar)%valid_max      = 6000.
    sl_typvar(ivar)%clong_name     = TRIM(cprefixlongnam)//'isopycnal_depth'
    sl_typvar(ivar)%cshort_name    = 'isodep'

    ivar=ivar+1
    sl_typvar(ivar)%cname          = 'bintrp'//TRIM(csuffixvarnam)
    sl_typvar(ivar)%cunits         = 'SV'
    sl_typvar(ivar)%valid_min      = -5.
    sl_typvar(ivar)%valid_max      = 5.
    sl_typvar(ivar)%clong_name     = TRIM(cprefixlongnam)//'Binned_transport'
    sl_typvar(ivar)%cshort_name    = 'bintrp'

    ivar=ivar+1
    sl_typvar(ivar)%cname          = 'sumtrp'//TRIM(csuffixvarnam)
    sl_typvar(ivar)%cunits         = 'SV'
    sl_typvar(ivar)%valid_min      = -20.
    sl_typvar(ivar)%valid_max      = 20.
    sl_typvar(ivar)%clong_name     = TRIM(cprefixlongnam)//'cumulated_transport'
    sl_typvar(ivar)%cshort_name    = 'sumtrp'

    icout = create      (cf_nc, 'none',    npts, 1, nbins, cdep='levels'                 )
    ierr  = createvar   (icout, sl_typvar, ivar, ipk, id_varout, cdglobal=TRIM(cglobal)  )
    ierr  = putheadervar(icout, cf_tfil,   npts, 1, nbins, &
         &   pnavlon=rlonlat, pnavlat=rlonlat, pdep=REAL(dsigma_lev), cdep='levels'      )

    PRINT *, 'NBINS', nbins, npts
    DO jk = 1, nbins-1
       zdum(:,1)=dhiso   (:,jk)      ; ierr = putvar ( icout, id_varout(1), zdum, jk, npts, 1 )
    END DO
    DO jk = 1, nbins
       zdum(:,1)=dwtrpbin(:,jk)/1.e6 ; ierr = putvar ( icout, id_varout(2), zdum, jk, npts, 1 )
       zdum(:,1)=dwtrp   (:,jk)/1.e6 ; ierr = putvar ( icout, id_varout(3), zdum, jk, npts, 1 )
    END DO
    ierr = closeout(icout)

    DEALLOCATE ( zdum )

  END SUBROUTINE cdf_writ

  SUBROUTINE print_out(ksec)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE print_out  ***
    !!
    !! ** Purpose :  Print results on standard output 
    !!
    !! ** Method  :  Most of the variables are global and already known 
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) :: ksec  ! number of the section

    INTEGER(KIND=4)             :: ji, jk, jiso, jbin
    !!----------------------------------------------------------------------
    WRITE(cfmt_9000,'(a,i4,a)') '(i7,  ',npts,'f8.3)'
    WRITE(cfmt_9001,'(a,i4,a)') '(i7,  ',npts,'f8.0)'
    WRITE(cfmt_9002,'(a,i4,a)') '(f7.3,',npts,'f8.0)'
    WRITE(cfmt_9003,'(a,i4,a)') '(f7.3,',npts,'f8.3)'
    PRINT *,' T (deg C)' 
    DO jk=1,nk
       PRINT cfmt_9000, jk,  (zt(ji,jk),ji=1,npts)
    END DO

    PRINT *,' S (PSU)'
    DO jk=1,nk
       PRINT cfmt_9000,  jk,  (zs(ji,jk),ji=1,npts)
    END DO

    PRINT *,' SIG (kg/m3 - 1000 )'
    DO jk=1,nk
       PRINT cfmt_9000, jk,  (dsig(ji,jk),ji=1,npts)
    END DO

    PRINT *,' VELOCITY (cm/s ) '
    DO jk=1,nk
       PRINT cfmt_9000, jk,  (zu(ji,jk)*100,ji=1,npts)
    END DO

    PRINT *,' GDEPU (m) '
    DO jk=1,nk
       PRINT cfmt_9001,jk,  (ddepu(ji,jk)*zmask(ji,jk),ji=1,npts)
    END DO

    PRINT *, 'E3 (m)'
    DO jk=1,nk
       PRINT cfmt_9001,jk,  (de3(ji,jk)*zmask(ji,jk),ji=1,npts)
    END DO

    PRINT *,' DEP ISO ( m )'
    DO  jiso =1, nbins+1
       PRINT cfmt_9002, dsigma_lev(jiso),(dhiso(ji,jiso),ji=1,npts)
    END DO

    PRINT *,' TRP SURF -->  ISO (SV)'
    DO  jiso =1, nbins+1
       PRINT  cfmt_9003, dsigma_lev(jiso),(dwtrp(ji,jiso)/1.d6,ji=1,npts)
    END DO

    PRINT *,' TRP bins (SV)'
    DO  jbin =1, nbins
       PRINT  cfmt_9003, dsigma_lev(jbin),(dwtrpbin(ji,jbin)/1.d6,ji=1,npts), dtrpbin(ksec,jbin)/1.d6
    END DO
  END SUBROUTINE print_out

  SUBROUTINE CreateOutput(ksec)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) :: ksec  ! section index
    !!----------------------------------------------------------------------

    IF ( cvarname(ksec) /= 'none' ) THEN
       csuffixvarname='_'//TRIM(cvarname(ksec))
    ELSE
       csuffixvarname=''
    ENDIF
    IF ( clongname(ksec) /= 'none' ) THEN
       cprefixlongname=TRIM(clongname(ksec))//'_'
    ELSE
       cprefixlongname=''
    ENDIF

    stypvar%rmissing_value    = 99999.
    stypvar%scale_factor      = 1.
    stypvar%add_offset        = 0.
    stypvar%savelog10         = 0.
    stypvar%iwght             = iweight
    stypvar%conline_operation = 'N/A'
    stypvar%caxis             = 'ZT'

    IF ( ltemp ) THEN
       stypvar(1)%cname          = 'temp_class'
       stypvar(1)%cunits         = '[]'
       stypvar(1)%valid_min      = 0.
       stypvar(1)%valid_max      = 100.
       stypvar(1)%clong_name     = 'class of potential temperature'
       stypvar(1)%cshort_name    = 'temp_class'

       stypvar(2)%cname          = 'temptrp'//TRIM(csuffixvarname)
       stypvar(2)%cunits         = 'Sv'
       stypvar(2)%valid_min      = -1000.
       stypvar(2)%valid_max      = 1000.
       stypvar(2)%clong_name     = TRIM(cprefixlongname)//'transport in temperature class'
       stypvar(2)%cshort_name    = 'temptrp'
    ELSE
       stypvar(1)%cname          = 'sigma_class'
       stypvar(1)%cunits         = '[]'
       stypvar(1)%valid_min      = 0.
       stypvar(1)%valid_max      = 100.
       stypvar(1)%clong_name     = 'class of potential density'
       stypvar(1)%cshort_name    = 'sigma_class'

       stypvar(2)%cname          = 'sigtrp'//TRIM(csuffixvarname)
       stypvar(2)%cunits         = 'Sv'
       stypvar(2)%valid_min      = -1000.
       stypvar(2)%valid_max      = 1000.
       stypvar(2)%clong_name     = TRIM(cprefixlongname)//'transport in sigma class'
       stypvar(2)%cshort_name    = 'sigtrp'
    ENDIF

    ! create output fileset
    IF (ltemp) THEN  ; cf_outnc = TRIM(csection(ksec))//'_trptemp.nc'
    ELSE             ; cf_outnc = TRIM(csection(ksec))//'_trpsig.nc'
    ENDIF

    ncout = create      (cf_outnc, 'none',  ikx,      iky, nbins, cdep=cv_dep               )
    ierr  = createvar   (ncout,    stypvar, nboutput, ipk, id_varout, cdglobal=TRIM(cglobal))
    ierr  = putheadervar(ncout,    cf_tfil, ikx,      iky, nbins, &
         &   pnavlon=rdumlon, pnavlat=rdumlat, pdep=REAL(dsigma_lev), cdep=cv_dep           )

    dtim = getvar1d(cf_tfil, cn_vtimec, 1     )
    ierr = putvar1d(ncout,   dtim,      1, 'T')

  END SUBROUTINE CreateOutput

  SUBROUTINE file_example
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE file_example  ***
    !!
    !! ** Purpose :  give an example for the dens_section.dat 
    !!
    !!----------------------------------------------------------------------
    PRINT *,'  '
    PRINT *,'  '
    PRINT *,'  '
    PRINT *,'   EXAMPLE of dens_section.dat file'
    PRINT *,'   --------------------------------'
    PRINT *,'  Each section is described by 2 lines :'
    PRINT *,'    line#1 : name_of_section [variable_suffix] [ long_name_prefix ]'
    PRINT *,'    line#2 : imin  imax jmin jmax '
    PRINT *,'     IMPORTANT : the points indicated by imin,jmin  imax,jmax are F-point.'
    PRINT *,'  '
    PRINT *,'  example (taken for ORCA12 configuration) :'
    PRINT *,'01_Denmark_strait denma Denmark_Strait_transport_in_sigma_classes'
    PRINT *,'3043 3162 2496 2496'
    PRINT *,'02_Faoes_Bank_Channel faroe Faroes_Bank_Channel_transport_in_sigma_classes'
    PRINT *,'3318 3318 2398 2420'
    PRINT *,'03_Gibraltar gibra Gibraltar_Strait_transport_in_sigma_classes'
    PRINT *,'3378 3378 1956 1961'
    PRINT *,'04_Gibbs_FZ gibbs Gibbs_Fracture_Zone_transport_in_sigma_classes'
    PRINT *,'3075 3075 2215 2235'
    PRINT *,'  '
    PRINT *,'   with this example there will be 4 output files named :'
    PRINT *,'   01_Denmark_strait_trpsig.nc'
    PRINT *,'   02_Faoes_Bank_Channel_trpsig.nc'
    PRINT *,'   03_Gibraltar.nc'
    PRINT *,'   04_Gibbs_FZ'
    PRINT *,' and taking 03_Gibraltar.nc as example, the variable will be (from ncdump):'
    PRINT *,'  '
    PRINT *,'....  '
    PRINT *,' float sigtrp_gibra(time_counter, levels, y, x) ;'
    PRINT *,' 	sigtrp_gibra:units = "Sv" ;'
    PRINT *,' 	sigtrp_gibra:_FillValue = 99999.f ;'
    PRINT *,' 	sigtrp_gibra:valid_min = -1000.f ;'
    PRINT *,' 	sigtrp_gibra:valid_max = 1000.f ;'
    PRINT *,' 	sigtrp_gibra:long_name = "Gibraltar_Strait_transport_in_sigma_classes_transport in sigma class" ;'
    PRINT *,' 	sigtrp_gibra:short_name = "sigtrp" ;'
    PRINT *,' 	sigtrp_gibra:iweight = 73 ;'
    PRINT *,' 	sigtrp_gibra:online_operation = "N/A" ;'
    PRINT *,' 	sigtrp_gibra:axis = "ZT" ;'
    PRINT *,' 	sigtrp_gibra:scale_factor = 1.f ;'
    PRINT *,' 	sigtrp_gibra:add_offset = 0.f ;'
    PRINT *,' 	sigtrp_gibra:savelog10 = 0.f ;'
    PRINT *,'....  '
    PRINT *, ' '
    PRINT *,'  If no variable name suffix nor longname prefix is given, variable names are '
    PRINT *,'  the same for all sections, preventing the concatenation of all the output '
    PRINT *,'  files into a single file.'
    PRINT *,'  Note that in the given example, the long name is not particularly well chosen,'
    PRINT *,'  as ''_transport in sigma class'' is automatically append to the prefix ! :) '
    PRINT *,'  '
    PRINT *,'  If you are using the DRAKKAR monitoring tools (DMONTOOLS), consider the use'
    PRINT *,'  of the create_sections_list.ksh script to build the dens_section.dat file '
    PRINT *,'  from a configuration data base called drakkar_sections_table.txt.'
    PRINT *,'  '
    PRINT *,'  '
    PRINT *,'  '

  END SUBROUTINE file_example

END PROGRAM cdfsigtrp
