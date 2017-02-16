PROGRAM cdfsigtrp_broken
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
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   routines      : description
   !!  section_init   : initialize section names and positions
   !!  print_out      : routine which performs standard output if required
   !!  bimg_writ      : routine which performs bimg output if required
   !!----------------------------------------------------------------------
   USE cdfio
   USE eos          ! for sigma0, sigmai
   USE modcdfnames  ! for ReadCdfNames
   USE modutils     ! for SetGlobalAtt
   !!----------------------------------------------------------------------
   !! CDFTOOLS_3.0 , MEOM 2011
   !! $Id: cdfsigtrp.f90 699 2013-06-24 14:17:21Z molines $
   !! Copyright (c) 2011, J.-M. Molines
   !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
   !!----------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER(KIND=4)                               :: ji, jk, jclass, jsec ! dummy loop index
   INTEGER(KIND=4)                               :: jiso, jbin, jarg     ! dummy loop index
   INTEGER(KIND=4)                               :: nbins                ! number of density classes
   INTEGER(KIND=4)                               :: ipos                 ! working variable
   INTEGER(KIND=4)                               :: narg, iargc          ! command line 
   INTEGER(KIND=4)                               :: ijarg, ireq          ! command line
   INTEGER(KIND=4)                               :: npk, nk              ! vertical size, number of wet layers
   INTEGER(KIND=4)                               :: numbimg=10           ! optional bimg logical unit
   INTEGER(KIND=4)                               :: numout=11            ! ascii output
   INTEGER(KIND=4)                               :: nsection             ! number of sections (overall)
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
   REAL(KIND=4), DIMENSION(1)                    :: tim                  ! time counter
   REAL(KIND=4), DIMENSION(1)                    :: rdummy1, rdummy2     ! working variable
   REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: gdept, gdepw         ! depth of T and W points 
   REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: eu                   ! either e1v or e2u
   REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: e3t1d, e3w1d         ! vertical metrics in case of full step
   REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: rlonlat              ! longitudes/latitudes if the section
   REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zs, zt               ! salinity and temperature from file 
   REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: rdumlon, rdumlat     ! dummy longitude and latitude for output
   REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zu                   ! velocity
   REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zmask                ! mask
   REAL(KIND=4), DIMENSION(:,:,:),   ALLOCATABLE :: tmpm, tmpz           ! temporary arrays
   REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: tmpz1d               ! temporary arrays
   ! double precision for cumulative variables and densities
   REAL(KIND=8)                                  :: dsigma_min           ! minimum density for bining
   REAL(KIND=8)                                  :: dsigma_max, dltsig   ! maximum density for bining, step
   REAL(KIND=8)                                  :: dsigma, dalfa        ! working sigma, interpolation coeff.
   REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dsigma_lev           ! built array with sigma levels
   REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: de3                  ! vertical metrics
   REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: ddepu                ! depth of vel points
   REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dsig                 ! density
   REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dhiso                ! depth of isopycns
   REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dwtrp, dwtrpbin      ! transport arrays
   REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dtrpbin              ! transport arrays

   TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvar              ! structure of output

   CHARACTER(LEN=2048)                            :: cf_tfil              ! temperature salinity velocity file
   !CHARACTER(LEN=2048)                            :: cf_ufil              ! zonal velocity file
   !CHARACTER(LEN=2048)                            :: cf_vfil              ! meridional velocity file
   CHARACTER(LEN=2048)                            :: cf_section='dens_section.dat'  ! input section file
   CHARACTER(LEN=2048)                            :: cf_out='trpsig.txt'  ! output  ascii file
   CHARACTER(LEN=2048)                            :: cf_bimg              ! output bimg file (2d)
   CHARACTER(LEN=2048)                            :: cf_nc                ! output netcdf file (2d)
   CHARACTER(LEN=2048)                            :: cf_outnc             ! output netcdf file (1d, 0d))
   CHARACTER(LEN=2048)                            :: cv_dep               ! depth variable
   CHARACTER(LEN=2048)                            :: cldum                ! dummy string
   CHARACTER(LEN=2048)                            :: cglobal              ! global attribute
   CHARACTER(LEN=80 )                            :: cfmt_9000            ! format string 
   CHARACTER(LEN=80 )                            :: cfmt_9001            ! format string
   CHARACTER(LEN=80 )                            :: cfmt_9002            ! format string
   CHARACTER(LEN=80 )                            :: cfmt_9003            ! format string
   CHARACTER(LEN=2048)                            :: cl_vnam, cl_lname    ! working variables
   CHARACTER(LEN=2048)                            :: csuffixvarname       !
   CHARACTER(LEN=2048)                            :: cprefixlongname      !
   CHARACTER(LEN=2048), DIMENSION(:), ALLOCATABLE :: cv_names             ! names of input variables
   CHARACTER(LEN=2048), DIMENSION(:), ALLOCATABLE :: csection             ! section name
   CHARACTER(LEN=2048), DIMENSION(:), ALLOCATABLE :: cvarname             ! output variable name (root)
   CHARACTER(LEN=2048), DIMENSION(:), ALLOCATABLE :: clongname            ! output long name (root)

   LOGICAL                                       :: l_merid              ! flag for meridional section
   LOGICAL                                       :: ltemp  =.FALSE.      ! flag for use of temperature
   LOGICAL                                       :: lprint =.FALSE.      ! flag for extra print
   LOGICAL                                       :: lbimg  =.FALSE.      ! flag for bimg output
   LOGICAL                                       :: lncdf  =.FALSE.      ! flag for extra netcdf output
   LOGICAL                                       :: lfull  =.FALSE.      ! flag for full step 
   LOGICAL                                       :: lneutral  =.FALSE.   ! flag for neutral density
   LOGICAL                                       :: lchk   =.FALSE.      ! flag for missing files
   !!----------------------------------------------------------------------
   CALL ReadCdfNames()

   narg= iargc()
   IF ( narg < 4 ) THEN
      PRINT *,' usage :  cdfsigtrp_broken TSV-file sigma_min sigma_max nbins ...'
      PRINT *,'              ... [-print ] [-bimg ] [-full ] [ -refdep ref_depth] ...'
      PRINT *,'              ... [-neutral ] [-section file ] [-temp ]'
      PRINT *,'      '
      PRINT *,'     PURPOSE :'
      PRINT *,'       Compute density class transports, according to the density class' 
      PRINT *,'       definition ( minimum, maximum and number of bins) given in arguments.'
      PRINT *,'       Section position are given in ',TRIM(cf_section),', an ASCII file '
      PRINT *,'       with pairs of lines giving section name and section location as'
      PRINT *,'       imin imax jmin jmax. Only zonal or meridional section are allowed.'
      PRINT *,'       The name of this file can be specified with the -section option, if'
      PRINT *,'       it differs from the standard name. Optionaly, a netcdf root variable '
      PRINT *,'       name and a netcdf root long-name can be provided on the line giving '
      PRINT *,'       the section name.'
      PRINT *,'       Pedro: The section used is the result of cdf_xtrac_broken. In this way,'
      PRINT *,'       it is possible to calculate the transport of density class in'
      PRINT *,'       oblicous sections (non dependance on zonal or meridional).'
      PRINT *,'      '
      PRINT *,'       This program can also be used to compute transport by class of '
      PRINT *,'       temperatures, provided the temperatures decrease monotonically '
      PRINT *,'       downward. In this case, use -temp option and of course specify'
      PRINT *,'       sigma_min, sigma_max as temperatures.'
      PRINT *,'      '
      PRINT *,'     ARGUMENTS :'
      PRINT *,'       TSV-file : netcdf_broken_line file with temperature, salinity' 
      PRINT *,'       and the normal velocity through the section'
      PRINT *,'       sigma_min : minimum density for binning'
      PRINT *,'       sigma_max : maximum density for binning'
      PRINT *,'       nbins : number of bins. This will fix the bin ''width'' '
      PRINT *,'      '
      PRINT *,'     OPTIONS :'
      PRINT *,'       [ -full ] : for full step configuration' 
      PRINT *,'       [ -bimg ] : produce extra bimg output file which shows the details'
      PRINT *,'               of the sections (normal velocity, density, temperature, '
      PRINT *,'               salinity, transports, isopycnal depths. (to be change to '
      PRINT *,'               netcdf files for more common use.'
      PRINT *,'       [ -ncdf ] : produce extra netcdf output file which shows the details'
      PRINT *,'               of the sections (normal velocity, density, temperature, '
      PRINT *,'               salinity, transports, isopycnal depths. '
      PRINT *,'       [ -print ]: write the binned transports on standard output, for each'
      PRINT *,'               sections.'
      PRINT *,'       [ -refdep ref_depth ]: give a reference depths for the computation of'
      PRINT *,'               potential density. Sigma_min, sigma_max must be adapted '
      PRINT *,'               accordingly.'
      PRINT *,'       [ -neutral ]: use neutral density instead of potential density '
      PRINT *,'       [ -section file] : give the name of section file.'
      PRINT *,'               Default is ', TRIM(cf_section)
      PRINT *,'       [ -temp ] : use temperature instead of density for binning'
      PRINT *,'      '
      PRINT *,'     REQUIRED FILES :'
      PRINT *,'       ', TRIM(cn_fhgr),', ', TRIM(cn_fzgr),' and ', TRIM(cf_section)
      PRINT *,'      '
      PRINT *,'     OUTPUT : '
      PRINT *,'       Netcdf file : There is 1 netcdf file per section. File name is build'
      PRINT *,'         from section name : Section_name_trpsig.nc'
      PRINT *,'         variables : sigma_class (upper limit of the bin)'
      PRINT *,'                     sigtrp : transport (Sv per bin)'
      PRINT *,'      '
      PRINT *,'       ascii file  : ', TRIM(cf_out) 
      PRINT *,'      '
      PRINT *,'       bimg  file  :  There are 2 bimg files whose name is build from section'
      PRINT *,'         name : section_name_trpdep.bimg and section_name_trpsig.bimg.'
      PRINT *,'         This file is written only if -bimg option is used.'
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
   ijarg = 1 ; ireq = 0
   DO WHILE ( ijarg <= narg )
      CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
      SELECT CASE ( cldum )
      CASE ( '-full' ) ; lfull  = .TRUE.
      CASE ( '-bimg' ) ; lbimg  = .TRUE.
      CASE ( '-ncdf' ) ; lncdf  = .TRUE.
      CASE ( '-print') ; lprint = .TRUE.
      CASE ( '-temp')  ; ltemp  = .TRUE. 
      CASE ( '-refdep' ) ; CALL getarg(ijarg, cldum      ) ; ijarg=ijarg+1 ; READ(cldum,*) refdep
      CASE ( '-section') ; CALL getarg(ijarg, cf_section ) ; ijarg=ijarg+1 
      CASE ( '-neutral') ; lneutral = .TRUE.
      CASE DEFAULT
         ireq=ireq+1
         SELECT CASE ( ireq)
         CASE ( 1 ) ; cf_tfil = cldum
         !CASE ( 2 ) ; cf_ufil = cldum
         !CASE ( 3 ) ; cf_vfil = cldum
         CASE ( 2 ) ; READ(cldum,*) dsigma_min
         CASE ( 3 ) ; READ(cldum,*) dsigma_max
         CASE ( 4 ) ; READ(cldum,*) nbins
         CASE DEFAULT 
            PRINT *,' Too many arguments ' ; STOP
         END SELECT
      END SELECT
   END DO
   !Commented Pedro
   ! check for file existence
   lchk = lchk .OR. chkfile( cn_fzgr    )
   lchk = lchk .OR. chkfile( cn_fhgr    )
   lchk = lchk .OR. chkfile( cf_section )
   lchk = lchk .OR. chkfile( cf_tfil    )
   !lchk = lchk .OR. chkfile( cf_ufil    )
   !lchk = lchk .OR. chkfile( cf_vfil    )
   IF ( lchk ) STOP ! missing file
   IF ( ltemp)  THEN  ! temperature decrease downward. Change sign and swap min/max
      refdep = -10. ! flag value
      dltsig     = dsigma_max  ! use dltsig as dummy variable for swapping
      dsigma_max = -dsigma_min
      dsigma_min = -dltsig
   ENDIF

   ! define global attribute with command line
   CALL SetGlobalAtt( cglobal)

   !! Commented Pedro   
   ! get the attribute iweight from vozocrtx
   !iweight = getatt(cf_ufil, cn_vozocrtx, 'iweight')
   iweight = 1
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
   nsection = 0 ; CALL section_init(cf_section, csection,cvarname,clongname,iimina, iimaxa, ijmina, ijmaxa, nsection)
   ALLOCATE ( csection(nsection), cvarname(nsection), clongname(nsection) )
   ALLOCATE ( iimina(nsection), iimaxa(nsection), ijmina(nsection),ijmaxa(nsection) )
   CALL section_init(cf_section, csection,cvarname,clongname, iimina,iimaxa,ijmina,ijmaxa, nsection)

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
   gdept(:) = getvare3(cn_fzgr, cn_gdept, npk)
   gdepw(:) = getvare3(cn_fzgr, cn_gdepw, npk)

   IF ( lfull )  THEN 
      e3t1d(:) = getvare3(cn_fzgr, cn_ve3t, npk)
      e3w1d(:) = getvare3(cn_fzgr, cn_ve3w, npk)
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

      ALLOCATE ( zu(npts,npk), zt(npts,npk), zs(npts,npk), dsig(npts,0:npk)  )
      ALLOCATE ( eu(npts), de3(npts,npk), ddepu(npts, 0:npk), zmask(npts,npk) )
      ALLOCATE ( tmpm(1,npts,2), tmpz(npts,1,2)                              )
      ALLOCATE ( dwtrp(npts, nbins+1), dhiso(npts,nbins+1), dwtrpbin(npts,nbins) )
      ALLOCATE ( rlonlat(npts,1) )


      IF (l_merid ) THEN   ! meridional section at i=iimin=iimax
         EXIT
      !   tmpm(:,:,1)   = getvar(cn_fhgr, cn_ve2u,   1, 1, npts, kimin=iimin, kjmin=ijmin+1)
      !   eu(:)         = tmpm(1,:,1)  ! metrics varies only horizontally
      !   tmpm(:,:,1)   = getvar(cn_fhgr, cn_vlat2d, 1, 1, npts, kimin=iimin, kjmin=ijmin+1)
      !   rlonlat(:,1)  = tmpm(1,:,1)  ! latitude in this case
      !   DO jk = 1,npk
      !      ! initiliaze ddepu to gdept()
      !      ddepu(:,jk) = gdept(jk)

      !      IF ( lfull ) THEN
      !         de3(:,jk)   = e3t1d(jk)
      !         tmpm(1,:,1) = e3w1d(jk)
      !         tmpm(1,:,2) = e3w1d(jk)
      !      ELSE
      !         ! vertical metrics (PS case)
      !         tmpm(:,:,1) = getvar(cf_ufil, 'e3u_ps', jk, 1, npts, kimin=iimin,   kjmin=ijmin+1, ldiom=.TRUE.)
      !         de3(:,jk)   = tmpm(1,:,1)
      !         tmpm(:,:,1) = getvar(cf_ufil, 'e3w_ps', jk, 1, npts, kimin=iimin,   kjmin=ijmin+1, ldiom=.TRUE.)
      !         tmpm(:,:,2) = getvar(cf_ufil, 'e3w_ps', jk, 1, npts, kimin=iimin+1, kjmin=ijmin+1, ldiom=.TRUE.)
      !      ENDIF

      !      IF (jk >= 2 ) THEN
      !         DO ji=1,npts
      !            ddepu(ji,jk)= ddepu(ji,jk-1) + MIN(tmpm(1,ji,1), tmpm(1,ji,2))
      !         END DO
      !      ENDIF

      !      ! Normal velocity
      !      tmpm(:,:,1) = getvar(cf_ufil,cn_vozocrtx,jk,1,npts, kimin=iimin, kjmin=ijmin+1)
      !      zu(:,jk)    = tmpm(1,:,1)

      !      ! salinity and deduce umask for the section
      !      tmpm(:,:,1) = getvar(cf_tfil,cn_vosaline,jk,1,npts, kimin=iimin  , kjmin=ijmin+1)
      !      tmpm(:,:,2) = getvar(cf_tfil,cn_vosaline,jk,1,npts, kimin=iimin+1, kjmin=ijmin+1)
      !      zmask(:,jk) = tmpm(1,:,1)*tmpm(1,:,2)
      !      WHERE ( zmask(:,jk) /= 0 ) zmask(:,jk)=1
            ! do not take special care for land value, as the corresponding velocity point is masked
      !      zs(:,jk) = 0.5 * ( tmpm(1,:,1) + tmpm(1,:,2) )

      !      ! limitation to 'wet' points
      !      IF ( SUM(zs(:,jk))  == 0 ) THEN
      !         nk=jk ! first vertical point of the section full on land
      !         EXIT  ! as soon as all the points are on land
      !      ENDIF

      !      ! temperature
      !      tmpm(:,:,1) = getvar(cf_ufil, cn_votemper, jk, 1, npts, kimin=iimin,   kjmin=ijmin+1)
      !      tmpm(:,:,2) = getvar(cf_ufil, cn_votemper, jk, 1, npts, kimin=iimin+1, kjmin=ijmin+1)
      !      zt(:,jk) = 0.5 * ( tmpm(1,:,1) + tmpm(1,:,2) )
      !   END DO

      ELSE                !!Commented Pedro
         ! zonal section at j=ijmin=ijmax
         !tmpz(:,:,1)  = getvar(cn_fhgr, cn_ve1v,   1, npts, 1, kimin=iimin, kjmin=ijmin)
         tmpz(:,:,1)  = getvar(cf_tfil, cn_ve1v,   1, npts, 1, kimin=iimin, kjmin=ijmin)
         eu(:)        = tmpz(:,1,1)
         tmpz(:,:,1)  = getvar(cf_tfil, cn_vlon2d, 1, npts, 1, kimin=iimin, kjmin=ijmin)
         rlonlat(:,1) = tmpz(:,1,1)  ! longitude in this case
         DO jk=1,npk
            ! initiliaze ddepu to gdept()
            ddepu(:,jk) = gdept(jk)

            IF ( lfull ) THEN
               de3(:,jk)   = e3t1d(jk)
               tmpm(:,1,1) = e3w1d(jk)
               tmpm(:,1,2) = e3w1d(jk)
            ELSE
               !!Commented Pedro
               ! vertical metrics (PS case)
               !tmpz(:,:,1)=getvar(cn_fzgr,'e3v_ps',jk, npts, 1, kimin=iimin+1, kjmin=ijmin, ldiom=.TRUE.)
               tmpz(:,:,1)=getvar(cf_tfil,'e3v_ps',jk, npts, 1, kimin=iimin+1, kjmin=ijmin, ldiom=.TRUE.)
               de3(:,jk) = tmpz(:,1,1)
               !tmpz(:,:,1)=getvar(cn_fzgr,'e3w_ps',jk,npts,1, kimin=iimin+1, kjmin=ijmin,   ldiom=.TRUE.)
               !tmpz(:,:,2)=getvar(cn_fzgr,'e3w_ps',jk,npts,1, kimin=iimin+1, kjmin=ijmin+1, ldiom=.TRUE.)
            ENDIF
            !!Commented Pedro
            IF (jk >= 2 ) THEN
               DO ji=1,npts
                  ddepu(ji,jk)= ddepu(ji,jk-1) + tmpz(ji,1,1)!MIN(tmpz(ji,1,1), tmpz(ji,1,2))
               END DO
            ENDIF

            ! Normal velocity
            tmpz(:,:,1)=getvar(cf_tfil,cn_vomecrty,jk,npts,1, kimin=iimin+1, kjmin=ijmin)
            zu(:,jk)=tmpz(:,1,1)
            !!Commented Pedro
            ! salinity and mask
            !tmpz(:,:,2)=getvar(cf_tfil,cn_vosaline,jk, npts, 1, kimin=iimin+1, kjmin=ijmin+1)
            tmpz(:,:,1)=getvar(cf_tfil,cn_vosaline,jk,npts,1, kimin=iimin+1, kjmin=ijmin)
            zmask(:,jk)=tmpz(:,1,1)
            WHERE ( zmask(:,jk) /= 0 ) zmask(:,jk)=1
            ! do not take special care for land value, as the corresponding velocity point is masked
            !zs(:,jk) = 0.5 * ( tmpz(:,1,1) + tmpz(:,1,2) )
            zs(:,jk) = tmpz(:,1,1)
            ! limitation to 'wet' points
            IF ( SUM(zs(:,jk))  == 0 ) THEN
               nk=jk ! first vertical point of the section full on land
               EXIT  ! as soon as all the points are on land
            ENDIF

            !!Commented Pedro
            ! temperature
            tmpz(:,:,1)=getvar(cf_tfil,cn_votemper,jk, npts, 1, kimin=iimin+1, kjmin=ijmin)
            !tmpz(:,:,2)=getvar(cf_tfil,cn_votemper,jk, npts, 1, kimin=iimin+1, kjmin=ijmin+1)
            !zt(:,jk) = 0.5 * ( tmpz(:,1,1) + tmpz(:,1,2) )
            zt(:,jk) = tmpz(:,1,1)
         END DO

      ENDIF
      !! Commented Pedro
      ! compute density only for wet points
      IF ( lneutral ) THEN 
         dsig(:,1:nk)=sigmantr( zt, zs,         npts, nk)*zmask(:,:)
      ELSE
        IF ( refdep == -10. ) THEN
           dsig(:,1:nk)= -zt(:,:)  ! change sign 
        ELSEIF ( refdep == 0. ) THEN
           dsig(:,1:nk)=sigma0( zt, zs,         npts, nk)*zmask(:,:)
        ELSE
           dsig(:,1:nk)=sigmai( zt, zs, refdep, npts, nk)*zmask(:,:)
        ENDIF
      ENDIF

      dsig(:,0)=dsig(:,1)-1.e-4   ! dummy layer for easy interpolation

      ! compute depth of isopynals (nbins+1 )
      DO  jiso =1, nbins+1
         dsigma=dsigma_lev(jiso)
!!!  REM : I and K loop can be inverted if necessary
         DO ji=1,npts
            dhiso(ji,jiso) = gdept(npk)
            DO jk=1,nk 
               IF ( dsig(ji,jk) < dsigma ) THEN
               ELSE
                  ! interpolate between jk-1 and jk
                  dalfa=(dsigma - dsig(ji,jk-1)) / ( dsig(ji,jk) -dsig(ji,jk-1) )
                  IF (ABS(dalfa) > 1.1 .OR. dalfa < 0 ) THEN   ! case dsig(0) = dsig(1)-1.e-4
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
      IF (lbimg ) CALL bimg_writ(jsec)
      IF (lncdf ) CALL cdf_writ(jsec)
      PRINT *,' Total transport in all bins :',TRIM(csection(jsec)),' ',SUM(dtrpbin(jsec,:) )/1.d6

      ! free memory for the next section
      DEALLOCATE ( zu, zt, zs, dsig, ddepu, dhiso, dwtrp, dwtrpbin )
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
   ! previously in cdf_writ(  case lncdf=true )
   IF (lncdf) THEN
      CALL section_init(cf_section, csection,cvarname,clongname, iimina,iimaxa,ijmina,ijmaxa, nsection)
   ENDIF

   DO jsec=1,nsection
      ! setup output variables (section dependant for adaptative variable name (if possible)
      ! define new variables for output 
      IF ( cvarname(jsec) /= 'none' ) THEN
         csuffixvarname='_'//TRIM(cvarname(jsec))
      ELSE
         csuffixvarname=''
      ENDIF
      IF ( clongname(jsec) /= 'none' ) THEN
         cprefixlongname=TRIM(clongname(jsec))//'_'
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
      IF (ltemp) THEN
         cf_outnc = TRIM(csection(jsec))//'_trptemp.nc'
      ELSE
         cf_outnc = TRIM(csection(jsec))//'_trpsig.nc'
      ENDIF

      ncout = create      (cf_outnc, 'none',  ikx,      iky, nbins, cdep=cv_dep               )
      ierr  = createvar   (ncout,    stypvar, nboutput, ipk, id_varout, cdglobal=TRIM(cglobal))
      ierr  = putheadervar(ncout,    cf_tfil, ikx,      iky, nbins, &
           &   pnavlon=rdumlon, pnavlat=rdumlat, pdep=REAL(dsigma_lev), cdep=cv_dep           )

      tim  = getvar1d(cf_tfil, cn_vtimec, 1     )
      ierr = putvar1d(ncout,   tim,       1, 'T')

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
      CHARACTER(LEN=2048), DIMENSION(knumber), INTENT(out  ) :: cdsection
      CHARACTER(LEN=2048), DIMENSION(knumber), INTENT(out  ) :: cdvarname
      CHARACTER(LEN=2048), DIMENSION(knumber), INTENT(out  ) :: cdlongname
      INTEGER(KIND=4),                        INTENT(inout) :: knumber
      INTEGER(KIND=4), DIMENSION(knumber),    INTENT(out  ) :: kimin, kimax, kjmin, kjmax

      ! Local variables
      INTEGER(KIND=4)                                       :: jsec
      INTEGER(KIND=4)                                       :: ii, inum=10
      INTEGER(KIND=4)                                       :: ipos  
      CHARACTER(LEN=2048)                                    :: cline
      CHARACTER(LEN=2048), DIMENSION(3)                       :: cldum
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
         ipos = index(cline,' ')
         DO WHILE ( ipos > 1 ) 
            ii = ii + 1
            cldum(ii) = cline(1:ipos - 1 )
            cline = TRIM ( cline(ipos+1:) )
            ipos  = index( cline,' ' ) 
            IF ( ii >= 3 ) EXIT
         END DO
         cdsection(jsec) = TRIM(cldum(1) )
         cdvarname(jsec) = TRIM(cldum(2) )
         cdlongname(jsec) = TRIM(cldum(3) )
         READ(inum,*    ) kimin(jsec), kimax(jsec), kjmin(jsec), kjmax(jsec)
      END DO

      CLOSE(inum)

   END SUBROUTINE section_init

   SUBROUTINE bimg_writ( ksec)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE bimg_writ  ***
      !!
      !! ** Purpose :  Write output bimg files if required 
      !!
      !! ** Method  :  Most of the variables are global 
      !!
      !!----------------------------------------------------------------------
      INTEGER(KIND=4), INTENT(in) :: ksec  ! number of the section

      INTEGER(KIND=4)             :: ji, jk
      !!----------------------------------------------------------------------
      ! (along section, depth ) 2D variables
      cf_bimg=TRIM(csection(ksec))//'_trpdep.bimg'
      OPEN(numbimg,FILE=cf_bimg,FORM='UNFORMATTED')
      cldum=' 4 dimensions in this isopycnal file '
      WRITE(numbimg) cldum

      cldum=' 1: T ;  2: S ; 3: sigma ; 4: Velocity '
      WRITE(numbimg) cldum

      WRITE(cldum,'(a,4i5.4)') ' from '//TRIM(csection(ksec)), iimin,iimax,ijmin,ijmax
      WRITE(numbimg) cldum

      cldum=' file '//TRIM(cf_tfil)
      WRITE(numbimg) cldum

      WRITE(numbimg) npts,nk,1,1,4,0
      WRITE(numbimg) 1.,-float(nk),1.,1., 0.
      WRITE(numbimg) 0.
      WRITE(numbimg) 0.

      WRITE(numbimg) (( REAL(zt(ji,jk)  ), ji=1,npts), jk=nk,1,-1 ) ! temperature 
      WRITE(numbimg) (( REAL(zs(ji,jk)  ), ji=1,npts), jk=nk,1,-1 ) ! salinity
      WRITE(numbimg) (( REAL(dsig(ji,jk)), ji=1,npts), jk=nk,1,-1 ) ! density
      WRITE(numbimg) (( REAL(zu(ji,jk)  ), ji=1,npts), jk=nk,1,-1 ) ! normal velocity
      CLOSE(numbimg)

      ! (along section, sigma ) 2D variables
      cf_bimg=TRIM(csection(ksec))//'_trpsig.bimg'
      OPEN(numbimg,FILE=cf_bimg,FORM='UNFORMATTED')
      cldum=' 3 dimensions in this isopycnal file '
      WRITE(numbimg) cldum
      cldum=' 1: hiso ;  2: bin trp ; 3: cumulated  trp '
      WRITE(numbimg) cldum
      WRITE(cldum,'(a,4i5.4)') ' from '//TRIM(csection(ksec)), iimin,iimax,ijmin,ijmax
      WRITE(numbimg) cldum
      cldum=' file '//TRIM(cf_tfil)
      WRITE(numbimg) cldum
      WRITE(numbimg) npts,nbins,1,1,3,0
      WRITE(numbimg) 1.,-REAL(dsigma_lev(nbins)),1.,REAL(dltsig), 0.
      WRITE(numbimg) 0.
      WRITE(numbimg) 0.
      WRITE(numbimg) (( REAL(dhiso(ji,jiso)   ),      ji=1,npts), jiso=nbins,1,-1) ! isopyc depth
      WRITE(numbimg) (( REAL(dwtrpbin(ji,jiso))/1.e6, ji=1,npts), jiso=nbins,1,-1) ! binned transport
      WRITE(numbimg) (( REAL(dwtrp(ji,jiso)   )/1.e6, ji=1,npts), jiso=nbins,1,-1) ! cumulated transport
      CLOSE(numbimg)

   END SUBROUTINE bimg_writ

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
      CHARACTER(LEN=2048)            :: csuffixvarnam
      CHARACTER(LEN=2048)            :: cprefixlongnam
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

      !    tim  = getvar1d(cf_tfil, cn_vtimec, 1     )
      !    ierr = putvar1d(icout,   tim,       1, 'T')

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

END PROGRAM cdfsigtrp_broken
