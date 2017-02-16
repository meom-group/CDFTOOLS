PROGRAM cdf_xtract_brokenline
   !!======================================================================
   !!                     ***  PROGRAM  cdf_xtract_brokenline  ***
   !!=====================================================================
   !!  ** Purpose : Extract temperature, Salinity and velocity components
   !!               along a broken line formed by various legs
   !!
   !!  ** Method  : A broken line is defined by various segments or leg.
   !!               Each leg is defined by its starting and endig point
   !!               given as geographical coordinates on standard input.
   !!
   !!
   !! History : 2.1  : 12/2009  : R. Dussin    : Original code
   !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
   !!           3.0  : 05/2013  : T. Penduff & R. Dussin  : Saving new variables
   !!           3.0  : 05/2013  : J.M. Molines : Code review, generalization 
   !!           3.0  : 06/2013  : J.M. Molines : allows multiple section files (optimization).
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   routines      : description
   !!----------------------------------------------------------------------
   USE cdfio
   USE cdftools
   USE modcdfnames
   !!----------------------------------------------------------------------
   !! CDFTOOLS_3.0 , MEOM 2011
   !! $Id$
   !! Copyright (c) 2011, J.-M. Molines
   !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
   !!----------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER(KIND=4) :: jsec, jleg, jt, jk,  jipt, jvar ! dummy loop index
   INTEGER(KIND=4) :: narg, iargc, ijarg, ifree       ! command line
   INTEGER(KIND=4) :: numin=10                        ! logical unit for input section file
   INTEGER(KIND=4) :: numout=11                       ! logical unit for output section.dat (used in cdftransport)
   INTEGER(KIND=4) :: npiglo, npjglo, npk, npt        ! size of the domain
   INTEGER(KIND=4) :: iimin, iimax, ijmin, ijmax      ! ending points of a leg in model I J
   INTEGER(KIND=4) :: ii, ij, ii1, ij1, ipoint        ! working integer
   INTEGER(KIND=4) :: ierr                            ! Netcdf error and ncid
   INTEGER(KIND=4) :: nvar = 16                       ! number of output variables (modified after if options)
   INTEGER(KIND=4) :: np_tem, np_sal, np_una, np_vna  ! index for output variable
   INTEGER(KIND=4) :: np_isec, np_jsec, np_e2vn       !  "
   INTEGER(KIND=4) :: np_e1vn, np_e3un, np_e3vn       !  "
   INTEGER(KIND=4) :: np_vmod, np_e1v,  np_e3v        !  "
   INTEGER(KIND=4) :: np_vmsk, np_baro, np_bat        !  "
   INTEGER(KIND=4) :: np_ssh,  np_mld,  np_vt, np_vs  !  "
   INTEGER(KIND=4) :: np_icethick, np_icefra          !  "
   INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ncout                     ! Netcdf error and ncid
   INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout  ! netcdf output stuff

   ! broken line definition
   INTEGER(KIND=4) :: nsec = 1                       ! number of sections
   INTEGER(KIND=4) :: nn                             ! working integer (number of points in a leg)
   INTEGER(KIND=4) :: nstamax                        ! maximum number of points per defined broken line
   INTEGER(KIND=4) :: npsecmax                       ! maximum number of points per defined model broken line
   INTEGER(KIND=4) :: ista                           ! working integer
   INTEGER(KIND=4), DIMENSION(:),     ALLOCATABLE :: npsec            ! number of points defining the model broken line
   INTEGER(KIND=4), DIMENSION(:),     ALLOCATABLE :: nsta             ! number of points defining the broken line
   INTEGER(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: normu_sec, normv_sec ! velocity normalization per section
   INTEGER(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: iisec, ijsec     ! F-index of points on the broken line
   INTEGER(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: norm_u, norm_v   ! velocity normalization per leg
   INTEGER(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: iista, ijsta     ! I,J position of the point on the broken line
   INTEGER(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: ikeepn           ! Number of points per leg
   INTEGER(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: iilegs, ijlegs   ! F-index of points on the broken line per leg

   REAL(KIND=4)                              :: xmin, xmax, ymin, ymax !
   REAL(KIND=4)                              :: ztmp

   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                 ! Model time array
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rxx, ryy            ! leg i j index of F points
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlonsta, rlatsta    ! Geographic position defining legs

   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1v, e3v            ! V point relevant metric
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2u, e3u            ! U point relevant metric
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: hdepw               ! model bathymetry
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rmbat               ! model bathymetry (levels)

   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlonu, rlatu        ! model long and lat of U points
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlonv, rlatv        ! model long and lat of U points
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlonf, rlatf        ! model long and lat of F points
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: temper, saline      ! model Temperature and salinity
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: uzonal, vmerid      ! model zonal and meridional velocity
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ssh, rmld           ! model SSH and MLD
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ricethick, ricefra  ! ice thickness and fraction
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zvmod               ! ice thickness and fraction
   ! along section array (dimension x,z or x,1 )
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: tempersec, salinesec, uzonalsec, vmeridsec
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: sshsec, rmldsec
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: ricethicksec, ricefrasec
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rlonsec, rlatsec, risec, rjsec
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e3usec, e3vsec
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: batsec
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: vmasksec
   REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: e1vsec, e2usec  ! 3rd dimension for sections

   REAL(KIND=8)                              :: dtmp    ! temporary cumulating variable
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE   :: dbarot  ! for barotropic transport computation

   CHARACTER(LEN=2048) :: cf_tfil , cf_ufil, cf_vfil   ! input T U V files
   CHARACTER(LEN=2048) :: cf_icefil                    ! input ice file
   CHARACTER(LEN=2048) :: cf_root=''                   ! root name used as prefix
   CHARACTER(LEN=2048) :: cf_out                       ! output file
   CHARACTER(LEN=2048) :: cf_secdat                    ! output section file (suitable for cdftransport or cdfsigtrp)
   CHARACTER(LEN=2048) :: cverb='n'                    ! verbose key for findij
   CHARACTER(LEN=5  ) :: cstar, cend                  ! dummy character variable
   CHARACTER(LEN=2048) :: cldum                        ! can handle a long list of section files ...
   CHARACTER(LEN=2048), DIMENSION(:), ALLOCATABLE :: cf_sec    ! input section file dim: nsec
   CHARACTER(LEN=2048), DIMENSION(:), ALLOCATABLE :: csection  ! section name

   LOGICAL  :: lchk                                  ! flag for missing files
   LOGICAL  :: lverbose = .FALSE.                    ! flag for verbosity
   LOGICAL  :: lsecfile = .FALSE.                    ! flag for input section file
   LOGICAL  :: lssh     = .FALSE.                    ! flag for saving ssh
   LOGICAL  :: lmld     = .FALSE.                    ! flag for saving mld
   LOGICAL  :: lice     = .FALSE.                    ! flag for saving ice*
   LOGICAL  :: lvt      = .FALSE.                    ! flag for saving products vt, vs
   LOGICAL  :: ll_ssh, ll_mld, ll_ice                ! working flag for jk =1

   TYPE (variable), DIMENSION(:), ALLOCATABLE :: stypvar  ! variable definition and attributes
   !!----------------------------------------------------------------------
   ! 1. : Initialization
   ! --------------------
   CALL ReadCdfNames()

   ! check argument number and show usage if necessary
   narg = iargc()
   IF ( narg < 3 ) THEN
      PRINT *,' usage :  cdf_xtrac_brokenline T-file U-file V-file [ice-file] ....'
      PRINT *,'    [-f section_filei,sec_file2, ... ] [-verbose] [-ssh ] [-mld] [-ice] '
      PRINT *,'    [-vt] [-o ROOT_name]'
      PRINT *,'      '
      PRINT *,'     PURPOSE :'
      PRINT *,'        This tool extracts model variables from model files for a geographical' 
      PRINT *,'      broken line, similar to an oceanographic campaign where an oceanic '
      PRINT *,'      section is formed by one or more legs.' 
      PRINT *,'        The broken line is specified by the position of ending points of each'
      PRINT *,'      leg, given in an ASCII file. OVIDE section is taken as default, when no'
      PRINT *,'      section file is provided.'
      PRINT *,'        This tool provides a netcdf file similar to a model file, but with a '
      PRINT *,'      degenerated y dimension (1). In order to be able to use standard CDFTOOLS'
      PRINT *,'      relevant metric variables are saved into the output file, such as pseudo'
      PRINT *,'      e1v and e3v_ps and vmask. Therefore the output file can be considered as'
      PRINT *,'      a mesh_hgr, mesh_zgr and mask file for any ''meridional'' computation.'
      PRINT *,'        This tools works with temperatures, salinities and normal velocities.'
      PRINT *,'      The broken line is approximated in the model, by a succession of segments'
      PRINT *,'      joining F-points. The velocity is taken as either U or V depending on the'
      PRINT *,'      orientation of the segment, temperatures and salinities are interpolated'
      PRINT *,'      on the velocity points. When progressing along the broken line, velocity'
      PRINT *,'      is positive when heading to the right of the progression.' 
      PRINT *,'        The barotropic transport across the broken line is computed, using the'
      PRINT *,'      same sign convention. On a closed broken line, the barotropic transport'
      PRINT *,'      should be very small.'
      PRINT *,'      ' 
      PRINT *,'     ARGUMENTS :'
      PRINT *,'      T-file   :  model gridT file '
      PRINT *,'      U-file   :  model gridU file '
      PRINT *,'      V-file   :  model gridV file '
      PRINT *,'      ice-file :  model ice file '
      PRINT *,'      ' 
      PRINT *,'     OPTIONS :'
      PRINT *,'      -f section_file1,section_file2,... : provide a comma separated list of'
      PRINT *,'              files for section definition. Section_file is an ascii file as '
      PRINT *,'              follows:'
      PRINT *,'             * line #1 : name of the section (e.g. ovide). '
      PRINT *,'                  Will be used for naming the output file.'
      PRINT *,'             * line #2 : number of points defining the broken line.'
      PRINT *,'             * line #3-end : a pair of Longitude latitude values defining'
      PRINT *,'                   the points. If not supplied, use hard-coded information'
      PRINT *,'                   for OVIDE section. A comment can be added at the end of'
      PRINT *,'                   of the lines, using a # as separator'
      PRINT *,'      -verbose : increase verbosity  ' 
      PRINT *,'      -ssh     : also save ssh along the broken line.'
      PRINT *,'      -mld     : also save mld along the broken line.'
      PRINT *,'      -ice     : also save ice properties along the broken line.'
      PRINT *,'      -vt      : also save products vt and vs along the broken line.'
      PRINT *,'      -o ROOT-name : specified the prefix to be used for the output file name.'
      PRINT *,'                 Note that it may be a good idea to include a separator '
      PRINT *,'                 character such as _ at the end of the ROOT_name.'
      PRINT *,'     '
      PRINT *,'     REQUIRED FILES :'
      PRINT *,'      ', TRIM(cn_fhgr),' and ',TRIM(cn_fzgr),' in the current directory ' 
      PRINT *,'      '
      PRINT *,'     OUTPUT : '
      PRINT *,'       netcdf file : <section_name>.nc (default). If -o option is used, the'
      PRINT *,'                     name will be <ROOT-name><section_name>.nc'
      PRINT *,'         variables : temperature, salinity, normal velocity, pseudo V metrics,'
      PRINT *,'                    mask, barotropic transport, bathymetry of velocity points.'
      PRINT *,'                    Additional variables can be set when using options.'
      PRINT *,'       ASCII file : <section_name>_section.dat usefull for cdftransport, gives'
      PRINT *,'                  the position in I,J of the geographical input points.'
      PRINT *,'      '
      PRINT *,'     SEE ALSO :'
      PRINT *,'        cdftransport, cdfmoc, cdfmocsig. This tool replaces cdfovide.' 
      PRINT *,'      '
      STOP
   ENDIF

   ! Decode command line
   ijarg = 1 ; ifree=0
   DO WHILE ( ijarg <= narg ) 
      CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1
      SELECT CASE  ( cldum   )
      CASE ( '-verbose' ) ; lverbose=.TRUE.  ; cverb='y'
      CASE ( '-ssh'     ) ; lssh    =.TRUE.  ; nvar = nvar + 1  ! 
      CASE ( '-mld'     ) ; lmld    =.TRUE.  ; nvar = nvar + 1  !
      CASE ( '-ice'     ) ; lice    =.TRUE.  ; nvar = nvar + 2  !
      CASE ( '-vt '     ) ; lvt     =.TRUE.  ; nvar = nvar + 2  !
      CASE ( '-o '      ) ;  CALL getarg(ijarg, cf_root) ; ijarg = ijarg + 1  !
      CASE ( '-f' )       ;  CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; lsecfile=.TRUE.
         CALL ParseFiles(cldum)        ! many section files can be given separated with comma
      CASE DEFAULT 
         ifree = ifree + 1
         SELECT CASE ( ifree )
         CASE ( 1 ) ;             cf_tfil   = cldum
         CASE ( 2 ) ;             cf_ufil   = cldum
         CASE ( 3 ) ;             cf_vfil   = cldum
         CASE ( 4 ) ; IF ( lice ) cf_icefil = cldum
         END SELECT
      END SELECT
   ENDDO

   ! check file existence
   lchk = chkfile(cn_fhgr )
   lchk = chkfile(cn_fzgr ) .OR. lchk
   lchk = chkfile(cf_tfil ) .OR. lchk
   lchk = chkfile(cf_ufil ) .OR. lchk
   lchk = chkfile(cf_vfil ) .OR. lchk
   IF ( lsecfile ) THEN 
      DO jsec = 1, nsec
         lchk = chkfile(cf_sec(jsec) ) .OR. lchk
      ENDDO
   ENDIF
   IF ( lchk     ) STOP ! missing files

   ! nvar and nsec are  now fixed
   ALLOCATE( stypvar(nvar), ipk(nvar), id_varout(nvar) )
   ALLOCATE( csection (nsec),nsta(nsec), ncout(nsec), dbarot(nsec)  )

   ! read section file if required
   IF ( lsecfile ) THEN
      nstamax = 0
      DO jsec =1, nsec
         OPEN(numin, file=cf_sec(jsec) )
         READ(numin,'(a)') csection(jsec)
         READ(numin,*    ) nsta(jsec)
         nstamax  = MAX( nstamax, nsta(jsec))
         CLOSE (numin)
      END DO
      IF ( lverbose ) PRINT *, 'NSTAMAX = ', nstamax
      ALLOCATE ( iista(nstamax,nsec), ijsta(nstamax,nsec), ikeepn(nstamax -1,nsec )  )
      ALLOCATE ( npsec(nsec) )
      ALLOCATE ( rlonsta(nstamax,nsec), rlatsta(nstamax,nsec) )

      DO jsec =1, nsec
         OPEN(numin, file=cf_sec(jsec) )
         READ(numin,'(a)') csection(jsec)
         READ(numin,*    ) ista
         DO jipt = 1, ista
            READ(numin, * ) rlonsta(jipt,jsec), rlatsta(jipt,jsec)
         ENDDO
         CLOSE (numin)
      ENDDO
   ELSE    ! default to OVIDE section
      nsec = 1
      nstamax = 1
      nsta(1) = 5
      csection(1) = 'ovide'
      ALLOCATE ( iista(nstamax,nsec), ijsta(nstamax,nsec), ikeepn(nstamax -1,nsec )  )
      ALLOCATE ( rlonsta(nstamax,nsec), rlatsta(nstamax,nsec) )

      ! D. Desbruyeres : Location of leg points that define the 4 legs of the OVIDE section
      rlonsta(1,1) = -43.70 ; rlatsta(1,1) = 59.90    ! 
      rlonsta(2,1) = -30.30 ; rlatsta(2,1) = 58.90    ! 
      rlonsta(3,1) = -19.40 ; rlatsta(3,1) = 44.90    ! 
      rlonsta(4,1) = -12.65 ; rlatsta(4,1) = 40.33    ! 
      rlonsta(5,1) = -08.70 ; rlatsta(5,1) = 40.33    ! 
   ENDIF

   ! 2. Find the model F-points along the legs of the section
   ! --------------------------------------------------------
   npiglo = getdim (cf_tfil, cn_x)
   npjglo = getdim (cf_tfil, cn_y)
   npk    = getdim (cf_tfil, cn_z)
   npt    = getdim (cf_tfil, cn_t)

   IF ( lverbose ) THEN
      PRINT *, 'NPIGLO = ', npiglo
      PRINT *, 'NPJGLO = ', npjglo
      PRINT *, 'NPK    = ', npk
      PRINT *, 'NPT    = ', npt
   ENDIF

   ! input  2D fields
   ALLOCATE(rlonu(npiglo,npjglo), rlatu(npiglo,npjglo))
   ALLOCATE(rlonv(npiglo,npjglo), rlatv(npiglo,npjglo))
   ALLOCATE(rlonf(npiglo,npjglo), rlatf(npiglo,npjglo))
   ALLOCATE(temper(npiglo,npjglo), saline(npiglo,npjglo))
   ALLOCATE(uzonal(npiglo,npjglo), vmerid(npiglo,npjglo))
   ALLOCATE(e1v(npiglo,npjglo))
   ALLOCATE(e2u(npiglo,npjglo))
   ALLOCATE(e3u(npiglo,npjglo), e3v(npiglo,npjglo))
   ALLOCATE(hdepw(npiglo,npjglo), rmbat(npiglo, npjglo))
   IF ( lssh ) ALLOCATE (ssh (npiglo, npjglo))
   IF ( lmld ) ALLOCATE (rmld(npiglo, npjglo))
   IF ( lice ) ALLOCATE(ricethick(npiglo,npjglo),ricefra(npiglo,npjglo))

   ! allocate section working arrays
   ALLOCATE ( iilegs(nstamax-1, npiglo+npjglo, nsec), ijlegs(nstamax-1, npiglo+npjglo, nsec) )
   ALLOCATE ( norm_u(nstamax-1, nsec) , norm_v(nstamax-1, nsec) )
   ALLOCATE ( rxx(npiglo+npjglo, nsec), ryy(npiglo+npjglo, nsec) )
   ALLOCATE ( tim (npt) )

   iilegs = 0  ; ijlegs = 0 ; npsec(:) = 0

   ! Loop on the sections 
   npsecmax = 0
   DO jsec  = 1, nsec
      ! loop on the legs
      DO jleg = 1, nsta(jsec)-1

         xmin = rlonsta(jleg,  jsec)
         ymin = rlatsta(jleg,  jsec)
         xmax = rlonsta(jleg+1,jsec)
         ymax = rlatsta(jleg+1,jsec)

         ! return ending points of a leg in I J model coordinates
         PRINT *,TRIM(csection(jsec)),' leg ',jleg
         CALL cdf_findij ( xmin, xmax, ymin, ymax, iimin, iimax, ijmin, ijmax, &
              &            cd_coord=cn_fhgr, cd_point='F', cd_verbose=cverb)

         ! save leg information
         iista(jleg  , jsec) = iimin
         ijsta(jleg  , jsec) = ijmin
         iista(jleg+1, jsec) = iimax
         ijsta(jleg+1, jsec) = ijmax

         ! find the broken line between P1 (iimin,ijmin) and P2 (iimax, ijmax)
         CALL broken_line( iimin, iimax, ijmin, ijmax, rxx(:,jsec), ryy(:,jsec), nn, npiglo, npjglo,  &
              &            norm_u(jleg,jsec), norm_v(jleg,jsec) )
         ikeepn(jleg,jsec) = nn  ! number of points (F) on leg jleg
         npsec(jsec)       = npsec(jsec) + nn   ! total number of points (F) on the broken line

         IF ( lverbose) PRINT *, 'Leg ', jleg,' : npoints : ', nn

         IF ( jleg == 1 ) THEN
            ! we want to ensure that the broken line start in the direction that we specify
            IF ( INT(rxx(1,jsec)) /= iimin .OR.  INT(ryy(1,jsec)) /= ijmin ) THEN 
              IF ( lverbose ) PRINT *,' First leg is to be reverse'
                iilegs(jleg,1:nn,jsec)=rxx(nn:1:-1,jsec)
                ijlegs(jleg,1:nn,jsec)=ryy(nn:1:-1,jsec)
              ELSE
                iilegs(jleg,1:nn,jsec)=rxx(1:nn,jsec)
                ijlegs(jleg,1:nn,jsec)=ryy(1:nn,jsec)
              ENDIF
         ELSE  ! check the continuity between legs
            IF ( iilegs(jleg-1, ikeepn(jleg-1,jsec),jsec) == rxx(1,jsec) .AND. &
             &   ijlegs(jleg-1, ikeepn(jleg-1,jsec),jsec) == ryy(1,jsec) ) THEN  ! continuity
               IF ( lverbose ) PRINT *,' Leg ',jleg ,' is continuous ...'
               iilegs(jleg,1:nn,jsec)=rxx(1:nn,jsec)
               ijlegs(jleg,1:nn,jsec)=ryy(1:nn,jsec)
            ELSE                          ! reverse sense
               IF ( lverbose ) PRINT *,' Leg ',jleg ,' requires inversion ...'
               iilegs(jleg,1:nn,jsec)=rxx(nn:1:-1,jsec)
               ijlegs(jleg,1:nn,jsec)=ryy(nn:1:-1,jsec)
            END IF
         ENDIF

         IF ( lverbose) THEN
            PRINT *, '          Leg      iileg        ijleg rxx          ryy '
            DO jipt = 1, ikeepn(jleg,jsec)  ! ( nn ! ) 
               PRINT *, jleg, iilegs(jleg,jipt,jsec), ijlegs(jleg,jipt,jsec) ,rxx(jipt,jsec), ryy(jipt,jsec)
            END DO
         ENDIF
      END DO !! loop on the legs
      npsecmax = MAX(npsecmax, npsec(jsec))  ! maximum number of point in any section
   END DO !! loop on the sections
   IF ( lverbose)  PRINT *,' NPSECMAX = ', npsecmax
   
   ! Now can allocate the section arrays 
   ALLOCATE( rlonsec(npsecmax,  1), rlatsec(npsecmax,1) )
   ALLOCATE( risec  (npsecmax,  1), rjsec  (npsecmax,1) )
   ALLOCATE( batsec   (npsecmax-1,1  ), vmasksec (npsecmax-1,npk) )
   ALLOCATE( tempersec(npsecmax-1,npk), salinesec(npsecmax-1,npk) )
   ALLOCATE( uzonalsec(npsecmax-1,npk), vmeridsec(npsecmax-1,npk) )
   ALLOCATE( zvmod (npsecmax-1,1) )  ! working array
   IF ( lssh ) ALLOCATE ( sshsec (npsecmax-1,1) )
   IF ( lmld ) ALLOCATE ( rmldsec(npsecmax-1,1) )
   IF ( lice ) ALLOCATE(ricethicksec(npsecmax-1,1),ricefrasec(npsecmax-1,1))

   ! Next arrays are initialized outside the vertical loop and thus require a section index
   ALLOCATE ( iisec    (npsecmax,nsec),     ijsec(npsecmax,nsec) ) 
   ALLOCATE ( normu_sec(npsecmax,nsec), normv_sec(npsecmax,nsec) ) 
   ALLOCATE ( e2usec   (npsecmax-1,1,nsec), e3usec(npsecmax-1,npk) )
   ALLOCATE ( e1vsec   (npsecmax-1,1,nsec), e3vsec(npsecmax-1,npk) )

   ! 3. : Extraction along the legs
   ! ------------------------------

   e1vsec = -9999.
   e2usec = -9999.

   risec(:,:) = 0.
   rjsec(:,:) = 0.

   rlonu(:,:) = getvar(cn_fhgr, cn_glamu, 1, npiglo, npjglo)
   rlatu(:,:) = getvar(cn_fhgr, cn_gphiu, 1, npiglo, npjglo)
   rlonv(:,:) = getvar(cn_fhgr, cn_glamv, 1, npiglo, npjglo)
   rlatv(:,:) = getvar(cn_fhgr, cn_gphiv, 1, npiglo, npjglo)
   rlonf(:,:) = getvar(cn_fhgr, cn_glamf, 1, npiglo, npjglo)
   rlatf(:,:) = getvar(cn_fhgr, cn_gphif, 1, npiglo, npjglo)
   e1v(:,:)   = getvar(cn_fhgr, cn_ve1v,  1, npiglo, npjglo)
   e2u(:,:)   = getvar(cn_fhgr, cn_ve2u,  1, npiglo, npjglo)
   rmbat(:,:) = getvar(cn_fzgr, cn_mbathy,1, npiglo, npjglo)
   hdepw(:,:) = getvar(cn_fzgr, cn_hdepw, 1, npiglo, npjglo)

   !  Loop on section for metrics and non z-depending variables
   DO jsec = 1, nsec   ! loop on sections
      cf_out    = TRIM(cf_root)//TRIM(csection(jsec))//'.nc'
      cf_secdat = TRIM(cf_root)//TRIM(csection(jsec))//'_section.dat'

      ipoint = 1
      DO jleg=1, nsta(jsec) -1      ! loop on legs 
         ipoint = ipoint -1  ! trick to avoid repetition of points in between legs
         DO jipt=1, ikeepn(jleg,jsec) 
            ipoint = ipoint + 1
            iisec    (ipoint,jsec) = iilegs(jleg,jipt,jsec)  ! i-index
            ijsec    (ipoint,jsec) = ijlegs(jleg,jipt,jsec)  ! j-index
            normu_sec(ipoint,jsec) = norm_u(jleg,jsec)
            normv_sec(ipoint,jsec) = norm_v(jleg,jsec)
         END DO
      END DO
      ! adjust npsec to its real value ( 2nd part of the trick)
      npsec(jsec) = ipoint

      ! now that we know the model grid and bathy do fancy print of the legs.
      CALL FancyPrint(jsec)

      ! loop on 2d arrays
      DO jipt = 1,npsec(jsec)
         ii = iisec(jipt,jsec)
         ij = ijsec(jipt,jsec)

         risec  (jipt,1) = ii
         rjsec  (jipt,1) = ij
      END DO

      DO jipt=1,npsec(jsec)-1
         ii  = iisec(jipt  ,jsec) ; ij  = ijsec(jipt  ,jsec)
         ii1 = iisec(jipt+1,jsec) ; ij1 = ijsec(jipt+1,jsec)
         IF ( ij1 == ij ) THEN ! horizontal segment
            e2usec(jipt,1,jsec) = 0.
            IF ( ii1 > ii ) THEN ! eastward
               e1vsec (jipt,1,jsec) = e1v  (ii+1,ij)
               rlonsec(jipt,1) = rlonv(ii+1,ij)
               rlatsec(jipt,1) = rlatv(ii+1,ij)
               batsec (jipt,1) = MIN( hdepw(ii+1,ij),  hdepw(ii+1,ij+1) )
            ELSE
               e1vsec (jipt,1,jsec) = e1v  (ii,ij)
               rlonsec(jipt,1) = rlonv(ii,ij)
               rlatsec(jipt,1) = rlatv(ii,ij)
               batsec (jipt,1) = MIN( hdepw(ii,ij),  hdepw(ii,ij+1) )
            ENDIF

         ELSEIF ( ii1 == ii ) THEN ! vertical segment
            e1vsec(jipt,1,jsec) = 0.
            IF ( ij1 < ij ) THEN ! southward
               e2usec (jipt,1,jsec) = e2u  (ii,ij)
               rlonsec(jipt,1) = rlonu(ii,ij)
               rlatsec(jipt,1) = rlatu(ii,ij)
               batsec (jipt,1) = MIN( hdepw(ii,ij),  hdepw(ii+1,ij) )
            ELSE
               e2usec (jipt,1,jsec) = e2u  (ii,ij+1)
               rlonsec(jipt,1) = rlonu(ii,ij+1)
               rlatsec(jipt,1) = rlatu(ii,ij+1)
               batsec (jipt,1) = MIN( hdepw(ii,ij+1),  hdepw(ii+1,ij+1) )
            ENDIF

         ELSE
            PRINT *, 'problem 1 for JIPT = ', jipt
            PRINT *, '             I(P2)=',ii1, 'J(P1)=', ii
            PRINT *, '             J(P2)=',ij1, 'J(P1)=', ij
            EXIT 
         ENDIF
      END DO

      ! Prepare output file ( here because rlonsec and rlatsec required )
      CALL CreateOutputFile (jsec )

      ierr = putvar (ncout(jsec), id_varout(np_isec), risec(:,1),                            1,  npsec(jsec)  , 1 )
      ierr = putvar (ncout(jsec), id_varout(np_jsec), rjsec(:,1),                            1,  npsec(jsec)  , 1 )
      ierr = putvar (ncout(jsec), id_varout(np_e2vn), e2usec(:,1,jsec),                      1,  npsec(jsec)-1, 1 )
      ierr = putvar (ncout(jsec), id_varout(np_e1vn), e1vsec(:,1,jsec),                      1,  npsec(jsec)-1, 1 )
      ierr = putvar (ncout(jsec), id_varout(np_e1v ), e2usec(:,1,jsec) + e1vsec(:,1,jsec),   1,  npsec(jsec)-1, 1 )
      ierr = putvar (ncout(jsec), id_varout(np_bat ), batsec(:,1),                           1,  npsec(jsec)-1, 1 )

   END DO  ! section for non depth dependent

   ! Temperature and salinity are interpolated on the respective U or V  point for better flux computation
   DO jt=1, npt  ! time loop
      dbarot(:) = 0.d0    ! reset barotropic transport  for all sections
      IF ( lssh ) ssh (:,:)     = getvar(cf_tfil  , cn_sossheig, 1, npiglo, npjglo, ktime = jt)
      IF ( lmld ) rmld(:,:)     = getvar(cf_tfil  , cn_somxl010, 1, npiglo, npjglo, ktime = jt)
      IF ( lice ) ricethick(:,:) = getvar(cf_icefil, cn_iicethic , 1, npiglo, npjglo, ktime = jt)
      IF ( lice ) ricefra(:,:)   = getvar(cf_icefil, cn_ileadfra , 1, npiglo, npjglo, ktime = jt)

      DO jk=1,npk   ! level loop , read only once the horizontal slab
         temper(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime = jt)
         saline(:,:) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime = jt)
         uzonal(:,:) = getvar(cf_ufil, cn_vozocrtx, jk, npiglo, npjglo, ktime = jt)
         vmerid(:,:) = getvar(cf_vfil, cn_vomecrty, jk, npiglo, npjglo, ktime = jt)
         e3u(:,:)    = getvar(cn_fzgr, 'e3u_ps',    jk, npiglo, npjglo, ldiom=.TRUE.)
         e3v(:,:)    = getvar(cn_fzgr, 'e3v_ps',    jk, npiglo, npjglo, ldiom=.TRUE.)

         ll_ssh = ( lssh .AND. jk == 1 )
         ll_mld = ( lmld .AND. jk == 1 )
         ll_ice = ( lice .AND. jk == 1 )
         tempersec(:,:) = 0.
         salinesec(:,:) = 0.
         uzonalsec(:,:) = 0.
         vmeridsec(:,:) = 0.

         DO jsec = 1, nsec  ! section loop at level jk
            DO jipt=1,npsec(jsec)-1
               ii  = iisec(jipt  ,jsec) ; ij  = ijsec(jipt  ,jsec)  ! F point  position
               ii1 = iisec(jipt+1,jsec) ; ij1 = ijsec(jipt+1,jsec)  ! Next F point position
               IF ( ij1  == ij ) THEN ! horizontal segment
                  uzonalsec(jipt,jk) = 0.
                  e3usec   (jipt,jk) = 0.
                  IF ( ii1 > ii ) THEN ! eastward

                     IF ( MIN( saline(ii+1,ij) , saline(ii+1,ij+1))  == 0. ) THEN
                        tempersec(jipt,jk) = 0. ; salinesec(jipt,jk) = 0.
                        IF ( ll_ssh ) sshsec(jipt,jk) = 0.
                        IF ( ll_mld ) rmldsec(jipt,jk) = 0.
                        IF ( ll_ice ) ricethicksec(jipt,jk) = 0.
                        IF ( ll_ice ) ricefrasec(jipt,jk) = 0.
                     ELSE
                        tempersec(jipt,jk) = 0.5 * ( temper(ii+1,ij) + temper(ii+1,ij+1) )
                        salinesec(jipt,jk) = 0.5 * ( saline(ii+1,ij) + saline(ii+1,ij+1) )
                        IF ( ll_ssh ) sshsec (jipt,jk) = 0.5 * ( ssh (ii+1,ij) + ssh (ii+1,ij+1) )
                        IF ( ll_mld ) rmldsec(jipt,jk) = 0.5 * ( rmld(ii+1,ij) + rmld(ii+1,ij+1) )
                        IF ( ll_ice ) ricethicksec(jipt,jk) = 0.5 * ( ricethick(ii+1,ij) + ricethick(ii+1,ij+1) )
                        IF ( ll_ice ) ricefrasec(jipt,jk)   = 0.5 * ( ricefra(ii+1,ij) + ricefra(ii+1,ij+1) )
                     ENDIF
                     vmeridsec(jipt,jk) = vmerid(ii+1,ij) * normv_sec(jipt,jsec)
                     e3vsec   (jipt,jk) = e3v   (ii+1,ij)

                  ELSE ! westward

                     IF ( MIN( saline(ii,ij) , saline(ii,ij+1) ) == 0. ) THEN
                        tempersec(jipt,jk) = 0. ; salinesec(jipt,jk) = 0.
                        IF ( ll_ssh ) sshsec (jipt,jk) = 0.
                        IF ( ll_mld ) rmldsec(jipt,jk) = 0.
                        IF ( ll_ice ) ricethicksec(jipt,jk) = 0.
                        IF ( ll_ice ) ricefrasec(jipt,jk) = 0.
                     ELSE
                        tempersec(jipt,jk) = 0.5 * ( temper(ii,ij) + temper(ii,ij+1) )
                        salinesec(jipt,jk) = 0.5 * ( saline(ii,ij) + saline(ii,ij+1) )
                        IF ( ll_ssh ) sshsec (jipt,jk) = 0.5 * ( ssh (ii,ij) + ssh (ii,ij+1) )
                        IF ( ll_mld ) rmldsec(jipt,jk) = 0.5 * ( rmld(ii,ij) + rmld(ii,ij+1) )
                        IF ( ll_ice ) ricethicksec(jipt,jk) = 0.5 * ( ricethick(ii,ij) + ricethick(ii,ij+1) )
                        IF ( ll_ice ) ricefrasec(jipt,jk) = 0.5 * ( ricefra(ii,ij) + ricefra(ii,ij+1) )
                     ENDIF
                     vmeridsec(jipt,jk) = vmerid(ii,ij) * normv_sec(jipt,jsec)
                     e3vsec   (jipt,jk) = e3v   (ii,ij)

                  ENDIF
               ELSEIF ( ii1 == ii ) THEN ! vertical segment
                  vmeridsec(jipt,jk) = 0.
                  e3vsec   (jipt,jk) = 0.
                  IF ( ij1 < ij ) THEN ! southward

                     IF ( MIN( saline(ii,ij) , saline(ii+1,ij) ) == 0. ) THEN
                        tempersec(jipt,jk) = 0. ; salinesec(jipt,jk) = 0.
                        IF ( ll_ssh ) sshsec (jipt,jk) = 0.
                        IF ( ll_mld ) rmldsec(jipt,jk) = 0.
                        IF ( ll_ice ) ricethicksec(jipt,jk) = 0.
                        IF ( ll_ice ) ricefrasec(jipt,jk) = 0.
                     ELSE
                        tempersec(jipt,jk) = 0.5 * ( temper(ii,ij) + temper(ii+1,ij) )
                        salinesec(jipt,jk) = 0.5 * ( saline(ii,ij) + saline(ii+1,ij) )
                        IF ( ll_ssh ) sshsec (jipt,jk) = 0.5 * ( ssh (ii,ij) + ssh (ii+1,ij) )
                        IF ( ll_mld ) rmldsec(jipt,jk) = 0.5 * ( rmld(ii,ij) + rmld(ii+1,ij) )
                        IF ( ll_ice ) ricethicksec(jipt,jk) = 0.5 * ( ricethick(ii,ij) + ricethick(ii+1,ij) )
                        IF ( ll_ice ) ricefrasec(jipt,jk) = 0.5 * ( ricefra(ii,ij) + ricefra(ii+1,ij) )
                     ENDIF
                     uzonalsec(jipt,jk) = uzonal(ii,ij) * normu_sec(jipt,jsec)
                     e3usec   (jipt,jk) = e3u   (ii,ij)

                  ELSE ! northward

                     IF ( MIN( saline(ii,ij+1) , saline(ii+1,ij+1) ) == 0. ) THEN
                        tempersec(jipt,jk) = 0. ; salinesec(jipt,jk) = 0.
                        IF ( ll_ssh ) sshsec (jipt,jk) = 0.
                        IF ( ll_mld ) rmldsec(jipt,jk) = 0.
                        IF ( ll_ice ) ricethicksec(jipt,jk) = 0.
                        IF ( ll_ice ) ricefrasec(jipt,jk) = 0.
                     ELSE
                        tempersec(jipt,jk) = 0.5 * ( temper(ii,ij+1) + temper(ii+1,ij+1) )
                        salinesec(jipt,jk) = 0.5 * ( saline(ii,ij+1) + saline(ii+1,ij+1) )
                        IF ( ll_ssh ) sshsec (jipt,jk) = 0.5 * ( ssh (ii,ij+1) + ssh (ii+1,ij+1) )
                        IF ( ll_mld ) rmldsec(jipt,jk) = 0.5 * ( rmld(ii,ij+1) + rmld(ii+1,ij+1) )
                        IF ( ll_ice ) ricethicksec(jipt,jk) = 0.5 * ( ricethick(ii,ij+1) + ricethick(ii+1,ij+1) )
                        IF ( ll_ice ) ricefrasec(jipt,jk) = 0.5 * ( ricefra(ii,ij+1) + ricefra(ii+1,ij+1) )
                     ENDIF
                     uzonalsec(jipt,jk) = uzonal(ii,ij+1) * normu_sec(jipt,jsec)
                     e3usec   (jipt,jk) = e3u   (ii,ij+1)

                  ENDIF

               ELSE
                  PRINT *, 'problem 2 for JIPT = ', jipt, 'JK=', jk
                  PRINT *, '             I(P2)=',ii1, 'J(P1)=', ii
                  PRINT *, '             J(P2)=',ij1, 'J(P1)=', ij
                  EXIT 
               ENDIF

               ! cumulate transport for barotropic calculation
               dtmp=1.d0* (uzonalsec(jipt,jk) + vmeridsec(jipt,jk))*    &
                   &   (e2usec(jipt,1,jsec )+ e1vsec(jipt,1,jsec ))*    &
                   &   (e3usec(jipt,jk)+ e3vsec(jipt,jk))
               dbarot(jsec)=dbarot(jsec)+dtmp
            END DO

            ! output section variable at level jk, in separated output section files
                        ierr = putvar (ncout(jsec), id_varout(np_tem), tempersec(:,jk), jk, npsec(jsec)-1, 1, ktime=jt )
                        ierr = putvar (ncout(jsec), id_varout(np_sal), salinesec(:,jk), jk, npsec(jsec)-1, 1, ktime=jt )
                        ierr = putvar (ncout(jsec), id_varout(np_una), uzonalsec(:,jk), jk, npsec(jsec)-1, 1, ktime=jt )
                        ierr = putvar (ncout(jsec), id_varout(np_vna), vmeridsec(:,jk), jk, npsec(jsec)-1, 1, ktime=jt )
            ! along-track normal velocity, horiz. and vert. resolution, and mask
            zvmod(:,1)= uzonalsec(:,jk) + vmeridsec(:,jk)
                        ierr = putvar (ncout(jsec), id_varout(np_vmod    ), zvmod(:,1),                jk, npsec(jsec)-1, 1, ktime=jt ) 
            IF (ll_ssh) ierr = putvar (ncout(jsec), id_varout(np_ssh     ), sshsec (:,jk),              1, npsec(jsec)-1, 1, ktime=jt )
            IF (ll_mld) ierr = putvar (ncout(jsec), id_varout(np_mld     ), rmldsec(:,jk),              1, npsec(jsec)-1, 1, ktime=jt )
            IF (ll_ice) ierr = putvar (ncout(jsec), id_varout(np_icethick), ricethicksec(:,jk),         1, npsec(jsec)-1, 1, ktime=jt )
            IF (ll_ice) ierr = putvar (ncout(jsec), id_varout(np_icefra  ), ricefrasec(:,jk),           1, npsec(jsec)-1, 1, ktime=jt )
            IF (lvt   ) ierr = putvar (ncout(jsec), id_varout(np_vt      ), zvmod(:,1)*tempersec(:,jk),jk, npsec(jsec)-1, 1, ktime=jt )
            IF (lvt   ) ierr = putvar (ncout(jsec), id_varout(np_vs      ), zvmod(:,1)*salinesec(:,jk),jk, npsec(jsec)-1, 1, ktime=jt )

            IF ( jt == 1 ) THEN   ! output of time independent variables at first time step only
               ! save a mask of the section
               vmasksec(:,:) = 1.
               WHERE( salinesec(:,:) == 0. ) vmasksec(:,:) = 0.

               ierr = putvar (ncout(jsec), id_varout(np_e3un), e3usec(:,jk),                jk, npsec(jsec)-1, 1 )
               ierr = putvar (ncout(jsec), id_varout(np_e3vn), e3vsec(:,jk),                jk, npsec(jsec)-1, 1 )
               ierr = putvar (ncout(jsec), id_varout(np_e3v ), e3usec(:,jk)+e3vsec(:,jk),   jk, npsec(jsec)-1, 1 )
               ierr = putvar (ncout(jsec), id_varout(np_vmsk), vmasksec(:,jk),              jk, npsec(jsec)-1, 1 )
            ENDIF
         END DO  ! section
      END DO  ! levels

      ! print and output barotropic transport once all levels have been processed
      DO jsec = 1, nsec
         PRINT 9010, TRIM(csection(jsec)),' BAROTROPIC TRANSPORT at time ',jt,' = ', dbarot(jsec)/1.d6, ' Sv.'
         ierr  = putvar0d ( ncout(jsec), id_varout(np_baro), REAL(dbarot(jsec)/1.d6), ktime = jt    )
      ENDDO

   END DO  ! time 

   ! close all output files
   DO jsec = 1, nsec
      ierr = closeout(ncout(jsec))
   ENDDO
9010 FORMAT ( "Section :",a15,a,i3,a3,f8.3,a4)

CONTAINS 
   SUBROUTINE CreateOutputFile(ksec)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE CreateOutputFile  ***
      !!
      !! ** Purpose :  Perform output file creation with all the variables 
      !!
      !! ** Method  :  Move this part of the code in a subroutine for clarity
      !!               All variables are global.  
      !!
      !!----------------------------------------------------------------------
      INTEGER(KIND=4), INTENT(in) :: ksec  ! section index
      INTEGER(KIND=4)             :: ivar  ! variable index
      !!----------------------------------------------------------------------
      ivar = 1

      stypvar%scale_factor= 1.
      stypvar%add_offset= 0.
      stypvar%savelog10= 0.
      stypvar%rmissing_value=0.
      stypvar%conline_operation='N/A'

      ! define new variables for output 
      np_tem = ivar
      stypvar(ivar)%cname       = cn_votemper
      stypvar(ivar)%cunits      = 'deg C'
      stypvar(ivar)%valid_min   = -2.
      stypvar(ivar)%valid_max   = 40.
      stypvar(ivar)%clong_name  = 'Temperature along '//TRIM(csection(ksec))//' section'
      stypvar(ivar)%cshort_name = cn_votemper
      stypvar(ivar)%caxis       = 'TZX'
      ipk(ivar)                 = npk
      ivar = ivar + 1

      np_sal = ivar
      stypvar(ivar)%cname       = cn_vosaline
      stypvar(ivar)%cunits      = 'PSU'
      stypvar(ivar)%valid_min   = 0.
      stypvar(ivar)%valid_max   = 50.
      stypvar(ivar)%clong_name  = 'Salinity along '//TRIM(csection(ksec))//' section'
      stypvar(ivar)%cshort_name = cn_vosaline
      stypvar(ivar)%caxis       = 'TZX'
      ipk(ivar)                 = npk
      ivar = ivar + 1

      np_una = ivar
      stypvar(ivar)%cname       = TRIM(cn_vozocrtx)//'_native'
      stypvar(ivar)%cunits      = 'm.s-1'
      stypvar(ivar)%valid_min   = -20.
      stypvar(ivar)%valid_max   = 20.
      stypvar(ivar)%clong_name  = 'Zonal velocity along '//TRIM(csection(ksec))//' section'
      stypvar(ivar)%cshort_name = TRIM(cn_vozocrtx)//'_native'
      stypvar(ivar)%caxis       = 'TZX'
      ipk(ivar)                 = npk
      ivar = ivar + 1

      np_vna = ivar
      stypvar(ivar)%cname       = TRIM(cn_vomecrty)//'_native'
      stypvar(ivar)%cunits      = 'm.s-1'
      stypvar(ivar)%valid_min   = -20.
      stypvar(ivar)%valid_max   = 20.
      stypvar(ivar)%clong_name  = 'Meridionnal velocity along '//TRIM(csection(ksec))//' section'
      stypvar(ivar)%cshort_name = TRIM(cn_vomecrty)//'_native'
      stypvar(ivar)%caxis       = 'TZX'
      ipk(ivar)                 = npk
      ivar = ivar + 1

      np_isec = ivar
      stypvar(ivar)%cname       = 'isec'
      stypvar(ivar)%valid_min   = 1.
      stypvar(ivar)%valid_max   = npiglo 
      stypvar(ivar)%caxis       = 'TX'
      ipk(ivar)                 = 1
      ivar = ivar + 1

      np_jsec = ivar
      stypvar(ivar)%cname       = 'jsec'
      stypvar(ivar)%valid_min   = 1.
      stypvar(ivar)%valid_max   = npjglo 
      stypvar(ivar)%caxis       = 'TX'
      ipk(ivar)                 = 1
      ivar = ivar + 1

      np_e2vn = ivar
      stypvar(ivar)%cname       = TRIM(cn_ve2u)//'_native'
      stypvar(ivar)%valid_min   = 1.
      stypvar(ivar)%valid_max   = 200000.
      stypvar(ivar)%caxis       = 'TX'
      ipk(ivar)                 = 1
      ivar = ivar + 1

      np_e1vn = ivar
      stypvar(ivar)%cname       = TRIM(cn_ve1v)//'_native'
      stypvar(ivar)%valid_min   = 1.
      stypvar(ivar)%valid_max   = 200000.
      stypvar(ivar)%caxis       = 'TX'
      ipk(ivar)                 = 1
      ivar = ivar + 1

      np_e3un = ivar
      stypvar(ivar)%cname       = 'e3u_native'
      stypvar(ivar)%valid_min   = 1.
      stypvar(ivar)%valid_max   = 200000.
      stypvar(ivar)%caxis       = 'TZX'
      ipk(ivar)                 =  npk
      ivar = ivar + 1

      np_e3vn = ivar
      stypvar(ivar)%cname      = 'e3v_native'
      stypvar(ivar)%valid_min   = 1.
      stypvar(ivar)%valid_max   = 200000.
      stypvar(ivar)%caxis      = 'TZX'
      ipk(ivar)                =  npk
      ivar = ivar + 1

      np_vmod = ivar
      stypvar(ivar)%cname       = cn_vomecrty
      stypvar(ivar)%cunits      = 'm.s-1'
      stypvar(ivar)%valid_min   = -20.
      stypvar(ivar)%valid_max   = 20.
      stypvar(ivar)%clong_name  = 'Normal velocity along '//TRIM(csection(ksec))//' section'
      stypvar(ivar)%cshort_name = cn_vomecrty
      stypvar(ivar)%caxis       = 'TZX'
      ipk(ivar)                 =  npk
      ivar = ivar + 1

      np_e1v = ivar
      stypvar(ivar)%cname       = cn_ve1v
      stypvar(ivar)%cunits      = 'm'
      stypvar(ivar)%valid_min   = 0.
      stypvar(ivar)%valid_max   = 1000000.
      stypvar(ivar)%clong_name  = 'Local horiz. resolution along '//TRIM(csection(ksec))//' section'
      stypvar(ivar)%cshort_name = cn_ve1v
      stypvar(ivar)%caxis       = 'TX'
      ipk(ivar)                 = 1
      ivar = ivar + 1

      np_e3v = ivar
      stypvar(ivar)%cname       = 'e3v_ps'
      stypvar(ivar)%cunits      = 'm'
      stypvar(ivar)%valid_min   = 0.
      stypvar(ivar)%valid_max   = 100000000.
      stypvar(ivar)%clong_name  = 'Local vert. resolution along '//TRIM(csection(ksec))//' section'
      stypvar(ivar)%cshort_name = 'e3v_ps'
      stypvar(ivar)%caxis       = 'TZX'
      ipk(ivar)                 =  npk
      ivar = ivar + 1

      np_vmsk = ivar
      stypvar(ivar)%cname       = 'vmask'
      stypvar(ivar)%cunits      ='1/0'
      stypvar(ivar)%valid_min   = 0.
      stypvar(ivar)%valid_max   = 1.
      stypvar(ivar)%rmissing_value = 9999.
      stypvar(ivar)%clong_name  ='Mask along '//TRIM(csection(ksec))//' section'
      stypvar(ivar)%cshort_name = 'vmask'
      stypvar(ivar)%caxis       = 'TZX'
      ipk(ivar)                 =  npk
      ivar = ivar + 1

      np_baro = ivar
      stypvar(ivar)%cname       = 'barotrop_'//TRIM(csection(ksec))
      stypvar(ivar)%cunits      ='Sv'
      stypvar(ivar)%valid_min   = -500.
      stypvar(ivar)%valid_max   = 500.
      stypvar(ivar)%rmissing_value = -99999.
      stypvar(ivar)%clong_name  = 'Barotropic_transport for '//TRIM(csection(ksec))//' section'
      stypvar(ivar)%cshort_name = 'barotrop_'//TRIM(csection(ksec))
      stypvar(ivar)%caxis       = 'T'
      ipk(ivar)                 =  -1
      ivar = ivar + 1

      np_bat = ivar
      stypvar(ivar)%cname       = 'Bathymetry'
      stypvar(ivar)%cunits      = 'm'
      stypvar(ivar)%valid_min   = 0.
      stypvar(ivar)%valid_max   = 1000000.
      stypvar(ivar)%clong_name  = 'Bathymetry along '//TRIM(csection(ksec))//' section'
      stypvar(ivar)%cshort_name = 'Bathymetry'
      stypvar(ivar)%caxis       = 'TX'
      ipk(ivar)                 = 1
      ivar = ivar + 1

      IF ( lssh ) THEN
         np_ssh = ivar
         stypvar(ivar)%cname       = cn_sossheig
         stypvar(ivar)%cunits      = 'm'
         stypvar(ivar)%valid_min   = 0.
         stypvar(ivar)%valid_max   = 1000000.
         stypvar(ivar)%clong_name  = 'SSH  along '//TRIM(csection(ksec))//' section'
         stypvar(ivar)%cshort_name = cn_sossheig
         stypvar(ivar)%caxis       = 'TX'
         ipk(ivar)                 = 1
         ivar = ivar + 1
      ENDIF

      IF ( lmld ) THEN
         np_mld = ivar
         stypvar(ivar)%cname       = cn_somxl010
         stypvar(ivar)%cunits      = 'm'
         stypvar(ivar)%valid_min   = 0.
         stypvar(ivar)%valid_max   = 100000.
         stypvar(ivar)%clong_name  = 'Mixed Layer Depth 0.01  along '//TRIM(csection(ksec))//' section'
         stypvar(ivar)%cshort_name = cn_somxl010
         stypvar(ivar)%caxis       = 'TX'
         ipk(ivar)                 = 1
         ivar = ivar + 1
      ENDIF

      IF ( lice ) THEN
         np_icethick = ivar
         stypvar(ivar)%cname       = cn_iicethic
         stypvar(ivar)%cunits      = 'm'
         stypvar(ivar)%valid_min   = -10000.
         stypvar(ivar)%valid_max   = 1000000.
         stypvar(ivar)%clong_name  = 'icethick along '//TRIM(csection(ksec))//' section'
         stypvar(ivar)%cshort_name = cn_iicethic
         stypvar(ivar)%caxis       = 'TX'
         ipk(ivar)                 = 1
         ivar = ivar + 1
   
         np_icefra = ivar
         stypvar(ivar)%cname       = cn_ileadfra
         stypvar(ivar)%cunits      = 'm'
         stypvar(ivar)%valid_min   = -10000.
         stypvar(ivar)%valid_max   = 1000000.
         stypvar(ivar)%clong_name  = 'icefra along '//TRIM(csection(ksec))//' section'
         stypvar(ivar)%cshort_name = cn_ileadfra
         stypvar(ivar)%caxis       = 'TX'
         ipk(ivar)                 = 1
         ivar = ivar + 1
      ENDIF
     
      IF ( lvt ) THEN
         np_vt = ivar
         stypvar(ivar)%cname       = cn_vomevt
         stypvar(ivar)%cunits      = 'C.m/s'
         stypvar(ivar)%valid_min   = -1000000.
         stypvar(ivar)%valid_max   = 1000000.
         stypvar(ivar)%clong_name  = 'VT product along '//TRIM(csection(ksec))//' section'
         stypvar(ivar)%cshort_name = cn_vomevt
         stypvar(ivar)%caxis       = 'TZX'
         ipk(ivar)                 = npk
         ivar = ivar + 1

         np_vs = ivar
         stypvar(ivar)%cname       = cn_vomevs
         stypvar(ivar)%cunits      = 'PSU.m/s'
         stypvar(ivar)%valid_min   = -1000000.
         stypvar(ivar)%valid_max   = 1000000.
         stypvar(ivar)%clong_name  = 'VS product along '//TRIM(csection(ksec))//' section'
         stypvar(ivar)%cshort_name = cn_vomevs
         stypvar(ivar)%caxis       = 'TZX'
         ipk(ivar)                 = npk
         ivar = ivar + 1
      ENDIF

      ! create output fileset
      ncout(ksec) = create      (cf_out, cf_tfil, npsec(ksec),  1, npk, cdep='deptht'                      )
      ierr  = createvar   (ncout(ksec),  stypvar, nvar,  ipk, id_varout                                    )
      ierr  = putheadervar(ncout(ksec),  cf_tfil, npsec(ksec)-1,  1, npk, pnavlon=rlonsec, pnavlat=rlatsec )
      tim   = getvar1d    (cf_tfil, 'time_counter', npt                                                    )
      ierr  = putvar1d    (ncout(ksec), tim, npt, 'T'                                                      )

   END SUBROUTINE CreateOutputFile

   SUBROUTINE ParseFiles (cdum)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE ParseFiles  ***
      !!
      !! ** Purpose :  Decode -f option from command line 
      !!
      !! ** Method  :  look for , in the argument string and set the number of 
      !!         sections (nsec), allocate cf_sec array and fill it with the 
      !!         decoded  names.
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: cdum

      CHARACTER(LEN=2048), DIMENSION(100) :: cl_dum  ! 100 is arbitrary
      INTEGER  :: ji
      INTEGER  :: inchar,  i1=1
      !!----------------------------------------------------------------------

      inchar= LEN(TRIM(cdum))
      ! scan the input string and look for ',' as separator
      DO ji=1,inchar
         IF ( cdum(ji:ji) == ',' ) THEN
            cl_dum(nsec) = cdum(i1:ji-1)
            i1=ji+1
            nsec=nsec+1
         ENDIF
      ENDDO

      ! last name of the list does not have a ','
      cl_dum(nsec) = cdum(i1:inchar)

      ALLOCATE ( cf_sec(nsec) )

      DO ji=1, nsec
         cf_sec(ji) = cl_dum(ji)
      ENDDO
   END SUBROUTINE ParseFiles

   SUBROUTINE FancyPrint(ksec)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE FancyPrint  ***
      !!
      !! ** Purpose :   perform Fancy Print for sections definitions 
      !!
      !!----------------------------------------------------------------------
      INTEGER(KIND=4), INTENT(in) :: ksec
      !!----------------------------------------------------------------------

      PRINT 9005
      PRINT 9006, TRIM(csection(ksec))
      PRINT 9005
      PRINT 9000
      PRINT 9005

      OPEN(numout, file=cf_secdat )   ! open section.dat file for output

      DO jleg = 1, nsta(ksec) -1
         ! start point 
         ii = iista(jleg,ksec)  ; ij = ijsta(jleg,ksec)
         ztmp = MIN (rmbat(ii,ij), rmbat(ii+1,ij), rmbat(ii+1,ij+1), rmbat(ii,ij+1) ) 
         IF ( ztmp == 0. ) THEN 
            cstar = 'LAND'
         ELSE
            cstar = ' SEA'
         ENDIF

         ! end  point 
         ii1 = iista(jleg+1,ksec)  ; ij1 = ijsta(jleg+1,ksec)
         ztmp = MIN (rmbat(ii1,ij1), rmbat(ii1+1,ij1), rmbat(ii1+1,ij1+1), rmbat(ii1,ij1+1) )
         IF ( ztmp == 0. ) THEN
            cend = 'LAND'
         ELSE
            cend = ' SEA'
         ENDIF
         PRINT 9001, jleg, rlatsta(jleg,ksec), rlonsta(jleg,ksec),    rlatsta(jleg+1,ksec), rlonsta(jleg+1,ksec)
         PRINT 9002,          ii,ij,                        ii1,ij1
         PRINT 9003, rlatf(ii,ij), rlonf(ii,ij),            rlatf(ii1,ij1), rlonf(ii1,ij1)
         PRINT 9004, TRIM(cstar),                           TRIM(cend)
         PRINT 9005
         WRITE(numout,'(i2.2,"_",a)') jleg, TRIM(csection(ksec))
         WRITE(numout,*) ii, ii1, ij, ij1
      ENDDO

      WRITE(numout,'("EOF")')
      CLOSE(numout)

9000  FORMAT ("  |  Leg #  |    Start point     |      End point     | ")
9001  FORMAT ("  |   ",i3,"   | ", f6.2," N ", f7.2, " E | ", f6.2," N ", f7.2, " E |" )
9002  FORMAT ("  |        F| I =", i5,", J =",i5 " | I =", i5,", J =",i5 " |" )
9003  FORMAT ("  |      mod| ", f6.2," N ", f7.2, " E | ", f6.2," N ", f7.2, " E |" )
9004  FORMAT ("  |         |      ",a4,"          |     ",a4,"           | ")
9005  FORMAT ("  |---------|--------------------|--------------------| ")
9006  FORMAT ("  |    Section : ",a33,"    |" )

   END SUBROUTINE FancyPrint

END PROGRAM cdf_xtract_brokenline
