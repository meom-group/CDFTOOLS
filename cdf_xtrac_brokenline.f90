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

   INTEGER(KIND=4) :: jleg, jt, jk,  jipt, jvar      ! dummy loop index
   INTEGER(KIND=4) :: narg, iargc, ijarg, ifree      ! command line
   INTEGER(KIND=4) :: numin=10                       ! logical unit for input section file
   INTEGER(KIND=4) :: numout=11                      ! logical unit for output section.dat (used in cdftransport)
   INTEGER(KIND=4) :: npiglo, npjglo, npk, npt       ! size of the domain
   INTEGER(KIND=4) :: iimin, iimax, ijmin, ijmax     ! ending points of a leg in model I J
   INTEGER(KIND=4) :: ii, ij, ii1, ij1, ipoint       ! working integer
   INTEGER(KIND=4) :: ierr, ncout                    ! Netcdf error and ncid
   INTEGER(KIND=4) :: nvar = 16                      ! number of output variables
   INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout  ! netcdf output stuff

   ! broken line definition
   INTEGER(KIND=4) :: nsta=5                         ! number of points defining the broken line
   INTEGER(KIND=4) :: nsec=0                         ! total number of points on the broken line
   INTEGER(KIND=4) :: nn                             ! working integer (number of points in a leg)
   INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: norm_u, norm_v   ! velocity normalization per leg
   INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: normu_sec, normv_sec ! velocity normalization per section
   INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: iista, ijsta     ! I,J position of the point on the broken line
   INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: ikeepn           ! Number of points per leg
   INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: iisec, ijsec     ! F-index of points on the broken line
   INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: iilegs, ijlegs   ! F-index of points on the broken line per leg

   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: rlonsta, rlatsta    ! Geographic position defining legs
   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: rxx, ryy            ! leg i j index of F points
   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                 ! Model time array

   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1v, e3v            ! V point relevant metric
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2u, e3u            ! U point relevant metric
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: hdepw               ! model bathymetry
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlonu, rlatu        ! model long and lat of U points
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlonv, rlatv        ! model long and lat of U points
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlonf, rlatf        ! model long and lat of F points
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: temper, saline      ! model Temperature and salinity
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: uzonal, vmerid      ! model zonal and meridional velocity
   ! along section array (dimension x,z or x,1 )
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tempersec, salinesec, uzonalsec, vmeridsec
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlonsec, rlatsec, risec, rjsec
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1vsec, e2usec
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e3usec, e3vsec
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: batsec
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: vmasksec
   REAL(KIND=4)                              :: xmin, xmax, ymin, ymax !
   REAL(KIND=4)                              :: ztmp

   REAL(KIND=8)                              :: dtmp, dbarot  ! for barotropic transport computation

   CHARACTER(LEN=80) :: cf_tfil , cf_ufil, cf_vfil   ! input T U V files
   CHARACTER(LEN=80) :: cf_out                       ! output file
   CHARACTER(LEN=80) :: cf_secdat                    ! output section file (suitable for cdftransport or cdfsigtrp)
   CHARACTER(LEN=80) :: cf_sec                       ! input section file
   CHARACTER(LEN=80) :: csection                     ! section name
   CHARACTER(LEN=80) :: cverb='n'                    ! verbose key for findij
   CHARACTER(LEN=80) :: cldum, cstar, cend           ! dummy character variable

   LOGICAL  :: lchk                                  ! flag for missing files
   LOGICAL  :: lverbose = .FALSE.                    ! flag for verbosity
   LOGICAL  :: lsecfile = .FALSE.                    ! flag for input section file

   TYPE (variable), DIMENSION(:), ALLOCATABLE :: stypvar  ! variable definition and attributes
   !!----------------------------------------------------------------------
   ! 1. : Initialization
   ! --------------------
   CALL ReadCdfNames()

   ! check argument number and show usage if necessary
   narg = iargc()
   IF ( narg < 3 ) THEN
      PRINT *,' usage :  cdf_xtrac_brokenline T-file U-file V-file [-f section_file ] ...'
      PRINT *,'                     ... [-verbose]'
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
      PRINT *,'      T-file :  model gridT file '
      PRINT *,'      U-file :  model gridU file '
      PRINT *,'      V-file :  model gridV file '
      PRINT *,'      ' 
      PRINT *,'     OPTIONS :'
      PRINT *,'      -f section_file : provide a file for section definition.'
      PRINT *,'             section_file is an ascii file as follows:'
      PRINT *,'             * line #1 : name of the section (e.g. ovide). '
      PRINT *,'                  Will be used for naming the output file.'
      PRINT *,'             * line #2 : number of points defining the broken line.'
      PRINT *,'             * line #3-end : a pair of Longitude latitude values defining'
      PRINT *,'                   the points. If not supplied, use hard-coded information'
      PRINT *,'                   for OVIDE section. A comment can be added at the end of'
      PRINT *,'                   of the lines, using a # as separator'
      PRINT *,'      -verbose : increase verbosity  ' 
      PRINT *,'     '
      PRINT *,'     REQUIRED FILES :'
      PRINT *,'      ', TRIM(cn_fhgr),' and ',TRIM(cn_fzgr),' must be in the current directory ' 
      PRINT *,'      '
      PRINT *,'     OUTPUT : '
      PRINT *,'       netcdf file : section_name.nc'
      PRINT *,'         variables : temperature, salinity, normal velocity, pseudo V metrics,'
      PRINT *,'                     mask, barotropic transport, bathymetry of velocity points.'
      PRINT *,'       ASCII file : section_name_section.dat usefull for cdftransport '
      PRINT *,'      '
      PRINT *,'     SEE ALSO :'
      PRINT *,'        cdftransport, cdfmoc, cdfmocsig. This tool replace cdfovide.' 
      PRINT *,'      '
      STOP
   ENDIF

   ! Decode command line
   ijarg = 1 ; ifree=0
   DO WHILE ( ijarg <= narg ) 
      CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1
      SELECT CASE  ( cldum   )
      CASE ( '-verbose' ) ; lverbose=.true.  ; cverb='y'
      CASE ( '-f' )       ;  CALL getarg(ijarg, cf_sec) ; ijarg = ijarg + 1 ; lsecfile=.TRUE.
      CASE DEFAULT 
         ifree = ifree + 1
         SELECT CASE ( ifree )
         CASE ( 1 ) ; cf_tfil = cldum
         CASE ( 2 ) ; cf_ufil = cldum
         CASE ( 3 ) ; cf_vfil = cldum
         END SELECT
      END SELECT
   ENDDO

   ! check file existence
   lchk = chkfile(cn_fhgr )
   lchk = chkfile(cn_fzgr ) .OR. lchk
   lchk = chkfile(cf_tfil ) .OR. lchk
   lchk = chkfile(cf_ufil ) .OR. lchk
   lchk = chkfile(cf_vfil ) .OR. lchk
   IF ( lsecfile ) lchk = chkfile(cf_sec ) .OR. lchk
   IF ( lchk     ) STOP ! missing files

   ! read section file if required
   IF ( lsecfile ) THEN
      OPEN(numin, file=cf_sec )
      READ(numin,'(a)') csection
      READ(numin,*    ) nsta
      ALLOCATE ( iista(nsta), ijsta(nsta), ikeepn(nsta -1 )  )
      ALLOCATE ( rlonsta(nsta), rlatsta(nsta) )
      DO jipt = 1, nsta
         READ(numin, * ) rlonsta(jipt), rlatsta(jipt)
      ENDDO
      CLOSE (numin)
   ELSE    ! default to OVIDE section
      nsta = 5
      csection = 'ovide'
      ALLOCATE ( iista(nsta), ijsta(nsta), ikeepn(nsta -1 )  )
      ALLOCATE ( rlonsta(nsta), rlatsta(nsta) )  
      ! R. Dussin : Location of leg points that define the 3 legs of OVIDE section
      !rlonsta(1) = -43.00 ; rlatsta(1) = 60.60    ! Greenland
      !rlonsta(2) = -31.30 ; rlatsta(2) = 58.90    ! Reykjanes Ridge
      !rlonsta(3) = -12.65 ; rlatsta(3) = 40.33    ! Off Portugal
      !rlonsta(4) =  -8.70 ; rlatsta(4) = 40.33    ! Lisboa

      ! D. Desbruyeres : Location of leg points that define the 4 legs of the OVIDE section
      rlonsta(1) = -43.70 ; rlatsta(1) = 59.90    ! 
      rlonsta(2) = -30.30 ; rlatsta(2) = 58.90    ! 
      rlonsta(3) = -19.40 ; rlatsta(3) = 44.90    ! 
      rlonsta(4) = -12.65 ; rlatsta(4) = 40.33    ! 
      rlonsta(5) = -08.70 ; rlatsta(5) = 40.33    ! 
   ENDIF

   cf_out=TRIM(csection)//'.nc'
   cf_secdat=TRIM(csection)//'_section.dat'

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

   ALLOCATE ( iilegs(nsta-1, npiglo+npjglo), ijlegs(nsta-1, npiglo+npjglo) )
   ALLOCATE ( norm_u(nsta-1) , norm_v(nsta-1) )
   ALLOCATE ( rxx(npiglo+npjglo), ryy(npiglo+npjglo) )
   ALLOCATE ( tim (npt) )

   iilegs = 0  ; ijlegs = 0

   ! loop on the legs
   DO jleg = 1, nsta-1

      xmin = rlonsta(jleg  )
      ymin = rlatsta(jleg  )
      xmax = rlonsta(jleg+1)
      ymax = rlatsta(jleg+1)

      ! return ending points of a leg in I J model coordinates
      CALL cdf_findij ( xmin, xmax, ymin, ymax, iimin, iimax, ijmin, ijmax, &
           &            cd_coord=cn_fhgr, cd_point='F', cd_verbose=cverb)

      ! save leg information
      iista(jleg  ) = iimin
      ijsta(jleg  ) = ijmin
      iista(jleg+1) = iimax
      ijsta(jleg+1) = ijmax

      ! find the broken line between P1 (iimin,ijmin) and P2 (iimax, ijmax)
      CALL broken_line( iimin, iimax, ijmin, ijmax, rxx, ryy, nn, npiglo, npjglo, norm_u(jleg), norm_v(jleg) )
      ikeepn(jleg) = nn  ! number of points (F) on leg jleg
      nsec = nsec + nn   ! total number of points (F) on the broken line

      IF ( lverbose) PRINT *, 'Leg ', jleg,' : npoints : ', nn

      IF ( jleg == 1 ) THEN
        IF (rxx(1) < rxx(nn) ) THEN ! leg is oriented eastward
           iilegs(jleg,1:nn)=rxx(1:nn)
           ijlegs(jleg,1:nn)=ryy(1:nn)
        ELSE                        ! leg is oriented westward
           iilegs(jleg,1:nn)=rxx(nn:1:-1)
           ijlegs(jleg,1:nn)=ryy(nn:1:-1)
        END IF
      ELSE  ! check the continuity between legs
      IF ( iilegs(jleg-1, ikeepn(jleg-1)) == rxx(1) ) THEN  ! continuity
         iilegs(jleg,1:nn)=rxx(1:nn)
         ijlegs(jleg,1:nn)=ryy(1:nn)
      ELSE                          ! reverse sense
         iilegs(jleg,1:nn)=rxx(nn:1:-1)
         ijlegs(jleg,1:nn)=ryy(nn:1:-1)
      END IF
      ENDIF

      IF ( lverbose) THEN
         PRINT *, 'Leg  rxx   ryy '
         DO jk = 1, nn
            PRINT *, jleg, iilegs(jleg,jk), ijlegs(jleg,jk) ,rxx(jk), ryy(jk)
         END DO
      ENDIF
   END DO !! loop on the legs

   ! 3. : Extraction along the legs
   ! ------------------------------
   ALLOCATE (iisec(nsec), ijsec(nsec)) 
   ALLOCATE (normu_sec(nsec), normv_sec(nsec)) 

   ipoint = 1
   DO jleg=1, nsta-1      ! loop on legs 
      ipoint = ipoint -1  ! trick to avoid repetition of points in between legs
      DO jipt=1, ikeepn(jleg) 
         ipoint = ipoint + 1
         iisec(ipoint)=iilegs(jleg,jipt)  ! i-index
         ijsec(ipoint)=ijlegs(jleg,jipt)  ! j-index
         normu_sec(ipoint) = norm_u(jleg)
         normv_sec(ipoint) = norm_v(jleg)
      END DO
   END DO

   ! adjust nsec to its real value ( 2nd part of the trick)
   nsec = ipoint

   ! input fields
   ALLOCATE(rlonu(npiglo,npjglo), rlatu(npiglo,npjglo))
   ALLOCATE(rlonv(npiglo,npjglo), rlatv(npiglo,npjglo))
   ALLOCATE(rlonf(npiglo,npjglo), rlatf(npiglo,npjglo))
   ALLOCATE(temper(npiglo,npjglo), saline(npiglo,npjglo))
   ALLOCATE(uzonal(npiglo,npjglo), vmerid(npiglo,npjglo))
   ALLOCATE(e1v(npiglo,npjglo))
   ALLOCATE(e2u(npiglo,npjglo))
   ALLOCATE(e3u(npiglo,npjglo), e3v(npiglo,npjglo))
   ALLOCATE(hdepw(npiglo,npjglo) )

   ! output fields
   ALLOCATE(rlonsec(nsec,1), rlatsec(nsec,1) )
   ALLOCATE(risec  (nsec,1), rjsec  (nsec,1) )
   ALLOCATE(e2usec(nsec-1,1), e3usec(nsec-1,npk) )
   ALLOCATE(e1vsec(nsec-1,1), e3vsec(nsec-1,npk) )
   ALLOCATE(batsec(nsec-1,1) )
   ALLOCATE(tempersec(nsec-1,npk), salinesec(nsec-1,npk) )
   ALLOCATE(uzonalsec(nsec-1,npk), vmeridsec(nsec-1,npk) )
   ALLOCATE(vmasksec (nsec-1,npk))

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
   hdepw(:,:) = getvar(cn_fzgr, cn_mbathy,1, npiglo, npjglo)

  ! now that we know the model grid and bathy do fancy print of the legs.
  PRINT 9005
  PRINT 9000
  PRINT 9005

  OPEN(numout, file=cf_secdat )   ! open section.dat file for output

  DO jleg = 1, nsta -1
     ! start point 
     ii = iista(jleg)  ; ij = ijsta(jleg)
     ztmp = MIN (hdepw(ii,ij), hdepw(ii+1,ij), hdepw(ii+1,ij+1), hdepw(ii,ij+1) ) 
     IF ( ztmp == 0. ) THEN 
       cstar = 'LAND'
     ELSE
       cstar = ' SEA'
     ENDIF

     ! end  point 
     ii1 = iista(jleg+1)  ; ij1 = ijsta(jleg+1)
     ztmp = MIN (hdepw(ii1,ij1), hdepw(ii1+1,ij1), hdepw(ii1+1,ij1+1), hdepw(ii1,ij1+1) )
     IF ( ztmp == 0. ) THEN
       cend = 'LAND'
     ELSE
       cend = ' SEA'
     ENDIF
     PRINT 9001, jleg, rlatsta(jleg), rlonsta(jleg),    rlatsta(jleg+1), rlonsta(jleg+1)
     PRINT 9002,          ii,ij,                        ii1,ij1
     PRINT 9003, rlatf(ii,ij), rlonf(ii,ij),            rlatf(ii1,ij1), rlonf(ii1,ij1)
     PRINT 9004, TRIM(cstar),                           TRIM(cend)
     PRINT 9005
     WRITE(numout,'("leg_",i2.2)') jleg
     WRITE(numout,*) ii, ii1, ij, ij1
  ENDDO

  WRITE(numout,'("EOF")')
  CLOSE(numout)

9000 FORMAT ("  |  Leg #  |    Start point     |      End point     | ")
9001 FORMAT ("  |   ",i3,"   | ", f6.2," N ", f7.2, " E | ", f6.2," N ", f7.2, " E |" )
9002 FORMAT ("  |        F| I =", i5,", J =",i5 " | I =", i5,", J =",i5 " |" )
9003 FORMAT ("  |      mod| ", f6.2," N ", f7.2, " E | ", f6.2," N ", f7.2, " E |" )
9004 FORMAT ("  |         |      ",a4,"          |     ",a4,"           | ")
9005 FORMAT ("  |---------|--------------------|--------------------| ")

   ! now set hdepw to its true value
   hdepw(:,:) = getvar(cn_fzgr, cn_hdepw, 1, npiglo, npjglo)

   ! loop on 2d arrays
   DO jipt = 1,nsec
      ii = iisec(jipt)
      ij = ijsec(jipt)

      risec  (jipt,1) = ii
      rjsec  (jipt,1) = ij
   END DO

   DO jipt=1,nsec-1
      ii  = iisec(jipt  ) ; ij  = ijsec(jipt  )
      ii1 = iisec(jipt+1) ; ij1 = ijsec(jipt+1)
      IF ( ij1 == ij ) THEN ! horizontal segment
         e2usec(jipt,1) = 0.
         IF ( ii1 > ii ) THEN ! eastward
            e1vsec (jipt,1) = e1v  (ii+1,ij)
            rlonsec(jipt,1) = rlonv(ii+1,ij)
            rlatsec(jipt,1) = rlatv(ii+1,ij)
            batsec (jipt,1) = MIN( hdepw(ii+1,ij),  hdepw(ii+1,ij+1) )
         ELSE
            e1vsec (jipt,1) = e1v  (ii,ij)
            rlonsec(jipt,1) = rlonv(ii,ij)
            rlatsec(jipt,1) = rlatv(ii,ij)
            batsec (jipt,1) = MIN( hdepw(ii,ij),  hdepw(ii,ij+1) )
         ENDIF

      ELSEIF ( ii1 == ii ) THEN ! vertical segment
         e1vsec(jipt,1) = 0.
         IF ( ij1 < ij ) THEN ! southward
            e2usec (jipt,1) = e2u  (ii,ij)
            rlonsec(jipt,1) = rlonu(ii,ij)
            rlatsec(jipt,1) = rlatu(ii,ij)
            batsec (jipt,1) = MIN( hdepw(ii,ij),  hdepw(ii+1,ij) )
         ELSE
            e2usec (jipt,1) = e2u  (ii,ij+1)
            rlonsec(jipt,1) = rlonu(ii,ij+1)
            rlatsec(jipt,1) = rlatu(ii,ij+1)
            batsec (jipt,1) = MIN( hdepw(ii,ij+1),  hdepw(ii+1,ij+1) )
         ENDIF

      ELSE
         PRINT *, 'problem 1 for JIPT = ', jipt
         PRINT *, '             I(P2)=',ii1, 'J(P1)=', ii
         PRINT *, '             J(P2)=',ij1, 'J(P1)=', ij
         exit 
      ENDIF
   END DO

   ! Prepare output file ( here because rlonsec and rlatsec required )
   ALLOCATE ( stypvar(nvar), ipk(nvar), id_varout(nvar) )
   CALL CreateOutputFile 

   ierr = putvar (ncout, id_varout(5), risec(:,1),                            1,  nsec  , 1 )
   ierr = putvar (ncout, id_varout(6), rjsec(:,1),                            1,  nsec  , 1 )
   ierr = putvar (ncout, id_varout(7), e2usec(:,1),                           1,  nsec-1, 1 )
   ierr = putvar (ncout, id_varout(8), e1vsec(:,1),                           1,  nsec-1, 1 )
   ierr = putvar (ncout, id_varout(12),e2usec(:,1) + e1vsec(:,1),             1,  nsec-1, 1 )
   ierr = putvar (ncout, id_varout(16),batsec(:,1),                           1,  nsec-1, 1 )

   ! Temperature and salinity are interpolated on the respective U or V  point for better flux computation
   DO jt=1, npt
      dbarot = 0.d0    ! reset barotropic transport 
      DO jk=1,npk
         temper(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime = jt)
         saline(:,:) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime = jt)
         uzonal(:,:) = getvar(cf_ufil, cn_vozocrtx, jk, npiglo, npjglo, ktime = jt)
         vmerid(:,:) = getvar(cf_vfil, cn_vomecrty, jk, npiglo, npjglo, ktime = jt)

         e3u(:,:)    = getvar(cn_fzgr, 'e3u_ps',    jk, npiglo, npjglo, ldiom=.true.)
         e3v(:,:)    = getvar(cn_fzgr, 'e3v_ps',    jk, npiglo, npjglo, ldiom=.true.)

         DO jipt=1,nsec-1
            ii  = iisec(jipt  ) ; ij  = ijsec(jipt  )  ! F point  position
            ii1 = iisec(jipt+1) ; ij1 = ijsec(jipt+1)  ! Next F point position
            IF ( ij1  == ij ) THEN ! horizontal segment
               uzonalsec(jipt,jk) = 0.
               e3usec   (jipt,jk) = 0.
               IF ( ii1 > ii ) THEN ! eastward

                  IF ( MIN( saline(ii+1,ij) , saline(ii+1,ij+1))  == 0. ) THEN
                     tempersec(jipt,jk) = 0. ; salinesec(jipt,jk) = 0.
                  ELSE
                     tempersec(jipt,jk) = 0.5 * ( temper(ii+1,ij) + temper(ii+1,ij+1) )
                     salinesec(jipt,jk) = 0.5 * ( saline(ii+1,ij) + saline(ii+1,ij+1) )
                  ENDIF
                  vmeridsec(jipt,jk) = vmerid(ii+1,ij) * normv_sec(jipt)
                  e3vsec   (jipt,jk) = e3v   (ii+1,ij)

               ELSE ! westward

                  IF ( MIN( saline(ii,ij) , saline(ii,ij+1) ) == 0. ) THEN
                     tempersec(jipt,jk) = 0. ; salinesec(jipt,jk) = 0.
                  ELSE
                     tempersec(jipt,jk) = 0.5 * ( temper(ii,ij) + temper(ii,ij+1) )
                     salinesec(jipt,jk) = 0.5 * ( saline(ii,ij) + saline(ii,ij+1) )
                  ENDIF
                  vmeridsec(jipt,jk) = vmerid(ii,ij) * normv_sec(jipt)
                  e3vsec   (jipt,jk) = e3v   (ii,ij)

               ENDIF
            ELSEIF ( ii1 == ii ) THEN ! vertical segment
               vmeridsec(jipt,jk) = 0.
               e3vsec   (jipt,jk) = 0.
               IF ( ij1 < ij ) THEN ! southward

                  IF ( MIN( saline(ii,ij) , saline(ii+1,ij) ) == 0. ) THEN
                     tempersec(jipt,jk) = 0. ; salinesec(jipt,jk) = 0.
                  ELSE
                     tempersec(jipt,jk) = 0.5 * ( temper(ii,ij) + temper(ii+1,ij) )
                     salinesec(jipt,jk) = 0.5 * ( saline(ii,ij) + saline(ii+1,ij) )
                  ENDIF
                  uzonalsec(jipt,jk) = uzonal(ii,ij) * normu_sec(jipt)
                  e3usec   (jipt,jk) = e3u   (ii,ij)

               ELSE ! northward

                  IF ( MIN( saline(ii,ij+1) , saline(ii+1,ij+1) ) == 0. ) THEN
                     tempersec(jipt,jk) = 0. ; salinesec(jipt,jk) = 0.
                  ELSE
                     tempersec(jipt,jk) = 0.5 * ( temper(ii,ij+1) + temper(ii+1,ij+1) )
                     salinesec(jipt,jk) = 0.5 * ( saline(ii,ij+1) + saline(ii+1,ij+1) )
                  ENDIF
                  uzonalsec(jipt,jk) = uzonal(ii,ij+1) * normu_sec(jipt)
                  e3usec   (jipt,jk) = e3u   (ii,ij+1)

               ENDIF

            ELSE
               PRINT *, 'problem 2 for JIPT = ', jipt, 'JK=', jk
               PRINT *, '             I(P2)=',ii1, 'J(P1)=', ii
               PRINT *, '             J(P2)=',ij1, 'J(P1)=', ij
               exit 
            ENDIF

            dtmp=1.d0* (uzonalsec(jipt,jk) + vmeridsec(jipt,jk))*    &
                 (e2usec(jipt,1 )+ e1vsec(jipt,1 ))*                 &
                 (e3usec(jipt,jk)+ e3vsec(jipt,jk))
            dbarot=dbarot+dtmp
         END DO

         ierr = putvar (ncout, id_varout(1), tempersec(:,jk), jk, nsec-1, 1, ktime=jt )
         ierr = putvar (ncout, id_varout(2), salinesec(:,jk), jk, nsec-1, 1, ktime=jt )
         ierr = putvar (ncout, id_varout(3), uzonalsec(:,jk), jk, nsec-1, 1, ktime=jt )
         ierr = putvar (ncout, id_varout(4), vmeridsec(:,jk), jk, nsec-1, 1, ktime=jt )
         ! along-track normal velocity, horiz. and vert. resolution, and mask
         ierr = putvar (ncout, id_varout(11),uzonalsec(:,jk) + vmeridsec(:,jk), &
           &                                                  jk, nsec-1, 1, ktime=jt ) 

         IF ( jt == 1 ) THEN 
            ! save a mask of the section
            vmasksec(:,:) = 1.
            WHERE( salinesec(:,:) == 0. ) vmasksec(:,:) = 0.

            ierr = putvar (ncout, id_varout(9), e3usec(:,jk),                jk, nsec-1, 1 )
            ierr = putvar (ncout, id_varout(10),e3vsec(:,jk),                jk, nsec-1, 1 )
            ierr = putvar (ncout, id_varout(13),e3usec  (:,jk)+e3vsec(:,jk), jk, nsec-1, 1 )
            ierr = putvar (ncout, id_varout(14),vmasksec(:,jk),              jk, nsec-1, 1 )
         ENDIF
      END DO
      PRINT *, 'BAROTROPIC TRANSPORT at time ',jt,' = ', dbarot/1.d6, ' Sv.'
      ierr  = putvar0d ( ncout, id_varout(15), REAL(dbarot/1.d6), ktime = jt    )

   END DO

   ierr = closeout(ncout)

CONTAINS 
   SUBROUTINE CreateOutputFile
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE CreateOutputFile  ***
      !!
      !! ** Purpose :  Perform output file creation with all the variables 
      !!
      !! ** Method  :  Move this part of the code in a subroutine for clarity
      !!               All variables are global.  
      !!
      !!----------------------------------------------------------------------

      stypvar%scale_factor= 1.
      stypvar%add_offset= 0.
      stypvar%savelog10= 0.
      stypvar%rmissing_value=0.
      stypvar%conline_operation='N/A'

      ! define new variables for output 
      stypvar(1)%cname       = cn_votemper
      stypvar(1)%cunits      = 'deg C'
      stypvar(1)%valid_min   = -2.
      stypvar(1)%valid_max   = 40.
      stypvar(1)%clong_name  = 'Temperature along '//TRIM(csection)//' section'
      stypvar(1)%cshort_name = cn_votemper
      stypvar(1)%caxis       = 'TZX'
      ipk(1)                 = npk

      stypvar(2)%cname       = cn_vosaline
      stypvar(2)%cunits      = 'PSU'
      stypvar(2)%valid_min   = 0.
      stypvar(2)%valid_max   = 50.
      stypvar(2)%clong_name  = 'Salinity along '//TRIM(csection)//' section'
      stypvar(2)%cshort_name = cn_vosaline
      stypvar(2)%caxis       = 'TZX'
      ipk(2)                 = npk

      stypvar(3)%cname       = TRIM(cn_vozocrtx)//'_native'
      stypvar(3)%cunits      = 'm.s-1'
      stypvar(3)%valid_min   = -20.
      stypvar(3)%valid_max   = 20.
      stypvar(3)%clong_name  = 'Zonal velocity along '//TRIM(csection)//' section'
      stypvar(3)%cshort_name = TRIM(cn_vozocrtx)//'_native'
      stypvar(3)%caxis       = 'TZX'
      ipk(3)                 = npk

      stypvar(4)%cname       = TRIM(cn_vomecrty)//'_native'
      stypvar(4)%cunits      = 'm.s-1'
      stypvar(4)%valid_min   = -20.
      stypvar(4)%valid_max   = 20.
      stypvar(4)%clong_name  = 'Meridionnal velocity along '//TRIM(csection)//' section'
      stypvar(4)%cshort_name = TRIM(cn_vomecrty)//'_native'
      stypvar(4)%caxis       = 'TZX'
      ipk(4)                 = npk

      stypvar(5)%cname       = 'isec'
      stypvar(5)%valid_min   = 1.
      stypvar(5)%valid_max   = npiglo 
      stypvar(5)%caxis       = 'TX'
      ipk(5)                 = 1

      stypvar(6)%cname       = 'jsec'
      stypvar(6)%valid_min   = 1.
      stypvar(6)%valid_max   = npjglo 
      stypvar(6)%caxis       = 'TX'
      ipk(6)                 = 1

      stypvar(7)%cname       = TRIM(cn_ve2u)//'_native'
      stypvar(7)%valid_min   = 1.
      stypvar(7)%valid_max   = 200000.
      stypvar(7)%caxis       = 'TX'
      ipk(7)                 = 1

      stypvar(8)%cname       = TRIM(cn_ve1v)//'_native'
      stypvar(8)%valid_min   = 1.
      stypvar(8)%valid_max   = 200000.
      stypvar(8)%caxis       = 'TX'
      ipk(8)                 = 1

      stypvar(9)%cname       = 'e3u_native'
      stypvar(9)%valid_min   = 1.
      stypvar(9)%valid_max   = 200000.
      stypvar(9)%caxis       = 'TZX'
      ipk(9)                 =  npk

      stypvar(10)%cname      = 'e3v_native'
      stypvar(10)%valid_min   = 1.
      stypvar(10)%valid_max   = 200000.
      stypvar(10)%caxis      = 'TZX'
      ipk(10)                =  npk

      stypvar(11)%cname       = cn_vomecrty
      stypvar(11)%cunits      = 'm.s-1'
      stypvar(11)%valid_min   = -20.
      stypvar(11)%valid_max   = 20.
      stypvar(11)%clong_name  = 'Normal velocity along '//TRIM(csection)//' section'
      stypvar(11)%cshort_name = cn_vomecrty
      stypvar(11)%caxis       = 'TZX'
      ipk(11)                 =  npk

      stypvar(12)%cname       = cn_ve1v
      stypvar(12)%cunits      = 'm'
      stypvar(12)%valid_min   = 0.
      stypvar(12)%valid_max   = 1000000.
      stypvar(12)%clong_name  = 'Local horiz. resolution along '//TRIM(csection)//' section'
      stypvar(12)%cshort_name = cn_ve1v
      stypvar(12)%caxis       = 'TX'
      ipk(12)                 = 1

      stypvar(13)%cname       = 'e3v_ps'
      stypvar(13)%cunits      = 'm'
      stypvar(13)%valid_min   = 0.
      stypvar(13)%valid_max   = 100000000.
      stypvar(13)%clong_name  = 'Local vert. resolution along '//TRIM(csection)//' section'
      stypvar(13)%cshort_name = 'e3v_ps'
      stypvar(13)%caxis       = 'TZX'
      ipk(13)                 =  npk

      stypvar(14)%cname       = 'vmask'
      stypvar(14)%cunits      ='1/0'
      stypvar(14)%valid_min   = 0.
      stypvar(14)%valid_max   = 1.
      stypvar(14)%rmissing_value = 9999.
      stypvar(14)%clong_name  ='Mask along '//TRIM(csection)//' section'
      stypvar(14)%cshort_name = 'vmask'
      stypvar(14)%caxis       = 'TZX'
      ipk(14)                 =  npk

      stypvar(15)%cname       = 'barotrop_'//TRIM(csection)
      stypvar(15)%cunits      ='Sv'
      stypvar(15)%valid_min   = -500.
      stypvar(15)%valid_max   = 500.
      stypvar(15)%rmissing_value = -99999.
      stypvar(15)%clong_name  = 'Barotropic_transport for '//TRIM(csection)//' section'
      stypvar(15)%cshort_name = 'barotrop_'//TRIM(csection)
      stypvar(15)%caxis       = 'T'
      ipk(15)                 =  -1

      stypvar(16)%cname       = 'Bathymetry'
      stypvar(16)%cunits      = 'm'
      stypvar(16)%valid_min   = 0.
      stypvar(16)%valid_max   = 1000000.
      stypvar(16)%clong_name  = 'Bathymetry along '//TRIM(csection)//' section'
      stypvar(16)%cshort_name = 'Bathymetry'
      stypvar(16)%caxis       = 'TX'
      ipk(16)                 = 1


      ! create output fileset
      ncout = create      (cf_out, cf_tfil, nsec,  1, npk, cdep='deptht'                   )
      ierr  = createvar   (ncout,  stypvar, nvar,  ipk, id_varout                          )
      ierr  = putheadervar(ncout,  cf_tfil, nsec-1,  1, npk, pnavlon=rlonsec, pnavlat=rlatsec )
      tim   = getvar1d    (cf_tfil, 'time_counter', npt                                    )
      ierr  = putvar1d    (ncout, tim, npt, 'T'                                            )

   END SUBROUTINE CreateOutputFile



END PROGRAM cdf_xtract_brokenline
