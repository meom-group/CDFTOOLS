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
!============================
! WORK IN PROGRESS DO NOT USE 
!============================
   USE cdfio
   USE cdftools
   USE modcdfnames
   !!----------------------------------------------------------------------
   !! CDFTOOLS_3.0 , MEOM 2011
   !! $Id: cdfovide.f90 539 2011-07-11 10:33:35Z molines $
   !! Copyright (c) 2011, J.-M. Molines
   !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
   !!----------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER(KIND=4) :: jleg,  jk,  jipt, jvar         ! dummy loop index
   INTEGER(KIND=4) :: narg, iargc, ijarg, ifree      ! command line
   INTEGER(KIND=4) :: numin=10                       ! logical unit for input section file
   INTEGER(KIND=4) :: npiglo, npjglo, npk, npt       ! size of the domain
   INTEGER(KIND=4) :: iimin, iimax, ijmin, ijmax     ! ending points of a leg in model I J
   INTEGER(KIND=4) :: ii, ij, ipoint                 ! working integer
   INTEGER(KIND=4) :: ierr, ncout                    ! Netcdf error and ncid
   INTEGER(KIND=4) :: nvar = 14                    ! number of output variables
   INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout  ! netcdf output stuff

   ! broken line definition
   INTEGER(KIND=4) :: nsta=5                         ! number of points defining the broken line
   INTEGER(KIND=4) :: nsec=0                         ! total number of points on the broken line
   INTEGER(KIND=4) :: nn                             ! working integer (number of points in a leg)
   INTEGER(KIND=4) :: norm_u, norm_v                 ! velocity normalization (further use)
   INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: iista, ijsta     ! I,J position of the point on the broken line
   INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: ikeepn           ! Number of points per leg
   INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: iisec, ijsec     ! F-index of points on the broken line
   INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: iilegs, ijlegs   ! F-index of points on the broken line per leg

   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: rlonsta, rlatsta    ! Geographic position defining legs
   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: rxx, ryy            ! leg i j index of F points
   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                 ! Model time array

   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1v, e3v            ! V point relevant metric
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2u, e3u            ! U point relevant metric
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlon, rlat          ! model long and lat of T points
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: temper, saline      ! model Temperature and salinity
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: uzonal, vmerid      ! model zonal and meridional velocity
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tempersec, salinesec, uzonalsec, vmeridsec
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlonsec, rlatsec, risec, rjsec
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1vsec, e2usec
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e3usec, e3vsec
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: vmasksec
   REAL(KIND=4)                              :: xmin, xmax, ymin, ymax !

   REAL(KIND=8)                              :: dtmp, dbarot  ! for barotropic transport computation

   CHARACTER(LEN=80) :: cf_tfil , cf_ufil, cf_vfil   ! input T U V files
   CHARACTER(LEN=80) :: cf_out                       ! output file
   CHARACTER(LEN=80) :: cf_sec                       ! input section file
   CHARACTER(LEN=80) :: csection                     ! section name
   CHARACTER(LEN=80) :: cverb='n'                    ! verbose key for findij
   CHARACTER(LEN=80) :: cldum                        ! dummy character variable

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
      PRINT *,'       ' 
      PRINT *,'      '
      PRINT *,'     ARGUMENTS :'
      PRINT *,'       T-file :  model gridT file '
      PRINT *,'       U-file :  model gridU file '
      PRINT *,'       V-file :  model gridV file '
      PRINT *,'      ' 
      PRINT *,'     OPTIONS :'
      PRINT *,'       -f section_file : provide a file for section definition.'
      PRINT *,'              section_file is an ascii file as follows:'
      PRINT *,'              * line #1 : name of the section (e.g. ovide). '
      PRINT *,'                   Will be used for naming the output file.'
      PRINT *,'              * line #2 : number of points defining the broken line.'
      PRINT *,'              * line #3-end : a pair of Longitude latitude values defining'
      PRINT *,'                    the points. If not supplied, use hard-coded information'
      PRINT *,'                    for OVIDE section. A comment can be added at the end of'
      PRINT *,'                    of the lines, using a # as separator'
      PRINT *,'       -verbose : increase verbosity  ' 
      PRINT *,'      '
      PRINT *,'     REQUIRED FILES :'
      PRINT *,'       ', TRIM(cn_fhgr),' and ',TRIM(cn_fzgr),' must be in the current directory ' 
      PRINT *,'      '
      PRINT *,'     OUTPUT : '
      PRINT *,'       netcdf file : section_name.nc'
!     PRINT *,'         variables : ', TRIM(cv_out),' (    )'
      PRINT *,'      '
      PRINT *,'     SEE ALSO :'
      PRINT *,'      ' 
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
      CALL broken_line( iimin, iimax, ijmin, ijmax, rxx, ryy, nn, npiglo, npjglo, norm_u, norm_v )

      IF ( lverbose) PRINT *, 'Leg ', jleg,' : npoints : ', nn

      IF (rxx(1) < rxx(nn) ) THEN ! leg is oriented eastward
         iilegs(jleg,1:nn)=rxx(1:nn)
         ijlegs(jleg,1:nn)=ryy(1:nn)
      ELSE                        ! leg is oriented westward
         iilegs(jleg,1:nn)=rxx(nn:1:-1)
         ijlegs(jleg,1:nn)=ryy(nn:1:-1)
      END IF

      IF ( lverbose) THEN
        PRINT *, 'Leg  rxx   ryy '
        DO jk = 1, nn
           PRINT *, jleg, iilegs(jleg,jk), ijlegs(jleg,jk) ,rxx(jk), ryy(jk)
        END DO
      ENDIF

      ikeepn(jleg) = nn  ! number of points (F) on leg jleg
      nsec = nsec + nn   ! total number of points (F) on the broken line
   END DO !! loop on the legs

   ! fancy control print
   DO jleg = 1, nsta -1
     WRITE(*,*) '------------------------------------------------------------'
     WRITE(*,9100) 'leg ',jleg,'  start at ', rlonsta(jleg) ,'N ', rlatsta(jleg), 'W and ends at ', rlonsta(jleg+1) ,'N ', rlatsta(jleg+1), 'W'
     WRITE(*,9101) 'corresponding to F-gridpoints(', iista(jleg),',',ijsta(jleg),') and (', iista(jleg+1),',',ijsta(jleg+1),')' 
     WRITE(*,*) '------------------------------------------------------------'
   ENDDO

9100 FORMAT(a,i3,a,f6.2,a,f6.2,a,f6.2,a,f6.2,a)
9101 FORMAT(a,i4,a,i4,a,i4,a,i4,a)

   ! 3. : Extraction along the legs
   ! ------------------------------
   ALLOCATE (iisec(nsec), ijsec(nsec)) 

   ipoint = 0
   DO jleg=1, nsta-1  ! loop on legs 
      DO jipt=1, ikeepn(jleg) 
         ipoint = ipoint + 1
         iisec(ipoint)=iilegs(jleg,jipt)  ! i-index
         ijsec(ipoint)=ijlegs(jleg,jipt)  ! j-index
      END DO
   END DO


   ! input fields
   ALLOCATE(rlon(npiglo,npjglo), rlat(npiglo,npjglo))
   ALLOCATE(temper(npiglo,npjglo), saline(npiglo,npjglo))
   ALLOCATE(uzonal(npiglo,npjglo), vmerid(npiglo,npjglo))
   ALLOCATE(e1v(npiglo,npjglo))
   ALLOCATE(e2u(npiglo,npjglo))
   ALLOCATE(e3u(npiglo,npjglo), e3v(npiglo,npjglo))

   ! output fields
   ALLOCATE(rlonsec(1,nsec), rlatsec(1,nsec) )
   ALLOCATE(risec(1,nsec), rjsec(1,nsec) )
   ALLOCATE(e2usec(1,nsec-1), e3usec(nsec-1,npk) )
   ALLOCATE(e1vsec(1,nsec-1), e3vsec(nsec-1,npk) )
   ALLOCATE(tempersec(nsec-1,npk), salinesec(nsec-1,npk) )
   ALLOCATE(uzonalsec(nsec-1,npk), vmeridsec(nsec-1,npk) )
   ALLOCATE(vmasksec(nsec-1,npk))

   e1vsec = -9999.
   e2usec = -9999.

   risec(:,:) = 0.
   rjsec(:,:) = 0.

   rlon(:,:) = getvar(cn_fhgr, cn_glamt, 1, npiglo, npjglo)
   rlat(:,:) = getvar(cn_fhgr, cn_gphit, 1, npiglo, npjglo)
   e1v(:,:)  = getvar(cn_fhgr, cn_ve1v,  1, npiglo, npjglo)
   e2u(:,:)  = getvar(cn_fhgr, cn_ve2u,  1, npiglo, npjglo)

   ! il faut faire un test sur la continuité des segments
   ! on va prendre T et S comme etant la moyenne du point
   ! en dessous et au-dessus du segment pour pouvoir calculer
   ! les fluxs de maniere optimales...

   ! loop on 2d arrays
   DO jipt = 1,nsec
      ii = iisec(jipt)
      ij = ijsec(jipt)

      risec  (1,jipt) = ii
      rjsec  (1,jipt) = ij
      rlonsec(1,jipt) = rlon(ii, ij)
      rlatsec(1,jipt) = rlat(ii, ij)
   END DO

   jk=1

   DO jipt=1,nsec-1
      !PRINT*, 'jipt=', jipt

      IF ( ijsec(jipt+1) == ijsec(jipt) ) THEN ! horizontal segment
         e2usec(jk,jipt) = 0.
         IF ( iisec(jipt+1) > iisec(jipt) ) THEN ! eastward
            e1vsec(jk,jipt) = e1v(iisec(jipt)+1,ijsec(jipt))
         ELSE
            e1vsec(jk,jipt) = e1v(iisec(jipt),ijsec(jipt))
         ENDIF

      ELSEIF ( iisec(jipt+1) == iisec(jipt) ) THEN ! vertical segment
         e1vsec(jk,jipt) = 0.
         IF ( ijsec(jipt+1) < ijsec(jipt) ) THEN ! southward
            e2usec(jk,jipt) = e2u(iisec(jipt),ijsec(jipt))
         ELSE
            e2usec(jk,jipt) = e2u(iisec(jipt),ijsec(jipt)+1)
         ENDIF

      ELSE
         PRINT *, 'problem'
         exit 
      ENDIF
   END DO

   ! loop on 3d arrays
   DO jk=1,npk
      temper(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo)
      saline(:,:) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo)
      uzonal(:,:) = getvar(cf_ufil, cn_vozocrtx, jk, npiglo, npjglo)
      vmerid(:,:) = getvar(cf_vfil, cn_vomecrty, jk, npiglo, npjglo)
      e3u(:,:)    = getvar(cn_fzgr, 'e3u_ps',    jk, npiglo, npjglo, ldiom=.true.)
      e3v(:,:)    = getvar(cn_fzgr, 'e3v_ps',    jk, npiglo, npjglo, ldiom=.true.)

      DO jipt=1,nsec-1
         IF ( ijsec(jipt+1) == ijsec(jipt) ) THEN ! horizontal segment
            uzonalsec(jipt,jk) = 0.
            e3usec(jipt,jk) = 0.
            IF ( iisec(jipt+1) > iisec(jipt) ) THEN ! eastward

               IF ( MIN( saline(iisec(jipt)+1,ijsec(jipt)) , saline(iisec(jipt)+1,ijsec(jipt)+1) ) == 0. ) THEN
                  tempersec(jipt,jk) = 0. ; salinesec(jipt,jk) = 0.
               ELSE
                  tempersec(jipt,jk) = 0.5 * ( temper(iisec(jipt)+1,ijsec(jipt)) + temper(iisec(jipt)+1,ijsec(jipt)+1) )
                  salinesec(jipt,jk) = 0.5 * ( saline(iisec(jipt)+1,ijsec(jipt)) + saline(iisec(jipt)+1,ijsec(jipt)+1) )
               ENDIF
               vmeridsec(jipt,jk) = vmerid(iisec(jipt)+1,ijsec(jipt))
               e3vsec(jipt,jk) = e3v(iisec(jipt)+1,ijsec(jipt))

            ELSE ! westward

               IF ( MIN( saline(iisec(jipt),ijsec(jipt)) , saline(iisec(jipt),ijsec(jipt)+1) ) == 0. ) THEN
                  tempersec(jipt,jk) = 0. ; salinesec(jipt,jk) = 0.
               ELSE
                  tempersec(jipt,jk) = 0.5 * ( temper(iisec(jipt),ijsec(jipt)) + temper(iisec(jipt),ijsec(jipt)+1) )
                  salinesec(jipt,jk) = 0.5 * ( saline(iisec(jipt),ijsec(jipt)) + saline(iisec(jipt),ijsec(jipt)+1) )
               ENDIF
               vmeridsec(jipt,jk) = vmerid(iisec(jipt),ijsec(jipt))
               e3vsec(jipt,jk) = e3v(iisec(jipt),ijsec(jipt))

            ENDIF
         ELSEIF ( iisec(jipt+1) == iisec(jipt) ) THEN ! vertical segment
            vmeridsec(jipt,jk) = 0.
            e3vsec(jipt,jk) = 0.
            IF ( ijsec(jipt+1) < ijsec(jipt) ) THEN ! southward

               IF ( min( saline(iisec(jipt),ijsec(jipt)) , saline(iisec(jipt)+1,ijsec(jipt)) ) == 0. ) THEN
                  tempersec(jipt,jk) = 0. ; salinesec(jipt,jk) = 0.
               ELSE
                  tempersec(jipt,jk) = 0.5 * ( temper(iisec(jipt),ijsec(jipt)) + temper(iisec(jipt)+1,ijsec(jipt)) )
                  salinesec(jipt,jk) = 0.5 * ( saline(iisec(jipt),ijsec(jipt)) + saline(iisec(jipt)+1,ijsec(jipt)) )
               ENDIF
               uzonalsec(jipt,jk) = uzonal(iisec(jipt),ijsec(jipt))
               e3usec(jipt,jk) = e3u(iisec(jipt),ijsec(jipt))

            ELSE ! northward

               IF ( min( saline(iisec(jipt),ijsec(jipt)+1) , saline(iisec(jipt)+1,ijsec(jipt)+1) ) == 0. ) THEN
                  tempersec(jipt,jk) = 0. ; salinesec(jipt,jk) = 0.
               ELSE
                  tempersec(jipt,jk) = 0.5 * ( temper(iisec(jipt),ijsec(jipt)+1) + temper(iisec(jipt)+1,ijsec(jipt)+1) )
                  salinesec(jipt,jk) = 0.5 * ( saline(iisec(jipt),ijsec(jipt)+1) + saline(iisec(jipt)+1,ijsec(jipt)+1) )
               ENDIF
               uzonalsec(jipt,jk) = uzonal(iisec(jipt),ijsec(jipt)+1)
               e3usec(jipt,jk) = e3u(iisec(jipt),ijsec(jipt)+1)

            ENDIF

         ELSE
            PRINT *, 'problem'
            exit 
         ENDIF

      END DO
   END DO

   ! save a mask of the section
   vmasksec(:,:) = 1.
   WHERE( salinesec(:,:) == 0. ) vmasksec(:,:) = 0.


   ! 4. Compute barotropic transport across the section for validation
   ! -----------------------------------------------------------------
   dbarot = 0.d0

   DO jipt=1,nsec-1
      DO jk=1,npk
         dtmp=1.d0* (uzonalsec(jipt,jk) + vmeridsec(jipt,jk))*&
              (e2usec(1,jipt )+ e1vsec(1,jipt ))*                 &
              (e3usec(jipt,jk)+ e3vsec(jipt,jk))*                 &
              vmasksec(jipt,jk)
         dbarot=dbarot+dtmp
      ENDDO
   ENDDO

   PRINT*, 'BAROTROPIC TRANSPORT = ', dbarot/1.e6, ' Sv.'

   ! 5. Output of the section in netcdf file
   ! ----------------------------------------
   ALLOCATE ( stypvar(nvar), ipk(nvar), id_varout(nvar) )

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
   stypvar(7)%valid_min   = MINVAL(e2usec(1,:))
   stypvar(7)%valid_max   = MAXVAL(e2usec(1,:)) 
   stypvar(7)%caxis       = 'TX'
   ipk(7)                 = 1

   stypvar(8)%cname       = TRIM(cn_ve1v)//'_native'
   stypvar(8)%valid_min   = MINVAL(e1vsec(1,:))
   stypvar(8)%valid_max   = MAXVAL(e1vsec(1,:))
   stypvar(8)%caxis       = 'TX'
   ipk(8)                 = 1

   stypvar(9)%cname       = 'e3u_native'
   stypvar(9)%valid_min   = MINVAL(e3usec(:,:))
   stypvar(9)%valid_max   = MAXVAL(e3usec(:,:)) 
   stypvar(9)%caxis       = 'TZX'
   ipk(9)                 =  npk

   stypvar(10)%cname      = 'e3v_native'
   stypvar(10)%valid_min  = MINVAL(e3vsec(:,:))
   stypvar(10)%valid_max  = MAXVAL(e3vsec(:,:)) 
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

   ! create output fileset
   ncout = create      (cf_out, cf_tfil, nsec,  1, npk, cdep='deptht'                   )
   ierr  = createvar   (ncout,  stypvar, nvar,  ipk, id_varout                          )
   ierr  = putheadervar(ncout,  cf_tfil, nsec,  1, npk, pnavlon=rlonsec, pnavlat=rlatsec )
   tim   = getvar1d    (cf_tfil, 'time_counter', npt                                    )
   ierr  = putvar1d    (ncout, tim, npt, 'T'                                            )



   ! netcdf output 
   DO jk =1, npk
      ierr = putvar (ncout, id_varout(1), tempersec(:,jk), jk, nsec-1, 1 )
      ierr = putvar (ncout, id_varout(2), salinesec(:,jk), jk, nsec-1, 1 )
      ierr = putvar (ncout, id_varout(3), uzonalsec(:,jk), jk, nsec-1, 1 )
      ierr = putvar (ncout, id_varout(4), vmeridsec(:,jk), jk, nsec-1, 1 )
      ierr = putvar (ncout, id_varout(9), e3usec(:,jk),      jk, nsec-1, 1 )
      ierr = putvar (ncout, id_varout(10),e3vsec(:,jk),      jk, nsec-1, 1 )

      ! along-track normal velocity, horiz. and vert. resolution, and mask
      ierr = putvar (ncout, id_varout(11),uzonalsec(:,jk) + vmeridsec(:,jk),     jk, nsec-1, 1 )
      ierr = putvar (ncout, id_varout(13),e3usec(:,jk) + e3vsec(:,jk),           jk, nsec-1, 1 )
      ierr = putvar (ncout, id_varout(14),vmasksec(:,jk),                        jk, nsec-1, 1 )
   END DO
      ierr = putvar (ncout, id_varout(5), risec(1,:),                            1,  nsec  , 1 )
      ierr = putvar (ncout, id_varout(6), rjsec(1,:),                            1,  nsec  , 1 )
      ierr = putvar (ncout, id_varout(7), e2usec(1,:),                           1,  nsec-1, 1 )
      ierr = putvar (ncout, id_varout(8), e1vsec(1,:),                           1,  nsec-1, 1 )
      ierr = putvar (ncout, id_varout(12),e2usec(1,:) + e1vsec(1,:),             1,  nsec-1, 1 )

   ierr = closeout(ncout)

END PROGRAM cdf_xtract_brokenline
