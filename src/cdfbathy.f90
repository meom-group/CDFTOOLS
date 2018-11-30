PROGRAM cdfbathy
  !!======================================================================
  !!                     ***  PROGRAM  cdfbathy  ***
  !!=====================================================================
  !!  ** Purpose : Utility to modify a bathymetric file according to 
  !!               specific option (eg : fill an area, modify points ...)
  !!               Using -var option and -lev, can also edit any file with
  !!               the same tool, except the specific actions dedicated to
  !!               the bathymetry (eg : zstep like ...)
  !!
  !!  ** Method  : All modifications are save in a fortran file ready to be
  !!               used to replay all the modif at once.
  !!
  !! History : 2.1  : 11/2007  : J.M. Molines : Original code
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!                : 04/1014  : P. Mathiot   : add fill_pool option
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   zgr_zps 
  !!   zgr_read
  !!   prlog
  !!   fillzone
  !!   raz_zone
  !!   raz_below
  !!   set_below
  !!   dumpzone
  !!   nicedumpzone
  !!   replacezone
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class data_transformation
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  !
  INTEGER(KIND=4)                              :: jk, jt              ! loop index
  INTEGER(KIND=4)                              :: narg, iargc, ijarg  ! browse command line
  INTEGER(KIND=4)                              :: iimin, iimax        ! selected area
  INTEGER(KIND=4)                              :: ijmin, ijmax        ! selected area
  INTEGER(KIND=4)                              :: ierr                ! error status
  INTEGER(KIND=4)                              :: icrit               ! maximal size of pool 
  INTEGER(KIND=4)                              :: iklev               ! selected level
  INTEGER(KIND=4)                              :: itime               ! selected time
  INTEGER(KIND=4)                              :: npiglo, npjglo      ! domain size
  INTEGER(KIND=4)                              :: npk, npt            ! domaine size
  INTEGER(KIND=4)                              :: nnpk, nnpt            ! domaine size
  INTEGER(KIND=4)                              :: iversion=1          ! version counter for working copy
  INTEGER(KIND=4)                              :: iostat, ipos        ! used for version control
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mbathy              ! map of wet model level

  REAL(KIND=4)                                 :: e3zps_min=25.       ! minimum thickness of bottom cell
  REAL(KIND=4)                                 :: e3zps_rat=0.2       ! minimum ratio e3bot/e3_0
  REAL(KIND=4)                                 :: rdepmin=600.        ! default value for depmin (full step like)
  REAL(KIND=4)                                 :: rdepfill=0.         ! default filling value
  REAL(KIND=4)                                 :: scale_factor=1.     ! divide by scale factor when reading
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: e3t, e3w            ! vertical metrics
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: gdept, gdepw        ! depth at T and W points
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e3_bot              ! bottom depth (partial steps)
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: bathyin             ! initial data value
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: bathy               ! modified data value
  !
  CHARACTER(LEN=256)                           :: cf_in               ! original input file name
  CHARACTER(LEN=256)                           :: cf_root             ! root part of the file name
  CHARACTER(LEN=256)                           :: cf_dump             ! dump txt file name (out)
  CHARACTER(LEN=256)                           :: cf_replace          ! replace txt file name (in)
  CHARACTER(LEN=80)                            :: cf_batfile = 'zgrbat.txt' ! txt file giving vertical mesh
  CHARACTER(LEN=80)                            :: cf_log = 'log.f90'  ! default log file
  CHARACTER(LEN=80)                            :: cv_in               ! variable name
  CHARACTER(LEN=256)                           :: cwkc                ! filename of working copy
  CHARACTER(LEN=256)                           :: cldum               ! dummy string

  LOGICAL :: lexist     = .TRUE.,  lfill      = .FALSE.       ! all required flags for options
  LOGICAL :: lfullstep  = .FALSE., lappend    = .FALSE.       ! all required flags for options
  LOGICAL :: lreplace   = .FALSE., ldump      = .FALSE.       ! all required flags for options
  LOGICAL :: lmodif     = .FALSE., loverwrite = .FALSE.       ! all required flags for options
  LOGICAL :: lraz       = .FALSE., ldumpn     = .FALSE.       ! all required flags for options
  LOGICAL :: lrazb      = .FALSE., lsetb      = .FALSE.       ! all required flags for options
  LOGICAL :: lsetz      = .FALSE., lseta      = .FALSE.       ! all required flags for options
  LOGICAL :: lchk       = .FALSE., lfillpool  = .FALSE.       ! all required flags for options
  LOGICAL :: ll_lev     = .FALSE., ll_time    = .FALSE.       ! all level, all time ==> true
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfbathy/cdfvar -f IN-file [options]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Allow manual modification of the input file. Very convenient for ' 
     PRINT *,'       bathymetric files, can also be used with any other model file.'
     PRINT *,'       Keep a log.f90 file of the modifications for automatic reprocessing'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file : original input file. The program works on a copy of the'
     PRINT *,'                original file (default)' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT 9999, '   -file (or -f )       : name of input file '
     PRINT 9999, '   -var  (or -v )       : name of cdf variable [default: Bathymetry]'
     PRINT 9999, '   -lev  (or -l )       : level to work with (0 means all levels)'
     PRINT 9999, '   -time (or -t )       : time to work with (0 means all times)'
     PRINT 9999, '   -scale  s            : use s as a scale factor (divide when read the file)'
     PRINT 9999, '   -zoom (or -z )       : sub area of the bathy file to work with (imin imax jmin jmax)'
     PRINT 9999, '   -fillzone (or -fz )  : sub area will be filled with 0 up to the first coast line '
     PRINT 9999, '   -fillpool (or -fp ) [ icrit ] : the whole file is check and fill all the pool smaller than (icrit) cell by 0'
     PRINT 9999, '   -raz_zone (or -raz ) : sub area will be filled with 0 up '
     PRINT 9999, '   -raz_below depmin    : any depth less than depmin in subarea will be replaced by 0 '
     PRINT 9999, '      (or -rb depmin )  '
     PRINT 9999, '   -set_below depmin    : any depth less than depmin in subarea will be replaced by depmin '
     PRINT 9999, '      (or -sb depmin ) '
     PRINT 9999, '   -set_above depmax    : any depth larger than depmax in subarea will be replaced by depmax '
     PRINT 9999, '      (or -sa depmax ) '
     PRINT 9999, '   -set_zone value      : all value in area will be set to value'
     PRINT 9999, '      (or -sz value ) '
     PRINT 9999, '   -fullstep depmin     : sub area will be reshaped as full-step, below depmin'
     PRINT 9999, '      (or -fs depmin )    requires the presence of the file zgr_bat.txt (from ocean.output, eg )'
     PRINT 9999, '   -dumpzone (or -d )   : sub area will be output to an ascii file, which can be used by -replace'
     PRINT 9999, '                          after manual editing '
     PRINT 9999, '   -nicedumpzone        : sub area will be output to an ascii file (nice output)'
     PRINT 9999, '            (or -nd )'
     PRINT 9999, '   -replace (or -r )    : sub area defined by the file will replace the original bathy'
     PRINT 9999, '   -append (or -a )     : fortran log file (log.f90) will be append with actual modif'
     PRINT 9999, '                          Standard behaviour is to overwrite/create log file'
     PRINT 9999, '   -overwrite (or -o )  : input bathy file will be used as output.'
     PRINT 9999, '                          Standard behaviour is to use a work copy of the original file'
     PRINT 9999, '                          (indexed from 01 to 99 if necessary ) '
     PRINT 9999, '   -log logfile         : log file for change (default is log.f90) '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT 9999, '      netcdf file : according to used options, if the original file is to be modified'
     PRINT 9999, '             a sequence number is added at the end of the input file name, to keep'
     PRINT 9999, '             modifications.'
     PRINT *,'            variables : same as input file'
     STOP 
  ENDIF
9999 FORMAT(5x,a)

  ijarg = 1
  iimin=-10 ; iimax=-10 ; ijmin=-10 ; ijmax=-10

  cv_in = cn_bathymet  ! default value
  iklev = 1
  itime = 1

  DO  WHILE (ijarg <=  narg)
     CALL getarg(ijarg, cldum) ;  ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-file' , '-f') ; CALL getarg(ijarg, cf_in) ; ijarg = ijarg + 1 
        ; lchk = ( lchk .OR. chkfile (cf_in) )
     CASE ( '-var' , '-v' ) ; CALL getarg(ijarg, cv_in) ; ijarg = ijarg + 1 
     CASE ( '-lev' , '-k' ) ; CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) iklev
     CASE ( '-time' , '-t') ; CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) itime
     CASE ( '-scale'      ) ; CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) scale_factor
     CASE ( '-zoom' , '-z') ; CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) iimin
        ;                     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) iimax
        ;                     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmin
        ;                     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmax
     CASE ('-fillzone','-fz' ) ; lfill =.TRUE. ; lmodif =.TRUE.
     CASE ('-raz_zone','-raz') ; lraz  =.TRUE. ; lmodif =.TRUE.
     CASE ('-raz_below','-rb') ; CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rdepfill
        ; lrazb =.TRUE. ; lmodif =.TRUE.
     CASE ('-set_below','-sb') ; CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rdepfill
        ; lsetb =.TRUE. ; lmodif =.TRUE.
     CASE ('-set_above','-sa') ; CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rdepfill
        ; lseta =.TRUE. ; lmodif =.TRUE.
     CASE ('-set_zone','-sz') ; CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rdepfill
        ; lsetz =.TRUE. ; lmodif =.TRUE.
     CASE ('-fullstep','-fs' ) ; lfullstep =.TRUE. ; lmodif=.TRUE.
        ; CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rdepmin
     CASE ('-append' , '-a'  ) ; lappend=.TRUE.
     CASE ('-overwrite' ,'-o') ; loverwrite=.TRUE.
     CASE ('-fillpool','-fp' ) ; lfillpool =.TRUE. ; lmodif =.TRUE.
        ; CALL getarg(ijarg, cldum) ; ijarg = ijarg +1 ; READ(cldum,*) icrit
     CASE ('-replace','-r'   ) ; lreplace =.TRUE. ; lmodif =.TRUE.
        ; CALL getarg(ijarg, cf_replace) ; ijarg = ijarg +1
        ; lchk = ( lchk .OR. chkfile (cf_replace) )
     CASE ( '-log'           ) ; CALL getarg(ijarg, cf_log) ; ijarg = ijarg +1
     CASE ( '-dumpzone','-d' ) ; ldump =.TRUE.
        ; CALL getarg(ijarg, cf_dump) ; ijarg = ijarg +1
     CASE ('-nicedumpzone','-nd'); ldumpn =.TRUE.
        ; CALL getarg(ijarg, cf_dump) ; ijarg = ijarg +1
     CASE DEFAULT              ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  IF ( lchk ) STOP 99  ! missing files

  IF ( lmodif .AND. .NOT. loverwrite) THEN
     ! creating a working copy of the file indexed by iversion
     ipos=INDEX(cf_in,'.',.TRUE.)
     READ(cf_in(ipos+1:),*,IOSTAT=iostat) iversion
     IF (iostat /=0 ) THEN 
        iversion=0
        cf_root=cf_in
     ELSE
        cf_root=cf_in(1:ipos-1)
     ENDIF
     iversion=iversion+1

     DO WHILE ( lexist )
        WRITE(cwkc,'(a,a,i2.2)') TRIM(cf_root),'.',iversion
        INQUIRE(FILE=cwkc,EXIST=lexist)
        iversion=iversion+1
     END DO
     PRINT *, 'Working copy will be : ' ,TRIM(cwkc)
     CALL system(' cp -f '//cf_in//' '//cwkc )
  ELSE
     cwkc=cf_in
  ENDIF

  npiglo = getdim(cwkc,cn_x)
  npjglo = getdim(cwkc,cn_y)
  npk    = getdim(cwkc,cn_z)
  npt    = getdim(cwkc,cn_t)
  IF (npk == 0 ) npk = 1
  IF (npt == 0 ) npt = 1

  nnpk=1
  nnpt=1

  IF ( iklev == 0 )  THEN
     ll_lev=.true.
     nnpk = npk
  ENDIF


  IF ( iklev > npk ) THEN
     PRINT *,' ERROR : not enough levels in input file ', TRIM(cwkc)
  ENDIF

  IF ( itime > npt ) THEN
     PRINT *,' ERROR : not enough times in input file ', TRIM(cwkc)
  ENDIF

  IF ( itime == 0 )  THEN
     ll_time=.true.
     nnpt = npt
  ENDIF

  IF ( iimin == -10 ) THEN  ! no zoom option passed
     iimin=1 ; iimax=npiglo
     ijmin=1 ; ijmax=npjglo
  END IF

  PRINT *, 'NPIGLO = ', npiglo
  PRINT *, 'NPJGLO = ', npjglo
  PRINT *, 'IMIN IMAX JMIN JMAX :', iimin, iimax,ijmin,ijmax

  ALLOCATE (mbathy(npiglo,npjglo), e3_bot( npiglo,npjglo))
  ALLOCATE (bathy( npiglo,npjglo), bathyin(npiglo,npjglo))

  DO jt = 1,nnpt
     IF (ll_time ) itime=jt
     DO jk =1,nnpk
        IF ( ll_lev ) iklev = jk

        ! we use bathy as variable name but it can be any field from cf_in
        bathy(:,:) = getvar(cwkc, cv_in, iklev, npiglo, npjglo, ktime=itime)
        bathy(:,:) = bathy(:,:)/scale_factor
        bathyin    = bathy  ! save original 

        IF (lfullstep ) THEN  ;  CALL zgr_read ; CALL zgr_zps(iimin, iimax, ijmin, ijmax); ENDIF
        IF (lfill     )       CALL fillzone  (iimin, iimax, ijmin, ijmax)
        IF (lraz      )       CALL raz_zone  (iimin, iimax, ijmin, ijmax)
        IF (lrazb     )       CALL raz_below (iimin, iimax, ijmin, ijmax, rdepfill)
        IF (lsetb     )       CALL set_below (iimin, iimax, ijmin, ijmax, rdepfill)
        IF (lseta     )       CALL set_above (iimin, iimax, ijmin, ijmax, rdepfill)
        IF (lsetz     )       CALL set_zone (iimin, iimax, ijmin, ijmax, rdepfill)
        IF (ldump     )       CALL dumpzone     (cf_dump, iimin, iimax, ijmin, ijmax)
        IF (ldumpn    )       CALL nicedumpzone (cf_dump, iimin, iimax, ijmin, ijmax)
        IF (lreplace  )       CALL replacezone  (cf_replace)
        IF (lfillpool )       CALL fillpool  (icrit, iimin, iimax, ijmin, ijmax)

        IF (lmodif ) THEN   ! save log 
           CALL prlog(bathyin, bathy, npiglo, npjglo, lappend)
           ierr = putvar(cwkc, cv_in, iklev, iimax-iimin+1, ijmax-ijmin+1, kimin=iimin, kjmin=ijmin, &
                &            ptab=bathy(iimin:iimax,ijmin:ijmax)*scale_factor, ktime=itime)
        ENDIF
     ENDDO
  ENDDO

CONTAINS 

  SUBROUTINE zgr_zps ( kimin, kimax ,kjmin, kjmax )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE zgr_zps  ***
    !!
    !! ** Purpose :   Build the partial steps
    !!
    !! ** Method  :  Use NEMO routine
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4) ,INTENT(in) :: kimin, kimax, kjmin, kjmax
    !! * Local declarations
    INTEGER(KIND=4)  ::   ji, jj, jk        ! dummy loop indices
    INTEGER(KIND=4)  ::   ik, it            ! temporary integers
    INTEGER(KIND=4), PARAMETER :: wp=4      ! working precision is 4 in the CDFTOOLS

    REAL(wp)           :: ze3tp, ze3wp  ! Last ocean level thickness at T- and W-points
    REAL(wp)           :: zdepwp        ! Ajusted ocean depth to avoid too small e3t
    REAL(wp)           :: zdepth        !    "         "
    REAL(wp)           :: zmax, zmin    ! Maximum and minimum depth
    REAL(wp)           :: zdiff         ! temporary scalar
    !!----------------------------------------------------------------------
    ! Initialization of constant
    zmax = gdepw(npk) + e3t(npk)
    zmin = gdepw(4)

    ! initialize mbathy to the maximum ocean level available
    mbathy(kimin:kimax,kjmin:kjmax) = npk-1

    ! storage of land and island's number (zero and negative values) in mbathy
    WHERE (bathy(kimin:kimax,kjmin:kjmax) <= 0. ) mbathy(kimin:kimax,kjmin:kjmax)=INT( bathy(kimin:kimax,kjmin:kjmax) )

    ! bounded value of bathy
    ! minimum depth == 3 levels
    ! maximum depth == gdepw(jpk)+e3t(jpk) 
    ! i.e. the last ocean level thickness cannot exceed e3t(jpkm1)+e3t(jpk)
    WHERE (bathy(kimin:kimax,kjmin:kjmax) <= 0 ) 
       bathy(kimin:kimax,kjmin:kjmax)=0.
    ELSEWHERE (bathy(kimin:kimax,kjmin:kjmax) < zmin ) 
       bathy(kimin:kimax,kjmin:kjmax) = zmin
    ELSEWHERE (bathy(kimin:kimax,kjmin:kjmax) >= zmax )
       bathy(kimin:kimax,kjmin:kjmax) = zmax
    END WHERE

    ! Compute mbathy for ocean points (i.e. the number of ocean levels)
    ! find the number of ocean levels such that the last level thickness
    ! is larger than the minimum of e3zps_min and e3zps_rat * e3t (where
    ! e3t is the reference level thickness
    DO jk = npk-1, 1, -1
       !      zdepth = gdepw(jk) + MIN( e3zps_min, e3t(jk)*e3zps_rat )
       zdepth = gdept(jk)
       WHERE ( bathy(kimin:kimax,kjmin:kjmax) > 0. .AND. bathy (kimin:kimax,kjmin:kjmax) <= zdepth )  
          mbathy(kimin:kimax,kjmin:kjmax)=jk-1
          e3_bot(kimin:kimax,kjmin:kjmax)= bathy(kimin:kimax,kjmin:kjmax) - gdepw(jk-1)
       END WHERE
    END DO

    DO ji=kimin,kimax
       DO jj=kjmin,kjmax
          jk=mbathy(ji,jj)
          IF (jk /= 0 ) THEN 
             IF (gdepw(jk+1) > rdepmin ) bathy(ji,jj)=gdepw(jk+1)-0.1
          ENDIF
       ENDDO
    END DO
  END SUBROUTINE zgr_zps


  SUBROUTINE zgr_read()
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE zgr_read  ***
    !!
    !! ** Purpose :  Read zgrbat.txt file (cf_batfile) to set the gdep[tw]_0
    !!               and e3[tw] 
    !!
    !! ** Method  :  Read the ocean output format ( ie, cf_batfile is just 
    !!               a copy of the ocean.output concerning zgrbat 
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4) :: inumzgr = 10, il, iostat, idum, ifoo
    CHARACTER(LEN=256) :: cline, clfile
    !!----------------------------------------------------------------------
    clfile = cf_batfile       ! defined in the main program
    il=0
    OPEN(inumzgr, FILE=clfile,IOSTAT=iostat)

    DO WHILE ( iostat == 0 )
       READ(inumzgr,'(a)',IOSTAT=iostat) cline
       READ(cline,*,IOSTAT=idum )il
       IF ( idum == 0 ) npk=il
    END DO

    ALLOCATE ( gdept(npk), gdepw(npk), e3t(npk), e3w(npk) )
    REWIND(inumzgr)

    il=0 ; iostat=0
    DO WHILE ( iostat == 0 )
       READ(inumzgr,'(a)', IOSTAT=iostat) cline
       READ(cline,*,IOSTAT=idum) il
       IF ( idum == 0  ) READ(cline,*) ifoo, gdept(il), gdepw(il), &
            &             e3t(il), e3w(il)
    END DO
  END SUBROUTINE zgr_read


  SUBROUTINE prlog (ptabold, ptab ,kpi, kpj, ldapp)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE prlog  ***
    !!
    !! ** Purpose :   Print a fortran 90 log file describing the modifications
    !!              done to the bathymetry
    !!
    !! ** Method  :  File is append instead of created if ldapp true
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(:,:), INTENT(in) :: ptabold  ! original array
    REAL(KIND=4), DIMENSION(:,:), INTENT(in) :: ptab     ! modified array
    INTEGER(KIND=4),              INTENT(in) :: kpi, kpj ! size of the array
    LOGICAL,                      INTENT(in) :: ldapp    ! append flag

    INTEGER(KIND=4)        :: ji, jj
    INTEGER(KIND=4)        :: inumlog=10
    CHARACTER(LEN=80)      :: clfile
    !!----------------------------------------------------------------------
    clfile = cf_log

    IF (ldapp ) THEN 
       OPEN (inumlog, FILE=clfile, POSITION='append')
    ELSE
       OPEN (inumlog, FILE=clfile)
    ENDIF

    WRITE(inumlog,'(a,a)') '! modification from original file : ', TRIM(cf_in)
    WRITE(inumlog,'(a,a)') '! written to : ', TRIM(cwkc)
    DO ji=1,kpi
       DO jj=1,kpj
          IF ( ABS( ptabold(ji,jj) - ptab(ji,jj)) > 0.02  ) THEN    ! allow a 2 cm tolerance for rounding purposes
             WRITE(inumlog,'(a,i4,a,i4,a,f8.2,a,f8.2)') ' bathy(',ji,',',jj,')=',ptab(ji,jj)*scale_factor,  &
                  & ' ! instead of ',ptabold(ji,jj)*scale_factor
          END IF
       END DO
    END DO

    CLOSE(inumlog)
  END SUBROUTINE prlog


  SUBROUTINE fillzone(kimin, kimax, kjmin, kjmax)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE fillzone  ***
    !!
    !! ** Purpose :  Fill a subarea with 0 up to encounter a coast on the East
    !!
    !! ** Method  :  Assume that first point is sea point. Mask it and do so with
    !!               all points to the east (j=cst) up to a land point. 
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) :: kimin, kimax, kjmin, kjmax ! position of the data windows

    INTEGER(KIND=4) :: jj 
    INTEGER(KIND=4) :: ii
    !!----------------------------------------------------------------------
    DO jj=kjmin,kjmax
       ii=kimin
       IF ( bathy(ii,jj)  /= 0 ) THEN
          DO WHILE ( bathy(ii,jj)  /= 0 .AND. ii <= kimax )
             bathy(ii,jj) = 0.
             ii=ii+1
          END DO
       END IF
    END DO
  END SUBROUTINE fillzone


  SUBROUTINE raz_zone(kimin, kimax, kjmin, kjmax)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE raz_zone  ***
    !!
    !! ** Purpose :  Fill a sub area of a bathy file with 0 
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) :: kimin, kimax, kjmin, kjmax ! position of the data windows
    !!----------------------------------------------------------------------
    bathy(kimin:kimax, kjmin:kjmax) = 0.

  END SUBROUTINE raz_zone


  SUBROUTINE raz_below(kimin, kimax, kjmin, kjmax, pdepmin)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE raz_below  ***
    !!
    !! ** Purpose : Fill point (set to 0) that are below pdepmin
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) :: kimin, kimax, kjmin, kjmax ! position of the data windows
    REAL(KIND=4),    INTENT(in) :: pdepmin                    ! threshold bathy value
    !!----------------------------------------------------------------------
    WHERE ( bathy(kimin:kimax, kjmin:kjmax) <= pdepmin)  bathy(kimin:kimax, kjmin:kjmax) = 0.

  END SUBROUTINE raz_below


  SUBROUTINE set_below(kimin, kimax, kjmin, kjmax, pdepmin)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE set_below  ***
    !!
    !! ** Purpose : Set bathy points to pdepmin if less than pdepmin in the
    !!              original bathy
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) :: kimin, kimax, kjmin, kjmax ! position of the data windows
    REAL(KIND=4),    INTENT(in) :: pdepmin                    ! threshold bathy value
    !!----------------------------------------------------------------------
    WHERE ( bathy(kimin:kimax, kjmin:kjmax) <= pdepmin .AND. bathy(kimin:kimax, kjmin:kjmax) > 0 ) &
         &                 bathy(kimin:kimax, kjmin:kjmax) = pdepmin

  END SUBROUTINE set_below

  SUBROUTINE set_above(kimin, kimax, kjmin, kjmax, pdepmax)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE set_above  ***
    !!
    !! ** Purpose : Set bathy points pdepmax if larger than pdepmax
    !!              original bathy
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) :: kimin, kimax, kjmin, kjmax ! position of the data windows
    REAL(KIND=4),    INTENT(in) :: pdepmax                    ! threshold bathy value
    !!----------------------------------------------------------------------
    WHERE ( bathy(kimin:kimax, kjmin:kjmax) >= pdepmax ) &
         &                 bathy(kimin:kimax, kjmin:kjmax) = pdepmax

  END SUBROUTINE set_above


  SUBROUTINE set_zone(kimin, kimax, kjmin, kjmax, pdepmin)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE set_zone  ***
    !!
    !! ** Purpose : Set bathy points to pdepmin 
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) :: kimin, kimax, kjmin, kjmax ! position of the data windows
    REAL(KIND=4),    INTENT(in) :: pdepmin                    ! threshold bathy value
    !!----------------------------------------------------------------------
    bathy(kimin:kimax, kjmin:kjmax) =pdepmin

  END SUBROUTINE set_zone



  SUBROUTINE dumpzone(cdumpf, kimin, kimax, kjmin, kjmax)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE dumpzone  ***
    !!
    !! ** Purpose :  Print subarea to cdumpf ascii file 
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdumpf                     ! name of the dump file
    INTEGER(KIND=4),  INTENT(in) :: kimin, kimax, kjmin, kjmax ! position of the data windows

    INTEGER(KIND=4)    :: ji, jj
    INTEGER(KIND=4)    :: inumdmp=20 , ini
    CHARACTER(LEN=256) :: cl_fmtr, cl_fmti
    !!----------------------------------------------------------------------
    ini = kimax - kimin + 1
    WRITE(cl_fmtr,99) ini
    WRITE(cl_fmti,98) ini
    OPEN(inumdmp,FILE=cdumpf)
    WRITE(inumdmp,*) kimin, kimax, kjmin, kjmax, TRIM(cl_fmtr)
99  FORMAT('(I5,',i4.4,'f8.2)')
98  FORMAT('(5x,',i4.4,'I8)')
    WRITE(inumdmp,cl_fmti)(ji,ji=kimin,kimax)
    DO jj= kjmax,kjmin,-1
       WRITE(inumdmp,cl_fmtr) jj, bathy(kimin:kimax,jj)
    ENDDO
    CLOSE(inumdmp)

  END SUBROUTINE dumpzone


  SUBROUTINE nicedumpzone(cdumpf, kimin, kimax, kjmin, kjmax)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE nicedumpzone  ***
    !!
    !! ** Purpose : Print subarea to cdumpf ascii file with a nice format  
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdumpf                     ! name of the dump file
    INTEGER(KIND=4),  INTENT(in) :: kimin, kimax, kjmin, kjmax ! position of the data windows

    INTEGER(KIND=4)    :: ji, jj
    INTEGER(KIND=4)    :: inumdmp=20 , ini
    CHARACTER(LEN=256) :: cl_fmtr, cl_fmti
    !!----------------------------------------------------------------------
    ini=kimax-kimin+1
    WRITE(cl_fmtr,99) ini
    WRITE(cl_fmti,98) ini
    OPEN(inumdmp,FILE=cdumpf)
    WRITE(inumdmp,*) kimin,kimax,kjmin,kjmax, TRIM(cl_fmtr)
99  FORMAT('(I5,',i4.4,'I5)')
98  FORMAT('(5x,',i4.4,'I5)')
    WRITE(inumdmp,cl_fmti)(ji,ji=kimin,kimax)
    DO jj= kjmax,kjmin,-1
       WRITE(inumdmp,cl_fmtr) jj, INT(bathy(kimin:kimax,jj))
       WRITE(inumdmp,*)
       WRITE(inumdmp,*)
    ENDDO
    CLOSE(inumdmp)

  END SUBROUTINE nicedumpzone


  SUBROUTINE replacezone(cdreplace)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE replacezone  ***
    !!
    !! ** Purpose :  Replace a bathy area by data read from an ascii input file
    !!               formely generated by -dump option (and manualy modified) 
    !!
    !! ** Method  :  Read format in the header part of the file 
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdreplace

    INTEGER(KIND=4) :: jj
    INTEGER(KIND=4) :: iimin, iimax, ijmin, ijmax
    INTEGER(KIND=4) :: inumrep=20, idum
    !!----------------------------------------------------------------------
    OPEN(inumrep,FILE=cdreplace)
    READ(inumrep,*) iimin, iimax, ijmin, ijmax
    READ(inumrep,*)  ! skip 1 line
    DO jj=ijmax,ijmin,-1
       READ(inumrep,*) idum, bathy(iimin:iimax,jj)
    END DO
    CLOSE(inumrep)

  END SUBROUTINE replacezone


  SUBROUTINE fillpool(kcrit, kimin, kimax, kjmin, kjmax) 
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE replacezone  ***
    !!
    !! ** Purpose :  Replace all area surrounding by mask value by mask value
    !!
    !! ** Method  :  flood fill algorithm
    !!
    !!----------------------------------------------------------------------
    INTEGER, INTENT(in) :: kcrit        ! maximal allowed pool 
    INTEGER(KIND=4),  INTENT(in) :: kimin, kimax, kjmin, kjmax ! position of the data windows

    INTEGER :: ik                       ! number of point change
    INTEGER :: ip                       ! size of the pile
    INTEGER :: ji, jj                   ! loop index
    INTEGER :: iip1, iim1, ii, ij       ! working integer
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ipile    ! pile variable
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ioptm    ! matrix to check already tested value 

    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zbathy   ! new bathymetry
    !!----------------------------------------------------------------------
    PRINT *, 'WARNING North fold case not coded'
    ! allocate variable
    ALLOCATE(ipile(((kimax-kimin)+1)*((kjmax-kjmin)+1),2))
    ALLOCATE(zbathy(npiglo,npjglo), ioptm(npiglo,npjglo))

    ioptm = bathy
    WHERE (ioptm /=  0)
       ioptm = 1
    END WHERE

    PRINT *, 'Filling area in progress ... (it can take a while)'    

    DO ji=kimin,kimax
       IF (MOD(ji,100) == 0) PRINT *, ji,'/',npiglo
       DO jj=kjmin,kjmax
          ! modify something only if seed point is a non 0 cell
          IF (ioptm(ji,jj) == 1) THEN
             ! initialise variables
             zbathy=bathy
             ipile(:,:)=0
             ipile(1,:)=[ji,jj]
             ip=1; ik=0

             ! loop until the pile size is 0 or if the pool is larger than the critical size
             DO WHILE ( ip /= 0 .AND. ik < kcrit);
                ik=ik+1 
                ii=ipile(ip,1); ij=ipile(ip,2)

                ! update bathy and update pile size
                zbathy(ii,ij)=0.0 
                ipile(ip,:)  =[0,0]; ip=ip-1

                ! check neighbour cells and update pile
                iip1=ii+1; IF ( iip1 == npiglo+1 ) iip1=2
                iim1=ii-1; IF ( iim1 == 0        ) iim1=npiglo-1
                IF (zbathy(ii, ij+1) /=  0.0) THEN
                   ip=ip+1; ipile(ip,:)=[ii  ,ij+1] 
                   ioptm (ii, ij+1) = 0
                END IF
                IF (zbathy(ii, ij-1) /= 0.0) THEN
                   ip=ip+1; ipile(ip,:)=[ii  ,ij-1]
                   ioptm(ii, ij-1) = 0
                END IF
                IF (zbathy(iip1, ij) /=  0.0) THEN
                   ip=ip+1; ipile(ip,:)=[iip1,ij  ]
                   ioptm(iip1, ij) = 0
                END IF
                IF (zbathy(iim1, ij) /=  0.0) THEN
                   ip=ip+1; ipile(ip,:)=[iim1,ij  ]
                   ioptm(iim1, ij) = 0
                END IF
             END DO
             IF (ik < kcrit) bathy=zbathy;
          END IF
       END DO
    END DO

    DEALLOCATE(ipile); DEALLOCATE(zbathy, ioptm)

  END SUBROUTINE fillpool


END PROGRAM cdfbathy
