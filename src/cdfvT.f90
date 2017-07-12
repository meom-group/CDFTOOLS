PROGRAM cdfvT
  !!======================================================================
  !!                     ***  PROGRAM  cdfvT  ***
  !!=====================================================================
  !!  ** Purpose : Compute the average values for the products 
  !!               V.T, V.S, U.T and U.S, used afterward for heat and
  !!               salt transport.
  !!
  !!  ** Method  : pass the CONFIG name and a series of tags as arguments.
  !!
  !! History : 2.1  : 11/2004  : J.M. Molines : Original code
  !!           2.1  : 02/2010  : J.M. Molines : handle multiframes input files.
  !!           3.0  : 04/2011  : J.M. Molines : Doctor norm + Lic.
  !!                : 10/2012  : M. Balmaseda : Split T and S file eventually
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils       ! SetFileName function
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class second_order_moments
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jk, jv, jtt, jt  ! dummy loop index
  INTEGER(KIND=4)                           :: ierr                 ! working integer
  INTEGER(KIND=4)                           :: narg, iargc          ! command line
  INTEGER(KIND=4)                           :: ijarg                ! command line
  INTEGER(KIND=4)                           :: npiglo,npjglo        ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt             ! size of the domain
  INTEGER(KIND=4)                           :: ntframe              ! Cumul of time frame
  INTEGER(KIND=4)                           :: ntags                ! number of tags in the list
  INTEGER(KIND=4)                           :: ncout                ! ncid of output file
  INTEGER(KIND=4)                           :: nvaro                ! Number of output variables
  INTEGER(KIND=4), DIMENSION(:),ALLOCATABLE :: ipk, id_varout       ! level and varid's of output vars

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztemp, zsal          ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zu, zv               ! Velocity component
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zworku, zworkv       ! working arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e3u, e3v             ! vertical metrics for vvl case

  REAL(KIND=8), DIMENSION(1)                :: dtimean              ! mean time
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim                 ! time counter of individual files
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulut, dcumulus   ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulvt, dcumulvs   ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumule3u, dcumule3v ! Arrays for cumulated vertical metrics (vvl case)
  REAL(KIND=8)                              :: dtotal_time          ! cumulated time

  CHARACTER(LEN=256)                        :: cf_tfil              ! T file name
  CHARACTER(LEN=256)                        :: cf_sfil              ! S file name [default: idem T file]
  CHARACTER(LEN=256)                        :: cf_ufil              ! zonal velocity file
  CHARACTER(LEN=256)                        :: cf_vfil              ! meridional velocity file
  CHARACTER(LEN=256)                        :: cf_out='vt.nc'       ! output file
  CHARACTER(LEN=256)                        :: config               ! configuration name
  CHARACTER(LEN=256)                        :: ctag                 ! current tag to work with               
  CHARACTER(LEN=256)                        :: cldum                ! dummy character argument
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: ctag_lst         ! dummy character argument

  TYPE (variable), DIMENSION(:), ALLOCATABLE:: stypvar          ! structure for attributes

  LOGICAL                                   :: lcaltmean            ! flag for mean time computation
  LOGICAL                                   :: lnc4=.FALSE.         ! flag for netcdf4 output with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfvT -c CONFIG-CASE -l LST-tags [-o OUT-file ] [-nc4 ] [-vvl]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the time average values for second order moments V.T, V.S, U.T' 
     PRINT *,'       and U.S used in heat and salt transport computation.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -c CONFIG-CASE is the config name of a given experiment (eg ORCA025-G70)'
     PRINT *,'            The program will look for gridT, gridU and gridV files for this ' 
     PRINT *,'            config ( grid_T, grid_U and grid_V are also accepted).'
     PRINT *,'            Additionaly, if gridS or grid_S file is found, it will be taken'
     PRINT *,'            in place of gridT for the salinity variable.'
     PRINT *,'       -l LST-tags : a blank-separated list of time tags that will be used for'
     PRINT *,'            time averaging, e.g.:  y2000m01d05 y2000m01d10 ...'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file] : specify output file name instead of ', TRIM(cf_out)
     PRINT *,'       [-nc4 ] use netcdf4 output with chunking and deflation 1'
     PRINT *,'       [-vvl ] use time varying vertical metrics.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option used.'
     PRINT *,'       variables : ',TRIM(cn_vozout),' time-mean ( U.T ) at U point.'
     PRINT *,'                   ',TRIM(cn_vozous),' time-mean ( U.S ) at U point.'
     PRINT *,'                   ',TRIM(cn_vomevt),' time-mean ( V.T ) at V point.'
     PRINT *,'                   ',TRIM(cn_vomevs),' time-mean ( V.S ) at V point.'
     PRINT *,'      '
     STOP 
  ENDIF

  !! Initialisation from 1st file (all file are assume to have the same geometry)
  ijarg = 1 
  DO WHILE ( ijarg <= narg ) 
     CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum ) 
     CASE ( '-c'   ) ; CALL getarg (ijarg, config) ; ijarg=ijarg+1
     CASE ( '-l'   ) ; CALL GetTagList
        ! options
     CASE ( '-o'   ) ; CALL getarg (ijarg, cf_out) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4   = .TRUE.
     CASE ( '-vvl' ) ; lg_vvl = .TRUE.
     CASE DEFAULT    ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  cf_tfil = SetFileName( config, ctag_lst(1), 'T')

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk

  IF ( lg_vvl ) THEN
     nvaro = 6  ! UT US VT VS e3u, e3v
     cn_ve3u = cn_ve3uvvl
     cn_ve3v = cn_ve3vvvl
  ELSE
     nvaro = 4  ! UT US VT VS
  ENDIF

  ALLOCATE( dcumulut(npiglo,npjglo), dcumulus(npiglo,npjglo) )
  ALLOCATE( dcumulvt(npiglo,npjglo), dcumulvs(npiglo,npjglo) )
  ALLOCATE( zu(npiglo,npjglo),    zv(npiglo,npjglo) )
  ALLOCATE( zworku(npiglo,npjglo),   zworkv(npiglo,npjglo) )
  ALLOCATE( ztemp(npiglo,npjglo),    zsal(npiglo,npjglo) )
  ALLOCATE( id_varout(nvaro), ipk(nvaro), stypvar(nvaro) )
  IF (lg_vvl ) THEN
     ALLOCATE(       e3u(npiglo,npjglo),       e3v(npiglo,npjglo))
     ALLOCATE( dcumule3u(npiglo,npjglo), dcumule3v(npiglo,npjglo) )
  ENDIF

  CALL CreateOutput

  lcaltmean=.TRUE.
  DO jk = 1, npk
     PRINT *,'level ',jk
     dcumulut(:,:) = 0.d0 ;  dcumulvt(:,:) = 0.d0 ; dtotal_time = 0.d0
     dcumulus(:,:) = 0.d0 ;  dcumulvs(:,:) = 0.d0 ; ntframe = 0
     IF ( lg_vvl ) THEN
        dcumule3u(:,:) = 0.d0 ;  dcumule3v(:,:) = 0.d0
     ENDIF

     DO jt = 1, ntags
        ctag=ctag_lst(jt) 
        cf_tfil = SetFileName( config, ctag, 'T', ld_stop=.TRUE. )
        cf_sfil = SetFileName( config, ctag, 'S', ld_stop=.FALSE.)      ! do not stop if gridS/grid_S not found !
        IF ( chkfile (cf_sfil, ld_verbose=.FALSE.) ) cf_sfil = cf_tfil  ! do not complain if not found

        npt = getdim (cf_tfil, cn_t)
        IF ( lcaltmean ) THEN
           ALLOCATE ( dtim(npt) )
           dtim        = getvar1d(cf_tfil, cn_vtimec, npt)
           dtotal_time = dtotal_time + SUM(dtim(1:npt) )
           DEALLOCATE( dtim )
        END IF

        ! assume U and V file have same time span ...
        cf_ufil = SetFileName( config, ctag, 'U' )
        cf_vfil = SetFileName( config, ctag, 'V' )

        DO jtt = 1, npt  ! loop on time frame in a single file
           ntframe = ntframe+1
           zu(:,:)    = getvar(cf_ufil, cn_vozocrtx, jk, npiglo, npjglo, ktime=jtt )
           zv(:,:)    = getvar(cf_vfil, cn_vomecrty, jk, npiglo, npjglo, ktime=jtt )
           ztemp(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jtt )
           zsal(:,:)  = getvar(cf_sfil, cn_vosaline, jk, npiglo, npjglo, ktime=jtt )

           IF ( lg_vvl ) THEN
              e3u(:,:) = getvar(cf_ufil, cn_ve3u, jk, npiglo, npjglo, ktime=jtt )
              e3v(:,:) = getvar(cf_vfil, cn_ve3v, jk, npiglo, npjglo, ktime=jtt )
           ENDIF

           ! temperature at u point, v points ( As in NEMO, with or without  vvl Temp a U/V points is 
           ! computed as follow ( second order) 
           zworku(:,:) = 0. ; zworkv(:,:) = 0.
           DO ji=1, npiglo-1
              DO jj = 1, npjglo -1
                 zworku(ji,jj) = 0.5 * ( ztemp(ji,jj) + ztemp(ji+1,jj) )  ! temper at Upoint
                 zworkv(ji,jj) = 0.5 * ( ztemp(ji,jj) + ztemp(ji,jj+1) )  ! temper at Vpoint
              END DO
           END DO

           IF (lg_vvl ) THEN
              dcumulut(:,:) = dcumulut(:,:) + zworku(:,:) * e3u(:,:) * zu(:,:)*1.d0
              dcumulvt(:,:) = dcumulvt(:,:) + zworkv(:,:) * e3v(:,:) * zv(:,:)*1.d0
              dcumule3u(:,:) = dcumule3u(:,:) +  e3u(:,:) *1.d0
              dcumule3v(:,:) = dcumule3v(:,:) +  e3v(:,:) *1.d0
           ELSE
              dcumulut(:,:) = dcumulut(:,:) + zworku(:,:) * zu(:,:)*1.d0
              dcumulvt(:,:) = dcumulvt(:,:) + zworkv(:,:) * zv(:,:)*1.d0
           ENDIF

           ! salinity at u points, v points
           zworku(:,:) = 0. ; zworkv(:,:) = 0.
           DO ji=1, npiglo-1
              DO jj = 1, npjglo -1
                 zworku(ji,jj) = 0.5 * ( zsal(ji,jj) + zsal(ji+1,jj) )  ! salinity  at Upoint
                 zworkv(ji,jj) = 0.5 * ( zsal(ji,jj) + zsal(ji,jj+1) )  ! salinity  at Vpoint
              END DO
           END DO

           IF (lg_vvl ) THEN
              dcumulus(:,:) = dcumulus(:,:) + zworku(:,:) * e3u(:,:) * zu(:,:)*1.d0
              dcumulvs(:,:) = dcumulvs(:,:) + zworkv(:,:) * e3v(:,:) * zv(:,:)*1.d0
           ELSE
              dcumulus(:,:) = dcumulus(:,:) + zworku(:,:) * zu(:,:)*1.d0
              dcumulvs(:,:) = dcumulvs(:,:) + zworkv(:,:) * zv(:,:)*1.d0
           ENDIF

        END DO  !jtt
     END DO  ! arg list
     ! finish with level jk ; compute mean (assume spval is 0 )
     ! JMM In the following line multiplication by 1.e0 is a trick to write real*4 values.
     IF ( lg_vvl) THEN
        ierr = putvar(ncout, id_varout(1), dcumulvt(:,:)/dcumule3v(:,:)*1.e0, jk,npiglo, npjglo, kwght=ntframe )
        ierr = putvar(ncout, id_varout(2), dcumulvs(:,:)/dcumule3v(:,:)*1.e0, jk,npiglo, npjglo, kwght=ntframe )
        ierr = putvar(ncout, id_varout(3), dcumulut(:,:)/dcumule3u(:,:)*1.e0, jk,npiglo, npjglo, kwght=ntframe )
        ierr = putvar(ncout, id_varout(4), dcumulus(:,:)/dcumule3u(:,:)*1.e0, jk,npiglo, npjglo, kwght=ntframe )
        ierr = putvar(ncout, id_varout(5), dcumule3u(:,:)/ntframe      *1.e0, jk,npiglo, npjglo, kwght=ntframe )
        ierr = putvar(ncout, id_varout(6), dcumule3v(:,:)/ntframe      *1.e0, jk,npiglo, npjglo, kwght=ntframe )
     ELSE
        ierr = putvar(ncout, id_varout(1), dcumulvt(:,:)/ntframe*1.d0, jk,npiglo, npjglo, kwght=ntframe )
        ierr = putvar(ncout, id_varout(2), dcumulvs(:,:)/ntframe*1.d0, jk,npiglo, npjglo, kwght=ntframe )
        ierr = putvar(ncout, id_varout(3), dcumulut(:,:)/ntframe*1.d0, jk,npiglo, npjglo, kwght=ntframe )
        ierr = putvar(ncout, id_varout(4), dcumulus(:,:)/ntframe*1.d0, jk,npiglo, npjglo, kwght=ntframe )
     ENDIF

     IF (lcaltmean )  THEN
        dtimean(1) = dtotal_time/ntframe
        ierr       = putvar1d(ncout, dtimean, 1, 'T')
     END IF
     lcaltmean=.FALSE. ! tmean already computed

  END DO  ! loop to next level

  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------

    ipk(:)= npk  ! all variables (input and output are 3D)
    ! define output variables
    stypvar%rmissing_value    = 0.
    stypvar%valid_min         = -100.
    stypvar%valid_max         = 100.
    stypvar%conline_operation = 'N/A'
    stypvar%caxis             = 'TZYX'

    DO jv = 1, nvaro
       stypvar(jv)%ichunk = (/npiglo,MAX(1,npjglo/30), 1, 1 /)
    ENDDO

    stypvar(1)%cname          = cn_vomevt       ; stypvar(1)%cunits        = 'm.DegC.s-1'
    stypvar(2)%cname          = cn_vomevs       ; stypvar(2)%cunits        = 'm.PSU.s-1'
    stypvar(3)%cname          = cn_vozout       ; stypvar(3)%cunits        = 'm.DegC.s-1'
    stypvar(4)%cname          = cn_vozous       ; stypvar(4)%cunits        = 'm.PSU.s-1'

    stypvar(1)%clong_name     = 'Meridional_VT' ; stypvar(1)%cshort_name   = cn_vomevt
    stypvar(2)%clong_name     = 'Meridional_VS' ; stypvar(2)%cshort_name   = cn_vomevs
    stypvar(3)%clong_name     = 'Zonal_UT'      ; stypvar(3)%cshort_name   = cn_vozout
    stypvar(4)%clong_name     = 'Zonal_US'      ; stypvar(4)%cshort_name   = cn_vozous

    IF ( lg_vvl ) THEN
       stypvar(5)%cname       = cn_ve3u         ; stypvar(5)%cunits        = 'm'
       stypvar(6)%cname       = cn_ve3v         ; stypvar(6)%cunits        = 'm'

       stypvar(5)%clong_name  = 'Mean_e3u'      ; stypvar(5)%cshort_name   = cn_ve3u
       stypvar(6)%clong_name  = 'Mean_e3v'      ; stypvar(6)%cshort_name   = cn_ve3v
    ENDIF

    ! create output fileset
    ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk, ld_xycoo=.TRUE. , ld_nc4=lnc4 )
    ierr  = createvar   (ncout , stypvar, nvaro,      ipk,    id_varout        , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk, ld_xycoo=.TRUE. )

  END SUBROUTINE CreateOutput

  SUBROUTINE GetTagList
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetTagList  ***
    !!
    !! ** Purpose :  Set up a tag list given on the command line as 
    !!               blank separated list
    !!
    !! ** Method  :  Scan the command line until a '-' is found
    !!----------------------------------------------------------------------
    INTEGER (KIND=4)  :: ji
    INTEGER (KIND=4)  :: icur
    !!----------------------------------------------------------------------
    !!
    ntags=0
    ! need to read a list of file ( number unknow ) 
    ! loop on argument till a '-' is found as first char
    icur=ijarg                          ! save current position of argument number
    DO ji = icur, narg                  ! scan arguments till - found
       CALL getarg ( ji, cldum )
       IF ( cldum(1:1) /= '-' ) THEN ; ntags = ntags+1
       ELSE                          ; EXIT
       ENDIF
    ENDDO
    ALLOCATE (ctag_lst(ntags) )
    DO ji = icur, icur + ntags -1
       CALL getarg(ji, ctag_lst( ji -icur +1 ) ) ; ijarg=ijarg+1
    END DO
  END SUBROUTINE GetTagList

END PROGRAM cdfvT
