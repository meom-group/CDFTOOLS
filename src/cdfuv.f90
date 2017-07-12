PROGRAM cdfuv
  !!======================================================================
  !!                     ***  PROGRAM  cdfuv  ***
  !!=====================================================================
  !!  ** Purpose : Compute the average values for the products 
  !!               u.v at T-point
  !!
  !!  ** Method  : pass the CONFIG name and a series of tags as arguments.
  !!
  !! History : 2.1  : 11/2004  : J.M. Molines : Original code
  !!           2.1  : 02/2010  : J.M. Molines : handle multiframes input files.
  !!           3.0  : 04/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!
  !! * Reference for mean and variance computing:
  !! Higham N. (2002) :"Accuracy and Stability of Numerical Algorithms." 
  !!         Second edition, Siam editor, 710pp. (see chapter 1.9 )
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

  INTEGER(KIND=4)                           :: ji, jj, jk, jt, jtt  ! dummy loop index
  INTEGER(KIND=4)                           :: ierr                 ! working integer
  INTEGER(KIND=4)                           :: narg, iargc, ijarg   ! command line
  INTEGER(KIND=4)                           :: npiglo,npjglo        ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt, npkf       ! size of the domain
  INTEGER(KIND=4)                           :: ntframe              ! Cumul of time frame
  INTEGER(KIND=4)                           :: ntags                ! number of tags to process
  INTEGER(KIND=4)                           :: ncout                ! ncid of output file
  INTEGER(KIND=4), DIMENSION(4)             :: ipk, id_varout       ! level and varid's of output vars

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zu, zv               ! Velocity component
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zworku, zworkv       ! working arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmean                ! temporary mean value for netcdf write
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zlon                 ! longitude of T points to check periodicity

  REAL(KIND=8), DIMENSION(1)                :: dtimean              ! mean time
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim                 ! time counter of individual files
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumuluv             ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulu              ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulv              ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dupvp                ! Arrays for U'.V'
  REAL(KIND=8)                              :: dtotal_time          ! cumulated time

  CHARACTER(LEN=256)                        :: cf_tfil              ! Temperature file for reference only
  CHARACTER(LEN=256)                        :: cf_ufil              ! zonal velocity file
  CHARACTER(LEN=256)                        :: cf_vfil              ! meridional velocity file
  CHARACTER(LEN=256)                        :: cf_out='uv.nc'       ! output file
  CHARACTER(LEN=256)                        :: config               ! configuration name
  CHARACTER(LEN=256)                        :: ctag                 ! current tag to work with               
  CHARACTER(LEN=256)                        :: cl_name              ! temporary variable name
  CHARACTER(LEN=256)                        :: cldum                ! working char variable
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: ctag_lst         ! tag list

  TYPE (variable), DIMENSION(4)             :: stypvar              ! structure for attributes

  LOGICAL                                   :: lcaltmean            ! flag for mean time computation
  LOGICAL                                   :: lperio = .FALSE.     ! flag for E-W periodicity
  LOGICAL                                   :: lnc4   = .FALSE.     ! Use nc4 with chunking and deflation
  LOGICAL                                   :: lopt   = .FALSE.     ! optimized algorithm flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfuv -c CONFIG-CASE -l LST-tags [-opt] [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the time average values for U.V product, at T point.' 
     PRINT *,'       Mean U, mean V  and mean U''.V''  at T point are also computed.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -c CONFIG-CASE is the config name of a given experiment (eg ORCA025-G70)'
     PRINT *,'            The program will look for gridU and gridV files for this config'
     PRINT *,'            grid_U and grid_V are also accepted).'
     PRINT *,'       -l LST-tags : a list of time tags that will be used for time averaging.'
     PRINT *,'            e.g. :  y2000m01d05 y2000m01d10 ...'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-opt :] use optimized algorithm, minimizing truncation errors in the'
     PRINT *,'             evaluation of mean U, meanV and mean (U''.V'').'
     PRINT *,'       [-o OUT-file] : specify output filename instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ] :  Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'             This option is effective only if cdftools are compiled with'
     PRINT *,'             a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'       variables : ',TRIM(cn_vouv), '  : Mean U.V at T point'
     PRINT *,'                   ',TRIM(cn_vozocrtx)//'_t : Mean U at T point'
     PRINT *,'                   ',TRIM(cn_vomecrty)//'_t : Mean V at T point'
     PRINT *,'                   ',TRIM(cn_vouv)//'_prime : Mean U''.V'' at T point'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg=1
  DO WHILE (ijarg <= narg) 
     CALL getarg(ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-c'  ) ; CALL getarg(ijarg, config) ; ijarg=ijarg+1
     CASE ( '-l'  ) ; CALL GetTagList
        ! options
     CASE ( '-opt') ; lopt = .TRUE.
     CASE ( '-o'  ) ; CALL getarg(ijarg, cf_out) ; ijarg=ijarg+1
     CASE ( '-nc4') ; lnc4 = .TRUE.
     CASE DEFAULT   ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  cf_tfil = SetFileName( config, ctag_lst(1), 'T')
  cf_ufil = SetFileName( config, ctag_lst(1), 'U')

  npiglo = getdim (cf_ufil,cn_x)
  npjglo = getdim (cf_ufil,cn_y)
  npkf   = getdim (cf_ufil,cn_z)

  IF ( npkf == 0 ) THEN ; npk = 1
  ELSE                  ; npk=npkf
  ENDIF

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk

  ALLOCATE( dcumuluv(npiglo,npjglo) )
  ALLOCATE( dcumulu(npiglo,npjglo) )
  ALLOCATE( dcumulv(npiglo,npjglo) )
  ALLOCATE( dupvp(npiglo,npjglo) )
  ALLOCATE( zu(npiglo,npjglo),     zv(npiglo,npjglo) )
  ALLOCATE( zworku(npiglo,npjglo), zworkv(npiglo,npjglo) )
  ALLOCATE( zmean(npiglo,npjglo))
  ALLOCATE( zlon(npiglo,npjglo))

  ! check for E_W periodicity
  zlon(:,:) = getvar(cf_tfil, cn_vlon2d, 1, npiglo, npjglo )
  IF ( zlon(1,1) ==  zlon(npiglo-1,1) ) THEN
     lperio = .TRUE.
     PRINT *,' E-W periodicity detected '
  ENDIF

  CALL CreateOutput

  lcaltmean=.TRUE.
  DO jk = 1, npk
     PRINT *,'level ',jk
     dcumuluv(:,:) = 0.d0 ;  dtotal_time  = 0.d0 ; ntframe = 0
     dcumulu(:,:)  = 0.d0 ;  dcumulv(:,:) = 0.d0

     DO jt = 1, ntags           ! loop on tags
        ctag=ctag_lst(jt)

        cf_ufil = SetFileName( config, ctag, 'U', ld_stop=.TRUE. )
        npt = getdim (cf_ufil, cn_t)
        IF ( lcaltmean ) THEN
           ALLOCATE ( dtim(npt) )
           dtim        = getvar1d(cf_ufil, cn_vtimec, npt)
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

           ! U and V velocities at T point
           zworku(:,:) = 0. ; zworkv(:,:) = 0.
           DO ji=2, npiglo
              DO jj = 1, npjglo
                 zworku(ji,jj) = 0.5 * ( zu(ji,jj) + zu(ji-1,jj) )  ! U at T point
              END DO
           END DO

           DO ji=1, npiglo
              DO jj = 2, npjglo 
                 zworkv(ji,jj) = 0.5 * ( zv(ji,jj) + zv(ji,jj-1) )  ! V at T point
              END DO
           END DO
           IF ( lopt ) THEN
              IF ( ntframe == 1 ) THEN  ! initialize recurence formula
                dcumulu(:,:) = zworku(:,:)
                dcumulv(:,:) = zworkv(:,:)
                dupvp(:,:)   = 0.d0
                dcumuluv(:,:) =zworku(:,:) * zworkv(:,:)*1.d0
              ELSE
                dupvp(:,:)= dupvp(:,:) + (1.*(ntframe -1))/ntframe*(zworku(:,:) - dcumulu(:,:))*(zworkv(:,:) - dcumulv(:,:))
                dcumulu(:,:) = dcumulu(:,:)  +  1./ntframe*(zworku(:,:)            *1.d0 - dcumulu(:,:))
                dcumulv(:,:) = dcumulv(:,:)  +  1./ntframe*(zworkv(:,:)            *1.d0 - dcumulv(:,:))
                dcumuluv(:,:) = dcumuluv(:,:)+  1./ntframe*(zworku(:,:)*zworkv(:,:)*1.d0 - dcumuluv(:,:))
              ENDIF
           ELSE
             dcumuluv(:,:) = dcumuluv(:,:) + zworku(:,:) * zworkv(:,:)*1.d0
             dcumulu(:,:)  = dcumulu(:,:)  + zworku(:,:)*1.d0
             dcumulv(:,:)  = dcumulv(:,:)  + zworkv(:,:)*1.d0
           ENDIF

        END DO  !jtt
     END DO  ! jt
     ! finish with level jk ; compute mean (assume spval is 0 ) 
     !  < > is temporal mean here
     ! < U.V> 
     IF ( lopt) THEN ; zmean(:,:) = dcumuluv(:,:)
     ELSE            ; zmean(:,:) = dcumuluv(:,:)/ntframe
     ENDIF
     IF ( lperio ) zmean(1, :) = zmean(npiglo-1,:)
     ierr = putvar(ncout, id_varout(1), zmean, jk,npiglo, npjglo, kwght=ntframe )

     ! < U >
     IF ( lopt ) THEN ; zmean(:,:) = dcumulu(:,:)
     ELSE             ; zmean(:,:) = dcumulu(:,:)/ntframe
     ENDIF
     IF ( lperio ) zmean(1, :) = zmean(npiglo-1,:)
     ierr = putvar(ncout, id_varout(2), zmean, jk,npiglo, npjglo, kwght=ntframe )

     ! < V >
     IF ( lopt ) THEN ; zmean(:,:) = dcumulv(:,:)
     ELSE             ; zmean(:,:) = dcumulv(:,:)/ntframe
     ENDIF
     IF ( lperio ) zmean(1, :) = zmean(npiglo-1,:)
     ierr = putvar(ncout, id_varout(3), zmean, jk,npiglo, npjglo, kwght=ntframe )

     ! < U'.V' >
     IF ( lopt ) THEN ; dupvp(:,:)= 1./(ntframe -1 ) * dupvp(:,:) ! use unbiased estimate 
     ELSE             ; dupvp(:,:) = dcumuluv(:,:)/ntframe - dcumulu(:,:)*dcumulv(:,:)/ntframe/ntframe
     ENDIF
     zmean(:,:) = dupvp(:,:)
     IF ( lperio ) zmean(1, :) = zmean(npiglo-1,:)
     ierr = putvar(ncout, id_varout(4), zmean, jk,npiglo, npjglo, kwght=ntframe )

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

    stypvar(1)%ichunk         = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname          = cn_vouv                  ; stypvar(1)%cunits        = 'm2/s2'
    stypvar(1)%clong_name     = 'U.V product at T point' ; stypvar(1)%cshort_name   = cn_vouv

    cl_name = TRIM(cn_vozocrtx)//'_t'
    stypvar(2)%ichunk         = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname          = cl_name                  ; stypvar(2)%cunits        = 'm/s'
    stypvar(2)%clong_name     = 'Mean U at T point '     ; stypvar(2)%cshort_name   = cl_name

    cl_name = TRIM(cn_vomecrty)//'_t'
    stypvar(3)%ichunk         = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(3)%cname          = cl_name                  ; stypvar(3)%cunits        = 'm/s'
    stypvar(3)%clong_name     = 'Mean V at T point '     ; stypvar(3)%cshort_name   = cl_name

    cl_name = TRIM(cn_vouv)//'_prime'
    stypvar(4)%ichunk         = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(4)%cname          = cl_name                     ; stypvar(3)%cunits        = 'm2/s2'
    stypvar(4)%clong_name     = 'Uprime .Vprime at T point' ; stypvar(3)%cshort_name   = cl_name

    ! create output fileset
    ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npkf, ld_xycoo=.TRUE., ld_nc4=lnc4 )
    ierr  = createvar   (ncout , stypvar, 4,      ipk,    id_varout           , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npkf, ld_xycoo=.TRUE. )

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
    PRINT *,' NTAGS', ntags
    ALLOCATE (ctag_lst(ntags) )
    DO ji = icur, icur + ntags -1
       CALL getarg(ji, ctag_lst( ji -icur +1 ) ) ; ijarg=ijarg+1
    END DO
  END SUBROUTINE GetTagList

END PROGRAM cdfuv
