PROGRAM cdfvsig
  !!======================================================================
  !!                     ***  PROGRAM  cdfvsig  ***
  !!=====================================================================
  !!  ** Purpose : Compute the average values for the products 
  !!               U.sig, V.sig, W.sig where sig is the potential density.
  !!
  !!  ** Method  : pass the CONFIG name and a series of tags as arguments.
  !!               Tracers are interpolated on velocity points. The product
  !!               is evaluated at velocity points.
  !!
  !! History : 2.1  : 11/2004  : J.M. Molines : Original code
  !!           2.1  : 02/2010  : J.M. Molines : handle multiframes input files.
  !!           3.0  : 04/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jk, jt, jtt  ! dummy loop index
  INTEGER(KIND=4)                           :: ierr                 ! working integer
  INTEGER(KIND=4)                           :: narg, iargc          ! command line
  INTEGER(KIND=4)                           :: npiglo,npjglo        ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt             ! size of the domain
  INTEGER(KIND=4)                           :: ntframe              ! Cumul of time frame
  INTEGER(KIND=4)                           :: ncoutu               ! ncid of output file
  INTEGER(KIND=4)                           :: ncoutv               ! ncid of output file
  INTEGER(KIND=4)                           :: ncoutw               ! ncid of output file
  INTEGER(KIND=4), DIMENSION(3)             :: ipku, id_varoutu     ! level and varid's of output vars
  INTEGER(KIND=4), DIMENSION(3)             :: ipkv, id_varoutv     ! level and varid's of output vars
  INTEGER(KIND=4), DIMENSION(3)             :: ipkw, id_varoutw     ! level and varid's of output vars

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztemp,  zsal         ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zu, zv, zw           ! Velocity component
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztempu, zsalu        ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztempv, zsalv        ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztempw, zsalw        ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: umask,  vmask, wmask ! masks
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                  ! time counter of individual files
  REAL(KIND=4), DIMENSION(1)                :: timean               ! mean time

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulus             ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulvs             ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulws             ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulsu             ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulsv             ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulsw             ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulu              ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulv              ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulw              ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dsigu                ! Array for sigma0 at u point
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dsigv                ! Array for sigma0 at v point
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dsigw                ! Array for sigma0 at w point
  REAL(KIND=8)                              :: dtotal_time          ! cumulated time

  CHARACTER(LEN=256)                        :: cf_tfil              ! TS file name
  CHARACTER(LEN=256)                        :: cf_ufil              ! zonal velocity file
  CHARACTER(LEN=256)                        :: cf_vfil              ! meridional velocity file
  CHARACTER(LEN=256)                        :: cf_wfil              ! vertical velocity file
  CHARACTER(LEN=256)                        :: cf_outu='usig.nc'    ! output file
  CHARACTER(LEN=256)                        :: cf_outv='vsig.nc'    ! output file
  CHARACTER(LEN=256)                        :: cf_outw='wsig.nc'    ! output file
  CHARACTER(LEN=256)                        :: config               ! configuration name
  CHARACTER(LEN=256)                        :: ctag                 ! current tag to work with               

  TYPE (variable), DIMENSION(3)             :: stypvaru             ! structure for attributes
  TYPE (variable), DIMENSION(3)             :: stypvarv             ! structure for attributes
  TYPE (variable), DIMENSION(3)             :: stypvarw             ! structure for attributes

  LOGICAL                                   :: lcaltmean            ! flag for mean time computation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfvsig CONFIG ''list_of_tags'' '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the time average values for second order products ' 
     PRINT *,'       U.sig,  V.sig and W.sig.  Also save mean sigma-0 interpolated at'
     PRINT *,'       velocity points, as well as mean velocity component, for further use.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       CONFIG is the config name of a given experiment (eg ORCA025-G70)'
     PRINT *,'            The program will look for gridT, gridU, gridV  and gridW files for' 
     PRINT *,'            this config ( grid_T, grid_U, grid_V and grid_W are also accepted).'
     PRINT *,'       list_of_tags : a list of time tags that will be used for time'
     PRINT *,'            averaging. e.g. y2000m01d05 y2000m01d10 ...'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_outu),', ',TRIM(cf_outv),' and ', TRIM(cf_outw)
     PRINT *,'       variables : vousig, vovsig, vowsig : mean product v x sigma-0 '
     PRINT *,'                                            at velocity point.'
     PRINT *,'                   vosigu, vosigv, vosigw : mean sigma-0 at velocity point.'
     PRINT *,'                   ',TRIM(cn_vozocrtx),', ',TRIM(cn_vomecrty),', ',TRIM(cn_vovecrtz),' : mean velocity components.'
     STOP
  ENDIF

  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, config)
  CALL getarg (2, ctag  )

  cf_tfil = SetFileName ( config, ctag, 'T')
  cf_ufil = SetFileName ( config, ctag, 'U')
  cf_vfil = SetFileName ( config, ctag, 'V')
  cf_wfil = SetFileName ( config, ctag, 'W')

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)

  ipku(:)= npk  ! all variables (input and output are 3D)
  ipkv(:)= npk  !   "                     "
  ipkw(:)= npk  !   "                     "

  ! define output variables  U points
  stypvaru%rmissing_value    = 0.
  stypvaru%valid_min         = -100.
  stypvaru%valid_max         = 100.
  stypvaru%conline_operation = 'N/A'
  stypvaru%caxis             = 'TZYX'

  stypvaru(1)%cname          = 'vousig'           ; stypvaru(1)%cunits        = 'kg.m-2.s-1'
  stypvaru(2)%cname          = 'vosigu'           ; stypvaru(2)%cunits        = 'kg.m-3'
  stypvaru(3)%cname          =  cn_vozocrtx       ; stypvaru(3)%cunits        = 'm/s'

  stypvaru(1)%clong_name     = 'Mean U x sigma0'  ; stypvaru(1)%cshort_name   = 'vousig'
  stypvaru(2)%clong_name     = 'Mean sigma0 at U' ; stypvaru(2)%cshort_name   = 'vosigu'
  stypvaru(3)%clong_name     = 'Mean zonal vel'   ; stypvaru(3)%cshort_name   = cn_vozocrtx

  ! define output variables  V points
  stypvarv%rmissing_value    = 0.
  stypvarv%valid_min         = -100.
  stypvarv%valid_max         = 100.
  stypvarv%conline_operation = 'N/A'
  stypvarv%caxis             = 'TZYX'

  stypvarv(1)%cname          = 'vovsig'           ; stypvarv(1)%cunits        = 'kg.m-2.s-1'
  stypvarv(2)%cname          = 'vosigv'           ; stypvarv(2)%cunits        = 'kg.m-3'
  stypvarv(3)%cname          =  cn_vomecrty       ; stypvarv(3)%cunits        = 'm/s'

  stypvarv(1)%clong_name     = 'Mean V x sigma0'  ; stypvarv(1)%cshort_name   = 'vovsig'
  stypvarv(2)%clong_name     = 'Mean sigma0 at V' ; stypvarv(2)%cshort_name   = 'vosigv'
  stypvarv(3)%clong_name     = 'Mean merid vel'   ; stypvarv(3)%cshort_name   = cn_vomecrty

  ! define output variables  W points
  stypvarw%rmissing_value    = 0.
  stypvarw%valid_min         = -100.
  stypvarw%valid_max         = 100.
  stypvarw%conline_operation = 'N/A'
  stypvarw%caxis             = 'TZYX'

  stypvarw(1)%cname          = 'vowsig'           ; stypvarw(1)%cunits        = 'kg.m-2.s-1'
  stypvarw(2)%cname          = 'vosigw'           ; stypvarw(2)%cunits        = 'kg.m-3'
  stypvarw(3)%cname          =  cn_vovecrtz       ; stypvarw(3)%cunits        = 'm/s'

  stypvarw(1)%clong_name     = 'Mean W x sigma0'  ; stypvarw(1)%cshort_name   = 'vowsig'
  stypvarw(2)%clong_name     = 'Mean sigma0 at W' ; stypvarw(2)%cshort_name   = 'vosigw'
  stypvarw(3)%clong_name     = 'Mean vert. vel'   ; stypvarw(3)%cshort_name   = cn_vovecrtz

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk

  ALLOCATE( dcumulus(npiglo,npjglo), dcumulvs(npiglo,npjglo), dcumulws(npiglo,npjglo) )
  ALLOCATE( dcumulsu(npiglo,npjglo), dcumulsv(npiglo,npjglo), dcumulsw(npiglo,npjglo) )
  ALLOCATE( dcumulu(npiglo,npjglo),  dcumulv(npiglo,npjglo),  dcumulw(npiglo,npjglo)  )
  ALLOCATE( ztemp(npiglo,npjglo),    zsal(npiglo,npjglo)                              )
  ALLOCATE( zu(npiglo,npjglo),       zv(npiglo,npjglo),       zw(npiglo,npjglo)       )
  ALLOCATE( dsigu(npiglo,npjglo),    dsigv(npiglo,npjglo),    dsigw(npiglo,npjglo)    )
  ALLOCATE( umask(npiglo,npjglo),    vmask(npiglo,npjglo),    wmask(npiglo,npjglo)    )

  ! create output fileset
  ncoutu = create      (cf_outu, cf_ufil,  npiglo, npjglo, npk        )
  ierr   = createvar   (ncoutu,  stypvaru, 3,      ipku,   id_varoutu )
  ierr   = putheadervar(ncoutu,  cf_ufil,  npiglo, npjglo, npk        )

  ncoutv = create      (cf_outv, cf_vfil,  npiglo, npjglo, npk        )
  ierr   = createvar   (ncoutv,  stypvarv, 3,      ipkv,   id_varoutv )
  ierr   = putheadervar(ncoutv,  cf_vfil,  npiglo, npjglo, npk        )

  ncoutw = create      (cf_outw, cf_wfil,  npiglo, npjglo, npk        )
  ierr   = createvar   (ncoutw,  stypvarw, 3,      ipku,   id_varoutw )
  ierr   = putheadervar(ncoutw,  cf_wfil,  npiglo, npjglo, npk        )

  lcaltmean=.TRUE.
  DO jk = 1, npk
     PRINT *,'level ',jk
     dcumulus(:,:) = 0.d0 ;  dcumulvs(:,:) = 0.d0 ; dcumulws(:,:) = 0.d0
     dcumulsu(:,:) = 0.d0 ;  dcumulsv(:,:) = 0.d0 ; dcumulsw(:,:) = 0.d0
     dcumulu(:,:)  = 0.d0 ;  dcumulv(:,:)  = 0.d0 ; dcumulw(:,:)  = 0.d0
     dtotal_time   = 0.d0 ;  ntframe       = 0
     
     umask(:,:) = getvar(cn_fmsk, 'umask' , jk, npiglo, npjglo )
     vmask(:,:) = getvar(cn_fmsk, 'vmask' , jk, npiglo, npjglo )
     wmask(:,:) = getvar(cn_fmsk, 'tmask' , jk, npiglo, npjglo )

     DO jt = 2, narg           ! loop on tags
        CALL getarg (jt, ctag)
        cf_tfil = SetFileName ( config, ctag, 'T' )
        cf_ufil = SetFileName ( config, ctag, 'U' )
        cf_vfil = SetFileName ( config, ctag, 'V' )
        cf_wfil = SetFileName ( config, ctag, 'W' )

        npt = getdim (cf_tfil, cn_t)
        IF ( lcaltmean ) THEN
           ALLOCATE ( tim(npt) )
           tim = getvar1d(cf_tfil, cn_vtimec, npt)
           dtotal_time = dtotal_time + SUM(tim(1:npt) )
           DEALLOCATE( tim )
        END IF

        DO jtt = 1, npt  ! loop on time frame in a single file
           ntframe    = ntframe+1
           zu(:,:)    = getvar(cf_ufil, cn_vozocrtx, jk, npiglo, npjglo, ktime=jtt )
           zv(:,:)    = getvar(cf_vfil, cn_vomecrty, jk, npiglo, npjglo, ktime=jtt )
           zw(:,:)    = getvar(cf_wfil, cn_vovecrtz, jk, npiglo, npjglo, ktime=jtt )
           ztemp(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jtt )
           zsal(:,:)  = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jtt )

           ! temperature at u point, v points
           dsigu(:,:) = 0.d0  ; dsigv(:,:) = 0.d0
           DO ji=1, npiglo-1
              DO jj = 1, npjglo -1
                 ztempu(ji,jj) = 0.5 * ( ztemp(ji,jj) + ztemp(ji+1,jj) )  ! temper at Upoint
                 ztempv(ji,jj) = 0.5 * ( ztemp(ji,jj) + ztemp(ji,jj+1) )  ! temper at Vpoint
                 zsalu(ji,jj)  = 0.5 * ( zsal(ji,jj)  +  zsal(ji+1,jj) )  ! sal at U point
                 zsalv(ji,jj)  = 0.5 * ( zsal(ji,jj)  +  zsal(ji,jj+1) )  ! sal at v point
              END DO
           END DO
           
           dsigu(:,:) = sigma0(ztempu, zsalu, npiglo, npjglo) * umask(:,:)
           dsigv(:,:) = sigma0(ztempv, zsalv, npiglo, npjglo) * vmask(:,:)

           dcumulus(:,:) = dcumulus(:,:) + dsigu(:,:) * zu(:,:) * 1.d0
           dcumulvs(:,:) = dcumulvs(:,:) + dsigv(:,:) * zv(:,:) * 1.d0
           dcumulsu(:,:) = dcumulsu(:,:) + dsigu(:,:) * 1.d0
           dcumulsv(:,:) = dcumulsv(:,:) + dsigv(:,:) * 1.d0
           dcumulu(:,:)  = dcumulu(:,:)  + zu(:,:)    * 1.d0
           dcumulv(:,:)  = dcumulv(:,:)  + zv(:,:)    * 1.d0
          
           IF ( jk > 1 ) THEN ! now wsig
               ztempw(:,:)   = 0.5 * ( ztemp(:,:) + getvar(cf_tfil, cn_votemper, jk-1, npiglo, npjglo, ktime=jtt ))
               zsalw(:,:)    = 0.5 * ( zsal(:,:)  + getvar(cf_tfil, cn_vosaline, jk-1, npiglo, npjglo, ktime=jtt ))
               dsigw(:,:)    = sigma0(ztempw, zsalw, npiglo, npjglo) * wmask(:,:)
               dcumulws(:,:) = dcumulws(:,:) + dsigw(:,:) * zw(:,:) * 1.d0
               dcumulsw(:,:) = dcumulsw(:,:) + dsigw(:,:) * 1.d0
               dcumulw(:,:)  = dcumulw(:,:)  + zw(:,:)    * 1.d0
           ENDIF

        END DO  !jtt
     END DO  ! jt
     ! finish with level jk ; compute mean (assume spval is 0 )
     ierr = putvar(ncoutu, id_varoutu(1), SNGL(dcumulus(:,:)/ntframe), jk, npiglo, npjglo, kwght=ntframe )
     ierr = putvar(ncoutu, id_varoutu(2), SNGL(dcumulsu(:,:)/ntframe), jk, npiglo, npjglo, kwght=ntframe )
     ierr = putvar(ncoutu, id_varoutu(3), SNGL(dcumulu(:,:) /ntframe), jk, npiglo, npjglo, kwght=ntframe )

     ierr = putvar(ncoutv, id_varoutv(1), SNGL(dcumulvs(:,:)/ntframe), jk, npiglo, npjglo, kwght=ntframe )
     ierr = putvar(ncoutv, id_varoutv(2), SNGL(dcumulsv(:,:)/ntframe), jk, npiglo, npjglo, kwght=ntframe )
     ierr = putvar(ncoutv, id_varoutv(3), SNGL(dcumulv(:,:) /ntframe), jk, npiglo, npjglo, kwght=ntframe )

     ierr = putvar(ncoutw, id_varoutw(1), SNGL(dcumulws(:,:)/ntframe), jk, npiglo, npjglo, kwght=ntframe )
     ierr = putvar(ncoutw, id_varoutw(2), SNGL(dcumulsw(:,:)/ntframe), jk, npiglo, npjglo, kwght=ntframe )
     ierr = putvar(ncoutw, id_varoutw(3), SNGL(dcumulw(:,:) /ntframe), jk, npiglo, npjglo, kwght=ntframe )

     IF ( lcaltmean )  THEN
        timean(1) = dtotal_time/ntframe
        ierr      = putvar1d(ncoutu, timean, 1, 'T')
        ierr      = putvar1d(ncoutv, timean, 1, 'T')
        ierr      = putvar1d(ncoutw, timean, 1, 'T')
     END IF

     lcaltmean = .FALSE. ! tmean already computed

  END DO  ! loop to next level

  ierr = closeout(ncoutu)
  ierr = closeout(ncoutv)
  ierr = closeout(ncoutw)

END PROGRAM cdfvsig
