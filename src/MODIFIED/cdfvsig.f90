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
  !!         :  4.0  : 03/2017  : J.M. Molines  
  !!           2.1  : 02/2010  : J.M. Molines : handle multiframes input files.
  !!           3.0  : 04/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class second_order_moments
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jk, jt, jtt  ! dummy loop index
  INTEGER(KIND=4)                           :: jsig                 ! dummy loop index
  INTEGER(KIND=4)                           :: ierr                 ! working integer
  INTEGER(KIND=4)                           :: narg, iargc          ! command line
  INTEGER(KIND=4)                           :: ijarg, iiarg         ! argument counter
  INTEGER(KIND=4)                           :: ndep                 ! number of reference depth to deal with
  INTEGER(KIND=4)                           :: npiglo,npjglo        ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt             ! size of the domain
  INTEGER(KIND=4)                           :: ntframe              ! Cumul of time frame
  INTEGER(KIND=4)                           :: nopt                 ! number of options
  INTEGER(KIND=4)                           :: ncoutu               ! ncid of output file
  INTEGER(KIND=4)                           :: ncoutv               ! ncid of output file
  INTEGER(KIND=4)                           :: ncoutw               ! ncid of output file
  INTEGER(KIND=4)                           :: nfieldu, nfieldv, nfieldw  ! ncid of output file
  INTEGER(KIND=4)                           :: ivaru, ivarv, ivarw  ! variable counter
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipku, id_varoutu     ! level and varid's of output vars
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipkv, id_varoutv     ! level and varid's of output vars
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipkw, id_varoutw     ! level and varid's of output vars

  REAL(KIND=4)                              :: zdepref              ! reference level for potential density
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztemp,  zsal         ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zu, zv, zw           ! Velocity component
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztempu, zsalu        ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztempv, zsalv        ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztempw, zsalw        ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: umask,  vmask, wmask ! masks
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: refdep               ! Reference depth table
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                  ! time counter of individual files
  REAL(KIND=4), DIMENSION(1)                :: timean               ! mean time

  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dcumulus             ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dcumulvs             ! Arrays for cumulated values

  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dcumulws             ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dcumulsu             ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dcumulsv             ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dcumulsw             ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dsigu                ! Array for sigmai at u point
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dsigv                ! Array for sigmai at v point
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dsigw                ! Array for sigmai at w point
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dcumulu, dcumulu2    ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dcumulv, dcumulv2    ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dcumulw, dcumulw2    ! Arrays for cumulated values
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
  CHARACTER(LEN=256)                        :: cldum                ! dummy character var for browsing

  TYPE (variable), DIMENSION(:), ALLOCATABLE :: stypvaru             ! structure for attributes
  TYPE (variable), DIMENSION(:), ALLOCATABLE :: stypvarv             ! structure for attributes
  TYPE (variable), DIMENSION(:), ALLOCATABLE :: stypvarw             ! structure for attributes

  LOGICAL                                   :: lcaltmean            ! flag for mean time computation
  LOGICAL                                   :: lwo=.true.           ! flag for -no-w option
  LOGICAL                                   :: lsigo=.true.         ! flag for -no-sig option
  LOGICAL                                   :: luvo=.true.          ! flag for -no-uv option
  LOGICAL                                   :: lTpt=.false.         ! flag for -no-uv option
  LOGICAL                                   :: lpref=.false.        ! flag for -pref option
  LOGICAL                                   :: lperio=.false.       ! checking E-W periodicity

  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfvsig CONFIG  [-no-w] [-no-sig]  [-no-uv] [-T ] [-pref pref1,pref2,...]'
     PRINT *,'        ... ''list_of_tags'' '
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
     PRINT *,'            averaging. e.g. y2000m01d05 y2000m01d10.'
     PRINT *,'            ! IMPORTANT : list_of_tag are at the end of the command line ! '
     PRINT *,'      '
     PRINT *,'     OPTIONS ( to be used before the list_of tags ):'
     PRINT *,'        -T : compute u and v at T points, so that usig, vsig will be at T point'
     PRINT *,'        -no-w : no computation of vertical products'
     PRINT *,'        -no-sig : no output of density on U V points'
     PRINT *,'        -no-uv : no output of mean velocity components'
     PRINT *,'        -pref pref1,pref2,..: give comma separated list of reference depths for'
     PRINT *,'             density computation. eg : -pref 0,2000,3000  If not specified '
     PRINT *,'             assumes pref=0.'
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
  ijarg = 1 ; iiarg=0 ; nopt=0
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-no-w'   ) ; lwo=.false.   ; nopt=nopt+1
     CASE ( '-no-sig' ) ; lsigo=.false. ; nopt=nopt+1
     CASE ( '-no-uv'  ) ; luvo=.false.  ; nopt=nopt+1
     CASE ( '-T'      ) ; lTpt=.true.   ; nopt=nopt+1
     CASE ( '-pref'   ) ; lpref=.true.  ; CALL getarg(ijarg, cldum) ; ijarg=ijarg+1 
        CALL ParseRefDep(cldum)
        nopt=nopt+2 ! this option count for 2 args
     CASE DEFAULT
        iiarg=iiarg+1
        SELECT CASE (iiarg)
        CASE ( 1 )  ; config=cldum ;  nopt=nopt+1
        CASE ( 2 )  ; ctag=cldum
        END SELECT
     END SELECT
  END DO

  ! initialize refdep if not done on command line
  IF ( .NOT. lpref ) THEN
     ndep = 1
     ALLOCATE(refdep(ndep) )
     refdep(1) =0.0
  ENDIF

  ! |  always            |   if lsigo        |   if luvo    |
  ! | ndep [ u x sigma ] +   ndep [ sigma ]  +  2 [ u, u2 ] |
  nfieldu = ndep + ndep * COUNT ( (/lsigo /)  ) + 2 * COUNT ( (/luvo/) )
  nfieldv = nfieldu
  nfieldw = 2 * COUNT ( (/lwo /)) * nfieldu

  cf_tfil = SetFileName ( config, ctag, 'T')
  cf_ufil = SetFileName ( config, ctag, 'U')
  cf_vfil = SetFileName ( config, ctag, 'V')
  IF ( lwo ) cf_wfil = SetFileName ( config, ctag, 'W')

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk

  ALLOCATE( dcumulus(npiglo,npjglo,ndep), dcumulvs(npiglo,npjglo,ndep) )
  ALLOCATE( ztemp(npiglo,npjglo),    zsal(npiglo,npjglo)               )
  ALLOCATE( ztempu(npiglo,npjglo),  zsalu(npiglo,npjglo)               )
  ALLOCATE( ztempv(npiglo,npjglo),  zsalv(npiglo,npjglo)               )
  ALLOCATE( zu(npiglo,npjglo),       zv(npiglo,npjglo)                 )
  ALLOCATE( dsigu(npiglo,npjglo,ndep),    dsigv(npiglo,npjglo,ndep)    )
  ALLOCATE( umask(npiglo,npjglo),    vmask(npiglo,npjglo),    wmask(npiglo,npjglo)    )

  IF ( lwo   )           ALLOCATE( dcumulws(npiglo,npjglo,ndep) )
  IF ( lsigo )           ALLOCATE( dcumulsu(npiglo,npjglo,ndep), dcumulsv(npiglo,npjglo,ndep) )
  IF ( lwo .AND. lsigo ) ALLOCATE( dcumulsw(npiglo,npjglo,ndep) )
  IF ( luvo  )           ALLOCATE( dcumulu(npiglo,npjglo),  dcumulv(npiglo,npjglo)  )
  IF ( luvo  )           ALLOCATE( dcumulu2(npiglo,npjglo), dcumulv2(npiglo,npjglo) )
  IF ( lwo .AND. luvo  ) ALLOCATE( dcumulw(npiglo,npjglo), dcumulw2(npiglo,npjglo)  )

  IF ( lwo  )  ALLOCATE( zw(npiglo,npjglo)                             )
  IF ( lwo  )  ALLOCATE( ztempw(npiglo,npjglo),  zsalw(npiglo,npjglo)  )
  IF ( lwo  )  ALLOCATE( dsigw(npiglo,npjglo,ndep)                     )

  ! check periodicity (use ztemp as dummy array )
  ztemp(:,:) = getvar(cf_tfil, cn_vlon2d, 1, npiglo, npjglo )
  IF ( ztemp(1,1) == ztemp(npiglo-1,1) )  THEN
     lperio = .TRUE.
     PRINT *,' E-W periodicity detected '
  ENDIF


  CALL CreateOutputFile

  lcaltmean=.TRUE.
  DO jk = 1, npk
     PRINT *,'level ',jk
     dcumulus(:,:,:) = 0.d0 ;  dcumulvs(:,:,:) = 0.d0 
     IF ( lwo   ) dcumulws(:,:,:) = 0.d0
     IF ( lsigo ) THEN  ; dcumulsu(:,:,:) = 0.d0 ;  dcumulsv(:,:,:) = 0.d0  ;  ENDIF
     IF ( lsigo .AND. lwo ) dcumulsw(:,:,:) = 0.d0
     IF ( luvo  ) THEN   ; dcumulu(:,:)    = 0.d0  ;  dcumulv(:,:)   = 0.d0 ;  ENDIF
     IF ( luvo  ) THEN   ; dcumulu2(:,:)   = 0.d0  ;  dcumulv2(:,:)  = 0.d0 ;  ENDIF
     IF (luvo .AND. lwo ) THEN ; dcumulw(:,:) = 0.d0 ; dcumulw2(:,:) = 0.d0 ;  ENDIF
     dtotal_time   = 0.d0 ;  ntframe     = 0

     umask(:,:) = getvar(cn_fmsk, 'umask' , jk, npiglo, npjglo )
     vmask(:,:) = getvar(cn_fmsk, 'vmask' , jk, npiglo, npjglo )
     IF ( lwo .OR. lTpt ) wmask(:,:) = getvar(cn_fmsk, 'tmask' , jk, npiglo, npjglo )

     DO jt = nopt+1, narg            ! loop on tags
        CALL getarg (jt, ctag)
        cf_tfil = SetFileName ( config, ctag, 'T' )
        cf_ufil = SetFileName ( config, ctag, 'U' )
        cf_vfil = SetFileName ( config, ctag, 'V' )
        IF ( lwo ) cf_wfil = SetFileName ( config, ctag, 'W' )

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
           ztemp(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jtt )
           zsal(:,:)  = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jtt )
           IF ( lwo ) zw(:,:)    = getvar(cf_wfil, cn_vovecrtz, jk, npiglo, npjglo, ktime=jtt )

           dsigu(:,:,:) = 0.d0  ; dsigv(:,:,:) = 0.d0
           IF ( lTpt )  THEN
              ! u,v at T point
              DO ji=npiglo,2,-1
                 zu(ji,:) = 0.5 * ( zu(ji-1,:) + zu(ji,:) )
              END DO
              DO jj=npjglo,2,-1
                 zv(:,jj) = 0.5 * ( zv(:,jj-1) + zv(:,jj) )
              END DO
              IF ( lperio ) THEN
                 zu(1,:) = zu(npiglo-1,:)
              ENDIF
              IF ( luvo ) THEN
                 dcumulu(:,:)   = dcumulu(:,:)  + zu(:,:)    * 1.d0
                 dcumulv(:,:)   = dcumulv(:,:)  + zv(:,:)    * 1.d0
                 dcumulu2(:,:)  = dcumulu2(:,:)  + zu(:,:)*zu(:,:) * 1.d0
                 dcumulv2(:,:)  = dcumulv2(:,:)  + zv(:,:)*zv(:,:) * 1.d0
              ENDIF
              DO jsig=1,ndep
                 zdepref=refdep(jsig)
                 ! rem : use dsigu for sig at T point, for the sake of simplicity ...
                 dsigu(:,:,jsig) = sigmai(ztemp, zsal, zdepref, npiglo, npjglo) * wmask(:,:)

                 dcumulus(:,:,jsig) = dcumulus(:,:,jsig) + dsigu(:,:,jsig) * zu(:,:) * 1.d0
                 dcumulvs(:,:,jsig) = dcumulvs(:,:,jsig) + dsigu(:,:,jsig) * zv(:,:) * 1.d0
                 IF ( lsigo ) THEN
                    dcumulsu(:,:,jsig) = dcumulsu(:,:,jsig) + dsigu(:,:,jsig) * 1.d0
                 ENDIF
              ENDDO
           ELSE
              ! temperature at u point, v points
              DO ji=1, npiglo-1
                 DO jj = 1, npjglo -1
                    ztempu(ji,jj) = 0.5 * ( ztemp(ji,jj) + ztemp(ji+1,jj) )  ! temper at Upoint
                    ztempv(ji,jj) = 0.5 * ( ztemp(ji,jj) + ztemp(ji,jj+1) )  ! temper at Vpoint
                    zsalu(ji,jj)  = 0.5 * ( zsal(ji,jj)  +  zsal(ji+1,jj) )  ! sal at U point
                    zsalv(ji,jj)  = 0.5 * ( zsal(ji,jj)  +  zsal(ji,jj+1) )  ! sal at v point
                 END DO
              END DO

              IF ( luvo ) THEN
                 dcumulu(:,:)   = dcumulu(:,:)  + zu(:,:)    * 1.d0
                 dcumulv(:,:)   = dcumulv(:,:)  + zv(:,:)    * 1.d0
                 dcumulu2(:,:)  = dcumulu2(:,:)  + zu(:,:)*zu(:,:) * 1.d0
                 dcumulv2(:,:)  = dcumulv2(:,:)  + zv(:,:)*zv(:,:) * 1.d0
              ENDIF
              DO jsig=1,ndep
                 zdepref=refdep(jsig)
                 dsigu(:,:,jsig) = sigmai(ztempu, zsalu, zdepref, npiglo, npjglo) * umask(:,:)
                 dsigv(:,:,jsig) = sigmai(ztempv, zsalv, zdepref, npiglo, npjglo) * vmask(:,:)

                 dcumulus(:,:,jsig) = dcumulus(:,:,jsig) + dsigu(:,:,jsig) * zu(:,:) * 1.d0
                 dcumulvs(:,:,jsig) = dcumulvs(:,:,jsig) + dsigv(:,:,jsig) * zv(:,:) * 1.d0
                 IF ( lsigo ) THEN
                    dcumulsu(:,:,jsig) = dcumulsu(:,:,jsig) + dsigu(:,:,jsig) * 1.d0
                    dcumulsv(:,:,jsig) = dcumulsv(:,:,jsig) + dsigv(:,:,jsig) * 1.d0
                 ENDIF

                 IF ( lwo ) THEN
                    IF ( jk > 1 ) THEN ! now wsig
                       ztempw(:,:)   = 0.5 * ( ztemp(:,:) + getvar(cf_tfil, cn_votemper, jk-1, npiglo, npjglo, ktime=jtt ))
                       zsalw(:,:)    = 0.5 * ( zsal(:,:)  + getvar(cf_tfil, cn_vosaline, jk-1, npiglo, npjglo, ktime=jtt ))
                       dsigw(:,:,jsig)    = sigmai(ztempw, zsalw, zdepref, npiglo, npjglo) * wmask(:,:)
                       dcumulws(:,:,jsig) = dcumulws(:,:,jsig) + dsigw(:,:,jsig) * zw(:,:) * 1.d0
                       IF ( lsigo ) dcumulsw(:,:,jsig) = dcumulsw(:,:,jsig) + dsigw(:,:,jsig) * 1.d0
                       IF ( luvo  ) dcumulw(:,:)   = dcumulw(:,:)   + zw(:,:)         * 1.d0
                       IF ( luvo  ) dcumulw2(:,:)  = dcumulw2(:,:)  + zw(:,:)*zw(:,:) * 1.d0
                    ENDIF
                 ENDIF
              ENDDO
           ENDIF  ! Tpoint

        END DO  !jtt
     END DO  ! jt
     ! finish with level jk ; compute mean (assume spval is 0 )
     ivaru=0 ; ivarv=0 ; ivarw=0
     !    U
     DO jsig=1,ndep
        ivaru=ivaru + 1
        ierr = putvar(ncoutu, id_varoutu(ivaru), SNGL(dcumulus(:,:,jsig)/ntframe), jk, npiglo, npjglo, kwght=ntframe )
        IF ( lsigo) THEN 
           ivaru=ivaru + 1
           ierr = putvar(ncoutu, id_varoutu(ivaru), SNGL(dcumulsu(:,:,jsig)/ntframe), jk, npiglo, npjglo, kwght=ntframe )
        ENDIF
     END DO

     IF ( luvo ) THEN 
        ivaru=ivaru + 1
        ierr = putvar(ncoutu, id_varoutu(ivaru), SNGL(dcumulu(:,:)      /ntframe), jk, npiglo, npjglo, kwght=ntframe )
        ivaru=ivaru + 1
        ierr = putvar(ncoutu, id_varoutu(ivaru), SNGL(dcumulu2(:,:)     /ntframe), jk, npiglo, npjglo, kwght=ntframe )
     ENDIF

     !    V
     DO jsig=1,ndep
        ivarv=ivarv + 1
        ierr = putvar(ncoutv, id_varoutv(ivarv), SNGL(dcumulvs(:,:,jsig)/ntframe), jk, npiglo, npjglo, kwght=ntframe )
        IF ( lsigo) THEN 
           ivarv=ivarv + 1
           ierr = putvar(ncoutv, id_varoutv(ivarv), SNGL(dcumulsv(:,:,jsig)/ntframe), jk, npiglo, npjglo, kwght=ntframe )
        ENDIF
     ENDDO
     IF ( luvo ) THEN 
        ivarv=ivarv + 1
        ierr = putvar(ncoutv, id_varoutv(ivarv), SNGL(dcumulv(:,:)      /ntframe), jk, npiglo, npjglo, kwght=ntframe )
        ivarv=ivarv + 1
        ierr = putvar(ncoutv, id_varoutv(ivarv), SNGL(dcumulv2(:,:)     /ntframe), jk, npiglo, npjglo, kwght=ntframe )
     ENDIF

     !    W
     IF ( lwo ) THEN
        DO jsig=1,ndep
           ivarw=ivarw + 1
           ierr = putvar(ncoutw, id_varoutw(ivarw), SNGL(dcumulws(:,:,jsig)/ntframe), jk, npiglo, npjglo, kwght=ntframe )
           IF (lsigo) THEN 
              ivarw=ivarw + 1
              ierr = putvar(ncoutw, id_varoutw(ivarw), SNGL(dcumulsw(:,:,jsig)/ntframe), jk, npiglo, npjglo, kwght=ntframe )
           ENDIF
        ENDDO
        IF ( luvo ) THEN 
           ivarw=ivarw + 1
           ierr = putvar(ncoutw, id_varoutw(ivarw), SNGL(dcumulw(:,:)      /ntframe), jk, npiglo, npjglo, kwght=ntframe )
           ivarw=ivarw + 1
           ierr = putvar(ncoutw, id_varoutw(ivarw), SNGL(dcumulw2(:,:)     /ntframe), jk, npiglo, npjglo, kwght=ntframe )
        ENDIF
     ENDIF

     IF ( lcaltmean )  THEN
        timean(1) = dtotal_time/ntframe
        ierr      = putvar1d(ncoutu, timean, 1, 'T')
        ierr      = putvar1d(ncoutv, timean, 1, 'T')
        IF ( lwo ) ierr = putvar1d(ncoutw, timean, 1, 'T')
     END IF

     lcaltmean = .FALSE. ! tmean already computed

  END DO  ! loop to next level

  ierr = closeout(ncoutu)
  ierr = closeout(ncoutv)
  IF ( lwo ) ierr = closeout(ncoutw)

CONTAINS

  SUBROUTINE CreateOutputFile
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutputFile  ***
    !!
    !! ** Purpose :  Create netcdf file according to flags
    !!
    !! ** Method  :  Flags are global variables and there dore known in
    !!               in this routine.
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4)    :: jsig, ivaru, ivarv, ivarw
    CHARACTER(LEN=1)   :: cldep
    CHARACTER(LEN=256) :: cl_global, cl_refu, cl_refv, cl_refw
    !!----------------------------------------------------------------------
    ALLOCATE ( stypvaru(nfieldu), ipku(nfieldu), id_varoutu(nfieldu)    )
    ALLOCATE ( stypvarv(nfieldv), ipkv(nfieldv), id_varoutv(nfieldv)     )
    IF ( lwo ) ALLOCATE ( stypvarw(nfieldw), ipkw(nfieldw), id_varoutw(nfieldw)  )

    IF ( lTpt ) THEN
       cl_global = ' All variables computed on T points'
       cl_refu= cf_tfil
       cl_refv= cf_tfil
       cl_refw= cf_tfil
    ELSE
       cl_global = ' All variables on native U or V points'
       cl_refu= cf_ufil
       cl_refv= cf_vfil
       cl_refw= cf_wfil
    ENDIF


    ! define output variables  U points
    stypvaru%rmissing_value    = 0.
    stypvaru%valid_min         = -100.
    stypvaru%valid_max         = 100.
    stypvaru%conline_operation = 'N/A'
    stypvaru%caxis             = 'TZYX'
    ipku(:)= npk  ! all variables (input and output are 3D)

    ivaru=0
    DO jsig = 1, ndep
       WRITE(cldep,'(I1)') INT(refdep(jsig)/1000)
       ivaru = ivaru + 1
       stypvaru(ivaru)%cname      = 'vousig'//cldep          ; stypvaru(ivaru)%cunits    = 'kg.m-2.s-1'
       stypvaru(ivaru)%clong_name = 'Mean U x sigma'//cldep  ; stypvaru(ivaru)%cshort_name   = 'vousig'//cldep

       IF ( lsigo ) THEN
          ivaru = ivaru + 1 
          stypvaru(ivaru)%cname      = 'vosigu'//cldep              ; stypvaru(ivaru)%cunits  = 'kg.m-3'
          stypvaru(ivaru)%clong_name = 'Mean sigma'//cldep//' at U' ; stypvaru(ivaru)%cshort_name = 'vosigu'//cldep
       ENDIF
    ENDDO

    IF ( luvo ) THEN
       ivaru = ivaru + 1 
       stypvaru(ivaru)%cname      = cn_vozocrtx      ; stypvaru(ivaru)%cunits        = 'm/s'
       stypvaru(ivaru)%clong_name = 'Mean zonal vel' ; stypvaru(ivaru)%cshort_name   = cn_vozocrtx
       ivaru = ivaru + 1 
       stypvaru(ivaru)%cname      = TRIM(cn_vozocrtx)//'_sqd'      ; stypvaru(ivaru)%cunits        = '(m/s)^2'
       stypvaru(ivaru)%clong_name = 'Mean zonal vel squared' ; stypvaru(ivaru)%cshort_name   = TRIM(cn_vozocrtx)//'_sqd' 
    ENDIF

    ! create output fileset
    ncoutu = create      (cf_outu, cl_refu,  npiglo, npjglo, npk        )
    ierr   = createvar   (ncoutu,  stypvaru, nfieldu,      ipku,   id_varoutu , cdglobal=cl_global)
    ierr   = putheadervar(ncoutu,  cl_refu,  npiglo, npjglo, npk        )


    ! define output variables  V points
    stypvarv%rmissing_value    = 0.
    stypvarv%valid_min         = -100.
    stypvarv%valid_max         = 100.
    stypvarv%conline_operation = 'N/A'
    stypvarv%caxis             = 'TZYX'
    ipkv(:)= npk  !   "                     "

    ivarv=0
    DO jsig = 1, ndep 
       WRITE(cldep,'(I1)') INT(refdep(jsig)/1000)
       ivarv = ivarv + 1
       stypvarv(ivarv)%cname      = 'vovsig'//cldep          ; stypvarv(ivarv)%cunits      = 'kg.m-2.s-1'
       stypvarv(ivarv)%clong_name = 'Mean V x sigma'//cldep  ; stypvarv(ivarv)%cshort_name = 'vovsig'//cldep

       IF ( lsigo ) THEN
          ivarv = ivarv + 1
          stypvarv(ivarv)%cname      = 'vosigv'//cldep              ; stypvarv(ivarv)%cunits      = 'kg.m-3'
          stypvarv(ivarv)%clong_name = 'Mean sigma'//cldep//' at V' ; stypvarv(ivarv)%cshort_name = 'vosigv'//cldep
       ENDIF
    ENDDO

    IF ( luvo ) THEN
       ivarv = ivarv + 1 
       stypvarv(ivarv)%cname      = cn_vomecrty      ; stypvarv(ivarv)%cunits        = 'm/s'
       stypvarv(ivarv)%clong_name = 'Mean merid vel' ; stypvarv(ivarv)%cshort_name   = cn_vomecrty
       ivarv = ivarv + 1 
       stypvarv(ivarv)%cname      = TRIM(cn_vomecrty)//'_sqd'      ; stypvarv(ivarv)%cunits  = '(m/s)^2'
       stypvarv(ivarv)%clong_name = 'Mean merid vel squared' ; stypvarv(ivarv)%cshort_name   = TRIM(cn_vomecrty)//'_sqd'
    ENDIF

    ncoutv = create      (cf_outv, cl_refv,  npiglo, npjglo, npk        )
    ierr   = createvar   (ncoutv,  stypvarv, nfieldv,      ipkv,   id_varoutv, cdglobal=cl_global )
    ierr   = putheadervar(ncoutv,  cl_refv,  npiglo, npjglo, npk        )

    IF ( lwo ) THEN
       ! define output variables  W points
       stypvarw%rmissing_value    = 0.
       stypvarw%valid_min         = -100.
       stypvarw%valid_max         = 100.
       stypvarw%conline_operation = 'N/A'
       stypvarw%caxis             = 'TZYX'
       ipkw(:)= npk  !   "                     "

       ivarw=0
       DO jsig = 1, ndep
          WRITE(cldep,'(I1)') INT(refdep(jsig)/1000)
          ivarw = ivarw + 1
          stypvarw(ivarw)%cname      = 'vowsig'//cldep          ; stypvarw(ivarw)%cunits      = 'kg.m-2.s-1'
          stypvarw(ivarw)%clong_name = 'Mean W x sigma'//cldep  ; stypvarw(ivarw)%cshort_name = 'vowsig'//cldep

          IF ( lsigo ) THEN
             ivarw = ivarw + 1
             stypvarw(ivarw)%cname      = 'vosigw'//cldep              ; stypvarw(ivarw)%cunits      = 'kg.m-3'
             stypvarw(ivarw)%clong_name = 'Mean sigma'//cldep//' at W' ; stypvarw(ivarw)%cshort_name = 'vosigw'//cldep
          ENDIF
       ENDDO

       IF ( luvo ) THEN
          ivarw = ivarw + 1
          stypvarw(ivarw)%cname      = cn_vovecrtz      ; stypvarw(ivarw)%cunits        = 'm/s'
          stypvarw(ivarw)%clong_name = 'Mean vert. vel' ; stypvarw(ivarw)%cshort_name   = cn_vovecrtz
          ivarw = ivarw + 1
          stypvarw(ivarw)%cname      = TRIM(cn_vovecrtz)//'_sqd' ; stypvarw(ivarw)%cunits       = '(m/s)^2'
          stypvarw(ivarw)%clong_name = 'Mean vert. vel squared' ; stypvarw(ivarw)%cshort_name   = TRIM(cn_vovecrtz)//'_sqd'
       ENDIF

       ncoutw = create      (cf_outw, cl_refw,  npiglo, npjglo, npk        )
       ierr   = createvar   (ncoutw,  stypvarw, nfieldw,      ipku,   id_varoutw, cdglobal=cl_global )
       ierr   = putheadervar(ncoutw,  cl_refw,  npiglo, npjglo, npk        )
    ENDIF


  END SUBROUTINE  CreateOutputFile

  SUBROUTINE ParseRefDep (cdum)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ParseRefDep  ***
    !!
    !! ** Purpose :  Decode variable name  option from command line
    !!
    !! ** Method  :  look for , in the argument string and set the number of
    !!         ref_dep (ndep), allocate  array and fill it with the
    !!         decoded  names.
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdum

    CHARACTER(LEN=80), DIMENSION(100) :: cl_dum  ! 100 is arbitrary
    INTEGER  :: ji
    INTEGER  :: inchar,  i1=1
    !!----------------------------------------------------------------------
    inchar= LEN(TRIM(cdum))
    ndep = 1
    ! scan the input string and look for ',' as separator
    DO ji=1,inchar
       IF ( cdum(ji:ji) == ',' ) THEN
          cl_dum(ndep) = cdum(i1:ji-1)
          i1=ji+1
          ndep=ndep+1
       ENDIF
    ENDDO

    ! last name of the list does not have a ','
    cl_dum(ndep) = cdum(i1:inchar)

    ALLOCATE ( refdep(ndep) )
    DO ji=1, ndep
       READ(cl_dum(ji),*) refdep(ji)
    ENDDO
    PRINT *, refdep
  END SUBROUTINE ParseRefDep


END PROGRAM cdfvsig
