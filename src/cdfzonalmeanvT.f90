PROGRAM cdfzonalmeanvT
  !!======================================================================
  !!                     ***  PROGRAM  cdfzonalmeanvT  ***
  !!=====================================================================
  !!  ** Purpose : Compute the mean product of zonal mean V by zonal mean
  !!               of tracer (T and S )
  !!
  !!  ** Method  : In this program the 'zonal' meanvT is in fact a meanvT 
  !!               along the I coordinate. 
  !! History : 3.0  : 06/2013  : J.M. Molines : from cdfzonalmean and cdfvT
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class integration
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                              :: ji, jj, jk ,jt      ! dummy loop index
  INTEGER(KIND=4)                              :: jbasin, jtag        ! dummy loop index
  INTEGER(KIND=4)                              :: npbasins=1          ! number of subbasin
  INTEGER(KIND=4)                              :: ivar = 0            ! output variable counter
  INTEGER(KIND=4)                              :: narg, iargc         ! command line 
  INTEGER(KIND=4)                              :: ijarg               ! command line 
  INTEGER(KIND=4)                              :: itag, ntags         ! arg index of 1rst tag, number of tags
  INTEGER(KIND=4)                              :: ntframe             ! time frame counter
  INTEGER(KIND=4)                              :: npiglo, npjglo      ! size of the domain
  INTEGER(KIND=4)                              :: npk, npt            ! size of the domain
  INTEGER(KIND=4)                              :: ncout               ! ncid of output file
  INTEGER(KIND=4)                              :: ierr                ! working integers
  INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: ipk, id_varout      ! jpbasin x nvar
  INTEGER(KIND=4), DIMENSION(2)                :: ijloc               ! for maxloc

  REAL(KIND=4)                                 :: zspval=99999.       ! missing value 
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE :: gdep                ! gdept or gdepw
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: e1, e2, gphi        ! metrics, latitude
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: ztem, zsal, zvel    ! temp, sal and velocity
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zdumlon             ! dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zdumlat             ! latitude for i = north pole
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zvmask              ! vmask
  REAL(KIND=4), DIMENSION (:,:,:), ALLOCATABLE :: zmask               ! basin mask jpbasins x npiglo x npjglo

  REAL(KIND=8)                                 :: dtotal_time         ! cumul of time in seconds
  REAL(KIND=8), DIMENSION (:),     ALLOCATABLE :: dtim                ! time counter
  REAL(KIND=8), DIMENSION (:),     ALLOCATABLE :: dzovel, dzotem, dzosal  ! npjglo
  REAL(KIND=8), DIMENSION (:),     ALLOCATABLE :: darea, dl_tmp       ! npjglo
  REAL(KIND=8), DIMENSION (:,:,:,:), ALLOCATABLE :: dzovt, dzovs       ! 1xnpjglo x npk x npbasins

  CHARACTER(LEN=256)                           :: cf_tfil, cf_sfil    ! input files names for T S
  CHARACTER(LEN=256)                           :: cf_vfil             ! input files names for V
  CHARACTER(LEN=256)                           :: cf_out='zonalmeanvt.nc' ! output file name
  CHARACTER(LEN=256)                           :: cf_basins='none'    ! sub basin file name
  CHARACTER(LEN=256)                           :: cldum               ! dummy character variable
  CHARACTER(LEN=256)                           :: confcase            ! confcase name
  CHARACTER(LEN=4  ), DIMENSION(5)             :: cbasin=(/'_glo','_atl','_inp','_ind','_pac'/) ! sub basin suffixes
  CHARACTER(LEN=80 ),DIMENSION(:), ALLOCATABLE :: ctag_lst


  TYPE(variable), DIMENSION(:),    ALLOCATABLE :: stypvar             ! structure for input variables

  LOGICAL            :: lpdep    =.FALSE.   ! flag for depth sign (default dep < 0)
  LOGICAL            :: lndep_in =.FALSE.   ! flag for depth sign (default dep < 0) in input file
  LOGICAL            :: ldebug   =.FALSE.   ! flag for activated debug print 
  LOGICAL            :: lchk     =.FALSE.   ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfzonalmeanvT -c CONFIG-CASE -l LST-tags [-b BASIN-file] [-pdep] ...'
     PRINT *,'                 ...  [-ndep_in] [-o OUT-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the time average of zonal-mean(V) x zonal-mean(T/S) for the'
     PRINT *,'       set of files corresponding to the list of tags, passed as arguments.'
     PRINT *,'       This quantity is the mean-flow contribution to the heat/salt transport'
     PRINT *,'       overturning component. < > being the zonal average, we have: '
     PRINT *,'              Total       =        mean-flow             +     eddy.'
     PRINT *,'        time_mean(<V><T>) = time_mean(<V>)*time_mean(<T>)+time_mean(<V>''<T>'')'
     PRINT *,'      '
     PRINT *,'       Zonal mean is in fact the mean value computed along the I coordinate.'
     PRINT *,'       The result is a vertical slice, in the meridional direction.'
     PRINT *,'      '
     PRINT *,'       REMARKS:  Partial steps are not handled properly (but probably minor'
     PRINT *,'           impact on results) nor vvl case !.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -c CONFIG-CASE is the config name of a given experiment (eg ORCA025-G70)'
     PRINT *,'            The program will look for gridT, gridU and gridV files for this'
     PRINT *,'            config (grid_T, grid_U and grid_V are also accepted). In addition,'
     PRINT *,'            if gridS or grid_S file is found, it will be taken in place of '
     PRINT *,'            gridT for the salinity variable.'
     PRINT *,'       -l LST-tags : a blank-separated list of time tags that will be used'
     PRINT *,'            for time averaging. e.g. y2000m01d05 y2000m01d10 ...'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-b BASIN-file] : netcdf file describing sub basins, similar to '
     PRINT *,'            ', TRIM(cn_fbasins),'. If this name is not given as option, only'
     PRINT *,'            the global zonal mean is computed.'
     PRINT *,'       [-pdep  ] : use positive depths in the output file.'
     PRINT *,'            Default behaviour is to have negative depths.'
     PRINT *,'       [-ndep_in ] : negative depths are used in the input file.'
     PRINT *,'            Default behaviour is to have positive depths.'
     PRINT *,'       [-o OUT-file] : specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-debug] : add some print for debug'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fhgr),', ', TRIM(cn_fzgr),' and ', TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : zovzot : mean product of zonal_mean(V) x zonal_mean(T)'
     PRINT *,'                     zovzot : mean product of zonal_mean(V) x zonal_mean(S)'
     PRINT *,'                       A suffix _bas is append to variable name oin order to'
     PRINT *,'                     indicate the basin (atl, inp, ind, pac) or glo for global'
     PRINT *,'         '
     STOP 
  ENDIF

  ! decode command line
  ijarg = 1  
  DO WHILE ( ijarg <= narg ) 
     CALL getarg( ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE (cldum)
     CASE ( '-c'      ) ; CALL getarg( ijarg, confcase ) ; ijarg=ijarg+1  
     CASE ( '-l'      ) ; CALL GetTagList
     CASE ( '-pdep'   ) ; lpdep    =.TRUE.
     CASE ( '-ndep_in') ; lndep_in =.TRUE.
     CASE ( '-debug'  ) ; ldebug   =.TRUE.
     CASE ( '-b'      ) ; CALL getarg( ijarg, cf_basins) ; ijarg=ijarg+1 ; npbasins = 5
     CASE ( '-o'      ) ; CALL getarg( ijarg, cf_out   ) ; ijarg=ijarg+1 
     CASE DEFAULT       ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  IF ( ldebug ) THEN
     PRINT *, ' CONFIG-CASE = ', TRIM(confcase)
     PRINT *, ' NTAGS       = ', ntags
     PRINT *, '  arg pos    = ', itag
     PRINT *, '  -o         = ', TRIM(cf_out)
  ENDIF

  ! check  files existence
  lchk = lchk .OR. chkfile (cn_fhgr)
  lchk = lchk .OR. chkfile (cn_fzgr)
  lchk = lchk .OR. chkfile (cn_fmsk)
  IF ( npbasins /=1 ) THEN
     lchk = lchk .OR. chkfile (cf_basins  )
  ENDIF
  IF ( lchk ) STOP 99 ! missing files

  cf_tfil = SetFileName( confcase, cldum, 'T')  ! look in first T file for dimensions
  IF ( chkfile (cf_tfil) ) STOP 99 

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)

  ! Allocation ...
  ALLOCATE ( ipk(2*npbasins), id_varout(2*npbasins) )
  ALLOCATE ( stypvar(2*npbasins) )

  ALLOCATE ( gdep (npk) )
  ALLOCATE ( e1(npiglo,npjglo), e2(npiglo,npjglo), gphi(npiglo,npjglo) )
  ALLOCATE ( zvmask(npiglo, npjglo))
  ALLOCATE ( zvel(npiglo,npjglo), ztem(npiglo,npjglo), zsal(npiglo,npjglo) )
  ALLOCATE ( zdumlon(1,npjglo), zdumlat(1,npjglo) )
  ALLOCATE ( zmask(npbasins,npiglo, npjglo))

  ALLOCATE ( dzovel(npjglo), dzotem(npjglo), dzosal(npjglo), darea(npjglo), dl_tmp(npjglo) )
  ALLOCATE ( dzovt(1,npjglo, npk, npbasins), dzovs(1,npjglo, npk, npbasins) )

  ! read config information previous to file creation
  gdep(:)   = getvare3(cn_fzgr, cn_gdept, npk)
  gphi(:,:) = getvar(cn_fhgr, cn_gphit, 1, npiglo, npjglo)

  ! Look for the i-index that go through the North Pole
  ijloc        = MAXLOC(gphi)
  zdumlat(1,:) = gphi(ijloc(1),:)  ! 
  zdumlon(:,:) = 0.                ! set the dummy longitude to 0
  IF ( .NOT. lpdep ) gdep(:)   = -1.*  gdep(:)     ! helps for plotting the results

  IF ( ldebug ) PRINT *, 'Create Output files ...'
  CALL CreateOutput
  IF ( ldebug ) PRINT *, 'done.'

  ! initialization of 2D time independant fields
  zmask(1,:,:) = getvar(cn_fmsk, cn_tmask, 1, npiglo, npjglo)
  IF ( cf_basins /= 'none' ) THEN
     zmask(2,:,:) = getvar(cf_basins, cn_tmaskatl, 1, npiglo, npjglo )
     zmask(4,:,:) = getvar(cf_basins, cn_tmaskind, 1, npiglo, npjglo )
     zmask(5,:,:) = getvar(cf_basins, cn_tmaskpac, 1, npiglo, npjglo )
     zmask(3,:,:) = zmask(5,:,:) + zmask(4,:,:)
     ! ensure that there are no overlapping on the masks
     WHERE(zmask(3,:,:) > 0 ) zmask(3,:,:) = 1
  ENDIF
  e1(:,:)   = getvar(cn_fhgr, cn_ve1v,  1, npiglo, npjglo) 
  e2(:,:)   = getvar(cn_fhgr, cn_ve2v,  1, npiglo, npjglo) 

  ! tag loop
  ijarg = itag
  ntframe  = 0  ! reset count index for time mean

  DO jtag = 1, ntags 
     cldum = ctag_lst(jtag)
     cf_tfil = SetFileName( confcase, cldum, 'T', ld_stop=.TRUE.  )
     cf_sfil = SetFileName( confcase, cldum, 'S', ld_stop=.FALSE. )
     cf_vfil = SetFileName( confcase, cldum, 'V'                  )
     IF ( chkfile (cf_sfil, ld_verbose=.FALSE.) ) cf_sfil = cf_tfil  ! do not complain if not found
     IF ( ldebug ) THEN
        PRINT *, ' T-FILE = ', TRIM(cf_tfil)
        PRINT *, ' S-FILE = ', TRIM(cf_sfil)
        PRINT *, ' V-FILE = ', TRIM(cf_vfil)
     ENDIF

     npt = getdim (cf_tfil,cn_t)   ! case of multiple time frames in a single file, assume identical of V file
     ALLOCATE( dtim(npt) ) 
     dtim = getvar1d(cf_tfil, cn_vtimec, npt)
     IF ( ldebug ) PRINT *, 'TIME : ', dtim(:)
     dtotal_time = dtotal_time + SUM(dtim(1:npt) )
     DEALLOCATE( dtim )
     dzovt = 0.d0
     dzovs = 0.d0

     DO jt = 1, npt  
        ntframe = ntframe + 1
        DO jk = 1, npk 
           IF ( ldebug) PRINT *,' JTAG JT JK', jtag, jt, jk
           ! read variables 
           zsal(:,:) = getvar(cf_sfil,  cn_vosaline, jk, npiglo, npjglo, ktime=jt )
           ztem(:,:) = getvar(cf_tfil,  cn_votemper, jk, npiglo, npjglo, ktime=jt )
           zvel(:,:) = getvar(cf_vfil,  cn_vomecrty, jk, npiglo, npjglo, ktime=jt )
           ! do not read e3 metrics at level jk (to do as in cdfzonal mean ... JMM : to be improved !
           zvmask(:,:) = getvar(cn_fmsk, cn_vmask,    jk ,npiglo, npjglo          )

           ! put T and S at V points
           ztem(:,1:npjglo-1) = 0.5 * ( ztem(:,1:npjglo-1) + ztem(:,2:npjglo) ) * zvmask(:,1:npjglo-1)
           zsal(:,1:npjglo-1) = 0.5 * ( zsal(:,1:npjglo-1) + zsal(:,2:npjglo) ) * zvmask(:,1:npjglo-1)

           ! For all basins 
           DO jbasin = 1, npbasins
              IF ( ldebug) PRINT *,'    JBASIN ', jbasin
              dzovel(:) = 0.d0
              dzotem(:) = 0.d0
              dzosal(:) = 0.d0
              darea (:) = 0.d0
              IF ( ldebug) PRINT *,'      reset done.'

              ! integrates V 'zonally' (along i-coordinate)
              DO ji=1,npiglo
                 dl_tmp(:) = 1.d0*e1(ji,:)*e2(ji,:)* zmask(jbasin,ji,:) 
                 dzovel(:) = dzovel(:) + dl_tmp(:)*zvel  (ji,:)
                 dzotem(:) = dzotem(:) + dl_tmp(:)*ztem  (ji,:)
                 dzosal(:) = dzosal(:) + dl_tmp(:)*zsal  (ji,:)
                 darea (:) = darea (:) + dl_tmp(:)*zvmask(ji,:)
              END DO
              IF ( ldebug) PRINT *,'      Zonal integral done.'

              ! compute the mean value if the darea is not 0, else assign spval
              WHERE (darea /= 0 )
                 dzovel=dzovel/darea
                 dzotem=dzotem/darea
                 dzosal=dzosal/darea
                 ! cumulate in time
                 dzovt(1,:,jk,jbasin) = dzovt(1,:,jk,jbasin) + dzovel(:)*dzotem(:)
                 dzovs(1,:,jk,jbasin) = dzovs(1,:,jk,jbasin) + dzovel(:)*dzosal(:)
              ELSEWHERE
                 dzovel=zspval
                 dzotem=zspval
                 dzosal=zspval
                 ! cumulate in time
                 dzovt(1,:,jk,jbasin) = zspval
                 dzovs(1,:,jk,jbasin) = zspval
              ENDWHERE
              IF ( ldebug) PRINT *,'      mean and masking done.'

           END DO  !next basin
        END DO ! next level
     ENDDO   ! next time in file
  ENDDO  ! next file

  ! normalize before output
  WHERE ( dzovt /= zspval ) 
     dzovt(:,:,:,:) = dzovt(:,:,:,:) / ntframe
     dzovs(:,:,:,:) = dzovs(:,:,:,:) / ntframe
  ELSEWHERE
     dzovt(:,:,:,:) = zspval
     dzovs(:,:,:,:) = zspval
  ENDWHERE

  ALLOCATE ( dtim(1) )
  dtim(1) = dtotal_time/ntframe
  IF ( ldebug ) PRINT *, ' mean time ', dtim(1), ntframe

  ! output file 
  ierr   = putvar1d(ncout, dtim, 1, 'T')
  ivar = 0
  DO jbasin = 1, npbasins
     ivar = ivar + 1
     DO jk = 1, npk
        ierr = putvar (ncout, id_varout(ivar  ), REAL(dzovt(:,:,jk,jbasin) ), jk, 1, npjglo, kwght=ntframe )
        ierr = putvar (ncout, id_varout(ivar+1), REAL(dzovs(:,:,jk,jbasin) ), jk, 1, npjglo, kwght=ntframe )
     ENDDO
     ivar = ivar + 1
  ENDDO

  ierr = closeout(ncout)

CONTAINS 

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Set up all things required for the output file, create
    !!               the file and write the header part.
    !!
    !! ** Method  :  Use global module variables
    !!
    !!----------------------------------------------------------------------
    ipk(:) = npk
    ivar = 0
    DO jbasin = 1, npbasins
       ivar = ivar + 1
       stypvar(ivar)%cname             = 'zovzot'//TRIM(cbasin(jbasin) )
       stypvar(ivar)%cunits            = 'm.DegC.s-1'
       stypvar(ivar)%rmissing_value    = zspval
       stypvar(ivar)%valid_min         = -50.
       stypvar(ivar)%valid_max         =  50.
       stypvar(ivar)%clong_name        = 'product of zonalmean V x zonalmean T for'//TRIM(cbasin(jbasin) )
       stypvar(ivar)%cshort_name       = 'zovzot'//TRIM(cbasin(jbasin) )
       stypvar(ivar)%conline_operation = 'N/A'
       stypvar(ivar)%caxis             = 'TZY'

       ivar = ivar + 1
       stypvar(ivar)%cname             = 'zovzos'//TRIM(cbasin(jbasin) )
       stypvar(ivar)%cunits            = 'm.PSU.s-1'
       stypvar(ivar)%rmissing_value    = zspval
       stypvar(ivar)%valid_min         = -50.
       stypvar(ivar)%valid_max         = 50.
       stypvar(ivar)%clong_name        = 'product of zonalmean V x zonalmean S for'//TRIM(cbasin(jbasin) )
       stypvar(ivar)%cshort_name       = 'zovzos'//TRIM(cbasin(jbasin) )
       stypvar(ivar)%conline_operation = 'N/A'
       stypvar(ivar)%caxis             = 'TZY'
    END DO

    ! create output fileset
    ncout = create      (cf_out, cf_tfil,          1, npjglo, npk    )
    ierr  = createvar   (ncout,  stypvar, 2*npbasins, ipk, id_varout )
    ierr  = putheadervar(ncout,  cf_tfil,          1, npjglo, npk, pnavlon=zdumlon, pnavlat=zdumlat, pdep=gdep )

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

END PROGRAM cdfzonalmeanvT
