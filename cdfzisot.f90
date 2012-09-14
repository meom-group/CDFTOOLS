PROGRAM cdfzisot
  !!======================================================================
  !!                     ***  PROGRAM  cdfzisot  ***
  !!=====================================================================
  !!  ** Purpose : Compute isothermal depth
  !!
  !!  ** Method  : - compute surface properties
  !!               - initialize depths and model levels number
  !!
  !! History : 3.0  : 07/2012  : F.Hernandez: Original code
  !!           
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id: cdfzisot.f90  $
  !! Copyright (c) 2012, F. Hernandez
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4),PARAMETER                    :: pnvarout = 2   ! number of output variables
  INTEGER(KIND=4)                              :: ji, jj, jk, jt ! dummy loop index
  INTEGER(KIND=4)                              :: jref           ! dummy loop index
  INTEGER(KIND=4)                              :: narg, iargc, ii ! browse line
  INTEGER(KIND=4)                              :: npiglo, npjglo ! domain size
  INTEGER(KIND=4)                              :: npk, npt       ! domain size
  INTEGER(KIND=4)                              :: ncout, ierr    ! ncid of output file, error status
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: ipk, id_varout ! levels and varid's of output vars
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mbathy         ! mbathy metric

  REAL(KIND=4)                                 :: rtref          ! reference temperature
  REAL(KIND=4)                                 :: rmisval        ! Missing value of temperature
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: gdept          ! depth of T levels
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: gdepw          ! depth of W levels
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: tim            ! time counter
  REAL(KIND=4), DIMENSION(1)                   :: rdep           ! dummy depth for output
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rtem, rtemxz   ! temperature
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: tmask          ! temperature mask
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: glam,gphi      ! lon/lat
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rzisot         ! depth of the isotherm
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rzisotup       ! depth of the isotherm above
                                                                 ! in case of inversion

  CHARACTER(LEN=256)                           :: cf_tfil        ! input T file
  CHARACTER(LEN=256)                           :: cf_out='zisot.nc'! defaults output file name
  CHARACTER(LEN=256)                           :: cdum           ! dummy value

  TYPE(variable), DIMENSION(pnvarout)          :: stypvar        ! structure for output var. attributes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg < 2 ) THEN
     PRINT *,' usage : cdfzisot T-file RefTemp [Output File]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute depth of an isotherm given as argument'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file  : input netcdf file (gridT)' 
     PRINT *,'       RefTemp : Temperature of the isotherm.'
     PRINT *,'       Output File : netCDF Optional (defaults: zisot.nc)' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cn_fzgr)
     PRINT *,'         In case of FULL STEP configuration, ',TRIM(cn_fbathylev),' is also required.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     STOP
  ENDIF

  ! get input gridT filename
  CALL getarg (1, cf_tfil)
  IF ( chkfile(cf_tfil) .OR. chkfile(cn_fzgr) ) STOP ! missing file

  ! get reference temperature
  CALL getarg (2, cdum)
  IF (INDEX(cdum,'.') == 0 ) THEN
     READ(cdum,'(I10)') ii
     rtref = float(ii)
  ELSE
     READ(cdum,'(f12.6)') rtref
     cdum=cdum( 1:INDEX(cdum,'.')-1 )//'_'//cdum( INDEX(cdum,'.')+1:LEN_TRIM(cdum) )
  ENDIF
  IF ( rtref > 50 .OR. rtref < -3. ) THEN
     PRINT*,'Sea Water temperature on Earth, not Pluton ! ',rtref
     STOP
  ENDIF

  ! get (optional) output file name
  IF ( narg == 3 ) CALL getarg (3, cf_out)
  PRINT*,TRIM(cf_out)

  ! read dimensions 
  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)
!!$  PRINT *, 'npiglo = ', npiglo
!!$  PRINT *, 'npjglo = ', npjglo
!!$  PRINT *, 'npk    = ', npk
!!$  PRINT *, 'npt    = ', npt


  ! define structure to write the computed isotherm depth
  rdep(1) = 0.
  ipk(:)                    = 1
  stypvar(1)%cname          = 'zisot'//TRIM(cdum)
  stypvar(2)%cname          = 'zisotup'//TRIM(cdum)
  stypvar%cunits            = 'm'
  stypvar%rmissing_value    = 32767.
  stypvar%valid_min         = 0.
  stypvar%valid_max         = 7000.
  stypvar(1)%clong_name     = 'Depth_of_'//TRIM(cdum)//'C_isotherm' 
  stypvar(2)%clong_name     = 'Depth_of_'//TRIM(cdum)//'C_upper_isotherm' 
  stypvar(1)%cshort_name    = 'D'//TRIM(cdum)
  stypvar(2)%cshort_name    = 'D'//TRIM(cdum)//'up'
  stypvar%conline_operation = 'N/A'
  stypvar%caxis             = 'TYX'


  ! dynamical allocations
  ALLOCATE (gdept(npk), gdepw(npk), tim(npt) )
  ALLOCATE (rtem(npiglo,npjglo), rtemxz(npiglo,npk) )
  ALLOCATE (tmask(npiglo,npjglo), glam(npiglo,npjglo), gphi(npiglo,npjglo) )
  ALLOCATE (rzisot(npiglo,npjglo) , rzisotup(npiglo,npjglo) )
  ALLOCATE (mbathy(npiglo,npjglo) )

  ! read metrics gdept and gdepw
  gdept(:)    = getvare3(cn_fzgr, cn_gdept, npk )
  gdepw(:)    = getvare3(cn_fzgr, cn_gdepw, npk)
  ! read "mbathy"
  mbathy(:,:) = getvar(cn_fzgr,      'mbathy',    1, npiglo, npjglo)

  ! get missing value of votemper
  rmisval = getatt(cf_tfil, cn_votemper, cn_missing_value )

  ! get longitude and latitude
  glam(:,:) = getvar(cf_tfil, cn_vlon2d, 1, npiglo, npjglo)
  gphi(:,:) = getvar(cf_tfil, cn_vlat2d, 1, npiglo, npjglo)


  ! initialize tmask: 1=valid / 0= no valid of SST
  tmask(:,:) = 1.
  rtem( :,:) = getvar(cf_tfil, cn_votemper, 1, npiglo, npjglo, ktime=1 )  
  WHERE ( rtem == rmisval ) tmask = 0.
  ! initialise matrix of results
  rzisot = 0. ; WHERE ( rtem == rmisval ) rzisot = rmisval
  rzisotup = 0. ; WHERE ( rtem == rmisval ) rzisotup = rmisval

  ! Create output file, based on existing input gridT file
  ncout = create      (cf_out, cf_tfil, npiglo, npjglo, 1           )
  ierr  = createvar   (ncout,  stypvar, pnvarout,      ipk,    id_varout   )
  ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, 1, pdep=rdep)

  tim  = getvar1d(cf_tfil, cn_vtimec, npt     )
  ierr = putvar1d(ncout,   tim,       npt, 'T')


  ! compute depth for the isotherm, loop first on time
  DO jt=1,npt

     ! then loop on Y axis
     DO jj = 1 , npjglo

        ! read temperature on x-z slab
        rtemxz(:,:) = getvarxz(cf_tfil, cn_votemper, jj, npiglo, npk, kimin=1, kkmin=1, ktime=jt )

        ! loop on X axis
        DO ji = 1, npiglo

           ! proceed to the depth of isotherm computation if mask = OK
           IF ( tmask(ji,jj) == 1 ) THEN
              IF ( COUNT( rtemxz(ji,:)>=rtref .AND. rtemxz(ji,:) .NE.rmisval ) > 0 ) THEN

                 jk = 1 ! count level down

                 ! take into account temperature inversion from the surface
                 IF ( rtemxz(ji,1)<rtref .AND. rtemxz(ji,1) .NE.rmisval ) THEN

                    ! search first level with T >= rtref
                    DO WHILE ( jk < npk .AND. rtemxz(ji,jk) < rtref .AND. rtemxz(ji,jk) .NE. rmisval )
                       jref = jk
                       jk = jk + 1
                    ENDDO

                    ! compute depth of the above isotherm
                    rzisotup(ji,jj) = ( gdept(jk-1)*( rtemxz(ji,jk)-rtref ) + &
                         & gdept(jk)*( rtref-rtemxz(ji,jk-1) ) ) / &
                         & ( rtemxz(ji,jk)-rtemxz(ji,jk-1) )

                    !write(12,*)ji,jj,glam(ji,jj),gphi(ji,jj),rzisotup(ji,jj)

                 ENDIF

                 ! then start from the first level with T >= rtref
                 ! and search first value below rtref
                 jref = 0
                 DO WHILE ( jk < npk .AND. rtemxz(ji,jk) > rtref .AND. rtemxz(ji,jk) .NE. rmisval )
                    jref = jk
                    jk = jk + 1
                 ENDDO

                 ! test if the level is the last "wet level" in model metrics
                 ! OR if next temperature value is missing
                 ! Or at the bottom
                 ! --> give value of the bottom of the layer: gdepw(k+1)
                 IF ( jref == mbathy(ji,jj) .OR. rtemxz(ji,jref+1) == rmisval .OR. jref == npk-1 ) THEN
                    rzisot(ji,jj) = gdepw(jref+1)
                 ELSE
                    rzisot(ji,jj) = ( gdept(jref)*( rtemxz(ji,jref+1)-rtref ) + &
                         & gdept(jref+1)*( rtref-rtemxz(ji,jref) ) ) / &
                         & ( rtemxz(ji,jref+1)-rtemxz(ji,jref) )
                 ENDIF
 
              ENDIF  ! COUNT( rtemxz(ji,:)>=rtref
           ENDIF  ! tmask(ji,jj) == 1

        ENDDO ! ji = 1, npiglo
     ENDDO ! jj = 1 , npjglo

     ! Store the zisot variable in output file
     ierr = putvar(ncout, id_varout(1), rzisot, 1, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(2), rzisotup, 1, npiglo, npjglo, ktime=jt)

  END DO ! time loop

  ierr = closeout(ncout)

END PROGRAM cdfzisot
