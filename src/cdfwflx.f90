PROGRAM cdfwflx
  !!======================================================================
  !!                     ***  PROGRAM  cdfwflx  ***
  !!=====================================================================
  !!  ** Purpose : Produce a file with the water flux separated into 
  !!               4 components: E (soevap), P (soprecip), R (sorunoff),
  !!               dmp (sowafldp).
  !!               The total water flux is E -P -R + dmp. Units in this 
  !!               program are mm/days.
  !!
  !!  ** Method  : Evap is computed from the latent heat flux : evap=-qla/Lv
  !!               Runoff is read from the climatological input file
  !!               dmp is read from the file (sowafldp)
  !!               Precip is then computed as the difference between the
  !!               total water flux (sowaflup) and the E-R+dmp. In the high 
  !!               latitudes this precip includes the effect of snow 
  !!               (storage/melting). Therefore it may differ slightly from
  !!               the input precip file.
  !!
  !! History : 2.1  : 01/2008  : J.M. Molines : Original code
  !!           3.0  : 04/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class forcing
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER                 :: jpvarout = 5     ! number of output variables
  INTEGER(KIND=4)                            :: jj, jk, ji, jt   ! dummy loop index
  INTEGER(KIND=4)                            :: npiglo, npjglo   ! size of the domain
  INTEGER(KIND=4)                            :: npt              ! size of the domain
  INTEGER(KIND=4)                            :: narg, iargc      ! command line 
  INTEGER(KIND=4)                            :: ijarg            !
  INTEGER(KIND=4)                            :: ncout, ierr      ! netcdf i/o
  INTEGER(KIND=4), DIMENSION(jpvarout)       :: ipk, id_varout   ! levels and varid of output vars

  REAL(KIND=4)                               :: Lv=2.5e6         ! latent HF <--> evap conversion
  REAL(KIND=4), DIMENSION(1)                 :: dep              ! dummy depth
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: zmask, zwk       !  work array
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: evap, precip     ! water flux components
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: runoff, wdmp     ! water flux components

  REAL(KIND=8), DIMENSION(1)                 :: dtim             ! time_counter

  CHARACTER(LEN=256)                         :: cf_tfil          ! input gridT file name
  CHARACTER(LEN=256)                         :: cf_ffil          ! input flxT file name
  CHARACTER(LEN=256)                         :: cf_rnf           ! input runoff file name
  CHARACTER(LEN=256)                         :: cf_out='wflx.nc' ! output file
  CHARACTER(LEN=256)                         :: cldum            ! working char variable

  TYPE(variable), DIMENSION(jpvarout)        :: stypvar          ! structure for attributes

  LOGICAL                                    :: lchk             ! flag for missing files
  LOGICAL                                    :: lnc4 = .FALSE.   ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfwflx -t T-file -r RNF-file [-f FLX-file] [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the water fluxes components. Suitable for annual means files.'
     PRINT *,'       All output variables are in mm/days.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file   : model output file with water fluxes (gridT). '
     PRINT *,'       -r RNF-file : file with the climatological runoff on the model grid.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-f FLX-file]: model output file with water fluxes if not in T-file.'
     PRINT *,'       [-o OUT-file]: specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ] : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'            This option is effective only if cdftools are compiled with'
     PRINT *,'            a netcdf library supporting chunking and deflation.'
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option is used.'
     PRINT *,'       variables : soevap, soprecip, sorunoff, sowadmp, sowaflux'
     PRINT *,'      '
     STOP 
  ENDIF

  cf_ffil='none'
  ijarg=1
  DO WHILE ( ijarg <= narg) 
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum)
     CASE( '-t'  ) ; CALL getarg(ijarg, cf_tfil) ; ijarg=ijarg+1
     CASE( '-r'  ) ; CALL getarg(ijarg, cf_rnf ) ; ijarg=ijarg+1
        ! options
     CASE( '-f'  ) ; CALL getarg(ijarg, cf_ffil) ; ijarg=ijarg+1
     CASE( '-o'  ) ; CALL getarg(ijarg, cf_out ) ; ijarg=ijarg+1
     CASE( '-nc4') ; lnc4 = .TRUE.
     CASE DEFAULT  ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  IF (cf_ffil =='none' ) cf_ffil=cf_tfil

  lchk = lchk .OR. chkfile ( cf_tfil)
  lchk = lchk .OR. chkfile ( cf_rnf )
  IF ( lchk ) STOP 99 ! missing file

  npiglo= getdim (cf_tfil, cn_x)
  npjglo= getdim (cf_tfil, cn_y)
  npt   = getdim (cf_tfil, cn_t)

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npt   =', npt

  ALLOCATE ( zmask(npiglo,npjglo), zwk(npiglo,npjglo))
  ALLOCATE ( evap(npiglo,npjglo), precip(npiglo,npjglo), runoff(npiglo,npjglo), wdmp(npiglo,npjglo) )

  CALL CreateOutput

  DO jt =1, npt
     ! read vosaline for masking purpose
     zwk(:,:)    =  getvar(cf_tfil, cn_vosaline,  1 ,npiglo,npjglo, ktime=jt )
     zmask       =  1. ; WHERE ( zwk == 0 ) zmask = 0.
     evap(:,:)   = -1.* getvar(cf_ffil, cn_solhflup, 1 ,npiglo, npjglo, ktime=jt )/Lv*86400. *zmask(:,:)  ! mm/days
     wdmp(:,:)   =      getvar(cf_ffil, cn_sowafldp, 1 ,npiglo, npjglo, ktime=jt )   *86400. *zmask(:,:)  ! mm/days
     runoff(:,:) =      getvar(cf_rnf,  cn_sorunoff, 1 ,npiglo, npjglo, ktime=jt )   *86400. *zmask(:,:)  ! mm/days
     zwk(:,:)    =      getvar(cf_ffil, cn_sowaflup, 1 ,npiglo, npjglo, ktime=jt )   *86400. *zmask(:,:)  ! mm/days
     precip(:,:)= evap(:,:) - runoff(:,:) + wdmp(:,:) - zwk(:,:)                                          ! mm/day

     ierr = putvar(ncout, id_varout(1), evap,   1, npiglo, npjglo, ktime=jt )
     ierr = putvar(ncout, id_varout(2), precip, 1, npiglo, npjglo, ktime=jt )
     ierr = putvar(ncout, id_varout(3), runoff, 1, npiglo, npjglo, ktime=jt )
     ierr = putvar(ncout, id_varout(4), wdmp,   1, npiglo, npjglo, ktime=jt )
     ierr = putvar(ncout, id_varout(5), zwk,    1, npiglo, npjglo, ktime=jt )
  ENDDO

  ierr=closeout(ncout)

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
    INTEGER(KIND=4) :: jv
    !!----------------------------------------------------------------------
    ! prepare output variables
    dep(1) = 0.
    ipk(:) = 1  ! all variables ( output are 2D)
    DO jv=1, jpvarout
       stypvar(jv)%ichunk     = (/npiglo,MAX(1,npjglo/30),1,1 /)
    ENDDO

    stypvar%rmissing_value    = 0.
    stypvar%valid_min         = -100.
    stypvar%valid_max         =  100.
    stypvar%cunits            = 'mm/day'
    stypvar%conline_operation = 'N/A'
    stypvar%caxis             = 'TYX'

    stypvar(1)%cname = 'soevap'   ; stypvar(1)%clong_name = 'Evaporation'      ; stypvar(1)%cshort_name = 'soevap'
    stypvar(2)%cname = 'soprecip' ; stypvar(2)%clong_name = 'Precipitation'    ; stypvar(2)%cshort_name = 'soprecip'
    stypvar(3)%cname = 'sorunoff' ; stypvar(3)%clong_name = 'Runoff'           ; stypvar(3)%cshort_name = 'sorunoff'
    stypvar(4)%cname = 'sowadmp'  ; stypvar(4)%clong_name = 'SSS damping'      ; stypvar(4)%cshort_name = 'sowadmp'
    stypvar(5)%cname = 'sowaflux' ; stypvar(5)%clong_name = 'Total water flux' ; stypvar(5)%cshort_name = 'sowaflux'

    ! Write output file
    ncout = create      (cf_out, cf_tfil, npiglo,   npjglo, 1          , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar, jpvarout, ipk,    id_varout  , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_tfil, npiglo,   npjglo, 1, pdep=dep)


    dtim  = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr  = putvar1d(ncout,   dtim,      npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfwflx
