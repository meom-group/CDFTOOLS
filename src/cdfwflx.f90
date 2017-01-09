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
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER                 :: jpvarout = 5     ! number of output variables
  INTEGER(KIND=4)                            :: jj, jk, ji       ! dummy loop index
  INTEGER(KIND=4)                            :: npiglo, npjglo   ! size of the domain
  INTEGER(KIND=4)                            :: narg, iargc      ! command line 
  INTEGER(KIND=4)                            :: ncout, ierr      ! netcdf i/o
  INTEGER(KIND=4), DIMENSION(jpvarout)       :: ipk, id_varout   ! levels and varid of output vars

  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: zmask, zwk       !  work array
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: evap, precip     ! water flux components
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: runoff, wdmp     ! water flux components
  REAL(KIND=4), DIMENSION(1)                 :: tim, dep         ! time_counter and dummy depth
  REAL(KIND=4)                               :: Lv=2.5e6         ! latent HF <--> evap conversion

  CHARACTER(LEN=256)                         :: cf_tfil          ! input gridT file name
  CHARACTER(LEN=256)                         :: cf_rnf           ! input runoff file name
  CHARACTER(LEN=256)                         :: cf_out='wflx.nc' ! output file

  TYPE(variable), DIMENSION(jpvarout)        :: stypvar          ! structure for attributes
  
  LOGICAL                                    :: lchk             ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfwflx T-file Runoff'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Computes the water fluxes components. Suitable for '
     PRINT *,'       annual means files. All output variables are in mm/days.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file  : model output file with water fluxes (gridT) '
     PRINT *,'       Runoff : file with the climatological runoff on the'
     PRINT *,'                model grid.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'       variables : soevap, soprecip, sorunoff, sowadmp, sowaflux'
     STOP
  ENDIF

  CALL getarg (1, cf_tfil)
  CALL getarg (2, cf_rnf )

  lchk = lchk .OR. chkfile ( cf_tfil)
  lchk = lchk .OR. chkfile ( cf_rnf )
  IF ( lchk ) STOP ! missing file

  npiglo= getdim (cf_tfil, cn_x)
  npjglo= getdim (cf_tfil, cn_y)

  ! prepare output variables
  dep(1) = 0.
  ipk(:) = 1  ! all variables ( output are 2D)

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

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo

  ALLOCATE ( zmask(npiglo,npjglo), zwk(npiglo,npjglo))
  ALLOCATE ( evap(npiglo,npjglo), precip(npiglo,npjglo), runoff(npiglo,npjglo), wdmp(npiglo,npjglo) )

 ! read vosaline for masking purpose
     zwk(:,:) =  getvar(cf_tfil, cn_vosaline,  1 ,npiglo,npjglo)
     zmask    = 1. ; WHERE ( zwk == 0 ) zmask = 0.
 ! Evap : 
     evap(:,:) = -1.* getvar(cf_tfil, cn_solhflup, 1 ,npiglo, npjglo)/Lv*86400. *zmask(:,:)  ! mm/days
      print *,'Evap done'
 ! Wdmp
     wdmp(:,:)   = getvar(cf_tfil, cn_sowafldp, 1 ,npiglo, npjglo) * 86400. * zmask(:,:)     ! mm/days
      print *,'Damping done'
 ! Runoff
     runoff(:,:) = getvar(cf_rnf,  'sorunoff',  1 ,npiglo, npjglo) * 86400. * zmask(:,:)     ! mm/days
      print *,'Runoff done'
 ! total water flux
     zwk(:,:)    = getvar(cf_tfil, cn_sowaflup, 1 ,npiglo, npjglo) * 86400. *zmask(:,:)      ! mm/days
      print *,'Total water flux done'
 ! Precip:
     precip(:,:)= evap(:,:) - runoff(:,:) + wdmp(:,:) - zwk(:,:)                            ! mm/day
      print *,'Precip done'

 ! Write output file
 ncout = create      (cf_out, cf_tfil, npiglo,   npjglo, 1          )
 ierr  = createvar   (ncout,  stypvar, jpvarout, ipk,    id_varout  )
 ierr  = putheadervar(ncout,  cf_tfil, npiglo,   npjglo, 1, pdep=dep)

 ierr = putvar(ncout, id_varout(1), evap,   1, npiglo, npjglo)
 ierr = putvar(ncout, id_varout(2), precip, 1, npiglo, npjglo)
 ierr = putvar(ncout, id_varout(3), runoff, 1, npiglo, npjglo)
 ierr = putvar(ncout, id_varout(4), wdmp,   1, npiglo, npjglo)
 ierr = putvar(ncout, id_varout(5), zwk,    1, npiglo, npjglo)

 tim   = getvar1d(cf_tfil, cn_vtimec, 1     )
 ierr  = putvar1d(ncout,   tim,       1, 'T')

 ierr=closeout(ncout)

END PROGRAM cdfwflx
