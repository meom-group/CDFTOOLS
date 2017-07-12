PROGRAM cdfmaskdmp
  !!======================================================================
  !!                     ***  PROGRAM  cdfmaskdmp  ***
  !!=====================================================================
  !!  ** Purpose : Compute 3D mask for AABW relaxation from T and S 
  !!               climatology.
  !!               Store the results on a cdf file.
  !!
  !!  **  Method: read temp and salinity, compute sigma-2
  !!              compute coefs, create mask
  !!
  !! History : 2.1  : 09/2010  : R. Dussin    : Original code from JLS Py version
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class preprocessing
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jk, jt     ! dummy loop index
  INTEGER(KIND=4)                           :: ierr               ! error status
  INTEGER(KIND=4)                           :: narg, iargc,ijarg  ! browse command line
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                           :: ncout              ! ncid of output file
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout     ! level and varid's

  REAL(KIND=4)                              :: ref_dep=2000.      ! reference depth in meters
  REAL(KIND=4)                              :: zsnmin=37.16       ! minimum density
  REAL(KIND=4)                              :: zswidth=0.025      ! tapering width
  REAL(KIND=4)                              :: hmin=1000.         ! depth limit
  REAL(KIND=4)                              :: hwidth=100.        ! depth tapering height
  REAL(KIND=4)                              :: rlatmax=-20.       ! max latitude
  REAL(KIND=4)                              :: rlatwidth=2.       ! latitude tapering width
  REAL(KIND=4)                              :: wdep, wsig, wlat   ! tapering function dep, sigma and lat
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztemp              ! temperature
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsal               ! salinity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsigi              ! sigma-i
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask              ! 2D mask at current level
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zwdmp              ! 2D build mask at current level
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zlat               ! latitudes
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: zdep               ! deptht

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim               ! time counter

  CHARACTER(LEN=256)                        :: cf_tfil            ! input filename for temperature
  CHARACTER(LEN=256)                        :: cf_sfil ='none'    ! input filename for salinity
  CHARACTER(LEN=256)                        :: cf_out='mask_dmp.nc' ! output file name
  CHARACTER(LEN=256)                        :: cldum              ! dummy string
  CHARACTER(LEN=256)                        :: cglobal            ! Global attribute with command name

  TYPE (variable), DIMENSION(1)             :: stypvar            ! structure for attributes

  LOGICAL                                   :: lnc4 = .FALSE.     ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmaskdmp -t T-file [-s S-file] [-refdep REF-depth] ...'
     PRINT *,'             ... [-dens smin width] [-dep hmin width] [-lat latmax width] ...'
     PRINT *,'             ... [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute a damping mask with smooth transition according to density,'
     PRINT *,'       depth and latitude criteria. For each kind of criterion, a minimum  '
     PRINT *,'       value and a width of transition is used. Both minimum and width can'
     PRINT *,'       be adjusted with corresponding options. Reasonable default values are'
     PRINT *,'       provided.'
     PRINT *,'      '
     PRINT *,'       This tool was designed for building a mask in the deep Southern Ocean,'
     PRINT *,'       for dense waters. This explains that we consider the limit as a minimum'
     PRINT *,'       depth and density, but a maximum latitude. It needs some adjustments for'
     PRINT *,'       other situations.'
     PRINT *,'       '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file : temperature/salinity file used to compute the potential ' 
     PRINT *,'             density relative to the reference depth.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-s S-file] : salinity file in case it differs from the temperature file.' 
     PRINT *,'       [-refdep REF-depth] : reference depth for potential density.'
     PRINT *,'       [-dens smin width] : set minimum density and width for density tapering'
     PRINT *,'       [-dep  hmin width] : set minimum depth and width for depth tapering'
     PRINT *,'       [-lat  latmax width]: set max latitude and width for latitude tapering'
     PRINT *,'       [-o OUT-file] : specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ]  : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'             This option is effective only if cdftools are compiled with'
     PRINT *,'             a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'       Actual default values are :'
     PRINT '("             -refdep ", f6.0)', ref_dep
     PRINT '("             -dens   ", 2f8.3)', zsnmin, zswidth
     PRINT '("             -dep    ", 2f5.0)', hmin, hwidth
     PRINT '("             -lat    ", 2f5.0)', rlatmax, rlatwidth
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fmsk) 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option is used.'
     PRINT *,'         variables : wdmp'
     STOP 
  ENDIF

  ijarg=1
  DO WHILE (ijarg <= narg )
     CALL getarg(ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum)
     CASE ( '-t'     ) ; CALL getarg(ijarg, cf_tfil ) ; ijarg=ijarg+1
        ! options
     CASE ( '-s'     ) ; CALL getarg(ijarg, cf_sfil ) ; ijarg=ijarg+1
     CASE ( '-refdep') ; CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) ref_dep
     CASE ( '-dens'  ) ; CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) zsnmin
        ;                CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) zswidth
     CASE ( '-dep'   ) ; CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) hmin
        ;                CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) hwidth
     CASE ( '-lat'   ) ; CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) rlatmax
        ;                CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) rlatwidth
     CASE ( '-o'     ) ; CALL getarg(ijarg, cf_out  ) ; ijarg=ijarg+1
     CASE ( '-nc4'   ) ; lnc4 = .TRUE. 
     CASE DEFAULT      ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO
  IF ( cf_sfil == 'none' ) cf_sfil =cf_tfil

  WRITE(cglobal,'(" cdfmaskdmp -t ",a," -s ",a," -refdep ",f6.0," -dens ",2f8.3," -dep ",2f5.0," -lat ",2f5.0," -o ",a)')  &
       & TRIM(cf_tfil), TRIM(cf_sfil), ref_dep, zsnmin,zswidth, hmin, hwidth, rlatmax, rlatwidth,TRIM(cf_out)

  IF ( chkfile(cf_tfil) .OR. chkfile(cf_sfil) .OR. chkfile(cn_fmsk) ) STOP 99 ! missing files

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE (ztemp(npiglo,npjglo), zsal( npiglo,npjglo)                      )
  ALLOCATE (zsigi(npiglo,npjglo), zmask(npiglo,npjglo), zlat(npiglo,npjglo) )
  ALLOCATE (zwdmp(npiglo,npjglo) )
  ALLOCATE (dtim(npt) , zdep(npk) )

  CALL CreateOutput

  zdep(:)   = getvar1d(cf_tfil, cn_vdeptht, npk              )
  zlat(:,:) = getvar  (cf_tfil, cn_vlat2d,  1, npiglo, npjglo)

  DO jt = 1, npt
     PRINT *,'time: ',jt
     DO jk = 1, npk
        PRINT *, 'jk = ', jk
        ztemp(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
        zsal( :,:) = getvar(cf_sfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)
        zmask(:,:) = getvar(cn_fmsk, cn_tmask,    jk, npiglo, npjglo          )

        zsigi(:,:) = sigmai( ztemp, zsal, ref_dep, npiglo, npjglo)* zmask(:,:)

        DO jj=1,npjglo
           DO ji=1,npiglo

              wdep = TANH( (zdep(jk    ) - hmin   ) / hwidth   ) / 2. + 0.5
              wsig = TANH( (zsigi(ji,jj) - zsnmin ) / zswidth  ) / 2. + 0.5
              wlat = TANH(-(zlat( ji,jj) - rlatmax) / rlatwidth) / 2. + 0.5

              zwdmp(ji,jj) = wdep * wsig * wlat

           ENDDO
        ENDDO

        zwdmp(:,:) = zwdmp(:,:) * zmask(:,:)

        ierr = putvar(ncout, id_varout(1), zwdmp, jk,npiglo, npjglo, ktime=jt)

     END DO  ! loop to next level
  END DO  ! loop on time

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

    ipk(:)                    = npk  
    stypvar(1)%ichunk         = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname          = 'wdmp'
    stypvar(1)%cunits         = '[0-1]'
    stypvar(1)%rmissing_value = 1.e+20
    stypvar(1)%caxis          = 'TZYX'
    stypvar(1)%valid_min      = 0.
    stypvar(1)%valid_max      = 1.
    stypvar(1)%clong_name     = 'Damping mask build on density criteria'
    stypvar(1)%cshort_name    = 'wdmp'

    ! create output fileset
    ncout = create      (cf_out, cf_tfil,  npiglo, npjglo, npk      , ld_nc4=lnc4                  )
    ierr  = createvar   (ncout,  stypvar,  1,     ipk,     id_varout, ld_nc4=lnc4, cdglobal=cglobal)
    ierr  = putheadervar(ncout,  cf_tfil,  npiglo, npjglo, npk                        )

    dtim = getvar1d(cf_tfil, cn_vtimec,  npt )
    ierr = putvar1d(ncout,dtim, npt   , 'T'  )

  END SUBROUTINE CreateOutput

END PROGRAM cdfmaskdmp
