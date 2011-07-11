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
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jk, jt     ! dummy loop index
  INTEGER(KIND=4)                           :: ierr               ! error status
  INTEGER(KIND=4)                           :: narg, iargc        ! browse command line
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                           :: ncout              ! ncid of output file
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout     ! level and varid's

  REAL(KIND=4)                              :: ref_dep=2000.      ! reference depth in meters
  REAL(KIND=4)                              :: zsnmin=37.16       ! minimum density
  REAL(KIND=4)                              :: zswidth=0.025      ! tapering width
  REAL(KIND=4)                              :: hmin=1000.         ! depth limit
  REAL(KIND=4)                              :: hwidth=100.        ! depth tapering height
  REAL(KIND=4)                              :: rlatmax=-20        ! max latitude
  REAL(KIND=4)                              :: rlatwidth=2        ! latitude tapering width
  REAL(KIND=4)                              :: wdep, wsig, wlat   ! tapering function dep, sigma and lat
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztemp              ! temperature
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsal               ! salinity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsigi              ! sigma-i
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask              ! 2D mask at current level
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zwdmp              ! 2D build mask at current level
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zlat               ! latitudes
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                ! time counter
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: zdep               ! deptht

  CHARACTER(LEN=256)                        :: cf_tfil            ! input filename for temperature
  CHARACTER(LEN=256)                        :: cf_sfil            ! input filename for salinity
  CHARACTER(LEN=256)                        :: cf_out='mask_dmp.nc' ! output file name
  CHARACTER(LEN=256)                        :: cldum              ! dummy string

  TYPE (variable), DIMENSION(1)             :: stypvar            ! structure for attributes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmaskdmp T-file S-file  ... '
     PRINT *,'               ... [ref_dep snmin swidth hmin hwidth latmax latwidth]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute a damping mask with smooth transition according to density,'
     PRINT *,'       depth and latitude criteria.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file : temperature file' 
     PRINT *,'       S-file : salinity file' 
     PRINT *,'        They can be the same file, but as many climatologied are provided'
     PRINT *,'        in separate files, we decided to put both in the command line.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        ** If used, they must all be provided in the correct order (!) **'
     PRINT *,'       ref_dep  : reference depth for potential density.'
     PRINT *,'       snmin    : density minimum for the mask.'
     PRINT *,'       swidth   : density width for tapering'
     PRINT *,'       hmin     : minimum depth'
     PRINT *,'       hwidth   : depth width  for tapering'
     PRINT *,'       latmax   : maximum latitude'
     PRINT *,'       latwidth : latitude width  for tapering'
     PRINT *,'      '
     PRINT *,'       Actual default values are :'
     PRINT *,'        ref_dep  = ', ref_dep
     PRINT *,'        snmin    = ', zsnmin
     PRINT *,'        swidth   = ', zswidth
     PRINT *,'        hmin     = ', hmin
     PRINT *,'        hwidth   = ', hwidth
     PRINT *,'        latmax   = ', rlatmax
     PRINT *,'        latwidth = ', rlatwidth
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fmsk) 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : wdmp'
     STOP
  ENDIF

  IF ( narg > 2 .AND. narg < 9 ) THEN
     PRINT *,'wrong number of arguments'
     STOP
  ENDIF

  CALL getarg (1, cf_tfil)
  CALL getarg (2, cf_sfil)

  IF ( chkfile(cf_tfil) .OR. chkfile(cf_sfil) ) STOP ! missing files

  IF ( narg == 9 ) THEN
     CALL getarg (3, cldum) ; READ(cldum,*) ref_dep
     CALL getarg (4, cldum) ; READ(cldum,*) zsnmin
     CALL getarg (5, cldum) ; READ(cldum,*) zswidth
     CALL getarg (6, cldum) ; READ(cldum,*) hmin
     CALL getarg (7, cldum) ; READ(cldum,*) hwidth
     CALL getarg (8, cldum) ; READ(cldum,*) rlatmax
     CALL getarg (9, cldum) ; READ(cldum,*) rlatwidth
  ENDIF

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  ipk(:)                    = npk  
  stypvar(1)%cname          = 'wdmp'
  stypvar(1)%rmissing_value = 1.e+20
  stypvar(1)%caxis          = 'TZYX'

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE (ztemp(npiglo,npjglo), zsal( npiglo,npjglo)                      )
  ALLOCATE (zsigi(npiglo,npjglo), zmask(npiglo,npjglo), zlat(npiglo,npjglo) )
  ALLOCATE (zwdmp(npiglo,npjglo) )
  ALLOCATE (tim(npt) , zdep(npk) )

  ! create output fileset
  ncout = create      (cf_out, cf_tfil,  npiglo, npjglo, npk       )
  ierr  = createvar   (ncout,  stypvar,  1,     ipk,     id_varout )
  ierr  = putheadervar(ncout,  cf_tfil,  npiglo, npjglo, npk       )

  tim(:)    = getvar1d(cf_tfil, cn_vtimec,  npt              )
  zdep(:)   = getvar1d(cf_tfil, cn_vdeptht, npk              )
  zlat(:,:) = getvar  (cf_tfil, cn_vlat2d,  1, npiglo, npjglo)

  ierr=putvar1d(ncout, tim, npt, 'T')

  DO jt = 1, npt
     PRINT *,'time: ',jt
     DO jk = 1, npk
        PRINT *, 'jk = ', jk
        ztemp(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
        zsal( :,:) = getvar(cf_sfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)
        zmask(:,:) = getvar(cn_fmsk, 'tmask',     jk, npiglo, npjglo          )

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

END PROGRAM cdfmaskdmp
