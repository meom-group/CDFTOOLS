PROGRAM cdfbottomsig
  !!======================================================================
  !!                     ***  PROGRAM  cdfbottomsig  ***
  !!=====================================================================
  !!  ** Purpose : Compute the bottom sigma from gridT file.
  !!               Store the results on a 'similar' cdf file.
  !!
  !!  **  Method:  Uses vosaline do determine the bottom points. A depth
  !!               reference can be specify to compute density refered to
  !!               this depth.
  !!
  !! History : 2.1  : 11/2005  : J.M. Molines : Original code
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: jk, jt         ! dummy loop index
  INTEGER(KIND=4)                            :: ierr           ! working integer
  INTEGER(KIND=4)                            :: narg, iargc    ! 
  INTEGER(KIND=4)                            :: npiglo, npjglo ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt       ! size of the domain
  INTEGER(KIND=4)                            :: ncout          ! ncid of output file
  INTEGER(KIND=4), DIMENSION(1)              :: ipk            ! outptut variables : number of levels,
  INTEGER(KIND=4), DIMENSION(1)              :: id_varout      ! ncdf varid's
  INTEGER(KIND=4), DIMENSION(2)              :: ismin, ismax   ! location of min and max sigmabot

  REAL(KIND=4)                               :: zsigmn, zsigmx ! value of min and max of sigmabot
  REAL(KIND=4)                               :: zref           ! value of min and max of sigmabot
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: ztemp, zsal    ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: ztemp0, zsal0  ! temporary array to read temp, sal
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zsig           ! potential density 
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zmask          ! 2D mask at surface
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: tim            ! time counter

  CHARACTER(LEN=256)                         :: cf_out='botsig.nc' ! Output file name
  CHARACTER(LEN=256)                         :: cf_tfil        ! input filename
  CHARACTER(LEN=256)                         :: cv_sig         ! output variable name
  CHARACTER(LEN=256)                         :: cref           ! message for depth reference
  CHARACTER(LEN=256)                         :: cldum          ! dummy char variable

  TYPE (variable), DIMENSION(1)              :: stypvar        ! structure for attributes

  LOGICAL                                    :: lsigi=.FALSE.  ! flag for sigma-i computation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfbottomsig  T-file [zref]' 
     PRINT *,'      '
     PRINT *,'     PURPOSE :' 
     PRINT *,'       Create a 2D file with bottom density. In case a depth reference' 
     PRINT *,'       is given, the density is refered to this depth. By default sigma-0'
     PRINT *,'       is used. Bottom most point is determined from the last non zero '
     PRINT *,'       salinity point in the water column.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file : input file with temperature and salinity '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [zref] : depth reference for potential density'
     PRINT *,'             If not given assume sigma-0'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : sobotsig0 or sobotsigi ( kg/m3 - 1000 )' 
     STOP
  ENDIF

  cv_sig = 'sobotsig0'
  cref=''
  CALL getarg (1, cf_tfil)
  IF ( chkfile(cf_tfil) ) STOP ! missing file

  IF ( narg == 2 ) THEN
     lsigi = .TRUE.
     CALL getarg (2, cldum) ; READ(cldum,*) zref
     cv_sig = 'sobotsigi'
     WRITE(cref,'("_refered_to_",i4.4,"_m")') NINT(zref)
  ENDIF

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  ipk(:)= 1  ! all variables (input and output are 3D)

  stypvar(1)%cname             = cv_sig
  stypvar(1)%cunits            = 'kg/m3'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = 0.001
  stypvar(1)%valid_max         = 40.
  stypvar(1)%clong_name        = 'Bottom_Potential_density'//TRIM(cref)
  stypvar(1)%cshort_name       = cv_sig
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TYX'

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  ALLOCATE (ztemp( npiglo,npjglo), zsal( npiglo,npjglo), zsig(npiglo,npjglo) ,zmask(npiglo,npjglo))
  ALLOCATE (ztemp0(npiglo,npjglo), zsal0(npiglo,npjglo) )
  ALLOCATE ( tim (npt) )

  ! create output fileset

  ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk       )
  ierr  = createvar   (ncout,  stypvar, 1     , ipk   , id_varout )
  ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk       )

  zsal  = 0.
  ztemp = 0.
  zmask = 1.

  DO jt = 1, npt
     DO jk = 1, npk
        PRINT *,'level ',jk
        zsal0(:,:) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)
        ztemp0(:,:)= getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
        IF (jk == 1  )  THEN
           WHERE( zsal0 == 0. ) zmask=0.
        END IF
        WHERE ( zsal0 /= 0 )
          zsal=zsal0 ; ztemp=ztemp0
        END WHERE
     ENDDO
     
     IF (lsigi ) THEN
        zsig(:,:) = sigmai ( ztemp, zsal, zref, npiglo, npjglo ) * zmask(:,:)
     ELSE
        zsig(:,:) = sigma0 ( ztemp, zsal,       npiglo, npjglo ) * zmask(:,:)
     ENDIF


     zsigmn=minval(zsig(2:npiglo-1,2:npjglo-1), zmask(2:npiglo-1,2:npjglo-1)==1)
     zsigmx=maxval(zsig(2:npiglo-1,2:npjglo-1), zmask(2:npiglo-1,2:npjglo-1)==1)
     ismin= minloc(zsig(2:npiglo-1,2:npjglo-1), zmask(2:npiglo-1,2:npjglo-1)==1)
     ismax= maxloc(zsig(2:npiglo-1,2:npjglo-1), zmask(2:npiglo-1,2:npjglo-1)==1)

     PRINT *,'Bottom density : min = ', zsigmn,' at ', ismin(1), ismin(2)
     PRINT *,'               : max = ', zsigmx,' at ', ismax(1), ismax(2)

     ierr = putvar(ncout, id_varout(1), zsig, 1, npiglo, npjglo, ktime=jt)
  ENDDO

  tim  = getvar1d(cf_tfil, cn_vtimec, npt     )
  ierr = putvar1d(ncout,   tim      , npt, 'T')
  ierr = closeout(ncout)

END PROGRAM cdfbottomsig
