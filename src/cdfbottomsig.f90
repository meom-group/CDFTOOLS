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
  !! @class bottom
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: jk, jt         ! dummy loop index
  INTEGER(KIND=4)                            :: ierr           ! working integer
  INTEGER(KIND=4)                            :: narg, iargc    ! 
  INTEGER(KIND=4)                            :: ijarg          ! argument counter
  INTEGER(KIND=4)                            :: npiglo, npjglo ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt       ! size of the domain
  INTEGER(KIND=4)                            :: ncout          ! ncid of output file
  INTEGER(KIND=4), DIMENSION(1)              :: ipk            ! outptut variables : number of levels,
  INTEGER(KIND=4), DIMENSION(1)              :: id_varout      ! ncdf varid's
  INTEGER(KIND=4), DIMENSION(2)              :: ismin, ismax   ! location of min and max sigmabot

  REAL(KIND=4)                               :: zsigmn, zsigmx ! value of min and max of sigmabot
  REAL(KIND=4)                               :: zref           ! value of min and max of sigmabot
  REAL(KIND=4)                               :: zsps           ! Missing value for salinity
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: ztemp, zsal    ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: ztemp0, zsal0  ! temporary array to read temp, sal
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zsig           ! potential density 
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zmask          ! 2D mask at surface

  REAL(KIND=8), DIMENSION(:),    ALLOCATABLE :: dtim           ! time counter

  CHARACTER(LEN=256)                         :: cf_out='botsig.nc' ! Output file name
  CHARACTER(LEN=256)                         :: cf_tfil        ! input filename
  CHARACTER(LEN=256)                         :: cv_sig         ! output variable name
  CHARACTER(LEN=256)                         :: cref           ! message for depth reference
  CHARACTER(LEN=256)                         :: cldum          ! dummy char variable

  TYPE (variable), DIMENSION(1)              :: stypvar        ! structure for attributes

  LOGICAL                                    :: lsigi  =.FALSE.! flag for sigma-i computation
  LOGICAL                                    :: lsigntr=.FALSE.! flag for sigma-Neutral computation
  LOGICAL                                    :: lnc4   =.FALSE.! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfbottomsig  -t T-file [-r REF-depth ] [-ntr] [-o OUT-file] [-nc4]' 
     PRINT *,'      '
     PRINT *,'     PURPOSE :' 
     PRINT *,'       Create a 2D file with bottom density. In case a depth reference is ' 
     PRINT *,'       given, the density is referred to this depth. By default sigma-0 is'
     PRINT *,'       used. Bottom most point is determined from the last non zero salinity'
     PRINT *,'       point in the water column.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file : input file with temperature and salinity.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-r REF-depth] : depth reference for potential density.'
     PRINT *,'             Without -r nor -ntr options sigma-0 is assumed.'
     PRINT *,'       [-ntr ]: Will use neutral density.'
     PRINT *,'             Without -r nor -ntr options sigma-0 is assumed.'
     PRINT *,'       [-o OUT-file] : Specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4]:  Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'             This option is effective only if cdftools are compiled with'
     PRINT *,'             a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     OPENMP SUPPORT : yes'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out),' unless -o option is used.'
     PRINT *,'         variables : sobotsig0 or sobotsigi ( kg/m3 - 1000 )' 
     PRINT *,'                     or sobotsigntr (kg/m3)'
     PRINT *,'      '
     STOP 
  ENDIF

  cv_sig = 'sobotsig0'
  cref=''
  ijarg=1
  DO WHILE ( ijarg <= narg )
     CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-t'  ) ; CALL getarg(ijarg, cf_tfil) ;  ijarg=ijarg+1
        ! options
     CASE ( '-r'  ) ; CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) zref
        ;             lsigi  = .TRUE. 
        ;             cv_sig = 'sobotsigi'
        ;             WRITE(cref,'("_refered_to_",i4.4,"_m")') NINT(zref)
     CASE ( '-ntr') ; lsigntr = .TRUE.
        ;             cv_sig  = 'sobotsigntr'
        ;             WRITE(cref,'("_Neutral")')
     CASE ( '-o'  ) ; CALL getarg(ijarg, cf_out ) ;  ijarg=ijarg+1
     CASE ( '-nc4') ; lnc4     = .TRUE.
     CASE DEFAULT   ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  IF ( chkfile(cf_tfil) ) STOP 99 ! missing file
  ! look for MissingValue for salinity
  zsps = getspval(cf_tfil, cn_vosaline)


  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  ALLOCATE (ztemp( npiglo,npjglo), zsal( npiglo,npjglo), zsig(npiglo,npjglo) ,zmask(npiglo,npjglo))
  ALLOCATE (ztemp0(npiglo,npjglo), zsal0(npiglo,npjglo) )
  ALLOCATE ( dtim (npt) )

  CALL CreateOutput

  zsal  = 0.
  ztemp = 0.
  zmask = 1.

  DO jt = 1, npt
     DO jk = 1, npk
        PRINT *,'level ',jk
        zsal0(:,:) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)
        ztemp0(:,:)= getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
        IF (jk == 1  )  THEN
           WHERE( zsal0 == zsps ) zmask=0.
        END IF
        WHERE ( zsal0 /= zsps )
           zsal=zsal0 ; ztemp=ztemp0
        END WHERE
     ENDDO

     IF (lsigi ) THEN
        zsig(:,:) = sigmai ( ztemp, zsal, zref, npiglo, npjglo ) * zmask(:,:)
     ELSE IF (lsigntr ) THEN
        zsig(:,:) = sigmantr (ztemp, zsal,      npiglo, npjglo )* zmask(:,:)
     ELSE
        zsig(:,:) = sigma0 ( ztemp, zsal,       npiglo, npjglo ) * zmask(:,:)
     ENDIF

     zsigmn=MINVAL(zsig(2:npiglo-1,2:npjglo-1), zmask(2:npiglo-1,2:npjglo-1)==1)
     zsigmx=MAXVAL(zsig(2:npiglo-1,2:npjglo-1), zmask(2:npiglo-1,2:npjglo-1)==1)
     ismin= MINLOC(zsig(2:npiglo-1,2:npjglo-1), zmask(2:npiglo-1,2:npjglo-1)==1)
     ismax= MAXLOC(zsig(2:npiglo-1,2:npjglo-1), zmask(2:npiglo-1,2:npjglo-1)==1)

     PRINT *,'Bottom density : min = ', zsigmn,' at ', ismin(1), ismin(2)
     PRINT *,'               : max = ', zsigmx,' at ', ismax(1), ismax(2)

     ierr = putvar(ncout, id_varout(1), zsig, 1, npiglo, npjglo, ktime=jt)
  ENDDO

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
    ipk(:)= 1  ! all variables (input and output are 3D)

    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = cv_sig
    stypvar(1)%cunits            = 'kg/m3'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = 0.001
    stypvar(1)%valid_max         = 40.
    stypvar(1)%clong_name        = 'Bottom_Potential_density'//TRIM(cref)
    stypvar(1)%cshort_name       = cv_sig
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TYX'

    ! create output fileset
    ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk      , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar, 1     , ipk   , id_varout, ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk       )

    dtim = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr = putvar1d(ncout,   dtim     , npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfbottomsig
