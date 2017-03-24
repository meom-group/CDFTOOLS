PROGRAM cdfeddyscale
  !!======================================================================
  !!                     ***  PROGRAM  cdfeddyscale  ***
  !!=====================================================================
  !!  ** Purpose : Compute: -the Taylor scale or large scale eddy (lambda1)
  !!                        -the small scale eddy (lambda2) 
  !!                        -and the inertial range (scar) on F-points
  !!               lambda1 = sqrt(mean Kinetic Energie / Enstrophy)
  !!               lambda2 = sqrt(Enstrophy / Palinstrophy)
  !!               scar    = lambda1 / lambda2
  !!
  !!               Enstrophy = 1/2 * ( mean((RV)^2) )
  !!               Palinstrophy = 1/2 * ( mean((dx(RV))^2 + (dy(RV))^2) ) 
  !!
  !!  ** Method  : Use the mean of the variables (vozocrtx2 + vomecrty2), socurl2, 
  !!               (sodxcurl2 + sodycurl2) of the cdfeddyscale_pass1.f90  
  !!               
  !!  ** Warning : - the square of curl socurl2 is on F-points, 
  !!               - vozocrtx2, vomecrty2, sodxcurl2 and sodycurl2 are on UV-points
  !!              
  !! History : 3.0  : 12/2013  : C.Q. Akuetevi (original)
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence  (Licence/CDFTOOLSCeCILL.txt)
  !! @class energy_diagnostics
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jt         ! dummy loop index
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                           :: narg, iargc, ijarg ! browse command line
  INTEGER(KIND=4)                           :: ncout, ierr        ! browse command line
  INTEGER(KIND=4), DIMENSION(3)             :: ipk, id_varout     ! output variable properties

  REAL(KIND=4), DIMENSION(1)                :: tim                ! time counter in output file
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zun, zvn           ! working arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rotn2              ! square curl
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zdxrotn2, zdyrotn2 ! square curl gradient components
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: vozocrtx2          ! square of velocity components
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: vomecrty2          ! square of velocity components
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmke               ! mean kinetic energy
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ens                ! enstrophy
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zpal               ! palinstrophy
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zlambda1           ! Taylor or large scale eddy
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zlambda2           ! small scale eddy
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: scar               ! inertial range

  CHARACTER(LEN=256)                        :: cf_meanfil   ! file names
  CHARACTER(LEN=256)                        :: cf_out = 'lambda.nc' ! output file name
  CHARACTER(LEN=256)                        :: cldum              ! dummy string

  TYPE (variable), DIMENSION(3)             :: stypvar            ! structure for attibutes

  LOGICAL                                   :: lchk     = .FALSE. ! flag for missing files
  LOGICAL                                   :: lperio   = .FALSE. ! flag for E-W periodicity
  LOGICAL                                   :: lnc4     = .FALSE. ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfeddyscale -f PASS1-meanfile [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        This program computes 3 quantities relative to the eddy scales,'
     PRINT *,'        from a file preprocessed by cdfeddyscale_pass1 :'
     PRINT *,'              -the Taylor scale or large scale eddy (lambda1)'
     PRINT *,'              -the small scale eddy (lambda2)'
     PRINT *,'              -and the inertial range (scar) on F-points'
     PRINT *,'      '
     PRINT *,'           lambda1 = sqrt(mean Kinetic Energie / Enstrophy)'
     PRINT *,'           lambda2 = sqrt(Enstrophy / Palinstrophy)'
     PRINT *,'           Inertial Range    = lambda1 / lambda2'
     PRINT *,'      ' 
     PRINT *,'           Enstrophy = 1/2 * ( mean((RV)^2) )'
     PRINT *,'           Palenstrophy = 1/2 * ( mean((dx(RV))^2 + (dy(RV))^2) )'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f  PASS1-meanfile : time average of cdfeddyscale_pass1 preprocessed'
     PRINT *,'           files.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file ]: specify the output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ]  : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'              This option is effective only if cdftools are compiled with'
     PRINT *,'              a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out),' unless -o option is used.'
     PRINT *,'         variables : solambda1 (m), solambda2 (m), soscar'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfeddyscale_pass1 '
     STOP
  ENDIF

  ijarg=1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_meanfil ) ; ijarg=ijarg+1
        ! option
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out     ) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP
     END SELECT
  ENDDO

  lchk = chkfile(cf_meanfil ) .OR. lchk
  IF ( lchk ) STOP ! missing files

  npiglo = getdim(cf_meanfil,cn_x)
  npjglo = getdim(cf_meanfil,cn_y)

  PRINT *, 'npiglo = ',npiglo
  PRINT *, 'npjglo = ',npjglo

  ! Allocate the memory
  ALLOCATE ( rotn2(npiglo,npjglo)                                )
  ALLOCATE ( zdxrotn2(npiglo,npjglo)  , zdyrotn2(npiglo,npjglo)  )
  ALLOCATE ( vozocrtx2(npiglo,npjglo) , vomecrty2(npiglo,npjglo) )
  ALLOCATE ( zmke(npiglo,npjglo)                                 )
  ALLOCATE ( ens(npiglo,npjglo)                                  )
  ALLOCATE ( zpal(npiglo,npjglo)                                 )
  ALLOCATE ( zlambda1(npiglo,npjglo)   , zlambda2(npiglo,npjglo) )
  ALLOCATE ( scar(npiglo,npjglo)                                 )

  CALL CreateOutputFile

  !load the mean variables
  rotn2(:,:)     =  getvar(cf_meanfil, 'socurl2'  , 1,npiglo,npjglo )
  zdxrotn2(:,:)  =  getvar(cf_meanfil, 'sodxcurl2', 1,npiglo,npjglo )
  zdyrotn2(:,:)  =  getvar(cf_meanfil, 'sodycurl2', 1,npiglo,npjglo )
  vozocrtx2(:,:) =  getvar(cf_meanfil, 'vozocrtx2', 1,npiglo,npjglo )
  vomecrty2(:,:) =  getvar(cf_meanfil, 'vomecrty2', 1,npiglo,npjglo )

  ! Enstrophy
  ens(:,:) = 0.5 * rotn2(:,:)

  ! compute the Kinetic Energy on F-points
  zmke = -9999.
  DO jj = 1, npjglo-1
     DO ji = 1, npiglo-1   ! vector opt.
        zmke(ji,jj) =  0.25 *( (vozocrtx2(ji,jj+1) + vozocrtx2(ji,jj))   &
             &                +(vomecrty2(ji+1,jj) + vomecrty2(ji,jj))   )
     END DO
  END DO

  ! compute the Palinstrophy on  F-points
  zpal = -9999.
  DO jj = 1, npjglo-1 
     DO ji = 1, npiglo-1   ! vector opt.
        zpal(ji,jj) =  0.25 * (zdxrotn2(ji+1,jj) + zdxrotn2(ji,jj))  &
             & + 0.25 * (zdyrotn2(ji,jj+1) + zdyrotn2(ji,jj))   
     END DO
  END DO

  ! compute the Taylor and small scale eddy
  zlambda1(:,:) = 0.
  zlambda2(:,:) = 0.
  WHERE( ens > 0.  ) zlambda1(:,:) = SQRT ( zmke(:,:)/ens(:,:))
  WHERE( zpal > 0. ) zlambda2(:,:) = SQRT ( ens(:,:)/zpal(:,:))

  ! compute the Inertial Range
  scar(:,:) = 0. 
  WHERE( zlambda2 > 0. ) scar(:,:) = zlambda1(:,:)/zlambda2(:,:) 

  ! write zlambda1 on file 
  ierr = putvar(ncout, id_varout(1), zlambda1, 1, npiglo, npjglo )

  ! write zlamdba2 on file
  ierr = putvar(ncout, id_varout(2), zlambda2, 1, npiglo, npjglo )

  ! write scar on file
  ierr = putvar(ncout, id_varout(3), scar,     1, npiglo, npjglo ) 

  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutputFile
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutputFile  ***
    !!
    !! ** Purpose :  Create Output file 
    !!
    !! ** Method  :  Use global variable 
    !!
    !!----------------------------------------------------------------------

    ! define new variables for output
    ! Taylor or large scale eddy F point
    ipk(1)                       = 1   !2D
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = 'solambda1'
    stypvar(1)%cunits            = 'm'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = 0.
    stypvar(1)%valid_max         =  100000.
    stypvar(1)%clong_name        = 'Taylor_large_eddy_scale (lambda1)'
    stypvar(1)%cshort_name       = 'solambda1'
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TYX'

    ! Small scale eddy F point
    ipk(2)                       = 1   !2D
    stypvar(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname             = 'solambda2'
    stypvar(2)%cunits            = 'm'
    stypvar(2)%rmissing_value    = 0.
    stypvar(2)%valid_min         = 0.
    stypvar(2)%valid_max         =  100000.
    stypvar(2)%clong_name        = 'Small scale eddy (lambda2)'
    stypvar(2)%cshort_name       = 'solambda2'
    stypvar(2)%conline_operation = 'N/A'
    stypvar(2)%caxis             = 'TYX'

    ! Inertial Range F point
    ipk(3)                       = 1   !2D
    stypvar(3)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(3)%cname             = 'soscar'
    stypvar(3)%cunits            = ''
    stypvar(3)%rmissing_value    = 0.
    stypvar(3)%valid_min         = 0.
    stypvar(3)%valid_max         =  20.
    stypvar(3)%clong_name        = 'Inertial range (scar)'
    stypvar(3)%cshort_name       = 'soscar'
    stypvar(3)%conline_operation = 'N/A'
    stypvar(3)%caxis             = 'TYX'

    ! create output fileset
    ncout = create      (cf_out, cf_meanfil, npiglo, npjglo, 0     , ld_nc4=lnc4 )
    ierr  = createvar   (ncout , stypvar, 3,      ipk,    id_varout, ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_meanfil, npiglo, npjglo, 0 )

    tim  = getvar1d(cf_meanfil, cn_vtimec, 1      )
    ierr = putvar1d(ncout,      tim,       1,  'T')

  END SUBROUTINE CreateOutputFile

END PROGRAM cdfeddyscale

