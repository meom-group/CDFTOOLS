PROGRAM cdfhdy
  !!======================================================================
  !!                     ***  PROGRAM  cdfhdy  ***
  !!=====================================================================
  !!  ** Purpose : Compute dynamical height anomaly field from gridT file
  !!               Store the results on a 2D cdf file.
  !!
  !!  ** Method  : the integral of (1/g) *10e4 * sum [ delta * dz ]
  !!               with delta = (1/rho - 1/rho0)
  !!               10e4 factor is conversion decibar/pascal
  !!
  !! History : 2.1  : 05/2010  : R. Dussin    : Original code
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos, ONLY : sigmai
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: jk, jt           ! dummy loop index
  INTEGER(KIND=4)                           :: ierr             ! working integer
  INTEGER(KIND=4)                           :: narg, iargc,ijarg! browse arguments
  INTEGER(KIND=4)                           :: npiglo, npjglo   ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt         !  "           "
  INTEGER(KIND=4)                           :: nlev1, nlev2     ! limit of vertical integration
  INTEGER(KIND=4)                           :: ncout            ! ncid of output fileset
  INTEGER(KIND=4), DIMENSION(1)             :: ipk              ! outptut variables : number of levels,
  INTEGER(KIND=4), DIMENSION(1)             :: id_varout        ! ncdf varid's

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: temp, sal        ! Temperature and salinity at current level
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: temp0, sal0      ! reference temperature and salinity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask            ! 2D mask at current level
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rdep, rdepth     ! depth at current level including SSH
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ssh              ! Sea Surface Heigh
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim, e3t_1d      ! time counter, vertical level spacing

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dhdy, dterm      ! dynamic height, working array
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dsig0, dsig      ! In situ density (reference, local)
  REAL(KIND=8)                              :: drau0 = 1000.d0  ! density of fresh water
  REAL(KIND=8)                              :: dgrav = 9.81d0   ! gravity

  CHARACTER(LEN=256)                        :: cf_tfil          ! input file name
  CHARACTER(LEN=256)                        :: cf_out='cdfhdy.nc' ! output file name
  CHARACTER(LEN=256)                        :: cldum            ! dummy string

  TYPE(variable) , DIMENSION(1)             :: stypvar          ! structure for attributes

  LOGICAL                                   :: lnc4 = .FALSE.   ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg <= 3 ) THEN
     PRINT *,' usage : cdfhdy -f T-file -limit level1 level2 [-o OUT-file] [-nc4] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute dynamical height anomaly field from gridT file.'
     PRINT *,'        It is computed as the integral of (1/g) *10e4 * sum [ delta * dz ]'
     PRINT *,'            where delta = (1/rho - 1/rho0)'
     PRINT *,'            10e4 factor is for the conversion decibar to pascal.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        -f T-file  : netcdf file with temperature and salinity'
     PRINT *,'        -limit level1  level2  : level1 is the upper limit for vertical '
     PRINT *,'                   integration (usually 1 = surface), while level2 is the'
     PRINT *,'                   lower limit for vertical integration.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        [-o OUT-file] : specify name of the output file instead of ', TRIM(cf_out) 
     PRINT *,'        [-nc4   ]     : Use netcdf4 output with chunking and deflation level 1'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fmsk),' and ', TRIM(cn_fzgr)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option is used.'
     PRINT *,'         variables : sohdy (m)'
     STOP
  ENDIF

  ijarg=1
  DO WHILE (ijarg <= narg )
     CALL getarg( ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'     ) ;  CALL getarg( ijarg, cf_tfil) ; ijarg=ijarg+1
     CASE ( '-limit' ) ;  CALL getarg( ijarg, cldum  ) ; ijarg=ijarg+1  ;  READ(cldum,*) nlev1
        ;  CALL getarg( ijarg, cldum  ) ; ijarg=ijarg+1  ;  READ(cldum,*) nlev2
        ! options
     CASE ( '-o'     ) ;  CALL getarg( ijarg, cf_out ) ; ijarg=ijarg+1
     CASE ( '-nc4'   ) ;  lnc4 = .TRUE.
     CASE DEFAULT      ;  PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP
     END SELECT
  ENDDO

  IF ( chkfile (cf_tfil) .OR. chkfile(cn_fmsk) .OR. chkfile(cn_fzgr)  ) STOP ! missing file

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'npt   =', npt

  ALLOCATE (temp0(npiglo,npjglo), sal0(npiglo,npjglo), dsig0(npiglo,npjglo) ,tmask(npiglo,npjglo))
  ALLOCATE (temp(npiglo,npjglo), sal(npiglo,npjglo), dsig(npiglo,npjglo) , dhdy(npiglo,npjglo), dterm(npiglo,npjglo))
  ALLOCATE (rdep(npiglo,npjglo), rdepth(npiglo,npjglo), ssh(npiglo,npjglo), e3t_1d(npk))
  ALLOCATE (tim(npt))

  CALL CreateOutput

  ! Temperature and salinity for reference profile
  temp0(:,:) =  0.
  sal0 (:,:) = 35.

  tmask(:,:) = getvar(cn_fmsk,  cn_tmask, nlev2, npiglo, npjglo)
  ssh(:,:)   = getvar(cf_tfil,  cn_sossheig, 1,  npiglo, npjglo)
  e3t_1d(:)  = getvare3(cn_fzgr, cn_ve3t, npk)

  DO jt=1,npt
     PRINT *,' TIME = ', jt, tim(jt)/86400.,' days'
     dhdy  (:,:) = 0.d0
     rdepth(:,:) = 0.

     DO jk = nlev1, nlev2

        !rdep(:,:)   = getvar(cn_fzgr, 'e3t_ps', jk,npiglo,npjglo,ldiom=.true.)
        ! we degrade the computation to smooth the results
        rdep(:,:) = e3t_1d(jk)

        IF ( jk == 1 ) THEN
           rdep(:,:) = rdep(:,:) + ssh(:,:)
        ENDIF

        ! depth at current level, including ssh (used for computation of rho in situ)
        rdepth(:,:) = rdepth(:,:) + rdep(:,:)

        temp(:,:) = getvar(cf_tfil, cn_votemper,  jk, npiglo, npjglo, ktime=jt)
        sal (:,:) = getvar(cf_tfil, cn_vosaline,  jk, npiglo, npjglo, ktime=jt)

        dsig0 = sigmai(temp0, sal0, rdepth, npiglo, npjglo) 
        dsig  = sigmai(temp , sal , rdepth, npiglo, npjglo) 

        ! we compute the term of the integral : (1/g) *10e4 * sum [ delta * dz ]
        ! with delta = (1/rho - 1/rho0)
        ! 10e4 factor is conversion decibar/pascal
        !
        dterm = ( ( 1.d0 / ( drau0 + dsig(:,:) ) ) - ( 1.d0 / ( drau0 + dsig0(:,:) ) ) ) * 10000.d0 * rdep / dgrav
        ! in land, it seems appropriate to stop the computation
        WHERE(sal == 0 ) dterm = 0

        dhdy(:,:) = dhdy(:,:) + dterm(:,:)

     END DO  ! loop to next level

     ! we mask with the last level of the integral
     dhdy(:,:) = dhdy(:,:) * tmask(:,:)

     ierr = putvar(ncout, id_varout(1) ,REAL(dhdy), 1, npiglo, npjglo, ktime=jt)

  END DO  ! next time frame

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
    ipk(:) = 1 
    stypvar(1)%ichunk           = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname            = 'sohdy'
    stypvar(1)%cunits           = 'm'
    stypvar(1)%rmissing_value   = 0.
    stypvar(1)%valid_min        = -100.
    stypvar(1)%valid_max        =  100.
    stypvar(1)%clong_name       = 'Dynamical height anomaly'
    stypvar(1)%cshort_name      = 'sohdy'
    stypvar(1)%conline_operation= 'N/A'
    stypvar(1)%caxis            = 'TYX'

    ! create output fileset
    ncout = create      (cf_out, cf_tfil,  npiglo, npjglo, npk       , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar,  1,      ipk,    id_varout , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_tfil,  npiglo, npjglo, npk       )

    tim   = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr  = putvar1d(ncout,   tim,       npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfhdy
