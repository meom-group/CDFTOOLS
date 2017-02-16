PROGRAM cdfmxlsaltc
  !!======================================================================
  !!                     ***  PROGRAM  cdfmxlsaltc  ***
  !!=====================================================================
  !!  ** Purpose : Compute the salt content in the mixed layer. Work for
  !!               partial steps (default) or full step (-full option)
  !!
  !!  ** Method  : compute the sum ( rho S  * e1 *e2 * e3 *mask )
  !!               for the mixed layer stored into gridT file
  !!
  !! History : 2.1  : 04/2006  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modutils
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jt              ! dummy loop index
  INTEGER(KIND=4)                               :: narg, iargc         ! command line 
  INTEGER(KIND=4)                               :: ijarg, ireq         ! command line 
  INTEGER(KIND=4)                               :: npiglo, npjglo      ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt            ! size of the domain,
  INTEGER(KIND=4)                               :: ncout, ierr         ! ncid and error status
  INTEGER(KIND=4), DIMENSION(1)                 :: ipk, id_varout      ! levels and varid's of output vars

  REAL(KIND=4), PARAMETER                       :: rprho0=1020.        ! rho reference density
  REAL(KIND=4), PARAMETER                       :: rpcp=4000.          ! calorific capacity
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: e3                  ! metrics
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zs                  ! temperature in the MXL
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zmxl                ! depth of the MXL
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zmask               ! mask
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: gdepw               ! vertical levels
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: e31d                ! vertical metric full
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: tim                 ! time counter
  REAL(KIND=4), DIMENSION(1)                    :: rdep                ! dummy depth output

  REAL(KIND=8)                                  :: dvol                ! total volume
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dmxlsaltc           ! heat content

  CHARACTER(LEN=2048)                            :: cf_tfil             ! input file name
  CHARACTER(LEN=2048)                            :: cf_out='mxlsaltc.nc'! output file
  CHARACTER(LEN=2048)                            :: cv_out='somxlsaltc' ! input file name
  CHARACTER(LEN=2048)                            :: cglobal             ! global attribute
  CHARACTER(LEN=2048)                            :: cldum               ! dummy string

  TYPE(variable), DIMENSION(1)                  :: stypvar             ! stucture for attributes (output)

  LOGICAL                                       :: lfull=.false.       ! full step flag
  LOGICAL                                       :: lchk                ! file existence flag (true if missing)
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmxlsaltc T-file [-full ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the salt content in the mixed layer.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file : netcdf file with salinity and mixed layer deptht.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-full ] : indicate a full step configuration.'
     PRINT *,'       [-o OUT-file ] : specify output file instead of ',TRIM(cf_out) 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fzgr),' and ', TRIM(cn_fmsk) 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ', TRIM(cv_out),' (kg/m2 )'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfmxl, cdfmxlhcsc, cdfmxlheatc ' 
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 ; ireq = 0
  
  DO WHILE ( ijarg <= narg )
    CALL getarg (ijarg, cldum   ) ; ijarg = ijarg + 1 
    SELECT CASE ( cldum )
    CASE ( '-full'    ) ; lfull = .true.
    CASE ( '-partial' ) ; lfull = .false.
    CASE ('-o') 
        CALL getarg (ijarg, cf_out) ; ijarg=ijarg+1
    CASE DEFAULT 
      ireq=ireq+1
      SELECT CASE ( ireq )
      CASE ( 1 )    ; cf_tfil=cldum    
      CASE DEFAULT  ; PRINT *,' Too many arguments'
      END SELECT
    END SELECT
  END DO

  lchk = chkfile (cn_fzgr)
  lchk = chkfile (cn_fmsk) .OR. lchk
  lchk = chkfile (cf_tfil  ) .OR. lchk
  IF ( lchk ) STOP ! missing files

  CALL SetGlobalAtt( cglobal )

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  rdep(1)                      = 0.
  ipk(:)                       = 1
  stypvar(1)%cname             = cv_out
  stypvar(1)%cunits            = 'kg/m2'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = 0
  stypvar(1)%valid_max         =  1.e9
  stypvar(1)%clong_name        = 'Mixed_Layer_Salt_Content'
  stypvar(1)%cshort_name       = cv_out
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TYX'

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo), dmxlsaltc(npiglo, npjglo) )
  ALLOCATE ( zs(npiglo,npjglo), zmxl(npiglo,npjglo)          )
  ALLOCATE ( e3(npiglo,npjglo)                               )
  ALLOCATE ( gdepw(npk), tim(npt)                            )

  IF ( lfull ) ALLOCATE ( e31d(npk)                          )

  ! Initialize output file
  ncout = create      (cf_out, cf_tfil, npiglo, npjglo, 1                           )
  ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout, cdglobal=cglobal )
  ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, 1, pdep=rdep                )

  tim  = getvar1d(cf_tfil, cn_vtimec, npt     )
  ierr = putvar1d(ncout,   tim,       npt, 'T')

               gdepw(:) = getvare3(cn_fzgr, cn_gdepw, npk)
  IF ( lfull ) e31d( :) = getvare3(cn_fzgr, cn_ve3t,  npk)



  DO jt=1,npt
     dvol= 0.d0
     dmxlsaltc(:,:) = 0.d0
     zmxl( :,:) = getvar(cf_tfil, cn_somxl010, 1,  npiglo, npjglo, ktime=jt)

     DO jk = 1, npk
        zs(   :,:) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)
        zmask(:,:) = getvar(cn_fmsk, 'tmask',     jk, npiglo, npjglo          )

        ! get e3 at level jk ( ps...)
        IF ( lfull ) THEN
           e3(:,:) = e31d(jk)
        ELSE
           e3(:,:) = getvar(cn_fzgr, 'e3t_ps', jk, npiglo, npjglo, ldiom=.TRUE.)
        ENDIF

        !  e3 is used as a flag for the mixed layer; It is 0 outside the mixed layer
        e3(:,:)=MAX ( 0., MIN(e3, zmxl-gdepw(jk) ) )
        WHERE ( e3 == 0 ) zmask = 0.

        dvol      = SUM( DBLE(e3 * zmask) )
        dmxlsaltc = zs * e3 * zmask * 1.d0 + dmxlsaltc

        IF (dvol /= 0 )THEN
           !   go on !
        ELSE
           !   no more layer below !
           EXIT   ! get out of the jk loop
        ENDIF

     END DO


     ! Output to netcdf file : Kg/m2
     dmxlsaltc = rprho0*dmxlsaltc
     ierr = putvar(ncout, id_varout(1), REAL(dmxlsaltc), 1, npiglo, npjglo, ktime=jt)
  END DO


  ierr = closeout(ncout)

END PROGRAM cdfmxlsaltc
