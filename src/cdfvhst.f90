PROGRAM cdfvhst
  !!======================================================================
  !!                     ***  PROGRAM  cdfvhst  ***
  !!=====================================================================
  !!  ** Purpose : Compute Verticaly integrated  Heat Salt Transport.
  !!
  !!  ** Method  : Take VT files computed by cdfvT.f90 and integrate
  !!               vertically to produce a 2D file
  !!
  !! History : 2.1  : 01/2005  : J.M. Molines : Original code
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
  !! @class transport
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: jk, jt          ! dummy loop index
  INTEGER(KIND=4)                           :: it              ! time index for vvl
  INTEGER(KIND=4)                           :: ierr            ! working integer
  INTEGER(KIND=4)                           :: narg, iargc     ! command line 
  INTEGER(KIND=4)                           :: ijarg           ! argument counter
  INTEGER(KIND=4)                           :: npiglo, npjglo  ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt        ! size of the domain
  INTEGER(KIND=4)                           :: ncout           ! ncdf id of output file
  INTEGER(KIND=4), DIMENSION(4)             :: ipk, id_varout  ! output variable levels and id's

  REAL(KIND=4), PARAMETER                   :: pp_rau0=1000.   ! fresh water density ( kg/m3)
  REAL(KIND=4), PARAMETER                   :: pp_rcp=4000.    ! heat capacity of water (J/kg/K)
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: e31d            ! vertical metrics when full step
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1v, e2u        ! horizontal metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e3u, e3v        ! vertical metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zut, zus        ! heat and salt zonal copmponents
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zvt, zvs        ! heat and salt meridional components

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim            ! time counter
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dtrput, dtrpus  ! zonal transport
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dtrpvt, dtrpvs  ! meridional transport

  TYPE (variable), DIMENSION(4)             :: stypvar         ! structure output variables

  CHARACTER(LEN=256)                        :: cf_vtfil        ! input file name (vt)
  CHARACTER(LEN=256)                        :: cf_out='trp.nc' ! output file name
  CHARACTER(LEN=256)                        :: cldum           ! dummy char variable

  LOGICAL                                   :: lfull   =.FALSE.! flag for full step
  LOGICAL                                   :: lnc4    =.FALSE.! flag for netcdf4
  LOGICAL                                   :: lchk    =.FALSE.! flag for checking files
  LOGICAL                                   :: lchkvar =.FALSE.! flag for missing variables
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfvhst  -f VT-file [-full] [-vvl] [-o OUT-file] [-nc4] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'         Compute the vertically integrated heat and salt transports '
     PRINT *,'         at each grid cell. In other words, it is the vertical integral'
     PRINT *,'         of the VT-file.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'         -f VT-file : file which contains UT, VT, US, VS quantities'
     PRINT *,'              (produced by cdfvT.f90)'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'         [-full ] : use full step computation (default is partial steps).'
     PRINT *,'         [-vvl  ] : use time-varying vertical metrics.'
     PRINT *,'         [-o OUT-file ] : specify output file name, instead of ', TRIM(cf_out)
     PRINT *,'         [-nc4  ] : use netcdf4 chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'         Files ',TRIM(cn_fhgr),', ',TRIM(cn_fzgr)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'         Netcdf file : ',TRIM(cf_out), ' unless -o option is used.'
     PRINT *,'         Variables : ', TRIM(cn_somevt),', ',TRIM(cn_somevs),', ',TRIM(cn_sozout),' and  ',TRIM(cn_sozous)
     STOP 
  ENDIF

  ijarg = 1
  DO WHILE (ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg = ijarg+1
     SELECT CASE (cldum )
     CASE ( '-f'    ) ;  CALL getarg(ijarg, cf_vtfil) ; ijarg = ijarg+1
        ! options
     CASE ( '-full' ) ;  lfull  = .TRUE.
     CASE ( '-vvl'  ) ;  lg_vvl = .TRUE.
     CASE ( '-o'    ) ;  CALL getarg(ijarg, cf_out  ) ; ijarg = ijarg+1
     CASE ( '-nc4'  ) ;  lnc4   = .TRUE.
     CASE DEFAULT     ;  PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.'  ; STOP 99
     END SELECT
  END DO

  lchk = lchk .AND. chkfile(cn_fhgr  )
  lchk = lchk .AND. chkfile(cn_fzgr  )
  lchk = lchk .AND. chkfile(cf_vtfil )

  IF ( lchk  ) STOP 99 ! missing file
  IF ( lg_vvl) THEN 
     cn_fe3u = cf_vtfil 
     cn_fe3v = cf_vtfil 
     cn_ve3u = cn_ve3uvvl
     cn_ve3v = cn_ve3vvvl
     lchkvar = lchkvar .AND.chkvar( cn_fe3u, cn_ve3u)
     lchkvar = lchkvar .AND.chkvar( cn_fe3v, cn_ve3v)
     IF ( lchkvar ) THEN ; PRINT *,'no vertical metrics for vvl' ; STOP 99 ! missing e3 metrics in VT file 
     ENDIF
  ENDIF

  npiglo= getdim (cf_vtfil,cn_x )
  npjglo= getdim (cf_vtfil,cn_y )
  npk   = getdim (cf_vtfil,cn_z )
  npt   = getdim (cf_vtfil,cn_t )

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ! Allocate arrays
  ALLOCATE ( zvt(npiglo,npjglo), zvs(npiglo,npjglo) )
  ALLOCATE ( zut(npiglo,npjglo), zus(npiglo,npjglo) )
  ALLOCATE ( e1v(npiglo,npjglo), e3v(npiglo,npjglo) )
  ALLOCATE ( e2u(npiglo,npjglo), e3u(npiglo,npjglo) )
  ALLOCATE ( dtrpvt(npiglo,npjglo), dtrpvs(npiglo,npjglo))
  ALLOCATE ( dtrput(npiglo,npjglo), dtrpus(npiglo,npjglo))
  ALLOCATE ( dtim(npt), e31d(npk) )

  CALL CreateOutput

  ! read level independent metrics
  e1v(:,:) = getvar(cn_fhgr,   cn_ve1v, 1, npiglo, npjglo)
  e2u(:,:) = getvar(cn_fhgr,   cn_ve2u, 1, npiglo, npjglo)
  e31d(:)  = getvare3(cn_fzgr, cn_ve3t1d, npk              ) ! used only for full step

  DO jt=1, npt
     IF ( lg_vvl ) THEN ; it =jt
     ELSE               ; it = 1
     ENDIF
     ! reset transport to 0
     dtrpvt(:,:) = 0.d0 ; dtrpvs(:,:) = 0.d0 ; dtrput(:,:) = 0.d0 ; dtrpus(:,:) = 0.d0

     DO jk = 1,npk
        PRINT *,'level ',jk, ' time ', jt
        ! Get heat/salt transport component at jk
        zvt(:,:)= getvar(cf_vtfil, cn_vomevt, jk ,npiglo, npjglo, ktime=jt)
        zvs(:,:)= getvar(cf_vtfil, cn_vomevs, jk ,npiglo, npjglo, ktime=jt)
        zut(:,:)= getvar(cf_vtfil, cn_vozout, jk ,npiglo, npjglo, ktime=jt)
        zus(:,:)= getvar(cf_vtfil, cn_vozous, jk ,npiglo, npjglo, ktime=jt)

        ! get e3v at level jk ( and multiply by respective horizontal metric)
        IF ( lfull ) THEN
           e3v(:,:) = e31d(jk) * e1v(:,:)
           e3u(:,:) = e31d(jk) * e2u(:,:)
        ELSE
           e3v(:,:) = getvar(cn_fe3v, cn_ve3v, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl) * e1v(:,:)
           e3u(:,:) = getvar(cn_fe3u, cn_ve3u, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl) * e2u(:,:)
        ENDIF

        ! integrates vertically 
        dtrpvt(:,:) = dtrpvt(:,:) + zvt(:,:) * e3v(:,:) * pp_rau0*pp_rcp * 1.d0
        dtrpvs(:,:) = dtrpvs(:,:) + zvs(:,:) * e3v(:,:) * 1.d0
        dtrput(:,:) = dtrput(:,:) + zut(:,:) * e3u(:,:) * pp_rau0*pp_rcp * 1.d0
        dtrpus(:,:) = dtrpus(:,:) + zus(:,:) * e3u(:,:) * 1.d0

     END DO  ! loop to next level

     ! output on file
     ierr = putvar(ncout, id_varout(1) ,SNGL(dtrpvt), 1, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(2) ,SNGL(dtrpvs), 1, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(3) ,SNGL(dtrput), 1, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(4) ,SNGL(dtrpus), 1, npiglo, npjglo, ktime=jt)
  END DO  ! loop on time step

  ierr = closeout (ncout)

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
    ipk(:)                    = 1
    stypvar%rmissing_value    = 0.
    stypvar%valid_min         = -100.
    stypvar%valid_max         = 100.
    stypvar%conline_operation = 'N/A'
    stypvar%caxis             = 'TYX'

    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(3)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(4)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)

    stypvar(1)%cname       = cn_somevt
    stypvar(2)%cname       = cn_somevs
    stypvar(3)%cname       = cn_sozout
    stypvar(4)%cname       = cn_sozous

    stypvar(1)%cunits      = 'W'
    stypvar(2)%cunits      = 'kg.s-1'
    stypvar(3)%cunits      = 'W'
    stypvar(4)%cunits      = 'kg.s-1'

    stypvar(1)%clong_name  = 'Meridional_heat_transport'
    stypvar(2)%clong_name  = 'Meridional_salt_transport'
    stypvar(3)%clong_name  = 'Zonal_heat_transport'
    stypvar(4)%clong_name  = 'Zonal_salt_transport'

    stypvar(1)%cshort_name = cn_somevt
    stypvar(2)%cshort_name = cn_somevs
    stypvar(3)%cshort_name = cn_sozout
    stypvar(4)%cshort_name = cn_sozous

    ! create output fileset
    ncout = create      (cf_out, cf_vtfil, npiglo, npjglo, 1        , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar,  4,      ipk,    id_varout, ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_vtfil, npiglo, npjglo, 1         )

    dtim  = getvar1d(cf_vtfil, cn_vtimec, npt     )
    ierr  = putvar1d(ncout,    dtim,      npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfvhst
