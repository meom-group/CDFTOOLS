PROGRAM cdfconvert
  !!======================================================================
  !!                     ***  PROGRAM  cdfconvert  ***
  !!=====================================================================
  !!  ** Purpose : Convert a set of dimgfile (Clipper like)
  !!               to a set of CDF files (Drakkar like )
  !!
  !!  ** Method  : Read tag then open the respective T S 2D U V files to create
  !!              gridT, gridU and gridV files.
  !!              Requires  mesh_hgr.nc and mesh_zgr.nc  files
  !!
  !! History : 2.1  : 01/2007  : J.M. Molines : Original code
  !!           3.0  : 05/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!  isdirect       : integer function which return the record length
  !!                   of the file in argument if a dimgfile, 0 else.
  !!                  
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

  INTEGER(KIND=4)                            :: ji, jj, jk      ! dummy loop index
  INTEGER(KIND=4)                            :: jt, jvar        ! dummy loop index
  INTEGER(KIND=4)                            :: narg, iargc     ! command line
  INTEGER(KIND=4)                            :: nvar            ! number of output variables
  INTEGER(KIND=4)                            :: npiglo, npjglo  ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt        ! size of the domain
  INTEGER(KIND=4)                            :: irecl, ii, ndim ! dimg stuff variables
  INTEGER(KIND=4)                            :: numu=10         ! logical id for input dimg file
  INTEGER(KIND=4)                            :: numv=11         !      "             "
  INTEGER(KIND=4)                            :: numt=12         !      "             "
  INTEGER(KIND=4)                            :: nums=14         !      "             "
  INTEGER(KIND=4)                            :: num2d=15        !      "             "
  INTEGER(KIND=4)                            :: numssh=16       !      "             "
  INTEGER(KIND=4)                            :: numuu=17        !      "             "
  INTEGER(KIND=4)                            :: numvv=18        !      "             "
  INTEGER(KIND=4)                            :: ncout           ! ncid of output netcdf file
  INTEGER(KIND=4)                            :: ierr            ! error status
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout  ! outpur variables levels and id's

  REAL(KIND=4)                               :: x1, y1          ! dimg header ( SW corner)
  REAL(KIND=4)                               :: dx, dy          ! dimg header ( x,y step)
  REAL(KIND=4)                               :: zspval          ! dimg header ( special value)
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: v2d, glam, gphi ! working arrays
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: zdep            ! depth
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: tim             ! time counter

  CHARACTER(LEN=256)                         :: cf_ufil         ! output gridU file
  CHARACTER(LEN=256)                         :: cf_vfil         ! output gridV file
  CHARACTER(LEN=256)                         :: cf_tfil         ! output gridT file
  CHARACTER(LEN=256)                         :: cf_bsfil        ! output BSF file
  CHARACTER(LEN=256)                         :: cf_dimgu        ! input dimg U file
  CHARACTER(LEN=256)                         :: cf_dimgv        ! input dimg V file
  CHARACTER(LEN=256)                         :: cf_dimgt        ! input dimg T file
  CHARACTER(LEN=256)                         :: cf_dimgs        ! input dimg S file
  CHARACTER(LEN=256)                         :: cf_dimg2d       ! input dimg 2D file
  CHARACTER(LEN=256)                         :: cf_dimguu       ! input dimg U2 file
  CHARACTER(LEN=256)                         :: cf_dimgvv       ! input dimg V2 file
  CHARACTER(LEN=256)                         :: cf_dimgssh      ! input dimg SSH file
  CHARACTER(LEN=256)                         :: ctag            ! time tag
  CHARACTER(LEN=256)                         :: confcase        ! config-case
  CHARACTER(LEN=80 )                         :: cheader         ! comment in header of dimg file
  CHARACTER(LEN=4  )                         :: cver            ! dimg version

  TYPE(variable), DIMENSION(:),  ALLOCATABLE :: stypvar         ! output data structure

  LOGICAL                                    :: lexist          ! flag for existing file
  LOGICAL                                    :: lchk = .FALSE.  ! flag for missing files
  !!----------------------------------------------------------------------

  !!  Read command line
  narg= iargc()
  IF ( narg /= 2 ) THEN
     PRINT *,' usage : cdfconvert CLIPPER_tag CLIPPER_Confcase'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Convert dimg files (CLIPPER like) to netcdf (DRAKKAR like).'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       CLIPPER_tag      : a string such as y2000m01d15 for time identification.' 
     PRINT *,'       CLIPPER_confcase : CONFIG-CASE of the files to be converted (eg ATL6-V6)'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cn_fhgr),' and ', TRIM(cn_fzgr)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : gridT, gridU, gridV files'
     PRINT *,'         variables : same as in standard NEMO output'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfflxconv, cdfsstconv, cdfstrconv'
     PRINT *,'      '
     STOP
  ENDIF
  !!
  CALL getarg (1, ctag)
  CALL getarg (2, confcase)

  lchk = lchk .OR. chkfile( cn_fhgr )
  lchk = lchk .OR. chkfile (cn_fzgr )

  !! Build dimg file names
  cf_dimgu   = TRIM(confcase)//'_U_'  //TRIM(ctag)//'.dimg' ; lchk = lchk .OR. chkfile(cf_dimgu )
  cf_dimgv   = TRIM(confcase)//'_V_'  //TRIM(ctag)//'.dimg' ; lchk = lchk .OR. chkfile(cf_dimgv )
  cf_dimgt   = TRIM(confcase)//'_T_'  //TRIM(ctag)//'.dimg' ; lchk = lchk .OR. chkfile(cf_dimgt )
  cf_dimgs   = TRIM(confcase)//'_S_'  //TRIM(ctag)//'.dimg' ; lchk = lchk .OR. chkfile(cf_dimgs )
  cf_dimg2d  = TRIM(confcase)//'_2D_' //TRIM(ctag)//'.dimg' ; lchk = lchk .OR. chkfile(cf_dimg2d)
  IF ( lchk ) STOP ! missing file

  cf_dimgssh = TRIM(confcase)//'_SSH_'//TRIM(ctag)//'.dimg'
  cf_dimguu  = TRIM(confcase)//'_UU_' //TRIM(ctag)//'.dimg'
  cf_dimgvv  = TRIM(confcase)//'_VV_' //TRIM(ctag)//'.dimg'

  cf_ufil    = TRIM(confcase)//'_'    //TRIM(ctag)//'_gridU.nc'
  cf_vfil    = TRIM(confcase)//'_'    //TRIM(ctag)//'_gridV.nc'
  cf_tfil    = TRIM(confcase)//'_'    //TRIM(ctag)//'_gridT.nc'
  cf_bsfil   = TRIM(confcase)//'_'    //TRIM(ctag)//'_PSI.nc'

  ! open (and check ?? if they exists )
  irecl=isdirect(cf_dimgu ) ; OPEN( numu,  FILE=cf_dimgu,  FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
  irecl=isdirect(cf_dimgv ) ; OPEN( numv,  FILE=cf_dimgv,  FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
  irecl=isdirect(cf_dimgt ) ; OPEN( numt,  FILE=cf_dimgt,  FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
  irecl=isdirect(cf_dimgs ) ; OPEN( nums,  FILE=cf_dimgs,  FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
  irecl=isdirect(cf_dimg2d) ; OPEN( num2d, FILE=cf_dimg2d, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )

  READ(numt,REC=1) cver, cheader, ii, npiglo, npjglo, npk, npt

  ALLOCATE (v2d(npiglo, npjglo), glam(npiglo,npjglo), gphi(npiglo,npjglo), zdep(npk), tim(npt) )

  READ(numt,REC=1) cver, cheader, ii, npiglo, npjglo, npk, npt, ndim, &
       &      x1,y1,dx,dy,zspval, &
       &    ( zdep(jk),jk=1,npk), &
       ( tim(jt), jt=1,npt)

  ! transform Clipper days to drakkar seconds ...
  tim(:)=tim(:)*86400.

  !###############
  !# GRID T FILE #
  !###############
  ! Build gridT file with votemper, vosaline, sossheig, ... fluxes ...
  INQUIRE(FILE=cf_dimgssh, EXIST=lexist)
  IF ( lexist ) THEN
     irecl = isdirect(cf_dimgssh) 
     OPEN( numssh,FILE=cf_dimgssh, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
     nvar = 10 
  ELSE
     nvar = 9
  ENDIF

  ALLOCATE ( stypvar(nvar), ipk(nvar), id_varout(nvar) )
  jvar=1 
  ipk(jvar)                       = npk       
  stypvar(jvar)%cname             = cn_votemper
  stypvar(jvar)%cunits            = 'C'      
  stypvar(jvar)%rmissing_value    = 0.      
  stypvar(jvar)%valid_min         = -2.    
  stypvar(jvar)%valid_max         = 40.   
  stypvar(jvar)%clong_name        = 'Potential Temperature'
  stypvar(jvar)%cshort_name       = cn_votemper
  stypvar(jvar)%conline_operation = 'N/A'                
  stypvar(jvar)%caxis             = 'TZYX'              
  jvar=jvar+1

  ipk(jvar)                       = npk
  stypvar(jvar)%cname             = cn_vosaline
  stypvar(jvar)%cunits            = 'PSU'
  stypvar(jvar)%rmissing_value    = 0.
  stypvar(jvar)%valid_min         = 0.
  stypvar(jvar)%valid_max         = 45.
  stypvar(jvar)%clong_name        = 'Salinity'
  stypvar(jvar)%cshort_name       = cn_vosaline
  stypvar(jvar)%conline_operation = 'N/A'
  stypvar(jvar)%caxis             = 'TZYX'
  jvar=jvar+1

  IF ( lexist ) THEN
     ipk(jvar)                       = 1
     stypvar(jvar)%cname             = cn_sossheig
     stypvar(jvar)%cunits            = 'm'
     stypvar(jvar)%rmissing_value    = 0.
     stypvar(jvar)%valid_min         = -10.
     stypvar(jvar)%valid_max         = 10.
     stypvar(jvar)%clong_name        = 'Sea_Surface_height'
     stypvar(jvar)%cshort_name       = cn_sossheig
     stypvar(jvar)%conline_operation = 'N/A'
     stypvar(jvar)%caxis             = 'TYX'
     jvar=jvar+1
  ENDIF

  ipk(jvar)                       = 1
  stypvar(jvar)%cname             = cn_somxl010         ! rec 12 of dimg file 2D
  stypvar(jvar)%cunits            = 'm'
  stypvar(jvar)%rmissing_value    = 0.
  stypvar(jvar)%valid_min         = 0.
  stypvar(jvar)%valid_max         = 7000.
  stypvar(jvar)%clong_name        = 'Mixed_Layer_Depth_on_0.01_rho_crit'
  stypvar(jvar)%cshort_name       = cn_somxl010
  stypvar(jvar)%conline_operation = 'N/A'
  stypvar(jvar)%caxis             = 'TYX'
  jvar=jvar+1

  ipk(jvar)                       = 1
  stypvar(jvar)%cname             = 'sohefldo'         ! rec 4 of dimg file 2D
  stypvar(jvar)%cunits            = 'W/m2'
  stypvar(jvar)%rmissing_value    = 0.
  stypvar(jvar)%valid_min         = -1000.
  stypvar(jvar)%valid_max         = 1000.
  stypvar(jvar)%clong_name        = 'Net_Downward_Heat_Flux'
  stypvar(jvar)%cshort_name       = 'sohefldo'
  stypvar(jvar)%conline_operation = 'N/A'
  stypvar(jvar)%caxis             = 'TYX'
  jvar=jvar+1

  ipk(jvar)                       = 1
  stypvar(jvar)%cname             = cn_soshfldo         ! rec 8 of dimg file 2D (qsr)
  stypvar(jvar)%cunits            = 'W/m2'
  stypvar(jvar)%rmissing_value    = 0.
  stypvar(jvar)%valid_min         = -1000.
  stypvar(jvar)%valid_max         = 1000.
  stypvar(jvar)%clong_name        = 'Short_Wave_Radiation'
  stypvar(jvar)%cshort_name       = cn_soshfldo
  stypvar(jvar)%conline_operation = 'N/A'
  stypvar(jvar)%caxis             = 'TYX'
  jvar=jvar+1

  ipk(jvar)                       = 1
  stypvar(jvar)%cname             = cn_sowaflup       ! rec 5 of dimg file 2D (emp)
  stypvar(jvar)%cunits            = 'kg/m2/s'         ! conversion required from CLIPPER /86400.
  stypvar(jvar)%rmissing_value    = 0.
  stypvar(jvar)%valid_min         = -1000.
  stypvar(jvar)%valid_max         = 1000.
  stypvar(jvar)%clong_name        = 'Net_Upward_Water_Flux'
  stypvar(jvar)%cshort_name       = cn_sowaflup
  stypvar(jvar)%conline_operation = 'N/A'
  stypvar(jvar)%caxis             = 'TYX'
  jvar=jvar+1

  ipk(jvar)                       = 1
  stypvar(jvar)%cname             = 'sowafldp'       ! rec 10 of dimg file 2D (erp)
  stypvar(jvar)%cunits            = 'kg/m2/s'         ! conversion required from CLIPPER /jvar.
  stypvar(jvar)%rmissing_value    = 0.
  stypvar(jvar)%valid_min         = -1000.
  stypvar(jvar)%valid_max         = 1000.
  stypvar(jvar)%clong_name        = 'Surface_Water_Flux:Damping'
  stypvar(jvar)%cshort_name       = 'sowafldp'
  stypvar(jvar)%conline_operation = 'N/A'
  stypvar(jvar)%caxis             = 'TYX'
  jvar=jvar+1

  ipk(jvar)                       = 1
  stypvar(jvar)%cname             = cn_soicecov         ! rec 13 of dimg file 2D (erp)
  stypvar(jvar)%cunits            = '%'         
  stypvar(jvar)%rmissing_value    = 0.
  stypvar(jvar)%valid_min         = 0.
  stypvar(jvar)%valid_max         = 1.
  stypvar(jvar)%clong_name        = 'Ice Cover'
  stypvar(jvar)%cshort_name       = cn_soicecov
  stypvar(jvar)%conline_operation = 'N/A'
  stypvar(jvar)%caxis             = 'TYX'
  jvar=jvar+1

  ipk(jvar)                       = 1
  stypvar(jvar)%cname             = 'sohefldp'         ! rec 9 of dimg file 2D (erp)
  stypvar(jvar)%cunits            = 'W/m2'         
  stypvar(jvar)%rmissing_value    = 0.
  stypvar(jvar)%valid_min         = -10.
  stypvar(jvar)%valid_max         = 10.
  stypvar(jvar)%clong_name        = 'Surface Heat Flux: Damping'
  stypvar(jvar)%cshort_name       = 'sohefldp'
  stypvar(jvar)%conline_operation = 'N/A'
  stypvar(jvar)%caxis             = 'TYX'

  glam = getvar  (cn_fhgr, cn_glamt, 1, npiglo, npjglo)
  gphi = getvar  (cn_fhgr, cn_gphit, 1, npiglo, npjglo)
  zdep = getvare3(cn_fzgr, cn_gdept, npk              )

  ncout = create      (cf_tfil, 'none',  npiglo, npjglo, npk, cdep=cn_vdeptht                       )
  ierr  = createvar   (ncout,   stypvar, nvar,   ipk,    id_varout                                  )
  ierr  = putheadervar(ncout,   'none',  npiglo, npjglo, npk, pnavlon=glam, pnavlat=gphi, pdep=zdep )

  jvar=1
  ! T
  DO jk=1, npk
     READ(numt,REC=jk+1) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
     ierr = putvar(ncout, id_varout(jvar), v2d, jk, npiglo, npjglo)
  END DO
  jvar  = jvar+1
  PRINT *, 'Done for T'

  ! S
  DO jk=1, npk
     READ(nums,REC=jk+1) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
     ierr = putvar(ncout, id_varout(jvar), v2d, jk, npiglo, npjglo)
  END DO
  jvar  = jvar+1
  PRINT *, 'Done for S'

  IF ( lexist ) THEN
     ! SSH
     READ(numssh,REC=2) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
     ierr = putvar(ncout, id_varout(jvar), v2d, 1, npiglo, npjglo)
     jvar = jvar+1
     PRINT *, 'Done for SSH'
  ENDIF

  ! MXL
  READ(num2d,REC=12) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  ierr = putvar(ncout, id_varout(jvar), v2d, 1, npiglo, npjglo)
  jvar = jvar+1
  PRINT *, 'Done for MXL'

  ! QNET
  READ(num2d,REC=4 ) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  ierr = putvar(ncout, id_varout(jvar), v2d, 1, npiglo, npjglo)
  jvar = jvar+1
  PRINT *, 'Done for QNET'

  ! QSR
  READ(num2d,REC=8)  (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  ierr = putvar(ncout, id_varout(jvar), v2d, 1, npiglo, npjglo)
  jvar = jvar+1
  PRINT *, 'Done for QSR'

  ! EMP
  READ(num2d,REC=5)  (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  v2d  = v2d/86400. ! to change units
  ierr = putvar(ncout, id_varout(jvar), v2d, 1, npiglo, npjglo)
  jvar = jvar+1
  PRINT *, 'Done for EMP'

  ! ERP
  READ(num2d,REC=10) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  v2d  = v2d/86400. ! to change units
  ierr = putvar(ncout, id_varout(jvar), v2d, 1, npiglo, npjglo)
  jvar = jvar+1
  PRINT *, 'Done for ERP'

  ! FREEZE
  READ(num2d,REC=13) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  ierr = putvar(ncout, id_varout(jvar), v2d, 1, npiglo, npjglo)
  jvar = jvar+1
  PRINT *, 'Done for FREEZE'

  ! QRP
  READ(num2d,REC=9)  (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  ierr = putvar(ncout, id_varout(jvar), v2d, 1, npiglo, npjglo)
  jvar = jvar+1
  PRINT *, 'Done for QRP'

  ierr = putvar1d(ncout, tim, npt, 'T')
  ierr = closeout(ncout)
  DEALLOCATE ( stypvar, ipk, id_varout )


  !###############
  !# GRID U FILE #
  !###############
  ! Build gridU file with vozocrtx, sozotaux
  INQUIRE(FILE=cf_dimguu, EXIST=lexist)
  IF ( lexist ) THEN
     irecl = isdirect(cf_dimguu)
     OPEN( numuu, FILE=cf_dimguu, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
     nvar=3 
  ELSE
     nvar=2
  ENDIF

  ALLOCATE ( stypvar(nvar), ipk(nvar), id_varout(nvar) )

  jvar = 1
  ipk(jvar)                       = npk
  stypvar(jvar)%cname             = cn_vozocrtx
  stypvar(jvar)%cunits            = 'm/s'
  stypvar(jvar)%rmissing_value    = 0.
  stypvar(jvar)%valid_min         = 0.
  stypvar(jvar)%valid_max         = 20.
  stypvar(jvar)%clong_name        = 'Zonal Velocity '
  stypvar(jvar)%cshort_name       = cn_vozocrtx
  stypvar(jvar)%conline_operation = 'N/A'
  stypvar(jvar)%caxis             = 'TZYX'
  jvar = jvar+1

  ipk(jvar)                       =  1
  stypvar(jvar)%cname             = 'sozotaux'
  stypvar(jvar)%cunits            = 'N/m2'
  stypvar(jvar)%rmissing_value    = 0.
  stypvar(jvar)%valid_min         = 0.
  stypvar(jvar)%valid_max         = 20.
  stypvar(jvar)%clong_name        = 'Zonal Wind Stress'
  stypvar(jvar)%cshort_name       = 'sozotaux'
  stypvar(jvar)%conline_operation = 'N/A'
  stypvar(jvar)%caxis             = 'TYX'
  jvar = jvar+1

  IF ( lexist ) THEN
     ipk(jvar)      = npk
     stypvar(jvar)%cname             = TRIM(cn_vozocrtx)//'_sqd'
     stypvar(jvar)%cunits            = 'm2/s2'
     stypvar(jvar)%rmissing_value    = 0.
     stypvar(jvar)%valid_min         = 0.
     stypvar(jvar)%valid_max         = 100.
     stypvar(jvar)%clong_name        = 'MS_Zonal_Velocity'
     stypvar(jvar)%cshort_name       = TRIM(cn_vozocrtx)//'_sqd'
     stypvar(jvar)%conline_operation = 'N/A'
     stypvar(jvar)%caxis             = 'TZYX'
  ENDIF

  glam = getvar  (cn_fhgr, cn_glamu, 1, npiglo, npjglo)
  gphi = getvar  (cn_fhgr, cn_gphiu, 1, npiglo, npjglo)
  zdep = getvare3(cn_fzgr, cn_gdept, npk              )

  ncout = create      (cf_ufil, 'none',  npiglo, npjglo, npk, cdep=cn_vdepthu                       )
  ierr  = createvar   (ncout,   stypvar, nvar,   ipk,    id_varout                                  )
  ierr  = putheadervar(ncout,   'none',  npiglo, npjglo, npk, pnavlon=glam, pnavlat=gphi, pdep=zdep )

  jvar=1
  DO jk=1, npk
     READ(numu,REC=jk+1) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
     ierr = putvar(ncout, id_varout(jvar), v2d, jk, npiglo, npjglo)
  END DO
  jvar  = jvar+1
  PRINT *, 'Done for U'

  READ(num2d, REC=2 ) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  ierr = putvar(ncout, id_varout(jvar), v2d, 1,  npiglo, npjglo)
  jvar = jvar+1
  PRINT *, 'Done for TAUX'

  IF ( lexist ) THEN
     DO jk=1, npk
        READ(numuu,REC=jk+1) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
        ierr = putvar(ncout, id_varout(jvar), v2d, jk,  npiglo, npjglo)
     END DO
     PRINT *, 'Done for UU'
  ENDIF

  ierr = putvar1d(ncout, tim, npt, 'T')
  ierr = closeout(ncout               )

  DEALLOCATE ( stypvar, ipk, id_varout )

  !###############
  !# GRID V FILE #
  !###############
  ! Build gridV file with vomecrty, sometauy
  INQUIRE(FILE=cf_dimgvv, EXIST=lexist)
  IF ( lexist ) THEN
     irecl = isdirect(cf_dimgvv)
     OPEN( numvv, FILE=cf_dimgvv, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
     nvar=3 
  ELSE
     nvar=2
  ENDIF
  ALLOCATE ( stypvar(nvar), ipk(nvar), id_varout(nvar) )

  jvar=1
  ipk(jvar)                       = npk
  stypvar(jvar)%cname             = cn_vomecrty
  stypvar(jvar)%cunits            = 'm/s'
  stypvar(jvar)%rmissing_value    = 0.
  stypvar(jvar)%valid_min         = 0.
  stypvar(jvar)%valid_max         = 20.
  stypvar(jvar)%clong_name        = 'Meridinal  Velocity '
  stypvar(jvar)%cshort_name       = cn_vomecrty
  stypvar(jvar)%conline_operation = 'N/A'
  stypvar(jvar)%caxis             = 'TZYX'
  jvar = jvar+1

  ipk(jvar)                       = 1
  stypvar(jvar)%cname             = 'sometauy'
  stypvar(jvar)%cunits            = 'N/m2'
  stypvar(jvar)%rmissing_value    = 0.
  stypvar(jvar)%valid_min         = 0.
  stypvar(jvar)%valid_max         = 20.
  stypvar(jvar)%clong_name        = 'Meridional Wind Stress'
  stypvar(jvar)%cshort_name       = 'sometauy'
  stypvar(jvar)%conline_operation = 'N/A'
  stypvar(jvar)%caxis             = 'TYX'
  jvar=jvar+1

  IF ( lexist ) THEN
     ipk(jvar)                       = npk
     stypvar(jvar)%cname             = TRIM(cn_vomecrty)//'_sqd'
     stypvar(jvar)%cunits            = 'm2/s2'
     stypvar(jvar)%rmissing_value    = 0.
     stypvar(jvar)%valid_min         = 0.
     stypvar(jvar)%valid_max         = 100.
     stypvar(jvar)%clong_name        = 'MS_Meridional_Velocity'
     stypvar(jvar)%cshort_name       = TRIM(cn_vomecrty)//'_sqd'
     stypvar(jvar)%conline_operation = 'N/A'
     stypvar(jvar)%caxis             = 'TZYX'
  ENDIF


  glam = getvar  (cn_fhgr, cn_glamv, 1,  npiglo, npjglo)
  gphi = getvar  (cn_fhgr, cn_gphiv, 1,  npiglo, npjglo)
  zdep = getvare3(cn_fzgr, cn_gdept, npk               )

  ncout = create      (cf_vfil, 'none',  npiglo, npjglo, npk, cdep=cn_vdepthv                       )
  ierr  = createvar   (ncout,   stypvar, nvar,   ipk,    id_varout                                  )
  ierr  = putheadervar(ncout,   'none',  npiglo, npjglo, npk, pnavlon=glam, pnavlat=gphi, pdep=zdep )

  jvar = 1
  DO jk=1, npk
     READ(numv,REC=jk+1)  (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
     ierr = putvar (ncout, id_varout(jvar), v2d, jk, npiglo, npjglo)
  END DO
  jvar  = jvar+1
  PRINT *, 'Done for V'

  READ(num2d, REC=3) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  ierr = putvar(ncout, id_varout(jvar), v2d, 1, npiglo, npjglo)
  jvar = jvar+1
  PRINT *, 'Done for TAUY'

  IF ( lexist ) THEN
     DO jk=1, npk
        READ(numvv,REC=jk+1) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
        ierr = putvar(ncout, id_varout(jvar), v2d, jk,  npiglo, npjglo)
     END DO
     PRINT *, 'Done for VV'
  ENDIF

  ierr = putvar1d(ncout, tim, npt, 'T')
  ierr = closeout(ncout               )

  DEALLOCATE ( stypvar, ipk, id_varout )

  !###############
  !# PSI FILE #
  !###############
  ! Build PSI file with sobarstf
  nvar=1  
  ALLOCATE ( stypvar(nvar), ipk(nvar), id_varout(nvar) )
  ipk(1)                       = 1
  stypvar(1)%cname             = 'sobarstf'
  stypvar(1)%cunits            = 'm3/s'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = -3.e8
  stypvar(1)%valid_max         = 3.e8
  stypvar(1)%clong_name        = 'Barotropic_Stream_Function'
  stypvar(1)%cshort_name       = 'sobarstf'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TYX'

  glam = getvar  (cn_fhgr, cn_glamf, 1, npiglo, npjglo)
  gphi = getvar  (cn_fhgr, cn_gphif, 1, npiglo, npjglo)
  zdep = getvare3(cn_fzgr, cn_gdept, 1                )

  ncout = create      (cf_bsfil, 'none',  npiglo, npjglo, 1, cdep=cn_vdepthu                       )
  ierr  = createvar   (ncout,    stypvar, nvar,   ipk,    id_varout                                )
  ierr  = putheadervar(ncout,    'none',  npiglo, npjglo, 1, pnavlon=glam, pnavlat=gphi, pdep=zdep )

  jvar = 1
  READ(num2d,REC=7) (( v2d(ji,jj), ji=1, npiglo), jj=1,npjglo)
  ierr = putvar(ncout, id_varout(jvar),v2d, 1, npiglo, npjglo)
  PRINT *, 'Done for PSI'

  ierr = putvar1d(ncout, tim, npt, 'T')
  ierr = closeout(ncout               )

  DEALLOCATE ( stypvar, ipk, id_varout )

CONTAINS

  INTEGER(KIND=4) FUNCTION isdirect(cdname)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION isdirect  ***
    !!
    !! ** Purpose :  This integer function returns the record length if cdname 
    !!               is a valid dimg file, it returns 0 either.
    !!
    !! ** Method  :  Open the file and look for the key characters (@!01) for
    !!               identification.
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdname

    ! --
    INTEGER(KIND=4)              :: irecl
    INTEGER(KIND=4)              :: inum = 100

    CHARACTER(LEN=4)             :: clver
    CHARACTER(LEN=80)            :: clheader
    !!----------------------------------------------------------------------

    !
    OPEN(inum,FILE=cdname, FORM = 'UNFORMATTED', ACCESS = 'DIRECT', RECL = 88)
    READ(inum,REC=1) clver ,clheader, irecl
    CLOSE(inum)
    !
    IF (clver ==  '@!01' ) THEN
       isdirect = irecl
    ELSE
       isdirect = 0
    END IF
    !
  END FUNCTION isdirect
END PROGRAM cdfconvert
