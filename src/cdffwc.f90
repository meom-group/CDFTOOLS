PROGRAM cdffwc
  !!======================================================================
  !!                      ***  PROGRAM  cdffwc   ***
  !!=====================================================================
  !!  ** Purpose : Computes the freshwater content in a given basin from top
  !!               to bottom for each layer. Can handle full step configuration
  !!               using the -full option.'
  !!
  !!  ** Method  : compute fwc = sum(e1*e2*e3) - sum(S*e1*e2*e3)/Sref * btmask
  !!                 with reference salinity Sref (=34.7) and 
  !!                 the sub-basin mask btmask
  !!
  !!               based on cdfvertmean routine
  !!
  !! History : 0.1  : 09/2016  : M. Scheinert : First adaption
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2016
  !! $Id$
  !! Copyright (c) 2016, J.-M. Molines & Markus Scheinert (GEOMAR)
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

   INTEGER(KIND=4)                            :: jk, jvar, jvarin, jt ! dummy loop index
   INTEGER(KIND=4)                            :: ierr, ij, iko        ! working integer
   INTEGER(KIND=4)                            :: narg, iargc, ijarg   ! command line 
   INTEGER(KIND=4)                            :: npiglo, npjglo       ! size of the domain
   INTEGER(KIND=4)                            :: npk, npt             ! size of the domain
   INTEGER(KIND=4)                            :: nvars, ivar          ! variables in input
   INTEGER(KIND=4)                            :: nvaro=1              ! variables for output
   INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout       ! levels and varid's of output vars
   INTEGER(KIND=4)                            :: ncout 

   REAL(KIND=4), PARAMETER                    :: pprho0 = 1020.       ! water density (kg/m3)
   REAL(KIND=4), PARAMETER                    :: ppcp   = 4000.       ! calorific capacity (J/kg/m3)
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: e3t,area,ssh         ! vertical metric
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: zt                   ! working input variable
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: tmask                ! npiglo x npjglo
   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE  :: tim                  ! time counter
   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE  :: e31d                 ! vertical metrics in case of full step
   REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE  :: zmask              ! npiglo x npjglo
   REAL(KIND=4)                               :: rdep1, rdep2         ! depth counters
   REAL(KIND=4)                               :: tol  = 1.0           ! tolerance 
   REAL(KIND=4)                               :: sclf = 1.0           ! scale factor
   REAL(KIND=4), PARAMETER                    :: ppspval= 9999.99     ! missing value

   REAL(KIND=8)                               :: ds0=34.7             ! reference salinity
   REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dfwc                ! fwc. as 2dim to be consistent with putvar()
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE  :: dl_vint1, dl_vol2d   ! verticall int quantity         

   CHARACTER(LEN=2048)                         :: cf_in                ! input file
   CHARACTER(LEN=2048)                         :: cf_out='fwc.nc'      ! output file 
   CHARACTER(LEN=2048)                         :: cldum                ! dummy string for command line browsing
   CHARACTER(LEN=2048)                         :: cf_subbas='subbasins.nc'     ! subbasins file
   CHARACTER(LEN=2048)                         :: cv_cur               ! variable name

  CHARACTER(LEN=2048), DIMENSION(:), ALLOCATABLE :: cv_names           ! name of input variables
  CHARACTER(LEN=2048), DIMENSION(:), ALLOCATABLE :: cv_in              ! name of output variables



   LOGICAL                                    :: lfull  =.FALSE.      ! flag for full step computation
   LOGICAL                                    :: lchk   =.FALSE.      ! flag for missing files
   LOGICAL                                    :: lsref  =.FALSE.      ! flag for user reference salinity
   LOGICAL                                    :: laccum =.FALSE.      ! flag for accumulated fwc
   LOGICAL                                    :: lssh   =.FALSE.      ! flag for using ssh for top layer thickness

   TYPE(variable), DIMENSION(:), ALLOCATABLE  :: stypvarin            ! stucture for attributes (input)
   TYPE(variable), DIMENSION(:), ALLOCATABLE  :: stypvar              ! extension for attributes
   !!----------------------------------------------------------------------
   CALL ReadCdfNames()

   narg= iargc()
   IF ( narg == 0 ) THEN
      PRINT *,' usage : cdffwc IN-file BASIN-var1,var2,.. [-o OUT-file] [-sref REFSAL]'
      PRINT *,'                [-full] [-accum] [-ssh]'
      PRINT *,'      '
      PRINT *,'     PURPOSE :'
      PRINT *,'       Computes the freshwater content in a given basin from top'
      PRINT *,'       to bottom for each layer. Can handle full step configuration'
      PRINT *,'       using the -full option.'
      PRINT *,'      '
      PRINT *,'     ARGUMENTS :'
      PRINT *,'        IN-file            : netcdf input file.' 
      PRINT *,'        BASIN-var1,var2,.. : Comma separated list of sub-basin variables'
      PRINT *,'                             to process.'
      PRINT *,'        OUT-file           : use specified output file instead of <IN-var>.nc'
      PRINT *,'      '
      PRINT *,'     OPTIONS :'
      PRINT *,'        -full  : for full step computation ' 
      PRINT *,'        -accum : compute accumulated content from top to bottom' 
      PRINT *,'        -ssh   : take ssh into account for surface layer' 
      PRINT *,'        -sref  : reference salinity (= 34.7 by deafult)'
      PRINT *,'      '
      PRINT *,'     REQUIRED FILES :'
      PRINT *,'       ', TRIM(cn_fzgr),', ',TRIM(cn_fhgr),' and ',TRIM(cf_subbas) ,' and ',TRIM(cn_fmsk)
      PRINT *,'      '
      PRINT *,'     OUTPUT : '
      PRINT *,'       netcdf file :  fwc.nc (or specified with -o option)'
      PRINT *,'         variables :  fwc_BASIN, where BASIN was set by argument BASIN-var*'
      PRINT *,'                      (cAsE sensitive !)'
      PRINT *,'      '
      STOP
   ENDIF


   ! browse command line
   ijarg = 1   ; ij = 0
   DO WHILE ( ijarg <= narg ) 
      CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
      SELECT CASE ( cldum)
      CASE ( '-full' ) ; lfull  = .TRUE. 
      CASE ( '-o'    ) ; CALL getarg (ijarg, cf_out) ; ijarg = ijarg + 1
      CASE ( '-sref' ) ; lsref  = .TRUE. ; CALL getarg (ijarg, cldum) ; READ(cldum,*) ds0; ijarg = ijarg + 1
      CASE ( '-accum') ; laccum = .TRUE. 
      CASE ( '-ssh'  ) ; lssh   = .TRUE. 
      CASE DEFAULT     
         ij = ij + 1
         SELECT CASE ( ij)
         CASE ( 1 ) ; cf_in = cldum
         CASE ( 2 ) ; CALL ParseVars(cldum)
         CASE DEFAULT ; PRINT *, ' ERROR: Too many arguments ! ' ; STOP
         END SELECT
      END SELECT
   END DO

   ! Security check
   lchk = chkfile ( cf_in   )
   lchk = chkfile ( cn_fmsk ) .OR. lchk
   lchk = chkfile ( cn_fhgr ) .OR. lchk
   lchk = chkfile ( cn_fzgr ) .OR. lchk
   lchk = chkfile (cf_subbas) .OR. lchk
   IF ( lchk ) STOP ! missing files

   sclf      = 1.

   ! log information so far
   PRINT *,' OUTPUT FILE     : ' , TRIM(cf_out)

   npiglo = getdim (cf_in, cn_x )
   npjglo = getdim (cf_in, cn_y )
   npk    = getdim (cf_in, cn_z )
   npt    = getdim (cf_in, cn_t )

   PRINT *, ' NPIGLO = ', npiglo
   PRINT *, ' NPJGLO = ', npjglo
   PRINT *, ' NPK    = ', npk
   PRINT *, ' NPT    = ', npt

  nvars       = getnvar(cf_subbas)
  
  ALLOCATE(  cv_names( nvars       ) )
  ALLOCATE( stypvarin( nvars       ) )
  
  ALLOCATE(   stypvar( nvaro       ) )
  ALLOCATE(       ipk( nvaro       ) )
  ALLOCATE( id_varout( nvaro       ) )
  ALLOCATE(       dfwc( nvaro, 1, 1 ) )
  
  cv_names(:) = getvarname(cf_subbas, nvars, stypvarin)

  ! just chck if var exist in file 
  DO jvar = 1, nvaro
     IF ( chkvar( cf_subbas, cv_in(jvar)) ) STOP  ! message is written in cdfio.chkvar
  ENDDO

   
   ! Allocate arrays
   ALLOCATE (      tim( npt                    ) )
   ALLOCATE (    zmask( nvaro,   npiglo,npjglo ) )
   ALLOCATE (    tmask(          npiglo,npjglo ) )
   ALLOCATE (       zt(          npiglo,npjglo ) )
   ALLOCATE (      e3t(          npiglo,npjglo ) )
   ALLOCATE (     area(          npiglo,npjglo ) )
   ALLOCATE ( dl_vint1(          npiglo,npjglo ) )
   ALLOCATE ( dl_vol2d(          npiglo,npjglo ) )
   IF ( lssh  ) ALLOCATE ( ssh(  npiglo,npjglo ) )
   IF ( lfull ) ALLOCATE ( e31d( npk ) )

   ! prepare output variable
   ipk(:)                       = npk
   write(cldum, *) ds0


   DO jvar=1, nvaro     !Go through sub-basins
     DO jvarin=1,nvars
        IF ( cv_in(jvar) == stypvarin(jvarin)%cname ) EXIT  ! cv_in match cv_varin.
     END DO
     stypvar(jvar)%cname             = 'fwc_'//TRIM(cv_in(jvar))
     stypvar(jvar)%cunits            = 'km3'
     stypvar(jvar)%rmissing_value    = ppspval
     stypvar(jvar)%valid_min         = 0.0
     stypvar(jvar)%valid_max         = 0.0
     IF ( laccum ) THEN 
        stypvar(jvar)%clong_name        = 'freshwater content accumulated from top to bottom for '//TRIM(cv_in(jvar))//' based on S0='//TRIM(cldum)
     ELSE
        stypvar(jvar)%clong_name        = 'freshwater content per layer for '//TRIM(cv_in(jvar))//' based on S0='//TRIM(cldum)
     ENDIF
     stypvar(jvar)%clong_name        = 'freshwater content for '//TRIM(cv_in(jvar))//' based on S0='//TRIM(cldum)
     stypvar(jvar)%cshort_name       = 'fwc_'//TRIM(cv_in(jvar))
     stypvar(jvar)%conline_operation = 'N/A'
     stypvar(jvar)%caxis             = 'T'
   END DO

 
   ! Read area field
   area(:,:)   = getvar(cn_fhgr  , cn_ve1t,  1, npiglo, npjglo)  ! Save memory: e1t is first read into area
   zt(  :,:)   = getvar(cn_fhgr  , cn_ve2t,  1, npiglo, npjglo)  ! ..and zv is used to read e2t
   area(:,:)   = area(:,:) * zt(:,:)                             ! ..

   ! Read vertical axis in full-step-case
   IF ( lfull ) e31d(:) = getvare3(cn_fzgr, cn_ve3t, npk)

   ! Read sub-basin masks (2D) and keep them for all layers
   DO jvar = 1,nvaro          ! Loop through sub-basins
        zmask(jvar,:,:) = getvar(cf_subbas, cv_in(jvar) ,  1, npiglo, npjglo )
   END DO
   
   ! Initialize output file
   ncout = create      (cf_out, 'none', 1, 1, npk, cdep=cn_vdepthw, ld_xycoo=.FALSE. )
   ierr  = createvar   (ncout, stypvar, nvaro, ipk, id_varout)
   ierr  = putheadervar(ncout, cf_in,   1, 1, npk,     ld_xycoo=.FALSE.)
   
   tim   = getvar1d    (cf_in, cn_vtimec, npt     )
   ierr  = putvar1d    (ncout, tim,       npt, 'T')

   
   
   PRINT *, 'Output files initialised ...'


   DO jt=1,npt
        dl_vol2d(  :,:) = 0.d0
        dl_vint1(  :,:) = 0.d0
             dfwc(:,:,:) = 0.d0
        IF ( lssh ) ssh(:,:) = 0.e0

        DO jk = 1, npk
           ! Get values at jk
           zt(   :,:) = getvar(cf_in,     cn_vosaline, jk, npiglo, npjglo, ktime=jt)     ! Read Salinity(2D) at level jk and time step jt
           tmask(:,:) = getvar(cn_fmsk,       'tmask', jk, npiglo, npjglo          )
           
           ! get e3t at level jk ( ps...)
           IF ( lfull ) THEN ; e3t(:,:) = e31d(jk)
           ELSE              ; e3t(:,:) = getvar(cn_fzgr, 'e3t_ps', jk, npiglo, npjglo, ldiom=.TRUE.)
           ENDIF
           
           IF ( jk == 1 .AND. lssh ) THEN
                ssh(:,:) = getvar(cf_in,     cn_sossheig, 1, npiglo, npjglo, ktime=jt)
                e3t(:,:) = e3t(:,:) + ssh(:,:)
           ENDIF
           
           DO jvar = 1,nvaro          ! Loop through sub-basins
            
                dl_vol2d     = area * e3t * tmask * zmask(jvar,:,:) * 1.d0

                dl_vint1   = ( ds0 - zt ) / ds0 * dl_vol2d
                IF ( laccum ) THEN
                    dfwc(jvar,1,1)   = SUM( dl_vint1 ) + dfwc(jvar,1,1)
                ELSE
                    dfwc(jvar,1,1)   = SUM( dl_vint1 )
                ENDIF

                ! Output to netcdf file 
                ierr = putvar(ncout, id_varout(jvar), REAL(dfwc(jvar,:,:)), jk, 1, 1, ktime=jt)
            
           END DO   ! jvar
        
        END DO
        
   END DO  ! loop on time


   ierr = closeout(ncout)

  CONTAINS
  
  SUBROUTINE ParseVars (cdum)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ParseVars  ***
    !!
    !! ** Purpose :  Decode variable name  option from command line
    !!
    !! ** Method  :  look for , in the argument string and set the number of
    !!         variable (nvaro), allocate cv_in array and fill it with the
    !!         decoded  names.
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdum

    CHARACTER(LEN=2048), DIMENSION(100) :: cl_dum  ! 100 is arbitrary
    INTEGER  :: ji
    INTEGER  :: inchar,  i1=1
    !!----------------------------------------------------------------------
    inchar= LEN(TRIM(cdum))
    ! scan the input string and look for ',' as separator
    DO ji=1,inchar
       IF ( cdum(ji:ji) == ',' ) THEN
          cl_dum(nvaro) = cdum(i1:ji-1)
          i1=ji+1
          nvaro=nvaro+1
       ENDIF
    ENDDO

    ! last name of the list does not have a ','
    cl_dum(nvaro) = cdum(i1:inchar)

    ALLOCATE ( cv_in(nvaro) )
    DO ji=1, nvaro
       cv_in(ji) = cl_dum(ji)
    ENDDO
  END SUBROUTINE ParseVars


END PROGRAM cdffwc
