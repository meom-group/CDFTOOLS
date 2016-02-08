PROGRAM cdfvint
   !!======================================================================
   !!                     ***  PROGRAM  cdfvint  ***
   !!=====================================================================
   !!  ** Purpose : Compute vertically integrated temperature or salinity.
   !!
   !!  ** Method  : Compute the integral from top to bottom and save 
   !!               cumulated values. For temperature, cumulated values are
   !!               transformed to heat content (J.K.m^-2). For salinity
   !!               they are saved as PSU.m
   !!
   !! History : 2.1  : 10/2012  : M.A. Balmaseda : Original code from cdfmxlhc
   !!           3.0  : 11/2012  : J.M. Molines   : Doctor norm + Lic + ...
   !!                  02/2016  : S. Leroux      : Add -OCCI option
   !!----------------------------------------------------------------------
   USE cdfio
   USE modcdfnames
   !!----------------------------------------------------------------------
   !! CDFTOOLS_3.0 , MEOM 2012
   !! $Id$
   !! Copyright (c) 2012, J.-M. Molines
   !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
   !!----------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER(KIND=4)                           :: jk, jt              ! dummy loop index
   INTEGER(KIND=4)                           :: ierr, ij, iko       ! working integer
   INTEGER(KIND=4)                           :: narg, iargc, ijarg  ! command line 
   INTEGER(KIND=4)                           :: npiglo, npjglo      ! size of the domain
   INTEGER(KIND=4)                           :: npk, npt            ! size of the domain
   INTEGER(KIND=4)                           :: npko                ! size of the domain
   INTEGER(KIND=4),             DIMENSION(1) :: ipk, id_varout      ! only one output variable
   INTEGER(KIND=4)                           :: ncout 

   REAL(KIND=4), PARAMETER                   :: pprho0 = 1020.     ! water density (kg/m3)
   REAL(KIND=4), PARAMETER                   :: ppcp   = 4000.     ! calorific capacity (J/kg/m3)
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e3t                 ! vertical metric
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zt                  ! working input variable
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask               ! npiglo x npjglo
   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: gdepw               ! depth
   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: gdepo               ! output depth 
   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                 ! time counter
   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: e31d                ! vertical metrics in case of full step
   REAL(KIND=4)                              :: rdep1, rdep2        ! depth counters
   REAL(KIND=4)                              :: tol  = 1.0             ! tolerance 
   REAL(KIND=4)                              :: sclf = 1.0            ! scale factor

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_vint1, dl_vint2  ! verticall int quantity         

   CHARACTER(LEN=256)                        :: cf_in, cf_out       ! input/output file
   CHARACTER(LEN=256)                        :: cv_in, cv_out       ! variable name in and out
   CHARACTER(LEN=256)                        :: cunits, clongname   ! variable attributes
   CHARACTER(LEN=256)                        :: cldum               ! dummy string for command line browsing

   LOGICAL                                   :: lfull =.FALSE.      ! flag for full step computation
   LOGICAL                                   :: lgsop =.FALSE.      ! selected depths gsop intercomparison
   LOGICAL                                   :: locci =.FALSE.      ! selected 3 depths for occiput 
   LOGICAL                                   :: lchk  =.FALSE.      ! flag for missing files
   LOGICAL                                   :: lout                ! check for output

   TYPE(variable), DIMENSION(1)              :: stypvar             ! extension for attributes
   !!----------------------------------------------------------------------
   CALL ReadCdfNames()

   narg= iargc()
   IF ( narg == 0 ) THEN
      PRINT *,' usage : cdfvint T-file [IN-var] [-GSOP] [-OCCI]  [-full] '
      PRINT *,'      '
      PRINT *,'     PURPOSE :'
      PRINT *,'          Compute the vertical integral of the variable from top '
      PRINT *,'       to bottom, and save the cumulated valued, level by level.'
      PRINT *,'       For temperature (default var), the integral is transformed'
      PRINT *,'       to Heat Content ( J.K. m^-2) hence for salinity, the integral'
      PRINT *,'       represents PSU.m '
      PRINT *,'      '
      PRINT *,'     ARGUMENTS :'
      PRINT *,'         T-file : gridT file holding either Temperature or salinity '
      PRINT *,'        [IN-var ] : name of input variable to process. Default is '
      PRINT *,'               ', TRIM(cn_votemper),'. Can also be ',TRIM(cn_vosaline)
      PRINT *,'      '
      PRINT *,'     OPTIONS :'
      PRINT *,'        -GSOP : Use 7 GSOP standard level for the output '
      PRINT *,'                Default is to take the model levels for the output'
      PRINT *,'        -OCCI : Use 3 levels for the output: 700m, 2000m and bottom'
      PRINT *,'                Default is to take the model levels for the output'
      PRINT *,'        -full : for full step computation ' 
      PRINT *,'      '
      PRINT *,'     REQUIRED FILES :'
      PRINT *,'       ', TRIM(cn_fmsk),', ',TRIM(cn_fhgr),' and ', TRIM(cn_fzgr) 
      PRINT *,'      '
      PRINT *,'     OUTPUT : '
      PRINT *,'       netcdf file :  VAR-name.nc'
      PRINT *,'         variables :  either voheatc or vohsalt'
      PRINT *,'      '
      PRINT *,'     SEE ALSO :'
      PRINT *,'        cdfvertmean, cdfheatc, cdfmxlhcsc and  cdfmxlheatc'
      PRINT *,'      '
      STOP
   ENDIF

   ! default values
   cv_in = cn_votemper

   ! browse command line
   ijarg = 1   ; ij = 0
   DO WHILE ( ijarg <= narg ) 
      CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
      SELECT CASE ( cldum)
      CASE ( '-GSOP' ) ; lgsop = .TRUE.
      CASE ( '-OCCI' ) ; locci = .TRUE.
      CASE ( '-full' ) ; lfull = .TRUE. 
      CASE DEFAULT     
         ij = ij + 1
         SELECT CASE ( ij)
         CASE ( 1 ) ; cf_in = cldum
         CASE ( 2 ) ; cv_in = cldum
         CASE DEFAULT ; PRINT *, ' ERROR: Too many arguments ! ' ; STOP
         END SELECT
      END SELECT
   END DO

   ! Security check
   lchk = chkfile ( cf_in   )
   lchk = chkfile ( cn_fmsk ) .OR. lchk
   lchk = chkfile ( cn_fhgr ) .OR. lchk
   lchk = chkfile ( cn_fzgr ) .OR. lchk
   IF ( lchk ) STOP ! missing files

   ! Set output information according to variable name
   IF ( cv_in == cn_votemper ) THEN
      cv_out    = 'voheatc'
      clongname = 'Heat Content per unit area'
      cunits    = '10^6 J/m2'
      sclf      = pprho0*ppcp/1.e6
   ELSEIF ( cv_in == cn_vosaline ) THEN
      cv_out    = 'vohsalt'
      clongname = 'Vertically Integrated Salinity'
      cunits    = 'psu.m'
      sclf      =1.
   ELSE
      PRINT *,'  ERROR: Variable ', TRIM(cv_in), ' not pre-defined ...'
      PRINT *,'     Accepted variables are ', TRIM(cn_votemper),' and ',TRIM(cn_vosaline)
      STOP
   ENDIF

   ! log information so far
   cf_out = TRIM(cv_out)//'.nc'
   PRINT *,' INPUT VARIABLE  : ' , TRIM(cv_in)
   PRINT *,' OUTPUT VARIABLE : ' , TRIM(cv_out)
   PRINT *,' OUTPUT FILE     : ' , TRIM(cf_out)

   npiglo = getdim (cf_in, cn_x )
   npjglo = getdim (cf_in, cn_y )
   npk    = getdim (cf_in, cn_z )
   npt    = getdim (cf_in, cn_t )

   IF ( lgsop ) THEN ; PRINT *,' using GSOP depths' ; npko = 7 ! SL: commented out
   ENDIF

   IF ( locci ) THEN ; PRINT *,' using OCCIPUT depths' ; npko = 3
   ENDIF
   
   IF ((.not.lgsop).and.(.not.locci)) THEN ; PRINT *,' using model depths'; npko = npk
   ENDIF

   PRINT *, ' NPIGLO = ', npiglo
   PRINT *, ' NPJGLO = ', npjglo
   PRINT *, ' NPK    = ', npk
   PRINT *, ' NPKO   = ', npko
   PRINT *, ' NPT    = ', npt

   ! Allocate arrays
   ALLOCATE ( tim(npt) )
   ALLOCATE ( tmask(npiglo,npjglo) )
   ALLOCATE ( zt(npiglo,npjglo)  )
   ALLOCATE ( e3t(npiglo,npjglo) )
   ALLOCATE ( e31d(npk)   )
   ALLOCATE ( gdepw(npk), gdepo(npko) ) 

   ALLOCATE ( dl_vint1(npiglo, npjglo), dl_vint2(npiglo,npjglo) )

   ! prepare output variable
   ipk(:)                       = npko
   stypvar(1)%cname             = TRIM(cv_out)
   stypvar(1)%cunits            = TRIM(cunits)
   stypvar(1)%rmissing_value    = 0.
   stypvar(1)%valid_min         = -1.e15
   stypvar(1)%valid_max         =  1.e15
   stypvar(1)%clong_name        = TRIM(clongname)
   stypvar(1)%cshort_name       = TRIM(cv_out)
   stypvar(1)%conline_operation = 'N/A'
   stypvar(1)%caxis             = 'TZYX'

   ! Initialize output file
   gdepw(:) = getvare3(cn_fzgr, cn_gdepw, npk )
   e31d(:)  = getvare3(cn_fzgr, cn_ve3t,  npk )

   IF ( lgsop ) THEN
   gdepo(:) = (/100.,300.,500.,700.,800.,2000.,6000./) ! GEOP standard levs
   ENDIF

   IF ( locci ) THEN 
   gdepo(:) = (/700.,2000.,6000./) ! SL: occiput levels
   ENDIF

   IF ((.not.lgsop).and.(.not.locci)) THEN 
      gdepo(1:npk-1) = gdepw(2:npk)
      gdepo(npk)     = 6000.
   ENDIF



   PRINT*,'OUTPUT DEPTHS ARE : ',gdepo
   ncout = create      (cf_out, 'none', npiglo, npjglo, npko, cdep=cn_vdepthw, ld_xycoo=.FALSE.)
   ierr  = createvar   (ncout, stypvar, 1, ipk, id_varout )
   ierr  = putheadervar(ncout, cf_in,   npiglo, npjglo, npko, pdep=gdepo,     ld_xycoo=.FALSE.)

   tim   = getvar1d    (cf_in, cn_vtimec, npt     )
   ierr  = putvar1d    (ncout, tim,       npt, 'T')

   PRINT *, 'Output files initialised ...'

   DO jt = 1, npt
      dl_vint1(:,:) = 0.d0
      iko  = 1
      rdep1 = 0.0 ; rdep2 = 0.0
      lout = .TRUE.
      DO jk = 1, npk
         IF ( lgsop .or. locci ) lout = .FALSE.
         rdep1          = rdep2
         dl_vint2(:,:) = dl_vint1 (:,:)

         tmask(:,:)= getvar(cn_fmsk, 'tmask', jk, npiglo, npjglo           )
         zt(:,:)   = getvar(cf_in,   cv_in,   jk, npiglo, npjglo, ktime=jt )
         IF ( lfull ) THEN ; e3t(:,:) = e31d(jk)
                      ELSE ; e3t(:,:) = getvar(cn_fzgr, 'e3t_ps', jk, npiglo, npjglo, ldiom=.TRUE.)
         ENDIF

         rdep2 = rdep1 + e31d(jk)
         dl_vint1(:,:) = dl_vint1(:,:)+ zt(:,:)*e3t(:,:)*tmask(:,:)*1.d0

         IF ( rdep2 >= (gdepo(iko) - tol ) ) THEN
            lout=.TRUE.
            !modify vertical thickness for output
            WHERE (e3t(:,:) >  gdepo(iko) - rdep1 ) e3t(:,:) = gdepo(iko)-rdep1
            dl_vint2(:,:) = dl_vint2(:,:)+ zt(:,:)*e3t(:,:)*tmask(:,:)*1.d0
         ENDIF

         IF ( lout ) THEN
            dl_vint2(:,:) = dl_vint2(:,:) * sclf
            IF (jt == 1 ) THEN
               PRINT *,'Output for level ',iko
               PRINT *,'rdep1, rdep2, depo ',rdep1,rdep2,gdepo(iko) 
            ENDIF
            ierr = putvar(ncout, id_varout(1) ,REAL(dl_vint2), iko, npiglo, npjglo, ktime=jt)
            iko = iko + 1
         ENDIF
      END DO  ! loop to next level
   END DO  ! next time frame

   ierr = closeout(ncout)

END PROGRAM cdfvint
