PROGRAM cdfbotpressure
   !!======================================================================
   !!                     ***  PROGRAM  cdfbotpressure  ***
   !!=====================================================================
   !!  ** Purpose : Compute bottom pressure from insitu density
   !!
   !!  ** Method  : Vertical integral of rho g dz
   !!               eventually takes into account the SSH, and full step
   !!
   !! History : 3.0  : 01/2013  : J.M. Molines   : Original code from cdfvint
   !!----------------------------------------------------------------------
   USE cdfio
   USE eos
   USE modcdfnames
   USE modutils
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
   INTEGER(KIND=4),             DIMENSION(1) :: ipk, id_varout      ! only one output variable
   INTEGER(KIND=4)                           :: ncout 

   REAL(KIND=4), PARAMETER                   :: pp_grav = 9.81      ! Gravity
   REAL(KIND=4), PARAMETER                   :: pp_rau0 = 1035.e0   ! Reference density (as in NEMO)
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e3t                 ! vertical metric
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: hdept               ! vertical metric
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zt                  ! Temperature
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zs                  ! Salinity
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask               ! npiglo x npjglo
   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: gdepw               ! depth
   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                 ! time counter
   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: e31d                ! vertical metrics in case of full step

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_psurf            ! Surface pressure
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_bpres            ! Bottom pressure
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_sigi             ! insitu density

   CHARACTER(LEN=256)                        :: cf_in, cf_out       ! input/output file
   CHARACTER(LEN=256)                        :: cldum               ! dummy string for command line browsing
   CHARACTER(LEN=256)                        :: cglobal            ! Global attribute

   LOGICAL                                   :: lfull =.FALSE.      ! flag for full step computation
   LOGICAL                                   :: lssh  =.FALSE.      ! Use ssh and cst surf. density in the bot pressure
   LOGICAL                                   :: lssh2 =.FALSE.      ! Use ssh and variable surf.density in the bot pressure
   LOGICAL                                   :: lchk  =.FALSE.      ! flag for missing files

   TYPE(variable), DIMENSION(1)              :: stypvar             ! extension for attributes
   !!----------------------------------------------------------------------
   CALL ReadCdfNames()

   narg= iargc()
   IF ( narg == 0 ) THEN
      PRINT *,' usage : cdfbotpressure T-file [-full] [-ssh] [-ssh2 ] '
      PRINT *,'      '
      PRINT *,'     PURPOSE :'
      PRINT *,'          Compute the vertical bottom pressure (pa) from in situ density'
      PRINT *,'      '
      PRINT *,'     ARGUMENTS :'
      PRINT *,'         T-file : gridT file holding either Temperature and  salinity '
      PRINT *,'      '
      PRINT *,'     OPTIONS :'
      PRINT *,'        -full : for full step computation ' 
      PRINT *,'        -ssh  : Also take SSH into account in the computation'
      PRINT *,'                In this case, use rau0=',pp_rau0,' kg/m3 for '
      PRINT *,'                surface density (as in NEMO)'
      PRINT *,'                If you want to use 2d surface density from '
      PRINT *,'                the model, use option -ssh2'
      PRINT *,'        -ssh2 : as option -ssh but surface density is taken from '
      PRINT *,'                the model instead of a constant'
      PRINT *,'      '
      PRINT *,'     REQUIRED FILES :'
      PRINT *,'       ', TRIM(cn_fmsk),' and ', TRIM(cn_fzgr) 
      PRINT *,'      '
      PRINT *,'     OUTPUT : '
      PRINT *,'       netcdf file :  botpressure.nc'
      PRINT *,'         variables :  sobotpres'
      PRINT *,'      '
      PRINT *,'     SEE ALSO :'
      PRINT *,'        cdfvint'
      PRINT *,'      '
      STOP
   ENDIF

   ! browse command line
   ijarg = 1   ; ij = 0
   DO WHILE ( ijarg <= narg ) 
      CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
      SELECT CASE ( cldum)
      CASE ( '-ssh'  ) ; lssh  = .TRUE. 
      CASE ( '-ssh2' ) ; lssh2 = .TRUE. 
      CASE ( '-full' ) ; lfull = .TRUE. 
      CASE DEFAULT     
         ij = ij + 1
         SELECT CASE ( ij)
         CASE ( 1 ) ; cf_in = cldum
         CASE DEFAULT ; PRINT *, ' ERROR: Too many arguments ! ' ; STOP
         END SELECT
      END SELECT
   END DO
   CALL SetGlobalAtt(cglobal)

   ! Security check
   lchk = chkfile ( cf_in   )
   lchk = chkfile ( cn_fmsk ) .OR. lchk
   lchk = chkfile ( cn_fzgr ) .OR. lchk
   IF ( lchk ) STOP ! missing files

   ! log information so far
   cf_out = 'botpressure.nc'

   npiglo = getdim (cf_in, cn_x )
   npjglo = getdim (cf_in, cn_y )
   npk    = getdim (cf_in, cn_z )
   npt    = getdim (cf_in, cn_t )

   PRINT *, ' NPIGLO = ', npiglo
   PRINT *, ' NPJGLO = ', npjglo
   PRINT *, ' NPK    = ', npk
   PRINT *, ' NPT    = ', npt

   ! Allocate arrays
   ALLOCATE ( tim(npt) )
   ALLOCATE ( tmask(npiglo,npjglo) )
   ALLOCATE ( zt(npiglo,npjglo)  )
   ALLOCATE ( zs(npiglo,npjglo)  )
   ALLOCATE ( e3t(npiglo,npjglo) )
   ALLOCATE ( hdept(npiglo,npjglo) )
   ALLOCATE ( e31d(npk)   )
   ALLOCATE ( gdepw(npk)) 

   ALLOCATE ( dl_bpres(npiglo, npjglo))
   ALLOCATE ( dl_psurf(npiglo, npjglo))
   ALLOCATE ( dl_sigi(npiglo, npjglo))

   ! prepare output variable
   ipk(:)                       = 1
   stypvar(1)%cname             = 'sobotpres'
   stypvar(1)%cunits            = 'Pascal'
   stypvar(1)%rmissing_value    = 0.
   stypvar(1)%valid_min         = 0
   stypvar(1)%valid_max         =  1.e15
   stypvar(1)%clong_name        = 'Bottom Pressure'
   stypvar(1)%cshort_name       = 'sobotpres'
   stypvar(1)%conline_operation = 'N/A'
   stypvar(1)%caxis             = 'TYX'

   ! Initialize output file
   gdepw(:) = getvare3(cn_fzgr, cn_gdepw, npk )
   e31d(:)  = getvare3(cn_fzgr, cn_ve3t,  npk )

   ncout = create      (cf_out, cf_in, npiglo, npjglo, 1                    )
   ierr  = createvar   (ncout, stypvar, 1, ipk, id_varout, cdglobal=cglobal )
   ierr  = putheadervar(ncout, cf_in,   npiglo, npjglo, 1                   )

   tim   = getvar1d    (cf_in, cn_vtimec, npt     )
   ierr  = putvar1d    (ncout, tim,       npt, 'T')

   PRINT *, 'Output files initialised ...'

   DO jt = 1, npt
      IF ( lssh ) THEN
        zt(:,:)       = getvar(cf_in, cn_sossheig, 1, npiglo, npjglo, ktime=jt )
        dl_psurf(:,:) = pp_grav * pp_rau0 * zt(:,:)
      ELSE IF ( lssh2 ) THEN 
         zt(:,:)    = getvar(cf_in,   cn_votemper, 1, npiglo, npjglo, ktime=jt )
         zs(:,:)    = getvar(cf_in,   cn_vosaline, 1, npiglo, npjglo, ktime=jt )
     
         dl_sigi(:,:) = 1000. + sigmai(zt, zs, 0., npiglo, npjglo)

         !  CAUTION : hdept is used for reading SSH in the next line
         hdept(:,:)   = getvar(cf_in, cn_sossheig, 1, npiglo, npjglo, ktime=jt )
        dl_psurf(:,:) = pp_grav * dl_sigi * hdept(:,:)
      ELSE
        dl_psurf(:,:)=0.d0
      ENDIF

      dl_bpres(:,:) = dl_psurf(:,:)

      DO jk = 1, npk
         tmask(:,:) = getvar(cn_fmsk, 'tmask',      jk, npiglo, npjglo           )
         hdept(:,:) = getvar(cn_fzgr, cn_hdept,     jk, npiglo, npjglo           )
         zt(:,:)    = getvar(cf_in,   cn_votemper,  jk, npiglo, npjglo, ktime=jt )
         zs(:,:)    = getvar(cf_in,   cn_vosaline,  jk, npiglo, npjglo, ktime=jt )
     
         dl_sigi(:,:) = 1000. + sigmai(zt, zs, hdept, npiglo, npjglo)
         IF ( lfull ) THEN ; e3t(:,:) = e31d(jk)
                      ELSE ; e3t(:,:) = getvar(cn_fzgr, 'e3t_ps', jk, npiglo, npjglo, ldiom=.TRUE.)
         ENDIF
         dl_bpres(:,:) = dl_bpres(:,:) + dl_sigi(:,:) * e3t(:,:) * pp_grav * 1.d0 * tmask(:,:) 

      END DO  ! loop to next level
         ierr = putvar(ncout, id_varout(1) ,REAL(dl_bpres), 1, npiglo, npjglo, ktime=jt)
   END DO  ! next time frame

   ierr = closeout(ncout)

END PROGRAM cdfbotpressure
