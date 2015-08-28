PROGRAM cdfhgradb
   !!======================================================================
   !!                     ***  PROGRAM  cdfhgradb ***
   !!=====================================================================
   !!  ** Purpose : Compute the norm of the lateral buoyancy gradient.
   !!
   !!  ** Method  :
   !!
   !! History : 3.0  : 08/2015  : J. Le Sommer (from cdfgradT)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   routines      : description
   !!   CreateOutput  :  perform all the stuff linked with output file 
   !!----------------------------------------------------------------------
   USE cdfio
   USE eos
   USE modcdfnames
   !!----------------------------------------------------------------------
   !! CDFTOOLS_3.0 , MEOM 2015
   !! $Id$
   !! Copyright (c) 2012, J.-M. Molines
   !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
   !!----------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER(KIND=4), PARAMETER                  :: jp_varout = 1     ! number of output variables
   INTEGER(KIND=4)                             :: jk, jt            ! dummy loop index
   INTEGER(KIND=4)                             :: ji, jj, jlev      ! dummy loop index
   INTEGER(KIND=4)                             :: narg, iargc       ! command line
   INTEGER(KIND=4)                             :: npiglo, npjglo    ! size of the domain
   INTEGER(KIND=4)                             :: npk, npt          ! size of the domain
   INTEGER(KIND=4)                             :: ncout             ! ncid of output variable
   INTEGER(KIND=4)                             :: ierr              ! error status
   INTEGER(KIND=4), DIMENSION(jp_varout)       :: ipk, id_varout    ! output variable

   REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: tim               ! time variable
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: umask, vmask, tmask ! relevant mask
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1u, e2v          ! metrics
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zt, zs            ! Temperature Salinity on 2 levels
   REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: gdept             ! depth of Tlevels
   REAL(KIND=4)                                :: zdep              ! temporary scalar

   REAL(KIND=8)                                :: zgrav=9.81

   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dgradt_xu, dgradt_yv! Temperature gradient, native grid
   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dgrads_xu, dgrads_yv! Salinity gradient, native grid

   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dgradt_xt, dgradt_yt! Temperature gradient, t-point
   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dgrads_xt, dgrads_yt! Salinity gradient, t-point

   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dgradb_xt, dgradb_yt ! Buyoyancy gradient, t-point
   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dgradb              ! norm of the buyoyancy gradient, t-point

   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zalbet, zbeta         ! for alpha and beta

   CHARACTER(LEN=256)                         :: cf_tfil             ! input file name for T and S
   CHARACTER(LEN=256)                         :: cf_sfil             ! input file name for S (optional)
   CHARACTER(LEN=256)                         :: cf_out = 'hgradb_gridT.nc' ! output file name

   TYPE(variable), DIMENSION(jp_varout)       :: stypvar             ! output data structure

   LOGICAL                                    :: lchk = .FALSE.      ! flag for missing files
   !!----------------------------------------------------------------------
   CALL ReadCdfNames()

   narg= iargc()
   IF ( narg == 0 ) THEN
      PRINT *,' usage : cdfhgradb T-file [S-file] '
      PRINT *,'      '
      PRINT *,'     PURPOSE :'
      PRINT *,'        Compute the norm of the horizontal buoyancy gradient.'
      PRINT *,'      Results are saved at T points.'
      PRINT *,'      '
      PRINT *,'     ARGUMENTS :'
      PRINT *,'       T-file : File with ',TRIM(cn_votemper),' and ',TRIM(cn_vosaline),' variables' 
      PRINT *,'           If ',TRIM(cn_vosaline),' not in T-file give a second name for S-file.'
      PRINT *,'      '
      PRINT *,'     OPTIONS :'
      PRINT *,'       S-file : File with ',TRIM(cn_vosaline),' variable if not in T file'
      PRINT *,'      '
      PRINT *,'     REQUIRED FILES :'
      PRINT *,'       ', TRIM(cn_fhgr),' ',TRIM(cn_fmsk),' and ',TRIM(cn_fzgr)
      PRINT *,'      '
      PRINT *,'     OUTPUT : '
      PRINT *,'       netcdf file : ', TRIM(cf_out) 
      PRINT *,'         ',jp_varout,' variables : '
      PRINT *,'              vohgradb: norm of the horizontal buoyancy gradient at t-point'
      PRINT *,'      '
      PRINT *,'     SEE ALSO :'
      PRINT *,'      '
      PRINT *,'      '
      STOP
   ENDIF

   CALL getarg (1, cf_tfil)
   IF ( narg == 2 ) THEN
      CALL getarg( 2, cf_sfil)
   ELSE
      cf_sfil = cf_tfil   ! if only one file provided, they are identical
   ENDIF

   lchk = ( lchk .OR. chkfile(cf_tfil) )
   lchk = ( lchk .OR. chkfile(cf_sfil) )
   lchk = ( lchk .OR. chkfile(cn_fhgr) )
   lchk = ( lchk .OR. chkfile(cn_fzgr) )
   lchk = ( lchk .OR. chkfile(cn_fmsk) )
   IF (lchk ) STOP ! missing file

   npiglo = getdim (cf_tfil, cn_x)
   npjglo = getdim (cf_tfil, cn_y)
   npk    = getdim (cf_tfil, cn_z)
   npt    = getdim (cf_tfil, cn_t)

   PRINT *, 'npiglo = ', npiglo
   PRINT *, 'npjglo = ', npjglo
   PRINT *, 'npk    = ', npk
   PRINT *, 'npt    = ', npt

   !!  Allocate arrays
   ALLOCATE (tim(npt) )
   ALLOCATE (gdept(npk) )
   ALLOCATE (e1u  (npiglo,npjglo), e2v  (npiglo,npjglo) )
   ALLOCATE (umask(npiglo,npjglo), vmask(npiglo,npjglo), tmask(npiglo,npjglo))
   ALLOCATE (zt   (npiglo,npjglo), zs (npiglo,npjglo))
   ALLOCATE ( zalbet(npiglo,npjglo), zbeta(npiglo, npjglo) )
   ALLOCATE (dgradt_xu(npiglo,npjglo), dgradt_yv(npiglo,npjglo))
   ALLOCATE (dgrads_xu(npiglo,npjglo), dgrads_yv(npiglo,npjglo))
   ALLOCATE (dgradt_xt(npiglo,npjglo), dgradt_yt(npiglo,npjglo))
   ALLOCATE (dgrads_xt(npiglo,npjglo), dgrads_yt(npiglo,npjglo))
   ALLOCATE (dgradb_xt(npiglo,npjglo), dgradb_yt(npiglo,npjglo))
   ALLOCATE (dgradb(npiglo,npjglo))

   CALL CreateOutput

   tim  = getvar1d(cf_tfil, cn_vtimec, npt     )
   ierr = putvar1d(ncout,  tim,        npt, 'T')

   e1u =  getvar(cn_fhgr, cn_ve1u, 1, npiglo, npjglo)
   e2v =  getvar(cn_fhgr, cn_ve2v, 1, npiglo, npjglo)
   gdept(:)    = getvare3(cn_fzgr, cn_gdept, npk )

   DO jt = 1,npt
      DO jk = npk, 1, -1  

         PRINT *,'level ',jk
         zdep = gdept(jk)

         ! read files
         zt(:,:) = getvar(cf_tfil, cn_votemper, jk,   npiglo, npjglo, ktime=jt)
         zs(:,:) = getvar(cf_sfil, cn_vosaline, jk,   npiglo, npjglo, ktime=jt)

         umask(:,:) = getvar(cn_fmsk, 'umask' , jk, npiglo, npjglo )
         vmask(:,:) = getvar(cn_fmsk, 'vmask' , jk, npiglo, npjglo )
         tmask(:,:) = getvar(cn_fmsk, 'tmask' , jk, npiglo, npjglo )

         ! zonal grad located at U point at current level
         dgradt_xu(:,:) = 0.d0
         dgradt_xu(1:npiglo-1,:) = 1.d0 / e1u(1:npiglo-1,:) * &
              &                 ( zt(2:npiglo,:) - zt(1:npiglo-1,:) ) * umask(1:npiglo-1,:)
         dgrads_xu(:,:) = 0.d0
         dgrads_xu(1:npiglo-1,:) = 1.d0 / e1u(1:npiglo-1,:) * &
              &                 ( zs(2:npiglo,:) - zs(1:npiglo-1,:) ) * umask(1:npiglo-1,:)

         ! meridional grad located at V point at current level
         dgradt_yv(:,:) = 0.d0
         dgradt_yv(:,1:npjglo-1) = 1.d0 / e2v(:,1:npjglo-1) * &
              &                 ( zt(:,2:npjglo) - zt(:,1:npjglo-1) ) * vmask(:,1:npjglo-1)
         dgrads_yv(:,:) = 0.d0
         dgrads_yv(:,1:npjglo-1) = 1.d0 / e2v(:,1:npjglo-1) * &
              &                 ( zs(:,2:npjglo) - zs(:,1:npjglo-1) ) * vmask(:,1:npjglo-1)

         ! get temperature and salinity gradients at t-point 
         dgradt_xt(:,:) = 0.d0
         dgradt_yt(:,:) = 0.d0
         dgrads_xt(:,:) = 0.d0
         dgrads_yt(:,:) = 0.d0

         DO ji=npiglo,2,-1
             DO jj=1,npjglo
               dgradt_xt(ji,jj) = 0.5*(dgradt_xu(ji-1,jj)+dgradt_xu(ji,jj))
               dgrads_xt(ji,jj) = 0.5*(dgrads_xu(ji-1,jj)+dgrads_xu(ji,jj))
               dgradt_yt(ji,jj) = 0.5*(dgradt_yv(ji,jj-1)+dgradt_yv(ji,jj))
               dgrads_yt(ji,jj) = 0.5*(dgrads_yv(ji,jj-1)+dgrads_yv(ji,jj))
             ENDDO
         ENDDO
         
         !DO ji=1,npiglo
         !    DO jj=npjglo,2 -1
         !      dgradt_yt(ji,jj) = 0.5*(dgradt_yv(ji,jj-1)+dgradt_yv(ji,jj))
         !      dgrads_yt(ji,jj) = 0.5*(dgrads_yv(ji,jj-1)+dgrads_yv(ji,jj))
         !    ENDDO
         ! ENDDO

         ! compute alpha and beta at t-point
         zalbet(:,:) = albet ( zt, zs, zdep, npiglo, npjglo) ! not exact for partial-step level
         zbeta (:,:) = beta  ( zt, zs, zdep, npiglo, npjglo) ! 

         ! compute the buoyancy gradients at t-point
         dgradb_xt(:,:) = zbeta(:,:) * ( zalbet(:,:) * dgradt_xt(:,:) - dgrads_xt(:,:) ) * tmask(:,:)
         dgradb_yt(:,:) = zbeta(:,:) * ( zalbet(:,:) * dgradt_yt(:,:) - dgrads_yt(:,:) ) * tmask(:,:)

         ! get the norm of the buoyancy gradients at t-point
         dgradb(:,:) = zgrav * SQRT( dgradb_xt(:,:) * dgradb_xt(:,:) + dgradb_yt(:,:) * dgradb_yt(:,:) )

         ! write
         ierr = putvar(ncout, id_varout(1), REAL(dgradb), jk, npiglo, npjglo, ktime=jt)
      END DO
   END DO

   ierr = closeout(ncout)

CONTAINS
   SUBROUTINE CreateOutput
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE CreateOutput  ***
      !!
      !! ** Purpose :  Set up all things required for the output file, create
      !!               the file and write the header part. 
      !!
      !! ** Method  :  Use global module variables 
      !!
      !!----------------------------------------------------------------------
      ipk(:) = npk  !  3D


      stypvar(1)%cname             = 'vohgradb'
      stypvar(1)%cunits            = 's^{-2}'
      stypvar(1)%rmissing_value    = -1000.
      stypvar(1)%valid_min         = -1.
      stypvar(1)%valid_max         = 1.
      stypvar(1)%clong_name        = 'Horizontal gradient of buoyancy.'
      stypvar(1)%cshort_name       = 'vohgradb'
      stypvar(1)%conline_operation = 'N/A'

      ! create output fileset
      ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk       )
      ierr  = createvar   (ncout,  stypvar, jp_varout, ipk, id_varout )
      ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk       )

   END SUBROUTINE CreateOutput

END PROGRAM cdfhgradb
