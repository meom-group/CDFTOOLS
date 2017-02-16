PROGRAM cdfgradT
   !!======================================================================
   !!                     ***  PROGRAM  cdfgradT ***
   !!=====================================================================
   !!  ** Purpose :
   !!
   !!  ** Method  :
   !!
   !! History : 3.0  : 05/2013  : N. Ducousso
   !!           3.0  : 06/2013  : J.M. Molines : Doctor norm + Lic.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   routines      : description
   !!   CreateOutput  :  perform all the stuff linked with output file 
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

   INTEGER(KIND=4), PARAMETER                  :: jp_varout = 6     ! number of output variables
   INTEGER(KIND=4)                             :: jk, jt            ! dummy loop index
   INTEGER(KIND=4)                             :: narg, iargc       ! command line
   INTEGER(KIND=4)                             :: npiglo, npjglo    ! size of the domain
   INTEGER(KIND=4)                             :: npk, npt          ! size of the domain
   INTEGER(KIND=4)                             :: ncout             ! ncid of output variable
   INTEGER(KIND=4)                             :: ierr              ! error status
   INTEGER(KIND=4)                             :: iup = 1, icur = 2 ! 2 working level
   INTEGER(KIND=4)                             :: itmp              ! working integer for swapping arrays
   INTEGER(KIND=4), DIMENSION(jp_varout)       :: ipk, id_varout    ! output variable

   REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: tim               ! time variable
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: umask, vmask, wmask ! relevant mask
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1u, e2v, e3w     ! metrics
   REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: zt, zs            ! Temperature Salinity on 2 levels

   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dgradt_x, dgradt_y, dgradt_z  ! Temperature gradient component
   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dgrads_x, dgrads_y, dgrads_z  ! Salinity gradient component

   CHARACTER(LEN=2048)                         :: cf_tfil             ! input file name for T and S
   CHARACTER(LEN=2048)                         :: cf_sfil             ! input file name for S (optional)
   CHARACTER(LEN=2048)                         :: cf_out = 'gradT.nc' ! output file name

   TYPE(variable), DIMENSION(jp_varout)       :: stypvar             ! output data structure

   LOGICAL                                    :: lchk = .FALSE.      ! flag for missing files
   !!----------------------------------------------------------------------
   CALL ReadCdfNames()

   narg= iargc()
   IF ( narg == 0 ) THEN
      PRINT *,' usage : cdfgradT T-file [S-file] '
      PRINT *,'      '
      PRINT *,'     PURPOSE :'
      PRINT *,'        Compute horizontal and vertical gradient of temperature and salinity.'
      PRINT *,'      Results are saved at U point for zonal gradient, V point for meridional'
      PRINT *,'      gradient and W for vertical gradient.'
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
      PRINT *,'              vozogradt, vomegradt, vovegradt : 3 component of the temperature'
      PRINT *,'                          located respectively at U, V and W points'
      PRINT *,'              vozograds, vomegrads, vovegrads : 3 component of the salinity'
      PRINT *,'                          located respectively at U, V and W points'
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
      cf_sfil = cf_tfil   ! ff only one file provided, they are identical
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
   ALLOCATE (e1u  (npiglo,npjglo), e2v  (npiglo,npjglo), e3w  (npiglo,npjglo))
   ALLOCATE (umask(npiglo,npjglo), vmask(npiglo,npjglo), wmask(npiglo,npjglo))
   ALLOCATE (zt   (npiglo,npjglo,2), zs (npiglo,npjglo,2))
   ALLOCATE (dgradt_x(npiglo,npjglo), dgradt_y(npiglo,npjglo), dgradt_z(npiglo,npjglo))
   ALLOCATE (dgrads_x(npiglo,npjglo), dgrads_y(npiglo,npjglo), dgrads_z(npiglo,npjglo))

   CALL CreateOutput

   tim  = getvar1d(cf_tfil, cn_vtimec, npt     )
   ierr = putvar1d(ncout,  tim,        npt, 'T')

   e1u =  getvar(cn_fhgr, cn_ve1u, 1, npiglo, npjglo)
   e2v =  getvar(cn_fhgr, cn_ve2v, 1, npiglo, npjglo)

   DO jt = 1,npt
      DO jk = npk, 1, -1  !! Main loop : (2 levels of T are required : iup, icur)

         PRINT *,'level ',jk
         ! read files
         IF ( jk == npk   ) THEN  ! first time need to read 2 levels
            zt(:,:,iup ) = getvar(cf_tfil, cn_votemper, jk-1, npiglo, npjglo, ktime=jt)
            zt(:,:,icur) = getvar(cf_tfil, cn_votemper, jk,   npiglo, npjglo, ktime=jt)
            zs(:,:,iup ) = getvar(cf_sfil, cn_vosaline, jk-1, npiglo, npjglo, ktime=jt)
            zs(:,:,icur) = getvar(cf_sfil, cn_vosaline, jk,   npiglo, npjglo, ktime=jt)
         ELSEIF ( jk == 1 ) THEN ! surface : vertical gradient is N/A
            iup = icur        
         ELSE                  ! only read iup, icur has been swapped
            zt(:,:,iup ) = getvar(cf_tfil, cn_votemper, jk-1, npiglo, npjglo, ktime=jt)
            zs(:,:,iup ) = getvar(cf_sfil, cn_vosaline, jk-1, npiglo, npjglo, ktime=jt)
         ENDIF

!        IF (jk == 1) THEN
!           zt(:,:,iup ) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
!           zt(:,:,icur) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
!           zs(:,:,iup ) = getvar(cf_sfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)
!           zs(:,:,icur) = getvar(cf_sfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)
!        ELSE        
!           zt(:,:,iup ) = getvar(cf_tfil, cn_votemper, jk-1, npiglo, npjglo, ktime=jt)
!           zt(:,:,icur) = getvar(cf_tfil, cn_votemper, jk,   npiglo, npjglo, ktime=jt)
!           zs(:,:,iup ) = getvar(cf_sfil, cn_vosaline, jk-1, npiglo, npjglo, ktime=jt)
!           zs(:,:,icur) = getvar(cf_sfil, cn_vosaline, jk,   npiglo, npjglo, ktime=jt)
!        END IF

         e3w(:,:)   = getvar(cn_fzgr, 'e3w_ps', jk, npiglo, npjglo, ldiom=.TRUE.)

         umask(:,:) = getvar(cn_fmsk, 'umask' , jk, npiglo, npjglo )
         vmask(:,:) = getvar(cn_fmsk, 'vmask' , jk, npiglo, npjglo )
         wmask(:,:) = getvar(cn_fmsk, 'tmask' , jk, npiglo, npjglo )

         ! zonal grad located at U point at current level
         dgradt_x(:,:) = 0.d0
         dgradt_x(1:npiglo-1,:) = 1.d0 / e1u(1:npiglo-1,:) * &
              &                 ( zt(2:npiglo,:,icur) - zt(1:npiglo-1,:,icur) ) * umask(1:npiglo-1,:)
         dgrads_x(:,:) = 0.d0
         dgrads_x(1:npiglo-1,:) = 1.d0 / e1u(1:npiglo-1,:) * &
              &                 ( zs(2:npiglo,:,icur) - zs(1:npiglo-1,:,icur) ) * umask(1:npiglo-1,:)

         ! meridional grad located at V point at current level
         dgradt_y(:,:) = 0.d0
         dgradt_y(:,1:npjglo-1) = 1.d0 / e2v(:,1:npjglo-1) * &
              &                 ( zt(:,2:npjglo,icur) - zt(:,1:npjglo-1,icur) ) * vmask(:,1:npjglo-1)
         dgrads_y(:,:) = 0.d0
         dgrads_y(:,1:npjglo-1) = 1.d0 / e2v(:,1:npjglo-1) * &
              &                 ( zs(:,2:npjglo,icur) - zs(:,1:npjglo-1,icur) ) * vmask(:,1:npjglo-1)

         ! vertical grad located at W point 
         dgradt_z(:,:) = 0.d0
         dgradt_z(:,:) = 1.d0 / e3w(:,:) * ( zt(:,:,iup) - zt(:,:,icur) ) * wmask(:,:)
         dgrads_z(:,:) = 0.d0
         dgrads_z(:,:) = 1.d0 / e3w(:,:) * ( zs(:,:,iup) - zs(:,:,icur) ) * wmask(:,:)

         ! write
         ierr = putvar(ncout, id_varout(1), REAL(dgradt_x), jk, npiglo, npjglo, ktime=jt)
         ierr = putvar(ncout, id_varout(2), REAL(dgradt_y), jk, npiglo, npjglo, ktime=jt)
         ierr = putvar(ncout, id_varout(3), REAL(dgradt_z), jk, npiglo, npjglo, ktime=jt)
         ierr = putvar(ncout, id_varout(4), REAL(dgrads_x), jk, npiglo, npjglo, ktime=jt)
         ierr = putvar(ncout, id_varout(5), REAL(dgrads_y), jk, npiglo, npjglo, ktime=jt)
         ierr = putvar(ncout, id_varout(6), REAL(dgrads_z), jk, npiglo, npjglo, ktime=jt)
        
         ! swap levels and read only up
         itmp = iup ; iup  = icur ; icur = itmp

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

      stypvar(1)%cname             = 'vozogradt'
      stypvar(1)%cunits            = ''
      stypvar(1)%rmissing_value    = -1000.
      stypvar(1)%valid_min         = -1.
      stypvar(1)%valid_max         = 1.
      stypvar(1)%clong_name        = 'zonal temper gradient'
      stypvar(1)%cshort_name       = 'vozogradt'
      stypvar(1)%conline_operation = 'N/A'
      stypvar(1)%caxis             = 'TZYX'

      stypvar(2)%cname             = 'vomegradt'
      stypvar(2)%cunits            = ''
      stypvar(2)%rmissing_value    = -1000.
      stypvar(2)%valid_min         = -1.
      stypvar(2)%valid_max         = 1.
      stypvar(2)%clong_name        = 'meridional temper gradient'
      stypvar(2)%cshort_name       = 'vomegradt'
      stypvar(2)%conline_operation = 'N/A'
      stypvar(2)%caxis             = 'TZYX'

      stypvar(3)%cname             = 'vovegradt'
      stypvar(3)%cunits            = ''
      stypvar(3)%rmissing_value    = -1000.
      stypvar(3)%valid_min         = -1.
      stypvar(3)%valid_max         = 1.
      stypvar(3)%clong_name        = 'vertical temper gradient'
      stypvar(3)%cshort_name       = 'vovegradt'
      stypvar(3)%conline_operation = 'N/A'
      stypvar(3)%caxis             = 'TZYX'

      stypvar(4)%cname             = 'vozograds'
      stypvar(4)%cunits            = ''
      stypvar(4)%rmissing_value    = -1000.
      stypvar(4)%valid_min         = -1.
      stypvar(4)%valid_max         = 1.
      stypvar(4)%clong_name        = 'zonal saline gradient'
      stypvar(4)%cshort_name       = 'vozograds'
      stypvar(4)%conline_operation = 'N/A'
      stypvar(4)%caxis             = 'TZYX'

      stypvar(5)%cname             = 'vomegrads'
      stypvar(5)%cunits            = ''
      stypvar(5)%rmissing_value    = -1000.
      stypvar(5)%valid_min         = -1.
      stypvar(5)%valid_max         = 1.
      stypvar(5)%clong_name        = 'meridional saline gradient'
      stypvar(5)%cshort_name       = 'vomegrads'
      stypvar(5)%conline_operation = 'N/A'
      stypvar(5)%caxis             = 'TZYX'

      stypvar(6)%cname             = 'vovegrads'
      stypvar(6)%cunits            = ''
      stypvar(6)%rmissing_value    = -1000.
      stypvar(6)%valid_min         = -1.
      stypvar(6)%valid_max         = 1.
      stypvar(6)%clong_name        = 'vertical saline gradient'
      stypvar(6)%cshort_name       = 'vovegrads'
      stypvar(6)%conline_operation = 'N/A'

      ! create output fileset
      ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk       )
      ierr  = createvar   (ncout,  stypvar, jp_varout, ipk, id_varout )
      ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk       )

   END SUBROUTINE CreateOutput

END PROGRAM cdfgradT
