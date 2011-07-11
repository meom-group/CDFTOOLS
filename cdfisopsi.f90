PROGRAM cdfisopsi
  !!======================================================================
  !!                     ***  PROGRAM  cdfisopsi  ***
  !!=====================================================================
  !!  ** Purpose : Compute a geostrophic streamfunction projected
  !!               on an isopycn (Ref: McDougall and ?, need reference)
  !!
  !!  ** Method  : read temp and salinity, compute sigmainsitu and sigma
  !!              at a reference level, projection of p,T,S on a given
  !!              isopycnal, compute specific volume anomaly and
  !!              integrates it.
  !!
  !! History : 2.1  : 12/2010  : R. Dussin 
  !!           3.0  : 03/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER(KIND=4), PARAMETER  :: jp_vars=7

  INTEGER(KIND=4)                             :: jj, ji, jk, jt            ! dummy loop index
  INTEGER(KIND=4)                             :: ierr                      ! working integer
  INTEGER(KIND=4)                             :: narg, iargc               ! command line arguments
  INTEGER(KIND=4)                             :: npiglo, npjglo  ! size of the domain
  INTEGER(KIND=4)                             :: npk, npt  ! size of the domain
  INTEGER(KIND=4)                             :: ik, ik0                        ! 
  INTEGER(KIND=4)                             :: ncout
  INTEGER(KIND=4), DIMENSION(jp_vars)         :: ipk              ! outptut variables : number of levels,
  INTEGER(KIND=4), DIMENSION(jp_vars)         :: id_varout        ! ncdf varid's

  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: v3d, ztemp3         ! 3d array
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: zsal3, zsva3        ! 3d array
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: ztemp, zsal , zssh  ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: ztemp0, zsal0       ! Arrays for reference profile
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zsiginsitu          ! in-situ density
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zsig0, zsigsurf     ! potential density of ref profile and surface
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zmask, zdep           ! 2D mask at current level, level depths
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: v2d             !  2d working arrays
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: ztempint ! 2d working arrays
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zsalint ! 2d working arrays
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zint ! 2d working arrays
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: pint    ! 2d working arrays
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: alpha ! 2d working arrays
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1t, e2t
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rdeltapsi1
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rdeltapsi2
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: psi0
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: psi
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zsva2
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: prof, tim              ! prof (m) and time (sec)
  REAL(KIND=4)                                :: P1, P2
  REAL(KIND=4)                                :: zspval            ! missing value
  REAL(KIND=4)                                :: refdepth
  REAL(KIND=4)                                :: zsigmaref      !
  REAL(KIND=4)                                :: ztmean, zsmean  ! mean temperature and salinity on isopycnal
  REAL(KIND=4)                                :: hmean, pmean  ! mean isopycnal depth and mean pressure

  CHARACTER(LEN=256) :: cf_tfil                              ! input gridT file
  CHARACTER(LEN=256) :: cf_out='isopsi.nc'                   ! output file name
  CHARACTER(LEN=256) :: cldum                                ! dummy character variable for reading

  TYPE(variable) , DIMENSION(jp_vars) :: stypvar         ! structure for attributes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfisopsi ref_level sigma_ref gridT '
     PRINT *,'      Compute  a geostrophic streamfunction'
     PRINT *,'      projected  on an isopycn.'
     PRINT *,'         ref_level = reference level for pot. density'
     PRINT *,'         sigma_ref = density level to project on'
     PRINT *,'         gridT     = input file for temperature and salinity'
     PRINT *,'  '
     PRINT *,'         Output on ',TRIM(cf_out),' variable soisopsi'
     PRINT *,'         Depths are taken from input file '
     PRINT *,'         requires ',TRIM(cn_fhgr),' and ',TRIM(cn_fzgr) 
     STOP
  ENDIF

  CALL getarg (1, cldum) ; READ (cldum,*) refdepth
  CALL getarg (2, cldum) ; READ (cldum,*) zsigmaref
  CALL getarg (3, cf_tfil)

  IF ( chkfile(cf_tfil) .OR. chkfile(cn_fzgr) .OR. chkfile(cn_fhgr) ) STOP  ! missing file

  PRINT *, 'Potential density referenced at ', refdepth , ' meters'
  PRINT *, 'Isopycn for projection is ', zsigmaref

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE ( prof(npk)         , tim(npt)           )
  ALLOCATE ( e1t(npiglo,npjglo), e2t(npiglo,npjglo) )

  e1t(:,:) = getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo)
  e2t(:,:) = getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo)

  !--------------------------------------------------------------------
  !!  Output file
  ipk(:)= 1  ! all variables are 2d
  stypvar(1)%cname             = 'votemper_interp'
  stypvar(1)%cunits            = 'DegC'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = -2.
  stypvar(1)%valid_max         = 45.
  stypvar(1)%clong_name        = 'Temperature interpolated on isopycnal layer'
  stypvar(1)%cshort_name       = 'votemper_interp'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TZYX'

  stypvar(2)%cname             = 'vosaline_interp'
  stypvar(2)%cunits            = 'PSU'
  stypvar(2)%rmissing_value    = 0.
  stypvar(2)%valid_min         = 0.
  stypvar(2)%valid_max         = 50.
  stypvar(2)%clong_name        = 'Salinity interpolated on isopycnal layer'
  stypvar(2)%cshort_name       = 'vosaline_interp'
  stypvar(2)%conline_operation = 'N/A'
  stypvar(2)%caxis             = 'TZYX'

  stypvar(3)%cname             = 'depth_interp'
  stypvar(3)%cunits            = 'meters'
  stypvar(3)%rmissing_value    = 0.
  stypvar(3)%valid_min         = 0.0
  stypvar(3)%valid_max         = 8000.
  stypvar(3)%clong_name        = 'Depth of the isopycnal layer'
  stypvar(3)%cshort_name       = 'depth_interp'
  stypvar(3)%conline_operation = 'N/A'
  stypvar(3)%caxis             = 'TZYX'

  stypvar(4)%cname             = 'soisopsi'
  stypvar(4)%cunits            = 'm2s-2 (to be verified)'
  stypvar(4)%rmissing_value    = 0.
  stypvar(4)%valid_min         = -500.
  stypvar(4)%valid_max         =  500.
  stypvar(4)%clong_name        = 'Total streamfunction on the isopycnal layer'
  stypvar(4)%cshort_name       = 'soisopsi'
  stypvar(4)%conline_operation = 'N/A'
  stypvar(4)%caxis             = 'TZYX'

  stypvar(5)%cname             = 'soisopsi0'
  stypvar(5)%cunits            = 'm2s-2 (to be verified)'
  stypvar(5)%rmissing_value    = 0.
  stypvar(5)%valid_min         = -500.
  stypvar(5)%valid_max         =  500.
  stypvar(5)%clong_name        = 'Contribution of the SSH'
  stypvar(5)%cshort_name       = 'soisopsi'
  stypvar(5)%conline_operation = 'N/A'
  stypvar(5)%caxis             = 'TZYX'

  stypvar(6)%cname             = 'soisopsi1'
  stypvar(6)%cunits            = 'm2s-2 (to be verified)'
  stypvar(6)%rmissing_value    = 0.
  stypvar(6)%valid_min         = -500.
  stypvar(6)%valid_max         =  500.
  stypvar(6)%clong_name        = 'Contribution of specific volume anomaly vertical integration'
  stypvar(6)%cshort_name       = 'soisopsi'
  stypvar(6)%conline_operation = 'N/A'
  stypvar(6)%caxis             = 'TZYX'

  stypvar(7)%cname             = 'soisopsi2'
  stypvar(7)%cunits            = 'm2s-2 (to be verified)'
  stypvar(7)%rmissing_value    = 0.
  stypvar(7)%valid_min         = -500.
  stypvar(7)%valid_max         =  500.
  stypvar(7)%clong_name        = 'Contribution of pressure term on the isopycnal layer'
  stypvar(7)%cshort_name       = 'soisopsi'
  stypvar(7)%conline_operation = 'N/A'
  stypvar(7)%caxis             = 'TZYX'

  ! create output fileset
  ncout = create      (cf_out, cf_tfil, npiglo,  npjglo, npk       )
  ierr  = createvar   (ncout,  stypvar, jp_vars, ipk,    id_varout ) 
  ierr  = putheadervar(ncout,  cf_tfil, npiglo,  npjglo, npk       )

  prof(:) = getvar1d(cf_tfil, cn_vdeptht, npk     )
  tim     = getvar1d(cf_tfil, cn_vtimec,  npt     )
  ierr    = putvar1d(ncout,   tim,        npt, 'T')

  zspval   = getatt(cf_tfil, cn_vosaline, cn_missing_value )

  !---------------------------------------------------------------------------
  DO jt=1,npt
     PRINT *,'time ',jt, tim(jt)/86400.,' days'

     !------------------------------------------------------------------------------
     ! 1. First we compute the potential density and store it into a 3d array
     ALLOCATE (ztemp(npiglo,npjglo), zsal(npiglo,npjglo), zmask(npiglo,npjglo))
     ALLOCATE (v3d(npiglo,npjglo,npk)                                         )

     DO jk = 1, npk
        zmask(:,:) = 1.

        ztemp(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
        zsal(:,:)  = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)

        WHERE(zsal == zspval ) zmask = 0

        v3d(:,:,jk) = sigmai(ztemp, zsal, refdepth, npiglo, npjglo ) * zmask(:,:)
     END DO  ! loop to next level

     DEALLOCATE ( ztemp, zsal, zmask )
     !------------------------------------------------------------------------------
     ! 2. Projection of T,S and p on the chosen isopycnal layer (from cdfrhoproj)
     ALLOCATE ( alpha(npiglo,npjglo) )

     !! Compute coefficients
     DO ji=1,npiglo
        DO jj = 1, npjglo
           ik = 1
           !  Assume that rho (z) is increasing downward (no inversion)
           !     Caution with sigma0 at great depth !
           DO WHILE (zsigmaref >=  v3d(ji,jj,ik) .AND. ik <= npk &
                &                .AND. v3d(ji,jj,ik) /=  zspval )
              ik=ik+1
           END DO
           ik=ik-1
           ik0=ik
           IF (ik == 0) THEN
              ik=1
              alpha(ji,jj) = 0.
           ELSE IF (v3d(ji,jj,ik+1) == zspval ) THEN
              ik0=0
              alpha(ji,jj) = 0.
           ELSE
              ! ... alpha is always in [0,1]. Adding ik0 ( >=1 ) for saving space for ik0
              alpha(ji,jj)= (zsigmaref-v3d(ji,jj,ik))/(v3d(ji,jj,ik+1)-v3d(ji,jj,ik)) + ik0
           ENDIF
        END DO
     END DO
     DEALLOCATE (v3d)

     ! Working on temperature first
     ALLOCATE( ztempint(npiglo, npjglo), zint(npiglo, npjglo), pint(npiglo, npjglo) )
     ALLOCATE( ztemp3(npiglo, npjglo,npk) )

     DO jk=1,npk
        ztemp3(:,:,jk) = getvar(cf_tfil, cn_votemper, jk ,npiglo, npjglo, ktime=jt)
     ENDDO

     DO ji=1,npiglo
        DO jj=1,npjglo
           ! ik0 is retrieved from alpha, taking the integer part.
           ! The remnant is alpha. 
           ik0=INT(alpha(ji,jj))
           alpha(ji,jj) =  alpha(ji,jj) - ik0
           IF (ik0 /= 0) THEN
              P1=ztemp3(ji,jj,ik0)
              P2=ztemp3(ji,jj,ik0+1)
              IF (P1 /= zspval .AND. P2 /= zspval) THEN
                 ztempint(ji,jj) = alpha(ji,jj) * P2          + (1-alpha(ji,jj)) * P1
                 zint(ji,jj)     = alpha(ji,jj) * prof(ik0+1) + (1-alpha(ji,jj)) * prof(ik0)
              ELSE
                 ztempint(ji,jj) = zspval
                 zint(ji,jj)     = zspval
              ENDIF
           ELSE
              ztempint(ji,jj) = zspval
              zint(ji,jj)     = zspval
           ENDIF
           ! re-add ik0 to alpha for the next computation
           alpha(ji,jj) =  alpha(ji,jj) + ik0
        END DO
     END DO

     pint = zint / 10. ! pressure on the isopycnal layer = depth / 10.

     ierr = putvar(ncout, id_varout(1), ztempint, 1, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(3), zint,     1, npiglo, npjglo, ktime=jt)
     DEALLOCATE( ztemp3 )  

     ! Working on salinity
     ALLOCATE( zsalint(npiglo, npjglo) )
     ALLOCATE( zsal3(npiglo, npjglo,npk) )

     DO jk=1,npk
        zsal3(:,:,jk) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)
     ENDDO

     DO ji=1,npiglo
        DO jj=1,npjglo
           ! ik0 is retrieved from alpha, taking the integer part.
           ! The remnant is alpha. 
           ik0=INT(alpha(ji,jj))
           alpha(ji,jj) =  alpha(ji,jj) - ik0
           IF (ik0 /= 0) THEN
              P1=zsal3(ji,jj,ik0)
              P2=zsal3(ji,jj,ik0+1)
              IF (P1 /= zspval .AND. P2 /= zspval) THEN
                 zsalint(ji,jj) = alpha(ji,jj) * P2 + (1-alpha(ji,jj)) * P1
              ELSE
                 zsalint(ji,jj) = zspval
              ENDIF
           ELSE
              zsalint(ji,jj) = zspval
           ENDIF
           ! re-add ik0 to alpha for the next computation
           alpha(ji,jj) =  alpha(ji,jj) + ik0
        END DO
     END DO

     ierr = putvar(ncout, id_varout(2), zsalint, 1, npiglo, npjglo, ktime=jt)
     DEALLOCATE( zsal3 )

     ! 3. Compute means for T,S and depth on the isopycnal layer
     ALLOCATE( zmask(npiglo, npjglo) )
     zmask=1. ! define a new mask which correspond to the isopycnal layer
     WHERE( zint == 0. ) zmask = 0.

     ztmean = SUM( ztempint * e1t * e2t * zmask ) / SUM( e1t * e2t * zmask )
     zsmean = SUM( zsalint  * e1t * e2t * zmask ) / SUM( e1t * e2t * zmask )
     ! JMM rem : hmean never used ...
     ! hmean = SUM( zint     * e1t * e2t * zmask ) / SUM( e1t * e2t * zmask )
     pmean = SUM( pint     * e1t * e2t * zmask ) / SUM( e1t * e2t * zmask )

     DEALLOCATE ( ztempint, zsalint )

     ! 4. Compute specific volume anomaly
     ALLOCATE( zsva3(npiglo,npjglo,npk) )
     ALLOCATE( zsiginsitu(npiglo,npjglo), zsig0(npiglo,npjglo) )
     ALLOCATE( ztemp(npiglo,npjglo),  zsal(npiglo,npjglo) )
     ALLOCATE( ztemp0(npiglo,npjglo), zsal0(npiglo,npjglo) )

     DO jk=1,npk
        ztemp(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
        zsal (:,:) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)

        ztemp0(:,:) = ztmean
        zsal0 (:,:) = zsmean

        ! again land/sea mask
        zmask (:,:) = 1.
        WHERE( zsal == zspval ) zmask = 0.

        zsiginsitu(:,:) = sigmai ( ztemp,  zsal,  prof(jk), npiglo, npjglo ) * zmask(:,:) ! in-situ density
        zsig0(:,:)      = sigmai ( ztemp0, zsal0, prof(jk), npiglo, npjglo ) * zmask(:,:) ! density of reference profile

        zsva3(:,:,jk)    = ( 1. / zsiginsitu(:,:) ) - ( 1. / zsig0(:,:) )
     ENDDO

     DEALLOCATE( zsiginsitu, zsig0, ztemp0, zsal0 )

     ! 5. Integrates from surface to depth of isopycnal layer
     ALLOCATE( zdep(npiglo, npjglo), rdeltapsi1(npiglo, npjglo) )

     rdeltapsi1(:,:) = 0.
     DO jk=1, npk
        zdep(:,:) = getvar(cn_fzgr, 'e3t_ps', jk,npiglo,npjglo,ldiom=.TRUE.)

        ! For each point we integrate from surface to zint(ji,jj) which is the depth
        ! of the isopycnal layer

        ! If isopycnal layer depth is below the current level
        WHERE( zint >= prof(jk) ) rdeltapsi1 = rdeltapsi1 - zsva3(:,:,jk) * zdep / 10.
        ! If isopycnal layer is between current level and previous level
        WHERE( zint < prof(jk) .AND. zint > prof(jk-1) ) rdeltapsi1 = rdeltapsi1 &
             & - zsva3(:,:,jk) * ( zint - prof(jk-1) ) / 10.
     ENDDO
     ierr = putvar(ncout, id_varout(6), rdeltapsi1, 1, npiglo, npjglo, ktime=jt)
     DEALLOCATE( zdep )

     ! 6. Projection of the specific volume anomaly on the isopycnal layer
     ALLOCATE( zsva2(npiglo,npjglo), rdeltapsi2(npiglo,npjglo) )

     DO ji=1,npiglo
        DO jj=1,npjglo
           ! ik0 is retrieved from alpha, taking the integer part.
           ! The remnant is alpha. 
           ik0=INT(alpha(ji,jj))
           alpha(ji,jj) =  alpha(ji,jj) - ik0
           IF (ik0 /= 0) THEN
              P1=zsva3(ji,jj,ik0)
              P2=zsva3(ji,jj,ik0+1)
              IF (P1 /= zspval .AND. P2 /= zspval) THEN
                 zsva2(ji,jj) = alpha(ji,jj) * P2 + ( 1-alpha(ji,jj) ) * P1
              ELSE
                 zsva2(ji,jj) = zspval
              ENDIF
           ELSE
              zsva2(ji,jj) = zspval
           ENDIF
           ! re-add ik0 to alpha for the next computation
           alpha(ji,jj) =  alpha(ji,jj) + ik0
        END DO
     END DO

     rdeltapsi2 = ( pint - pmean ) * zsva2
     ierr       = putvar(ncout, id_varout(7), rdeltapsi2, 1, npiglo, npjglo, ktime=jt)
     DEALLOCATE ( zsva3, zsva2, alpha, zint, pint )

     ! 6. Finally we compute the surface streamfunction
     ALLOCATE(zssh(npiglo,npjglo) , zsigsurf(npiglo,npjglo), psi0(npiglo,npjglo) )

     ztemp   (:,:) = getvar(cf_tfil, cn_votemper, 1, npiglo, npjglo, ktime=jt)
     zsal    (:,:) = getvar(cf_tfil, cn_vosaline, 1, npiglo, npjglo, ktime=jt)
     zssh    (:,:) = getvar(cf_tfil, cn_sossheig, 1, npiglo, npjglo, ktime=jt)

     ! land/sea mask at surface
     zmask (:,:) = 1.
     WHERE( zsal == zspval ) zmask = 0.

     zsigsurf(:,:) = sigmai ( ztemp, zsal, prof(1), npiglo, npjglo ) * zmask(:,:)

     psi0 = zsigsurf * zssh * (9.81 / 1020. )
     ierr = putvar(ncout, id_varout(5), psi0, 1, npiglo, npjglo, ktime=jt)
     DEALLOCATE(zssh, zsigsurf, ztemp, zsal )

     ! 7. At least we are done with the computations
     ALLOCATE( psi(npiglo,npjglo) )
     ! final mask for output : mask the contribution of SSH where isopycn outcrops
     zmask=1.
     WHERE(rdeltapsi1 == zspval ) zmask = 0.
     psi = ( psi0 * zmask ) + rdeltapsi1 + rdeltapsi2

     ierr = putvar(ncout, id_varout(4), psi, 1, npiglo, npjglo, ktime=jt)
     DEALLOCATE( psi, psi0, rdeltapsi1, rdeltapsi2, zmask )
  END DO  ! loop to next time

  ierr = closeout(ncout)

END PROGRAM cdfisopsi
