PROGRAM cdfpsi
  !!======================================================================
  !!                     ***  PROGRAM  cdfpsi  ***
  !!=====================================================================
  !!  ** Purpose : Compute Barotropic Stream Function
  !!
  !!  ** Method  : Compute the 2D fields dtrpu, dtrpv as the integral on 
  !!               the vertical of u, v on their respective points.
  !!               Then integrate from south to north : ==> dpsiu
  !!               Then integrate from West to East   : ==> dpsiv
  !!                  (should be almost the same (if no error ))
  !!               Default (appropriate for global model): output dpsiu;
  !!               normalizes the values setting psi (jpi,jpj) = 0
  !!               If option "V" is given as last argument, output dpsiv,
  !!               normalizes values setting psi(jpi,1) = 0.
  !!               This is appropriate for North Atlantic
  !!
  !! History : 2.1  : 05/2005  : J.M. Molines : Original code
  !!           3.0  : 05/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jk, jt  ! dummy loop index
  INTEGER(KIND=4)                           :: ierr            ! working integer
  INTEGER(KIND=4)                           :: narg, iargc     ! command line 
  INTEGER(KIND=4)                           :: ijarg, ireq     ! command line
  INTEGER(KIND=4)                           :: npiglo, npjglo  ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt        ! size of the domain
  INTEGER(KIND=4)                           :: ncout           ! ncid of output file
  INTEGER(KIND=4)                           :: iiref, ijref    ! reference i j point
  INTEGER(KIND=4)                           :: nvout=1         ! number of output variables
  INTEGER(KIND=4), DIMENSION(:),ALLOCATABLE :: ipk, id_varout  ! levels and id's of output vars

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask           ! mask
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1v, e3v        ! v metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2u, e3u        ! u metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zu, zv          ! velocity components
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: glamf, gphif    ! longitude/latitude
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsshu, zsshv    ! ssh at u and v point
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zssh            ! temporary array for ssh
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim             ! time counter
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: e31d            ! 1d vertical metrics, full step case

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dtrpu, dtrpv    ! transport working arrays
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dtrpsshu        ! transport working arrays
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dtrpsshv        ! transport working arrays
  REAL(KIND=8), TARGET, DIMENSION(:,:), ALLOCATABLE :: dpsiu   ! BSF ( U computation
  REAL(KIND=8), TARGET, DIMENSION(:,:), ALLOCATABLE :: dpsiv   ! BSF (V computation )
  REAL(KIND=8), TARGET, DIMENSION(:,:), ALLOCATABLE :: dpsisshu ! BSF ( SSHU computation
  REAL(KIND=8), TARGET, DIMENSION(:,:), ALLOCATABLE :: dpsisshv ! BSF ( SSHV computation )
  REAL(KIND=8), POINTER, DIMENSION(:,:)     :: dpsi            ! point to dpsiu or dpsiv
  REAL(KIND=8), POINTER, DIMENSION(:,:)     :: dpsissh         ! point to dpsisshu or dpsisshv

  CHARACTER(LEN=256)                        :: cf_ufil         ! gridU netcdf file name
  CHARACTER(LEN=256)                        :: cf_vfil         ! gridV netcdf file name
  CHARACTER(LEN=256)                        :: cf_tfil         ! gridT netcdf file name (-ssh option)
  CHARACTER(LEN=256)                        :: cf_out='psi.nc' ! output file name
  CHARACTER(LEN=256)                        :: cv_out='sobarstf' ! output variable name
  CHARACTER(LEN=256)                        :: cv_outssh='sobarstfssh'    ! output variable name
  CHARACTER(LEN=256)                        :: cv_outotal='sobarstftotal' ! output variable name
  CHARACTER(LEN=256)                        :: cldum           ! dummy character variable
  CHARACTER(LEN=256)                        :: cglobal         ! global attribute

  TYPE(variable), DIMENSION(:), ALLOCATABLE :: stypvar         ! structure for attributes

  LOGICAL                                   :: lchk  = .FALSE. ! flag for missing files
  LOGICAL                                   :: ll_u  = .TRUE.  ! flag for U integration
  LOGICAL                                   :: ll_v  = .FALSE. ! flag for V integration
  LOGICAL                                   :: lfull = .FALSE. ! flag for full step config
  LOGICAL                                   :: lmask = .FALSE. ! flag for masking output
  LOGICAL                                   :: lmean = .FALSE. ! flag for mean U,V calculation
  LOGICAL                                   :: lopen = .FALSE. ! flag for open calculation
  LOGICAL                                   :: lssh  = .FALSE. ! flag for ssh computation
  LOGICAL                                   :: lnc4  = .FALSE. ! flag for netcdf4 cchunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfpsi U-file V-file [V] [-full ] [-mask ] [-mean] [-nc4 ] ...'
     PRINT *,'          ... [-ssh T-file ] [-open ] [-ref iref jref ] [-o OUT-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Computes the barotropic stream function (a proxy ) as the integral of '
     PRINT *,'       the transport.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       U-file  : netcdf file of zonal velocity.' 
     PRINT *,'       V-file  : netcdf file of meridional velocity.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [V] : use V field instead of U field for integration.' 
     PRINT *,'       [ -full ] : indicates a full step case. Default is partial steps.'
     PRINT *,'       [ -mask ] : mask output fields. Note that the land value is significant.'
     PRINT *,'                   It correspond to the potential on this continent.'
     PRINT *,'       [ -mean ] : save the average of the computations done with U and V.'
     PRINT *,'       [ -nc4  ] : use netcdf4 output files with chunking and deflation'
     PRINT *,'       [ -ssh T-file ] : compute the transport in the ''ssh'' layer, using '
     PRINT *,'                  surface velocities. Take the ssh from T-file specified in '
     PRINT *,'                  this option. This is a experimental option, not certified ...'
     PRINT *,'       [ -open ] : for open domain configuration. See also -ref to set '
     PRINT *,'                   reference point.'
     PRINT *,'       [ -ref iref jref ] : Set the reference point in i,j coordinates.'
     PRINT *,'                   BSF at reference point is arbitrarly set to zero.'
     PRINT *,'       [ -o  OUT-file ] : specify output file name instead of default ',TRIM(cf_out)
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fhgr),' and ', TRIM(cn_fzgr),'.'
     PRINT *,'       ', TRIM(cn_fmsk),' is required only if -mask option used.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ', TRIM(cv_out),' (m3/s )'
     PRINT *,'       If option -ssh is used, 2 additional variables are added to the file :'
     PRINT *,'                     ', TRIM(cv_outssh),' (m3/s ) : contribution of SSH'
     PRINT *,'                     ', TRIM(cv_outotal),' (m3/s ) : total BSF'
     PRINT *,'      '
     STOP
  ENDIF

  CALL SetGlobalAtt (cglobal)
  iiref = -1 ; ijref= -1

  ijarg = 1 ; ireq = 0
  DO WHILE ( ijarg <= narg )
     CALL getarg( ijarg, cldum ) ; ijarg=ijarg + 1
     SELECT CASE ( cldum )
     CASE ('-full') ; lfull = .TRUE.
     CASE ('-mask') ; lmask = .TRUE.
     CASE ('-mean') ; lmean = .TRUE.  ; ll_v=.TRUE. ; ll_u=.TRUE.
     CASE ('-ssh' ) ; lssh  = .TRUE.  ; nvout=3
        CALL getarg( ijarg, cf_tfil ) ; ijarg=ijarg + 1 
     CASE ('-nc4' ) ; lnc4  = .TRUE. 
     CASE ('-open') ; lopen = .TRUE.  ; ll_v=.TRUE. ; ll_u=.TRUE.
     CASE ('-o'   ) ; CALL getarg( ijarg, cf_out )   ; ijarg=ijarg + 1 
     CASE ('-ref') 
        CALL getarg( ijarg, cldum )   ; ijarg=ijarg + 1 ; READ(cldum,*) iiref
        CALL getarg( ijarg, cldum )   ; ijarg=ijarg + 1 ; READ(cldum,*) ijref

     CASE DEFAULT
        ireq = ireq + 1
        SELECT CASE ( ireq)
        CASE ( 1 ) ; cf_ufil = cldum
        CASE ( 2 ) ; cf_vfil = cldum
        CASE ( 3 ) ; ll_v = .TRUE. ; ll_u = .FALSE.
        CASE DEFAULT
           PRINT *, ' Too many arguments !' ; STOP
        END SELECT
     END SELECT
  ENDDO

  lchk = lchk .OR. chkfile( cn_fhgr )
  lchk = lchk .OR. chkfile( cn_fzgr )
  IF ( lmask) lchk = lchk .OR. chkfile( cn_fmsk )
  IF ( lssh ) lchk = lchk .OR. chkfile( cf_tfil )
  lchk = lchk .OR. chkfile( cf_ufil )
  lchk = lchk .OR. chkfile( cf_vfil )

  IF ( lchk ) STOP ! missing file

  npiglo = getdim (cf_ufil, cn_x)
  npjglo = getdim (cf_ufil, cn_y)
  npk    = getdim (cf_ufil, cn_z)
  npt    = getdim (cf_ufil, cn_t)

  IF ( iiref == -1 .OR. ijref == -1 ) THEN
     iiref=npiglo
     ijref=npjglo
  ENDIF

  ALLOCATE (stypvar(nvout), ipk(nvout), id_varout(nvout))
  ! define new variables for output ( must update att.txt)
  ipk(:) = 1  !  2D ( X, Y , T )
  stypvar(:)%cunits            = 'm3/s'
  stypvar(:)%valid_min         = -300.e6
  stypvar(:)%valid_max         = 300.e6
  stypvar(:)%conline_operation = 'N/A'
  stypvar(:)%caxis             = 'TYX'
  DO ji=1,nvout
    stypvar(ji)%ichunk = (/npiglo,MAX(1,npjglo/30),1,1 /)
  ENDDO

  stypvar(1)%cname             = cv_out
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%clong_name        = 'Barotropic_Stream_Function'
  stypvar(1)%cshort_name       = cv_out

  IF ( lssh ) THEN
     stypvar(2)%cname             = cv_outssh
     stypvar(2)%rmissing_value    = 0.
     stypvar(2)%clong_name        = 'Barotropic_Stream_Function SSH contribution'
     stypvar(2)%cshort_name       = cv_outssh

     stypvar(3)%cname             = cv_outotal
     stypvar(3)%rmissing_value    = 0.
     stypvar(3)%clong_name        = 'Barotropic_Stream_Function SSH total'
     stypvar(3)%cshort_name       = cv_outotal
  ENDIF

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt
  PRINT *, ' Option is use :'
  PRINT *, '    -full :', lfull
  PRINT *, '    -mask :', lmask
  PRINT *, '    -mean :', lmean
  PRINT *, '    -ssh  :', lssh
  PRINT *, '    -open :', lopen
  PRINT *, '    -ref  :', iiref, ijref
  PRINT *, '  U-comp  :', ll_u
  PRINT *, '  V-comp  :', ll_v

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo)                 )
  ALLOCATE ( e1v(npiglo,npjglo),e3v(npiglo,npjglo))
  ALLOCATE ( e2u(npiglo,npjglo),e3u(npiglo,npjglo))
  ALLOCATE ( zu(npiglo,npjglo),dtrpu(npiglo,npjglo), dpsiu(npiglo,npjglo) )
  ALLOCATE ( zv(npiglo,npjglo),dtrpv(npiglo,npjglo), dpsiv(npiglo,npjglo) )
  ALLOCATE ( glamf(npiglo,npjglo), gphif(npiglo,npjglo))
  ALLOCATE ( tim(npt))
  IF ( lfull ) ALLOCATE ( e31d(npk))
  IF ( lssh  ) ALLOCATE ( zssh(npiglo,npjglo), zsshu(npiglo,npjglo), zsshv(npiglo,npjglo))
  IF ( lssh  ) ALLOCATE ( dpsisshu(npiglo,npjglo), dpsisshv(npiglo,npjglo) )
  IF ( lssh  ) ALLOCATE ( dtrpsshu(npiglo,npjglo), dtrpsshv(npiglo,npjglo) )

  glamf(:,:) = getvar(cn_fhgr, cn_glamf, 1, npiglo, npjglo)
  gphif(:,:) = getvar(cn_fhgr, cn_gphif, 1, npiglo, npjglo)

  ! create output fileset
  ncout = create      (cf_out, cf_ufil, npiglo, npjglo, 1                                , ld_nc4=lnc4 )
  ierr  = createvar   (ncout,  stypvar, nvout,  ipk,    id_varout, cdglobal=TRIM(cglobal), ld_nc4=lnc4 )
  ierr  = putheadervar(ncout,  cf_ufil, npiglo, npjglo, 1, glamf, gphif)

  tim  = getvar1d(cf_ufil, cn_vtimec, npt     )
  ierr = putvar1d(ncout,   tim,       npt, 'T')

  e1v(:,:)   = getvar(cn_fhgr, cn_ve1v, 1, npiglo, npjglo)
  e2u(:,:)   = getvar(cn_fhgr, cn_ve2u, 1, npiglo, npjglo)
  IF ( lmask) THEN
     zmask(:,:) = getvar(cn_fmsk, 'fmask', 1, npiglo, npjglo)
     WHERE ( zmask >= 2 ) zmask = 1
  ENDIF

  IF ( lfull) e31d(:)    = getvare3(cn_fzgr, cn_ve3t, npk )
  ! get rid of the free-slip/no-slip condition

  DO jt = 1, npt
     dtrpu(:,:)= 0.d0
     dtrpv(:,:)= 0.d0
     dpsiu(:,:)= 0.d0
     dpsiv(:,:)= 0.d0
     IF ( lssh ) THEN
       zsshu(:,:) = 0.0
       zsshv(:,:) = 0.0
       dpsisshu(:,:) = 0.d0
       dpsisshv(:,:) = 0.d0
       dtrpsshu(:,:) = 0.d0
       dtrpsshv(:,:) = 0.d0
       zssh(:,:) = getvar(cf_tfil, cn_sossheig, 1, npiglo, npjglo, ktime=jt)
       zsshu(1:npiglo-1,    :     ) = 0.5*( zssh(2:npiglo,:       ) + zssh(1:npiglo-1,:         ))
       zsshv(   :      ,1:npjglo-1) = 0.5*( zssh(:       ,2:npjglo) + zssh(:         ,1:npjglo-1))
     ENDIF


     DO jk = 1,npk
        IF ( ll_v ) THEN
           zv(:,:) = getvar(cf_vfil, cn_vomecrty, jk, npiglo, npjglo, ktime=jt )
           IF ( lfull ) THEN ; e3v(:,:) = e31d(jk)
           ELSE              ; e3v(:,:) = getvar(cn_fzgr, 'e3v_ps', jk, npiglo, npjglo, ldiom=.TRUE.)
           ENDIF
           dtrpv(:,:) = dtrpv(:,:) + zv(:,:)*e1v(:,:)*e3v(:,:)*1.d0  ! meridional transport of each grid cell
           IF ( lssh .AND. (jk == 1 ) ) THEN
             dtrpsshv(:,:) = dtrpsshv(:,:) + zv(:,:)*e1v(:,:)*zsshv(:,:)*1.d0  ! meridional transport of each grid cell
           ENDIF
        ENDIF

        IF ( ll_u) THEN
           zu(:,:) = getvar(cf_ufil, cn_vozocrtx, jk, npiglo, npjglo, ktime=jt )
           IF ( lfull ) THEN ; e3u(:,:) = e31d(jk)
           ELSE              ; e3u(:,:) = getvar(cn_fzgr, 'e3u_ps', jk, npiglo, npjglo, ldiom=.TRUE.)
           ENDIF
           dtrpu(:,:) = dtrpu(:,:) + zu(:,:)*e2u(:,:)*e3u(:,:)*1.d0  ! zonal transport of each grid cell
           IF ( lssh .AND. (jk == 1 ) ) THEN
             dtrpsshu(:,:) = dtrpsshu(:,:) + zv(:,:)*e2u(:,:)*zsshu(:,:)*1.d0  ! meridional transport of each grid cell
           ENDIF
        ENDIF
     END DO  ! loop to next level

     IF ( lopen ) THEN
        ! This case corresponds to arbitrary configuration: we chose to compute the transport
        ! across a first line ( eg, ji= 2 or jj= npjglo-1 ), assuming that this starting line is
        ! in the true ocean. If it is on true land, it is not a problem. But it cannot be on
        ! arbitrary masked points....
        IF ( lssh ) THEN
          dpsisshu(1,npjglo-2) = dtrpsshv(1, npjglo-2)
          DO ji = 2, npiglo
             dpsisshu(ji,npjglo-2) = dpsisshu(ji-1,npjglo-2) + dtrpsshv(ji,npjglo-2)
          END DO
          ! Then compute the transport with along U starting from this line
          DO jj= npjglo-3,1,-1
             DO ji = 1, npiglo
                dpsisshu(ji,jj) = dpsisshu(ji,jj+1) + dtrpsshu(ji,jj+1)
             END DO
          END DO
        ENDIF

        dpsiu(1,npjglo-2) = dtrpv(1, npjglo-2)
        DO ji = 2, npiglo
           dpsiu(ji,npjglo-2) = dpsiu(ji-1,npjglo-2) + dtrpv(ji,npjglo-2)
        END DO
        ! Then compute the transport with along U starting from this line
        DO jj= npjglo-3,1,-1
           DO ji = 1, npiglo
              dpsiu(ji,jj) = dpsiu(ji,jj+1) + dtrpu(ji,jj+1)
           END DO
        END DO

        IF ( lmean ) THEN  ! we need also the other estimate
           dpsiv(npiglo-2, npjglo) = dtrpu(npiglo-2, npjglo)
           DO jj= npjglo - 1, 1, -1
              dpsiv(npiglo-2,jj) = dpsiv(npiglo-2, jj+1) + dtrpu(npiglo-2, jj+1)
           END DO
           DO jj=npjglo,1,-1
              DO ji = npiglo -3,1,-1
                 dpsiv(ji,jj) = dpsiv(ji+1,jj) - dtrpv(ji+1,jj)
              END DO
           END DO
           dpsiu(:,:) = 0.5*(dpsiu(:,:) + dpsiv(:,:))

           IF ( lssh ) THEN
              dpsisshv(npiglo-2, npjglo) = dtrpsshu(npiglo-2, npjglo)
              DO jj= npjglo - 1, 1, -1
                 dpsisshv(npiglo-2,jj) = dpsisshv(npiglo-2, jj+1) + dtrpsshu(npiglo-2, jj+1)
              END DO
              DO jj=npjglo,1,-1
                 DO ji = npiglo -3,1,-1
                    dpsisshv(ji,jj) = dpsisshv(ji+1,jj) - dtrpsshv(ji+1,jj)
                 END DO
              END DO
              dpsisshu(:,:) = 0.5*(dpsisshu(:,:) + dpsisshv(:,:))
           ENDIF
        ENDIF

        dpsi => dpsiu
        IF ( lssh ) dpsissh => dpsisshu

     ELSE
        ! now perform zonal integration if requested
        IF ( ll_v ) THEN
           ! integrate zonally from east to west 
           ! This comfortable with NATL configurations as the eastern most points are land points.
           dpsiv(npiglo,:)= 0.d0
           DO ji=npiglo-1,1,-1
              dpsiv(ji,:) = dpsiv(ji+1,:) - dtrpv(ji,:)  ! psi at f point
           END DO
           dpsi => dpsiv
           IF ( lssh ) THEN
              dpsisshv(npiglo,:)= 0.d0
              DO ji=npiglo-1,1,-1
                 dpsisshv(ji,:) = dpsisshv(ji+1,:) - dtrpsshv(ji,:)  ! psissh at f point
              END DO
              dpsissh => dpsisshv
           ENDIF
        ENDIF

        ! now perform meridional integration if requested
        IF ( ll_u ) THEN
           ! integrate from the south to the north with zonal transport
           ! This is because on global configuration, line jj=1 is always land (Antarctic)
           dpsiu(:,:) = 0.d0
           DO jj = 2, npjglo
              dpsiu(:,jj) = dpsiu(:,jj-1) - dtrpu(:,jj)   ! psi at f point
           END DO
           dpsi => dpsiu
           IF ( lssh ) THEN
              dpsisshu(:,:) = 0.d0
              DO jj = 2, npjglo
                 dpsisshu(:,jj) = dpsisshu(:,jj-1) - dtrpsshu(:,jj)   ! psissh at f point
              END DO
                 dpsissh => dpsisshu
           ENDIF
        ENDIF

        IF ( lmean) THEN 
           dpsiu(:,:) = 0.5 * ( dpsiu(:,:) + dpsiv(:,:) )
           dpsi => dpsiu
           IF ( lssh ) THEN
              dpsisshu(:,:) = 0.5 * ( dpsisshu(:,:) + dpsisshv(:,:) )
              dpsissh => dpsisshu
           ENDIF
        ENDIF
     ENDIF

     ! output results after normalization
     dpsi = dpsi - dpsi(iiref,ijref)
     IF ( lmask ) THEN
        PRINT *,' Write masked BSF'
        ierr = putvar(ncout, id_varout(1), SNGL(dpsi)*zmask(:,:), 1, npiglo, npjglo, ktime=jt)
        IF ( lssh ) THEN
           ierr = putvar(ncout, id_varout(2), SNGL(dpsissh     )*zmask(:,:), 1, npiglo, npjglo, ktime=jt)
           ierr = putvar(ncout, id_varout(3), SNGL(dpsissh+dpsi)*zmask(:,:), 1, npiglo, npjglo, ktime=jt)
        ENDIF
     ELSE
        PRINT *,' Write BSF'
        ierr = putvar(ncout, id_varout(1), SNGL(dpsi)           , 1, npiglo, npjglo, ktime=jt)
        IF ( lssh ) THEN
           ierr = putvar(ncout, id_varout(2), SNGL(dpsissh     ), 1, npiglo, npjglo, ktime=jt)
           ierr = putvar(ncout, id_varout(3), SNGL(dpsissh+dpsi), 1, npiglo, npjglo, ktime=jt)
        ENDIF
     ENDIF
  ENDDO

  ierr = closeout (ncout)

END PROGRAM cdfpsi
