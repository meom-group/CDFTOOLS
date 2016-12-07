PROGRAM cdftransig_xy3d
  !!======================================================================
  !!                     ***  PROGRAM  cdftransig_xy3d  ***
  !!=====================================================================
  !!  ** Purpose : Calculates u and v transports at each grid cell
  !!               in rho coordinates.  produces a 3D field.
  !!
  !!  ** Method  : allow two 3D arrays for more efficient reading
  !!
  !! History : 2.1  : 02/2006  : A.M. Treguier : Original code
  !!           2.1  : 02/2011  : A.M. Treguier : Allow increased resolution in density
  !!                                             in deeper layers
  !!           3.0  : 06/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  USE modutils      ! for SetFileName, SetGlobalAtt
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                              :: ji, jj, jk         ! dummy loop index
  INTEGER(KIND=4)                              :: jt, jtag           ! dummy loop index
  INTEGER(KIND=4)                              :: ijb   
  INTEGER(KIND=4)                              :: ierr               ! working integer
  INTEGER(KIND=4)                              :: narg, iargc        ! command line 
  INTEGER(KIND=4)                              :: ijarg, ireq, istag
  INTEGER(KIND=4)                              :: iset
  INTEGER(KIND=4)                              :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                              :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                              :: nbins              ! density  bins
  INTEGER(KIND=4)                              :: ncout
  INTEGER(KIND=4)                              :: ntags, nframes
  INTEGER(KIND=4)                              :: nsigmax , ijtrans  ! dimension for itab, intermediate index
  INTEGER(KIND=4), DIMENSION(2)                :: id_varout , ipk    !
  INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: itab               ! look up table for density intervals
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ibinu, ibinv       ! integer value corresponding to density for binning

  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: zmasku,zmaskv      !  masks x,1,nbins
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e1v,  gphiv        !  2D x,y metrics, velocity
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e2u                !  metrics, velocity
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: zt, zs, zv, e3v    !  x,1,z arrays metrics, velocity
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: zu, e3u            !  metrics, velocity
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: zdensu, zdensv     ! density on u and v points 
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: gdept              ! array for depth of T points  
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: e31d               ! vertical metric in case of full step
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: tim
  REAL(KIND=4), DIMENSION(1)                   :: timean
  REAL(KIND=4)                                 :: pref               !  reference for density 

  REAL(KIND=8), DIMENSION(:,:,:),  ALLOCATABLE :: dusigsig,dvsigsig  ! cumulated transports,   
  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: dens2d 
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dsigma             ! density coordinate, center of bins
  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dsig_edge          ! density coordinate, edge of bins.
  REAL(KIND=8)                                 :: dsigtest
  REAL(KIND=8)                                 :: ds1min, ds1scal    !  min sigma and delta_sigma
  REAL(KIND=8)                                 :: ds1zoom = 999., ds1scalmin  !  min sigma for increased resolution
  REAL(KIND=8)                                 :: dtotal_time

  CHARACTER(LEN=80 )                           :: cf_out='uvxysig.nc'
  CHARACTER(LEN=80 )                           :: cf_tfil
  CHARACTER(LEN=80 )                           :: cf_ufil
  CHARACTER(LEN=80 )                           :: cf_vfil
  CHARACTER(LEN=80 )                           :: cv_outu='vouxysig'
  CHARACTER(LEN=80 )                           :: cv_outv='vovxysig'
  CHARACTER(LEN=80 )                           :: config
  CHARACTER(LEN=80 )                           :: ctag 
  CHARACTER(LEN=80 )                           :: cldum
  CHARACTER(LEN=80 )                           :: cldepcode='1000'
  CHARACTER(LEN=256)                           :: cglobal
  CHARACTER(LEN=7  )                           :: clsigma 

  TYPE (variable), DIMENSION(2)                :: stypvar     ! structure for attributes

  LOGICAL                                      :: lprint  = .FALSE.
  LOGICAL                                      :: lfull   = .FALSE.
  LOGICAL                                      :: lnotset = .FALSE.
  LOGICAL                                      :: lchk    = .FALSE.  ! flag for missing files
  LOGICAL                                      :: lperio  = .FALSE.  ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdftransig_xyz CONFCASE ''list_of_tags'' [-depref depcode ] ...'
     PRINT *,'                    ... [-depref depref ] [ -nbins nbins ] ... ' 
     PRINT *,'                    ... [-sigmin smin s-scal] [-sigzoom sminr s-scalr ] ...'
     PRINT *,'                    ... [-full ] [-v ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Computes the volume transport at each grid cell in density space ' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       CONFCASE  : a DRAKKAR CONFIG-CASE name '
     PRINT *,'       list_of_tags : a list of time tags to be processed'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-depcode depcode ] : depcode corresponds to pre-defined parameter '
     PRINT *,'              setting, in term of reference depths, density limits for '
     PRINT *,'              binning, number of bins, deeper layer refinement.'
     PRINT *,'          AVAILABLE depcode are :'
     PRINT *,'              _______________________________________________________________'
     PRINT *,'              depcode  |  depth_ref nbins   smin  s-scal  szoommin szoom-scal '
     PRINT *,'              ---------------------------------------------------------------'
     PRINT *,'                0      |    0      101     23.0    0.05                      '
     PRINT *,'              1000     |   1000     93     24.2    0.10   32.3     0.05      ' 
     PRINT *,'              1000-acc |   1000     88     24.5    0.10                      '
     PRINT *,'              2000     |   2000    174     29.0    0.05                      '
     PRINT *,'              none     |  parameters must be set individually                '
     PRINT *,'              ---------------------------------------------------------------'
     PRINT *,'          DEFAULT depcode is : ',TRIM(cldepcode)
     PRINT *,'              For other setting use the options to specify the settings'
     PRINT *,'              individually.'
     PRINT *,'       [-depref depref ] : give the depth reference for potential density'
     PRINT *,'       [-nbins nbins ] : give the number of density bins.'
     PRINT *,'       [-sigmin smin s-scal ] : give the minimum of density for binning and'
     PRINT *,'              the bin width.  ( take care of the reference depth).'
     PRINT *,'       [-sigzoom sminr s-scalr ] : allow density refinement from sminr, with'
     PRINT *,'              s-scalr bin width.'
     PRINT *,'       [-full ] : indicate a full step configuration.'
     PRINT *,'       [-v ] : verbose mode : extra print are performed during execution.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fhgr),' and ',TRIM(cn_fzgr) 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ',TRIM(cv_outu),' and ', TRIM(cv_outv),' in m3/s.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      cdfrhoproj, cdfsigtrp' 
     PRINT *,'      '
     STOP
  ENDIF

  ! browse command line according to options
  ijarg = 1 ; ireq = 0 ; ntags = 0 ; iset = 0
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-depcode' ) ; CALL getarg(ijarg, cldepcode ) ; ijarg=ijarg+1
     CASE ( '-depref'  ) ; CALL getarg(ijarg, cldum     ) ; ijarg=ijarg+1 ; READ(cldum,*) pref       ; iset = iset+1
                           WRITE(clsigma,'("sigma_",I1)'), NINT(pref/1000.)
     CASE ( '-nbins'   ) ; CALL getarg(ijarg, cldum     ) ; ijarg=ijarg+1 ; READ(cldum,*) nbins      ; iset = iset+1
     CASE ( '-sigmin'  ) ; CALL getarg(ijarg, cldum     ) ; ijarg=ijarg+1 ; READ(cldum,*) ds1min
     CASE ( '-sigzoom' ) ; CALL getarg(ijarg, cldum     ) ; ijarg=ijarg+1 ; READ(cldum,*) ds1zoom    ; iset = iset+1
                           CALL getarg(ijarg, cldum     ) ; ijarg=ijarg+1 ; READ(cldum,*) ds1scalmin
     CASE ( '-full'    ) ; lfull  = .TRUE.
     CASE ( '-v'       ) ; lprint = .TRUE.
     CASE DEFAULT    ! mandatory arguments 
        ireq=ireq+1      
        SELECT CASE (ireq)
        CASE ( 1 ) ;  config=cldum
        CASE DEFAULT
         IF ( ntags == 0 ) istag = ijarg - 1 ! remember the argument number corresponding to 1rst tag
         ntags=ntags + 1
        END SELECT
     END SELECT
  ENDDO

  ! set parameters for pre-defined depcode
  SELECT CASE ( cldepcode )
  CASE ( '0'                    ) 
      pref = 0.    ; nbins = 101 ; ds1min = 23.0 ; ds1scal = 0.03 ; ds1zoom = 999. ; ds1scalmin = 999. ; clsigma='sigma_0'
  CASE ( '1000'                 ) 
      pref = 1000. ; nbins =  93 ; ds1min = 24.2d0 ; ds1scal = 0.10d0 ; ds1zoom = 32.3d0 ; ds1scalmin = 0.05d0 ; clsigma='sigma_1'
  CASE ( '1000-acc', '1000-ACC' ) 
      pref = 1000. ; nbins =  88 ; ds1min = 24.5 ; ds1scal = 0.10 ; ds1zoom = 999. ; ds1scalmin = 999. ; clsigma='sigma_1'
  CASE ( '2000'                 ) 
      pref = 2000. ; nbins = 174 ; ds1min = 29.0 ; ds1scal = 0.05 ; ds1zoom = 999. ; ds1scalmin = 999. ; clsigma='sigma_2'
  CASE ( 'none'                 ) 
      ! in this case check that all parameters are set individually
      IF ( iset /= 3  ) THEN 
         PRINT *, ' You must set depref, nbins, sigmin  individually' ; STOP
      ENDIF
  CASE DEFAULT 
      PRINT *, ' this depcode :',TRIM(cldepcode),' is not available.' ; STOP
  END SELECT

  ds1scalmin = MIN ( ds1scalmin, ds1scal )
  IF ( lprint ) THEN
     PRINT *,' DEP REF  : ', pref, ' m'
     PRINT *,' NBINS    : ', nbins
     PRINT *,' SIGMIN   : ', ds1min
     PRINT *,' SIGSTP   : ', ds1scal
     PRINT *,' SIGIN R  : ', ds1zoom
     PRINT *,' SIGSTP R : ', ds1scalmin
  ENDIF
  ! use first tag to look for file dimension
  CALL getarg (istag, ctag)
  cf_vfil = SetFileName (config, ctag, 'V' )
  IF ( chkfile(cf_vfil) ) STOP ! missing file

  npiglo = getdim (cf_vfil, cn_x)
  npjglo = getdim (cf_vfil, cn_y)
  npk    = getdim (cf_vfil, cn_z)

  ALLOCATE ( dsigma(nbins), dsig_edge(nbins+1) )
  ! define densities at middle of bins and edges of bins
  ijtrans = 0
  DO ji=1,nbins
     dsigtest  = ds1min +(ji-0.5)*ds1scal
     IF ( dsigtest > ds1zoom ) THEN
        IF ( ijtrans == 0 ) ijtrans = ji
        dsigma(ji) = ds1zoom + (ji-ijtrans+0.5)*ds1scalmin
     ELSE
        dsigma(ji) = dsigtest
     ENDIF
  ENDDO

  IF (lprint) PRINT *, ' min density:',dsigma(1), ' max density:', dsigma(nbins)
  IF (lprint) PRINT *, ' verify sigma:', dsigma

  dsig_edge(1) = ds1min
  DO ji=2,nbins 
     dsig_edge(ji) = 0.5* (dsigma(ji)+dsigma(ji-1))
  END DO
  dsig_edge(nbins+1) = dsig_edge(nbins) + ds1scalmin
  IF (lprint) PRINT *, ' sig_edge : ', dsig_edge
  !
  !  define a lookup table array so that the density can be binned according to 
  !  the smallest interval ds1scalmin
  nsigmax = NINT( (dsig_edge(nbins+1)-dsig_edge(1))/ds1scalmin ) !+1
  ALLOCATE ( itab(nsigmax))
  itab(:) = 0
  DO ji=1,nsigmax
     dsigtest = ds1min+ (ji-0.5) * ds1scalmin
     DO jj=1,nbins
        IF ( dsigtest > dsig_edge(jj) .AND. dsigtest <= dsig_edge(jj+1) ) THEN
           itab(ji) = jj
        ENDIF
     END DO
  ENDDO
  IF (lprint) PRINT *, ' nsigmax=' , nsigmax
  IF (lprint) PRINT *, ' verify itab:', itab

  ! define new variables for output ( must update att.txt)
  ! define output variables
  CALL SetGlobalAtt(cglobal)

  ipk(:)                    = nbins   ! output file has  nbins sigma values
  stypvar%cunits            = 'm3/s'  ! transports
  stypvar%rmissing_value    = 0.
  stypvar%valid_min         = -10.    ! seem to be small
  stypvar%valid_max         = 10.
  stypvar%conline_operation = 'N/A'
  stypvar%caxis             = 'TSYX'

  stypvar(1)%cname          = cv_outu                ;    stypvar(2)%cname       = cv_outv
  stypvar(1)%clong_name     = 'Zonal_trsp_sig_coord' ;    stypvar(2)%clong_name  = 'Meridional_trsp_sig_coord'
  stypvar(1)%cshort_name    = cv_outu                ;    stypvar(2)%cshort_name = cv_outv

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'nbins  = ', nbins

  ! Allocate arrays
  ALLOCATE ( zv (npiglo,npjglo), zu (npiglo,npjglo) )
  ALLOCATE ( zt (npiglo,npjglo), zs (npiglo,npjglo) )
  ALLOCATE ( e3v(npiglo,npjglo), e3u(npiglo,npjglo) )
  ALLOCATE ( ibinu(npiglo, npjglo), ibinv(npiglo, npjglo) )
  ALLOCATE ( e1v(npiglo,npjglo), gphiv(npiglo,npjglo), gdept(npk) )
  ALLOCATE ( e2u(npiglo,npjglo) )
  ALLOCATE ( zdensu(npiglo,npjglo), zdensv(npiglo,npjglo) )
  ALLOCATE ( zmasku(npiglo,npjglo), zmaskv(npiglo,npjglo))
  ALLOCATE ( dusigsig(npiglo,npjglo,nbins), dvsigsig(npiglo,npjglo,nbins))  ! huge as nbins can be > 100
  ALLOCATE ( dens2d(npiglo,npjglo) )

  e1v(:,:)   = getvar  (cn_fhgr, cn_ve1v,  1, npiglo, npjglo) 
  e2u(:,:)   = getvar  (cn_fhgr, cn_ve2u,  1, npiglo, npjglo) 
  gphiv(:,:) = getvar  (cn_fhgr, cn_gphiv, 1, npiglo, npjglo)
  gdept(:)   = getvare3(cn_fzgr, cn_gdept, npk              )

  ! look for  E-W periodicity (using zu for temporary array
  zu(:,:)    = getvar  (cn_fhgr, cn_glamv, 1, npiglo, npjglo)
  IF ( zu(1,1) == zu(npiglo-1,1) ) lperio = .TRUE.

  IF ( lfull ) THEN 
     ALLOCATE ( e31d(npk) )
     e31d(:)  = getvare3(cn_fzgr, cn_ve3t, npk              )
  ENDIF

  ! create output fileset
  IF (lprint) PRINT *, ' ready to create file:',TRIM( cf_out), ' from reference:',TRIM(cf_vfil )
  ncout = create      (cf_out, cf_vfil, npiglo, npjglo, nbins,    cdep=clsigma            )
  ierr  = createvar   (ncout,  stypvar, 2,      ipk,    id_varout, cdglobal=TRIM(cglobal) )
  ierr  = putheadervar(ncout,  cf_vfil, npiglo, npjglo, nbins,     pdep=REAL(dsigma)      )

  dtotal_time = 0.d0

  ! initialize transport to 0
  dusigsig (:,:,:) = 0.d0 ; dvsigsig (:,:,:) = 0.d0;
  !    loop on time and depth ---------------------------------------------------
  ! 
  DO jk= 1, npk-1
     IF ( lprint ) PRINT *, ' working on depth jk=',jk
     IF ( lfull ) THEN 
       e3v(:,:) = e31d(jk)
       e3u(:,:) = e31d(jk)
     ELSE
       e3v(:,:) = getvar(cn_fzgr, 'e3v_ps', jk, npiglo,npjglo )
       e3u(:,:) = getvar(cn_fzgr, 'e3u_ps', jk, npiglo,npjglo )
     ENDIF

     ijarg = istag ; nframes = 0
     DO jtag = 1, ntags

        CALL getarg (ijarg, ctag) ; ijarg = ijarg + 1
        IF (lprint   ) PRINT *, ' working on  ctag=',TRIM(ctag)
        cf_tfil = SetFileName(config, ctag, 'T') 
        cf_ufil = SetFileName(config, ctag, 'U')
        cf_vfil = SetFileName(config, ctag, 'V')

        ! check existence of files
        lchk =           chkfile ( cf_tfil) 
        lchk = lchk .OR. chkfile ( cf_ufil) 
        lchk = lchk .OR. chkfile ( cf_vfil) 
        IF ( lchk ) STOP ! missing file
        
        IF (jk== 1 ) THEN
           npt = getdim (cf_tfil, cn_t)  ! assuming all files (U V ) contains same number of time frame
           ALLOCATE ( tim(npt) )
           tim         = getvar1d(cf_tfil, cn_vtimec, npt )
           dtotal_time = dtotal_time + SUM( DBLE(tim) )
           DEALLOCATE ( tim )
        ENDIF

        DO jt = 1, npt
           nframes = nframes + 1
           ! Get velocities u, v  and mask   if first time slot only 
           zv(:,:)= getvar ( cf_vfil, cn_vomecrty, jk, npiglo, npjglo, ktime=jt )
           zu(:,:)= getvar ( cf_ufil, cn_vozocrtx, jk, npiglo, npjglo, ktime=jt)
           IF (jtag == 1) THEN 
              zmasku(:,:)= 1; zmaskv(:,:)= 1.0 ;
              WHERE( zu == 0) zmasku(:,:)= 0.0 ;
              WHERE( zv == 0) zmaskv(:,:)= 0.0;
              IF (lprint  ) PRINT *, ' min,max u:',MINVAL(zu),MAXVAL(zu)
           ENDIF
           !                     density  
           zt(:,:)= getvar ( cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt )
           zs(:,:)= getvar ( cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt )

           IF ( pref == 0. ) THEN
              dens2d = sigma0(zt, zs,       npiglo, npjglo)
           ELSE
              dens2d = sigmai(zt, zs, pref, npiglo, npjglo)
           ENDIF

           !  density on u points masked by u  , single precision 
           zdensu(1:npiglo-1,:) = 0.5*( dens2d(1:npiglo-1,:) + dens2d(2:npiglo,:))

           ! check for periodic EW condition
           IF ( lperio ) THEN ; zdensu(npiglo,:)     = zdensu(2,:)
           ELSE               ; zdensu(npiglo,:) = 0.
           ENDIF

           zdensu(:,:) = zdensu(:,:) * zmasku(:,:)

           !  density on v points masked by v  , single precision
           zdensv(:,1:npjglo-1) = 0.5*( dens2d(:,1:npjglo-1) + dens2d(:,2:npjglo) )
           zdensv(:,:) = zdensv(:,:) * zmaskv(:,:)

           !  bins density - bins based on dens2d 
           DO jj=1,npjglo
              DO ji=1,npiglo
                 ijb   = INT( (zdensu(ji,jj) - ds1min)/ds1scalmin )+1
                 ijb   = MAX( ijb ,1   )
                 ijb   = MIN( ijb,nsigmax)
                 ibinu(ji,jj) = itab (ijb)
                 ijb   = INT( (zdensv(ji,jj) - ds1min)/ds1scalmin )+1
                 ijb   = MAX( ijb ,1   )
                 ijb   = MIN( ijb,nsigmax)
                 ibinv(ji,jj) =  itab(ijb)
              ENDDO
           ENDDO
           zu(:,:) = zu(:,:)*e3u(:,:)
           zv(:,:) = zv(:,:)*e3v(:,:)
           DO jj=1,npjglo
              DO ji=1,npiglo
                 dusigsig(ji,jj,ibinu(ji,jj)) = dusigsig(ji,jj,ibinu(ji,jj))+ e2u(ji,jj)*zu(ji,jj)*1.d0
                 dvsigsig(ji,jj,ibinv(ji,jj)) = dvsigsig(ji,jj,ibinv(ji,jj))+ e1v(ji,jj)*zv(ji,jj)*1.d0            
              END DO
           END DO
        END DO  ! end of loop on file time frame
        !  -----------------------------------------end of loop on ctags
     END DO
     !           
     ! -----------------  end of loop on jk
  END DO

  timean(1) = dtotal_time/nframes
  ierr      = putvar1d(ncout, timean, 1, 'T')
  DO jk=1, nbins
     zt   = dusigsig(:,:,jk) / nframes
     ierr = putvar (ncout, id_varout(1), zt, jk, npiglo, npjglo, kwght=nframes)
  ENDDO
  DO jk=1, nbins
     zt   = dvsigsig(:,:,jk) / nframes
     ierr = putvar (ncout, id_varout(2), zt, jk, npiglo, npjglo, kwght=nframes)
  ENDDO
  ierr = closeout(ncout)

END PROGRAM cdftransig_xy3d

   
