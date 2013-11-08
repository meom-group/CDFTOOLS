PROGRAM cdfcurl
  !!======================================================================
  !!                     ***  PROGRAM  cdfcurl  ***
  !!=====================================================================
  !!  ** Purpose : Compute the curl on F-points for given gridU gridV 
  !!               files and variables
  !!
  !!  ** Method  : Use the same algorithm than NEMO
  !!
  !! History : 2.1  : 05/2005  : J.M. Molines : Original code
  !!         : 2.1  : 06/2007  : P. Mathiot   : for use with forcing fields
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jt         ! dummy loop index
  INTEGER(KIND=4)                           :: ilev               ! level to be processed
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                           :: narg, iargc        ! browse command line
  INTEGER(KIND=4)                           :: ncout, ierr        ! browse command line
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout     ! output variable properties

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2v, e1u, e1f, e2f ! horizontql metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: un, vn             ! velocity field
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zun, zvn           ! working arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rotn, fmask        ! curl and fmask
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zrotn              ! curl at T point 
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                ! time counter
  REAL(KIND=4)                              :: zmask              ! mask at T point for -T option

  CHARACTER(LEN=256)                        :: cf_ufil, cf_vfil   ! file names
  CHARACTER(LEN=256)                        :: cf_out = 'curl.nc' ! output file name
  CHARACTER(LEN=256)                        :: cv_u, cv_v         ! variable names
  CHARACTER(LEN=256)                        :: cldum              ! dummy string

  TYPE (variable), DIMENSION(1)             :: stypvar            ! structure for attibutes

  LOGICAL                                   :: lforcing = .FALSE. ! forcing flag
  LOGICAL                                   :: lchk     = .FALSE. ! flag for missing files
  LOGICAL                                   :: lperio   = .FALSE. ! flag for E-W periodicity
  LOGICAL                                   :: ltpoint  = .FALSE. ! flag for T-point output
  !!----------------------------------------------------------------------
  CALL ReadCdfNames() 

  narg = iargc()
  IF ( narg < 5 ) THEN
     PRINT *,' usage : cdfcurl U-file V-file U-var V-var lev [-T]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the curl of a vector field, at a specified level.'  
     PRINT *,'       If level is specified as 0, assume that the input files are'
     PRINT *,'       forcing files, presumably on A-grid. In this latter case, the'
     PRINT *,'       vector field is interpolated on the C-grid. In any case, the'
     PRINT *,'       curl is computed on the F-point (unless -T option is used).'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       U-file : zonal component of the vector field.'
     PRINT *,'       V-file : meridional component of the vector field.'
     PRINT *,'       U-var  : zonal component variable name'
     PRINT *,'       V-var  : meridional component variable name.'
     PRINT *,'       lev    : level to be processed. If set to 0, assume forcing file '
     PRINT *,'                in input.'
     PRINT * 
     PRINT *,'     OPTIONS :'
     PRINT *,'       -T : compute curl at T point instead of default F-point'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ', TRIM(cn_fhgr)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : socurl or socurlt (if -T option), units : s^-1'
     STOP
  ENDIF

  CALL getarg(1, cf_ufil)
  CALL getarg(2, cf_vfil)
  CALL getarg(3, cv_u   )
  CALL getarg(4, cv_v   )
  CALL getarg(5, cldum  ) ;  READ(cldum,*) ilev
  IF ( narg == 6 ) THEN
  CALL getarg(6, cldum ) 
    IF ( cldum /= '-T' ) THEN
      PRINT *, TRIM(cldum),' : unknown option ' ; STOP
    ELSE
      ltpoint=.true.
    ENDIF
  ENDIF

  lchk = chkfile(cn_fhgr ) .OR. lchk
  lchk = chkfile(cf_ufil ) .OR. lchk
  lchk = chkfile(cf_vfil ) .OR. lchk
  IF ( lchk ) STOP ! missing files

  ! define new variables for output
  stypvar(1)%cname             = 'socurl'
  IF (ltpoint) stypvar(1)%cname             = 'socurlt'
  stypvar(1)%cunits            = 's-1'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = -1000.
  stypvar(1)%valid_max         =  1000.
  stypvar(1)%clong_name        = 'Relative_Vorticity (curl)'
  stypvar(1)%cshort_name       = 'socurl'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TYX'

  ipk(1) = 1  !  2D

  npiglo = getdim(cf_ufil,cn_x)
  npjglo = getdim(cf_ufil,cn_y)
  npk    = getdim(cf_ufil,cn_z)
  npt    = getdim(cf_ufil,cn_t) 

  PRINT *, 'npiglo = ',npiglo
  PRINT *, 'npjglo = ',npjglo
  PRINT *, 'npk    = ',npk
  PRINT *, 'npt    = ',npt
  PRINT *, 'ilev   = ',ilev

  !test if lev exists
  IF ( (npk==0) .AND. (ilev > 0) ) THEN
     PRINT *, 'Problem : npk = 0 and lev > 0 STOP'
     STOP
  END IF

  ! if forcing field 
  IF ( ilev==0 .AND. npk==0 ) THEN
     lforcing=.true.
     npk = 1 ; ilev=1
     PRINT *, 'npk =0, assume 1'
  END IF

  IF ( npt==0 ) THEN
     PRINT *, 'npt=0, assume 1'
     npt=1
  END IF
  ! check files and determines if the curl will be 2D of 3D
  ! ????????????

  ! Allocate the memory
  ALLOCATE ( e1u(npiglo,npjglo) , e1f(npiglo,npjglo) )
  ALLOCATE ( e2v(npiglo,npjglo) , e2f(npiglo,npjglo) )
  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo)  )
  ALLOCATE ( zun(npiglo,npjglo) , zvn(npiglo,npjglo) )
  ALLOCATE ( rotn(npiglo,npjglo) , fmask(npiglo,npjglo) )
  ALLOCATE ( tim(npt) )
  IF ( ltpoint) ALLOCATE (zrotn(npiglo,npjglo) )

  e1u =  getvar(cn_fhgr, cn_ve1u, 1, npiglo, npjglo)
  e1f =  getvar(cn_fhgr, cn_ve1f, 1, npiglo, npjglo)
  e2v =  getvar(cn_fhgr, cn_ve2v, 1, npiglo, npjglo)
  e2f =  getvar(cn_fhgr, cn_ve2f, 1, npiglo, npjglo)

  ! use zun and zvn to store f latitude and longitude for output
  zun = getvar(cn_fhgr, cn_glamf, 1, npiglo, npjglo)
  zvn = getvar(cn_fhgr, cn_gphif, 1, npiglo, npjglo)

  ! look for  E-W periodicity
  IF ( zun(1,1) == zun(npiglo-1,1) ) lperio = .TRUE.

  ! create output fileset
  ncout = create      (cf_out, cf_ufil, npiglo, npjglo, 0                           )
  ierr  = createvar   (ncout , stypvar, 1,      ipk,    id_varout                   )
  ierr  = putheadervar(ncout,  cf_ufil, npiglo, npjglo, 0, pnavlon=zun, pnavlat=zvn )

  tim  = getvar1d(cf_ufil, cn_vtimec, npt      )
  ierr = putvar1d(ncout,   tim,       npt,  'T')
  
  DO jt=1,npt
     IF (MOD(jt,100)==0 ) PRINT *, jt,'/',npt
        ! if files are forcing fields
        zun(:,:) =  getvar(cf_ufil, cv_u, ilev ,npiglo,npjglo, ktime=jt)
        zvn(:,:) =  getvar(cf_vfil, cv_v, ilev ,npiglo,npjglo, ktime=jt)

     IF ( lforcing ) THEN ! for forcing file u and v are on the A grid
        DO ji=1, npiglo-1
           un(ji,:) = 0.5*(zun(ji,:) + zun(ji+1,:))
        END DO
        !
        DO jj=1, npjglo-1
           vn(:,jj) = 0.5*(zvn(:,jj) + zvn(:,jj+1))
        END DO
        ! end compute u and v on U and V point
     ELSE
       un(:,:) = zun(:,:)
       vn(:,:) = zvn(:,:)
     END IF

     ! compute the mask
     IF ( jt==1 ) THEN
        DO jj = 1, npjglo - 1
           DO ji = 1, npiglo - 1
              fmask(ji,jj)=0.
              fmask(ji,jj)= un(ji,jj)*un(ji,jj+1) * vn(ji,jj)*vn(ji+1,jj)
              IF (fmask(ji,jj) /= 0.) fmask(ji,jj)=1.
           ENDDO
        ENDDO
     END IF

     rotn(:,:) = 0.
     DO jj = 1, npjglo -1 
        DO ji = 1, npiglo -1   ! vector opt.
           rotn(ji,jj) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ) - e2v(ji,jj) * vn(ji,jj)    &
                &         - e1u(ji  ,jj+1) * un(ji  ,jj+1) + e1u(ji,jj) * un(ji,jj)  ) &
                &         * fmask(ji,jj) / ( e1f(ji,jj) * e2f(ji,jj) )
        END DO
     END DO

     IF ( lperio ) rotn(npiglo,:) = rotn(2, :)
     IF ( ltpoint ) THEN
       zrotn(:,:) = 0.
       DO ji = 2, npiglo
         DO jj = 2, npjglo
          zmask = fmask(ji,jj)*fmask(ji,jj-1)*fmask(ji-1,jj)*fmask(ji-1,jj-1)
          zrotn(ji,jj) = 0.25*( rotn(ji,jj) + rotn(ji,jj-1) + rotn(ji-1,jj) + rotn(ji-1,jj-1) ) * zmask
         ENDDO
       ENDDO
       IF ( lperio ) zrotn(1,:) = zrotn(npiglo, :)
       rotn(:,:) = zrotn(:,:)
       
     ENDIF
     ! write rotn on file at level k and at time jt
     ierr = putvar(ncout, id_varout(1), rotn, 1, npiglo, npjglo, ktime=jt)
  END DO
  ierr = closeout(ncout)

END PROGRAM cdfcurl

