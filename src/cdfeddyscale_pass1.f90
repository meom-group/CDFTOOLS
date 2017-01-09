PROGRAM cdfeddyscale_pass1
   !!======================================================================
   !!                     ***  PROGRAM  cdfeddyscale_pass1  ***
   !!=====================================================================
   !!  ** Purpose : Compute: - the curl and the square of curl on F-points, 
   !!                        - the gradient components of the curl and the
   !!                          square of the gradient components on UV-points
   !!                        - the square of velocity components on UV-points 
   !!               for given gridU gridV files and variables.
   !!               These terms will used to compute the Taylor scale or large
   !!               scale eddy (lambda1) and the small scale eddy (lambda2) 
   !!               in the program cdflambda.f90.
   !!
   !!  ** Method  : Use the same algorithm than NEMO
   !!
   !!----------------------------------------------------------------------
   USE cdfio
   USE modcdfnames
   !!----------------------------------------------------------------------
   !! CDFTOOLS_3.0 , MEOM 2013
   !! $Id$
   !! Copyright (c) 2013, C. Q. C. Akuetevi & J.-M. Molines
   !! Software governed by the CeCILL licence
   !! (Licence/CDFTOOLSCeCILL.txt)
   !!----------------------------------------------------------------------
   IMPLICIT NONE
   ! index of output variables
   INTEGER(KIND=4),                PARAMETER :: jp_nvar=8
   INTEGER(KIND=4),                PARAMETER :: jp_curl=1, jp_curl2=2, jp_dxcurl=3, jp_dycurl=4
   INTEGER(KIND=4),                PARAMETER :: jp_dxcurl2=5, jp_dycurl2=6, jp_u2=7, jp_v2=8

   INTEGER(KIND=4)                           :: ji, jj, jt         ! dummy loop index
   INTEGER(KIND=4)                           :: ilev               ! level to be processed
   INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
   INTEGER(KIND=4)                           :: npk, npt           ! size of the domain
   INTEGER(KIND=4)                           :: narg, iargc        ! browse command line
   INTEGER(KIND=4)                           :: ncout, ierr        ! browse command line
   INTEGER(KIND=4), DIMENSION(jp_nvar)       :: ipk, id_varout     ! output variable properties

   REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                ! time counter
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1f, e2f           ! F-grid metrics
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1u, e2u           ! zonal horizontal metrics
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1v, e2v           ! meridional horizontal metrics
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: un, vn             ! velocity field
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zun, zvn           ! working arrays
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: fmask              ! fmask

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_rotn, dl_rotn2  ! curl and square curl
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dxrotn, dyrotn     ! curl gradient components
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dxrotn2, dyrotn2   ! square curl gradient components
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_vozocrtx2       ! square of velocity components
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_vomecrty2       ! square of velocity components


   CHARACTER(LEN=256)                        :: cf_ufil, cf_vfil   ! file names
   CHARACTER(LEN=256)                        :: cf_out = 'lambda_int.nc' ! output file name
   CHARACTER(LEN=256)                        :: cv_u, cv_v         ! variable names
   CHARACTER(LEN=256)                        :: cldum              ! dummy string

   TYPE (variable), DIMENSION(jp_nvar)       :: stypvar            ! structure for attibutes

   LOGICAL                                   :: lforcing = .FALSE. ! forcing flag
   LOGICAL                                   :: lchk     = .FALSE. ! flag for missing files
   LOGICAL                                   :: lperio   = .FALSE. ! flag for E-W periodicity
   !!----------------------------------------------------------------------
   CALL ReadCdfNames()

   narg = iargc()
   IF ( narg /= 5 ) THEN
      PRINT *,' usage : cdfeddyscale_pass1 U-file V-file U-var V-var lev'
      PRINT *,'      '
      PRINT *,'     PURPOSE :'
      PRINT *,'     Compute: - the curl and the square of curl on F-points,' 
      PRINT *,'              - the gradient components of the curl and the'
      PRINT *,'                square of the gradient components on UV-points,'
      PRINT *,'              - the square of velocity components on UV-points,' 
      PRINT *,'     for given gridU gridV files and variables. These variables are required'
      PRINT *,'     for computing eddy scales with cdfeddyscale. Therefore this program is'
      PRINT *,'     the first step in computing the eddy scales.'
      PRINT *,'     '
      PRINT *,'        These terms will used to compute the Taylor scale or large'
      PRINT *,'     scale eddy (lambda1) and the small scale eddy (lambda2) in'
      PRINT *,'     the program cdfeddyscale.'
      PRINT *,'      '
      PRINT *,'     ARGUMENTS :'
      PRINT *,'       U-file : zonal component of the vector field.'
      PRINT *,'       V-file : meridional component of the vector field.'
      PRINT *,'       U-var  : zonal component variable name'
      PRINT *,'       V-var  : meridional component variable name.'
      PRINT *,'       lev    : level to be processed. If set to 0, assume forcing file '
      PRINT *,'                in input.'
      PRINT *,'      '
      PRINT *,'     REQUIRED FILES :'
      PRINT *,'        ', TRIM(cn_fhgr)
      PRINT *,'      '
      PRINT *,'     OUTPUT : '
      PRINT *,'       netcdf file : ', TRIM(cf_out)
      PRINT *,'         variables : socurl (s^-1), socurl2 (s^-2)'
      PRINT *,'         variables : sodxcurl, sodycurl (s^-1.m^-1)'
      PRINT *,'         variables : sodxcurl2, sodycurl2 (s^-2.m^-2)'
      PRINT *,'         variables : vozocrtx2, vomecrty2 (m^2.s^-2)'
      PRINT *,'         WARNING : variables in the output file are not located at the same'
      PRINT *,'                 C-grid point.'
      PRINT *,'      '
      PRINT *,'     SEE ALSO : '
      PRINT *,'        cdfeddyscale'
      STOP
   ENDIF


   CALL getarg(1, cf_ufil)
   CALL getarg(2, cf_vfil)
   CALL getarg(3, cv_u   )
   CALL getarg(4, cv_v   )
   CALL getarg(5, cldum  ) ;  READ(cldum,*) ilev

   lchk = chkfile(cn_fhgr ) .OR. lchk
   lchk = chkfile(cf_ufil ) .OR. lchk
   lchk = chkfile(cf_vfil ) .OR. lchk
   IF ( lchk ) STOP ! missing files

   ! load the dimension
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
   IF ( (npk==0) .AND. (ilev > 1) ) THEN
      PRINT *, 'Problem : npk = 0 and lev > 0 STOP'
      STOP
   ELSE
      npk = 1
   END IF

   ! if forcing field 
   IF ( ilev==0 .AND. npk==0 ) THEN
      lforcing=.TRUE.
      npk = 1 ; ilev=1
      PRINT *, 'npk =0, assume 1'
   END IF

   IF ( npt==0 ) THEN
      PRINT *, 'npt=0, assume 1'
      npt=1
   END IF
   ! check files and determines if the curl will be 2D of 3D

   ! Allocate the memory
   ALLOCATE ( e1u(npiglo,npjglo)       , e2u(npiglo,npjglo)       )  
   ALLOCATE ( e1v(npiglo,npjglo)       , e2v(npiglo,npjglo)       )
   ALLOCATE ( e1f(npiglo,npjglo)       , e2f(npiglo,npjglo)       )
   ALLOCATE ( un(npiglo,npjglo)        , vn(npiglo,npjglo)        )
   ALLOCATE ( zun(npiglo,npjglo)       , zvn(npiglo,npjglo)       )
   ALLOCATE ( fmask(npiglo,npjglo)                                )
   ALLOCATE ( tim(npt)                                            )
   ALLOCATE ( dl_rotn(npiglo,npjglo)   , dl_rotn2(npiglo,npjglo)  )
   ALLOCATE ( dxrotn(npiglo,npjglo)    , dyrotn(npiglo,npjglo)    )
   ALLOCATE ( dxrotn2(npiglo,npjglo)   , dyrotn2(npiglo,npjglo)   )
   ALLOCATE ( dl_vozocrtx2(npiglo,npjglo)                         )
   ALLOCATE ( dl_vomecrty2(npiglo,npjglo)                         )


   ! Read the metrics from the mesh_hgr file
   e1u =  getvar(cn_fhgr, cn_ve1u, 1, npiglo, npjglo)
   e2u =  getvar(cn_fhgr, cn_ve2u, 1, npiglo, npjglo)
   e1v =  getvar(cn_fhgr, cn_ve1v, 1, npiglo, npjglo)
   e2v =  getvar(cn_fhgr, cn_ve2v, 1, npiglo, npjglo)
   e1f =  getvar(cn_fhgr, cn_ve1f, 1, npiglo, npjglo)
   e2f =  getvar(cn_fhgr, cn_ve2f, 1, npiglo, npjglo)

   CALL CreateOutputFile 

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
      ! compute the curl
      dl_rotn(:,:) = 0.d0
      DO jj = 1, npjglo -1
         DO ji = 1, npiglo -1   ! vector opt.
            dl_rotn(ji,jj) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ) - e2v(ji,jj) *vn(ji,jj) &
                 &         - e1u(ji  ,jj+1) * un(ji  ,jj+1) + e1u(ji,jj) *un(ji,jj)  ) &
                 &         * fmask(ji,jj) / ( e1f(ji,jj) * e2f(ji,jj) )
         END DO
      END DO

      IF ( lperio ) dl_rotn(npiglo,:) = dl_rotn(2, :)

      ! compute the square curl
      dl_rotn2(:,:) = dl_rotn(:,:) * dl_rotn(:,:)

      ! compute the gradient components
      dxrotn(:,:) = 0.d0
      dyrotn(:,:) = 0.d0
      DO jj = 2, npjglo 
         DO ji = 2, npiglo    ! vector opt.
            dxrotn(ji,jj) = (dl_rotn(ji,jj) - dl_rotn(ji-1,jj))/e1v(ji,jj)
            dyrotn(ji,jj) = (dl_rotn(ji,jj) - dl_rotn(ji,jj-1))/e2u(ji,jj)
         END DO
      END DO

      IF ( lperio ) dxrotn(1,:) = dxrotn(npiglo-1, :)
      IF ( lperio ) dyrotn(1,:) = dyrotn(npiglo-1, :)

      ! compute the square module of the gradient
      dxrotn2(:,:) = dxrotn(:,:) * dxrotn(:,:)
      dyrotn2(:,:) = dyrotn(:,:) * dyrotn(:,:)

      ! compute the square of the velocity components
      dl_vozocrtx2(:,:) = zun(:,:) * zun(:,:) 
      dl_vomecrty2(:,:) = zvn(:,:) * zvn(:,:)

      ! write dl_rotn on file at level k and at time jt
      ierr = putvar(ncout, id_varout(jp_curl),    REAL(dl_rotn),  1, npiglo, npjglo, ktime=jt)

      ! write dl_rotn2 on file at level k and at time jt
      ierr = putvar(ncout, id_varout(jp_curl2),   REAL(dl_rotn2), 1, npiglo, npjglo, ktime=jt)

      ! write dxrotn on file at level k and at time jt
      ierr = putvar(ncout, id_varout(jp_dxcurl),  REAL(dxrotn),   1, npiglo, npjglo, ktime=jt)

      ! write dyrotn on file at level k and at time jt
      ierr = putvar(ncout, id_varout(jp_dycurl),  REAL(dyrotn),   1, npiglo, npjglo, ktime=jt)

      ! write dxrotn2 on file at level k and at time jt
      ierr = putvar(ncout, id_varout(jp_dxcurl2), REAL(dxrotn2),  1, npiglo, npjglo, ktime=jt)

      ! write dyrotn2 on file at level k and at time jt
      ierr = putvar(ncout, id_varout(jp_dycurl2), REAL(dyrotn2),  1, npiglo, npjglo, ktime=jt)

      ! write dl_vozocrtx2 on file at level k and at time jt
      ierr = putvar(ncout, id_varout(jp_u2),  REAL(dl_vozocrtx2), 1, npiglo, npjglo, ktime=jt)

      ! write dl_vozocrtx2 on file at level k and at time jt
      ierr = putvar(ncout, id_varout(jp_v2),  REAL(dl_vomecrty2), 1, npiglo, npjglo, ktime=jt)

   END DO
   ierr = closeout(ncout)
CONTAINS

   SUBROUTINE CreateOutputFile
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE CreateOutputFile  ***
      !!
      !! ** Purpose :  Create output file 
      !!
      !! ** Method  :   Use global program variables 
      !!
      !!----------------------------------------------------------------------

      ! define new variables for output
      ! Relative Vorticity F point
      ipk(jp_curl)                       = 1   !2D
      stypvar(jp_curl)%cname             = 'socurl'
      stypvar(jp_curl)%cunits            = 's-1'
      stypvar(jp_curl)%rmissing_value    = 0.
      stypvar(jp_curl)%valid_min         = -1000.
      stypvar(jp_curl)%valid_max         =  1000.
      stypvar(jp_curl)%clong_name        = 'Relative_Vorticity (curl)'
      stypvar(jp_curl)%cshort_name       = 'socurl'
      stypvar(jp_curl)%conline_operation = 'N/A'
      stypvar(jp_curl)%caxis             = 'TYX'

      ! Square of Relative Vorticity F point
      ipk(jp_curl2)                       = 1   !2D
      stypvar(jp_curl2)%cname             = 'socurl2'
      stypvar(jp_curl2)%cunits            = 's-2'
      stypvar(jp_curl2)%rmissing_value    = 0.
      stypvar(jp_curl2)%valid_min         = -1000.
      stypvar(jp_curl2)%valid_max         =  1000.
      stypvar(jp_curl2)%clong_name        = 'Square of Relative_Vorticity (curl2)'
      stypvar(jp_curl2)%cshort_name       = 'socurl2'
      stypvar(jp_curl2)%conline_operation = 'N/A'
      stypvar(jp_curl2)%caxis             = 'TYX'


      ! Relative Vorticity zonal gradient V point
      ipk(jp_dxcurl)                       = 1
      stypvar(jp_dxcurl)%cname             = 'sodxcurl'
      stypvar(jp_dxcurl)%cunits            = 'm-1*s-1'
      stypvar(jp_dxcurl)%rmissing_value    = 0.
      stypvar(jp_dxcurl)%valid_min         = -1000.
      stypvar(jp_dxcurl)%valid_max         =  1000.
      stypvar(jp_dxcurl)%clong_name        = 'Relative_Vorticity zonal gradient (dx_curl)'
      stypvar(jp_dxcurl)%cshort_name       = 'sodxcurl'
      stypvar(jp_dxcurl)%conline_operation = 'N/A'
      stypvar(jp_dxcurl)%caxis             = 'TYX'

      ! Relative Vorticity meridional gradient U point
      ipk(jp_dycurl)                       = 1
      stypvar(jp_dycurl)%cname             = 'sodycurl'
      stypvar(jp_dycurl)%cunits            = 'm-1*s-1'
      stypvar(jp_dycurl)%rmissing_value    = 0.
      stypvar(jp_dycurl)%valid_min         = -1000.
      stypvar(jp_dycurl)%valid_max         =  1000.
      stypvar(jp_dycurl)%clong_name        = 'Relative Vorticity meridional gradient (dy_curl)'
      stypvar(jp_dycurl)%cshort_name       = 'sodycurl'
      stypvar(jp_dycurl)%conline_operation = 'N/A'
      stypvar(jp_dycurl)%caxis             = 'TYX'

      ! Square of Relative Vorticity zonal gradient V point
      ipk(jp_dxcurl2)                       = 1
      stypvar(jp_dxcurl2)%cname             = 'sodxcurl2'
      stypvar(jp_dxcurl2)%cunits            = 'm-2*s-2'
      stypvar(jp_dxcurl2)%rmissing_value    = 0.
      stypvar(jp_dxcurl2)%valid_min         = -1000.
      stypvar(jp_dxcurl2)%valid_max         =  1000.
      stypvar(jp_dxcurl2)%clong_name        = 'Square Relative Vorticity zonal gradient (dx_curl2)'
      stypvar(jp_dxcurl2)%cshort_name       = 'sodxcurl2'
      stypvar(jp_dxcurl2)%conline_operation = 'N/A'
      stypvar(jp_dxcurl2)%caxis             = 'TYX'

      ! Square of Relative Vorticity meridional gradient U point
      ipk(jp_dycurl2)                       = 1
      stypvar(jp_dycurl2)%cname             = 'sodycurl2'
      stypvar(jp_dycurl2)%cunits            = 'm-2*s-2'
      stypvar(jp_dycurl2)%rmissing_value    = 0.
      stypvar(jp_dycurl2)%valid_min         = -1000.
      stypvar(jp_dycurl2)%valid_max         =  1000.
      stypvar(jp_dycurl2)%clong_name        = 'Square of Relative Vorticity meridional gradient (dy_curl2)'
      stypvar(jp_dycurl2)%cshort_name       = 'sodycurl2'
      stypvar(jp_dycurl2)%conline_operation = 'N/A'
      stypvar(jp_dycurl2)%caxis             = 'TYX'

      ! Square of Zonal Velocity V point
      ipk(jp_u2)                       = 1
      stypvar(jp_u2)%cname             = 'vozocrtx2'
      stypvar(jp_u2)%cunits            = 'm^2/s^2'
      stypvar(jp_u2)%rmissing_value    = 0.
      stypvar(jp_u2)%valid_min         = -10.
      stypvar(jp_u2)%valid_max         =  10.
      stypvar(jp_u2)%clong_name        = 'Square Zonal Velocity (vozocrtx2)'
      stypvar(jp_u2)%cshort_name       = 'vozocrtx2'
      stypvar(jp_u2)%conline_operation = 'N/A'
      stypvar(jp_u2)%caxis             = 'TYX'

      ! Square of Zonal Velocity V point
      ipk(jp_v2)                       = 1
      stypvar(jp_v2)%cname             = 'vomecrty2'
      stypvar(jp_v2)%cunits            = 'm^2/s^2'
      stypvar(jp_v2)%rmissing_value    = 0.
      stypvar(jp_v2)%valid_min         = -10.
      stypvar(jp_v2)%valid_max         =  10.
      stypvar(jp_v2)%clong_name        = 'Square Meridional Velocity (vomecrty2)'
      stypvar(jp_v2)%cshort_name       = 'vomecrty2'
      stypvar(jp_v2)%conline_operation = 'N/A'
      stypvar(jp_v2)%caxis             = 'TYX'

      ! use zun and zvn to store f latitude and longitude for output
      zun = getvar(cn_fhgr, cn_glamf, 1, npiglo, npjglo)
      zvn = getvar(cn_fhgr, cn_gphif, 1, npiglo, npjglo)

      ! look for  E-W periodicity
      IF ( zun(1,1) == zun(npiglo-1,1) ) lperio = .TRUE.

      ! create output fileset
      ncout = create      (cf_out, cf_ufil, npiglo,  npjglo, 0)
      ierr  = createvar   (ncout , stypvar, jp_nvar, ipk,    id_varout)
      ierr  = putheadervar(ncout,  cf_ufil, npiglo,  npjglo, 0, pnavlon=zun, pnavlat=zvn )

      tim  = getvar1d(cf_ufil, cn_vtimec, npt      )
      ierr = putvar1d(ncout,   tim,       npt,  'T')

   END SUBROUTINE CreateOutputFile
END PROGRAM cdfeddyscale_pass1
