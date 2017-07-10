PROGRAM cdfsigintegr_pedro
  !!======================================================================
  !!                     ***  PROGRAM  cdfsigintegr  ***
  !!=====================================================================
  !!  ** Purpose : This program is used to integrate quantities between 
  !!               isopycnals
  !!
  !!  ** Method  : Linear interpolation is used on the vertical to define
  !!               the depth of the given isopycn.
  !!               Then, the integral is performed from the top of the ocean
  !!               down to the given isopycnal. Finaly, by making the 
  !!               difference between 2 isopycnals we obtain the required 
  !!               quantity.
  !!
  !! History : 2.1  : 12/2007  : J.M. Molines : Original code
  !!           3.0  : 06/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id: cdfsigintegr.f90 657 2013-05-21 16:36:14Z molines $
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class integration
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: ji, jj, jk, jt   ! dummy loop index
  INTEGER(KIND=4)                               :: jiso, jfich      ! dummy loop index
  INTEGER(KIND=4)                               :: jvar             ! dummy loop index
  INTEGER(KIND=4)                               :: npiglo, npjglo   ! domain size
  INTEGER(KIND=4)                               :: npk, npt         ! domain size
  INTEGER(KIND=4)                               :: npiso, nvars     ! number of isopycnals, variables
  INTEGER(KIND=4)                               :: narg, iargc      ! command line
  INTEGER(KIND=4)                               :: ijarg, ireq      ! command line
  INTEGER(KIND=4)                               :: nfiles           ! number of input files
  INTEGER(KIND=4)                               :: istrt_arg        ! argument number of first input file
  INTEGER(KIND=4)                               :: ik0              ! layer index
  INTEGER(KIND=4)                               :: ijk              ! layer index
  INTEGER(KIND=4)                               :: numin=10         ! logical unit for ascii input file
  INTEGER(KIND=4)                               :: ncout, ierr      ! ncid and status variable
  INTEGER(KIND=4), DIMENSION(9)                 :: ipk, id_varout   ! levels and id's of output variables
  !
  REAL(KIND=4), DIMENSION(:,:,:),   ALLOCATABLE :: v3d              ! 3D v working array (npk)
  REAL(KIND=4), DIMENSION(:,:,:),   ALLOCATABLE :: u3d              ! 3D u working array (npk)
  REAL(KIND=4), DIMENSION(:,:,:),   ALLOCATABLE :: m3d              ! 3D modulus working array (npk)
  REAL(KIND=4), DIMENSION(:,:,:),   ALLOCATABLE :: zint             ! pseudo 3D working array (2)
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d              ! 2D v working array
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: u2d              ! 2D u working array
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: e3               ! vertical metrics
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: e1t              ! U/V Horizontal metric !!Pedro
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: e2t              ! U/V Horizontal metric !!Pedro
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: tmask            ! mask of t points from rho
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zdum             ! dummy array for I/O
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zdum_u           ! dummy array for I/O for u
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zdum_v           ! dummy array for I/O for v
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zdum_m           ! dummy array for I/O for m
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zdum_a           ! dummy array for I/O for angle
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: rho_lev          ! value of isopycnals
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: tim              ! time counter
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: h1d              ! depth of rho points
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: gdepw            ! depth of W points
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: e31d             ! vertical metrics in full step
  REAL(KIND=4)                                  :: zspval=999999.   ! output missing value
  REAL(KIND=4)                                  :: zspvalz          ! missing value from rho file      
  REAL(KIND=4)                                  :: maskvalue=0.     ! masked values for ocean velocity output
  INTEGER(KIND=4)                               :: nboutput=9       ! number of values to write in cdf output
  REAL(KIND=4)                                  :: pi               !pi 

  REAL(KIND=8), DIMENSION(:,:,:),   ALLOCATABLE :: dv2dint          ! interpolated value / int of v 
  REAL(KIND=8), DIMENSION(:,:,:),   ALLOCATABLE :: du2dint          ! int of u
  REAL(KIND=8), DIMENSION(:,:,:),   ALLOCATABLE :: dm2dint          ! int of the modulus
  REAL(KIND=8), DIMENSION(:,:,:),   ALLOCATABLE :: dalpha           ! 3D coefficient (npiso)

  CHARACTER(LEN=256)                            :: cf_rholev = 'rho_lev' ! input file for rho surfaces
  CHARACTER(LEN=256)                            :: cf_ufil          ! input file for u data
  CHARACTER(LEN=256)                            :: cf_vfil          ! input file for v data
  CHARACTER(LEN=256)                            :: cf_rho           ! input file for density
  CHARACTER(LEN=256)                            :: cf_out           ! output file
  CHARACTER(LEN=256)                            :: cv_in            ! name of input variable
  CHARACTER(LEN=256)                            :: cldum            ! dummy string variable
  CHARACTER(LEN=256)                            :: cluni            ! dummy string variable for variable units
  CHARACTER(LEN=256)                            :: cglobal          ! global attribute
  CHARACTER(LEN=256)                            :: ctype='T'        ! position of variable on C grid
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names         ! temporary arry for variable name in file

  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvar          ! structure for attributes
  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypzvar         ! structure for attributes

  LOGICAL                                       :: lfull = .FALSE.  ! flag for full step
  LOGICAL                                       :: lchk  = .FALSE.  ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfsigintegr IN-var RHO-file list_of_files [ VAR-type ] ...'
     PRINT *,'              ... [ -sig sigma_name] [ -full ] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :' 
     PRINT *,'       Take a list of input files with specific IN-var variable, associated'
     PRINT *,'       with a reference density file. A set of isopycnal surfaces is defined'
     PRINT *,'       in an ASCII file (rho_lev by default), using same depth reference than'
     PRINT *,'       the input reference density file. This program computes the integral of'
     PRINT *,'       IN-var between the isopycnals defined in rho_lev. It also gives the '
     PRINT *,'       isopycnal depth and thickness of density layers.'
     PRINT *,'      '
     PRINT *,'       Rho_lev file first line indicates the number of following isopycnals.'
     PRINT *,'       Then a list of the densities is given, one per line.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       IN-var : input variable to be integrated' 
     PRINT *,'       RHO-file : netcdf file with already computed density' 
     PRINT *,'       list_of_files : a list of model netcdf files containing IN-var.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [ VAR-type ] : one of T U V F W which defined the position on' 
     PRINT *,'               IN-var in the model C-grid. Default is ', TRIM(ctype)
     PRINT *,'       [ -sig sigma_name ] : give the name of sigma variable in RHO-file.'
     PRINT *,'               Default is ',TRIM(cn_vosigma0)
     PRINT *,'       [ -full ] : indicate a full step configuration.'
     PRINT *,'       [ -rholev  file] : indicates name of file defining the limits for '
     PRINT *,'               integration. Default is ', TRIM(cf_rholev)
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fzgr),',',TRIM(cn_fhgr),' and ',TRIM(cf_rholev)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : IN-file.integr'
     PRINT *,'         variables : inv_IN-var  : inventory of IN-var from input file.'
     PRINT *,'                     ', TRIM(cn_vodepiso),' (m) : depth of isopycnal.'
     PRINT *,'                     ', TRIM(cn_isothick),' (m) : thickness of isopycnal layer.'
     PRINT *,'                     mean_IN-var (same unit as IN-var) : mean IN-var in the isopycnal'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      cdfrhoproj, cdfsigtrp, cdfisopycdep'
     PRINT *,'      '
     STOP
  ENDIF
  pi=ACOS(-1.)
  !ijarg = 1 ; ireq = 0 ; nfiles = 0
  ijarg = 1 ; ireq = 0
  DO WHILE ( ijarg <= narg ) 
     CALL getarg( ijarg, cldum ) ; ijarg = ijarg+1
     SELECT CASE ( cldum )
     !CASE ( 'T','t','U','u','V','v','F','f','W','w' )
        !ctype=cldum
     CASE ( '-sig '   ) ; CALL getarg( ijarg, cn_vosigma0) ; ijarg = ijarg+1
     CASE ( '-rholev ') ; CALL getarg( ijarg, cf_rholev  ) ; ijarg = ijarg+1
     CASE ( '-full '  ) ; lfull = .TRUE.
     CASE DEFAULT
        ireq=ireq+1
        SELECT CASE ( ireq )
        CASE ( 1 ) ; cf_rho = cldum
        CASE ( 2 ) ; cf_ufil  = cldum
        CASE ( 3 ) ; cf_vfil = cldum
        CASE DEFAULT 
           !nfiles=nfiles+1
           !IF ( nfiles == 1 ) istrt_arg = ijarg - 1
           PRINT *,' Too many arguments ' ; STOP
        END SELECT
     END SELECT
  END DO

  CALL SetGlobalAtt( cglobal )

  ! check for files existence
  lchk = lchk .OR. chkfile (cn_fzgr   )
  lchk = lchk .OR. chkfile (cn_fhgr   ) !!Modified Pedro (Added)
  lchk = lchk .OR. chkfile (cf_rholev )
  lchk = lchk .OR. chkfile (cf_rho    )
  lchk = lchk .OR. chkfile (cf_ufil   )
  lchk = lchk .OR. chkfile (cf_vfil   )
  IF ( lchk ) STOP ! missing file

  ! Read rho level between which the integral is being performed
  OPEN(numin,file=cf_rholev)
  READ(numin,*) npiso
  ALLOCATE (rho_lev(npiso) )
  PRINT *,' Density limits read in ',TRIM(cf_rholev)
  DO jiso=1,npiso
     READ(numin,*) rho_lev(jiso)
     PRINT *,rho_lev(jiso)
  END DO
  CLOSE(numin)

  npiglo = getdim(cf_rho, cn_x)
  npjglo = getdim(cf_rho, cn_y)
  npk    = getdim(cf_rho, cn_z)
  
  zspvalz=getspval(cf_rho, cn_vosigma0)

  !CALL getarg(istrt_arg, cf_ufil)
  !IF ( chkfile ( cf_ufil ) ) STOP ! missing file

  nvars=getnvar(cf_vfil)
  ALLOCATE(cv_names(nvars), stypzvar(nvars))

  cv_names(:)=getvarname(cf_vfil,nvars,stypzvar)
  ALLOCATE( e1t(npiglo,npjglo), e2t(npiglo,npjglo) ) !!Modified Pedro (added)
  ALLOCATE( v3d(npiglo,npjglo,npk), u3d(npiglo,npjglo,npk), m3d(npiglo,npjglo,npk), u2d(npiglo,npjglo))
  ALLOCATE(dalpha(npiglo,npjglo,npiso), e3(npiglo,npjglo) )
  ALLOCATE( dv2dint(npiglo,npjglo,2), v2d(npiglo,npjglo), zint(npiglo,npjglo,2)  )
  ALLOCATE( du2dint(npiglo,npjglo,2), dm2dint(npiglo,npjglo,2))
  ALLOCATE( h1d(npk) ,gdepw(npk) ,tmask(npiglo,npjglo), zdum_u(npiglo,npjglo) )
  ALLOCATE( zdum_v(npiglo,npjglo), zdum_m(npiglo,npjglo), zdum(npiglo,npjglo) )
  ALLOCATE( stypvar(nboutput) , zdum_a(npiglo,npjglo) )
  IF ( lfull ) ALLOCATE ( e31d(npk) )

  gdepw(:) = getvare3(cn_fzgr, cn_gdepw, npk)

  IF (lfull ) e31d(:) = getvare3(cn_fzgr, cn_ve3t1d, npk)

  h1d(:) = getvar1d(cf_rho, cn_vdeptht, npk)

  ! Note, if working with vertical slabs, one may avoid 3D array, but may be slow ...
  tmask=1.
  DO jk=1,npk
     v3d(:,:,jk) = getvar(cf_rho, cn_vosigma0, jk, npiglo, npjglo)
     IF ( jk == 1 ) THEN
        WHERE (v3d(:,:,jk) == zspvalz ) tmask=0.
     ENDIF
  END DO

  !! Modified Pedro
  !! Get horizontal grid
  !! Obtain dx on t point
  !IF (ctype == 'U') THEN
  e2t(:,:)   = getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo)
  !! Obtain dy on t point
  !ELSE IF (ctype == 'V') THEN
  e1t(:,:)   = getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo)
  !! end modifying
  !END IF

  !! ** Compute interpolation coefficients as well as the level used
  !!    to interpolate between
  !! Comented Pedro
  DO ji=1,npiglo
     DO jj = 1, npjglo
        ijk = 1
        !!DO jiso=1,npiso 
        jiso=1 !! La primera definicion de isopyc la dejamos como esta
           !  Assume that rho (z) is increasing downward (no inversion)
           !     Caution with sigma0 at great depth !
           DO WHILE (rho_lev(jiso) >=  v3d(ji,jj,ijk) .AND. ijk <= npk &
                &                .AND. v3d(ji,jj,ijk) /=  zspvalz )
              ijk = ijk+1
           END DO
           ijk = ijk-1
           ik0 = ijk
           IF (ijk == 0) THEN
              ijk = 1
              dalpha(ji,jj,jiso) = 0.d0
           ELSE IF (v3d(ji,jj,ijk+1) == zspvalz ) THEN !!Masked values
              ik0 = 0
              dalpha(ji,jj,jiso) = 0.d0
           ELSE 
              ! ... dalpha is always in [0,1]. Adding ik0 ( >=1 ) for saving space for ik0
              dalpha(ji,jj,jiso)= (rho_lev(jiso)-v3d(ji,jj,ijk))/(v3d(ji,jj,ijk+1)-v3d(ji,jj,ijk)) + ik0
           ENDIF
        !!END DO
     END DO
  END DO

      !! Commented Pedro
   DO ji=1,npiglo
      DO jj = 1, npjglo
         ijk = 1
         !!DO jiso=1,npiso
         jiso=2  !La segunda definicion de isopyc la modificamos para extraer la superficie inferior
            !  Assume that rho (z) is increasing downward (no inversion)
            !     Caution with sigma0 at great depth ! 
            maskvalue=v3d(ji,jj,npk) !We pick a value that we know is masked (last value)
            DO WHILE (maskvalue /=  v3d(ji,jj,ijk) .AND. ijk <= npk ) !While the value is not masked, +1
               ijk = ijk+1
            END DO
            ijk = ijk-1
            dalpha(ji,jj,jiso) = ijk !Equal the first value before the mask appears
         !!END DO
      END DO
  END DO

  ! define header of all files
  ipk(1)=npiso-1 ; ipk(2)=npiso-1 ; ipk(3)=npiso ; ipk(4)=npiso-1
  ipk(5)=npiso-1 ; ipk(6)=npiso-1 ; ipk(7)=npiso-1 ; ipk(8)=npiso-1 
  ipk(9)=npiso-1
  !DO jvar=1,nvars
  !   IF ( cv_in == stypzvar(jvar)%cname ) THEN 
  !      stypvar(1)=stypzvar(jvar)
  !      EXIT
  !   ENDIF
  !END DO
  ! save original long name for further process
  !cldum = TRIM(stypvar(1)%clong_name)
  !cluni = TRIM(cv_vomecrty%cunits)

  stypvar(1)%cname             = 'int_vomecrty'
  stypvar(1)%clong_name        = 'vomecrty integrated on sigma bin'
  stypvar(1)%cshort_name       = stypvar(1)%cname
  stypvar(1)%cunits            = 'Sv'
  stypvar(1)%rmissing_value    = zspval
  stypvar(1)%caxis             = 'TRYX'

  stypvar(2)%cname             = TRIM(cn_isothick)
  stypvar(2)%cunits            = 'm'
  stypvar(2)%rmissing_value    = zspval
  stypvar(2)%valid_min         = 0.
  stypvar(2)%valid_max         = 7000.
  stypvar(2)%clong_name        = 'Thickness_of_Isopycnals'
  stypvar(2)%cshort_name       = TRIM(cn_isothick)
  stypvar(2)%conline_operation = 'N/A'
  stypvar(2)%caxis             = 'TRYX'

  stypvar(3)%cname             = TRIM(cn_vodepiso)
  stypvar(3)%cunits            = 'm'
  stypvar(3)%rmissing_value    = zspval
  stypvar(3)%valid_min         = 0.
  stypvar(3)%valid_max         = 7000.
  stypvar(3)%clong_name        = 'Depth_of_Isopycnals'
  stypvar(3)%cshort_name       = TRIM(cn_vodepiso)
  stypvar(3)%conline_operation = 'N/A'
  stypvar(3)%caxis             = 'TRYX'

  stypvar(4)%cname             = 'mean_vomecrty'
  stypvar(4)%cunits            = 'm/s'
  stypvar(4)%rmissing_value    = zspval
  stypvar(4)%valid_min         = stypvar(1)%valid_min
  stypvar(4)%valid_max         = stypvar(1)%valid_min
  stypvar(4)%clong_name        = ' mean value in sigma layer'
  stypvar(4)%cshort_name       = stypvar(4)%cname
  stypvar(4)%conline_operation = 'N/A'
  stypvar(4)%caxis             = 'TRYX'

  stypvar(5)%cname             = 'int_vozocrtx'
  stypvar(5)%clong_name        = ' vozocrtx integrated on sigma bin'
  stypvar(5)%cshort_name       = stypvar(5)%cname
  stypvar(5)%cunits            = 'Sv'
  stypvar(5)%rmissing_value    = zspval
  stypvar(5)%caxis             = 'TRYX'

  stypvar(6)%cname             = 'mean_vozocrtx'
  stypvar(6)%cunits            = 'm/s'
  stypvar(6)%rmissing_value    = zspval
  stypvar(6)%valid_min         = stypvar(5)%valid_min
  stypvar(6)%valid_max         = stypvar(5)%valid_min
  stypvar(6)%clong_name        = 'mean value in sigma layer'
  stypvar(6)%cshort_name       = stypvar(6)%cname
  stypvar(6)%conline_operation = 'N/A'
  stypvar(6)%caxis             = 'TRYX'
  
  stypvar(7)%cname             = 'int_modulus'
  stypvar(7)%clong_name        = 'integration of the modulus'
  stypvar(7)%cshort_name       = stypvar(7)%cname
  stypvar(7)%cunits            = 'Sv'
  stypvar(7)%rmissing_value    = zspval
  stypvar(7)%caxis             = 'TRYX'

  stypvar(8)%cname             = 'modulus_mean' 
  stypvar(8)%cunits            = 'm/s'
  stypvar(8)%rmissing_value    = zspval
  stypvar(8)%valid_min         = stypvar(7)%valid_min
  stypvar(8)%valid_max         = stypvar(7)%valid_min
  stypvar(8)%clong_name        = 'mean value of the modulus in sigma layer'
  stypvar(8)%cshort_name       = stypvar(8)%cname
  stypvar(8)%conline_operation = 'N/A'
  stypvar(8)%caxis             = 'TRYX'
  
  stypvar(9)%cname             = 'angle_modulus'
  stypvar(9)%clong_name        = 'Modulus angle'
  stypvar(9)%cshort_name       = stypvar(9)%cname
  stypvar(9)%cunits            = '°'
  stypvar(9)%rmissing_value    = zspval
  stypvar(9)%caxis             = 'TRYX'

!! ** Loop on the scalar files to project on choosen isopycnics surfaces
  !DO jfich=1, nfiles

     !CALL getarg(jfich+istrt_arg-1, cf_ufil)
     !IF ( chkfile (cf_ufil) ) STOP ! missing file
     !IF ( chkfile (cf_vfil) ) STOP ! missing file
     PRINT *,'working with ', TRIM(cf_ufil)
     PRINT *,'working with ', TRIM(cf_vfil)

     ! create output file
     !! Modified Pedro
     !cf_out=TRIM(cf_in)//'.integr'
     cf_out='cdfsigintegr_file'
     ncout = create      (cf_out, cf_rho,  npiglo, npjglo, npiso                       )
     ierr  = createvar   (ncout,  stypvar, 9,      ipk,    id_varout, cdglobal=cglobal )
     ierr  = putheadervar(ncout,  cf_rho,  npiglo, npjglo, npiso, pdep=rho_lev         )
     ! copy time arrays in output file
     npt = getdim ( cf_ufil, cn_t)
     ALLOCATE ( tim(npt) )
     tim(:) = getvar1d(cf_ufil, cn_vtimec, npt     )
     ierr   = putvar1d(ncout, tim,       npt, 'T')
     DEALLOCATE ( tim )
     DO jt =1, npt
        DO jk=1,npk
           v2d(:,:) = getvar(cf_vfil, cn_vomecrty, jk, npiglo, npjglo, ktime = jt )
           u2d(:,:) = getvar(cf_ufil, cn_vozocrtx, jk, npiglo, npjglo, ktime = jt )
           !SELECT CASE ( ctype )
           !CASE ('T', 't' )
           !   v3d(:,:,jk) = v2d(:,:)
           !CASE ('U','u' )
              DO jj=1,npjglo
                 DO ji=2, npiglo
                    u3d(ji,jj,jk)=0.5*( u2d(ji,jj) + u2d(ji-1,jj) )  ! put variable on T point
                    u3d(1,jj,jk)=0.5*( u2d(2,jj) + u2d(1,jj) )  ! put variable on T point
                 END DO
              END DO
           !CASE ('V','v' )
              DO jj=2,npjglo
                 DO ji=1, npiglo
                    v3d(ji,jj,jk)=0.5*( v2d(ji,jj) + v2d(ji,jj-1) )  ! put variable on T point
                    v3d(ji,1,jk)=0.5*( v2d(ji,2) + v2d(ji,1) )  ! put variable on T point
                 END DO
              END DO              
           !CASE('W','w' )
              !v3d(:,:,jk) = v2d(:,:)
              !v2d(:,:) = getvar(cf_ufil, cv_in, jk+1, npiglo, npjglo, ktime = jt )
              !v3d(:,:,jk) = 0.5 * ( v3d(:,:,jk) + v2d(:,:) )
           !CASE('F','f' )
              !DO jj = 2, npjglo
                 !DO ji = 2, npiglo
                    !v3d(:,:,jk) = 0.25*( v2d(ji,jj) + v2d( ji, jj-1) + v2d (ji-1,jj-1) + v2d(ji-1, jj) )
                 !END DO
              !END DO
           !END SELECT
        END DO

        ! Compute integral from surface to isopycnal
        !!DO jiso=1,npiso
        jiso=1 !! La primera integral la dejamos como esta
           ! determine isopycnal surface
           DO ji=1,npiglo
              DO jj=1,npjglo
                 ! ik0 is retrieved from dalpha, taking the integer part.
                 ik0=INT(dalpha(ji,jj,jiso)) ; dalpha(ji,jj,jiso) =  dalpha(ji,jj,jiso) - ik0
                 IF (ik0 /= 0) THEN
                    zint (ji,jj,1)=dalpha(ji,jj,jiso)*h1d(ik0+1) + (1.d0-dalpha(ji,jj,jiso))*h1d(ik0)
                 ELSE 
                    zint  (ji,jj,1)=0. !zspval
                 ENDIF
              END DO
           END DO

        ! Compute integral from surface to isopycnal
        jiso=2 !! En la segunda utilizamos la batimetria
           ! determine isopycnal surface
           DO ji=1,npiglo
              DO jj=1,npjglo
                 ! ik0 is retrieved from dalpha, taking the integer part.
                 ik0=INT(dalpha(ji,jj,jiso)) ; dalpha(ji,jj,jiso) =  dalpha(ji,jj,jiso) - ik0
                 IF (zint(ji,jj,1) /= 0) THEN
                    zint (ji,jj,2)=dalpha(ji,jj,jiso)*h1d(ik0+1) + (1.d0-dalpha(ji,jj,jiso))*h1d(ik0)
                 ELSE
                    zint  (ji,jj,2)=0. !zspval
                 ENDIF
              END DO
           END DO
           ! integrate from jk=1 to zint
           du2dint(:,:,1) = 0.d0
           du2dint(:,:,2) = 0.d0
           dv2dint(:,:,1) = 0.d0
           dv2dint(:,:,2) = 0.d0
           dm2dint(:,:,1) = 0.d0
           dm2dint(:,:,2) = 0.d0

           DO jk=1,npk-1
              ! get metrixs at level jk
              IF ( lfull ) THEN 
                 e3(:,:) = e31d(jk)
              ELSE
                 e3(:,:)=getvar(cn_fzgr,'e3t_ps',jk,npiglo,npjglo,ldiom=.TRUE.)
              ENDIF

              DO ji=1,npiglo
                 DO jj=1,npjglo
                    IF ( gdepw(jk)+e3(ji,jj) < zint(ji,jj,1) ) THEN  ! full cell
                       !!dv2dint(ji,jj,1)=dv2dint(ji,jj,1) + e3(ji,jj)* v3d(ji,jj,jk)
                       du2dint(ji,jj,1)=du2dint(ji,jj,1) + e3(ji,jj)* u3d(ji,jj,jk) * e2t(ji,jj) / 1000000 !!Modified Pedro
                       dv2dint(ji,jj,1)=dv2dint(ji,jj,1) + e3(ji,jj)* v3d(ji,jj,jk) * e1t(ji,jj) / 1000000!!Modified Pedro
                       !dm2dint(ji,jj,1)=dm2dint(ji,jj,1) + SQRT(du2dint(ji,jj,1)*du2dint(ji,jj,1) + &
                       !                                   & dv2dint(ji,jj,1)*dv2dint(ji,jj,1))
                    ELSE IF (( zint(ji,jj,1) <= gdepw(jk)+e3(ji,jj) ) .AND. (zint(ji,jj,1) > gdepw(jk)) ) THEN
                       !!dv2dint(ji,jj,1)=dv2dint(ji,jj,1)+ (zint(ji,jj,1) - gdepw(jk) )* v3d(ji,jj,jk)
                       du2dint(ji,jj,1)=du2dint(ji,jj,1)+ (zint(ji,jj,1) - gdepw(jk) )* u3d(ji,jj,jk) * e2t(ji,jj) / 1000000!! Modified Pedro
                       dv2dint(ji,jj,1)=dv2dint(ji,jj,1)+ (zint(ji,jj,1) - gdepw(jk) )* v3d(ji,jj,jk) * e1t(ji,jj) / 1000000!! Modified Pedro
                       !dm2dint(ji,jj,1)=dm2dint(ji,jj,1)+ (zint(ji,jj,1) - gdepw(jk) )* SQRT(du2dint(ji,jj,1)* &
                       !                 & du2dint(ji,jj,1) + dv2dint(ji,jj,1)*dv2dint(ji,jj,1)) !! Modified Pedro
                    ELSE   ! below the isopycnal 
                       ! do nothing for this i j point
                    ENDIF
                 END DO
              END DO
              
              DO ji=1,npiglo
                 DO jj=1,npjglo
                    IF ( gdepw(jk)+e3(ji,jj) < zint(ji,jj,2) ) THEN  ! full cell
                       !!dv2dint(ji,jj,1)=dv2dint(ji,jj,1) + e3(ji,jj)* v3d(ji,jj,jk)
                       du2dint(ji,jj,2)=du2dint(ji,jj,2) + e3(ji,jj)* u3d(ji,jj,jk) * e2t(ji,jj) / 1000000
                       dv2dint(ji,jj,2)=dv2dint(ji,jj,2) + e3(ji,jj)* v3d(ji,jj,jk) * e1t(ji,jj) / 1000000!!Modified Pedro
                       !dm2dint(ji,jj,2)=dm2dint(ji,jj,2) + SQRT(du2dint(ji,jj,2)*du2dint(ji,jj,2) + &
                       !                                    & dv2dint(ji,jj,2)*dv2dint(ji,jj,2))
                    ELSE IF (( zint(ji,jj,2) <= gdepw(jk)+e3(ji,jj) ) .AND. (zint(ji,jj,2) > gdepw(jk)) ) THEN
                       !!dv2dint(ji,jj,1)=dv2dint(ji,jj,1)+ (zint(ji,jj,1) - gdepw(jk) )* v3d(ji,jj,jk)
                       du2dint(ji,jj,2)=du2dint(ji,jj,2)+ (zint(ji,jj,2) - gdepw(jk) )* u3d(ji,jj,jk) * e2t(ji,jj) / 1000000!! Modified Pedro
                       dv2dint(ji,jj,2)=dv2dint(ji,jj,2)+ (zint(ji,jj,2) - gdepw(jk) )* v3d(ji,jj,jk) * e1t(ji,jj) / 1000000!! Modified Pedro
                       !dm2dint(ji,jj,2)=dm2dint(ji,jj,2)+ (zint(ji,jj,2) - gdepw(jk) )* SQRT(du2dint(ji,jj,2)* &
                       !                 & du2dint(ji,jj,2) + dv2dint(ji,jj,2)*dv2dint(ji,jj,2))!! Modified Pedro
                    ELSE   ! below the isopycnal 
                       ! do nothing for this i j point
                    ENDIF
                 END DO
              END DO
           END DO   ! end on vertical integral for isopynal jiso

           zdum=zint(:,:,1)

           WHERE (tmask == 0. ) zint(:,:,1)=zspval
           WHERE (tmask == 0. ) zint(:,:,2)=zspval
           !ierr = putvar(ncout,id_varout(3), zdum, jiso, npiglo, npjglo, ktime=jt ) !Saves depth of isopycnals
           ierr = putvar(ncout,id_varout(3), zint(:,:,1), 1, npiglo, npjglo, ktime=jt )
           ierr = putvar(ncout,id_varout(3), zint(:,:,2), 2, npiglo, npjglo, ktime=jt )

           !IF (jiso > 1  ) THEN  ! compute the difference ie the inventory in the layer between 2 isopycnals
              !zdum=dv2dint(:,:,1) - dv2dint(:,:,2) ; WHERE ((tmask == 0.)  .OR. (zdum < 0 ) ) zdum = zspval
           zdum_u=du2dint(:,:,2) - du2dint(:,:,1) ; WHERE (tmask == 0. ) zdum_u = zspval !! Modified Pedro
           zdum_v=dv2dint(:,:,2) - dv2dint(:,:,1) ; WHERE (tmask == 0. ) zdum_v = zspval !! Modified Pedro
           zdum_m=dm2dint(:,:,2) - dm2dint(:,:,1) ; WHERE (tmask == 0. ) zdum_m = zspval !! Modified Pedro
              !ierr = putvar(ncout, id_varout(1), zdum, jiso-1, npiglo, npjglo, ktime=jt) !Saves integral
           ierr = putvar(ncout, id_varout(5), zdum_u, 1, npiglo, npjglo, ktime=jt)
           ierr = putvar(ncout, id_varout(1), zdum_v, 1, npiglo, npjglo, ktime=jt)
           !ierr = putvar(ncout, id_varout(7), zdum_m, 1, npiglo, npjglo, ktime=jt)
     
           DO ji=1,npiglo
              DO jj=1,npjglo
                 zdum_m(ji,jj) = SQRT(zdum_v(ji,jj) * zdum_v(ji,jj) + zdum_u(ji,jj) * zdum_u(ji,jj))
                 zdum_a(ji,jj)= atan2(zdum_v(ji,jj), zdum_u(ji,jj)) * 180./pi 
                 IF ( zdum_v(ji,jj) < 0. ) zdum_a(ji,jj) = 360.+zdum_a(ji,jj)
              END DO
           END DO
           WHERE (tmask == 0. ) zdum_m=zspval
           WHERE (tmask == 0. ) zdum_a = zspval
           ierr = putvar(ncout, id_varout(7), zdum_m, 1, npiglo, npjglo, ktime=jt)
           ierr = putvar(ncout, id_varout(9), zdum_a , 1, npiglo, npjglo, ktime=jt)

              !zdum=zint  (:,:,1) - zint  (:,:,2) ; WHERE ((tmask == 0.)  .OR. (zdum < 0 ) ) zdum = zspval
           zdum=zint  (:,:,2) - zint  (:,:,1) ; WHERE (tmask == 0.) zdum = zspval
           !ierr = putvar(ncout, id_varout(2), zdum, jiso-1, npiglo, npjglo, ktime=jt) !Saves thickness
           ierr = putvar(ncout, id_varout(2), zdum, 1, npiglo, npjglo, ktime=jt)

              WHERE ( zdum /= zspval .AND. zdum /= 0.) 
                 zdum_u=(du2dint(:,:,2) - du2dint(:,:,1))/ zdum
                 zdum_v=(dv2dint(:,:,2) - dv2dint(:,:,1))/ zdum
                 zdum_m = zdum_m / zdum
              !ELSEWHERE
                 !zdum_u=zspval
                 !zdum_v=zspval
                 !zdum_m=zspval
              ENDWHERE
              !ierr = putvar(ncout, id_varout(4), zdum, jiso-1, npiglo, npjglo, ktime=jt) !Saves average
              ierr = putvar(ncout, id_varout(6), zdum_u, 1, npiglo, npjglo, ktime=jt)
              ierr = putvar(ncout, id_varout(4), zdum_v, 1, npiglo, npjglo, ktime=jt)
              ierr = putvar(ncout, id_varout(8), zdum_m, 1, npiglo, npjglo, ktime=jt)

           !ENDIF
           !dv2dint(:,:,2) = dv2dint(:,:,1)
           !zint   (:,:,2) = zint   (:,:,1)

        END DO
     !END DO
     ierr = closeout(ncout)
  !END DO  ! loop on scalar files
  PRINT *,' integral between isopycnals completed successfully'
END  PROGRAM cdfsigintegr_pedro
