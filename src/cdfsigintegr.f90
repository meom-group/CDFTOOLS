PROGRAM cdfsigintegr
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
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class integration
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: ji, jj, jk, jt   ! dummy loop index
  INTEGER(KIND=4)                               :: jiso, jfich      ! dummy loop index
  INTEGER(KIND=4)                               :: jvar             ! dummy loop index
  INTEGER(KIND=4)                               :: it               ! time index for vvl
  INTEGER(KIND=4)                               :: npiglo, npjglo   ! domain size
  INTEGER(KIND=4)                               :: npk, npt         ! domain size
  INTEGER(KIND=4)                               :: npiso, nvars     ! number of isopycnals, variables
  INTEGER(KIND=4)                               :: narg, iargc      ! command line
  INTEGER(KIND=4)                               :: ijarg, ireq      ! command line
  INTEGER(KIND=4)                               :: nfiles           ! number of input files
  INTEGER(KIND=4)                               :: nsfiles          ! number of density input files
  INTEGER(KIND=4)                               :: istrt_arg        ! argument number of first input file
  INTEGER(KIND=4)                               :: ijk              ! layer index
  INTEGER(KIND=4)                               :: ik               ! working integer
  INTEGER(KIND=4)                               :: ik0              ! working integer
  INTEGER(KIND=4)                               :: numin=10         ! logical unit for ascii input file
  INTEGER(KIND=4)                               :: ncout, ierr      ! ncid and status variable
  INTEGER(KIND=4), DIMENSION(4)                 :: ipk, id_varout   ! levels and id's of output variables
  !
  REAL(KIND=4)                                  :: zspval=999999.   ! output missing value
  REAL(KIND=4)                                  :: spvalz          ! missing value from rho file      
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: rho_lev          ! value of isopycnals
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: h1d              ! depth of rho points
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: gdepw            ! depth of W points
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: e31d             ! vertical metrics in full step
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d              ! 2D working array
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: e3               ! vertical metrics
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: tmask            ! mask of t points from rho
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zdum             ! dummy array for I/O
  REAL(KIND=4), DIMENSION(:,:,:),   ALLOCATABLE :: v3d              ! 3D working array (npk)
  REAL(KIND=4), DIMENSION(:,:,:),   ALLOCATABLE :: zint             ! pseudo 3D working array (2)

  REAL(KIND=8)                                  :: dl_wght
  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dtim             ! time counter
  REAL(KIND=8), DIMENSION(:,:,:),   ALLOCATABLE :: dv2dint          ! interpolated value 
  REAL(KIND=8), DIMENSION(:,:,:),   ALLOCATABLE :: dalpha           ! 3D coefficient (npiso)

  CHARACTER(LEN=256)                            :: cf_rholev = 'rho_lev' ! input file for rho surfaces
  CHARACTER(LEN=256)                            :: cf_in            ! input file for data
  CHARACTER(LEN=256)                            :: cf_tsig          ! text file with density filename
  CHARACTER(LEN=256)                            :: cf_tdta          ! text file with data filename
  CHARACTER(LEN=256)                            :: cf_wk            ! working character variable
  CHARACTER(LEN=256)                            :: cf_rho           ! input file for density
  CHARACTER(LEN=256)                            :: cf_out           ! output file
  CHARACTER(LEN=256)                            :: cf_e3            ! e3 filename ( for vvl )
  CHARACTER(LEN=256)                            :: cv_in            ! name of input variable
  CHARACTER(LEN=256)                            :: cldum            ! dummy string variable
  CHARACTER(LEN=256)                            :: cluni            ! dummy string variable for variable units
  CHARACTER(LEN=256)                            :: cglobal          ! global attribute
  CHARACTER(LEN=256)                            :: ctype='T'        ! position of variable on C grid
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names         ! temporary arry for variable name in file
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cf_lst           ! list of input files
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cf_slst          ! list of input density files (option -sl)

  TYPE(variable), DIMENSION(4)                  :: stypvar          ! structure for attributes
  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypzvar         ! structure for attributes

  LOGICAL                                       :: lscli = .FALSE.  ! flag for climatological  density
  LOGICAL                                       :: lfull = .FALSE.  ! flag for full step
  LOGICAL                                       :: lchk  = .FALSE.  ! flag for missing files
  LOGICAL                                       :: lnc4  = .FALSE.  ! flag for netcdf4 output
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfsigintegr -s RHO-file | -sl LST-RHO-files | -st LST-RHO-txt ...'
     PRINT *,'              ... -l LST-files | -lt LST-DATA-txt  -v IN-var [-p C-type ] ...'
     PRINT *,'              ... [-sig sigma_name] [-full] [-nc4] [-vvl]'
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
     PRINT *,'       -v IN-var : input variable to be integrated' 
     PRINT *,'       -s RHO-file : netcdf file with already computed density' 
     PRINT *,'          or'
     PRINT *,'       -sl LST-RHO-file : a blank separated list of density files '
     PRINT *,'           It is mandatory that, when using this option, the number of time '
     PRINT *,'           frames is fully coherent with the data files.'
     PRINT *,'          or'
     PRINT *,'       -sl LST-RHO-txt : Pass a text file with the name of density files (1 per line).'
     PRINT *,'       -l LST-files : a blank separated list of model netcdf files '
     PRINT *,'              containing IN-var.'
     PRINT *,'          or'
     PRINT *,'       -lt LST-DATA-txt : Pass a text file with the name of data files (1 per line).'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-p  C-type ] : one of T U V F W which defined the position of' 
     PRINT *,'               IN-var in the model C-grid. Default is ', TRIM(ctype)
     PRINT *,'       [-sig sigma_name ] : give the name of sigma variable in RHO-file.'
     PRINT *,'               Default is ',TRIM(cn_vosigma0)
     PRINT *,'       [-full ] : indicate a full step configuration.'
     PRINT *,'       [-rholev  file] : indicates name of file defining the limits for '
     PRINT *,'               integration. Default is ', TRIM(cf_rholev)
     PRINT *,'       [-vvl ] : use time-varying vertical metrics.'
     PRINT *,'       [-nc4 ] : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'            This option is effective only if cdftools are compiled with'
     PRINT *,'            a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fzgr),' and ',TRIM(cf_rholev)
     PRINT *,'      '
     PRINT *,'     OPENMP SUPPORT : yes'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : IN-file.integr'
     PRINT *,'         variables : inv_<IN-var>  : inventory of IN-var from input file.'
     PRINT *,'            ', TRIM(cn_vodepiso),' (m) : depth of isopycnal.'
     PRINT *,'            ', TRIM(cn_isothick),' (m) : thickness of isopycnal layer.'
     PRINT *,'            mean_<IN-var> (same unit as IN-var) : mean IN-var in the isopycnal'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      cdfrhoproj, cdfsigtrp, cdfisopycdep'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1 ; ireq = 0 ; nfiles = 0 ; nsfiles = 0
  cf_tsig='none' ; cf_tdta='none'
  DO WHILE ( ijarg <= narg ) 
     CALL getarg( ijarg, cldum ) ; ijarg = ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-v'      ) ; CALL getarg( ijarg, cv_in      ) ; ijarg = ijarg+1 ; ireq=ireq+1
     CASE ( '-s'      ) ; CALL getarg( ijarg, cf_rho     ) ; ijarg = ijarg+1 ; ireq=ireq+1
      ! or 
     CASE ( '-sl'     ) ; CALL GetFileList(cf_slst,nsfiles)                  ; ireq=ireq+1
      ! or 
     CASE ( '-st'     ) ; CALL getarg( ijarg, cf_tsig    );  ijarg = ijarg+1 ; ireq=ireq+1
     CASE ( '-l'      ) ; CALL GetFileList(cf_lst,nfiles)                    ; ireq=ireq+1
      ! or 
     CASE ( '-lt'     ) ; CALL getarg( ijarg, cf_tdta    );  ijarg = ijarg+1 ; ireq=ireq+1
        ! options
     CASE ( '-p'      ) ; CALL getarg( ijarg, ctype      ) ; ijarg = ijarg+1 
     CASE ( '-sig'    ) ; CALL getarg( ijarg, cn_vosigma0) ; ijarg = ijarg+1 
     CASE ( '-rholev ') ; CALL getarg( ijarg, cf_rholev  ) ; ijarg = ijarg+1
     CASE ( '-full '  ) ; lfull  = .TRUE.
     CASE ( '-vvl'    ) ; lg_vvl = .TRUE.
     CASE ( '-nc4'    ) ; lnc4   = .TRUE.
     CASE DEFAULT       ; PRINT *, ' ERROR : ', TRIM(cldum), 'unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( ireq /= 3 ) THEN 
      PRINT *,' Exactly 3 options are mandatory. Look to usage message !' 
      PRINT *,'   In particular you cannot have both -s AND -sl  or -st options !' ; STOP 99;
  ENDIF

  IF ( cf_tdta /= 'none' ) THEN
     PRINT *,' Setting data file list from ',TRIM(cf_tdta)
     CALL GetTxtFileList(cf_tdta,cf_lst,nfiles)
  ENDIF

  IF ( cf_tsig /= 'none' ) THEN
     PRINT *,' Setting density file list from ',TRIM(cf_tsig)
     CALL GetTxtFileList(cf_tsig,cf_slst,nsfiles)
  ENDIF

  IF ( nsfiles == 0 ) THEN  ! means there is a climatological density file
     lscli = .true.
     PRINT *, " Use a single climatological  density files."
     ALLOCATE(cf_slst(1) )
     cf_slst = cf_rho
  ELSE
     PRINT *, " Use ",nsfiles,"  density files."
     PRINT *, " Use ",nfiles,"  data files."
     ! nsfiles should be equal to nfiles in case of multiple files
     lscli = .false.
     IF ( nsfiles /= nfiles) STOP ' ERROR: Number of density files differs form number of data file'
     cf_rho=cf_slst(1)
  ENDIF

  CALL SetGlobalAtt( cglobal )

  ! check for files
  lchk = lchk .OR. chkfile (cn_fzgr   )
  lchk = lchk .OR. chkfile (cf_rholev )
  IF ( lscli) lchk = lchk .OR. chkfile (cf_rho    )
  IF ( lchk ) STOP 99 ! missing file

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

  ! assume here that all files have the same geometry (as well as density files) use the 1rst of the list
  cf_wk = cf_slst(1)
  IF ( chkfile (cf_wk) ) STOP 'missing density file 1' 
  npiglo = getdim(cf_wk, cn_x)
  npjglo = getdim(cf_wk, cn_y)
  npk    = getdim(cf_wk, cn_z)

  spvalz=getspval(cf_wk, cn_vosigma0)

  cf_in =  cf_lst(1)
  IF ( chkfile ( cf_in) ) STOP 'missing data file 1'  ! missing file

  nvars=getnvar(cf_in )
  ALLOCATE(cv_names(nvars), stypzvar(nvars))

  cv_names(:)=getvarname(cf_in,nvars,stypzvar)

  ALLOCATE( v3d(npiglo,npjglo,npk), dalpha(npiglo,npjglo,npiso), e3(npiglo,npjglo) )
  ALLOCATE( dv2dint(npiglo,npjglo,2), v2d(npiglo,npjglo), zint(npiglo,npjglo,2)  )
  ALLOCATE( h1d(npk) ,gdepw(npk) ,tmask(npiglo,npjglo), zdum(npiglo,npjglo) )
  IF ( lfull ) ALLOCATE ( e31d(npk) )

  gdepw(:) = getvare3(cn_fzgr, cn_gdepw, npk)

  IF (lfull ) e31d(:) = getvare3(cn_fzgr, cn_ve3t1d, npk)

  h1d(:) = getvar1d(cf_rho, cn_vdeptht, npk)

  ! Note, if working with vertical slabs, one may avoid 3D array, but may be slow ...
  IF ( lscli ) THEN
    CALL Weight( 1 )
  ENDIF

  ! Create output variables
  CALL CreateOutputVar

  !! ** Loop on the scalar files to project on choosen isopycnics surfaces
  DO jfich=1, nfiles
     cf_in = cf_lst( jfich)
     IF ( .not. lscli ) THEN  ! need to process every density files
        cf_rho=cf_slst(jfich) 
        IF ( chkfile (cf_rho) ) STOP 99 ! missing file
     ENDIF

     IF ( chkfile (cf_in) ) STOP 99 ! missing file
     PRINT *,'working with ', TRIM(cf_in)
     PRINT *,'    and with ', TRIM(cf_rho)
     ! JMM : not obvious to find file wirh correct e3t
     IF (lg_vvl )  THEN
        cn_fe3t = cf_rho
        cn_ve3t = cn_ve3tvvl
     ENDIF

     ! create output file
     cf_out=TRIM(cf_in)//'.integr'

     ! creation of output file is done within the file loop, but do not interfere with 
     ! possible parallelization as output files name are different.
     ncout = create      (cf_out, cf_rho,  npiglo, npjglo, npiso                      , ld_nc4=lnc4 )
     ierr  = createvar   (ncout,  stypvar, 4,      ipk,    id_varout, cdglobal=cglobal, ld_nc4=lnc4 )
     ierr  = putheadervar(ncout,  cf_rho,  npiglo, npjglo, npiso, pdep=rho_lev                      )

     ! copy time arrays in output file
     npt = getdim ( cf_in, cn_t)
     ALLOCATE ( dtim(npt) )
     dtim = getvar1d(cf_in, cn_vtimec, npt     )
     ierr = putvar1d(ncout, dtim,      npt, 'T')
     DEALLOCATE ( dtim )

     DO jt =1, npt
        IF ( lg_vvl) THEN ; it=jt
        ELSE              ; it=1
        ENDIF
        IF ( .not. lscli ) THEN  ! need to process every density files
           CALL Weight( jt )
        ENDIF
        DO jk=1,npk
           v2d(:,:) = getvar(cf_in, cv_in, jk, npiglo, npjglo, ktime = jt )
           SELECT CASE ( ctype )
           CASE ('T', 't' ) ; v3d(:,:,jk) = v2d(:,:)
           CASE ('U','u'  ) ; 
              DO jj=1,npjglo
                 DO ji=2, npiglo
                    v3d(ji,jj,jk)=0.5*( v2d(ji,jj) + v2d(ji-1,jj) )  ! put variable on T point
                 END DO
              END DO
           CASE ('V','v' )
              DO jj=2,npjglo
                 DO ji=1, npiglo
                    v3d(ji,jj,jk)=0.5*( v2d(ji,jj) + v2d(ji,jj-1) )  ! put variable on T point
                 END DO
              END DO
           CASE('W','w' )
              v3d(:,:,jk) = v2d(:,:)
              v2d(:,:) = getvar(cf_in, cv_in, jk+1, npiglo, npjglo, ktime = jt )
              v3d(:,:,jk) = 0.5 * ( v3d(:,:,jk) + v2d(:,:) )
           CASE('F','f' )
              DO jj = 2, npjglo
                 DO ji = 2, npiglo
                    v3d(:,:,jk) = 0.25*( v2d(ji,jj) + v2d( ji, jj-1) + v2d (ji-1,jj-1) + v2d(ji-1, jj) )
                 END DO
              END DO
           END SELECT
        END DO

        ! Compute integral from surface to isopycnal
        DO jiso=1,npiso
           ! determine isopycnal surface
           DO ji=1,npiglo
              DO jj=1,npjglo
                 ! 
                 ik      = INT(dalpha(ji,jj,jiso))
                 dl_wght = dalpha(ji,jj,jiso) - ik
                 IF (ik /= 0) THEN
                    zint (ji,jj,1)=dl_wght*h1d(ik+1) + (1.d0-dl_wght)*h1d(ik)
                 ELSE 
                    zint  (ji,jj,1)=0.  !zspval  
                 ENDIF
              END DO
           END DO
           ! integrate from jk=1 to zint
           dv2dint(:,:,1) = 0.d0

           DO jk=1,npk-1
              ! get metrixs at level jk
              IF ( lfull ) THEN  ; e3(:,:) = e31d(jk)
              ELSE               ; e3(:,:) = getvar(cn_fe3t, cn_ve3t, jk,npiglo,npjglo, ktime=it, ldiom=.NOT.lg_vvl )
              ENDIF

              DO ji=1,npiglo
                 DO jj=1,npjglo
                    IF ( gdepw(jk)+e3(ji,jj) < zint(ji,jj,1) ) THEN  ! full cell
                       dv2dint(ji,jj,1)=dv2dint(ji,jj,1) + e3(ji,jj)* v3d(ji,jj,jk)
                    ELSE IF (( zint(ji,jj,1) <= gdepw(jk)+e3(ji,jj) ) .AND. (zint(ji,jj,1) > gdepw(jk)) ) THEN
                       dv2dint(ji,jj,1)=dv2dint(ji,jj,1)+ (zint(ji,jj,1) - gdepw(jk) )* v3d(ji,jj,jk)
                    ELSE   ! below the isopycnal 
                       ! do nothing for this i j point
                    ENDIF
                 END DO
              END DO
           END DO   ! end on vertical integral for isopynal jiso

           zdum=zint(:,:,1)

           WHERE (tmask == 0. ) zdum=zspval
           ierr = putvar(ncout,id_varout(3), zdum, jiso, npiglo, npjglo, ktime=jt )

           IF (jiso > 1  ) THEN  ! compute the difference ie the inventory in the layer between 2 isopycnals
              zdum=dv2dint(:,:,1) - dv2dint(:,:,2) ; WHERE ((tmask == 0.)  .OR. (zdum < 0 ) ) zdum = zspval
              ierr = putvar(ncout, id_varout(1), zdum, jiso-1, npiglo, npjglo, ktime=jt)

              zdum=zint  (:,:,1) - zint  (:,:,2) ; WHERE ((tmask == 0.)  .OR. (zdum < 0 ) ) zdum = zspval
              ierr = putvar(ncout, id_varout(2), zdum, jiso-1, npiglo, npjglo, ktime=jt)

              WHERE ( zdum /= zspval .AND. zdum /= 0.)  ; zdum=(dv2dint(:,:,1) - dv2dint(:,:,2))/ zdum
              ELSEWHERE                                 ; zdum=zspval
              ENDWHERE

              ierr = putvar(ncout, id_varout(4), zdum, jiso-1, npiglo, npjglo, ktime=jt)
           ENDIF
           dv2dint(:,:,2) = dv2dint(:,:,1)
           zint   (:,:,2) = zint   (:,:,1)

        END DO
     END DO
     ierr = closeout(ncout)
  END DO  ! loop on scalar files
  PRINT *,' integral between isopycnals completed successfully'

CONTAINS
  SUBROUTINE Weight ( kt) 
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE Weight  ***
    !!
    !! ** Purpose :  Compute weights for vertical interpoaltion
    !!
    !! ** Method  :  linear interpolation
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) :: kt
    ! local variables
    INTEGER(KIND=4 )  :: ji,jj,jk, jiso
    INTEGER(KIND=4)   :: ijk
    REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: z3d
    !!----------------------------------------------------------------------
    ALLOCATE (z3d(npiglo,npjglo,npk) )
    ! in this case, this is to be done once only !
    DO jk=1,npk
       z3d(:,:,jk) = getvar(cf_rho, cn_vosigma0, jk, npiglo, npjglo,ktime=kt)
       IF ( jk == 1 .AND. kt == 1 ) THEN
          tmask=1.
          WHERE (z3d(:,:,jk) == spvalz ) tmask=0.
       ENDIF
    END DO
    !$OMP PARALLEL DO SCHEDULE(RUNTIME)
    DO ji=1,npiglo
       DO jj = 1, npjglo
          ijk = 1
          DO jiso=1,npiso
             !  Assume that rho (z) is increasing downward (no inversion)
             !     Caution with sigma0 at great depth !
             DO WHILE (rho_lev(jiso) >=  z3d(ji,jj,ijk) .AND. ijk <= npk &
                  &                .AND. z3d(ji,jj,ijk) /=  spvalz )
                ijk = ijk+1
             END DO
             ijk = ijk-1
             ik0 = ijk
             IF (ijk == 0) THEN
                ijk = 1
                dalpha(ji,jj,jiso) = 0.d0
             ELSE IF (z3d(ji,jj,ijk+1) == spvalz ) THEN
                ik0 = 0
                dalpha(ji,jj,jiso) = 0.d0
             ELSE 
                ! ... dalpha is always in [0,1] we add ik0
                dalpha(ji,jj,jiso)= (rho_lev(jiso)-z3d(ji,jj,ijk))/(z3d(ji,jj,ijk+1)-z3d(ji,jj,ijk)) +ik0
             ENDIF
          END DO
       END DO
    END DO
    !$OMP END PARALLEL DO
    DEALLOCATE (z3d )

  END SUBROUTINE Weight

  SUBROUTINE CreateOutputVar
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutputVar  ***
    !!
    !! ** Purpose :  Create netcdf output  variables
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ! define header of all files
    ipk(1)=npiso-1 ; ipk(2)=npiso-1 ; ipk(3)=npiso ; ipk(4)=npiso-1

    DO jvar=1,nvars
       IF ( cv_in == stypzvar(jvar)%cname ) THEN 
          stypvar(1)=stypzvar(jvar)
          EXIT
       ENDIF
    END DO
    ! save original long name for further process
    cldum = TRIM(stypvar(1)%clong_name)
    cluni = TRIM(stypvar(1)%cunits)

    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = 'inv'//TRIM(cv_in)
    stypvar(1)%clong_name        = TRIM(cldum)//' integrated on sigma bin'
    stypvar(1)%cshort_name       = stypvar(1)%cname
    stypvar(1)%cunits            = TRIM(cluni)//'.m'
    stypvar(1)%rmissing_value    = zspval
    stypvar(1)%caxis             = 'TRYX'

    stypvar(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname             = TRIM(cn_isothick)
    stypvar(2)%cunits            = 'm'
    stypvar(2)%rmissing_value    = zspval
    stypvar(2)%valid_min         = 0.
    stypvar(2)%valid_max         = 7000.
    stypvar(2)%clong_name        = 'Thickness_of_Isopycnals'
    stypvar(2)%cshort_name       = TRIM(cn_isothick)
    stypvar(2)%conline_operation = 'N/A'
    stypvar(2)%caxis             = 'TRYX'

    stypvar(3)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(3)%cname             = TRIM(cn_vodepiso)
    stypvar(3)%cunits            = 'm'
    stypvar(3)%rmissing_value    = zspval
    stypvar(3)%valid_min         = 0.
    stypvar(3)%valid_max         = 7000.
    stypvar(3)%clong_name        = 'Depth_of_Isopycnals'
    stypvar(3)%cshort_name       = TRIM(cn_vodepiso)
    stypvar(3)%conline_operation = 'N/A'
    stypvar(3)%caxis             = 'TRYX'

    stypvar(4)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(4)%cname             = 'mean'//TRIM(cv_in)
    stypvar(4)%cunits            = TRIM(cluni)
    stypvar(4)%rmissing_value    = zspval
    stypvar(4)%valid_min         = stypvar(1)%valid_min
    stypvar(4)%valid_max         = stypvar(1)%valid_max
    stypvar(4)%clong_name        = TRIM(cldum)//' mean value in sigma layer'
    stypvar(4)%cshort_name       = stypvar(4)%cname
    stypvar(4)%conline_operation = 'N/A'
    stypvar(4)%caxis             = 'TRYX'

  END SUBROUTINE CreateOutputVar

  SUBROUTINE GetFileList(cd_lst, kfiles)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetFileList  ***
    !!
    !! ** Purpose :  Set up a file list given on the command line as 
    !!               blank separated list
    !!
    !! ** Method  :  Scan the command line until a '-' is found
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: cd_lst
    INTEGER(KIND=4),                             INTENT(out)   :: kfiles
    INTEGER (KIND=4)  :: icur 
    !!----------------------------------------------------------------------
    !!
    kfiles=0
    ! need to read a list of file ( number unknow ) 
    ! loop on argument till a '-' is found as first char
    icur=ijarg                          ! save current position of argument number
    DO ji = icur, narg                  ! scan arguments till - found
       CALL getarg ( ji, cldum )
       IF ( cldum(1:1) /= '-' ) THEN ; kfiles = kfiles+1
       ELSE                          ; EXIT
       ENDIF
    ENDDO
    ALLOCATE (cd_lst(kfiles) )
    DO ji = icur, icur + kfiles -1
       CALL getarg(ji, cd_lst( ji -icur +1 ) ) ; ijarg=ijarg+1
    END DO
  END SUBROUTINE GetFileList
  
  SUBROUTINE GetTxtFileList(cd_txt, cd_lst, kfiles)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetTxtFileList  ***
    !!
    !! ** Purpose :  Set up a file list read in the input text file
    !!
    !! ** Method  :  read an count the lines of input file
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),                            INTENT(in   ) :: cd_txt
    CHARACTER(LEN=*), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: cd_lst
    INTEGER(KIND=4),                             INTENT(  out) :: kfiles

    INTEGER(KIND=4)                                            :: jl
    INTEGER(KIND=4)                                            :: inum = 11
    !!----------------------------------------------------------------------
    kfiles=0
    OPEN(inum,file=cd_txt)
    DO 
      READ(inum,*,END=999)
      kfiles=kfiles+1
    ENDDO
999 CONTINUE
    ALLOCATE( cd_lst(kfiles) )
    REWIND(inum)
    DO jl = 1, kfiles
      READ(inum,'(a)') cd_lst(jl)
    ENDDO
    CLOSE(inum)
    
  END SUBROUTINE GetTxtFileList

END  PROGRAM cdfsigintegr
