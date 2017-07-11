PROGRAM cdfrhoproj
  !!======================================================================
  !!                     ***  PROGRAM  cdfrhoproj  ***
  !!=====================================================================
  !!  ** Purpose : This program is used to project any scalar on the A grid
  !!               onto given isopycnic surfaces.
  !!
  !!  ** Method  : Linear interpolation is used on the vertical to define
  !!               the depth of the given isopycn and linear interpolation
  !!               is also performed on the scalar to determine its value at
  !!               this depth.
  !!
  !! History :      :  1996    : J.M. Molines for SPEM in Dynamo
  !!                :  1999    : J.O. Beismann for OPA
  !!                :  2000    : J.M. Molines for normalization
  !!           2.1  : 11/2005  : J.M. Molines : netcdf 
  !!           3.0  : 05/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class derived_fields
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: ji, jj, jk, jsig, jt ! dummy loop index
  INTEGER(KIND=4)                               :: jfich, jvar          ! dummy loop index
  INTEGER(KIND=4)                               :: npiglo, npjglo       ! domain size
  INTEGER(KIND=4)                               :: npk, npsig=1, npt    ! domain size
  INTEGER(KIND=4)                               :: nvars, nvout=2       ! number of variables in/out
  INTEGER(KIND=4)                               :: narg, iargc          ! command line parser
  INTEGER(KIND=4)                               :: ijarg                ! command line parser
  INTEGER(KIND=4)                               :: ik0, ijk             ! working integer
  INTEGER(KIND=4)                               :: nfiles               ! file counter
  INTEGER(KIND=4)                               :: numlev=10            ! logical unit for rholev file
  INTEGER(KIND=4)                               :: ncout, ierr          ! netcdf working variables
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk, id_varout       ! for output variables
  !
  REAL(KIND=4), DIMENSION(:,:,:),   ALLOCATABLE :: zsig, alpha          ! data and interp coef 3D
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2dint, zint, v2d    ! working 2D arrays
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: zi, tim, h1d         ! working 1D arrays
  REAL(KIND=4)                                  :: P1, P2               ! values used in the vertical interpolation
  REAL(KIND=4)                                  :: zalpha               ! working real
  REAL(KIND=4)                                  :: zspvalo=999999.      ! output special value
  REAL(KIND=4)                                  :: zspvali=0.           ! input special values
  REAL(KIND=4)                                  :: sigmin, sigstp       ! definition of sigma layers

  CHARACTER(LEN=256)                            :: cf_rholev='rho_lev'  ! default name of rho_lev file
  CHARACTER(LEN=256)                            :: cf_dta               ! working input file
  CHARACTER(LEN=256)                            :: cf_rhofil            ! density reference file
  CHARACTER(LEN=256)                            :: cf_out               ! output file name
  CHARACTER(LEN=256)                            :: cf_iso = 'isopycdep.nc' ! default isodep file name
  CHARACTER(LEN=256)                            :: cv_in                ! input variable
  CHARACTER(LEN=256)                            :: cv_sig               ! default density variable name
  CHARACTER(LEN=256)                            :: ctype='T'            ! default C-grid type
  CHARACTER(LEN=256)                            :: cldum                ! working char variable
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names             ! temporary array for variable name in file
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cf_lst               ! list of input files

  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvar              ! structure for attributes (out)
  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypzvar             ! structure for attributes (in)
  !
  LOGICAL                                       :: lsingle = .FALSE.    ! flag for use of -s0 option
  LOGICAL                                       :: lchk    = .FALSE.    ! flag for file existence
  LOGICAL                                       :: lisodep = .FALSE.    ! flag for isodep only computation
  LOGICAL                                       :: liso    = .TRUE.     ! flag for isodep computation
  LOGICAL                                       :: ldebug  = .FALSE.    ! flag for extra debugging print
  LOGICAL                                       :: lnc4    = .FALSE.    ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  cv_sig = cn_vosigma0

  narg=iargc()
  IF ( narg == 0  ) THEN
     PRINT *,' usage : cdfrhoproj-v IN-var -s RHO-file -l LST-files [-p C-type] [-debug]...'
     PRINT *,'       ... [-isodep] [-s0 sig0 | -s0 sigmin,sigstp,nsig] [-sig sigma_name]..'
     PRINT *,'       ... [-noiso] [-rholev TXT-file] [-o OUT-isodep] [-nc4] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       This program aims at projecting the model variable IN-var, from a list'
     PRINT *,'       of model files (LST-files) on some isopycnic surfaces, inferred from'
     PRINT *,'       a 3D density file, passed as one the arguments of the program. '
     PRINT *,'      '
     PRINT *,'       The density values corresponding to the isopycnic surfaces can be'
     PRINT *,'       specified in three ways :'
     PRINT *,'          1. Using a predefined text file ',TRIM(cf_rholev),' with density'
     PRINT *,'             values. The format is straightforward: one value per line, first'
     PRINT *,'             line giving the number of isopycnic to consider. (The default'
     PRINT *,'             name of this text file can be changed using -rholev option).'
     PRINT *,'          2. Using the -s0 option, with 3 parameters, defining equally spaced'
     PRINT *,'             isopycnic surfaces.'
     PRINT *,'          3. Using the -s0 option with only one parameter, defining then a'
     PRINT *,'             single surface.'
     PRINT *,'       '
     PRINT *,'     WARNING: This cdftool is one of the few using 3D arrays. Additional'
     PRINT *,'         development is required to work with vertical slabs instead.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -v IN-var   : name of the input variable to be projected' 
     PRINT *,'       -s RHO-file : netcdf file with potential density field. If not a sigma0'
     PRINT *,'             file, use -sig option to indicate the name of the density '
     PRINT *,'             variable.'
     PRINT *,'       -l LST-files: List of netcdf file with variable IN-var to process.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-s0 sigma  | -s0 sigmin,sigstp,nsig] : In the first form define a '
     PRINT *,'            single sigma surface, while in the 2nd form, it defines a set of'
     PRINT *,'            ''nsig'' density values, starting from ''sigmin'' and spaced every'
     PRINT *,'            ''sigstp''. This option prevails the use of ',TRIM(cf_rholev),' file.'
     PRINT *,'       [-rholev TXT-file] : Specify the name of the ''rholev'' text file, '
     PRINT *,'            instead of ',TRIM(cf_rholev),'.'
     PRINT *,'       [-p C-type] : position of IN-var on the C-grid ( either T U V F W S ),'
     PRINT *,'            default is ''T''. ''S'' special point is used in case of section'
     PRINT *,'            files created by cdf_xtract_brokenline.'
     PRINT *,'       [-sig sigma_name] : name of the density variable in RHO_file. Default is'
     PRINT *,'           ',TRIM(cv_sig),'.'
     PRINT *,'       [-isodep ] : Only computes the isopycnic depth, then stops. '
     PRINT *,'       [-noiso]   : Does not save isopycnic depth (suitable for big files).'
     PRINT *,'       [-debug]   : Produces extra prints.'
     PRINT *,'       [-o OUT-isodep]: specify the name of isodep file (-isodep option),'
     PRINT *,'           instead of ',TRIM(cf_iso),'.'
     PRINT *,'       [ -nc4 ]   : Use netcdf4 output with chunking and deflation level 1..'
     PRINT *,'            This option is effective only if cdftools are compiled with'
     PRINT *,'            a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       If not using -s0 option ',TRIM(cf_rholev), 'is required, unless '
     PRINT *,'       the default name is changed by -rholev option.'
     PRINT *,'     ' 
     PRINT *,'     OPENMP SUPPORT : yes'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       There are as many output files than input files, with ''.interp'' '
     PRINT *,'       suffix added to the original name.'
     PRINT *,'       netcdf files : <IN-file>.interp'
     PRINT *,'         variables : VAR-in (unit is the same as input var)'
     PRINT *,'                     ', TRIM(cn_vodepiso),' (m) : depth of isopycnic.'
     PRINT *,'      '
     PRINT *,'       If option -isodep is used, only isopycnic depth is output on '
     PRINT *,'       netcdf file : ',TRIM(cf_iso),' (unless -o option is used).'
     PRINT *,'         variables : ',TRIM(cn_vodepiso),' (m) '
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfmocsig'
     PRINT *,'       '
     STOP
  ENDIF

  ijarg = 1 

  DO WHILE ( ijarg <= narg )
     CALL getarg( ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-v'    ) ; CALL getarg( ijarg, cv_in    ) ; ijarg=ijarg+1
     CASE ( '-s'    ) ; CALL getarg( ijarg, cf_rhofil) ; ijarg=ijarg+1
     CASE ( '-l'    ) ; CALL GetFileList
        ! options
     CASE ( '-s0'   ) ; CALL ParseS0Val
     CASE ( '-p'    ) ; CALL getarg( ijarg, ctype    ) ; ijarg=ijarg+1
     CASE ('-sig'   ) ; CALL getarg( ijarg, cv_sig   ) ; ijarg=ijarg+1 
     CASE ('-isodep') ; lisodep = .TRUE.  ; nvout=1 
     CASE ('-noiso' ) ; liso    = .FALSE. ; nvout=1
     CASE ('-debug' ) ; ldebug  = .TRUE.
     CASE ('-rholev') ; CALL getarg( ijarg, cf_rholev) ; ijarg=ijarg+1
     CASE ( '-o'    ) ; CALL getarg( ijarg, cf_iso   ) ; ijarg=ijarg+1
     CASE ( '-nc4'  ) ; lnc4    = .TRUE.
     CASE DEFAULT     ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP
     END SELECT
  END DO

  lchk = chkfile(cf_rhofil)
  IF ( .NOT. lsingle ) lchk = lchk .OR. chkfile(cf_rholev)
  IF ( lchk ) STOP ! missing file

  IF ( .NOT.  lsingle ) THEN
     OPEN(numlev,FILE=cf_rholev)
     READ(numlev,*) npsig
     IF (ldebug) PRINT *, TRIM(cf_rholev),' contains :'
     IF (ldebug) PRINT *, npsig
     ALLOCATE ( zi(npsig) )
     DO jsig=1,npsig
        READ(numlev,*)      zi(jsig)
        IF (ldebug) PRINT *,zi(jsig)
     END DO
     CLOSE(numlev)
  ENDIF

  ! Read reference density  file
  npiglo = getdim(cf_rhofil,cn_x)
  npjglo = getdim(cf_rhofil,cn_y)
  npk    = getdim(cf_rhofil,cn_z)
  npt    = getdim(cf_rhofil,cn_t)

  IF ( npt > 1 ) THEN
     PRINT *,'  WARNING : more than 1 time step in ',TRIM(cf_rhofil)
     PRINT *,'            Only the first one will be used.'
  ENDIF

  cf_dta = cf_lst(1)   ! work with the first file of the list for dimension determination
  nvars=getnvar(cf_dta)

  ALLOCATE(cv_names(nvars), stypzvar(nvars))
  ALLOCATE(ipk(nvout), id_varout(nvout), stypvar(nvout) )

  cv_names(:)=getvarname(cf_dta, nvars, stypzvar)

  ALLOCATE( zsig(npiglo,npjglo,npk), alpha(npiglo, npjglo, npsig)            )  ! 3D arrays !!
  ALLOCATE( v2dint(npiglo, npjglo), v2d(npiglo,npjglo), zint(npiglo,npjglo)  )
  ALLOCATE( tim(npt), h1d(npk)                                               )
alpha=0.

  tim(:)=getvar1d(cf_rhofil, cn_vtimec,  npt)
  h1d(:)=getvar1d(cf_rhofil, cn_vdeptht, npk)

  DO jk=1,npk
     zsig(:,:,jk) = getvar(cf_rhofil, cv_sig, jk, npiglo, npjglo)
  END DO

  !! Compute interpolation coefficients as well as the level used
  !! to interpolate between
  !!  This work is done once, on the RHO file.
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)
  DO ji=1,npiglo
     DO jj = 1, npjglo
        ijk = 1
        DO jsig=1,npsig
           !  Assume that rho (z) is increasing downward (no inversion)
           !     Caution with sigma0 at great depth !
           DO WHILE (zi(jsig) >=  zsig(ji,jj,ijk) .AND. ijk < npk &
                &   .AND. zsig(ji,jj,ijk) /=  zspvali )
              ijk=ijk+1
           END DO
           ijk=ijk-1
           ik0=ijk
           IF (ijk == 0) THEN
              ijk=1
              alpha(ji,jj,jsig) = 0.
           ELSE IF (zsig(ji,jj,ijk+1) == zspvali ) THEN
              ik0=0
              alpha(ji,jj,jsig) = 0.
           ELSE 
              ! ... alpha is always in [0,1]. Adding ik0 ( >=1 ) for saving space for ik0
              alpha(ji,jj,jsig)= &
                   &  (zi(jsig)-zsig(ji,jj,ijk))/(zsig(ji,jj,ijk+1)-zsig(ji,jj,ijk)) + ik0
           ENDIF
        END DO
     END DO
  END DO
  !$OMP END PARALLEL DO

  IF ( lisodep ) THEN
     CALL CreateOutputIsodep
     DO jsig=1,npsig
        !$OMP PARALLEL DO SCHEDULE(RUNTIME)
        DO ji=1,npiglo
           DO jj=1,npjglo
              ! ik0 is retrieved from alpha, taking the integer part.
              ! The remnant is alpha.
              ik0 = INT(alpha(ji,jj,jsig))
              zalpha =  alpha(ji,jj,jsig) - ik0
              IF (ik0 /= 0) THEN
                 P1 = zsig(ji,jj,ik0  )
                 P2 = zsig(ji,jj,ik0+1)
                 IF (P1 /= zspvali .AND. P2 /= zspvali) THEN
                    zint (ji,jj) = zalpha *h1d(ik0+1) &
                         &         +(1-zalpha)*h1d(ik0  )
                 ELSE
                    zint  (ji,jj)=zspvalo
                 ENDIF
              ELSE
                 zint  (ji,jj)=zspvalo
              ENDIF
           END DO
        END DO
        !$OMP END PARALLEL DO
        ierr = putvar(ncout, id_varout(1), zint , jsig, npiglo, npjglo)
     END DO
     ierr = closeout(ncout    )
     PRINT *,' -isodep option in use: only compute depth of isopycnic surfaces.'
     STOP 
  ENDIF
  DEALLOCATE ( tim) 

  !! ** Loop on the scalar files to project on choosen isopycnic surfaces
  DO jfich= 1, nfiles
     cf_dta = cf_lst(jfich)
     PRINT *,'working with ', TRIM(cf_dta)
     npt    = getdim(cf_dta, cn_t)
     CALL CreateOutput
     DO jt =1, npt

        ! fill in zsig with 3D input data interpolated at T-point
        DO jk=1,npk
           v2d(:,:) = getvar(cf_dta, cv_in, jk, npiglo, npjglo, ktime=jt)
           SELECT CASE ( ctype )
           CASE ('T', 't', 'S', 's' )
              zsig(:,:,jk) = v2d(:,:)
           CASE ('U','u' )
              DO ji=2,npiglo
                 DO jj=1, npjglo
                    zsig(ji,jj,jk)=0.5*( v2d(ji,jj) + v2d(ji-1,jj) )  ! put variable on T point
                 END DO
              END DO
           CASE ('V','v' )
              DO jj=2,npjglo
                 DO ji=1, npiglo
                    zsig(ji,jj,jk)=0.5*( v2d(ji,jj) + v2d(ji,jj-1) )  ! put variable on T point
                 END DO
              END DO
           CASE('W','w' )
              zint(:,:) = getvar(cf_dta, cv_in, jk+1, npiglo, npjglo)
              DO jj=1,npjglo
                 DO ji=1, npiglo
                    zsig(ji,jj,jk)=0.5*( v2d(ji,jj) + zint(ji,jj) )  ! put variable on T point
                 END DO
              END DO
           CASE('F','f' )
              DO jj=2,npjglo
                 DO ji=2, npiglo
                    zsig(ji,jj,jk)=0.25*( v2d(ji,jj) + v2d(ji,jj-1) + v2d(ji-1,jj) + v2d(ji-1,jj-1 )) ! put variable on T point
                 END DO
              END DO
           END SELECT
        END DO

        DO jsig=1,npsig
           !$OMP PARALLEL DO SCHEDULE(RUNTIME)
           DO ji=1,npiglo
              DO jj=1,npjglo
                 ! ik0 is retrieved from alpha, taking the integer part.
                 ! The remnant is alpha. 
                 ik0    = INT(alpha(ji,jj,jsig))
                 zalpha =     alpha(ji,jj,jsig) - ik0
                 IF (ik0 > 0 .AND. ik0 < npk ) THEN
                    P1 = zsig(ji,jj,ik0  )
                    P2 = zsig(ji,jj,ik0+1)
                    IF (P1 /= zspvali .AND. P2 /= zspvali) THEN
                       v2dint(ji,jj) = zalpha *P2  &
                            &         +(1-zalpha)*P1
                       IF( liso) zint (ji,jj) = zalpha *h1d(ik0+1) &
                            &         +(1-zalpha)*h1d(ik0  )
                    ELSE 
                       v2dint(ji,jj)=zspvalo
                       IF( liso )zint  (ji,jj)=zspvalo
                    ENDIF
                 ELSE 
                    v2dint(ji,jj)=zspvalo
                    IF ( liso ) zint  (ji,jj)=zspvalo
                 ENDIF
              END DO
           END DO
           !$OMP END PARALLEL DO
           ierr = putvar(ncout, id_varout(1), v2dint, jsig, npiglo, npjglo, ktime=jt)
           IF (liso) ierr = putvar(ncout, id_varout(2), zint  , jsig, npiglo, npjglo, ktime=jt)
        END DO
     ENDDO   ! loop on time
     ierr = closeout(ncout             )
  END DO  ! loop on scalar files
  PRINT *,'Projection on isopycns completed successfully'

CONTAINS

  SUBROUTINE ParseS0Val
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ParseS0Val  ***
    !!
    !! ** Purpose :  Parse -s0 option if used to set up the equally spaced
    !!               isopycnic surfaces to use for projection 
    !!
    !! ** Method  :  Assume cldum is a comma separated list, use global module
    !!               variables.
    !!----------------------------------------------------------------------
    CHARACTER(LEN=80), DIMENSION(100) :: cl_dum  ! 100 is arbitrary
    INTEGER  :: ji
    INTEGER  :: inchar,  i1=1, ival
    !!----------------------------------------------------------------------
    CALL getarg(ijarg, cldum); ijarg=ijarg+1
    inchar = LEN(TRIM(cldum))
    ival  = 1
    lsingle = .TRUE.  ! this flag tells the program that no rholev file is
    ! required.
    ! scan the input string and look for ',' as separator
    DO ji=1,inchar
       IF ( cldum(ji:ji) == ',' ) THEN
          cl_dum(ival) = cldum(i1:ji-1)
          i1=ji+1
          ival = ival + 1
       ENDIF
    ENDDO
    ! last name of the list does not have a ','
    cl_dum(ival) = cldum(i1:inchar)

    IF ( ival == 3 ) THEN  ! sigmin,sigstp,npsig
       READ(cl_dum(1),*) sigmin
       READ(cl_dum(2),*) sigstp
       READ(cl_dum(3),*) npsig
    ELSE IF (ival == 1 ) THEN ! single value
       READ(cl_dum(1),*) sigmin
       npsig  = 1
       sigstp = 0.
    ELSE
       PRINT *,' Error in -s0 option : either -s0 val  or -s0 sigmin,sigstp,nsig'
       STOP
    ENDIF

    ALLOCATE ( zi(npsig) )
    zi(1) = sigmin
    DO ji=2, npsig
       zi(ji) = zi(ji-1) + sigstp
    ENDDO
    IF ( ldebug ) THEN
       PRINT *, TRIM(cf_rholev),' like output '
       PRINT *, '---------------------'
       PRINT *, npsig
       DO ji = 1, npsig
          PRINT *, zi(ji)
       END DO
    ENDIF

  END SUBROUTINE ParseS0Val
  
  SUBROUTINE GetFileList
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetFileList  ***
    !!
    !! ** Purpose :  Set up a file list given on the command line as 
    !!               blank separated list
    !!
    !! ** Method  :  Scan the command line until a '-' is found
    !!----------------------------------------------------------------------
    INTEGER (KIND=4)  :: ji
    INTEGER (KIND=4)  :: icur
    !!----------------------------------------------------------------------
    !!
    nfiles=0
    ! need to read a list of file ( number unknow ) 
    ! loop on argument till a '-' is found as first char
    icur=ijarg                          ! save current position of argument number
    DO ji = icur, narg                  ! scan arguments till - found
       CALL getarg ( ji, cldum )
       IF ( cldum(1:1) /= '-' ) THEN ; nfiles = nfiles+1
       ELSE                          ; EXIT
       ENDIF
    ENDDO
    ALLOCATE (cf_lst(nfiles) )
    DO ji = icur, icur + nfiles -1
       CALL getarg(ji, cf_lst( ji -icur +1 ) ) ; ijarg=ijarg+1
    END DO
  END SUBROUTINE GetFileList

  SUBROUTINE CreateOutputIsodep
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutputIsodep  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ipk(1)                       = npsig

    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = cn_vodepiso
    stypvar(1)%cunits            = 'm'
    stypvar(1)%rmissing_value    = 999999.
    stypvar(1)%valid_min         = 0.
    stypvar(1)%valid_max         = 7000.
    stypvar(1)%clong_name        = 'Depth_of_Isopycnals'
    stypvar(1)%cshort_name       = cn_vodepiso
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TRYX'

    ncout = create      (cf_out, cf_rhofil, npiglo, npjglo, npsig         , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar,   nvout,  ipk,    id_varout     , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout , cf_rhofil, npiglo, npjglo, npsig, pdep=zi )

    ierr  = putvar1d(ncout, tim, 1, 'T')

  END SUBROUTINE CreateOutputIsodep

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ! ... open output file and write header
    cf_out=TRIM(cf_dta)//'.interp'

    ipk(:)=npsig
    DO jvar=1,nvars
       IF ( cv_in == stypzvar(jvar)%cname ) THEN 
          stypvar(1)=stypzvar(jvar)
          EXIT
       ENDIF
    END DO
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%clong_name        = TRIM(stypvar(2)%clong_name)//' on iso sigma'
    stypvar(1)%rmissing_value    = zspvalo
    stypvar(1)%caxis             = 'TRYX'

    IF ( liso ) THEN
       stypvar(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
       stypvar(2)%cname             = cn_vodepiso
       stypvar(2)%cunits            = 'm'
       stypvar(2)%rmissing_value    = 999999.
       stypvar(2)%valid_min         = 0.
       stypvar(2)%valid_max         = 7000.
       stypvar(2)%clong_name        = 'Depth_of_Isopycnes'
       stypvar(2)%cshort_name       = cn_vodepiso
       stypvar(2)%conline_operation = 'N/A'
       stypvar(2)%caxis             = 'TRYX'
    ENDIF

    ncout = create      (cf_out, cf_rhofil, npiglo, npjglo, npsig         , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar,   nvout,  ipk,    id_varout     , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout , cf_rhofil, npiglo, npjglo, npsig, pdep=zi              )

    ALLOCATE ( tim(npt) )
    tim(:)=getvar1d(cf_dta, cn_vtimec, npt)
    ierr = putvar1d(ncout, tim, npt, 'T')
    DEALLOCATE ( tim )

  END SUBROUTINE CreateOutput

END  PROGRAM cdfrhoproj
