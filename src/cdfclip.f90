PROGRAM cdfclip
  !!======================================================================
  !!                     ***  PROGRAM  cdfclip  ***
  !!=====================================================================
  !!  ** Purpose : An alternative to ncks to clip model file. It is 
  !!               usefull when the clipping area cross the E-W
  !!               periodic folding line. Additionally it does not
  !!               mess up the order of the dimensions and variables, 
  !!               which was a problem for coordinates.nc files with
  !!               IOIPSL
  !!
  !! History : 2.1  : 02/2007  : J.M. Molines : Original code
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class file_operations
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jt, jvar, jv        ! dummy loop index
  INTEGER(KIND=4)                               :: ik1, ik2, ik            !
  INTEGER(KIND=4)                               :: ierr                    ! working integer
  INTEGER(KIND=4)                               :: iimin, iimax            !
  INTEGER(KIND=4)                               :: ijmin, ijmax            !
  INTEGER(KIND=4)                               :: ikmin=-9999, ikmax=-9999  !
  INTEGER(KIND=4)                               :: narg, iargc, ijarg      ! 
  INTEGER(KIND=4)                               :: npiglo, npjglo, npk     !
  INTEGER(KIND=4)                               :: npkk, npt               ! size of the domain
  INTEGER(KIND=4)                               :: nvars                   ! Number of variables in a file
  INTEGER(KIND=4)                               :: ncout                   !
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var                  ! arrays of var id's
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk, ipkk               ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout , ndim        !

  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d, rlon, rlat         !
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2dxz, v2dyz, zxz, zyz  !
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: rdepg, rdep             !
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: tim                     !

  CHARACTER(LEN=256)                            :: cf_in                   ! input file name
  CHARACTER(LEN=256)                            :: cf_out='cdfclip.nc'     ! output file name
  CHARACTER(LEN=256)                            :: cv_dep, cv_tim          ! depth and time variable names
  CHARACTER(LEN=255)                            :: cglobal                 ! global attribute to write on output file
  CHARACTER(LEN=256)                            :: cldum                   ! dummy character variable
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names                ! array of var name

  TYPE (variable), DIMENSION(:),   ALLOCATABLE  :: stypvar                 !

  LOGICAL                                       :: lzonal=.FALSE.          !
  LOGICAL                                       :: lmeridian=.FALSE.       !
  LOGICAL                                       :: lnc4     = .FALSE.      ! Use nc4 with chunking and deflation
  !!-------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfclip -f IN-file  -zoom imin imax jmin jmax ...'
     PRINT *,'            ...[-klim kmin kmax] [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Clips the input file according to the indices given in the zoom option.'
     PRINT *,'       If no vertical zoomed area is indicated, the whole water column is '
     PRINT *,'       considered.'
     PRINT *,'      '
     PRINT *,'       This program is able to extract data for a region crossing the E-W '
     PRINT *,'       periodic boundary of a global configuration. It does so if imax < imin.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file : specify the input file to be clipped' 
     PRINT *,'       -zoom imin imax jmin jmax : specify the domain to be extracted.'
     PRINT *,'             If imin=imax, or jmin = jmax assume a vertical section either '
     PRINT *,'             meridional or zonal.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file ] : use OUT-file instead of ',TRIM(cf_out),' for output file'
     PRINT *,'       [-nc4 ] : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'             This option is effective only if cdftools are compiled with'
     PRINT *,'             a netcdf library supporting chunking and deflation.'
     PRINT *,'       [-klim kmin kmax ] : specify vertical limits for the zoom, in order to '
     PRINT *,'             reduce the extracted area to some levels. Default is to take the' 
     PRINT *,'             while water column.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out),' This can be changed using -o option'
     PRINT *,'         variables : same as input variables.'
     STOP
  ENDIF
  !!
  ijarg=1 
  DO WHILE (ijarg <= narg )
     CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum)
     CASE ('-f'   ) ; CALL getarg(ijarg,cf_in ) ; ijarg = ijarg + 1 
     CASE ('-zoom') ; CALL getarg(ijarg,cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) iimin
        ;             CALL getarg(ijarg,cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) iimax
        ;             CALL getarg(ijarg,cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmin
        ;             CALL getarg(ijarg,cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmax
     CASE ('-klim') ; CALL getarg(ijarg,cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) ikmin
        ;             CALL getarg(ijarg,cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) ikmax
     CASE ('-o'   ) ; CALL getarg(ijarg,cf_out) ; ijarg = ijarg + 1 
     CASE ('-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT   ; PRINT *,' ERROR : ', TRIM(cldum) ,' : unknwn option.';  STOP
     END SELECT
  ENDDO

  IF ( chkfile (cf_in ) ) STOP ! missing file

  ! set global attribute for output file
  IF ( ikmin > 0 ) THEN
     WRITE(cglobal,'(a,a,a,6i5)') 'cdfclip -f ',TRIM(cf_in),' -zoom ',iimin,iimax,ijmin,ijmax, ikmin, ikmax
  ELSE
     WRITE(cglobal,'(a,a,a,4i5)') 'cdfclip -f ',TRIM(cf_in),' -zoom ',iimin,iimax,ijmin,ijmax
  ENDIF

  IF ( iimin == iimax ) THEN ; lmeridian=.TRUE.; PRINT *,' Meridional section ' ;
  ENDIF
  IF ( ijmin == ijmax ) THEN ; lzonal=.TRUE.   ; PRINT *,' Zonal section '      ;
  ENDIF

  IF (iimax < iimin ) THEN ! we assume that this is the case when we cross the periodic line in orca (Indian ocean)
     npiglo= getdim (cf_in,cn_x)
     npiglo=iimax+(npiglo-iimin) -1
  ELSE
     npiglo= iimax-iimin+1
  ENDIF

  npjglo= ijmax-ijmin+1

  ! look for possible name for vertical dim                       :
  npk   = getdim (cf_in,cn_z,cdtrue=cv_dep, kstatus=ierr)     ! depthxxx
  PRINT *,'ist',ierr,TRIM(cn_z)
  IF (ierr /= 0 ) THEN
     npk   = getdim (cf_in,'z',cdtrue=cv_dep,kstatus=ierr)       ! zxxx
     PRINT *,'ist',ierr,'z'
     IF (ierr /= 0 ) THEN
        npk   = getdim (cf_in,'sigma',cdtrue=cv_dep,kstatus=ierr) ! sigmaxxx
        PRINT *,'ist',ierr,'sigma'
        IF (ierr /= 0 ) THEN
           PRINT *,' assume file with no depth'
           IF ( ikmin > 0 ) THEN
              PRINT *,' You cannot specify limits on k level !' ; STOP
           ENDIF
           npk=0  ! means no dim level in file (implicitly 1 level)
        ENDIF
     ENDIF
  ENDIF

  ! replace flag value (-9999) by standard value (no ikmin ikmax specified = whole column)
  IF ( ikmin < 0 ) ikmin = 1
  IF ( ikmax < 0 ) ikmax = npk
  npkk = ikmax - ikmin +1   ! number of extracted levels. If no level in file, it is 0: 0 -1 + 1 !
  IF (npk == 0 ) ikmax = 1

  ! look for possible name for time dimension
  npt     = getdim(cf_in,cn_t, cdtrue=cv_tim, kstatus=ierr)
  IF ( ierr /= 0 ) THEN
     npt     = getdim(cf_in,'time', cdtrue=cv_tim, kstatus=ierr)
     IF ( ierr /= 0 ) THEN
        npt     = getdim(cf_in,'t', cdtrue=cv_tim, kstatus=ierr)
        IF ( ierr /= 0 ) THEN
           PRINT *, 'no time dimension found'
           npt=1
        ENDIF
     ENDIF
  ENDIF

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk ,' npkk  =', npkk
  PRINT *, 'npt    = ', npt

  IF (npkk > npk ) THEN
     PRINT *,' It seems that you want levels that are not represented '
     PRINT *,' in any of the variables that are in the file ',TRIM(cf_in)
     STOP
  ENDIF

  ALLOCATE( v2d(npiglo,npjglo),rlon(npiglo,npjglo), rlat(npiglo,npjglo), rdepg(npk) , rdep(npkk))
  ALLOCATE( zxz(npiglo,1), zyz(1,npjglo) )
  ALLOCATE( tim(npt) )

  nvars = getnvar(cf_in)
  PRINT *,' nvars =', nvars

  ALLOCATE (cv_names(nvars), ndim(nvars) )
  ALLOCATE (stypvar(nvars))
  ALLOCATE (id_var(nvars), ipk(nvars), id_varout(nvars), ipkk(nvars))

  rlon =getvar(cf_in, cn_vlon2d, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin)  ! nav_lon
  rlat =getvar(cf_in, cn_vlat2d, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin)  ! nav_lat

  IF ( npk /= 0 ) THEN
     rdepg   = getvar1d(cf_in, cv_dep, npk)
     rdep(:) = rdepg(ikmin:ikmax)  
  ENDIF

  CALL CreateOutput

  DO jvar = 1,nvars
     ! skip dimension variables (already done when creating the output file)
     ik1=MAX(1,ikmin) ; ik2=ipk(jvar)
     SELECT CASE (cv_names(jvar) )
        !
     CASE ('none' )
        ! skip
     CASE DEFAULT
        IF ( lzonal ) THEN
           ALLOCATE( v2dxz(npiglo,ipk(jvar)) )
           DO jt=1,npt
              v2dxz=getvarxz(cf_in, cv_names(jvar), ijmin, npiglo, ipk(jvar), kimin=iimin, kkmin=1, ktime=jt)
              DO jk=ik1,ik2
                 ik = jk - ik1 + 1 
                 zxz(:,1) = v2dxz(:,jk)
                 ierr=putvar(ncout, id_varout(jvar), zxz, ik, npiglo, 1, ktime=jt)
              ENDDO
           ENDDO
           DEALLOCATE ( v2dxz )
        ELSEIF (lmeridian) THEN
           ALLOCATE(  v2dyz(npjglo,ipk(jvar)) )
           DO jt=1,npt
              v2dyz=getvaryz(cf_in, cv_names(jvar), iimin, npjglo, ipk(jvar), kjmin=ijmin, kkmin=1, ktime=jt)
              DO jk=ik1, ik2
                 ik = jk - ik1 + 1 
                 zyz(1,:) = v2dyz(:,jk)
                 ierr=putvar(ncout, id_varout(jvar), zyz, ik, 1, npjglo, ktime=jt)
              ENDDO
           ENDDO
           DEALLOCATE ( v2dyz )
        ELSE
           DO jt = 1, npt
              DO jk=ik1,ik2
                 ik = jk - ik1 + 1
                 v2d  = getvar(cf_in, cv_names(jvar),  jk,      npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=jt)
                 ierr = putvar(ncout, id_varout(jvar), v2d, ik, npiglo, npjglo,                           ktime=jt)
              ENDDO
           ENDDO
        ENDIF
     END SELECT
  END DO ! loop to next var in file

  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ! get list of variable names and collect attributes in stypvar (optional)
    cv_names(:)=getvarname(cf_in, nvars, stypvar)

    ! save variable dimension in ndim
    !  1 = either time or depth : noclip
    !  2 = nav_lon, nav_lat
    !  3 = X,Y,T  or X,Y,Z   <-- need to fix the ambiguity ...
    !  4 = X,Y,Z,T
    DO jvar=1,nvars
       ndim(jvar) = getvdim(cf_in, cv_names(jvar)) + 1   !  we add 1 because vdim is dim - 1 ...
    END DO

    id_var(:)  = (/(jv, jv=1,nvars)/)

    ! ipk gives the number of level or 0 if not a T[Z]YX  variable
    ipk(:) = getipk (cf_in,nvars,cdep=cv_dep)
    ipk(:) = MIN ( ipk , ikmax )            ! reduce max depth to the required maximum
    ipkk(:)= MAX( 0 , ipk(:) - ikmin + 1 )  ! for output variable. For 2D input var, 
    ! ipkk is set to 0 if ikmin > 1 ... OK ? 
    WHERE( ipkk == 0 ) cv_names='none'
    stypvar(:)%cname = cv_names

    DO jv = 1, nvars
       stypvar(jv)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    ENDDO

    ! create output fileset
    ! create output file taking the sizes in cf_in
    ncout = create      (cf_out, cf_in,   npiglo, npjglo, npkk, cdep=cv_dep            , ld_nc4=lnc4)
    ierr  = createvar   (ncout,  stypvar, nvars,  ipkk,   id_varout, cdglobal=cglobal,   ld_nc4=lnc4)
    ierr  = putheadervar(ncout,  cf_in,   npiglo, npjglo, npkk, pnavlon=rlon, pnavlat=rlat, pdep=rdep, cdep=cv_dep)

    tim  = getvar1d(cf_in, cn_vtimec, npt     )
    ierr = putvar1d(ncout, tim,       npt, 'T')
  END SUBROUTINE CreateOutput

END PROGRAM cdfclip
