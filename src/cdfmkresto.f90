PROGRAM cdfmkresto
  !!======================================================================
  !!                     ***  PROGRAM  cdfmkresto  ***
  !!=====================================================================
  !!  ** Purpose : implement DRAKKAR restoring strategy off-line
  !!
  !!  ** Method  : port DRAKKAR resto_patch routine off-line
  !!
  !! History :  4.0  : 05/2017  : J.M. Molines :  Original code
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   PrintCfgInfo  : Give details on the format for input cfg file
  !!   resto_patch   : Drakkar/Nemo routine used for building a local patch
  !!   ReadCfg       : Read Config file to feed the patch structure
  !!   GetCoord      : Read Config file to feed the patch structure
  !!   CreateOutput  : Create output file
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE cdftools
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class preprocessing
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER                 :: wp=4
  INTEGER(KIND=4)                            :: jjpat, jk
  INTEGER(KIND=4)                            :: npiglo, npjglo, npk
  INTEGER(KIND=4)                            :: nvar=1
  INTEGER(KIND=4)                            :: ijarg, narg, iargc
  INTEGER(KIND=4)                            :: npatch
  INTEGER(KIND=4)                            :: ncout, ierr
  INTEGER(KIND=4), DIMENSION(1)              :: id_varout, ipk

  REAL(KIND=4)                               :: ra    = 6371229.   !: earth radius
  REAL(KIND=4)                               :: rad = 3.141592653589793 / 180.
  REAL(KIND=4)                               :: rvalue
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE:: gdept_1d
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE:: gphi, glam
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE:: resto

  REAL(KIND=8)                               :: dlon1, dlon2, dlat1, dlat2
  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE:: dtim

  CHARACTER(LEN=255)                         :: cf_coord  ! coordinate file
  CHARACTER(LEN=255)                         :: cf_cfg    ! config file (txt)
  CHARACTER(LEN=255)                         :: cf_out = 'damping_coef.nc' ! default output name
  CHARACTER(LEN=255)                         :: cf_dep    ! depth file (txt)
  CHARACTER(LEN=255)                         :: cf_resto  ! Previous resto file for appending new zones
  CHARACTER(LEN=255)                         :: cv_resto  ! Previous resto variable in cf_resto
  CHARACTER(LEN=255)                         :: cv_out = 'resto' ! output variable name
  CHARACTER(LEN=255)                         :: ctype='T' ! C-point type for the output file
  CHARACTER(LEN=255)                         :: cldum     ! dummy character variable

  TYPE                                       :: patch     ! Structure handling a patch
     CHARACTER(1) :: ctyp         ! type ( R or C)
     REAL(KIND=4) :: rlon1 
     REAL(KIND=4) :: rlon2 
     REAL(KIND=4) :: rlat1 
     REAL(KIND=4) :: rlat2 
     REAL(KIND=4) :: rim
     REAL(KIND=4) :: radius
     REAL(KIND=4) :: tresto
     REAL(KIND=4) :: rdep1
     REAL(KIND=4) :: rdep2
  END TYPE patch

  TYPE (patch )  , DIMENSION(:), ALLOCATABLE :: spatch
  TYPE (variable), DIMENSION(1)              :: stypvar            ! structure for attributes

  LOGICAL                                    :: lnc4      = .FALSE.     ! Use nc4 with chunking and deflation
  LOGICAL                                    :: lfdep     = .TRUE.      ! flag for ascii depth file
  LOGICAL                                    :: lprev     = .FALSE.     ! flag for previous restoring file
  LOGICAL                                    :: ltime     = .TRUE.      ! flag for specific value (not time)
  LOGICAL                                    :: l2d       = .FALSE.     ! flag for 2d fields
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg=iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfmkresto -c COORD-file -i CFG-file [-d DEP-file] [-o DMP-file]...'
     PRINT *,'                     ...[-ov VAR-out] [-2d] [-prev RESTO-file RESTO-var ] ...'
     PRINT *,'                     ...[-p C-TYPE] [-val VALUE] [-nc4] [-h]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Create a file with a 3D damping coefficient suitable for the NEMO'
     PRINT *,'       TS restoring performed by tradmp. Units for the damping_coefficient '
     PRINT *,'       are s^-1.' 
     PRINT *,'      '
     PRINT *,'       The restoring zone is defined by a series of patches defined in the'
     PRINT *,'       configuration file, and that can have either a rectangular or a circular'
     PRINT *,'       shape. Overlapping patches are not added.'
     PRINT *,'       '
     PRINT *,'       This tool has been improved with new options (-val, -2d) so that it can'
     PRINT *,'       be used for building ''bfr2d_coef'' as well as ''shlat2d''.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -c COORD-file : pass the name of the file with horizontal coordinates.'
     PRINT *,'       -i CFG-file : pass the name of an ascii configuration file used to '
     PRINT *,'               to define the restoring properties. (use cdfmkresto -h for '
     PRINT *,'               a detailed description of this CFG-file, with some examples.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-h ]: print a detailed description of the configuration file.'
     PRINT *,'       [-d DEP-file]: name on an ASCII file with gdept_1d, if ',TRIM(cn_fzgr)
     PRINT *,'                is not available. This file is just a list of the deptht, in'
     PRINT *,'                one column.'
     PRINT *,'       [-prev RESTO-file RESTO-var] : Use RESTO-file and RESTO-var for the'
     PRINT *,'                initialization of the restoring coefficient. Units MUST be'
     PRINT *,'                s^-1 !!! '
     PRINT *,'       [-p C-type] : indicate on which grid point (T or F -so far-) the '
     PRINT *,'               variable is computed in the output file.'
     PRINT *,'       [-val VALUE ] : with this option, the ''restoring'' coefficient is'
     PRINT *,'                  set to VALUE , instead of a time scale.'
     PRINT *,'       [-2d ] : Create a 2D file instead of a default 3D file.'
     PRINT *,'       [-o DMP-file]: name of the output file instead of ',TRIM(cf_out),'.'
     PRINT *,'       [-ov VAR-out]: name of the output variable instead of ',TRIM(cv_out),'.'
     PRINT *,'       [-nc4]  : Use netcdf4 output with chunking and deflation level 1'
     PRINT *,'            This option is effective only if cdftools are compiled with'
     PRINT *,'            a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'         ', TRIM(cn_fzgr),'. If not available, use the -d option.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option is used.'
     PRINT *,'         variables : ', TRIM(cv_out),' (s^-1 ) (unless -ov option used)'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'         cdfmkdmp'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-c'   ) ; CALL getarg(ijarg, cf_coord ) ; ijarg=ijarg+1
     CASE ( '-i'   ) ; CALL getarg(ijarg, cf_cfg   ) ; ijarg=ijarg+1
        ! option
     CASE ( '-h'   ) ; CALL PrintCfgInfo             ; STOP 0
     CASE ( '-d'   ) ; CALL getarg(ijarg, cf_dep   ) ; ijarg=ijarg+1
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out   ) ; ijarg=ijarg+1
     CASE ( '-ov'  ) ; CALL getarg(ijarg, cv_out   ) ; ijarg=ijarg+1
     CASE ( '-prev') ; CALL getarg(ijarg, cf_resto ) ; ijarg=ijarg+1
        ;            ; CALL getarg(ijarg, cv_resto ) ; ijarg=ijarg+1
        ;            ; lprev = .TRUE.
     CASE ( '-p'   ) ; CALL getarg(ijarg, ctype    ) ; ijarg=ijarg+1 
     CASE ( '-val' ) ; CALL getarg(ijarg, cldum    ) ; ijarg=ijarg+1 ; READ(cldum,*) rvalue
        ;            ; ltime = .FALSE.
     CASE ( '-2d'  ) ; l2d  = .TRUE.
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 99 
     END SELECT
  ENDDO

  ! Sanity check
  IF ( lprev ) THEN
        IF ( chkfile(cf_resto) ) STOP 99
  ENDIF
  IF ( .NOT. l2d ) THEN
     IF ( chkfile(cf_dep, ld_verbose=.FALSE.) ) THEN
        ! look for cn_fzgr file
        lfdep=.FALSE.
        IF ( chkfile(cn_fzgr) ) STOP 99
     ENDIF
  ENDIF
  IF ( chkfile(cf_coord) .OR. chkfile(cf_cfg) ) STOP 99 ! missing file
  IF ( ctype /= 'T' .AND. ctype /='F' ) THEN
     PRINT *,' ERROR : C-TYPE can be only T or F (so far...)'
     STOP 99
  ENDIF

  CALL ReadCfg  ! read configuration file and set variables

  CALL GetCoord ! read model horizontal coordinates and vertical levels
  ALLOCATE ( resto(npiglo, npjglo, npk) )

  IF ( lprev ) THEN 
    DO jk=1,npk
      resto(:,:,jk) = getvar(cf_resto,cv_resto, jk, npiglo, npjglo) 
    ENDDO
  ELSE
    resto(:,:,:) = 0.
  ENDIF

  CALL CreateOutput ! prepare netcdf output file 

  DO jjpat =1, npatch
     CALL resto_patch ( spatch(jjpat), resto  )
  ENDDO
  DO jk = 1, npk
     ierr = putvar( ncout, id_varout(1), resto(:,:,jk), jk, npiglo, npjglo)
  ENDDO
  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE PrintCfgInfo
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE PrintCfgInfo  ***
    !!
    !! ** Purpose :  Display detailed information about the configuration
    !!               file.
    !!
    !! ** Method  :  Just print text.
    !!
    !!----------------------------------------------------------------------
    PRINT *,'       '
    PRINT *,'#   FORMAT OF THE CONFIGURATION FILE FOR cdfmkresto TOOL.'
    PRINT *,'       '
    PRINT *,'## CONTEXT:      '
    PRINT *,'      The restoring zone is defined by a series of patches described in the'
    PRINT *,'    configuration file, and that can have either a rectangular, a circular'
    PRINT *,'    or a disk shape. Overlapping patches uses the maximum between the patches.'
    PRINT *,'       '
    PRINT *,'    The configuration file is used to describe the series of patches.'
    PRINT *,'       '
    PRINT *,'## FILE FORMAT:   '
    PRINT *,'      There are as many lines as patches in the configuration files. No blank'
    PRINT *,'    lines are allowed but lines starting with a # are skipped.'
    PRINT *,'      Each line starts with either R, C, D or I to indicate either a Rectangular'
    PRINT *,'    patch, a Circular patch a Disk patch or IJ patch. Remaining fields on the '
    PRINT *,'    line  depend on the type of patch:'
    PRINT *,'###   Case of rectangular patches: the line looks like:'
    PRINT *,'       R lon1 lon2 lat1 lat2 band_width tresto z1 z2'
    PRINT *,'     In the rectangular case, the rectangle is defined by its geographical '
    PRINT *,'     limits (lon1, lon2, lat1, lat2 in degrees).'
    PRINT *,'     tresto represents the restoring time scale for the maximum restoring.'
    PRINT *,'     The restoring coefficient decay to zero outside the patch. The decay is'
    PRINT *,'     applied across a rim which width is given by the rim_width parameter'
    PRINT *,'     given in degrees.'
    PRINT *,'     Additionaly, z1 and z2 indicates a depth range (m) limiting the restoring.'
    PRINT *,'     If z1=z2, the restoring will be applied on the full water column.'
    PRINT *,'       '
    PRINT *,'###   Case of circulat patches: the line looks like:'
    PRINT *,'       C lon1 lat1 radius tresto z1 z2'
    PRINT *,'     In the circular case, the circle is defined by its center (lon1, lat1) in '
    PRINT *,'     degrees  and its radius in km. ''tresto'' represents the restoring time'
    PRINT *,'     scale at the center of the circle. The damping coefficient decays to zero'
    PRINT *,'     outside the defined circle. As in the rectangular case, z1 and z2 define a'
    PRINT *,'     depth range for  limiting the restoring. If z1=z2, the restoring is '
    PRINT *,'     applied to the full water column.'
    PRINT *,'       '
    PRINT *,'###   Case of disk patches: the line looks like:'
    PRINT *,'       D lon1 lat1 radius rim tresto z1 z2'
    PRINT *,'      In the disk case, a circle is defined with its center (lon1, lat1) in '
    PRINT *,'      degrees, and its radius in km.  In the disk case, the damping coefficient'
    PRINT *,'      is constant over the disk, and a there is a linear decay in a ring '
    PRINT *,'      around the disk, which width is the rim value (in km).'
    PRINT *,'      '
    PRINT *,'###   Case of a IJ patch : In this case, a rectangular patch is constructed '
    PRINT *,'      using directly the I,J coordinates (config dependent) passed on the line.'
    PRINT *,'      line looks like:'
    PRINT *,'      I  imin  imax  jmin  jmax  tresto  z1 z2 '
    PRINT *,'      In this case, if -val option is used, the value is exactly tresto'
    PRINT *,'       '
    PRINT *,'##  EXAMPLE:  '
    PRINT *,'       The standard DRAKKAR restoring procedure corresponds to the following '
    PRINT *,'     configuration file:  '
    PRINT *,'    '
    PRINT *,'# DRAKKAR restoring configuration file'
    PRINT *,'# Black Sea'
    PRINT *,'# type lon1  lon2 lat1 lat2 rim_width tresto   z1   z2'
    PRINT *,'    R   27.4 42.0 41.0 47.5   0.       180.     0    0'
    PRINT *,'# Red Sea'
    PRINT *,'    R   20.4 43.6 12.9 30.3   0.       180.     0    0'
    PRINT *,'# Persian Gulf'
    PRINT *,'    R   46.5 57.0 23.0 31.5   1.       180.     0    0'
    PRINT *,'# Restoring in the overflow regions'
    PRINT *,'# Gulf of Cadix (Gibraltar overflow)'
    PRINT *,'# type lon1  lat1 radius tresto  z1    z2'
    PRINT *,'    C  -7.0  36.0  80.      6.  600. 1300.'
    PRINT *,'# Gulf of Aden (Bab-el-Mandeb overflow)'
    PRINT *,'    C  44.75 11.5 100.      6.  0.   0.'
    PRINT *,'# Arabian Gulf  (Ormuz Strait  overflow)'
    PRINT *,'    C  57.75 25.0 100.      6.  0.   0.'
    PRINT *,'# Disk example ( not used in drakkar so far)'
    PRINT *,'# type lon1  lat1 radius  rim tresto  z1    z2'
    PRINT *,'   D    -30  -50   200    10   720    0     500'
    PRINT *,'# IJ example ( dummy , just an example)'
    PRINT *,'   I    250  252   300 301  10.  0 500 '
    PRINT *,'    '
    PRINT *,'       This file defines 8 patches, (3 rectangular, 3 circular,1 disk and 1 IJ).'
    PRINT *,'     Lines starting by # are helpfull comment for documenting the'
    PRINT *,'     restoring strategy. Note that as far as the position of the patches are '
    PRINT *,'     given in geographical coordinates, the same file can be used for different'
    PRINT *,'     model configuration (and resolution) !'
    PRINT *,'       Note that the time scale is not used if a value is passed with the '
    PRINT *,'     -val option. This latter option, combined with -2d option can be used'
    PRINT *,'     for producing 2D bottom friction enhancement, of 2D shlat coefficient.'
    PRINT *,'    '


  END SUBROUTINE PrintCfgInfo

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ipk(1)                    = npk
    stypvar(1)%ichunk         = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname          = cv_out
    IF ( ltime ) THEN
      stypvar(1)%cunits         = '[s^1]'
    ELSE
      stypvar(1)%cunits         = '[ ]'
    ENDIF
    stypvar(1)%rmissing_value = 1.e+20
    IF ( l2d ) THEN
      stypvar(1)%caxis          = 'TYX'
    ELSE
      stypvar(1)%caxis          = 'TZYX'
    ENDIF
    stypvar(1)%valid_min      = 0.
    stypvar(1)%valid_max      = 500.
    IF ( ltime )  THEN
      stypvar(1)%clong_name     = 'Restoring coefficent'
    ELSE
      stypvar(1)%clong_name     = 'Mask coefficent'
    ENDIF
    stypvar(1)%cshort_name    = cv_out

    IF ( l2d ) THEN
      ncout = create      (cf_out, 'none',  npiglo, npjglo, 0, ld_nc4=lnc4 )
      ierr  = createvar   (ncout,  stypvar,  1,     ipk,        id_varout,       ld_nc4=lnc4 )
      ierr  = putheadervar(ncout,  cf_coord,  npiglo, npjglo, 0 ,  &
           pnavlon=glam, pnavlat=gphi)
    ELSE
      ncout = create      (cf_out, 'none',  npiglo, npjglo, npk  ,cdep='deptht', ld_nc4=lnc4 )
      ierr  = createvar   (ncout,  stypvar,  1,     ipk,        id_varout,       ld_nc4=lnc4 )
      ierr  = putheadervar(ncout,  cf_coord,  npiglo, npjglo, npk ,  &
           pnavlon=glam, pnavlat=gphi, pdep=gdept_1d, cdep=cn_vdeptht )
    ENDIF

    ALLOCATE (dtim(1) )
    dtim = 0.d0
    ierr = putvar1d(ncout, dtim, 1   , 'T'  )

  END SUBROUTINE CreateOutput

  SUBROUTINE ReadCfg
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ReadCfg  ***
    !!
    !! ** Purpose :  Read config file and set patch parameters
    !!
    !! ** Method  :  read and fill in corresponding variables
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4)    :: inum = 10! logical unit for configuration file
    INTEGER(KIND=4)    :: ierr ! error flag
    CHARACTER(LEN=255) :: cline, cltyp

    OPEN(inum, FILE=cf_cfg)
    ! count the number of patches and allocate spatch
    ierr=0
    npatch=0
    DO WHILE ( ierr == 0 )
       READ(inum,'(a)', iostat=ierr) cline
       IF ( ierr == 0 ) THEN
          IF ( cline(1:1) /= '#' ) THEN
             npatch=npatch+1
          ENDIF
       ENDIF
    ENDDO
    PRINT *,' NPATCH = ', npatch
    ALLOCATE ( spatch(npatch) )

    ! read again to fill the patch structure
    REWIND(inum)
    ierr = 0
    npatch = 0
    DO WHILE ( ierr == 0 )
       READ(inum,'(a)', iostat=ierr) cline
       IF ( ierr == 0 ) THEN
          IF ( cline(1:1) /= '#' ) THEN
             npatch=npatch+1
             PRINT *,'   Patch ',npatch,' :  ',TRIM(cline)
             READ(cline,*) cltyp
             SELECT CASE (cltyp)
             CASE ( 'R', 'r' ) ; READ(cline,*) spatch(npatch)%ctyp, &  ! Rectangle
                  & spatch(npatch)%rlon1, spatch(npatch)%rlon2,     &
                  & spatch(npatch)%rlat1, spatch(npatch)%rlat2,     &
                  & spatch(npatch)%rim, spatch(npatch)%tresto,      &
                  & spatch(npatch)%rdep1, spatch(npatch)%rdep2 
             CASE ( 'C', 'c' ) ; READ(cline,*) spatch(npatch)%ctyp, &   ! Circle
                  & spatch(npatch)%rlon1, spatch(npatch)%rlat1,     &
                  & spatch(npatch)%radius, spatch(npatch)%tresto,   &
                  & spatch(npatch)%rdep1, spatch(npatch)%rdep2 
             CASE ( 'D', 'd' ) ; READ(cline,*) spatch(npatch)%ctyp, &   ! Disque
                  & spatch(npatch)%rlon1, spatch(npatch)%rlat1,     &
                  & spatch(npatch)%radius,     &
                  & spatch(npatch)%rim, spatch(npatch)%tresto,   &
                  & spatch(npatch)%rdep1, spatch(npatch)%rdep2 
             CASE ( 'I', 'i' ) ; READ(cline,*) spatch(npatch)%ctyp, &   ! I
                  & spatch(npatch)%rlon1, spatch(npatch)%rlon2,     &
                  & spatch(npatch)%rlat1, spatch(npatch)%rlat2,     &
                  & spatch(npatch)%tresto,                          &
                  & spatch(npatch)%rdep1, spatch(npatch)%rdep2 
                  PRINT *, ' W A R N I N G : You are using a patch defined with I,J index'
                  PRINT *, ' =============   It depends on your configuration !'
             END SELECT
          ENDIF
       ENDIF
    ENDDO


  END SUBROUTINE ReadCfg

  SUBROUTINE GetCoord
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetCoord  ***
    !!
    !! ** Purpose :  Read glam, gphi from horizontal grid information
    !!               Read vertical grid information  
    !!
    !! ** Method  :  Open cf_coord file for horizontal grid.
    !!               Read vertical information from eiher mesh_zgr.nc
    !!               or cf_gdep ascii file
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4) :: inum=11, ieof=0, ilev=0
    !!----------------------------------------------------------------------
    npiglo = getdim(cf_coord,cn_x)
    npjglo = getdim(cf_coord,cn_y)
    ALLOCATE (glam(npiglo,npjglo), gphi(npiglo,npjglo) )
    IF ( ctype == 'T' ) THEN
      glam(:,:)=getvar(cf_coord,cn_glamt, 1,npiglo,npjglo)
      gphi(:,:)=getvar(cf_coord,cn_gphit, 1,npiglo,npjglo)
    ELSE IF  ( ctype == 'F' ) THEN
      glam(:,:)=getvar(cf_coord,cn_glamf, 1,npiglo,npjglo)
      gphi(:,:)=getvar(cf_coord,cn_gphif, 1,npiglo,npjglo)
    ELSE
      PRINT *, "ERROR : C-Type can be only T or F" ; STOP 99
    ENDIF
    ! now deal with vertical levels ( suppose z or zps ! )
    IF ( l2d ) THEN
       npk=1
       ALLOCATE(gdept_1d(npk))
       gdept_1d(1) = 0.
    ELSE
      IF ( lfdep ) THEN
         OPEN(inum, FILE=cf_dep)
         ! first read to look for number of levels
         DO WHILE ( ieof == 0 )
            READ(inum,*,iostat=ieof)
            ilev=ilev+1
         ENDDO
         npk=ilev - 1
         ALLOCATE(gdept_1d(npk))
         REWIND(inum)
         DO jk=1,npk
            READ(inum,*) gdept_1d(jk)
         ENDDO
         CLOSE(inum)
      ELSE
         npk = getdim(cn_fzgr,'z', kstatus=ierr)   ! depth dimension in mesh_zgr is 'z' 
         IF ( ierr /= 0 ) THEN
            npk   = getdim (cn_fzgr,'nav_lev', kstatus=ierr)
         ENDIF
  
         ALLOCATE( gdept_1d(npk) )
         gdept_1d(:) =  getvare3(cn_fzgr, cn_gdept, npk)
      ENDIF
    ENDIF

  END SUBROUTINE GetCoord


  SUBROUTINE resto_patch ( sd_patch, presto )
    !!------------------------------------------------------------------------
    !!                 ***  Routine resto_patch  ***
    !!
    !! ** Purpose :   modify resto array on a geographically defined zone.
    !!                The restoring can be limited on the vertical.
    !!
    !! ** Method  :  Use glam, gphi arrays. If the defined zone is outside 
    !!              the domain, resto is unchanged. Reactangular and Circulat
    !!              patches can be used.
    !!
    !!------------------------------------------------------------------------
    TYPE (patch),               INTENT(in )   :: sd_patch
    REAL(wp), DIMENSION(:,:,:), INTENT(inout) :: presto 
    !
    REAL(wp)          :: zlon1, zlon2, zlat1, zlat2, zbw, ztmax
    REAL(wp)          :: zz1, zz2
    !!
    INTEGER :: ji,jj, jk            ! dummy loop index
    INTEGER :: ii1, ii2, ij1, ij2   ! limiting horizontal index corresponding (I case)
    INTEGER :: ik1, ik2             ! limiting vertical index corresponding to zz1,zz2
    INTEGER, DIMENSION(1)  :: iloc 

    REAL(wp) :: zv1, zv2, zv3, zv4, zcoef, ztmp, zdist, zradius, zradius2, zcoef2
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zpatch
    REAL(wp), DIMENSION(:)  , ALLOCATABLE :: zmask

    CHARACTER(LEN=1) :: cl_typ
    !!------------------------------------------------------------------------
    ALLOCATE ( zpatch(npiglo,npjglo), zmask(npk) )

    cl_typ   = sd_patch%ctyp
    zlon1    = sd_patch%rlon1
    zlon2    = sd_patch%rlon2
    zlat1    = sd_patch%rlat1
    zlat2    = sd_patch%rlat2
    zbw      = sd_patch%rim
    zradius  = sd_patch%radius
    zradius2 = zradius * zradius
    ztmax    = sd_patch%tresto
    zz1      = sd_patch%rdep1
    zz2      = sd_patch%rdep2

    zpatch = 0._wp
    IF ( ltime ) THEN
       zcoef  = 1._wp/ztmax/86400._wp
    ELSE
       zcoef  = rvalue
    ENDIF

    SELECT CASE ( cl_typ )
    CASE ( 'C', 'c' )   ! Circular patch 
       !  mask for horizontal extent
       DO jj = 1, npjglo
          DO ji = 1 , npiglo
             zpatch(ji,jj) =  SIN(gphi(ji,jj)*rad)* SIN(zlat1*rad)  &
                  &         + COS(gphi(ji,jj)*rad)* COS(zlat1*rad)  &
                  &         * COS(rad*(zlon1-glam(ji,jj)))
          ENDDO
       ENDDO

       WHERE ( ABS (zpatch ) > 1 ) zpatch = 1.
       ! applying spatial horizontal variation
       DO jj = 1, npjglo
          DO ji= 1, npiglo 
             ztmp = zpatch(ji,jj)
             zdist = ATAN(SQRT( (1.-ztmp)/(1+ztmp)) )*2.*ra/1000.
             zpatch(ji,jj) = EXP( - zdist*zdist/zradius2 )
          ENDDO
       ENDDO
       ! clean cut off
       WHERE (ABS(zpatch) < 0.01 ) zpatch = 0.

    CASE ( 'D', 'd' )   ! Circular patch 
       !  mask for horizontal extent
       DO jj = 1, npjglo
          DO ji = 1 , npiglo
            dlon1=zlon1*1.d0         ; dlat1=zlat1*1.d0
            dlon2=glam(ji,jj)*1.d0  ; dlat2=gphi(ji,jj)*1.d0
            zpatch(ji,jj) =  dist( dlon1, dlon2, dlat1, dlat2 )
          ENDDO
       ENDDO

       WHERE ( ABS (zpatch ) < zradius     ) zpatch = 1.
       WHERE ( ABS (zpatch ) > zradius+zbw ) zpatch = 0.
       ! applying spatial horizontal variation
       DO jj = 1, npjglo
          DO ji= 1, npiglo
           IF (zpatch(ji,jj) >= zradius ) THEN
              zpatch(ji,jj) = (zradius+zbw - zpatch(ji,jj) )/ zbw
           ENDIF
          ENDDO
       ENDDO
       ! clean cut off
       WHERE (ABS(zpatch) < 0.01 ) zpatch = 0.

       ! JMM : eventually add some checking to avoid locally large resto.


    CASE ( 'R','r' )
       ! horizontal extent
       zcoef2=1./(zbw +1.e-20 ) ! to avoid division by 0
       DO jj=1,npjglo
          DO ji=1,npiglo
             zv1=MAX(0., zcoef2*( glam(ji,jj) - zlon1)  )
             zv2=MAX(0., zcoef2*( zlon2 - glam(ji,jj))  )
             zv3=MAX(0., zcoef2*( gphi(ji,jj) - zlat1)  )
             zv4=MAX(0., zcoef2*( zlat2 - gphi(ji,jj))  )
             zpatch(ji,jj)= MIN( 1., MIN( 1., zv1,zv2,zv3,zv4 ) )
          ENDDO
       ENDDO
     CASE ( 'I','i')
       ii1=NINT(zlon1)  ; ii2 = NINT(zlon2)
       ij1=NINT(zlat1)  ; ij2 = NINT(zlat2)
       DO jj= ij1, ij2
          DO ji=ii1, ii2
             PRINT *, ' I-Patch at ',ji,jj, ztmax
             zpatch(ji,jj) = 1.
          ENDDO
       ENDDO
       ! Alter zcoef 
       IF ( .NOT. ltime  ) zcoef = ztmax
    END SELECT

    ! Vertical limitation same treatment for both types
    zmask(:) = 1.
    IF ( zz1 /= zz2 ) THEN
       WHERE ( gdept_1d < zz1 .OR. gdept_1d > zz2 ) zmask = 0.
       ! look for first 1
       iloc=MAXLOC(zmask) ; ik1 = iloc(1)
       ! now look for first 0
       zmask(1:ik1) = 1.
       iloc=MINLOC(zmask) ; ik2 = iloc(1) - 1

       zmask = 0._wp
       IF (ik2-ik1 > 4 ) THEN ! vertical ramp
          zmask(ik1        ) = 0.25_wp
          zmask(ik1+1      ) = 0.75_wp
          zmask(ik1+2:ik2-2) = 1.0_wp
          zmask(ik2-1      ) = 0.75_wp
          zmask(ik2        ) = 0.25_wp
       ELSE
          zmask(ik1:ik2) = 1. ! all the water column in zz1-zz2 is restored the same
       ENDIF
    ENDIF

    DO jk=1, npk
       presto(:,:,jk)= MAX(presto(:,:,jk),zpatch * zcoef * zmask(jk) )
    ENDDO
    !
    DEALLOCATE ( zmask, zpatch )
    !
  END SUBROUTINE resto_patch
END PROGRAM cdfmkresto
