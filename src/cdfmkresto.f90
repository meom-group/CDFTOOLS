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

  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE:: resto
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: gphit, glamt
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE    :: gdept_1d, tim
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE    :: rlon1, rlon2, rlat1, rlat2
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE    :: rbw, rtmax, rz1, rz2
  REAL(KIND=4)                               :: ra    = 6371229.   !: earth radius
  REAL(KIND=4)                               :: rad = 3.141592653589793 / 180.

  CHARACTER(LEN=255)                         :: cf_coord
  CHARACTER(LEN=255)                         :: cf_cfg
  CHARACTER(LEN=255)                         :: cf_out = 'damping_coef.nc'
  CHARACTER(LEN=255)                         :: cf_dep 
  CHARACTER(LEN=255)                         :: cv_out = 'resto'
  CHARACTER(LEN=255)                         :: cldum

  TYPE                                       :: patch
    CHARACTER(1) :: ctyp
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
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg=iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfmkresto -c COORD-file -i CFG-file [-d DEP-file] [-o DMP-file]...'
     PRINT *,'                  ...  [-nc4] [-h]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Create a file with a 3D damping coefficient suitable for the NEMO'
     PRINT *,'       TS restoring performed by tradmp. Units for the damping_coefficient '
     PRINT *,'       is days^-1.' 
     PRINT *,'      '
     PRINT *,'       The restoring zone is defined by a series of patches defined in the'
     PRINT *,'       configuration file, and that can have either a rectangular or a circular'
     PRINT *,'       shape. Overlapping patches are not added.'
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
     PRINT *,'       [-o DMP-file]: name of the output file instead of ',TRIM(cf_out),'.'
     PRINT *,'       [ -nc4] : Use netcdf4 output with chunking and deflation level 1'
     PRINT *,'            This option is effective only if cdftools are compiled with'
     PRINT *,'            a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'         ', TRIM(cn_fzgr),'. If not available, use the -d option.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option is used.'
     PRINT *,'         variables : ', TRIM(cv_out),' (days^-1 )'
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
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  IF ( chkfile(cf_dep, ld_verbose=.FALSE.) ) THEN
     ! look for cn_fzgr file
     lfdep=.FALSE.
     IF ( chkfile(cn_fzgr) ) STOP
  ENDIF
  IF ( chkfile(cf_coord) .OR. chkfile(cf_cfg) ) STOP ! missing file

  CALL ReadCfg  ! read configuration file and set variables
  
  CALL GetCoord ! read model horizontal coordinates and vertical levels
  ALLOCATE ( resto(npiglo, npjglo, npk) )
  resto(:,:,:) = 0.

  CALL CreateOutput ! prepare netcdf output file 

  DO jjpat =1, npatch
!     CALL resto_patch ( rlon1(jjpat), rlon2(jjpat), rlat1(jjpat), rlat2(jjpat), &
!          &             rbw(jjpat), rtmax(jjpat), resto, rz1(jjpat), rz2(jjpat) )
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
    PRINT *,'    configuration file, and that can have either a rectangular or a circular'
    PRINT *,'    shape. Overlapping patches uses the maximum between the patches.'
    PRINT *,'       '
    PRINT *,'    The configuration file is used to describe the series of patches.'
    PRINT *,'       '
    PRINT *,'## FILE FORMAT:   '
    PRINT *,'      There are as many lines as patches in the configuration files. No blank'
    PRINT *,'    lines are allowed but lines starting with a # are skipped.'
    PRINT *,'      Each line starts with either R or C to indicate either a Rectangular'
    PRINT *,'    patch or a Circular patch. Remaining fields on the line depend on the type'
    PRINT *,'    of patch:'
    PRINT *,'###   Case of rectangular patches: the line looks like:'
    PRINT *,'       R lon1 lon2 lat1 lat2 band_width tresto z1 z2'
    PRINT *,'     In the rectangular case, the rectangle is defined by its geographical '
    PRINT *,'     limits (lon1, lon2, lat1, lat2 in degrees).'
    PRINT *,'     tresto represents the restoring time scale for the maximum restoring.'
    PRINT *,'     The restoring coefficient decay to zero outside the patch. The decay is'
    PRINT *,'     applied across a rim which width is given by the rim_width parameter'
    PRINT *,'     given in kilometers.'
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
    PRINT *,'    '
    PRINT *,'       This file defines 6 patches, (3 rectangular and 3 circular).'
    PRINT *,'     Lines starting by # are helpfull comment for documenting the'
    PRINT *,'     restoring strategy. Note that as far as the position of the patches are '
    PRINT *,'     given in geographical coordinates, the same file can be used for different'
    PRINT *,'     model configuration (and resolution) !'
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
    stypvar(1)%cunits         = '[days^1]'
    stypvar(1)%rmissing_value = 1.e+20
    stypvar(1)%caxis          = 'TZYX'
    stypvar(1)%valid_min      = 0.
    stypvar(1)%valid_max      = 500.
    stypvar(1)%clong_name     = 'Restoring coefficent'
    stypvar(1)%cshort_name    = cv_out

    ncout = create      (cf_out, 'none',  npiglo, npjglo, npk  ,cdep='deptht', ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar,  1,     ipk,        id_varout,       ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_coord,  npiglo, npjglo, npk ,  &
                              pnavlon=glamt, pnavlat=gphit, pdep=gdept_1d, cdep=cn_vdeptht )

    ALLOCATE (tim(1) )
    tim(1)= 0.
    ierr   = putvar1d(ncout, tim, 1   , 'T'  )

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
          CASE ( 'R', 'r' ) ; READ(cline,*) spatch(npatch)%ctyp ,                &
                                & spatch(npatch)%rlon1, spatch(npatch)%rlon2,    &
                                & spatch(npatch)%rlat1, spatch(npatch)%rlat2,    &
                                & spatch(npatch)%rim, spatch(npatch)%tresto,     &
                                & spatch(npatch)%rdep1, spatch(npatch)%rdep2 
          CASE ( 'C', 'c' ) ; READ(cline,*) spatch(npatch)%ctyp ,                &
                                & spatch(npatch)%rlon1, spatch(npatch)%rlat1,    &
                                & spatch(npatch)%radius, spatch(npatch)%tresto,  &
                                & spatch(npatch)%rdep1, spatch(npatch)%rdep2 
          END SELECT
        ENDIF
      ENDIF
    ENDDO


  END SUBROUTINE ReadCfg
  
  SUBROUTINE GetCoord
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetCoord  ***
    !!
    !! ** Purpose :  Read glamt, gphit from horizontal grid information
    !!               Read vertical grid information  
    !!
    !! ** Method  :  Open cf_coord file for horizontal grid.
    !!               Read vertical information from eiher mesh_zgr.nc
    !!               or cf_gdep ascii file
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4) :: inum=11, ieof=0, ilev=0
    
     npiglo = getdim(cf_coord,cn_x)
     npjglo = getdim(cf_coord,cn_y)
     ALLOCATE (glamt(npiglo,npjglo), gphit(npiglo,npjglo) )
     glamt(:,:)=getvar(cf_coord,cn_glamt, 1,npiglo,npjglo)
     gphit(:,:)=getvar(cf_coord,cn_gphit, 1,npiglo,npjglo)
     ! now deal with vertical levels ( suppose z or zps ! )
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
       npk = getdim(cn_fzgr,'z')   ! depth dimension in mesh_zgr is 'z' 
       ALLOCATE( gdept_1d(npk) )
       gdept_1d(:) =  getvare3(cn_fzgr, cn_gdept, npk)
     ENDIF

  END SUBROUTINE GetCoord
  

!   SUBROUTINE resto_patch ( plon1, plon2, plat1, plat2, pbw ,ptmax, presto, pz1, pz2 )
   SUBROUTINE resto_patch ( sd_patch, presto )
      !!------------------------------------------------------------------------
      !!                 ***  Routine resto_patch  ***
      !!
      !! ** Purpose :   modify resto array on a geographically defined zone.
      !!
      !! ** Method  :  Use glamt, gphit arrays. If the defined zone is outside 
      !!              the domain, resto is unchanged. If pz1 and pz2 are provided
      !!              then plon1, plat1 is taken as the position of a the center
      !!              of a circle with decay radius is pbw (in km) 
      !!
      !! ** Action  : IF not present pz1, pz2 : 
      !!              - plon1, plon2 : min and max longitude of the zone (Deg)
      !!              - plat1, plat2 : min and max latitude of the zone (Deg)
      !!              - pbw : band width of the linear decaying restoring (Deg)
      !!              - ptmax : restoring time scale for the inner zone (days)
      !!              IF present pz1 pz2
      !!              - plon1, plat1 : position of the center of the circle
      !!              - pbw = radius (km) of the restoring circle
      !!              - ptmax = time scale at maximum restoring
      !!              - pz1, pz2 : optional: if used, define the depth range (m)
      !!                          for restoring. If not all depths are considered
      !!------------------------------------------------------------------------
      TYPE (patch),               INTENT(in )   :: sd_patch
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) :: presto 
      !
      REAL(wp)          :: plon1, plon2, plat1, plat2, pbw, ptmax
      REAL(wp)          :: pz1, pz2
      !!
      INTEGER :: ji,jj, jk    ! dummy loop index
      INTEGER :: ik1, ik2     ! limiting vertical index corresponding to pz1,pz2
      INTEGER :: ij0, ij1, iiO, ii1
      INTEGER, DIMENSION(1)           :: iloc 

      REAL(wp) :: zv1, zv2, zv3, zv4, zcoef, ztmp, zdist, zradius2, zcoef2
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zpatch
      REAL(wp), DIMENSION(:)  , ALLOCATABLE :: zmask

      CHARACTER(LEN=1) :: cl_typ
      !!------------------------------------------------------------------------
      ALLOCATE ( zpatch(npiglo,npjglo), zmask(npk) )

      cl_typ = sd_patch%ctyp
      plon1=sd_patch%rlon1
      plon2=sd_patch%rlon2
      plat1=sd_patch%rlat1
      plat2=sd_patch%rlat2
      pbw  =sd_patch%rim
      zradius2  =sd_patch%radius*sd_patch%radius
      ptmax=sd_patch%tresto
      pz1  =sd_patch%rdep1
      pz2  =sd_patch%rdep2

      zpatch = 0._wp
      zcoef  = 1._wp/ptmax/86400._wp

      SELECT CASE ( cl_typ )
      CASE ( 'C', 'c' ) 
        ! horizontal extent
!        zradius2 = pbw * pbw !  radius squared
        DO jj = 1, npjglo
          DO ji = 1 , npiglo
            zpatch(ji,jj) =  sin(gphit(ji,jj)*rad)* sin(plat1*rad)  &
       &                         + cos(gphit(ji,jj)*rad)* cos(plat1*rad)  &
       &                         * cos(rad*(plon1-glamt(ji,jj)))
          ENDDO
        ENDDO

        WHERE ( abs (zpatch ) > 1 ) zpatch = 1.
        DO jj = 1, npjglo
          DO ji= 1, npiglo 
             ztmp = zpatch(ji,jj)
             zdist = atan(sqrt( (1.-ztmp)/(1+ztmp)) )*2.*ra/1000.
             zpatch(ji,jj) = exp( - zdist*zdist/zradius2 )
          ENDDO
        ENDDO
        ! clean cut off
        WHERE (ABS(zpatch) < 0.01 ) zpatch = 0.
        ! Vertical limitation
        zmask = 1.
        WHERE ( gdept_1d < pz1 .OR. gdept_1d > pz2 ) zmask = 0.
        iloc=MAXLOC(zmask) ; ik1 = iloc(1)
        zmask(1:ik1) = 1.
        iloc=MAXLOC(zmask) ; ik2 = iloc(1) - 1
        zmask(1:ik1) = 1.
        iloc=MAXLOC(zmask) ; ik2 = iloc(1) - 1
        IF (ik2 > 2 ) THEN
          zmask = 0._wp
          zmask(ik1       ) = 0.25_wp
          zmask(ik1+1     ) = 0.75_wp
          zmask(ik1+2:ik2-2) = 1.0_wp
          zmask(ik2-1     ) = 0.75_wp
          zmask(ik2       ) = 0.25_wp
        ELSE
          zmask = 1.   ! all the water column is restored the same
        ENDIF

        DO jk=1, npk
          presto(:,:,jk)= presto(:,:,jk) + zpatch * zcoef * zmask(jk) 
        ENDDO
        ! JMM : eventually add some checking to avoid locally large resto.

      CASE ( 'R','r' )
        ! horizontal extent
        zcoef2=1./(pbw +1.e-20 ) ! to avoid division by 0
        DO jj=1,npjglo
          DO ji=1,npiglo
             zv1=MAX(0., zcoef2*( glamt(ji,jj) - plon1)  )
             zv2=MAX(0., zcoef2*( plon2 - glamt(ji,jj))  )
             zv3=MAX(0., zcoef2*( gphit(ji,jj) - plat1)  )
             zv4=MAX(0., zcoef2*( plat2 - gphit(ji,jj))  )
             zpatch(ji,jj)= MIN( 1., MIN( 1., zv1,zv2,zv3,zv4 ) )
          ENDDO
        ENDDO
          ! resto all over the water column
          presto(:,:,1)= presto(:,:,1) + zpatch *zcoef
          WHERE (zpatch /= 0 ) presto(:,:,1) = MIN( presto(:,:,1), zcoef )
          DO jk=2, npk
             presto(:,:,jk)=presto(:,:,1)
          ENDDO
      END SELECT
      !
      DEALLOCATE ( zmask, zpatch )
      !
   END SUBROUTINE resto_patch
END PROGRAM cdfmkresto
