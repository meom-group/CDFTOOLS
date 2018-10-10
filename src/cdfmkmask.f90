PROGRAM cdfmkmask
  !!======================================================================
  !!                     ***  PROGRAM  cdfmkmask  ***
  !!=====================================================================
  !!  ** Purpose : Build mask file from a salinity output
  !!
  !!  ** Method  : Read vosaline and set tmask to 1 where sal is not 0
  !!               then umask, vmask and fmask are deduced from tmask
  !!               REM: the result may be locally different for fmask than
  !!                   fmask produced online as there are computed on line
  !!               merged with cdfmkmask-zone by adding a zoom option. When
  !!               used with -zoom option, the mask is 0 outside the zoom
  !!               area.
  !!
  !! History : 2.1  : 11/2005  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !! Modified  3.0  : 08/2011  : P.   Mathiot : Add zoomij, zoombat, zoomvar and time option
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!         : 4.0  : 08/2017  : P.   Mathiot : Add flood filling option
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils
  USE cdftools
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class mask
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jk, jt           ! dummy loop index
  INTEGER(KIND=4)                           :: ierr                     ! working integer
  INTEGER(KIND=4)                           :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                           :: npiglo, npjglo, npk, npt ! size of the domain
  INTEGER(KIND=4)                           :: npkk                     ! handle case without vertical dim
  INTEGER(KIND=4)                           :: iimin, iimax, iipts      ! limit in i
  INTEGER(KIND=4)                           :: ijmin, ijmax, ijpts      ! limit in j
  INTEGER(KIND=4)                           :: ncout                    ! ncid of output file
  INTEGER(KIND=4), DIMENSION(4)             :: ipk, id_varout           ! outptut variables : number of levels,
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mbathy                ! bathymetry  in levels

  REAL(KIND=4)                              :: rlonmin, rlonmax         ! limit in longitude
  REAL(KIND=4)                              :: rlatmin, rlatmax         ! limit in latitude
  REAL(KIND=4)                              :: rbatmin, rbatmax         ! limit in latitude
  REAL(KIND=4)                              :: rlonpts, rlatpts         ! seed point for lfilllonlat
  REAL(KIND=4)                              :: rvarmin, rvarmax         ! limit in variable
  REAL(KIND=4), DIMENSION(:)  , ALLOCATABLE :: rdep                     ! depth 
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask, rmask, ssmask     ! 2D masks at current level and non depth and time dependent mask
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlon, rlat               ! latitude and longitude
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rbat                     ! bathymetry 

  REAL(KIND=8), DIMENSION(:)  , ALLOCATABLE :: dtim                     ! time counter

  CHARACTER(LEN=256)                        :: cf_sfil                  ! file name
  CHARACTER(LEN=256)                        :: cf_out = 'mask_sal.nc'   ! output file
  CHARACTER(LEN=256)                        :: cf_boundary = 'boundary.txt' ! default boundary input file
  CHARACTER(LEN=256)                        :: cv_mask                  ! variable name
  CHARACTER(LEN=256)                        :: cv_dep                   ! variable name
  CHARACTER(LEN=256)                        :: cldum                    ! dummy string

  TYPE (variable), DIMENSION(4)             :: stypvar                  ! output attribute

  LOGICAL                                   :: lzoom    = .FALSE.       ! zoom flag lat/lon
  LOGICAL                                   :: lzoomij  = .FALSE.       ! zoom flag i/j
  LOGICAL                                   :: lzoombat = .FALSE.       ! zoom flag bat
  LOGICAL                                   :: lzoomvar = .FALSE.       ! zoom flag var
  LOGICAL                                   :: lfill    = .FALSE.       ! flood fill algo flag    
  LOGICAL                                   :: lfilllonlat = .FALSE.    ! flood fill algo flag    
  LOGICAL                                   :: lboundf  = .FALSE.       ! section flag var
  LOGICAL                                   :: lboundflonlat = .FALSE.  ! section flag var
  LOGICAL                                   :: ltime    = .FALSE.       ! time flag    
  LOGICAL                                   :: lmbathy  = .FALSE.       ! mbathy flag    
  LOGICAL                                   :: l2dmask  = .FALSE.       ! 2d mask flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmkmask -s S-file [-zoom lonmin lonmax latmin latmax] ...'
     PRINT *,'                   ... [-zoomij iimin iimax ijmin ijmax] ...'
     PRINT *,'                   ... [-zoombat bathymin bathymax]  ...'
     PRINT *,'                   ... [-zoomvar varname varmin varmax]  ...'
     PRINT *,'                   ... [-fill iiseed jjseed] ...'
     PRINT *,'                   ... [-bfij BOUND_IJ-file.txt] ...'
     PRINT *,'                   ... [-bflonlat BOUND_LONLAT-file.txt] ...'
     PRINT *,'                   ... [-time ] [-o OUT-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       In its simplest use, this tool builds a mask file from the salinity '
     PRINT *,'       field read from the input file, assuming that land points are set to'
     PRINT *,'       _Fill_Value.'
     PRINT *,'       '
     PRINT *,'       A growing set of options, extend this tool to a much more sophisticated'
     PRINT *,'       mask builder: use of other variables than salinity, range criteria on'
     PRINT *,'       the variables, geographical or model limits, and the possibility of'
     PRINT *,'       using any combination of the criteria.'
     PRINT *,'       '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -s S-file : netcdf file with salinity.' 
     PRINT *,'               * if S-file = -maskfile, we assume a reference file named ',TRIM(cn_fmsk)
     PRINT *,'                   with tmask variable.' 
     PRINT *,'               * if S-file = -2dmaskfile, we assume a reference file named '
     PRINT *,'                  ',TRIM(cn_fmsk),' with tmaskutil variable.' 
     PRINT *,'               * if S-file = -mbathy, we assume a reference file named '
     PRINT *,'                   bathylevel.nc with mbathy variable, giving the number of '
     PRINT *,'                   levels in the ocean.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-zoom lonmin lonmax latmin latmax] : geographical windows used to'
     PRINT *,'                        limit the area where the mask is builded. Outside'
     PRINT *,'                        this area, the mask is set to 0.'
     PRINT *,'       [-zoomij iimin iimax ijmin ijmax] : model grid windows used to'
     PRINT *,'                        limit the area where the mask is builded. Outside'
     PRINT *,'                        this area, the mask is set to 0.'
     PRINT *,'       [-zoombat bathymin bathymax] : depth windows used to'
     PRINT *,'                        limit the area where the mask is builded. Outside'
     PRINT *,'                        this area, the mask is set to 0.' 
     PRINT *,'                        Need mesh_zgr.nc'
     PRINT *,'       [-zoomvar varname varmin varmax] : range of varname variable used to'
     PRINT *,'                        limit the area where the mask is builded. Outside'
     PRINT *,'                        this area, the mask is set to 0.'
     PRINT *,'       [-fill iiseed jjseed] : mask everything except the cells into the'
     PRINT *,'                        non mask area where the point (iiseed,jjseed) is.'
     PRINT *,'       [-filllonlat lon lat] : mask everything except the cells into the'
     PRINT *,'                        non mask area where the point (lon,lat) is.'
     PRINT *,'       [-bfij BOUND_IJ-file.txt] : Specify a boundary file expressed in'
     PRINT *,'                        model grid point. It is used in conjunction with'
     PRINT *,'                        -fill series of options, to define an arbitrary'
     PRINT *,'                        closed area where mask values will be 1.'
     PRINT *,'                        See format of boundary files in dedicated paragraph.'
     PRINT *,'       [-bflonlat BOUND_LONLAT-file.txt] :  Specify a boundary file expressed '
     PRINT *,'                        in geographical coordinates. It is used in conjunction'
     PRINT *,'                        with the -fill series of options, to define an '
     PRINT *,'                        arbitrary closed area, where mask values will be 1.'
     PRINT *,'                        See format of boundary files in dedicated paragraph.'
     PRINT *,'       [-time ] : Build a mask corresponding to each time step of the '
     PRINT *,'                       input file'
     PRINT *,'       [-o OUT-file ] : output file name to be used in place of standard'
     PRINT *,'                        name [ ',TRIM(cf_out),' ]'
     PRINT *,'      '
     PRINT *,'     FORMAT OF BOUNDARY FILES :'
     PRINT *,'       Boundary files describes a set of model sections forming a closed line.'
     PRINT *,'       Each section (segment of the boundary) is defined by the position '
     PRINT *,'       (either in model I-J, or in geographical LON-LAT) of the ending points'
     PRINT *,'       of the respective section. Therefore, each section in the boundary file'
     PRINT *,'       corresponds to 2 lines :'
     PRINT *,'       NAME'
     PRINT *,'       xmin xmax ymin ymax logical_flag'
     PRINT *,'       NAME2 '
     PRINT *,'       .......'
     PRINT *,'       EOF'
     PRINT *,'         where xmin xmax ymin ymax specify the ending points of the section,'
     PRINT *,'         and   logical_flag (either T ior F), indicates if the corresponding'
     PRINT *,'         section is to be taken into account (T). '
     PRINT *,'       Note that last line of the file should be just ''EOF'' '
     PRINT *,'      '
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       If option -zoombat is used, file ', TRIM(cn_fzgr),' is required.'
     PRINT *,'       If option T-file is -maskfile then ', TRIM(cn_fmsk), ' is required.'
     PRINT *,'       If option T-file is -mbathy then bathylevel.nc and ', TRIM(cn_fzgr) 
     PRINT *,'        are required.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out), ' or OUT-file.'
     PRINT *,'         variables : tmask, umask, vmask, fmask'
     PRINT *,'                fmask can differ from standard fmask because it does not'
     PRINT *,'                reflect the slip/noslip lateral condition.'
     STOP 
  ENDIF

  ijarg = 1
  DO WHILE ( ijarg <= narg ) 
     CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ('-s','-f') ; CALL getarg (ijarg, cf_sfil) ; ijarg = ijarg + 1
     CASE ( '-zoom' )  ! read a zoom lat/lon area
        lzoom = .TRUE.
        CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rlonmin
        CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rlonmax
        CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rlatmin
        CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rlatmax
        !
     CASE ( '-zoomij' )  ! read a zoom i/j area
        lzoomij = .TRUE.
        CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) iimin
        CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) iimax
        CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmin
        CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmax
        !
     CASE ( '-zoombat' )  ! read a zoom bathy area 
        lzoombat = .TRUE.
        CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rbatmin
        CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rbatmax
        !
     CASE ( '-zoomvar' ) ! read a zoom variable area
        lzoomvar = .TRUE.
        CALL getarg (ijarg, cv_mask) ; ijarg = ijarg + 1 ;
        CALL getarg (ijarg, cldum)   ; ijarg = ijarg + 1 ; READ(cldum,*) rvarmin 
        CALL getarg (ijarg, cldum)   ; ijarg = ijarg + 1 ; READ(cldum,*) rvarmax 
     CASE ( '-fill' )  ! read a seed point and a boundary file
        lfill = .TRUE.
        CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) iipts
        CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ijpts
     CASE ( '-filllonlat' )  ! read a seed point and a boundary file
        lfilllonlat = .TRUE.
        CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rlonpts
        CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rlatpts
     CASE ( '-bf' )    ! read boundary file name
        lboundf=.TRUE.
        CALL getarg (ijarg, cf_boundary) ; ijarg = ijarg + 1
     CASE ( '-bflonlat' )    ! read boundary file name
        lboundflonlat=.TRUE.
        CALL getarg (ijarg, cf_boundary) ; ijarg = ijarg + 1
     CASE ( '-time' )  ! create a mask for each time step of the file
        ltime=.TRUE.
     CASE ( '-o'    )  ! change output file name
        CALL getarg (ijarg, cf_out) ; ijarg = ijarg + 1
        !
     CASE DEFAULT  ; PRINT *, ' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( lfill .AND. lfilllonlat ) THEN
     PRINT *, 'E R R O R: 2 seeds for the filling specified (-fill and -filllonlat), STOP'; STOP 99
  END IF
  IF ( lboundf .AND. lboundflonlat ) THEN
     PRINT *, 'E R R O R: 2 boundary files for the filling specified (-bf and -bflonlat), STOP'; STOP 99
  END IF

  IF ( lzoom .AND. lzoomij ) PRINT *, 'WARNING 2 spatial condition for mask'

  IF (.NOT. lzoomvar) cv_mask = cn_vosaline
  IF (TRIM(cf_sfil)=='-maskfile') THEN
     cv_mask = cn_tmask
     cf_sfil = cn_fmsk
     cn_z    = 'z'
  END IF

  IF (TRIM(cf_sfil)=='-2dmaskfile') THEN
     cv_mask = 'tmaskutil'
     cf_sfil = cn_fmsk
     cn_z    = 'z'
     l2dmask=.TRUE.
  END IF

  IF (TRIM(cf_sfil)=='-mbathy') THEN
     cv_mask = cn_mbathy
     cv_dep  = 'nav_lev'
     cf_sfil = 'bathylevel.nc'
     cn_z    = 'z'
     lmbathy = .TRUE.
     IF ( chkfile(cn_fzgr) ) STOP 99 ! missing file
  END IF

  IF ( chkfile(cf_sfil) ) STOP 99 ! missing file

  npiglo = getdim (cf_sfil,cn_x)
  npjglo = getdim (cf_sfil,cn_y)
  IF ( lmbathy ) THEN
     npk  = getdim (cn_fzgr,cn_z)
     ALLOCATE ( rdep(npk) )
  ELSE IF ( l2dmask ) THEN
     PRINT *,' npk is forced to 1'
     npk  = 0
  ELSE
     npk  = getdim (cf_sfil,cn_z)
  ENDIF
  npt    = getdim (cf_sfil,cn_t)

  PRINT *,' npiglo = ', npiglo
  PRINT *,' npjglo = ', npjglo
  PRINT *,' npk    = ', npk
  PRINT *,' npt    = ', npt

  IF ( npt == 0 ) THEN
     PRINT *,' npt is forced to 1'
     npt = 1
  ENDIF

  IF ((npt > 1) .AND. (.NOT. ltime)) THEN 
     PRINT *, "WARNING npt > 1"
     PRINT *, "we used only the first time step"
     npt=1
  END IF

  IF ( npk == 0 ) THEN ; npkk = 1
  ELSE                 ; npkk = npk
  ENDIF

  ALLOCATE (dtim(npt))
  CALL CreateOutput

  !! Allocate only usefull variable and read only usefull variable
  ALLOCATE (tmask(npiglo,npjglo), rmask(npiglo,npjglo), ssmask(npiglo,npjglo))
  ssmask(:,:) = 1.

  !! apply constraint constant over time and depth
  !! mbathy constrain
  IF ( lmbathy ) THEN
     ALLOCATE (mbathy(npiglo,npjglo))
     mbathy(:,:) = getvar(cf_sfil, cv_mask, 1, npiglo, npjglo)
     WHERE (mbathy < jk ) ssmask = 0.
  ENDIF

  !! bathy constraint
  IF ( lzoombat ) THEN
     IF ( chkfile(cn_fzgr) ) STOP 99 ! missing file
     ALLOCATE ( rbat  (npiglo,npjglo) )
     rbat(:,:)= getvar(cn_fbathymet, cn_bathymet,  1 ,npiglo, npjglo)
     WHERE (rbat < rbatmin .OR. rbat > rbatmax) ssmask = 0
  ENDIF

  !! lat/lon constrain
  IF ( lzoom ) THEN
     ALLOCATE (rlon(npiglo,npjglo), rlat(npiglo,npjglo))
     rlon(:,:) = getvar(cf_sfil, cn_vlon2d, 1, npiglo, npjglo)
     rlat(:,:) = getvar(cf_sfil, cn_vlat2d, 1, npiglo, npjglo)
     IF (rlonmax > rlonmin) THEN
        WHERE (rlon > rlonmax ) ssmask = 0
        WHERE (rlon < rlonmin ) ssmask = 0
     ELSE
        WHERE (rlon < rlonmin .AND. rlon > rlonmax ) ssmask = 0
     END IF

     WHERE (rlat > rlatmax ) ssmask = 0
     WHERE (rlat < rlatmin ) ssmask = 0
  ENDIF

  !! i/j constrain
  IF ( lzoomij ) THEN
     ssmask(1:iimin-1,:     ) = 0   ! West
     ssmask(iimax+1:npiglo,:) = 0   ! East
     ssmask(:,ijmax+1:npjglo) = 0   ! North
     ssmask(:,1:ijmin-1     ) = 0   ! South
  ENDIF

  !! Now compute the mask 
  DO jt=1, npt
     IF (MOD(jt,10)==0) PRINT *,jt,'/',npt,' ...'
     DO jk=1, npkk
        PRINT *, jk,'/',npkk

        tmask(:,:) = getvar(cf_sfil, cv_mask,  jk, npiglo, npjglo, ktime=jt)
        tmask(:,:) = tmask(:,:) * ssmask(:,:)

        !! variable constrain
        IF ( lzoomvar ) THEN
           rmask=tmask
           WHERE ((tmask >= rvarmin) .AND. (tmask <= rvarmax)) rmask = 1
           WHERE ((tmask <  rvarmin) .OR.  (tmask >  rvarmax)) rmask = 0
           tmask=rmask
        ELSE
           WHERE (tmask > 0 ) tmask = 1
           WHERE (tmask <=0 ) tmask = 0
        ENDIF

        !! fill constrain
        IF ( lfill .OR. lfilllonlat ) THEN
           CALL FillMask(tmask,iipts,ijpts,rlonpts,rlatpts)
        ENDIF

        !! write t- u- v- mask
        ierr       = putvar(ncout, id_varout(1), tmask, jk ,npiglo, npjglo, ktime=jt)
        ! umask
        rmask = 0.
        DO ji=1,npiglo-1
           DO jj=1,npjglo
              rmask(ji,jj) = tmask(ji,jj)*tmask(ji+1,jj)
           END DO
        END DO
        ierr       = putvar(ncout, id_varout(2), rmask, jk ,npiglo, npjglo, ktime=jt)
        ! vmask
        rmask=0.
        DO ji=1,npiglo
           DO jj=1,npjglo-1
              rmask(ji,jj) = tmask(ji,jj)*tmask(ji,jj+1)
           END DO
        END DO
        ierr       = putvar(ncout, id_varout(3), rmask, jk, npiglo, npjglo, ktime=jt)
        !fmask
        rmask=0.
        DO ji=1,npiglo-1
           DO jj=1,npjglo-1
              rmask(ji,jj) = tmask(ji,jj)*tmask(ji,jj+1)*tmask(ji+1,jj)*tmask(ji+1,jj+1)
           END DO
        END DO
        ierr       = putvar(ncout, id_varout(4), rmask, jk, npiglo, npjglo, ktime=jt)
     END DO  ! loop to next level
  END DO ! loop to next time

  ierr   = closeout(ncout)

  PRINT *,''
  PRINT *,'Mask file ',TRIM(cf_out),' has been created' 

CONTAINS

  SUBROUTINE FillMask(pmask, kipts, kjpts, plonpts, platpts)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE FillMask  ***
    !!
    !! ** Purpose :  build mask when using flood filling option
    !!
    !! ** Method  : Use FillPool2D routine after eventually setting up the
    !!              boundary of the domain to fill.
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(:,:),  INTENT(inout) :: pmask
    REAL(KIND=4),                  INTENT(in   ) :: plonpts, platpts
    INTEGER(KIND=4),               INTENT(inout) :: kipts, kjpts ! seeding point coordinates

    INTEGER(KIND=4)                              :: js           ! dummy loop index
    INTEGER(KIND=4)                              :: ii, ij       ! working integer
    INTEGER(KIND=4)                              :: iunit=10
    INTEGER(KIND=4)                              :: ipts         ! local integer with number of points per section
    INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: imask

    REAL(KIND=4)                                 :: zlonmin, zlonmax, zlatmin, zlatmax
    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE    :: zmskline
    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE    :: zmsk_bck
    REAL(KIND=4), DIMENSION(:)  , ALLOCATABLE    :: zxx, zyy       ! working variables

    CHARACTER(LEN=256)                           :: cldum          ! dummy char variable
    CHARACTER(LEN=256)                           :: cline          ! dummy char variable
    CHARACTER(LEN=256)                           :: cl_section     ! section names

    LOGICAL :: ll_inc, ll_section
    !!----------------------------------------------------------------------
    ALLOCATE(zmskline(npiglo, npjglo), zmsk_bck(npiglo, npjglo))
    ALLOCATE(imask(npiglo, npjglo))

    zmskline(:,:) = 0.0
    zmsk_bck(:,:) = pmask(:,:)

    IF (lboundf .OR. lboundflonlat) THEN
       PRINT *,''
       PRINT *,'Boundary file: ',TRIM(cf_boundary),' is used to close the basin'
       PRINT *,''

       IF ( chkfile(cf_boundary) ) STOP 99 ! missing file
       ! allocate arrays which are only used in this case ( lboundf OR lboundflonlat )
       ALLOCATE(zxx(npiglo+npjglo), zyy(npiglo+npjglo)) ! dimension specified in broken_line
       ! optimal dimension could be ABS(imax-imin +1) + ABS(jmax-jmin+1) - 1

       OPEN(iunit, FILE=cf_boundary)
       ll_section = .TRUE.
       DO WHILE (ll_section)
          zxx(:)=1; zyy(:)=1

          ! read section name
          READ(iunit,'(a)') cl_section
          IF (TRIM(cl_section) == 'EOF' ) THEN
             ll_section = .FALSE.
          ELSE
             ! read section coordinates
             IF (lboundflonlat) THEN
                READ(iunit,*) zlonmin, zlonmax, zlatmin, zlatmax, ll_inc
                CALL cdf_findij ( zlonmin, zlonmax, zlatmin, zlatmax, iimin, iimax, ijmin, ijmax, &
                     &            cd_coord=cn_fhgr, cd_point='T', cd_verbose='T')
             ELSE
                READ(iunit,*) iimin, iimax, ijmin, ijmax, ll_inc
             ENDIF

             ! get index of cell included into the section
             CALL broken_line(iimin, iimax, ijmin, ijmax, zxx, zyy, ipts, npiglo, npjglo)

             ! mask boundary and keep location in zmskline
             DO js=1,ipts
                ii=zxx(js)
                ij=zyy(js)
                IF (ll_inc)  zmskline(ii,ij)=1.0 * zmsk_bck(ii,ij)
                pmask(ii,ij)=0.0
             END DO
          ENDIF
       END DO
       CLOSE(iunit)
    ELSE
       PRINT *,''
       PRINT *, 'NO BOUNDARIES ARE ADDED TO THE INPUT FILE TO CLOSE THE BASIN'
       PRINT *,''
    END IF

    ! fill area
    ! find ij if lon/lat given
    IF (lfilllonlat) THEN
       CALL cdf_findij ( plonpts, plonpts, platpts, platpts, kipts, kipts, kjpts, kjpts, &
            &            cd_coord=cn_fhgr, cd_point='T', cd_verbose='F')   
    ENDIF
    imask = NINT(pmask,2)
    CALL FillPool2D(kipts, kjpts, imask, -1) ! fill pool (use -1 to flag the
    ! area and avoid infinit loop in the algo

    ! keep only the point selected by the flood filling algo
    WHERE ( imask == -1 )
       pmask = 1.0
    ELSE WHERE
       pmask = 0.0
    END WHERE

    ! apply mskline condition (ll_inc)
    WHERE (zmskline == 1.0)
       pmask = zmsk_bck
    END WHERE

  END SUBROUTINE FillMask

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ipk(1:4)                       = npkk

    stypvar(1)%cname               = cn_tmask
    stypvar(2)%cname               = cn_umask
    stypvar(3)%cname               = cn_vmask
    stypvar(4)%cname               = cn_fmask

    stypvar(1:4)%cunits            = '1/0'
    stypvar(1:4)%rmissing_value    = 9999.
    stypvar(1:4)%valid_min         = 0.
    stypvar(1:4)%valid_max         = 1.

    stypvar(1)%clong_name          = cn_tmask
    stypvar(2)%clong_name          = cn_umask
    stypvar(3)%clong_name          = cn_vmask
    stypvar(4)%clong_name          = cn_fmask

    stypvar(1)%cshort_name         = cn_tmask
    stypvar(2)%cshort_name         = cn_umask
    stypvar(3)%cshort_name         = cn_vmask
    stypvar(4)%cshort_name         = cn_fmask

    stypvar(1:4)%conline_operation = 'N/A'
    stypvar(1:4)%caxis             = 'TZYX'
    stypvar(1:4)%cprecision        = 'i2'

    ncout = create      (cf_out, cf_sfil,  npiglo, npjglo, npk)
    ierr  = createvar   (ncout,    stypvar, 4,      ipk,    id_varout )

    IF ( lmbathy ) THEN ; rdep(:) = getvare3(cn_fzgr, cv_dep ,npk)
       ; ierr  = putheadervar(ncout,    cf_sfil,  npiglo, npjglo, npk, pdep=rdep, cdep='nav_lev')
    ELSE                ; ierr  = putheadervar(ncout,    cf_sfil,  npiglo, npjglo, npk)
    ENDIF
    dtim = 0.d0
    ierr = putvar1d(ncout, dtim, npt,'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfmkmask
