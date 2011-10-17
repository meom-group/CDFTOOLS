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
  !! Modified : 3.0 : 08/2011  : P.   Mathiot : Add zoomij, zoombat, zoomvar and time option
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

  INTEGER(KIND=4)                           :: ji, jj, jk, jt           ! dummy loop index
  INTEGER(KIND=4)                           :: ierr                     ! working integer
  INTEGER(KIND=4)                           :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                           :: npiglo, npjglo, npk, nt  ! size of the domain
  INTEGER(KIND=4)                           :: iimin, iimax             ! limit in i
  INTEGER(KIND=4)                           :: ijmin, ijmax             ! limit in j
  INTEGER(KIND=4)                           :: ncout                    ! ncid of output file
  INTEGER(KIND=4), DIMENSION(4)             :: ipk, id_varout           ! outptut variables : number of levels,

  REAL(KIND=4)                              :: rlonmin, rlonmax         ! limit in longitude
  REAL(KIND=4)                              :: rlatmin, rlatmax         ! limit in latitude
  REAL(KIND=4)                              :: rbatmin, rbatmax         ! limit in latitude
  REAL(KIND=4)                              :: rvarmin, rvarmax         ! limit in variable
  REAL(KIND=4), DIMENSION(:)  , ALLOCATABLE :: tim                      ! time counter
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask, zmask             ! 2D mask at current level
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlon, rlat               ! latitude and longitude
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rbat                     ! bathymetry 

  CHARACTER(LEN=256)                        :: cf_tfil                  ! file name
  CHARACTER(LEN=256)                        :: cf_out = 'mask_sal.nc'   ! output file
  CHARACTER(LEN=256)                        :: cv_mask                  ! variable name
  CHARACTER(LEN=256)                        :: cldum                    ! dummy string

  TYPE (variable), DIMENSION(4)             :: stypvar                  ! output attribute
 
  LOGICAL                                   :: lzoom    = .false.       ! zoom flag lat/lon
  LOGICAL                                   :: lzoomij  = .false.       ! zoom flag i/j
  LOGICAL                                   :: lzoombat = .false.       ! zoom flag bat
  LOGICAL                                   :: lzoomvar = .false.       ! zoom flag var
  LOGICAL                                   :: ltime    = .false.       ! time flag    
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmkmask T-file [-zoom lonmin lonmax latmin latmax] ...'
     PRINT *,'                   ... [-zoomij iimin iimax ijmin ijmax] ...'
     PRINT *,'                   ... [-zoombat bathymin bathymax]  ...'
     PRINT *,'                   ... [-o OUT-file ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Build a mask file from vosaline array read from the input file.' 
     PRINT *,'       It assumes that land salinity values are set to 0.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file : netcdf file with salinity.' 
     PRINT *,'                if T-file = -maskfile, we assume a reference file named ',TRIM(cn_fmsk)
     PRINT *,'                with tmask variable.' 
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
     PRINT *,'       [-zoomvar varname varmin varmax] : range of varname used to'
     PRINT *,'                        limit the area where the mask is builded. Outside'
     PRINT *,'                        this area, the mask is set to 0.'
     PRINT *,'       [-time ] : If further time step is available'
     PRINT *,'                        a mask for each time step is done'
     PRINT *,'       [-o OUT-file ] : output file name to be used in place of standard'
     PRINT *,'                        name [ ',TRIM(cf_out),' ]'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       If option -zoombat is used, file ', TRIM(cn_fzgr),' is required.'
     PRINT *,'       If option T-file is -maskfile then ', TRIM(cn_fmsk), ' is required.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out), ' or OUT-file.'
     PRINT *,'         variables : tmask, umask, vmask, fmask'
     PRINT *,'                fmask can differ from standard fmask because it does not'
     PRINT *,'                reflect the slip/noslip lateral condition.'
     STOP
  ENDIF

  ijarg = 1
  CALL getarg (ijarg, cf_tfil) ; ijarg = ijarg + 1

  DO WHILE ( ijarg <= narg ) 
    CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
    SELECT CASE ( cldum )
    !
    CASE ( '-zoom' )  ! read a zoom lat/lon area
       lzoom = .true.
       CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rlonmin
       CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rlonmax
       CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rlatmin
       CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rlatmax
    !
    CASE ( '-zoomij' )  ! read a zoom i/j area
       lzoomij = .true.
       CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) iimin
       CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) iimax
       CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmin
       CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmax
    !
    CASE ( '-zoombat' )  ! read a zoom bathy area 
       lzoombat = .true.
       CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rbatmin
       CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rbatmax
    !
    CASE ( '-zoomvar' ) ! read a zoom variable area
       lzoomvar = .true.
       CALL getarg (ijarg, cv_mask) ; ijarg = ijarg + 1 ;
       CALL getarg (ijarg, cldum)   ; ijarg = ijarg + 1 ; READ(cldum,*) rvarmin 
       CALL getarg (ijarg, cldum)   ; ijarg = ijarg + 1 ; READ(cldum,*) rvarmax 
    CASE ( '-time' )  ! create a mask for each time step of the file
       ltime=.true.
    CASE ( '-o'    )  ! change output file name
       CALL getarg (ijarg, cf_out) ; ijarg = ijarg + 1
    !
    CASE DEFAULT
        PRINT *, 'ERROR : unknown option :', TRIM(cldum)
        STOP
    END SELECT
  ENDDO

  IF ( lzoom .AND. lzoomij ) PRINT *, 'WARNING 2 spatial condition for mask'
  
  IF (.NOT. lzoomvar) cv_mask = cn_vosaline
  IF (TRIM(cf_tfil)=='-maskfile') THEN
     cv_mask = 'tmask'
     cf_tfil = cn_fmsk
     cn_z    = 'z'
  END IF    

  IF ( chkfile(cf_tfil) ) STOP ! missing file

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  nt     = getdim (cf_tfil,cn_t)

  PRINT *,' npiglo = ', npiglo
  PRINT *,' npjglo = ', npjglo
  PRINT *,' npk    = ', npk
  PRINT *,' nt     = ', nt 
  
  IF ((nt .GT. 1) .AND. (.NOT. ltime)) THEN 
     PRINT *, "WARNING nt > 1"
     PRINT *, "we used only the first time step"
     nt=1
  END IF


  ipk(1:4)                       = npk
  stypvar(1)%cname               = 'tmask'
  stypvar(2)%cname               = 'umask'
  stypvar(3)%cname               = 'vmask'
  stypvar(4)%cname               = 'fmask'

  stypvar(1:4)%cunits            = '1/0'
  stypvar(1:4)%rmissing_value    = 9999.
  stypvar(1:4)%valid_min         = 0.
  stypvar(1:4)%valid_max         = 1.

  stypvar(1)%clong_name          = 'tmask'
  stypvar(2)%clong_name          = 'umask'
  stypvar(3)%clong_name          = 'vmask'
  stypvar(4)%clong_name          = 'fmask'

  stypvar(1)%cshort_name         = 'tmask'
  stypvar(2)%cshort_name         = 'umask'
  stypvar(3)%cshort_name         = 'vmask'
  stypvar(4)%cshort_name         = 'fmask'

  stypvar(1:4)%conline_operation = 'N/A'
  stypvar(1:4)%caxis             = 'TZYX'
  stypvar(1:4)%cprecision        = 'i2'

  ncout = create      (cf_out, cf_tfil,  npiglo, npjglo, npk)
  ierr  = createvar   (ncout,    stypvar, 4,      ipk,    id_varout )
  ierr  = putheadervar(ncout,    cf_tfil,  npiglo, npjglo, npk)

  !! Allocate only usefull variable and read only usefull variable
  ALLOCATE (tmask(npiglo,npjglo), zmask(npiglo,npjglo), tim(nt))

  IF ( lzoom ) THEN
    ALLOCATE (rlon(npiglo,npjglo), rlat(npiglo,npjglo))
    rlon(:,:) = getvar(cf_tfil, cn_vlon2d, 1, npiglo, npjglo)
    rlat(:,:) = getvar(cf_tfil, cn_vlat2d, 1, npiglo, npjglo)
  ENDIF

  IF ( lzoombat ) THEN
     ALLOCATE ( rbat  (npiglo,npjglo) )
     rbat(:,:)= getvar(cn_fzgr, cn_hdepw,  1 ,npiglo, npjglo)
  END IF

  !! Now compute the mask 
  DO jt=1, nt
    IF (MOD(jt,10)==0) PRINT *,jt,'/',nt,' ...'
  DO jk=1, npk
     ! tmask
      tmask(:,:) = getvar(cf_tfil, cv_mask,  jk, npiglo, npjglo, ktime=jt)
      
      IF ( lzoomvar ) THEN
	    zmask=tmask
	    WHERE ((tmask .GE. rvarmin) .AND. (tmask .LE. rvarmax)) zmask = 1
	    WHERE ((tmask .LT. rvarmin) .OR.  (tmask .GT. rvarmax)) zmask = 0
	    tmask=zmask
      ELSE
        WHERE (tmask > 0 ) tmask = 1
        WHERE (tmask <=0 ) tmask = 0
      ENDIF

     IF ( lzoom ) THEN
        IF (rlonmax > rlonmin) THEN
           WHERE (rlon > rlonmax ) tmask = 0
           WHERE (rlon < rlonmin ) tmask = 0
        ELSE
           WHERE (rlon < rlonmin .AND. rlon > rlonmax ) tmask = 0
        END IF

        WHERE (rlat > rlatmax ) tmask = 0
        WHERE (rlat < rlatmin ) tmask = 0
     ENDIF

     IF ( lzoomij ) THEN
        tmask(1:iimin-1,:)      = 0
        tmask(iimax+1:npiglo,:) = 0
        tmask(:,ijmax+1:npjglo) = 0
        tmask(:,ijmax+1:npjglo) = 0   
     ENDIF

     IF ( lzoombat ) THEN
        WHERE (rbat < rbatmin .OR. rbat > rbatmax) tmask = 0
     ENDIF

      ierr       = putvar(ncout, id_varout(1), tmask, jk ,npiglo, npjglo, ktime=jt)
     ! umask
     zmask = 0.
     DO ji=1,npiglo-1
       DO jj=1,npjglo
        zmask(ji,jj) = tmask(ji,jj)*tmask(ji+1,jj)
       END DO
     END DO
      ierr       = putvar(ncout, id_varout(2), zmask, jk ,npiglo, npjglo, ktime=jt)
    ! vmask
     zmask=0.
     DO ji=1,npiglo
       DO jj=1,npjglo-1
        zmask(ji,jj) = tmask(ji,jj)*tmask(ji,jj+1)
       END DO
     END DO
      ierr       = putvar(ncout, id_varout(3), zmask, jk, npiglo, npjglo, ktime=jt)
     !fmask
     zmask=0.
     DO ji=1,npiglo-1
       DO jj=1,npjglo-1
        zmask(ji,jj) = tmask(ji,jj)*tmask(ji,jj+1)*tmask(ji+1,jj)*tmask(ji+1,jj+1)
       END DO
     END DO
      ierr       = putvar(ncout, id_varout(4), zmask, jk, npiglo, npjglo, ktime=jt)
  END DO  ! loop to next level
  END DO

  tim(:) = 0.
  ierr   = putvar1d(ncout, tim, nt,'T')
  ierr   = closeout(ncout              )

  PRINT *,''
  PRINT *,'Mask file ',TRIM(cf_out),' has been created' 
 
END PROGRAM cdfmkmask
