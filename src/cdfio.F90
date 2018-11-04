 MODULE cdfio
  !!======================================================================
  !!                     ***  MODULE  cdfio  ***
  !! Implement all I/O related to netcdf in CDFTOOLS
  !!=====================================================================
  !! History : 2.1 : 2005  : J.M. Molines   : Original code
  !!               : 2009  : R. Dussin      : add putvar_0d function
  !!           3.0 : 12/2010 : J.M. Molines : Doctor + Licence     
  !! Modified: 3.0 : 08/2011 : P.   Mathiot : Add chkvar function           
  !!         : 4.0 : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------

  !!----------------------------------------------------------------------
  !!   routines      : description
  !! .............................
  !!   ERR_HDL       : Error Handler routine to catch netcdf errors
  !!   gettimeseries : print a 2 column array (time, variable) for a given
  !!                   file, variable and depth
  !!
  !!   functions     : description
  !! .............................
  !!   chkfile       : check the existence of a file
  !!   chkvar        : check the existence of a variable in a file
  !!   closeout      : close output file
  !!   copyatt       : copy attributes from a file taken as model
  !!   create        : create a netcdf data set
  !!   createvar     : create netcdf variables in a new data set
  !!   cvaratt       : change some var attributes
  !!   edatt_char    : edit attributes of char type
  !!   edatt_r4      : edit attributes of float type
  !!   getatt        : get attributes of a variable
  !!   getdim        : return the value of the dimension passed as argument
  !!   getipk        : get the vertical dimension of the variable
  !!   getnvar       : get the number of variable in a file
  !!   getspval      : get spval of a given variable
  !!   getvar1d      : read 1D variable (eg depth, time_counter) from a file
  !!   getvar3d      : read 3D variable  at once
  !!   getvaratt     : read variable attributes
  !!   gettimeatt    : get time attributes
  !!   getvar        : read the variable
  !!   getvare3      : read e3 type variable
  !!   getvarid      : get the varid of a variable in a file
  !!   getvarname    : get the name of a variable, according to its varid
  !!   getvarxz      : get a x-z slice of 3D data
  !!   getvaryz      : get a y-z slice of 3D data
  !!   getvdim       : get the number of dim of a variable
  !!   getvardim     : get the values of the vertical coordinates variable
  !!   ncopen        : open a netcdf file and return its ncid
  !!   putatt        : write variable attribute
  !!   puttimeatt    : write time variable attribute
  !!   putheadervar  : write header variables such as nav_lon, nav_lat etc ... from a file taken
  !!                 : as template
  !!   putvar0d      : write a 0d variable (constant)
  !!   putvar1d4     : write a 1d variable
  !!   putvari2      : write a 2d Integer*2 variable
  !!   putvarr4      : write a 2d Real*4 variable
  !!   putvarr8      : write a 2d Real*8 variable
  !!   putvarzo      : write a zonally integrated/mean field
  !!   reputvarr4    : re-write a real*4 variable
  !!   reputvar1d4   : re-write a real*4 1d variable 
  !!------------------------------------------------------------------------------------------------------
  USE netcdf        
  USE modcdfnames

  IMPLICIT NONE

  ! Global var defining the mesh_mask_version (not in namelist)
  CHARACTER(LEN=256), PUBLIC  :: cg_zgr_ver='v3.6'
  LOGICAL, PUBLIC             :: lg_vvl=.FALSE.  ! flag (global and public) for vvl ( set to T when vvl )

  PRIVATE 
  INTEGER(KIND=4) :: nid_x, nid_y, nid_z, nid_t, nid_lat, nid_lon, nid_dep, nid_tim
  INTEGER(KIND=4) :: nid_lon1d, nid_lat1d
  LOGICAL         :: l_mbathy = .false.
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mbathy         !: for reading e3._ps in nemo3.x
  REAL(KIND=4),    DIMENSION(:,:), ALLOCATABLE :: e3t_ps, e3w_ps !: for reading e3._ps in nemo3.x
  REAL(KIND=4),    DIMENSION(:),   ALLOCATABLE :: e3t_0, e3w_0   !: for readinf e3._ps in nemo3.x
  
  INTEGER(KIND=4)    :: nstart_date  = -1                      !# from global file attribute
  CHARACTER(LEN=256) :: ctime_units  = 'seconds since 0001-01-01 00:00:00'
  CHARACTER(LEN=256) :: ctime_origin = 'N/A'                   !# 
  CHARACTER(LEN=256) :: calendar     = 'N/A'                   !# gregorian noleap xxxd
  CHARACTER(LEN=256) :: config                                 !# config name associated with file
  CHARACTER(LEN=256) :: ccase                                  !# case name associated with  file
  CHARACTER(LEN=10 ) :: cfreq                                  !# raw model output frequency (5d, 30d 1h ..)

  TYPE, PUBLIC ::   variable 
     CHARACTER(LEN=256) :: cname             !# variable name
     CHARACTER(LEN=256) :: cunits            !# variable unit
     REAL(KIND=4)       :: rmissing_value    !# variable missing value or spval
     REAL(KIND=4)       :: valid_min         !# valid minimum
     REAL(KIND=4)       :: valid_max         !# valid maximum
     REAL(KIND=4)       :: scale_factor=1.   !# scale factor
     REAL(KIND=4)       :: add_offset=0.     !# add offset
     REAL(KIND=4)       :: savelog10=0.      !# flag for log10 transform
     INTEGER(KIND=4)    :: iwght=1           !# weight of the variable for cdfmoy_weighted
     INTEGER(KIND=4), DIMENSION(4) :: ichunk !# chunksize
     CHARACTER(LEN=256) :: clong_name        !# Long Name of the variable
     CHARACTER(LEN=256) :: cshort_name       !# short name of the variable
     CHARACTER(LEN=256) :: conline_operation !# ???
     CHARACTER(LEN=256) :: caxis             !# string defining the dim of the variable
     CHARACTER(LEN=256) :: cprecision='r4'   !# possible values are i2, r4, r8
  END TYPE variable

  TYPE, PUBLIC  :: ncfile                         ! logical structure reflecting file structure
     INTEGER(KIND=4)                              :: ncid    ! file ncid
     INTEGER(KIND=4)                              :: ndims   ! number of dims
     INTEGER(KIND=4)                              :: nvars   ! number of vars
     INTEGER(KIND=4)                              :: natts   ! number of global attributes
     INTEGER(KIND=4)                              :: iunlim  ! ID of unlimited dimension
     INTEGER(KIND=4)                              :: npi     ! i-size of file
     INTEGER(KIND=4)                              :: npj     ! j-size of file
     INTEGER(KIND=4)                              :: npk     ! k-size of file
     INTEGER(KIND=4)                              :: npt     ! t-size of file
     INTEGER(KIND=4)                              :: npb     ! time_bound size
     INTEGER(KIND=4)                              :: idx     ! x dimid
     INTEGER(KIND=4)                              :: idy     ! y dimid
     INTEGER(KIND=4)                              :: idz     ! z dimid
     INTEGER(KIND=4)                              :: idt     ! t dimid
     INTEGER(KIND=4)                              :: idb     ! time-bound  dimid
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: ideflat ! deflate level (nvar)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nvatt   ! number of att of each variable (var)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nvid    ! varid of each variable (var)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nvdim   ! dimension of each variable (var)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: itype   ! type of each variable (var)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nlen    ! len of each dimension ( ndims)
     INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ichunk  ! size of chunk (nvar, ndims)
     INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: idimids ! dimids of each variable (nvar, ndims) 
     CHARACTER(LEN=255)                            :: c_fnam ! name of working file
     CHARACTER(LEN=80 ), DIMENSION(:), ALLOCATABLE :: c_vnam ! name of each variable (var)
     CHARACTER(LEN=80 ), DIMENSION(:), ALLOCATABLE :: c_dnam ! name of each dimension (ndims)
     LOGICAL,           DIMENSION(:), ALLOCATABLE :: lconti  ! contiguous flag (nvar)
     !   extra information for global attribute 
     INTEGER(KIND=4)                :: number_total          ! DOMAIN_number_total
     INTEGER(KIND=4)                :: number                ! DOMAIN_number
     INTEGER(KIND=4), DIMENSION(2)  :: idimensions_ids       ! DOMAIN_dimensions_ids
     INTEGER(KIND=4), DIMENSION(2)  :: isize_global          ! DOMAIN_size_global
     INTEGER(KIND=4), DIMENSION(2)  :: isize_local           ! DOMAIN_size_local
     INTEGER(KIND=4), DIMENSION(2)  :: iposition_first       ! DOMAIN_position_first
     INTEGER(KIND=4), DIMENSION(2)  :: iposition_last        ! DOMAIN_position_last 
     INTEGER(KIND=4), DIMENSION(2)  :: ihalo_size_start      ! DOMAIN_halo_size_start
     INTEGER(KIND=4), DIMENSION(2)  :: ihalo_size_end        ! DOMAIN_halo_size_end
     CHARACTER(LEN=80)              :: c_type                ! DOMAIN_type
  END TYPE ncfile


  INTEGER(KIND=4), PARAMETER :: jp_missing_nm = 3
  
  CHARACTER(LEN=256), DIMENSION(jp_missing_nm) :: & ! take care of same length for each element
        & cl_missing_nm = (/'missing_value','Fillvalue    ','_Fillvalue   '/)
  CHARACTER(LEN=256 ) :: cl_dum              !# dummy char argument

  INTERFACE putvar
     MODULE PROCEDURE putvarr8, putvarr4, putvari2, putvarzo, reputvarr4, putvare3
  END INTERFACE

  INTERFACE putvar1d   
     MODULE PROCEDURE putvar1d4, putvar1d8, reputvar1d4, reputvar1d8
  END INTERFACE

  INTERFACE putvar0d   
     MODULE PROCEDURE putvar0dt, putvar0ds
  END INTERFACE

  INTERFACE atted
     MODULE PROCEDURE atted_char, atted_r4
  END INTERFACE

  INTERFACE getatt
     MODULE PROCEDURE  getattr, getatti, getattc
  END INTERFACE

  PUBLIC :: chkfile, chkvar
  PUBLIC :: copyatt, create, createvar, getvaratt, cvaratt, gettimeatt
  PUBLIC :: putatt, putheadervar, putvar, putvar1d, putvar0d, atted, puttimeatt
  PUBLIC :: getatt, getdim, getvdim, getdimvar,getipk, getnvar, getvarname, getvarid
  PUBLIC :: getvar, getvarxz, getvaryz, getvar1d, getvare3, getvar3d, getvar3dt, getvar4d, getspval
  PUBLIC :: gettimeseries
  PUBLIC :: closeout, ncopen
  PUBLIC :: ERR_HDL


  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class system
  !!----------------------------------------------------------------------

CONTAINS

  INTEGER(KIND=4) FUNCTION copyatt (cdvar, kidvar, kcin, kcout)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION copyatt  ***
    !!
    !! ** Purpose :   Copy attributes for variable cdvar, which have id 
    !!                kidvar in kcout, from file id kcin
    !!
    !! ** Method  :   Use NF90_COPY_ATT
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdvar
    INTEGER(KIND=4),  INTENT(in) :: kidvar, kcin, kcout

    INTEGER(KIND=4)    :: ja
    INTEGER(KIND=4)    :: istatus, idvar, iatt
    CHARACTER(LEN=256) :: clatt
    !!----------------------------------------------------------------------
    IF ( kcin /= -9999) THEN    ! there is a reference file open
       istatus = NF90_INQ_VARID(kcin, cdvar, idvar)
       istatus = NF90_INQUIRE_VARIABLE(kcin, idvar, natts=iatt)
       DO ja = 1, iatt
          istatus = NF90_INQ_ATTNAME(kcin,idvar,ja,clatt)
          istatus = NF90_COPY_ATT(kcin,idvar,clatt,kcout,kidvar)
       END DO
    ELSE                        ! no reference file
       SELECT CASE (TRIM(cdvar) )
       CASE ('nav_lon', 'lon', 'x', 'longitude' )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'units',     'degrees_east')
          istatus=NF90_PUT_ATT(kcout, kidvar, 'valid_min', -180.         )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'valid_max',  180.         )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'long_name', 'Longitude'   )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'nav_model', 'Default grid')
       CASE ('nav_lat' ,'lat', 'y', 'latitude' )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'units',    'degrees_north')
          istatus=NF90_PUT_ATT(kcout, kidvar, 'valid_min', -90.          )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'valid_max',  90.          )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'long_name','Latitude'     )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'nav_model', 'Default grid')
       CASE ('time_counter', 'time', 't' )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'calendar',   calendar     )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'units',      ctime_units  )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'time_origin',ctime_origin )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'title',      'Time'       )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'long_name',  'Time axis'  )
       CASE ('deptht', 'depthu' ,'depthv' , 'depthw', 'dep', 'gdept'     )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'units',      'm'          )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'positive',  'unknown'     )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'valid_min',  0.           )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'valid_max',  5875.        )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'title',     TRIM(cdvar)   )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'long_name', 'Vertical Levels')
       CASE ('sigma', 'levels')
          istatus=NF90_PUT_ATT(kcout, kidvar, 'units', 'kg/m3'           )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'positive', 'unknown'      )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'valid_min', 0.            )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'valid_max', 40.           )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'title', TRIM(cdvar)       )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'long_name', 'Sigma bin limits')
       END SELECT
    ENDIF

    copyatt = istatus
  END FUNCTION copyatt


  INTEGER(KIND=4) FUNCTION create( cdfile, cdfilref ,kx,ky,kz ,cdep, cdepvar, &
       &                           cdlonvar, cdlatvar,  ld_xycoo, ld_nc4 )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION create  ***
    !!
    !! ** Purpose : Create the file, and creates dimensions, and copy attributes 
    !!              from a cdilref reference file (for the nav_lon, nav_lat etc ...)
    !!              If optional cdep given : take as depth variable name instead of 
    !!              cdfilref. Return the ncid of the created file, and leave it open
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),           INTENT(in) :: cdfile, cdfilref ! input file and reference file
    INTEGER(KIND=4),            INTENT(in) :: kx, ky, kz       ! dimension of the variable
    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: cdep     ! name of vertical dim name if not standard
    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: cdepvar  ! name of vertical var name if it differs
                                                       ! from vertical dimension name
    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: cdlonvar ! name of 1D longitude
    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: cdlatvar ! name of 1D latitude
    LOGICAL,          OPTIONAL, INTENT(in) :: ld_xycoo ! if false then DO NOT read nav_lat nav_lat from input file
    LOGICAL,          OPTIONAL, INTENT(in) :: ld_nc4   ! create NETCDF4 file with chunking and deflation

    INTEGER(KIND=4)               :: istatus, icout, incid, idum
    INTEGER(KIND=4) ,DIMENSION(4) :: invdim
    CHARACTER(LEN=256)            :: cldep, cldepref, cldepvar, clonvar, clatvar
    LOGICAL                       :: ll_xycoo, ll_nc4
    !!----------------------------------------------------------------------
    IF ( PRESENT (ld_nc4 ) ) THEN 
       ll_nc4 = ld_nc4
    ELSE
       ll_nc4 = .false. 
    ENDIF
#if defined key_netcdf4
    IF ( ll_nc4 ) THEN
      istatus = NF90_CREATE(cdfile,cmode=or(NF90_CLOBBER,NF90_NETCDF4     ), ncid=icout)
    ELSE
      istatus = NF90_CREATE(cdfile,cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=icout)
    ENDIF
#else
    istatus = NF90_CREATE(cdfile,cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=icout)
#endif
    istatus = NF90_DEF_DIM(icout, cn_x, kx, nid_x)
    istatus = NF90_DEF_DIM(icout, cn_y, ky, nid_y)

    IF ( kz /= 0 ) THEN
       ! try to find out the name I will use for depth dimension in the new file ...
       IF (PRESENT (cdep) ) THEN
          cldep = cdep
          idum=getdim(cdfilref,cldep,cldepref)   ! look for depth dimension name in ref file
         IF (cldepref =='unknown' ) cldepref=cdep
       ELSE 
          idum=getdim(cdfilref,cn_z,cldep   )   ! look for depth dimension name in ref file
          cldepref=cldep
       ENDIF
       cldepvar=cldep
       istatus = NF90_DEF_DIM(icout,TRIM(cldep),kz, nid_z)
       IF (PRESENT (cdepvar) ) THEN
         cldepvar=cdepvar
       ENDIF
    ENDIF


    istatus = NF90_DEF_DIM(icout,cn_t,NF90_UNLIMITED, nid_t)

    invdim(1) = nid_x ; invdim(2) = nid_y ; invdim(3) = nid_z ; invdim(4) = nid_t

    ! Open reference file if any,  otherwise set ncid to flag value (for copy att)
    IF ( TRIM(cdfilref) /= 'none' ) THEN
       istatus = NF90_OPEN(cdfilref,NF90_NOWRITE,incid)
    ELSE
       incid = -9999
    ENDIF

    IF (PRESENT (ld_xycoo) ) THEN
      ll_xycoo = ld_xycoo
    ELSE
      ll_xycoo = .true.
    ENDIF 
    ! define variables and copy attributes
    IF ( ll_xycoo ) THEN
       istatus = NF90_DEF_VAR(icout,cn_vlon2d,NF90_FLOAT,(/nid_x, nid_y/), nid_lon)
       istatus = copyatt(cn_vlon2d, nid_lon,incid,icout)
       istatus = NF90_DEF_VAR(icout,cn_vlat2d,NF90_FLOAT,(/nid_x, nid_y/), nid_lat)
       istatus = copyatt(cn_vlat2d, nid_lat,incid,icout)
    ENDIF
    IF ( PRESENT(cdlonvar) ) THEN
       istatus = NF90_DEF_VAR(icout,cdlonvar,NF90_FLOAT,(/nid_x/), nid_lon1d)
    ENDIF
    IF ( PRESENT(cdlatvar) ) THEN
       istatus = NF90_DEF_VAR(icout,cdlatvar,NF90_FLOAT,(/nid_y/), nid_lat1d)
    ENDIF
    IF ( kz /= 0 ) THEN
       istatus = NF90_DEF_VAR(icout,TRIM(cldepvar),NF90_FLOAT,(/nid_z/), nid_dep)
       ! JMM bug fix : if cdep present, then chose attribute from cldepref
       istatus = copyatt(TRIM(cldepvar), nid_dep,incid,icout)
    ENDIF

    istatus = NF90_DEF_VAR(icout,cn_vtimec,NF90_DOUBLE,(/nid_t/), nid_tim)
    istatus = copyatt(cn_vtimec, nid_tim,incid,icout)

    ! Add Global General attribute at first call
    istatus=NF90_PUT_ATT(icout,NF90_GLOBAL,'start_date', nstart_date )
    istatus=NF90_PUT_ATT(icout,NF90_GLOBAL,'output_frequency', cfreq )
    istatus=NF90_PUT_ATT(icout,NF90_GLOBAL,'CONFIG',          config )
    istatus=NF90_PUT_ATT(icout,NF90_GLOBAL,'CASE',             ccase )

    istatus = NF90_CLOSE(incid)

    create=icout
  END FUNCTION create


  INTEGER(KIND=4) FUNCTION createvar(kout, sdtyvar, kvar, kpk, kidvo, cdglobal, ld_nc4)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION createvar  ***
    !!
    !! ** Purpose :  Create kvar  variables cdvar(:), in file id kout,
    !!
    !! ** Method  : INPUT:
    !!                 kout  = ncid of output file
    !!                 cdvar = array of name of variables
    !!                 kvar  = number of variables to create
    !!                 kpk   = number of vertical dimensions foreach variable
    !!             OUTPUT:
    !!                 kidvo = arrays with the varid of the variables just created.
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),                  INTENT(in   ) :: kout    ! ncid of output file
    TYPE (variable), DIMENSION(kvar) ,INTENT(inout) :: sdtyvar ! variable structure
    INTEGER(KIND=4),                  INTENT(in   ) :: kvar    ! number of variable
    INTEGER(KIND=4), DIMENSION(kvar), INTENT(in   ) :: kpk     ! number of level/var
    INTEGER(KIND=4), DIMENSION(kvar), INTENT(out  ) :: kidvo   ! varid's of output var
    CHARACTER(LEN=*), OPTIONAL,       INTENT(in   ) :: cdglobal! Global Attribute
    LOGICAL         , OPTIONAL,       INTENT(in   ) :: ld_nc4  ! user chunking and deflation

    INTEGER(KIND=4)               :: jv             ! dummy loop index
    INTEGER(KIND=4)               :: idims, istatus 
    INTEGER(KIND=4), DIMENSION(4) :: iidims, ichunk
    INTEGER(KIND=4)               :: iprecision
    LOGICAL                       :: ll_nc4
    !!----------------------------------------------------------------------
    IF ( PRESENT (ld_nc4 ) ) THEN 
       ll_nc4 = ld_nc4
    ELSE
       ll_nc4 = .false. 
    ENDIF

    DO jv = 1, kvar
       ! Create variables whose name is not 'none'
       IF ( sdtyvar(jv)%cname /= 'none' ) THEN
          IF (kpk(jv) == 1 ) THEN
             idims=3
             iidims(1) = nid_x ; iidims(2) = nid_y ; iidims(3) = nid_t
          ELSE IF (kpk(jv) > 1 ) THEN
             idims=4
             iidims(1) = nid_x ; iidims(2) = nid_y ; iidims(3) = nid_z ; iidims(4) = nid_t
          ELSE IF (kpk(jv) < 0 ) THEN
             idims=1
             iidims(1) = nid_t 
          ELSE
             PRINT *,' ERROR: ipk = ',kpk(jv), jv , sdtyvar(jv)%cname
             STOP 98
          ENDIF
    
          SELECT CASE ( sdtyvar(jv)%cprecision ) ! check the precision of the variable to create
          !
          CASE ( 'r8' ) ; iprecision = NF90_DOUBLE
          !
          CASE ( 'i2' ) ; iprecision = NF90_SHORT
          !
          CASE ( 'by' ) ; iprecision = NF90_BYTE
          !
          CASE DEFAULT  ! r4
                          iprecision = NF90_FLOAT
             IF ( sdtyvar(jv)%scale_factor /= 1. .OR. sdtyvar(jv)%add_offset /= 0. ) THEN
                          iprecision = NF90_SHORT
                          sdtyvar(jv)%cprecision='i2' 
             ENDIF
          END SELECT

#if defined key_netcdf4
         IF ( ll_nc4 ) THEN
          istatus = NF90_DEF_VAR(kout, sdtyvar(jv)%cname, iprecision, iidims(1:idims) ,kidvo(jv), & 
                  &               chunksizes=sdtyvar(jv)%ichunk(1:idims), deflate_level=1 )
         ELSE
          istatus = NF90_DEF_VAR(kout, sdtyvar(jv)%cname, iprecision, iidims(1:idims) ,kidvo(jv) )
         ENDIF
#else
          istatus = NF90_DEF_VAR(kout, sdtyvar(jv)%cname, iprecision, iidims(1:idims) ,kidvo(jv) )
#endif
          ! add attributes
          istatus = putatt(sdtyvar(jv), kout, kidvo(jv), cdglobal=cdglobal)
          createvar=istatus
       ENDIF
    END DO
    istatus = NF90_ENDDEF(kout)

  END FUNCTION createvar


  FUNCTION getvarid( cdfile, knvars )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION getvarid  ***
    !!
    !! ** Purpose :  return a real array with the nvar variable id 
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),       INTENT(in) :: cdfile
    INTEGER(KIND=4),        INTENT(in) :: knvars    ! Number of variables in cdfile
    INTEGER(KIND=4), DIMENSION(knvars) :: getvarid  ! return function

    INTEGER(KIND=4)                       :: jv     ! dummy loop index
    CHARACTER(LEN=256), DIMENSION(knvars) :: cdvar
    INTEGER(KIND=4)                       :: incid
    INTEGER(KIND=4)                       :: istatus
    !!----------------------------------------------------------------------
    istatus = NF90_OPEN(cdfile, NF90_NOWRITE, incid)
    DO jv = 1, knvars
       istatus = NF90_INQUIRE_VARIABLE(incid, jv, cdvar(jv) )
       istatus = NF90_INQ_VARID(incid, cdvar(jv), getvarid(jv))
    ENDDO
    istatus=NF90_CLOSE(incid)

  END FUNCTION getvarid


  INTEGER(KIND=4) FUNCTION getvaratt (cdfile, cdvar, cdunits, pmissing_value, cdlong_name, cdshort_name)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION getvaratt  ***
    !!
    !! ** Purpose : Get specific attributes for a variable (units, missing_value, 
    !!              long_name, short_name
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in)  :: cdfile, cdvar
    REAL(KIND=4), INTENT(out)       :: pmissing_value
    CHARACTER(LEN=*), INTENT(out) :: cdunits, cdlong_name, cdshort_name

    INTEGER(KIND=4) :: istatus
    INTEGER(KIND=4) :: incid, ivarid
    !!----------------------------------------------------------------------
    istatus = NF90_OPEN(cdfile, NF90_NOWRITE, incid)
    istatus = NF90_INQ_VARID(incid, cdvar, ivarid)

    istatus = NF90_GET_ATT(incid, ivarid, 'units',           cdunits        )
    pmissing_value = getspval ( cdfile, cdvar )
    istatus = NF90_GET_ATT(incid, ivarid, 'long_name',       cdlong_name    )
    istatus = NF90_GET_ATT(incid, ivarid, 'short_name',      cdshort_name   )

    getvaratt = istatus
    istatus   = NF90_CLOSE(incid)

  END FUNCTION getvaratt


  INTEGER(KIND=4) FUNCTION gettimeatt (cdfile, cdvartime, ctcalendar, cttitle, &
        &                            ctlong_name, ctaxis, ctunits, cttime_origin )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION gettimeatt  ***
    !!
    !! ** Purpose : Get specific attributes for time variable
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in)  :: cdfile
    CHARACTER(LEN=*), INTENT(in)  :: cdvartime
    CHARACTER(LEN=*), INTENT(out) :: ctcalendar, cttitle, ctlong_name, ctaxis, ctunits, cttime_origin

    INTEGER(KIND=4) :: istatus
    INTEGER(KIND=4) :: incid, ivarid
    !!----------------------------------------------------------------------
    ctcalendar    = 'unknown'
    cttitle       = 'unknown'
    ctlong_name   = 'unknown'
    ctaxis        = 'unknown'
    ctunits       = 'unknown'
    cttime_origin = 'unknown'

    istatus = NF90_OPEN(cdfile, NF90_NOWRITE, incid)
    istatus = NF90_INQ_VARID(incid, cdvartime, ivarid)

    istatus = NF90_GET_ATT(incid, ivarid, 'calendar',         ctcalendar    )
    istatus = NF90_GET_ATT(incid, ivarid, 'title',            cttitle       )
    istatus = NF90_GET_ATT(incid, ivarid, 'long_name',        ctlong_name   )
    istatus = NF90_GET_ATT(incid, ivarid, 'axis',             ctaxis        )
    istatus = NF90_GET_ATT(incid, ivarid, 'units',            ctunits       )
    istatus = NF90_GET_ATT(incid, ivarid, 'time_origin',      cttime_origin )

    gettimeatt = istatus
    istatus   = NF90_CLOSE(incid)

  END FUNCTION gettimeatt

  INTEGER(KIND=4) FUNCTION puttimeatt (kout, cdvartime, ctcalendar, cttitle, &
        &                          ctlong_name, ctaxis, ctunits, cttime_origin )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION puttimeatt  ***
    !!
    !! ** Purpose : Put specific attributes for time variable
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),     INTENT(in) :: kout
    CHARACTER(LEN=20),  INTENT(in)  :: cdvartime
    CHARACTER(LEN=256), INTENT(in) :: ctcalendar, cttitle, ctlong_name, ctaxis, ctunits, cttime_origin

    INTEGER(KIND=4) :: ivarid
    !!----------------------------------------------------------------------

    puttimeatt=NF90_INQ_VARID(kout, cdvartime, ivarid)
    IF (puttimeatt /= 0 ) THEN 
      PRINT *, TRIM(NF90_STRERROR(puttimeatt) )
      PRINT *, 'puttimeatt : ',TRIM(cdvartime),' does not exist'
      STOP 98
    ENDIF
   
    puttimeatt = NF90_REDEF(kout)
    puttimeatt=NF90_PUT_ATT(kout,ivarid,'calendar',ctcalendar)
    IF (puttimeatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(puttimeatt); PRINT *,' in puttimeatt calendar'     ; STOP 98 ; ENDIF
    puttimeatt=NF90_PUT_ATT(kout,ivarid,'title',cttitle)
    IF (puttimeatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(puttimeatt); PRINT *,' in puttimeatt title'        ; STOP 98 ; ENDIF
    puttimeatt=NF90_PUT_ATT(kout,ivarid,'long_name',ctlong_name)
    IF (puttimeatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(puttimeatt); PRINT *,' in puttimeatt long_name'    ; STOP 98 ; ENDIF
    puttimeatt=NF90_PUT_ATT(kout,ivarid,'axis',ctaxis)
    IF (puttimeatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(puttimeatt); PRINT *,' in puttimeatt axis'         ; STOP 98 ; ENDIF
    puttimeatt=NF90_PUT_ATT(kout,ivarid,'units',ctunits)
    IF (puttimeatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(puttimeatt); PRINT *,' in puttimeatt units'        ; STOP 98 ; ENDIF
    puttimeatt=NF90_PUT_ATT(kout,ivarid,'time_origin',cttime_origin)
    IF (puttimeatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(puttimeatt); PRINT *,' in puttimeatt time_origin'  ; STOP 98 ; ENDIF

    puttimeatt=NF90_ENDDEF(kout)

  END FUNCTION puttimeatt


  INTEGER(KIND=4) FUNCTION cvaratt (cdfile, cdvar, cdunits, pmissing_value, cdlong_name, cdshort_name)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION cvaratt  ***
    !!
    !! ** Purpose : Change variable attributs in an existing variable
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=256), INTENT(in) :: cdfile, cdvar
    CHARACTER(LEN=256), INTENT(in) :: cdunits, cdlong_name, cdshort_name
    REAL(KIND=4),       INTENT(in) :: pmissing_value

    INTEGER(KIND=4) :: istatus
    INTEGER(KIND=4) :: incid, ivarid
    REAL(KIND=4)    :: zspval
    CHARACTER(LEN=256) :: clmissing   ! get the actual missing_value attribute name
    !!----------------------------------------------------------------------
    istatus = NF90_OPEN (cdfile, NF90_WRITE, incid)
    istatus = NF90_REDEF(incid)
    istatus = NF90_INQ_VARID(incid, cdvar, ivarid)

!   istatus=NF90_RENAME_ATT(incid, ivarid, 'units',         cdunits        )
    istatus=NF90_PUT_ATT(incid, ivarid, 'units',         cdunits        )
    zspval = getspval   ( cdfile, cdvar, clmissing                       )
    istatus=NF90_PUT_ATT(incid, ivarid, clmissing, pmissing_value        )
    istatus=NF90_PUT_ATT(incid, ivarid, 'long_name',     cdlong_name    )
    istatus=NF90_PUT_ATT(incid, ivarid, 'short_name',    cdshort_name   )

    istatus=NF90_ENDDEF(incid)
    cvaratt=istatus
    istatus=NF90_CLOSE(incid)

  END FUNCTION cvaratt


  INTEGER(KIND=4) FUNCTION putatt (sdtyvar, kout, kid, cdglobal)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION putatt  ***
    !!
    !! ** Purpose : Put attribute for variable defined in the data structure 
    !!
    !!----------------------------------------------------------------------
    TYPE (variable),            INTENT(in) :: sdtyvar
    INTEGER(KIND=4),            INTENT(in) :: kout, kid
    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: cdglobal   !: global attribute
    !!----------------------------------------------------------------------
    putatt=NF90_PUT_ATT(kout,kid,'units',sdtyvar%cunits) 
    IF (putatt /= NF90_NOERR ) THEN ;PRINT *, NF90_STRERROR(putatt); PRINT *,' putatt for units'  ; STOP 98; ENDIF

    ! With netcdf4, missing value must have the same precision than the variable. Need to convert
    ! to sdtyvar%cprecision previous PUT_ATT
    SELECT CASE (sdtyvar%cprecision )
    CASE ( 'r8' ) ; putatt=NF90_PUT_ATT(kout,kid,cn_missing_value,REAL(sdtyvar%rmissing_value,8) )  
    CASE ( 'i2' ) ; putatt=NF90_PUT_ATT(kout,kid,cn_missing_value, INT(sdtyvar%rmissing_value,2) )  
    CASE ( 'by' ) ; putatt=NF90_PUT_ATT(kout,kid,cn_missing_value, INT(sdtyvar%rmissing_value,1) )  
    CASE DEFAULT  ; putatt=NF90_PUT_ATT(kout,kid,cn_missing_value,REAL(sdtyvar%rmissing_value,4) )
    END SELECT
    IF (putatt /= NF90_NOERR ) THEN ;PRINT *, NF90_STRERROR(putatt); PRINT *,' putatt for missing value'    ; STOP 98; ENDIF

    putatt=NF90_PUT_ATT(kout,kid,'valid_min',sdtyvar%valid_min) 
    IF (putatt /= NF90_NOERR ) THEN ;PRINT *, NF90_STRERROR(putatt); PRINT *,' putatt for valid_min'        ; STOP 98; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'valid_max',sdtyvar%valid_max)
    IF (putatt /= NF90_NOERR ) THEN ;PRINT *, NF90_STRERROR(putatt); PRINT *,' putatt for valid_max'        ; STOP 98; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'long_name',sdtyvar%clong_name)
    IF (putatt /= NF90_NOERR ) THEN ;PRINT *, NF90_STRERROR(putatt); PRINT *,' putatt for long_name'        ; STOP 98; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'short_name',sdtyvar%cshort_name) 
    IF (putatt /= NF90_NOERR ) THEN ;PRINT *, NF90_STRERROR(putatt); PRINT *,' putatt for short_name'       ; STOP 98; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'iweight',sdtyvar%iwght) 
    IF (putatt /= NF90_NOERR ) THEN ;PRINT *, NF90_STRERROR(putatt); PRINT *,' putatt for iweight'          ; STOP 98; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'online_operation',sdtyvar%conline_operation) 
    IF (putatt /= NF90_NOERR ) THEN ;PRINT *, NF90_STRERROR(putatt); PRINT *,' putatt for online_operation' ; STOP 98; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'axis',sdtyvar%caxis) 
    IF (putatt /= NF90_NOERR ) THEN ;PRINT *, NF90_STRERROR(putatt); PRINT *,' putatt for axis'             ; STOP 98; ENDIF

    ! Optional attributes (scale_factor, add_offset )
    putatt=NF90_PUT_ATT(kout,kid,'scale_factor',sdtyvar%scale_factor) 
    IF (putatt /= NF90_NOERR ) THEN ;PRINT *, NF90_STRERROR(putatt); PRINT *,' putatt for scale_factor'     ; STOP 98; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'add_offset',sdtyvar%add_offset) 
    IF (putatt /= NF90_NOERR ) THEN ;PRINT *, NF90_STRERROR(putatt); PRINT *,' putatt for add_offset'       ; STOP 98; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'savelog10',sdtyvar%savelog10) 
    IF (putatt /= NF90_NOERR ) THEN ;PRINT *, NF90_STRERROR(putatt); PRINT *,' putatt for savelog10'        ; STOP 98; ENDIF

    ! Global attribute
    IF ( PRESENT(cdglobal) ) THEN
      putatt=NF90_PUT_ATT(kout,NF90_GLOBAL,'history',cdglobal)
      IF (putatt /= NF90_NOERR ) THEN ;PRINT *, NF90_STRERROR(putatt); PRINT *,' putatt for history'        ; STOP 98; ENDIF
    ENDIF

  END FUNCTION putatt


  REAL(KIND=4) FUNCTION getattr (cdfile, cdvar, cdatt)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION getatt  ***
    !!
    !! ** Purpose : return a REAL value with the values of the
    !!              attribute cdatt for all the variable cdvar in cdfile  
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdfile  ! file name
    CHARACTER(LEN=*), INTENT(in) :: cdvar   ! var name
    CHARACTER(LEN=*), INTENT(in) :: cdatt   ! attribute name to look for

    INTEGER(KIND=4) :: istatus, jv, incid, idum
    !!----------------------------------------------------------------------
    istatus = NF90_OPEN  (cdfile, NF90_NOWRITE, incid)
    istatus = NF90_INQ_VARID(incid, cdvar, idum)

    IF ( istatus /= NF90_NOERR) PRINT *, TRIM(NF90_STRERROR(istatus)),' when looking for ',TRIM(cdvar),' in getatt.'

    istatus = NF90_GET_ATT(incid, idum, cdatt, getattr)
    IF ( istatus /= NF90_NOERR ) THEN
       PRINT *,' getatt problem :',NF90_STRERROR(istatus)
       PRINT *,' attribute :', TRIM(cdatt)
       PRINT *,' variable  :', TRIM(cdvar)
       PRINT *,' file      :', TRIM(cdfile)
       PRINT *,' return default 0 '
       getattr=0.
    ENDIF

    istatus=NF90_CLOSE(incid)

  END FUNCTION getattr

  INTEGER(KIND=4) FUNCTION getatti (cdfile, cdvar, cdatt, kidum)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION getatti  ***
    !!
    !! ** Purpose : return a INTEGER value with the values of the
    !!              attribute cdatt for all the variable cdvar in cdfile  
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdfile  ! file name
    CHARACTER(LEN=*), INTENT(in) :: cdvar   ! var name
    CHARACTER(LEN=*), INTENT(in) :: cdatt   ! attribute name to look for
    INTEGER(KIND=4)     :: kidum

    INTEGER(KIND=4) :: istatus, jv, incid, idum
    !!----------------------------------------------------------------------
    istatus = NF90_OPEN  (cdfile, NF90_NOWRITE, incid)
    istatus = NF90_INQ_VARID(incid, cdvar, idum)

    IF ( istatus /= NF90_NOERR) PRINT *, TRIM(NF90_STRERROR(istatus)),' when looking for ',TRIM(cdvar),' in getatt.'

    istatus = NF90_GET_ATT(incid, idum, cdatt, getatti)
    IF ( istatus /= NF90_NOERR ) THEN
       PRINT *,' getatt problem :',NF90_STRERROR(istatus)
       PRINT *,' attribute :', TRIM(cdatt)
       PRINT *,' variable  :', TRIM(cdvar)
       PRINT *,' file      :', TRIM(cdfile)
       PRINT *,' return default 0 '
       getatti=0
    ENDIF

    istatus=NF90_CLOSE(incid)

  END FUNCTION getatti


  CHARACTER(LEN=256) FUNCTION getattc (cdfile, cdvar, cdatt, cdum)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION getattr  ***
    !!
    !! ** Purpose : return a string value with the values of the
    !!              attribute cdatt for all the variable cdvar in cdfile  
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdfile  ! file name
    CHARACTER(LEN=*), INTENT(in) :: cdvar   ! var name
    CHARACTER(LEN=*), INTENT(in) :: cdatt   ! attribute name to look for

    CHARACTER(LEN=*)   :: cdum

    INTEGER(KIND=4) :: istatus, jv, incid, idum
    !!----------------------------------------------------------------------
    istatus = NF90_OPEN  (cdfile, NF90_NOWRITE, incid)
    istatus = NF90_INQ_VARID(incid, cdvar, idum)

    IF ( istatus /= NF90_NOERR) PRINT *, TRIM(NF90_STRERROR(istatus)),' when looking for ',TRIM(cdvar),' in getatt.'

    istatus = NF90_GET_ATT(incid, idum, cdatt, getattc)
    IF ( istatus /= NF90_NOERR ) THEN
       PRINT *,' getatt problem :',NF90_STRERROR(istatus)
       PRINT *,' attribute :', TRIM(cdatt)
       PRINT *,' variable  :', TRIM(cdvar)
       PRINT *,' file      :', TRIM(cdfile)
       PRINT *,' return default N/A '
       getattc="N/A"
    ENDIF

    istatus=NF90_CLOSE(incid)

  END FUNCTION getattc



  INTEGER(KIND=4) FUNCTION atted_char ( cdfile, cdvar, cdatt, cdvalue )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION atted_char  ***
    !!
    !! ** Purpose : attribute editor : modify existing attribute or create
    !!              new attribute for variable cdvar in cdfile 
    !!
    !! ** Method  : just put_att after some check.
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),  INTENT(in) :: cdfile  ! input file
    CHARACTER(LEN=*),  INTENT(in) :: cdvar   ! variable name
    CHARACTER(LEN=*),  INTENT(in) :: cdatt   ! attribute  name
    CHARACTER(LEN=*),  INTENT(in) :: cdvalue ! attribute value

    INTEGER(KIND=4)               :: incid,  istatus, idvar, ityp
    !!-------------------------------------------------------------------------
    istatus = NF90_OPEN(cdfile, NF90_WRITE, incid)
    istatus = NF90_INQ_VARID(incid, cdvar, idvar)
    IF ( istatus /= NF90_NOERR ) THEN
       PRINT *, NF90_STRERROR(istatus),' in atted ( inq_varid)'
       STOP 98
    ENDIF
    istatus = NF90_INQUIRE_ATTRIBUTE(incid, idvar, cdatt, xtype=ityp )
    IF ( istatus /= NF90_NOERR ) THEN
       PRINT *, ' Attribute does not exist. Create it'
       istatus = NF90_REDEF(incid)
       istatus = NF90_PUT_ATT(incid, idvar, cdatt, cdvalue)
       atted_char = istatus
    ELSE
       IF ( ityp == NF90_CHAR ) THEN
         istatus = NF90_REDEF(incid)
         istatus = NF90_PUT_ATT(incid, idvar, cdatt, cdvalue)
         atted_char = istatus
       ELSE
         PRINT *, ' Mismatch in attribute type in atted_char'
         STOP 98
       ENDIF
    ENDIF
    istatus=NF90_CLOSE(incid)

  END FUNCTION atted_char


  INTEGER(KIND=4) FUNCTION atted_r4 ( cdfile, cdvar, cdatt, pvalue )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION atted_r4  ***
    !!
    !! ** Purpose : attribute editor : modify existing attribute or create
    !!              new attribute for variable cdvar in cdfile
    !!
    !! ** Method  : just put_att after some check.
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),  INTENT(in) :: cdfile  ! input file
    CHARACTER(LEN=*),  INTENT(in) :: cdvar   ! variable name
    CHARACTER(LEN=*),  INTENT(in) :: cdatt   ! attribute  name
    REAL(KIND=4),      INTENT(in) :: pvalue  ! attribute value

    INTEGER(KIND=4)               :: incid,  istatus, idvar, ityp
    !!-------------------------------------------------------------------------
    istatus = NF90_OPEN(cdfile, NF90_WRITE, incid)
    istatus = NF90_INQ_VARID(incid, cdvar, idvar)
    IF ( istatus /= NF90_NOERR ) THEN
       PRINT *, NF90_STRERROR(istatus),' in atted ( inq_varid)'
       STOP 98
    ENDIF
    istatus = NF90_INQUIRE_ATTRIBUTE(incid, idvar, cdatt, xtype=ityp )
    IF ( istatus /= NF90_NOERR ) THEN
       PRINT *, ' Attribute does not exist. Create it'
       istatus = NF90_REDEF(incid)
       istatus = NF90_PUT_ATT(incid, idvar, cdatt, pvalue)
       atted_r4 = istatus
    ELSE
       IF ( ityp == NF90_FLOAT ) THEN
         istatus = NF90_REDEF(incid)
         istatus = NF90_PUT_ATT(incid, idvar, cdatt, pvalue)
         atted_r4 = istatus
       ELSE
         PRINT *, ' Mismatch in attribute type in atted_r4'
         STOP 98
       ENDIF
    ENDIF
    istatus=NF90_CLOSE(incid)

  END FUNCTION atted_r4


  INTEGER(KIND=4)  FUNCTION  getdim (cdfile, cdim_name, cdtrue, kstatus, ldexact)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION getdim  ***
    !!
    !! ** Purpose : Return the INTEGER value of the dimension
    !!              identified with cdim_name in cdfile 
    !!
    !! ** Method  : This function look for a dimension name that contains 
    !!              cdim_name, in cdfile. In option it returns the error 
    !!              status which can be used to make another intent, changing 
    !!              the dim name. Finally, with the last optional argument 
    !!              ldexact, exact match to cdim_name can be required.
    !!              
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),             INTENT(in ) :: cdfile     ! File name to look at
    CHARACTER(LEN=*),             INTENT(in ) :: cdim_name  ! File name to look at
    CHARACTER(LEN=256), OPTIONAL, INTENT(out) :: cdtrue     ! full name of the read dimension
    INTEGER(KIND=4),    OPTIONAL, INTENT(out) :: kstatus    ! status of the nf inquire
    LOGICAL,            OPTIONAL, INTENT(in ) :: ldexact    ! when true look for exact cdim_name

    INTEGER(KIND=4)    :: jdim
    INTEGER(KIND=4)    :: incid, id_dim, idv
    INTEGER(KIND=4)    :: istatus
    INTEGER(KIND=4)    :: idims
    CHARACTER(LEN=256) :: clnam
    LOGICAL            :: lexact   = .false.
    LOGICAL, SAVE      :: ll_first = .true.
    !!-----------------------------------------------------------
    clnam = '-------------'

    IF ( PRESENT(kstatus)  ) kstatus=0
    IF ( PRESENT(ldexact)  ) lexact=ldexact
    IF ( cdim_name == cn_x ) lexact=.true.  ! fix for XIOS files having now a new dimension xaxis_bound which match getdim ('x') ....
                                            ! more clever fix must be found for identification of the dimensions in the input files
    istatus=NF90_OPEN(cdfile, NF90_NOWRITE, incid)
    IF ( istatus == NF90_NOERR ) THEN
       istatus=NF90_INQUIRE(incid, ndimensions=idims)

       IF ( lexact ) THEN
          istatus=NF90_INQ_DIMID(incid, cdim_name, id_dim)
          IF (istatus /= NF90_NOERR ) THEN
            PRINT *,NF90_STRERROR(istatus)
            PRINT *,' Exact dimension name ', TRIM(cdim_name),' not found in ',TRIM(cdfile) ; STOP 98
          ENDIF
          istatus=NF90_INQUIRE_DIMENSION(incid, id_dim, len=getdim)
          IF ( PRESENT(cdtrue) ) cdtrue=cdim_name
          jdim = 0
       ELSE  ! scann all dims to look for a partial match
         DO jdim = 1, idims
            istatus=NF90_INQUIRE_DIMENSION(incid, jdim, name=clnam, len=getdim)
            IF ( INDEX(clnam, TRIM(cdim_name)) /= 0 ) THEN
               IF ( PRESENT(cdtrue) ) cdtrue=clnam
               EXIT
            ENDIF
         ENDDO
       ENDIF

       IF ( jdim > idims ) THEN   ! dimension not found
          IF ( PRESENT(kstatus) ) kstatus=1    ! error send optionally to the calling program
          getdim=0
          IF ( PRESENT(cdtrue) ) cdtrue='unknown'
       ENDIF
       ! first call 
       IF ( ll_first ) THEN
          ll_first = .false. 
          ! take the opportunity to initialize time_counter attributes :
          istatus = NF90_INQ_VARID(incid, cn_vtimec, idv )
          IF ( istatus == NF90_NOERR ) THEN 
             istatus = NF90_GET_ATT(incid, idv, 'units',  ctime_units  )
             istatus = NF90_GET_ATT(incid, idv, 'time_origin', ctime_origin )
             istatus = NF90_GET_ATT(incid, idv, 'calendar',    calendar     )
          ENDIF
 
          ! read global attributes 

          ! start_date
          cl_dum = Get_Env ( 'start_date' )  ! look for environment variable
          istatus = NF90_INQUIRE_ATTRIBUTE(incid, NF90_GLOBAL, 'start_date')
          IF ( istatus == NF90_NOERR ) THEN
             istatus = NF90_GET_ATT(incid, NF90_GLOBAL, 'start_date', nstart_date )
          ELSE IF ( cl_dum /= '' ) THEN 
             READ(cl_dum, * ) nstart_date
          ELSE
             nstart_date = -1
          ENDIF

          ! output_frequency 
          cl_dum = Get_Env ( 'output_frequency' )  ! look for environment variable
          istatus = NF90_INQUIRE_ATTRIBUTE(incid, NF90_GLOBAL, 'output_frequency')
          IF ( istatus == NF90_NOERR ) THEN
             istatus = NF90_GET_ATT(incid, NF90_GLOBAL, 'output_frequency', cfreq )
          ELSE IF ( cl_dum /= '' ) THEN
             cfreq = TRIM(cl_dum)
          ELSE
             cfreq = 'N/A'
          ENDIF

          ! CONFIG
          cl_dum = Get_Env ( 'CONFIG' )  ! look for environment variable
          istatus = NF90_INQUIRE_ATTRIBUTE(incid, NF90_GLOBAL, 'CONFIG')
          IF ( istatus == NF90_NOERR ) THEN
             istatus = NF90_GET_ATT(incid, NF90_GLOBAL, 'CONFIG', config )
          ELSE IF ( cl_dum /= '' ) THEN
             config = TRIM(cl_dum)
          ELSE
             config = 'N/A'
          ENDIF

          ! CASE
          cl_dum = Get_Env ( 'CASE' )  ! look for environment variable
          istatus = NF90_INQUIRE_ATTRIBUTE(incid, NF90_GLOBAL, 'CASE')
          IF ( istatus == NF90_NOERR ) THEN
             istatus = NF90_GET_ATT(incid, NF90_GLOBAL, 'CASE', ccase )
          ELSE IF ( cl_dum /= '' ) THEN
             ccase = TRIM(cl_dum)
          ELSE
             ccase = 'N/A'
          ENDIF
       ENDIF
       istatus=NF90_CLOSE(incid)
    ELSE              ! problem with the file
       IF ( PRESENT(cdtrue) ) cdtrue='unknown'
       IF ( PRESENT(kstatus) ) kstatus=1 
    ENDIF
    ! reset lexact to false for next call 
    lexact=.false.

  END FUNCTION getdim


  REAL(KIND=4) FUNCTION  getspval (cdfile, cdvar, cdmissing )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION getspval  ***
    !!
    !! ** Purpose : return the SPVAL value of the variable cdvar in cdfile
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),           INTENT(in ) :: cdfile    ! File name to look at
    CHARACTER(LEN=*),           INTENT(in ) :: cdvar     ! variable name
    CHARACTER(LEN=*), OPTIONAL, INTENT(out) :: cdmissing ! missing att. name

    INTEGER(KIND=4) :: incid, id_var
    INTEGER(KIND=4) :: istatus
    INTEGER(KIND=4) :: jtry
    !!----------------------------------------------------------------------

    IF ( PRESENT (cdmissing) ) cdmissing = cn_missing_value

    istatus=NF90_OPEN      (cdfile, NF90_NOWRITE, incid               )
    istatus=NF90_INQ_VARID (incid, cdvar, id_var                      )
    istatus=NF90_GET_ATT   (incid, id_var, cn_missing_value, getspval )

    IF ( istatus /= NF90_NOERR ) THEN 
      DO jtry = 1, jp_missing_nm
         IF ( PRESENT (cdmissing) ) cdmissing = TRIM(cl_missing_nm(jtry))
         istatus = NF90_GET_ATT (incid, id_var, cl_missing_nm(jtry) , getspval )
         IF ( istatus == NF90_NOERR ) EXIT
         IF ( PRESENT (cdmissing) ) cdmissing = cn_missing_value
         getspval = 0.
      ENDDO
    ENDIF
    istatus=NF90_CLOSE     (incid )

  END FUNCTION getspval


  INTEGER(KIND=4) FUNCTION getvdim (cdfile, cdvar)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION getvdim  ***
    !!
    !! ** Purpose : Return the number of dimensions for variable cdvar in cdfile 
    !!
    !! ** Method  : Inquire for variable cdvar in cdfile. If found,
    !!              determines the number of dimensions , assuming that variables
    !!              are either (x,y,dep,time) or (x,y,time)
    !!              If cdvar is not found, give an interactive choice for an existing
    !!              variable, cdvar is then updated to this correct name.  
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in)    :: cdfile   ! File name to look at
    CHARACTER(LEN=*), INTENT(inout) :: cdvar    ! variable name to look at.

    INTEGER(KIND=4)    :: jvar
    INTEGER(KIND=4)    :: istatus, incid, id_var, ivar, idi, istatus0
    CHARACTER(LEN=256) :: clongname='long_name', clongn
    !!----------------------------------------------------------------------
    CALL ERR_HDL(NF90_OPEN(cdfile,NF90_NOWRITE,incid))

    istatus0 = NF90_INQ_VARID ( incid,cdvar,id_var)
    DO WHILE  ( istatus0 == NF90_ENOTVAR ) 
       ivar=getnvar(cdfile)
       PRINT *, 'Give the number corresponding to the variable you want to work with '
       DO jvar = 1, ivar
          clongn=''
          istatus=NF90_INQUIRE_VARIABLE (incid, jvar, cdvar, ndims=idi)
          istatus=NF90_GET_ATT (incid, jvar, clongname, clongn)
          IF (istatus /= NF90_NOERR ) clongn='unknown'
          PRINT *, jvar, ' ',TRIM(cdvar),' ',TRIM(clongn)
       ENDDO
       READ *,id_var
       istatus0=NF90_INQUIRE_VARIABLE (incid, id_var, cdvar, ndims=idi)
    ENDDO
    ! 
    CALL ERR_HDL(NF90_INQUIRE_VARIABLE (incid, id_var, cdvar, ndims=idi))
    getvdim = idi - 1
    CALL ERR_HDL (NF90_CLOSE(incid))

  END FUNCTION getvdim

  FUNCTION getdimvar (cdfile, kpk, cd_depnam)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION getdimvar  ***
    !!
    !! ** Purpose :  Try to infer the name of the depth variable associated
    !!               with z dimension, and return the values as a 1d array 
    !!
    !! ** Method  :  Trial and error method ... 
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdfile
    INTEGER(KIND=4),  INTENT(in) :: kpk
    REAL(KIND=4), DIMENSION(kpk) :: getdimvar
    CHARACTER(LEN=*), OPTIONAL, INTENT(out) :: cd_depnam

    INTEGER(KIND=4) :: ji
    INTEGER(KIND=4) :: incid, idims, iuldid, idimv, ivars, istatus
    INTEGER(KIND=4), DIMENSION(4) :: idimt
    CHARACTER(LEN=80) :: clvar ='none'
    
    !!----------------------------------------------------------------------
    istatus = NF90_OPEN   (cdfile,NF90_NOWRITE,incid)
    istatus = NF90_INQUIRE(incid, ndimensions=idims, unlimiteddimid=iuldid,&
                  &              nvariables=ivars)
    ! look for variables with only one dim, not unlimited dim ...
    DO ji = 1, ivars
       istatus = NF90_INQUIRE_VARIABLE(incid, ji, name=clvar, ndims=idimv, &
                  & dimids=idimt )
       IF ( idimv == 1 ) THEN ! candidate
          IF ( idimt(1) /= iuldid) THEN ! this is The Variable
             PRINT *, ' Found vertical variable :', TRIM(clvar)
             istatus = NF90_GET_VAR(incid, ji, getdimvar(:) )
             EXIT
          ENDIF
       ENDIF
    ENDDO
    IF ( present(cd_depnam) ) cd_depnam=TRIM(clvar)
    IF ( ji == ivars +1 ) THEN
       PRINT *,'Sorry, no vertical dim variables inferred ...'
       getdimvar=0.
    ENDIF

  END FUNCTION

  INTEGER(KIND=4) FUNCTION  getnvar (cdfile)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION getnvar  ***
    !!
    !! ** Purpose :  return the number of variables in cdfile 
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) ::  cdfile   ! file to look at

    INTEGER(KIND=4) :: incid
    INTEGER(KIND=4) :: istatus
    !!----------------------------------------------------------------------
    istatus = NF90_OPEN    (cdfile, NF90_NOWRITE, incid )
    istatus = NF90_INQUIRE (incid, nvariables = getnvar )
    istatus = NF90_CLOSE   (incid )

  END FUNCTION getnvar


  FUNCTION  getipk (cdfile,knvars,cdep)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION getipk  ***
    !!
    !! ** Purpose : Return the number of levels for all the variables
    !!              in cdfile. Return 0 if the variable in 1d.
    !!
    !! ** Method  : returns npk when 4D variables ( x,y,z,t )
    !!              returns  1  when 3D variables ( x,y,  t )
    !!              returns  0  when other ( vectors )
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),           INTENT(in) :: cdfile   ! File to look at
    INTEGER(KIND=4),            INTENT(in) ::  knvars  ! Number of variables in cdfile
    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: cdep     ! optional depth dim name
    INTEGER(KIND=4), DIMENSION(knvars)     :: getipk   ! array (variables ) of levels

    INTEGER(KIND=4)    :: incid, ipk, jv, iipk
    INTEGER(KIND=4)    :: istatus
    CHARACTER(LEN=256) :: cldep='dep'
    !!----------------------------------------------------------------------
    istatus=NF90_OPEN(cdfile,NF90_NOWRITE,incid)

    IF (  PRESENT (cdep) ) cldep = cdep

    ! Note the very important TRIM below : if not, getdim crashes as it never find the correct dim !
    iipk = getdim(cdfile, TRIM(cldep), kstatus=istatus)

    IF ( istatus /= 0 ) THEN
       PRINT *,' getipk : vertical dim not found ...assume 1'
       iipk=1
    ENDIF

    DO jv = 1, knvars
       istatus=NF90_INQUIRE_VARIABLE(incid, jv, ndims=ipk)
       IF (ipk == 4 ) THEN
          getipk(jv) = iipk
       ELSE IF (ipk == 3 ) THEN
          getipk(jv) = 1
       ELSE
          getipk(jv) = 0
       ENDIF
    END DO

    istatus=NF90_CLOSE(incid)

  END FUNCTION getipk


  FUNCTION  getvarname (cdfile, knvars, sdtypvar)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION getvarname  ***
    !!
    !! ** Purpose : return a character array with the knvars variable
    !!              name corresponding to cdfile 
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),          INTENT(in) :: cdfile
    INTEGER(KIND=4),           INTENT(in) :: knvars    ! Number of variables in cdfile
    TYPE (variable),   DIMENSION (knvars) :: sdtypvar  ! Retrieve variables attribute
    CHARACTER(LEN=256), DIMENSION(knvars) :: getvarname

    INTEGER(KIND=4)    :: incid,  jv, ilen
    INTEGER(KIND=4)    :: istatus
    INTEGER(KIND=4)    :: iatt
    REAL(KIND=4)       :: zatt
    CHARACTER(LEN=256) :: cldum=''
    !!----------------------------------------------------------------------
    istatus=NF90_OPEN(cdfile,NF90_NOWRITE,incid)

    DO jv = 1, knvars
       istatus=NF90_INQUIRE_VARIABLE(incid, jv, name=getvarname(jv) )
       sdtypvar(jv)%cname=getvarname(jv)


       ! look for standard attibutes
       IF ( NF90_INQUIRE_ATTRIBUTE(incid, jv, 'units', len=ilen) == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(incid, jv, 'units', cldum(1:ilen))
          sdtypvar(jv)%cunits = TRIM(cldum)
          cldum = ''
       ELSE 
          sdtypvar(jv)%cunits = 'N/A'
       ENDIF

       sdtypvar(jv)%rmissing_value = getspval ( cdfile, sdtypvar(jv)%cname )

       IF ( NF90_INQUIRE_ATTRIBUTE(incid, jv, 'valid_min') == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(incid, jv, 'valid_min', zatt)
          sdtypvar(jv)%valid_min = zatt
       ELSE
          sdtypvar(jv)%valid_min = 0.
       ENDIF

       IF ( NF90_INQUIRE_ATTRIBUTE(incid, jv, 'valid_max') == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(incid, jv, 'valid_max', zatt)
          sdtypvar(jv)%valid_max = zatt
       ELSE
          sdtypvar(jv)%valid_max = 0.
       ENDIF

       IF ( NF90_INQUIRE_ATTRIBUTE(incid, jv, 'iweight') == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(incid, jv, 'iweight', iatt)
          sdtypvar(jv)%iwght = iatt
       ELSE
          sdtypvar(jv)%iwght = 1
       ENDIF

       IF ( NF90_INQUIRE_ATTRIBUTE(incid, jv, 'long_name', len=ilen) == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(incid, jv, 'long_name', cldum(1:ilen))
          sdtypvar(jv)%clong_name = TRIM(cldum)
          cldum = ''
       ELSE
          sdtypvar(jv)%clong_name = 'N/A'
       ENDIF

       IF ( NF90_INQUIRE_ATTRIBUTE(incid, jv, 'short_name', len=ilen) == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(incid, jv, 'short_name', cldum(1:ilen))
          sdtypvar(jv)%cshort_name = TRIM(cldum)
          cldum = ''
       ELSE
          sdtypvar(jv)%cshort_name = 'N/A'
       ENDIF

       IF ( NF90_INQUIRE_ATTRIBUTE(incid, jv, 'online_operation', len=ilen) == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(incid, jv, 'online_operation', cldum(1:ilen))
          sdtypvar(jv)%conline_operation = TRIM(cldum)
          cldum = ''
       ELSE
          sdtypvar(jv)%conline_operation = 'N/A'
       ENDIF

       IF ( NF90_INQUIRE_ATTRIBUTE(incid, jv, 'axis', len=ilen) == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(incid, jv, 'axis', cldum(1:ilen))
          sdtypvar(jv)%caxis = TRIM(cldum)
          cldum = ''
       ELSE
          sdtypvar(jv)%caxis = 'N/A'
       ENDIF

    END DO
    istatus=NF90_CLOSE(incid)

  END FUNCTION getvarname

  FUNCTION  getvar (cdfile,cdvar,klev,kpi,kpj,kimin,kjmin, ktime, ldiom)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION  getvar  ***
    !!
    !! ** Purpose : Return the 2D REAL variable cvar, from cdfile at level klev.
    !!              kpi,kpj are the horizontal size of the 2D variable
    !!
    !! ** Method  : Initially a quite straigth forward function. But with the
    !!              NEMO variation about the e3t in partial steps, I try to adapt
    !!              the code to all existing mesh_zgr format, which reduces the
    !!              readibility of the code. One my think of specific routine for
    !!              getvar (e3._ps ...)
    !!
    !!---------------------------------------------------------------------
    CHARACTER(LEN=*),          INTENT(in) :: cdfile       ! file name to work with
    CHARACTER(LEN=*),          INTENT(in) :: cdvar        ! variable name to work with
    INTEGER(KIND=4), OPTIONAL, INTENT(in) :: klev         ! Optional variable. If missing 1 is assumed
    INTEGER(KIND=4),           INTENT(in) :: kpi, kpj     ! horizontal size of the 2D variable
    INTEGER(KIND=4), OPTIONAL, INTENT(in) :: kimin, kjmin ! Optional variable. If missing 1 is assumed
    INTEGER(KIND=4), OPTIONAL, INTENT(in) :: ktime        ! Optional variable. If missing 1 is assumed
    LOGICAL,         OPTIONAL, INTENT(in) :: ldiom        ! Optional variable. If missing false is assumed
    REAL(KIND=4), DIMENSION(kpi,kpj) :: getvar            ! 2D REAL 4 holding variable field at klev

    INTEGER(KIND=4), DIMENSION(4)               :: istart, icount, inldim
    INTEGER(KIND=4)                             :: incid, id_var, id_dimunlim, inbdim, inbdim2
    INTEGER(KIND=4)                             :: istatus, ilev, imin, jmin
    INTEGER(KIND=4)                             :: itime, ilog, ipiglo, imax
    INTEGER(KIND=4), SAVE                       :: ii, ij, ik0, ji, jj, ik1, ik
    REAL(KIND=4)                                :: sf=1., ao=0.        !: Scale factor and add_offset
    REAL(KIND=4)                                :: spval  !: missing value
    REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: zend, zstart
    CHARACTER(LEN=256)                          :: clvar
    LOGICAL                                     :: lliom=.false., llperio=.false.
    LOGICAL                                     :: llog=.FALSE. , lsf=.FALSE. , lao=.FALSE.
    !!
    INTEGER(KIND=4)                :: ityp
    !INTEGER(KIND=4), DIMENSION(:)  :: dimids
    !INTEGER(KIND=4)                :: nAtts
    !!---------------------------------------------------------------------
    llperio=.false.
    IF (PRESENT(klev) ) THEN
       ilev=klev
    ELSE
       ilev=1
    ENDIF
    ! Optionall arguments

    IF (PRESENT(kimin) ) THEN
       imin=kimin

       ipiglo=getdim(cdfile, cn_x, ldexact=.true.)
       IF (imin+kpi-1 > ipiglo ) THEN 
         llperio=.true.
         imax=kpi+1 +imin -ipiglo
       ENDIF
    ELSE
       imin=1
    ENDIF

    IF (PRESENT(kjmin) ) THEN
       jmin=kjmin
    ELSE
       jmin=1
    ENDIF

    IF (PRESENT(ktime) ) THEN
       itime=ktime
    ELSE
       itime=1
    ENDIF

    IF (PRESENT(ldiom) ) THEN
       lliom=ldiom
    ELSE
       lliom=.false.
    ENDIF

    ! Must reset the flags to false for every call to getvar
    clvar=cdvar
    llog = .FALSE.
    lsf  = .FALSE.
    lao  = .FALSE.

    CALL ERR_HDL(NF90_OPEN(cdfile,NF90_NOWRITE,incid) )

    IF ( lliom) THEN  !
      IF ( clvar == cn_ve3t ) THEN
        SELECT CASE ( cg_zgr_ver )
        CASE ( 'v2.0' ) ; clvar = 'e3t_ps'
        CASE ( 'v3.0' ) ; clvar = 'e3t'
        CASE ( 'v3.6' ) ; clvar = 'e3t_0'
        END SELECT
      ELSE IF ( clvar == cn_ve3u ) THEN
        SELECT CASE ( cg_zgr_ver )
        CASE ( 'v2.0' ) ; clvar = 'e3u_ps'
        CASE ( 'v3.0' ) ; clvar = 'e3u'
        CASE ( 'v3.6' ) ; clvar = 'e3u_0'
        END SELECT
      ELSE IF ( clvar == cn_ve3v ) THEN
        SELECT CASE ( cg_zgr_ver )
        CASE ( 'v2.0' ) ; clvar = 'e3v_ps'
        CASE ( 'v3.0' ) ; clvar = 'e3v'
        CASE ( 'v3.6' ) ; clvar = 'e3v_0'
        END SELECT
      ELSE IF ( clvar == cn_ve3w ) THEN
        SELECT CASE ( cg_zgr_ver )
        CASE ( 'v2.0' ) ; clvar = 'e3w_ps'
        CASE ( 'v3.0' ) ; clvar = 'e3w'
        CASE ( 'v3.6' ) ; clvar = 'e3w_0'
        END SELECT
      ENDIF
    ENDIF

    istatus=NF90_INQUIRE(incid, unlimitedDimId=id_dimunlim)
    CALL ERR_HDL(NF90_INQ_VARID ( incid,clvar,id_var))

    ! look for time dim in variable
    inldim=0
    istatus=NF90_INQUIRE_VARIABLE(incid, id_var, ndims=inbdim,dimids=inldim(:) )

    istart(1) = imin
    istart(2) = jmin
    ! JMM ! it workd for X Y Z T file,   not for X Y T .... try to found a fix !
    IF ( inldim(3) == id_dimunlim ) THEN
      istart(3) = itime
      istart(4) = 1
    ELSE
      istart(3) = ilev
      istart(4) = itime
    ENDIF

    icount(1)=kpi
    icount(2)=kpj
    icount(3)=1
    icount(4)=1

    spval = getspval ( cdfile, cdvar)  ! try many kind of missing_value (eg _FillValue _Fillvalue Fillvalue ...)

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'savelog10')
    IF (istatus == NF90_NOERR ) THEN
       ! The file is saved a the log of the field
       istatus=NF90_GET_ATT(incid,id_var,'savelog10',ilog)
       IF ( ilog /= 0 ) llog=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'scale_factor')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'scale_factor',sf)
       IF ( sf /= 1. ) lsf=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'add_offset')
    IF (istatus == NF90_NOERR ) THEN
       ! there is an add_offset for this variable
       istatus=NF90_GET_ATT(incid,id_var,'add_offset', ao)
       IF ( ao /= 0.) lao=.TRUE.
    ENDIF


    IF (llperio ) THEN ! Deal with E-W periodic conditions ( used when reading across the folding line)
      ALLOCATE (zend (ipiglo-imin,kpj), zstart(imax-1,kpj) )
       istatus=NF90_GET_VAR(incid,id_var,zend,   start=(/imin,jmin,ilev,itime/),count=(/ipiglo-imin,kpj,1,1/))
       istatus=NF90_GET_VAR(incid,id_var,zstart, start=(/2   ,jmin,ilev,itime/),count=(/imax-1,     kpj,1,1/))
       getvar(1:ipiglo-imin    ,:) = zend
       getvar(ipiglo-imin+1:kpi,:) = zstart
      DEALLOCATE(zstart, zend )
    ELSE
      istatus=NF90_GET_VAR(incid,id_var,getvar, start=istart,count=icount)
    ENDIF

    IF ( istatus /= 0 ) THEN
       PRINT *,' Problem in getvar for ', TRIM(clvar)
       CALL ERR_HDL(istatus)
       STOP 98
    ENDIF

    ! Caution : order does matter !
    IF (lsf )  WHERE (getvar /= spval )  getvar=getvar*sf
    IF (lao )  WHERE (getvar /= spval )  getvar=getvar + ao
    IF (llog)  WHERE (getvar /= spval )  getvar=10**getvar

    istatus=NF90_CLOSE(incid)

  END FUNCTION getvar

  FUNCTION  getvar3d (cdfile,cdvar,kpi,kpj,kpz, kimin, kjmin, kkmin, ktime )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION  getvar3d  ***
    !!
    !! ** Purpose : Return the 3D REAL variable cdvar(x,y,z), from cdfile 
    !!              kpi,kpj,kpz are the horizontal size of the 3D variable
    !!              Sub domain can be optionaly specified with its starting
    !!              point, kimin, kjmin, kkmin
    !!
    !! ** Method  : Use NF90 primitive to read the block of 3D data
    !!
    !!---------------------------------------------------------------------
    CHARACTER(LEN=*),          INTENT(in) :: cdfile
    CHARACTER(LEN=*),          INTENT(in) :: cdvar
    INTEGER(KIND=4),           INTENT(in) :: kpi,   kpj,   kpz
    INTEGER(KIND=4), OPTIONAL, INTENT(in) :: kimin, kjmin, kkmin
    INTEGER(KIND=4), OPTIONAL, INTENT(in) :: ktime         ! Optional variable. If missing 1 is assumed 
    REAL(KIND=4), DIMENSION(kpi,kpj,kpz)  :: getvar3d      ! 3D REAL 

    INTEGER(KIND=4), DIMENSION(4) :: istart, icount
    INTEGER(KIND=4)               :: incid, id_var
    INTEGER(KIND=4)               :: istatus
    INTEGER(KIND=4)               :: iimin, ijmin, ikmin
    INTEGER(KIND=4)               :: itime, ilog
    INTEGER(KIND=4)               :: idum
    REAL(KIND=4)                  :: sf=1., ao=0.       !  Scale factor and add_offset
    REAL(KIND=4)                  :: spval              !  Missing values
    LOGICAL                       :: llog=.FALSE. , lsf=.FALSE. , lao=.FALSE.
    !!---------------------------------------------------------------------
    IF (PRESENT(kimin) ) THEN
       iimin=kimin
    ELSE
       iimin=1
    ENDIF

    IF (PRESENT(kjmin) ) THEN
       ijmin=kjmin
    ELSE
       ijmin=1
    ENDIF

    IF (PRESENT(kkmin) ) THEN
       ikmin=kkmin
    ELSE
       ikmin=1
    ENDIF

    IF (PRESENT(ktime) ) THEN
       itime=ktime
    ELSE
       itime=1
    ENDIF
   
    ! Must reset the flags to false for every call to getvar
    llog=.FALSE.
    lsf=.FALSE.
    lao=.FALSE.
    PRINT *,' GETVAR3D '
    PRINT *,'  ', TRIM(cdfile)
    PRINT *,'  ', TRIM(cdvar )
    PRINT *,'    KPI KPJ KPZ ',kpi, kpj, kpz

    CALL ERR_HDL(NF90_OPEN(cdfile,NF90_NOWRITE,incid) )
    CALL ERR_HDL(NF90_INQ_VARID ( incid, cdvar, id_var) )
    istart=(/iimin, ijmin, ikmin, itime/)
    icount=(/kpi,   kpj,   kpz,   1    /)

    spval = getspval ( cdfile, cdvar )

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'savelog10')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'savelog10',ilog)
       IF ( ilog /= 0 ) llog=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'scale_factor')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'scale_factor',sf)
       IF ( sf /= 1. ) lsf=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'add_offset')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'add_offset',ao)
       IF ( ao /= 0.) lao=.TRUE.
    ENDIF

    istatus=NF90_GET_VAR(incid,id_var,getvar3d, start=istart,count=icount)
    IF ( istatus /= 0 ) THEN
       PRINT *,' Problem in getvar3d for ', TRIM(cdvar)
       CALL ERR_HDL(istatus)
       STOP 98
    ENDIF

    ! Caution : order does matter !
    IF (lsf )  WHERE (getvar3d /= spval )  getvar3d=getvar3d*sf
    IF (lao )  WHERE (getvar3d /= spval )  getvar3d=getvar3d + ao
    IF (llog)  WHERE (getvar3d /= spval )  getvar3d=10**getvar3d

    istatus=NF90_CLOSE(incid)

  END FUNCTION getvar3d

  FUNCTION  getvar3dt (cdfile,cdvar,kk, kpi,kpj,kpt, kimin, kjmin, ktmin )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION  getvar3d  ***
    !!
    !! ** Purpose : Return the 3D REAL variable cdvar(x,y,t), at level kk
    !!              from cdfile.  kpi,kpj,kpt are the size of the variable
    !!              subdomain starting point can be specified by kimin, kjmin, ktmin
    !!
    !! ** Method  : Use NF90 primitive to read the block of 3D data
    !!
    !!---------------------------------------------------------------------
    CHARACTER(LEN=*),          INTENT(in) :: cdfile
    CHARACTER(LEN=*),          INTENT(in) :: cdvar
    INTEGER(KIND=4),           INTENT(in) :: kpi,   kpj,   kk, kpt
    INTEGER(KIND=4), OPTIONAL, INTENT(in) :: kimin, kjmin, ktmin
    REAL(KIND=4), DIMENSION(kpi,kpj,kpt)  :: getvar3dt      ! 3D REAL 

    INTEGER(KIND=4), DIMENSION(4) :: istart, icount
    INTEGER(KIND=4)               :: jt
    INTEGER(KIND=4)               :: incid, id_var, iid
    INTEGER(KIND=4)               :: istatus
    INTEGER(KIND=4)               :: iimin, ijmin, itmin
    INTEGER(KIND=4)               :: itime, ilog
    INTEGER(KIND=4)               :: idum
    REAL(KIND=4)                  :: sf=1., ao=0.       !  Scale factor and add_offset
    REAL(KIND=4)                  :: spval              !  Missing values
    LOGICAL                       :: llog=.FALSE. , lsf=.FALSE. , lao=.FALSE.
    !!---------------------------------------------------------------------
    IF (PRESENT(kimin) ) THEN
       iimin=kimin
    ELSE
       iimin=1
    ENDIF

    IF (PRESENT(kjmin) ) THEN
       ijmin=kjmin
    ELSE
       ijmin=1
    ENDIF

    IF (PRESENT(ktmin) ) THEN
       itmin=ktmin
    ELSE
       itmin=1
    ENDIF

    ! Must reset the flags to false for every call to getvar
    llog=.FALSE.
    lsf=.FALSE.
    lao=.FALSE.
    PRINT *,' GETVAR3DT '
    PRINT *,'  ', TRIM(cdfile)
    PRINT *,'  ', TRIM(cdvar )
    PRINT *,'    KPI KPJ KPT, KK ',kpi, kpj, kpt, kk

    CALL ERR_HDL(NF90_OPEN(cdfile,NF90_NOWRITE,incid) )
    CALL ERR_HDL(NF90_INQ_VARID ( incid, cdvar, id_var) )
    istatus= NF90_INQUIRE_VARIABLE(incid, id_var, ndims=iid)
    IF ( iid == 4 ) THEN
      istart=(/iimin, ijmin, kk, itmin/)
      icount=(/kpi,   kpj,    1, kpt  /)
    ELSEIF  ( iid == 3 ) THEN ! assume X Y T ( cannot be X, Z, T nor Y Z T )
      istart=(/iimin, ijmin, itmin, 1/)
      icount=(/kpi,   kpj,   kpt  , 1/)
    ENDIF

    spval = getspval ( cdfile, cdvar )

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'savelog10')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'savelog10',ilog)
       IF ( ilog /= 0 ) llog=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'scale_factor')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'scale_factor',sf)
       IF ( sf /= 1. ) lsf=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'add_offset')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'add_offset',ao)
       IF ( ao /= 0.) lao=.TRUE.
    ENDIF

    istatus=NF90_GET_VAR(incid,id_var, getvar3dt(:,:,:), start=istart,count=icount)
    IF ( istatus /= 0 ) THEN
       PRINT *,' Problem in getvar3dt for ', TRIM(cdvar)
       CALL ERR_HDL(istatus)
       STOP 98
    ENDIF

    ! Caution : order does matter !
    IF (lsf )  WHERE (getvar3dt /= spval )  getvar3dt=getvar3dt*sf
    IF (lao )  WHERE (getvar3dt /= spval )  getvar3dt=getvar3dt + ao
    IF (llog)  WHERE (getvar3dt /= spval )  getvar3dt=10**getvar3dt

    istatus=NF90_CLOSE(incid)

  END FUNCTION getvar3dt

  FUNCTION  getvar4d (cdfile,cdvar,kpi,kpj,kpz,kpt, kimin, kjmin, kkmin, ktmin )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION  getvar4d  ***
    !!
    !! ** Purpose : Return the 4D REAL variable cdvar(x,y,z,t).
    !!              kpi,kpj,kpz,kpt are the dimensions of the 4D variable.
    !!              A subdomain can be specified using kimin, kjmin, kkmin, ktmin
    !!
    !! ** Method  : Use NF90 primitive to read the block of 3D data
    !!
    !!---------------------------------------------------------------------
    CHARACTER(LEN=*),          INTENT(in) :: cdfile
    CHARACTER(LEN=*),          INTENT(in) :: cdvar
    INTEGER(KIND=4),           INTENT(in) :: kpi,   kpj,   kpz, kpt
    INTEGER(KIND=4), OPTIONAL, INTENT(in) :: kimin, kjmin, kkmin, ktmin
    REAL(KIND=4), DIMENSION(kpi,kpj,kpz,kpt) :: getvar4d      ! 3D REAL 

    INTEGER(KIND=4), DIMENSION(4) :: istart, icount
    INTEGER(KIND=4)               :: incid, id_var, iid
    INTEGER(KIND=4)               :: istatus
    INTEGER(KIND=4)               :: iimin, ijmin, ikmin, itmin
    INTEGER(KIND=4)               :: ilog
    INTEGER(KIND=4)               :: idum
    REAL(KIND=4)                  :: sf=1., ao=0.       !  Scale factor and add_offset
    REAL(KIND=4)                  :: spval              !  Missing values
    LOGICAL                       :: llog=.FALSE. , lsf=.FALSE. , lao=.FALSE.
    !!---------------------------------------------------------------------
    IF (PRESENT(kimin) ) THEN
       iimin=kimin
    ELSE
       iimin=1
    ENDIF

    IF (PRESENT(kjmin) ) THEN
       ijmin=kjmin
    ELSE
       ijmin=1
    ENDIF

    IF (PRESENT(kkmin) ) THEN
       ikmin=kkmin
    ELSE
       ikmin=1
    ENDIF

    IF (PRESENT(ktmin) ) THEN
       itmin=ktmin
    ELSE
       itmin=1
    ENDIF
   
    ! Must reset the flags to false for every call to getvar
    llog=.FALSE.
    lsf=.FALSE.
    lao=.FALSE.
    PRINT *,' GETVAR4D '
    PRINT *,'  ', TRIM(cdfile)
    PRINT *,'  ', TRIM(cdvar )
    PRINT *,'    KPI KPJ KPZ, KPT ',kpi, kpj, kpz, kpt

    CALL ERR_HDL(NF90_OPEN(cdfile,NF90_NOWRITE,incid) )
    CALL ERR_HDL(NF90_INQ_VARID ( incid, cdvar, id_var) )
    istatus= NF90_INQUIRE_VARIABLE(incid, id_var, ndims=iid)
    IF ( iid == 4 ) THEN
      istart=(/iimin, ijmin, ikmin, itmin/)
      icount=(/kpi,   kpj,   kpz,   kpt  /)
    ELSEIF  ( iid == 3 ) THEN ! assume X Y T ( cannot be X, Z, T nor Y Z T )
      istart=(/iimin, ijmin, itmin, 1/)
      icount=(/kpi,   kpj,   kpt,   1/)
    ENDIF

    spval = getspval ( cdfile, cdvar )

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'savelog10')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'savelog10',ilog)
       IF ( ilog /= 0 ) llog=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'scale_factor')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'scale_factor',sf)
       IF ( sf /= 1. ) lsf=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'add_offset')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'add_offset',ao)
       IF ( ao /= 0.) lao=.TRUE.
    ENDIF

    istatus=NF90_GET_VAR(incid,id_var,getvar4d, start=istart,count=icount)
    IF ( istatus /= 0 ) THEN
       PRINT *,' Problem in getvar4d for ', TRIM(cdvar)
       CALL ERR_HDL(istatus)
       STOP 98
    ENDIF

    ! Caution : order does matter !
    IF (lsf )  WHERE (getvar4d /= spval )  getvar4d=getvar4d*sf
    IF (lao )  WHERE (getvar4d /= spval )  getvar4d=getvar4d + ao
    IF (llog)  WHERE (getvar4d /= spval )  getvar4d=10**getvar4d

    istatus=NF90_CLOSE(incid)

  END FUNCTION getvar4d


  FUNCTION  getvarxz (cdfile, cdvar, kj, kpi, kpz, kimin, kkmin, ktime)
    !!-------------------------------------------------------------------------
    !!                  ***  FUNCTION  getvarxz  ***
    !!
    !! ** Purpose : Return the 2D REAL variable x-z slab cvar, from cdfile at j=kj
    !!              kpi,kpz are the  size of the 2D variable
    !!
    !!-------------------------------------------------------------------------
    CHARACTER(LEN=*),          INTENT(in) :: cdfile        ! file name to work with
    CHARACTER(LEN=*),          INTENT(in) :: cdvar         ! variable name to work with
    INTEGER(KIND=4),           INTENT(in) :: kj            ! Optional variable. If missing 1 is assumed
    INTEGER(KIND=4),           INTENT(in) :: kpi, kpz      ! size of the 2D variable
    INTEGER(KIND=4), OPTIONAL, INTENT(in) :: kimin, kkmin  ! Optional variable. If missing 1 is assumed
    INTEGER(KIND=4), OPTIONAL, INTENT(in) :: ktime         ! Optional variable. If missing 1 is assumed 
    REAL(KIND=4), DIMENSION(kpi,kpz)      :: getvarxz      ! 2D REAL 4 holding variable x-z slab at kj

    INTEGER(KIND=4), DIMENSION(4) :: istart, icount
    INTEGER(KIND=4)               :: incid, id_var
    INTEGER(KIND=4)               :: istatus, ilev, imin, kmin
    INTEGER(KIND=4)               :: itime, ilog
    INTEGER(KIND=4)               :: idum
    REAL(KIND=4)                  :: sf=1., ao=0.       !  Scale factor and add_offset
    REAL(KIND=4)                  :: spval              !  Missing values
    LOGICAL                       :: llog=.FALSE. , lsf=.FALSE. , lao=.FALSE.
    !!-------------------------------------------------------------------------

    IF (PRESENT(kimin) ) THEN
       imin=kimin
    ELSE
       imin=1
    ENDIF

    IF (PRESENT(kkmin) ) THEN
       kmin=kkmin
    ELSE
       kmin=1
    ENDIF

    IF (PRESENT(ktime) ) THEN
       itime=ktime
    ELSE
       itime=1
    ENDIF

    ! Must reset the flags to false for every call to getvar
    llog=.FALSE.
    lsf=.FALSE.
    lao=.FALSE.


    CALL ERR_HDL(NF90_OPEN(cdfile,NF90_NOWRITE,incid) )
    CALL ERR_HDL(NF90_INQ_VARID ( incid,cdvar,id_var))

    spval = getspval ( cdfile, cdvar )

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'savelog10')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'savelog10',ilog)
       IF ( ilog /= 0 ) llog=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'scale_factor')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'scale_factor',sf)
       IF ( sf /= 1. ) lsf=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'add_offset')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'add_offset',ao)
       IF ( ao /= 0.) lao=.TRUE.
    ENDIF

    ! detect if there is a y dimension in cdfile
    istatus=NF90_INQ_DIMID(incid,'y',idum)
    IF ( istatus == NF90_NOERR ) THEN  ! the file has a 'y' dimension
      istart=(/imin,kj,kmin,itime/)
      ! JMM ! it workd for X Y Z T file,   not for X Y T .... try to found a fix !
      icount=(/kpi,1,kpz,1/)
    ELSE    ! no y dimension
      istart=(/imin,kmin,itime,1/)
      icount=(/kpi,kpz,1,1/)
    ENDIF

    istatus=NF90_GET_VAR(incid,id_var,getvarxz, start=istart,count=icount)
    IF ( istatus /= 0 ) THEN
       PRINT *,' Problem in getvarxz for ', TRIM(cdvar)
       CALL ERR_HDL(istatus)
       STOP 98
    ENDIF

    ! Caution : order does matter !
    IF (lsf )  WHERE (getvarxz /= spval )  getvarxz=getvarxz*sf
    IF (lao )  WHERE (getvarxz /= spval )  getvarxz=getvarxz + ao
    IF (llog)  WHERE (getvarxz /= spval )  getvarxz=10**getvarxz

    istatus=NF90_CLOSE(incid)

  END FUNCTION getvarxz


  FUNCTION  getvaryz (cdfile, cdvar, ki, kpj, kpz, kjmin, kkmin, ktime)
    !!-------------------------------------------------------------------------
    !!                  ***  FUNCTION  getvar  ***
    !!
    !! ** Purpose : Return the 2D REAL variable y-z slab cvar, from cdfile at i=ki
    !!              kpj,kpz are the  size of the 2D variable
    !!
    !!-------------------------------------------------------------------------
    CHARACTER(LEN=*),          INTENT(in) :: cdfile       ! file name to work with
    CHARACTER(LEN=*),          INTENT(in) :: cdvar        ! variable name to work with
    INTEGER(KIND=4),           INTENT(in) :: ki           ! 
    INTEGER(KIND=4),           INTENT(in) :: kpj,kpz      ! size of the 2D variable
    INTEGER(KIND=4), OPTIONAL, INTENT(in) :: kjmin, kkmin ! Optional variable. If missing 1 is assumed
    INTEGER(KIND=4), OPTIONAL, INTENT(in) :: ktime        ! Optional variable. If missing 1 is assumed
    REAL(KIND=4), DIMENSION(kpj,kpz)      :: getvaryz     ! 2D REAL 4 holding variable x-z slab at kj

    INTEGER(KIND=4), DIMENSION(4)       :: istart, icount
    INTEGER(KIND=4)                     :: incid, id_var
    INTEGER(KIND=4)                     :: istatus, ilev, jmin, kmin
    INTEGER(KIND=4)                     :: itime, ilog
    INTEGER(KIND=4)                     :: idum

    REAL(KIND=4)                        :: sf=1., ao=0.   !  Scale factor and add_offset
    REAL(KIND=4)                        :: spval          !  Missing values
    LOGICAL                             :: llog=.FALSE. , lsf=.FALSE. , lao=.FALSE.
    !!-------------------------------------------------------------------------

    IF (PRESENT(kjmin) ) THEN
       jmin=kjmin
    ELSE
       jmin=1
    ENDIF

    IF (PRESENT(kkmin) ) THEN
       kmin=kkmin
    ELSE
       kmin=1
    ENDIF

    IF (PRESENT(ktime) ) THEN
       itime=ktime
    ELSE
       itime=1
    ENDIF

    ! Must reset the flags to false for every call to getvar
    llog=.FALSE.
    lsf=.FALSE.
    lao=.FALSE.


    CALL ERR_HDL(NF90_OPEN(cdfile,NF90_NOWRITE,incid) )
    CALL ERR_HDL(NF90_INQ_VARID ( incid,cdvar,id_var))

    spval = getspval ( cdfile, cdvar )

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'savelog10')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'savelog10',ilog)
       IF ( ilog /= 0 ) llog=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'scale_factor')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'scale_factor',sf)
       IF ( sf /= 1. ) lsf=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(incid,id_var,'add_offset')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(incid,id_var,'add_offset', ao)
       IF ( ao /= 0.) lao=.TRUE.
    ENDIF

    ! detect if there is a x dimension in cdfile
    istatus=NF90_INQ_DIMID(incid,'x',idum)
    IF ( istatus == NF90_NOERR ) THEN  ! the file has a 'x' dimension
      istart=(/ki,jmin,kmin,itime/)
      ! JMM ! it workd for X Y Z T file,   not for X Y T .... try to found a fix !
      icount=(/1,kpj,kpz,1/)
    ELSE    ! no x dimension
      istart=(/jmin,kmin,itime,1/)
      icount=(/kpj,kpz,1,1/)
    ENDIF

    istatus=NF90_GET_VAR(incid,id_var,getvaryz, start=istart,count=icount)
    IF ( istatus /= 0 ) THEN
       PRINT *,' Problem in getvaryz for ', TRIM(cdvar)
       CALL ERR_HDL(istatus)
       STOP 98
    ENDIF

    ! Caution : order does matter !
    IF (lsf )  WHERE (getvaryz /= spval )  getvaryz=getvaryz*sf
    IF (lao )  WHERE (getvaryz /= spval )  getvaryz=getvaryz + ao
    IF (llog)  WHERE (getvaryz /= spval )  getvaryz=10**getvaryz

    istatus=NF90_CLOSE(incid)

  END FUNCTION getvaryz


  FUNCTION  getvar1d (cdfile, cdvar, kk, kstatus)
    !!-------------------------------------------------------------------------
    !!                  ***  FUNCTION  getvar1d  ***
    !!
    !! ** Purpose :  return 1D variable cdvar from cdfile, of size kk
    !!
    !!-------------------------------------------------------------------------
    CHARACTER(LEN=*),           INTENT(in) :: cdfile   ! file name to work with
    CHARACTER(LEN=*),           INTENT(in) :: cdvar    ! variable name to work with
    INTEGER(KIND=4),            INTENT(in) :: kk       ! size of 1D vector to be returned
    INTEGER(KIND=4), OPTIONAL, INTENT(out) :: kstatus  ! return status concerning the variable existence
    REAL(KIND=8), DIMENSION(kk)            :: getvar1d ! returned vector (double precision --for time --)
                                                       ! 

    INTEGER(KIND=4), DIMENSION(1) :: istart, icount
    INTEGER(KIND=4) :: incid, id_var
    INTEGER(KIND=4) :: istatus
    !!-------------------------------------------------------------------------
    istart(:) = 1
    icount(1)=kk
    IF ( PRESENT(kstatus) ) kstatus = 0

    istatus=NF90_OPEN(cdfile,NF90_NOWRITE,incid)
    istatus=NF90_INQ_VARID ( incid,cdvar,id_var)
    IF ( istatus == NF90_NOERR ) THEN
       istatus=NF90_GET_VAR(incid,id_var,getvar1d,start=istart,count=icount)
    ELSE
       IF ( PRESENT(kstatus) ) kstatus= istatus
       getvar1d=99999999999.d0
    ENDIF

    istatus=NF90_CLOSE(incid)

  END FUNCTION getvar1d


  FUNCTION  getvare3 (cdfile,cdvar,kk)
    !!-------------------------------------------------------------------------
    !!                  ***  FUNCTION  getvare3  ***
    !!
    !! ** Purpose :  Special routine for e3, which in fact is a 1D variable
    !!               but defined as e3 (1,1,npk,1) in coordinates.nc (!!)
    !!
    !!-------------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdfile   ! file name to work with
    CHARACTER(LEN=*), INTENT(in) :: cdvar    ! variable name to work with
    INTEGER(KIND=4),  INTENT(in) :: kk       ! size of 1D vector to be returned
    REAL(KIND=4),  DIMENSION(kk) :: getvare3 ! return e3 variable form the coordinate file

    INTEGER(KIND=4), DIMENSION(4) :: istart, icount
    INTEGER(KIND=4)               :: incid, id_var
    INTEGER(KIND=4)               :: istatus
    CHARACTER(LEN=256)            :: clvar   ! local name for cdf var (modified)
    !!-------------------------------------------------------------------------
    istart(:) = 1
    icount(:) = 1
    icount(3) = kk
    clvar=cdvar

    istatus=NF90_OPEN(cdfile,NF90_NOWRITE,incid)
    ! check for IOM style mesh_zgr or coordinates :
    ! IOIPSL (x_a=y_a=1)(2.0)          IOM(3.0)           3.6
    ! gdept(time,z,y_a,x_a)            gdept_0(t,z)    gdept_1d(t,z)
    ! gdepw(time,z,y_a,x_a)            gdepw_0(t,z)    gdepw_1d(t,z)
    !   e3t(time,z,y_a,x_a)            e3t_0(t,z)      e3t_1d(t,z)
    !   e3w(time,z,y_a,x_a)            e3w_0(t,z)      e3w_1d(t,z)

    ! change icount for 1D variables in mesh_zgr only
    IF (    clvar == cn_gdept  .OR. &
          & clvar == cn_gdepw  .OR. &
          & clvar == cn_ve3t1d .OR. &
          & clvar == cn_ve3w1d      ) THEN
       SELECT CASE ( cg_zgr_ver ) 
       CASE ( 'v2.0') ; icount(1)=1  ; icount(3)=kk
       CASE ( 'v3.0') ; icount(1)=kk ; icount(3)=1
       CASE ( 'v3.6') ; icount(1)=kk ; icount(3)=1
       END SELECT
    ENDIF
   
    IF ( clvar == cn_gdept) THEN
    SELECT CASE ( cg_zgr_ver )
    CASE ( 'v2.0')
      clvar = 'gdept'
    CASE ( 'v3.0')
      clvar = 'gdept_0'
    CASE ( 'v3.6')
      clvar = 'gdept_1d'
    END SELECT
    ENDIF

    IF ( clvar == cn_gdepw) THEN
    SELECT CASE ( cg_zgr_ver )
    CASE ( 'v2.0')
      clvar = 'gdepw'
    CASE ( 'v3.0')
      clvar = 'gdepw_0'
    CASE ( 'v3.6')
      clvar = 'gdepw_1d'
    END SELECT
    ENDIF

    IF ( clvar == cn_ve3t1d) THEN
    SELECT CASE ( cg_zgr_ver )
    CASE ( 'v2.0')
      clvar = 'e3t'
    CASE ( 'v3.0')
      clvar = 'e3t_0'
    CASE ( 'v3.6')
      clvar = 'e3t_1d'
    END SELECT
    ENDIF

    IF ( clvar == cn_ve3w1d) THEN
    SELECT CASE ( cg_zgr_ver )
    CASE ( 'v2.0')
      clvar = 'e3w'
    CASE ( 'v3.0')
      clvar = 'e3w_0'
    CASE ( 'v3.6')
      clvar = 'e3w_1d'
    END SELECT
    ENDIF
    
    istatus=NF90_INQ_VARID ( incid,clvar,id_var)
    istatus=NF90_GET_VAR(incid,id_var,getvare3,start=istart,count=icount)
    IF ( istatus /= 0 ) THEN
       PRINT *,' Problem in getvare3 for ', TRIM(cdvar)
       PRINT *,TRIM(cdfile), kk
       CALL ERR_HDL(istatus)
       STOP 98
    ENDIF

    istatus=NF90_CLOSE(incid)

  END FUNCTION getvare3


  INTEGER(KIND=4) FUNCTION putheadervar(kout, cdfile, kpi, kpj, kpk, pnavlon, pnavlat , pdep, cdep, ld_xycoo)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION  putheadervar  ***
    !!
    !! ** Purpose :  copy header variables from cdfile to the already open ncfile (ncid=kout)
    !!
    !! ** Method  :  header variables are nav_lat, nav_lon and either (deptht, depthu, or depthv )
    !!               Even if the use of different variable name for deptht, depthu depthv is
    !!               one of the many non sense of IOIPSL, we are forced to stick with !
    !!               (Note that these 3 depth are identical in OPA. On the other hand, nav_lon, nav_lat
    !!               differ for U and V and T points but have the same variable name).
    !!               If pnavlon and pnavlat are provided as arguments, they are used for nav_lon, nav_lat
    !!               instead of the nav_lon,nav_lat read on the file cdfile.
    !!
    !! ** Action  : header variables for file kout is copied from cdfile
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),                            INTENT(in) :: kout     ! ncid of the outputfile (already open )
    CHARACTER(LEN=*),                           INTENT(in) :: cdfile   ! file from where the headers will be copied
    INTEGER(KIND=4),                            INTENT(in) :: kpi, kpj ! dimension of nav_lon (kpi,kpj)
    INTEGER(KIND=4),                            INTENT(in) :: kpk      ! dimension of depht(kpk)
    LOGICAL,      OPTIONAL,                     INTENT(in) :: ld_xycoo ! option to put yx info
    REAL(KIND=4), OPTIONAL, DIMENSION(kpi,kpj), INTENT(in) :: pnavlon  ! array provided optionaly to overrid the
    REAL(KIND=4), OPTIONAL, DIMENSION(kpi,kpj), INTENT(in) :: pnavlat  ! corresponding arrays in cdfile 
    REAL(KIND=4), OPTIONAL, DIMENSION(kpk),     INTENT(in) :: pdep     ! dep array if not on cdfile
    CHARACTER(LEN=*), OPTIONAL,                 INTENT(in) :: cdep     ! optional name of vertical variable

    INTEGER(KIND=4), PARAMETER                :: jpdep=6
    INTEGER(KIND=4)                           :: istatus, idep, jj
    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: z2d
    REAL(KIND=4), DIMENSION(kpk)              :: z1d
    CHARACTER(LEN=256), DIMENSION(jpdep )     :: cldept= (/'deptht ','depthu ','depthv ','depthw ','nav_lev','z      '/)
    CHARACTER(LEN=256)                        :: cldep
    LOGICAL                                   :: ll_xycoo
    !!----------------------------------------------------------------------
    IF (PRESENT(ld_xycoo) ) THEN 
      ll_xycoo = ld_xycoo
    ELSE
      ll_xycoo = .true.
    ENDIF

    cldept = (/cn_vdeptht, cn_vdepthu, cn_vdepthv, cn_vdepthw,'nav_lev','z      '/)

    IF ( ll_xycoo  ) THEN   
       ALLOCATE ( z2d (kpi,kpj) )
   
       IF (PRESENT(pnavlon) ) THEN 
          z2d = pnavlon
       ELSE
          IF ( chkvar ( cdfile, cn_vlon2d )) THEN
              IF (chkvar (cdfile, cn_vlon1d) ) THEN
                PRINT *, '... dummy value used!'
                z2d = 0.
              ELSE
                z2d(:,1) = getvar1d(cdfile, cn_vlon1d, kpi, istatus ) 
                DO jj=2,kpj
                  z2d(:,jj) = z2d(:,1)
                ENDDO
              ENDIF
          ELSE
            z2d=getvar(cdfile,cn_vlon2d, 1,kpi,kpj)
          ENDIF
       ENDIF
       istatus = putvar(kout, nid_lon,z2d,1,kpi,kpj)
   
       IF (PRESENT(pnavlat) ) THEN
          z2d = pnavlat
       ELSE
          IF ( chkvar ( cdfile, cn_vlat2d )) THEN
              IF (chkvar (cdfile, cn_vlat1d) ) THEN
                PRINT *, '... dummy value used!'
                z2d = 0.
              ELSE
                z2d(1,:) = getvar1d(cdfile, cn_vlat1d, kpj, istatus ) 
                DO jj=2,kpi
                  z2d(jj,:) = z2d(1,:)
                ENDDO
              ENDIF
          ELSE
            z2d=getvar(cdfile,cn_vlat2d, 1,kpi,kpj)
          ENDIF
       ENDIF
   
       istatus = putvar(kout, nid_lat,z2d,1,kpi,kpj)
       DEALLOCATE (z2d)
    ENDIF

    IF (kpk /= 0 ) THEN
       IF (PRESENT(pdep) ) THEN
          z1d = pdep
       ELSE
          idep = NF90_NOERR

          IF ( PRESENT (cdep)) THEN
             z1d=getvar1d(cdfile,cdep,kpk,idep)
          ENDIF

          IF ( .NOT. PRESENT(cdep) .OR. idep /= NF90_NOERR ) THEN  ! look for standard dep name
             DO jj = 1,jpdep
                cldep=cldept(jj)
                z1d=getvar1d(cdfile,cldep,kpk,idep)
                IF ( idep == NF90_NOERR )  EXIT
             END DO
             IF (jj == jpdep +1 ) THEN
                PRINT *,' No depth variable found in ', TRIM(cdfile)
                STOP 98
             ENDIF
          ENDIF
       ENDIF

       istatus = putvar1d(kout,z1d,kpk,'D')
    ENDIF

    putheadervar=istatus

  END FUNCTION putheadervar


  INTEGER(KIND=4) FUNCTION putvarr8(kout, kid, ptab, klev, kpi, kpj, ktime, kwght)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION putvarr8  ***
    !!            
    !! ** Purpose : copy a 2D level of ptab in already open file kout, 
    !!              using variable kid
    !!
    !! ** Method  : this corresponds to the generic function putvar with r8 arg.
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),                  INTENT(in) :: kout     ! ncid of output file
    INTEGER(KIND=4),                  INTENT(in) :: kid      ! varid of output variable
    REAL(KIND=8), DIMENSION(kpi,kpj), INTENT(in) :: ptab     ! 2D array to write in file
    INTEGER(KIND=4),                  INTENT(in) :: klev     ! level at which ptab will be written
    INTEGER(KIND=4),                  INTENT(in) :: kpi, kpj ! dimension of ptab
    INTEGER(KIND=4), OPTIONAL,        INTENT(in) :: ktime    ! dimension of ptab
    INTEGER(KIND=4), OPTIONAL,        INTENT(in) :: kwght    ! weight of this variable

    INTEGER(KIND=4)               :: istatus, itime, id_dimunlim, inbdim
    INTEGER(KIND=4), DIMENSION(4) :: istart, icount, inldim
    !!----------------------------------------------------------------------
    IF (PRESENT(ktime) ) THEN
       itime=ktime
    ELSE
       itime=1
    ENDIF

     ! look for unlimited dim (time_counter)
    istatus=NF90_INQUIRE         (kout, unlimitedDimId=id_dimunlim       )
    istatus=NF90_INQUIRE_VARIABLE(kout,kid,ndims=inbdim,dimids=inldim(:) )

    !  if the last dim of id_var is time, then adjust the starting point
    istart(:) = 1    ; icount(:) = 1    ! default
    icount(1) = kpi  ; icount(2) = kpj  ! in any case
    IF ( inldim(inbdim) == id_dimunlim ) istart(inbdim) = itime ! assume than last dim is UNLIM
    IF ( inbdim == 4                   ) istart(3     ) = klev

    istatus=NF90_PUT_VAR(kout,kid, ptab, start=istart,count=icount)

    IF (PRESENT(kwght) ) THEN
      istatus=NF90_PUT_ATT(kout, kid, 'iweight', kwght)
    ENDIF
    putvarr8=istatus

  END FUNCTION putvarr8

  INTEGER(KIND=4) FUNCTION putvarr4(kout, kid, ptab, klev, kpi, kpj, ktime, kwght)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION putvarr4  ***
    !!            
    !! ** Purpose : copy a 2D level of ptab in already open file kout, 
    !!              using variable kid
    !!
    !! ** Method  : this corresponds to the generic function putvar with r4 arg.
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),                  INTENT(in) :: kout     ! ncid of output file
    INTEGER(KIND=4),                  INTENT(in) :: kid      ! varid of output variable
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptab     ! 2D array to write in file
    INTEGER(KIND=4),                  INTENT(in) :: klev     ! level at which ptab will be written
    INTEGER(KIND=4),                  INTENT(in) :: kpi, kpj ! dimension of ptab
    INTEGER(KIND=4), OPTIONAL,        INTENT(in) :: ktime    ! dimension of ptab
    INTEGER(KIND=4), OPTIONAL,        INTENT(in) :: kwght    ! weight of this variable

    INTEGER(KIND=4)               :: istatus, itime, id_dimunlim, inbdim
    INTEGER(KIND=4), DIMENSION(4) :: istart, icount, inldim
    !!----------------------------------------------------------------------
    IF (PRESENT(ktime) ) THEN
       itime=ktime
    ELSE
       itime=1
    ENDIF

     ! look for unlimited dim (time_counter)
    istatus=NF90_INQUIRE         (kout, unlimitedDimId=id_dimunlim       )
    istatus=NF90_INQUIRE_VARIABLE(kout,kid,ndims=inbdim,dimids=inldim(:) )

    !  if the last dim of id_var is time, then adjust the starting point
    istart(:) = 1    ; icount(:) = 1    ! default
    icount(1) = kpi  ; icount(2) = kpj  ! in any case
    IF ( inldim(inbdim) == id_dimunlim ) istart(inbdim) = itime ! assume than last dim is UNLIM
    IF ( inbdim == 4                   ) istart(3     ) = klev

    istatus=NF90_PUT_VAR(kout,kid, ptab, start=istart,count=icount)

    IF (PRESENT(kwght) ) THEN
      istatus=NF90_PUT_ATT(kout, kid, 'iweight', kwght)
    ENDIF
    putvarr4=istatus

  END FUNCTION putvarr4


  INTEGER(KIND=4) FUNCTION putvare3(kout, kid, ptab, kmax, cd_vert, ktime, kwght)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION putvarr4  ***
    !!            
    !! ** Purpose : copy a 2D level of ptab in already open file kout, 
    !!              using variable kid
    !!
    !! ** Method  : this corresponds to the generic function putvar with r4 arg.
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),                  INTENT(in) :: kout     ! ncid of output file
    INTEGER(KIND=4),                  INTENT(in) :: kid      ! varid of output variable
    REAL(KIND=4), DIMENSION(kmax),    INTENT(in) :: ptab     ! 1D array to write in file 
    INTEGER(KIND=4),                  INTENT(in) :: kmax     ! number of level in file
    CHARACTER(LEN=*),                 INTENT(in) :: cd_vert  ! dummy var to indicate vertical profile
    INTEGER(KIND=4), OPTIONAL,        INTENT(in) :: ktime    ! 
    INTEGER(KIND=4), OPTIONAL,        INTENT(in) :: kwght    ! weight of this variable

    INTEGER(KIND=4)               :: istatus, itime, id_dimunlim, inbdim
    INTEGER(KIND=4), DIMENSION(4) :: istart, icount, inldim
    !!----------------------------------------------------------------------
    IF (PRESENT(ktime) ) THEN
       itime=ktime
    ELSE
       itime=1
    ENDIF

     ! look for unlimited dim (time_counter)
    istatus=NF90_INQUIRE         (kout, unlimitedDimId=id_dimunlim       )
    istatus=NF90_INQUIRE_VARIABLE(kout,kid,ndims=inbdim,dimids=inldim(:) )

    !  if the last dim of id_var is time, then adjust the starting point
    istart(:) = 1    ; icount(:) = 1    ! default
    icount(1) = 1  ; icount(2) = 1  ; icount(3) = kmax ! in any case
    IF ( inldim(inbdim) == id_dimunlim ) istart(inbdim) = itime ! assume than last dim is UNLIM

    istatus=NF90_PUT_VAR(kout,kid, ptab, start=istart,count=icount)

    IF (PRESENT(kwght) ) THEN
      istatus=NF90_PUT_ATT(kout, kid, 'iweight', kwght)
    ENDIF
    putvare3=istatus

  END FUNCTION putvare3

  INTEGER(KIND=4) FUNCTION putvari2(kout, kid, ktab, klev, kpi, kpj, ktime, kwght)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION putvari2  ***
    !!            
    !! ** Purpose : copy a 2D level of ptab in already open file kout, 
    !!              using variable kid
    !!
    !! ** Method  : this corresponds to the generic function putvar with i2 arg.
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),                     INTENT(in) :: kout     ! ncid of output file
    INTEGER(KIND=4),                     INTENT(in) :: kid      ! varid of output variable
    INTEGER(KIND=2), DIMENSION(kpi,kpj), INTENT(in) :: ktab     ! 2D array to write in file
    INTEGER(KIND=4),                     INTENT(in) :: klev     ! level at which ktab will be written
    INTEGER(KIND=4),                     INTENT(in) :: kpi, kpj ! dimension of ktab
    INTEGER(KIND=4), OPTIONAL,           INTENT(in) :: ktime    ! dimension of ktab
    INTEGER(KIND=4), OPTIONAL,           INTENT(in) :: kwght    ! weight of this variable

    INTEGER(KIND=4)               :: istatus, itime, id_dimunlim, inbdim
    INTEGER(KIND=4), DIMENSION(4) :: istart, icount, inldim
    !!----------------------------------------------------------------------
    IF (PRESENT(ktime) ) THEN
       itime=ktime
    ELSE
       itime=1
    ENDIF

     ! look for unlimited dim (time_counter)
    istatus=NF90_INQUIRE         (kout, unlimitedDimId=id_dimunlim       )
    istatus=NF90_INQUIRE_VARIABLE(kout,kid,ndims=inbdim,dimids=inldim(:) )

    !  if the last dim of id_var is time, then adjust the starting point
    istart(:) = 1    ; icount(:) = 1    ! default
    icount(1) = kpi  ; icount(2) = kpj  ! in any case
    IF ( inldim(inbdim) == id_dimunlim ) istart(inbdim) = itime ! assume than last dim is UNLIM
    IF ( inbdim == 4                   ) istart(3     ) = klev

    istatus=NF90_PUT_VAR(kout,kid, ktab, start=istart,count=icount)

    IF (PRESENT(kwght) ) THEN
      istatus=NF90_PUT_ATT(kout, kid, 'iweight', kwght)
    ENDIF
    putvari2=istatus

  END FUNCTION putvari2


  INTEGER(KIND=4) FUNCTION reputvarr4 (cdfile, cdvar, klev, kpi, kpj, kimin, kjmin, ktime, ptab, kwght)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION reputvarr4  ***
    !!
    !! ** Purpose :  Change an existing variable in inputfile 
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),                 INTENT(in) :: cdfile       ! file name to work with
    CHARACTER(LEN=*),                 INTENT(in) :: cdvar        ! variable name to work with
    INTEGER(KIND=4), OPTIONAL,        INTENT(in) :: klev         ! Optional variable. If missing 1 is assumed
    INTEGER(KIND=4),                  INTENT(in) :: kpi, kpj     ! horizontal size of the 2D variable
    INTEGER(KIND=4), OPTIONAL,        INTENT(in) :: kimin, kjmin ! Optional variable. If missing 1 is assumed
    INTEGER(KIND=4), OPTIONAL,        INTENT(in) :: ktime        ! Optional variable. If missing 1 is assumed
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptab        ! 2D REAL 4 holding variable field at klev
    INTEGER(KIND=4), OPTIONAL,        INTENT(in) :: kwght        ! weight of this variable

    INTEGER(KIND=4) :: incid, id_var, id_dimunlim, inbdim
    INTEGER(KIND=4) :: istatus, ilev, iimin, ijmin, itime
    INTEGER(KIND=4), DIMENSION(4) :: istart, icount, inldim
    !!----------------------------------------------------------------------
    ilev  = 1 ; IF (PRESENT(klev ) ) ilev  = klev
    iimin = 1 ; IF (PRESENT(kimin) ) iimin = kimin
    ijmin = 1 ; IF (PRESENT(kjmin) ) ijmin = kjmin
    itime = 1 ; IF (PRESENT(ktime) ) itime = ktime

    istatus=NF90_OPEN(cdfile,NF90_WRITE,incid)
    istatus=NF90_INQ_VARID(incid,cdvar,id_var)
     ! look for unlimited dim (time_counter)
    istatus=NF90_INQUIRE         (incid,    unlimitedDimId=id_dimunlim       )
    istatus=NF90_INQUIRE_VARIABLE(incid,id_var,ndims=inbdim,dimids=inldim(:) )

    !  if the last dim of id_var is time, then adjust the starting point
    istart(:) = 1     ; icount(:) = 1    ! default
    istart(1) = iimin ; istart(2) = ijmin
    icount(1) = kpi   ; icount(2) = kpj  ! in any case
    IF ( inldim(inbdim) == id_dimunlim ) istart(inbdim) = itime ! assume than last dim is UNLIM
    IF ( inbdim == 4                   ) istart(3     ) = ilev

    istatus=NF90_PUT_VAR(incid,id_var, ptab,start=istart, count=icount )

    IF (PRESENT(kwght)) THEN
      istatus=NF90_PUT_ATT(incid,id_var,'iweight',kwght)
    ENDIF

    reputvarr4=istatus

    istatus=NF90_CLOSE(incid)

  END FUNCTION reputvarr4


  INTEGER(KIND=4) FUNCTION putvarzo(kout, kid, ptab, klev, kpi, kpj, ktime)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION putvarzo  ***
    !!
    !! ** Purpose : Copy a 2D level of ptab in already open file kout, using variable kid
    !!              This variant deals with degenerated 2D (1 x jpj) zonal files
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),              INTENT(in) :: kout             ! ncid of output file
    INTEGER(KIND=4),              INTENT(in) :: kid              ! varid of output variable
    REAL(KIND=4), DIMENSION(:),   INTENT(in) :: ptab             ! 1D array to write in file (x-z or y-z )
    INTEGER(KIND=4),              INTENT(in) :: klev             ! level at which ptab will be written
    INTEGER(KIND=4),              INTENT(in) :: kpi, kpj         ! dimension of ptab
    INTEGER(KIND=4), OPTIONAL,    INTENT(in) :: ktime            ! time to write

    INTEGER(KIND=4)               :: istatus, itime, ilev, id_dimunlim, inbdim
    INTEGER(KIND=4), DIMENSION(4) :: istart, icount, inldim
    !!----------------------------------------------------------------------
    ilev=klev
    IF (PRESENT(ktime) ) THEN
       itime=ktime
    ELSE
       itime=1
    ENDIF

     ! look for unlimited dim (time_counter)
    istatus=NF90_INQUIRE         (kout, unlimitedDimId=id_dimunlim       )
    istatus=NF90_INQUIRE_VARIABLE(kout,kid,ndims=inbdim,dimids=inldim(:) )

    !  if the third dim of id_var is time, then adjust the starting point 
    !  to take ktime into account (case XYT file)
    istart(:) = 1    ; icount(:) = 1    ! default 
    icount(1) = kpi  ; icount(2) = kpj  ! in any case
    IF ( inldim(inbdim) == id_dimunlim ) istart(inbdim) = itime ! assume than last dim is UNLIM
    IF ( inbdim == 4                   ) istart(3) = klev

    istatus=NF90_PUT_VAR(kout,kid, ptab, start=istart,count=icount)
    putvarzo=istatus

  END FUNCTION putvarzo


  INTEGER(KIND=4) FUNCTION putvar1d4(kout, ptab, kk, cdtype)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION putvar1d4  ***
    !!
    !! ** Purpose :  Copy 1D variable (size kk) hold in ptab,  with id 
    !!               kid, into file id kout 
    !!
    !! ** Method  : cdtype is either T (time_counter) or D (depth.)
    !!                         LON (1D longitude) or LAT (1D latitude)
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),            INTENT(in) :: kout   ! ncid of output file
    REAL(KIND=4), DIMENSION(kk),INTENT(in) :: ptab   ! 1D array to write in file
    INTEGER(KIND=4),            INTENT(in) :: kk     ! number of elements in ptab
    CHARACTER(LEN=1),           INTENT(in) :: cdtype ! either T or D LON or LAT

    INTEGER(KIND=4)               :: istatus, iid
    INTEGER(KIND=4), DIMENSION(1) :: istart, icount
    !!----------------------------------------------------------------------
    SELECT CASE ( cdtype )
    CASE ('T', 't' ) 
       iid = nid_tim
    CASE ('D', 'd' )
       iid = nid_dep
    CASE ('X', 'x' )
       iid = nid_lon1d
    CASE ('Y', 'y' )
       iid = nid_lat1d
    END SELECT

    istart(:) = 1
    icount(:) = kk
    istatus=NF90_PUT_VAR(kout,iid, ptab, start=istart,count=icount)
    putvar1d4=istatus

  END FUNCTION putvar1d4

  INTEGER(KIND=4) FUNCTION putvar1d8(kout, ddtab, kk, cdtype)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION putvar1d8  ***
    !!
    !! ** Purpose :  Copy 1D variable (size kk) hold in ptab,  with id 
    !!               kid, into file id kout (double precision)
    !!
    !! ** Method  : cdtype is either T (time_counter) or D (depth.)
    !!                         LON (1D longitude) or LAT (1D latitude)
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),            INTENT(in) :: kout   ! ncid of output file
    REAL(KIND=8), DIMENSION(kk),INTENT(in) :: ddtab  ! 1D array to write in file
    INTEGER(KIND=4),            INTENT(in) :: kk     ! number of elements in ptab
    CHARACTER(LEN=1),           INTENT(in) :: cdtype ! either T or D LON or LAT

    INTEGER(KIND=4)               :: istatus, iid
    INTEGER(KIND=4), DIMENSION(1) :: istart, icount
    !!----------------------------------------------------------------------
    SELECT CASE ( cdtype )
    CASE ('T', 't' ) 
       iid = nid_tim
    CASE ('D', 'd' )
       iid = nid_dep
    CASE ('X', 'x' )
       iid = nid_lon1d
    CASE ('Y', 'y' )
       iid = nid_lat1d
    END SELECT

    istart(:) = 1
    icount(:) = kk
    istatus=NF90_PUT_VAR(kout,iid, ddtab, start=istart,count=icount)
    putvar1d8=istatus

  END FUNCTION putvar1d8

  INTEGER(KIND=4) FUNCTION reputvar1d4(cdfile, cdvar, ptab, kk )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION reputvar1d4  ***
    !!
    !! ** Purpose : Copy 1d variable cdfvar in cdfile, an already existing file
    !!              ptab is the 1d array to write and kk the size of ptab
    !!
    !! ** Method  :   
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),            INTENT(in) :: cdfile      ! filename
    CHARACTER(LEN=*),            INTENT(in) :: cdvar       ! variable name
    REAL(KIND=4), DIMENSION(kk), INTENT(in) :: ptab        ! 1D array to write in file
    INTEGER(KIND=4),             INTENT(in) :: kk          ! number of elements in ptab

    INTEGER                                 :: istatus, incid, id
    !!-----------------------------------------------------------
    istatus = NF90_OPEN(cdfile, NF90_WRITE, incid)
    istatus = NF90_INQ_VARID(incid, cdvar, id )
    istatus = NF90_PUT_VAR(incid, id, ptab, start=(/1/), count=(/kk/) )
    reputvar1d4 = istatus
    istatus = NF90_CLOSE(incid)

  END FUNCTION reputvar1d4

  INTEGER(KIND=4) FUNCTION reputvar1d8(cdfile, cdvar, ddtab, kk )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION reputvar1d8  ***
    !!
    !! ** Purpose : Copy 1d variable cdfvar in cdfile, an already existing file
    !!              tab is the 1d array to write and kk the size of ptab
    !!
    !! ** Method  :   double precision
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),            INTENT(in) :: cdfile      ! filename
    CHARACTER(LEN=*),            INTENT(in) :: cdvar       ! variable name
    REAL(KIND=8), DIMENSION(kk), INTENT(in) :: ddtab       ! 1D array to write in file
    INTEGER(KIND=4),             INTENT(in) :: kk          ! number of elements in ptab

    INTEGER                                 :: istatus, incid, id
    !!-----------------------------------------------------------
    istatus = NF90_OPEN(cdfile, NF90_WRITE, incid)
    istatus = NF90_INQ_VARID(incid, cdvar, id )
    istatus = NF90_PUT_VAR(incid, id, ddtab, start=(/1/), count=(/kk/) )
    reputvar1d8 = istatus
    istatus = NF90_CLOSE(incid)

  END FUNCTION reputvar1d8


  INTEGER(KIND=4) FUNCTION putvar0dt(kout, kid, pvalue, ktime)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION putvar0dt  ***
    !!
    !! ** Purpose : Copy single value, with id varid, into file id kout
    !!
    !! ** Method  :  use argument as dummy array(1,1) 
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),              INTENT(in) :: kout   ! ncid of output file
    INTEGER(KIND=4),              INTENT(in) :: kid    ! id of the variable
    REAL(KIND=4), DIMENSION(1,1), INTENT(in) :: pvalue ! single value to write in file
    INTEGER(KIND=4), OPTIONAL,    INTENT(in) :: ktime  ! time frame to write

    INTEGER(KIND=4) :: istatus
    INTEGER(KIND=4) :: itime
    INTEGER(KIND=4) :: id_dimunlim, inbdim
    INTEGER(KIND=4), DIMENSION(4) :: inldim, istart, icount
    !!----------------------------------------------------------------------
    IF (PRESENT(ktime) ) THEN
      itime = ktime
    ELSE
      itime = 1
    ENDIF 
    
     ! look for unlimited dim (time_counter)
    istatus=NF90_INQUIRE(kout, unlimitedDimId=id_dimunlim)
    istatus=NF90_INQUIRE_VARIABLE(kout,kid,ndims=inbdim,dimids=inldim(:) )

    istart(:) = 1
    icount(:) = 1
    IF ( inldim(inbdim) == id_dimunlim ) THEN ! var has an unlimited dim
       istart(inbdim) = itime ! assume than last dim is UNLIM 
    ENDIF

    istatus=NF90_PUT_VAR(kout, kid, pvalue, start=istart, count=icount )

    putvar0dt=istatus

  END FUNCTION putvar0dt

  INTEGER(KIND=4) FUNCTION putvar0ds(kout, kid, pvalue, ktime)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION putvar0ds  ***
    !!
    !! ** Purpose : Copy single value, with id varid, into file id kout
    !!
    !! ** Method  : use argument as scalar
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),              INTENT(in) :: kout   ! ncid of output file
    INTEGER(KIND=4),              INTENT(in) :: kid    ! id of the variable
    REAL(KIND=4),                 INTENT(in) :: pvalue ! single value to write in file
    INTEGER(KIND=4), OPTIONAL,    INTENT(in) :: ktime  ! time frame to write

    INTEGER(KIND=4) :: istatus
    INTEGER(KIND=4) :: itime
    INTEGER(KIND=4) :: id_dimunlim, inbdim
    INTEGER(KIND=4), DIMENSION(4) :: inldim, istart, icount

    REAL(KIND=4), DIMENSION(1,1)             :: ztab   ! dummy array for PUT_VAR
    !!----------------------------------------------------------------------
    IF (PRESENT(ktime) ) THEN
      itime = ktime
    ELSE
      itime = 1
    ENDIF

     ! look for unlimited dim (time_counter)
    istatus=NF90_INQUIRE(kout, unlimitedDimId=id_dimunlim)
    istatus=NF90_INQUIRE_VARIABLE(kout,kid,ndims=inbdim,dimids=inldim(:) )

    istart(:) = 1
    icount(:) = 1
    IF ( inldim(inbdim) == id_dimunlim ) THEN ! var has an unlimited dim
       istart(inbdim) = itime ! assume than last dim is UNLIM 
    ENDIF
    ztab = pvalue

    istatus=NF90_PUT_VAR(kout, kid, ztab, start=istart, count=icount )

    putvar0ds=istatus

  END FUNCTION putvar0ds

  INTEGER(KIND=4) FUNCTION closeout(kout)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION closeout  ***
    !!
    !! ** Purpose :  close opened output files 
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) :: kout   ! ncid of file to be closed
    !!----------------------------------------------------------------------
    closeout=NF90_CLOSE(kout)

  END FUNCTION closeout

  INTEGER(KIND=4) FUNCTION ncopen(cdfile)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION  ncopen  ***
    !!
    !! ** Purpose : open file cdfile and return file ID
    !!
    !!---------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: cdfile ! file name

      INTEGER(KIND=4) :: istatus, incid
    !!---------------------------------------------------------------------
      istatus = NF90_OPEN(cdfile,NF90_WRITE,incid)

      ncopen=incid

  END FUNCTION ncopen

  SUBROUTINE ERR_HDL(kstatus)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ERR_HDL  ***
    !!
    !! ** Purpose : Error handle for NetCDF routine.
    !!              Stop if kstatus indicates error conditions.  
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) ::  kstatus
    !!----------------------------------------------------------------------
    IF (kstatus /=  NF90_NOERR ) THEN
       PRINT *, 'ERROR in NETCDF routine, status=',kstatus
       PRINT *,NF90_STRERROR(kstatus)
       STOP 98
    END IF

  END SUBROUTINE ERR_HDL


  SUBROUTINE gettimeseries (cdfile, cdvar, kilook, kjlook, klev, ldncdf)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE gettimeseries  ***
    !!
    !! ** Purpose : Display a 2 columns output ( time, variable) for
    !!              a given variable of a given file at a given point 
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),          INTENT(in) :: cdfile, cdvar
    INTEGER(KIND=4),           INTENT(in) :: kilook,kjlook
    INTEGER(KIND=4), OPTIONAL, INTENT(in) :: klev
    LOGICAL,         OPTIONAL, INTENT(in) :: ldncdf

    INTEGER(KIND=4)                         :: jt, jk
    INTEGER(KIND=4)                         :: iint
    INTEGER(KIND=4)                         :: istatus
    INTEGER(KIND=4)                         :: incid, id_t, id_var
    INTEGER(KIND=4)                         :: indim
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: ztime, zval
    REAL(KIND=4)                            :: ztmp  
    REAL(KIND=4)                            :: zao=0., zsf=1.0   !: add_offset, scale_factor
    CHARACTER(LEN=256)                      :: clname
    LOGICAL                                 :: ll_netcdf=.false.
    !!----------------------------------------------------------------------
    ! Klev can be used to give the model level we want to look at
    IF ( PRESENT(klev) ) THEN
       jk=klev
    ELSE
       jk=1
    ENDIF

    IF ( PRESENT(ldncdf) ) THEN
       ll_netcdf=ldncdf
    ELSE
       ll_netcdf=.false.
    ENDIF


    ! Open cdf dataset
    istatus=NF90_OPEN(cdfile,NF90_NOWRITE,incid)

    ! read time dimension
    istatus=NF90_INQ_DIMID(incid, cn_t, id_t)
    istatus=NF90_INQUIRE_DIMENSION(incid,id_t, len=iint)

    ! Allocate space
    ALLOCATE (ztime(iint), zval(iint) )

    ! gettime
    istatus=NF90_INQ_VARID(incid,cn_vtimec,id_var)
    istatus=NF90_GET_VAR(incid,id_var,ztime,(/1/),(/iint/) )

    ! read variable
    istatus=NF90_INQ_VARID(incid,cdvar,id_var)

    ! look for scale_factor and add_offset attribute:
    istatus=NF90_GET_ATT(incid,id_var,'add_offset',ztmp)
    IF ( istatus == NF90_NOERR ) zao = ztmp
    istatus=NF90_GET_ATT(incid,id_var,'scale_factor',ztmp)
    IF ( istatus == NF90_NOERR ) zsf = ztmp

    ! get number of dimension of the variable ( either x,y,t or x,y,z,t )
    istatus=NF90_INQUIRE_VARIABLE(incid,id_var, ndims=indim)
    IF ( indim == 3 ) THEN
       istatus=NF90_GET_VAR(incid,id_var,zval,(/kilook,kjlook,1/),(/1,1,iint/) )
    ELSE IF ( indim == 4 ) THEN
       istatus=NF90_GET_VAR(incid,id_var,zval,(/kilook,kjlook,jk,1/),(/1,1,1,iint/) )
    ELSE 
       PRINT *,'  ERROR : variable ',TRIM(cdvar),' has ', indim, &
            &       ' dimensions !. Only 3 or 4 supported'
       STOP 98
    ENDIF

    ! convert to physical values
    zval=zval*zsf + zao

    ! display results :
    DO jt=1,iint
       PRINT *,ztime(jt)/86400., zval(jt)
    ENDDO

    istatus=NF90_CLOSE(incid)
    IF ( ll_netcdf )  THEN
      WRITE (clname,'("probe_",i4.4,"_",i4.4,"_",a,".nc")') kilook, kjlook, TRIM(cdvar)
    ENDIF
    ! not finished ...

  END SUBROUTINE gettimeseries

  LOGICAL FUNCTION chkfile (cd_file, ld_verbose ) 
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION chkfile  ***
    !!
    !! ** Purpose :  Check if cd_file exists.
    !!               Return false if it exists, true if it does not
    !!               Do nothing is filename is 'none'
    !!
    !! ** Method  : Doing it this way allow statements such as
    !!              IF ( chkfile( cf_toto) ) STOP 99  ! missing file
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),  INTENT(in) :: cd_file
    LOGICAL, OPTIONAL, INTENT(in) :: ld_verbose

    INTEGER(KIND=4)               :: ierr
    LOGICAL                       :: ll_exist, ll_verbose
    !!----------------------------------------------------------------------
    IF ( PRESENT(ld_verbose) ) THEN
       ll_verbose = ld_verbose
    ELSE
       ll_verbose = .TRUE.
    ENDIF

    IF ( TRIM(cd_file) /= 'none')  THEN 
       INQUIRE (file = TRIM(cd_file), EXIST=ll_exist)

       IF (ll_exist) THEN
          chkfile = .false.
          IF ( cd_file == cn_fzgr ) ierr = SetMeshZgrVersion ()
       ELSE
          IF ( ll_verbose ) PRINT *, ' File ',TRIM(cd_file),' is missing '
          chkfile = .true.
       ENDIF
    ELSE  
       chkfile = .false.  ! 'none' file is not checked
    ENDIF

  END FUNCTION chkfile

  LOGICAL FUNCTION chkvar (cd_file, cd_var, ld_verbose )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION chkvar  ***
    !!
    !! ** Purpose :  Check if cd_var exists in file cd_file.
    !!               Return false if it exists, true if it does not
    !!               Do nothing is varname is 'none'
    !!
    !! ** Method  : Doing it this way allow statements such as
    !!              IF ( chkvar( cf_toto, cv_toto) ) STOP 99  ! missing var
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),  INTENT(in) :: cd_file
    CHARACTER(LEN=*),  INTENT(in) :: cd_var
    LOGICAL, OPTIONAL, INTENT(in) :: ld_verbose

    INTEGER(KIND=4)              :: istatus
    INTEGER(KIND=4)              :: incid, id_t, id_var
    LOGICAL                      :: ll_verbose
    !!----------------------------------------------------------------------
    IF ( TRIM(cd_var) /= 'none')  THEN
       IF ( PRESENT(ld_verbose) ) THEN
          ll_verbose = ld_verbose
       ELSE
          ll_verbose = .TRUE.
       ENDIF
    
       ! Open cdf dataset
       istatus = NF90_OPEN(cd_file, NF90_NOWRITE,incid)
       ! Read variable
       istatus = NF90_INQ_VARID(incid, cd_var, id_var)

       IF ( istatus == NF90_NOERR ) THEN
          chkvar = .false.
       ELSE
          IF ( ll_verbose ) PRINT *, ' Var ',TRIM(cd_var),' is missing in file ',TRIM(cd_file)
          chkvar = .true.
       ENDIF
       
       ! Close file
       istatus = NF90_CLOSE(incid) 
    ELSE
       chkvar = .false.  ! 'none' file is not checked
    ENDIF

  END FUNCTION chkvar

  CHARACTER(LEN=256) FUNCTION Get_Env ( cd_env )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION Get_Env  ***
    !!
    !! ** Purpose :  A wrapper for system routine getenv
    !!
    !! ** Method  :  Call getenv
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_env
    !!----------------------------------------------------------------------
    CALL getenv( TRIM(cd_env), Get_Env )
    IF ( TRIM(Get_Env) /= '' ) THEN
      PRINT *,'Environment found : ',TRIM(cd_env),' = ', TRIM(Get_Env)
    ENDIF

  END FUNCTION Get_Env

  FUNCTION GetNcFile (cd_file)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION GetNcFile  ***
    !!
    !! ** Purpose :  fills in the ncfile  structure corresponding to the file
    !!               given in argument 
    !!
    !! ** Method  :  Use NF90 function to get the ad-hoc information 
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_file
    TYPE(ncfile)                 :: GetNcFile

    INTEGER(KIND=4) :: jvar, jdim                      ! loop index
    INTEGER(KIND=4) :: ierr, idx, idy, idz, idt, idb   ! error status and dimids
    !!----------------------------------------------------------------------
    GetNcFile%c_fnam = cd_file
    ierr = NF90_OPEN(cd_file, NF90_NOWRITE, GetNcFile%ncid )
    ierr = NF90_INQUIRE(GetNcFile%ncid, nDimensions    = GetNcFile%ndims,  &
         &                              nVariables     = GetNcFile%nvars,  &
         &                              nAttributes    = GetNcFile%natts,  &
         &                              unlimitedDimId = GetNcFile%iunlim  )
    ALLOCATE (GetNcFile%nvdim   (GetNcFile%nvars) )
    ALLOCATE (GetNcFile%nvid    (GetNcFile%nvars) )
    ALLOCATE (GetNcFile%c_vnam  (GetNcFile%nvars) )
    ALLOCATE (GetNcFile%nvatt   (GetNcFile%nvars) )
    ALLOCATE (GetNcFile%itype   (GetNcFile%nvars) )
    ALLOCATE (GetNcFile%c_dnam  (GetNcFile%ndims) )
    ALLOCATE (GetNcFile%nlen    (GetNcFile%ndims) )
    ALLOCATE (GetNcFile%idimids (GetNcFile%nvars,GetNcFile%ndims) )

    ! Look for dimensions
    DO jdim = 1, GetNcFile%ndims
       ierr = NF90_INQUIRE_DIMENSION(GetNcFile%ncid,jdim,                 &
            &                         name   = GetNcFile%c_dnam(jdim),    &
            &                         len    = GetNcFile%nlen  (jdim) )
    ENDDO
    ! Look for variables
    DO jvar = 1, GetNcFile%nvars
       ierr = NF90_INQUIRE_VARIABLE (GetNcFile%ncid, jvar,                &
            &                         name   = GetNcFile%c_vnam(jvar),    &
            &                         xtype  = GetNcFile%itype(jvar),     &
            &                         ndims  = GetNcFile%nvdim(jvar),     &
            &                         dimids = GetNcFile%idimids(jvar,:), &
            &                         nAtts  = GetNcFile%nvatt(jvar),     &
            &                         contiguous = GetNcFile%lconti(jvar),   &
            &                         chunksizes = GetNcFile%ichunk(jvar,:), &
            &                         deflate_level = GetNcFile%ideflat(jvar)      )
    END DO
    ! Look for attributes
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_number_total'   , GetNcFile%number_total       )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_number      '   , GetNcFile%number             )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_dimensions_ids' , GetNcFile%idimensions_ids(:) )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_size_global'    , GetNcFile%isize_global(:)    )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_size_local'     , GetNcFile%isize_local(:)     )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_position_first' , GetNcFile%iposition_first(:) )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_position_last'  , GetNcFile%iposition_last(:)  )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_halo_size_start', GetNcFile%ihalo_size_start(:))
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_halo_size_end'  , GetNcFile%ihalo_size_end(:)  )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_type'           , GetNcFile%c_type             )

    ! NOTE : for recombined files, no more DOMAIN attributes !!!
    ! DOMAIN_dimensions_ids gives ids for x, y 
    idx = GetNcFile%idimensions_ids(1)
    idy = GetNcFile%idimensions_ids(2)
    ! time is unlimited dim
    idt = GetNcFile%iunlim

    ! try to infer size of the domain assuming some basis:
    ! (1) 2D var are (x,y)
    ! (2) time dim is unlimited
    ! (3) allowed shape of var : x,y ; x,y,t ; x,y,z,t  ; x,y,z   ; t   ; z 

    ! Look for x y z t tbound dim id
    idx=-1 ; idy=-1 ; idz=-1 ; idt=-1 ; idb=-1

    DO jvar = 1, GetNcFile%nvars
       GetNcFile%nvid(jvar) = jvar
       IF ( GetNcFile%nvdim(jvar) == 2 ) THEN  
          IF ( GetNcFile%idimids(jvar,2) == GetNcFile%iunlim ) THEN  ! catch a time_bounds, time_counter var
            idb = GetNcFile%idimids(jvar,1)
          ELSE
          idx = GetNcFile%idimids(jvar,1) 
          idy = GetNcFile%idimids(jvar,2) 
          ENDIF
       ELSE IF ( GetNcFile%nvdim(jvar) == 3 ) THEN 
          idx = GetNcFile%idimids(jvar,1) 
          idy = GetNcFile%idimids(jvar,2) 
          IF ( GetNcFile%idimids(jvar,3) == GetNcFile%iunlim ) THEN
             idt = GetNcFile%idimids(jvar,3) 
          ELSE
             idz = GetNcFile%idimids(jvar,3) 
          ENDIF
       ELSE IF ( GetNcFile%nvdim(jvar) == 4 ) THEN
          idx = GetNcFile%idimids(jvar,1)  
          idy = GetNcFile%idimids(jvar,2)  
          idz = GetNcFile%idimids(jvar,3)  
          IF ( GetNcFile%idimids(jvar,4) /= GetNcFile%iunlim ) THEN
             PRINT *, ' 4D variables must have an unlimited time dimension ...'
             PRINT *, ' Cannot process this file :', TRIM(cd_file)
             STOP 98
          ENDIF
          idt = GetNcFile%idimids(jvar,4) 
       ENDIF
    END DO

    GetNcFile%idx=idx
    GetNcFile%idy=idy
    GetNcFile%idz=idz
    GetNcFile%idt=idt
    GetNcFile%idb=idb

    IF ( idx == -1 .OR. idy == -1 ) THEN 
       PRINT *, ' ERROR : no x, y dimensions found'
       STOP 98
    ENDIF

    ! get dimensions
    GetNcFile%npi=0 ; GetNcFile%npj=0 ; GetNcFile%npk=0 ; GetNcFile%npt=0 ; GetNcFile%npb=0

    ierr = NF90_INQUIRE_DIMENSION( GetNcFile%ncid, idx,         &
         &                        name = GetNcFile%c_dnam(idx), &
         &                        len  = GetNcFile%npi          )
    GetNcFile%nlen(idx) = GetNcFile%npi
    ierr = NF90_INQUIRE_DIMENSION( GetNcFile%ncid, idy,         &
         &                        name = GetNcFile%c_dnam(idy), &
         &                        len  = GetNcFile%npj          )
    GetNcFile%nlen(idy) = GetNcFile%npj

    IF ( idz /= -1 ) THEN
       ierr = NF90_INQUIRE_DIMENSION( GetNcFile%ncid, idz,         &
            &                        name = GetNcFile%c_dnam(idz), &
            &                        len  = GetNcFile%npk          )
       GetNcFile%nlen(idz) = GetNcFile%npk
    ENDIF

    IF ( idt /= -1 ) THEN
       ierr = NF90_INQUIRE_DIMENSION( GetNcFile%ncid, idt,         &
            &                        name = GetNcFile%c_dnam(idt), &
            &                        len  = GetNcFile%npt          )
       GetNcFile%nlen(idt) = GetNcFile%npt
    ENDIF

    IF ( idb /= -1 ) THEN
       ierr = NF90_INQUIRE_DIMENSION( GetNcFile%ncid, idb,         &
            &                        name = GetNcFile%c_dnam(idb), &
            &                        len  = GetNcFile%npb          )
       GetNcFile%nlen(idb) = GetNcFile%npb
    ENDIF
    ! fill in DOMAIN attributes 
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_number_total'   , GetNcFile%number_total       )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_number      '   , GetNcFile%number             )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_dimensions_ids' , GetNcFile%idimensions_ids(:) )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_size_global'    , GetNcFile%isize_global(:)    )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_size_local'     , GetNcFile%isize_local(:)     )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_position_first' , GetNcFile%iposition_first(:) )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_position_last'  , GetNcFile%iposition_last(:)  )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_halo_size_start', GetNcFile%ihalo_size_start(:))
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_halo_size_end'  , GetNcFile%ihalo_size_end(:)  )
    ierr = NF90_GET_ATT (GetNcFile%ncid, NF90_GLOBAL, 'DOMAIN_type'           , GetNcFile%c_type             )

    ! CORRECT for Halo :
    GetNcfile%isize_local(:)     = GetNcfile%isize_local(:)     - GetNcFile%ihalo_size_start(:) - GetNcFile%ihalo_size_end(:)
    GetNcFile%iposition_first(:) = GetNcFile%iposition_first(:) + GetNcFile%ihalo_size_start(:)
    GetNcFile%iposition_last(:)  = GetNcFile%iposition_last(:)  - GetNcFile%ihalo_size_end(:)

    GetNcFile%npi       = GetNcfile%isize_local(1)
    GetNcFile%npj       = GetNcfile%isize_local(2)
    GetNcFile%nlen(idx) = GetNcFile%npi
    GetNcFile%nlen(idy) = GetNcFile%npj

  END FUNCTION GetNcFile

  INTEGER FUNCTION  SetMeshZgrVersion()
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE SetMeshZgrVersion  ***
    !!
    !! ** Purpose :  This routine set the global variable cg_zgr_ver 
    !!              according to the format of the set of mesh_mask file
    !!
    !! ** Method  : Use the algorithm formelly in getvar with option ldiom=true
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=10)  :: clvar='e3t_0'
    !!----------------------------------------------------------------------
   
    IF ( chkvar (cn_fzgr, clvar, .false.) ) THEN  ! use quiet mode
         cg_zgr_ver='v2.0'
    ELSE
      IF ( getvdim (cn_fzgr, clvar) == 1 ) THEN
         cg_zgr_ver='v3.0'
      ELSE
         cg_zgr_ver='v3.6'
      ENDIF
    ENDIF
    PRINT *,' mesh_zgr version is ', TRIM( cg_zgr_ver )
    SetMeshZgrVersion=1

  END FUNCTION SetMeshZgrVersion


END MODULE cdfio

