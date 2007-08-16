MODULE cdfio
  !!---------------------------------------------------------------------------------------------------
  !!                     ***  MODULE  cdfio  ***
  !!
  !!    ** Purpose : this module manage all the IO with Netcdf Library
  !!
  !!    ** Method : provide functions that are used in the different
  !!                 subprograms for performing dedicated tasks
  !!  
  !!   history:
  !!         Original : J.M. Molines (2005 )
  !!------------------------------------------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  USE netcdf 

  IMPLICIT NONE
  INTEGER :: id_x, id_y, id_z, id_t, id_lat, id_lon, id_dep, id_tim

  TYPE, PUBLIC ::   variable 
     CHARACTER(LEN=80)::  name
     CHARACTER(LEN=80):: units
     REAL(kind=4)    :: missing_value
     REAL(kind=4)    :: valid_min
     REAL(kind=4)    :: valid_max
     REAL(kind=4)    :: scale_factor=1.
     REAL(kind=4)    :: add_offset=0.
     REAL(kind=4)    :: savelog10=0.
     CHARACTER(LEN=80):: long_name
     CHARACTER(LEN=80):: short_name
     CHARACTER(LEN=80):: online_operation
     CHARACTER(LEN=80):: axis
     CHARACTER(LEN=80):: PRECISION='r4'  ! possible values are i2, r4, r8
  END TYPE variable

  INTERFACE putvar
     MODULE PROCEDURE putvarr4, putvari2, putvarzo, reputvarr4
  END INTERFACE


  PRIVATE 
  PUBLIC  copyatt, create, createvar, getvaratt,cvaratt
  PUBLIC  putatt, putheadervar, putvar, putvar1d
  PUBLIC  getatt, getdim, getvdim, getipk, getnvar, getvarname, getvarid, getspval
  PUBLIC  getvar, getvarxz, getvaryz, getvar1d, getvare3
  PUBLIC gettimeseries
  PUBLIC closeout, ncopen
  PUBLIC ERR_HDL


CONTAINS
  FUNCTION copyatt(cdvar,kidvar,kcin,kcout)
    !! ----------------------------------------------------------------------------------------------------
    !!  *** Copy attributes for variable cdvar, which have id kidvar in kcout, from file id kcin
    !!
    !! ----------------------------------------------------------------------------------------------------
    ! * Arguments
    INTEGER, INTENT(in) :: kidvar, kcout
    INTEGER,  INTENT(in) :: kcin
    CHARACTER(LEN=*), INTENT(in) :: cdvar
    INTEGER :: copyatt

    ! * Local variable
    INTEGER :: istatus, idvar, iatt, ja
    CHARACTER(LEN=80) :: clatt

    IF ( kcin /= -9999) THEN
       istatus = NF90_INQ_VARID(kcin,cdvar,idvar)
       istatus = NF90_INQUIRE_VARIABLE(kcin,idvar,natts=iatt)
       DO ja = 1, iatt
          istatus = NF90_INQ_ATTNAME(kcin,idvar,ja,clatt)
          istatus = NF90_COPY_ATT(kcin,idvar,clatt,kcout,kidvar)
       END DO
    ELSE
       SELECT CASE (cdvar )
       CASE ('nav_lon' )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'units', 'degrees_east')
          istatus=NF90_PUT_ATT(kcout, kidvar, 'valid_min', -180.)
          istatus=NF90_PUT_ATT(kcout, kidvar, 'valid_max', 180.)
          istatus=NF90_PUT_ATT(kcout, kidvar, 'long_name', 'Longitude')
          istatus=NF90_PUT_ATT(kcout, kidvar, 'nav_model', 'Default grid')
       CASE ('nav_lat' )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'units', 'degrees_north')
          istatus=NF90_PUT_ATT(kcout, kidvar, 'valid_min', -90.)
          istatus=NF90_PUT_ATT(kcout, kidvar, 'valid_max', 90.)
          istatus=NF90_PUT_ATT(kcout, kidvar, 'long_name', 'Latitude')
          istatus=NF90_PUT_ATT(kcout, kidvar, 'nav_model', 'Default grid')
       CASE ('time_counter' )
          istatus=NF90_PUT_ATT(kcout, kidvar, 'calendar', 'gregorian')
          istatus=NF90_PUT_ATT(kcout, kidvar, 'units', 'seconds since 0006-01-01 00:00:00')
          istatus=NF90_PUT_ATT(kcout, kidvar, 'time_origin', '0001-JAN-01 00:00:00')
          istatus=NF90_PUT_ATT(kcout, kidvar, 'title', 'Time')
          istatus=NF90_PUT_ATT(kcout, kidvar, 'long_name', 'Time axis')
       CASE ('deptht', 'depthu' ,'depthv' , 'depthw', 'dep')
          istatus=NF90_PUT_ATT(kcout, kidvar, 'units', 'm')
          istatus=NF90_PUT_ATT(kcout, kidvar, 'positive', 'unknown')
          istatus=NF90_PUT_ATT(kcout, kidvar, 'valid_min', 0.)
          istatus=NF90_PUT_ATT(kcout, kidvar, 'valid_max', 5875.)
          istatus=NF90_PUT_ATT(kcout, kidvar, 'title', TRIM(cdvar))
          istatus=NF90_PUT_ATT(kcout, kidvar, 'long_name', 'Vertical Levels')
       END SELECT
    ENDIF

    copyatt = istatus
  END FUNCTION copyatt


  FUNCTION create( cdfile, cdfilref ,kx,ky,kz ,cdep)
    !! ------------------------------------------------------------------------------------------
    !! ***  Create the file, and creates dimensions, and copy attributes from a cdilref
    !!      reference file ( for the nav_lon, nav_lat etc ...)
    !!      If optional cdep given : take as depth variable name instead of cdfilref
    !!      Return the nc id of the created file, and leave it open
    !!
    !! ------------------------------------------------------------------------------------------
    ! * Arguments
    CHARACTER(LEN=*), INTENT(in) :: cdfile,cdfilref
    INTEGER, INTENT(in) :: kx,ky,kz
    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: cdep 
    INTEGER :: create

    ! * Local Variable
    INTEGER   :: istatus, icout,ncid, idum
    INTEGER ,DIMENSION(4) :: nvdim
    CHARACTER(LEN=80) :: cldep, cldepref

    istatus = NF90_CREATE(cdfile,NF90_CLOBBER, icout)
    ! Open reference file if any,  otherwise set ncid to flag value (for copy att)
    IF ( TRIM(cdfilref) /= 'none' ) THEN
       istatus = NF90_OPEN(cdfilref,NF90_NOWRITE,ncid)
    ELSE
       ncid = -9999
    ENDIF

    istatus = NF90_DEF_DIM(icout,'x',kx, id_x)
    istatus = NF90_DEF_DIM(icout,'y',ky, id_y)

    IF ( kz /= 0 ) THEN
       ! try to find out the name I will use for depth dimension in the new file ...
       IF (PRESENT (cdep) ) THEN
          cldep = cdep
          idum=getdim(cdfilref,cldep,cldepref)   ! look for depth dimension name in ref file
       ELSE 
          idum=getdim(cdfilref,'depth',cldep   )   ! look for depth dimension name in ref file
          cldepref=cldep
       ENDIF
       istatus = NF90_DEF_DIM(icout,TRIM(cldep),kz, id_z)
    ENDIF

    istatus = NF90_DEF_DIM(icout,'time_counter',NF90_UNLIMITED, id_t)

    nvdim(1) = id_x ; nvdim(2) = id_y ; nvdim(3) = id_z ; nvdim(4) = id_t

    ! define variables and copy attributes
    istatus = NF90_DEF_VAR(icout,'nav_lon',NF90_FLOAT,(/id_x,id_y/),id_lon)
    istatus = copyatt('nav_lon',id_lon,ncid,icout)
    istatus = NF90_DEF_VAR(icout,'nav_lat',NF90_FLOAT,(/id_x,id_y/),id_lat)
    istatus = copyatt('nav_lat',id_lat,ncid,icout)
    IF ( kz /= 0 ) THEN
       ! here we assume that dep variable has same name as dep dim .... jmm
       istatus = NF90_DEF_VAR(icout,TRIM(cldep),NF90_FLOAT,(/id_z/),id_dep)
       ! JMM bug fix : if cdep present, then chose attribute from cldepref
       istatus = copyatt(TRIM(cldepref),id_dep,ncid,icout)
    ENDIF

    istatus = NF90_DEF_VAR(icout,'time_counter',NF90_FLOAT,(/id_t/),id_tim)
    istatus = copyatt('time_counter',id_tim,ncid,icout)

    istatus = NF90_CLOSE(ncid)

    create=icout
  END FUNCTION create

  FUNCTION createvar(kout,ptyvar,kvar,kpk, kidvo)
    !! ----------------------------------------------------------------------------------------------------
    !!  *** Create kvar n-2D variables cdvar(:), in file id kout, kpk gives the number of vertical levels
    !!      idvo(:) contains the id of the crated variables.
    !!     INPUT:
    !!       kout = ncid of output file
    !!       cdvar= array of name of variables
    !!       kvar = number of variables to create
    !!       kpk  = number of vertical dimensions foreach variable
    !!
    !!     OUTPUT:
    !!       kidvo = arrays with the varid of the variables just created.
    !!
    !! ----------------------------------------------------------------------------------------------------
    ! * Arguments
    INTEGER, INTENT(in) :: kout, kvar
    INTEGER, DIMENSION(kvar), INTENT(in) :: kpk
    INTEGER, DIMENSION(kvar), INTENT(out) :: kidvo
    INTEGER :: createvar
    TYPE (variable), DIMENSION(kvar) ,INTENT(in) :: ptyvar

    ! * Local variables
    INTEGER :: jv,idims, istatus
    INTEGER, DIMENSION(4):: iidims

    DO jv = 1, kvar

       ! Create variables whose name is not 'none'
       IF ( ptyvar(jv)%name /= 'none' ) THEN
          IF (kpk(jv) == 1 ) THEN
             idims=3
             iidims(1) = id_x ; iidims(2) = id_y ; iidims(3) = id_t
          ELSE IF (kpk(jv) > 1 ) THEN
             idims=4
             iidims(1) = id_x ; iidims(2) = id_y ; iidims(3) = id_z ; iidims(4) = id_t
          ELSE
             PRINT *,' ERROR: ipk = ',kpk(jv), jv , ptyvar(jv)%name
             STOP
          ENDIF

          IF ( ptyvar(jv)%precision == 'r8' ) THEN
             istatus = NF90_DEF_VAR(kout,ptyvar(jv)%name,NF90_DOUBLE,iidims(1:idims) ,kidvo(jv) )
          ELSE
             IF ( ptyvar(jv)%scale_factor == 1. .AND. ptyvar(jv)%add_offset == 0. ) THEN
                istatus = NF90_DEF_VAR(kout,ptyvar(jv)%name,NF90_FLOAT,iidims(1:idims) ,kidvo(jv) )
             ELSE
                istatus = NF90_DEF_VAR(kout,ptyvar(jv)%name,NF90_SHORT,iidims(1:idims) ,kidvo(jv) )
             ENDIF
          ENDIF

          ! add attributes
          istatus = putatt(ptyvar(jv), kout,kidvo(jv))
          createvar=istatus
       ENDIF
    END DO
    istatus = NF90_ENDDEF(kout)

  END FUNCTION createvar

  FUNCTION getvarid( cdfile, knvars )
    !! ------------------------------------------------------------------------------------------
    !! ***  return a real array with the nvar variable id
    !!
    !! ------------------------------------------------------------------------------------------
    ! * Arguments
    CHARACTER(LEN=*), INTENT(in) :: cdfile
    INTEGER, INTENT(in)  ::  knvars                  ! Number of variables in cdfile
    INTEGER, DIMENSION(knvars) :: getvarid

    !! * local declarations
    CHARACTER(LEN=80), DIMENSION(knvars) :: cdvar
    INTEGER :: ncid, jv
    INTEGER :: istatus


    istatus = NF90_OPEN(cdfile,NF90_NOWRITE,ncid)
    DO jv = 1, knvars
       istatus = NF90_INQUIRE_VARIABLE(ncid,jv,cdvar(jv) )
       istatus = NF90_INQ_VARID(ncid,cdvar(jv),getvarid(jv))
    ENDDO
    istatus=NF90_CLOSE(ncid)

  END FUNCTION getvarid

  FUNCTION getvaratt (cdfile,cdvar,cdunits, pmissing_value, cdlong_name, cdshort_name)
    !! ----------------------------------------------------------------------------------------------------
    !!  ***  Change variable attributs in an existing variable
    !!
    !! ----------------------------------------------------------------------------------------------------
    ! * Arguments
    CHARACTER(LEN=80), INTENT(in) :: cdfile, cdvar
    CHARACTER(LEN=80), INTENT(out) :: cdunits, cdlong_name, cdshort_name
    REAL(KIND=4), INTENT(out) :: pmissing_value
    INTEGER :: getvaratt

    !! * local declarations
    INTEGER :: istatus
    INTEGER :: ncid, varid

    istatus = NF90_OPEN(cdfile,NF90_NOWRITE,ncid)
    istatus = NF90_INQ_VARID(ncid,cdvar,varid)

    istatus=NF90_GET_ATT(ncid, varid, 'units', cdunits)
    istatus=NF90_GET_ATT(ncid, varid, 'missing_value', pmissing_value)
    istatus=NF90_GET_ATT(ncid, varid, 'long_name', cdlong_name)
    istatus=NF90_GET_ATT(ncid, varid, 'short_name', cdshort_name)

!   istatus = NF90_ENDDEF(ncid)
    getvaratt=istatus
    istatus=NF90_CLOSE(ncid)

  END FUNCTION getvaratt


  FUNCTION cvaratt (cdfile,cdvar,cdunits,pmissing_value, cdlong_name, cdshort_name)
    !! ----------------------------------------------------------------------------------------------------
    !!  ***  Change variable attributs in an existing variable
    !!
    !! ----------------------------------------------------------------------------------------------------
    ! * Arguments
    CHARACTER(LEN=80), INTENT(in) :: cdfile, cdvar
    CHARACTER(LEN=80), INTENT(in) :: cdunits, cdlong_name, cdshort_name
    INTEGER :: cvaratt
    REAL(KIND=4) :: pmissing_value

    !! * local declarations
    INTEGER :: istatus
    INTEGER :: ncid, varid

    istatus = NF90_OPEN(cdfile,NF90_WRITE,ncid)
    istatus = NF90_REDEF(ncid)
    istatus = NF90_INQ_VARID(ncid,cdvar,varid)

    istatus=NF90_RENAME_ATT(ncid, varid, 'units', cdunits)
    istatus=NF90_PUT_ATT(ncid, varid, 'missing_value', pmissing_value)
    istatus=NF90_RENAME_ATT(ncid, varid, 'long_name', cdlong_name)
    istatus=NF90_RENAME_ATT(ncid, varid, 'short_name', cdshort_name)

    istatus=NF90_ENDDEF(ncid)
    cvaratt=istatus
    istatus=NF90_CLOSE(ncid)

  END FUNCTION cvaratt


  FUNCTION putatt (tyvar,kout,kid)
    !! ----------------------------------------------------------------------------------------------------
    !!  ***  Scan file att.txt for finding the line corresponding to cdvar, then read the attributes
    !!       for this variables ,whose id is kid and  write them in file id kout
    !!
    !! ----------------------------------------------------------------------------------------------------
    ! * Arguments
    INTEGER :: putatt
    INTEGER, INTENT(in) :: kout, kid
    TYPE (variable) ,INTENT(in) :: tyvar
    putatt=NF90_PUT_ATT(kout,kid,'units',tyvar%units) 
    IF (putatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(putatt)  ; STOP 'putatt'; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'missing_value',tyvar%missing_value)  
    IF (putatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(putatt)  ; STOP 'putatt'; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'valid_min',tyvar%valid_min) 
    IF (putatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(putatt)  ; STOP 'putatt'; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'valid_max',tyvar%valid_max)
    IF (putatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(putatt)  ; STOP 'putatt'; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'long_name',tyvar%long_name)
    IF (putatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(putatt)  ; STOP 'putatt'; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'short_name',tyvar%short_name) 
    IF (putatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(putatt)  ; STOP 'putatt'; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'online_operation',tyvar%online_operation) 
    IF (putatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(putatt)  ; STOP 'putatt'; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'axis',tyvar%axis) 
    IF (putatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(putatt)  ; STOP 'putatt'; ENDIF
    ! Optional attributes (scale_factor, add_offset )
    putatt=NF90_PUT_ATT(kout,kid,'scale_factor',tyvar%scale_factor) 
    IF (putatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(putatt)  ; STOP 'putatt'; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'add_offset',tyvar%add_offset) 
    IF (putatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(putatt)  ; STOP 'putatt'; ENDIF
    putatt=NF90_PUT_ATT(kout,kid,'savelog10',tyvar%savelog10) 
    IF (putatt /= 0 ) THEN ;PRINT *, NF90_STRERROR(putatt)  ; STOP 'putatt'; ENDIF

  END FUNCTION putatt

  FUNCTION getatt(cdfile,cdvar,cdatt)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  getatt  ***
    !!
    !! ** Purpose : return a REAL value with the values of the 
    !!              attribute cdatt for all the variable cdvar  in cdfile
    !!  
    !! ** Method : open, read attribute close
    !!
    !! history:
    !!    27/04/2005 : Jean-Marc Molines : Original Code 
    !!    12/03/2007 :  J.M. Molines : modif 
    !!-----------------------------------------------------------
    !! * Arguments declarations

    CHARACTER(LEN=*), INTENT(in) :: cdatt,   &   ! attribute name to look for
         &                          cdfile,  &   ! file to look at
         &                          cdvar

    REAL(KIND=4) :: getatt

    !! * Local declarations

    INTEGER :: istatus, jv, ncid, idum
    !! ----------------------------------------------------------
    istatus=NF90_OPEN(cdfile,NF90_NOWRITE,ncid)
    istatus=NF90_INQ_VARID(ncid,cdvar,idum)
    IF ( istatus /= NF90_NOERR) PRINT *, TRIM(NF90_STRERROR(istatus)),' when looking for ',TRIM(cdvar),' in getatt.'
    istatus = NF90_GET_ATT(ncid, idum,cdatt, getatt)
    IF ( istatus /= NF90_NOERR ) THEN
       PRINT *,' getatt problem :',NF90_STRERROR(istatus)
       PRINT *,' attribute :', TRIM(cdatt)
       PRINT *,' return default 0 '
       getatt=0.
    ENDIF
    istatus=NF90_CLOSE(ncid)

  END FUNCTION getatt

  FUNCTION  getdim (cdfile,cdim_name,cdtrue,kstatus)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  getdim  ***
    !!
    !! ** Purpose : return the INTEGER value of the dimension
    !!              identified with cdim_name in cdfile
    !!  
    !! ** Method  : Scan all the dimension name in cdfile and
    !!              select the one which match cdim_name.
    !!              cdim_name can be only a fraction of the total name
    !!              (eg: depth  will be ok for depht, or dephu, or dephv )
    !!
    !! history:
    !!    27/04/2005 : Jean-Marc Molines : Original Code 
    !!-----------------------------------------------------------
    !! * Arguments declarations
    CHARACTER(LEN=*), INTENT(in) :: cdfile , &  ! File name to look at
         &                           cdim_name   ! dimension name to look at
    CHARACTER(LEN=80),OPTIONAL, INTENT(out) ::  cdtrue ! full name of the read dimension
    INTEGER, OPTIONAL, INTENT(out) :: kstatus   ! status of the nf inquire
    INTEGER :: getdim                           ! the value for dim cdim_name, in file cdfile

    ! * Local variables
    INTEGER :: ncid, id_var
    INTEGER :: istatus
    INTEGER :: idims
    CHARACTER(LEN=80) :: clnam
    clnam = '-------------'

    IF ( PRESENT(kstatus) ) kstatus=0
    istatus=NF90_OPEN(cdfile,NF90_NOWRITE,ncid)
    IF ( istatus == NF90_NOERR ) THEN
       istatus=NF90_INQUIRE(ncid,ndimensions=idims)

       id_var = 1
       ! Look for dim name containing at least 'cdim_name'
       !     DO WHILE ( INDEX(clnam,cdim_name) == 0 .AND. id_var <= idims )
       !      istatus=NF90_INQUIRE_DIMENSION(ncid,id_var,name=clnam,len=getdim)
       !      id_var = id_var + 1
       !     END DO

       DO id_var = 1,idims
          istatus=NF90_INQUIRE_DIMENSION(ncid,id_var,name=clnam,len=getdim)
          IF ( INDEX(clnam,cdim_name) /= 0 ) THEN
             IF ( PRESENT(cdtrue) ) cdtrue=clnam
             EXIT
          ENDIF
       ENDDO

       IF ( id_var > idims ) THEN
          !      PRINT *,' warning: problem in getdim for ', TRIM(cdim_name),' in ', TRIM(cdfile)
          IF ( PRESENT(kstatus) ) kstatus=1    ! error send optionally to the calling program
          getdim=0
          IF ( PRESENT(cdtrue) ) cdtrue='unknown'
       ENDIF
       istatus=NF90_CLOSE(ncid)
    ELSE
       IF ( PRESENT(cdtrue) ) cdtrue='unknown'
       IF ( PRESENT(kstatus) ) kstatus=1 
    ENDIF

  END FUNCTION getdim

  FUNCTION  getspval (cdfile,cdvar)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  getspval  ***
    !!
    !! ** Purpose : return the SPVAL value of the variable
    !!              cdvar  in cdfile
    !!
    !!-----------------------------------------------------------
    !! * Arguments declarations
    CHARACTER(LEN=*), INTENT(in) :: cdfile , &  ! File name to look at
         &                           cdvar      ! variable name
    REAL(KIND=4) :: getspval                               ! the missing value for cdvar

    ! * Local variables
    INTEGER :: ncid, id_var
    INTEGER :: istatus

    istatus=NF90_OPEN(cdfile,NF90_NOWRITE,ncid)
    istatus=NF90_INQ_VARID ( ncid,cdvar,id_var)
    istatus=NF90_GET_ATT(ncid,id_var,"missing_value",getspval)
    istatus=NF90_CLOSE(ncid)

  END FUNCTION getspval

  FUNCTION getvdim (cdfile, cdvar)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  getvdim  ***
    !!
    !! ** Purpose : return the number of dimensions for variable cdvar in cdfile
    !!
    !! ** Method  : Inquire for variable cdvar in cdfile. If found,
    !!              determines the number of dimensions , assuming that variables
    !!              are either (x,y,dep,time) or (x,y,time)
    !!              If cdvar is not found, give an interactive choice for an existing
    !!              variable, cdvar is then updated to this correct name.
    !!
    !! history:
    !!    31/10/2005 : Jean-Marc Molines : Original Code
    !!-----------------------------------------------------------
    !! * Arguments declarations
    CHARACTER(LEN=*), INTENT(in) :: cdfile      ! File name to look at
    CHARACTER(LEN=*), INTENT(inout) :: cdvar    ! variable name to look at.
    INTEGER :: getvdim                          ! number of lebvels for cdvar

    !! * Local variables
    INTEGER :: istatus, ncid, id_var, ivar, idi, istatus0
    INTEGER :: jvar
    CHARACTER(LEN=80) :: clongname='long_name', clongn

    CALL ERR_HDL(NF90_OPEN(cdfile,NF90_NOWRITE,ncid))
    istatus0 = NF90_INQ_VARID ( ncid,cdvar,id_var)
    DO WHILE  ( istatus0 == NF90_ENOTVAR ) 
       ivar=getnvar(cdfile)
       PRINT *, 'Give the number corresponding to the variable you want to work with '
       DO jvar = 1, ivar
          clongn=''
          istatus=NF90_INQUIRE_VARIABLE (ncid, jvar, cdvar,ndims=idi)
          istatus=NF90_GET_ATT (ncid,jvar,clongname,clongn)
          IF (istatus /= NF90_NOERR ) clongn='unknown'
          PRINT *, jvar, ' ',TRIM(cdvar),' ',TRIM(clongn)
       ENDDO
       READ *,id_var
       istatus0=NF90_INQUIRE_VARIABLE (ncid, id_var, cdvar,ndims=idi)
    ENDDO
    ! 
    CALL ERR_HDL(NF90_INQUIRE_VARIABLE (ncid, id_var, cdvar,ndims=idi))
    getvdim=idi-1
    CALL ERR_HDL (NF90_CLOSE(ncid))
  END FUNCTION getvdim

  FUNCTION  getnvar (cdfile)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  getnvar  ***
    !!
    !! ** Purpose : return the number of variables in cdfile
    !!
    !! ** Method  :
    !!
    !! history:
    !!    27/04/2005 : Jean-Marc Molines : Original Code
    !!-----------------------------------------------------------
    !! * Arguments declarations
    CHARACTER(LEN=*), INTENT(in) ::  cdfile   ! file to look at
    INTEGER :: getnvar                        ! return the number of variables

    !! * Local variables
    INTEGER :: ncid
    INTEGER :: istatus

    istatus=NF90_OPEN(cdfile,NF90_NOWRITE,ncid)
    istatus=NF90_INQUIRE(ncid,nvariables= getnvar)
    istatus=NF90_CLOSE(ncid)

  END FUNCTION getnvar

  FUNCTION  getipk (cdfile,knvars,cdep)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  getipk  ***
    !!
    !! ** Purpose : return the number of levels for all the variables
    !!              in cdfile. Return 0 if the variable in a vector.
    !!
    !! ** Method  : returns npk when 4D variables ( x,y,z,t )
    !!              returns  1  when 3D variables ( x,y,  t )
    !!              returns  0  when other ( vectors )
    !!
    !! history:
    !!    27/04/2005 : Jean-Marc Molines : Original Code
    !!-----------------------------------------------------------
    !! * Arguments declarations
    CHARACTER(LEN=*), INTENT(in) :: cdfile   ! File to look at
    INTEGER, INTENT(in)  ::  knvars          ! Number of variables in cdfile
    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: cdep ! optional depth dim name
    INTEGER, DIMENSION(knvars) :: getipk     ! array (variables ) of levels

    !! * local declarations
    INTEGER :: ncid, ipk, jv, iipk
    INTEGER :: istatus
    CHARACTER(LEN=80) :: cldep='dep'


    istatus=NF90_OPEN(cdfile,NF90_NOWRITE,ncid)
    IF (  PRESENT (cdep) ) cldep = cdep
    ! Note the very important TRIM below : if not, getdim crashes as it never find the correct dim !
    iipk = getdim(cdfile,TRIM(cldep),kstatus=istatus)
    IF ( istatus /= 0 ) THEN
       PRINT *,' getipk : vertical dim not found ...assume 1'
       iipk=1
    ENDIF
    DO jv = 1, knvars
       istatus=NF90_INQUIRE_VARIABLE(ncid,jv, ndims=ipk)
       IF (ipk == 4 ) THEN
          getipk(jv) = iipk
       ELSE IF (ipk == 3 ) THEN
          getipk(jv) = 1
       ELSE
          getipk(jv) = 0
       ENDIF
    END DO
    istatus=NF90_CLOSE(ncid)

  END FUNCTION getipk

  FUNCTION  getvarname (cdfile, knvars, ptypvar)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  getvarname  ***
    !!
    !! ** Purpose : return a character array with the knvars variable
    !!              name corresponding to cdfile
    !!
    !! ** Method  : 
    !!
    !!
    !! history:
    !!    27/04/2005 : Jean-Marc Molines : Original Code
    !!-----------------------------------------------------------
    !! * Arguments declarations
    CHARACTER(LEN=*), INTENT(in) :: cdfile
    INTEGER, INTENT(in)  ::  knvars                  ! Number of variables in cdfile
    CHARACTER(LEN=80), DIMENSION(knvars) :: getvarname
    TYPE (variable), DIMENSION (knvars) :: ptypvar  ! Retrieve variables attribute

    !! * local declarations
    INTEGER :: ncid,  jv, ILEN
    INTEGER :: istatus
    CHARACTER(LEN=80) :: cldum=''
    REAL(KIND=4) :: zatt

    istatus=NF90_OPEN(cdfile,NF90_NOWRITE,ncid)
    DO jv = 1, knvars
       istatus=NF90_INQUIRE_VARIABLE(ncid,jv,name=getvarname(jv) )
       ptypvar(jv)%name=getvarname(jv)
       ! look for standard attibutes
       IF ( NF90_INQUIRE_ATTRIBUTE(ncid,jv,'units',len=ILEN) == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(ncid,jv,'units',cldum(1:ILEN))
          ptypvar(jv)%units=TRIM(cldum)
          cldum =''
       ELSE 
          ptypvar(jv)%units='N/A'
       ENDIF

       IF ( NF90_INQUIRE_ATTRIBUTE(ncid,jv,'missing_value') == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(ncid,jv,'missing_value',zatt)
          ptypvar(jv)%missing_value=zatt
       ELSE 
          ptypvar(jv)%missing_value=0.
       ENDIF

       IF ( NF90_INQUIRE_ATTRIBUTE(ncid,jv,'valid_min') == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(ncid,jv,'valid_min',zatt)
          ptypvar(jv)%valid_min=zatt
       ELSE
          ptypvar(jv)%valid_min=0.
       ENDIF

       IF ( NF90_INQUIRE_ATTRIBUTE(ncid,jv,'valid_max') == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(ncid,jv,'valid_max',zatt)
          ptypvar(jv)%valid_max=zatt
       ELSE
          ptypvar(jv)%valid_max=0.
       ENDIF

       IF ( NF90_INQUIRE_ATTRIBUTE(ncid,jv,'long_name',len=ILEN) == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(ncid,jv,'long_name',cldum(1:ILEN))
          ptypvar(jv)%long_name=TRIM(cldum)
          cldum=''
       ELSE
          ptypvar(jv)%long_name='N/A'
       ENDIF

       IF ( NF90_INQUIRE_ATTRIBUTE(ncid,jv,'short_name',len=ILEN) == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(ncid,jv,'short_name',cldum(1:ILEN))
          ptypvar(jv)%short_name=TRIM(cldum)
          cldum=''
       ELSE
          ptypvar(jv)%short_name='N/A'
       ENDIF

       IF ( NF90_INQUIRE_ATTRIBUTE(ncid,jv,'online_operation',len=ILEN) == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(ncid,jv,'online_operation',cldum(1:ILEN))
          ptypvar(jv)%online_operation=TRIM(cldum)
          cldum=''
       ELSE
          ptypvar(jv)%online_operation='N/A'
       ENDIF

       IF ( NF90_INQUIRE_ATTRIBUTE(ncid,jv,'axis',len=ILEN) == NF90_NOERR ) THEN
          istatus=NF90_GET_ATT(ncid,jv,'axis',cldum(1:ILEN))
          ptypvar(jv)%axis=TRIM(cldum)
          cldum=''
       ELSE
          ptypvar(jv)%axis='N/A'
       ENDIF

    END DO
    istatus=NF90_CLOSE(ncid)

  END FUNCTION getvarname

  FUNCTION  getvar (cdfile,cdvar,klev,kpi,kpj,kimin,kjmin, ktime)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  getvar  ***
    !!
    !! ** Purpose : Return the 2D REAL variable cvar, from cdfile at level klev.
    !!              kpi,kpj are the horizontal size of the 2D variable
    !!
    !! ** Method  :
    !!
    !!
    !! history:
    !!    27/04/2005 : Jean-Marc Molines : Original Code
    !!-----------------------------------------------------------
    !! * Arguments declarations
    CHARACTER(LEN=*), INTENT(in) :: cdfile,     &   ! file name to work with
         &                          cdvar           ! variable name to work with
    INTEGER, INTENT(in) :: kpi,kpj                  ! horizontal size of the 2D variable
    INTEGER, OPTIONAL, INTENT(in) :: klev           ! Optional variable. If missing 1 is assumed
    INTEGER, OPTIONAL, INTENT(in) :: kimin,kjmin    ! Optional variable. If missing 1 is assumed
    INTEGER, OPTIONAL, INTENT(in) :: ktime          ! Optional variable. If missing 1 is assumed
    REAL(KIND=4), DIMENSION(kpi,kpj) :: getvar      ! 2D REAL 4 holding variable field at klev

    !! * Local variables
    INTEGER, DIMENSION(4) :: istart, icount
    INTEGER :: ncid, id_var
    INTEGER :: istatus, ilev, imin, jmin, itime, ilog

    LOGICAL :: llog=.FALSE. , lsf=.FALSE. , lao=.FALSE.
    REAL(KIND=4) :: sf=1., ao=0.        !: Scale factor and add_offset
    REAL(KIND=4) :: spval               !: missing value

    IF (PRESENT(klev) ) THEN
       ilev=klev
    ELSE
       ilev=1
    ENDIF

    IF (PRESENT(kimin) ) THEN
       imin=kimin
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

    ! Must reset the flags to false for every call to getvar
    llog=.FALSE.
    lsf=.FALSE.
    lao=.FALSE.

    istart(1) = imin
    istart(2) = jmin
    ! JMM ! it workd for X Y Z T file,   not for X Y T .... try to found a fix !
    istart(3) = ilev
    istart(4) = itime

    icount(1)=kpi
    icount(2)=kpj
    icount(3)=1
    icount(4)=1

    CALL ERR_HDL(NF90_OPEN(cdfile,NF90_NOWRITE,ncid) )
    CALL ERR_HDL(NF90_INQ_VARID ( ncid,cdvar,id_var))

    istatus=NF90_INQUIRE_ATTRIBUTE(ncid,id_var,'missing_value')
    IF (istatus == NF90_NOERR ) THEN
       istatus=NF90_GET_ATT(ncid,id_var,'missing_value',spval)
    ELSE
       ! assume spval is 0 ?
       spval = 0.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(ncid,id_var,'savelog10')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(ncid,id_var,'savelog10',ilog)
       IF ( ilog /= 0 ) llog=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(ncid,id_var,'scale_factor')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(ncid,id_var,'scale_factor',sf)
       IF ( sf /= 1. ) lsf=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(ncid,id_var,'add_offset')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(ncid,id_var,'add_offset',ao)
       IF ( ao /= 0.) lao=.TRUE.
    ENDIF


    istatus=NF90_GET_VAR(ncid,id_var,getvar, start=istart,count=icount)
    IF ( istatus /= 0 ) THEN
       PRINT *,' Problem in getvar for ', TRIM(cdvar)
       CALL ERR_HDL(istatus)
       STOP
    ENDIF

    ! Caution : order does matter !
    IF (lsf )  WHERE (getvar /= spval )  getvar=getvar*sf
    IF (lao )  WHERE (getvar /= spval )  getvar=getvar + ao
    IF (llog)  WHERE (getvar /= spval )  getvar=10**getvar

    istatus=NF90_CLOSE(ncid)

  END FUNCTION getvar

  FUNCTION  getvarxz (cdfile,cdvar,kj,kpi,kpz,kimin,kkmin,ktime)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  getvar  ***
    !!
    !! ** Purpose : Return the 2D REAL variable x-z slab cvar, from cdfile at j=kj
    !!              kpi,kpz are the  size of the 2D variable
    !!
    !! ** Method  :
    !!
    !! history:
    !!    03/03/2006 : Jean-Marc Molines : Original code
    !!-----------------------------------------------------------
    !! * Arguments declarations
    CHARACTER(LEN=*), INTENT(in) :: cdfile,     &   ! file name to work with
         &                          cdvar           ! variable name to work with
    INTEGER, INTENT(in) :: kpi,kpz                  ! size of the 2D variable
    INTEGER, INTENT(in) :: kj                       ! Optional variable. If missing 1 is assumed
    INTEGER, OPTIONAL, INTENT(in) :: kimin,kkmin    ! Optional variable. If missing 1 is assumed
    INTEGER, OPTIONAL, INTENT(in) :: ktime          ! Optional variable. If missing 1 is assumed 
    REAL(KIND=4), DIMENSION(kpi,kpz) :: getvarxz    ! 2D REAL 4 holding variable x-z slab at kj

    !! * Local variables
    INTEGER, DIMENSION(4) :: istart, icount
    INTEGER :: ncid, id_var
    INTEGER :: istatus, ilev, imin, kmin,  itime, ilog

    LOGICAL :: llog=.FALSE. , lsf=.FALSE. , lao=.FALSE.
    REAL(KIND=4) :: sf=1., ao=0.        !: Scale factor and add_offset
    REAL(KIND=4)         :: spval       !: Missing values


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

    istart(1) = imin
    istart(2) = kj
    istart(3) = kmin
    ! JMM ! it workd for X Y Z T file,   not for X Y T .... try to found a fix !
    istart(4) = itime

    icount(1)=kpi
    icount(2)=1
    icount(3)=kpz
    icount(4)=1

    CALL ERR_HDL(NF90_OPEN(cdfile,NF90_NOWRITE,ncid) )
    CALL ERR_HDL(NF90_INQ_VARID ( ncid,cdvar,id_var))

    istatus=NF90_INQUIRE_ATTRIBUTE(ncid,id_var,'missing_value')
    IF (istatus == NF90_NOERR ) THEN
       istatus=NF90_GET_ATT(ncid,id_var,'missing_value',spval)
    ELSE
       ! assume spval is 0 ?
       spval = 0.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(ncid,id_var,'savelog10')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(ncid,id_var,'savelog10',ilog)
       IF ( ilog /= 0 ) llog=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(ncid,id_var,'scale_factor')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(ncid,id_var,'scale_factor',sf)
       IF ( sf /= 1. ) lsf=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(ncid,id_var,'add_offset')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(ncid,id_var,'add_offset',ao)
       IF ( ao /= 0.) lao=.TRUE.
    ENDIF

    istatus=NF90_GET_VAR(ncid,id_var,getvarxz, start=istart,count=icount)
    IF ( istatus /= 0 ) THEN
       PRINT *,' Problem in getvarxz for ', TRIM(cdvar)
       CALL ERR_HDL(istatus)
       STOP
    ENDIF

    ! Caution : order does matter !
    IF (lsf )  WHERE (getvarxz /= spval )  getvarxz=getvarxz*sf
    IF (lao )  WHERE (getvarxz /= spval )  getvarxz=getvarxz + ao
    IF (llog)  WHERE (getvarxz /= spval )  getvarxz=10**getvarxz

    istatus=NF90_CLOSE(ncid)

  END FUNCTION getvarxz

  FUNCTION  getvaryz (cdfile,cdvar,ki,kpj,kpz,kjmin,kkmin,ktime)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  getvar  ***
    !!
    !! ** Purpose : Return the 2D REAL variable y-z slab cvar, from cdfile at i=ki
    !!              kpj,kpz are the  size of the 2D variable
    !!
    !! ** Method  :
    !!
    !! history:
    !!    03/03/2006 : Jean-Marc Molines : Original code
    !!-----------------------------------------------------------
    !! * Arguments declarations
    CHARACTER(LEN=*), INTENT(in) :: cdfile,     &   ! file name to work with
         &                          cdvar           ! variable name to work with
    INTEGER, INTENT(in) :: kpj,kpz                  ! size of the 2D variable
    INTEGER, INTENT(in) :: ki                       ! 
    INTEGER, OPTIONAL, INTENT(in) :: kjmin,kkmin    ! Optional variable. If missing 1 is assumed
    INTEGER, OPTIONAL, INTENT(in) :: ktime          ! Optional variable. If missing 1 is assumed
    REAL(KIND=4), DIMENSION(kpj,kpz) :: getvaryz    ! 2D REAL 4 holding variable x-z slab at kj

    !! * Local variables
    INTEGER, DIMENSION(4) :: istart, icount
    INTEGER :: ncid, id_var
    INTEGER :: istatus, ilev, jmin, kmin, itime, ilog

    LOGICAL :: llog=.FALSE. , lsf=.FALSE. , lao=.FALSE.
    REAL(KIND=4) :: sf=1., ao=0.        !: Scale factor and add_offset
    REAL(KIND=4)         :: spval       !: Missing values

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

    istart(1) = ki
    istart(2) = jmin
    istart(3) = kmin
    istart(4) = 1

    icount(1)=1
    icount(2)=kpj
    icount(3)=kpz
    ! JMM ! it workd for X Y Z T file,   not for X Y T .... try to found a fix !
    icount(4)=itime

    CALL ERR_HDL(NF90_OPEN(cdfile,NF90_NOWRITE,ncid) )
    CALL ERR_HDL(NF90_INQ_VARID ( ncid,cdvar,id_var))

    istatus=NF90_INQUIRE_ATTRIBUTE(ncid,id_var,'missing_value')
    IF (istatus == NF90_NOERR ) THEN
       istatus=NF90_GET_ATT(ncid,id_var,'missing_value',spval)
    ELSE
       ! assume spval is 0 ?
       spval = 0.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(ncid,id_var,'savelog10')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(ncid,id_var,'savelog10',ilog)
       IF ( ilog /= 0 ) llog=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(ncid,id_var,'scale_factor')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(ncid,id_var,'scale_factor',sf)
       IF ( sf /= 1. ) lsf=.TRUE.
    ENDIF

    istatus=NF90_INQUIRE_ATTRIBUTE(ncid,id_var,'add_offset')
    IF (istatus == NF90_NOERR ) THEN
       ! there is a scale factor for this variable
       istatus=NF90_GET_ATT(ncid,id_var,'add_offset',ao)
       IF ( ao /= 0.) lao=.TRUE.
    ENDIF

    istatus=NF90_GET_VAR(ncid,id_var,getvaryz, start=istart,count=icount)
    IF ( istatus /= 0 ) THEN
       PRINT *,' Problem in getvaryz for ', TRIM(cdvar)
       CALL ERR_HDL(istatus)
       STOP
    ENDIF

    ! Caution : order does matter !
    IF (lsf )  WHERE (getvaryz /= spval )  getvaryz=getvaryz*sf
    IF (lao )  WHERE (getvaryz /= spval )  getvaryz=getvaryz + ao
    IF (llog)  WHERE (getvaryz /= spval )  getvaryz=10**getvaryz

    istatus=NF90_CLOSE(ncid)

  END FUNCTION getvaryz

  FUNCTION  getvar1d (cdfile,cdvar,kk,kstatus)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  getvar1d  ***
    !!
    !! ** Purpose :  return 1D variable cdvar from cdfile, of size kk
    !!
    !! ** Method  :
    !!
    !!
    !! history:
    !!    27/04/2005 : Jean-Marc Molines : Original Code
    !!-----------------------------------------------------------
    !! * Arguments declarations
    CHARACTER(LEN=*), INTENT(in) :: cdfile,     &   ! file name to work with
         &                          cdvar           ! variable name to work with
    INTEGER, INTENT(in) :: kk                       ! size of 1D vector to be returned
    INTEGER, OPTIONAL, INTENT(out) :: kstatus       ! return status concerning the variable existence
    REAL(KIND=4), DIMENSION(kk) :: getvar1d         ! real returned vector

    !! * Local variables
    INTEGER, DIMENSION(1) :: istart, icount
    INTEGER :: ncid, id_var
    INTEGER :: istatus

    istart(:) = 1
    icount(1)=kk
    IF ( PRESENT(kstatus) ) kstatus = 0

    istatus=NF90_OPEN(cdfile,NF90_NOWRITE,ncid)
    istatus=NF90_INQ_VARID ( ncid,cdvar,id_var)
    IF ( istatus == NF90_NOERR ) THEN
       istatus=NF90_GET_VAR(ncid,id_var,getvar1d,start=istart,count=icount)
    ELSE
       IF ( PRESENT(kstatus) ) kstatus= istatus
       getvar1d=99999999999.
    ENDIF

    istatus=NF90_CLOSE(ncid)

  END FUNCTION getvar1d

  FUNCTION  getvare3 (cdfile,cdvar,kk)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  getvare3  ***
    !!
    !! ** Purpose :  Special routine for e3, which in fact is a 1D variable
    !!               but defined as e3 (1,1,npk,1) in coordinates.nc (!!)
    !!
    !! ** Method  :
    !!
    !!
    !! history:
    !!    27/04/2005 : Jean-Marc Molines : Original Code
    !!-----------------------------------------------------------
    !! * Arguments declarations
    CHARACTER(LEN=*), INTENT(in) :: cdfile,     &   ! file name to work with
         &                          cdvar           ! variable name to work with
    INTEGER, INTENT(in) :: kk                       ! size of 1D vector to be returned
    REAL(KIND=4), DIMENSION(kk) :: getvare3         ! return e3 variable form the coordinate file

    !! * Local variables   
    INTEGER, DIMENSION(4) :: istart, icount
    INTEGER :: ncid, id_var
    INTEGER :: istatus

    istart(:) = 1
    icount(:) = 1
    icount(3)=kk

    istatus=NF90_OPEN(cdfile,NF90_NOWRITE,ncid)
    istatus=NF90_INQ_VARID ( ncid,cdvar,id_var)
    istatus=NF90_GET_VAR(ncid,id_var,getvare3,start=istart,count=icount)
    IF ( istatus /= 0 ) THEN
       PRINT *,' Problem in getvare3 for ', TRIM(cdvar)
       PRINT *,TRIM(cdfile), kk
       CALL ERR_HDL(istatus)
       STOP
    ENDIF

    istatus=NF90_CLOSE(ncid)
  END FUNCTION getvare3


  FUNCTION putheadervar(kout, cdfile, kpi,kpj,kpk, pnavlon, pnavlat ,pdep,cdep)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  putheadervar  ***
    !!
    !! ** Purpose :  copy header variables from cdfile to the already open ncfile (ncid=kout)
    !!
    !! ** Method  :  header variables are nav_lat, nav_lon and either (deptht, depthu, or depthv )
    !!               Even if the use of different variable name for deptht, depthu depthv is 
    !!               one of the many non sense of IOIPSL, we are forced to stick with !
    !!               (Note that these 3 depth are identical in OPA. On the other hand, nav_lon, nav_lat
    !!                differ for U and V and T points but have the same variable name). 
    !!               If pnavlon and pnavlat are provided as arguments, they are used for nav_lon, nav_lat
    !!               instead of the nav_lon,nav_lat read on the file cdfile.
    !! 
    !! ** Action  : header variables for file kout is copied from cdfile
    !!
    !! history:
    !!    27/04/2005 : Jean-Marc Molines : Original Code
    !!-----------------------------------------------------------
    !! * Arguments declarations
    INTEGER, INTENT(in) :: kout             ! ncid of the outputfile (already open )
    CHARACTER(LEN=*), INTENT(in) :: cdfile  ! file from where the headers will be copied
    INTEGER, INTENT(in) :: kpi,kpj,kpk      ! dimension of nav_lon,nav_lat (kpi,kpj), and depht(kpk)
    REAL(KIND=4), OPTIONAL, DIMENSION(kpi,kpj), INTENT(in) :: pnavlon, pnavlat  ! array provided optionaly to overrid the
    !                                   !   corresponding arrays in cdfile
    REAL(KIND=4), OPTIONAL,DIMENSION(kpk), INTENT(in) :: pdep   ! dep array if not on cdfile
    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: cdep     ! optional name of vertical variable
    INTEGER :: putheadervar                 ! return status

    !! * Local variables
    INTEGER , PARAMETER :: jpdep=5
    INTEGER :: istatus, idep, jj
    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: z2d
    REAL(KIND=4), DIMENSION(kpk) :: z1d
    CHARACTER(LEN=80),DIMENSION(jpdep ) :: cldept=(/'deptht ','depthu ','depthv ','depthw ','nav_lev'/)
    CHARACTER(LEN=80) :: cldep

    ALLOCATE ( z2d (kpi,kpj) )
    IF (PRESENT(pnavlon) ) THEN 
       z2d = pnavlon
    ELSE
       z2d=getvar(cdfile,'nav_lon', 1,kpi,kpj)
    ENDIF
    istatus = putvar(kout,id_lon,z2d,1,kpi,kpj)

    IF (PRESENT(pnavlat) ) THEN
       z2d = pnavlat
    ELSE
       z2d=getvar(cdfile,'nav_lat', 1,kpi,kpj)
    ENDIF

    istatus = putvar(kout,id_lat,z2d,1,kpi,kpj)

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
                STOP
             ENDIF
          ENDIF
       ENDIF

       istatus = putvar1d(kout,z1d,kpk,'D')
    ENDIF
    putheadervar=istatus
    DEALLOCATE (z2d)

  END FUNCTION putheadervar

  FUNCTION putvarr4(kout, kid,ptab, klev, kpi, kpj,ktime)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  putvar  ***
    !!
    !! ** Purpose :  copy a 2D level of ptab in already open file kout, using variable kid
    !!
    !! ** Method  :  
    !!
    !! ** Action  : variable level written
    !!
    !! history:
    !!    27/04/2005 : Jean-Marc Molines : Original Code
    !!-----------------------------------------------------------
    !! * Arguments declarations
    INTEGER, INTENT(in) :: kout  ,  &       ! ncid of output file
         &                  kid              ! varid of output variable
    INTEGER, INTENT(in) :: klev             ! level at which ptab will be written
    INTEGER, INTENT(in) :: kpi,kpj          ! dimension of ptab
    INTEGER, OPTIONAL, INTENT(in) :: ktime  ! dimension of ptab
    REAL(KIND=4), DIMENSION(kpi,kpj),INTENT(in) :: ptab ! 2D array to write in file
    INTEGER :: putvarr4                       ! return status

    !! * Local variables    
    INTEGER :: istatus, itime
    INTEGER, DIMENSION(4) :: istart, icount

    IF (PRESENT(ktime) ) THEN
       itime=ktime
    ELSE
       itime=1
    ENDIF

    istart(:) = 1 ; istart(3)=klev ; istart(4)=itime
    icount(:) = 1 ; icount(1) = kpi ; icount(2) = kpj
    istatus=NF90_PUT_VAR(kout,kid, ptab, start=istart,count=icount)
    putvarr4=istatus

  END FUNCTION putvarr4

  FUNCTION reputvarr4 (cdfile,cdvar,klev,kpi,kpj,kimin,kjmin, ktime,ptab)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  putvar  ***
    !!
    !! ** Purpose :  Change an existing variable in inputfile
    !!
    !! ** Method  :  
    !!
    !! ** Action  : variable level written
    !!
    !! history:
    !!    27/04/2005 : Jean-Marc Molines : Original Code
    !!-----------------------------------------------------------
    !! * Arguments declarations

    CHARACTER(LEN=*), INTENT(in) :: cdfile,     &   ! file name to work with
         &                          cdvar           ! variable name to work with
    INTEGER, INTENT(in) :: kpi,kpj                  ! horizontal size of the 2D variable
    INTEGER, OPTIONAL, INTENT(in) :: klev           ! Optional variable. If missing 1 is assumed
    INTEGER, OPTIONAL, INTENT(in) :: kimin,kjmin    ! Optional variable. If missing 1 is assumed
    INTEGER, OPTIONAL, INTENT(in) :: ktime          ! Optional variable. If missing 1 is assumed
    REAL(KIND=4), DIMENSION(kpi,kpj) ::  ptab     ! 2D REAL 4 holding variable field at klev
    INTEGER :: reputvarr4

    !! * Local variables
    INTEGER, DIMENSION(4) :: istart, icount
    INTEGER :: ncid, id_var
    INTEGER :: istatus, ilev, imin, jmin, itime

    ilev=1  ; IF (PRESENT(klev)  ) ilev=klev
    imin=1  ; IF (PRESENT(kimin) ) imin=kimin
    jmin=1  ; IF (PRESENT(kjmin) ) jmin=kjmin
    itime=1 ; IF (PRESENT(ktime) ) itime=ktime

    istatus=NF90_OPEN(cdfile,NF90_WRITE,ncid)
    istatus=NF90_INQ_VARID(ncid,cdvar,id_var)
    istatus=NF90_PUT_VAR(ncid,id_var,ptab,start=(/imin,jmin,klev,itime/), count=(/kpi,kpj,1,1/) )
    !PRINT *,TRIM(NF90_STRERROR(istatus)),' in reputvar'
    reputvarr4=istatus
    istatus=NF90_CLOSE(ncid)

  END FUNCTION reputvarr4

  FUNCTION putvarzo(kout, kid,ptab, klev, kpi, kpj,ktime)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  putvarzo  ***
    !!
    !! ** Purpose :  copy a 2D level of ptab in already open file kout, using variable kid
    !!               This variant deals with degenerated 2D (1 x jpj) zonal files
    !!
    !! ** Method  :
    !!
    !! ** Action  : variable level written
    !!
    !! history:
    !!    27/04/2005 : Jean-Marc Molines : Original Code
    !!-----------------------------------------------------------
    !! * Arguments declarations
    INTEGER, INTENT(in) :: kout  ,  &       ! ncid of output file
         &                  kid              ! varid of output variable
    INTEGER, INTENT(in) :: klev             ! level at which ptab will be written
    INTEGER, INTENT(in) :: kpi,kpj          ! dimension of ptab
    INTEGER, OPTIONAL, INTENT(in) :: ktime  ! dimension of ptab
    REAL(KIND=4), DIMENSION(kpj),INTENT(in) :: ptab ! 2D array to write in file
    INTEGER :: putvarzo                       ! return status

    !! * Local variables
    INTEGER :: istatus, itime
    INTEGER, DIMENSION(4) :: istart, icount

    IF (PRESENT(ktime) ) THEN
       itime=ktime
    ELSE
       itime=1
    ENDIF

    istart(:) = 1 ; istart(3)=klev ; istart(4)=itime
    icount(:) = 1 ; icount(1) = kpi ; icount(2) = kpj
    istatus=NF90_PUT_VAR(kout,kid, ptab, start=istart,count=icount)
    putvarzo=istatus

  END FUNCTION putvarzo


  FUNCTION putvari2(kout, kid,ptab, klev, kpi, kpj,ktime)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  putvar  ***
    !!
    !! ** Purpose :  copy a 2D level of ptab in already open file kout, using variable kid
    !!
    !! ** Method  :
    !!
    !! ** Action  : variable level written
    !!
    !! history:
    !!    27/04/2005 : Jean-Marc Molines : Original Code
    !!-----------------------------------------------------------
    !! * Arguments declarations
    INTEGER, INTENT(in) :: kout  ,  &       ! ncid of output file
         &                  kid              ! varid of output variable
    INTEGER, INTENT(in) :: klev             ! level at which ptab will be written
    INTEGER, INTENT(in) :: kpi,kpj          ! dimension of ptab
    INTEGER, OPTIONAL, INTENT(in) :: ktime  ! dimension of ptab
    INTEGER(KIND=2), DIMENSION(kpi,kpj),INTENT(in) :: ptab ! 2D array to write in file
    INTEGER :: putvari2                       ! return status

    !! * Local variables
    INTEGER :: istatus, itime
    INTEGER, DIMENSION(4) :: istart, icount

    IF (PRESENT(ktime) ) THEN
       itime=ktime
    ELSE
       itime=1
    ENDIF

    istart(:) = 1 ; istart(3)=klev ; istart(4)=itime
    icount(:) = 1 ; icount(1) = kpi ; icount(2) = kpj
    istatus=NF90_PUT_VAR(kout,kid, ptab, start=istart,count=icount)
    putvari2=istatus

  END FUNCTION putvari2


  FUNCTION putvar1d(kout,ptab,kk,cdtype)
    !!-----------------------------------------------------------
    !!                       ***  FUNCTION  putvar1d  ***
    !!
    !! ** Purpose : Copy 1D variable (size kk) hold in ptab,  with id kid, into file id kout
    !!
    !! ** Method  :  cdtype is either T (time_counter) or D (depth?)
    !!
    !! ** Action  : 1D variable  written
    !!
    !! history:
    !!    27/04/2005 : Jean-Marc Molines : Original Code
    !!-----------------------------------------------------------
    !! * Arguments declarations
    INTEGER, INTENT(in) :: kout             ! ncid of output file
    INTEGER, INTENT(in) :: kk               ! number of elements in ptab
    REAL(KIND=4), DIMENSION(kk),INTENT(in) :: ptab ! 1D array to write in file
    CHARACTER(LEN=1), INTENT(in)  :: cdtype ! either T or D
    INTEGER :: putvar1d                     ! return status

    !! * Local variables     
    INTEGER :: istatus, iid
    INTEGER, DIMENSION(1) :: istart, icount

    SELECT CASE ( cdtype )
    CASE ('T', 't' ) 
       iid = id_tim
    CASE ('D', 'd' )
       iid = id_dep
    END SELECT

    istart(:) = 1
    icount(:) = kk
    istatus=NF90_PUT_VAR(kout,iid, ptab, start=istart,count=icount)
    putvar1d=istatus

  END FUNCTION putvar1d

  FUNCTION closeout(kout)
    !!----------------------------------------------------------
    !!                       ***  FUNCTION  closeout  ***
    !!
    !! ** Purpose : close open output files
    !!
    !! history:
    !!    27/04/2005 : Jean-Marc Molines : Original Code
    !!-----------------------------------------------------------
    INTEGER,INTENT(in) :: kout   ! ncid of file to be closed
    INTEGER :: closeout          ! return status
    closeout=NF90_CLOSE(kout)
  END FUNCTION closeout

  FUNCTION ncopen(cdfile)
    !!----------------------------------------------------------
    !!                       ***  FUNCTION  ncopen  ***
    !!
    !! ** Purpose : open file cdfile and return file ID
    !!
    !!-----------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: cdfile ! file name
      INTEGER :: ncopen                      ! return status
    ! * Local variables
      INTEGER :: istatus, ncid
      istatus = NF90_OPEN(cdfile,NF90_WRITE,ncid)
      ncopen=ncid
  END FUNCTION ncopen

  SUBROUTINE ERR_HDL(kstatus)
    !! ----------------------------------------------------------
    !!   ***  SUBROUTINE err_hdl
    !!
    !!   ** Purpose :  Error handle for NetCDF routine.
    !!          Stop if kstatus indicates error conditions.
    !!
    !! History :
    !!     Original: J.M. Molines (01/99)
    !!
    !! -----------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(in) ::  kstatus
    IF (kstatus /=  NF90_NOERR ) THEN
       PRINT *, 'ERROR in NETCDF routine, status=',kstatus
       PRINT *,NF90_STRERROR(kstatus)
       STOP
    END IF
  END SUBROUTINE ERR_HDL

  SUBROUTINE gettimeseries (cdfile, cdvar, kilook, kjlook,klev)
    !! ----------------------------------------------------------
    !!   ***  SUBROUTINE gettimeseries  ***
    !!
    !!   ** Purpose : Display a 2 column output ( time, variable) for
    !!           a given variable of a given file at a given point
    !!
    !! History :
    !!     Original: J.M. Molines (03/2007)
    !!
    !! -----------------------------------------------------------
    !* Arguments
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(in) :: cdfile, cdvar
    INTEGER,INTENT(in) :: kilook,kjlook
    INTEGER, OPTIONAL, INTENT(in) :: klev
    !* Local variables
    INTEGER ::  nt, jt
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: ztime, zval
    REAL(KIND=4) :: ztmp, zao=0., zsf=1.0   !: add_offset, scale_factor

    ! netcdf stuff
    INTEGER :: istatus
    INTEGER :: ncid, id_t, id_var, ndim, jk

    ! Klev can be used to give the model level we want to look at
    IF ( PRESENT(klev) ) THEN
       jk=klev
    ELSE
       jk=1
    ENDIF

    ! Open cdf dataset
    istatus=NF90_OPEN(cdfile,NF90_NOWRITE,ncid)

    ! read time dimension
    istatus=NF90_INQ_DIMID(ncid,'time_counter',id_t)
    istatus=NF90_INQUIRE_DIMENSION(ncid,id_t,len=nt)

    ! Allocate space
    ALLOCATE (ztime(nt), zval(nt) )

    ! gettime
    istatus=NF90_INQ_VARID(ncid,'time_counter',id_var)
    istatus=NF90_GET_VAR(ncid,id_var,ztime,(/1/),(/nt/) )

    ! read variable
    istatus=NF90_INQ_VARID(ncid,cdvar,id_var)
    !   look for scale_factor and add_offset attribute:
    istatus=NF90_GET_ATT(ncid,id_var,'add_offset',ztmp)
    IF ( istatus == NF90_NOERR ) zao = ztmp
    istatus=NF90_GET_ATT(ncid,id_var,'scale_factor',ztmp)
    IF ( istatus == NF90_NOERR ) zsf = ztmp

    ! get number of dimension of the variable ( either x,y,t or x,y,z,t )
    istatus=NF90_INQUIRE_VARIABLE(ncid,id_var,ndims=ndim)
    IF ( ndim == 3 ) THEN
       istatus=NF90_GET_VAR(ncid,id_var,zval,(/kilook,kjlook,1/),(/1,1,nt/) )
    ELSE IF ( ndim == 4 ) THEN
       istatus=NF90_GET_VAR(ncid,id_var,zval,(/kilook,kjlook,jk,1/),(/1,1,1,nt/) )
    ELSE 
       PRINT *,'  ERROR : variable ',TRIM(cdvar),' has ', ndim, &
            &       ' dimensions !. Only 3 or 4 supported'
       STOP
    ENDIF

    ! convert to physical values
    zval=zval*zsf + zao

    ! display results :
    DO jt=1,nt
       PRINT *,ztime(jt)/86400., zval(jt)
    ENDDO

    istatus=NF90_CLOSE(ncid)

  END SUBROUTINE gettimeseries

END MODULE cdfio

