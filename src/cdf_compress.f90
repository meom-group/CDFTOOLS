PROGRAM cdf_compress
  !!======================================================================
  !!                     ***  PROGRAM  cdf_compress  ***
  !!=====================================================================
  !!  ** Purpose : a low memory consuming tool to be used in place of
  !!               ncks for building netcdf4/Hdf5 deflated files
  !!
  !!  ** Method  :  Read and write deflated
  !!
  !! History :  4.0  : 01/2020  : J.M. Molines :  (Energetics)
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------

  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class file_operations
  !!----------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE

  INTEGER(KIND=4)    :: jd, jv, j3,j4
  INTEGER(KIND=4)    :: n1 , n2, n3, n4
  INTEGER(KIND=4)    :: ncin, ncou
  INTEGER(KIND=4)    :: narg, iargc, ijarg
  INTEGER(KIND=4)    :: idef_lev, nvertdim, ilen, ndim, ierr
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: icnk_dmn

  REAL(KIND=4) , DIMENSION(:,:,:), ALLOCATABLE :: r3d

  REAL(KIND=8) , DIMENSION(:), ALLOCATABLE    ::  d1d
  REAL(KIND=8) , DIMENSION(:,:), ALLOCATABLE  ::  d2d

  CHARACTER(LEN=255) :: cf_in, cf_ou, cdum
  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: c_vertdim
  
  LOGICAL :: l_3d=.FALSE.

  TYPE :: ncfile
     INTEGER(KIND=4)                              :: ncid   ! file id
     INTEGER(KIND=4)                              :: ndims  ! number of variables
     INTEGER(KIND=4)                              :: nvars   ! number of vars
     INTEGER(KIND=4)                              :: natts   ! number of global attributes
     INTEGER(KIND=4)                              :: iunlim  ! ID of unlimited dimension
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: ideflat ! deflate level (nvar)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nvatt   ! number of att of each variable (var)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nvid    ! varid of each variable (var)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nvdim   ! dimension of each variable (var)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: itype   ! type of each variable (var)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nlen    ! len of each dimension ( ndims)
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: idimid  ! dimid of each dimension ( ndims)
     INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ichunk  ! size of chunk (nvar, ndims)
     INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: idimids ! dimids of each variable (nvar, ndims) 
     CHARACTER(LEN=255)                            :: c_fnam ! name of working file
     CHARACTER(LEN=80 ), DIMENSION(:), ALLOCATABLE :: c_vnam ! name of each variable (var)
     CHARACTER(LEN=80 ), DIMENSION(:), ALLOCATABLE :: c_dnam ! name of each dimension (ndims)

  END TYPE ncfile

  TYPE(ncfile) :: scdf_in, scdf_ou
  !!----------------------------------------------------------------------
  ! default values
  idef_lev = 1
  nvertdim = 0
  cf_ou='none'

  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdf_compress -f FILE-in [-v LIST-vertical_Dimensions ] [-d DEF-lev ]'
     PRINT *,'                     [-o FILE-out] [-3D]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compress input file using netcdf4/HDF5 with deflation level <DEF_lev>.'
     PRINT *,'       Need to specify vertical dimension names if any. Note that this tool'
     PRINT *,'       is designed to work with any netcdf file, not only NEMO. It is a low'
     PRINT *,'       memory alternative of ncks for huge big files. Sizes of chunks are '
     PRINT *,'       predermined, but always 1 for vertical and time axis.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f FILE-in : name of netcdfile to be compressed.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -v LIST-vertical_dimensions : comma separated list of the names of'
     PRINT *,'             the vertical dimensions (in order to set chunk size of 1). '
     PRINT *,'       -d DEF-lev  : Specify deflation level ( default is 1).'
     PRINT *,'       -o FILE-out : Specify output file instead of default <FILE-in>4 '
     PRINT *,'       -3D : use 3D arrays for eventual increased performances (not sure)'
     PRINT *,'            Use more in core memory.' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : <FILE_in>4 or file specified with -o option.'
     PRINT *,'         variables : same as input ( with same attributes).'
     PRINT *,'      '
     STOP
  ENDIF

  ijarg=1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg,cdum) ; ijarg=ijarg+1
     SELECT CASE ( cdum)
     CASE ( '-f' ) ; CALL getarg(ijarg,cf_in) ; ijarg=ijarg+1
        ! options
     CASE ( '-d'  ) ; CALL getarg(ijarg,cdum ) ; ijarg=ijarg+1 ; READ(cdum,*) idef_lev
     CASE ( '-v'  ) ; CALL getarg(ijarg,cdum ) ; ijarg=ijarg+1 ; CALL getlist(cdum)
     CASE ( '-o'  ) ; CALL getarg(ijarg,cf_ou) ; ijarg=ijarg+1
     CASE ( '-3D' ) ; l_3d = .TRUE.
     CASE DEFAULT  ; PRINT *,'  E R R O R : Unknown option :', TRIM(cdum) ; STOP 99
     END SELECT
  END DO

  IF ( cf_ou == 'none') THEN
     cf_ou=TRIM(cf_in)//'4'
  ENDIF

  scdf_in%c_fnam = cf_in

  ! retrieve file structure 
  CALL getncfile(scdf_in)

  ! copy input structure on output structure
  scdf_ou = scdf_in
  ! weird but nice ... this copy of the structure also allocate the target ???
  ! print *, size( scdf_ou%nvdim )
  ! print *, scdf_in%nvdim(:)
  ! print *, scdf_ou%nvdim(:)

  ! Change the output structure
  ! 1) file name
  scdf_ou%c_fnam = cf_ou
  ! 2) set deflation level  for each variables ( at this stage, all var with same def level ...
  scdf_ou%ideflat(:) = idef_lev
  ! 3) set the chunk size for each dimension according to some criteria
  ALLOCATE( icnk_dmn( scdf_ou%ndims ) )
  icnk_dmn(:) = 0
  DO jd = 1, scdf_ou%ndims
     DO jv = 1, nvertdim
        IF ( scdf_ou%c_dnam(jd) == c_vertdim(jv) ) THEN
           icnk_dmn(jd) = 1
        ENDIF
     ENDDO
     IF ( jd == scdf_ou%iunlim ) THEN
        icnk_dmn(jd) = 1
     ENDIF
     IF ( INDEX( scdf_ou%c_dnam(jd),'x') /= 0 ) THEN
        ilen=scdf_ou%nlen(jd)
        icnk_dmn(jd)=ilen !! 
     ENDIF
     IF ( icnk_dmn(jd) == 0 ) THEN  ! not treated above
        ilen=scdf_ou%nlen(jd)
        icnk_dmn(jd)=MIN(ilen, 200) !!  500 hard coded !! 
     ENDIF
  ENDDO

  ! Create the skeleton of the output file
  CALL makncfile(scdf_ou) 

  !  Read in/ write out ... Loop on variables.
  ncin = scdf_in%ncid
  ncou = scdf_ou%ncid

  DO jv = 1, scdf_in%nvars
     PRINT * , 'compress and copy ', TRIM(scdf_in%c_vnam(jv) )

     ndim=scdf_in%nvdim(jv)
     SELECT CASE ( ndim )
     CASE ( 1 )    ! can be t or z
        ALLOCATE ( d1d (scdf_in%nlen(scdf_in%idimids(jv,1)) ) )
        ierr = NF90_GET_VAR(ncin,jv,d1d )
        ierr = NF90_PUT_VAR(ncou,jv,d1d )
        ierr = NF90_SYNC(ncou)
        DEALLOCATE ( d1d)
     CASE ( 2 )   ! can be t bound or x,y
        n1 = scdf_in%nlen(scdf_in%idimids(jv,1))
        n2 = scdf_in%nlen(scdf_in%idimids(jv,2))
        ALLOCATE ( d2d (n1,n2)  )
        ierr = NF90_GET_VAR(ncin,jv,d2d)
        ierr = NF90_PUT_VAR(ncou,jv,d2d )
        ierr = NF90_SYNC(ncou)
        DEALLOCATE ( d2d)
     CASE ( 3 )   !  can be : x,y,t   x,y,z    loop on 3rd dim
        n1 = scdf_in%nlen(scdf_in%idimids(jv,1))
        n2 = scdf_in%nlen(scdf_in%idimids(jv,2))
        n3 = scdf_in%nlen(scdf_in%idimids(jv,3))
        ALLOCATE ( d2d (n1,n2) )
        DO j3 = 1, n3
           ierr = NF90_GET_VAR(ncin,jv,d2d, start=(/1,1,j3/),count=(/n1,n2,1/) )
           ierr = NF90_PUT_VAR(ncou,jv,d2d, start=(/1,1,j3/),count=(/n1,n2,1/) )
        ENDDO
        ierr = NF90_SYNC(ncou)
        DEALLOCATE ( d2d)
     CASE ( 4 )  !  can be   x,y,z,t or ???  loop on 3 rd and 4th dims
        n1 = scdf_in%nlen(scdf_in%idimids(jv,1))
        n2 = scdf_in%nlen(scdf_in%idimids(jv,2))
        n3 = scdf_in%nlen(scdf_in%idimids(jv,3))
        n4 = scdf_in%nlen(scdf_in%idimids(jv,4))
        IF ( l_3d ) THEN
          ALLOCATE ( r3d (n1,n2,n3) )
          DO j4 =   1, n4
             ierr = NF90_GET_VAR(ncin,jv,r3d, start=(/1,1,1,j4/),count=(/n1,n2,n3,1/) )
             ierr = NF90_PUT_VAR(ncou,jv,r3d, start=(/1,1,1,j4/),count=(/n1,n2,n3,1/) )
          ENDDO
          DEALLOCATE ( r3d)
        ELSE
          ALLOCATE ( d2d (n1,n2) )
          DO j4 =   1, n4
             DO j3 = 1, n3
                ierr = NF90_GET_VAR(ncin,jv,d2d, start=(/1,1,j3,j4/),count=(/n1,n2,1,1/) )
                ierr = NF90_PUT_VAR(ncou,jv,d2d, start=(/1,1,j3,j4/),count=(/n1,n2,1,1/) )
             ENDDO
          ENDDO
          DEALLOCATE ( d2d)
        ENDIF
        ierr = NF90_SYNC(ncou)

     CASE DEFAULT
        PRINT *, 'WARNING : found variable with ', ndim,' dimension : not yet supported : skip it !'
     END SELECT
  ENDDO

  ierr = NF90_CLOSE(ncin)
  ierr = NF90_CLOSE(ncou)

CONTAINS 
  SUBROUTINE getncfile (sd_ncfile)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE getncfile  ***
    !!
    !! ** Purpose :  Fill a ncfile structure with information taken from
    !!               the input file.
    !!
    !! ** Method  :  Use NF90 primitive functions
    !!
    !! References :  netcdf documentation
    !!----------------------------------------------------------------------
    TYPE(ncfile), INTENT(inout) :: sd_ncfile
    ! local vars
    INTEGER(KIND=4) :: jd, jv   ! loop index
    INTEGER(KIND=4) :: ierr, ncid, nvars, ndims  ! short names for some structure elements
    !!----------------------------------------------------------------------

    ierr = NF90_OPEN(sd_ncfile%c_fnam,NF90_NOWRITE,sd_ncfile%ncid)
    ncid =  sd_ncfile%ncid
    ! now fill in ncfile structure ...
    ierr = NF90_INQUIRE(ncid, nDimensions  = sd_ncfile%ndims  &
         &                , nVariables     = sd_ncfile%nvars  &
         &                , nAttributes    = sd_ncfile%natts  &
         &                , unlimitedDimid = sd_ncfile%iunlim  )
    ! look for dimensions 
    ndims = sd_ncfile%ndims
    ALLOCATE ( sd_ncfile%nlen ( ndims) , sd_ncfile%c_dnam ( ndims), sd_ncfile%idimid ( ndims) )
    DO jd = 1 , ndims
       ierr = NF90_INQUIRE_DIMENSION(ncid,jd, name=sd_ncfile%c_dnam(jd) , len=sd_ncfile%nlen (jd) ) 
       PRINT *, " Dimension ",TRIM(sd_ncfile%c_dnam(jd) )," : ", sd_ncfile%nlen (jd)
    ENDDO
    PRINT *, "   Unlimited dimension is ", TRIM(sd_ncfile%c_dnam(sd_ncfile%iunlim))
    ! look for variables
    nvars=sd_ncfile%nvars
    ALLOCATE ( sd_ncfile%ideflat(nvars),  sd_ncfile%nvatt(nvars),  sd_ncfile%nvid(nvars) )
    ALLOCATE ( sd_ncfile%nvdim(nvars),  sd_ncfile%itype(nvars),  sd_ncfile%c_vnam(nvars) )
    ALLOCATE ( sd_ncfile%ichunk(nvars,ndims),  sd_ncfile%idimids(nvars,ndims)            )
    DO  jv = 1, nvars
       sd_ncfile%nvid(jv) = jv
       ierr = NF90_INQUIRE_VARIABLE(ncid, jv, name        = sd_ncfile%c_vnam(jv)    &
            &                             , xtype         = sd_ncfile%itype(jv)     &
            &                             , ndims         = sd_ncfile%nvdim(jv)     &
            &                             , nAtts         = sd_ncfile%nvatt(jv)     &
            &                             , deflate_level = sd_ncfile%ideflat(jv)   &
            &                             , dimids        = sd_ncfile%idimids(jv,:) &
            &                             , chunksizes    = sd_ncfile%ichunk(jv,:)   )
    ENDDO

    ! leave input file open !!! 

  END SUBROUTINE getncfile

  SUBROUTINE makncfile (sd_ncfile )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE <routine>  ***
    !!
    !! ** Purpose :    Create a new file from the information stored in the
    !!               ncfile structure passed as argument.
    !!
    !! ** Method  :   use NF90 primitives
    !!
    !! References :  netcdf documentation
    !!----------------------------------------------------------------------
    TYPE(ncfile), INTENT(inout) :: sd_ncfile
    !
    INTEGER(KIND=4) :: ja   ! loop index
    INTEGER(KIND=4) :: ierr, ncid, ilen
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ichunksz
    CHARACTER(LEN=80) :: cl_anam
    !!----------------------------------------------------------------------
    ! Create output data set
    ierr = NF90_CREATE(sd_ncfile%c_fnam,NF90_NETCDF4,sd_ncfile%ncid)
    ncid = sd_ncfile%ncid

    ! Create dimensions
    DO jd = 1, sd_ncfile%ndims
       IF ( jd == sd_ncfile%iunlim ) THEN
          ilen = NF90_UNLIMITED
       ELSE        
          ilen =sd_ncfile%nlen(jd)
       ENDIF
       ierr = NF90_DEF_DIM(ncid,sd_ncfile%c_dnam(jd),ilen,sd_ncfile%idimid(jd))
    ENDDO

    ! Create variables
    DO jv = 1, sd_ncfile%nvars
       ALLOCATE ( ichunksz(sd_ncfile%nvdim(jv)) )
       DO jd = 1, sd_ncfile%nvdim(jv)
          ichunksz(jd)=icnk_dmn(sd_ncfile%idimids(jv,jd) )
       ENDDO
       ierr=NF90_DEF_VAR(ncid,    sd_ncfile%c_vnam(jv)                         &
            &                   , sd_ncfile%itype(jv)                          &
            &                   , sd_ncfile%idimids(jv,1:sd_ncfile%nvdim(jv))  &
            &                   , sd_ncfile%nvid(jv)                           &
            &                   , chunksizes    = ichunksz(:)                  &
            &                   , deflate_level = sd_ncfile%ideflat(jv)  )
       DEALLOCATE ( ichunksz )
    ENDDO

    ! Create attribute : copy from scdf_in  ( global variable )
    ncin=scdf_in%ncid
    ncou=ncid
    ! Global attributes
    DO ja = 1, scdf_in%natts
       ierr = NF90_INQ_ATTNAME(ncin,NF90_GLOBAL,ja,cl_anam)
       ierr = NF90_COPY_ATT(ncin,NF90_GLOBAL,cl_anam,ncou,NF90_GLOBAL)
    ENDDO
    ! variable attributes
    DO jv = 1, sd_ncfile%nvars
       DO ja=1,  sd_ncfile%nvatt(jv)
          ierr = NF90_INQ_ATTNAME(ncin,jv,ja,cl_anam)
          ierr = NF90_COPY_ATT(ncin,jv,cl_anam,ncou,jv)
       ENDDO
    ENDDO

    ierr = NF90_ENDDEF(ncid)
  END SUBROUTINE makncfile


  SUBROUTINE getlist (cd_list)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE getlist  ***
    !!
    !! ** Purpose :  decode  a comma separated list and put each element
    !!               in the global character array variable  c_vertdim
    !!
    !! ** Method  :  Count the number of ',' in the input list and infer
    !!               the number of elements .. etc 
    !!
    !!----------------------------------------------------------------------

    CHARACTER(LEN=*) , INTENT(in) :: cd_list
    !
    CHARACTER(LEN=80)  ::cl_list
    INTEGER(KIND=4) :: ilen, ipos, ii
    INTEGER(KIND=4) :: jc
    !!----------------------------------------------------------------------
    cl_list=TRIM(cd_list)

    ! count number of ','
    nvertdim=0  ! number of vertical dimensions
    ilen=LEN(TRIM(cd_list))
    DO jc=1,ilen
       IF ( cd_list(jc:jc) == ',' ) THEN
          nvertdim=nvertdim+1
       ENDIF
    ENDDO
    nvertdim=nvertdim + 1  ! number of vertical variables
    PRINT *, ' Vert dim = ', nvertdim
    ALLOCATE (c_vertdim(nvertdim) )
    ipos=0 ; ii=0
    IF ( nvertdim == 1 ) THEN
       c_vertdim(1) = cl_list
    ELSE
       DO jc=1, nvertdim
          ipos = INDEX(cl_list,',' )
          IF ( ipos /= 0 ) THEN
             ii = ii + 1
             c_vertdim(ii) = cl_list(1:ipos-1)
             cl_list=TRIM(cl_list(ipos+1:))
          ELSE
             c_vertdim(nvertdim)= cl_list
          ENDIF
       ENDDO
    ENDIF

    DO jc = 1, nvertdim
       PRINT *, TRIM( c_vertdim(jc) )
    ENDDO

  END SUBROUTINE getlist

END PROGRAM cdf_compress
