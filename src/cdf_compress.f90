PROGRAM cdf_compress

   USE netcdf
   IMPLICIT NONE

   INTEGER(KIND=4)    :: jd, jv, j3,j4
   INTEGER(KIND=4)    :: n1 , n2, n3, n4
   INTEGER(KIND=4)    :: ncin, ncou
   INTEGER(KIND=4)    :: narg, iargc, ijarg
   INTEGER(KIND=4)    :: idef_lev, nvertdim, ilen, ndim, ierr
   INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: icnk_dmn

   REAL(KIND=8) , DIMENSION(:), ALLOCATABLE    ::  d1d
   REAL(KIND=8) , DIMENSION(:,:), ALLOCATABLE  ::  d2d
   CHARACTER(LEN=255) :: cf_in, cf_ou, cdum
   CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: c_vertdim

   TYPE :: ncfile
     INTEGER(KIND=4)                              :: ncid   ! file id
     INTEGER(KIND=4)                              :: ndims  ! number of variables
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
     INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: idimid  ! dimid of each dimension ( ndims)
     INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ichunk  ! size of chunk (nvar, ndims)
     INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: idimids ! dimids of each variable (nvar, ndims) 
     CHARACTER(LEN=255)                            :: c_fnam ! name of working file
     CHARACTER(LEN=80 ), DIMENSION(:), ALLOCATABLE :: c_vnam ! name of each variable (var)
     CHARACTER(LEN=80 ), DIMENSION(:), ALLOCATABLE :: c_dnam ! name of each dimension (ndims)
     LOGICAL,           DIMENSION(:), ALLOCATABLE :: lconti  ! contiguous flag (nvar)

   END TYPE ncfile

   TYPE(ncfile) :: scdf_in, scdf_ou

! --------------------
! default values
  idef_lev = 1
  nvertdim = 0

  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdf_compress -f FILE_in [-v LIST_vertical_Dimensions ] [-d DEF_lev ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compress input file using netcdf4/HDF5 with deflation level <DEF_lev>.'
     PRINT *,'       Need to specify vertical dimension names if any.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f FILE_in : name of netcdfile to be compressed.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -v LIST_vertical_dimensions : comma separated list of the names of'
     PRINT *,'             the vertical dimensions (in order to set chunk size of 1). '
     PRINT *,'       -d DEF_lev : Specify deflation level ( default is 1).'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : <FILE_in>4 '
     PRINT *,'         variables : same as input ( with same attributes).'
     PRINT *,'      '
     PRINT *,'      '
     STOP
  ENDIF

   ijarg=1
   DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg,cdum) ; ijarg=ijarg+1
     SELECT CASE ( cdum)
     CASE ( '-f' ) ; CALL getarg(ijarg,cf_in) ; ijarg=ijarg+1
     CASE ( '-d' ) ; CALL getarg(ijarg,cdum ) ; ijarg=ijarg+1 ; READ(cdum,*) idef_lev
     CASE ( '-v' ) ; CALL getarg(ijarg,cdum ) ; ijarg=ijarg+1 ; CALL getlist(cdum)
     CASE DEFAULT  ; PRINT *,'  E R R O R : Unknown option :', TRIM(cdum) ; STOP 99
     END SELECT
   END DO
   
   cf_ou=TRIM(cf_in)//'4'

   scdf_in%c_fnam = cf_in

   ! retrieve file structure 
   CALL getncfile(scdf_in)

   scdf_ou = scdf_in
   ! weird but nice ... this copy of the structure also allocate the target ???
 ! print *, size( scdf_ou%nvdim )
 ! print *, scdf_in%nvdim(:)
 ! print *, scdf_ou%nvdim(:)

   ! now change the output structure
   scdf_ou%c_fnam = cf_ou
   ! set deflation level  for each variables ( at this stage, all var with same def level ...
   scdf_ou%ideflat(:) = idef_lev
   ! set the chunk size for each dimension ... 
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
     IF ( icnk_dmn(jd) == 0 ) THEN  ! not treated above
        ilen=scdf_ou%nlen(jd)
        icnk_dmn(jd)=MIN(ilen, 500) !!  500 hard coded !! 
     ENDIF
   ENDDO

    ! Now create the skeleton of the output file
    CALL makncfile(scdf_ou) 

    ! now Read in/ write out ... Loop on variables.
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
           ALLOCATE ( d2d (n1,n2) )
           DO j4 =   1, n4
              DO j3 = 1, n3
               ierr = NF90_GET_VAR(ncin,jv,d2d, start=(/1,1,j3,j4/),count=(/n1,n2,1,1/) )
               ierr = NF90_PUT_VAR(ncou,jv,d2d, start=(/1,1,j3,j4/),count=(/n1,n2,1,1/) )
              ENDDO
           ENDDO
           ierr = NF90_SYNC(ncou)
           DEALLOCATE ( d2d)

       CASE DEFAULT
           PRINT *, 'WARNING : found variable with ', ndim,' dimension : not yet supported : skip it !'
       END SELECT
    ENDDO

    ierr = NF90_CLOSE(ncin)
    ierr = NF90_CLOSE(ncou)
   

   CONTAINS 
    SUBROUTINE getncfile (sd_ncfile)
      TYPE(ncfile), INTENT(inout) :: sd_ncfile
      INTEGER(KIND=4) :: ierr, ncid, nvars, ndims
      INTEGER(KIND=4) :: jd, jv

      ierr = NF90_OPEN(sd_ncfile%c_fnam,NF90_NOWRITE,sd_ncfile%ncid)
      ncid =  sd_ncfile%ncid
      ! now fill in ncfile structure ...
      ierr = NF90_INQUIRE(ncid, nDimensions    = sd_ncfile%ndims  &
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
         ierr = NF90_INQUIRE_VARIABLE(ncid, jv, name          = sd_ncfile%c_vnam(jv)    &
                &                             , xtype         = sd_ncfile%itype(jv)     &
                &                             , ndims         = sd_ncfile%nvdim(jv)     &
                &                             , nAtts         = sd_ncfile%nvatt(jv)     &
                &                             , deflate_level = sd_ncfile%ideflat(jv)   &
                &                             , dimids        = sd_ncfile%idimids(jv,:) &
                &                             , chunksizes    = sd_ncfile%ichunk(jv,:)   )
      ENDDO

      ! leave input file open
      
    END SUBROUTINE getncfile

    SUBROUTINE makncfile (sd_ncfile )
       TYPE(ncfile), INTENT(inout) :: sd_ncfile
       !
       INTEGER(KIND=4) :: ja
       INTEGER(KIND=4) :: ierr, ncid, ilen
       INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ichunksz
       CHARACTER(LEN=80) :: cl_anam

       ! Create output data set
       ierr = NF90_CREATE(sd_ncfile%c_fnam,NF90_NETCDF4,sd_ncfile%ncid)
       ncid = sd_ncfile%ncid

       ! Create dimensions
       DO jd = 1, sd_ncfile%ndims
         IF ( jd == sd_ncfile%iunlim ) THEN
           ilen = NF90_UNLIMITED
         ELSE        
           ilen=sd_ncfile%nlen(jd)
         ENDIF
         ierr = NF90_DEF_DIM(ncid,sd_ncfile%c_dnam(jd),ilen,sd_ncfile%idimid(jd))
       ENDDO

       ! Create variables
       DO jv = 1, sd_ncfile%nvars
         ALLOCATE ( ichunksz(sd_ncfile%nvdim(jv)) )
         DO jd = 1, sd_ncfile%nvdim(jv)
            ichunksz(jd)=icnk_dmn(sd_ncfile%idimids(jv,jd) )
         ENDDO
         ierr=NF90_DEF_VAR(ncid, sd_ncfile%c_vnam(jv) &
              &                   , sd_ncfile%itype(jv)  &
              &                   , sd_ncfile%idimids(jv,1:sd_ncfile%nvdim(jv))  &
              &                   , sd_ncfile%nvid(jv)  &
              &                   , chunksizes=ichunksz(:)  &
              &                   , deflate_level = sd_ncfile%ideflat(jv)  )
         DEALLOCATE ( ichunksz )
       ENDDO
       
       ! Create attribute : copy from scdf_in
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
      CHARACTER(LEN=*) , INTENT(in) :: cd_list
      !
      CHARACTER(LEN=80)  ::cl_list
      INTEGER(KIND=4) :: ilen, ipos, ii
      INTEGER(KIND=4) :: jc

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
