PROGRAM cdfmoy
  !!======================================================================
  !!                     ***  PROGRAM  cdfmoy  ***
  !!=====================================================================
  !!  ** Purpose : Compute mean values for all the variables in a bunch
  !!               of cdf files given as argument
  !!               Store the results on a 'similar' cdf file.
  !!
  !!  ** Method  : Also store the mean squared values for the nn_sqdvar
  !!               variables belonging to cn_sqdvar(:), than can be changed 
  !!               in the nam_cdf_names namelist if wished.
  !!               Optionally order 3 moments for some variables can be
  !!               computed.
  !!
  !! History : 2.0  : 11/2004  : J.M. Molines : Original code
  !!         : 2.1  : 06/2007  : P. Mathiot   : Modif for forcing fields
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!                  04/2015  : S. Leroux    : add nomissincl option
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   varchk2       : check if variable is candidate for square mean
  !!   varchk3       : check if variable is candidate for cubic mean
  !!   zeromean      : substract mean value from input field
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!-----------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jfil,jdep      ! dummy loop index
  INTEGER(KIND=4)                               :: jvar, jv, jt       ! dummy loop index
  INTEGER(KIND=4)                               :: ierr               ! working integer
  INTEGER(KIND=4)                               :: idep, idep_max     ! possible depth index, maximum
  INTEGER(KIND=4)                               :: narg, iargc, ijarg ! browsing command line
  INTEGER(KIND=4)                               :: nfil               ! number of files to average
  INTEGER(KIND=4)                               :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                               :: nvars              ! number of variables in a file
  INTEGER(KIND=4)                               :: ntframe            ! cumul of time frame
  INTEGER(KIND=4)                               :: ncout              ! ncid of output files
  INTEGER(KIND=4)                               :: ncout2             ! ncid of output files
  INTEGER(KIND=4)                               :: ncout3             ! ncid of output files
  INTEGER(KIND=4)                               :: ncout4             ! ncid of output files
  INTEGER(KIND=4)                               :: nperio=4           ! periodic flag
  INTEGER(KIND=4)                               :: iwght              ! weight of variable
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var             ! arrays of var id's
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk4               ! arrays of vertical level for min/max
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout          ! varid's of average vars
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout2         ! varid's of sqd average vars
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout3         ! varid's of cub average vars
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout4         ! varid's of cub average vars

  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: rmask2d             ![from SL]
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d                ! array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: rmax               ! array for maximum value
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: rmin               ! array for minimum value
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: rmean              ! average
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: rmean2             ! squared average
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: rmean3             ! cubic average
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: tim                ! time counter
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: zspval_in          ! time counter
  REAL(KIND=4), DIMENSION(1)                    :: timean             ! mean time

  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dtab, dtab2        ! arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dtab3              ! arrays for cumulated values
  REAL(KIND=8)                                  :: dtotal_time        ! to compute mean time

  CHARACTER(LEN=256)                            :: cf_in              ! input file names
  CHARACTER(LEN=256)                            :: cf_root='cdfmoy'       ! optional root of output files 
  CHARACTER(LEN=256)                            :: cf_out  = 'cdfmoy.nc'  ! output file for average
  CHARACTER(LEN=256)                            :: cf_out2 = 'cdfmoy2.nc' ! output file for squared average
  CHARACTER(LEN=256)                            :: cf_out3 = 'cdfmoy3.nc' ! output file for squared average
  CHARACTER(LEN=256)                            :: cf_out4 = 'cdfmoy_minmax.nc'  ! output file for min/max
  CHARACTER(LEN=256)                            :: cv_dep             ! depth dimension name
  CHARACTER(LEN=256)                            :: cldum              ! dummy string argument
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cf_list            ! list of input files
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_nam             ! array of var name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_nam2            ! array of var2 name for output
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_nam3            ! array of var3 name for output
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_nam4            ! array of var3 name for output
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: clv_dep            ! array of possible depth name (or 3rd dimension)
  
  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar            ! attributes for average values
  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar2           ! attributes for square averaged values
  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar3           ! attributes for cubic averaged values
  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar4           ! attributes for min/max

  LOGICAL                                       :: lcaltmean          ! mean time computation flag
  LOGICAL                                       :: lspval0 = .false.  ! cdfmoy_chsp flag
  LOGICAL                                       :: lcubic  = .false.  ! 3rd momment computation
  LOGICAL                                       :: lzermean = .false. ! flag for zero-mean process
  LOGICAL                                       :: lmax    = .false.  ! flag for min/max computation
  LOGICAL                                       :: lchk    = .false.  ! flag for missing files
  LOGICAL                                       :: lnc4    = .false.  ! flag for netcdf4 output
  LOGICAL                                       :: lnomissincl =.false.! [from SL] flag for excuding gridpoints where some values are missing
  !!----------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmoy list_of_model_files [-spval0] [-cub ] [-zeromean] [-max]'
     PRINT *,'               [-nomissincl] [-nc4 ] [-o output_file_root ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Computes the time average of a list of files given as arguments.' 
     PRINT *,'       The program assumes that all files in the list are of same'
     PRINT *,'       type (shape, variables etc...). '
     PRINT *,'       For some variables, the program also computes the time average '
     PRINT *,'       of the squared variables, which is used in other cdftools '
     PRINT *,'       (cdfeke, cdfrmsssh, cdfstdevw, cdfstddevts ... The actual variables'
     PRINT *,'       selected for squared average are :'
     PRINT '(10x,"- ",a)' , (TRIM(cn_sqdvar(jv)), jv=1, nn_sqdvar)
     PRINT *,'       This selection can be adapted with the nam_cdf_namelist process.'
     PRINT *,'       (See cdfnamelist -i for details).'
     PRINT *,'       If you want to compute the average of already averaged files,'
     PRINT *,'       consider using cdfmoy_weighted instead, in order to take into'
     PRINT *,'       account a particular weight for each file in the list.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       A list of similar model output files. '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [ -spval0 ] :  set missing_value attribute to 0 for all output'
     PRINT *,'               variables and take care of the input missing_value.'
     PRINT *,'               This option is usefull if missing_values differ from files '
     PRINT *,'               to files; it was formely done by cdfmoy_chsp).'
     PRINT *,'       [ -cub ] :  use this option if you want to compute third order moments'
     PRINT *,'               for the eligible variables, which are at present :'
     PRINT '(15x,"- ",a)' , (TRIM(cn_cubvar(jv)), jv=1, nn_cubvar)
     PRINT *,'              This selection can be adapted with the nam_cdf_namelist process.'
     PRINT *,'              (See cdfnamelist -i for details).'
     PRINT *,'       [ -zeromean ] : with this option, the spatial mean value for each '
     PRINT *,'              time frame is substracted from the original field before '
     PRINT *,'              averaging, square averaging and eventually cubic averaging.'
     PRINT *,'       [-max ] : with this option, a file with the minimum and maximum values'
     PRINT *,'              of the variables is created.'
     PRINT *,'       [-nomissincl ] : with this option, the output mean is set to missing' 
     PRINT *,'              value at any gridpoint where the variable contains a  missing'
     PRINT *,'              value for at least one timestep. You should combine with option'
     PRINT *,'              -spval0 if missing values are not 0 in all  the input files.'
     PRINT *,'       [ -nc4 ] Use netcdf4 output with chunking and deflation level 1'
     PRINT *,'               This option is effective only if cdftools are compiled with'
     PRINT *,'               a netcdf library supporting chunking and deflation.'
     PRINT *,'       [ -o output file root ] Default is ', TRIM(cf_root) 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       If -zeromean option is used, need ', TRIM(cn_fhgr),' and ',TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out),' and ',TRIM(cf_out2)
     PRINT *,'       variables : are the same than in the input files. For squared averages' 
     PRINT *,'       _sqd is append to the original variable name.'
     PRINT *,'       If -cub option is used, the file ', TRIM(cf_out3),' is also created'
     PRINT *,'       with _cub append to the original variable name.'
     PRINT *,'       If -max option is used, file ',TRIM(cf_out4),' is also created, with '
     PRINT *,'       same variable names.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfmoy_weighted, cdfstdev'
     PRINT *,'      '
     STOP
  ENDIF


  ALLOCATE ( cf_list(narg) )
  ! look for -spval0 option and set up cf_list, nfil 
  ijarg = 1 
  nfil = 0
  DO WHILE ( ijarg <= narg )
     CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-spval0' )   ! option to reset spval to 0 in the output files
        lspval0 = .true.
     CASE ( '-cub' )   ! option to reset spval to 0 in the output files
        lcubic = .true.
     CASE ( '-zeromean' )   ! option to reset spval to 0 in the output files
        lzermean = .true.
     CASE ( '-max' )   ! option to reset spval to 0 in the output files
        lmax = .true.
     CASE ( '-nomissincl' )   ! [from SL] option to mask the output at gridpoints where some values are missing
        lnomissincl = .true.   ! [from SL]
     CASE ( '-nc4' )   !  allow chunking and deflation on output
        lnc4 = .true.
     CASE ( '-o' )     ! specify root of output files
        CALL getarg (ijarg, cf_root) ; ijarg = ijarg + 1

     CASE DEFAULT         ! then the argument is a file
        nfil          = nfil + 1
        cf_list(nfil) = TRIM(cldum)
     END SELECT
  END DO
  cf_out=TRIM(cf_root)//'.nc'
  cf_out2=TRIM(cf_root)//'2.nc'
  cf_out3=TRIM(cf_root)//'3.nc'
  cf_out4=TRIM(cf_root)//'_minmax.nc'

  IF ( lzermean ) THEN
    lchk = lchk .OR. chkfile ( cn_fhgr )
    lchk = lchk .OR. chkfile ( cn_fmsk )
    IF ( lchk ) STOP ! missing files
  ENDIF

  ! Initialisation from  1rst file (all file are assume to have the same geometry)
  ! time counter can be different for each file in the list. It is read in the
  ! loop for files

  cf_in = cf_list(1)
  IF ( chkfile (cf_in) ) STOP ! missing file
!
  npiglo = getdim (cf_in, cn_x)
  npjglo = getdim (cf_in, cn_y)
  
  ! looking for npk among various possible name
  idep_max=8
  ALLOCATE ( clv_dep(idep_max) )
  clv_dep(:) = (/cn_z,'z','sigma','nav_lev','levels','ncatice','icbcla','icbsect'/)
  idep=1  ; ierr=1000
  DO WHILE ( ierr /= 0 .AND. idep <= idep_max )
     npk  = getdim (cf_in, clv_dep(idep), cdtrue=cv_dep, kstatus=ierr)
     idep = idep + 1
  ENDDO

  IF ( ierr /= 0 ) THEN  ! none of the dim name was found
      PRINT *,' assume file with no depth'
      npk=0
  ENDIF

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk , 'Dep name :' , TRIM(cv_dep)

  ALLOCATE( dtab(npiglo,npjglo), dtab2(npiglo,npjglo), v2d(npiglo,npjglo),rmask2d(npiglo,npjglo)) ! [from SL]
  ALLOCATE( rmean(npiglo,npjglo), rmean2(npiglo,npjglo) )
  IF ( lcubic ) THEN
     ALLOCATE( dtab3(npiglo,npjglo), rmean3(npiglo,npjglo) )
  ENDIF

  nvars = getnvar(cf_in)
  PRINT *,' nvars = ', nvars

  ALLOCATE (cv_nam(nvars), cv_nam2(nvars) )
  ALLOCATE (stypvar(nvars), stypvar2(nvars) )
  ALLOCATE (id_var(nvars),ipk(nvars),id_varout(nvars), id_varout2(nvars)  )
  IF ( lcubic ) THEN
     ALLOCATE (cv_nam3(nvars), stypvar3(nvars), id_varout3(nvars)  )
  ENDIF
  IF ( lmax ) THEN
     ALLOCATE ( ipk4(2*nvars) )
     ALLOCATE ( cv_nam4(2*nvars), stypvar4(2*nvars), id_varout4(2*nvars)  )
     ALLOCATE ( rmin(npiglo, npjglo), rmax(npiglo,npjglo) )
  ENDIF

  ! get list of variable names and collect attributes in stypvar (optional)
  cv_nam(:) = getvarname(cf_in,nvars,stypvar)

  ! choose chunk size for output ... not easy not used if lnc4=.false. but anyway ..
  DO jv = 1, nvars
     stypvar(jv)%ichunk=(/npiglo,MAX(1,npjglo/30),1,1 /)
     stypvar2(jv)%ichunk=(/npiglo,MAX(1,npjglo/30),1,1 /)
  ENDDO

  IF ( lspval0 ) THEN 
     ALLOCATE ( zspval_in(nvars) )
     zspval_in(:) = stypvar(:)%rmissing_value
     stypvar(:)%rmissing_value = 0.
  ENDIF

  IF ( lcubic) THEN
     ! force votemper to be squared saved
     nn_sqdvar = nn_sqdvar + 1
     cn_sqdvar(nn_sqdvar) = TRIM(cn_votemper)  
  ENDIF

  DO jvar = 1, nvars
     ! variables that will not be computed or stored are named 'none'
     IF ( varchk2 ( cv_nam(jvar) ) ) THEN 
        cv_nam2(jvar)                    = TRIM(cv_nam(jvar))//'_sqd'
        stypvar2(jvar)%cname             = TRIM(stypvar(jvar)%cname)//'_sqd'         ! name
        stypvar2(jvar)%cunits            = '('//TRIM(stypvar(jvar)%cunits)//')^2'    ! unit
        stypvar2(jvar)%rmissing_value    = stypvar(jvar)%rmissing_value              ! missing_value
        stypvar2(jvar)%valid_min         = 0.                                        ! valid_min = zero
        stypvar2(jvar)%valid_max         = stypvar(jvar)%valid_max**2                ! valid_max *valid_max
        stypvar2(jvar)%scale_factor      = 1.
        stypvar2(jvar)%add_offset        = 0.
        stypvar2(jvar)%savelog10         = 0.
        stypvar2(jvar)%clong_name        = TRIM(stypvar(jvar)%clong_name)//'_Squared'  ! 
        stypvar2(jvar)%cshort_name       = TRIM(stypvar(jvar)%cshort_name)//'_sqd'     !
        stypvar2(jvar)%conline_operation = TRIM(stypvar(jvar)%conline_operation) 
        stypvar2(jvar)%caxis             = TRIM(stypvar(jvar)%caxis) 
     ELSE
         cv_nam2(jvar) = 'none'
     END IF

     ! check for cubic average
     IF ( lcubic ) THEN
       IF ( varchk3 ( cv_nam(jvar) ) ) THEN 
          cv_nam3(jvar)                    = TRIM(cv_nam(jvar))//'_cub'
          stypvar3(jvar)%cname             = TRIM(stypvar(jvar)%cname)//'_cub'         ! name
          stypvar3(jvar)%cunits            = '('//TRIM(stypvar(jvar)%cunits)//')^3'    ! unit
          stypvar3(jvar)%rmissing_value    = stypvar(jvar)%rmissing_value              ! missing_value
          stypvar3(jvar)%valid_min         = 0.                                        ! valid_min = zero
          stypvar3(jvar)%valid_max         = stypvar(jvar)%valid_max**3                ! valid_max *valid_max
          stypvar3(jvar)%scale_factor      = 1.
          stypvar3(jvar)%add_offset        = 0.
          stypvar3(jvar)%savelog10         = 0.
          stypvar3(jvar)%clong_name        = TRIM(stypvar(jvar)%clong_name)//'_Cubed'   ! 
          stypvar3(jvar)%cshort_name       = TRIM(stypvar(jvar)%cshort_name)//'_cub'    !
          stypvar3(jvar)%conline_operation = TRIM(stypvar(jvar)%conline_operation) 
          stypvar3(jvar)%caxis             = TRIM(stypvar(jvar)%caxis) 

          stypvar3(jvar)%ichunk=(/npiglo,MAX(1,npjglo/30),1,1 /)
       ELSE
          cv_nam3(jvar) = 'none'
       END IF
     ENDIF  

     IF ( lmax ) THEN
          cv_nam4(jvar)                    = TRIM(cv_nam(jvar))//'_max'
          stypvar4(jvar)%cname             = TRIM(stypvar(jvar)%cname)//'_max'         ! name
          stypvar4(jvar)%cunits            = '('//TRIM(stypvar(jvar)%cunits)//')'      ! unit
          stypvar4(jvar)%rmissing_value    = stypvar(jvar)%rmissing_value              ! missing_value
          stypvar4(jvar)%valid_min         = 0.                                        ! valid_min = zero
          stypvar4(jvar)%valid_max         = stypvar(jvar)%valid_max                   ! valid_max *valid_max
          stypvar4(jvar)%scale_factor      = 1.
          stypvar4(jvar)%add_offset        = 0.
          stypvar4(jvar)%savelog10         = 0.
          stypvar4(jvar)%clong_name        = TRIM(stypvar(jvar)%clong_name)//'_max'   ! 
          stypvar4(jvar)%cshort_name       = TRIM(stypvar(jvar)%cshort_name)//'_max'  !
          stypvar4(jvar)%conline_operation = TRIM(stypvar(jvar)%conline_operation)
          stypvar4(jvar)%caxis             = TRIM(stypvar(jvar)%caxis)

          stypvar4(jvar)%ichunk=(/npiglo,MAX(1,npjglo/30),1,1 /)
      
          cv_nam4(nvars+jvar)                    = TRIM(cv_nam(jvar))//'_min'
          stypvar4(nvars+jvar)%cname             = TRIM(stypvar(jvar)%cname)//'_min'         ! name
          stypvar4(nvars+jvar)%cunits            = '('//TRIM(stypvar(jvar)%cunits)//')'      ! unit
          stypvar4(nvars+jvar)%rmissing_value    = stypvar(jvar)%rmissing_value              ! missing_value
          stypvar4(nvars+jvar)%valid_min         = 0.                                        ! valid_min = zero
          stypvar4(nvars+jvar)%valid_max         = stypvar(jvar)%valid_max                   ! valid_max *valid_max
          stypvar4(nvars+jvar)%scale_factor      = 1.
          stypvar4(nvars+jvar)%add_offset        = 0.
          stypvar4(nvars+jvar)%savelog10         = 0.
          stypvar4(nvars+jvar)%clong_name        = TRIM(stypvar(jvar)%clong_name)//'_min'   ! 
          stypvar4(nvars+jvar)%cshort_name       = TRIM(stypvar(jvar)%cshort_name)//'_min'  !
          stypvar4(nvars+jvar)%conline_operation = TRIM(stypvar(jvar)%conline_operation)
          stypvar4(nvars+jvar)%caxis             = TRIM(stypvar(jvar)%caxis)

          stypvar4(nvars+jvar)%ichunk=(/npiglo,MAX(1,npjglo/30),1,1 /)
     ENDIF


  END DO

  id_var(:)  = (/(jv, jv=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cf_in,nvars,cdep=cv_dep)
  WHERE( ipk == 0 ) cv_nam='none'
  IF ( lmax ) THEN 
      ipk4(1      :nvars  ) = ipk(1:nvars)
      ipk4(nvars+1:2*nvars) = ipk(1:nvars)
      WHERE( ipk4 == 0 ) cv_nam4='none'
  ENDIF
                stypvar (:)%cname = cv_nam
                stypvar2(:)%cname = cv_nam2
  IF ( lcubic ) stypvar3(:)%cname = cv_nam3
  IF ( lmax   ) stypvar4(:)%cname = cv_nam4

  ! create output file taking the sizes in cf_in
  ncout  = create      (cf_out,  cf_in,    npiglo, npjglo, npk, cdep=cv_dep, ld_nc4=lnc4)
  ierr   = createvar   (ncout ,  stypvar,  nvars,  ipk,    id_varout       , ld_nc4=lnc4)
  ierr   = putheadervar(ncout,   cf_in,    npiglo, npjglo, npk, cdep=cv_dep      )

  ncout2 = create      (cf_out2, cf_in,    npiglo, npjglo, npk, cdep=cv_dep, ld_nc4=lnc4)
  ierr   = createvar   (ncout2,  stypvar2, nvars,  ipk,    id_varout2      , ld_nc4=lnc4)
  ierr   = putheadervar(ncout2,  cf_in,    npiglo, npjglo, npk, cdep=cv_dep      )

  IF ( lcubic) THEN
     ncout3 = create      (cf_out3, cf_in,    npiglo, npjglo, npk, cdep=cv_dep, ld_nc4=lnc4)
     ierr   = createvar   (ncout3,  stypvar3, nvars,  ipk,    id_varout3      , ld_nc4=lnc4)
     ierr   = putheadervar(ncout3,  cf_in,    npiglo, npjglo, npk, cdep=cv_dep      )
  ENDIF

  IF ( lmax ) THEN
     ncout4 = create      (cf_out4, cf_in,    npiglo, npjglo, npk, cdep=cv_dep, ld_nc4=lnc4)
     ierr   = createvar   (ncout4,  stypvar4, 2*nvars,  ipk4,    id_varout4   , ld_nc4=lnc4)
     ierr   = putheadervar(ncout4,  cf_in,    npiglo, npjglo, npk, cdep=cv_dep      )
  ENDIF

  lcaltmean=.TRUE.
  DO jvar = 1,nvars
     iwght=0
     IF ( cv_nam(jvar) == cn_vlon2d .OR. &     ! nav_lon
          cv_nam(jvar) == cn_vlat2d ) THEN     ! nav_lat
        ! skip these variable
     ELSE
        PRINT *,' Working with ', TRIM(cv_nam(jvar)), ipk(jvar)
        DO jk = 1, ipk(jvar)
           PRINT *,'level ',jk
           dtab(:,:) = 0.d0 ; dtab2(:,:) = 0.d0 ; dtotal_time = 0.
           rmask2d(:,:) = 1.
           IF ( lcubic ) THEN  ; dtab3(:,:) = 0.d0                       ; ENDIF
           IF ( lmax   ) THEN  ; rmin (:,:) = 1.e20 ; rmax(:,:) = -1.e20 ; ENDIF
           ntframe = 0
           DO jfil = 1, nfil
              cf_in = cf_list(jfil)
              IF ( jk == 1 ) THEN
                  IF ( chkfile (cf_in) ) STOP ! missing file
                  iwght=iwght+MAX(1,INT(getatt( cf_in, cv_nam(jvar), 'iweight')))
              ENDIF

              npt = getdim (cf_in, cn_t)
              IF ( lcaltmean )  THEN
                 ALLOCATE ( tim(npt) )
                 tim         = getvar1d(cf_in, cn_vtimec, npt)
                 dtotal_time = dtotal_time + SUM(DBLE(tim(:)))
                 DEALLOCATE (tim )
              END IF
              DO jt=1,npt
                ntframe = ntframe + 1
                v2d(:,:)  = getvar(cf_in, cv_nam(jvar), jk ,npiglo, npjglo,ktime=jt )
                IF ( lspval0  )  WHERE (v2d == zspval_in(jvar))  v2d = 0.  ! change missing values to 0
                WHERE (v2d == 0.) rmask2d = 0.                              ! [from SL]
                IF ( lzermean ) CALL zeromean (jk, v2d )
                dtab(:,:) = dtab(:,:) + v2d(:,:)*1.d0
                IF (cv_nam2(jvar) /= 'none' ) dtab2(:,:) = dtab2(:,:) + v2d(:,:)*v2d(:,:)*1.d0
                IF ( lcubic ) THEN
                   IF (cv_nam3(jvar) /= 'none' ) dtab3(:,:) = dtab3(:,:) + v2d(:,:)*v2d(:,:)*v2d(:,:) *1.d0
                ENDIF
                IF ( lmax ) THEN
                  rmax(:,:) = MAX(v2d(:,:),rmax(:,:))
                  rmin(:,:) = MIN(v2d(:,:),rmin(:,:))
                ENDIF
              ENDDO
           END DO
           ! finish with level jk ; compute mean (assume spval is 0 )
           rmean(:,:) = dtab(:,:)/ntframe
           IF ( lnomissincl ) rmean(:,:) = rmean(:,:)*(rmask2d(:,:)*1.d0)    ! [from SL]
           IF (cv_nam2(jvar) /= 'none' ) THEN
                rmean2(:,:) = dtab2(:,:)/ntframe
                IF ( lnomissincl ) rmean2(:,:) = rmean2(:,:)*(rmask2d(:,:)*1.d0)    ! [from SL]
           ENDIF
           IF ( lcubic ) THEN
              IF (cv_nam3(jvar) /= 'none' ) THEN 
                  rmean3(:,:) = dtab3(:,:)/ntframe
                  IF ( lnomissincl ) rmean3(:,:) = rmean3(:,:)*(rmask2d(:,:)*1.d0)    ! [from SL]
              ENDIF
           ENDIF

           ! store variable on outputfile
           ierr = putvar(ncout, id_varout(jvar), rmean, jk, npiglo, npjglo, kwght=iwght)
           IF (cv_nam2(jvar) /= 'none' ) THEN 
               ierr = putvar(ncout2, id_varout2(jvar), rmean2, jk, npiglo, npjglo, kwght=iwght)
           ENDIF

           IF ( lcubic) THEN
              IF (cv_nam3(jvar) /= 'none' ) THEN 
                 ierr = putvar(ncout3, id_varout3(jvar), rmean3, jk, npiglo, npjglo, kwght=iwght)
              ENDIF
           ENDIF
           IF ( lmax  ) THEN
                 ierr = putvar(ncout4, id_varout4(      jvar), rmax, jk, npiglo, npjglo, kwght=iwght)
                 ierr = putvar(ncout4, id_varout4(nvars+jvar), rmin, jk, npiglo, npjglo, kwght=iwght)
           ENDIF

           IF (lcaltmean )  THEN
              timean(1) = dtotal_time/ntframe
                          ierr = putvar1d(ncout,  timean, 1, 'T')
                          ierr = putvar1d(ncout2, timean, 1, 'T')
              IF (lcubic) ierr = putvar1d(ncout3, timean, 1, 'T')
              IF (lmax  ) ierr = putvar1d(ncout4, timean, 1, 'T')
           END IF

           lcaltmean=.FALSE. ! tmean already computed
        END DO  ! loop to next level
     END IF
  END DO ! loop to next var in file

                ierr = closeout(ncout)
                ierr = closeout(ncout2)
  IF ( lcubic ) ierr = closeout(ncout3 ) 
  IF ( lmax   ) ierr = closeout(ncout4 ) 

CONTAINS 

  LOGICAL FUNCTION varchk2 ( cd_var ) 
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION varchk2  ***
    !!
    !! ** Purpose : Return true if cd_var is candidate for mean squared value  
    !!
    !! ** Method  : List of candidate is established in modcdfnames, and
    !!              can be changed via the nam_cdf_names namelist   
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_var

    INTEGER(KIND=4)              :: jv
    !!----------------------------------------------------------------------
    varchk2 = .FALSE.
    DO jv = 1, nn_sqdvar 
      IF ( cd_var == cn_sqdvar(jv) ) THEN
       varchk2 = .TRUE.
       exit
      ENDIF
    ENDDO 
    
  END FUNCTION varchk2

  LOGICAL FUNCTION varchk3 ( cd_var )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION varchk3  ***
    !!
    !! ** Purpose : Return true if cd_var is candidate for cubic mean average
    !!
    !! ** Method  : List of candidate is established in modcdfnames, and
    !!              can be changed via the nam_cdf_names namelist
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_var

    INTEGER(KIND=4)              :: jv
    !!----------------------------------------------------------------------
    varchk3 = .FALSE.
    DO jv = 1, nn_cubvar
      IF ( cd_var == cn_cubvar(jv) ) THEN
       varchk3 = .TRUE.
       exit
      ENDIF
    ENDDO

  END FUNCTION varchk3

  SUBROUTINE zeromean(kk, ptab)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE zeromean  ***
    !!
    !! ** Purpose :  Computes the spatial average of argument and
    !!               and substract it from the field 
    !!
    !! ** Method  :  requires the horizontal metrics 
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),              INTENT(   in) :: kk
    REAL(KIND=4), DIMENSION(:,:), INTENT(inout) :: ptab

    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE, SAVE :: ze2, ze1, tmask, tmask0

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: dareas
    REAL(KIND=8),                              SAVE :: darea
    REAL(KIND=8)                                    :: dmean

    LOGICAL, SAVE                                   :: lfirst=.true.
    !!----------------------------------------------------------------------

   IF (lfirst) THEN
     lfirst=.false.
     ! read e1 e2 and tmask ( assuming this prog only deal with T-points)
     ALLOCATE ( ze1(npiglo, npjglo), ze2(npiglo,npjglo)       )
     ALLOCATE ( tmask(npiglo,npjglo), tmask0(npiglo,npjglo)   )
     ALLOCATE ( dareas(npiglo,npjglo)                         )

     ze1(:,:)     = getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo)
     ze2(:,:)     = getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo)
     dareas(:,:)  = ze1(:,:) * ze2(:,:) *1.d0
   ENDIF
     tmask0(:,:)  = getvar(cn_fmsk, 'tmask', kk, npiglo, npjglo)
     tmask = tmask0
     tmask(1,:)=0 ; tmask(npiglo,:)=0 ; tmask(:,1) = 0.; tmask(:,npjglo) = 0 

     IF ( nperio == 3 .OR. nperio == 4 ) THEN
       tmask(npiglo/2+1:npiglo,npjglo-1) = 0.
     ENDIF


     darea = SUM( dareas * tmask )

     IF ( darea /= 0.d0 ) THEN
        dmean = SUM( ptab * dareas ) / darea
     ELSE
        dmean = 0.d0
     ENDIF

     WHERE ( ptab /= 0 )  ptab = ( ptab - dmean ) * tmask0
     
  END SUBROUTINE zeromean


END PROGRAM cdfmoy
