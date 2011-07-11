PROGRAM cdfbottom
  !!======================================================================
  !!                     ***  PROGRAM  cdfbottom  ***
  !!=====================================================================
  !!  ** Purpose : Extract the bottom value for the 3D variables
  !!               which are in the input file. Store the results
  !!               on a similar file, with the same variable name.
  !!
  !!  **  Method:  Uses the corresponding mask file to determine the bottom.
  !!               If no mask found it assumes that 0.0000 values corresponds
  !!               to masked values.
  !!
  !! History : 2.1  : 11/2005  : J.M. Molines : Original code
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER(KIND=4)                            :: jk , jv, jvar, jt        ! dummy loop index
  INTEGER(KIND=4)                            :: ierr                     ! working integer
  INTEGER(KIND=4)                            :: narg, iargc, ijarg       ! argument on line
  INTEGER(KIND=4)                            :: npiglo, npjglo, npk, npt ! size of the domain
  INTEGER(KIND=4)                            :: nvars                    ! number of variables in the input file
  INTEGER(KIND=4)                            :: ncout                    ! ncid of output ncdf file
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, ipko                ! outptut variables : number of levels,
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: id_var, id_varout        ! ncdf varid's

  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zfield                   ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zbot                     ! array to store the bottom value
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zmask                    ! 2D mask at current level
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: tim                      ! time counter of the file

  CHARACTER(LEN=256)                         :: cf_out='bottom.nc'       ! output file name
  CHARACTER(LEN=256)                         :: cf_in, cldum             ! working strings
  CHARACTER(LEN=256)                         :: cv_dep                   ! true name of dep dimension
  CHARACTER(LEN=5)                           :: cv_msk=' '               ! name of the mask variable
  CHARACTER(LEN=1)                           :: ctype=' '                ! point type (T U V ..)
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names              ! array of var name

  TYPE (variable), DIMENSION(:), ALLOCATABLE :: stypvar                  ! structure for variable attribute
  !!--------------------------------------------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfbottom  IN-file [ T | U | V | F]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Create a 2D file with bottom most values for all the variables'
     PRINT *,'       which are in the input 3D file.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       IN-file : input netcdf 3D file.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [ T | U | V | F] : specify the type of grid point on the C-grid' 
     PRINT *,'            if not given, assume that land points are values with 0.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fmsk),' file is required if the grid point is specified' 
     PRINT *,'                  or if the land value is not 0.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables :  same names than input file, long_name attribute is'
     PRINT *,'               prefixed by Bottom '
     STOP
  ENDIF

  ijarg = 1
  CALL getarg (ijarg, cf_in) ; ijarg = ijarg + 1

  IF ( chkfile(cf_in) /= 0 ) STOP  ! missing files

  npiglo = getdim (cf_in,cn_x)
  npjglo = getdim (cf_in,cn_y)
  npk    = getdim (cf_in,cn_z, cdtrue=cv_dep, kstatus=ierr)  ! defautl cn_z is depth

  IF (ierr /= 0 ) THEN
     npk   = getdim (cf_in,'z', cdtrue=cv_dep, kstatus=ierr)
     IF (ierr /= 0 ) THEN
       npk   = getdim (cf_in,'sigma', cdtrue=cv_dep, kstatus=ierr)
        IF ( ierr /= 0 ) THEN
          npk = getdim (cf_in,'nav_lev', cdtrue=cv_dep, kstatus=ierr)
            IF ( ierr /= 0 ) THEN
              PRINT *,' assume file with no depth'
              npk=0
            ENDIF
        ENDIF
     ENDIF
  ENDIF
  npt   = getdim (cf_in,cn_t)


  ALLOCATE (zfield(npiglo,npjglo), zbot(npiglo,npjglo), zmask(npiglo,npjglo))
  ALLOCATE (tim(npt) )

  DO WHILE ( ijarg <= narg )
     CALL getarg (ijarg, ctype ) ; ijarg = ijarg + 1
     IF ( chkfile (cn_fmsk ) ) STOP  ! missing mask file

     SELECT CASE ( ctype )
     CASE ( 'T', 't', 'S', 's' )
        cv_msk='tmask'
     CASE ( 'U', 'u' )
        cv_msk='umask'
     CASE ( 'V', 'v' )
        cv_msk='vmask'
     CASE ( 'F', 'f' )
        cv_msk='fmask'
        PRINT *, 'Be carefull with fmask ... !!!'
     CASE DEFAULT
        PRINT *, ' ERROR : This type of point ', ctype,' is not known !'
        STOP
     END SELECT

  END DO

  ! look for the number of variables in the input file
  nvars = getnvar(cf_in)
  PRINT *,' nvars =', nvars

  ALLOCATE (cv_names(nvars) ,stypvar(nvars))
  ALLOCATE (id_var(nvars), ipk(nvars), id_varout(nvars), ipko(nvars) )

  cv_names(:)=getvarname(cf_in,nvars,stypvar)

  id_var(:)  = (/(jv, jv=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cf_in,nvars,cdep=cv_dep)
  ipko(:)    = 1  ! all variables output are 2D

  WHERE( ipk <= 1 ) cv_names='none'
  DO jvar=1,nvars
    stypvar(jvar)%cname      = cv_names(jvar)
    stypvar(jvar)%caxis      = 'TYX'
    cldum=stypvar(jvar)%clong_name
    stypvar(jvar)%clong_name = 'Bottom '//TRIM(cldum)
  END DO
  ! create output fileset
  ! create output file taking the sizes in cf_in

  ncout = create      (cf_out,   cf_in  , npiglo, npjglo, 1         ) ! 1 level file
  ierr  = createvar   (ncout   , stypvar, nvars , ipko  , id_varout )
  ierr  = putheadervar(ncout   , cf_in  , npiglo, npjglo, 1         )

  DO jvar = 1,nvars
     zfield = 0.
     zbot   = 0.

     IF (cv_names(jvar) == 'none' ) THEN
        ! skip these variable
     ELSE
        PRINT *, ' WORKING with ', TRIM( cv_names(jvar) ), ipk(jvar)
        DO jt = 1, npt
          DO jk = 1, ipk(jvar)
             zmask = 1.
             zfield(:,:) = getvar(cf_in, cv_names(jvar),  jk, npiglo, npjglo, ktime=jt)
             IF ( cv_msk == ' ' ) THEN
                WHERE ( zfield /= 0 )
                   zbot = zfield
                END WHERE
             ELSE
                zmask(:,:) = getvar(cn_fmsk, cv_msk, jk, npiglo, npjglo)
                WHERE ( zmask /= 0 )
                   zbot = zfield
                END WHERE
             ENDIF
          END DO
          ierr = putvar(ncout, id_varout(jvar), zbot, 1, npiglo, npjglo, ktime=jt)
        ENDDO
     ENDIF
  END DO

  tim     = getvar1d(cf_in, cn_vtimec, npt     )
  ierr    = putvar1d(ncout, tim,       npt, 'T')
  ierr    = closeout(ncout)

END PROGRAM cdfbottom
