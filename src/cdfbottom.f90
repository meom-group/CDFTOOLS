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
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class bottom
  !!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER(KIND=4)                            :: jk , jv, jvar, jt        ! dummy loop index
  INTEGER(KIND=4)                            :: ierr                     ! working integer
  INTEGER(KIND=4)                            :: idep, idep_max           ! possible depth index, maximum
  INTEGER(KIND=4)                            :: narg, iargc, ijarg       ! argument on line
  INTEGER(KIND=4)                            :: npiglo, npjglo, npk, npt ! size of the domain
  INTEGER(KIND=4)                            :: nvars                    ! number of variables in the input file
  INTEGER(KIND=4)                            :: ncout                    ! ncid of output ncdf file
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, ipko                ! outptut variables : number of levels,
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: id_var, id_varout        ! ncdf varid's

  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zfield                   ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zbot                     ! array to store the bottom value
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zmask                    ! 2D mask at current level

  REAL(KIND=8), DIMENSION(:),    ALLOCATABLE :: dtim                     ! time counter of the file

  CHARACTER(LEN=256)                         :: cf_out='bottom.nc'       ! output file name
  CHARACTER(LEN=256)                         :: cf_in, cldum             ! working strings
  CHARACTER(LEN=256)                         :: cv_dep                   ! true name of dep dimension
  CHARACTER(LEN=5)                           :: cv_msk=' '               ! name of the mask variable
  CHARACTER(LEN=1)                           :: ctype=' '                ! point type (T U V ..)
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names              ! array of var name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: clv_dep               ! array of possible depth name (or 3rd dimension)

  TYPE (variable), DIMENSION(:), ALLOCATABLE :: stypvar                  ! structure for variable attribute

  LOGICAL                                    :: lnc4      = .FALSE.     ! Use nc4 with chunking and deflation
  !!--------------------------------------------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfbottom  -f IN-file [-p C-type] [-o OUT-file] [-nc4] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Create a 2D file with bottom most values for all the variables which '
     PRINT *,'       are in the input 3D file.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file : input netcdf 3D file.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-p C-type] : specify the type of grid point on the C-grid (T|U|V|F).' 
     PRINT *,'               If not given, assume that land points are values with 0.'
     PRINT *,'       [-o OUT-file ]: specify output filename instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ]     : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'               This option is effective only if cdftools are compiled with'
     PRINT *,'               a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fmsk),' file is required if the grid point is specified or if' 
     PRINT *,'               the land value is not 0.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option is used.'
     PRINT *,'         variables :  same names than input file, long_name attribute is'
     PRINT *,'               prefixed by Bottom '
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE (cldum )
     CASE ( '-f'  ) ;  CALL getarg (ijarg, cf_in) ; ijarg = ijarg + 1
        ! options
     CASE ( '-p'  ) ;  CALL getarg (ijarg, ctype) ; ijarg = ijarg + 1
        IF ( chkfile (cn_fmsk )) STOP 99  ! missing files
        SELECT CASE ( ctype )
        CASE ( 'T', 't', 'S', 's' ) ; cv_msk=cn_tmask
        CASE ( 'U', 'u'           ) ; cv_msk=cn_umask
        CASE ( 'V', 'v'           ) ; cv_msk=cn_vmask
        CASE ( 'F', 'f'           ) ; cv_msk=cn_fmask
           ; PRINT *, 'Be carefull with fmask (think of shlat)... !!!'
        CASE DEFAULT                ; PRINT *, ' ERROR : This type of point ', TRIM(ctype),' is not known !' ; STOP 99
        END SELECT
     CASE ( '-o'  ) ;  CALL getarg (ijarg, cf_out) ; ijarg = ijarg + 1
     CASE ( '-nc4') ;  lnc4 = .TRUE.
     CASE DEFAULT   ; PRINT *, ' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( chkfile(cf_in) ) STOP 99  ! missing files

  npiglo = getdim (cf_in,cn_x)
  npjglo = getdim (cf_in,cn_y)

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

  npt   = getdim (cf_in,cn_t)

  ALLOCATE (zfield(npiglo,npjglo), zbot(npiglo,npjglo), zmask(npiglo,npjglo))
  ALLOCATE (dtim(npt) )

  ! look for the number of variables in the input file
  nvars = getnvar(cf_in)
  PRINT *,' nvars =', nvars

  ALLOCATE (cv_names(nvars) ,stypvar(nvars))
  ALLOCATE (id_var(nvars), ipk(nvars), id_varout(nvars), ipko(nvars) )

  cv_names(:)=getvarname(cf_in,nvars,stypvar)
  id_var(:)  = (/(jv, jv=1,nvars)/)

  CALL CreateOutput

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
           END DO ! level-loop
           ierr = putvar(ncout, id_varout(jvar), zbot, 1, npiglo, npjglo, ktime=jt)
        ENDDO  ! time-loop
     ENDIF
  END DO       ! variable-loop

  ierr    = closeout(ncout)
CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ! ipk gives the number of level or 0 if not a T[Z]YX  variable
    ipk(:)     = getipk (cf_in,nvars,cdep=cv_dep)
    ipko(:)    = 1  ! all variables output are 2D

    WHERE( ipk <= 1 ) cv_names='none'
    DO jvar=1,nvars
       stypvar(jvar)%ichunk     = (/npiglo,MAX(1,npjglo/30),1,1 /)
       stypvar(jvar)%cname      = cv_names(jvar)
       stypvar(jvar)%caxis      = 'TYX'
       cldum=stypvar(jvar)%clong_name
       stypvar(jvar)%clong_name = 'Bottom '//TRIM(cldum)
    END DO
    ! create output fileset
    ! create output file taking the sizes in cf_in

    ncout = create      (cf_out,   cf_in  , npiglo, npjglo, 1         , ld_nc4=lnc4 ) ! 1 level file
    ierr  = createvar   (ncout   , stypvar, nvars , ipko  , id_varout , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout   , cf_in  , npiglo, npjglo, 1         )
    dtim  = getvar1d(cf_in, cn_vtimec, npt     )
    ierr  = putvar1d(ncout, dtim,      npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfbottom
