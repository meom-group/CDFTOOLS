PROGRAM cdfchgrid
  !!======================================================================
  !!                     ***  PROGRAM  cdfchgrid  ***
  !!======================================================================
  !!  ** Purpose : Transform an 1442x1021 ORCA025 grid variable into an 
  !!               4322x3059 ORCA12 grid variable.
  !!               No interpolation, only copying one grid cell into 9 grid cells.
  !!
  !!  ** Method  : Store the result on a 'cdfchgrid.nc' file similar to the input file
  !!               (except x and y dimension)
  !!
  !!  ** Restriction  : Caution for mask coherence !
  !!                    This tool is only adapted for drowned field
  !!
  !! History : 3.0 !  08/2012    A. Lecointre   : Original code with Full Doctor form + Lic.
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   chgrid        : Convert 1442x1021 ORCA025 2D var into 4322x3059 ORCA12 var
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id: cdfchgrid.f90 XXX YYYY-MM-DD MM:MM:SSZ molines $
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jt                 ! dummy loop index
  INTEGER(KIND=4)                               :: jvar,jjvar             ! dummy loop index
  INTEGER(KIND=4)                               :: ierr                   ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg     ! argument on line
  INTEGER(KIND=4)                               :: npiglo, npjglo         ! size of the input domain
  INTEGER(KIND=4), PARAMETER                    :: npigloout=4322, npjgloout=3059 ! size of the output domain
  INTEGER(KIND=4)                               :: npk, npt               ! size of the domain
  INTEGER(KIND=4)                               :: nvars                  ! number of variables in the input file
  INTEGER(KIND=4)                               :: ncout                  ! ncid of output ncdf file
  INTEGER(KIND=4), DIMENSION(1)                 :: ipk                    ! output variable : number of levels
  INTEGER(KIND=4), DIMENSION(1)                 :: id_varout              ! ncdf varid's

  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d                    ! array to read a layer of data 
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: u2d                    ! array onto ORCA12-grid
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: tim                    ! time counter of the file
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: dep                    ! depth of the file

  CHARACTER(LEN=256)                            :: cf_out='cdfchgrid.nc'  ! output file name
  CHARACTER(LEN=256)                            :: cf_in                  ! input file name
  CHARACTER(LEN=256)                            :: cv_in                  ! variable name
  CHARACTER(LEN=256)                            :: cldum                  ! working string
  CHARACTER(LEN=256)                            :: cv_dep                 ! true name of dep dimension
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names               ! array of var name

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar                ! structure for variable attribute
  !!--------------------------------------------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfchgrid -f IN-file -var IN-var'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Convert ORCA025-grid variable into ORCA12-grid variable'
     PRINT *,'       No interpolation, only copying one grid cell into nine grid cells'
     PRINT *,'      '
     PRINT *,'     RESTRICTION :'
     PRINT *,'       Caution for mask coherence !'
     PRINT *,'       This tool is only adapted for drowned field'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-var : input ORCA025-grid file'
     PRINT *,'       -var IN-var : input ORCA025-grid variable to be converted'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out)
     PRINT *,'         variable : same name as in input file'
     STOP
  ENDIF
  !!
  ijarg = 1
  ! Read command line
  DO  WHILE (ijarg <=  narg)
     CALL getarg(ijarg,cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-f' )
        CALL getarg(ijarg, cf_in) ; ijarg = ijarg + 1
     CASE ( '-var' )
        CALL getarg(ijarg,cv_in) ; ijarg = ijarg + 1
     CASE DEFAULT
        PRINT *, TRIM(cldum),' : unknown option '
        STOP
     END SELECT
  END DO

  IF ( chkfile(cf_in) ) STOP  ! missing files

  ! get domain dimension from input file
  npiglo = getdim (cf_in, cn_x)
  npjglo = getdim (cf_in, cn_y)
  npk    = getdim (cf_in, cn_z, cdtrue=cv_dep, kstatus=ierr)  ! defautl cn_z is depth
  npt    = getdim (cf_in, cn_t)

  IF ( npk == 0 )  npk = 1  ! assume a 2D variable
  IF ( npt == 0 )  npt = 1  ! assume a 1 time frame file

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE ( v2d(npiglo,npjglo) )
  ALLOCATE ( u2d(npigloout,npjgloout) )
  ALLOCATE ( tim(npt) )
  ALLOCATE ( dep(npk) )

  ! look for the number of variables in the input file
  nvars = getnvar(cf_in)
  ALLOCATE (cv_names(nvars) ,stypvar(nvars))
  cv_names(:)=getvarname(cf_in,nvars,stypvar)

  ! find the number of variable we are interested in
  jvar=0
  DO  WHILE (jvar <=  nvars)
     jvar=jvar+1
     IF ( cv_names(jvar) == cv_in ) jjvar=jvar
  END DO
  ipk(1)=npk
  ncout = create      (cf_out,   cf_in  , npigloout, npjgloout, npk         )
  ierr  = createvar   (ncout   , stypvar(jjvar), 1 , ipk  , id_varout )

  ! get time and write time and get deptht and write deptht
  tim=getvar1d(cf_in,cn_t,npt) ; ierr=putvar1d(ncout,tim,npt,'T')
  dep=getvar1d(cf_in,cv_dep,npk) ; ierr=putvar1d(ncout,dep,npk,'D')

  PRINT *,' Working with ', TRIM(cv_in), npk
  DO jt = 1, npt
     DO jk = 1, npk
        v2d(:,:) = getvar(cf_in, cv_in,  jk, npiglo, npjglo, ktime=jt)
        PRINT *,'level ',jk, 'time ',jt
        CALL chgrid(v2d,u2d,'025to12')
        ierr = putvar ( ncout , id_varout(1), REAL(u2d), jk, npigloout, npjgloout, ktime=jt)
     ENDDO
  ENDDO 

  ierr    = closeout(ncout)

CONTAINS

  SUBROUTINE chgrid (invar,outvar,cc)

  REAL(KIND=4), DIMENSION(npiglo,npjglo), INTENT(in)        :: invar
  REAL(KIND=4), DIMENSION(npigloout,npjgloout), INTENT(out) :: outvar
  CHARACTER(LEN=*), INTENT(in)                              :: cc
  INTEGER(KIND=4)                                           :: iin,jin,iout,jout ! dummy loop index

  SELECT CASE (cc)
  CASE ('025to12')
     DO iin = 2, 1441
        iout=3*iin-4
        jin=1   ! Fill only NORTH and EAST and WEST
        jout=3*jin-2
        outvar(iout  ,jout  ) = invar(iin,jin)
        outvar(iout  ,jout+1) = invar(iin,jin)
        outvar(iout-1,jout  ) = invar(iin,jin)
        outvar(iout+1,jout  ) = invar(iin,jin)
        outvar(iout+1,jout+1) = invar(iin,jin)
        outvar(iout-1,jout+1) = invar(iin,jin)
        DO jin = 2, 1020 ! Fill all: NORTH and SOUTH and EAST and WEST
           jout=3*jin-2
           outvar(iout  ,jout  ) = invar(iin,jin)
           outvar(iout+1,jout  ) = invar(iin,jin)
           outvar(iout-1,jout  ) = invar(iin,jin)
           outvar(iout  ,jout-1) = invar(iin,jin)
           outvar(iout  ,jout+1) = invar(iin,jin)
           outvar(iout+1,jout+1) = invar(iin,jin)
           outvar(iout+1,jout-1) = invar(iin,jin)
           outvar(iout-1,jout+1) = invar(iin,jin)
           outvar(iout-1,jout-1) = invar(iin,jin)
        ENDDO
     ENDDO
     iin=1442
     iout=3*iin-4
     jin=1 ! Fill only NORTH and WEST
     jout=3*jin-2
     outvar(iout  ,jout  ) = invar(iin,jin)
     outvar(iout  ,jout+1) = invar(iin,jin)
     outvar(iout-1,jout  ) = invar(iin,jin)
     outvar(iout-1,jout+1) = invar(iin,jin)
     DO jin = 2, 1020 ! Fill only NORTH and SOUTH and WEST
        jout=3*jin-2
        outvar(iout  ,jout  ) = invar(iin,jin)
        outvar(iout  ,jout-1) = invar(iin,jin)
        outvar(iout  ,jout+1) = invar(iin,jin)
        outvar(iout-1,jout  ) = invar(iin,jin)
        outvar(iout-1,jout-1) = invar(iin,jin)
        outvar(iout-1,jout+1) = invar(iin,jin)
     ENDDO
  CASE ('05to025')
    ! to do ...
  CASE ('05to12')
    ! to do ...
  CASE DEFAULT
     PRINT *, TRIM(cc),'  is not recognized !'
     PRINT *, 'No conversion will be performed'
  END SELECT

  END SUBROUTINE chgrid

END PROGRAM cdfchgrid
