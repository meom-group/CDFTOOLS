PROGRAM cdfpolymask
  !!======================================================================
  !!                     ***  PROGRAM  cdfpolymask  ***
  !!=====================================================================
  !!  ** Purpose : Create a nc file with 1 into subareas defined as a 
  !!               polygone.
  !!
  !!  ** Method  : Use polylib routine (from finite element mesh generator
  !!               Trigrid)
  !!               Read vertices of polygone in an ascii file an produce a
  !!               resulting file the same shape as file givent in argumment
  !!               (used only for size and header )
  !!
  !! History : 2.1  : 07/2007  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   polymask 
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class mask
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: narg, iargc, ijarg     ! browse line
  INTEGER(KIND=4)                           :: npiglo, npjglo, npk    ! size of the domain
  INTEGER(KIND=4)                           :: ncout                  ! ncid of output file
  INTEGER(KIND=4)                           :: ierr                   ! error status
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout         ! output var levels and varid 

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rpmask                 ! mask array

  REAL(KIND=8), DIMENSION(1)                :: dtim                   ! dummy time counter

  CHARACTER(LEN=256)                        :: cf_ref                 ! name of reference file
  CHARACTER(LEN=256)                        :: cf_poly                ! name of ascii poly file
  CHARACTER(LEN=256)                        :: cf_out='polymask.nc'   ! output file name
  CHARACTER(LEN=256)                        :: cldum                  ! dummy arguments

  TYPE(variable), DIMENSION(1)              :: stypvar                ! output attribute

  LOGICAL                                   :: lreverse=.FALSE.       ! reverse flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfpolymask -p POLY-file -ref REF-file [-r] [-o OUT_file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Create a maskfile with polymask variable having 1 inside the'
     PRINT *,'       polygon, and 0 outside. Option -r revert the behaviour (0 inside,'
     PRINT *,'       1 outside).'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -p POLY-file : input ASCII file describing a polyline in I J grid.'
     PRINT *,'            This file is structured by block, one block corresponding '
     PRINT *,'            to a polygon:'
     PRINT *,'              1rst line of the block gives a polygon name'
     PRINT *,'              2nd line gives the number of vertices (nvert) and a dummy 0'
     PRINT *,'              the block finishes  with nvert pairs of (I,J) describing '
     PRINT *,'              the polygon vertices.'
     PRINT *,'       -ref REF-file  : reference netcdf file for header of polymask file.'
     PRINT *,'             This file will be used to look for domain dimensions, and '
     PRINT *,'             in order to build the output file (nav_lon, nav_lat etc ...)'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        [-r ] : revert option. When used, 0 is inside the polygon,1 outside.'
     PRINT *,'        [-o OUT-file ] : spefify the name of the output mask file instead'
     PRINT *,'                 of ',TRIM(cf_out)
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ', TRIM(cn_polymask)
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1 

  DO WHILE ( ijarg <= narg ) 
     CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum ) 
     CASE ( '-p'   ) ; CALL getarg (ijarg, cf_poly ) ; ijarg = ijarg + 1
     CASE ( '-ref' ) ; CALL getarg (ijarg, cf_ref  ) ; ijarg = ijarg + 1
     CASE ( '-o'   ) ; CALL getarg (ijarg, cf_out  ) ; ijarg = ijarg + 1
     CASE ( '-r'   ) ; lreverse = .TRUE.
     CASE DEFAULT    ; PRINT *,' unknown optional argument (', TRIM(cldum),' )' ; STOP 99
     END SELECT
  END DO

  IF ( chkfile(cf_poly) .OR. chkfile(cf_ref) ) STOP 99 ! missing files

  npiglo = getdim (cf_ref, cn_x)
  npjglo = getdim (cf_ref, cn_y)
  npk    = 1

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo

  CALL CreateOutput

  ALLOCATE( rpmask(npiglo,npjglo) )

  CALL polymask(cf_poly, rpmask) 

  ierr   = putvar(ncout, id_varout(1), rpmask, 1, npiglo, npjglo)
  ierr   = closeout(ncout)

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
    ipk(1)                       = 1
    stypvar(1)%cname             = cn_polymask
    stypvar(1)%cunits            = '1/0'
    stypvar(1)%rmissing_value    = 999.
    stypvar(1)%valid_min         = 0.
    stypvar(1)%valid_max         = 1.
    stypvar(1)%clong_name        = 'Polymask'
    stypvar(1)%cshort_name       = 'polymask'
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TYX'


    ncout = create      (cf_out, cf_ref,  npiglo, npjglo, npk       )
    ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout )
    ierr  = putheadervar(ncout,  cf_ref,  npiglo, npjglo, npk       )
    dtim(:) = 0.d0
    ierr   = putvar1d(ncout, dtim, 1, 'T')

  END SUBROUTINE CreateOutput

  SUBROUTINE polymask( cdpoly, pmask)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE polymask  ***
    !!
    !! ** Purpose : Build polymask from asci polygon file
    !!
    !! ** Method  : Use Poly routines and functions from modpoly module  
    !!
    !!----------------------------------------------------------------------
    USE modpoly

    CHARACTER(LEN=*),             INTENT(in ) :: cdpoly         ! polygon file name
    REAL(KIND=4), DIMENSION(:,:), INTENT(out) :: pmask          ! mask array

    INTEGER(KIND=4)                           :: ji, jj, jjpoly ! dummy loop index
    INTEGER(KIND=4)                           :: infront        ! number of
    REAL(KIND=4)                              :: zin, zout      ! 
    CHARACTER(LEN=256), DIMENSION(jpolys)     :: cl_area        ! name of the areas 
    LOGICAL                                   :: ll_in          ! flag for in/out poly
    !!----------------------------------------------------------------------
    IF ( lreverse ) THEN
       zin = 0. ; zout = 1.
    ELSE
       zin = 1. ; zout = 0.
    ENDIF

    pmask(:,:) = zout
    CALL ReadPoly(cdpoly, infront, cl_area)
    DO jjpoly=1, infront
       CALL PrepPoly(jjpoly)
       DO jj=npjglo, 1, -1
          DO ji=1,npiglo
             CALL InPoly(jjpoly,float(ji), float(jj), ll_in)
             IF (ll_in ) pmask(ji,jj) = zin
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE polymask

END PROGRAM cdfpolymask
