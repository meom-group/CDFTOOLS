PROGRAM cdffindij
  !!======================================================================
  !!                     ***  PROGRAM  cdffindij  ***
  !!=====================================================================
  !!  ** Purpose : Return the window index (imin imax jmin jmax )
  !!               for the geographical windows given on input 
  !!               (longmin longmax latmin matmax)
  !!
  !!  ** Method  : Read the coordinate/mesh_hgr file and look for the glam,
  !!               gphi variables.
  !!               Then use a search algorithm to find the corresponding I J
  !!               The point type ( T U V F ) is specified on the command 
  !!               line as well as the name of the coordinate/mesh hgr file.
  !!
  !! History : 2.1  : 11/2005  : J.M. Molines : Original code
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!           3.0  : 02/2016  : J.M. Molines : add -f, -d, -o options
  !!----------------------------------------------------------------------
  USE cdfio
  USE cdftools
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)    :: ji                         ! dummy loop index
  INTEGER(KIND=4)    :: narg, iargc                ! command line
  INTEGER(KIND=4)    :: ijarg, ireq                ! command line
  INTEGER(KIND=4)    :: iimin, iimax, ijmin, ijmax ! model grid window
  INTEGER(KIND=4)    :: inum=10, iout=6            ! logical unit of assci files
  INTEGER(KIND=4)    :: nfields                    ! number of fields in file_list
  INTEGER(KIND=4)    :: ipx, ipy                   ! field number for X and Y in file_list

  REAL(KIND=4)       :: xmin, xmax, ymin, ymax     ! geographical window
  REAL(KIND=4)       :: zlon, zlat                 ! position of model point

  CHARACTER(LEN=256) :: cltype='F'                 ! point type to search for
  CHARACTER(LEN=256) :: cldum                      ! dummy character variable
  CHARACTER(LEN=256) :: clcoo                      ! dummy character variable
  CHARACTER(LEN=256) :: cf_list                    ! list_file name
  CHARACTER(LEN=256) :: cf_out                     ! output file name
  CHARACTER(LEN=256) :: cldes='XY'                 ! descriptor for input file
  CHARACTER(LEN=50), DIMENSION(:), ALLOCATABLE  :: cfields ! string array to receive
                                                   ! fields of the list_file

  LOGICAL            :: l_file_in=.false.          ! flag for input file
  LOGICAL            :: l_file_ou=.false.          ! flag for output file
  LOGICAL            :: l_append =.false.          ! flag for appending x,y to existing data on line
  LOGICAL            :: l_lonlat =.false.          ! flag for adding lon lat to the output
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  clcoo = cn_fcoo

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg < 4 ) THEN
     PRINT *,' usage :   cdffindij  xmin xmax ymin ymax  [-c COOR-file] [-p point_type]...'
     PRINT *,'                    [-f list_file ] [-d decriptor] [-o output_file] [-a] [-l]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Return the model limit (i,j space) of the geographical window ' 
     PRINT *,'       given on the input line. If using -f list_file option, then the output'
     PRINT *,'       is just a single point, not a window, and xmin, xmax, ymin ymax are not'
     PRINT *,'       used at all.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       xmin xmax ymin ymax : geographical limits of the window, in lon/lat' 
     PRINT *,'       (relevant only if -f option not used.)'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-c COOR-file ] : specify a particular coordinate file' 
     PRINT *,'                     default is ',TRIM(cn_fcoo)
     PRINT *,'       [-p point type] : specify the point on the C-grid (T U V F)'
     PRINT *,'                     default is ',TRIM(cltype)
     PRINT *,'       [-f list_file ] : list_file is an ascii file describing the location'
     PRINT *,'                (one per line) of geographical points to be translated to '
     PRINT *,'                model (i,j) point. Unless specified with -d option, this list'
     PRINT *,'                file contains Longitude (X) Latitudes (Y) information.'
     PRINT *,'       [-d descriptor] : descriptor is a string indicating the position of'
     PRINT *,'                X and Y coordinates for the lines of list_file. Default value'
     PRINT *,'                of the descriptor is ''XY''. Any other field on the line is '
     PRINT *,'                indicated with any characterm except X or Y. Example of valid'
     PRINT *,'                descriptor : ''oXYooo'' or ''ooYabcdfXooo'' '
     PRINT *,'       [-a  ] : With this option, output is similar to input with I,J appended'
     PRINT *,'                to the corresponding line.'
     PRINT *,'       [-l  ] : With this option, also output the exact model longitude and '
     PRINT *,'                latitude of the I,J point.'
     PRINT *,'       [-o output_file] : write output in ascii output_file instead of standard'
     PRINT *,'                output.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fcoo),' or the specified coordinates file.' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       Output is done on standard output.'
     STOP
  ENDIF

  ijarg = 1 ; ireq = 0
  DO WHILE ( ijarg <= narg ) 
    CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
    SELECT CASE ( cldum )
    CASE ( '-c' ) ; CALL getarg(ijarg, clcoo  ) ; ijarg=ijarg+1
    CASE ( '-p' ) ; CALL getarg(ijarg, cltype ) ; ijarg=ijarg+1
    CASE ( '-f' ) ; CALL getarg(ijarg, cf_list) ; ijarg=ijarg+1 ;  l_file_in=.true.
    CASE ( '-d' ) ; CALL getarg(ijarg, cldes  ) ; ijarg=ijarg+1
    CASE ( '-o' ) ; CALL getarg(ijarg, cf_out ) ; ijarg=ijarg+1 ;  l_file_ou=.true.
    CASE ( '-a' ) ;                                                l_append =.true.
    CASE ( '-l' ) ;                                                l_lonlat =.true.
    CASE DEFAULT
       ireq=ireq+1
       SELECT CASE (ireq)
       CASE ( 1 ) ; READ(cldum,*) xmin
       CASE ( 2 ) ; READ(cldum,*) xmax
       CASE ( 3 ) ; READ(cldum,*) ymin
       CASE ( 4 ) ; READ(cldum,*) ymax
       CASE DEFAULT 
         PRINT *,' Too many arguments !' ; STOP 99
       END SELECT
    END SELECT
  END DO
  IF ( l_file_in) THEN
     ! interpret descriptor
     nfields=LEN(TRIM(cldes) )
     ipx=INDEX(TRIM(cldes),'X')
     ipy=INDEX(TRIM(cldes),'Y')
     ALLOCATE( cfields(nfields))
     ! open list_file and loop over lines
     OPEN(inum, FILE=cf_list)
     IF ( l_file_ou ) OPEN(iout, FILE=cf_out)
     DO 
       READ(inum,*,END=999) cfields
       READ(cfields(ipx),*) xmin
       READ(cfields(ipy),*) ymin
       CALL cdf_findij ( xmin, xmin, ymin, ymin, iimin, iimax, ijmin, ijmax, cd_coord=clcoo, & 
           &cd_point=cltype, cd_verbose='n', plonmin=zlon, platmin=zlat)
       IF ( l_append ) THEN
         DO ji = 1, nfields
           WRITE(iout,'(a,x)',advance="no") TRIM(cfields(ji))
         ENDDO
       ENDIF
       IF ( l_lonlat ) THEN
           WRITE(iout,'(g20.9,x)',advance="no")  zlon
           WRITE(iout,'(g20.9,x)',advance="no")  zlat
       ENDIF
         WRITE(iout,'(i10,x)',advance="no") iimin
         WRITE(iout,'(i10,x)'             ) ijmin
     ENDDO
 999 CONTINUE

  ELSE
    CALL cdf_findij ( xmin, xmax, ymin, ymax, iimin, iimax, ijmin, ijmax, cd_coord=clcoo, cd_point=cltype, cd_verbose='y')
  ENDIF

END PROGRAM cdffindij
