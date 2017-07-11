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
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE cdftools
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class file_informations
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)    :: ji                         ! dummy loop index
  INTEGER(KIND=4)    :: narg, iargc                ! command line
  INTEGER(KIND=4)    :: ijarg                      ! command line
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
  CHARACTER(LEN=50), DIMENSION(:), ALLOCATABLE :: cfields ! string array to receive
                                                   ! fields of the list_file

  LOGICAL            :: l_file_in=.false.          ! flag for input file
  LOGICAL            :: l_file_ou=.false.          ! flag for output file
  LOGICAL            :: l_append =.false.          ! flag for appending x,y to existing data on line
  LOGICAL            :: l_lonlat =.false.          ! flag for adding lon lat to the output
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  clcoo = cn_fcoo

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :   cdffindij  -w xmin xmax ymin ymax  [-c COOR-file] [-p C-type]...'
     PRINT *,'                 ... [-f LST-file] [-d descriptor] [-o OUT-file] [-A] [-l]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Return the model limit (i,j space) of the geographical window given on' 
     PRINT *,'       the input line. If using -f list_file option, then the output is just'
     PRINT *,'       a single point, not a window, and there are no need to define the '
     PRINT *,'       window with -w.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -w xmin xmax ymin ymax : geographical limits of the window, in lon/lat' 
     PRINT *,'       (relevant only if -f option not used.)'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-c COOR-file ] : specify a particular coordinate file' 
     PRINT *,'                     default is ',TRIM(cn_fcoo)
     PRINT *,'       [-p C-type] : specify the point on the C-grid (T U V F). Default is ',TRIM(cltype),'.'
     PRINT *,'       [-f LST-file ] : LST-file is an ascii file describing the location'
     PRINT *,'                (one per line) of geographical points to be translated to '
     PRINT *,'                model (i,j) point. Unless specified with -d option, this list'
     PRINT *,'                file contains Longitude (X) Latitudes (Y) information.'
     PRINT *,'       [-d descriptor] : descriptor is a string indicating the position of'
     PRINT *,'                X and Y coordinates for the lines of list_file. Default value'
     PRINT *,'                of the descriptor is ''XY''. Any other field on the line is '
     PRINT *,'                indicated with any characterm except X or Y. Example of valid'
     PRINT *,'                descriptor : ''oXYooo'' or ''ooYabcdfXooo'' '
     PRINT *,'       [-A  ] : With this option, output is similar to input with I,J appended'
     PRINT *,'                to the corresponding line.'
     PRINT *,'       [-l  ] : With this option, also output the exact model longitude and '
     PRINT *,'                latitude of the I,J point.'
     PRINT *,'       [-o OUT-file]: write output in text OUT-file instead of standard output.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fcoo),' or the specified coordinates file.' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       Output is done on standard output.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO : '
     PRINT *,'       cdfwhereij'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-w' ) ; CALL getarg(ijarg, cldum  ) ; ijarg=ijarg+1 ; READ(cldum,*) xmin
        ;            CALL getarg(ijarg, cldum  ) ; ijarg=ijarg+1 ; READ(cldum,*) xmax
        ;            CALL getarg(ijarg, cldum  ) ; ijarg=ijarg+1 ; READ(cldum,*) ymin
        ;            CALL getarg(ijarg, cldum  ) ; ijarg=ijarg+1 ; READ(cldum,*) ymax
     CASE ( '-c' ) ; CALL getarg(ijarg, clcoo  ) ; ijarg=ijarg+1
     CASE ( '-p' ) ; CALL getarg(ijarg, cltype ) ; ijarg=ijarg+1
     CASE ( '-f' ) ; CALL getarg(ijarg, cf_list) ; ijarg=ijarg+1 ;  l_file_in=.true.
     CASE ( '-d' ) ; CALL getarg(ijarg, cldes  ) ; ijarg=ijarg+1
     CASE ( '-o' ) ; CALL getarg(ijarg, cf_out ) ; ijarg=ijarg+1 ;  l_file_ou=.true.
     CASE ( '-A' ) ;                                                l_append =.true.
     CASE ( '-l' ) ;                                                l_lonlat =.true.
     CASE DEFAULT  ; PRINT *,' ERROR : ',TRIM(cldum),' unknown option.' ; STOP 99
     END SELECT
  END DO

  IF ( l_file_in) THEN
     ! interpret descriptor
     nfields = LEN(TRIM(cldes) )
     ipx     = INDEX(TRIM(cldes),'X')
     ipy     = INDEX(TRIM(cldes),'Y')
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
999  CONTINUE

  ELSE
     CALL cdf_findij ( xmin, xmax, ymin, ymax, iimin, iimax, ijmin, ijmax, cd_coord=clcoo, cd_point=cltype, cd_verbose='y')
  ENDIF

END PROGRAM cdffindij
