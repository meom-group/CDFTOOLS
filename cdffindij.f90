PROGRAM cdffindij
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdffindij  ***
  !!
  !!  **  Purpose  :  return the window index (imin imax jmin jmax )
  !!          for the geographical windows given on input (longmin longmax latmin matmax)
  !!  
  !!  **  Method   :  Read the coordinate/mesh_hgr file and look
  !!                  for the glam, gphi variables
  !!                  Then use a seach algorithm to find the corresponding I J
  !!                 The point type ( T U V F ) is specified on the command line
  !!                 as well as the name of the coordinate/mesh hgr file.
  !!
  !! history ;
  !!  Original :  J.M. Molines (November 2005 )
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdftools

  !! * Local variables
  IMPLICIT NONE
  INTEGER :: narg, iargc, niter
  INTEGER :: imin, imax, jmin, jmax
  REAL(KIND=4)                              :: xmin, xmax, ymin, ymax
  CHARACTER(LEN=256) :: cdum, coord='coordinates.nc', ctype='F'

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg < 4 ) THEN
     PRINT *,' Usage : cdffindij  xmin xmax ymin ymax  [coord_file] [point_type]'
     PRINT *,' return the i,j  position for the zoomed area (nearest point ) '
     PRINT *,' as read in coord_file for the point type specified by point_type'
     PRINT *,' Example : cdffindij  -70 15 -20 25  coordinate_ORCA025.nc F '
     STOP
  ENDIF

  CALL getarg (1, cdum ) ; READ(cdum,*) xmin
  CALL getarg (2, cdum ) ; READ(cdum,*) xmax
  CALL getarg (3, cdum ) ; READ(cdum,*) ymin
  CALL getarg (4, cdum ) ; READ(cdum,*) ymax

  ! if 5th argument not given coordinates.nc is assumed
  IF ( narg > 4 ) THEN
     CALL getarg (5, coord )
  ENDIF
  ! if 6th argument not given, assume F point
  IF ( narg == 6 ) THEN
     CALL getarg (6, ctype )
  ENDIF

   CALL cdf_findij ( xmin, xmax, ymin, ymax, imin, imax, jmin, jmax, cd_coord=coord, cd_point=ctype)
END PROGRAM cdffindij
