PROGRAM cdfwhereij
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfwhereij  ***
  !!
  !!  **  Purpose  :  Give the values of longitude latitude for a given i, j
  !!  
  !!  **  Method   :  Read the coordinate/mesh_hgr file and look
  !!                  for the glam, gphi variables
  !!                 The point type ( T U V F ) is specified on the command line
  !!
  !! history ;
  !!  Original :  J.M. Molines (May 2005 )
  !!-------------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER :: narg, iargc
  INTEGER :: imin, imax, jmin, jmax
  INTEGER :: npiglo, npjglo

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: glam, gphi
 
  CHARACTER(LEN=80) :: cdum, coord, ctype
  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg /= 6 ) THEN
     PRINT *,' Usage : cdfwhereij  imin imax jmin jmax  coord_file point_type'
     PRINT *,' return the geographical position for the zoomed area '
     PRINT *,' as read in coord_file for the point type specified by point_type'
     PRINT *,' Example : cdfwhereij  200 400 600 750  coordinate_ORCA025.nc F '
     STOP
  ENDIF

  CALL getarg (1, cdum )
  READ(cdum,*) imin
  CALL getarg (2, cdum )
  READ(cdum,*) imax
  CALL getarg (3, cdum )
  READ(cdum,*) jmin
  CALL getarg (4, cdum )
  READ(cdum,*) jmax
  CALL getarg (5, coord )
  CALL getarg (6, ctype )

  npiglo= getdim (coord,'x')
  npjglo= getdim (coord,'y')
  IF ( imax > npiglo ) THEN
    PRINT *,' ERROR : imax is greater than the maximum size ', imax, npiglo
    STOP
  ENDIF

  IF ( jmax > npjglo ) THEN
    PRINT *,' ERROR : jmax is greater than the maximum size ', jmax, npjglo
    STOP
  END IF
  
  ALLOCATE (glam(npiglo,npjglo), gphi(npiglo,npjglo) )

  SELECT CASE ( ctype )
  CASE ('T' , 't' )
     glam(:,:) = getvar(coord, 'glamt',1,npiglo,npjglo)
     gphi(:,:) = getvar(coord, 'gphit',1,npiglo,npjglo)
  CASE ('U','u' )
     glam(:,:) = getvar(coord, 'glamu',1,npiglo,npjglo)
     gphi(:,:) = getvar(coord, 'gphiu',1,npiglo,npjglo)
  CASE ('V','v' )
     glam(:,:) = getvar(coord, 'glamv',1,npiglo,npjglo)
     gphi(:,:) = getvar(coord, 'gphiv',1,npiglo,npjglo)
  CASE ('F','f' )
     glam(:,:) = getvar(coord, 'glamf',1,npiglo,npjglo)
     gphi(:,:) = getvar(coord, 'gphif',1,npiglo,npjglo)
  CASE DEFAULT
     PRINT *,' ERROR : type of point not known: ', TRIM(ctype)
  END SELECT

     PRINT '(2a)'     ,' Type of point   : ', TRIM(ctype)
     PRINT '(a,4i6)'  ,'   I J zoom      : ', imin, imax, jmin, jmax
     PRINT '(a,4f9.3)','   LON LAT zoom  : ', glam(imin,jmin), glam(imax,jmax), gphi(imin,jmin), gphi(imax,jmax)

   END PROGRAM cdfwhereij
