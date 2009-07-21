PROGRAM cdfstatcoord
  !!-------------------------------------------------------------------
  !!                 ***  PROGRAM cdfstatcoord ***
  !!
  !!  **  Purpose: Compute statistics about the grid metric versus latitude
  !!  
  !!  **  Method: bins e1 and e2 by latitudes and takes the mean value of each bin
  !!
  !! history:
  !!     Original : J.M. Molines 07/07 for T. Penduff
  !!--------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !!
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: narg, iargc                         !: 
  INTEGER   :: npiglo,npjglo, npk ,npt             !: size of the domain
  INTEGER   :: ngood
  REAL(kind=4) , DIMENSION(:,:), ALLOCATABLE  :: e1, e2, gphi, tmask
  LOGICAL, DIMENSION(:,:), ALLOCATABLE  :: lgood
  REAL(KIND=4) :: binsize=2., rlatmin=-80., rlatmax=90. , rlat, rlat1, rlat2
  REAL(KIND=8) :: e1mean, e2mean

  CHARACTER(LEN=256) :: coord   ='mesh_hgr.nc' , cmask='mask.nc', cvmask='tmask' !:
  TYPE(variable), DIMENSION(3) :: typvar          !: structure for attribute
 
  !!  Read command line
  narg= iargc()
  IF ( narg < 2 ) THEN
     PRINT *,' Usage : cdfstatcoord coordinate-file mask [mask variable name]'
     PRINT *,'    coordinate file is the file where e1t e2t anf gphit can be found'
     PRINT *,'    if mask variable is not tmask, give it as optional argument'
     PRINT *,'    results is given on standard output '
     STOP
  ENDIF

  CALL getarg (1, coord)
  CALL getarg (2, cmask)
  IF ( narg == 3 ) CALL getarg(3,cvmask)

  npiglo= getdim (coord,'x')
  npjglo= getdim (coord,'y')
 
  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo

  ALLOCATE ( e1(npiglo,npjglo)  , e2(npiglo,npjglo) )
  ALLOCATE ( gphi(npiglo,npjglo), tmask(npiglo,npjglo), lgood(npiglo,npjglo)  )

  ! read grid metrics and latitude
  e1=  getvar(coord, 'e1t', 1,npiglo,npjglo)
  e2=  getvar(coord, 'e2t', 1,npiglo,npjglo)
  gphi=  getvar(coord, 'gphit', 1,npiglo,npjglo)
  ! read tmask (1)
  tmask=  getvar(cmask, cvmask, 1,npiglo,npjglo)

  rlat=rlatmin+binsize/2.
  DO WHILE ( rlat <= rlatmax ) 
    rlat1= rlat -binsize/2. ; rlat2 = rlat+binsize/2.
    lgood=.false.
    WHERE ( rlat1 <= gphi .AND. gphi < rlat2  .AND. tmask /= 0 )  lgood=.true.
      ngood=count(lgood)
      IF ( ngood /= 0 ) THEN
        e1mean=SUM( e1, mask=lgood)/ngood
        e2mean=SUM( e2, mask=lgood)/ngood
      ELSE
        e1mean=-999.
        e2mean=-999.
      ENDIF
     PRINT '(f8.3, 3f15.3,i8)', rlat, e1mean, e2mean ,e1mean/e2mean, ngood
     rlat=rlat+binsize
  ENDDO
END PROGRAM cdfstatcoord
