PROGRAM cdfmsk
  !!======================================================================
  !!                     ***  PROGRAM  cdfmsk  ***
  !!=====================================================================
  !!  ** Purpose : Computes the number of land points from the mask
  !!
  !! History : 2.1  : 05/2005  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: jk                   ! dummy loop index
  INTEGER(KIND=4)                           :: npoint               ! number of points
  INTEGER(KIND=4)                           :: narg, iargc          ! browse line
  INTEGER(KIND=4)                           :: npiglo, npjglo, npk  ! size of the domain

  REAL(KIND=4)                              :: zss                  !
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask                ! 2D mask at current level

  CHARACTER(LEN=2048)                        :: cf_msk               ! file name
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmsk MSK-file'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the number of ocean points, land points and display' 
     PRINT *,'       some statistics.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       MSK-file : input mask file (which contains tmask)' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none apart the mask file passed as argument.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       Standard output'
     STOP
  ENDIF
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmsk  maskfile '
     STOP
  ENDIF
  
  CALL getarg (1, cf_msk)

  IF ( chkfile(cf_msk) ) STOP ! missing file

  npiglo = getdim (cf_msk, cn_x)
  npjglo = getdim (cf_msk, cn_y)
  npk    = getdim (cf_msk, cn_z)

  ALLOCATE (zmask(npiglo,npjglo))

  npoint = 0
  DO jk=1, npk
     zmask(:,:) = getvar(cf_msk, 'tmask', jk ,npiglo, npjglo)
     zss        = SUM(zmask)
     npoint     = npoint + zss
  END DO  

  PRINT *, ' Number of Ocean points :', npoint ,'  ',(1.*npoint )/npiglo/npjglo/npk*100,' %'
  PRINT *, ' Number of Land points :', npiglo*npjglo*npk - npoint ,'  ',(npiglo*npjglo*npk -1.*npoint )/npiglo/npjglo/npk*100,' %'

END PROGRAM cdfmsk
