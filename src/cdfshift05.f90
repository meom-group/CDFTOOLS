PROGRAM cdfshift05
  !!======================================================================
  !!                     ***  PROGRAM  cdfshift05  ***
  !!=====================================================================
  !!  ** Purpose : Shift the periodic folding line for ORCA05
  !!               (It may work, after verification for other ORCA grid,
  !!               but I checked in detail the case jperio=6 in this program)
  !!
  !!  ** Method  :  Work on a copy of the input file
  !!                Scan the variables and apply shift on 'map' variable
  !!                Shift and take care of the periodicity afterward
  !!  ** Remarks : As the files to be shifted are not (in general) NEMO
  !!               output, I decided to use netcdf function within this
  !!               program instead of the wrappers in cdfio. Therefore
  !!               this tool is not fully following CDFTOOLS spirit :) !
  !!
  !! History :  4.0  : 03/2021  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  !!   shift05       : perform the shift of 180 degre
  !!----------------------------------------------------------------------
  USE netcdf
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2021
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class file_operations
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER(KIND=4) :: npiglo, npjglo, npk, nt
  INTEGER(KIND=4) :: ji,jj,jk, jt, jvar, jdim
  INTEGER(KIND=4) :: ipivot, ii
  INTEGER(KIND=4) :: narg, ijarg
  INTEGER(KIND=4) :: nvar=5, nvartotal,ndims
  INTEGER(KIND=4) :: ierr, ncid, id, idunlim
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: idims
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: idimids

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: v2di
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: v2do

  CHARACTER(LEN=255) :: cf_in
  CHARACTER(LEN=255) :: cf_out
  CHARACTER(LEN=255) :: cvar
  CHARACTER(LEN=255) :: cldum
  CHARACTER(LEN=5) :: cdimx='x', cdimy='y'
  CHARACTER(LEN=255), DIMENSION(:), ALLOCATABLE :: cvart
  CHARACTER(LEN=255), DIMENSION(:), ALLOCATABLE :: cdimt
  !!----------------------------------------------------------------------
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfshift05 -f IN-file [ -x X-dmn] [-y Y-dmn]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Shift ORCA05 files in order to put the periodic line in the pacific' 
     PRINT *,'       The names of the ''x'' and ''y'' dimension need to be known.'
     PRINT *,'        Time dimension (if any) is supposed to be unlimited. '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file : pass the name of the file to shift' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-x X-dmn] : pass the name of the x dimension if not ''x'' ' 
     PRINT *,'       [-y Y-dmn] : pass the name of the y dimension if not ''y'' ' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       None ' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file :  IN-file.shifted'
     PRINT *,'         variables :  same as input file'
     PRINT *,'      '
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_in ) ; ijarg=ijarg+1
        ! option
     CASE ( '-x'   ) ; CALL getarg(ijarg, cdimx ) ; ijarg=ijarg+1
     CASE ( '-y'   ) ; CALL getarg(ijarg, cdimy ) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO
  cf_out=TRIM(cf_in)//'.shifted'

  !  
  CALL SYSTEM ("cp "//TRIM(cf_in)//" "//TRIM(cf_out) )
  ierr = NF90_OPEN(cf_out,NF90_WRITE,ncid)
  ierr = NF90_INQUIRE( ncid, nVariables=nvartotal, nDimensions=ndims,unlimitedDimId=idunlim)
  PRINT *, 'NDIMS     ', ndims
  PRINT *, 'NVARTOTAL ', nvartotal
  ALLOCATE(cvart(nvartotal), idims(nvartotal), idimids(ndims,nvartotal) )
  ALLOCATE(cdimt(ndims))
  DO jdim = 1 , ndims
    ierr = NF90_INQUIRE_DIMENSION(ncid, jdim, name=cdimt(jdim) )
!   print *, jdim, TRIM(cdimt(jdim))
  ENDDO

  DO jvar = 1, nvartotal
     ierr = NF90_INQUIRE_VARIABLE(ncid, jvar, name=cvart(jvar), ndims=idims(jvar),dimids=idimids(:,jvar))
!    print *, TRIM(cvart(jvar) )
!    DO jdim = 1, idims(jvar)
!      print *, '   dim ',jdim," ",TRIM(cdimt( idimids(jdim,jvar) ))
!   ENDDO
  ENDDO
  ierr = NF90_INQ_DIMID( ncid, cdimx,id) ; ierr=NF90_INQUIRE_DIMENSION(ncid,id, len=npiglo)
  ierr = NF90_INQ_DIMID( ncid, cdimy,id) ; ierr=NF90_INQUIRE_DIMENSION(ncid,id, len=npjglo)
  PRINT *, 'npiglo, npjglo', npiglo, npjglo
  ALLOCATE(v2di(npiglo,npjglo), v2do(npiglo,npjglo))

  DO jvar=1,nvartotal
     cvar=cvart(jvar)
     SELECT CASE (idims(jvar))
     CASE ( 2 )
        IF (cdimt(idimids(1,jvar)) == cdimx .AND. cdimt(idimids(2,jvar)) == cdimy) THEN  ! 2D x,y
           PRINT *,' SHIFTING 2D x y ', TRIM(cvar)
           ierr=NF90_INQ_VARID(ncid,cvar,id)
           ierr=NF90_GET_VAR(ncid,id,v2di,start=(/1,1/), count=(/npiglo,npjglo/) )
           CALL shift05
           ierr=NF90_PUT_VAR(ncid,id,v2do,start=(/1,1/), count=(/npiglo,npjglo/) )
        ENDIF
     CASE ( 3 )
        IF (cdimt(idimids(1,jvar)) == cdimx .AND. cdimt(idimids(2,jvar)) == cdimy) THEN  ! 2D x,y ,?
           IF (idimids(3,jvar) == idunlim ) THEN          ! x,y,t 
             ierr=NF90_INQUIRE_DIMENSION(ncid,idunlim, len=nt)
             PRINT *,' SHIFTING 3D x y t', TRIM(cvar)
             DO jt = 1, nt
               ierr=NF90_INQ_VARID(ncid,cvar,id)
               ierr=NF90_GET_VAR(ncid,id,v2di,start=(/1,1,jt/), count=(/npiglo,npjglo,1/) )
               CALL shift05
               ierr=NF90_PUT_VAR(ncid,id,v2do,start=(/1,1,jt/), count=(/npiglo,npjglo,1/) )
             ENDDO
           ELSE                                           ! x,y,z
             ierr=NF90_INQUIRE_DIMENSION(ncid,idimids(3,jvar), len=npk )
             PRINT *,' SHIFTING 3D x y z', TRIM(cvar)
             DO jk = 1, npk
               ierr=NF90_INQ_VARID(ncid,cvar,id)
               ierr=NF90_GET_VAR(ncid,id,v2di,start=(/1,1,jk/), count=(/npiglo,npjglo,1/) )
               CALL shift05
               ierr=NF90_PUT_VAR(ncid,id,v2do,start=(/1,1,jk/), count=(/npiglo,npjglo,1/) )
             ENDDO
           ENDIF
        ENDIF
     CASE ( 4 )
        ierr=NF90_INQUIRE_DIMENSION(ncid,idunlim, len=nt)
        ierr=NF90_INQUIRE_DIMENSION(ncid,idimids(3,jvar), len=npk )
        PRINT *,' SHIFTING 4D x y z t', TRIM(cvar)
        DO jt = 1, nt
          DO jk=1,npk
           ierr=NF90_INQ_VARID(ncid,cvar,id)
           ierr=NF90_GET_VAR(ncid,id,v2di,start=(/1,1,jk,jt/), count=(/npiglo,npjglo,1,1/) )
           CALL shift05
           ierr=NF90_PUT_VAR(ncid,id,v2do,start=(/1,1,jk,jt/), count=(/npiglo,npjglo,1,1/) )
          ENDDO
        ENDDO
     END SELECT
  ENDDO

  ierr = NF90_CLOSE(ncid)

CONTAINS

  SUBROUTINE shift05 
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE shift05  ***
    !!
    !! ** Purpose :  perform the shifting of v2di into v2do 
    !!
    !! ** Method  :  use all variables in the global space 
    !!               checked for the case jperio = 6 (F-Pivot)
    !!               Although it may work for other cases.
    !!----------------------------------------------------------------------

    ipivot=npiglo/2+1
    DO jj=1,npjglo
       DO ji=2,ipivot-1
          ii=ipivot+ji-2
          v2do(ji,jj) = v2di(ii,jj)
       ENDDO
       DO ji=ipivot,npiglo-1
          ii=ji-ipivot+2
          v2do(ji,jj) = v2di(ii,jj)
       ENDDO
       ! Periodic condition
       v2do(1,jj)=v2do(npiglo-1,jj)
       v2do(npiglo,jj) = v2do(2,jj)

    ENDDO
  END SUBROUTINE shift05

END PROGRAM cdfshift05
