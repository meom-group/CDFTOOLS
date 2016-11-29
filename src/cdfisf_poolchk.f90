PROGRAM cdfisf_poolchk
  !!======================================================================
  !!                     ***  PROGRAM  cdfisf_poolchk  ***
  !!=====================================================================
  !!  ** Purpose : Build a mask-like file that marks the un-connected point
  !!              in a 3D mask files. Un-connected points are points which
  !!              do not communicate with the open ocean. This case is
  !!              frequent for the ocean cavity below the ice shelves.
  !!
  !!  ** Method  : Use a fillpool3D algorithm
  !!
  !! History :   3.0  : 11/2016  : J.M. Molines. P. Mathiot (original code)
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------


  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2012
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE

  INTEGER(KIND=4) :: ji, jj, jk
  INTEGER(KIND=4) :: npiglo, npjglo, npk
  INTEGER(KIND=4) :: iiseed, ijseed, ikseed
  INTEGER(KIND=4) :: ifill = 2
  INTEGER(KIND=4) :: narg
  INTEGER(KIND=4) :: ncid, id, ierr
  INTEGER(KIND=2), DIMENSION(:,:),   ALLOCATABLE :: itab
  INTEGER(KIND=2), DIMENSION(:,:,:), ALLOCATABLE :: itab3d

  CHARACTER(LEN=255) :: cf_in

  !============================
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *, 'Usage : cdfisf_poolchk mask_file '
     STOP
  ENDIF
  CALL getarg(1, cf_in )

  CALL system ( 'cp '//TRIM(cf_in)//' copymask.nc' )
  cf_in='copymask.nc'
  ierr = NF90_OPEN(cf_in, NF90_WRITE, ncid)
  ierr = NF90_INQ_DIMID(ncid, 'x',id ) ; ierr= NF90_INQUIRE_DIMENSION( ncid, id, len=npiglo)
  ierr = NF90_INQ_DIMID(ncid, 'y',id ) ; ierr= NF90_INQUIRE_DIMENSION( ncid, id, len=npjglo)
  ierr = NF90_INQ_DIMID(ncid, 'z',id ) ; ierr= NF90_INQUIRE_DIMENSION( ncid, id, len=npk   )

  PRINT *, ' NPIGLO = ', npiglo
  PRINT *, ' NPJGLO = ', npjglo

  ALLOCATE ( itab(npiglo, npjglo), itab3d( npiglo, npjglo, npk ) )
  PRINT *, 'WORK WITH tmaskutil'
  ierr = NF90_INQ_VARID(ncid,'tmaskutil', id)
  ierr = NF90_GET_VAR( ncid, id, itab, start=(/1,1,1/), count=(/npiglo,npjglo,1/) )

  ! put a wall at Northern boundary, and on upper level 
  itab(:,npjglo) = 0

  ! iiseed=npiglo/2 ; ijseed = npjglo-2 ; ikseed = 3
  ! iiseed=npiglo/2 ; ijseed = npjglo/2 ; ikseed = 3
  iiseed=400 ; ijseed =490 ; ikseed = 3
  PRINT *, itab( iiseed,ijseed )
  CALL fillpool2d( iiseed, ijseed,        itab,   -ifill )
  PRINT *, '  number of disconected points : ', COUNT(  (itab == 1) )
  ierr = NF90_PUT_VAR(ncid, id, itab, start=(/1,1,1/), count=(/npiglo,npjglo,1/) )
  PRINT *, NF90_STRERROR(ierr)

  PRINT *, 'NOW WORK WITH tmask'
  ierr = NF90_INQ_VARID(ncid,'tmask', id)
  ierr = NF90_GET_VAR( ncid, id, itab3d, start=(/1,1,1,1/), count=(/npiglo,npjglo,npk,1/) )
  itab3d(:,npjglo,:) = 0
  itab3d(:,:,1) = 0
  PRINT *, itab3d( iiseed,ijseed,ikseed)
  CALL fillpool3d( iiseed, ijseed,ikseed, itab3d, -ifill )
  ierr = NF90_PUT_VAR(ncid, id, itab3d, start=(/1,1,1,1/), count=(/npiglo,npjglo,npk,1/) )
  PRINT *, NF90_STRERROR(ierr)
  ierr = NF90_CLOSE( ncid )
  PRINT *, '  number of disconected points : ', COUNT(  (itab3d == 1) )

CONTAINS

  SUBROUTINE fillpool2d(kiseed, kjseed, kdta, kifill)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE fillpool  ***
    !!
    !! ** Purpose :  Replace all area surrounding by mask value by mask value
    !!
    !! ** Method  :  flood fill algorithm
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),                 INTENT(in)    :: kiseed, kjseed
    INTEGER(KIND=4),                 INTENT(in)    :: kifill         ! pool value
    INTEGER(KIND=2), DIMENSION(:,:), INTENT(inout) :: kdta           ! mask

    INTEGER :: ik                       ! number of point change
    INTEGER :: ip                       ! size of the pile
    INTEGER :: ji, jj                   ! loop index
    INTEGER :: iip1, iim1, ii, ij       ! working integer
    INTEGER :: ipiglo, ipjglo           ! size of the domain, infered from kdta size

    INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: ipile    ! pile variable
    INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: idata   ! new bathymetry
    !!----------------------------------------------------------------------
    ! infer domain size from input array
    ipiglo = SIZE(kdta,1)
    ipjglo = SIZE(kdta,2)

    ! allocate variable
    ALLOCATE(ipile(2*ipiglo*ipjglo,2))
    ALLOCATE(idata(ipiglo,ipjglo))

    ! initialise variables
    idata=kdta
    ipile(:,:)=0
    ipile(1,:)=[kiseed,kjseed]
    ip=1; ik=0

    ! loop until the pile size is 0 or if the pool is larger than the critical size
    DO WHILE ( ip /= 0 ) ! .AND. ik < 600000);
       ik=ik+1
       ii=ipile(ip,1); ij=ipile(ip,2)
       IF ( MOD(ik, 10000) == 0 ) PRINT *, 'IP =', ip, ik, ii,ij

       ! update bathy and update pile size
       idata(ii,ij) =kifill
       ipile(ip,:)  =[0,0]; ip=ip-1

       ! check neighbour cells and update pile ( assume E-W periodicity )
       iip1=ii+1; IF ( iip1 == ipiglo+1 ) iip1=2
       iim1=ii-1; IF ( iim1 == 0        ) iim1=ipiglo-1

       IF (idata(ii, ij+1) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[ii  ,ij+1]
       END IF
       IF (idata(ii, ij-1) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[ii  ,ij-1]
       END IF
       IF (idata(iip1, ij) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[iip1,ij  ]
       END IF
       IF (idata(iim1, ij) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[iim1,ij  ]
       END IF

    END DO
    kdta=idata;

    DEALLOCATE(ipile); DEALLOCATE(idata)

  END SUBROUTINE fillpool2d

  SUBROUTINE fillpool3d(kiseed, kjseed, kkseed, kdta, kifill)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE fillpool  ***
    !!
    !! ** Purpose :  Replace all area surrounding by mask value by mask value
    !!
    !! ** Method  :  flood fill algorithm
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),                   INTENT(in)    :: kiseed, kjseed, kkseed
    INTEGER(KIND=4),                   INTENT(in)    :: kifill   ! new bathymetry
    INTEGER(KIND=2), DIMENSION(:,:,:), INTENT(inout) :: kdta     ! new bathymetry

    INTEGER :: ik,iik                   ! number of point change
    INTEGER :: ip                       ! size of the pile
    INTEGER :: ji, jj                   ! loop index
    INTEGER :: ipiglo, ipjglo, ipk      ! size of the domain inferred from kdta
    INTEGER :: iip1, iim1, ii, ij       ! working integer
    INTEGER :: ijp1, ijm1, ikp1, ikm1
    INTEGER(KIND=2), DIMENSION(:,:),   ALLOCATABLE :: ipile    ! pile variable
    INTEGER(KIND=2), DIMENSION(:,:,:), ALLOCATABLE :: idata   ! new bathymetry
    !!----------------------------------------------------------------------
    ! infer domain size from input array
    ipiglo = SIZE(kdta,1)
    ipjglo = SIZE(kdta,2)
    ipk    = SIZE(kdta,3)

    ! allocate variable
    ALLOCATE(ipile(2*ipiglo*ipjglo*ipk,3))
    ALLOCATE(idata(ipiglo,ipjglo,ipk))

    ! initialise variables
    idata=kdta
    ipile(:,:)=0
    ipile(1,:)=[kiseed,kjseed,kkseed]
    ip=1; iik=0

    ! loop until the pile size is 0 or if the pool is larger than the critical size
    DO WHILE ( ip /= 0 ) !.AND. iik < 600000);
       iik=iik+1
       ii=ipile(ip,1); ij=ipile(ip,2) ; ik=ipile(ip,3)
       IF ( MOD( ip, 1000000) == 0 ) PRINT *, 'IP =', ip, iik, ii,ij, ik

       ! update bathy and update pile size
       idata(ii,ij,ik) = kifill
       ipile(ip,:)  =[0,0,0]; ip=ip-1

       ! check neighbour cells and update pile ( assume E-W periodicity )
       iip1=ii+1; IF ( iip1 == ipiglo+1 ) iip1=2
       iim1=ii-1; IF ( iim1 == 0        ) iim1=ipiglo-1
       ijp1=ij+1 ; ijm1 =ij-1 
       ikp1=ik+1 ; ikm1 =ik-1 ; IF (ikm1 == 0 ) ikm1=ik


       IF (idata(ii, ijp1,ik) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[ii  ,ijp1,ik]
       END IF
       IF (idata(ii, ijm1,ik) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[ii  ,ijm1,ik]
       END IF

       IF (idata(iip1, ij,ik) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[iip1,ij  ,ik]
       END IF
       IF (idata(iim1, ij,ik) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[iim1,ij  ,ik]
       END IF
       IF (idata(ii, ij,ikp1) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[ii  ,ij,ikp1]
       END IF

       IF (idata(ii, ij,ikm1) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[ii  ,ij,ikm1]
       END IF

    END DO

    kdta=idata;

    DEALLOCATE(ipile); DEALLOCATE(idata)

  END SUBROUTINE fillpool3d

END PROGRAM cdfisf_poolchk
