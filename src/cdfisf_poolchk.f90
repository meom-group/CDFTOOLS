PROGRAM cdfisf_poolchk
  !!======================================================================
  !!                     ***  PROGRAM  cdfisf_poolchk  ***
  !!=====================================================================
  !!  ** Purpose : Build a mask-like file that marks the un-connected point
  !!              in a 3D mask files. Un-connected points are points which
  !!              do not communicate with the open ocean. This case is
  !!              frequent for the ocean cavity below the ice shelves.
  !!
  !!  ** Method  : Use a fillpool3D algorithm. This program use 3D arrays
  !!               and may be very memory consuming for big domains.
  !!
  !! History :   3.0  : 11/2016  : J.M. Molines. P. Mathiot (original code)
  !!----------------------------------------------------------------------
  !!   routines                                         : description
  !!----------------------------------------------------------------------
  !!  fillpool2d(kiseed, kjseed,         kdta, kifill)  : 2D fill algo
  !!  fillpool3d(kiseed, kjseed, kkseed, kdta, kifill)  : 3D fill algo
  !!----------------------------------------------------------------------
  USE netcdf
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2012
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4) :: ji, jj, jk
  INTEGER(KIND=4) :: ierr, ncout
  INTEGER(KIND=4) :: npiglo, npjglo, npk
  INTEGER(KIND=4) :: iiseed, ijseed, ikseed
  INTEGER(KIND=4) :: ifill = 2
  INTEGER(KIND=4) :: narg, ijarg
  INTEGER(KIND=4) :: ncid, id, ierr
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout          ! varid's of average vars

  INTEGER(KIND=2), DIMENSION(:,:),   ALLOCATABLE :: itab
  INTEGER(KIND=2), DIMENSION(:,:,:), ALLOCATABLE :: itab3d

  CHARACTER(LEN=255) :: cf_in
  CHARACTER(LEN=255) :: cf_out='poolmask.nc'
  CHARACTER(LEN=255) :: cdum

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar            ! attributes for average values

  LOGICAL            :: lnc4=.FALSE.
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfisf_poolchk -m MASK-file [-nc4] [-o OUT-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Produce a netcdf mask file with 1 everywhere, except for points '
     PRINT *,'       not connected to the open ocean (Frequent for cavities below '
     PRINT *,'       ice-shelves), which have 0 value. Both 3D and 2D variables are'
     PRINT *,'       created, the 2D variables beiing used for cdfisf_forcing.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -m MASK-file : name of the input mask file, holding the tmask variable.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -nc4 : use netcdf4 with chunking and deflation for the output.'
     PRINT *,'       -o OUT-file : name of the output file. [Default : ',TRIM(cf_out),' ]' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       Only the mask file given as argument' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out),' unless -o option is used.'
     PRINT *,'         variables : tmask_pool3d, tmask_pool2d'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      cdfisf_fill, cdfisf_forcing, cdfisf_rnf' 
     PRINT *,'      '
     STOP
  ENDIF
  
  ijarg=1
  DO WHILE ( ijarg <= narg )
     CALL getarg(1, cdum) ; ijarg = ijarg+1
     SELECT CASE ( cdum )
     CASE ( '-m'  ) ; CALL getarg(1, cf_in ) ; ijarg = ijarg+1
     CASE ( '-o'  ) ; CALL getarg(1, cf_out) ; ijarg = ijarg+1
     CASE ( '-nc4') ; lnc4=.TRUE.
     CASE DEFAULT 
        PRINT *,' Unknown option : ', TRIM(cdum)
        STOP
     END SELECT
  ENDDO

  npiglo = getdim (cf_in,cn_x )
  npjglo = getdim (cf_in,cn_y )
  npk    = getdim (cf_in,cn_z )

  PRINT *, ' NPIGLO = ', npiglo
  PRINT *, ' NPJGLO = ', npjglo
  PRINT *, ' NPK    = ', npk

  ALLOCATE ( itab(npiglo, npjglo), itab3d( npiglo, npjglo, npk ) )
  ALLOCATE ( ipk(2), id_varout(2), stypvar(2))

  CALL CreateOutput

  PRINT *, 'NOW WORK WITH tmask'
  itab3d(:,:,:) = getvar3d (cf_in, cn_tmask, npiglo, npjglo, npk)
  ! set limits  for fillpool algo
  itab3d(:,npjglo,:) = 0
  itab3d(:,:,     1) = 0

  CALL fillpool3d( iiseed, ijseed,ikseed, itab3d, -ifill )
  
  DO jk = 1, npk 
    ierr = putvar( ncout, id_varout(2), itab3d(:,:,jk), jk, npiglo, npjglo)
  ENDDO
  PRINT *, '  Number of disconected points : ', COUNT(  (itab3d == 1) )

CONTAINS
  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create the output file. This is done outside the main
    !!               in order to increase readability of the code. 
    !!
    !! ** Method  :  Use global variables, defined in mail 
    !!----------------------------------------------------------------------
  ! define new variables for output
  stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
  stypvar(1)%cname             = 'tmask_pool2d'
  stypvar(1)%rmissing_value    =  -99.
  stypvar(1)%valid_min         =  0.
  stypvar(1)%valid_max         =  1.
  stypvar(1)%clong_name        = '2d isf pool mask'
  stypvar(1)%cshort_name       = 'tmask_pool2d'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TYX'
  ipk(1) = 1  !  2D
  ! define new variables for output
  stypvar(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
  stypvar(2)%cname             = 'tmask_pool3d'
  stypvar(2)%rmissing_value    =  -99.
  stypvar(2)%valid_min         =  0.
  stypvar(2)%valid_max         =  1
  stypvar(2)%clong_name        = '3d isf pool mask'
  stypvar(2)%cshort_name       = 'tmask_pool3d'
  stypvar(2)%conline_operation = 'N/A'
  stypvar(2)%caxis             = 'TYX'
  ipk(2) = npk  !  3D


  END SUBROUTINE CreateOutput

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
