PROGRAM cdffracinv
  !!-------------------------------------------------------------------
  !!              PROGRAM CDFFRACINV
  !!              ******************
  !!
  !!  **  Purpose: Computes fraction of inventory for passive tracers 
  !!                output. This is the ratio between inventory at a 
  !!                grid point and total inventory
  !!  
  !!  **  Method: takes TRC files as input
  !!
  !! history:
  !!    Original:  C.O. Dufour (Jul. 2010)
  !!-------------------------------------------------------------------
  !!--------------------------------------------------------------
  !!
  !! * Modules used
  USE cdfio 

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk, jarg
  INTEGER   :: narg, iargc                                  !: 
  INTEGER   :: npiglo,npjglo, npk                                !: size of the domain
  INTEGER, DIMENSION(1) ::  ipk, id_varout
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: trcinvij, fracinv 
  REAL(KIND=4) ,DIMENSION(1)                  :: timean

  CHARACTER(LEN=256) :: cfiletrc, cfileout='fracinv.nc'            !: file name
  CHARACTER(LEN=256) :: cinv='invcfc' , cdum
  TYPE(variable), DIMENSION(1) :: typvar

  INTEGER    :: ncout
  INTEGER    :: istatus, ierr

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdffracinv ''TRC file'' [-inv inventory_name ]'
     PRINT *,' if not given, inventory name is invcfc '
     PRINT *,' Output on fracinv.nc ,variable fracinv (no unit) '
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cfiletrc)
  IF ( narg > 1 ) THEN
    jarg=2
    DO WHILE (jarg <= narg )
      CALL getarg(jarg,cdum)
      SELECT CASE (cdum)
       CASE ('-inv') ; jarg=jarg+1 ; CALL getarg(jarg,cinv) ; jarg=jarg+1
       CASE DEFAULT ; PRINT *, 'option ', TRIM(cdum),' not understood' ; STOP
      END SELECT
    END DO
   ENDIF

  npiglo = getdim (cfiletrc,'x')
  npjglo = getdim (cfiletrc,'y')
  npk    = getdim (cfiletrc,'depth')

  ipk(1)      = 1
  typvar(1)%name='fracinv'
  typvar(1)%units=''
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= 0.
  typvar(1)%valid_max= 10000.
  typvar(1)%long_name='Fraction of inventory'
  typvar(1)%short_name='fracinv'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( trcinvij(npiglo,npjglo), fracinv(npiglo,npjglo) )

  ncout =create(cfileout, cfiletrc,npiglo,npjglo,1)

  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfiletrc, npiglo, npjglo,1)

    fracinv(:,:)=0.
    trcinvij(:,:) = getvar(cfiletrc,cinv,1 ,npiglo, npjglo)
    fracinv(:,:)=trcinvij(:,:)/SUM(trcinvij(:,:))
    ierr=putvar(ncout,id_varout(1), fracinv, 1 ,npiglo, npjglo)

    timean=getvar1d(cfiletrc,'time_counter',1)
    ierr=putvar1d(ncout,timean,1,'T')
    istatus = closeout(ncout)

END PROGRAM cdffracinv
