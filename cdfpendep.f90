PROGRAM cdfpendep_new
  !!-------------------------------------------------------------------
  !!              PROGRAM CDFPENDEP
  !!              *****************
  !!
  !!  **  Purpose: Computes penetration depth for passive tracer
  !!                output. This is the ratio between inventory
  !!                and surface concentration (2D) field
  !!  
  !!  **  Method: takes TRC files as input
  !!
  !! history:
  !!    Original:  J.M. Molines (Feb. 2008(
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
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
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: trcinv, trcsurf, pendep
  REAL(KIND=4) ,DIMENSION(1)                  :: timean

  CHARACTER(LEN=256) :: cfiletrc, cfiledia, cfileout='pendep.nc'            !: file name
  CHARACTER(LEN=256) :: cinv='invcfc' , ctrc='cfc11', cdum
  TYPE(variable), DIMENSION(1) :: typvar

  INTEGER    :: ncout
  INTEGER    :: istatus, ierr

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfpendep ''TRC file'' ''DIA file'' [-inv inventory_name  -trc trc_name ]'
     PRINT *,' if not given, inventory name is invcfc, and trc name is cfc11 '
     PRINT *,'   Output on pendep.nc ,variable pendep (m) '
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cfiletrc)
  CALL getarg (2, cfiledia)
  IF ( narg > 2 ) THEN
    jarg=3
    DO WHILE (jarg <= narg )
      CALL getarg(jarg,cdum)
      SELECT CASE (cdum)
       CASE ('-inv') ; jarg=jarg+1 ; CALL getarg(jarg,cinv) ; jarg=jarg+1
       CASE ('-trc') ; jarg=jarg+1 ; CALL getarg(jarg,ctrc) ; jarg=jarg+1
       CASE DEFAULT ; PRINT *, 'option ', TRIM(cdum),' not understood' ; STOP
      END SELECT
    END DO
   ENDIF

  npiglo = getdim (cfiletrc,'x')
  npjglo = getdim (cfiletrc,'y')
  npk    = getdim (cfiletrc,'deptht')

  ipk(1)      = 1
  typvar(1)%name='pendep'
  typvar(1)%units='m'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= 0.
  typvar(1)%valid_max= 10000.
  typvar(1)%long_name='Penetration depth'
  typvar(1)%short_name='pendep'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( trcinv(npiglo,npjglo), trcsurf(npiglo,npjglo), pendep(npiglo,npjglo) )

  ncout =create(cfileout, cfiletrc,npiglo,npjglo,1)

  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfiletrc, npiglo, npjglo,1)

    pendep(:,:)=0.
    trcinv(:,:) = getvar(cfiledia,cinv,1 ,npiglo, npjglo)
    trcsurf(:,:) = getvar(cfiletrc,ctrc,1 ,npiglo, npjglo)
    WHERE( trcsurf /= 0 ) pendep=trcinv/trcsurf
    ierr=putvar(ncout,id_varout(1), pendep, 1 ,npiglo, npjglo)

    timean=getvar1d(cfiletrc,'time_counter',1)
    ierr=putvar1d(ncout,timean,1,'T')
    istatus = closeout(ncout)

END PROGRAM cdfpendep_new
