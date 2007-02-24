PROGRAM cdfimprovechk
  !!-------------------------------------------------------------------
  !!             ***  PROGRAM cdfimprovechk  ***
  !!
  !!  **  Purpose:  Estimate the improvement/deterioration
  !!                of a test run, compared with a reference run
  !!                relative to some observations
  !!                given zobs (observed field), zref (reference run field)
  !!                and ztst (test run field)
  !!               compute zchk as the ratio : zchk=(zref - ztst) / (zref - zobs )
  !!               Where  0 < zchk <=1  the correction act in the right direction
  !!               Where  1 < zchk      the correction is too strong, in the right way
  !!               Where  zchk < 0      the correction is in the wrong way (deterioration)
  !!               or deterioration (-1)
  !!               store results on file
  !!  
  !!  **  Method: Try to avoid 3 d arrays 
  !!
  !! history: 
  !!     Original :   J.M. Molines (Nov. 2005)
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk                                  !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: nvpk                                !: dim of the working variable
  INTEGER, DIMENSION(1) ::  ipk, &                 !: outptut variables : number of levels,
       &                    id_varout              !: ncdf varid's
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: zobs, zref, ztst ,& !: Array to read a layer of data
       &                                         zmask ,&         !: 2D mask at surface
       &                                         zchk             !: check index output
  REAL(KIND=4),DIMENSION(1)                   ::  tim

  CHARACTER(LEN=80) :: cfilobs, cfilref, cfiltst ,cvar ,cfileout='chk.nc' !:
  TYPE (variable), DIMENSION(1) :: typvar         !: structure for attributes

  INTEGER    :: ncout
  INTEGER    :: istatus

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfimprovechk  cdfvar obs.nc ref.nc  tst.nc '
     PRINT *,' Output on chk.nc, same variable '
     STOP
  ENDIF

  CALL getarg (1, cvar)
  CALL getarg (2, cfilobs)
  CALL getarg (3, cfilref)
  CALL getarg (4, cfiltst)

  npiglo= getdim (cfilref,'x')
  npjglo= getdim (cfilref,'y')
  npk   = getdim (cfilref,'depth')

  nvpk  = getvdim (cfilref,cvar)
  IF (nvpk == 2 ) nvpk = 1
  IF (nvpk == 3 ) nvpk = npk

  ipk(:)= nvpk  ! all variables 
  typvar(1)%name=TRIM(cvar)
  typvar(1)%units='%'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= 0.
  typvar(1)%valid_max= 100.
  typvar(1)%long_name='Checking ratio for'//TRIM(cvar)
  typvar(1)%short_name=cvar
  typvar(1)%online_operation='N/A'
  IF (nvpk == npk ) typvar(1)%axis='TZYX'
  IF (nvpk == 1   ) typvar(1)%axis='TYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE (zobs(npiglo,npjglo), zref(npiglo,npjglo), ztst(npiglo,npjglo) ,zmask(npiglo,npjglo))
  ALLOCATE (zchk(npiglo,npjglo) )

  ! create output fileset

  ncout =create(cfileout, cfilref, npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfilref,npiglo, npjglo,npk)

  zref = 0.
  zobs = 0.
  zmask = 1.
  DO jk = 1, npk
     PRINT *,'level ',jk
     zchk = 0.
     zobs(:,:) = getvar(cfilobs, cvar,  jk ,npiglo, npjglo)
     zref(:,:) = getvar(cfilref, cvar,  jk ,npiglo, npjglo)
     ztst(:,:) = getvar(cfiltst, cvar,  jk ,npiglo, npjglo)
     IF (jk == 1  )  THEN
        tim=getvar1d(cfilref,'time_counter',1)
        WHERE( zref == 0. ) zmask=0.
     END IF
       WHERE ( (zref -zobs) /= 0 ) 
         zchk= (zref - ztst ) / ( zref - zobs) * zmask
       END WHERE
     ierr = putvar(ncout, id_varout(1) ,zchk, jk,npiglo, npjglo)

  ENDDO

        ierr=putvar1d(ncout,tim,1,'T')

  istatus = closeout(ncout)
END PROGRAM cdfimprovechk
