PROGRAM cdfnan
  !!-----------------------------------------------------------------------
  !!                 ***  PROGRAM cdfnan  ***
  !!
  !!  **  Purpose: Replace the nan values by spval or another value given in argument
  !!  
  !! history :
  !!   Original :  J.M. Molines from cdfcsp
  !!-----------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !!
  !! * Modules used
  USE cdfio 

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jf,jk,jvar, jt, jkk                          !: dummy loop index
  INTEGER   :: ierr                                         !: working integer
  INTEGER   :: ncid, narg, iargc                            !: 
  INTEGER   :: npiglo,npjglo, npk , nt                      !: size of the domain
  INTEGER   ::  nvars                                       !: Number of variables in a file
  INTEGER , DIMENSION(:), ALLOCATABLE :: ipk                !: arrays of vertical level for each var
  INTEGER , DIMENSION(:), ALLOCATABLE :: id_var             !: arrays of var id
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: tab        !: Arrays for cumulated values
  REAL(KIND=4)                                :: spval, replace
  CHARACTER(LEN=256) :: cfile                                !: file name
  CHARACTER(LEN=256) :: cunits, clname, csname
  CHARACTER(LEN=256) ,DIMENSION(:), ALLOCATABLE:: cvarname   !: array of var name

  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvar       !: type for attributes
  LOGICAL     :: l_replace=.false.

  INTEGER    :: istatus


  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfnan ''list_of_ioipsl_model_output_files''[-value replace] '
     PRINT *,'  Option -value must be the last on the line'
     PRINT *,'  When used, it replace Nan by this value instead of the spval'
     STOP
  ENDIF
  PRINT *, 'narg=', narg
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  IF ( narg >= 2 ) THEN
    CALL getarg ( narg - 1, cfile)
    IF (TRIM(cfile) == '-value' ) THEN
      CALL getarg(narg,cfile) ; READ(cfile,*) replace ; l_replace=.true.
    ENDIF
  ENDIF
  CALL getarg (1, cfile)

  npiglo= getdim (cfile,'x')
  npjglo= getdim (cfile,'y')
  npk   = getdim (cfile,'depth',kstatus=istatus)
  nt    = getdim (cfile,'time_counter')
  IF (istatus /= 0 ) THEN
     npk   = getdim (cfile,'z',kstatus=istatus)
     IF (istatus /= 0 ) THEN
       PRINT *, "ASSUME NO VERTICAL DIMENSIONS !"
       npk=0
     ENDIF
  ENDIF

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( tab(npiglo,npjglo) )

  nvars = getnvar(cfile)

  ALLOCATE (cvarname(nvars), id_var(nvars),ipk(nvars), typvar(nvars))

  print *,' in getvarname'
  cvarname(:)=getvarname(cfile,nvars,typvar)
  print *,' in getipk'
  ipk(:)      = getipk(cfile,nvars)
  print *,' done'
  id_var(:)    = getvarid(cfile,nvars)

  DO jf = 1, narg
     CALL getarg (jf, cfile)
     PRINT *, 'Change spval on file ', cfile
     ncid = ncopen(cfile)
     nt    = getdim (cfile,'time_counter')
     DO jvar = 1,nvars
        IF (cvarname(jvar) == 'nav_lon' .OR. &
             cvarname(jvar) == 'nav_lat'  .OR. &
             cvarname(jvar) == 'time_counter'  .OR. &
             cvarname(jvar) == 'deptht'  .OR.  &
             cvarname(jvar) == 'depthu'  .OR.  &
             cvarname(jvar) == 'depthv' )  THEN
           ! skip these variable
        ELSE
           IF ( l_replace ) THEN
              spval=replace
           ELSE
              ierr = getvaratt (cfile,cvarname(jvar),cunits,spval,clname,csname)
           ENDIF

          DO jt=1,nt
           DO jk = 1, ipk(jvar) 
              jkk=jk
              IF (npk == 0 ) jkk=jt
              tab(:,:) = getvar(cfile, cvarname(jvar), jkk ,npiglo, npjglo, ktime=jt )
              WHERE( isnan(tab(:,:)) ) tab(:,:) = spval
              ierr = putvar(ncid, id_var(jvar) ,tab, jkk, npiglo, npjglo, ktime=jt)
           ENDDO
          END DO
        ENDIF 
     ENDDO
  ENDDO

  istatus = closeout(ncid)

END PROGRAM cdfnan
