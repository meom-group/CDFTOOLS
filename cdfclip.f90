PROGRAM cdfclip
  !!-----------------------------------------------------------------------
  !!                 ***  PROGRAM cdfclip  ***
  !!
  !!  **  Purpose: same functionality than ncks but without changing the order of the dim/variables
  !!  
  !!  **  Method: read zoomed area on the command line ( imin imax jmin jmax)
  !!              read the sub zome
  !!              write the subzone
  !!
  !! history :
  !!     Original code :   J.M. Molines (Feb 2007)
  !!-----------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !!
  USE cdfio 

  IMPLICIT NONE
  INTEGER   :: jk,jt,jvar, jv                               !: dummy loop index
  INTEGER   :: ierr                                         !: working integer
  INTEGER   :: imin, imax, jmin, jmax
  INTEGER   :: narg, iargc , jarg                           !: 
  INTEGER   :: npiglo,npjglo, npk                           !: size of the domain
  INTEGER   ::  nvars                                       !: Number of variables in a file
  INTEGER , DIMENSION(:), ALLOCATABLE :: id_var , &         !: arrays of var id's
       &                             ipk    , &         !: arrays of vertical level for each var
       &                             id_varout , ndim
  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE :: tab  !: Arrays for cumulated values
  REAL(KIND=8)                               :: total_time
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: v2d ,rlon, rlat
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE    :: dep
  REAL(KIND=4), DIMENSION(1)                 :: timean, tim

  CHARACTER(LEN=80) :: cfile ,cfileout                         !: file name
  CHARACTER(LEN=80) ::  cdep, cdum
  CHARACTER(LEN=80) ,DIMENSION(:), ALLOCATABLE:: cvarname   !: array of var name
  
  TYPE (variable), DIMENSION(:), ALLOCATABLE :: typvar

  INTEGER    :: ncout
  INTEGER    :: istatus
  !!

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfclip -f file -zoom imin imax jmin jmax '
     STOP
  ENDIF
  !!
  jarg=1
  DO WHILE (jarg < 7 )
    CALL getarg (jarg, cdum)
    SELECT CASE ( cdum)
    CASE ('-f' )
       jarg=jarg+1 ; CALL getarg(jarg,cdum) ; cfile=cdum
    CASE ('-zoom')
       jarg=jarg+1 ; CALL getarg(jarg,cdum) ; READ(cdum,*) imin
       jarg=jarg+1 ; CALL getarg(jarg,cdum) ; READ(cdum,*) imax
       jarg=jarg+1 ; CALL getarg(jarg,cdum) ; READ(cdum,*) jmin
       jarg=jarg+1 ; CALL getarg(jarg,cdum) ; READ(cdum,*) jmax
    CASE DEFAULT
       PRINT *,' Unknown option :', TRIM(cdum) ; STOP
    END SELECT
    jarg=jarg+1
  ENDDO
   

! npiglo= getdim (cfile,'x')
! npjglo= getdim (cfile,'y')
  npiglo= imax-imin+1
  npjglo= jmax-jmin+1
  npk   = getdim (cfile,'depth',cdtrue=cdep, kstatus=istatus)

  IF (istatus /= 0 ) THEN
     npk   = getdim (cfile,'z',cdtrue=cdep,kstatus=istatus)
     IF (istatus /= 0 ) THEN
       npk   = getdim (cfile,'sigma',cdtrue=cdep,kstatus=istatus)
        ELSE
          IF (istatus /= 0 ) THEN
! STOP 'depth dimension name not suported'
          PRINT *,' assume file with no depth'
          npk=0
        ENDIF
     ENDIF
  ENDIF
  

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( v2d(npiglo,npjglo),rlon(npiglo,npjglo), rlat(npiglo,npjglo), dep(npk) )

  nvars = getnvar(cfile)
  PRINT *,' nvars =', nvars

  ALLOCATE (cvarname(nvars),ndim(nvars) )
  ALLOCATE (typvar(nvars))
  ALLOCATE (id_var(nvars),ipk(nvars),id_varout(nvars))

  rlon=getvar(cfile,'nav_lon',1,npiglo,npjglo,kimin=imin,kjmin=jmin)
  rlat=getvar(cfile,'nav_lat',1,npiglo,npjglo,kimin=imin,kjmin=jmin)
  dep=getvar1d(cfile,cdep,npk)

  ! get list of variable names and collect attributes in typvar (optional)
  cvarname(:)=getvarname(cfile,nvars,typvar)

  ! save variable dimension in ndim
  !  1 = either time or depth : noclip
  !  2 = nav_lon, nav_lat
  !  3 = X,Y,T  or X,Y,Z   <-- need to fix the ambiguity
  !  4 = X,Y,Z,T
  DO jvar=1,nvars
      ndim(jvar) = getvdim(cfile,cvarname(jvar)) + 1   !  we add 1 because vdim is dim - 1 ...
  END DO

  id_var(:)  = (/(jv, jv=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  IF ( npk /= 0 ) THEN
  ipk(:)     = getipk (cfile,nvars,cdep=cdep)
  WHERE( ipk == 0 ) cvarname='none'
  ENDIF
  typvar(:)%name=cvarname

  ! create output fileset
  cfileout='cdfclip.nc'
  ! create output file taking the sizes in cfile

  ncout =create(cfileout, cfile,npiglo,npjglo,npk,cdep=cdep)
  ierr= createvar(ncout , typvar,  nvars, ipk, id_varout )
  ierr= putheadervar(ncout , cfile, npiglo, npjglo, npk,pnavlon=rlon, pnavlat=rlat,pdep=dep,cdep=cdep)

  DO jvar = 1,nvars
      ! skip dimension variables (already done when creating the output file)
     SELECT CASE (cvarname(jvar) )
      CASE ('none' )
          ! skip
      CASE DEFAULT
            DO jk=1,ipk(jvar)
             v2d=getvar(cfile,cvarname(jvar),jk,npiglo,npjglo,kimin=imin,kjmin=jmin)
             ierr=putvar(ncout,id_varout(jvar),v2d,jk,npiglo,npjglo)
            END DO
      END SELECT
  END DO ! loop to next var in file
  ierr=putvar1d(ncout,timean,1,'T')

  istatus = closeout(ncout)


END PROGRAM cdfclip
