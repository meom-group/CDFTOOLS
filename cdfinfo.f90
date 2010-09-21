PROGRAM cdfinfo
  !!-----------------------------------------------------------------------
  !!                 ***  PROGRAM cdfinfo  ***
  !!
  !!  **  Purpose: Give very basic informations for Netcdf File
  !!  
  !!  **  Method: 
  !!
  !! history :
  !!     Original code :   J.M. Molines (Sep. 2010) 
  !!-----------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !!
  USE cdfio 

  IMPLICIT NONE
  INTEGER   :: jk,jt,jvar, jv , jtt,jkk                     !: dummy loop index
  INTEGER   :: ierr                                         !: working integer
  INTEGER   :: narg, iargc                                  !: 
  INTEGER   :: npiglo,npjglo, npk ,nt                       !: size of the domain
  INTEGER   :: nvars                                        !: Number of variables in a file
  INTEGER   :: ntframe                                      !: Cumul of time frame
  INTEGER , DIMENSION(:), ALLOCATABLE :: id_var , &         !: arrays of var id's
       &                             ipk    , &             !: arrays of vertical level for each var
       &                             id_varout

  CHARACTER(LEN=256) :: cfile                                !: file name
  CHARACTER(LEN=256) ::  cdep
  CHARACTER(LEN=256) ,DIMENSION(:), ALLOCATABLE:: cvarname   !: array of var name
  TYPE (variable), DIMENSION(:), ALLOCATABLE :: typvar
  
  INTEGER    :: istatus


  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfinfo ''model cdf file'' '
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cfile)

  npiglo= getdim (cfile,'x')
  npjglo= getdim (cfile,'y')
  npk   = getdim (cfile,'depth',cdtrue=cdep, kstatus=istatus)

  IF (istatus /= 0 ) THEN
     npk   = getdim (cfile,'z',cdtrue=cdep,kstatus=istatus)
     IF (istatus /= 0 ) THEN
       npk   = getdim (cfile,'sigma',cdtrue=cdep,kstatus=istatus)
        IF ( istatus /= 0 ) THEN 
          npk = getdim (cfile,'nav_lev',cdtrue=cdep,kstatus=istatus)
            IF ( istatus /= 0 ) THEN 
              npk = getdim (cfile,'levels',cdtrue=cdep,kstatus=istatus)
              IF ( istatus /= 0 ) THEN 
                PRINT *,' assume file with no depth'
                npk=0
              ENDIF
            ENDIF
        ENDIF
     ENDIF
  ENDIF
  

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk

  PRINT *,' Depth dimension name is ', TRIM(cdep)

  nvars = getnvar(cfile)
  PRINT *,' nvars =', nvars

  ALLOCATE (cvarname(nvars)  )
  ALLOCATE (typvar(nvars)  )

  ! get list of variable names and collect attributes in typvar (optional)
  cvarname(:)=getvarname(cfile,nvars,typvar)

  DO jvar = 1, nvars
   PRINT *, 'variable# ',jvar,' is : ',TRIM(cvarname(jvar))
  END DO

END PROGRAM cdfinfo
