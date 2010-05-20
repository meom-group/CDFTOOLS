PROGRAM cdfmoy_weighted
  !!-----------------------------------------------------------------------
  !!                 ***  PROGRAM cdfmoy_weighted  ***
  !!
  !!  **  Purpose: Compute weighted mean values from monthly mean
  !!  
  !!  **  Method: monthly mean were computed (cdfmoy) with all dumps that fall within a montn
  !!              thus, all month have different weigth : Feb = 5. March, Dec. = 7 other = 6
  !!
  !! history :
  !!     Original code :   J.M. Molines (Nov 2004 ) for ORCA025
  !!                       J.M. Molines (Apr 2005 ) put all NCF stuff in module
  !!                              now valid for grid T U V W icemod
  !!-----------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !!
  USE cdfio 

  IMPLICIT NONE
  INTEGER   :: jk,jt,jvar, jv                               !: dummy loop index
  INTEGER   :: ierr,idum                                    !: working integer
  INTEGER   :: narg, iargc                                  !: 
  INTEGER   :: npiglo,npjglo, npk                           !: size of the domain
  INTEGER   ::  nvars                                       !: Number of variables in a file
  INTEGER , DIMENSION(:), ALLOCATABLE :: id_var , &         !: arrays of var id's
       &                             ipk    , &         !: arrays of vertical level for each var
       &                             id_varout 
  INTEGER, DIMENSION(:), ALLOCATABLE :: iweight
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: tab  !: Arrays for cumulated values
  REAL(KIND=8)                                :: total_time, sumw
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: v2d ,&       !: Array to read a layer of data
       &                                   rmean
  REAL(KIND=4),DIMENSION(1)                   :: timean, tim

  CHARACTER(LEN=256) :: cfile ,cfileout           !: file name
  CHARACTER(LEN=256) ::  cdep
  CHARACTER(LEN=256) ,DIMENSION(:), ALLOCATABLE:: cvarname   !: array of var name
  CHARACTER(LEN=256) ,DIMENSION(:), ALLOCATABLE:: cdummy     !: array of var name
  
  TYPE (variable), DIMENSION(:), ALLOCATABLE :: typvar, typvardum

  INTEGER    :: ncout
  INTEGER    :: istatus

  !!

  !!  Read command line
  narg= iargc()
  IF ( narg ==  0 ) THEN
     PRINT *,' Usage : cdfmoy_weighted ''list of files'' '
     STOP
  ENDIF
  ALLOCATE (iweight(narg) )  ! as maby weights as files
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

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( tab(npiglo,npjglo), v2d(npiglo,npjglo) )
  ALLOCATE( rmean(npiglo,npjglo) )

  nvars = getnvar(cfile)
  PRINT *,' nvars =', nvars

  ALLOCATE (cvarname(nvars),cdummy(nvars) )
  ALLOCATE (typvar(nvars), typvardum(nvars)  )
  ALLOCATE (id_var(nvars),ipk(nvars),id_varout(nvars) )

  ! get list of variable names and collect attributes in typvar (optional)
  cvarname(:)=getvarname(cfile,nvars,typvar)


  id_var(:)  = (/(jv, jv=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cfile,nvars,cdep=cdep)
  WHERE( ipk == 0 ) cvarname='none'
  typvar(:)%name=cvarname

  ! create output fileset
  cfileout='cdfmoy_weighted.nc'
  ! create output file taking the sizes in cfile

  ncout =create(cfileout, cfile,npiglo,npjglo,npk,cdep=cdep)

  ierr= createvar(ncout , typvar,  nvars, ipk, id_varout )

  ierr= putheadervar(ncout , cfile, npiglo, npjglo, npk,cdep=cdep)

  DO jvar = 1,nvars
     ! fill iweight for each variables: need to scan all the input files
     DO jt=1,narg  ! this is far from optimal : think about a special function
                   ! for retrieving an attribute of a variable
       CALL getarg(jt,cfile)
       cdummy(:)=getvarname(cfile,nvars,typvardum)
       iweight(jt)=typvardum(jvar)%iwght
     ENDDO
     PRINT *, iweight
     IF (cvarname(jvar) == 'nav_lon' .OR. &
          cvarname(jvar) == 'nav_lat' ) THEN
        ! skip these variable
     ELSE
        PRINT *,' Working with ', TRIM(cvarname(jvar)), ipk(jvar)
        DO jk = 1, ipk(jvar)
           PRINT *,'level ',jk
           tab(:,:) = 0.d0 ; total_time = 0. ; sumw=0.
           DO jt = 1, narg
              sumw = sumw + iweight(jt)
              IF (jk == 1 .AND. jvar == nvars )  THEN
                 tim=getvar1d(cfile,'time_counter',1)
                 total_time = total_time + tim(1)
              END IF
              CALL getarg (jt, cfile)
              v2d(:,:)= getvar(cfile, cvarname(jvar), jk ,npiglo, npjglo )
              tab(:,:) = tab(:,:) + iweight(jt)* v2d(:,:)
           END DO
           ! finish with level jk ; compute mean (assume spval is 0 )
           rmean(:,:) = tab(:,:)/sumw
           ! store variable on outputfile
           ierr = putvar(ncout, id_varout(jvar) ,rmean, jk, npiglo, npjglo,kwght=INT(sumw) )
           IF (jk == 1 .AND. jvar == nvars )  THEN
              timean(1)= total_time/narg
              ierr=putvar1d(ncout,timean,1,'T')
           END IF
        END DO  ! loop to next level
     END IF
  END DO ! loop to next var in file

  istatus = closeout(ncout)


END PROGRAM cdfmoy_weighted
