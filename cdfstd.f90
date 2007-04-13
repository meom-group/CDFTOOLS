PROGRAM cdfstd
  !!-----------------------------------------------------------------------
  !!                 ***  PROGRAM cdfrms  ***
  !!
  !!  **  Purpose: Compute Standard deviation values for all the variables in a bunch
  !!                of cdf files given as argument
  !!                Store the results on a 'similar' cdf file.
  !!  
  !!  **  Method: Try to avoid 3 d arrays 
  !!
  !! history :
  !!     Original code :   F. Castruccio (2.0, from cdfmoy) 04/2007
  !!                       J.M. Molines for 2.1 (04/07)
  !!-----------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !!
  USE cdfio 

  IMPLICIT NONE
  INTEGER   :: jk,jt,jtt,jvar, jv                           !: dummy loop index
  INTEGER   :: ierr                                         !: working integer
  INTEGER   :: narg, iargc                                  !: 
  INTEGER   :: npiglo,npjglo, npk,nt                        !: size of the domain
  INTEGER   :: nvars                                        !: Number of variables in a file
  INTEGER   :: ntframe                                      !: Cumul of time frame
  INTEGER, DIMENSION(:), ALLOCATABLE :: id_var , &          !: arrays of var id's
       &                             ipk    , &             !: arrays of vertical level for each var
       &                             id_varout
  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE :: tab, tab2  !: Arrays for cumulated values
  REAL(KIND=8)                                :: total_time
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: v2d ,&       !: Array to read a layer of data
       &                                   rmean, rmean2, std
  REAL(KIND=4), DIMENSION(1)                   :: timean
  REAL(KIND=4), DIMENSION(2000)                :: tim

  CHARACTER(LEN=80) :: cfile ,cfileout                       !: file name
  CHARACTER(LEN=80) :: cdep
  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE:: cvarname    !: array of var name
  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE:: cvarnameo   !: array of var name for output

  TYPE ( variable ), DIMENSION(:), ALLOCATABLE :: typvar, typvaro

  INTEGER    :: ncout
  INTEGER    :: istatus
  LOGICAL    :: lcaltmean

  !!

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfstd ''list_of_ioipsl_model_output_files'' '
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cfile)

  npiglo= getdim (cfile,'x')
  npjglo= getdim (cfile,'y')
  npk   = getdim (cfile,'depth',cdtrue=cdep, kstatus=istatus)

  IF (istatus /= 0 ) THEN
     npk   = getdim (cfile,'z',kstatus=istatus)
     IF (istatus /= 0 ) THEN
        npk   = getdim (cfile,'sigma',cdtrue=cdep,kstatus=istatus)
     ELSE
        PRINT *,' assume file with no depth'
        npk=0
     ENDIF
  ENDIF

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( tab(npiglo,npjglo), tab2(npiglo,npjglo), v2d(npiglo,npjglo) )
  ALLOCATE( rmean(npiglo,npjglo), rmean2(npiglo,npjglo), std(npiglo,npjglo) )

  nvars = getnvar(cfile)
  PRINT *,' nvars =', nvars

  ALLOCATE (cvarname(nvars), cvarnameo(nvars) )
  ALLOCATE (typvar(nvars), typvaro(nvars) )
  ALLOCATE (id_var(nvars),ipk(nvars),id_varout(nvars) )

  cvarname(:)=getvarname(cfile,nvars,typvar)

  id_var(:)  = (/(jv, jv=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cfile,nvars,cdep=cdep)
  DO jvar = 1, nvars 
     cvarnameo(jvar)=TRIM(cvarname(jvar))//'_std' 
  ENDDO

  WHERE( ipk == 0 ) cvarnameo='none'

  DO jvar = 1, nvars
     typvaro(jvar)=typvar(jvar)
     typvaro(jvar)%name=cvarnameo(jvar)
     typvaro(jvar)%long_name='Std Deviation of '//TRIM(cvarname(jvar))
     typvaro(jvar)%short_name=cvarnameo(jvar)
  END DO

  ! create output fileset
  cfileout='cdfstd.nc'
  ! create output file taking the sizes in cfile

  ncout =create(cfileout, cfile, npiglo, npjglo, npk,cdep=cdep )
  ierr= createvar(ncout, typvaro,  nvars, ipk, id_varout )
  ierr= putheadervar(ncout, cfile, npiglo, npjglo, npk,cdep=cdep )

  lcaltmean=.TRUE.
  DO jvar = 1,nvars
     IF (cvarname(jvar) == 'nav_lon' .OR. &
          cvarname(jvar) == 'nav_lat' .OR. &
          cvarnameo(jvar) == 'none'  ) THEN
        ! skip these variable
     ELSE
        PRINT *,' Working with ', TRIM(cvarname(jvar)), ipk(jvar)
        DO jk = 1, ipk(jvar)
           PRINT *,'level ',jk
           tab(:,:) = 0.d0 ; tab2(:,:) = 0.d0 ; total_time = 0.; ntframe=0
           DO jt = 1, narg
              CALL getarg (jt, cfile)
              nt = getdim (cfile,'time_counter')
              IF ( lcaltmean )  THEN
                 tim(1:nt)=getvar1d(cfile,'time_counter',nt)
                 total_time = total_time + SUM(tim(1:nt) )
              END IF
              DO jtt=1,nt
                 ntframe=ntframe+1
                 v2d(:,:)= getvar(cfile, cvarname(jvar), jk ,npiglo, npjglo ,ktime=jtt)
                 tab(:,:) = tab(:,:) + v2d(:,:)
                 tab2(:,:) = tab2(:,:) + v2d(:,:)*v2d(:,:)
              END DO
           END DO
           ! finish with level jk ; compute mean (assume spval is 0 )
           rmean(:,:) = tab(:,:)/ntframe
           rmean2(:,:) = tab2(:,:)/ntframe
           std(:,:) = SQRT(rmean2(:,:) - (rmean(:,:)*rmean(:,:)))
           ! store variable on outputfile
           ierr = putvar(ncout, id_varout(jvar) ,std, jk, npiglo, npjglo)
           IF (lcaltmean )  THEN
              timean(1)= total_time/ntframe
              ierr=putvar1d(ncout,timean,1,'T')
           END IF
           lcaltmean=.FALSE. ! tmean already computed
        END DO  ! loop to next level
     END IF

  END DO ! loop to next var in file

  istatus = closeout(ncout)

END PROGRAM cdfstd
