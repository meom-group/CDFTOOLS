PROGRAM cdfmoy_freq
  !!-----------------------------------------------------------------------
  !!                 ***  PROGRAM cdfmoy  ***
  !!
  !!  **  Purpose: Compute mean values for all the variables in a bunch
  !!                of cdf files given as argument
  !!                Store the results on a 'similar' cdf file.
  !!  
  !!  **  Method: Try to avoid 3 d arrays 
  !!
  !! history :
  !!     Original code :   J.M. Molines (Nov 2004 ) for ORCA025
  !!                       J.M. Molines (Apr 2005 ) put all NCF stuff in module
  !!                              now valid for grid T U V W icemod
  !!     Modified      :   P. Mathiot (June 2007) update for forcing fields
  !!-----------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !!
  USE cdfio 

  IMPLICIT NONE
  INTEGER   :: nt_in, nt_out, nmois
  INTEGER   :: jk,jt,jvar, jv , jtt,jkk                     !: dummy loop index
  INTEGER   :: ierr                                         !: working integer
  INTEGER   :: narg, iargc                                  !: 
  INTEGER   :: npiglo,npjglo, npk ,nt                       !: size of the domain
  INTEGER   :: nvars                                        !: Number of variables in a file
  INTEGER   :: ntframe                                      !: Cumul of time frame
  INTEGER, DIMENSION(12) :: njm

  INTEGER , DIMENSION(:), ALLOCATABLE :: id_var , &         !: arrays of var id's
       &                             ipk    , &         !: arrays of vertical level for each var
       &                             id_varout
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: tab         !: Arrays for cumulated values
  REAL(KIND=8)                                :: total_time
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: v2d ,&       !: Array to read a layer of data
       &                                   rmean
  REAL(KIND=4),DIMENSION(1)                   :: time
  REAL(KIND=4),DIMENSION(365)                   ::  tim

  CHARACTER(LEN=256) :: cfile ,cfileout                      !: file name
  CHARACTER(LEN=256) ::  cdep, cfreq_out, cfreq_in
  CHARACTER(LEN=256) ,DIMENSION(:), ALLOCATABLE:: cvarname   !: array of var nam
  
  TYPE (variable), DIMENSION(:), ALLOCATABLE :: typvar

  INTEGER    :: ncout, ncout2
  INTEGER    :: istatus
  LOGICAL    :: lcaltmean

  !!

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmoy_freq forcing_field frequency (monthly or daily or annual)'
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cfile)
  CALL getarg (2, cfreq_out)

  IF (TRIM(cfreq_out) .EQ. 'daily')  nt_out=365 ! 
  IF (TRIM(cfreq_out) .EQ. 'monthly') nt_out=12  ! 
  IF (TRIM(cfreq_out) .EQ. 'annual') nt_out=1   ! 
  IF ((TRIM(cfreq_out) .NE. 'annual') .AND. (TRIM(cfreq_out) .NE. 'daily') .AND. (TRIM(cfreq_out) .NE. 'monthly')) THEN
     PRINT *, 'Pb : this frequency is not allowed, used please daily, monthly or annual'
     STOP
  END IF
  

  npiglo= getdim (cfile,'x')
  npjglo= getdim (cfile,'y')
  npk   = getdim (cfile,'depth',cdtrue=cdep, kstatus=istatus)

  IF (istatus /= 0 ) THEN
     npk   = getdim (cfile,'z',cdtrue=cdep,kstatus=istatus)
     IF (istatus /= 0 ) THEN
       npk   = getdim (cfile,'sigma',cdtrue=cdep,kstatus=istatus)
        IF ( istatus /= 0 ) THEN 
          PRINT *,' assume file with no depth'
          npk=0
        ENDIF
     ENDIF
  ENDIF
  

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( tab(npiglo,npjglo), v2d(npiglo,npjglo) )
  ALLOCATE( rmean(npiglo,npjglo))

  nvars = getnvar(cfile)
  PRINT *,' nvars =', nvars

  ALLOCATE (cvarname(nvars))
  ALLOCATE (typvar(nvars))
  ALLOCATE (id_var(nvars),ipk(nvars),id_varout(nvars))

  ! get list of variable names and collect attributes in typvar (optional)
  cvarname(:)=getvarname(cfile,nvars,typvar)

  id_var(:)  = (/(jv, jv=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cfile,nvars,cdep=cdep)
  !
  WHERE( ipk == 0 ) cvarname='none'
  typvar(:)%name=cvarname

  PRINT *, '',cvarname

  ! create output fileset
  cfileout='cdfmoy_'//TRIM(cfreq_out)//'.nc'
  ! create output file taking the sizes in cfile

  ncout =create(cfileout, cfile,npiglo,npjglo, 0)
  ierr= createvar(ncout , typvar,  nvars, ipk, id_varout )
  ierr= putheadervar(ncout , cfile, npiglo, npjglo, 0)
  time=getvar1d(cfile,'time_counter',1)
  ierr=putvar1d(ncout,time,1,'T')


  nt = getdim (cfile,'time_counter')
  nt_in=nt
  IF (nt==1460) THEN
     PRINT *, 'Frequency of this file : 6h '
!     cfreq_in='6h'
!     nfreq_in=6
  END IF
  IF (nt==365) THEN
     PRINT *, 'Frequency of this file : daily '
!     cfreq_in='daily'
!     nfreq_in=24
  END IF
  IF (nt==12) THEN
     PRINT *, 'Frequency of this file : monthly '
!     cfreq_in='monthly'
!     nfreq_in=720
  END IF

  IF (nt .LE. nt_out) THEN
     PRINT *, 'You don''t need to use it, or it is impossible'
     STOP
  END IF
  jt=0
  njm= (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  DO jvar = 1,nvars
     IF (cvarname(jvar) == 'nav_lon' .OR. &
          cvarname(jvar) == 'nav_lat' .OR. cvarname(jvar) == 'none') THEN
        ! skip these variable
     ELSE
        PRINT *,' Working with ', TRIM(cvarname(jvar))
        tab(:,:) = 0.d0 ; total_time = 0.;  ntframe=0; nmois=1
        DO jtt=1,nt_in
           ntframe=ntframe+1
           ! If forcing fields is without depth dimension
           v2d(:,:)= getvar(cfile, cvarname(jvar), jtt ,npiglo, npjglo,ktime=jtt )
           tab(:,:) = tab(:,:) + v2d(:,:)
           !PRINT *, '',v2d(100,100), tab(100,100), ntframe
           IF (nt_out==12) THEN
              IF (ntframe .EQ. njm(nmois)*nt_in/365) THEN
                 PRINT *, nmois, jtt,'/',nt
                 jt=jt+1
                 ! finish with level jk ; compute mean (assume spval is 0 )
                 rmean(:,:) = tab(:,:)/ntframe
                 ! store variable on outputfile
                 ierr = putvar(ncout, id_varout(jvar) ,rmean, jt, npiglo, npjglo, jt)
                 tab(:,:) = 0.d0 ; total_time = 0.;  ntframe=0; nmois=nmois+1
              END IF
           ELSE
              !PRINT *, jtt,'/',nt,' et on enregistre tous les ',nt_in/nt_out
              IF (MOD(jtt,nt_in/nt_out)==0) THEN
                 jt=jt+1
                 PRINT *, jtt,'/',nt,' dumping every ',nt_in/nt_out
                 ! finish with level jk ; compute mean (assume spval is 0 )
                 rmean(:,:) = tab(:,:)/ntframe
                 PRINT *, '',rmean(100,100)
                 ! store variable on outputfile
                 ierr = putvar(ncout, id_varout(jvar) ,rmean, jt, npiglo, npjglo, jt)
                 tab(:,:) = 0.d0 ; total_time = 0.;  ntframe=0
              END IF
           END IF
        ENDDO
     END IF
  END DO ! loop to next var in file

  istatus = closeout(ncout)


END PROGRAM cdfmoy_freq
