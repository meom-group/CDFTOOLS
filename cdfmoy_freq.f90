PROGRAM cdfmoy_freq
  !!======================================================================
  !!                     ***  PROGRAM  cdfmoy_freq  ***
  !!=====================================================================
  !!  ** Purpose : Mainly in case of forcing file (gathered as yearly file)
  !!               compute annual mean, monthl mean or diurnal means.
  !!
  !!  ** Method  : Detect the frequency of the input file according to the
  !!               number of fields in the file.
  !!
  !! History : 2.1  : 06/2007  : P. Mathiot   : Original code from cdfmoy
  !!           3.0  : 06/2011  : J.M. Molines : Doctor norm + Lic.
  !!           3.0  : 10/2011  : P. Mathiot   : Add seasonal option and 
  !!                                            allow file with 73 time steps
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: nt_in, nt_out
  INTEGER(KIND=4)                               :: jk, jvar    ! dummy loop index
  INTEGER(KIND=4)                               :: jv, jtt     ! dummy loop index
  INTEGER(KIND=4)                               :: ierr            ! working integer
  INTEGER(KIND=4)                               :: itime    ! dummy loop index
  INTEGER(KIND=4)                               :: narg, iargc     ! 
  INTEGER(KIND=4)                               :: npiglo, npjglo  ! size of the domain
  INTEGER(KIND=4)                               :: npk ,npt        ! size of the domain
  INTEGER(KIND=4)                               :: nvars           ! Number of variables in a file
  INTEGER(KIND=4)                               :: ntframe         ! Cumul of time frame
  INTEGER(KIND=4)                               :: ncout, ncout2
  INTEGER(KIND=4), DIMENSION(365)               :: njd  ! day vector
  INTEGER(KIND=4), DIMENSION( 12)               :: njm  ! month vector
  INTEGER(KIND=4), DIMENSION(  4)               :: njs  !season vector
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var, ipk, id_varout

  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d, rmean
  REAL(KIND=4), DIMENSION(1)                    :: time
  REAL(KIND=4), DIMENSION(365)                  :: tim

  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dtab         ! Arrays for cumulated values
  REAL(KIND=8)                                  :: dtotal_time

  CHARACTER(LEN=256)                            :: cf_in               !
  CHARACTER(LEN=256)                            :: cf_out              ! file name
  CHARACTER(LEN=256)                            :: cv_dep
  CHARACTER(LEN=256)                            :: cfreq_out, cfreq_in
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names            ! array of var nam

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar

  LOGICAL                                       :: lcaltmean, lprint
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmoy_freq IN-file output_frequency'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute annual mean or monthly mean or daily mean from a yearly'
     PRINT *,'       input forcing file given on input.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       IN-file : netcdf input file corresponding to 1 year of forcing variable (1460, 365, 73 or 12 time step) '
     PRINT *,'       output_frequency : either one of monthly, daily, seasonal or annual.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file :  cdfmoy_outputFreaquency.nc'
     PRINT *,'         variables :  same as variables in input file.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      cdfmoy, cdfmoy_weighted'
     PRINT *,'      '
     STOP
  ENDIF

  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cf_in    )
  CALL getarg (2, cfreq_out)

  IF ( chkfile ( cf_in ) ) STOP ! missing file

  SELECT CASE ( cfreq_out )
  CASE ('daily'   ) ; nt_out = 365
  CASE ('monthly' ) ; nt_out =  12
  CASE ('seasonal') ; nt_out =   4
  CASE ('annual'  ) ; nt_out =   1
  CASE DEFAULT 
     PRINT *, 'Pb : this output_frequency is not allowed, please use daily, monthly, seasonal or annual'
     STOP
  END SELECT

  npiglo= getdim (cf_in, cn_x                             )
  npjglo= getdim (cf_in, cn_y                             )
  npk   = getdim (cf_in, cn_z, cdtrue=cv_dep, kstatus=ierr)

  IF (ierr /= 0 ) THEN
     npk   = getdim (cf_in,'z',cdtrue=cv_dep,kstatus=ierr)
     IF (ierr /= 0 ) THEN
        npk   = getdim (cf_in,'sigma',cdtrue=cv_dep,kstatus=ierr)
        IF ( ierr /= 0 ) THEN 
           PRINT *,' assume file with no depth'
           npk=0
        ENDIF
     ENDIF
  ENDIF

  npt   = getdim (cf_in, cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE( dtab(npiglo,npjglo), v2d(npiglo,npjglo) )
  ALLOCATE( rmean(npiglo,npjglo)                    )

  nvars = getnvar(cf_in)
  PRINT *,' nvars =', nvars

  ALLOCATE (cv_names(nvars)                            )
  ALLOCATE (stypvar(nvars)                             )
  ALLOCATE (id_var(nvars), ipk(nvars), id_varout(nvars))

  ! get list of variable names and collect attributes in stypvar (optional)
  cv_names(:)=getvarname(cf_in, nvars, stypvar)

  id_var(:)  = (/(jv, jv=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cf_in, nvars, cdep=cv_dep)
  !
  WHERE( ipk == 0 ) cv_names='none'
  stypvar(:)%cname = cv_names

  !
  ! create output file taking the sizes in cf_in

  cf_out = 'cdfmoy_'//TRIM(cfreq_out)//'.nc'
  ncout = create      (cf_out, cf_in,   npiglo, npjglo, npk, cdep=cv_dep )
  ierr  = createvar   (ncout,  stypvar, nvars,  ipk,    id_varout )
  ierr  = putheadervar(ncout,  cf_in,   npiglo, npjglo, npk, cdep=cv_dep )

  time=getvar1d(cf_in, cn_vtimec, 1)
  ierr=putvar1d(ncout, time, 1, 'T')

  npt   = getdim (cf_in, cn_t)
  nt_in = npt

  SELECT CASE ( npt )
  CASE ( 1460 ) ; PRINT *, 'Frequency of this file : 6h '
  CASE (  365 ) ; PRINT *, 'Frequency of this file : daily '
  CASE (   73 ) ; PRINT *, 'Frequency of this file : 5 day '
  CASE (   12 ) ; PRINT *, 'Frequency of this file : monthly '
  END SELECT

  IF (npt <= nt_out) THEN
     PRINT *, 'You don''t need to use it, or it is impossible (npt_in <= npt_out)'
     STOP
  END IF

  ! number of day by month, season and day
  njd(:)= 1
  njm= (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  njs= (/ 90, 91, 92, 92 /)

  DO jvar = 1,nvars
     IF ( cv_names(jvar) == cn_vlon2d .OR.                            &
          cv_names(jvar) == cn_vlat2d .OR. cv_names(jvar) == 'none') THEN
        ! skip these variable
     ELSE
        PRINT *,' Working with ', TRIM(cv_names(jvar))
	DO jk=1,ipk(jvar)

    ! initialisation
           dtab(:,:) = 0.d0 ; dtotal_time = 0.d0;  ntframe=0; itime=1; 

	   ! time loop
           DO jtt=1, nt_in
	      lprint=.FALSE.
              ntframe=ntframe+1
              ! load data
              v2d(:,:)  = getvar(cf_in, cv_names(jvar), jk, npiglo, npjglo, ktime=jtt )
              dtab(:,:) = dtab(:,:) + v2d(:,:)*1.d0

              ! detection of time when you have to print the average
              IF (nt_out==365) THEN
	         IF (MOD(jtt,FLOOR(SUM(njd(1:itime)) * nt_in/365.))==0) lprint=.TRUE.
              ELSE IF (nt_out==12) THEN
	         IF (MOD(jtt,FLOOR(SUM(njm(1:itime)) * nt_in/365.))==0) lprint=.TRUE.
              ELSE IF (nt_out==4) THEN   
	         IF (MOD(jtt,FLOOR(SUM(njs(1:itime)) * nt_in/365.))==0) lprint=.TRUE.
              ELSE IF (nt_out==1) THEN
	         IF (MOD(jtt,                          nt_in      )==0) lprint=.TRUE.
              END IF
              !
              ! Compute and print the average at the right time              
              IF ( lprint ) THEN
                 IF (jk==1) PRINT *, itime, jtt,'/',npt 
                 ! compute mean
                 rmean(:,:) = dtab(:,:)/ntframe
                 ! store variable on outputfile
                 ierr = putvar(ncout, id_varout(jvar) ,rmean, jk, npiglo, npjglo, ktime=itime)
                 dtab(:,:) = 0.d0 ; dtotal_time = 0.;  ntframe=0; itime=itime+1
              END IF
	      !
           ENDDO ! loop to next time

	ENDDO ! loop to next level

     END IF
  END DO ! loop to next var in file

  ierr = closeout(ncout)


END PROGRAM cdfmoy_freq
