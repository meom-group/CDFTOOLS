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
  !!                : 05/2015  : J.M. Molines : Rewrite to be compliant with XIOS
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class time_averaging
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jvar    ! dummy loop index
  INTEGER(KIND=4)                               :: jv, jtt     ! dummy loop index
  INTEGER(KIND=4)                               :: jframe      ! dummy loop index
  INTEGER(KIND=4)                               :: it1, it2    ! box limits
  INTEGER(KIND=4)                               :: ierr        ! working integer
  INTEGER(KIND=4)                               :: narg, iargc     ! 
  INTEGER(KIND=4)                               :: ijarg, ijm
  INTEGER(KIND=4)                               :: npiglo, npjglo  ! size of the domain
  INTEGER(KIND=4)                               :: npk ,npt        ! size of the domain
  INTEGER(KIND=4)                               :: ip, nf          ! working integer
  INTEGER(KIND=4)                               :: nvars           ! Number of variables in a file
  INTEGER(KIND=4)                               :: nframes         ! Number of frames in the output file
  INTEGER(KIND=4)                               :: ndyr, nhyr      ! Days and hours per year
  INTEGER(KIND=4)                               :: nhfri           ! input freq in hours
  INTEGER(KIND=4)                               :: ncout, ncout2
  INTEGER(KIND=4), DIMENSION( 12)               :: njm             ! month vector
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var, ipk, id_varout, ibox

  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d, rmean
  REAL(KIND=4), DIMENSION(:,:,:),   ALLOCATABLE :: v3d
  REAL(KIND=4), DIMENSION(:,:,:,:), ALLOCATABLE :: v4d

  REAL(KIND=8)                                  :: dtotal_time
  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dtim, dtimean
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dtab         ! Arrays for cumulated values

  CHARACTER(LEN=256)                            :: cf_in               !
  CHARACTER(LEN=256)                            :: cf_out='cdfmoy_'    ! file name
  CHARACTER(LEN=256)                            :: cv_dep
  CHARACTER(LEN=256)                            :: cfreq_o             ! output frequency
  CHARACTER(LEN=256)                            :: cldum               ! dummy character arguments
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names            ! array of var nam
  CHARACTER(LEN=2  ), DIMENSION(4)              :: cfreq_a=(/'h ','d ','mo','y '/) ! authorized keys for frequency specif
  CHARACTER(LEN=2  )                            :: cfr_id              ! current output freq id

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar

  LOGICAL                                       :: lcaltmean=.TRUE.
  LOGICAL                                       :: lleap
  LOGICAL                                       :: lerr
  LOGICAL                                       :: lnc4 = .FALSE.
  LOGICAL                                       :: lv3d = .FALSE.
  LOGICAL                                       :: lv4d = .FALSE.
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmoy_freq -f IN-file -avg AVG-length [-v3d] [-v4d] [-o OUT-rootname]'
     PRINT *,'            ... [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       This program takes a file covering 1 year of data (evenly spaced) and'
     PRINT *,'       sub-samples the data by performing box averages, which span is given as'
     PRINT *,'       argument. The original data sampling can be hours, days, monthes or'
     PRINT *,'       even seasons.'
     PRINT *,'       The program recognizes leap years, and when feb. 29 is found, it is '
     PRINT *,'       included in the current ''box'' (averaging length is thus increased by'
     PRINT *,'       1 day.)'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file : gives the name of the yearly file containing either 365 '
     PRINT *,'              or 366 days'
     PRINT *,'       -avg AVG-length : Set the time size of the averaging box. Averaging '
     PRINT *,'              length is specified using XIOS convention (e.g. 1d,5d, 1mo, 1y ;'
     PRINT *,'              4mo stands for seasonal means )'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-v3d] : use 3d variable (x,y,t) : save execution time, increase memory'
     PRINT *,'       [-v4d] : use 4d variable (x,y,z,t): save execution time, increase memory'
     PRINT *,'       [-o OUT-rootname] : specify the root of the output file name instead '
     PRINT *,'                   of ',TRIM(cf_out),'. Final name will have <freq> appened'
     PRINT *,'                   to the root.'
     PRINT *,'       [-nc4] : use netcdf4 with chunking and deflation for the output file'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file :  cdfmoy_output<freq>.nc'
     PRINT *,'         variables :  same as variables in input file.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      cdfmoy, cdfmoy_weighted'
     PRINT *,'      '
     STOP 
  ENDIF

  ! parse command line
  ijarg = 1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum) ; ijarg = ijarg +1
     SELECT CASE ( cldum ) 
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_in   ) ; ijarg = ijarg + 1
     CASE ( '-avg' ) ; CALL getarg(ijarg, cfreq_o ) ; ijarg = ijarg + 1
        ! options
     CASE ( '-v3d' ) ; lv3d=.TRUE.
     CASE ( '-v4d' ) ; lv4d=.TRUE.
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out  ) ; ijarg = ijarg + 1
     CASE ( '-nc4' ) ; lnc4=.TRUE.
     CASE DEFAULT    ; PRINT *,' ERROR : ', TRIM(cldum) ,' : unknown option.' ; STOP 99
     END SELECT
  END DO

  IF ( chkfile ( cf_in ) ) STOP 99 ! missing file

  ! parse the cfreqo to determine the output frequency. 
  ! Allowed syntax is nf<cfr_id> where nf is an integer >0, <cfr_id> is h, d, mo or y
  cfr_id='--'
  DO jk =1, 4
     ip=INDEX(cfreq_o, cfreq_a(jk) )
     IF ( ip /= 0 ) THEN
        cfr_id=cfreq_o(ip:)
        READ(cfreq_o(1:ip-1), *) nf
        EXIT
     ENDIF
  ENDDO

  SELECT CASE ( cfr_id )
  CASE ( '--' ) 
     PRINT *, ' +++ ERROR : Cannot determine the output frequency'
     PRINT *, '            You should use a character string such as 6h, 5d, 1mo, 1y '
  CASE ( 'd' ) 
     IF ( nf /= 1 .AND. nf/=5 ) THEN
        PRINT *, ' +++ ERROR : only 1d or 5d are acceptable !'
        STOP 99
     ENDIF
  CASE ( 'y' )
     IF ( nf > 1 ) THEN
        PRINT *, ' +++ ERROR : Cannot have output freq > 1 y !'
        STOP 99
     ENDIF
  END SELECT

  ! get domain size from the input file
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

  ! Now look at input file and try to look for input frequency
  ! Note that we know that the file contains either 365d or 366 days of data
  ! We suppose that input file freq is a multiple of 1 hour.
  lleap = .FALSE.
  ndyr  = 365     ! number of days per year
  nhyr  = ndyr*24 ! number of hours per year

  IF ( MOD( nhyr, npt ) /= 0 ) THEN
     ndyr = 366     ! try leap year
     nhyr = ndyr*24 ! number of hours per leap year
     IF ( MOD( nhyr, npt ) /= 0 ) THEN
        PRINT *," +++ ERROR : npt do not fit in 365 nor 366 days "
        STOP 99
     ELSE
        lleap=.TRUE.
     ENDIF
  ENDIF
  nhfri = 24*ndyr/npt  ! input frequency in hours

  PRINT *, 'INPUT FILE : ', TRIM(cf_in)
  PRINT *, 'NPIGLO = ', npiglo
  PRINT *, 'NPJGLO = ', npjglo
  PRINT *, 'NPK    = ', npk
  PRINT *, 'NPT    = ', npt
  PRINT *, ' '
  PRINT *, ' LEAP YEAR  : ', lleap
  PRINT *, ' INPUT FREQ : ', nhfri,' hours '

  ! Now determines the number of frames in the output file and detect impossible case
  !  Also determines the number of input frames in boxes. (can be variable)
  ! number of day by month
  njm(:)= (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  IF ( lleap ) THEN 
     njm(:) = (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  ENDIF

  lerr = .TRUE.
  SELECT CASE ( cfr_id )
  CASE ( 'h' ) 
     IF ( MOD( nf, nhfri ) == 0 ) THEN
        nframes = ndyr*24/nf 
        ALLOCATE (ibox( nframes) )
        ibox(:)=nf/nhfri ; lerr=.FALSE.
     ENDIF
  CASE ( 'd' ) ! note : 365 = 73 *5  ( 366 = 73 *5 + 1 ) 
     IF ( MOD ( nf*24, nhfri ) == 0 ) THEN
        IF ( nf == 1 ) nframes = ndyr
        IF ( nf == 5 ) nframes = 73  ! in case of leap year frame#12 will have 6day average
        ALLOCATE (ibox( nframes) )
        lerr=.FALSE.
        ibox(:)=nf*24/nhfri 
        IF ( lleap .AND. nf == 5 ) ibox(12)=6*24/nhfri
     ENDIF
  CASE ( 'mo' )
     !##################################
     !JM : do not work if nhfri > 24 !!! 
     !##################################
     IF ( MOD( 12, nf ) == 0 ) THEN
        nframes=12/nf ; lerr=.FALSE.
        ALLOCATE (ibox( nframes) )  ! 12 6 4 3 2 
        SELECT CASE ( nframes )
        CASE ( 12 )
           ibox(:) = njm(:)*24/nhfri 
        CASE ( 6 )
           DO jframe= 1, nframes
              ijm=jframe*2-1
              ibox(jframe) = (njm(ijm)+njm(ijm+1)) *24/nhfri 
           ENDDO
        CASE ( 4 )
           DO jframe= 1, nframes
              ijm=jframe*3-2
              ibox(jframe) = (njm(ijm)+njm(ijm+1)+njm(ijm+2)) *24/nhfri 
           ENDDO
        CASE ( 3 )
           DO jframe= 1, nframes
              ijm=jframe*4-3
              ibox(jframe) = (njm(ijm)+njm(ijm+1)+njm(ijm+2)+njm(ijm+3)) *24/nhfri 
           ENDDO
        CASE ( 2 )
           DO jframe= 1, nframes
              ijm=jframe*6-5
              ibox(jframe) = (njm(ijm)+njm(ijm+1)+njm(ijm+2)+njm(ijm+3)+njm(ijm+4)+njm(ijm+5)) *24/nhfri 
           ENDDO
        END SELECT
     ENDIF
  CASE ( 'y' )  ! 1y average  all is to be taken !
     nframes = 1
     ALLOCATE (ibox( nframes) )
     ibox(:) = npt
     lerr = .FALSE.
  END SELECT

  IF ( lerr ) THEN
     PRINT *, ' +++ ERROR : Input and output frequency incompatible.'
     PRINT *, '         Input  : ',  nhfri,' hours '
     PRINT *, '         Output : ',  nf,' hours '
     STOP 99
  ENDIF

  ALLOCATE( dtab(npiglo,npjglo), v2d(npiglo,npjglo) )
  ALLOCATE( rmean(npiglo,npjglo)                    )
  IF (lv3d)   ALLOCATE( v3d(npiglo,npjglo,npt)      )

  nvars = getnvar(cf_in)
  PRINT *,' nvars =', nvars

  ALLOCATE (cv_names(nvars)                            )
  ALLOCATE (stypvar(nvars)                             )
  ALLOCATE (id_var(nvars), ipk(nvars), id_varout(nvars))
  ALLOCATE( dtim( npt), dtimean(nframes)             )

  ! get list of variable names and collect attributes in stypvar (optional)
  cv_names(:)=getvarname(cf_in, nvars, stypvar)

  CALL CreateOutput

  dtim(:)=getvar1d(cf_in, cn_vtimec, npt)
  DO jvar = 1,nvars
     IF ( cv_names(jvar) == cn_vlon2d .OR.                            &
          cv_names(jvar) == cn_vlat2d .OR. cv_names(jvar) == 'none') THEN
        ! skip these variable
     ELSE
        PRINT *,' Working with ', TRIM(cv_names(jvar))
        IF (lv4d)   ALLOCATE( v4d(npiglo, npjglo,ipk(jvar),npt)  )

        IF (lv4d)  v4d(:,:,:,:) = getvar4d(cf_in, cv_names(jvar),npiglo, npjglo,ipk(jvar),npt)
        DO jk=1,ipk(jvar)

           ! initialisation
           dtab(:,:) = 0.d0 ; dtotal_time = 0.d0

           ! time loop
           it1=1
           IF ( lv3d )  v3d(:,:,:)=getvar3dt(cf_in, cv_names(jvar),jk,npiglo, npjglo, npt)
           DO jframe = 1, nframes
              it2=it1+ibox(jframe)-1
              DO jtt=it1, it2
                 ! load data
                  IF ( lv4d) THEN 
                       v2d(:,:) = v4d(:,:,jk,jtt)
                  ELSE
                    IF ( lv3d ) THEN ; v2d(:,:)  = v3d(:,:,jtt)
                    ELSE             ; v2d(:,:)  = getvar(cf_in, cv_names(jvar), jk, npiglo, npjglo, ktime=jtt )
                    ENDIF
                  ENDIF
                 dtab(:,:) = dtab(:,:) + v2d(:,:)*1.d0
                 IF ( lcaltmean ) THEN
                    dtotal_time = dtotal_time + dtim(jtt) 
                 ENDIF
              ENDDO
              rmean(:,:) = dtab(:,:)/ibox(jframe)
              ierr = putvar(ncout, id_varout(jvar) ,rmean, jk, npiglo, npjglo, ktime=jframe)
              IF ( lcaltmean ) THEN
                 dtimean(jframe) = dtotal_time/ibox(jframe)
              ENDIF
              dtab(:,:) = 0.d0 ; dtotal_time = 0.d0
              it1 = it2 + 1
              !
           ENDDO ! loop to next time
           lcaltmean=.FALSE.

        ENDDO ! loop to next level
        IF (lv4d)   DEALLOCATE( v4d )
     END IF
  END DO ! loop to next var in file

  ierr = putvar1d(ncout,   dtimean,  nframes  , 'T')
  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
  id_var(:)  = (/(jv, jv=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cf_in, nvars, cdep=cv_dep)
  !
  WHERE( ipk == 0 ) cv_names='none'
  stypvar(:)%cname = cv_names
  DO jv=1, nvars
    stypvar(jv)%ichunk=(/npiglo,MAX(1,npjglo/30), 1, 1 /)
  ENDDO

  ! create output file taking the sizes in cf_in
  cf_out = TRIM(cf_out)//'_'//TRIM(cfreq_o)//'.nc'
  ncout = create      (cf_out, cf_in,   npiglo, npjglo, npk, cdep=cv_dep , ld_nc4=lnc4 )
  ierr  = createvar   (ncout,  stypvar, nvars,  ipk,    id_varout        , ld_nc4=lnc4 )
  ierr  = putheadervar(ncout,  cf_in,   npiglo, npjglo, npk, cdep=cv_dep               )


  END SUBROUTINE CreateOutput

END PROGRAM cdfmoy_freq
