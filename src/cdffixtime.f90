PROGRAM cdffixtime
  !!======================================================================
  !!                     ***  PROGRAM  cdffixtime  ***
  !!=====================================================================
  !!  ** Purpose : Correct time inconsistency in model output file or
  !!               mean fields.
  !!
  !!  ** Method  : Adjust the values of time_counters in order to be 
  !!               coherent with the time_origin and units attribute.
  !!               According to drakkar the time in seconds represents
  !!               the time of the model at the moment of output, ie at
  !!               the end of the averaging period. The time origin is
  !!               shifted back half the averaging period in order to
  !!               indicate the center of the averaging period.
  !!                  This program is intended to manage both leap year
  !!               and noleap year calendars.
  !!
  !! History :  3.0  : 04/2011  : J.M. Molines : Doctor norm + Lic.
  !!         :  4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   jcnes         : return the jcnes Julian day from time tag
  !!   julday        : return the true Julian day 
  !!   caldatjm      : Return the calendar date from the input jcnes day
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class data_transformation
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)            :: narg            ! number of arguments
  INTEGER(KIND=4)            :: iargc           ! f90 function
  INTEGER(KIND=4)            :: ijarg           ! argument counter
  INTEGER(KIND=4)            :: is, ie          ! starting and ending position of the tag in file name
  INTEGER(KIND=4)            :: iyr_init        ! initial date (year)
  INTEGER(KIND=4)            :: imm_init        ! initial date (month)
  INTEGER(KIND=4)            :: idd_init        ! initial date (day)
  INTEGER(KIND=4)            :: ihr_init        ! ititial time (hour)
  INTEGER(KIND=4)            :: imn_init        ! ititial time (minutes)
  INTEGER(KIND=4)            :: isec_init       ! ititial time (seconds)
  INTEGER(KIND=4)            :: ierr            ! error status for i/o

  REAL(KIND=4)               :: rpp_one_year = 365 ! 365.2425
  REAL(KIND=4)               :: rdt_obs = 5.    ! time interval between file fields (days)
  REAL(KIND=4)               :: rday0           ! CNES julian day corresponding to tag of initial date
  REAL(KIND=4)               :: rday_origin     ! CNES julian day corresponding to origin date
  REAL(KIND=4), DIMENSION(1) :: rdaycnes        ! CNES julian day corresponding to current tag

  REAL(KIND=8), DIMENSION(1) :: dseconds        ! seconds since rday0

  CHARACTER(LEN=80)          :: cf_in           ! input file 
  CHARACTER(LEN=80)          :: cldum           ! dummy character variable
  CHARACTER(LEN=80)          :: ctag='none'     ! tag default. Interpreted from file name if possible
  CHARACTER(LEN=80)          :: cldate, ctim    ! date and time as string
  CHARACTER(LEN=80)          :: ctag0           ! time tag from input initial date/time
  CHARACTER(LEN=80)          :: ctim_unit       ! attribute value for time_counter unit
  CHARACTER(LEN=80)          :: ctim_origin     ! attribute value for time_counter time_origin
  CHARACTER(LEN=3)           :: cmm             ! month in character

  LOGICAL                    :: lnoleap=.TRUE.  ! flag for noleap years
  LOGICAL                    :: lagrif=.FALSE.  ! flag for agrif files
  LOGICAL                    :: lkeep=.FALSE.  ! flag for agrif files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdffixtime  -f IN-file -i initial date [-t tag] [-dt freq] ... '
     PRINT *,'               ...  [-keep ] [-leap] [-noleap]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Change time_counter in file to set it according to drakkar rule,' 
     PRINT *,'        time_counter attibutes ''units'' and ''time_origin'' are ajusted.'
     PRINT *,'          * units are ''seconds since yyyy-mm-dd hh:mm:ss'' '
     PRINT *,'          * time_origin is set to ''yyyy-MMM-dd hh:mm:ss'', MMM represents a'
     PRINT *,'            literal abbreviation for the month (eg: JAN FEB MAR ...)'
     PRINT *,'        Once fixed, the time_counter indicates the middle of the output '
     PRINT *,'        interval (in case of averaged output, of course).'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file     : specify the file whose time_counter need adjustment' 
     PRINT *,'       -i inital date : indicate the time origin in a fixed 2 words format'
     PRINT *,'                   yyyy-mm-dd hh:mm:ss ( eg: 1956-05-16 04:30:00 )'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-t tag]  : supply a time tag corresponding to the file. If not'
     PRINT *,'                   supplied, tag is taken from the name of the input file'
     PRINT *,'                   assuming DRAKKAR convention ( CONFIG-CASE_tag_xxxx.nc )'
     PRINT *,'       [-dt freq]: number of days between model output [ 5d ]'
     PRINT *,'       [-leap]   : assume a calendar with leap years'
     PRINT *,'       [-noleap] : assume a calendar without leap years (default)'
     PRINT *,'       [-keep]   : keep the actual value of time_counter, adjust time_counter'
     PRINT *,'                   attributes only;'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none ' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : Input file is modified (only attributes)'
     PRINT *,'      '
     STOP 
  ENDIF

  ! parse command line 
  ijarg=1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum) ; ijarg=ijarg + 1 
     SELECT CASE (cldum)
     CASE ( '-f'      ) ; CALL getarg(ijarg,cf_in ) ; ijarg=ijarg +1
     CASE ( '-t'      ) ; CALL getarg(ijarg,ctag  ) ; ijarg=ijarg +1
     CASE ( '-dt'     ) ; CALL getarg(ijarg,cldum ) ; ijarg=ijarg +1 ; READ(cldum,*) rdt_obs
     CASE ( '-i'      ) ; CALL getarg(ijarg,cldate) ; ijarg=ijarg +1
                        ; CALL getarg(ijarg,ctim  ) ; ijarg=ijarg +1
     CASE ( '-leap'   ) ; rpp_one_year = 365.2425
                        ; lnoleap      = .FALSE.
     CASE ( '-noleap' ) ; rpp_one_year = 365
                        ; lnoleap      = .TRUE.
     CASE ( '-keep'   ) ; lkeep        = .TRUE.
     CASE DEFAULT       ; PRINT *,' ERROR : ',TRIM(cldum),' :  unknown options.' ; STOP 99
     END SELECT
  END DO

  IF ( chkfile(cf_in) ) STOP 99 ! missing file
  PRINT *,' Changing time on file :', TRIM(cf_in)

  ! if ctag = none, try to find it from the file name.
  IF ( TRIM(ctag) == 'none' ) THEN  ! no tag given as arguments
     is = INDEX(cf_in,'_')
     IF ( is == 2 ) THEN 
        PRINT *,' ASSUME AGRIF file for ', TRIM(cf_in)
        lagrif = .TRUE.
     ENDIF

     IF (lagrif) THEN
        is=INDEX(cf_in(3:),'_' )+2
        ie=INDEX(cf_in(is+1:),'_' )
        ctag=cf_in(is+1:is+ie-1) 
     ELSE
        is=INDEX(cf_in,'_')
        ie=INDEX(cf_in(is+1:),'_' )
        ctag=cf_in(is+1:is+ie-1) 
     ENDIF

     is=INDEX(ctag,'d')
     IF ( is == 0 ) THEN ! not a model output but a mean value
        is=INDEX(ctag,'m') 
        IF ( is == 0 ) THEN ! annual mean set pseudo date to 01/07
           ctag=ctag(1:5)//"m07d01"
        ELSE  ! monthly mean set pseudo date to the 15 of month
           ctag=ctag(1:8)//"d15"
        ENDIF
     ENDIF
  ENDIF

  PRINT *,'            Using tag = ', TRIM(ctag)

  ! interpret ctim and cldate
  READ(cldate,'(i4,1x,i2,1x,i2)'            ) iyr_init, imm_init, idd_init
  READ(ctim,  '(i2,1x,i2,1x,i2)'            ) ihr_init, imn_init, isec_init
  WRITE(ctag0,'("y",i4.4,"m",i2.2,"d",i2.2)') iyr_init, imm_init, idd_init

  ! jcnes of initial date including time as fraction of days
  rday0 = jcnes(ctag0) + ihr_init/24.0 + imn_init/60./24. + isec_init/3600./24.

  ! compute the pseudo time_origin and set up variable attributes
  rday_origin = rday0 - rdt_obs/2.   ! offset of -1/2 of time interval
  CALL caldatjm(rday_origin, iyr_init, imm_init, idd_init, ihr_init, imn_init, isec_init)

  WRITE(cldate,'(i4.4,"-",i2.2,"-",i2.2)') iyr_init, imm_init, idd_init
  WRITE(ctim,  '(i2.2,":",i2.2,":",i2.2)') ihr_init, imn_init, isec_init

  ! Compute initial julian day
  SELECT CASE ( imm_init )
  CASE (  1 ) ; cmm='JAN' 
  CASE (  2 ) ; cmm='FEB'
  CASE (  3 ) ; cmm='MAR'
  CASE (  4 ) ; cmm='APR'
  CASE (  5 ) ; cmm='MAY'
  CASE (  6 ) ; cmm='JUN'
  CASE (  7 ) ; cmm='JUL'
  CASE (  8 ) ; cmm='AUG'
  CASE (  9 ) ; cmm='SEP'
  CASE ( 10 ) ; cmm='OCT'
  CASE ( 11 ) ; cmm='NOV'
  CASE ( 12 ) ; cmm='DEC'
  END SELECT

  WRITE(ctim_unit,  '("seconds since ",a,i3.2,":",i2.2,":",i2.2   )') TRIM(cldate), ihr_init, imn_init, isec_init
  WRITE(ctim_origin,'(i5,"-",a,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)') iyr_init,cmm, idd_init, ihr_init, imn_init, isec_init

  PRINT *, "          ",TRIM(cn_vtimec)," units set to : ", TRIM(ctim_unit)
  PRINT *, "          ",TRIM(cn_vtimec)," time origin set to : ", TRIM(ctim_origin)

  ! Compute corresponding jcnes
  rdaycnes=jcnes(ctag)
  dseconds=(rdaycnes - rday0 +1 ) * 86400.d0 

  ! Modify cdfile !! CAUTION : Original file will be modified  !!
  IF ( .NOT. lkeep )  THEN
     ierr = putvar1d( cf_in, cn_vtimec, dseconds, 1 )
  ENDIF
  ierr = atted   ( cf_in, cn_vtimec, 'units',       ctim_unit  )
  ierr = atted   ( cf_in, cn_vtimec, 'time_origin', ctim_origin)

CONTAINS

  REAL(KIND=4) FUNCTION jcnes(cdtag)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION jcnes  ***
    !!
    !! ** Purpose : return the JCNES corresponding to time tag. JCNES is a julian
    !!              day refered from 1950-01-01
    !!
    !! ** Method  : Interface with function julday 
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),INTENT(in) :: cdtag

    INTEGER(KIND=4)             :: iyear, imon, iday
    REAL(KIND=4)                :: zsec = 0.
    REAL(KIND=4)                :: zjuldeb, zjulfin, zjulday

    READ(cdtag,'(1x,i4.4,1x,i2.2,1x,i2.2)') iyear, imon, iday
    zsec=0.
    !---------------------------------------------------------------------
    zjulfin = julday(iyear, imon, iday, zsec)
    zjuldeb = julday(1950, 01, 01, 0.)
    jcnes   = zjulfin - zjuldeb
  END FUNCTION jcnes

  REAL(KIND=4) FUNCTION julday(kyear, kmonth, kday, psec)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION julday  ***
    !!
    !! ** Purpose : Converts year, month, day and seconds into a julian day
    !!
    !! ** Method  : In 1968 in a letter to the editor of Communications of 
    !!              the ACM (CACM, volume 11, number 10, October 1968, p.657)
    !!              Henry F. Fliegel and Thomas C. Van Flandern presented 
    !!              such an algorithm.
    !!                  In the case of the Gregorian calendar we have chosen 
    !!              to use the Lilian day numbers. This is the day counter 
    !!              which starts on the 15th October 1582.
    !!              This is the day at which Pope Gregory XIII introduced the
    !!              Gregorian calendar.
    !!              Compared to the true Julian calendar, which starts some
    !!              7980 years ago, the Lilian days are smaller and are dealt 
    !!              with easily on 32 bit machines. With the true Julian days
    !!              you can only the fraction of the day in the real part to 
    !!              a precision of a 1/4 of a day with 32 bits.
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) :: kyear, kmonth, kday  ! input date
    REAL(KIND=4),    INTENT(in) :: psec                 ! input seconds

    REAL(KIND=4),     PARAMETER :: pp_one_day = 86400.0
    INTEGER(KIND=4),  PARAMETER :: jp_mon_len(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)

    INTEGER(KIND=4)             :: in_m, in_y, in_d
    INTEGER(KIND=4)             :: ijd, iml
    !---------------------------------------------------------------------

    in_m = kmonth
    in_y = kyear
    in_d = kday
    !- We deduce the calendar from the length of the year as it
    !- is faster than an INDEX on the calendar variable.
    !-
    !- Gregorian
    IF ( (rpp_one_year > 365.0) .AND. (rpp_one_year < 366.0) ) THEN
       ijd = (1461*(in_y+4800+INT(( in_m-14 )/12)))/4 &
            &      +(367*(in_m-2-12*(INT(( in_m-14 )/12))))/12 &
            &      -(3*((in_y+4900+INT((in_m-14)/12))/100))/4 &
            &      +in_d-32075
       ijd = ijd-2299160
       !- No leap or All leap
    ELSE IF (ABS(rpp_one_year-365.0) <= EPSILON(rpp_one_year) .OR. &
         &   ABS(rpp_one_year-366.0) <= EPSILON(rpp_one_year)) THEN
       iml = SUM(jp_mon_len(1:in_m-1))
       ijd = in_y*INT(rpp_one_year)+iml+(in_d-1)
       !- Calendar with regular month
       !  ELSE
       !    iml = INT(one_year)/12
       !    ijd = y*INT(one_year)+(m-1)*iml+(d-1)
    ENDIF
    !-
    julday = ijd + psec / pp_one_day
  END FUNCTION julday


  SUBROUTINE caldatjm( pjcnes, ky, km, kd, kh, kmn, ksec )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE caldatjm  ***
    !!
    !! ** Purpose : Compute the calendar date from the julian CNES day
    !!              given as input  
    !!
    !! ** Method  : Take care of the leap/noleap calendar. That's why we
    !!              cannot use the standard caldat from numerical recipe
    !!----------------------------------------------------------------------
    REAL(KIND=4),       INTENT(in) :: pjcnes
    INTEGER(KIND=4),   INTENT(out) :: ky, km, kd, kh, kmn, ksec

    INTEGER(KIND=4)                :: jd, jm      ! dummy loop index
    INTEGER(KIND=4)                :: isec, idays
    INTEGER(KIND=4), DIMENSION(12) :: indays=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    INTEGER(KIND=4), DIMENSION(12) :: icumul
    !!--------------------------------------------------------------------------------
    ! initialize the cumulated time
    icumul(1)    = indays(1)
    DO jm=2,12
       icumul(jm) = icumul(jm-1) + indays(jm)
    ENDDO

    ! look for time part of pjcnes
    isec = (pjcnes-INT(pjcnes) ) * 86400.
    kh   = isec/3600
    kmn  = (isec - kh * 3600 )/60
    ksec =  isec - kh * 3600 - kmn * 60

    ! number of years since 1950
    IF ( lnoleap ) THEN ! no leap years
       ky=1950 + INT(pjcnes)/365
       idays= ( INT(pjcnes)/ 365. - INT(pjcnes)/365 )* 365 
       km=1 ; kd=0
       DO jd=1, idays
          IF ( jd > icumul(km) ) THEN
             km=km+1
             kd=1
          ELSE
             kd=kd+1
          ENDIF
       ENDDO
    ELSE
       ! use caldat from Numerical Recipe
       CALL caldat_nr ( pjcnes, ky, km, kd, kh, kmn, ksec )
    ENDIF
  END SUBROUTINE caldatjm

  SUBROUTINE caldat_nr( pjcnes, kiyyy, kmm, kid, kh, kmn, ksec ) 
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE caldat_nr  ***
    !!
    !! ** Purpose : This routine convert a julian day in calendar date.
    !!
    !! ** Method  :  This routine comes directly from the Numerical Recipe Book,
    !!
    !!   Arguments
    !!     kjulian : input julian day number
    !!     kmm     : output, corresponding month
    !!     kid     : output, corresponding day
    !!     kiyyy   : output, corresponding year, positive IF a.d, negative b.c.
    !!
    !! References  : Numerical Recipe Book,  Press et al., numerical recipes,
    !!               cambridge univ. press, 1986.
    !!----------------------------------------------------------------------
    IMPLICIT NONE

    REAL(KIND=4),    INTENT(in)  :: pjcnes
    INTEGER(KIND=4), INTENT(out) :: kiyyy, kmm, kid, kh, kmn, ksec
    ! * Local
    INTEGER(KIND=4), PARAMETER  :: jpgreg = 2299161
    INTEGER(KIND=4)             :: ijulian
    INTEGER(KIND=4)             :: ia, ialpha, ib, ic, id, ie, isec
    REAL(KIND=4)                :: zjul1950
    !!----------------------------------------------------------------------
    ! look for time part of pjcnes
    isec = (pjcnes-INT(pjcnes) ) * 86400.
    kh   = isec/3600
    kmn  = (isec - kh * 3600 )/60
    ksec =  isec - kh * 3600 - kmn * 60
    zjul1950 = julday_nr( 01, 01, 1950)
    ijulian  = INT(pjcnes + zjul1950)
    !
    IF ( ijulian >= jpgreg) THEN
       ialpha = INT ((( ijulian - 1867216) - 0.25)/36524.25 )
       ia     = ijulian +1 + ialpha -INT (0.25*ialpha)
    ELSE
       ia = ijulian
    END IF
    !
    ib = ia + 1524
    ic = INT (6680. + (( ib -2439870) - 122.1)/365.25 )
    id = 365* ic + INT (0.25*ic)
    ie = INT (( ib - id )/30.6001)
    !
    kid = ib - id - INT (30.6001*ie)
    kmm = ie -1
    IF ( kmm > 12 ) kmm = kmm - 12
    kiyyy = ic - 4715
    IF ( kmm   >  2 ) kiyyy = kiyyy - 1
    IF ( kiyyy <= 0 ) kiyyy = kiyyy - 1
  END SUBROUTINE caldat_nr

  INTEGER(KIND=4) FUNCTION julday_nr(kmm,kid,kiyyy)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION julday_nr  ***
    !!
    !! ** Purpose : his routine returns the julian day number which begins at noon
    !!         of the calendar date specified by month kmm, day kid, and year kiyyy.
    !!         positive year signifies a.d.; negative, b.c.  (remember that the
    !!         year after 1 b.c. was 1 a.d.)
    !!         routine handles changeover to gregorian calendar on oct. 15, 1582.
    !!
    !! ** Method:  This routine comes directly from the Numerical Recipe Book,
    !!
    !!----------------------------------------------------------------------
    INTEGER, INTENT(in) :: kiyyy
    INTEGER, INTENT(in) :: kmm, kid
    !  * Local
    INTEGER, PARAMETER ::jpgreg=15+31*(10+12*1582)
    INTEGER  ::ky, iy, im, ia
    !!----------------------------------------------------------------------
    ky = kiyyy
    ! ... Year 0 never existed ...
    IF (ky == 0) STOP 99
    !
    IF (ky < 0) ky = ky + 1
    IF (kmm > 2) THEN
       iy = ky
       im = kmm + 1
    ELSE
       iy = ky - 1
       im = kmm + 13
    END IF
    !
    julday_nr = INT(365.25*iy) + INT(30.6001*im) + kid + 1720995
    IF (kid+31*(kmm+12*ky).GE.jpgreg) THEN
       ia = INT(0.01*iy)
       julday_nr = julday_nr + 2 - ia + INT(0.25*ia)
    END IF
  END FUNCTION julday_nr

END PROGRAM cdffixtime
