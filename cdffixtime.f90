PROGRAM cdffixtime
  !--------------------------------------------------------------------------------------
  !                            *** PROGRAM cdffixtime  ***
  !
  !        ** Purpose: change time variable to jcness  deduce from time tag given in arguments
  !
  !   History:
  !          Jean-Marc Molines (March 2007) from old jcness
  !---------------------------------------------------------------------------------------
  USE cdfio
!  USE netcdf
  IMPLICIT NONE
  ! parameter to set the behaviour of the calendar : 365= no leap year
  !                                                  365.2425 = leap year
  REAL :: rpp_un_an = 365 !365.2425

  INTEGER :: narg, iargc, jarg
  INTEGER :: is, ie            !: starting and ending position of the tag in file name
  INTEGER :: iyear, imon, iday
  INTEGER :: iyr_init, imm_init, idd_init
  INTEGER :: ihr_init, imn_init, isec_init
  REAL(KIND=4) :: rdt_obs=5.   !: time interval between the observations (jcness will be offset by -rdt_obs/2
  REAL(KIND=4) :: rday0, rday_origin
  !  with respect to time tag
  REAL(KIND=4),DIMENSION(1) :: rdaycnes, rseconds
  CHARACTER(LEN=80) :: cfile, cdum, ctag='none', cdate, ctim
  CHARACTER(LEN=80) :: ctag0, ctim_unit, ctim_origin
  CHARACTER(LEN=3) :: cmm
  LOGICAL :: lnoleap=.true., lagrif=.false.

  ! Netcdf Stuff
  INTEGER :: istatus, ncid, id_time
  !---------------------------------------------------------------------------------------
  ! * 
  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdffixtime  -f file -i initial date [-t tag] [-leap] [ -noleap]'
     PRINT *,'        Change time_counter in file to set it according to drakkar rule'
     PRINT *,'     -i initial_date : to indicate time origin (yyyy-mm-dd hh:mm:ss) (2 words)'
     PRINT *,'     [-t tag ] : if not supplied, tag is taken from the name''s file'
     PRINT *,'            (assuming Drakkar convention ( CONFIG-CASE_tag_xxxxx.nc )'
     PRINT *,'     [-dt freq ] : number of days between model output [ 5 ]'
     PRINT *,'     [-leap ]   : assume a calendar with leap years'
     PRINT *,'     [-noleap ] : assume a calendar without leap years (default)'
     STOP
  ENDIF

  jarg=1
  DO WHILE ( jarg <= narg )
     CALL getarg(jarg, cdum) ; jarg=jarg + 1 
     SELECT CASE (cdum)
     CASE ( '-f' )
        CALL getarg(jarg,cfile) ;jarg=jarg +1
     CASE ( '-t' )
        CALL getarg(jarg,ctag) ; jarg=jarg +1
     CASE ( '-dt' )
        CALL getarg(jarg,cdum) ; jarg=jarg +1
        READ(cdum,*) rdt_obs
     CASE ( '-i' )
        CALL getarg(jarg,cdate) ; jarg=jarg +1
        CALL getarg(jarg,ctim)  ; jarg=jarg +1
     CASE ( '-leap' )
        rpp_un_an=365.2425
        lnoleap=.false.
     CASE ( '-noleap' )
        rpp_un_an=365
        lnoleap=.true.
     CASE DEFAULT 
         PRINT *,' Option ',TRIM(cdum),' unknown'
         STOP
     END SELECT
  END DO

  ! if ctag = none, try to find it from the file name.

  IF ( TRIM(ctag) == 'none' ) THEN
    is = INDEX(cfile,'_')
    IF ( is == 2 ) THEN 
      PRINT *,' ASSUME AGRIF file for ', TRIM(cfile)
      lagrif = .TRUE.
    ENDIF
    IF (lagrif) THEN
     is=INDEX(cfile(3:),'_' )+2
     ie=INDEX(cfile(is+1:),'_' )
     ctag=cfile(is+1:is+ie-1) 
    ELSE
     is=INDEX(cfile,'_')
     ie=INDEX(cfile(is+1:),'_' )
     ctag=cfile(is+1:is+ie-1) 
    ENDIF
  ENDIF

  PRINT *,' Changing time on file :', TRIM(cfile)
  is=INDEX(ctag,'d')
  IF ( is == 0 ) THEN ! not a model output but a mean value
    is=INDEX(ctag,'m') 
    IF ( is == 0 ) THEN ! annual mean set pseudo date to 01/07
      ctag=ctag(1:5)//"m07d01"
    ELSE  ! monthly mean
      ctag=ctag(1:8)//"d15"
    ENDIF
  ENDIF
  PRINT *,'            Using tag = ', TRIM(ctag)
 
  ! interpret ctim and cdate
  READ(cdate,'(i4,1x,i2,1x,i2)') iyr_init, imm_init, idd_init
  READ(ctim,'(i2,1x,i2,1x,i2)') ihr_init, imn_init, isec_init
  WRITE(ctag0,'("y",i4.4,"m",i2.2,"d",i2.2)') iyr_init, imm_init, idd_init
  rday0=jcnes(ctag0)+ihr_init/24.0 + imn_init/60./24. + isec_init/3600./24.
  rday_origin = rday0 - rdt_obs/2.
  CALL caldatjm( rday_origin, iyr_init, imm_init, idd_init, ihr_init, imn_init, isec_init)

  WRITE(cdate,'(i4.4,"-",i2.2,"-",i2.2)') iyr_init, imm_init, idd_init
  WRITE(ctim, '(i2.2,":",i2.2,":",i2.2)') ihr_init, imn_init, isec_init

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


  WRITE(ctim_unit,'("seconds since ",a,i3.2,":",i2.2,":",i2.2)') TRIM(cdate), ihr_init, imn_init, isec_init
  WRITE(ctim_origin,'(i5,"-",a,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)') iyr_init,cmm,idd_init, ihr_init, imn_init, isec_init
  PRINT *, iyr_init, imm_init, idd_init, ihr_init, imn_init, isec_init
  PRINT *, TRIM(ctim_unit)
  PRINT *, TRIM(ctim_origin)

  ! Compute corresponding jcnes
  rdaycnes=jcnes(ctag)
  rseconds=(rdaycnes - rday0 +1 ) * 86400. 

  ! Modify cdfile !! CAUTION : Original file will be modified  !!
  istatus = putvar1d( cfile, 'time_counter', rseconds, 1 )
  istatus = atted(cfile,'time_counter','units',ctim_unit)
  istatus = atted(cfile,'time_counter','time_origin',ctim_unit)

CONTAINS

  FUNCTION jcnes(cdtag)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(in) :: cdtag
    REAL(KIND=4) :: jcnes 
    ! local variables
    INTEGER :: iyear,imon,iday
    REAL(KIND=4)    :: sec=0.
    REAL(KIND=4)    :: rjuldeb, rjulfin, rjulday

    READ(cdtag,'(1x,i4.4,1x,i2.2,1x,i2.2)') iyear, imon, iday
    sec=0.
    !---------------------------------------------------------------------
    rjulfin = julday(iyear,imon,iday,sec)
    rjuldeb = julday(1950,01,01,0.)
    jcnes = rjulfin - rjuldeb
  END FUNCTION jcnes

  FUNCTION julday(kyear,kmonth,kday,rsec)
    !---------------------------------------------------------------------
    !- Converts year, month, day and seconds into a julian day
    !-
    !- In 1968 in a letter to the editor of Communications of the ACM
    !- (CACM, volume 11, number 10, October 1968, p.657) Henry F. Fliegel
    !- and Thomas C. Van Flandern presented such an algorithm.
    !-
    !- See also : http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm
    !-
    !- In the case of the Gregorian calendar we have chosen to use
    !- the Lilian day numbers. This is the day counter which starts
    !- on the 15th October 1582.
    !- This is the day at which Pope Gregory XIII introduced the
    !- Gregorian calendar.
    !- Compared to the true Julian calendar, which starts some
    !- 7980 years ago, the Lilian days are smaler and are dealt with
    !- easily on 32 bit machines. With the true Julian days you can only
    !- the fraction of the day in the real part to a precision of
    !- a 1/4 of a day with 32 bits.
    !---------------------------------------------------------------------
    IMPLICIT NONE
    !-
    INTEGER, INTENT(in)      :: kyear,kmonth,kday
    REAL(KIND=4),INTENT(in)  :: rsec
    REAL(KIND=4)             :: julday

    ! Local variables
    REAL,PARAMETER :: pp_un_jour = 86400.0
    INTEGER,PARAMETER :: jp_mon_len(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)

    INTEGER :: mm, iy, id, jd, ml
    INTEGER :: julian_day
    REAL    :: rjulian_sec

    CHARACTER(LEN=3),PARAMETER :: &
         &  cal(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
         &              'JUL','AUG','SEP','OCT','NOV','DEC'/)
    !---------------------------------------------------------------------

    mm = kmonth
    iy = kyear
    id = kday
    !-
    !- We deduce the calendar from the length of the year as it
    !- is faster than an INDEX on the calendar variable.
    !-
    !- Gregorian
    IF ( (rpp_un_an > 365.0).AND.(rpp_un_an < 366.0) ) THEN
       jd = (1461*(iy+4800+INT(( mm-14 )/12)))/4 &
            &      +(367*(mm-2-12*(INT(( mm-14 )/12))))/12 &
            &      -(3*((iy+4900+INT((mm-14)/12))/100))/4 &
            &      +id-32075
       jd = jd-2299160
       !- No leap or All leap
    ELSE IF (ABS(rpp_un_an-365.0) <= EPSILON(rpp_un_an) .OR. &
         &   ABS(rpp_un_an-366.0) <= EPSILON(rpp_un_an)) THEN
       ml = SUM(jp_mon_len(1:mm-1))
       jd = iy*INT(rpp_un_an)+ml+(id-1)
       !- Calendar with regular month
       !  ELSE
       !    ml = INT(un_an)/12
       !    jd = y*INT(un_an)+(m-1)*ml+(d-1)
    ENDIF
    !-
    julian_day = jd
    rjulian_sec = rsec
    julday = julian_day+rjulian_sec / pp_un_jour
  END FUNCTION julday

  SUBROUTINE caldatjm( pjcnes, ky, km, kd, kh, kmn, ksec )
  !!--------------------------------------------------------------------------------
  !!                  *** ROUTINE caldatjm ***
  !!
  !!   Purpose : return the calendar date from the jcnes given in argument
  !!
  !!   Method  : jcnes= 0 is 1950/01/01 00:00:00
  !!             Take care of leap/noleap year 
  !!--------------------------------------------------------------------------------
  REAL(KIND=4), INTENT(in) :: pjcnes
  INTEGER,      INTENT(out) :: ky, km, kd, kh, kmn, ksec
  
  INTEGER  :: isec, idays
  INTEGER  :: jd, jm
  INTEGER, DIMENSION(12) :: indays=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  INTEGER, DIMENSION(12) :: icumul
  !!--------------------------------------------------------------------------------
  icumul(1) = indays(1)
  DO jm=2,12
    icumul(jm)=icumul(jm-1)+indays(jm)
  ENDDO
  
  ! look for time 
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
    PRINT *, 'Not done yet for leap years'
  ENDIF
  END SUBROUTINE caldatjm


  
END PROGRAM cdffixtime
