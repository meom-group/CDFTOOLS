PROGRAM cdfvar
  !!----------------------------------------------------------------------------
  !!                ***   PROGRAM cdfvar ***
  !!
  !!  ** Purpose:    Locally transform a data file .... ???
  !!
  !!  ** Method:  Use OPA9 routine to look for zps. Locally force the depth to give
  !!              full depth. Save the modifs as source fortran code.
  !!
  !!  ** Usage :  cdfvar -f file -zoom imin imax jmin jmax 
  !!
  !!   History: 
  !!       2007 : J-M Molines : Original
  !!       2008 : P. Mathiot : Adaptation from cdfbathy for any variable of a file
  !!
  !!----------------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  ! * Module used
  USE cdfio

  ! * Local Variable
  IMPLICIT NONE
  !
  INTEGER :: numin,jk,ji,jj,jt,jl, jd, jarg
  INTEGER :: narg, iargc
  INTEGER :: imin, imax, jmin, jmax, klev, istatus, jtime
  INTEGER :: npiglo, npjglo, npk
  INTEGER, DIMENSION(:), ALLOCATABLE :: level
  INTEGER, DIMENSION (:,:), ALLOCATABLE :: mbathy, mask
  ! REAL(KIND=4) :: e3zps_min=25, e3zps_rat=0.2
  REAL(KIND=4) :: e3zps_min=1000, e3zps_rat=1, depmin=600.
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: gdept, gdepw, e3t, e3w
  !
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE     :: h, rtime
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: bathyin,bathy, e3_bot
  !
  CHARACTER(LEN=100) ::  cfilein, cline1, cline2, ctmp, cfileroot, creplace, cdump
  CHARACTER(LEN=80) :: cdim, cvar

  LOGICAL :: lexist=.TRUE., lfill=.FALSE., lfullstep=.FALSE., lappend=.FALSE., lreplace=.FALSE.
  LOGICAL :: ldump = .FALSE., lmodif=.FALSE., loverwrite=.false., lraz=.false., ldumpn=.false.
  INTEGER :: iversion=1, iostat, ipos
  !!
  !! 1. Initializations:
  !! -------------------
  !!
  narg = iargc()
  IF (narg == 0) THEN
     PRINT 9999,'USAGE :cdfvar -f file '// &
          '-zoom imin imax jmin jmax klev jtime -fillzone -fullstep depmin'
     PRINT 9999,'      -replace ''file'' -dumpzone ''file'' -a -o ' 
     PRINT 9999
     PRINT 9999, ' DESCRIPTION OF OPTIONS '
     PRINT 9999, ' ---------------------- '
     PRINT 9999, '   -file (or -f ) : name of var file '
     PRINT 9999, '   -var  (or -v ) : name of variable used '
     PRINT 9999, '   -zoom (or -z ) : sub area of the var file to work with (imin imax jmin jmax klev jtime)'
     PRINT 9999, '   -fillzone (or -fz ) : sub area will be filled with 0 up to the first coast line '
     PRINT 9999, '   -raz_zone (or -raz ) : sub area will be filled with 0 up '
     PRINT 9999, '   -fullstep (or -fs ) : sub area will be reshaped as full-step, below depmin'
     PRINT 9999, '               requires the presence of the file zgr_bat.txt (from ocean.output, eg )'
     PRINT 9999, '   -dumpzone (or -d ): sub area will be output to an ascii file, which can be used by -replace'
     PRINT 9999, '               after manual editing '
     PRINT 9999, '   -nicedumpzone (or -nd ): sub area will be output to an ascii file (nice output)'
     PRINT 9999, '   -replace (or -r ) : sub area defined by the file will replace the original var'
     PRINT 9999, '   -append (or -a )  : fortran log file (log.f90) will be append with actual modif'
     PRINT 9999, '              Standard behaviour is to overwrite/create log file'
     PRINT 9999, '   -overwrite (or -o ): input var file will be used as output.'
     PRINT 9999, '              Standard behaviour is to use a work copy of the original file'
     PRINT 9999, '               (indexed from 01 to 99 if necessary ) '
     STOP
  END IF
9999 FORMAT(a)
  ! Read command line
  jarg=1
  imin=-10 ; imax=-10 ; jmin=-10 ; jmax=-10
  DO  WHILE (jarg <=  narg)
     CALL getarg(jarg,cline1) ;  jarg = jarg + 1
     IF (cline1 == '-file ' .OR. cline1 == '-f') THEN
        CALL getarg(jarg,cline2) ; jarg = jarg + 1
        cfilein=cline2
     ELSE IF (cline1 == '-zoom' .OR. cline1 == '-z') THEN
        CALL getarg(jarg,cline2) ; jarg = jarg + 1
        READ(cline2,*) imin
        CALL getarg(jarg,cline2) ; jarg = jarg + 1
        READ(cline2,*) imax
        CALL getarg(jarg,cline2) ; jarg = jarg + 1
        READ(cline2,*) jmin
        CALL getarg(jarg,cline2) ; jarg = jarg + 1
        READ(cline2,*) jmax
       CALL getarg(jarg,cline2) ; jarg = jarg + 1
        READ(cline2,*) klev
       CALL getarg(jarg,cline2) ; jarg = jarg + 1
        READ(cline2,*) jtime
     ELSE IF (cline1 == '-var' .OR. cline1 == '-v' ) THEN
        CALL getarg(jarg,cline2) ; jarg = jarg + 1
        cvar=cline2
     ELSE IF (cline1 == '-fillzone' .OR. cline1 == '-fz' ) THEN
        lfill=.TRUE. ; lmodif=.TRUE.
     ELSE IF (cline1 == '-raz_zone' .OR. cline1 == '-raz' ) THEN
        lraz=.TRUE. ; lmodif=.TRUE.
     ELSE IF (cline1 == '-fullstep' .OR. cline1 == '-fs' ) THEN
        lfullstep=.TRUE. ; lmodif=.TRUE.
        CALL getarg(jarg,cline2) ; jarg = jarg + 1
        READ(cline2,*) depmin
     ELSE IF (cline1 == '-append' .OR. cline1 == '-a' ) THEN
        lappend=.TRUE.
     ELSE IF (cline1 == '-overwrite' .OR. cline1 == '-o' ) THEN
        loverwrite=.TRUE.
     ELSE IF (cline1 == '-replace' .OR. cline1 == '-r') THEN
        lreplace=.TRUE. ; lmodif=.TRUE.
        CALL getarg(jarg,creplace) ; jarg = jarg +1
     ELSE IF (cline1 == '-dumpzone' .OR. cline1 == '-d') THEN
        ldump=.TRUE.
        CALL getarg(jarg,cdump) ; jarg = jarg +1
     ELSE IF (cline1 == '-nicedumpzone' .OR. cline1 == '-nd') THEN
        ldumpn=.TRUE.
        CALL getarg(jarg,cdump) ; jarg = jarg +1
     ELSE
        PRINT *, cline1,' : unknown option '
        STOP
     END IF
  END DO

  IF ( lmodif .AND. .NOT. loverwrite) THEN
     ipos=INDEX(cfilein,'.',.TRUE.)
     READ(cfilein(ipos+1:),*,IOSTAT=iostat) iversion
     IF (iostat /=0 ) THEN 
        iversion=0
        cfileroot=cfilein
     ELSE
        cfileroot=cfilein(1:ipos-1)
     ENDIF
     iversion=iversion+1

     DO WHILE ( lexist )
        WRITE(ctmp,'(a,a,i2.2)') TRIM(cfileroot),'.',iversion
        INQUIRE(FILE=ctmp,EXIST=lexist)
        iversion=iversion+1
     END DO
     PRINT *, 'Working copy will be : ' ,TRIM(ctmp)
     CALL system(' cp -f '//cfilein//' '//ctmp )
  ELSE
     ctmp=cfilein
  ENDIF
  npiglo=getdim(ctmp,'x')
  npjglo=getdim(ctmp,'y')
  IF ( imin == -10 ) THEN  ! no zoom option passed
     imin=1 ; imax=npiglo
     jmin=1 ; jmax=npjglo
  END IF
  PRINT *, 'NPIGLO = ', npiglo
  PRINT *, 'NPJGLO = ', npjglo
  PRINT *, 'IMIN IMAX JMIN JMAX :', imin, imax,jmin,jmax

  ALLOCATE (mbathy(npiglo,npjglo), bathy(npiglo,npjglo),bathyin(npiglo,npjglo),e3_bot(npiglo,npjglo))
  ALLOCATE (mask(npiglo,npjglo))
  mask = 0
  bathy(:,:)=getvar(ctmp,cvar,klev, npiglo, npjglo, ktime=jtime)
  bathyin=bathy  ! save original 

  IF (lfullstep ) THEN 
     CALL zgr_read ; CALL zgr_zps(imin, imax, jmin, jmax)
  ENDIF
  IF (lfill ) CALL fillzone( imin, imax, jmin, jmax)
  IF (lraz )  CALL raz_zone( imin, imax, jmin, jmax)
  IF (ldump)    CALL dumpzone(cdump,imin, imax, jmin, jmax)
  IF (ldumpn)    CALL nicedumpzone(cdump,imin, imax, jmin, jmax)
  IF (lreplace) CALL replacezone(creplace)

  IF (lmodif ) THEN
     CALL prlog(bathyin,bathy,npiglo,npjglo,lappend)
     istatus=putvar(ctmp,cvar,klev,imax-imin+1,jmax-jmin+1,kimin=imin,kjmin=jmin,ktime=jtime,ptab=bathy(imin:imax,jmin:jmax))
  ENDIF

CONTAINS 
  SUBROUTINE zgr_zps ( kimin,kimax ,kjmin, kjmax )
    INTEGER ,INTENT(in) :: kimin,kimax, kjmin,kjmax
    !! * Local declarations
    INTEGER  ::   ji, jj, jk    ! dummy loop indices
    INTEGER  ::   ik, it            ! temporary integers
    INTEGER, PARAMETER :: wp=4

    REAL(wp) ::   &  
         ze3tp, ze3wp,    &  ! Last ocean level thickness at T- and W-points
         zdepwp,          &  ! Ajusted ocean depth to avoid too small e3t
         zdepth,          &  !    "         "
         zmax, zmin,      &  ! Maximum and minimum depth
         zdiff               ! temporary scalar

    ! Initialization of constant
    zmax = gdepw(npk) + e3t(npk)
    zmin = gdepw(4)

    ! initialize mbathy to the maximum ocean level available
    mbathy(kimin:kimax,kjmin:kjmax) = npk-1

    ! storage of land and island's number (zero and negative values) in mbathy
    WHERE (bathy(kimin:kimax,kjmin:kjmax) <= 0. ) mbathy(kimin:kimax,kjmin:kjmax)=INT( bathy(kimin:kimax,kjmin:kjmax) )

    ! bounded value of bathy
    ! minimum depth == 3 levels
    ! maximum depth == gdepw(jpk)+e3t(jpk) 
    ! i.e. the last ocean level thickness cannot exceed e3t(jpkm1)+e3t(jpk)
    WHERE (bathy(kimin:kimax,kjmin:kjmax) <= 0 ) 
       bathy(kimin:kimax,kjmin:kjmax)=0.
    ELSEWHERE (bathy(kimin:kimax,kjmin:kjmax) < zmin ) 
       bathy(kimin:kimax,kjmin:kjmax) = zmin
    ELSEWHERE (bathy(kimin:kimax,kjmin:kjmax) >= zmax )
       bathy(kimin:kimax,kjmin:kjmax) = zmax
    END WHERE

    ! Compute mbathy for ocean points (i.e. the number of ocean levels)
    ! find the number of ocean levels such that the last level thickness
    ! is larger than the minimum of e3zps_min and e3zps_rat * e3t (where
    ! e3t is the reference level thickness
    DO jk = npk-1, 1, -1
!      zdepth = gdepw(jk) + MIN( e3zps_min, e3t(jk)*e3zps_rat )
       zdepth = gdept(jk)
       WHERE ( bathy(kimin:kimax,kjmin:kjmax) > 0. .AND. bathy (kimin:kimax,kjmin:kjmax) <= zdepth )  
          mbathy(kimin:kimax,kjmin:kjmax)=jk-1
          e3_bot(kimin:kimax,kjmin:kjmax)= bathy(kimin:kimax,kjmin:kjmax) - gdepw(jk-1)
       END WHERE
    END DO

    DO ji=kimin,kimax
      DO jj=kjmin,kjmax
        jk=mbathy(ji,jj)
        IF (jk /= 0 ) THEN 
           IF (gdepw(jk+1) > depmin ) bathy(ji,jj)=gdepw(jk+1)-0.1
        ENDIF
      ENDDO
    END DO
  END SUBROUTINE zgr_zps

  SUBROUTINE zgr_read
    INTEGER :: numzgr = 10, il, iostat, idum
    CHARACTER(LEN=80) :: cline, cfile='zgrbat.txt'
    il=0
    OPEN(numzgr, FILE=cfile,IOSTAT=iostat)

    DO WHILE ( iostat == 0 )
       READ(numzgr,'(a)',IOSTAT=iostat) cline
       READ(cline,*,IOSTAT=idum )il
       IF ( idum == 0 )npk=il
    END DO

    ALLOCATE ( level(npk), gdept(npk), gdepw(npk), e3t(npk), e3w(npk) )
    REWIND(numzgr)

    il=0 ; iostat=0
    DO WHILE ( iostat == 0 )
       READ(numzgr,'(a)', IOSTAT=iostat) cline
       READ(cline,*,IOSTAT=idum) il
       IF ( idum == 0  ) READ(cline,*) level(il), gdept(il), gdepw(il), &
            &             e3t(il), e3w(il)
    END DO
  END SUBROUTINE zgr_read

  SUBROUTINE prlog (ptabold, ptab ,kpi,kpj,ldapp)
    ! * save differences in a log fill
    ! * if  ldapp results are append to the logfile
    INTEGER :: kpi,kpj
    REAL(KIND=4), DIMENSION(kpi,kpj) :: ptabold, ptab
    LOGICAL :: ldapp
    ! * Local variables
    INTEGER :: numlog=10

    IF (ldapp ) THEN 
       OPEN (numlog, FILE='log.f90', POSITION='append')
    ELSE
       OPEN (numlog, FILE='log.f90')
    ENDIF

    WRITE(numlog,'(a,a)') '! modification from original file : ', TRIM(cfilein)
    WRITE(numlog,'(a,a)') '! written to : ', TRIM(ctmp)
    DO ji=1,kpi
       DO jj=1,kpj
          IF ( ABS( ptabold(ji,jj) - ptab(ji,jj)) > 0.02  ) THEN    ! allow a 2 cm tolerance for rounding purposes
             WRITE(numlog,'(a,i4,a,i4,a,f8.2,a,f8.2)') ' bathy(',ji,',',jj,')=',ptab(ji,jj),' ! instead of ',ptabold(ji,jj)
          END IF
       END DO
    END DO
    CLOSE(numlog)
  END SUBROUTINE prlog

  SUBROUTINE fillzone(kimin,kimax,kjmin,kjmax)
    ! * Fill subzone of the bathy file
    INTEGER :: kimin, kimax, kjmin,kjmax
    INTEGER :: ji,jj
    DO jj=kjmin,kjmax
       ji=kimin
       IF ( bathy(ji,jj)  /= 0 ) THEN
          DO WHILE ( bathy(ji,jj)  /= 0 .AND. ji <= kimax )
             bathy(ji,jj) = 0.
             ji=ji+1
          END DO
       END IF
    END DO
  END SUBROUTINE fillzone

  SUBROUTINE raz_zone(kimin,kimax,kjmin,kjmax)
    ! * Fill subzone of the bathy file
    INTEGER :: kimin, kimax, kjmin,kjmax
    bathy(kimin:kimax, kjmin:kjmax) = 0.
  END SUBROUTINE raz_zone


  SUBROUTINE dumpzone(cdumpf,kimin,kimax,kjmin,kjmax)
    CHARACTER(LEN=*), INTENT(in) :: cdumpf
    INTEGER, INTENT(in) :: kimin,kimax,kjmin,kjmax
    INTEGER :: ji,jj
    INTEGER :: numdmp=20 , ni
    CHARACTER(LEN=80) :: cfmtr, cfmti
    !   PRINT *,' Dumpzone not yet operational' ; STOP
    ni=kimax-kimin+1
    WRITE(cfmtr,99) ni
    WRITE(cfmti,98) ni
    OPEN(numdmp,FILE=cdumpf)
    WRITE(numdmp,*) kimin,kimax,kjmin,kjmax, TRIM(cfmtr)
99  FORMAT('(I5,',i4.4,'f8.2)')
98  FORMAT('(5x,',i4.4,'I8)')
    WRITE(numdmp,cfmti)(ji,ji=kimin,kimax)
    DO jj= kjmax,kjmin,-1
       WRITE(numdmp,cfmtr) jj, bathy(kimin:kimax,jj)
    ENDDO
    CLOSE(numdmp)
  END SUBROUTINE dumpzone

  SUBROUTINE nicedumpzone(cdumpf,kimin,kimax,kjmin,kjmax)
    CHARACTER(LEN=*), INTENT(in) :: cdumpf
    INTEGER, INTENT(in) :: kimin,kimax,kjmin,kjmax
    INTEGER :: ji,jj
    INTEGER :: numdmp=20 , ni
    CHARACTER(LEN=80) :: cfmtr, cfmti
    ni=kimax-kimin+1
    WRITE(cfmtr,99) ni
    WRITE(cfmti,98) ni
    OPEN(numdmp,FILE=cdumpf)
    WRITE(numdmp,*) kimin,kimax,kjmin,kjmax, TRIM(cfmtr)
99  FORMAT('(I5,',i4.4,'I5)')
98  FORMAT('(5x,',i4.4,'I5)')
    WRITE(numdmp,cfmti)(ji,ji=kimin,kimax)
    DO jj= kjmax,kjmin,-1
       WRITE(numdmp,cfmtr) jj, INT(bathy(kimin:kimax,jj))
       WRITE(numdmp,*)
       WRITE(numdmp,*)
    ENDDO
    CLOSE(numdmp)
  END SUBROUTINE nicedumpzone


  SUBROUTINE replacezone(cdreplace)
    CHARACTER(LEN=*), INTENT(in) :: cdreplace
    INTEGER :: jj
    INTEGER :: iimin,iimax,ijmin,ijmax
    INTEGER :: numrep=20, idum
    !   PRINT *,' replacezone not yet operational' ; STOP
    OPEN(numrep,FILE=cdreplace)
    READ(numrep,*) iimin, iimax, ijmin, ijmax
    READ(numrep,*)  ! skip 1 line
    DO jj=ijmax,ijmin,-1
       READ(numrep,*) idum, bathy(iimin:iimax,jj)
    END DO
    CLOSE(numrep)
  END SUBROUTINE replacezone


END PROGRAM cdfvar
