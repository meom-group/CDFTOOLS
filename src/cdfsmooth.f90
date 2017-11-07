PROGRAM cdfsmooth
  !!======================================================================
  !!                     ***  PROGRAM  cdfsmooth  ***
  !!=====================================================================
  !!  ** Purpose :  perform a spatial filtering on input file.
  !!               - various filters are available :
  !!               1: Lanczos (default)
  !!               2: hanning
  !!               3: shapiro
  !!
  !!  ** Method  : read file level by level and perform a x direction 
  !!               filter, then y direction filter
  !!
  !! History : --   : 1995     : J.M. Molines : Original code for spem
  !!         : 2.1  : 07/2007  : J.M. Molines : port in cdftools
  !!         : 2.1  : 05/2010  : R. Dussin    : Add shapiro filter
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!           3.0  : 07/2011  : R. Dussin    : Add anisotropic box 
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!  filterinit   : initialise weight
  !!  filter       : main routine for filter computation
  !!  initlanc     : initialise lanczos weights
  !!  inithann     : initialise hanning weights
  !!  initshap     : initialise shapiro routine
  !!  initbox      : initialize weight for box car average
  !!  lislanczos2d : Lanczos filter
  !!  lishan2d     : hanning 2d filter
  !!  lisshapiro1d : shapiro filter
  !!  lisbox       : box car filter
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class data_transformation
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  !
  INTEGER(KIND=4), PARAMETER                    :: jp_lanc=1         ! lancszos id
  INTEGER(KIND=4), PARAMETER                    :: jp_hann=2         ! hanning id
  INTEGER(KIND=4), PARAMETER                    :: jp_shap=3         ! shapiro id
  INTEGER(KIND=4), PARAMETER                    :: jp_boxc=4         ! box car id
  INTEGER(KIND=4)                               :: jk, jt, jvar      ! dummy loop index
  INTEGER(KIND=4)                               :: npiglo, npjglo    ! size of the domain
  INTEGER(KIND=4)                               :: npk, npkf, npt    ! size of the domain
  INTEGER(KIND=4)                               :: narg, iargc       ! browse arguments
  INTEGER(KIND=4)                               :: ijarg             ! argument index for browsing line
  INTEGER(KIND=4)                               :: ncut, nband       ! cut period/ length, bandwidth
  INTEGER(KIND=4)                               :: nfilter = jp_lanc ! default value
  INTEGER(KIND=4)                               :: nvars, ierr       ! number of vars
  INTEGER(KIND=4)                               :: ncout             ! ncid of output file
  INTEGER(KIND=4)                               :: ilev              ! level to process if not 0
  INTEGER(KIND=4)                               :: ijk               ! indirect level addressing
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var            ! arrays of var id's
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk               ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout         ! id of output variables
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: iklist            ! list of k-level to process
  INTEGER(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: iw                ! flag for bad values (or land masked )

  REAL(KIND=4)                                  :: fn, rspval        ! cutoff freq/wavelength, spval
  REAL(KIND=4)                                  :: ranis             ! anistropy
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: gdep, gdeptmp     ! depth array 
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d, w2d          ! raw data,  filtered result

  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dtim              ! time array
  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dec, de           ! weight in r8, starting index 0:nband
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dec2d             ! working array

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar           ! struture for attribute

  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names          ! array of var name
  CHARACTER(LEN=256)                            :: cf_in, cf_out     ! file names
  CHARACTER(LEN=256)                            :: cv_dep, cv_tim    ! variable name for depth and time
  CHARACTER(LEN=256)                            :: ctyp              ! filter type
  CHARACTER(LEN=256)                            :: cldum             ! dummy character variable
  CHARACTER(LEN=256)                            :: clklist           ! ciphered k-list of level

  LOGICAL                                       :: lnc4 = .FALSE.    ! flag for netcdf4 output with chinking and deflation

  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfsmooth -f IN-file -c ncut [-t FLT-type] [-k LST-level] ...'
     PRINT *,'       [-anis ratio ] [-nc4 ] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Perform a spatial smoothing on the file using a particular filter as'
     PRINT *,'       specified in the ''-t'' option. Available filters are : Lanczos, Hanning,' 
     PRINT *,'       Shapiro and Box car average. Default is Lanczos filter.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f  IN-file  : input data file. All variables will be filtered'
     PRINT *,'       -c  ncut     : number of grid step to be filtered, or number'
     PRINT *,'                    of iteration of the Shapiro filter.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-t FLT-type] : Lanczos      , L, l  (default)'
     PRINT *,'                       Hanning      , H, h'
     PRINT *,'                       Shapiro      , S, s'
     PRINT *,'                       Box          , B, b'
     PRINT *,'       [-anis ratio ] : Specify an anisotropic ratio in case of Box-car filter.'
     PRINT *,'               With ratio=1, the box is a square 2.ncut x 2.ncut grid points.'
     PRINT *,'               In general, the box is then a rectangle 2.ncut*ratio x 2.ncut.'
     PRINT *,'       [-k LST-level ] : levels to be filtered (default = all levels)'
     PRINT *,'               LST-level is a comma-separated list of levels. For example,'
     PRINT *,'               the syntax 1-3,6,9-12 will select 1 2 3 6 9 10 11 12'
     PRINT *,'       [-nc4] : produce netcdf4 output file with chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       Output file name is build from input file name with indication'
     PRINT *,'       of the filter type (1 letter) and of ncut.'
     PRINT *,'       netcdf file :   IN-file[LHSB]ncut'
     PRINT *,'         variables : same as input variables.'
     PRINT *,'      '
     STOP 
  ENDIF
  !
  ijarg = 1
  ilev  = 0
  ranis = 1   ! anisotropic ratio for Box car filter
  ctyp  = 'L'
  ncut  = 0   ! hence program exit if none specified on command line
  DO WHILE (ijarg <= narg )
     CALL getarg ( ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE (cldum)
     CASE ( '-f'  ) ; CALL getarg ( ijarg, cf_in   ) ; ijarg=ijarg+1
     CASE ( '-c'  ) ; CALL getarg ( ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) ncut
     CASE ( '-t'  ) ; CALL getarg ( ijarg, ctyp    ) ; ijarg=ijarg+1 
     CASE ( '-k'  ) ; CALL getarg ( ijarg, clklist ) ; ijarg=ijarg+1 
                    ; CALL GetList (clklist, iklist, ilev )
     CASE ('-anis') ; CALL getarg ( ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) ranis
     CASE ( '-nc4') ; lnc4 = .TRUE.
     CASE DEFAULT   ; PRINT *,' ERROR :' ,TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( ncut == 0 ) THEN ; PRINT *, ' cdfsmooth : ncut = 0 --> nothing to do !' ; STOP 99
  ENDIF

  IF ( chkfile(cf_in) ) STOP 99 ! missing file

  !  remark: for a spatial filter, fn=dx/lambda where dx is spatial step, lamda is cutting wavelength
  fn    = 1./ncut
  nband = 2*ncut    ! Bandwidth of filter is twice the filter span

  ALLOCATE ( dec(0:nband) , de(0:nband) )

  WRITE(cf_out,'(a,a,i3.3)') TRIM(cf_in),'L',ncut   ! default name

  SELECT CASE ( ctyp)
  CASE ( 'Lanczos','L','l') 
     nfilter=jp_lanc
     WRITE(cf_out,'(a,a,i3.3)') TRIM(cf_in),'L',ncut
     PRINT *,' Working with Lanczos filter'
  CASE ( 'Hanning','H','h')
     nfilter=jp_hann
     ALLOCATE ( dec2d(0:2,0:2) )
     WRITE(cf_out,'(a,a,i3.3)') TRIM(cf_in),'H',ncut
     PRINT *,' Working with Hanning filter'
  CASE ( 'Shapiro','S','s')
     nfilter=jp_shap
     WRITE(cf_out,'(a,a,i3.3)') TRIM(cf_in),'S',ncut
     PRINT *,' Working with Shapiro filter'
  CASE ( 'Box','B','b')
     nfilter=jp_boxc
     WRITE(cf_out,'(a,a,i3.3)') TRIM(cf_in),'B',ncut
     IF ( ranis /=1. ) THEN
        PRINT *, 'Anisotropic box car with ratio Lx = ', ranis, 'x Ly'
     ELSE
        PRINT *,' Working with Box filter'
     ENDIF
  CASE DEFAULT
     PRINT *, TRIM(ctyp),' : undefined filter ' ; STOP 99
  END SELECT

  CALL filterinit (nfilter, fn, nband)
  ! Look for input file and create outputfile
  npiglo = getdim (cf_in,cn_x)
  npjglo = getdim (cf_in,cn_y)
  npk    = getdim (cf_in,cn_z, cdtrue=cv_dep, kstatus=ierr)
  npkf   = npk
  IF ( ierr /= 0 ) THEN
     npk   = getdim (cf_in,'z', cdtrue=cv_dep, kstatus=ierr)
     npkf  = npk
     IF ( ierr /= 0 ) THEN
        npk   = getdim (cf_in, 'sigma', cdtrue=cv_dep, kstatus=ierr)
        npkf  = npk
        IF ( ierr /= 0 ) THEN 
           PRINT *,' assume file with no depth'
           npk  = 1  ! Data have 1 level...
           npkf = 0  ! file have no deptht
        ENDIF
     ENDIF
  ENDIF
  npt    = getdim (cf_in,cn_t, cdtrue=cv_tim)

  PRINT *, 'npiglo = ',npiglo
  PRINT *, 'npjglo = ',npjglo
  PRINT *, 'npk    = ',npk
  PRINT *, 'npt    = ',npt

  ALLOCATE ( v2d(npiglo,npjglo),iw(npiglo,npjglo), w2d(npiglo,npjglo), dtim(npt) )
  nvars = getnvar(cf_in)
  PRINT *, 'nvars = ', nvars
  ALLOCATE (cv_names(nvars) )
  ALLOCATE (stypvar(nvars) )
  ALLOCATE (id_var(nvars),ipk(nvars),id_varout(nvars) )

  ALLOCATE ( gdeptmp(npk)  )
  IF (npkf /= 0 ) THEN ; gdeptmp(:) = getvar1d(cf_in, cv_dep, npk )  
  ELSE                 ; gdeptmp(:) = 0.  ! dummy value
  ENDIF

  ! get list of variable names and collect attributes in stypvar (optional)
  cv_names(:) = getvarname(cf_in, nvars, stypvar)

  DO jvar=1,nvars
     ! choose chunk size for output ... not easy not used if lnc4=.false. but anyway ..
     stypvar(jvar)%ichunk=(/npiglo,MAX(1,npjglo/30),1,1 /)
  ENDDO

  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cf_in, nvars, cdep=cv_dep)
  WHERE( ipk == 0 ) cv_names='none'
  stypvar(:)%cname=cv_names

  IF ( ilev /= 0 ) THEN   ! selected level on the command line
     WHERE (ipk(:) == npk ) ipk = ilev 
     npk = ilev
  ELSE                    ! all levels
     ilev = npk
     ALLOCATE(iklist(ilev) )
     iklist(:)=(/ (jk,jk=1,npk) /)
  ENDIF

  ALLOCATE ( gdep(ilev ) )
  gdep(:) = (/ (gdeptmp(iklist(jk)), jk=1,ilev) /)

  ! create output file taking the sizes in cf_in
  PRINT *, 'Output file name : ', TRIM(cf_out)
  ncout = create      (cf_out, cf_in,   npiglo, npjglo, npkf, cdep=cv_dep, ld_nc4=lnc4)
  ierr  = createvar   (ncout , stypvar, nvars,  ipk,    id_varout        , ld_nc4=lnc4)
  ierr  = putheadervar(ncout , cf_in,   npiglo, npjglo, npkf, pdep=gdep, cdep=cv_dep)
  dtim  = getvar1d(cf_in, cv_tim, npt)
  !
  DO jvar = 1,nvars
     IF ( cv_names(jvar) == cn_vlon2d .OR.                     &
          cv_names(jvar) == cn_vlat2d .OR. cv_names(jvar) == 'none' ) THEN
        ! skip these variables
     ELSE
        rspval=stypvar(jvar)%rmissing_value
        DO jt=1,npt
           DO jk=1,ipk(jvar)
              PRINT *, jt,'/',npt,' and ',jk,'/',ipk(jvar)
              ijk = iklist(jk) 
              v2d(:,:) = getvar(cf_in,cv_names(jvar),ijk,npiglo,npjglo,ktime=jt)
              iw(:,:) = 1
              WHERE ( v2d == rspval ) iw =0
              IF ( ncut /= 0 ) CALL filter( nfilter, v2d, iw, w2d)
              IF ( ncut == 0 ) w2d = v2d
              w2d  = w2d *iw  ! mask filtered data
              ierr = putvar(ncout, id_varout(jvar), w2d, jk, npiglo, npjglo, ktime=jt)
              !
           END DO
        END DO
     ENDIF
  END DO
  ierr = putvar1d(ncout, dtim, npt, 'T')
  ierr = closeout(ncout                )

CONTAINS

  SUBROUTINE filterinit(kfilter, pfn, kband)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE filterinit  ***
    !!
    !! ** Purpose :   initialise weight according to filter type
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) :: kfilter  ! filter number
    REAL(KIND=4),    INTENT(in) :: pfn      ! filter cutoff frequency/wavelength
    INTEGER(KIND=4), INTENT(in) :: kband    ! filter bandwidth
    !!----------------------------------------------------------------------
    SELECT CASE ( kfilter)
    CASE ( jp_lanc ) ; CALL initlanc (pfn, kband)
    CASE ( jp_hann ) ; CALL inithann (pfn, kband)
    CASE ( jp_shap ) ; CALL initshap (pfn, kband)
    CASE ( jp_boxc ) ; CALL initbox  (pfn, kband)
    END SELECT

  END SUBROUTINE filterinit

  SUBROUTINE filter (kfilter, px, kpx, py)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE filter  ***
    !!
    !! ** Purpose :  Call the proper filter routine according to filter type 
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),                 INTENT(in) :: kfilter ! filter number
    REAL(KIND=4), DIMENSION(:,:),    INTENT(in) :: px      ! input data
    INTEGER(KIND=4), DIMENSION(:,:), INTENT(in) :: kpx     ! validity flag
    REAL(KIND=4), DIMENSION(:,:),   INTENT(out) :: py      ! output data
    !!----------------------------------------------------------------------
    SELECT CASE ( kfilter)
    CASE ( jp_lanc ) ; CALL lislanczos2d (px, kpx, py, npiglo, npjglo, fn, nband)
    CASE ( jp_hann ) ; CALL lishan2d     (px, kpx, py, ncut, npiglo, npjglo)
    CASE ( jp_shap ) ; CALL lisshapiro1d (px, kpx, py, ncut, npiglo, npjglo)
    CASE ( jp_boxc ) ; CALL lisbox       (px, kpx, py, npiglo, npjglo, fn, nband, ranis)
    END SELECT

  END SUBROUTINE filter

  SUBROUTINE initlanc(pfn, knj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE initlanc  ***
    !!
    !! ** Purpose : initialize lanczos weights
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4),    INTENT(in) :: pfn  ! cutoff freq/wavelength
    INTEGER(KIND=4), INTENT(in) :: knj  ! bandwidth

    INTEGER(KIND=4)             :: ji   ! dummy loop index
    REAL(KIND=8)                :: dl_pi, dl_ey, dl_coef
    !!----------------------------------------------------------------------
    dl_pi   = ACOS(-1.d0)
    dl_coef = 2*dl_pi*pfn

    de(0) = 2.d0*pfn
    DO  ji=1,knj
       de(ji) = SIN(dl_coef*ji)/(dl_pi*ji)
    END DO
    !
    dec(0) = 2.d0*pfn
    DO ji=1,knj
       dl_ey   = dl_pi*ji/knj
       dec(ji) = de(ji)*SIN(dl_ey)/dl_ey
    END DO

  END SUBROUTINE initlanc

  SUBROUTINE inithann(pfn, knj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE inithann  ***
    !!
    !! ** Purpose : Initialize hanning weight 
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4),    INTENT(in) :: pfn  ! cutoff freq/wavelength
    INTEGER(KIND=4), INTENT(in) :: knj  ! bandwidth

    REAL(KIND=8)                :: dl_sum
    !!----------------------------------------------------------------------
    dec2d(:,:) = 0.d0 
    ! central point
    dec2d(1,1) = 4.d0
    ! along one direction
    dec2d(1,0) = 1.d0 ;  dec2d(1,2) = 1.d0
    ! and the other 
    dec2d(0,1) = 1.d0 ;  dec2d(2,1) = 1.d0

  END SUBROUTINE inithann

  SUBROUTINE initshap(pfn, knj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE initshap  ***
    !!
    !! ** Purpose :  Dummy routine to respect program structure 
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4),    INTENT(in) :: pfn  ! cutoff freq/wavelength
    INTEGER(KIND=4), INTENT(in) :: knj  ! bandwidth
    !!----------------------------------------------------------------------
    !   nothing to do 

  END SUBROUTINE initshap

  SUBROUTINE initbox(pfn, knj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE initbox  ***
    !!
    !! ** Purpose :  Init weights for box car 
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4),    INTENT(in) :: pfn  ! cutoff freq/wavelength
    INTEGER(KIND=4), INTENT(in) :: knj  ! bandwidth
    !!----------------------------------------------------------------------
    dec(:) = 1.d0

  END SUBROUTINE initbox


  SUBROUTINE lislanczos2d(px, kiw, py, kpi, kpj, pfn, knj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE lislanczos2d  ***
    !!
    !! ** Purpose : Perform lanczos filter
    !!
    !! ** Method  :   px      = input data
    !!                kiw     = validity of input data
    !!                py      = output filter
    !!                kpi,kpj = number of input/output data
    !!                pfn     = cutoff frequency
    !!                knj     = bandwith of the filter
    !!
    !! References : E. Blayo (1992) from CLS source and huge optimization 
    !!----------------------------------------------------------------------
    REAL(KIND=4),    DIMENSION(:,:), INTENT(in ) :: px               ! input array
    INTEGER(KIND=4), DIMENSION(:,:), INTENT(in ) :: kiw              ! flag input array
    REAL(KIND=4),    DIMENSION(:,:), INTENT(out) :: py               ! output array
    INTEGER(KIND=4),                 INTENT(in ) :: kpi, kpj         ! size of input/output
    REAL(KIND=4),                    INTENT(in ) :: pfn              ! cutoff frequency/wavelength
    INTEGER(KIND=4),                 INTENT(in ) :: knj              ! filter bandwidth

    INTEGER(KIND=4)                              :: ji, jj, jmx, jkx !  dummy loop index
    INTEGER(KIND=4)                              :: ik1x, ik2x, ikkx
    INTEGER(KIND=4)                              :: ifrst=0
    INTEGER(KIND=4)                              :: inxmin, inxmaxi
    INTEGER(KIND=4)                              :: inymin, inymaxi
    REAL(KIND=8), DIMENSION(kpi,kpj)             :: dl_tmpx, dl_tmpy
    REAL(KIND=8)                                 :: dl_yy, dl_den
    !!----------------------------------------------------------------------
    inxmin   =  knj
    inxmaxi  =  kpi-knj+1
    inymin   =  knj
    inymaxi  =  kpj-knj+1

    PRINT *,' filtering parameters'
    PRINT *,'    nx    = ', kpi
    PRINT *,'    nband = ', knj
    PRINT *,'    fn    = ', pfn

    DO jj=1,kpj
       DO  jmx=1,kpi
          ik1x = -knj
          ik2x =  knj
          !
          IF (jmx <= inxmin ) ik1x = 1-jmx
          IF (jmx >= inxmaxi) ik2x = kpi-jmx
          !
          dl_yy  = 0.d0
          dl_den = 0.d0
          !
          DO jkx=ik1x,ik2x
             ikkx=ABS(jkx)
             IF (kiw(jkx+jmx,jj)  ==  1) THEN
                dl_den = dl_den + dec(ikkx)
                dl_yy  = dl_yy  + dec(ikkx)*px(jkx+jmx,jj)
             END IF
          END DO
          !
          dl_tmpx(jmx,jj)=dl_yy/dl_den
       END DO
    END DO

    DO ji=1,kpi
       DO  jmx=1,kpj
          ik1x = -knj
          ik2x =  knj
          !
          IF (jmx <= inymin ) ik1x = 1-jmx
          IF (jmx >= inymaxi) ik2x = kpj-jmx
          !
          dl_yy  = 0.d0
          dl_den = 0.d0
          !
          DO jkx=ik1x,ik2x
             ikkx=ABS(jkx)
             IF (kiw(ji,jkx+jmx)  ==  1) THEN
                dl_den = dl_den + dec(ikkx)
                dl_yy  = dl_yy  + dec(ikkx)*dl_tmpx(ji,jkx+jmx)
             END IF
          END DO
          py(ji,jmx)=0.
          IF (dl_den /=  0.) py(ji,jmx) = dl_yy/dl_den
       END DO
    END DO
    !
  END SUBROUTINE lislanczos2d

  SUBROUTINE lishan2d(px, kiw, py, korder, kpi, kpj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE lishan2d  ***
    !!
    !! ** Purpose : compute hanning filter at order korder
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4),    DIMENSION(:,:), INTENT(in ) :: px      ! input data
    INTEGER(KIND=4), DIMENSION(:,:), INTENT(in ) :: kiw     ! validity flags
    REAL(KIND=4),    DIMENSION(:,:), INTENT(out) :: py      ! output data
    INTEGER(KIND=4),                 INTENT(in ) :: korder  ! order of the filter
    INTEGER(KIND=4),                 INTENT(in ) :: kpi, kpj ! size of the data

    INTEGER(KIND=4)                              :: jj, ji, jorder  ! loop indexes
    INTEGER(KIND=4)                              :: iiplus1, iiminus1
    INTEGER(KIND=4)                              :: ijplus1, ijminus1
    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE    :: ztmp
    !!----------------------------------------------------------------------
    ALLOCATE( ztmp(kpi,kpj) )

    py(:,:)   = 0.
    ztmp(:,:) = px(:,:)

    DO jorder = 1, korder
       DO jj   = 2, kpj-1
          DO ji = 2, kpi-1
             !treatment of the domain frontiers
             iiplus1 = MIN(ji+1,kpi) ; iiminus1 = MAX(ji-1,1) 
             ijplus1 = MIN(jj+1,kpj) ; ijminus1 = MAX(jj-1,1) 

             ! we don't compute in land
             IF ( kiw(ji,jj) == 1 ) THEN
                py(ji,jj) = SUM( dec2d(:,:) * ztmp(iiminus1:iiplus1,ijminus1:ijplus1) * kiw(iiminus1:iiplus1,ijminus1:ijplus1) )
                py(ji,jj) = py(ji,jj) / SUM( dec2d(:,:) * kiw(iiminus1:iiplus1,ijminus1:ijplus1) )   ! normalisation
             ENDIF
          ENDDO
       ENDDO
       ! update the ztmp array
       ztmp(:,:) = py(:,:)
    ENDDO

  END SUBROUTINE lishan2d

  SUBROUTINE lisshapiro1d(px, kiw, py, korder, kpi, kpj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE lisshapiro1d  ***
    !!
    !! ** Purpose :  compute shapiro filter 
    !!
    !! References :  adapted from Mercator code
    !!----------------------------------------------------------------------
    REAL(KIND=4),    DIMENSION(:,:), INTENT(in ) :: px      ! input data
    INTEGER(KIND=4), DIMENSION(:,:), INTENT(in ) :: kiw     ! validity flags
    REAL(KIND=4),    DIMENSION(:,:), INTENT(out) :: py      ! output data
    INTEGER(KIND=4),                 INTENT(in ) :: korder  ! order of the filter
    INTEGER(KIND=4),                 INTENT(in ) :: kpi, kpj ! size of the data

    INTEGER(KIND=4)                              :: jj, ji, jorder  ! loop indexes
    INTEGER(KIND=4)                              :: imin, imax, ihalo=0
    REAL(KIND=4), PARAMETER                      :: rp_aniso_diff_XY = 2.25 !  anisotrope case
    REAL(KIND=4)                                 :: zalphax, zalphay, znum
    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE    :: ztmp , zpx , zpy, zkiw
    LOGICAL                                      :: ll_cycl = .TRUE.
    !!----------------------------------------------------------------------

    IF(ll_cycl) ihalo=1
    ! we allocate with an ihalo
    ALLOCATE( ztmp(0:kpi+ihalo,kpj) , zkiw(0:kpi+ihalo,kpj) )
    ALLOCATE( zpx (0:kpi+ihalo,kpj) , zpy (0:kpi+ihalo,kpj) )

    IF(ll_cycl) THEN
       zpx(1:kpi,:) = px(:  ,:) ;  zkiw(1:kpi,:) = kiw(:  ,:)
       zpx(0    ,:) = px(kpi,:) ;  zkiw(0    ,:) = kiw(kpi,:)
       zpx(kpi+1,:) = px(1  ,:) ;  zkiw(kpi+1,:) = kiw(1  ,:)
    ELSE
       zpx(:    ,:) = px(:  ,:)
    ENDIF

    zpy (:,:) = zpx(:,:)  ! init?
    ztmp(:,:) = zpx(:,:)  ! init

    zalphax=1./2.
    zalphay=1./2.

    !  Dx/Dy=rp_aniso_diff_XY  , D_ = vitesse de diffusion
    !  140 passes du fitre, Lx/Ly=1.5, le rp_aniso_diff_XY correspondant est:

    IF ( rp_aniso_diff_XY >=  1. ) zalphay=zalphay/rp_aniso_diff_XY
    IF ( rp_aniso_diff_XY <   1. ) zalphax=zalphax*rp_aniso_diff_XY

    DO jorder=1,korder
       imin = 2     - ihalo
       imax = kpi-1 + ihalo
       DO ji = imin,imax
          DO jj = 2,kpj-1
             ! We crop on the coast
             znum =    ztmp(ji,jj)                                                  &
                  &    + 0.25*zalphax*(ztmp(ji-1,jj  )-ztmp(ji,jj))*zkiw(ji-1,jj  ) &
                  &    + 0.25*zalphax*(ztmp(ji+1,jj  )-ztmp(ji,jj))*zkiw(ji+1,jj  ) &
                  &    + 0.25*zalphay*(ztmp(ji  ,jj-1)-ztmp(ji,jj))*zkiw(ji  ,jj-1) &
                  &    + 0.25*zalphay*(ztmp(ji  ,jj+1)-ztmp(ji,jj))*zkiw(ji  ,jj+1)
             zpy(ji,jj) = znum*zkiw(ji,jj)+zpx(ji,jj)*(1.-zkiw(ji,jj))
          ENDDO  ! end loop ji
       ENDDO  ! end loop jj

       IF ( ll_cycl ) THEN
          zpy(0    ,:) = zpy(kpi,:) 
          zpy(kpi+1,:) = zpy(1  ,:) 
       ENDIF

       ! update the tmp array
       ztmp(:,:) = zpy(:,:)

    ENDDO

    ! return this array
    IF( ll_cycl ) THEN
       py(:,:) = zpy(1:kpi,:)
    ELSE
       py(:,:) = zpy(:    ,:)
    ENDIF

  END SUBROUTINE lisshapiro1d

  SUBROUTINE lisbox(px, kiw, py, kpi, kpj, pfn, knj,panis)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE lisbox  ***
    !!
    !! ** Purpose :  Perform box car filtering 
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4),    DIMENSION(:,:), INTENT(in ) :: px               ! input array
    INTEGER(KIND=4), DIMENSION(:,:), INTENT(in ) :: kiw              ! flag input array
    REAL(KIND=4),    DIMENSION(:,:), INTENT(out) :: py               ! output array
    INTEGER(KIND=4),                 INTENT(in ) :: kpi, kpj         ! size of input/output
    REAL(KIND=4),                    INTENT(in ) :: pfn              ! cutoff frequency/wavelength
    INTEGER(KIND=4),                 INTENT(in ) :: knj              ! filter bandwidth
    REAL(KIND=4),                    INTENT(in ) :: panis            ! anisotrop

    INTEGER(KIND=4)                              :: ji, jj
    INTEGER(KIND=4)                              :: ik1x, ik2x, ik1y, ik2y
    REAL(KIND=8)                                 :: dl_den
    LOGICAL, DIMENSION(kpi,kpj)                  :: ll_mask
    !!----------------------------------------------------------------------
    ll_mask=.TRUE.
    WHERE (kiw == 0 ) ll_mask=.FALSE.
    DO ji=1,kpi
       ik1x = ji-NINT( panis * knj)  ; ik2x = ji+NINT( panis * knj)
       ik1x = MAX(1,ik1x)            ; ik2x = MIN(kpi,ik2x)
       DO jj=1,kpj
          ik1y = jj-knj       ; ik2y = jj+knj
          ik1y = MAX(1,ik1y)  ; ik2y = MIN(kpj,ik2y)
          dl_den = SUM(kiw(ik1x:ik2x,ik1y:ik2y) )
          IF ( dl_den /= 0 ) THEN
             py(ji,jj) = SUM(px(ik1x:ik2x,ik1y:ik2y), mask=ll_mask(ik1x:ik2x,ik1y:ik2y) )/dl_den
          ELSE
             py(ji,jj) = rspval
          ENDIF
       END DO
    END DO

  END SUBROUTINE lisbox

END PROGRAM cdfsmooth
