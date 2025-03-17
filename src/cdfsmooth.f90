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
  !!                     2024  : M.J. Bell    : Shapiro filter routine rewritten
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
  !!  lisshapiro2d : shapiro filter
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
  INTEGER(KIND=4)                               :: npass             ! number of passes of Shapiro filter
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
  LOGICAL                                       :: lap               ! choice 1D or 2D (Laplacian) for Shapiro filter 
  LOGICAL                                       :: lsm               ! T => land points not used (for Shapiro filter only)  
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: gdep, gdeptmp     ! depth array 
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d, w2d          ! raw data,  filtered result

  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dtim              ! time array
  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dec, de           ! weight in r8, starting index 0:nband
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dec2d             ! working array

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar           ! struture for attribute

  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names          ! array of var name
  CHARACTER(LEN=256)                            :: cf_in, cf_out, cf_nam   ! file names
  CHARACTER(LEN=256)                            :: cf_in_stem, cf_suffix   ! file name parts
  CHARACTER(LEN=256)                            :: cv_dep, cv_tim    ! variable name for depth and time
  CHARACTER(LEN=256)                            :: ctyp              ! filter type
  CHARACTER(LEN=256)                            :: cldum             ! dummy character variable
  CHARACTER(LEN=256)                            :: clklist           ! ciphered k-list of level

  LOGICAL                                       :: lnc4 = .FALSE.    ! flag for netcdf4 output with chinking and deflation

  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfsmooth -f IN-file -c ncut [-t FLT-type] [-npass npass] [-k LST-level] ...'
     PRINT *,'       [-anis ratio ] [-lap lap] [-lsm lsm] [-nc4 ] '
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
     PRINT *,'       [-n NAM-file] : namelist file. Only required for Shapiro filter.'
     PRINT *,'       [-t FLT-type] : Lanczos      , L, l  (default)'
     PRINT *,'                       Hanning      , H, h'
     PRINT *,'                       Shapiro      , S, s'
     PRINT *,'                       Box          , B, b'
     PRINT *,'       [-npass npass ] : Number of passes of filter; only used for Shapiro filter'
     PRINT *,'       [-k LST-level ] : levels to be filtered (default = all levels)'
     PRINT *,'               LST-level is a comma-separated list of levels. For example,'
     PRINT *,'               the syntax 1-3,6,9-12 will select 1 2 3 6 9 10 11 12'
     PRINT *,'       [-anis ratio ] : Specify an anisotropic ratio in case of Box-car filter.'
     PRINT *,'               With ratio=1, the box is a square 2.ncut x 2.ncut grid points.'
     PRINT *,'               In general, the box is then a rectangle 2.ncut*ratio x 2.ncut.'
     PRINT *,'       [-lap lap ] : .TRUE. implies Laplacian; .FALSE. implies 1D lines; only used for Shapiro filter'
     PRINT *,'       [-lsm lsm ] : .TRUE. implies land pts not used; .FALSE. implies full 2D field (no land/sea mask); only used for Shapiro filter'
     PRINT *,'       [-nc4] : produce netcdf4 output file with chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       Output file name is build from input file name with indication'
     PRINT *,'       of the filter type (1 letter) and of ncut .'
     PRINT *,'       netcdf file :   IN-file[LHSB]ncut'
     PRINT *,'         variables : same as input variables.'
     PRINT *,'      '
     STOP 
  ENDIF
  !
  cf_nam = ""
  ijarg = 1
  ilev  = 0
  ranis = 1   ! anisotropic ratio for Box car filter
  lap = .TRUE.
  lsm = .TRUE.
  ctyp  = 'L'
  ncut  = 0   ! hence program exit if none specified on command line
  npass = 1   ! just one pass
  DO WHILE (ijarg <= narg )
     CALL getarg ( ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE (cldum)
     CASE ( '-f'  ) ; CALL getarg ( ijarg, cf_in   ) ; ijarg=ijarg+1
     CASE ( '-n'  ) ; CALL getarg ( ijarg, cf_nam  ) ; ijarg=ijarg+1
     CASE ( '-c'  ) ; CALL getarg ( ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) ncut
     CASE ( '-t'  ) ; CALL getarg ( ijarg, ctyp    ) ; ijarg=ijarg+1 
     CASE ( '-npass') ; CALL getarg ( ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) npass
     CASE ( '-k'  ) ; CALL getarg ( ijarg, clklist ) ; ijarg=ijarg+1 
                    ; CALL GetList (clklist, iklist, ilev )
     CASE ('-anis') ; CALL getarg ( ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) ranis
     CASE ('-lap') ; CALL getarg ( ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) lap
     CASE ('-lsm') ; CALL getarg ( ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) lsm
     CASE ( '-nc4') ; lnc4 = .TRUE.
     CASE DEFAULT   ; PRINT *,' ERROR :' ,TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( ncut == 0 ) THEN ; PRINT *, ' cdfsmooth : ncut = 0 --> nothing to do !' ; STOP 99
  ENDIF

  IF ( chkfile(cf_in) ) STOP 99 ! missing file
  IF ( index( cf_in, ".nc" ) > 0 ) THEN
     cf_in_stem = cf_in(1:len(trim(cf_in))-3)
     cf_suffix = ".nc"
  ELSE
     cf_in_stem = cf_in
     cf_suffix = ""
  ENDIF

  !  remark: for a spatial filter, fn=dx/lambda where dx is spatial step, lamda is cutting wavelength
  fn    = 1./ncut
  nband = 2*ncut    ! Bandwidth of filter is twice the filter span

  ALLOCATE ( dec(0:nband) , de(0:nband) )

  WRITE(cf_out,'(a,a,i3.3,a)') TRIM(cf_in_stem),'_L',ncut,TRIM(cf_suffix)   ! default name

  SELECT CASE ( ctyp)
  CASE ( 'Lanczos','L','l') 
     nfilter=jp_lanc
     WRITE(cf_out,'(a,a,i3.3,a)') TRIM(cf_in_stem),'_L',ncut,TRIM(cf_suffix)
     PRINT *,' Working with Lanczos filter'
  CASE ( 'Hanning','H','h')
     nfilter=jp_hann
     ALLOCATE ( dec2d(0:2,0:2) )
     WRITE(cf_out,'(a,a,i3.3,a)') TRIM(cf_in_stem),'_H',ncut,TRIM(cf_suffix)
     PRINT *,' Working with Hanning filter'
  CASE ( 'Shapiro','S','s')
     nfilter=jp_shap
     WRITE(cf_out,'(a,a,2i1.1,2l1,a)') TRIM(cf_in_stem),'_S',ncut,npass,lap,lsm,TRIM(cf_suffix)
     PRINT *,' Working with Shapiro filter'
  CASE ( 'Box','B','b')
     nfilter=jp_boxc
     WRITE(cf_out,'(a,a,i3.3,a)') TRIM(cf_in_stem),'_B',ncut,TRIM(cf_suffix)
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
    CASE ( jp_shap ) ; CALL lisshapiro2d (px, kpx, py, ncut, npass, lap, lsm, npiglo, npjglo)
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

  SUBROUTINE lisshapiro2d(px, kiw, py, korder, kpass, lap, lsm, kpi, kpj)   
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE lisshapiro2d  ***
    !!
    !! ** Purpose :  apply a korder 2D shapiro filter kpass times. The land/sea mask and 
    !!               values at selected points may be forced to their initial values.  
    !!
    !! References :  Originally adapted from Mercator code.
    !!               This version rewritten by Mike Bell to implement the algorithm from
    !!               Francis, P.E., Quart. J. Roy. Met. Soc 101, pp 567-582 (1975)
    !!----------------------------------------------------------------------
    REAL(KIND=4),    DIMENSION(:,:), INTENT(in ) :: px      ! input data
    INTEGER(KIND=4), DIMENSION(:,:), INTENT(in ) :: kiw     ! validity flags
    REAL(KIND=4),    DIMENSION(:,:), INTENT(out) :: py      ! output data
    INTEGER(KIND=4),                 INTENT(in ) :: korder  ! order of the filter
    INTEGER(KIND=4),                 INTENT(in ) :: kpass   ! number of passes of the filter 
    LOGICAL,                         INTENT(in ) :: lap     ! True => Laplacian; False => lines
    LOGICAL,                         INTENT(in ) :: lsm     ! True => land points not used; False => no land/sea mask used (full 2D field) 
    INTEGER(KIND=4),                 INTENT(in ) :: kpi, kpj ! size of the data

    INTEGER(KIND=4)                              :: max_iterations = 10   ! max number of iterations allowed to enforce original land/sea mask and original values at selected points

    INTEGER(KIND=4)                              :: jn_fix_pts ! number of points where bathymetry is to remain fixed
    INTEGER(KIND=4)                              :: jj, ji, jorder, jpass, jiteration, jpt  ! loop indexes
    INTEGER(KIND=4)                              :: ijt     ! transposed ji index used for north pole fold
    INTEGER(KIND=4)                              :: ji_min, jj_min ! indices of min field values
    INTEGER(KIND=4)                              :: jcount_shallow, jcount_fixed  ! temporary indices for print out 
    REAL(KIND=4)                                 :: znum
    REAL(KIND=4)                                 :: rms_int, znpts
    REAL(KIND=4),  DIMENSION(:,:), ALLOCATABLE   :: ztmp , zpx , zpx_iteration, zpy, zkiw, zones
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE   :: ji_fix     ! i indices of bathymetry pts to hold fixed
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE   :: jj_fix     ! j indices of bathymetry pts to hold fixed
    LOGICAL                                      :: l_test_shallow, l_test_fixed ! temporary logicals 

!-----------------------------------------------------------------------------------------------
!! Variables that can be set through the namelist 

    LOGICAL                     :: ln_npol_fold      ! north pole fold in grid? 
    LOGICAL                     :: ln_ew_cycl        ! this filter has only been tested for EW cyclic grids 
!! MJB: for limited area models the bathymetry should be smoothed over a wider domain than that of the model; the margin should be at least korder*kpass
!! points and preferably 2*korder*kpass points (the bathymetry in nested models needs to match that of the outer domain in the nesting zone; that 
!! might be done by including the nesting zone in the list of fixed points - but simpler solutions might be adequate) 
  
    LOGICAL                     :: ln_single_point_response ! => study response to single non-zero value
    LOGICAL                     :: ln_pass_shallow_updates  ! => update values where bathy < zmin_val between passes
    LOGICAL                     :: ln_pass_fixed_pt_updates ! => update values where bathymetry should remain fixed

    CHARACTER(LEN=256)          :: cn_fixed_points          ! file containing list of fixed points
                                                            ! for ln_pass_fixed_pt_updates = T

    REAL(KIND=4)                :: rn_min_val, rn_factor_shallow, rn_tol_fixed, rn_tol_shallow ! see below

    INTEGER (KIND=4)            :: ji_single_pt, jj_single_pt ! indices of single point (see ln_single_point_response)
    INTEGER(KIND=4)             :: jst_prt, jend_prt ! 
 
!! The points that are printed out are set by ji_min_prt, jj_min_prt, ji_max_prt, jj_max_prt - these are local to prt_summary 

!-----------------------------------------------------------------------------------------------

    NAMELIST / nam_shapiro / ln_npol_fold, ln_ew_cycl, ln_single_point_response, ln_pass_shallow_updates, &
   &                         ln_pass_fixed_pt_updates, cn_fixed_points, ji_single_pt, jj_single_pt, jst_prt, jend_prt, &
   &                         rn_min_val, rn_tol_shallow, rn_tol_fixed, rn_factor_shallow 
!! Namelist default values 
    ln_npol_fold = .FALSE.
    ln_ew_cycl = .FALSE. 
    ln_single_point_response = .FALSE.
    ln_pass_shallow_updates = .FALSE.
    ln_pass_fixed_pt_updates = .FALSE.
    cn_fixed_points   = ""
    rn_min_val        = 5.0    ! minimum depth (e.g. 10.0 metres) 
    rn_tol_shallow    = 1.0    ! tolerance in metres of minimum shallow values (rn_min_val - rn_tol_shallow) 
    rn_tol_fixed      = 1.0    ! tolerance in metres for fit to bathymetry at point specified to be fixed  
    rn_factor_shallow = 1.5    ! rn_factor_shallow needs to be slightly greater than 1.0   ! 1.1 to 1.5 are reasonable values

    ji_single_pt = 622 ; jj_single_pt = 779
    jst_prt = 400 ;      jend_prt = 405

    IF( len(trim(cf_nam)) == 0 ) cf_nam='namelist_shapiro.txt'
    OPEN(UNIT=20, file = trim(cf_nam), form='formatted', status='old' )
    READ(NML=nam_shapiro, UNIT = 20) 
    WRITE(NML=nam_shapiro, UNIT=6)
    CLOSE(20)

    PRINT*, 'korder, kpass, lap, lsm = ', korder, kpass, lap, lsm
    
!-----------------------------------------------------------------------------------------------

! MJB 2021/02/10: This code has been re-written for files in which column 1 is the same as column kpi - 1; col 2 is same as col kpi.
!                 In other words the halo columns are included in the input file. 
!                 px  

    ! we do NOT allocate with an additional ihalo
    ALLOCATE( ztmp(kpi,kpj) , zkiw(kpi,kpj) )
    ALLOCATE( zpx (kpi,kpj) , zpy (kpi,kpj) )
    ALLOCATE( zones(kpi,kpj), zpx_iteration(kpi,kpj)  )
       
 ! print values from top and bottom rows to check they are OK 
    PRINT *, ' col 1 = ', (px(1,jj), jj = jst_prt, jend_prt)  
    PRINT *, ' col 2 = ', (px(2,jj), jj = jst_prt, jend_prt)  
    PRINT *, ' col kpi-1 = ', (px(kpi-1,jj), jj = jst_prt, jend_prt)  
    PRINT *, ' col kpi = ', (px(kpi,jj), jj = jst_prt, jend_prt)  

! option for testing out response to a single non-zero value 
    IF ( ln_single_point_response ) THEN 
       zpx(:,:) = 0.0 
       zpx(ji_single_pt,jj_single_pt) = 1.0     ! choose a point somewhere in the domain  
    END IF 

! read the indices of points to remain fixed 
    IF ( ln_pass_fixed_pt_updates ) THEN 
       OPEN (unit=20, file = trim(cn_fixed_points), form='formatted', status='old' ) 
       READ (20, *) jn_fix_pts

       ALLOCATE( ji_fix(jn_fix_pts), jj_fix(jn_fix_pts) )

       DO jpt = 1, jn_fix_pts
          READ (20, *) ji_fix(jpt), jj_fix(jpt)
       ENDDO 
       CLOSE (20)        
    END IF 

! write out land/sea mask and initial bathymetry values

    zpx(:,:)  = px(:  ,:)     ! px is only used at interior points and kiw is not used hereafter
    zkiw(:,:) = kiw(:,:)      ! halos could be inserted in zpx and zkiw here if not present in px and kiw 

    IF ( .NOT. lsm ) zkiw(:,:) = 1.0    !  used to test whether filter is stable when all points are included in the filter    
    
    zones(:,:) = 1.0                       
    PRINT *, ' point 1 zkiw' 
    CALL prt_summary( zkiw, zones, kpi, kpj) 

    IF ( rn_min_val > 0.0 ) THEN 
       DO jj = 2, kpj-1
         DO ji = 2,kpi-1
           IF ( zkiw(ji,jj) > 0.0 .AND. zpx(ji,jj) < rn_min_val ) zpx(ji,jj) = rn_min_val        
         ENDDO
       ENDDO 
    ENDIF 

    PRINT *, ' point 1 zpx' 
    CALL prt_summary( zpx, zkiw, kpi, kpj) 

! main calculations start 

! enforce cyclic conditions and north pole fold (southern boundary is assumed to be land) on zpx and zkiw
    CALL impose_bcs( zpx,  kpi, kpj, ln_ew_cycl, ln_npol_fold)
    CALL impose_bcs( zkiw, kpi, kpj, ln_ew_cycl, ln_npol_fold) 

    zpy (:,:) = zpx(:,:)  ! initialisation of zpy is necessary for row 1 (and other outer rows/columns if non-periodic)

    zpx_iteration(:,:) = zpx(:,:)    ! zpx_iteration  is only updated outside the jpass and jorder loops 
                                     ! zpx is a working array updated within the jpass loop    

    jiterationloop: DO jiteration=1,max_iterations

       zpx(:,:) = zpx_iteration(:,:)  ! initial values for zpx for this value of jiteration  

       jpassloop: DO jpass=1,kpass    
     
          ztmp(:,:) = zpx(:,:)  ! initialision of jorder loop 

          CALL impose_bcs( ztmp,  kpi, kpj, ln_ew_cycl, ln_npol_fold) 

          DO jorder=1,korder

             IF ( lap ) THEN 
                DO jj = 2,kpj-1
                   DO ji = 2,kpi-1
                      znum =      0.25*(ztmp(ji-1,jj  )-ztmp(ji,jj))*zkiw(ji-1,jj  ) &
                           &    + 0.25*(ztmp(ji+1,jj  )-ztmp(ji,jj))*zkiw(ji+1,jj  ) &
                           &    + 0.25*(ztmp(ji  ,jj-1)-ztmp(ji,jj))*zkiw(ji  ,jj-1) &
                           &    + 0.25*(ztmp(ji  ,jj+1)-ztmp(ji,jj))*zkiw(ji  ,jj+1)

                      zpy(ji,jj) = - 0.5*znum*zkiw(ji,jj)

                   ENDDO  ! end loop ji
                ENDDO  ! end loop jj
             ELSE 
                DO jj = 1,kpj
                   DO ji = 2,kpi-1
                      znum =      0.25*(ztmp(ji-1,jj  )-ztmp(ji,jj))*zkiw(ji-1,jj  ) &
                           &    + 0.25*(ztmp(ji+1,jj  )-ztmp(ji,jj))*zkiw(ji+1,jj  ) 

                      zpy(ji,jj) = - znum*zkiw(ji,jj)

                   ENDDO  ! end loop ji
                ENDDO  ! end loop jj

                ztmp(:,:) = zpy(:,:)

                DO jj = 2,kpj-1
                   DO ji = 2,kpi-1
                      znum =    + 0.25*(ztmp(ji  ,jj-1)-ztmp(ji,jj))*zkiw(ji  ,jj-1) &
                           &    + 0.25*(ztmp(ji  ,jj+1)-ztmp(ji,jj))*zkiw(ji  ,jj+1)

                      zpy(ji,jj) = - znum*zkiw(ji,jj)
 
                   ENDDO  ! end loop ji
                ENDDO  ! end loop jj
             ENDIF

             CALL impose_bcs( zpy, kpi, kpj, ln_ew_cycl, ln_npol_fold) 
       
             PRINT *, 'jorder, point 2 zpy = ', jorder 
             CALL prt_summary( zpy, zkiw, kpi, kpj) 

             ztmp(:,:) = zpy(:,:)  ! update ztmp for use with the next value of jorder

          ENDDO  ! jorder 

          zpy(:,:) = zpx(:,:) -  zpy(:,:)    !  zpy stores k-order filter after jpass iterations   
          zpx(:,:) = zpy(:,:)                !  update zpx for use with the next value of jpass

          PRINT *, 'jpass, point 3 zpx = ', jpass 
          CALL prt_summary( zpx, zkiw, kpi, kpj) 

       END DO jpassloop 

       IF ( .NOT. ( ln_pass_fixed_pt_updates .OR. ln_pass_shallow_updates ) ) EXIT jiterationloop  

! Find the number of points where the filtered bathymetry (zpy) is less than rn_min_val (to within tolerance rn_tol_shallow)  
     
       IF ( ln_pass_shallow_updates ) THEN 
          jcount_shallow = 0
          rms_int = 0.0 
          PRINT *, ' zpy(ji,jj), ji, jj, jcount where zpy < rn_min_val - rn_tol_shallow'  
          DO jj = 2, kpj-1
             DO ji = 2,kpi-1
                IF ( zkiw(ji,jj) > 0.0 .AND. zpy(ji,jj) < rn_min_val - rn_tol_shallow ) THEN  
                   jcount_shallow = jcount_shallow + 1
                   IF ( jcount_shallow < 50 ) PRINT *, zpy(ji,jj), ji, jj, jcount_shallow 
                ENDIF
                rms_int = rms_int + ( zpy(ji,jj) - px(ji,jj) )*( zpy(ji,jj) - px(ji,jj) ) 
             ENDDO
          ENDDO
          znpts = (kpj-2)*(kpi-2)
          rms_int = SQRT( rms_int / znpts ) 
          PRINT *, 'jcount_shallow = ', jcount_shallow
          PRINT *, 'rms_int = ', rms_int
       ENDIF ! ln_pass_shallow_updates

       IF ( ln_pass_fixed_pt_updates ) THEN 
          jcount_fixed = 0 
          PRINT *, ' px(ji,jj), zpy(ji,jj), ji, jj, jcount where ABS( zpy(ji,jj) - px(ji,jj) ) > rn_tol_fixed '  
          DO jpt = 1, jn_fix_pts 
             ji = ji_fix(jpt)
             jj = jj_fix(jpt)
             IF ( zkiw(ji,jj) > 0.0 .AND.  ABS( zpy(ji,jj) - px(ji,jj) ) > rn_tol_fixed ) THEN     ! zpy is updated value; px is original value 
                jcount_fixed = jcount_fixed + 1
                IF ( jcount_fixed < 50 ) PRINT *, px(ji,jj), zpy(ji,jj), ji, jj, jcount_fixed 
             ENDIF
          ENDDO
          PRINT *, 'jcount_fixed = ', jcount_fixed
       ENDIF ! ln_pass_fixed_pt_updates

! If output bathymetry is too small at some sea points or not close enough to the original at the selected points, increment zpx 

       l_test_shallow = jcount_shallow == 0 .OR. .NOT. ln_pass_shallow_updates 
       l_test_fixed   = jcount_fixed == 0   .OR. .NOT. ln_pass_fixed_pt_updates 
       IF ( l_test_shallow .AND. l_test_fixed  ) EXIT jiterationloop  

       IF ( ln_pass_shallow_updates ) THEN ! increment zpx before next pass at points where zpy < rn_min_val
          DO jj = 2, kpj-1
            DO ji = 2,kpi-1
              IF ( zkiw(ji,jj) > 0.0 .AND. zpy(ji,jj) < rn_min_val) THEN  
                zpx_iteration(ji,jj) = zpx_iteration(ji,jj) + MAX(px(ji,jj), rn_factor_shallow*rn_min_val) -  zpy(ji,jj)        
              ENDIF 
            ENDDO
          ENDDO 
       END IF 
       
       IF ( ln_pass_fixed_pt_updates ) THEN 
          DO jpt = 1, jn_fix_pts
             ji = ji_fix(jpt)
             jj = jj_fix(jpt)
             IF ( zkiw(ji,jj) > 0.0 ) zpx_iteration(ji,jj) = zpx_iteration(ji,jj) +  px(ji,jj) - zpy(ji,jj)  
          ENDDO 
       END IF 
       
    END DO jiterationloop 

    IF ( ln_pass_shallow_updates .AND. jcount_shallow == 0) THEN 
      DO jj = 1, kpj
        DO ji = 1,kpi
          IF ( zkiw(ji,jj) > 0.0 ) zpy(ji,jj) = MAX( zpy(ji,jj), rn_min_val )         
        ENDDO
      ENDDO
    ENDIF 
    
    py(:,:) = zpy(:    ,:)    ! first use of py (halos could easily be removed here) 
    
    DEALLOCATE( ztmp, zkiw )
    DEALLOCATE( zpx, zpy )
    DEALLOCATE( zones  )

  END SUBROUTINE lisshapiro2d

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

  SUBROUTINE prt_summary( pa, pkiw, kpi, kpj) 

  REAL(KIND=4),         DIMENSION(:,:), INTENT(in ) :: pa, pkiw   
  INTEGER(KIND=4),                      INTENT(in ) :: kpi,kpj

  INTEGER(KIND=4)    ::  ji, jj
  INTEGER(KIND=4)    ::  ji_min_prt, jj_min_prt
  INTEGER(KIND=4)    ::  ji_max_prt, jj_max_prt
  REAL(KIND=4)       ::  zmin, zmax

! User sets these values 
  ji_min_prt =  620 ;  jj_min_prt = 777 
  ji_max_prt =  624 ;  jj_max_prt = 781


  IF ( ji_max_prt > kpi )  ji_max_prt = kpi 
  IF ( jj_max_prt > kpj )  jj_max_prt = kpj 

  DO jj = jj_min_prt, jj_max_prt
    PRINT*, jj, (ji, pa(ji,jj), ji = ji_min_prt, ji_max_prt)     
  END DO 
  
  zmin = 1.E20  ; zmax = -1.E20 
  DO jj = 1, kpj
    DO ji = 1, kpi
       IF ( pkiw(ji,jj) .NE. 0 .AND. pa(ji,jj) < zmin ) THEN 
          ji_min_prt =  ji  ; jj_min_prt =  jj ; zmin = pa(ji,jj)  
       END IF 
       IF ( pkiw(ji,jj) .NE. 0 .AND. pa(ji,jj) > zmax ) THEN 
          ji_max_prt =  ji  ; jj_max_prt =  jj ; zmax = pa(ji,jj)  
       END IF 
    ENDDO 
  ENDDO

  PRINT*, 'ji, jj, min value = ', ji_min_prt, jj_min_prt, zmin 
  PRINT*, 'ji, jj, max value = ', ji_max_prt, jj_max_prt, zmax 

  RETURN
  END SUBROUTINE prt_summary
  
  SUBROUTINE impose_bcs( pfld, kpi, kpj, ln_ew_cycl, ln_npol_fold) 
  REAL(KIND=4),    DIMENSION(:,:), INTENT(inout) :: pfld      
  INTEGER,                         INTENT(in   ) :: kpi, kpj
  LOGICAL,                         INTENT(in   ) :: ln_ew_cycl, ln_npol_fold

  INTEGER ji, ijt 
    
! cyclic points in ji 
  IF ( ln_ew_cycl ) THEN
     pfld(1  ,:) = pfld(kpi-1,:) 
     pfld(kpi,:) = pfld(2    ,:) 
  ENDIF

! north pole fold;  (southern boundary is assumed to be land) 
  IF ( ln_npol_fold ) THEN 
     DO ji = 2, kpi 
        ijt = kpi-ji+2
        pfld(ji,kpj) = pfld(ijt,kpj-2)
     END DO 
     pfld(1,kpj) = pfld(3,kpj-2)
     DO ji = kpi/2+1, kpi
        ijt = kpi-ji+2
        pfld(ji,kpj-1) = pfld(ijt,kpj-1)
     END DO
  END IF 

  RETURN
  END SUBROUTINE impose_bcs

END PROGRAM cdfsmooth
