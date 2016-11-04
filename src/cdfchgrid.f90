PROGRAM cdfchgrid
  !!======================================================================
  !!                     ***  PROGRAM  cdfchgrid  ***
  !!======================================================================
  !!  ** Purpose : Transform an 1442x1021 ORCA025 grid variable into an 
  !!               4322x3059 ORCA12 grid variable.
  !!               No interpolation, only copying one grid cell into 9 grid cells.
  !!
  !!  ** Method  : Store the result on a 'cdfchgrid.nc' file similar to the input file
  !!               (except x and y dimension)
  !!
  !!  ** Restriction  : Caution for mask coherence !
  !!                    This tool is only adapted for drowned field
  !!
  !! History : 3.0 !  08/2012    A. Lecointre   : Original code with Full Doctor form + Lic.
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   chgrid        : Convert coarser grid into refined grid
  !!----------------------------------------------------------------------
  USE cdfio
  USE modutils
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id: cdfchgrid.f90 XXX YYYY-MM-DD MM:MM:SSZ molines $
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jt                 ! dummy loop index
  INTEGER(KIND=4)                               :: ivar,iivar             ! dummy loop index
  INTEGER(KIND=4)                               :: ierr                   ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg     ! argument on line
  INTEGER(KIND=4)                               :: npiglo, npjglo         ! size of the input domain
  INTEGER(KIND=4)                               :: npigloout              ! I-size of the output domain
  INTEGER(KIND=4)                               :: npjgloout              ! J-size of the output domain
  INTEGER(KIND=4)                               :: npk, npkk, npt          ! size of the domain
  INTEGER(KIND=4)                               :: nvars                  ! number of variables in the input file
  INTEGER(KIND=4)                               :: ncout                  ! ncid of output ncdf file
  INTEGER(KIND=4), DIMENSION(1)                 :: ipk                    ! output variable : number of levels
  INTEGER(KIND=4), DIMENSION(1)                 :: id_varout              ! ncdf varid's

  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d                    ! array to read a layer of data 
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: u2d                    ! array onto ORCA12-grid
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: tim                    ! time counter of the file
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: rdep                   ! depth of the file

  CHARACTER(LEN=256)                            :: cf_out='cdfchgrid.nc'  ! output file name
  CHARACTER(LEN=256)                            :: cf_in                  ! input file name
  CHARACTER(LEN=256)                            :: cf_ref                 ! reference file for output file
  CHARACTER(LEN=256)                            :: cv_in                  ! variable name
  CHARACTER(LEN=256)                            :: cldum                  ! working string
  CHARACTER(LEN=256)                            :: cv_dep                 ! true name of dep dimension
  CHARACTER(LEN=256)                            :: cl_trf                 ! conversion key
  CHARACTER(LEN=256)                            :: cglobal                ! Global attribute with command line
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names               ! array of var name

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar                ! structure for variable attribute
  LOGICAL                                       :: lnc4=.FALSE.           ! flag for nc4 output with chinking and deflation
  LOGICAL                                       :: ldbg=.FALSE.           ! flag for nc4 output with chinking and deflation
  !!--------------------------------------------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfchgrid -f IN-file -r REF-file -var IN-var [-nc4] [-o OUT-file] [-d]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Build a new file on a refined grid, from a coarser grid, assuming that'
     PRINT *,'       the two grids are embedded, with common points (hence an odd scaling '
     PRINT *,'       factor). Grid characteristics are hard wired in the code. Support for'
     PRINT *,'       ORCA025 --> ORCA12, eORCA025 --> eORCA12 is actually provided. Hooks '
     PRINT *,'       are ready in the code for adding new conversion.'
     PRINT *,'       No interpolation, only copying value of a coarse grid cell, onto '
     PRINT *,'       scale x scale cells of the output grid (scale is the refinement factor)'
     PRINT *,'      '
     PRINT *,'     RESTRICTION :'
     PRINT *,'       Caution for mask coherence !'
     PRINT *,'       This tool is only adapted for drowned field'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file  : input Coarser-grid file'
     PRINT *,'       -r REF-file : Reference file used for identification of the output grid'
     PRINT *,'               should be of same geometry than the output file.'
     PRINT *,'       -var IN-var : input coarser-grid variable to be converted'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -nc4        : use netcdf4 chunking and deflation for the output file'
     PRINT *,'       -o OUT-file : specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       -d          : Display some debugging information '
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out)
     PRINT *,'         variable : same name as in input file'
     STOP
  ENDIF
  !!
  ijarg = 1
  ! Read command line
  DO  WHILE (ijarg <=  narg)
     CALL getarg(ijarg,cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-f' )
        CALL getarg(ijarg, cf_in) ; ijarg = ijarg + 1
     CASE ( '-r' )
        CALL getarg(ijarg, cf_ref) ; ijarg = ijarg + 1
     CASE ( '-var' )
        CALL getarg(ijarg,cv_in ) ; ijarg = ijarg + 1
     CASE ( '-o' )
        CALL getarg(ijarg,cf_out) ; ijarg = ijarg + 1
     CASE ( '-nc4' )
        lnc4 = .TRUE.
     CASE ( '-d' )
        ldbg = .TRUE.
     CASE DEFAULT
        PRINT *, TRIM(cldum),' : unknown option '
        STOP
     END SELECT
  END DO

  IF ( chkfile(cf_in) .OR. chkfile(cf_ref) ) STOP  ! missing files

  ! get domain dimension from input file
  npiglo = getdim (cf_in, cn_x)
  npjglo = getdim (cf_in, cn_y)
  npk    = getdim (cf_in, cn_z, cdtrue=cv_dep, kstatus=ierr)  ! defautl cn_z is depth
  npt    = getdim (cf_in, cn_t)

  IF ( npt == 0 ) THEN
     PRINT *, 'npt is forced to 1'
     npt = 1
  ENDIF

  IF ( npk == 0 ) THEN
    npkk = 1
  ELSE
    npkk = npk
  ENDIF


  ! get output domain dimension from reference file
  npigloout = getdim (cf_ref, cn_x)
  npjgloout = getdim (cf_ref, cn_y)
  ! infer convertion key from input/output sizes 
  cl_trf = 'none'

  ! check input size
  IF ( npiglo ==  722 .AND. npjglo ==  511 ) cl_trf='05to'   ! ORCA05
  IF ( npiglo == 1442 .AND. npjglo == 1021 ) cl_trf='025to'  ! ORCA025
  IF ( npiglo == 1442 .AND. npjglo == 1207 ) cl_trf='e025to' ! eORCA025
  IF ( npiglo == 4322 .AND. npjglo == 3059 ) cl_trf='12to'   ! ORCA12
  IF ( npiglo == 4322 .AND. npjglo == 3606 ) cl_trf='e12to'  ! eORCA12
  ! check output size
  IF ( npigloout ==  722 .AND. npjgloout ==  511 ) cl_trf=TRIM(cl_trf)//'05'
  IF ( npigloout == 1442 .AND. npjgloout == 1021 ) cl_trf=TRIM(cl_trf)//'025'
  IF ( npigloout == 1442 .AND. npjgloout == 1207 ) cl_trf=TRIM(cl_trf)//'e025'
  IF ( npigloout == 4322 .AND. npjgloout == 3059 ) cl_trf=TRIM(cl_trf)//'12'
  IF ( npigloout == 4322 .AND. npjgloout == 3606 ) cl_trf=TRIM(cl_trf)//'e12'


  PRINT *,' INPUT GRID '
  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt
  PRINT *
  PRINT *,' OUTPUT GRID '
  PRINT *, 'npigloout = ', npigloout
  PRINT *, 'npjgloout = ', npjgloout

  ALLOCATE ( v2d(npiglo+1,npjglo+1) )
  ALLOCATE ( u2d(npigloout,npjgloout) )
  ALLOCATE ( tim(npt) )
  ALLOCATE ( rdep(npkk) )

  ! look for the number of variables in the input file
  nvars = getnvar(cf_in)
  ALLOCATE (cv_names(nvars) ,stypvar(nvars))
  cv_names(:)=getvarname(cf_in,nvars,stypvar)

  ! find the number of variable we are interested in
  ivar=0
  DO  WHILE (ivar <   nvars)
     ivar=ivar+1
     IF ( cv_names(ivar) == cv_in ) iivar=ivar
  END DO
  cglobal="File produced with cdfchgrid "
  CALL SetGlobalAtt( cglobal, "A" )
  rdep=getvar1d(cf_in,cv_dep,npkk) 

  ipk(1)=npkk
  stypvar(iivar)%ichunk = (/npigloout,MAX(1,npjgloout/30),1,1 /)

  ncout = create      (cf_out,   cf_in  , npigloout, npjgloout, npk,   ld_nc4=lnc4  )
  ierr  = createvar   (ncout   , stypvar(iivar), 1 , ipk  , id_varout, ld_nc4=lnc4, cdglobal=cglobal  )
  ierr  = putheadervar(ncout,    cf_ref,  npigloout, npjgloout, npk , cdep='deptht', pdep=rdep     )

  ! get time and write time and get deptht and write deptht
  tim=getvar1d(cf_in,cn_t,npt)    ; ierr=putvar1d(ncout,tim,npt,'T')
                                    ierr=putvar1d(ncout,rdep,npk,'D')

  PRINT *,' Working with ', TRIM(cv_in), npk
  DO jt = 1, npt
     DO jk = 1, npk
        v2d(1:npiglo,1:npjglo) = getvar(cf_in, cv_in,  jk, npiglo, npjglo, ktime=jt)
        ! duplicate last  row and column 
        v2d(npiglo+1, :) = v2d(npiglo,:)
        v2d(:,npjglo+1 ) = v2d(:,npjglo)
        PRINT *,'level ',jk, 'time ',jt
        CALL chgrid(v2d, u2d, cl_trf)
        ierr = putvar ( ncout , id_varout(1), REAL(u2d), jk, npigloout, npjgloout, ktime=jt)
     ENDDO
  ENDDO

  ierr    = closeout(ncout)

CONTAINS

  SUBROUTINE chgrid (pinvar,poutvar,cd_trf)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE chgrid  ***
    !!
    !! ** Purpose :  This routine just copy input coarser grid
    !!              onto a finer inbedded grid.
    !!
    !! ** Method  :  Each case ( defined by the key cd_trf) is treated
    !!              separatly 
    !!----------------------------------------------------------------------

    REAL(KIND=4), DIMENSION(npiglo+1,npjglo+1),   INTENT(in ) :: pinvar
    REAL(KIND=4), DIMENSION(npigloout,npjgloout), INTENT(out) :: poutvar
    CHARACTER(LEN=*),                             INTENT(in ) :: cd_trf

    INTEGER(KIND=4)                                           :: jiin,jjin ! dummy loop index
    INTEGER(KIND=4)                                           :: iiin,ijin ! 
    INTEGER(KIND=4)                                           :: iiout,ijout !
    INTEGER(KIND=4)                                           :: ii_offset, ij_offset, iscal
    INTEGER(KIND=4)                                           :: iicomc,ijcomc,  iicomf,ijcomf

    SELECT CASE (cd_trf)
    CASE ('025to12')
       iscal  = 3
       iicomc = 2
       ijcomc = 499
       iicomf = 2
       ijcomf = 1495
       CALL filltab( pinvar, poutvar, iscal, iicomc, ijcomc, iicomf,ijcomf)
       ! force E-W periodicity
       poutvar(npigloout-1,:)= poutvar(1,:)
       poutvar(npigloout,  :)= poutvar(2,:)

    CASE ('05to025')
       ! to do ...
       PRINT * ,' Conversion ', TRIM(cd_trf), ' not supported yet!'
    CASE ('05to12')
       ! to do ...
       PRINT * ,' Conversion ', TRIM(cd_trf), ' not supported yet!'
    CASE ('e025toe12')
       iscal  = 3
       iicomc = 2
       ijcomc = 685
       iicomf = 2
       ijcomf = 2042
       CALL filltab( pinvar, poutvar, iscal, iicomc, ijcomc, iicomf,ijcomf)
       ! force E-W periodicity
       poutvar(npigloout-1,:)= poutvar(1,:)
       poutvar(npigloout,  :)= poutvar(2,:)
       ! 
    CASE DEFAULT
       PRINT *, TRIM(cd_trf),'  is not recognized !'
       PRINT *, 'No conversion will be performed'
    END SELECT

  END SUBROUTINE chgrid

  SUBROUTINE filltab(  pinvar, poutvar,  kscal, &
       &         kicomc, kjcomc, kicomf, kjcomf )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE filltab  ***
    !!
    !! ** Purpose :  fill higher resolution poutvar with values from pinvar
    !!               repeating the values over the iscal x iscale square  
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(npiglo+1,npjglo+1)  , INTENT(in )  :: pinvar
    REAL(KIND=4), DIMENSION(npigloout,npjgloout), INTENT(out)  :: poutvar
    INTEGER(KIND=4),                              INTENT(in )  :: kscal
    INTEGER(KIND=4),                           INTENT(inout )  :: kicomc, kjcomc
    INTEGER(KIND=4),                           INTENT(inout )  :: kicomf, kjcomf
    !
    INTEGER(KIND=4)                           :: jic, jjc, jif, jjf, ji,jj
    INTEGER(KIND=4)                           :: iiglof, ijglof, ipi,ipj
    INTEGER(KIND=4)                           :: iiofset, ijofset , iifc, ijfc, iiff, ijff
    INTEGER(KIND=4)                           :: ii,ij, ii0,ij0,  ii1,ij1

    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zwvar  ! working array to be cropped out after

    !------------------------------------------------------------------------
    iiglof=(npiglo+1) * kscal
    ijglof=(npjglo+1) * kscal
    ALLOCATE(zwvar(iiglof,ijglof))
    IF (ldbg) THEN
       PRINT *,' Size of working array: ', iiglof,ijglof
       PRINT *,'   last NE corner point ', kscal*(npiglo+1) -kscal/2 ,  kscal*(npjglo+1) -kscal/2
    ENDIF

    ! fill working fine array with all possible values from coarse array
    DO jjc = 1, npjglo+1
       ijfc=kscal * jjc - kscal/2  ! take care : integer division of an odd number
       DO jic = 1, npiglo+1
          iifc=kscal * jic - kscal/2  ! take care : integer division of an odd number
          DO jjf = ijfc - kscal/2, ijfc + kscal/2
             DO jif = iifc - kscal/2, iifc + kscal/2
                zwvar(jif,jjf) = pinvar(jic,jjc)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! crop working array to fit final fine array
    ! use common point information 
    ! Compute off set of grid from the inverce of formula : Pf = kscal * Pc + offset 
    iiofset = kicomf - kscal*kicomc
    ijofset = kjcomf - kscal*kjcomc
    ! find the most SW common point for the 2 grids ( brute force!)
    DO ji=1, npiglo+1
       ii= ji * kscal + iiofset
       IF ( ii >= 1) THEN
          kicomc=ji ; kicomf=ii
          EXIT
       ENDIF
    ENDDO

    DO jj=1, npjglo+1
       ij= jj * kscal + ijofset
       IF ( ij >= 1) THEN
          kjcomc=jj ; kjcomf=ij
          EXIT
       ENDIF
    ENDDO
    IF ( ldbg ) THEN
      PRINT *,'   Matching points : '
      PRINT *,'      coarse : ', kicomc, kjcomc
      PRINT *,'      fine   : ', kicomf, kjcomf
   ENDIF
    ! now do crop zwvar to fit poutvar  iiff,ijff are the index of SW most common point in zwvar
    iiff = kscal * kicomc - kscal/2 
    ijff = kscal * kjcomc - kscal/2 

!   ii0=iiff-kicomf +1 ; ii1=MIN(npigloout , ii0+npigloout-1) ; ipi=ii1-ii0+1
!   ij0=ijff-kjcomf +1 ; ij1=MIN(npjgloout , ij0+npjgloout-1) ; ipj=ij1-ij0+1
    ii0=iiff-kicomf +1 ; ii1=MIN(iiglof , ii0+npigloout-1) ; ipi=ii1-ii0+1
    ij0=ijff-kjcomf +1 ; ij1=MIN(ijglof , ij0+npjgloout-1) ; ipj=ij1-ij0+1
    IF ( ldbg) THEN
      PRINT *,'  Index limit to crop the array'
      PRINT *,'     ', ii0,ii1,ij0,ij1
      PRINT *,'   crop size : ', ipi,ipj
    ENDIF

    poutvar(:,:) = 999.
    poutvar(1:ipi,1:ipj) = zwvar(ii0:ii1,ij0:ij1)

    DEALLOCATE (zwvar )

  END SUBROUTINE filltab

END PROGRAM cdfchgrid
