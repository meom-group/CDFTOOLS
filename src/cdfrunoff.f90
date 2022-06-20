PROGRAM cdfrunoff
  !!======================================================================
  !!                     ***  PROGRAM  cdfrunoff  ***
  !!=====================================================================
  !!  ** Purpose : create a runoff file from surface tmask an gridded
  !!               input runoff file (such as ISBA)
  !!
  !! History :  4.0  : 06/2022  : J.M. Molines   from python version by J.Jouano
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE cdftools
  USE netcdf
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2022
  !! $Id$
  !! Copyright (c) 2022, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class mask
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  REAL(KIND=4), PARAMETER                      :: pp_rho=1000.       ! volumic mass of fresh water [kg/m3]

  INTEGER(KIND=4)                              :: ierr               ! error status
  INTEGER(KIND=4)                              :: ji, jj, jt,jc      ! dummy loop index
  INTEGER(KIND=4)                              :: ii, ij , inn, inm, ntot, ic
  INTEGER(KIND=4)                              :: ii1, ij1
  INTEGER(KIND=4)                              :: narg, iargc, ijarg ! command line 
  INTEGER(KIND=4)                              :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                              :: isum               ! sum of the 9 tmask point
  INTEGER(KIND=4)                              :: ncout              ! ncid of output file
  INTEGER(KIND=4)                              :: iwidth=1           ! width of coastal band
  INTEGER(KIND=4), DIMENSION(1)                :: icmin
  INTEGER(KIND=4), DIMENSION(2)                :: ipk, id_varout     ! outptut variables : number of levels,
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE   :: iimin, ijmin

  ! ISBA rvdis file related variables
  INTEGER(KIND=4)                              :: ncid_rdis, id_rdis
  INTEGER(KIND=4)                              :: ni_rdis, nj_rdis, nt_rdis

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE      :: dlon_rdis, dlat_rdis
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE      :: dtim_rdis
  REAL(KIND=8)                                 :: dspval_rdis, dmin, dradius=300  !km
  REAL(KIND=8)                                 :: disin, disout
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE    :: rivdis             ! input river discharge [m3/s]

  CHARACTER(LEN=20)                            :: cvlon_rdis='longitude'
  CHARACTER(LEN=20)                            :: cvlat_rdis='latitude'
  CHARACTER(LEN=20)                            :: cvtim_rdis='time_counter'
  CHARACTER(LEN=256)                           :: cv_rdis='rivdis'   ! input runoff variable (in m3/s)

  ! NEMO grid
  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: itmask             ! mask at jk level 
  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: icoast             ! masked cv_in at jk level
  INTEGER(KIND=4)                              :: ncoas              ! number of coastal point
  INTEGER(KIND=2), DIMENSION(:), ALLOCATABLE   :: iicoas, ijcoas     ! masked cv_in at jk level

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE    :: runoff             ! output runoff variable [kg/m2/s]
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE    :: rnfmsk             ! output variable socoefr [0-0.5]
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE    :: e1t, e2t           ! horizontal metrics [ m ] 

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE    :: dlon, dlat         ! horizontal metrics [ m ] 
  REAL(KIND=8), DIMENSION(:)  , ALLOCATABLE    :: distkm

  CHARACTER(LEN=256)                           :: cf_rdis            ! input file name for river discharge
  CHARACTER(LEN=256)                           :: cf_out='runoff.nc' ! output file name
  CHARACTER(LEN=256)                           :: cf_msk             ! input mask file name
  CHARACTER(LEN=256)                           :: cldum              ! dummy char var

  TYPE (variable), DIMENSION(2)                :: stypvar                  ! output attribute

  LOGICAL                                      :: lnc4 =.FALSE.      ! use Netcdf4 chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfrunoff -r RNF-file -f MASK-file [-v MASK-var] [-nc4] '
     PRINT *,'         [-vr RNF-var] [-radius RADIUS]  [-o OUT-file] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       This program is used to create a NEMO runoff file from a gridded'
     PRINT *,'       runoff data set (such as ISBA). We assume the runoff units in the'
     PRINT *,'       imput file is m3/s (like ISBA).'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -r RNF-file  : name of the gridded runoff  file ' 
     PRINT *,'       -f MASK-file  : name of the mask file ' 
     PRINT *,'      '
     PRINT *,'      OPTIONS : '
     PRINT *,'       -v MASK-var : input netcdf mask variable. ['//TRIM(cn_tmask)//'] '
     PRINT *,'       -vr RNF-var : input netcdf runoff variable.['//TRIM(cv_rdis)//'] '
     PRINT *,'       -nc4 : use netcdf4/Hdf5 chunking and deflation.'
     PRINT *,'       -o OUT-file     : name of coastal_mask file.['//TRIM(cf_out)//']'
     PRINT *,'       -radius RADIUS  : Distance threshold for runoff distribution on NEMO'
     PRINT *,'               (km)  [',dradius,'].'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        mesh_hgr.nc'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       runoff.nc file unless -o option is used.'
     PRINT *,'       Variables : sorunoff, socoefr'
     STOP 
  ENDIF

  !  Get parameters from command line
  ijarg = 1
  DO WHILE (ijarg <= narg)
     CALL getarg (ijarg, cldum ) ; ijarg = ijarg + 1
     SELECT CASE ( cldum)
     CASE ( '-r'   ) ; CALL getarg(ijarg, cf_rdis  )  ; ijarg = ijarg + 1
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_msk   )  ; ijarg = ijarg + 1
     CASE ( '-v'   ) ; CALL getarg(ijarg, cn_tmask )  ; ijarg = ijarg + 1
     CASE ( '-vr'  ) ; CALL getarg(ijarg, cv_rdis  )  ; ijarg = ijarg + 1
     CASE ( '-radius') ; CALL getarg(ijarg, cldum  )  ; ijarg = ijarg + 1 ; READ(cldum,*) dradius
     CASE ( '-nc4' ) ; lnc4 = .TRUE. 
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out   )  ; ijarg = ijarg + 1 
     CASE DEFAULT    ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF (  chkfile(cf_msk) .OR. chkfile(cn_fhgr) .OR. chkfile(cf_rdis)  ) STOP 99 ! missing files

  ! Open and allocate rdis files and data
  CALL OpenDisFile( cf_rdis )

  ! Look at NEMO grid
  npiglo = getdim (cf_msk,cn_x)
  npjglo = getdim (cf_msk,cn_y)

  PRINT *, 'NPIGLO = ', npiglo
  PRINT *, 'NPJGLO = ', npjglo

  ! Allocate arrays
  ALLOCATE( itmask(npiglo,npjglo) )
  ALLOCATE( icoast(npiglo,npjglo) )
  ALLOCATE( e1t(npiglo,npjglo), e2t(npiglo,npjglo) )
  ALLOCATE( dlon(npiglo,npjglo), dlat(npiglo,npjglo) )
  ALLOCATE( runoff(npiglo,npjglo), rnfmsk(npiglo,npjglo) )
  rnfmsk=0.

  ! create icoast mask
  itmask(:,:) = getvar(cf_msk, cn_tmask, 1, npiglo, npjglo)
  icoast(:,:) = 0
  DO jj = 2, npjglo-1
     DO ji = 2, npiglo-1
        IF ( itmask(ji,jj)  == 1 ) THEN
           isum=SUM(itmask(ji-1:ji+1, jj)) +  SUM(itmask(ji, jj-1:jj+1))
           ! check against 6 instead of 5 because central point (i,j) is
           ! counted twice
           IF ( isum /= 6 ) icoast(ji,jj) = icoast(ji,jj) + 1
        ENDIF
     ENDDO
  ENDDO

  ! read NEMO metric
  e1t(:,:)  = getvar(cn_fhgr,cn_ve1t, 1,npiglo,npjglo)
  e2t(:,:)  = getvar(cn_fhgr,cn_ve2t, 1,npiglo,npjglo)
  dlon(:,:) = getvar(cn_fhgr,cn_glamt, 1,npiglo,npjglo)
  dlat(:,:) = getvar(cn_fhgr,cn_gphit, 1,npiglo,npjglo)

  ! put coast points into a 1d vector
  ncoas=COUNT( (icoast /= 0 ) )
  ALLOCATE( iicoas(ncoas), ijcoas(ncoas), distkm(ncoas) )
  ALLOCATE( iimin(ncoas), ijmin(ncoas) )
  PRINT *, ' NCOAS = ', ncoas
  ic = 0 
  DO jj = 1, npjglo
     DO ji = 1 , npiglo
        IF ( icoast(ji,jj) /= 0  ) THEN
           ic = ic  + 1 
           iicoas(ic) = ji
           ijcoas(ic) = jj
        ENDIF
     ENDDO
  ENDDO

  CALL CreateOutput

  DO jt = 1,nt_rdis
     ierr = NF90_GET_VAR(ncid_rdis,id_rdis, rivdis, start=(/1,1,jt/), count=(/ni_rdis,nj_rdis,1/) )
     runoff(:,:) = 0.
     IF ( jt == 1 ) THEN
        ! at first time step determine on which NEMO grid point  the river discharge is to be
        ! attributed.
        inn = 0
        inm = 0
        disin=0.d0
        disout=0.d0
        ntot=COUNT( (rivdis /=  dspval_rdis ) )
        DO jj = 1, nj_rdis
           DO ji = 1, ni_rdis
              IF ( rivdis( ji, jj) /= dspval_rdis ) THEN
                 ! we are at (ji,jj) with a river discharge
                 DO jc = 1, ncoas
                    ! scan the coastal points and compute distance to local point
                    ii=iicoas(jc)
                    ij=ijcoas(jc)
                    distkm(jc)= dist(dlon_rdis(ji), dlon(ii,ij), dlat_rdis(jj),dlat(ii,ij) )
                 ENDDO
                 dmin=MINVAL( distkm )   ! look for minimum distance
                 icmin=MINLOC(distkm)    ! look for which coast point
                 IF ( dmin <= dradius )  THEN  ! river point should be near enough the NEMO point
                    inn = inn + 1 
                    ic = icmin(1)
                    ii=iicoas(ic)
                    ij=ijcoas(ic)
                    iimin(ic) = ji
                    ijmin(ic) = jj
                    !            print *, ji,jj,dmin, rivdis(ji,jj), ii,  ij
                    disin= disin + rivdis(ji,jj)
                    runoff(ii,ij) = rivdis(ji,jj)/e1t(ii,ij)/e2t(ii,ij)*pp_rho
                    rnfmsk(ii,ij) = 0.5
                 ELSE
                    disout= disout + rivdis(ji,jj)
                    inm = inm + 1 
                    PRINT *,' MISS ', dlon_rdis(ji), dlat_rdis(jj), rivdis(ji,jj) ,dmin
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
        PRINT *, 'INN : ', inn
        PRINT *, 'INM : ', inm
        PRINT *, 'NTOT : ', ntot
        PRINT *, 'DISIN : ', disin
        PRINT *, 'DISOUT : ', disout
     ELSE 
        DO jc = 1, ncoas
           ii=iicoas(jc)
           ij=ijcoas(jc)
           ii1=iimin(jc)
           ij1=ijmin(jc)
           runoff(ii,ij) = rivdis(ii1,ij1)/e1t(ii,ij)/e2t(ii,ij)*pp_rho
        ENDDO
     ENDIF
     PRINT *, 'JT done ', jt
     ierr = putvar(ncout, id_varout(1), runoff, 1, npiglo, npjglo, ktime=jt)
  ENDDO
  ierr = putvar(ncout, id_varout(2), rnfmsk, 1, npiglo, npjglo, ktime=1)

  ierr =closeout(ncout)

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
    REAL(KIND=8), DIMENSION(1) :: dltim
    !!----------------------------------------------------------------------
    ipk                          = 1

    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = 'sorunoff'

    stypvar(1)%cunits            = 'kg/m2/s'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = 0.
    stypvar(1)%valid_max         = 1.

    stypvar(1)%clong_name        = 'river runoff'

    stypvar(1)%cshort_name       = 'sorunoff'

    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TYX'
    stypvar(1)%cprecision        = 'r4'

    stypvar(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname             = 'socoefr'

    stypvar(2)%cunits            = '0-0.5'
    stypvar(2)%rmissing_value    = 0.
    stypvar(2)%valid_min         = 0.
    stypvar(2)%valid_max         = 0.5

    stypvar(2)%clong_name        = 'runoff mask'

    stypvar(2)%cshort_name       = 'socoefr'

    stypvar(2)%conline_operation = 'N/A'
    stypvar(2)%caxis             = 'TYX'
    stypvar(2)%cprecision        = 'r4'


    ncout = create      (cf_out, cf_msk,  npiglo, npjglo, 0, ld_nc4=lnc4 )
    ierr  = createvar   (ncout, stypvar, 2, ipk, id_varout,  ld_nc4=lnc4 )

    ierr  = putheadervar(ncout,  cf_msk,  npiglo, npjglo, 0)
    dltim = 0.d0
    ierr = putvar1d(ncout, dtim_rdis, nt_rdis,'T')


  END SUBROUTINE CreateOutput

  SUBROUTINE OpenDisFile (cd_dis) 
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE OpenDisFile  ***
    !!
    !! ** Purpose :    Open Discharge file (regular) and allocate disch. Arrays
    !!
    !! ** Method  :    Use netcdf primitives
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),    INTENT(in) :: cd_dis

    INTEGER(KIND=4) :: ierr, id

    ierr = NF90_OPEN(cd_dis,NF90_NOWRITE,ncid_rdis)
    ierr = NF90_INQ_DIMID(ncid_rdis,cvlon_rdis ,id) ; ierr = NF90_INQUIRE_DIMENSION(ncid_rdis, id, len=ni_rdis)
    ierr = NF90_INQ_DIMID(ncid_rdis,cvlat_rdis ,id) ; ierr = NF90_INQUIRE_DIMENSION(ncid_rdis, id, len=nj_rdis)
    ierr = NF90_INQ_DIMID(ncid_rdis,cvtim_rdis ,id) ; ierr = NF90_INQUIRE_DIMENSION(ncid_rdis, id, len=nt_rdis)

    PRINT *,' NI_RDIS = ', ni_rdis
    PRINT *,' NJ_RDIS = ', nj_rdis
    PRINT *,' NT_RDIS = ', nt_rdis

    ierr = NF90_INQ_VARID(ncid_rdis,cv_rdis, id_rdis)
    ierr = NF90_GET_ATT(ncid_rdis, id_rdis,'_FillValue', dspval_rdis )

    ALLOCATE(dlon_rdis(ni_rdis), dlat_rdis(nj_rdis), dtim_rdis(nt_rdis), rivdis(ni_rdis,nj_rdis) )
    ! READ lon lat time
    ierr = NF90_INQ_VARID(ncid_rdis, cvlon_rdis, id ) ; ierr = NF90_GET_VAR(ncid_rdis, id, dlon_rdis )
    ierr = NF90_INQ_VARID(ncid_rdis, cvlat_rdis, id ) ; ierr = NF90_GET_VAR(ncid_rdis, id, dlat_rdis )
    ierr = NF90_INQ_VARID(ncid_rdis, cvtim_rdis, id ) ; ierr = NF90_GET_VAR(ncid_rdis, id, dtim_rdis )


  END SUBROUTINE OpenDisFile

END PROGRAM cdfrunoff
