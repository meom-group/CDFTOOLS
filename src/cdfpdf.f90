PROGRAM cdfpdf
  !!======================================================================
  !!                     ***  PROGRAM  cdfpdf  ***
  !!=====================================================================
  !!  ** Purpose : Build the pdf for a given variable on a given area.
  !!
  !!  ** Method  : Establish a data binning and count the number of element
  !!               in each bin.  Binning is defined from given  minimum 
  !!               value, maximum value, and number of bins
  !!
  !! History : 3.0  : 10/2015  : J.M. Molines : Original code
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class statistics
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jk, jt
  INTEGER(KIND=4)                           :: npiglo, npjglo            ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt                  ! size of the domain
  INTEGER(KIND=4)                           :: npi, npj, ilev            !
  INTEGER(KIND=4)                           :: narg, iargc, ijarg
  INTEGER(KIND=4)                           :: numout=10
  INTEGER(KIND=4)                           :: imin, imax 
  INTEGER(KIND=4)                           :: jmin, jmax
  INTEGER(KIND=4)                           :: nbin, nlim, ibin
  INTEGER(KIND=4)                           :: ncout             ! ncid of output variable
  INTEGER(KIND=4)                           :: ierr              ! error status
  INTEGER(KIND=4), DIMENSION(1)             :: ipk , id_varout   ! output variable

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zvar, zcount, ztimed, zlon
  REAL(KIND=4), DIMENSION(:)  , ALLOCATABLE :: vlim
  REAL(KIND=4)                              :: vmin, vmax, bin_siz
  REAL(KIND=4)                              :: below, above, spval

  REAL(KIND=8), DIMENSION(:)  , ALLOCATABLE :: dtim

  CHARACTER(LEN=256)                        :: cf_ifil
  CHARACTER(LEN=256)                        :: cf_asc='pdf.txt'
  CHARACTER(LEN=256)                        :: cf_out='pdf.nc'
  CHARACTER(LEN=256)                        :: cv_nam
  CHARACTER(LEN=256)                        :: cldum
  CHARACTER(LEN=256)                        :: cglobal

  TYPE(variable), DIMENSION(1)              :: stypvar           ! output data structure

  LOGICAL                                   :: lchk
  LOGICAL                                   :: l_nozoom  =.TRUE.
  LOGICAL                                   :: l_norange =.TRUE.
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfpdf -f IN-file -v IN-var [-zoom imin imax jmin jmax] ..'
     PRINT *,'       [-lev level] [-range vmin vmax nbin ] [-o OUT-ncfile] [-a OUT-ascfile]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Build the pdf of a given variable, on a given area, according'
     PRINT *,'       to bin specifications passed as argument of the program. If no' 
     PRINT *,'       particular specification is passed to the program, build 100 '
     PRINT *,'       bins between minimum and maximum value of the variable.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file : input file '
     PRINT *,'       -v IN-var  : variable name '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-zoom imin imax jmin jmax] : define a sub-area, in model '
     PRINT *,'                                     coordinates' 
     PRINT *,'       [-lev level ] : choose a level for pdf computation '
     PRINT *,'              If not specified, takes level 1. '
     PRINT *,'       [-range vmin vmax nbin  ] : define the limit for binning '
     PRINT *,'                               and number of bins.'
     PRINT *,'       [-o OUT-file] : specify name for netcdf output file'
     PRINT *,'       [-a ASC-file] : specify name for ascii output file'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'          (1) ascii file with bin number, value and mean field in bin'
     PRINT *,'          (2) netcdf file for 2d array where x dimension corresponds to bins'
     PRINT *,'              y dimension corresponds to time, thus the field value being '
     PRINT *,'              an array count(bin,time). The output file follows the nemo '
     PRINT *,'              standards, even, if nav_lon, nav_lat are no more longitude or'
     PRINT *,'              latitude.'
     PRINT *,'              netdf variable is <IN-var>_pdf'
     PRINT *,'      '
     PRINT *,'     SEE ALSO : '
     STOP 
  ENDIF

  ijarg = 1
  ilev  = 1
  l_nozoom = .TRUE.
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_ifil) ; ijarg=ijarg+1
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out ) ; ijarg=ijarg+1
     CASE ( '-a'   ) ; CALL getarg(ijarg, cf_asc ) ; ijarg=ijarg+1
     CASE ( '-v'   ) ; CALL getarg(ijarg, cv_nam ) ; ijarg=ijarg+1
     CASE ( '-zoom') ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) imin
        ;              CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) imax
        ;              CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) jmin
        ;              CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) jmax
        ;              l_nozoom = .FALSE.
     CASE ( '-lev' ) ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) ilev
     CASE ('-range') ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) vmin
        ;              CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) vmax
        ;              CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) nbin
        ;              l_norange = .FALSE.
     CASE DEFAULT    ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF (  chkfile ( cf_ifil) ) STOP 99  ! some compulsory files are missing

  ! set domain size from input ile
  npiglo = getdim (cf_ifil,cn_x)
  npjglo = getdim (cf_ifil,cn_y)
  npk    = getdim (cf_ifil,cn_z)
  npt    = getdim (cf_ifil,cn_t)

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  ! if no zoom specified, take the full domain
  IF ( l_nozoom ) THEN
     imin=1 ; imax=npiglo
     jmin=1 ; jmax=npjglo
  ENDIF

  npi=imax - imin + 1
  npj=jmax - jmin + 1
  ALLOCATE ( zvar(npi,npj) )

  ! takle the case when no range specified for binning
  IF ( l_norange ) THEN
     nbin =100
     spval = getspval( cf_ifil, cv_nam) 
     ! scan the file for min and max
     vmin=1.e10
     vmax=-1.e10
     DO jt=1,npt
        zvar(:,:) = getvar( cf_ifil, cv_nam, ilev, npi, npj, ktime=jt, kimin=imin, kjmin=jmin )
        zvar(:,:) = getvar( cf_ifil, cv_nam, ilev, npi, npj,  ktime=jt, kimin=imin, kjmin=jmin )
        vmin=MIN(vmin, MINVAL(zvar, (zvar /= spval) ) )
        vmax=MAX(vmax, MAXVAL(zvar, (zvar /= spval) ) )
     ENDDO
  ENDIF

  bin_siz=(vmax - vmin ) / nbin
  nlim = nbin+1

  PRINT *, ' NPI  = ', npi
  PRINT *, ' NPJ  = ', npj
  PRINT *, ' IMIN = ', imin
  PRINT *, ' IMAX = ', imax
  PRINT *, ' JMIN = ', jmin
  PRINT *, ' JMAX = ', jmax
  PRINT *, ' ILEV = ', ilev
  PRINT *, ' VMIN = ', vmin
  PRINT *, ' VMAX = ', vmax
  PRINT *, ' NBIN = ', nbin
  PRINT *, ' NLIM = ', nlim
  PRINT *, ' BINS = ', bin_siz


  ! Allocate memory
  ALLOCATE ( zcount(nbin,npt), vlim(nlim) )
  ALLOCATE ( dtim(npt), ztimed(nbin,npt) , zlon(nbin,npt))

  DO ji=1, nlim 
     vlim(ji)= vmin + (ji-1)*bin_siz
  ENDDO

  OPEN (numout, FILE=cf_asc)   ! this file can be plotted easily with graph
  ! time in seconds read from file 
  dtim  = getvar1d(cf_ifil, cn_vtimec, npt )

  ! convert in time in days since the begining of the file
  ! this will be the dummy 'latitude' for the output file

  DO jt=1, npt
     ztimed(:,jt) = (dtim(jt) - dtim(1) ) / 86400.
  ENDDO

  ! dummy longitude for the output file in the mean value of the bin
  DO ji = 1, nbin
     zlon(ji,:) =  (vlim(ji) + vlim(ji+1)) /2.
  ENDDO

  CALL CreateOutput 

  DO jt = 1, npt
     zcount(:,jt) = 0. ; below=0. ; above=0.
     zvar(:,:) = getvar( cf_ifil, cv_nam, ilev, npi, npj,  ktime=jt, kimin=imin, kjmin=jmin )
     DO jj=1,npj
        DO ji=1,npi
           ibin=INT((zvar(ji,jj)-vmin)/bin_siz) +1
           IF ( ibin < 1 ) THEN 
              below=below+1
           ELSE IF ( ibin > nbin ) THEN
              above=above+1
           ELSE
              zcount(ibin,jt) = zcount(ibin,jt) + 1
           ENDIF
        ENDDO
     ENDDO

     WRITE(numout,*) 
     WRITE(numout,*) vlim(1),  below
     DO ji=1, nbin
        WRITE(numout,*)  (vlim(ji) + vlim(ji+1)) /2.,  zcount (ji,jt )
     ENDDO
     WRITE(numout,*)  vlim(nlim),  above
  ENDDO
  ierr = putvar( ncout, id_varout(1), zcount, 1, nbin, npt)
  ierr = closeout(ncout)
  CLOSE(numout)

CONTAINS
  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Set up all things required for the output file, create
    !!               the file and write the header part. 
    !!
    !! ** Method  :  Use global module variables 
    !!
    !!----------------------------------------------------------------------
    ipk(:) = 1

    stypvar(1)%cname             = 'pdf_'//TRIM(cv_nam)
    stypvar(1)%cunits            = 'N/A'
    stypvar(1)%rmissing_value    = -1000.
    stypvar(1)%valid_min         = 0.
    stypvar(1)%valid_max         = npi*npj
    stypvar(1)%clong_name        = 'PDF of '//TRIM(cv_nam)
    stypvar(1)%cshort_name       = 'pdf_'//TRIM(cv_nam)
    stypvar(1)%conline_operation = 'N/A'

    CALL SetGlobalAtt (cglobal)

    ! create output fileset
    ncout = create      (cf_out, 'none', nbin, npt, 0      )
    ierr  = createvar   (ncout,  stypvar, 1, ipk, id_varout , cdglobal=cglobal        )
    ierr  = putheadervar(ncout,  cf_ifil, nbin, npt, 0 , pnavlon=zlon, pnavlat=ztimed )

  END SUBROUTINE CreateOutput
END PROGRAM cdfpdf
