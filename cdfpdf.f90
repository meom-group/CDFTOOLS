PROGRAM cdfpdf
  !!======================================================================
  !!                     ***  PROGRAM  cdfpdf  ***
  !!=====================================================================
  !!  ** Purpose : Build the pdf for a given variable on a given area.
  !!
  !!  ** Method  : Establish a data binning and count the number of element
  !!               in each bin.  Binning can be defined either using a
  !!               minimum value + the size of a bin and the number of bin
  !!               or using a minimum value, maximum value, and bin size
  !!
  !! History : 3.0  : 10/2015  : J.M. Molines : Original code
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
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
  REAL(KIND=4), DIMENSION(:)  , ALLOCATABLE :: vlim, ztimes
  REAL(KIND=4)                              :: vmin, vmax, bin_siz
  REAL(KIND=4)                              :: below, above

  CHARACTER(LEN=256)                        :: cf_ifil
  CHARACTER(LEN=256)                        :: cf_asc='pdf.txt'
  CHARACTER(LEN=256)                        :: cf_out='pdf.nc'
  CHARACTER(LEN=256)                        :: cv_nam
  CHARACTER(LEN=256)                        :: cldum
  CHARACTER(LEN=256)                        :: cglobal

  TYPE(variable), DIMENSION(1)              :: stypvar           ! output data structure

  LOGICAL                                   :: lchk
  LOGICAL                                   :: l_nozoom=.true.
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfpdf -f IN-file -v IN-var [-zoom imin imax jmin jmax] ..'
     PRINT *,'       [-lev level ] [-range vmin vmax nbin ] [-o OUT-ncfile] [-a OUT-ascfile]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file : input file '
     PRINT *,'       -v IN-var  : variable name '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-zoom imin imax jmin jmax] : define a sub-area, in model '
     PRINT *,'                                     coordinates' 
     PRINT *,'       [ -lev level ] : choose a level for pdf computation '
     PRINT *,'       [-range vmin vmax nbin  ] : define the limit for binning '
     PRINT *,'                               and number of bins.'
     PRINT *,'       [-o OUT-file] : specify name for netcdf output file'
     PRINT *,'       [-a ASC-file] : specify name for ascii output file'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     STOP
  ENDIF

  ijarg = 1
  ilev  = 0
  l_nozoom = .true.
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_ifil) ; ijarg=ijarg+1
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out ) ; ijarg=ijarg+1
     CASE ( '-a'   ) ; CALL getarg(ijarg, cf_asc ) ; ijarg=ijarg+1
     CASE ( '-v'   ) ; CALL getarg(ijarg, cv_nam ) ; ijarg=ijarg+1
     CASE ( '-zoom') ; l_nozoom = .false.
             CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) imin
             CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) imax
             CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) jmin
             CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) jmax
     CASE ( '-lev') 
             CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) ilev
     CASE ( '-range') 
             CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) vmin
             CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) vmax
             CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) nbin
     CASE DEFAULT
             PRINT *,' option ',TRIM(cldum),' not understood'
     END SELECT
  ENDDO

  IF (  chkfile ( cf_ifil) ) STOP  ! some compulsory files are missing

  ! set domain size from input ile
  npiglo = getdim (cf_ifil,cn_x)
  npjglo = getdim (cf_ifil,cn_y)
  npk    = getdim (cf_ifil,cn_z)
  npt    = getdim (cf_ifil,cn_t)

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt
  
  IF ( l_nozoom ) THEN
    imin=1 ; imax=npiglo
    jmin=1 ; jmax=npjglo
  ENDIF

  npi=imax - imin + 1
  npj=jmax - jmin + 1

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
  ALLOCATE ( zvar(npi,npj), zcount(nbin,npt), vlim(nlim) )
  ALLOCATE ( ztimes(npt), ztimed(nbin,npt) , zlon(nbin,npt))

  DO ji=1, nlim 
     vlim(ji)= vmin + (ji-1)*bin_siz
  ENDDO

  OPEN (numout, FILE=cf_asc)   ! this file can be plotted easily with graph
  ! time in seconds read from file 
  ztimes  = getvar1d(cf_ifil, cn_vtimec, npt )

  ! convert in time in days since the begining of the file
  ! this will be the dummy 'latitude' for the output file
 
  DO jt=1, npt
    ztimed(:,jt) = (ztimes(jt) - ztimes(1) ) / 86400.
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

      PRINT *, ' JT = ', jt

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
