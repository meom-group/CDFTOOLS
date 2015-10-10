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
  INTEGER(KIND=4)                           :: imin, imax 
  INTEGER(KIND=4)                           :: jmin, jmax
  INTEGER(KIND=4)                           :: nbin, nlim, ibin

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zvar
  REAL(KIND=4), DIMENSION(:)  , ALLOCATABLE :: vlim, zcount
  REAL(KIND=4)                              :: vmin, vmax, bin_siz
  REAL(KIND=4)                              :: below, above


  CHARACTER(LEN=256)                        :: cf_ifil
  CHARACTER(LEN=256)                        :: cf_out
  CHARACTER(LEN=256)                        :: cv_nam
  CHARACTER(LEN=256)                        :: cldum

  LOGICAL                                   :: lchk
  LOGICAL                                   :: l_nozoom=.true.
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfpdf -f IN-file -v IN-var [-zoom imin imax jmin jmax] ..'
     PRINT *,'            [-lev level ] [-range vmin vmax nbin ] '
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
     PRINT *,'       [-range vmin vmax nbin  ] : define the limit for binning '
     PRINT *,'                               and number of bins.'
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

  ! Allocate memory
  ALLOCATE ( zvar(npi,npj), zcount(nbin), vlim(nlim) )

  DO ji=1, nlim 
     vlim(ji)= vmin + (ji-1)*bin_siz
  ENDDO

  DO jt = 1, npt
      zcount(:) = 0. ; below=0. ; above=0.
      zvar(:,:) = getvar( cf_ifil, cv_nam, npi, npj, ilev, ktime=jt, kimin=imin, kjmin=jmin )
      DO jj=1,npj
        DO ji=1,npi
           ibin=INT((zvar(ji,ji)-vmin)/bin_siz) +1
           IF ( ibin < 1 ) THEN 
              below=below+1
           ELSE IF ( ibin > nbin ) THEN
              above=above+1
           ELSE
              zcount(ibin) = zcount(ibin) + 1
           ENDIF
        ENDDO
      ENDDO
      PRINT *, ' JT = ', jt
      print *, 'below', below
      DO ji=1, nbin
        print *, zcount (ji )
      ENDDO
      PRINT *, 'above', above
  ENDDO

END PROGRAM cdfpdf
