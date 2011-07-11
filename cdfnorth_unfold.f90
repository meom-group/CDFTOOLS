PROGRAM cdfnorth_unfold
  !!======================================================================
  !!                     ***  PROGRAM  cdfnorth_unfold  ***
  !!=====================================================================
  !!  ** Purpose : Unfold the arctic ocean in an ORCA like configuration
  !!               for all the variables of the file given in the arguments
  !!
  !!  ** Method  : read the filename, the limit of the extracted zone, and
  !!              the type of pivot to use and the C-grid point of variables
  !!
  !! History : 2.1  : 04/2010  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  !!
  !!   unfold     unfold the north pole of orca grid
  !!   chkisig    function to determine if the variable changes sign
  !!              when folded
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jt, jvar, jv         ! dummy loop index
  INTEGER(KIND=4)                               :: ierr                     ! working integer
  INTEGER(KIND=4)                               :: narg, iargc              ! browse line
  INTEGER(KIND=4)                               :: npiglo, npjglo, npk, npt ! size of the domain
  INTEGER(KIND=4)                               :: nvars                    ! Number of variables in a file
  INTEGER(KIND=4)                               :: ijatl, ijpacif           ! j starting position in atl and pacif
  INTEGER(KIND=4)                               :: npiarctic, npjarctic     ! size of the output arrays file
  INTEGER(KIND=4)                               :: isig                     ! change sign indicator
  INTEGER(KIND=4)                               :: nipivot                  ! i position of pivot
  INTEGER(KIND=4)                               :: ncout                    ! ncid of output file
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var                   ! arrays of var id's (input)
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk, id_varout           ! level and varid of output var

  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: tab                      ! output array
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: tablon, tablat           ! output longitude and latitude
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d                      ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: tim                      ! time counter

  CHARACTER(LEN=256)                            :: cf_in                    ! input file name
  CHARACTER(LEN=256)                            :: cf_out='unfold.nc'       ! output file names 
  CHARACTER(LEN=256)                            :: cv_dep                   ! depth name
  CHARACTER(LEN=256)                            :: cpivot                   ! pivot position
  CHARACTER(LEN=256)                            :: ctype                    ! variable position
  CHARACTER(LEN=256)                            :: cglobal                  ! variable position
  CHARACTER(LEN=256)                            :: cldum                    ! dummy string
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names                 ! array of var name

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar                  ! output var attribute

  LOGICAL                                       :: lchk=.false.             ! flag for consistency check
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg /=5 ) THEN
     PRINT *,' usage : cdfnorth_unfold IN-file jatl jpacif pivot Cgrid_point'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Unfold the Artic Ocean in an ORCA configuration. Produce a netcdf' 
     PRINT *,'       file with the Artic ocean as a whole. The area can be adjusted on'
     PRINT *,'       both Atlantic and Pacific sides.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       IN-file     : netcdf file to be unfolded.' 
     PRINT *,'       jatl        : J index to start the unfold process in the Atlantic.'
     PRINT *,'       jpacif      : J index to start the unfold process in the Pacific.'
     PRINT *,'       pivot       : type of pivot for the north fold condition ( T or F )'
     PRINT *,'       Cgrid_point : grid point where the variables in the input file are'
     PRINT *,'                     located. If all variables in a single file are not on'
     PRINT *,'                     the same C-grid location, there might be a problem ...'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : same name and units than in the input file.'
     STOP
  ENDIF

  CALL getarg (1, cf_in )
  CALL getarg (2, cldum ) ; READ(cldum,*) ijatl
  CALL getarg (3, cldum ) ; READ(cldum,*) ijpacif
  CALL getarg (4, cpivot) 
  CALL getarg (5, ctype )
  
  IF ( chkfile(cf_in) ) STOP ! missing file

  WRITE(cglobal,9000) 'cdfnorth_unfold ',TRIM(cf_in), ijatl, ijpacif, TRIM(cpivot), TRIM(ctype)
9000 FORMAT(a,a,2i5,a,1x,a)

  npiglo = getdim (cf_in, cn_x                             )
  npjglo = getdim (cf_in, cn_y                             )
  npt    = getdim (cf_in, cn_t                             )
  npk    = getdim (cf_in, cn_z, cdtrue=cv_dep, kstatus=ierr)

  IF (ierr /= 0 ) THEN
     npk   = getdim (cf_in,'z',cdtrue=cv_dep,kstatus=ierr)
     IF (ierr /= 0 ) THEN
        npk   = getdim (cf_in,'sigma',cdtrue=cv_dep,kstatus=ierr)
        IF ( ierr /= 0 ) THEN
           npk = getdim (cf_in,'nav_lev',cdtrue=cv_dep,kstatus=ierr)
           IF ( ierr /= 0 ) THEN
              PRINT *,' assume file with no depth'
              npk=0
           ENDIF
        ENDIF
     ENDIF
  ENDIF

  ! to be improved
  npiarctic = npiglo/2
  nipivot   = npiglo/2

  SELECT CASE ( cpivot )
  CASE ( 'T','t') ; npjarctic=(npjglo-ijatl+1)  + (npjglo -ijpacif +1) -3
  CASE ( 'F','f') ; npjarctic=(npjglo-ijatl+1)  + (npjglo -ijpacif +1) -2
  END SELECT

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE( tab(npiarctic, npjarctic),  v2d(npiglo,npjglo), tim(npt)   )
  ALLOCATE( tablon(npiarctic, npjarctic), tablat(npiarctic, npjarctic) )

  nvars = getnvar(cf_in)
  PRINT *,' nvars = ', nvars

  ALLOCATE (cv_names(nvars) )
  ALLOCATE (stypvar(nvars) )
  ALLOCATE (id_var(nvars), ipk(nvars), id_varout(nvars) )

  ! get list of variable names and collect attributes in stypvar (optional)
  cv_names(:) = getvarname(cf_in, nvars, stypvar)

  id_var(:)  = (/(jv, jv=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cf_in, nvars, cdep=cv_dep)

  WHERE( ipk == 0 ) cv_names = 'none'
  stypvar(:)%cname = cv_names

  v2d=getvar(cf_in, cn_vlon2d, 1, npiglo, npjglo)
  CALL unfold(v2d ,tablon, ijatl, ijpacif, cpivot, ctype, 1)

  v2d=getvar(cf_in, cn_vlat2d, 1, npiglo, npjglo)
  CALL unfold(v2d ,tablat, ijatl, ijpacif, cpivot, ctype, 1)

  ! create output file taking the sizes in cf_in
  ncout = create      (cf_out, cf_in,   npiarctic, npjarctic, npk,                                 cdep=cv_dep)
  ierr  = createvar   (ncout,  stypvar, nvars,     ipk,       id_varout, cdglobal=TRIM(cglobal)               )
  ierr  = putheadervar(ncout,  cf_in,   npiarctic, npjarctic, npk, pnavlon=tablon, pnavlat=tablat, cdep=cv_dep)

  tim  = getvar1d(cf_in, cn_vtimec, npt     )
  ierr = putvar1d(ncout, tim,       npt, 'T')

  DO jvar = 1,nvars
     PRINT *,' Working with ', TRIM(cv_names(jvar)), ipk(jvar)
     DO jk = 1, ipk(jvar)
        PRINT *,'level ',jk
        tab(:,:) = 0.
        isig =  1
        DO jt=1,npt
           v2d(:,:) = getvar(cf_in, cv_names(jvar), jk, npiglo, npjglo, ktime=jt )

           IF ( jk == 1 .AND. jt == 1) THEN ! look for correct isig
              isig=chkisig( cpivot, ctype, v2d, lchk)
              PRINT *,'ISIG=', isig
           ENDIF

           CALL unfold(v2d, tab, ijatl, ijpacif, cpivot, ctype, isig)
           ierr = putvar(ncout, id_varout(jvar), tab, jk, npiarctic, npjarctic)
        ENDDO
     END DO  ! loop to next level
  END DO ! loop to next var in file

  ierr = closeout(ncout)

CONTAINS

 
  INTEGER(KIND=4) FUNCTION chkisig (cdpivot, cdtype, ptab, ldchk)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION chkisig  ***
    !!
    !! ** Purpose : from the input data determine if the field is to be 
    !!              multiplied by -1 in the unfolding process  or not.
    !!              if ldchk is true, proceed to an extended check of the 
    !!              overlaping area (not written yet)
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),             INTENT(in) :: cdpivot, cdtype
    REAL(KIND=4), DIMENSION(:,:), INTENT(in) :: ptab
    LOGICAL,                      INTENT(in) :: ldchk
    ! 
    INTEGER(KIND=4)                          :: ii, ij
    REAL(KIND=4)                             :: zrat
    !!----------------------------------------------------------------------

    IF ( ldchk ) THEN
      PRINT *,' Full check not written yet ' ; STOP
    ELSE
    SELECT CASE ( cdpivot)
    CASE ( 'T','t')
       SELECT CASE (cdtype )
       CASE ( 'T','t') 
          ii = 1
          DO WHILE ( ptab(ii,npjglo-1)  == 0 .AND. ii < npiglo ) 
             ii = ii+1
          ENDDO
          IF ( ii /=  npiglo )  THEN
             ij = 2*nipivot - ii +2
             zrat = ptab(ij,npjglo-1) / ptab(ii,npjglo-1)
             IF ( ABS(zrat) /= 1. ) THEN
                PRINT *, 'INCOHERENT value in T point ', TRIM(cv_names(jvar)), zrat
                ierr = closeout(ncout)
                STOP
             ELSE
                chkisig = zrat
             ENDIF
          ENDIF
       CASE ( 'U','u') 
          ii = 1
          DO WHILE ( ptab(ii,npjglo-1) == 0  .AND. ii < npiglo )
             ii = ii+1
          ENDDO
          ij = 2*nipivot - ii + 1
          zrat = ptab(ij,npjglo-1) / ptab(ii,npjglo-1)
          IF ( ABS(zrat) /= 1. ) THEN
             PRINT *, 'INCOHERENT value in U point ', TRIM(cv_names(jvar)), zrat
             ierr = closeout(ncout)
             STOP
          ELSE
             chkisig=zrat
          ENDIF
       CASE ( 'V','v') 
          ii = 1
          DO WHILE ( ptab(ii,npjglo-1) == 0 .AND. ii < npiglo )
             ii = ii+1
          ENDDO
          ij = 2*nipivot - ii + 2
          zrat = ptab(ij,npjglo-2) / ptab(ii,npjglo-1)
          IF ( ABS(zrat) /= 1. ) THEN
             PRINT *, 'INCOHERENT value in V point ', TRIM(cv_names(jvar)), zrat
             ierr = closeout(ncout)
             STOP
          ELSE
             chkisig=zrat
          ENDIF
       END SELECT
    CASE ( 'F','f')
       PRINT *, 'F pivot not done yet ' ; STOP
    END SELECT
    ENDIF

  END FUNCTION chkisig


  SUBROUTINE unfold( ptabin, ptabout, kjatl, kjpacif, cdpivot, cdtype, ksig)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE unfold  ***
    !!
    !! ** Purpose : unfold the north pole
    !!
    !! -----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(npiglo,npjglo),       INTENT(in ) :: ptabin
    REAL(KIND=4), DIMENSION(npiarctic,npjarctic), INTENT(out) :: ptabout
    INTEGER(KIND=4),                              INTENT(in ) :: kjatl
    INTEGER(KIND=4),                              INTENT(in ) :: kjpacif
    CHARACTER(LEN=*),                             INTENT(in ) :: cdpivot
    CHARACTER(LEN=*),                             INTENT(in ) :: cdtype
    INTEGER(KIND=4),                              INTENT(in ) :: ksig

    INTEGER(KIND=4)                                           :: ji, jj
    INTEGER(KIND=4)                                           :: ipivot
    INTEGER(KIND=4)                                           :: ijn, ii, ij
    !! -----------------------------------------------------------------------
    !
    ipivot=npiglo/2
    DO jj=kjatl, npjglo
       ij = jj-kjatl+1
       ptabout(:,ij) = ptabin(ipivot:npiglo,jj)
    ENDDO

    ijn=ij

    SELECT CASE ( cdpivot )
    CASE ('T','t')   ! pivot
       SELECT CASE ( cdtype )
       CASE ('T','t')
          DO jj=npjglo-3,kjpacif, -1
             ij=  ijn + ( npjglo - 3  - jj ) +1 !  2 *npjglo - kjatl -1 -jj
             DO ji = 2, npiarctic
                ii = ipivot - ji + 3 
                ptabout(ji,ij) = ksig * ptabin(ii, jj)
             ENDDO
          ENDDO
       CASE ('V','v')
          DO jj=npjglo-4,kjpacif-1, -1
             ij=  ijn + ( npjglo - 4  - jj ) +1 !  2 *npjglo - kjatl -1 -jj
             DO ji = 2, npiarctic
                ii = ipivot - ji + 3
                ptabout(ji,ij) = ksig * ptabin(ii, jj)
             ENDDO
          ENDDO
       CASE ('U','u')
          DO jj=npjglo-3,kjpacif, -1
             ij=  ijn + ( npjglo - 3  - jj ) +1 !  2 *npjglo - kjatl -1 -jj
             DO ji = 1, npiarctic
                ii = ipivot -ji + 2
                ptabout(ji,ij) = ksig * ptabin(ii, jj)
             ENDDO
          ENDDO
       END SELECT
    CASE ('F','f')   ! pivot
       PRINT * , ' Not yet done for F pivot ' ; STOP
    END SELECT

  END SUBROUTINE unfold

END PROGRAM cdfnorth_unfold
