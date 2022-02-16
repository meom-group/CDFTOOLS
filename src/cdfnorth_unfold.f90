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
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!
  !!   unfold     unfold the north pole of orca grid
  !!   chkisig    function to determine if the variable changes sign
  !!              when folded
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class data_transformation
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jt, jvar, jv         ! dummy loop index
  INTEGER(KIND=4)                               :: ierr                     ! working integer
  INTEGER(KIND=4)                               :: idep, idep_max           ! possible depth index, maximum
  INTEGER(KIND=4)                               :: narg, iargc, ijarg       ! browse line
  INTEGER(KIND=4)                               :: npiglo, npjglo, npk, npt ! size of the domain
  INTEGER(KIND=4)                               :: nvars                    ! Number of variables in a file
  INTEGER(KIND=4)                               :: ijatl, ijpacif           ! j starting position in atl and pacif
  INTEGER(KIND=4)                               :: npiarctic, npjarctic     ! size of the output arrays file
  INTEGER(KIND=4)                               :: isig                     ! change sign indicator
  INTEGER(KIND=4)                               :: nipivot                  ! i position of pivot
  INTEGER(KIND=4)                               :: ncout                    ! ncid of output file
  INTEGER(KIND=4)                               :: nvaro                    ! ncid of output file
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var                   ! arrays of var id's (input)
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk, id_varout           ! level and varid of output var

  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: tab                      ! output array
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: tablon, tablat           ! output longitude and latitude
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d                      ! Array to read a layer of data

  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dtim                     ! time counter

  CHARACTER(LEN=256)                            :: cf_in                    ! input file name
  CHARACTER(LEN=256)                            :: cf_out='unfold.nc'       ! output file names 
  CHARACTER(LEN=256)                            :: cv_dep                   ! depth name
  CHARACTER(LEN=256)                            :: cpivot                   ! pivot position
  CHARACTER(LEN=256)                            :: ctype                    ! variable position
  CHARACTER(LEN=256)                            :: cglobal                  ! variable position
  CHARACTER(LEN=256)                            :: cldum                    ! dummy string
  CHARACTER(LEN=80 ), DIMENSION(:), ALLOCATABLE :: cv_in                    ! array of var name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names                 ! array of var name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: clv_dep                  ! array of possible depth name (or 3rd dimension)

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar                  ! output var attribute

  LOGICAL, DIMENSION(:)           , ALLOCATABLE :: ll_keep
  LOGICAL                                       :: lchk = .false.           ! flag for consistency check
  LOGICAL                                       :: lnc4 = .FALSE.           ! Use nc4 with chunking and deflation
  LOGICAL                                       :: llst = .FALSE.           ! Use selected list of variable
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfnorth_unfold -f IN-file -jatl jatl -jpacif jpacif -piv pivot ...'
     PRINT *,'              ... -p C-type [-v VAR-list] [-tdim TIME-dim] [-zdim Z-dim] ...'
     PRINT *,'      [-tvar TIME-var] [-zvar Z-var] [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Unfolds the Artic Ocean in an ORCA configuration. Produce a netcdf' 
     PRINT *,'       file with the Artic ocean as a whole. The area can be adjusted on'
     PRINT *,'       both Atlantic and Pacific sides.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file     : Input netcdf file to be unfolded.' 
     PRINT *,'       -jatl jatl     : J index to start the unfold process in the Atlantic.'
     PRINT *,'       -jpacif jpacif : J index to start the unfold process in the Pacific.'
     PRINT *,'       -piv pivot     : type of pivot for the north fold condition ( T or F )'
     PRINT *,'             ORCA1, ORCA05 use F-pivot, ORCA2, ORCA025, ORCA12 use T-pivot.'
     PRINT *,'       -p C-type : one of T|U|V|W|F , indicating the grid point where the'
     PRINT *,'             variables in the input file are located. If all variables in a '
     PRINT *,'             single file are not on the same C-grid location, there might be'
     PRINT *,'             a problem ...'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-v VAR-list] : give comma separated list of variable to unfold'
     PRINT *,'              Default : all variables in file.'
     PRINT *,'       [-tdim TIME-dim] : give the name of the time dimension if not '
     PRINT *,'               NEMO standard (time_counter)'
     PRINT *,'       [-zdim Z-dim] : give the name of the Z  dimension if not '
     PRINT *,'               NEMO standard (deptht, depthu ... ))'
     PRINT *,'       [-tvar TIME-var] : give the name of the time variable if not '
     PRINT *,'               NEMO standard (time_counter)'
     PRINT *,'       [-zvar TIME-var] : give the name of the Z variable if not '
     PRINT *,'               NEMO standard (deptht, depthu ...)'
     PRINT *,'       [-o OUT-file] : Specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ]       : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option is used.'
     PRINT *,'         variables : same name and units than in the input file.'
     PRINT *,'      '
     STOP 
  ENDIF
  
  ijarg=1
  DO WHILE ( ijarg <= narg )
      CALL getarg (ijarg, cldum )  ; ijarg=ijarg+1
      SELECT CASE ( cldum )
      CASE ( '-f'     ) ; CALL getarg (ijarg, cf_in ) ; ijarg=ijarg+1
      CASE ( '-jatl'  ) ; CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1 ;  READ(cldum,*) ijatl
      CASE ( '-jpacif') ; CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1 ;  READ(cldum,*) ijpacif
      CASE ( '-piv'   ) ; CALL getarg (ijarg, cpivot) ; ijarg=ijarg+1
      CASE ( '-p'     ) ; CALL getarg (ijarg, ctype ) ; ijarg=ijarg+1
      ! options
      CASE ( '-v'     ) ; CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1 ; CALL ParseVars(cldum) ; llst = .true.
      CASE ( '-tdim'  ) ; CALL getarg (ijarg, cn_t  )     ; ijarg=ijarg+1 ;
      CASE ( '-zdim'  ) ; CALL getarg (ijarg, cn_z  )     ; ijarg=ijarg+1 ;
      CASE ( '-tvar'  ) ; CALL getarg (ijarg, cn_vtimec ) ; ijarg=ijarg+1 ;
      CASE ( '-zvar'  ) ; CALL getarg (ijarg, cn_vdeptht) ; ijarg=ijarg+1 ;
      CASE ( '-o'     ) ; CALL getarg (ijarg, cf_out) ; ijarg=ijarg+1
      CASE ( '-nc4'   ) ; lnc4 = .TRUE.
      CASE DEFAULT      ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
      END SELECT
  ENDDO

  IF ( chkfile(cf_in) ) STOP 99 ! missing file

  WRITE(cglobal,9000) 'cdfnorth_unfold ',TRIM(cf_in), ijatl, ijpacif, TRIM(cpivot), TRIM(ctype)
9000 FORMAT(a,a,2i5,a,1x,a)

  npiglo = getdim (cf_in, cn_x                             )
  npjglo = getdim (cf_in, cn_y                             )
  npt    = getdim (cf_in, cn_t                             )

  ! looking for npk among various possible name
  idep_max=8
  ALLOCATE ( clv_dep(idep_max) )
  clv_dep(:) = (/cn_z,'z','sigma','nav_lev','levels','ncatice','icbcla','icbsect'/)
  idep=1  ; ierr=1000
  DO WHILE ( ierr /= 0 .AND. idep <= idep_max )
     npk  = getdim (cf_in, clv_dep(idep), cdtrue=cv_dep, kstatus=ierr)
     idep = idep + 1
  ENDDO

  IF ( ierr /= 0 ) THEN  ! none of the dim name was found
      PRINT *,' assume file with no depth'
      npk=0
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
  PRINT *, 'npt    = ', npt
  PRINT *, 'npk    = ', npk , 'Dep name :' , TRIM(cv_dep)

  ALLOCATE( tab(npiarctic, npjarctic),  v2d(npiglo,npjglo), dtim(npt)  )
  ALLOCATE( tablon(npiarctic, npjarctic), tablat(npiarctic, npjarctic) )

  nvars = getnvar(cf_in)
  PRINT *,' nvars = ', nvars

  ALLOCATE (cv_names(nvars) )
  ALLOCATE (stypvar(nvars) )
  ALLOCATE (id_var(nvars), ipk(nvars), id_varout(nvars) )
  ALLOCATE (ll_keep(nvars) )

  ! get list of variable names and collect attributes in stypvar (optional)
  cv_names(:) = getvarname(cf_in, nvars, stypvar)
  id_var(:)  = (/(jv, jv=1,nvars)/)
  IF ( llst ) THEN
    ll_keep(:)=.false.
    ! set cvnames to 'none' if not in cv_in
    DO jv=1,nvaro
      DO jvar=1,nvars
         IF ( cv_names(jvar) == cv_in(jv) ) THEN
            ll_keep(jvar) = .true.
         ENDIF
      ENDDO
    ENDDO
    WHERE (.not. ll_keep) cv_names='none'
  ENDIF
           

  CALL CreateOutput

  ! default value
  isig =  1

  DO jvar = 1,nvars
     IF (cv_names(jvar) /= 'none' ) THEN
     PRINT *,' Working with ', TRIM(cv_names(jvar)), ipk(jvar)
     DO jk = 1, ipk(jvar)
        PRINT *,'level ',jk
        tab(:,:) = 0.
        tab(:,:) = -9999.
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
     ENDIF
  END DO ! loop to next var in file

  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE ParseVars (cdum)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ParseVars  ***
    !!
    !! ** Purpose :  Decode variable name  option from command line
    !!
    !! ** Method  :  look for , in the argument string and set the number of
    !!         variable (nvaro), allocate cv_in array and fill it with the
    !!         decoded  names.
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdum

    CHARACTER(LEN=80), DIMENSION(100) :: cl_dum  ! 100 is arbitrary
    INTEGER  :: ji
    INTEGER  :: inchar,  i1=1
    !!----------------------------------------------------------------------
    PRINT *, "PARSE : ", trim(cdum)
    inchar= LEN(TRIM(cdum))
    ! scan the input string and look for ',' as separator
    nvaro=1
    DO ji=1,inchar
       IF ( cdum(ji:ji) == ',' ) THEN
          cl_dum(nvaro) = cdum(i1:ji-1)
          i1=ji+1
          nvaro=nvaro+1
       ENDIF
    ENDDO

    ! last name of the list does not have a ','
    cl_dum(nvaro) = cdum(i1:inchar)

    ALLOCATE ( cv_in(nvaro) )
    DO ji=1, nvaro
       cv_in(ji) = cl_dum(ji)
       PRINT *, ji,TRIM(cv_in(ji))
    ENDDO
  END SUBROUTINE ParseVars

 
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
      PRINT *,' Full check not written yet ' ; STOP 99
    ELSE
    SELECT CASE ( cdpivot)
    CASE ( 'T','t')
       SELECT CASE (cdtype )
       CASE ( 'T','t') 
          ii = 10  
          DO WHILE ( ptab(ii,npjglo-1)  == 0 .AND. ii < npiglo ) 
             ii = ii+1
          ENDDO
          IF ( ii /=  npiglo )  THEN
             ij = 2*nipivot - ii +2
             zrat = ptab(ij,npjglo-1) / ptab(ii,npjglo-1)
             IF ( ABS(zrat) /= 1. ) THEN
                PRINT *, 'INCOHERENT value in T point ', TRIM(cv_names(jvar)), zrat
                ierr = closeout(ncout)
                STOP 99
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
             STOP 99
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
             STOP 99
          ELSE
             chkisig=zrat
          ENDIF
       CASE ( 'F','f' ) 
          ii = 1
          DO WHILE ( ptab(ii,npjglo-1) == 0 .AND. ii < npiglo )
             ii = ii+1
          ENDDO
          ij = 2*nipivot - ii + 1
          zrat = ptab(ij,npjglo-2) / ptab(ii,npjglo-1)
          IF ( ABS(zrat) /= 1. ) THEN
             PRINT *, 'INCOHERENT value in V point ', TRIM(cv_names(jvar)), zrat
             ierr = closeout(ncout)
             STOP 99
          ELSE
             chkisig=zrat
          ENDIF
       END SELECT
    CASE ( 'F','f')
       PRINT *, 'F pivot not done yet ' ; STOP 99
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
             DO ji = 1, npiarctic
                ii = ipivot - ji + 3 
                ptabout(ji,ij) = ksig * ptabin(ii, jj)
             ENDDO
          ENDDO
       CASE ('V','v')
          DO jj=npjglo-4,kjpacif-1, -1
             ij=  ijn + ( npjglo - 4  - jj ) +1 !  2 *npjglo - kjatl -1 -jj
             DO ji = 1, npiarctic
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
       CASE ('F','f' )
          DO jj=npjglo-4,kjpacif-1, -1
             ij=  ijn + ( npjglo - 4  - jj ) +1 
             DO ji = 1, npiarctic
                ii = ipivot - ji + 2
                ptabout(ji,ij) = ksig * ptabin(ii, jj)
             ENDDO
          ENDDO
       END SELECT
    CASE ('F','f')   ! pivot
       PRINT * , ' Not yet done for F pivot ' ; STOP 99
    END SELECT

  END SUBROUTINE unfold

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cf_in, nvars, cdep=cv_dep)

  WHERE( ipk == 0 ) cv_names = 'none'
  stypvar(:)%cname = cv_names

  DO jvar = 1, nvars
     stypvar(jvar)%ichunk = (/npiarctic,MAX(1,npjarctic/30),1,1 /)
  ENDDO

  v2d=getvar(cf_in, cn_vlon2d, 1, npiglo, npjglo)
  CALL unfold(v2d ,tablon, ijatl, ijpacif, cpivot, ctype, 1)

  v2d=getvar(cf_in, cn_vlat2d, 1, npiglo, npjglo)
  CALL unfold(v2d ,tablat, ijatl, ijpacif, cpivot, ctype, 1)

  ! create output file taking the sizes in cf_in
  ncout = create      (cf_out, cf_in,   npiarctic, npjarctic, npk,                                 cdep=cv_dep, ld_nc4=lnc4)
  ierr  = createvar   (ncout,  stypvar, nvars,     ipk,       id_varout, cdglobal=TRIM(cglobal)               , ld_nc4=lnc4)
  ierr  = putheadervar(ncout,  cf_in,   npiarctic, npjarctic, npk, pnavlon=tablon, pnavlat=tablat, cdep=cv_dep)

  dtim = getvar1d(cf_in, cn_vtimec, npt     )
  ierr = putvar1d(ncout, dtim,      npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfnorth_unfold
