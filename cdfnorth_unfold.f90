PROGRAM cdfnorth_unfold
  !!-----------------------------------------------------------------------
  !!                 ***  PROGRAM cdfnorth_unfold  ***
  !!
  !!  **  Purpose: Unfold the arctic ocean in an ORCA like configuration
  !!               for all the variables of the file given in the arguments
  !!  
  !!  **  Method: read the filename, the limit of the extracted zone, and
  !!              the type of pivot to use and the C-grid point of variables
  !!
  !! history :
  !!     Original code :   J.M. Molines (Apr. 2010 )
  !!-------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !!
  USE cdfio 

  IMPLICIT NONE
  INTEGER   :: jk,jt,jvar, jv , jtt,jkk                     !: dummy loop index
  INTEGER   :: ji, ij                                       !: dummy loop index
  INTEGER   :: ierr                                         !: working integer
  INTEGER   :: narg, iargc                                  !: 
  INTEGER   :: npiglo,npjglo, npk ,nt                       !: size of the domain
  INTEGER   :: nvars                                        !: Number of variables in a file
  INTEGER   :: ijatl, ijpacif, npiarctic, npjarctic, isig
  INTEGER   :: ipivot
  INTEGER , DIMENSION(:), ALLOCATABLE :: id_var , &         !: arrays of var id's
       &                             ipk    , &         !: arrays of vertical level for each var
       &                             id_varout
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: tab 
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: tablon,tablat 
  REAL(KIND=4)                                :: zrat
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: v2d        !: Array to read a layer of data
  REAL(KIND=4),DIMENSION(:), ALLOCATABLE      :: tim, gdep

  CHARACTER(LEN=256) :: cfile ,cfileout                      !: file name
  CHARACTER(LEN=256) ::  cdep, cdum, cpivot, ctype
  CHARACTER(LEN=256) ,DIMENSION(:), ALLOCATABLE:: cvarname   !: array of var name
  
  TYPE (variable), DIMENSION(:), ALLOCATABLE :: typvar

  INTEGER    :: ncout
  INTEGER    :: istatus
  LOGICAL    :: lcaltmean

  !!
  !!  Read command line
  narg= iargc()
  IF ( narg /= 5 ) THEN
     PRINT *,' Usage : cdfnorth_unfold filename jatl jpacif pivot Cgrid_point'
     PRINT *, '    example: cdfnorth_unfold ORCA025-G70_y2000m10d02_gridT.nc 766 766 T T'
     PRINT *, '    a file named unfold.nc will be created '
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cfile)
  CALL getarg (2, cdum) ; READ(cdum,*) ijatl
  CALL getarg (3, cdum) ; READ(cdum,*) ijpacif
  CALL getarg (4, cpivot) 
  CALL getarg (5, ctype )


  npiglo= getdim (cfile,'x')
  npjglo= getdim (cfile,'y')
  npk   = getdim (cfile,'depth',cdtrue=cdep, kstatus=istatus)
  nt   = getdim (cfile,'time', kstatus=istatus)

  IF (istatus /= 0 ) THEN
     npk   = getdim (cfile,'z',cdtrue=cdep,kstatus=istatus)
     IF (istatus /= 0 ) THEN
       npk   = getdim (cfile,'sigma',cdtrue=cdep,kstatus=istatus)
        IF ( istatus /= 0 ) THEN 
          npk = getdim (cfile,'nav_lev',cdtrue=cdep,kstatus=istatus)
            IF ( istatus /= 0 ) THEN 
              PRINT *,' assume file with no depth'
              npk=0
            ENDIF
        ENDIF
     ENDIF
  ENDIF

  ! to be improved
  npiarctic=npiglo/2
  ipivot=npiglo/2

  SELECT CASE ( cpivot )
  CASE ( 'T','t') ; npjarctic=(npjglo-ijatl+1)  + (npjglo -ijpacif +1) -3 
  CASE ( 'F','f') ; npjarctic=(npjglo-ijatl+1)  + (npjglo -ijpacif +1) -2 
  END SELECT
  

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'nt    =', nt

  ALLOCATE( tab(npiarctic, npjarctic),  v2d(npiglo,npjglo), tim(nt), gdep(npk) )
  ALLOCATE( tablon(npiarctic, npjarctic), tablat(npiarctic, npjarctic) )

  nvars = getnvar(cfile)
  PRINT *,' nvars =', nvars

  ALLOCATE (cvarname(nvars) )
  ALLOCATE (typvar(nvars) )
  ALLOCATE (id_var(nvars),ipk(nvars),id_varout(nvars) )

  ! get list of variable names and collect attributes in typvar (optional)
  cvarname(:)=getvarname(cfile,nvars,typvar)

  id_var(:)  = (/(jv, jv=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cfile,nvars,cdep=cdep)
  WHERE( ipk == 0 ) cvarname='none'
  typvar(:)%name=cvarname

  v2d=getvar(cfile, 'nav_lon',1, npiglo,npjglo)
  CALL unfold(v2d ,tablon, ijatl, ijpacif, cpivot, ctype, 1)
  v2d=getvar(cfile, 'nav_lat',1, npiglo,npjglo)
  CALL unfold(v2d ,tablat, ijatl, ijpacif, cpivot, ctype, 1)

  ! create output fileset
  cfileout='unfold.nc'
  ! create output file taking the sizes in cfile

  ncout =create(cfileout, cfile,npiarctic,npjarctic,npk,cdep=cdep)
  ierr= createvar(ncout , typvar,  nvars, ipk, id_varout )
  tim=getvar1d(cfile,'time_counter',nt)
! gdep=getvar1d(cfile,cdep,npk)
  
  ierr= putheadervar(ncout , cfile, npiarctic,npjarctic, npk,pnavlon=tablon, pnavlat=tablat, cdep=cdep)
  ierr=putvar1d(ncout,tim,nt,'T')
! ierr=putvar1d(ncout,gdep,npk,'D')

  DO jvar = 1,nvars
        PRINT *,' Working with ', TRIM(cvarname(jvar)), ipk(jvar)
        DO jk = 1, ipk(jvar)
           PRINT *,'level ',jk
           tab(:,:) = 0.
           isig =  1
              DO jtt=1,nt
                jkk=jk
                ! If forcing fields is without depth dimension
                IF (npk==0) jkk=jtt 
                v2d(:,:)= getvar(cfile, cvarname(jvar), jkk ,npiglo, npjglo,ktime=jtt )
                IF ( jk == 1 ) THEN ! look for correct isig
                 SELECT CASE ( cpivot)
                 CASE ( 'T','t')
                   SELECT CASE (ctype )
                   CASE ( 'T','t') 
                     ji=1
                     DO WHILE ( v2d(ji,npjglo-1)  == 0 .AND. ji < npiglo ) 
                       ji=ji+1
                     ENDDO
                      IF ( ji /=  npiglo )  THEN
                      ij=2*ipivot - ji +2
                      zrat= v2d(ij,npjglo-1) / v2d(ji,npjglo-1)
                      IF ( ABS(zrat) /= 1. ) THEN
                        PRINT *, 'INCOHERENT value in T point '; stop
                      ELSE
                       isig=zrat
                      ENDIF
                      ENDIF
                   CASE ( 'U','u') 
                     ji=1
                     DO WHILE ( v2d(ji,npjglo-1) == 0  .AND. ji < npiglo )
                       ji=ji+1
                     ENDDO
                      ij=2*ipivot - ji + 1
                      zrat= v2d(ij,npjglo-1) / v2d(ji,npjglo-1)
                      IF ( ABS(zrat) /= 1. ) THEN
                        PRINT *, 'INCOHERENT value in U point '; stop
                      ELSE
                       isig=zrat
                      ENDIF
                   CASE ( 'V','v') 
                     ji=1
                     DO WHILE ( v2d(ji,npjglo-1) == 0 .AND. ji < npiglo )
                       ji=ji+1
                     ENDDO
                      ij=2*ipivot - ji + 2
                      zrat= v2d(ij,npjglo-2) / v2d(ji,npjglo-1)
                      IF ( ABS(zrat) /= 1. ) THEN
                        PRINT *, 'INCOHERENT value in V point '; stop
                      ELSE
                       isig=zrat
                      ENDIF
                   END SELECT
                 CASE ( 'F','f')
                 END SELECT
                PRINT *,'ISIG=', isig
                ENDIF
                
                CALL unfold(v2d, tab, ijatl, ijpacif, cpivot, ctype, isig)
                ierr = putvar(ncout, id_varout(jvar) ,tab, jkk, npiarctic, npjarctic)
              ENDDO
        END DO  ! loop to next level
  END DO ! loop to next var in file

  istatus = closeout(ncout)

CONTAINS
  SUBROUTINE unfold( ptabin, ptabout, kjatl, kjpacif, cdpivot, cdtype, ksig)
    !!------------------------------------------------------------------------
    !!            ** SUBROUTINE unfol **
    !!
    !!   Purpose : unfold the north pole 
    !! -----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(npiglo,npjglo)      , INTENT(in)  :: ptabin
    REAL(KIND=4), DIMENSION(npiarctic,npjarctic), INTENT(out) :: ptabout
    INTEGER, INTENT(in)   ::  kjatl
    INTEGER, INTENT(in)   ::  kjpacif
    INTEGER, INTENT(in)   ::  ksig
    CHARACTER(LEN=*), INTENT(in) :: cdpivot
    CHARACTER(LEN=*), INTENT(in) :: cdtype
    !!
    ! local variables :
    INTEGER :: jj,  ipivot, ij, ijn, ji, ii
    !
    ipivot=npiglo/2
    DO jj=kjatl, npjglo
      ij=jj-kjatl+1
      ptabout(:,ij) = ptabin (ipivot:npiglo,jj)
    ENDDO
    ijn=ij
    SELECT CASE ( cdpivot )
    CASE ('T','t')   ! pivot
    SELECT CASE ( cdtype)
    CASE ('T','t')
      DO jj=npjglo-3,kjpacif, -1
        ij=  ijn + ( npjglo - 3  - jj ) +1 !  2 *npjglo - kjatl -1 -jj
        DO ji = 2, npiarctic
!        ii = 2*ipivot -ji +2 -ipivot +1 
         ii = ipivot - ji + 3 
         ptabout(ji,ij)= ksig * ptabin(ii, jj)
        ENDDO
      ENDDO
    CASE ('V','v')
      DO jj=npjglo-4,kjpacif-1, -1
        ij=  ijn + ( npjglo - 4  - jj ) +1 !  2 *npjglo - kjatl -1 -jj
        DO ji = 2, npiarctic
!        ii = 2*ipivot -ji +2 -ipivot +1
         ii = ipivot - ji + 3
         ptabout(ji,ij)= ksig * ptabin(ii, jj)
        ENDDO
      ENDDO
    CASE ('U','u')
      DO jj=npjglo-3,kjpacif, -1
        ij=  ijn + ( npjglo - 3  - jj ) +1 !  2 *npjglo - kjatl -1 -jj
        DO ji = 1, npiarctic
!        ii = 2*ipivot -ji + 1 -ipivot + 1
         ii = ipivot -ji + 2
         ptabout(ji,ij)= ksig * ptabin(ii, jj)
        ENDDO
      ENDDO
    END SELECT
    CASE ('F','f')   ! pivot
     PRINT * , ' Not yet done for F pivot ' ; stop
    END SELECT

  END SUBROUTINE unfold

END PROGRAM cdfnorth_unfold
