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
  INTEGER   :: ierr                                         !: working integer
  INTEGER   :: narg, iargc                                  !: 
  INTEGER   :: npiglo,npjglo, npk ,nt                       !: size of the domain
  INTEGER   :: nvars                                        !: Number of variables in a file
  INTEGER   :: ijatl, ijpacif, npiarctic, npjarctic, isig
  INTEGER , DIMENSION(:), ALLOCATABLE :: id_var , &         !: arrays of var id's
       &                             ipk    , &         !: arrays of vertical level for each var
       &                             id_varout
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: tab  !: Arrays for cumulated values
  REAL(KIND=4)                                :: total_time
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: v2d        !: Array to read a layer of data
  REAL(KIND=4),DIMENSION(1)                   :: timean
  REAL(KIND=4),DIMENSION(365)                 ::  tim

  CHARACTER(LEN=256) :: cfile ,cfileout                      !: file name
  CHARACTER(LEN=256) ::  cdep, cdum, cpivot, ctype
  CHARACTER(LEN=256) ,DIMENSION(:), ALLOCATABLE:: cvarname   !: array of var name
  
  TYPE (variable), DIMENSION(:), ALLOCATABLE :: typvar

  INTEGER    :: ncout, ncout2
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
  CALL getarg (4, ctype )

  ! to be improved
  SELECT CASE ( ctype ) 
   CASE ( 'T','t') 
     isig=1
   CASE ('U','u','V','v')
     isig=-1
  END SELECT

  npiglo= getdim (cfile,'x')
  npjglo= getdim (cfile,'y')
  npk   = getdim (cfile,'depth',cdtrue=cdep, kstatus=istatus)

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
  npjarctic=npjglo-ijatl + npjglo -ijpacif 
  

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( tab(npiarctic, npjarctic),  v2d(npiglo,npjglo) )

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

  ! create output fileset
  cfileout='unfold.nc'
  ! create output file taking the sizes in cfile

  ncout =create(cfileout, cfile,npiarctic,npjarctic,npk,cdep=cdep)
  ierr= createvar(ncout , typvar,  nvars, ipk, id_varout )
  
!  ierr= putheadervar(ncout , cfile, npiarctic,npjarctic, npk,cdep=cdep)
  ierr=putvar1d(ncout,timean,1,'T')

  DO jvar = 1,nvars
        PRINT *,' Working with ', TRIM(cvarname(jvar)), ipk(jvar)
        DO jk = 1, ipk(jvar)
           PRINT *,'level ',jk
           tab(:,:) = 0.
              DO jtt=1,nt
                jkk=jk
                ! If forcing fields is without depth dimension
                IF (npk==0) jkk=jtt 
                v2d(:,:)= getvar(cfile, cvarname(jvar), jkk ,npiglo, npjglo,ktime=jtt )
                CALL unfold(v2d, tab, ijatl, ijpacif, cpivot, ctype, isig)
                ierr = putvar(ncout, id_varout(jvar) ,tab, jkk, npiarctic, npjarctic, ktime=jtt)
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
    INTEGER :: jj, jnorth, ipivot, ij
    !
    jnorth=npjglo -1

    ipivot=npiglo/2
    SELECT CASE ( cdtype)
    CASE ('T','t','V','v')
      DO jj=kjatl,npjglo
        jout=jj-kjatl+1
        ptabout(:,jout)=ptabin(ipivot:npiglo, jj)
      ENDDO

    DO jj=jnorth-1, kjpacif,-1
       ij= 1
       ptabout(:,ij)= ksig * ptabin(1:npiglo/2, jj)
    ENDDO

  END SUBROUTINE unfold

    
!   SUBROUTINE lbc_nfd_2d( pt2d, cd_type, psgn )
  SUBROUTINE lbc_nfd_2d( ptabin, ptabout, kjatl, kjpacif, cdpivot, cdtype, psgn)
      !!----------------------------------------------------------------------
      !!                  ***  routine lbc_nfd_2d  ***
      !!
      !! ** Purpose :   2D lateral boundary condition : North fold treatment
      !!       without processor exchanges. 
      !!
      !! ** Method  :   
      !!
      !! ** Action  :   pt2d with update value at its periphery
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER(len=1) , INTENT( in ) ::   &
         cdtype       ! define the nature of ptab array grid-points
      !             ! = T , U , V , F , W points
      !             ! = S : T-point, north fold treatment ???
      !             ! = G : F-point, north fold treatment ???
      REAL(wp), INTENT( in ) ::   &
         psgn          ! control of the sign change
      !             !   = -1. , the sign is changed if north fold boundary
      !             !   =  1. , the sign is kept  if north fold boundary
      REAL(wp), DIMENSION(:,:), INTENT( inout ) ::   &
         pt2d          ! 3D array on which the boundary condition is applied

      !! * Local declarations
      INTEGER  ::   ji, jl
      INTEGER  ::   ijt, iju, ijpj, ijpjm1

      ijpj = 4

      ijpjm1 = ijpj-1



      SELECT CASE ( cdpivot )

      CASE ( 'T','t' )                       ! *  North fold  T-point pivot
      jnorth= npjglo -1
      DO jj=kjatl,jnorth
        ptabout(:,jj)=ptabin(npiglo/2:npiglo, jj)
      ENDDO
      DO jj=jnorth, kjpacif,-1
        ij= 1
        ptabout(:,ij)= ksig * ptabin(1:npiglo/2, jj)
      ENDDO

         SELECT CASE ( cdtype )

         CASE ( 'T', 'S', 'W' )
            DO ji = 2, npiglo
               ijt=npiglo-ji+2
               pt2d(ji,ijpj) = psgn * pt2d(ijt,ijpj-2)
            END DO
            DO ji = npiglo/2+1, npiglo
               ijt=npiglo-ji+2
               pt2d(ji,ijpj-1) = psgn * pt2d(ijt,ijpj-1)
            END DO
         CASE ( 'U' )                                     ! U-point
            DO ji = 1, npiglo-1
               iju = npiglo-ji+1
               pt2d(ji,ijpj) = psgn * pt2d(iju,ijpj-2)
            END DO
            DO ji = npiglo/2, npiglo-1
               iju = npiglo-ji+1
               pt2d(ji,ijpjm1) = psgn * pt2d(iju,ijpjm1)
            END DO
         CASE ( 'V' )                                     ! V-point
            DO jl =-1, 0
               DO ji = 2, npiglo
                  ijt = npiglo-ji+2
                  pt2d(ji,ijpj+jl) = psgn * pt2d(ijt,ijpj-3-jl)
               END DO
            END DO
         CASE ( 'F' , 'G' )                               ! F-point
            DO jl =-1, 0
               DO ji = 1, npiglo-1
                  iju = npiglo-ji+1
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-3-jl)
               END DO
            END DO
         CASE ( 'I' )                                     ! ice U-V point
               pt2d(2,ijpj) = psgn * pt2d(3,ijpj-1)
               DO ji = 3, npiglo
                  iju = npiglo - ji + 3
                  pt2d(ji,ijpj) = psgn * pt2d(iju,ijpj-1)
               END DO
         END SELECT

      CASE ( 'F','f' )                        ! *  North fold  F-point pivot

         SELECT CASE ( cdtype )
         CASE ( 'T' , 'W' ,'S' )                          ! T-, W-point
            DO ji = 1, npiglo
               ijt = npiglo-ji+1
               pt2d(ji,ijpj) = psgn * pt2d(ijt,ijpj-1)
            END DO
         CASE ( 'U' )                                     ! U-point
            DO ji = 1, npiglo-1
               iju = npiglo-ji
               pt2d(ji,ijpj) = psgn * pt2d(iju,ijpj-1)
            END DO
         CASE ( 'V' )                                     ! V-point
            DO ji = 1, npiglo
               ijt = npiglo-ji+1
               pt2d(ji,ijpj) = psgn * pt2d(ijt,ijpj-2)
            END DO
            DO ji = npiglo/2+1, npiglo
               ijt = npiglo-ji+1
               pt2d(ji,ijpjm1) = psgn * pt2d(ijt,ijpjm1)
            END DO
         CASE ( 'F' , 'G' )                               ! F-point
               DO ji = 1, npiglo-1
                  iju = npiglo-ji
                  pt2d(ji,ijpj) = psgn * pt2d(iju,ijpj-2)
               END DO
            DO ji = npiglo/2+1, npiglo-1
               iju = npiglo-ji
               pt2d(ji,ijpjm1) = psgn * pt2d(iju,ijpjm1)
            END DO
         CASE ( 'I' )                                  ! ice U-V point
            pt2d( 2 ,ijpj:ijpj) = 0.e0
               DO ji = 2 , npiglo-1
                  ijt = npiglo - ji + 2
                  pt2d(ji,ijpj)= 0.5 * ( pt2d(ji,ijpj-1) + psgn * pt2d(ijt,ijpj-1) )
               END DO
         END SELECT

      CASE DEFAULT                           ! *  closed : the code probably never go through
         PRINT * ,' ERROR : not a North condition '

      END SELECT

   END SUBROUTINE lbc_nfd_2d

END PROGRAM cdfnorth_unfold
