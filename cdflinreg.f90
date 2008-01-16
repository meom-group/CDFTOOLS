PROGRAM cdflinreg
  !!-----------------------------------------------------------------------
  !!                 ***  PROGRAM cdflinreg  ***
  !!
  !!  **  Purpose: Compute linear regression coef from a bunch of input files.
  !!                of cdf files given as argument
  !!                Store the results on a 'similar' cdf file.
  !!  
  !!  **  Method: compute a and b such as yr = a . t + b 
  !!              yr is the estimation of the field value, t is the time (in days ).
  !!              a= ( moy(y.t) - moy(y).moy(t) ) / (moy(t2) - moy(t).moy(t) )
  !!              b= moy(y) - a . moy(t) 
  !!              R2 pearson value [0,1], giving the quality of the adjustment is also given
  !!              R2= ( a.a.moy(t.t) -2a.b.moy(t) +b.b -moy(y).moy(y) ) )/( moy(y.y) -moy(y).moy(y) )
  !!
  !! history :
  !!     Original code :   J.M. Molines (Jan 2008 ) from cdfmoy
  !!
  !!-----------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------

  USE cdfio 

  IMPLICIT NONE
  INTEGER, PARAMETER ::  jptmax=365                         !: maximum number of time frame
  INTEGER   :: jk,jt,jvar, jv , jtt,jkk                     !: dummy loop index
  INTEGER   :: ierr, ijvar                                  !: working integer
  INTEGER   :: narg, iargc                                  !: 
  INTEGER   :: npiglo,npjglo, npk ,nt                       !: size of the domain
  INTEGER   :: nvars                                        !: Number of variables in a file
  INTEGER   :: ntframe                                      !: Cumul of time frame
  INTEGER , DIMENSION(:), ALLOCATABLE :: id_var , &         !: arrays of var id's
       &                             ipk    , &             !: arrays of vertical level for each var
       &                             ipk2   , &             !: arrays of vertical level for each var
       &                             id_varout,& 
       &                             id_varout2
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: zy, zyy, zyt    !: Arrays for cumulated values
  REAL(KIND=8)                                :: zt, zt2    !: variables for cumulated time values
  REAL(KIND=8)                                :: total_time
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: v2d ,&      !: Array to read a layer of data
       &                                   rmean, rmean2,rmean3,  &
       &                                   areg, breg, rpear !: slope, origin ordinate, pearson coef
  REAL(KIND=4),DIMENSION(2)                   :: timean     !: trick : timean(1) hold moy(t) (days)
                                                            !:         timean(2) hold moy(t2) (days)**2
  REAL(KIND=4),DIMENSION(365)                 ::  tim
  REAL(KIND=4)                                :: spval = -99999.

  CHARACTER(LEN=80) :: cfile ,cfileout, cfileout2           !: file name
  CHARACTER(LEN=80) ::  cdep
  CHARACTER(LEN=80) ,DIMENSION(:), ALLOCATABLE:: cvarname   !: array of var name
  CHARACTER(LEN=80) ,DIMENSION(:), ALLOCATABLE:: cvarname2  !: array of var22 name for output
  
  TYPE (variable), DIMENSION(:), ALLOCATABLE :: typvar, typvar2

  INTEGER    :: ncout, ncout2
  INTEGER    :: istatus
  LOGICAL    :: lcaltmean

  !!

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdflinreg ''list_of_ioipsl_model_output_files'' '
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cfile)

  npiglo= getdim (cfile,'x')
  npjglo= getdim (cfile,'y')
  npk   = getdim (cfile,'depth',cdtrue=cdep, kstatus=istatus)

  IF (istatus /= 0 ) THEN
     npk   = getdim (cfile,'z',cdtrue=cdep,kstatus=istatus)
     IF (istatus /= 0 ) THEN
       npk   = getdim (cfile,'sigma',cdtrue=cdep,kstatus=istatus)
        IF ( istatus /= 0 ) THEN 
          PRINT *,' assume file with no depth'
          npk=0
        ENDIF
     ENDIF
  ENDIF
  

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( zy(npiglo,npjglo), zyt(npiglo,npjglo), zyy(npiglo,npjglo),v2d(npiglo,npjglo) )
  ALLOCATE( rmean(npiglo,npjglo), rmean2(npiglo,npjglo),  rmean3(npiglo,npjglo) )
  ALLOCATE( areg(npiglo,npjglo), breg(npiglo,npjglo) , rpear(npiglo,npjglo) )

  nvars = getnvar(cfile)
  PRINT *,' nvars =', nvars

  ALLOCATE (cvarname(nvars), cvarname2(3*nvars) )
  ALLOCATE (typvar(nvars), typvar2(3*nvars) )
  ALLOCATE (id_var(nvars),ipk(nvars),id_varout(nvars), id_varout2(3*nvars),ipk2(3*nvars)  )

  ! get list of variable names and collect attributes in typvar (optional)
  cvarname(:)=getvarname(cfile,nvars,typvar)

  DO jvar = 1, nvars
        ijvar=(jvar -1)*3 +1
        ! AREG
        cvarname2(ijvar)=TRIM(cvarname(jvar))//'_areg'
        typvar2(ijvar)%name =  TRIM(typvar(jvar)%name)//'_areg'      ! name
        typvar2(ijvar)%units = TRIM(typvar(jvar)%units)//'/year'      ! unit
        typvar2(ijvar)%missing_value = spval                         ! missing_value
        typvar2(ijvar)%valid_min = -100.                                    ! valid_min = zero
        typvar2(ijvar)%valid_max =  100.            ! valid_max *valid_max
        typvar2(ijvar)%scale_factor= 1.
        typvar2(ijvar)%add_offset= 0.
        typvar2(ijvar)%savelog10= 0.
        typvar2(ijvar)%long_name =TRIM(typvar(jvar)%long_name)//'_linear_slope'   ! 
        typvar2(ijvar)%short_name = TRIM(typvar(jvar)%short_name)//'_areg'     !
        typvar2(ijvar)%online_operation = TRIM(typvar(jvar)%online_operation) 
        typvar2(ijvar)%axis = TRIM(typvar(jvar)%axis) 
        ! BREG
        cvarname2(ijvar+1)=TRIM(cvarname(jvar))//'_breg'
        typvar2(ijvar+1)%name =  TRIM(typvar(jvar)%name)//'_breg'      ! name
        typvar2(ijvar+1)%units = TRIM(typvar(jvar)%units)              ! unit
        typvar2(ijvar+1)%missing_value = spval                         ! missing_value
        typvar2(ijvar+1)%valid_min = -100.                                    ! valid_min = zero
        typvar2(ijvar+1)%valid_max =  100.            ! valid_max *valid_max
        typvar2(ijvar+1)%scale_factor= 1.
        typvar2(ijvar+1)%add_offset= 0.
        typvar2(ijvar+1)%savelog10= 0.
        typvar2(ijvar+1)%long_name =TRIM(typvar(jvar)%long_name)//'_b'   !
        typvar2(ijvar+1)%short_name = TRIM(typvar(jvar)%short_name)//'_breg'     !
        typvar2(ijvar+1)%online_operation = TRIM(typvar(jvar)%online_operation)
        typvar2(ijvar+1)%axis = TRIM(typvar(jvar)%axis)
        ! R2 pearson
        cvarname2(ijvar+2)=TRIM(cvarname(jvar))//'_r2'
        typvar2(ijvar+2)%name =  TRIM(typvar(jvar)%name)//'_r2'      ! name
        typvar2(ijvar+2)%units = 'no unit'              ! unit
        typvar2(ijvar+2)%missing_value = spval                         ! missing_value
        typvar2(ijvar+2)%valid_min = 0.                                    ! valid_min = zero
        typvar2(ijvar+2)%valid_max =  1.            ! valid_max *valid_max
        typvar2(ijvar+2)%scale_factor= 1.
        typvar2(ijvar+2)%add_offset= 0.
        typvar2(ijvar+2)%savelog10= 0.
        typvar2(ijvar+2)%long_name =TRIM(typvar(jvar)%long_name)//'_r2_Pearson'   !
        typvar2(ijvar+2)%short_name = TRIM(typvar(jvar)%short_name)//'_r2'     !
        typvar2(ijvar+2)%online_operation = TRIM(typvar(jvar)%online_operation)
        typvar2(ijvar+2)%axis = TRIM(typvar(jvar)%axis)
  END DO

  id_var(:)  = (/(jv, jv=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cfile,nvars,cdep=cdep)
  DO jvar=1,nvars
    ipk2( (jvar-1)*3 +1 ) = ipk(jvar)
    ipk2( (jvar-1)*3 +2 ) = ipk(jvar)
    ipk2( (jvar-1)*3 +3 ) = ipk(jvar)
  ENDDO
  WHERE( ipk == 0 ) cvarname='none'
  WHERE( ipk2 == 0 ) cvarname2='none'
  typvar(:)%name=cvarname
  typvar2(:)%name=cvarname2

  ! create output fileset
  cfileout='cdflinreg.nc'
  cfileout2='linreg.nc'
  ! create output file taking the sizes in cfile

! ncout =create(cfileout, cfile,npiglo,npjglo,npk,cdep=cdep)
  ncout2=create(cfileout2,cfile,npiglo,npjglo,npk,cdep=cdep)

! ierr= createvar(ncout , typvar,  nvars, ipk, id_varout )
  ierr= createvar(ncout2, typvar2, 3*nvars, ipk2, id_varout2)

! ierr= putheadervar(ncout , cfile, npiglo, npjglo, npk,cdep=cdep)
  ierr= putheadervar(ncout2, cfile, npiglo, npjglo, npk,cdep=cdep)

  lcaltmean=.TRUE. ; zt=0.d0 ; zt2=0.d0
  DO jvar = 1,nvars
     ijvar=(jvar-1)*3 +1
     IF (cvarname(jvar) == 'nav_lon' .OR. &
          cvarname(jvar) == 'nav_lat' ) THEN
        ! skip these variable
     ELSE
        PRINT *,' Working with ', TRIM(cvarname(jvar)), ipk(jvar), jvar
        DO jk = 1, ipk(jvar)
!          PRINT *,'level ',jk
           zy(:,:) = 0.d0 ; zyt(:,:) = 0.d0 ; zyy(:,:) =0.d0 ; total_time = 0.;  ntframe=0
           DO jt = 1, narg
              CALL getarg (jt, cfile)
              nt = getdim (cfile,'time_counter')
              ntframe=ntframe+nt
              IF ( lcaltmean )  THEN
                tim(ntframe-nt+1:ntframe)=getvar1d(cfile,'time_counter',nt)/86400.d0/365.
!               tim(ntframe-nt+1:ntframe)=(/(ntframe-nt+jtt,jtt=1,nt)/)
              END IF
              DO jtt=1,nt
                jkk=jk
                ! If forcing fields is without depth dimension
                IF (npk==0) jkk=jtt 
                v2d(:,:)= getvar(cfile, cvarname(jvar), jkk ,npiglo, npjglo,ktime=jtt )
                zy(:,:)  = zy(:,:)  + v2d(:,:)
                zyy(:,:) = zyy(:,:) + v2d(:,:)*v2d(:,:)
                zyt(:,:) = zyt(:,:) + v2d(:,:)*tim(ntframe-nt+jtt)
              ENDDO
           END DO
           ! finish with level jk ; compute mean (assume spval is 0 )
           zt=sum(tim(1:ntframe))
           zt2=sum(tim(1:ntframe)*tim(1:ntframe) )
           rmean(:,:) = zy(:,:)/ntframe
           rmean2(:,:) = zyt(:,:)/ntframe
           rmean3(:,:) = zyy(:,:)/ntframe
           ! store variable on outputfile
!          ierr = putvar(ncout2,id_varout2(jvar) ,rmean, jk, npiglo, npjglo)
!          ierr = putvar(ncout2,id_varout2(jvar),rmean2, jk,npiglo, npjglo)
           IF (lcaltmean )  THEN
              timean(1)= zt/ntframe
              timean(2)= zt2/ntframe
!             ierr=putvar1d(ncout,timean,2,'T')
              ierr=putvar1d(ncout2,timean,1,'T')
           END IF
           !compute areg, breg, rpear
           WHERE (rmean /= 0 ) 
           areg(:,:)=( rmean2(:,:) - rmean(:,:) *timean(1) ) / ( timean(2) -timean(1)*timean(1) )
           breg(:,:)=rmean(:,:) - areg(:,:)*timean(1) 
  !!              R2= ( a.a.moy(t.t) -2a.b.moy(t) +b.b -moy(y).moy(y) ) )/( moy(y.y) -moy(y).moy(y) )
!          rpear(:,:) = (areg(:,:)*areg(:,:)*timean(2) -2*areg(:,:)*breg(:,:)*timean(1) -rmean(:,:)*rmean(:,:) ) / &
!            &          ( rmean2(:,:) -rmean(:,:)*rmean(:,:) )
           rpear(:,:) = areg(:,:)*areg(:,:)*( timean(2) -timean(1)*timean(1))/( rmean3(:,:) -rmean(:,:)*rmean(:,:) )
           WHERE (rpear < 0 ) rpear=0 ; WHERE (rpear > 1 ) rpear=1
           ELSEWHERE
            areg=spval ; breg=spval ; rpear=spval
           ENDWHERE
           
           ierr = putvar(ncout2,id_varout2(ijvar) ,REAL(areg), jk, npiglo, npjglo)
           ierr = putvar(ncout2,id_varout2(ijvar+1),REAL(breg), jk,npiglo, npjglo)
           ierr = putvar(ncout2,id_varout2(ijvar+2),REAL(rpear), jk,npiglo, npjglo)
           lcaltmean=.FALSE. ! tmean already computed
        END DO  ! loop to next level
     END IF
  END DO ! loop to next var in file

! istatus = closeout(ncout)
  istatus = closeout(ncout2)


END PROGRAM cdflinreg
