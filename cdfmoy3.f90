PROGRAM cdfmoy3
  !!-----------------------------------------------------------------------
  !!                 ***  PROGRAM cdfmoy3  ***
  !!
  !!  **  Purpose: Compute mean values of ssh ssh^2 and ssh^3
  !!                Store the results on a 'similar' cdf file.
  !!  
  !!  **  Method: Try to avoid 3 d arrays 
  !!
  !! history :
  !!     Original code :   J.M. Molines (Nov 2004 ) for ORCA025
  !!                       J.M. Molines (Apr 2005 ) put all NCF stuff in module
  !!                              now valid for grid T U V W icemod
  !!     Modified      :   P. Mathiot (June 2007) update for forcing fields
  !!-----------------------------------------------------------------------
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
  INTEGER   :: ntframe                                      !: Cumul of time frame
  INTEGER , DIMENSION(:), ALLOCATABLE :: id_var , &         !: arrays of var id's
       &                             ipk    , &         !: arrays of vertical level for each var
       &                             id_varout,& 
       &                             id_varout2,&
       &                             id_varout3
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: tab, tab2,tab3  !: Arrays for cumulated values
  REAL(KIND=8)                                :: total_time
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: v2d ,&       !: Array to read a layer of data
       &                                   rmean, rmean2, rmean3
  REAL(KIND=4),DIMENSION(1)                   :: timean
  REAL(KIND=4),DIMENSION(365)                   ::  tim

  CHARACTER(LEN=80) :: cfile ,cfileout, cfileout2, cfileout3  !: file name
  CHARACTER(LEN=80) ::  cdep
  CHARACTER(LEN=80) ,DIMENSION(:), ALLOCATABLE:: cvarname   !: array of var name
  CHARACTER(LEN=80) ,DIMENSION(:), ALLOCATABLE:: cvarname2   !: array of var22 name for output
  CHARACTER(LEN=80) ,DIMENSION(:), ALLOCATABLE:: cvarname3   !: array of var3 name for output
  
  TYPE (variable), DIMENSION(:), ALLOCATABLE :: typvar, typvar2,typvar3

  INTEGER    :: ncout, ncout2, ncout3
  INTEGER    :: istatus
  LOGICAL    :: lcaltmean

  !!

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmoy3 ''list_of_ioipsl_model_output_files'' '
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
  npk=1

  ALLOCATE( tab(npiglo,npjglo), tab2(npiglo,npjglo), v2d(npiglo,npjglo) )
  ALLOCATE( rmean(npiglo,npjglo), rmean2(npiglo,npjglo) )

  nvars = getnvar(cfile)
  PRINT *,' nvars =', nvars

  ALLOCATE (cvarname(nvars), cvarname2(nvars) , cvarname3(nvars) )
  ALLOCATE (typvar(nvars), typvar2(nvars) , typvar3(nvars) )
  ALLOCATE (id_var(nvars),ipk(nvars),id_varout(nvars), id_varout2(nvars) ,id_varout3(nvars) )

  ! get list of variable names and collect attributes in typvar (optional)
  cvarname(:)=getvarname(cfile,nvars,typvar)

  DO jvar = 1, nvars
     ! variables that will not be computed or stored are named 'none'
     IF (cvarname(jvar)  /= 'sossheig' ) THEN
          cvarname2(jvar) ='none'
     ELSE
        cvarname2(jvar)=TRIM(cvarname(jvar))//'_sqd'
        typvar2(jvar)%name =  TRIM(typvar(jvar)%name)//'_sqd'           ! name
        typvar2(jvar)%units = '('//TRIM(typvar(jvar)%units)//')^2'      ! unit
        typvar2(jvar)%missing_value = typvar(jvar)%missing_value        ! missing_value
        typvar2(jvar)%valid_min = 0.                                    ! valid_min = zero
        typvar2(jvar)%valid_max =  typvar(jvar)%valid_max**2            ! valid_max *valid_max
        typvar2(jvar)%scale_factor= 1.
        typvar2(jvar)%add_offset= 0.
        typvar2(jvar)%savelog10= 0.
        typvar2(jvar)%long_name =TRIM(typvar(jvar)%long_name)//'_Squared'   ! 
        typvar2(jvar)%short_name = TRIM(typvar(jvar)%short_name)//'_sqd'     !
        typvar2(jvar)%online_operation = TRIM(typvar(jvar)%online_operation) 
        typvar2(jvar)%axis = TRIM(typvar(jvar)%axis) 

        cvarname3(jvar)=TRIM(cvarname(jvar))//'_cub'
        typvar3(jvar)%name =  TRIM(typvar(jvar)%name)//'_cub'           ! name
        typvar3(jvar)%units = '('//TRIM(typvar(jvar)%units)//')^3'      ! unit
        typvar3(jvar)%missing_value = typvar(jvar)%missing_value        ! missing_value
        typvar3(jvar)%valid_min = 0.                                    ! valid_min = zero
        typvar3(jvar)%valid_max =  typvar(jvar)%valid_max**3            ! valid_max *valid_max
        typvar3(jvar)%scale_factor= 1.
        typvar3(jvar)%add_offset= 0.
        typvar3(jvar)%savelog10= 0.
        typvar3(jvar)%long_name =TRIM(typvar(jvar)%long_name)//'_Cubic'   !
        typvar3(jvar)%short_name = TRIM(typvar(jvar)%short_name)//'_cub'     !
        typvar3(jvar)%online_operation = TRIM(typvar(jvar)%online_operation)
        typvar3(jvar)%axis = TRIM(typvar(jvar)%axis)


     END IF
  END DO

  id_var(:)  = (/(jv, jv=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cfile,nvars,cdep=cdep)
  WHERE( ipk == 0 .OR. ipk > 1 ) cvarname='none'
  typvar(:)%name=cvarname
  typvar2(:)%name=cvarname2
  typvar3(:)%name=cvarname3

  ! create output fileset
  cfileout='cdfmoy.nc'
  cfileout2='cdfmoy2.nc'
  cfileout3='cdfmoy3.nc'
  ! create output file taking the sizes in cfile

  ncout =create(cfileout, cfile,npiglo,npjglo,npk,cdep=cdep)
  ncout2=create(cfileout2,cfile,npiglo,npjglo,npk,cdep=cdep)
  ncout3=create(cfileout3,cfile,npiglo,npjglo,npk,cdep=cdep)

  ierr= createvar(ncout , typvar,  nvars, ipk, id_varout )
  ierr= createvar(ncout2, typvar2, nvars, ipk, id_varout2)
  ierr= createvar(ncout3, typvar3, nvars, ipk, id_varout3)

  ierr= putheadervar(ncout , cfile, npiglo, npjglo, npk,cdep=cdep)
  ierr= putheadervar(ncout2, cfile, npiglo, npjglo, npk,cdep=cdep)
  ierr= putheadervar(ncout3, cfile, npiglo, npjglo, npk,cdep=cdep)

  lcaltmean=.TRUE.
  DO jvar = 1,nvars
     IF (cvarname(jvar) == 'nav_lon' .OR. &
          cvarname(jvar) == 'nav_lat' .OR. cvarname2(jvar) == 'none' ) THEN
        ! skip these variable
     ELSE
        PRINT *,' Working with ', TRIM(cvarname(jvar)), ipk(jvar)
        DO jk = 1, ipk(jvar)
           PRINT *,'level ',jk
           tab(:,:) = 0.d0 ; tab2(:,:) = 0.d0 ; tab3(:,:) = 0.d0 
           total_time = 0.;  ntframe=0
           DO jt = 1, narg
              CALL getarg (jt, cfile)
              nt = getdim (cfile,'time_counter')
              IF ( lcaltmean )  THEN
                 tim=getvar1d(cfile,'time_counter',nt)
                 total_time = total_time + SUM(tim(1:nt) )
              END IF
              DO jtt=1,nt
                ntframe=ntframe+1
                jkk=jk
                ! If forcing fields is without depth dimension
                IF (npk==0) jkk=jtt 
                v2d(:,:)= getvar(cfile, cvarname(jvar), jkk ,npiglo, npjglo,ktime=jtt )
                tab(:,:) = tab(:,:) + v2d(:,:)
                IF (cvarname2(jvar) /= 'none' ) tab2(:,:) = tab2(:,:) + v2d(:,:)*v2d(:,:)
                IF (cvarname3(jvar) /= 'none' ) tab3(:,:) = tab3(:,:) + v2d(:,:)*v2d(:,:)*v2d(:,:)
              ENDDO
           END DO
           ! finish with level jk ; compute mean (assume spval is 0 )
           rmean(:,:) = tab(:,:)/ntframe
           IF (cvarname2(jvar) /= 'none' ) rmean2(:,:) = tab2(:,:)/ntframe
           IF (cvarname3(jvar) /= 'none' ) rmean3(:,:) = tab3(:,:)/ntframe
           ! store variable on outputfile
           ierr = putvar(ncout, id_varout(jvar) ,rmean, jk, npiglo, npjglo)
           IF (cvarname2(jvar) /= 'none' ) ierr = putvar(ncout2,id_varout2(jvar),rmean2, jk,npiglo, npjglo)
           IF (cvarname3(jvar) /= 'none' ) ierr = putvar(ncout3,id_varout3(jvar),rmean3, jk,npiglo, npjglo)
           IF (lcaltmean )  THEN
              timean(1)= total_time/ntframe
              ierr=putvar1d(ncout,timean,1,'T')
              ierr=putvar1d(ncout2,timean,1,'T')
              ierr=putvar1d(ncout3,timean,1,'T')
           END IF
           lcaltmean=.FALSE. ! tmean already computed
        END DO  ! loop to next level
     END IF
  END DO ! loop to next var in file

  istatus = closeout(ncout)
  istatus = closeout(ncout2)
  istatus = closeout(ncout3)


END PROGRAM cdfmoy3
