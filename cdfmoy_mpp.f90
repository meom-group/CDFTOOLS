PROGRAM cdfmoy_mpp
  !!-----------------------------------------------------------------------
  !!                 ***  PROGRAM cdfmoy_mpp  ***
  !!
  !!  **  Purpose: Compute mean values for all the variables in a bunch
  !!                of cdf files given as argument
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
  USE mpi

  IMPLICIT NONE
  ! MPI stuff
  INTEGER :: jt,jjp, ji, ii                                    ! loop counters
  INTEGER :: ierror, iproc, nproc, narea
  INTEGER, DIMENSION(:), ALLOCATABLE :: nptag
  INTEGER :: irest, ntag, nused_proc, ndimtag , ntask
  CHARACTER(LEN=256) :: cdum
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: ctags         ! tag list
  CHARACTER(LEN=256), DIMENSION(:,:) , ALLOCATABLE:: cptag       ! processor tag list
  LOGICAL :: lwp
  !
  !
  INTEGER   :: jk,jvar, jv , jtt,jkk                     !: dummy loop index
  INTEGER   :: ierr                                         !: working integer
  INTEGER   :: narg, iargc                                  !: 
  INTEGER   :: npiglo,npjglo, npk ,nt                       !: size of the domain
  INTEGER   :: nvars                                        !: Number of variables in a file
  INTEGER   :: ntframe, ntframe_tot                         !: Cumul of time frame
  INTEGER , DIMENSION(:), ALLOCATABLE :: id_var , &         !: arrays of var id's
       &                             ipk    , &         !: arrays of vertical level for each var
       &                             id_varout,& 
       &                             id_varout2
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: tab, tab2  !: Arrays for cumulated values
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: tabtot, tab2tot  !: Arrays for cumulated values
  REAL(KIND=8)                                :: total_time, total_timetot
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: v2d ,&       !: Array to read a layer of data
       &                                   rmean, rmean2
  REAL(KIND=4),DIMENSION(1)                   :: timean
  REAL(KIND=4),DIMENSION(365)                   ::  tim

  CHARACTER(LEN=256) :: cfile ,cfileout, cfileout2           !: file name
  CHARACTER(LEN=256) ::  cdep
  CHARACTER(LEN=256) ,DIMENSION(:), ALLOCATABLE:: cvarname   !: array of var name
  CHARACTER(LEN=256) ,DIMENSION(:), ALLOCATABLE:: cvarname2   !: array of var22 name for output

  TYPE (variable), DIMENSION(:), ALLOCATABLE :: typvar, typvar2

  INTEGER    :: ncout, ncout2
  INTEGER    :: istatus
  LOGICAL    :: lcaltmean

  !! * Initialization

  ! Initialize MPI
  CALL mpi_init(ierror)
  CALL mpi_comm_rank(mpi_comm_world,iproc,ierror)
  CALL mpi_comm_size(mpi_comm_world,nproc,ierror)
  narea = iproc + 1
  lwp=( narea == 1 )
  !!

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     IF ( lwp ) THEN
        PRINT *,' Usage : cdfmoy_mpp ''list_of_ioipsl_model_output_files'' '
     ENDIF
     CALL mpi_finalize(ierror)
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


  IF ( lwp ) THEN
     PRINT *, 'npiglo=', npiglo
     PRINT *, 'npjglo=', npjglo
     PRINT *, 'npk   =', npk
  ENDIF

     ALLOCATE( tab(npiglo,npjglo), tab2(npiglo,npjglo), v2d(npiglo,npjglo) )
     ALLOCATE( tabtot(npiglo,npjglo), tab2tot(npiglo,npjglo))
     ALLOCATE( rmean(npiglo,npjglo), rmean2(npiglo,npjglo) )

     nvars = getnvar(cfile)
     PRINT *,' nvars =', nvars

     ALLOCATE (cvarname(nvars), cvarname2(nvars) )
     ALLOCATE (typvar(nvars), typvar2(nvars) )
     ALLOCATE (id_var(nvars),ipk(nvars),id_varout(nvars), id_varout2(nvars)  )

     ! get list of variable names and collect attributes in typvar (optional)
     cvarname(:)=getvarname(cfile,nvars,typvar)

     DO jvar = 1, nvars
        ! variables that will not be computed or stored are named 'none'
        IF (cvarname(jvar)  /= 'vozocrtx' .AND. &
             cvarname(jvar) /= 'vomecrty' .AND. &
             cvarname(jvar) /= 'vovecrtz' .AND. &
             cvarname(jvar) /= 'sossheig' ) THEN
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

        END IF
     END DO

     id_var(:)  = (/(jv, jv=1,nvars)/)
     ! ipk gives the number of level or 0 if not a T[Z]YX  variable
     ipk(:)     = getipk (cfile,nvars,cdep=cdep)
     WHERE( ipk == 0 ) cvarname='none'
     typvar(:)%name=cvarname
     typvar2(:)%name=cvarname2

     ! create output fileset
     cfileout='cdfmoy.nc'
     cfileout2='cdfmoy2.nc'
     ! create output file taking the sizes in cfile

     IF ( lwp) THEN
        ncout =create(cfileout, cfile,npiglo,npjglo,npk,cdep=cdep)
        ncout2=create(cfileout2,cfile,npiglo,npjglo,npk,cdep=cdep)

        ierr= createvar(ncout , typvar,  nvars, ipk, id_varout )
        ierr= createvar(ncout2, typvar2, nvars, ipk, id_varout2)

        ierr= putheadervar(ncout , cfile, npiglo, npjglo, npk,cdep=cdep)
        ierr= putheadervar(ncout2, cfile, npiglo, npjglo, npk,cdep=cdep)
     END IF

     ! Allocate space for taglist
     ALLOCATE ( ctags(narg) )

     DO jt=1,narg
        CALL getarg(jt,ctags(jt))
     END DO

     !! * Dispatch the tags  among the processors

     ! Max number of tags per processors
     irest = MOD(narg, nproc)
     ntag   = narg / nproc
     nused_proc = nproc

     IF ( ntag == 0 ) THEN   ! when there are more proc than tags
        ntag = 1             ! each working proc takes 1 tag
        irest = 0            ! no tags left
        nused_proc= narg     ! number of used proc is less than nproc
     END IF

     ! maximum possible tags per procs
     ndimtag = ntag  ;  IF ( irest /= 0 ) ndimtag = ndimtag + 1  ! irest task will be reparted

     ! Allocate space
     ALLOCATE (cptag(nproc,ndimtag) )  ! list of tags for each proc
     ALLOCATE (nptag(nproc) )          ! number of tags per proc

     nptag(:) = 0
     nptag(1:nused_proc) = ntag

     ! reparts the remaining tags on the first irest proc
     DO jjp= 1, irest
        nptag(jjp) = nptag(jjp) + 1
     END DO

     ! build processor tag list
     ii= 0
     DO jjp = 1, nused_proc
        DO jt = 1, nptag(jjp)
           ii = ii + 1
           cptag(jjp,jt) = ctags(ii)
        END DO
     END DO
!!!
     !! * Dispatch the work ..
     ntask = nptag(narea)


     lcaltmean=.TRUE.
     DO jvar = 1,nvars
        IF (cvarname(jvar) == 'nav_lon' .OR. &
             cvarname(jvar) == 'nav_lat' ) THEN
           ! skip these variable
        ELSE
           IF (lwp) PRINT *,' Working with ', TRIM(cvarname(jvar)), ipk(jvar)
           DO jk = 1, ipk(jvar)
              IF (lwp)    PRINT *,'level ',jk
              tab(:,:) = 0.d0 ; tab2(:,:) = 0.d0 ; total_time = 0.;  ntframe=0
              DO jt = 1, ntask
                 cfile=cptag(narea,jt)
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
                 ENDDO
              END DO
              !          ! finish with level jk ; compute mean (assume spval is 0 )
              !          rmean(:,:) = tab(:,:)/ntframe
              !          IF (cvarname2(jvar) /= 'none' ) rmean2(:,:) = tab2(:,:)/ntframe
              !          ! store variable on outputfile
              ! collect total number of frames 
              CALL MPI_REDUCE(ntframe, ntframe_tot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
              CALL MPI_REDUCE(total_time, total_timetot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
              ! Collect sum of tab
              CALL MPI_REDUCE(tab,tabtot,npiglo*npjglo,MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,ierr)
              IF (cvarname2(jvar)  /= 'none' ) THEN
                 CALL MPI_REDUCE(tab2,tab2tot,npiglo*npjglo,MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,ierr)
              ENDIF
              IF (lwp) THEN
                 rmean(:,:)=tabtot/ntframe_tot
                 IF (cvarname2(jvar)  /= 'none' ) rmean2(:,:)=tab2tot/ntframe_tot
                 ierr = putvar(ncout, id_varout(jvar) ,rmean, jk, npiglo, npjglo)
                 IF (cvarname2(jvar) /= 'none' ) ierr = putvar(ncout2,id_varout2(jvar),rmean2, jk,npiglo, npjglo)
                 IF (lcaltmean )  THEN
                    timean(1)= total_timetot/ntframe_tot
                    ierr=putvar1d(ncout,timean,1,'T')
                    ierr=putvar1d(ncout2,timean,1,'T')
                 END IF
              ENDIF
              lcaltmean=.FALSE. ! tmean already computed
           END DO  ! loop to next level
        END IF
     END DO ! loop to next var in file

     IF (lwp) THEN
        istatus = closeout(ncout)
        istatus = closeout(ncout2)
     ENDIF
     CALL mpi_finalize(ierror)

   END PROGRAM cdfmoy_mpp
