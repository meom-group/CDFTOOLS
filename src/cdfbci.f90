PROGRAM cdfbci
  !!======================================================================
  !!                     ***  PROGRAM  cdfbci  ***
  !!=====================================================================
  !!  ** Purpose : Compute the term of energetic transfert BCI
  !!               for the baroclinic instability
  !!
  !!  ** Method  : take  an input file which is the result of a preprocessing
  !!               tool cdfmoyuvwt.
  !!
  !! History : 2.1  : 02/2008  : A. Melet     : Original code
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
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

  INTEGER(KIND=4)                           :: ji, jj, jk
  INTEGER(KIND=4)                           :: ilev
  INTEGER(KIND=4)                           :: npiglo, npjglo, npk, npt
  INTEGER(KIND=4)                           :: narg, iargc
  INTEGER(KIND=4)                           :: ncout, ierr
  INTEGER(KIND=4), DIMENSION(5)             :: ipk, id_varout         ! 

  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2t, e1t 
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: anout, anovt, un, vn, tn
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: utn, vtn
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask, umask, vmask
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: bci
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rdtdx, rdtdy

  CHARACTER(LEN=256)                        :: cf_in
  CHARACTER(LEN=256)                        :: cf_out='bci.nc'

  TYPE (variable), DIMENSION(5)             :: stypvar         !: structure for attibutes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()   ! load cdf variable name

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfbci UVWT-file'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute elements for analysing the baroclinic instability' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       UVWT-file : input file is produced by cdfmoyuvwt, and the mean'
     PRINT *,'              must be computed on a long-enough period for the ' 
     PRINT *,'              statistics to be meaningful. Points are on T grid.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       Need ', TRIM(cn_fhgr) ,' file'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : 5 output variables'
     PRINT *,'             dTdx : zonal derivative of Tbar on T point (*1000)'
     PRINT *,'             dTdy : meridional derivative of Tbar on T point (*1000)'
     PRINT *,'             uT   : anomaly of u times anomaly of T on T point'
     PRINT *,'             vT   : anomaly of v times anomaly of T on T point'
     PRINT *,'             bci  : transfert of energy for the baroclinic instability (*1000)'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfmoyuvwt '
     STOP
  ENDIF

  CALL getarg(1, cf_in)
  IF (chkfile(cf_in) .OR. chkfile (cn_fhgr) ) STOP 99

  npiglo = getdim(cf_in, cn_x)
  npjglo = getdim(cf_in, cn_y)
  npk    = getdim(cf_in, cn_z)
  npt    = getdim(cf_in, cn_t)

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  ! define new variables for output ( must update att.txt)
  stypvar(1)%cname       = 'dTdx'
  stypvar(1)%clong_name  = 'zonal derivate of Tbar on T point (*1000)'
  stypvar(1)%cshort_name = 'dTdx'

  stypvar(2)%cname       = 'dTdy'
  stypvar(2)%clong_name  = 'meridional derivate of Tbar on T point (*1000)'
  stypvar(2)%cshort_name = 'dTdy'
 
  stypvar(3)%cname       = 'uT'
  stypvar(3)%clong_name  = 'anomaly of u times anomaly of T on T point'
  stypvar(3)%cshort_name = 'uT'
  
  stypvar(4)%cname       = 'vT'
  stypvar(4)%clong_name  = 'anomaly of v times anomaly of T on T point'
  stypvar(4)%cshort_name = 'vT'

  stypvar(5)%cname       = 'bci'
  stypvar(5)%clong_name  = 'transfert of energy for the baroclinic instability (*1000)'
  stypvar(5)%cshort_name = 'bci'

  stypvar%cunits            = '1000 (u"T" dT/dx + v"T" dT/dy)'
  stypvar%rmissing_value    = 0.
  stypvar%valid_min         = -1000.
  stypvar%valid_max         = 1000.
  stypvar%conline_operation = 'N/A'
  stypvar%caxis             = 'TYX'

  ipk(:) = npk  
  
  !test if lev exists
  IF ((npk==0) .AND. (ilev > 0) ) THEN
     PRINT *, 'Problem : npk = 0 and lev > 0 STOP'
     STOP 99
  END IF

  ! create output fileset
  ncout = create      (cf_out,   cf_in,   npiglo, npjglo, npk       )
  ierr  = createvar   (ncout ,   stypvar, 5,      ipk,    id_varout )
  ierr  = putheadervar(ncout,    cf_in,   npiglo, npjglo, npk       )

  ! Allocate the memory
  ALLOCATE ( e1t(npiglo,npjglo) , e2t(npiglo,npjglo) )
  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo)  )
  ALLOCATE ( tn(npiglo,npjglo) )
  ALLOCATE ( utn(npiglo,npjglo) , vtn(npiglo,npjglo) )
  ALLOCATE ( umask(npiglo,npjglo) , vmask(npiglo,npjglo) )
  ALLOCATE ( tmask(npiglo,npjglo) )
  ALLOCATE ( rdtdx(npiglo,npjglo)  , rdtdy(npiglo,npjglo)  )
  ALLOCATE ( anout(npiglo,npjglo)  , anovt(npiglo,npjglo)  )
  ALLOCATE ( bci(npiglo,npjglo), tim(npt) )

  e1t =  getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo)
  e2t =  getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo)

  tim  = getvar1d(cf_in, cn_vtimec, npt)
  ierr = putvar1d(ncout, tim, npt, 'T')
     
  DO jk=1, npk
     PRINT *,'            level ',jk
           
     rdtdx(:,:) = 0.0
     rdtdy(:,:) = 0.0
        
     anovt(:,:) = 0.0      
     anout(:,:) = 0.0

     un(:,:)  =  getvar(cf_in, 'ubar',  jk ,npiglo, npjglo, ktime=1)
     vn(:,:)  =  getvar(cf_in, 'vbar',  jk ,npiglo, npjglo, ktime=1)
     tn(:,:)  =  getvar(cf_in, 'tbar',  jk ,npiglo, npjglo, ktime=1)
     utn(:,:) =  getvar(cf_in, 'utbar', jk ,npiglo, npjglo, ktime=1)
     vtn(:,:) =  getvar(cf_in, 'vtbar', jk ,npiglo, npjglo, ktime=1)

     ! compute the mask
     DO jj = 2, npjglo
        DO ji = 2, npiglo
           umask(ji,jj)= un(ji,jj)*un(ji-1,jj) 
           vmask(ji,jj)= vn(ji,jj)*vn(ji,jj-1)
           tmask(ji,jj)= tn(ji,jj)
           IF (umask(ji,jj) /= 0.) umask(ji,jj)=1.
           IF (vmask(ji,jj) /= 0.) vmask(ji,jj)=1.
           IF (tmask(ji,jj) /= 0.) tmask(ji,jj)=1.
        END DO
     END DO

     DO jj = 2, npjglo  
        DO ji = 2, npiglo    ! vector opt.
           ! compute derivatives at T point
           rdtdx(ji,jj) = 1000/2. *( ( tn(ji,jj ) - tn(ji-1,jj) )      &
                &           * tmask(ji,jj)*tmask(ji-1,jj)              &
                &           / ( 0.5* ( e1t(ji,jj) + e1t(ji-1,jj) ))    &
                &           +( tn(ji+1,jj ) - tn(ji,jj) )              &
                &           * tmask(ji+1,jj)*tmask(ji,jj)              &
                &           / ( 0.5* ( e1t(ji+1,jj) + e1t(ji,jj) )))

           rdtdy(ji,jj) = 1000/2. *( ( tn(ji,jj) - tn(ji,jj-1) )       &
                &           * tmask(ji,jj)*tmask(ji,jj-1)              &
                &           / ( 0.5* ( e2t(ji,jj) + e2t(ji,jj-1) ))    &
                &           +( tn(ji,jj+1 ) - tn(ji,jj) )              &
                &           * tmask(ji,jj+1)*tmask(ji,jj)              &
                &           / ( 0.5* ( e2t(ji,jj+1) + e2t(ji,jj) )) )

           anout(ji,jj)    = ( utn(ji,jj)                                            &
                &                 -   1/2 * umask(ji,jj)*( un(ji,jj) + un(ji-1,jj) ) &
                &                   * tmask(ji,jj) * tn(ji,jj) )

           anovt(ji,jj)    = ( vtn(ji,jj)                                            &
                &                 -   1/2 * vmask(ji,jj)*( vn(ji,jj) + vn(ji,jj-1) ) &
                &                   * tmask(ji,jj) * tn(ji,jj) )

           ! compute bci term
           bci(ji,jj) = ( anout(ji,jj) * rdtdx(ji,jj) + anovt(ji,jj) * rdtdy(ji,jj) )

        END DO
     END DO
     ! 
     ierr = putvar(ncout, id_varout(1), rdtdx, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(2), rdtdy, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(3), anout, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(4), anovt, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(5), bci,   jk, npiglo, npjglo, ktime=1)
  END DO
  ierr = closeout(ncout)

END PROGRAM cdfbci

