PROGRAM cdfnrj_bci
  !!======================================================================
  !!                     ***  PROGRAM  cdfnrj_bci  ***
  !!=====================================================================
  !!  ** Purpose : Compute the term of energetic transfert BCI
  !!               for the baroclinic instability
  !!
  !!  ** Method  : take  an input file which is the result of a preprocessing
  !!               tool cdfuvwt.
  !!
  !! History : 2.1  : 02/2008  : A. Melet     : Original code
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class energy_diagnostics
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jk, jt, jv       ! dummy loop index
  INTEGER(KIND=4)                           :: npiglo, npjglo, npk, npt ! domain size
  INTEGER(KIND=4)                           :: narg, iargc, ijarg       ! line parser
  INTEGER(KIND=4)                           :: ncout, ierr              ! netcdf i/o
  INTEGER(KIND=4), DIMENSION(5)             :: ipk, id_varout           ! idem

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2t, e1t                 ! horizontal metric
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: anout, anovt, un, vn, tn ! working variables
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: utn, vtn                 !  idem
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask, umask, vmask      ! mask variable
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: bci                      ! baroclinic energy transfert
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rdtdx, rdtdy             ! horizontal gradients

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim                     ! time variable

  CHARACTER(LEN=256)                        :: cf_in                    ! input filename
  CHARACTER(LEN=256)                        :: cf_out='bci.nc'          ! output file name
  CHARACTER(LEN=256)                        :: cldum                    ! working char variable

  TYPE (variable), DIMENSION(5)             :: stypvar                  ! structure for attibutes

  LOGICAL                                   :: lnc4 = .FALSE.           ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()   ! load cdf variable name

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfnrj_bci -f UVWT-file [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute elements (see variable list below) for analysing the baroclinic'
     PRINT *,'       instability.' 
     PRINT *,'       Note : this program was formerly named cdfbci.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f UVWT-file : input file is produced by cdfuvwt, and the mean'
     PRINT *,'              must be computed on a long-enough period for the statistics' 
     PRINT *,'              to be meaningful. Points are on T grid.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file] : specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ]   : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'              This option is effective only if cdftools are compiled with'
     PRINT *,'              a netcdf library supporting chunking and deflation.'
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
     PRINT *,'       cdfuvwt, cdfnrj_bti, cdfnrj_components, cdfnrj_transfert '
     STOP 
  ENDIF

  ijarg = 1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_in ) ; ijarg=ijarg+1
        ! option
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 99
     END SELECT
  ENDDO

  IF (chkfile(cf_in) .OR. chkfile (cn_fhgr) ) STOP 99

  npiglo = getdim(cf_in, cn_x)
  npjglo = getdim(cf_in, cn_y)
  npk    = getdim(cf_in, cn_z)
  npt    = getdim(cf_in, cn_t)

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  ! Allocate the memory
  ALLOCATE ( e1t(npiglo,npjglo) , e2t(npiglo,npjglo) )
  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo)  )
  ALLOCATE ( tn(npiglo,npjglo) )
  ALLOCATE ( utn(npiglo,npjglo) , vtn(npiglo,npjglo) )
  ALLOCATE ( umask(npiglo,npjglo) , vmask(npiglo,npjglo) )
  ALLOCATE ( tmask(npiglo,npjglo) )
  ALLOCATE ( rdtdx(npiglo,npjglo)  , rdtdy(npiglo,npjglo)  )
  ALLOCATE ( anout(npiglo,npjglo)  , anovt(npiglo,npjglo)  )
  ALLOCATE ( bci(npiglo,npjglo), dtim(npt) )

  CALL CreateOutput

  e1t =  getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo)
  e2t =  getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo)

  DO jt=1,npt
     DO jk=1, npk
        PRINT *,'            level ',jk

        rdtdx(:,:) = 0.0
        rdtdy(:,:) = 0.0

        anovt(:,:) = 0.0      
        anout(:,:) = 0.0

        un(:,:)  =  getvar(cf_in, 'ubar',  jk ,npiglo, npjglo, ktime=jt)
        vn(:,:)  =  getvar(cf_in, 'vbar',  jk ,npiglo, npjglo, ktime=jt)
        tn(:,:)  =  getvar(cf_in, 'tbar',  jk ,npiglo, npjglo, ktime=jt)
        utn(:,:) =  getvar(cf_in, 'utbar', jk ,npiglo, npjglo, ktime=jt)
        vtn(:,:) =  getvar(cf_in, 'vtbar', jk ,npiglo, npjglo, ktime=jt)

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
                   &                 -   0.5 * umask(ji,jj)*( un(ji,jj) + un(ji-1,jj) ) &
                   &                   * tmask(ji,jj) * tn(ji,jj) )

              anovt(ji,jj)    = ( vtn(ji,jj)                                            &
                   &                 -   0.5 * vmask(ji,jj)*( vn(ji,jj) + vn(ji,jj-1) ) &
                   &                   * tmask(ji,jj) * tn(ji,jj) )

              ! compute bci term
              bci(ji,jj) = ( anout(ji,jj) * rdtdx(ji,jj) + anovt(ji,jj) * rdtdy(ji,jj) )

           END DO
        END DO
        ! 
        ierr = putvar(ncout, id_varout(1), rdtdx, jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(2), rdtdy, jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(3), anout, jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(4), anovt, jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(5), bci,   jk, npiglo, npjglo, ktime=jt)
     END DO !level-loop
  ENDDO ! time-loop
  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    DO jv = 1,5
       stypvar(jv)%ichunk  = (/npiglo,MAX(1,npjglo/30),1,1 /)
    ENDDO

    ipk(:) = npk  

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


    ! create output fileset
    ncout = create      (cf_out,   cf_in,   npiglo, npjglo, npk       , ld_nc4=lnc4 )
    ierr  = createvar   (ncout ,   stypvar, 5,      ipk,    id_varout , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,    cf_in,   npiglo, npjglo, npk       )

    dtim = getvar1d(cf_in, cn_vtimec, npt)
    ierr = putvar1d(ncout, dtim, npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfnrj_bci

