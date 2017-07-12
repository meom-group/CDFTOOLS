PROGRAM cdfnrj_bti
  !!======================================================================
  !!                     ***  PROGRAM  cdfnrj_bti  ***
  !!=====================================================================
  !!  ** Purpose : Compute the term of energetic transfert BTI
  !!               for the barotropic instability for given gridU 
  !!               gridV gridU2 gridV2 files and variables
  !!
  !!  ** Method  : Take an input file which is preprocessed by
  !!               cdfuvwt. See also cdfbci
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

  INTEGER(KIND=4), PARAMETER                :: jp_varout  = 8
  INTEGER(KIND=4), PARAMETER                :: jp_dudx    = 1, jp_dvdx    = 2
  INTEGER(KIND=4), PARAMETER                :: jp_dudy    = 3, jp_dvdy    = 4
  INTEGER(KIND=4), PARAMETER                :: jp_anousqrt= 5, jp_anovsqrt= 6
  INTEGER(KIND=4), PARAMETER                :: jp_anouv   = 7, jp_bti     = 8
  INTEGER(KIND=4)                           :: ji, jj, jk, jt,jv      ! dummy loop index
  INTEGER(KIND=4)                           :: npiglo, npjglo         ! domain size
  INTEGER(KIND=4)                           :: npk, npt               ! vertical and time
  INTEGER(KIND=4)                           :: narg, iargc, ijarg     ! command line parser
  INTEGER(KIND=4)                           :: ncout, ierr            ! ncid of output file, error status
  INTEGER(KIND=4), DIMENSION(jp_varout)     :: ipk, id_varout         ! 

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2t, e1t, e1f, e2f     !
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: un, vn, u2n, v2n, uvn  !
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: fmask, umask, vmask    !
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: anousqrt, anovsqrt     !
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: anouv, bti             !
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: dudx, dudy, dvdx, dvdy !

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim

  CHARACTER(LEN=256)                        :: cf_out='bti.nc'        ! output file name
  CHARACTER(LEN=256)                        :: cf_uvwtfil             ! input file name
  CHARACTER(LEN=256)                        :: cldum                  ! working char variable

  TYPE (variable), DIMENSION(jp_varout)     :: stypvar                ! structure for attibutes

  LOGICAL                                   :: lchk                   ! flag for missing files
  LOGICAL                                   :: lnc4 = .FALSE.         ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  !!
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfnrj_bti -f UVWT-file [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute  the terms in the barotropic energy tranfert equation.'
     PRINT *,'       The transfert of energy for the barotropic instability is '
     PRINT *,'       bti= -[(u''bar)^2*dubar/dx ...'
     PRINT *,'             +(v''bar)^2*dvbar/dy ...'
     PRINT *,'             +(u''v''*(dubar/dy +dvbar/dx))]'
     PRINT *,'       Note : This program was formerly named cdfbti.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f UVWT-file : netcdf file produced by cdfuvwt'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file] : specify the output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ]  : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'             This option is effective only if cdftools are compiled with'
     PRINT *,'             a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fhgr) 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : '
     PRINT *,'              dudx : zonal derivate of ubar on T point'
     PRINT *,'              dvdx : zonal derivate of vbar on T point'
     PRINT *,'              dudy : meridional derivate of ubar on T point'
     PRINT *,'              dvdy : meridional derivate of vbar on T point'
     PRINT *,'              anousqrt : mean of (u-ubar)^2 on T point'
     PRINT *,'              anovsqrt : mean of (v-vbar)^2 on T point'
     PRINT *,'              anouv : mean of (u-ubar)*(v-vbar) on T point'
     PRINT *,'              bti  : transfert of energy for the barotropic instability.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      cdfuvwt, cdfnrj_bci, cdfnrj_components, cdfnrj_transfert'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_uvwtfil ) ; ijarg=ijarg+1
        ! option
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out     ) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 99
     END SELECT
  ENDDO

  lchk = chkfile (cn_fhgr )
  lchk = lchk .OR. chkfile (cf_uvwtfil )
  IF ( lchk ) STOP 99 ! missing file

  npiglo  = getdim(cf_uvwtfil,cn_x)
  npjglo  = getdim(cf_uvwtfil,cn_y)
  npk     = getdim(cf_uvwtfil,cn_z)
  npt     = getdim(cf_uvwtfil,cn_t)

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  ! Allocate the memory
  ALLOCATE ( e1t(npiglo,npjglo) , e1f(npiglo,npjglo) )
  ALLOCATE ( e2t(npiglo,npjglo) , e2f(npiglo,npjglo) )
  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo)  )
  ALLOCATE ( fmask(npiglo,npjglo) )
  ALLOCATE ( umask(npiglo,npjglo) , vmask(npiglo,npjglo) )
  ALLOCATE ( dudx(npiglo,npjglo)  , dudy(npiglo,npjglo)  )
  ALLOCATE ( dvdx(npiglo,npjglo)  , dvdy(npiglo,npjglo)  )
  ALLOCATE ( u2n(npiglo,npjglo)  , v2n(npiglo,npjglo)  )
  ALLOCATE ( uvn(npiglo,npjglo) )
  ALLOCATE ( anousqrt(npiglo,npjglo) , anovsqrt(npiglo,npjglo)  )
  ALLOCATE ( anouv(npiglo,npjglo), bti(npiglo,npjglo) )
  ALLOCATE ( dtim(npt) )

  CALL CreateOutput

  e1t = getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo)
  e1f = getvar(cn_fhgr, cn_ve1f, 1, npiglo, npjglo)
  e2t = getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo)
  e2f = getvar(cn_fhgr, cn_ve2f, 1, npiglo, npjglo)


  DO jt = 1, npt
     DO jk=1, npk
        PRINT *,'            level ',jk
        dudx(:,:) = 0.d0
        dvdx(:,:) = 0.d0
        dudy(:,:) = 0.d0
        dvdy(:,:) = 0.d0

        anousqrt(:,:) = 0.d0
        anovsqrt(:,:) = 0.d0      
        anouv(:,:)    = 0.d0

        un(:,:)  =  getvar(cf_uvwtfil, 'ubar',  jk ,npiglo,npjglo, ktime=jt)
        vn(:,:)  =  getvar(cf_uvwtfil, 'vbar',  jk ,npiglo,npjglo, ktime=jt)
        u2n(:,:) =  getvar(cf_uvwtfil, 'u2bar', jk ,npiglo,npjglo, ktime=jt)
        v2n(:,:) =  getvar(cf_uvwtfil, 'v2bar', jk ,npiglo,npjglo, ktime=jt)
        uvn(:,:) =  getvar(cf_uvwtfil, 'uvbar', jk ,npiglo,npjglo, ktime=jt)

        ! compute the masks
        umask(:,:) = 0. ; vmask(:,:) = 0. ; fmask(:,:) = 0.
        DO jj = 2, npjglo
           DO ji = 2, npiglo
              umask(ji,jj)= un(ji,jj)*un(ji-1,jj  ) 
              vmask(ji,jj)= vn(ji,jj)*vn(ji  ,jj-1)
           ENDDO
        ENDDO

        WHERE ( umask /= 0. ) umask = 1.
        WHERE ( vmask /= 0. ) vmask = 1.

        DO jj = 1, npjglo-1
           DO ji = 1, npiglo-1
              fmask(ji,jj)= un(ji,jj)*un(ji,jj+1) * vn(ji,jj)*vn(ji+1,jj)
           ENDDO
        ENDDO

        WHERE ( fmask /= 0. ) fmask = 1.

        DO jj = 2, npjglo  
           DO ji = 2, npiglo    ! vector opt.
              ! compute derivates at T points
              dudx(ji,jj) = 100000. * ( un(ji,jj ) - un(ji-1,jj) )   &
                   &               * umask(ji,jj) / e1t(ji,jj) 

              dvdy(ji,jj) = 100000. * ( vn(ji,jj ) - vn(ji,jj-1) )   &
                   &               * vmask(ji,jj) / e2t(ji,jj)             

              dudy(ji,jj) = 100000./4. *( ( un(ji,jj+1 ) - un(ji,jj) )   &
                   &           * fmask(ji,jj) / e2f(ji,jj)               &
                   &       + (un(ji,jj ) - un(ji,jj-1) )                 &
                   &           * fmask(ji,jj-1) / e2f(ji,jj-1)           &
                   &       + (un(ji-1,jj+1 ) - un(ji-1,jj) )             &
                   &           * fmask(ji-1,jj) / e2f(ji-1,jj)           &         
                   &       + (un(ji-1,jj ) - un(ji-1,jj-1) )             &    
                   &           * fmask(ji-1,jj-1) / e2f(ji-1,jj-1) )            

              dvdx(ji,jj) = 100000./4. *( ( vn(ji,jj ) - vn(ji-1,jj) )   &
                   &           * fmask(ji-1,jj) / e1f(ji-1,jj)           &
                   &       + (vn(ji+1,jj ) - vn(ji,jj) )                 &
                   &           * fmask(ji,jj) / e1f(ji,jj)               &
                   &       + (vn(ji-1,jj-1 ) - vn(ji,jj-1) )             &
                   &           * fmask(ji-1,jj-1) / e1f(ji-1,jj-1)       &
                   &       + (vn(ji+1,jj-1 ) - vn(ji,jj-1) )             &
                   &           * fmask(ji,jj-1) / e1f(ji,jj-1) )         

              ! Compute Reynolds terms
              anousqrt(ji,jj) = 1000./2. * umask(ji,jj)*( ( u2n(ji,jj) - un(ji,jj)*un(ji,jj) ) &
                   &                     + ( u2n(ji-1,jj) - un(ji-1,jj)*un(ji-1,jj) ) )       

              anovsqrt(ji,jj) = 1000./2. * vmask(ji,jj)*( ( v2n(ji,jj) - vn(ji,jj)*vn(ji,jj) ) &
                   &                     + ( v2n(ji,jj-1) - vn(ji,jj)*vn(ji,jj-1) ) ) 

              anouv(ji,jj)    = 1000. * ( uvn(ji,jj) &
                   &                 -   0.5 * umask(ji,jj)*( un(ji,jj) + un(ji-1,jj) ) &
                   &                   * 0.5 * vmask(ji,jj)*( vn(ji,jj) + vn(ji,jj-1) ) )

              ! Compute bti
              bti(ji,jj) = -1. * ( anousqrt(ji,jj) * dudx(ji,jj)                &
                   &           + anovsqrt(ji,jj) * dvdy(ji,jj)                  &
                   &           + anouv(ji,jj) * ( dvdx(ji,jj) + dudy(ji,jj) ))

           END DO
        END DO
        ! 
        ierr = putvar(ncout, id_varout(jp_dudx),     dudx,     jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(jp_dvdx),     dvdx,     jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(jp_dudy),     dudy,     jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(jp_dvdy),     dvdy,     jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(jp_anousqrt), anousqrt, jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(jp_anovsqrt), anovsqrt, jk, npiglo, npjglo, ktime=jt) 
        ierr = putvar(ncout, id_varout(jp_anouv),    anouv,    jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(jp_bti),      bti,      jk, npiglo, npjglo, ktime=jt)
     END DO
  END DO  ! time loop

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
    ! define new variables for output ( must update att.txt)
    DO jv = 1, jp_varout
       stypvar(jv)%ichunk        = (/npiglo,MAX(1,npjglo/30),1,1 /)
    ENDDO
    ipk(:) = npk  

    stypvar(jp_dudx)%cname       = 'dudx'
    stypvar(jp_dudx)%clong_name  = 'zonal derivate of u on T point'
    stypvar(jp_dudx)%cshort_name = 'dudx'

    stypvar(jp_dvdx)%cname       = 'dvdx'
    stypvar(jp_dvdx)%clong_name  = 'zonal derivate of v on T point'
    stypvar(jp_dvdx)%cshort_name = 'dvdx'

    stypvar(jp_dudy)%cname       = 'dudy'
    stypvar(jp_dudy)%clong_name  = 'meridional derivate of u on T point'
    stypvar(jp_dudy)%cshort_name = 'dudy'

    stypvar(jp_dvdy)%cname       = 'dvdy'
    stypvar(jp_dvdy)%clong_name  = 'meridional derivate of v on T point'
    stypvar(jp_dvdy)%cshort_name = 'dvdy'

    stypvar(jp_anousqrt)%cname       = 'anousqrt'
    stypvar(jp_anousqrt)%clong_name  = 'temporal mean of the square of the zonal speed anomaly'
    stypvar(jp_anousqrt)%cshort_name = 'anousqrt'

    stypvar(jp_anovsqrt)%cname       = 'anovsqrt'
    stypvar(jp_anovsqrt)%clong_name  = 'temporal mean of the square of the meridional speed anomaly'
    stypvar(jp_anovsqrt)%cshort_name = 'anovsqrt'

    stypvar(jp_anouv)%cname       = 'anouv'
    stypvar(jp_anouv)%clong_name  = 'temporal mean of the Reynolds term'
    stypvar(jp_anouv)%cshort_name = 'anouanov'

    stypvar(jp_bti)%cname         = 'bti'
    stypvar(jp_bti)%clong_name    = 'transfert of energy for the barotropic instability'
    stypvar(jp_bti)%cshort_name   = 'bti'

    stypvar%cunits            = '100000 s-1'
    stypvar%rmissing_value    = 0.
    stypvar%valid_min         = -1000.
    stypvar%valid_max         = 1000.
    stypvar%conline_operation = 'N/A'
    stypvar%caxis             = 'TYX'

    ! create output fileset
    ncout = create      (cf_out, cf_uvwtfil, npiglo,    npjglo, npk      , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar,    jp_varout, ipk,    id_varout, ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_uvwtfil, npiglo,    npjglo, npk       )

    dtim = getvar1d(cf_uvwtfil, cn_vtimec, npt      )
    ierr = putvar1d(ncout, dtim,           npt, 'T' )

  END SUBROUTINE CreateOutput

END PROGRAM cdfnrj_bti

