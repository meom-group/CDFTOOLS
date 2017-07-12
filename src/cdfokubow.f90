PROGRAM cdfokubow
  !!======================================================================
  !!                     ***  PROGRAM  cdfokubow  ***
  !!=====================================================================
  !!  ** Purpose :  Compute the okubow weiss parameter on F-points for 
  !!                given gridU gridV files and variables
  !!
  !!  ** Method  :  ?
  !!
  !! History :  3.0  : 08/2012  : N. Djath     : original code from cdfcurl
  !!         :  4.0  : 03/2017  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !! * Modules used
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class energy_diagnostics
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jk, jt                   ! dummy loop index
  INTEGER(KIND=4)                           :: ik1, ik2, iko                    ! limits of jk loop, k-index for output
  INTEGER(KIND=4)                           :: ilev=-1                          ! level to be processed
  INTEGER(KIND=4)                           :: npiglo, npjglo                   ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt                         ! size of the domain
  INTEGER(KIND=4)                           :: narg, iargc, ijarg               ! browse command line
  INTEGER(KIND=4)                           :: ncout, ierr                      ! browse command line
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout                   ! output variable properties

  REAL(KIND=4)                              :: zcoef
  REAL(KIND=4)                              :: zstrnsym2                        ! square of the Nsymetrical Strain T-pt.
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2v, e1u, e1f, e2f, e1t, e2t     ! horizontal metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: un, vn                           ! velocity field
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: okubow, fmask, tmask             ! obuko-weiss and mask
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rotn, strsym, strnsy             ! curl, symetric and non symetric strain
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: gdepu                            ! depth for output

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim                             ! time counter

  CHARACTER(LEN=256)                        :: cf_ufil, cf_vfil                 ! file names
  CHARACTER(LEN=256)                        :: cf_out = 'okubow.nc'             ! output file name
  CHARACTER(LEN=256)                        :: cv_u, cv_v                       ! variable names
  CHARACTER(LEN=256)                        :: cldum                            ! dummy string

  TYPE (variable), DIMENSION(1)             :: stypvar                          ! structure for attibutes

  LOGICAL                                   :: lchk     = .FALSE.               ! flag for missing files
  LOGICAL                                   :: lperio   = .FALSE.               ! flag for E-W periodicity
  LOGICAL                                   :: lnc4     = .FALSE.               ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------------
  CALL ReadCdfNames() 

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfokubow -u U-file U-var -v V-file V-var [-l lev] ...'
     PRINT *,'                    ... [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute Okubo-Weiss parameter of a vector field at all level or at'  
     PRINT *,'       specified level (-l option).'
     PRINT *,'       This parameter represents the balance between strain and vorticity.'
     PRINT *,'       W = Sn^2 +Ss^2 - curl(V)^2. Sn and Ss are the non-symetrical and  '
     PRINT *,'       symetrical components of the strain, respectively.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -u U-file U-var: zonal component of the vector field:'
     PRINT *,'                 filename and variable name'
     PRINT *,'       -v V-file V-var: meridional component of the vector field:'
     PRINT *,'                 filename and variable name'
     PRINT *,'     '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-l lev]: level to be processed. Process all level by default.'
     PRINT *,'       [-o OUT-file] : specify output file instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ] : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'            This option is effective only if cdftools are compiled with'
     PRINT *,'            a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ', TRIM(cn_fhgr),' and ', TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : sokubow (s^-2)'
     STOP 
  ENDIF

  ijarg=1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-u'   ) ; CALL getarg(ijarg, cf_ufil ) ; ijarg=ijarg+1
        ;              CALL getarg(ijarg, cv_u    ) ; ijarg=ijarg+1
     CASE ( '-v'   ) ; CALL getarg(ijarg, cf_vfil ) ; ijarg=ijarg+1
        ;              CALL getarg(ijarg, cv_v    ) ; ijarg=ijarg+1
        ! option
     CASE ( '-l'   ) ; CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,* ) ilev
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out  ) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 99
     END SELECT
  ENDDO

  lchk = chkfile(cn_fhgr ) .OR. lchk
  lchk = chkfile(cn_fmsk ) .OR. lchk
  lchk = chkfile(cf_ufil ) .OR. lchk
  lchk = chkfile(cf_vfil ) .OR. lchk
  IF ( lchk ) STOP 99 ! missing files

  npiglo = getdim(cf_ufil,cn_x)
  npjglo = getdim(cf_ufil,cn_y)
  npk    = getdim(cf_ufil,cn_z)
  npt    = getdim(cf_ufil,cn_t) 

  IF ( npk==0 ) THEN
     PRINT *, 'npk=0, assume 1'
     npk=1
  END IF

  IF ( npt==0 ) THEN
     PRINT *, 'npt=0, assume 1'
     npt=1
  END IF

  PRINT *, 'npiglo = ',npiglo
  PRINT *, 'npjglo = ',npjglo
  PRINT *, 'npk    = ',npk
  PRINT *, 'npt    = ',npt
  PRINT *, 'ilev   = ',ilev

  ! Allocate the memory
  ALLOCATE ( e1u(npiglo,npjglo) , e1f(npiglo,npjglo) )
  ALLOCATE ( e2v(npiglo,npjglo) , e2f(npiglo,npjglo) )
  ALLOCATE ( e1t(npiglo,npjglo) , e2t(npiglo,npjglo) )
  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo)  )
  ALLOCATE ( strsym(npiglo,npjglo)  )
  ALLOCATE ( strnsy(npiglo,npjglo) , tmask(npiglo,npjglo) )
  ALLOCATE ( okubow(npiglo,npjglo) , fmask(npiglo,npjglo) )
  ALLOCATE ( rotn(npiglo,npjglo) , dtim(npt) , gdepu(npk) )

  e1u =  getvar(cn_fhgr, cn_ve1u, 1, npiglo, npjglo)
  e1f =  getvar(cn_fhgr, cn_ve1f, 1, npiglo, npjglo)
  e2v =  getvar(cn_fhgr, cn_ve2v, 1, npiglo, npjglo)
  e2f =  getvar(cn_fhgr, cn_ve2f, 1, npiglo, npjglo)
  e1t =  getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo)
  e2t =  getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo)

  ! use un and vn to store f latitude and longitude for CreateOutput
  un    = getvar(cn_fhgr, cn_glamf, 1, npiglo, npjglo)
  vn    = getvar(cn_fhgr, cn_gphif, 1, npiglo, npjglo)
  gdepu = getvar1d (cf_ufil, cn_vdepthu, npk )

  ! look for  E-W periodicity
  IF ( un(1,1) == un(npiglo-1,1) ) lperio = .TRUE.

  IF ( ilev < 0 ) THEN  ; ik1 = 1    ; ik2 = npk
  ELSE                  ; ik1 = ilev ; ik2 = ilev ; npk = 1 ; gdepu(1) = gdepu(ilev)
  ENDIF

  iko=0

  CALL CreateOutput

  DO jt=1,npt
     IF (MOD(jt,100)==0 ) PRINT *, jt,'/',npt
     DO jk = ik1,ik2
        iko=iko+1
        un(:,:)    = getvar(cf_ufil, cv_u    , jk, npiglo, npjglo, ktime=jt)
        vn(:,:)    = getvar(cf_vfil, cv_v    , jk, npiglo, npjglo, ktime=jt)
        tmask(:,:) = getvar(cn_fmsk, cn_tmask, jk, npiglo, npjglo          )

        ! compute the mask at level jk
        DO jj = 1, npjglo - 1
           DO ji = 1, npiglo - 1
              fmask(ji,jj)=0.
              fmask(ji,jj)= un(ji,jj)*un(ji,jj+1) * vn(ji,jj)*vn(ji+1,jj)
              IF (fmask(ji,jj) /= 0.) fmask(ji,jj)=1.
           ENDDO
        ENDDO

        rotn(:,:) = 0. ; strsym(:,:) = 0. ; strnsy(:,:) = 0. ; okubow(:,:) = 0.
        DO jj = 1, npjglo -1 
           DO ji = 1, npiglo -1   
              rotn(ji,jj) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ) - e2v(ji,jj) * vn(ji,jj)    &
                   &         - e1u(ji  ,jj+1) * un(ji  ,jj+1) + e1u(ji,jj) * un(ji,jj)  ) &
                   &         * fmask(ji,jj) / ( e1f(ji,jj) * e2f(ji,jj) )      ! quantity on f grid

              strsym(ji,jj) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ) - e2v(ji,jj) * vn(ji,jj)    &
                   &           + e1u(ji  ,jj+1) * un(ji  ,jj+1) - e1u(ji,jj) * un(ji,jj)  ) &
                   &           * fmask(ji,jj) / ( e1f(ji,jj) * e2f(ji,jj) )    ! quantity on f grid

              strnsy(ji,jj) = (  e1u(ji+1,jj  ) * un(ji+1,jj  ) - e1u(ji,jj) * un(ji,jj)    &
                   &           - e2v(ji  ,jj+1) * vn(ji  ,jj+1) + e2v(ji,jj) * vn(ji,jj)  ) &
                   &           * tmask(ji,jj) / ( e1t(ji,jj) * e2t(ji,jj) )    ! quantity on T grid
           ENDDO
        ENDDO

        DO jj = 1, npjglo -1 
           DO ji = 1, npiglo -1   
              ! compute non-symetrical squared strain on f-point
              zcoef      =  0.25 * fmask(ji,jj)
              zstrnsym2  =  zcoef * ( strnsy(ji  ,jj  ) * strnsy(ji  ,jj  ) &
                   &                + strnsy(ji+1,jj  ) * strnsy(ji+1,jj  ) &
                   &                + strnsy(ji  ,jj+1) * strnsy(ji  ,jj+1) &
                   &                + strnsy(ji+1,jj+1) * strnsy(ji+1,jj+1) )  ! quantity computed on f grid

              okubow(ji,jj) = strsym(ji,jj) * strsym(ji,jj) + zstrnsym2 - rotn(ji,jj)*rotn(ji,jj)

           END DO
        END DO

        IF ( lperio ) okubow(npiglo,:) = okubow(2, :)
        ! write okubow on file at level k and at time jt
        ierr = putvar(ncout, id_varout(1), okubow, iko,  npiglo, npjglo, ktime=jt)
     ENDDO
  END DO
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
    ! define new variables for output
    ipk(1) = npk  !  2D or 3D
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = 'sokubow'
    stypvar(1)%cunits            = 's-2'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = -1000.
    stypvar(1)%valid_max         =  1000.
    stypvar(1)%clong_name        = 'Okubo_Weiss_param (okubow)'
    stypvar(1)%cshort_name       = 'sokubow'
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TYX'

    ! create output fileset
    ncout = create      (cf_out, cf_ufil, npiglo, npjglo, npk                         , ld_nc4=lnc4)
    ierr  = createvar   (ncout , stypvar, 1,      ipk,    id_varout                   , ld_nc4=lnc4)
    ierr  = putheadervar(ncout,  cf_ufil, npiglo, npjglo, npk, pdep=gdepu,  pnavlon=un, pnavlat=vn )

    dtim = getvar1d(cf_ufil, cn_vtimec, npt      )
    ierr = putvar1d(ncout,   dtim,      npt,  'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfokubow

