PROGRAM cdfstrain
  !!======================================================================
  !!                     ***  PROGRAM  cdfstrain  ***
  !!=====================================================================
  !!  ** Purpose :  Compute the 2 component (symetric and antisymetric)
  !!                of the strain using the velocty field
  !!
  !!  ** Method  :  Just apply formula
  !!
  !! History :  4.0  : 04/2020  : F. Le Guillou, J.M. Molines from cdfokbw
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
  INTEGER(KIND=4), DIMENSION(2)             :: ipk, id_varout                   ! output variable properties

  REAL(KIND=4)                              :: zcoef                            ! working real
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2v, e1u, e1f, e2f, e1t, e2t     ! horizontal metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: un, vn                           ! velocity field
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: fmask, tmask                     ! mask
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: strsym, strnsy                   ! symetric and non symetric strain
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zstrsym                          ! Symetrical strain for output
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zstrnsy                          ! Anti-Symetrical strain for output
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: gdepu                            ! depth for output

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim                             ! time counter

  CHARACTER(LEN=256)                        :: cf_ufil, cf_vfil                 ! file names
  CHARACTER(LEN=256)                        :: cf_out = 'strain.nc'             ! output file name
  CHARACTER(LEN=256)                        :: cv_u, cv_v                       ! variable names
  CHARACTER(LEN=256)                        :: cldum                            ! dummy string

  TYPE (variable), DIMENSION(2)             :: stypvar                          ! structure for attibutes

  LOGICAL                                   :: lchk     = .FALSE.               ! flag for missing files
  LOGICAL                                   :: lperio   = .FALSE.               ! flag for E-W periodicity
  LOGICAL                                   :: lnc4     = .FALSE.               ! Use nc4 with chunking and deflation
  LOGICAL                                   :: lTpt     = .FALSE.               ! Output on T points
  LOGICAL                                   :: lFpt     = .FALSE.               ! Output on F points
  !!----------------------------------------------------------------------------
  CALL ReadCdfNames() 

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfstrain -u U-file U-var -v V-file V-var [-l lev] ...'
     PRINT *,'                    ... [-T] [-F] [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the strain of a vector field at all level or at'  
     PRINT *,'       specified level (-l option).'
     PRINT *,'       The symetric (Ss)  and antisymetric (Sn) component of the strain are '
     PRINT *,'       computed, on the native C-grid (respectively at F point and T points).'
     PRINT *,'       Options are provided for saving the 2 component of the strain at the'
     PRINT *,'       same grid point. Default behaviour is to save on native grid.'
     PRINT *,'       '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -u U-file U-var: zonal component of the vector field:'
     PRINT *,'                 filename and variable name'
     PRINT *,'       -v V-file V-var: meridional component of the vector field:'
     PRINT *,'                 filename and variable name'
     PRINT *,'     '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-l lev]: level to be processed. Process all level by default.'
     PRINT *,'       [ -T   ]: output both components on T points.'
     PRINT *,'       [ -F   ]: output both components on F points.'
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
     CASE ( '-T'   ) ; lTpt = .TRUE.
     CASE ( '-F'   ) ; lFpt = .TRUE.
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
  ALLOCATE ( strsym(npiglo,npjglo) , zstrsym(npiglo,npjglo)  )
  ALLOCATE ( strnsy(npiglo,npjglo) , zstrnsy(npiglo,npjglo)  )
  ALLOCATE ( tmask(npiglo, npjglo) ,fmask(npiglo,npjglo) )
  ALLOCATE ( dtim(npt) , gdepu(npk) )

  e1u =  getvar(cn_fhgr, cn_ve1u, 1, npiglo, npjglo)
  e1f =  getvar(cn_fhgr, cn_ve1f, 1, npiglo, npjglo)
  e1t =  getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo)
  e2v =  getvar(cn_fhgr, cn_ve2v, 1, npiglo, npjglo)
  e2f =  getvar(cn_fhgr, cn_ve2f, 1, npiglo, npjglo)
  e2t =  getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo)

  ! use un and vn to store f latitude and longitude for CreateOutput
  ! note : when computing on native grid, the 2 component of the strain are not on the same grid
  !  and lon, lat in the file (F point) are just indicative ...
  IF (lTpt ) THEN
     un    = getvar(cn_fhgr, cn_glamt, 1, npiglo, npjglo)
     vn    = getvar(cn_fhgr, cn_gphif, 1, npiglo, npjglo)
     gdepu = getvar1d (cf_ufil, cn_vdepthu, npk )

  ELSE ! lFpt or std output
     un    = getvar(cn_fhgr, cn_glamf, 1, npiglo, npjglo)
     vn    = getvar(cn_fhgr, cn_gphif, 1, npiglo, npjglo)
     gdepu = getvar1d (cf_ufil, cn_vdepthu, npk )
  ENDIF

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

        strsym(:,:) = 0. ; strnsy(:,:) = 0. ;  zstrsym(:,:) = 0. ; zstrnsy(:,:) = 0.
        DO jj = 1, npjglo -1 
           DO ji = 1, npiglo -1   
              strsym(ji,jj) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ) - e2v(ji,jj) * vn(ji,jj)    &
                   &           + e1u(ji  ,jj+1) * un(ji  ,jj+1) - e1u(ji,jj) * un(ji,jj)  ) &
                   &           * fmask(ji,jj) / ( e1f(ji,jj) * e2f(ji,jj) )    ! quantity on f grid

              strnsy(ji,jj) = (  e1u(ji+1,jj  ) * un(ji+1,jj  ) - e1u(ji,jj) * un(ji,jj)    &
                   &           - e2v(ji  ,jj+1) * vn(ji  ,jj+1) + e2v(ji,jj) * vn(ji,jj)  ) &
                   &           * tmask(ji,jj) / ( e1t(ji,jj) * e2t(ji,jj) )    ! quantity on T grid
           ENDDO
        ENDDO

        IF ( lTpt) THEN
           ! compute symetrical strain T-point
           DO jj=2, npjglo
              DO ji=2, npiglo
                 zcoef      =  0.25 * tmask(ji,jj)
                 zstrsym(ji,jj)   =  zcoef * ( strsym(ji  ,jj  ) + strsym(ji-1,jj  ) &
                      &                      + strsym(ji-1,jj-1) + strsym(ji  ,jj-1))
                 zstrnsy(ji,jj)  = strnsy(ji,jj)
              ENDDO
           ENDDO

        ELSEIF (lFpt ) THEN
           DO jj = 1, npjglo -1 
              DO ji = 1, npiglo -1   
                 ! compute non-symetrical strain on f-point
                 zcoef      =  0.25 * fmask(ji,jj)
                 zstrnsy(ji,jj)  =  zcoef * ( strnsy(ji  ,jj  ) + strnsy(ji+1,jj   ) &
                      &                     + strnsy(ji+1,jj+1) + strnsy(ji  ,jj+1 ) )
                 zstrsym(ji,jj)  = strsym(ji,jj)
              ENDDO
           ENDDO
        ELSE   ! just do nothing (consider later the use of pointer for perf)
           zstrsym(ji,jj)   = strsym(ji,jj)
           zstrnsy(ji,jj)   = strnsy(ji,jj)
        ENDIF

        IF ( lperio ) THEN
           zstrsym(npiglo,:)  = zstrsym(2, :)
           zstrnsy(npiglo,:)  = zstrnsy(2, :)
        ENDIF
        ! write okubow on file at level k and at time jt
        ierr = putvar(ncout, id_varout(1), zstrsym,  iko,  npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(2), zstrnsy,  iko,  npiglo, npjglo, ktime=jt)
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
    CHARACTER(LEN=30) :: clv_strsym, cln_strsym
    CHARACTER(LEN=30) :: clv_strnsy, cln_strnsy
    CHARACTER(LEN=30) :: cl_axis

    IF ( lTpt ) THEN
       clv_strsym='strsym_T'
       clv_strnsy='strnsy_T'
       cln_strsym='symetrical strain component on T point'
       cln_strnsy='anti-symetrical strain component on T point'
    ELSE IF ( lFpt ) THEN
       clv_strsym='strsym_F'
       clv_strnsy='strnsy_F'
       cln_strsym='symetrical strain component on F point'
       cln_strnsy='anti-symetrical strain component on F point'
    ELSE
       clv_strsym='strsym_F'
       clv_strnsy='strnsy_T'
       cln_strsym='symetrical strain component on F point'
       cln_strnsy='anti-symetrical strain component on T point'
    ENDIF

    IF (npk > 1 ) THEN
       cl_axis='TZYX'
    ELSE
       cl_axis='TYX'
    ENDIF
    ! define new variables for output
    ipk(1) = npk  !  2D or 3D
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = TRIM(clv_strsym)
    stypvar(1)%cunits            = 's-2'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = -1000.
    stypvar(1)%valid_max         =  1000.
    stypvar(1)%clong_name        = TRIM(cln_strsym)
    stypvar(1)%cshort_name       = TRIM(clv_strsym)
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = TRIM(cl_axis)

    ipk(2) = npk  !  2D or 3D
    stypvar(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname             = TRIM(clv_strnsy)
    stypvar(2)%cunits            = 's-2'
    stypvar(2)%rmissing_value    = 0.
    stypvar(2)%valid_min         = -1000.
    stypvar(2)%valid_max         =  1000.
    stypvar(2)%clong_name        = TRIM(cln_strnsy)
    stypvar(2)%cshort_name       = TRIM(clv_strnsy)
    stypvar(2)%conline_operation = 'N/A'
    stypvar(2)%caxis             = TRIM(cl_axis)

    ! create output fileset
    ncout = create      (cf_out, cf_ufil, npiglo, npjglo, npk                         , ld_nc4=lnc4)
    ierr  = createvar   (ncout , stypvar, 2,      ipk,    id_varout                   , ld_nc4=lnc4)
    ierr  = putheadervar(ncout,  cf_ufil, npiglo, npjglo, npk, pdep=gdepu,  pnavlon=un, pnavlat=vn )

    dtim = getvar1d(cf_ufil, cn_vtimec, npt      )
    ierr = putvar1d(ncout,   dtim,      npt,  'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfstrain

