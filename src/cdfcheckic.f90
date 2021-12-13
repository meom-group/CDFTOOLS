PROGRAM cdfcheckic
  !!======================================================================
  !!                     ***  PROGRAM  cdfcheckic  ***
  !!=====================================================================
  !!  ** Purpose : apply a non penetrative convection scheme
  !!               using same algoritm than NEMO
  !!
  !!  ** Method  : when static instability detected, mix the cell above with the cell down
  !!               the process is repeated until no more static unstability
  !!
  !!  ** Comments: 3d array are used as the NEMO algo works column-wise and most of the 
  !!               netcdf used in NEMO are chunk 2d along x,y (map access)
  !!
  !! History : 4.0  : 09/2020 : P. Mathiot : initial code
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames   ! for cdf variable names
  USE eos
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class Equation_of_state
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                              :: jk, jt                   ! dummy loop index
  INTEGER(KIND=4)                              :: jpi, jpj, jpk
  INTEGER(KIND=4)                              :: ierr                     ! working integer
  INTEGER(KIND=4)                              :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                              :: npiglo, npjglo, npk, npt ! size of the domain
  INTEGER(KIND=4)                              :: ncout                    ! ncid of output file
  INTEGER(KIND=4), DIMENSION(2)                :: ipk, id_varout           ! level and id of output variables
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mbkt, mikt               ! top and bottom level

  REAL(KIND=4), DIMENSION (:,:,:), ALLOCATABLE :: ztem, zsal, e3w,e3t, wmask, tmask   ! 3d array
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: tmskutil
  REAL(KIND=4), DIMENSION (:,:,:), ALLOCATABLE :: gdept, gdepw             ! depth
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE :: gdept_1d                 ! depth

  REAL(KIND=8), DIMENSION (:),     ALLOCATABLE :: dtim                     ! time

  CHARACTER(LEN=256)                           :: cf_tfil, cldum, cv_dep   ! input file name, ...
  CHARACTER(LEN=256)                           :: cf_sfil                  ! salinity file if different from temperature file
  CHARACTER(LEN=256)                           :: cf_out = 'cdfic.nc'      ! output file name
  CHARACTER(LEN=256)                           :: cglobal                  ! global attribute
  CHARACTER(LEN=80)                            :: cf_e3w                   ! file with e3w in case of vvl

  TYPE(variable), DIMENSION(2)                 :: stypvar                  ! variable attribute

  LOGICAL                                      :: lchk                     ! check missing files
  LOGICAL                                      :: lnc4  =.FALSE.           ! full step flag
  LOGICAL                                      :: ll_teos10  = .FALSE.     ! teos10 flag

  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfcheckic -t T-file [-s S-file] [-o OUT-file] [-nc4] [-teos10]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Apply the non penetrative convection scheme to the initial condition to' 
     PRINT *,'       avoid static instability'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file : netcdf input gridT file for temperature and salinity.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-s S-file] : netcdf input for salinity if not in T-file.'
     PRINT *,'       [-o OUT-file ] : specify output file name instead of ',TRIM(cf_out),'.'
     PRINT *,'       [-nc4 ]  : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
     PRINT *,'       [-teos10] : use TEOS10 equation of state instead of default EOS80'
     PRINT *,'                 Temperature should be conservative temperature (CT) in deg C.'
     PRINT *,'                 Salinity should be absolute salinity (SA) in g/kg.'

     PRINT *,'      '
     PRINT *,'     OPENMP SUPPORT : yes'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fzgr),' is needed for this program.' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option specified'
     PRINT *,'      '
     PRINT *,'    SEE ALSO :'
     PRINT *,'       cdfbn2'
     STOP 
  ENDIF

  cglobal = 'Partial step computation'
  cf_sfil = 'none'

  ijarg = 1
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE (cldum)
     CASE ( '-t'      ) ; CALL getarg(ijarg, cf_tfil) ; ijarg = ijarg + 1
        ! options
     CASE ( '-s'      ) ; CALL getarg(ijarg, cf_sfil) ; ijarg = ijarg + 1
     CASE ( '-o'      ) ; CALL getarg(ijarg, cf_out ) ; ijarg = ijarg + 1
     CASE ( '-nc4'    ) ; lnc4    = .TRUE.
     CASE ( '-teos10' ) ; ll_teos10 = .TRUE. 
     CASE DEFAULT       ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  CALL eos_init( ll_teos10 )

  IF ( cf_sfil == 'none') cf_sfil = cf_tfil


  lchk = chkfile (cn_fzgr )
  lchk = lchk .OR. chkfile (cn_fhgr  )
  lchk = lchk .OR. chkfile (cf_tfil  )
  lchk = lchk .OR. chkfile (cf_sfil  )
  IF ( lchk  ) STOP 99  ! missing files 

  npiglo = getdim (cf_tfil, cn_x); jpi=npiglo
  npjglo = getdim (cf_tfil, cn_y); jpj=npjglo
  npk    = getdim (cf_tfil, cn_z); jpk=npk
  npt    = getdim (cf_tfil, cn_t)

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  ALLOCATE (ztem (npiglo,npjglo,npk), zsal(npiglo,npjglo,npk) )
  ALLOCATE (wmask(npiglo,npjglo,npk), tmask(npiglo,npjglo,npk))
  ALLOCATE (e3w(npiglo,npjglo,npk),e3t(npiglo,npjglo,npk))
  ALLOCATE (tmskutil(npiglo,npjglo) )
  ALLOCATE (gdepw(npiglo,npjglo,npk), gdept(npiglo,npjglo,npk), gdept_1d(npk), dtim(npt)        )
  ALLOCATE (mikt(npiglo,npjglo), mbkt(npiglo,npjglo) )

  gdept_1d = getvare3(cn_fzgr, cn_gdept, npk)
  CALL CreateOutput

  tmskutil(:,:) = getvar(cn_fmsk, cn_tmaskutil, 1, npiglo, npjglo)
  mbkt(:,:) = getvar(cn_fzgr, 'mbathy', 1, npiglo, npjglo)
  mikt(:,:) = getvar(cn_fzgr, 'misf'  , 1, npiglo, npjglo)

  ! load needed static data (wmask used as temporary array here)
  PRINT *, 'load scale factor, depth ...'
  PRINT *, cg_zgr_ver
  IF ( cg_zgr_ver == 'v3.6') cn_ve3w = 'e3w_0'
  IF ( cg_zgr_ver == 'v3.6') cn_ve3t = 'e3t_0'
  e3w(:,:,:)   = getvar3d(cn_fzgr, cn_ve3w  , npiglo, npjglo, npk)
  e3t(:,:,:)   = getvar3d(cn_fzgr, cn_ve3t  , npiglo, npjglo, npk)
  tmask(:,:,:) = getvar3d(cn_fmsk, cn_tmask , npiglo, npjglo, npk)
  gdepw(:,:,1) = 0.
  gdept(:,:,1) = e3t(:,:,1)/2.
  DO jk = 2, npk
    gdepw(:,:,jk) = gdepw(:,:,jk-1) + e3t(:,:,jk-1)
    gdept(:,:,jk) = gdept(:,:,jk-1) + e3w(:,:,jk)
  ENDDO
!  gdepw(:,:,:) = getvar3d(cn_fzgr, cn_depw3d, npiglo, npjglo, npk)
! gdept(:,:,:) = getvar3d(cn_fzgr, cn_dept3d, npiglo, npjglo, npk)

  ! wmask is computed from tmask
  DO jk = 2,npk 
     wmask(:,:,jk) = tmask(:,:,jk) * tmask(:,:,jk-1)
  END DO

  ! start temporal loop
  DO jt=1,npt
     PRINT *,jt,'/',npt

     !  2  compute 3d t/s/wmask
     ztem (:,:,:) = getvar3d(cf_tfil, cn_votemper, npiglo, npjglo, npk, ktime=jt)
     zsal (:,:,:) = getvar3d(cf_sfil, cn_vosaline, npiglo, npjglo, npk, ktime=jt)

     ! 3 compute alpha, beta and bn2
     CALL tra_npc( ztem, zsal )

     ! 4 write data
     DO jk = 1,npk 
        ierr = putvar(ncout, id_varout(1), ztem(:,:,jk)*tmask(:,:,jk), jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout, id_varout(2), zsal(:,:,jk)*tmask(:,:,jk), jk, npiglo, npjglo, ktime=jt )
     END DO

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
    ipk(1)                       = npk  !  3D
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = 'votemper'
    stypvar(1)%cunits            = 'C'
    stypvar(1)%rmissing_value    = -9999.
    stypvar(1)%valid_min         = -10.
    stypvar(1)%valid_max         = 50.
    stypvar(1)%clong_name        = 'ocean temperature'
    stypvar(1)%cshort_name       = 'ocean temperature'
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TZYX'

    ipk(2)                       = npk  !  3D
    stypvar(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname             = 'vosaline'
    stypvar(2)%cunits            = 'C'
    stypvar(2)%rmissing_value    = -9999.
    stypvar(2)%valid_min         = 0.
    stypvar(2)%valid_max         = 60.
    stypvar(2)%clong_name        = 'ocean salinity'
    stypvar(2)%cshort_name       = 'ocean salinity'
    stypvar(2)%conline_operation = 'N/A'
    stypvar(2)%caxis             = 'TZYX'


    ! create output fileset
    ncout = create      (cf_out,   cf_tfil,  npiglo, npjglo, npk,                               ld_nc4=lnc4 )
    ierr  = createvar   (ncout ,   stypvar,  2,      ipk,    id_varout, cdglobal=TRIM(cglobal), ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,    cf_tfil,  npiglo, npjglo, npk, pdep=gdept_1d)

    dtim = getvar1d(cf_tfil, cn_vtimec, npt    )
    ierr = putvar1d(ncout,  dtim,       npt,'T')

  END SUBROUTINE CreateOutput

   SUBROUTINE tra_npc( pt, ps )
      !!----------------------------------------------------------------------
      !!  FROM NEMO NEMO/OCE 4.0
      !!                  ***  ROUTINE tranpc  ***
      !!
      !! ** Purpose : Non-penetrative convective adjustment scheme. solve
      !!      the static instability of the water column on after fields
      !!      while conserving heat and salt contents.
      !!
      !! ** Method  : updated algorithm able to deal with non-linear equation of state
      !!              (i.e. static stability computed locally)
      !!
      !! ** Action  : - tsa: after tracers with the application of the npc scheme
      !!              - send the associated trends for on-line diagnostics (l_trdtra=T)
      !!
      !! References :     Madec, et al., 1991, JPO, 21, 9, 1349-1371.
      !!----------------------------------------------------------------------
      !
      REAL(4 ), DIMENSION(jpi,jpj,jpk), INTENT(inout) :: pt, ps          ! active tracers and RHS of tracer equation
      REAL(4 ), DIMENSION(jpk) :: zgdept, zgdepw
      !
      INTEGER  ::   ji, jj, jk                                               ! dummy loop indices
      INTEGER  ::   jiter, iktop, ikbot, ikp, ikup, ikdown, ilayer, ik_low   ! local integers
      INTEGER  ::   jp_tem=1, jp_sal=2                                       ! temp and salinity index
      LOGICAL  ::   l_bottom_reached, l_column_treated                       ! flag to detect if the work is done
      REAL(4 ) ::   zta, zalfa, zsum_temp, zsum_alfa, zaw, zdz, zsum_z
      REAL(4 ) ::   zsa, zbeta, zsum_sali, zsum_beta, zbw, zrw, z1_rDt
      REAL(4 ), PARAMETER ::   zn2_zero = 1.e-14                ! acceptance criteria for neutrality (N2==0)
      REAL(4 ), PARAMETER ::   grav = 9.80665                   ! physics cst
      REAL(4 ), DIMENSION(        jpk     )   ::   zvn2         ! vertical profile of N2 at 1 given point...
      REAL(4 ), DIMENSION(        jpk, 2  )   ::   zvts, zvab   ! vertical profile of T & S , and  alpha & betaat 1 given point
      !!----------------------------------------------------------------------
      !
      DO jj = 1,jpj
         DO ji = 1,jpi
            !
            ! 1 compute alpha/beta
            !
            IF( tmskutil(ji,jj) == 1 ) THEN      ! At least 2 ocean points
               !                                     ! consider one ocean column 
               !
               ! load T/S
               zvts(:,jp_tem) = pt(ji,jj,:)      ! temperature
               zvts(:,jp_sal) = ps(ji,jj,:)      ! salinity
               !
               ! load depth
               zgdept(:) = gdept(ji,jj,:)
               zgdepw(:) = gdepw(ji,jj,:)
               !
               ! compute beta and alpha at T point
               zvab(:,jp_sal)= beta ( zvts(:,jp_tem), zvts(:,jp_sal), zgdept(:), jpk )
               zvab(:,jp_tem)= albet( zvts(:,jp_tem), zvts(:,jp_sal), zgdept(:), jpk ) * zvab(:,jp_sal)
               !
               ! compute N2 at W point, doing exactly as in eosbn2.F90:
               DO jk = 2,jpk
                  ! 
                  ! Interpolating alfa and beta at W point:
                  zrw =  (zgdepw(jk  ) - zgdept(jk)) &
                     & / (zgdept(jk-1) - zgdept(jk))
                  zaw = zvab(jk,jp_tem) * (1.    - zrw) + zvab(jk-1,jp_tem) * zrw * wmask(ji,jj,jk)
                  zbw = zvab(jk,jp_sal) * (1.    - zrw) + zvab(jk-1,jp_sal) * zrw * wmask(ji,jj,jk)
                  !
                  ! compute N2
                  zvn2(jk) = grav*( zaw * ( zvts(jk-1,jp_tem) - zvts(jk,jp_tem) )     &
                              &   - zbw * ( zvts(jk-1,jp_sal) - zvts(jk,jp_sal) )  )  &
                              &   / e3w(ji,jj,jk) * wmask(ji,jj,jk)
               END DO
               !
               ikbot = mbkt(ji,jj)   ! ikbot: ocean bottom T-level
               iktop = mikt(ji,jj)   ! iktop: ocean top    T-level
               ikp = iktop           ! because N2 is irrelevant at the surface level (will start at ikp=2)
               ilayer = 0
               jiter  = 0
               l_column_treated = .FALSE.
               !
               ! Iteration until no more static instability
               DO WHILE ( .NOT. l_column_treated )
                  !
                  jiter = jiter + 1
                  ! 
                  IF( jiter >= 400 ) STOP 'tra_npc :  PROBLEM #0 (tto much iteration)'
                  !
                  l_bottom_reached = .FALSE.
                  !
                  DO WHILE ( .NOT. l_bottom_reached )
                     !
                     ikp = ikp + 1
                     !
                     !! Testing level ikp for instability
                     !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                     IF( zvn2(ikp) <  -zn2_zero ) THEN ! Instability found!
                        !
                        ilayer = ilayer + 1    ! yet another instable portion of the water column found....
                        !
                        !! ikup is the uppermost point where mixing will start:
                        ikup = ikp - 1 ! ikup is always "at most at ikp-1", less if neutral levels overlying
                        !
                        !! If the points above ikp-1 have N2 == 0 they must also be mixed:
                        IF( ikp > iktop+1 ) THEN
                           DO jk = ikp-1, iktop+1, -1
                              IF( ABS(zvn2(jk)) < zn2_zero ) THEN
                                 ikup = ikup - 1  ! 1 more upper level has N2=0 and must be added for the mixing
                              ELSE
                                 EXIT
                              ENDIF
                           END DO
                        ENDIF
                        
                        IF( iktop < iktop )  STOP 'tra_npc :  PROBLEM #1'
                        !
                        zsum_temp = 0.   
                        zsum_sali = 0.   
                        zsum_alfa = 0.   
                        zsum_beta = 0.   
                        zsum_z    = 0.   
                        DO jk = ikup, ikbot     ! Inside the instable (and overlying neutral) portion of the column
                           !
                           zdz       = e3t(ji,jj,jk)
                           zsum_temp = zsum_temp + zvts(jk,jp_tem)*zdz
                           zsum_sali = zsum_sali + zvts(jk,jp_sal)*zdz
                           zsum_alfa = zsum_alfa + zvab(jk,jp_tem)*zdz
                           zsum_beta = zsum_beta + zvab(jk,jp_sal)*zdz
                           zsum_z    = zsum_z    + zdz
                           !                              
                           IF( jk == ikbot ) EXIT ! avoid array-index overshoot in case ikbot = jpk, cause we're calling jk+1 next line
                           !! EXIT when we have reached the last layer that is instable (N2<0) or neutral (N2=0):
                           IF( zvn2(jk+1) > zn2_zero ) EXIT
                        END DO
                       
                        ikdown = jk ! for the current unstable layer, ikdown is the deepest point with a negative or neutral N2
                        IF( ikup == ikdown )   STOP 'tra_npc :  PROBLEM #2'
                        !
                        ! Mixing Temperature, salinity, alpha and beta from ikup to ikdown included:
                        zta   = zsum_temp/zsum_z
                        zsa   = zsum_sali/zsum_z
                        zalfa = zsum_alfa/zsum_z
                        zbeta = zsum_beta/zsum_z
   
                        !! Homogenaizing the temperature, salinity, alpha and beta in this portion of the column
                        DO jk = ikup, ikdown
                           zvts(jk,jp_tem) = zta
                           zvts(jk,jp_sal) = zsa
                           zvab(jk,jp_tem) = zalfa
                           zvab(jk,jp_sal) = zbeta
                        END DO

                        !! Updating N2 in the relvant portion of the water column
                        !! Temperature, Salinity, Alpha and Beta have been homogenized in the unstable portion
                        !! => Need to re-compute N2! will use Alpha and Beta!
                        
                        ikup   = MAX(iktop+1,ikup)         ! ikup can never be 1 !
                        ik_low = MIN(ikdown+1,ikbot) ! we must go 1 point deeper than ikdown!
                        
                        DO jk = ikup, ik_low              ! we must go 1 point deeper than ikdown!
   
                           !! Interpolating alfa and beta at W point:
                           zrw =  (zgdepw(jk  ) - zgdept(jk)) &
                              & / (zgdept(jk-1) - zgdept(jk))
                           zaw = zvab(jk,jp_tem) * (1.    - zrw) + zvab(jk-1,jp_tem) * zrw * wmask(ji,jj,jk)
                           zbw = zvab(jk,jp_sal) * (1.    - zrw) + zvab(jk-1,jp_sal) * zrw * wmask(ji,jj,jk)
   
                           !! N2 at W point, doing exactly as in eosbn2.F90:
                           zvn2(jk) = grav*( zaw * ( zvts(jk-1,jp_tem) - zvts(jk,jp_tem) )     &
                              &            - zbw * ( zvts(jk-1,jp_sal) - zvts(jk,jp_sal) )  )  &
                              &       / e3w(ji,jj,jk) * wmask(ji,jj,jk)
   
                        END DO
                     
                        ikp = MIN(ikdown+1,ikbot)
                        
   
                     ENDIF  !IF( zvn2(ikp) < 0. )
   
                     IF( ikp == ikbot ) l_bottom_reached = .TRUE.
                     !
                  END DO ! DO WHILE ( .NOT. l_bottom_reached )
   
                  IF( ikp /= ikbot )   STOP 'tra_npc :  PROBLEM #3'
                 
                  ! ******* At this stage ikp == ikbot ! *******
                 
                  IF( ilayer > 0 ) THEN      !! least an unstable layer has been found
                     !
                     ikp    = 1     ! starting again at the surface for the next iteration
                     ilayer = 0
                     !
                  ENDIF
                  !
                  IF( ikp >= ikbot )   l_column_treated = .TRUE.
                  !
               END DO ! DO WHILE ( .NOT. l_column_treated )
   
               !! Updating pts:
               pt(ji,jj,:) = zvts(:,jp_tem)
               ps(ji,jj,:) = zvts(:,jp_sal)
   
            ENDIF ! IF( tmskutil == 1 ) THEN
   
         END DO
      END DO
      !
   END SUBROUTINE tra_npc

END PROGRAM
