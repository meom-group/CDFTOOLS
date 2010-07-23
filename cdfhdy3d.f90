PROGRAM cdfhdy3d
  !!-------------------------------------------------------------------
  !!             ***  PROGRAM cdfhdy3d  ***
  !!
  !!  **  Purpose: Compute dynamical height anomaly field from gridT file
  !!                Store the results on a 3D cdf file.
  !!  
  !!
  !! history: 
  !!     Original :   J.M. Molines (Nov 2004 ) for ORCA025
  !!                  J.M. Molines Apr 2005 : use modules
  !!-------------------------------------------------------------------
  !!  $Rev: 256 $
  !!  $Date: 2009-07-21 17:49:27 +0200 (mar. 21 juil. 2009) $
  !!  $Id: cdfsig0.f90 256 2009-07-21 15:49:27Z molines $
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio
  USE eos

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk, jt                              !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: 
  INTEGER   :: npiglo,npjglo, npk, npt             !: size of the domain
  INTEGER, DIMENSION(1) ::  ipk, &                 !: outptut variables : number of levels,
       &                    id_varout              !: ncdf varid's
  real(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: ztemp, zsal ,&   !: Array to read a layer of data
       &                                         ztemp0, zsal0 ,&   !: reference density
       &                                         zsig0 , &        !: potential density (sig-0)
       &                                         zsig  , &        !: potential density (sig-0)
       &                                         zmask , &        !: 2D mask at current level
       &                                         zhdy, zterm, zdep, zdepth, zssh
  REAL(KIND=4),DIMENSION(:),ALLOCATABLE   ::  tim, ze3t_1d

  CHARACTER(LEN=256) :: cfilet , cdum, cfileout='cdfhdy3d.nc', cmask='mask.nc' !:
  CHARACTER(LEN=256) :: coordzgr='mesh_zgr.nc'

  TYPE(variable) , DIMENSION(1) :: typvar         !: structure for attributes

  INTEGER    :: ncout
  INTEGER    :: istatus
  INTEGER    :: zlev1, zlev2
  INTEGER, DIMENSION (2) :: ismin, ismax
  REAL(KIND=4)   :: sigmin, sigmax, rau0=1000.


  !!  Read command line
  narg= iargc()
  IF ( narg .NE. 1 ) THEN
     PRINT *,' Usage : cdfhdy3d gridT '
     PRINT *,' reference is the sea surface, mask.nc and mesh_zgr.nc must be in your directory'
     PRINT *,' Output on cdfhdy3d.nc, variable vohdy'
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')
  npt   = getdim (cfilet,'time')

  ipk(:)= npk 
  typvar(1)%name= 'vohdy'
  typvar(1)%units='m'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= -100.
  typvar(1)%valid_max= 100.
  typvar(1)%long_name='Dynamical height anomaly'
  typvar(1)%short_name='vohdy'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'npt   =', npt

  ALLOCATE (ztemp0(npiglo,npjglo), zsal0(npiglo,npjglo), zsig0(npiglo,npjglo) ,zmask(npiglo,npjglo))
  ALLOCATE (ztemp(npiglo,npjglo), zsal(npiglo,npjglo), zsig(npiglo,npjglo) , zhdy(npiglo,npjglo), zterm(npiglo,npjglo))
  ALLOCATE (zdep(npiglo,npjglo), zdepth(npiglo,npjglo), zssh(npiglo,npjglo), ze3t_1d(npk))
  ALLOCATE (tim(npt))

  ! create output fileset

  ncout =create(cfileout, cfilet, npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfilet,npiglo, npjglo,npk)
  tim=getvar1d(cfilet,'time_counter',npt)
  ierr=putvar1d(ncout,tim,npt,'T')

  ! Temperature and salinity for reference profile
  ztemp0(:,:)=0.
  zsal0(:,:)=35.

  zssh(:,:)  = getvar(cfilet, 'sossheig', 1, npiglo, npjglo)
  ze3t_1d(:) = getvare3(coordzgr, 'e3t',npk)

  DO jt=1,npt
    PRINT *,' TIME = ', jt, tim(jt)/86400.,' days'

     zhdy(:,:) = 0.
     zdepth(:,:) = 0.

  DO jk = 1, npk

     !zdep(:,:)   = getvar(coordzgr, 'e3t_ps', jk,npiglo,npjglo,ldiom=.true.)
     ! we degrade the computation to smooth the results
     zdep(:,:) = ze3t_1d(jk)

     zmask(:,:) = getvar(cmask, 'tmask', jk, npiglo, npjglo)

     IF ( jk == 1) THEN
        zdep(:,:) = zdep(:,:) + zssh(:,:)
     ENDIF

     ! total depth at current level (used for computation of rho in situ)
     zdepth(:,:) = zdepth(:,:) + zdep(:,:)

     ztemp(:,:)= getvar(cfilet, 'votemper',  jk ,npiglo, npjglo,ktime=jt)
     zsal(:,:) = getvar(cfilet, 'vosaline',  jk ,npiglo, npjglo,ktime=jt)

     CALL eos_insitu( ztemp0, zsal0, zdepth, npiglo, npjglo, zsig0 ) 
     CALL eos_insitu( ztemp, zsal, zdepth, npiglo, npjglo, zsig ) 

     ! we compute the term of the integral : (1/g) *10e4 * sum [ delta * dz ]
     ! with delta = (1/rho - 1/rho0)
     ! 10e4 factor is conversion decibar/pascal
     !
     zterm = ( ( 1. / ( rau0 + zsig(:,:) ) ) - ( 1. / ( rau0 + zsig0(:,:) ) ) ) * 10000. * zdep / 9.81
     ! in land, it seems appropriate to stop the computation
     WHERE(zsal == 0 ) zterm = 0

     zhdy(:,:) = zhdy(:,:) + zterm(:,:)
     ! we mask with the mask of the level
     zhdy(:,:) = zhdy(:,:) * zmask(:,:)

     ierr = putvar(ncout, id_varout(1) ,zhdy, jk,npiglo, npjglo,ktime=jt)


  END DO  ! loop to next level
     
  END DO  ! next time frame

  istatus = closeout(ncout)

CONTAINS

SUBROUTINE eos_insitu( ptem, psal, pdepth, jpiglo, jpjglo, prd )
     !!----------------------------------------------------------------------
     !!                   ***  ROUTINE eos_insitu  ***
     !! 
     !! ** Purpose :   Compute the in situ density (ratio rho/rau0) from 
     !!       potential temperature and salinity using an equation of state
     !!       defined through the namelist parameter nn_eos.
     !!
     !! ** Method  : 
     !!      nn_eos = 0 : Jackett and McDougall (1994) equation of state.
     !!         the in situ density is computed directly as a function of
     !!         potential temperature relative to the surface (the opa t
     !!         variable), salt and pressure (assuming no pressure variation
     !!         along geopotential surfaces, i.e. the pressure p in decibars
     !!         is approximated by the depth in meters.
     !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
     !!         with pressure                      p        decibars
     !!              potential temperature         t        deg celsius
     !!              salinity                      s        psu
     !!              reference volumic mass        rau0     kg/m**3
     !!              in situ volumic mass          rho      kg/m**3
     !!              in situ density anomalie      prd      no units
     !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
     !!          t = 40 deg celcius, s=40 psu
     !!              prd(t,s) = rn_beta * s - rn_alpha * tn - 1.
     !!      Note that no boundary condition problem occurs in this routine
     !!      as (ptem,psal) are defined over the whole domain.
     !!
     !! ** Action  :   compute prd , the in situ density (no units)
     !!
     !! References :   Jackett and McDougall, J. Atmos. Ocean. Tech., 1994
     !!----------------------------------------------------------------------
     INTEGER, INTENT(in   )                           ::   jpiglo, jpjglo
     REAL(8), DIMENSION(jpiglo,jpjglo), INTENT(in   ) ::   ptem   ! potential temperature  [Celcius]
     REAL(8), DIMENSION(jpiglo,jpjglo), INTENT(in   ) ::   psal   ! salinity               [psu]
     REAL(8), DIMENSION(jpiglo,jpjglo), INTENT(in   ) ::   pdepth ! depth                  [m]
     REAL(8), DIMENSION(jpiglo,jpjglo), INTENT(  out) ::   prd    ! in situ density 
     !!
     INTEGER  ::   ji, jj, jk           ! dummy loop indices
     INTEGER  ::   jpkm1
     REAL(8) ::   zt , zs , zh , zsr   ! temporary scalars
     REAL(8) ::   zr1, zr2, zr3, zr4   !    -         -
     REAL(8) ::   zrhop, ze, zbw, zb   !    -         -
     REAL(8) ::   zd , zc , zaw, za    !    -         -
     REAL(8) ::   zb1, za1, zkw, zk0   !    -         -
     REAL(8) ::   zrau0r               !    -         -
     REAL(8), DIMENSION(jpiglo,jpjglo) ::   zws   ! temporary workspace
     INTEGER  ::   nn_eos   = 0        !: = 0/1/2 type of eq. of state and Brunt-Vaisala frequ.
     REAL(8) ::   rn_alpha = 2.0e-4   !: thermal expension coeff. (linear equation of state)
     REAL(8) ::   rn_beta  = 7.7e-4   !: saline  expension coeff. (linear equation of state)

     REAL(8) ::   ralpbet           !: alpha / beta ratio
      !!----------------------------------------------------------------------

     zrau0r = 1.e0 / rau0
     zws(:,:) = SQRT( ABS( psal(:,:) ) )
     !  
        DO jj = 1, jpjglo
           DO ji = 1, jpiglo
              zt = ptem  (ji,jj)
              zs = psal  (ji,jj)
              zh = pdepth(ji,jj)        ! depth
              zsr= zws   (ji,jj)        ! square root salinity
              !
              ! compute volumic mass pure water at atm pressure
              zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt   &
                 &      -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
              ! seawater volumic mass atm pressure
              zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt        &
                 &                   -4.0899e-3 ) *zt+0.824493
              zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
              zr4= 4.8314e-4
              !
              ! potential volumic mass (reference to the surface)
              zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
              !
              ! add the compression terms
              ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
              zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
              zb = zbw + ze * zs
              !
              zd = -2.042967e-2
              zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
              zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt - 4.721788
              za = ( zd*zsr + zc ) *zs + zaw
              !
              zb1=   (-0.1909078*zt+7.390729 ) *zt-55.87545
              za1= ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077
              zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt + 2098.925 ) *zt+190925.6
              zk0= ( zb1*zsr + za1 )*zs + zkw
              !
              ! masked in situ density anomaly
              prd(ji,jj) = (  zrhop / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  )    &
                 &             - rau0  ) ! * zrau0r ! * tmask(ji,jj)
           END DO
        END DO
END SUBROUTINE eos_insitu

END PROGRAM cdfhdy3d
