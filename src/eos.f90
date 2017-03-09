MODULE eos
  !!======================================================================
  !!                     ***  MODULE  eos  ***
  !! All routines dealing with the Equation Of State of sea water
  !!=====================================================================
  !! History : 2.1  !  2004   : J.M. Molines : Original code ported
  !!                                           from NEMO
  !!           3.0    12/2010 : J.M. Molines : Doctor norm + Lic.
  !!           4.0  : 03/2017 : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   sigma0     : compute sigma-0 
  !!   eosbn2     : compute Brunt Vaissala Frequency
  !!   sigmai     : compute sigma-i ( refered to a depth given in argument
  !!   albet      : Compute the ratio alpha/beta ( Thermal/haline exapnsion)
  !!   beta       : compute beta (haline expension)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sigma0
  PUBLIC :: sigmai
  PUBLIC :: spice
  PUBLIC :: sigmantr
  PUBLIC :: eosbn2
  PUBLIC :: albet
  PUBLIC :: beta

  INTERFACE sigmai
     MODULE PROCEDURE sigmai_dep, sigmai_dep2d
  END INTERFACE


  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class system
  !!----------------------------------------------------------------------

CONTAINS


  FUNCTION sigma0 ( ptem, psal, kpi, kpj)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION sigma0  ***
    !!
    !! ** Purpose : Compute the in situ density (ratio rho/rau0) and the
    !!              potential volumic mass (Kg/m3) from potential temperature 
    !!              and salinity fields using an equation of state defined 
    !!              through the namelist parameter neos.
    !!
    !! ** Method  : Jackett and McDougall (1994) equation of state.
    !!              The in situ density is computed directly as a function of
    !!              potential temperature relative to the surface (the opa t
    !!              variable), salt and pressure (assuming no pressure variation
    !!              along geopotential surfaces, i.e. the pressure p in decibars
    !!              is approximated by the depth in meters.
    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
    !!              rhop(t,s)  = rho(t,s,0)
    !!              with pressure                      p        decibars
    !!                   potential temperature         t        deg celsius
    !!                   salinity                      s        psu
    !!                   reference volumic mass        rau0     kg/m**3
    !!                   in situ volumic mass          rho      kg/m**3
    !!                   in situ density anomalie      prd      no units
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptem, psal  ! temperature and salinity
    INTEGER(KIND=4),                  INTENT(in) :: kpi, kpj    ! dimension of 2D arrays
    REAL(KIND=8), DIMENSION(kpi,kpj)             :: sigma0      ! returned value

    INTEGER(KIND=4)                   :: ji, jj 
    REAL(KIND=8), DIMENSION (kpi,kpj) :: zws
    REAL(KIND=8)                      :: zt, zs, zsr, zrau0=1000.
    REAL(KIND=8)                      :: zr1, zr2, zr3, zr4
    !!----------------------------------------------------------------------
    zws = 0.d0
    sigma0 = 0.d0
    DO jj = 1, kpj
       DO ji = 1, kpi
          zws(ji,jj) = SQRT( ABS( psal(ji,jj) ) )
       END DO
    END DO

    !$OMP PARALLEL DO SCHEDULE(RUNTIME)
    DO jj = 1, kpj
       DO ji = 1, kpi

          zt  = ptem (ji,jj)          ! interpolated T
          zs  = psal (ji,jj)          ! interpolated S
          zsr = zws  (ji,jj)          ! square root of interpolated S

          ! compute volumic mass pure water at atm pressure
          zr1 = ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt   &
               -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
          ! seawater volumic mass atm pressure
          zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 )*zt+7.6438e-5 ) *zt   &
               -4.0899e-3 ) *zt+0.824493
          zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
          zr4= 4.8314e-4

          ! potential volumic mass (reference to the surface)
          sigma0(ji,jj) = ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1 - zrau0
       END DO
    END DO
    !$OMP END PARALLEL DO

  END FUNCTION sigma0

  FUNCTION sigmantr( ptem, psal, kpi,kpj) 
    !! --------------------------------------------------------------------
    !! ** Purpose :   Compute the  neutral volumic mass (kg/m3) from known
    !!      potential temperature and salinity fields using an equation of
    !!      state. 
    !!
    !! ** Method  :
    !!       McDougall and Jackett (2005) equation of state.
    !!              potential temperature         t        deg celsius
    !!              salinity                      s        psu
    !!              neutral density               rho      kg/m**3
    !!       result is not masked at this stage.
    !!       Check value: rho(20,35) = 1024.59416751197 kg/m**3  -1000.
    !!       t = 20 deg celcius, s=35 psu
    !!
    !! ** References : McDougall and Jackett, The material derivative of neutral density
    !!        Journal of Marine Research, 63, 159-185, 2005
    !! --------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptem, psal ! temperature salinity
    INTEGER(KIND=4),                  INTENT(in) :: kpi,kpj    ! dimension of 2D arrays
    REAL(KIND=8), DIMENSION(kpi,kpj) :: sigmantr          ! return value

    INTEGER(KIND=4)                   :: ji, jj
    REAL(KIND=8), DIMENSION (kpi,kpj) :: dl_ws
    REAL(KIND=8)                      :: dl_t, dl_s, dl_sr
    REAL(KIND=8)                      :: dl_r1, dl_r2, dl_r3, dl_r4, dl_r5
    !! --------------------------------------------------------------------
    !$OMP PARALLEL DO SCHEDULE(RUNTIME)
    DO jj = 1, kpj
       DO ji = 1, kpi
          dl_t = ptem(ji,jj)
          dl_s = psal(ji,jj)
          dl_sr= SQRT( ABS(dl_s) )
       ! Numerator
          ! T-Polynome
          dl_r1= ( ( -4.3159255086706703d-4*dl_t+8.1157118782170051d-2 )*dl_t+2.2280832068441331d-1 )*dl_t+1002.3063688892480d0
          ! S-T Polynome
          dl_r2= ( -1.7052298331414675d-7*dl_s-3.1710675488863952d-3*dl_t-1.0304537539692924d-4 )*dl_s
       ! Denominator
          ! T-Polynome
          dl_r3= ( ( (-2.3850178558212048d-9*dl_t -1.6212552470310961d-7 )*dl_t+7.8717799560577725d-5 )*dl_t+4.3907692647825900d-5 )*dl_t+     1.0d0
          ! S-T Polynome
          dl_r4= ( ( -2.2744455733317707d-9*dl_t*dl_t+6.0399864718597388d-6)*dl_t-5.1268124398160734d-4 )*dl_s
          ! S-T Polynome
          dl_r5= ( -1.3409379420216683d-9*dl_t*dl_t-3.6138532339703262d-5)*dl_s*dl_sr

          ! Neutral density
          sigmantr(ji,jj) = ( dl_r1 + dl_r2 ) / ( dl_r3 + dl_r4 + dl_r5 ) -1000.d0
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO


  END FUNCTION sigmantr

  FUNCTION spice ( ptem, psal, kpi, kpj )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION spice  ***
    !!
    !! ** Purpose :  Compute spiciness from T S fields. 
    !!
    !! ** Method  :  spiciness = sum(i=0,5)[sum(j=0,4)[b(i,j)*theta^i*(s-35)^j]]
    !!                   with:  b     -> coefficients
    !!                          theta -> potential temperature
    !!                          s     -> salinity
    !!
    !!  **  Example:
    !!       spice(15,33)=   0.5445863      0.544586321373410  calcul en double
    !!       spice(15,33)=   0.5445864      (calcul en simple precision) 
    !!
    !!  ** References : Flament (2002) "A state variable for characterizing
    !!              water masses and their diffusive stability: spiciness."
    !!              Progress in Oceanography Volume 54, 2002, Pages 493-501.
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptem, psal ! temperature salinity
    INTEGER(KIND=4),                  INTENT(in) :: kpi,kpj    ! dimension of 2D arrays
    REAL(KIND=8), DIMENSION(kpi,kpj)             :: spice      ! return value

    INTEGER(KIND=4)                        :: ji,jj     ! dummy loop index
    INTEGER(KIND=4)                        :: jig,jjg     ! dummy loop index
    REAL(KIND=8), SAVE, DIMENSION(6,5)     :: dl_bet    ! coefficients of spiciness formula
    REAL(KIND=8),       DIMENSION(kpi,kpj) :: dl_spi    ! spiciness
    REAL(KIND=8),       DIMENSION(kpi,kpj) :: dl_salref ! reference salinity
    REAL(KIND=8),       DIMENSION(kpi,kpj) :: dl_tempt  ! working array
    REAL(KIND=8),       DIMENSION(kpi,kpj) :: dl_salt   ! working array

    LOGICAL, SAVE                        :: lfrst=.TRUE.!
    !!----------------------------------------------------------------------
    IF ( lfrst ) THEN
       lfrst = .false.
       ! Define coefficients to compute spiciness (R*8)
       dl_bet(1,:) = (/       0.d0,   7.7442d-01,    -5.85d-03,    -9.84d-04,    -2.06d-04/)
       dl_bet(2,:) = (/ 5.1655d-02,    2.034d-03,   -2.742d-04,     -8.5d-06,     1.36d-05/)
       dl_bet(3,:) = (/6.64783d-03,  -2.4681d-04,   -1.428d-05,    3.337d-05,    7.894d-06/)
       dl_bet(4,:) = (/-5.4023d-05,    7.326d-06,   7.0036d-06,  -3.0412d-06,  -1.0853d-06/)
       dl_bet(5,:) = (/  3.949d-07,   -3.029d-08,  -3.8209d-07,   1.0012d-07,   4.7133d-08/)
       dl_bet(6,:) = (/  -6.36d-10,   -1.309d-09,    6.048d-09,  -1.1409d-09,   -6.676d-10/)
    ENDIF
    ! spiciness 
    dl_spi(:,:)    = 0.d0
    dl_salref(:,:) = psal(:,:) - 35.d0
    dl_tempt(:,:)  = 1.d0
!$OMP PARALLEL DO SCHEDULE(RUNTIME)
    DO jjg=1, kpj
    DO ji=1,6
       dl_salt(:,jjg) = 1.d0
       DO jj=1,5
           dl_spi( :,jjg) = dl_spi (:,jjg) +   dl_bet (ji,jj) * dl_tempt(:,jjg) * dl_salt(:,jjg)
           dl_salt(:,jjg) = dl_salt(:,jjg) * dl_salref( :,jjg )
       END DO
       dl_tempt(:,jjg) = dl_tempt(:,jjg) * ptem(:,jjg)
    END DO
    END DO
!$OMP END PARALLEL DO

    spice(:,:) = dl_spi( :,:)

  END FUNCTION spice

  FUNCTION sigmai_dep ( ptem, psal, pref, kpi,kpj)
    !! --------------------------------------------------------------------
    !! ** Purpose :   Compute the  density referenced to pref (ratio rho/rau0) 
    !!       from potential temperature and
    !!      salinity fields using an equation of state defined through the
    !!     namelist parameter neos.
    !!
    !! ** Method  :
    !!       Jackett and McDougall (1994) equation of state.
    !!         the in situ density is computed directly as a function of
    !!         potential temperature relative to the surface (the opa t
    !!         variable), salt and pressure (assuming no pressure variation
    !!         along geopotential surfaces, i.e. the pressure p in decibars
    !!         is approximated by the depth in meters.
    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
    !!              rhop(t,s)  = rho(t,s,0)
    !!         with pressure                      p        decibars
    !!              potential temperature         t        deg celsius
    !!              salinity                      s        psu
    !!              reference volumic mass        rau0     kg/m**3
    !!              in situ volumic mass          rho      kg/m**3
    !!              in situ density anomalie      prd      no units
    !! --------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptem, psal ! temperature salinity
    INTEGER(KIND=4),                  INTENT(in) :: kpi,kpj    ! dimension of 2D arrays
    REAL(KIND=4),                     INTENT(in) :: pref       ! reference pressure (meters or db)
    REAL(KIND=8), DIMENSION(kpi,kpj) :: sigmai_dep             ! return value

    REAL(kind=8), PARAMETER :: dpr4=4.8314d-4, dpd=-2.042967d-2 , dprau0 = 1000.d0

    INTEGER(KIND=4) :: ji, jj 
    REAL(KIND=8), DIMENSION (kpi,kpj) ::  dlrs
    REAL(KIND=8) :: dlt, dls      
    REAL(KIND=8) :: dla, dla1, dlaw, dlb, dlb1, dlbw, dlc, dle, dlk0, dlkw 
    REAL(kind=8) :: dlrhop, dlr1, dlr2, dlr3, dlref
    !! --------------------------------------------------------------------
    dlref      = pref
    sigmai_dep = 0.d0
    DO jj = 1, kpj
       DO ji = 1, kpi
          dlrs(ji,jj) = SQRT( ABS( psal(ji,jj) ) )
       END DO
    END DO

    !$OMP PARALLEL DO SCHEDULE(RUNTIME)
    DO jj=1,kpj
       DO ji=1,kpi

          ! Convert T and S to double precision.
          dlt = DBLE(ptem(ji,jj))
          dls = DBLE(psal(ji,jj))

          ! Compute the volumic mass of pure water at atmospheric pressure.
          dlr1=((((6.536332d-9*dlt-1.120083d-6)&
               *dlt+1.001685d-4)&
               *dlt-9.095290d-3)&
               *dlt+6.793952d-2)&
               *dlt+999.842594d0

          ! Compute the seawater volumic mass at atmospheric pressure.
          dlr2=(((5.3875d-9*dlt-8.2467d-7)&
               *dlt+7.6438d-5)&
               *dlt-4.0899d-3)&
               *dlt+0.824493d0

          dlr3=(-1.6546d-6*dlt+1.0227d-4)&
               *dlt-5.72466d-3

          ! Compute the potential volumic mass (referenced to the surface).
          dlrhop=(dpr4*dls+dlr3*dlrs(ji,jj)+dlr2)*dls+dlr1

          ! Compute the compression terms.
          dle=(-3.508914d-8*dlt-1.248266d-8)&
               *dlt-2.595994d-6

          dlbw=(1.296821d-6*dlt-5.782165d-9)&
               *dlt+1.045941d-4

          dlb=dlbw+dle*dls

          dlc=(-7.267926d-5*dlt+2.598241d-3 )&
               *dlt+0.1571896d0

          dlaw=((5.939910d-6*dlt+2.512549d-3)&
               *dlt-0.1028859d0)&
               *dlt-4.721788d0

          dla=(dpd*dlrs(ji,jj)+dlc)*dls+dlaw

          dlb1=(-0.1909078d0*dlt+7.390729d0)&
               *dlt-55.87545d0

          dla1=((2.326469d-3*dlt+1.553190d0)&
               *dlt-65.00517d0)&
               *dlt+1044.077d0

          dlkw=(((-1.361629d-4*dlt-1.852732d-2)&
               *dlt-30.41638d0)&
               *dlt+2098.925d0)&
               *dlt+190925.6d0

          dlk0=(dlb1*dlrs(ji,jj)+dla1)*dls+dlkw

          ! Compute the potential density anomaly.
          sigmai_dep(ji,jj)=dlrhop/(1.0d0-dlref/(dlk0-dlref*(dla-dlref*dlb)))&
               -dprau0

       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END FUNCTION sigmai_dep

  FUNCTION sigmai_dep2d ( ptem, psal, pref, kpi,kpj)
    !! --------------------------------------------------------------------
    !! ** Purpose :   Compute the  density referenced to pref (ratio rho/rau0) 
    !!       from potential temperature and
    !!      salinity fields using an equation of state defined through the
    !!     namelist parameter neos.
    !!
    !! ** Method  :
    !!       Jackett and McDougall (1994) equation of state.
    !!         the in situ density is computed directly as a function of
    !!         potential temperature relative to the surface (the opa t
    !!         variable), salt and pressure (assuming no pressure variation
    !!         along geopotential surfaces, i.e. the pressure p in decibars
    !!         is approximated by the depth in meters.
    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
    !!              rhop(t,s)  = rho(t,s,0)
    !!         with pressure                      p        decibars
    !!              potential temperature         t        deg celsius
    !!              salinity                      s        psu
    !!              reference volumic mass        rau0     kg/m**3
    !!              in situ volumic mass          rho      kg/m**3
    !!              in situ density anomalie      prd      no units
    !! --------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptem, psal ! temperature salinity
    INTEGER(KIND=4),                  INTENT(in) :: kpi,kpj    ! dimension of 2D arrays
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: pref       ! reference pressure (meters or db) (2d Array)
    REAL(KIND=8), DIMENSION(kpi,kpj) :: sigmai_dep2d           ! return value

    REAL(kind=8), PARAMETER :: dpr4=4.8314d-4, dpd=-2.042967d-2 , dprau0 = 1000.d0

    INTEGER(KIND=4) :: ji, jj 
    REAL(KIND=8), DIMENSION (kpi,kpj) ::  dlrs
    REAL(KIND=8) :: dlt, dls      
    REAL(KIND=8) :: dla, dla1, dlaw, dlb, dlb1, dlbw, dlc, dle, dlk0, dlkw 
    REAL(kind=8) :: dlrhop, dlr1, dlr2, dlr3, dlref

    sigmai_dep2d = 0.d0
    DO jj = 1, kpj
       DO ji = 1, kpi
          dlrs(ji,jj) = SQRT( ABS( psal(ji,jj) ) )
       END DO
    END DO

    !$OMP PARALLEL DO SCHEDULE(RUNTIME)
    DO jj=1,kpj
       DO ji=1,kpi

          ! Convert T and S to double precision.
          dlt   = DBLE(ptem(ji,jj))
          dls   = DBLE(psal(ji,jj))
          dlref = DBLE(pref(ji,jj))

          ! Compute the volumic mass of pure water at atmospheric pressure.
          dlr1=((((6.536332d-9*dlt-1.120083d-6)&
               *dlt+1.001685d-4)&
               *dlt-9.095290d-3)&
               *dlt+6.793952d-2)&
               *dlt+999.842594d0

          ! Compute the seawater volumic mass at atmospheric pressure.
          dlr2=(((5.3875d-9*dlt-8.2467d-7)&
               *dlt+7.6438d-5)&
               *dlt-4.0899d-3)&
               *dlt+0.824493d0

          dlr3=(-1.6546d-6*dlt+1.0227d-4)&
               *dlt-5.72466d-3

          ! Compute the potential volumic mass (referenced to the surface).
          dlrhop=(dpr4*dls+dlr3*dlrs(ji,jj)+dlr2)*dls+dlr1

          ! Compute the compression terms.
          dle=(-3.508914d-8*dlt-1.248266d-8)&
               *dlt-2.595994d-6

          dlbw=(1.296821d-6*dlt-5.782165d-9)&
               *dlt+1.045941d-4

          dlb=dlbw+dle*dls

          dlc=(-7.267926d-5*dlt+2.598241d-3 )&
               *dlt+0.1571896d0

          dlaw=((5.939910d-6*dlt+2.512549d-3)&
               *dlt-0.1028859d0)&
               *dlt-4.721788d0

          dla=(dpd*dlrs(ji,jj)+dlc)*dls+dlaw

          dlb1=(-0.1909078d0*dlt+7.390729d0)&
               *dlt-55.87545d0

          dla1=((2.326469d-3*dlt+1.553190d0)&
               *dlt-65.00517d0)&
               *dlt+1044.077d0

          dlkw=(((-1.361629d-4*dlt-1.852732d-2)&
               *dlt-30.41638d0)&
               *dlt+2098.925d0)&
               *dlt+190925.6d0

          dlk0=(dlb1*dlrs(ji,jj)+dla1)*dls+dlkw

          ! Compute the potential density anomaly.
          sigmai_dep2d(ji,jj)=dlrhop/(1.0d0-dlref/(dlk0-dlref*(dla-dlref*dlb)))&
               -dprau0

       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END FUNCTION sigmai_dep2d


  FUNCTION eosbn2 ( ptem, psal, pdep, pe3w, kpi, kpj, kup, kdown)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION eosbn2  ***
    !!
    !! ** Purpose :  Compute the local Brunt-Vaisala frequency at the time-
    !!               step of the input arguments
    !!
    !! ** Method :  UNESCO sea water properties
    !!             The brunt-vaisala frequency is computed using the
    !!              polynomial expression of McDougall (1987):
    !!              N^2 = grav * beta * ( alpha/beta*dk[ t ] - dk[ s ] )/e3w
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj,2), INTENT(in) :: ptem, psal ! temperaature salinity
    REAL(KIND=4)                                   :: pdep       ! reference depth
    REAL(KIND=4), DIMENSION(kpi,kpj),   INTENT(in) :: pe3w       ! e3w of the current layer
    INTEGER(KIND=4),                    INTENT(in) :: kpi, kpj   ! size of the array
    INTEGER(KIND=4),                    INTENT(in) :: kup, kdown ! index of levels up and down
    REAL(KIND=4), DIMENSION(kpi,kpj)               :: eosbn2     ! returned values

    INTEGER(KIND=4) :: ji, jj         ! dummy loop indices
    REAL(KIND=8)    :: zgde3w, zt, zs, zh
    REAL(KIND=8)    :: zalbet, zbeta
    REAL(KIND=8)    :: zgrav=9.81
    !!----------------------------------------------------------------------

    zh = pdep
    !$OMP PARALLEL DO SCHEDULE(RUNTIME)
    DO jj = 1, kpj
       DO ji = 1, kpi
          zgde3w = zgrav / pe3w(ji,jj)
          zt = 0.5 * ( ptem(ji,jj,kup) + ptem(ji,jj,kdown) )          ! potential temperature at w-point
          zs = 0.5 * ( psal(ji,jj,kup) + psal(ji,jj,kdown) ) - 35.0   ! salinity anomaly (s-35) at w-point

          zalbet = ( ( ( - 0.255019e-07 * zt + 0.298357e-05 ) * zt     &   ! ratio alpha/beta
               &                               - 0.203814e-03 ) * zt   &
               &                               + 0.170907e-01 ) * zt   &
               &   + 0.665157e-01                                      &
               &   +     ( - 0.678662e-05 * zs                         &
               &           - 0.846960e-04 * zt + 0.378110e-02 ) * zs   &
               &   +   ( ( - 0.302285e-13 * zh                         &
               &           - 0.251520e-11 * zs                         &
               &           + 0.512857e-12 * zt * zt           ) * zh   &
               &           - 0.164759e-06 * zs                         &
               &        +(   0.791325e-08 * zt - 0.933746e-06 ) * zt   &
               &                               + 0.380374e-04 ) * zh

          zbeta  = ( ( -0.415613e-09 * zt + 0.555579e-07 ) * zt        &   ! beta
               &                            - 0.301985e-05 ) * zt      &
               &   + 0.785567e-03                                      &
               &   + (     0.515032e-08 * zs                           &
               &         + 0.788212e-08 * zt - 0.356603e-06 ) * zs     &
               &   +(  (   0.121551e-17 * zh                           &
               &         - 0.602281e-15 * zs                           &
               &         - 0.175379e-14 * zt + 0.176621e-12 ) * zh     &
               &                             + 0.408195e-10   * zs     &
               &     + ( - 0.213127e-11 * zt + 0.192867e-09 ) * zt     &
               &                             - 0.121555e-07 ) * zh

          eosbn2(ji,jj) = zgde3w * zbeta                                         &   ! N^2
               &          * ( zalbet * ( ptem(ji,jj,kup) - ptem(ji,jj,kdown) )   &
               &                     - ( psal(ji,jj,kup) - psal(ji,jj,kdown) ) )
       END DO
    END DO
    !$OMP END PARALLEL DO

  END FUNCTION eosbn2


  FUNCTION albet(  ptem, psal, pdep, kpi, kpj)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION  albet  ***
    !!
    !! ** Purpose :  Compute the ratio alpha/beta 
    !!
    !! ** Method  :  Follow Mc Dougal et al as in other functions
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptem, psal ! temperature salinity
    REAL(KIND=4),                     INTENT(in) :: pdep       ! refererence depth
    INTEGER(KIND=4),                  INTENT(in) :: kpi, kpj   ! size of the arrays

    REAL(KIND=8), DIMENSION(kpi,kpj)             :: albet      ! returned value

    INTEGER(KIND=4) :: ji, jj      ! dummy loop index
    REAL(KIND=8)    :: zt, zs, zh  ! working local variables
    !!----------------------------------------------------------------------
    zh = pdep
    !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     DO jj=1,kpj
       DO ji=1,kpi
          zt =  ptem(ji,jj)         ! potential temperature
          zs =  psal(ji,jj)- 35.0   ! salinity anomaly (s-35)

          albet(ji,jj) = ( ( ( - 0.255019e-07 * zt + 0.298357e-05 ) * zt   &   ! ratio alpha/beta
               &                               - 0.203814e-03 ) * zt   &
               &                               + 0.170907e-01 ) * zt   &
               &   + 0.665157e-01                                      &
               &   +     ( - 0.678662e-05 * zs                         &
               &           - 0.846960e-04 * zt + 0.378110e-02 ) * zs   &
               &   +   ( ( - 0.302285e-13 * zh                         &
               &           - 0.251520e-11 * zs                         &
               &           + 0.512857e-12 * zt * zt           ) * zh   &
               &           - 0.164759e-06 * zs                         &
               &        +(   0.791325e-08 * zt - 0.933746e-06 ) * zt   &
               &                               + 0.380374e-04 ) * zh
       END DO
    END DO
    !$OMP END PARALLEL DO

  END FUNCTION albet


  FUNCTION beta (  ptem, psal, pdep, kpi, kpj)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION beta  ***
    !!
    !! ** Purpose :  Compute the beta 
    !!
    !! ** Method  :  Follow Mc Dougal et al as in other functions 
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj),INTENT(in) :: ptem, psal ! temperature salinity
    REAL(KIND=4),                    INTENT(in) :: pdep       ! reference depth
    INTEGER(KIND=4),                 INTENT(in) :: kpi, kpj   ! size of the array
    REAL(KIND=8), DIMENSION(kpi,kpj)            :: beta       ! returned values

    INTEGER(KIND=4) :: ji, jj      ! dummy loop index
    REAL(KIND=8)    :: zt, zs, zh  ! working variables
    !!----------------------------------------------------------------------
    zh = pdep
    !$OMP PARALLEL DO SCHEDULE(RUNTIME)
    DO jj=1,kpj
       DO ji=1,kpi
          zt =  ptem(ji,jj)         ! potential temperature
          zs =  psal(ji,jj)- 35.0   ! salinity anomaly (s-35)

          beta(ji,jj)  = ( ( -0.415613e-09 * zt + 0.555579e-07 ) * zt      &   ! beta
               &                            - 0.301985e-05 ) * zt      &
               &   + 0.785567e-03                                      &
               &   + (     0.515032e-08 * zs                           &
               &         + 0.788212e-08 * zt - 0.356603e-06 ) * zs     &
               &   +(  (   0.121551e-17 * zh                           &
               &         - 0.602281e-15 * zs                           &
               &         - 0.175379e-14 * zt + 0.176621e-12 ) * zh     &
               &                             + 0.408195e-10   * zs     &
               &     + ( - 0.213127e-11 * zt + 0.192867e-09 ) * zt     &
               &                             - 0.121555e-07 ) * zh
       END DO
    END DO
    !$OMP END PARALLEL DO

  END FUNCTION beta

END MODULE eos
