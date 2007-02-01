MODULE eos

IMPLICIT NONE
PRIVATE
PUBLIC sigma0, eosbn2, sigmai

CONTAINS
     FUNCTION sigma0 ( ptem, psal, kpi,kpj)
       !! --------------------------------------------------------------------
       !! ** Purpose :   Compute the in situ density (ratio rho/rau0) and the
       !!      potential volumic mass (Kg/m3) from potential temperature and
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
 
       !! * Arguments
       INTEGER,INTENT(in) :: kpi,kpj  !: dimension of 2D arrays
       REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptem, psal
       REAL(KIND=8), DIMENSION(kpi,kpj) :: sigma0
 
       !! * local variables
       INTEGER :: ji,jj 
       REAL(KIND=8), DIMENSION (kpi,kpj) :: zws
       REAL(KIND=8) :: zt, zs, zsr, zr1, zr2, zr3, zr4, zrau0=1000.
 
 
       zws = 0.d0
       sigma0 = 0.d0
       DO jj = 1, kpj
          DO ji = 1, kpi
             zws(ji,jj) = SQRT( ABS( psal(ji,jj) ) )
          END DO
       END DO
 
       DO jj = 1, kpj
          !
          DO ji = 1, kpi
 
             zt = ptem (ji,jj)            ! interpolated T
             zs = psal (ji,jj)            ! interpolated S
             zsr= zws(ji,jj)            ! square root of interpolated S
 
             ! compute volumic mass pure water at atm pressure
             zr1 = ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt   &
                  -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
             ! seawater volumic mass atm pressure
             zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 )*zt+7.6438e-5 ) *zt   &
                  -4.0899e-3 ) *zt+0.824493
             zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
             zr4= 4.8314e-4
 
             ! potential volumic mass (reference to the surface)
             sigma0(ji,jj) = ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1 -zrau0
          END DO
       END DO
 
     END FUNCTION sigma0
     
     FUNCTION sigmai ( ptem, psal, pref, kpi,kpj)
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
 
       !! * Arguments
       INTEGER,INTENT(in) :: kpi,kpj  !: dimension of 2D arrays
       REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptem, psal
       REAL(KIND=4),                     INTENT(in) :: pref      !: reference pressure (meters or db)
       REAL(KIND=8), DIMENSION(kpi,kpj) :: sigmai
 
       REAL(kind=8),PARAMETER :: dpr4=4.8314d-4,dpd=-2.042967d-2 , dprau0 = 1000;

       !! * local variables
       INTEGER :: ji,jj 
       REAL(KIND=8), DIMENSION (kpi,kpj)  ::  dlrs
       REAL(KIND=8) ::  dlt, dls      
       REAL(KIND=8) :: dla,dla1,dlaw,dlb,dlb1,dlbw,dlc,dle,dlk0,dlkw 
       REAL(kind=8) :: dlrhop, dlr1,dlr2,dlr3, dlref
 
       dlref = pref
       sigmai = 0.d0
       DO jj = 1, kpj
          DO ji = 1, kpi
             dlrs(ji,jj) = SQRT( ABS( psal(ji,jj) ) )
          END DO
       END DO
 
       DO jj=1,kpj
         DO ji=1,kpi

! Convert T and S to double precision.
        dlt=DBLE(ptem(ji,jj))
        dls=DBLE(psal(ji,jj))


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
        sigmai(ji,jj)=dlrhop/(1.0d0-dlref/(dlk0-dlref*(dla-dlref*dlb)))&
                       -dprau0

      ENDDO
    ENDDO
 
     END FUNCTION sigmai

     FUNCTION eosbn2 ( ptem, psal, pdep,pe3w, kpi,kpj,kup,kdown)
       !! ----------------------------------------------------------------------
       !! ** Purpose :   Compute the local Brunt-Vaisala frequency at the time-
       !!      step of the input arguments
       !!
       !! ** Method :
       !!       *  UNESCO sea water properties
       !!         The brunt-vaisala frequency is computed using the
       !!      polynomial expression of McDougall (1987):
       !!            N^2 = grav * beta * ( alpha/beta*dk[ t ] - dk[ s ] )/e3w
       !!---------------------------------------------------------------------
       ! * Arguments
       INTEGER, INTENT(in)    :: kpi,kpj
       INTEGER, INTENT(in)    :: kup,kdown
       REAL(KIND=4), DIMENSION(kpi,kpj,2), INTENT(in) :: ptem, psal
       REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) ::  pe3w
       REAL(KIND=4)  :: pdep
       REAL(KIND=4), DIMENSION(kpi,kpj) :: eosbn2
 
       ! * Local variables
       INTEGER  ::   ji, jj         ! dummy loop indices
       REAL(KIND=8) ::   &
            zgde3w, zt, zs, zh,  &  ! temporary scalars
            zalbet, zbeta           !    "
       REAL(KIND=8) :: grav=9.81
 
       zh = pdep
       DO jj = 1, kpj
          DO ji = 1, kpi
             zgde3w = grav / pe3w(ji,jj)
             zt = 0.5 * ( ptem(ji,jj,kup) + ptem(ji,jj,kdown) )          ! potential temperature at w-point
             zs = 0.5 * ( psal(ji,jj,kup) + psal(ji,jj,kdown) ) - 35.0   ! salinity anomaly (s-35) at w-point
 
             zalbet = ( ( ( - 0.255019e-07 * zt + 0.298357e-05 ) * zt   &   ! ratio alpha/beta
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
 
             zbeta  = ( ( -0.415613e-09 * zt + 0.555579e-07 ) * zt      &   ! beta
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
 
             eosbn2(ji,jj) = zgde3w * zbeta            &   ! N^2
                  &          * ( zalbet * ( ptem(ji,jj,kup) - ptem(ji,jj,kdown) )   &
                  &                     - ( psal(ji,jj,kup) - psal(ji,jj,kdown) ) )
          END DO
       END DO
 
 
     END FUNCTION eosbn2

     FUNCTION albet(  ptem, psal, pdep, kpi,kpj)
       !!-------------------------------------------------------------------------------------------
       !!                 *** FUNCTION  albet ***
       !!
       !!     * Purpose: Compute the ratio alpha/beta 
       !! -----------------------------------------------------------------------------------------
       !! * Arguments
       INTEGER, INTENT(in) :: kpi, kpj
       REAL(KIND=4), DIMENSION(kpi,kpj),INTENT(in) :: ptem, psal
       REAL(KIND=4), INTENT(in) :: pdep

       REAL(KIND=8), DIMENSION(kpi,kpj) :: albet

       !! * Local variables
       INTEGER :: ji,jj
       REAL(KIND=8)  :: zt, zs, zh
       
       zh = pdep
       DO ji=1,kpi
          DO jj=1,kpj
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

     END FUNCTION albet
     
     FUNCTION beta (  ptem, psal, pdep, kpi,kpj)
       !!-------------------------------------------------------------------------------------------
       !!                 *** FUNCTION  beta ***
       !!
       !!     * Purpose: Compute the beta
       !! -----------------------------------------------------------------------------------------
       !! * Arguments
       INTEGER, INTENT(in) :: kpi, kpj
       REAL(KIND=4), DIMENSION(kpi,kpj),INTENT(in) :: ptem, psal
       REAL(KIND=4), INTENT(in) :: pdep

       REAL(KIND=8), DIMENSION(kpi,kpj) :: beta

       !! * Local variables
       INTEGER :: ji,jj
       REAL(KIND=8)  :: zt, zs, zh

       zh = pdep
       DO ji=1,kpi
          DO jj=1,kpj
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

     END FUNCTION beta

END MODULE eos
