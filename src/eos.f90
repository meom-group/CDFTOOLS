MODULE eos
  !!======================================================================
  !!                     ***  MODULE  eos  ***
  !! All routines dealing with the Equation Of State of sea water
  !!=====================================================================
  !! History : 2.1  !  2004   : J.M. Molines : Original code ported
  !!                                           from NEMO
  !!           3.0    12/2010 : J.M. Molines : Doctor norm + Lic.
  !!           4.0  : 03/2017 : J.M. Molines  
  !!           4.2  : 10/2021 : J.M. Molines  Code TEOS10 from NEMO 4
  !!                               ( keep same API )
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   sigma0     : compute sigma-0 
  !!   eosbn2     : compute Brunt Vaissala Frequency
  !!   sigmai     : compute sigma-i ( refered to a depth given in argument
  !!   albet      : Compute the ratio alpha/beta ( Thermal/haline expansion)
  !!   beta       : compute beta (haline expansion)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: eos_init
  PUBLIC :: sigma0
  PUBLIC :: sigmai
  PUBLIC :: spice
  PUBLIC :: sigmantr  ! no TEOS10 so far
  PUBLIC :: eosbn2
  PUBLIC :: albet
  PUBLIC :: alpha
  PUBLIC :: beta

  INTERFACE sigmai
     MODULE PROCEDURE sigmai_dep, sigmai_dep2d, sigmai_dep_80, sigmai_dep2d_80
  END INTERFACE

  INTERFACE spice
     MODULE PROCEDURE spice_pol,  spice_80
  END INTERFACE

  INTERFACE eosbn2
     MODULE PROCEDURE eosbn2_pol,  eosbn2_80
  END INTERFACE

  INTERFACE albet
     MODULE PROCEDURE  albet_2d, albet_z
  END INTERFACE

  INTERFACE alpha
     MODULE PROCEDURE  alpha_0d, alpha_2d, alpha_z
  END INTERFACE

  INTERFACE beta
     MODULE PROCEDURE  beta_0d, beta_2d, beta_z, beta_2d_80, beta_z_80
  END INTERFACE

   INTEGER(KIND=4), PARAMETER :: wp=8
   ! constant
   REAL(KIND=8) :: r1_rau0, rau0, dlU
   ! TEOS10/EOS80 parameters
   REAL(wp) ::   r1_S0, r1_T0, r1_Z0, rdeltaS

   ! REF peofile parameters
   REAL(wp) ::   R00, R01, R02, R03, R04, R05
   ! EOS parameters
   REAL(wp) ::   EOS000 , EOS100 , EOS200 , EOS300 , EOS400 , EOS500 , EOS600
   REAL(wp) ::   EOS010 , EOS110 , EOS210 , EOS310 , EOS410 , EOS510
   REAL(wp) ::   EOS020 , EOS120 , EOS220 , EOS320 , EOS420
   REAL(wp) ::   EOS030 , EOS130 , EOS230 , EOS330
   REAL(wp) ::   EOS040 , EOS140 , EOS240
   REAL(wp) ::   EOS050 , EOS150
   REAL(wp) ::   EOS060
   REAL(wp) ::   EOS001 , EOS101 , EOS201 , EOS301 , EOS401
   REAL(wp) ::   EOS011 , EOS111 , EOS211 , EOS311
   REAL(wp) ::   EOS021 , EOS121 , EOS221
   REAL(wp) ::   EOS031 , EOS131
   REAL(wp) ::   EOS041
   REAL(wp) ::   EOS002 , EOS102 , EOS202
   REAL(wp) ::   EOS012 , EOS112
   REAL(wp) ::   EOS022
   REAL(wp) ::   EOS003 , EOS103
   REAL(wp) ::   EOS013

   ! ALPHA parameters
   REAL(wp) ::   ALP000 , ALP100 , ALP200 , ALP300 , ALP400 , ALP500
   REAL(wp) ::   ALP010 , ALP110 , ALP210 , ALP310 , ALP410
   REAL(wp) ::   ALP020 , ALP120 , ALP220 , ALP320
   REAL(wp) ::   ALP030 , ALP130 , ALP230
   REAL(wp) ::   ALP040 , ALP140
   REAL(wp) ::   ALP050
   REAL(wp) ::   ALP001 , ALP101 , ALP201 , ALP301
   REAL(wp) ::   ALP011 , ALP111 , ALP211
   REAL(wp) ::   ALP021 , ALP121
   REAL(wp) ::   ALP031
   REAL(wp) ::   ALP002 , ALP102
   REAL(wp) ::   ALP012
   REAL(wp) ::   ALP003

   ! BETA parameters
   REAL(wp) ::   BET000 , BET100 , BET200 , BET300 , BET400 , BET500
   REAL(wp) ::   BET010 , BET110 , BET210 , BET310 , BET410
   REAL(wp) ::   BET020 , BET120 , BET220 , BET320
   REAL(wp) ::   BET030 , BET130 , BET230
   REAL(wp) ::   BET040 , BET140
   REAL(wp) ::   BET050
   REAL(wp) ::   BET001 , BET101 , BET201 , BET301
   REAL(wp) ::   BET011 , BET111 , BET211
   REAL(wp) ::   BET021 , BET121
   REAL(wp) ::   BET031
   REAL(wp) ::   BET002 , BET102
   REAL(wp) ::   BET012
   REAL(wp) ::   BET003

   ! PEN parameters
   REAL(wp) ::   PEN000 , PEN100 , PEN200 , PEN300 , PEN400
   REAL(wp) ::   PEN010 , PEN110 , PEN210 , PEN310
   REAL(wp) ::   PEN020 , PEN120 , PEN220
   REAL(wp) ::   PEN030 , PEN130
   REAL(wp) ::   PEN040
   REAL(wp) ::   PEN001 , PEN101 , PEN201
   REAL(wp) ::   PEN011 , PEN111
   REAL(wp) ::   PEN021
   REAL(wp) ::   PEN002 , PEN102
   REAL(wp) ::   PEN012

   ! ALPHA_PEN parameters
   REAL(wp) ::   APE000 , APE100 , APE200 , APE300
   REAL(wp) ::   APE010 , APE110 , APE210
   REAL(wp) ::   APE020 , APE120
   REAL(wp) ::   APE030
   REAL(wp) ::   APE001 , APE101
   REAL(wp) ::   APE011
   REAL(wp) ::   APE002

   ! BETA_PEN parameters
   REAL(wp) ::   BPE000 , BPE100 , BPE200 , BPE300
   REAL(wp) ::   BPE010 , BPE110 , BPE210
   REAL(wp) ::   BPE020 , BPE120
   REAL(wp) ::   BPE030
   REAL(wp) ::   BPE001 , BPE101
   REAL(wp) ::   BPE011
   REAL(wp) ::   BPE002

   ! Spiciness-0 from GSW
INTEGER(KIND=4), PARAMETER :: r8=8
real (r8), parameter :: s01 = -9.22982898371678e1_r8
real (r8), parameter :: s02 = -1.35727873628866e1_r8
real (r8), parameter :: s03 =  1.87353650994010e1_r8
real (r8), parameter :: s04 = -1.61360047373455e1_r8
real (r8), parameter :: s05 =  3.76112762286425e1_r8
real (r8), parameter :: s06 = -4.27086671461257e1_r8
real (r8), parameter :: s07 =  2.00820111041594e1_r8
real (r8), parameter :: s08 =  2.87969717584045e2_r8
real (r8), parameter :: s09 =  1.13747111959674e1_r8
real (r8), parameter :: s10 =  6.07377192990680e1_r8
real (r8), parameter :: s11 = -7.37514033570187e1_r8
real (r8), parameter :: s12 = -7.51171878953574e1_r8
real (r8), parameter :: s13 =  1.63310989721504e2_r8
real (r8), parameter :: s14 = -8.83222751638095e1_r8
real (r8), parameter :: s15 = -6.41725302237048e2_r8
real (r8), parameter :: s16 =  2.79732530789261e1_r8
real (r8), parameter :: s17 = -2.49466901993728e2_r8
real (r8), parameter :: s18 =  3.26691295035416e2_r8
real (r8), parameter :: s19 =  2.66389243708181e1_r8
real (r8), parameter :: s20 = -2.93170905757579e2_r8
real (r8), parameter :: s21 =  1.76053907144524e2_r8
real (r8), parameter :: s22 =  8.27634318120224e2_r8
real (r8), parameter :: s23 = -7.02156220126926e1_r8
real (r8), parameter :: s24 =  3.82973336590803e2_r8
real (r8), parameter :: s25 = -5.06206828083959e2_r8
real (r8), parameter :: s26 =  6.69626565169529e1_r8
real (r8), parameter :: s27 =  3.02851235050766e2_r8
real (r8), parameter :: s28 = -1.96345285604621e2_r8
real (r8), parameter :: s29 = -5.74040806713526e2_r8
real (r8), parameter :: s30 =  7.03285905478333e1_r8
real (r8), parameter :: s31 = -2.97870298879716e2_r8
real (r8), parameter :: s32 =  3.88340373735118e2_r8
real (r8), parameter :: s33 = -8.29188936089122e1_r8
real (r8), parameter :: s34 = -1.87602137195354e2_r8
real (r8), parameter :: s35 =  1.27096944425793e2_r8
real (r8), parameter :: s36 =  2.11671167892147e2_r8
real (r8), parameter :: s37 = -3.15140919876285e1_r8
real (r8), parameter :: s38 =  1.16458864953602e2_r8
real (r8), parameter :: s39 = -1.50029730802344e2_r8
real (r8), parameter :: s40 =  3.76293848660589e1_r8
real (r8), parameter :: s41 =  6.47247424373200e1_r8
real (r8), parameter :: s42 = -4.47159994408867e1_r8
real (r8), parameter :: s43 = -3.23533339449055e1_r8
real (r8), parameter :: s44 =  5.30648562097667_r8
real (r8), parameter :: s45 = -1.82051249177948e1_r8
real (r8), parameter :: s46 =  2.33184351090495e1_r8
real (r8), parameter :: s47 = -6.22909903460368_r8
real (r8), parameter :: s48 = -9.55975464301446_r8
real (r8), parameter :: s49 =  6.61877073960113_r8


  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class system
  !!----------------------------------------------------------------------

CONTAINS
  SUBROUTINE eos_init (ld_teos10 )
   LOGICAL,  INTENT(in) :: ld_teos10

!  rau0        = 1026._wp
   rau0        = 1000._wp
   r1_rau0     = 1._wp / rau0
  ! Reference profile coef, common to TEOS10 and EOS80 ( F. Roquet pers. comm )
  R00 =  4.6494977072d+01
  R01 = -5.2099962525d+00
  R02 =  2.2601900708d-01
  R03 =  6.4326772569d-02
  R04 =  1.5616995503d-02
  R05 = -1.7243708991d-03

   IF ( ld_teos10 ) THEN ! TEOS10 equation 
         PRINT *, ' Equation of state for sea-water : TEOS-10 '
         rdeltaS = 32._wp
         r1_S0  = 0.875_wp/35.16504_wp
         r1_T0  = 1._wp/40._wp
         r1_Z0  = 1.e-4_wp
         !
         EOS000 = 8.0189615746e+02_wp
         EOS100 = 8.6672408165e+02_wp
         EOS200 = -1.7864682637e+03_wp
         EOS300 = 2.0375295546e+03_wp
         EOS400 = -1.2849161071e+03_wp
         EOS500 = 4.3227585684e+02_wp
         EOS600 = -6.0579916612e+01_wp
         EOS010 = 2.6010145068e+01_wp
         EOS110 = -6.5281885265e+01_wp
         EOS210 = 8.1770425108e+01_wp
         EOS310 = -5.6888046321e+01_wp
         EOS410 = 1.7681814114e+01_wp
         EOS510 = -1.9193502195_wp
         EOS020 = -3.7074170417e+01_wp
         EOS120 = 6.1548258127e+01_wp
         EOS220 = -6.0362551501e+01_wp
         EOS320 = 2.9130021253e+01_wp
         EOS420 = -5.4723692739_wp
         EOS030 = 2.1661789529e+01_wp
         EOS130 = -3.3449108469e+01_wp
         EOS230 = 1.9717078466e+01_wp
         EOS330 = -3.1742946532_wp
         EOS040 = -8.3627885467_wp
         EOS140 = 1.1311538584e+01_wp
         EOS240 = -5.3563304045_wp
         EOS050 = 5.4048723791e-01_wp
         EOS150 = 4.8169980163e-01_wp
         EOS060 = -1.9083568888e-01_wp
         EOS001 = 1.9681925209e+01_wp
         EOS101 = -4.2549998214e+01_wp
         EOS201 = 5.0774768218e+01_wp
         EOS301 = -3.0938076334e+01_wp
         EOS401 = 6.6051753097_wp
         EOS011 = -1.3336301113e+01_wp
         EOS111 = -4.4870114575_wp
         EOS211 = 5.0042598061_wp
         EOS311 = -6.5399043664e-01_wp
         EOS021 = 6.7080479603_wp
         EOS121 = 3.5063081279_wp
         EOS221 = -1.8795372996_wp
         EOS031 = -2.4649669534_wp
         EOS131 = -5.5077101279e-01_wp
         EOS041 = 5.5927935970e-01_wp
         EOS002 = 2.0660924175_wp
         EOS102 = -4.9527603989_wp
         EOS202 = 2.5019633244_wp
         EOS012 = 2.0564311499_wp
         EOS112 = -2.1311365518e-01_wp
         EOS022 = -1.2419983026_wp
         EOS003 = -2.3342758797e-02_wp
         EOS103 = -1.8507636718e-02_wp
         EOS013 = 3.7969820455e-01_wp
         !
         ALP000 = -6.5025362670e-01_wp
         ALP100 = 1.6320471316_wp
         ALP200 = -2.0442606277_wp
         ALP300 = 1.4222011580_wp
         ALP400 = -4.4204535284e-01_wp
         ALP500 = 4.7983755487e-02_wp
         ALP010 = 1.8537085209_wp
         ALP110 = -3.0774129064_wp
         ALP210 = 3.0181275751_wp
         ALP310 = -1.4565010626_wp
         ALP410 = 2.7361846370e-01_wp
         ALP020 = -1.6246342147_wp
         ALP120 = 2.5086831352_wp
         ALP220 = -1.4787808849_wp
         ALP320 = 2.3807209899e-01_wp
         ALP030 = 8.3627885467e-01_wp
         ALP130 = -1.1311538584_wp
         ALP230 = 5.3563304045e-01_wp
         ALP040 = -6.7560904739e-02_wp
         ALP140 = -6.0212475204e-02_wp
         ALP050 = 2.8625353333e-02_wp
         ALP001 = 3.3340752782e-01_wp
         ALP101 = 1.1217528644e-01_wp
         ALP201 = -1.2510649515e-01_wp
         ALP301 = 1.6349760916e-02_wp
         ALP011 = -3.3540239802e-01_wp
         ALP111 = -1.7531540640e-01_wp
         ALP211 = 9.3976864981e-02_wp
         ALP021 = 1.8487252150e-01_wp
         ALP121 = 4.1307825959e-02_wp
         ALP031 = -5.5927935970e-02_wp
         ALP002 = -5.1410778748e-02_wp
         ALP102 = 5.3278413794e-03_wp
         ALP012 = 6.2099915132e-02_wp
         ALP003 = -9.4924551138e-03_wp
         !
         BET000 = 1.0783203594e+01_wp
         BET100 = -4.4452095908e+01_wp
         BET200 = 7.6048755820e+01_wp
         BET300 = -6.3944280668e+01_wp
         BET400 = 2.6890441098e+01_wp
         BET500 = -4.5221697773_wp
         BET010 = -8.1219372432e-01_wp
         BET110 = 2.0346663041_wp
         BET210 = -2.1232895170_wp
         BET310 = 8.7994140485e-01_wp
         BET410 = -1.1939638360e-01_wp
         BET020 = 7.6574242289e-01_wp
         BET120 = -1.5019813020_wp
         BET220 = 1.0872489522_wp
         BET320 = -2.7233429080e-01_wp
         BET030 = -4.1615152308e-01_wp
         BET130 = 4.9061350869e-01_wp
         BET230 = -1.1847737788e-01_wp
         BET040 = 1.4073062708e-01_wp
         BET140 = -1.3327978879e-01_wp
         BET050 = 5.9929880134e-03_wp
         BET001 = -5.2937873009e-01_wp
         BET101 = 1.2634116779_wp
         BET201 = -1.1547328025_wp
         BET301 = 3.2870876279e-01_wp
         BET011 = -5.5824407214e-02_wp
         BET111 = 1.2451933313e-01_wp
         BET211 = -2.4409539932e-02_wp
         BET021 = 4.3623149752e-02_wp
         BET121 = -4.6767901790e-02_wp
         BET031 = -6.8523260060e-03_wp
         BET002 = -6.1618945251e-02_wp
         BET102 = 6.2255521644e-02_wp
         BET012 = -2.6514181169e-03_wp
         BET003 = -2.3025968587e-04_wp
         !
         PEN000 = -9.8409626043_wp
         PEN100 = 2.1274999107e+01_wp
         PEN200 = -2.5387384109e+01_wp
         PEN300 = 1.5469038167e+01_wp
         PEN400 = -3.3025876549_wp
         PEN010 = 6.6681505563_wp
         PEN110 = 2.2435057288_wp
         PEN210 = -2.5021299030_wp
         PEN310 = 3.2699521832e-01_wp
         PEN020 = -3.3540239802_wp
         PEN120 = -1.7531540640_wp
         PEN220 = 9.3976864981e-01_wp
         PEN030 = 1.2324834767_wp
         PEN130 = 2.7538550639e-01_wp
         PEN040 = -2.7963967985e-01_wp
         PEN001 = -1.3773949450_wp
         PEN101 = 3.3018402659_wp
         PEN201 = -1.6679755496_wp
         PEN011 = -1.3709540999_wp
         PEN111 = 1.4207577012e-01_wp
         PEN021 = 8.2799886843e-01_wp
         PEN002 = 1.7507069098e-02_wp
         PEN102 = 1.3880727538e-02_wp
         PEN012 = -2.8477365341e-01_wp
         !
         APE000 = -1.6670376391e-01_wp
         APE100 = -5.6087643219e-02_wp
         APE200 = 6.2553247576e-02_wp
         APE300 = -8.1748804580e-03_wp
         APE010 = 1.6770119901e-01_wp
         APE110 = 8.7657703198e-02_wp
         APE210 = -4.6988432490e-02_wp
         APE020 = -9.2436260751e-02_wp
         APE120 = -2.0653912979e-02_wp
         APE030 = 2.7963967985e-02_wp
         APE001 = 3.4273852498e-02_wp
         APE101 = -3.5518942529e-03_wp
         APE011 = -4.1399943421e-02_wp
         APE002 = 7.1193413354e-03_wp
         !
         BPE000 = 2.6468936504e-01_wp
         BPE100 = -6.3170583896e-01_wp
         BPE200 = 5.7736640125e-01_wp
         BPE300 = -1.6435438140e-01_wp
         BPE010 = 2.7912203607e-02_wp
         BPE110 = -6.2259666565e-02_wp
         BPE210 = 1.2204769966e-02_wp
         BPE020 = -2.1811574876e-02_wp
         BPE120 = 2.3383950895e-02_wp
         BPE030 = 3.4261630030e-03_wp
         BPE001 = 4.1079296834e-02_wp
         BPE101 = -4.1503681096e-02_wp
         BPE011 = 1.7676120780e-03_wp
         BPE002 = 1.7269476440e-04_wp
         !
   ELSE                  ! EOS80 default equation
         PRINT *, ' Equation of state for sea-water : EOS-80 '
         rdeltaS = 20._wp
         r1_S0  = 1._wp/40._wp
         r1_T0  = 1._wp/40._wp
         r1_Z0  = 1.e-4_wp
         !
         EOS000 = 9.5356891948e+02_wp
         EOS100 = 1.7136499189e+02_wp
         EOS200 = -3.7501039454e+02_wp
         EOS300 = 5.1856810420e+02_wp
         EOS400 = -3.7264470465e+02_wp
         EOS500 = 1.4302533998e+02_wp
         EOS600 = -2.2856621162e+01_wp
         EOS010 = 1.0087518651e+01_wp
         EOS110 = -1.3647741861e+01_wp
         EOS210 = 8.8478359933_wp
         EOS310 = -7.2329388377_wp
         EOS410 = 1.4774410611_wp
         EOS510 = 2.0036720553e-01_wp
         EOS020 = -2.5579830599e+01_wp
         EOS120 = 2.4043512327e+01_wp
         EOS220 = -1.6807503990e+01_wp
         EOS320 = 8.3811577084_wp
         EOS420 = -1.9771060192_wp
         EOS030 = 1.6846451198e+01_wp
         EOS130 = -2.1482926901e+01_wp
         EOS230 = 1.0108954054e+01_wp
         EOS330 = -6.2675951440e-01_wp
         EOS040 = -8.0812310102_wp
         EOS140 = 1.0102374985e+01_wp
         EOS240 = -4.8340368631_wp
         EOS050 = 1.2079167803_wp
         EOS150 = 1.1515380987e-01_wp
         EOS060 = -2.4520288837e-01_wp
         EOS001 = 1.0748601068e+01_wp
         EOS101 = -1.7817043500e+01_wp
         EOS201 = 2.2181366768e+01_wp
         EOS301 = -1.6750916338e+01_wp
         EOS401 = 4.1202230403_wp
         EOS011 = -1.5852644587e+01_wp
         EOS111 = -7.6639383522e-01_wp
         EOS211 = 4.1144627302_wp
         EOS311 = -6.6955877448e-01_wp
         EOS021 = 9.9994861860_wp
         EOS121 = -1.9467067787e-01_wp
         EOS221 = -1.2177554330_wp
         EOS031 = -3.4866102017_wp
         EOS131 = 2.2229155620e-01_wp
         EOS041 = 5.9503008642e-01_wp
         EOS002 = 1.0375676547_wp
         EOS102 = -3.4249470629_wp
         EOS202 = 2.0542026429_wp
         EOS012 = 2.1836324814_wp
         EOS112 = -3.4453674320e-01_wp
         EOS022 = -1.2548163097_wp
         EOS003 = 1.8729078427e-02_wp
         EOS103 = -5.7238495240e-02_wp
         EOS013 = 3.8306136687e-01_wp
         !
         ALP000 = -2.5218796628e-01_wp
         ALP100 = 3.4119354654e-01_wp
         ALP200 = -2.2119589983e-01_wp
         ALP300 = 1.8082347094e-01_wp
         ALP400 = -3.6936026529e-02_wp
         ALP500 = -5.0091801383e-03_wp
         ALP010 = 1.2789915300_wp
         ALP110 = -1.2021756164_wp
         ALP210 = 8.4037519952e-01_wp
         ALP310 = -4.1905788542e-01_wp
         ALP410 = 9.8855300959e-02_wp
         ALP020 = -1.2634838399_wp
         ALP120 = 1.6112195176_wp
         ALP220 = -7.5817155402e-01_wp
         ALP320 = 4.7006963580e-02_wp
         ALP030 = 8.0812310102e-01_wp
         ALP130 = -1.0102374985_wp
         ALP230 = 4.8340368631e-01_wp
         ALP040 = -1.5098959754e-01_wp
         ALP140 = -1.4394226233e-02_wp
         ALP050 = 3.6780433255e-02_wp
         ALP001 = 3.9631611467e-01_wp
         ALP101 = 1.9159845880e-02_wp
         ALP201 = -1.0286156825e-01_wp
         ALP301 = 1.6738969362e-02_wp
         ALP011 = -4.9997430930e-01_wp
         ALP111 = 9.7335338937e-03_wp
         ALP211 = 6.0887771651e-02_wp
         ALP021 = 2.6149576513e-01_wp
         ALP121 = -1.6671866715e-02_wp
         ALP031 = -5.9503008642e-02_wp
         ALP002 = -5.4590812035e-02_wp
         ALP102 = 8.6134185799e-03_wp
         ALP012 = 6.2740815484e-02_wp
         ALP003 = -9.5765341718e-03_wp
         !
         BET000 = 2.1420623987_wp
         BET100 = -9.3752598635_wp
         BET200 = 1.9446303907e+01_wp
         BET300 = -1.8632235232e+01_wp
         BET400 = 8.9390837485_wp
         BET500 = -1.7142465871_wp
         BET010 = -1.7059677327e-01_wp
         BET110 = 2.2119589983e-01_wp
         BET210 = -2.7123520642e-01_wp
         BET310 = 7.3872053057e-02_wp
         BET410 = 1.2522950346e-02_wp
         BET020 = 3.0054390409e-01_wp
         BET120 = -4.2018759976e-01_wp
         BET220 = 3.1429341406e-01_wp
         BET320 = -9.8855300959e-02_wp
         BET030 = -2.6853658626e-01_wp
         BET130 = 2.5272385134e-01_wp
         BET230 = -2.3503481790e-02_wp
         BET040 = 1.2627968731e-01_wp
         BET140 = -1.2085092158e-01_wp
         BET050 = 1.4394226233e-03_wp
         BET001 = -2.2271304375e-01_wp
         BET101 = 5.5453416919e-01_wp
         BET201 = -6.2815936268e-01_wp
         BET301 = 2.0601115202e-01_wp
         BET011 = -9.5799229402e-03_wp
         BET111 = 1.0286156825e-01_wp
         BET211 = -2.5108454043e-02_wp
         BET021 = -2.4333834734e-03_wp
         BET121 = -3.0443885826e-02_wp
         BET031 = 2.7786444526e-03_wp
         BET002 = -4.2811838287e-02_wp
         BET102 = 5.1355066072e-02_wp
         BET012 = -4.3067092900e-03_wp
         BET003 = -7.1548119050e-04_wp
         !
         !
         PEN000 = -5.3743005340_wp
         PEN100 = 8.9085217499_wp
         PEN200 = -1.1090683384e+01_wp
         PEN300 = 8.3754581690_wp
         PEN400 = -2.0601115202_wp
         PEN010 = 7.9263222935_wp
         PEN110 = 3.8319691761e-01_wp
         PEN210 = -2.0572313651_wp
         PEN310 = 3.3477938724e-01_wp
         PEN020 = -4.9997430930_wp
         PEN120 = 9.7335338937e-02_wp
         PEN220 = 6.0887771651e-01_wp
         PEN030 = 1.7433051009_wp
         PEN130 = -1.1114577810e-01_wp
         PEN040 = -2.9751504321e-01_wp
         PEN001 = -6.9171176978e-01_wp
         PEN101 = 2.2832980419_wp
         PEN201 = -1.3694684286_wp
         PEN011 = -1.4557549876_wp
         PEN111 = 2.2969116213e-01_wp
         PEN021 = 8.3654420645e-01_wp
         PEN002 = -1.4046808820e-02_wp
         PEN102 = 4.2928871430e-02_wp
         PEN012 = -2.8729602515e-01_wp
         !
         APE000 = -1.9815805734e-01_wp
         APE100 = -9.5799229402e-03_wp
         APE200 = 5.1430784127e-02_wp
         APE300 = -8.3694846809e-03_wp
         APE010 = 2.4998715465e-01_wp
         APE110 = -4.8667669469e-03_wp
         APE210 = -3.0443885826e-02_wp
         APE020 = -1.3074788257e-01_wp
         APE120 = 8.3359333577e-03_wp
         APE030 = 2.9751504321e-02_wp
         APE001 = 3.6393874690e-02_wp
         APE101 = -5.7422790533e-03_wp
         APE011 = -4.1827210323e-02_wp
         APE002 = 7.1824006288e-03_wp
         !
         BPE000 = 1.1135652187e-01_wp
         BPE100 = -2.7726708459e-01_wp
         BPE200 = 3.1407968134e-01_wp
         BPE300 = -1.0300557601e-01_wp
         BPE010 = 4.7899614701e-03_wp
         BPE110 = -5.1430784127e-02_wp
         BPE210 = 1.2554227021e-02_wp
         BPE020 = 1.2166917367e-03_wp
         BPE120 = 1.5221942913e-02_wp
         BPE030 = -1.3893222263e-03_wp
         BPE001 = 2.8541225524e-02_wp
         BPE101 = -3.4236710714e-02_wp
         BPE011 = 2.8711395266e-03_wp
         BPE002 = 5.3661089288e-04_wp
         !
   ENDIF

  END SUBROUTINE eos_init


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

    sigma0(:,:) = sigmai (ptem,psal, 0.0, kpi, kpj )

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

  FUNCTION spice_80 ( ptem, psal, kpi, kpj )
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

    REAL(KIND=8), DIMENSION(kpi,kpj)             :: spice_80   ! return value

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

    spice_80(:,:) = dl_spi( :,:)

  END FUNCTION spice_80

  FUNCTION spice_pol ( ptem, psal, kpi, kpj, ld_teos10 )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION spice  ***
    !!
    !! ** Purpose :  Compute spiciness from T S fields. 
    !!
    !! ** Method  :  use GSW formula from GSW toolbox
    !!
    !!  **  Example:
    !!       spice(15,33)=   0.5445863      0.544586321373410  calcul en double
    !!       spice(15,33)=   0.5445864      (calcul en simple precision) 
    !!
    !!  ** References : Flament (2002) "A state variable for characterizing
    !!              water masses and their diffusive stability: spiciness."
    !!              Progress in Oceanography Volume 54, 2002, Pages 493-501.
    !!              gsw toolbox
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptem, psal ! temperature salinity
    INTEGER(KIND=4),                  INTENT(in) :: kpi,kpj    ! dimension of 2D arrays
    LOGICAL,                          INTENT(in) :: ld_teos10  ! flag for teos10 formula

    REAL(KIND=8), DIMENSION(kpi,kpj)             :: spice_pol  ! return value

    INTEGER(KIND=4)                        :: ji,jj     ! dummy loop index
    REAL(KIND=8) :: xs, ys, offset, sfac
    !!----------------------------------------------------------------------
    sfac= 0.0248826675584615_8
    offset=5.971840214030754e-1_8
    DO jj = 1, kpj
       DO ji=1,kpi
         xs = sqrt(sfac*psal(ji,jj) + offset)
         ys = ptem(ji,jj)*0.025

         spice_pol(ji,jj) =  &
          s01 + ys*(s02 + ys*(s03 + ys*(s04 + ys*(s05 + ys*(s06 + s07*ys))))) &
    + xs*(s08 + ys*(s09 + ys*(s10 + ys*(s11 + ys*(s12 + ys*(s13 + s14*ys))))) &
    + xs*(s15 + ys*(s16 + ys*(s17 + ys*(s18 + ys*(s19 + ys*(s20 + s21*ys))))) &
    + xs*(s22 + ys*(s23 + ys*(s24 + ys*(s25 + ys*(s26 + ys*(s27 + s28*ys))))) &
    + xs*(s29 + ys*(s30 + ys*(s31 + ys*(s32 + ys*(s33 + ys*(s34 + s35*ys))))) &
    + xs*(s36 + ys*(s37 + ys*(s38 + ys*(s39 + ys*(s40 + ys*(s41 + s42*ys))))) &
    + xs*(s43 + ys*(s44 + ys*(s45 + ys*(s46 + ys*(s47 + ys*(s48 + s49*ys)))))))))))
       ENDDO
    ENDDO

  END FUNCTION spice_pol

  FUNCTION sigmai_dep ( ptem, psal, pref, kpi, kpj )
    !! --------------------------------------------------------------------
    !! ** Purpose :   Compute the  density referenced to pref (ratio rho/rau0) 
    !!       from temperature and salinity fields using an equation of state 
    !!       defined through the ll_teos10 flag defined in eos_init
    !!
    !! ** Method  :   prd(t,s,z) = ( rho(t,s,z) - rau0 ) / rau0
    !!         with   prd    in situ density anomaly      no units
    !!                t      TEOS10: CT or EOS80: PT      Celsius
    !!                s      TEOS10: SA or EOS80: SP      TEOS10: g/kg or EOS80: psu
    !!                z      depth                        meters
    !!                rho    in situ density              kg/m^3
    !!                rau0   reference density            kg/m^3
    !!
    !!     ln_teos10 : polynomial TEOS-10 equation of state is used for rho(t,s,z).
    !!         Check value: rho = 1028.21993233072 kg/m^3 for z=3000 dbar, ct=3 Celsius, sa=35.5 g/kg
    !!
    !!     ln_eos80 : polynomial EOS-80 equation of state is used for rho(t,s,z).
    !!         Check value: rho = 1028.35011066567 kg/m^3 for z=3000 dbar, pt=3 Celsius, sp=35.5 psu
    !!
    !!     ln_seos : simplified equation of state
    !!              prd(t,s,z) = ( -a0*(1+lambda/2*(T-T0)+mu*z+nu*(S-S0))*(T-T0) + b0*(S-S0) ) / rau0
    !!              linear case function of T only: rn_alpha<>0, other coefficients = 0
    !!              linear eos function of T and S: rn_alpha and rn_beta<>0, other coefficients=0
    !!              Vallis like equation: use default values of coefficients
    !!
    !! ** Action  :   compute prd , the in situ density (no units)
    !!
    !! References :   Roquet et al, Ocean Modelling, in preparation (2014)
    !!                Vallis, Atmospheric and Oceanic Fluid Dynamics, 2006
    !!                TEOS-10 Manual, 2010
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptem, psal ! temperature salinity
    INTEGER(KIND=4),                  INTENT(in) :: kpi,kpj    ! dimension of 2D arrays
    REAL(KIND=4),                     INTENT(in) :: pref       ! reference pressure (meters or db)
    REAL(KIND=8), DIMENSION(kpi,kpj) :: sigmai_dep             ! return value

    INTEGER(KIND=4)  :: ji, jj
    REAL(KIND=8)     :: dlt, dlh,dls,dltm, dlr, dlr0, dlr1, dlr2, dlr3, dlref
    !!----------------------------------------------------------------------
            dlh = pref *r1_Z0
            sigmai_dep = 0.d0
!  Reference density at pref
            dlref = ( ( ( ( ( R05*dlh + R04 )*dlh + R03 )*dlh + R02 )*dlh + R01 )*dlh + R00 )*dlh
!  Compute anomaly with Polynomial_bsq formula
            DO jj = 1, kpj
               DO ji = 1, kpi
                  !
                  dlt  = ptem (ji,jj) * r1_T0                           ! temperature
                  dls  = SQRT( ABS( psal(ji,jj) + rdeltaS ) * r1_S0 )   ! square root salinity
                  IF ( psal(ji,jj) == 0 ) THEN
                     dltm=0.d0
                  ELSE
                     dltm=1.d0
                  ENDIF
                  !
                  dlr3 = EOS013*dlt   &
                     &   + EOS103*dls+EOS003
                     !
                  dlr2 = (EOS022*dlt   &
                     &   + EOS112*dls+EOS012)*dlt   &
                     &   + (EOS202*dls+EOS102)*dls+EOS002
                     !
                  dlr1 = (((EOS041*dlt   &
                     &   + EOS131*dls+EOS031)*dlt   &
                     &   + (EOS221*dls+EOS121)*dls+EOS021)*dlt   &
                     &   + ((EOS311*dls+EOS211)*dls+EOS111)*dls+EOS011)*dlt   &
                     &   + (((EOS401*dls+EOS301)*dls+EOS201)*dls+EOS101)*dls+EOS001
                     !
                  dlr0 = (((((EOS060*dlt   &
                     &   + EOS150*dls+EOS050)*dlt   &
                     &   + (EOS240*dls+EOS140)*dls+EOS040)*dlt   &
                     &   + ((EOS330*dls+EOS230)*dls+EOS130)*dls+EOS030)*dlt   &
                     &   + (((EOS420*dls+EOS320)*dls+EOS220)*dls+EOS120)*dls+EOS020)*dlt   &
                     &   + ((((EOS510*dls+EOS410)*dls+EOS310)*dls+EOS210)*dls+EOS110)*dls+EOS010)*dlt   &
                     &   + (((((EOS600*dls+EOS500)*dls+EOS400)*dls+EOS300)*dls+EOS200)*dls+EOS100)*dls+EOS000
                     !
                  dlr  = ( ( dlr3 * dlh + dlr2 ) * dlh + dlr1 ) * dlh+ dlr0
                  !
!                  sigmai_dep(ji,jj) = (  dlr * r1_rau0 - 1._wp  ) * dltm  ! density anomaly (masked)
                  sigmai_dep(ji,jj) = (  dlr +dlref - rau0  ) * dltm  ! density anomaly (masked)
                  !
               END DO
            END DO

  END FUNCTION sigmai_dep

  FUNCTION sigmai_dep_80 ( ptem, psal, pref, kpi,kpj,ld_eos80)
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
    LOGICAL,                          INTENT(in) :: ld_eos80   ! flog use to differentiate eos80 style in the INTERFACE
    REAL(KIND=8), DIMENSION(kpi,kpj) :: sigmai_dep_80          ! return value

    REAL(kind=8), PARAMETER :: dpr4=4.8314d-4, dpd=-2.042967d-2 , dprau0 = 1000.d0

    INTEGER(KIND=4) :: ji, jj 
    REAL(KIND=8), DIMENSION (kpi,kpj) ::  dlrs
    REAL(KIND=8) :: dlt, dls      
    REAL(KIND=8) :: dla, dla1, dlaw, dlb, dlb1, dlbw, dlc, dle, dlk0, dlkw 
    REAL(kind=8) :: dlrhop, dlr1, dlr2, dlr3, dlref
    !! --------------------------------------------------------------------
    dlref      = pref
    sigmai_dep_80 = 0.d0
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
          sigmai_dep_80(ji,jj)=dlrhop/(1.0d0-dlref/(dlk0-dlref*(dla-dlref*dlb)))-dprau0

       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END FUNCTION sigmai_dep_80

  FUNCTION sigmai_dep2d ( ptem, psal, pref, kpi, kpj )
    !! --------------------------------------------------------------------
    !! ** Purpose :   Compute the  density referenced to pref (ratio rho/rau0) 
    !!       from temperature and salinity fields using an equation of state 
    !!       defined through the ll_teos10 flag defined in eos_init
    !!
    !! ** Method  :   prd(t,s,z) = ( rho(t,s,z) - rau0 ) / rau0
    !!         with   prd    in situ density anomaly      no units
    !!                t      TEOS10: CT or EOS80: PT      Celsius
    !!                s      TEOS10: SA or EOS80: SP      TEOS10: g/kg or EOS80: psu
    !!                z      depth                        meters
    !!                rho    in situ density              kg/m^3
    !!                rau0   reference density            kg/m^3
    !!
    !!     ln_teos10 : polynomial TEOS-10 equation of state is used for rho(t,s,z).
    !!         Check value: rho = 1028.21993233072 kg/m^3 for z=3000 dbar, ct=3 Celsius, sa=35.5 g/kg
    !!
    !!     ln_eos80 : polynomial EOS-80 equation of state is used for rho(t,s,z).
    !!         Check value: rho = 1028.35011066567 kg/m^3 for z=3000 dbar, pt=3 Celsius, sp=35.5 psu
    !!
    !!     ln_seos : simplified equation of state
    !!              prd(t,s,z) = ( -a0*(1+lambda/2*(T-T0)+mu*z+nu*(S-S0))*(T-T0) + b0*(S-S0) ) / rau0
    !!              linear case function of T only: rn_alpha<>0, other coefficients = 0
    !!              linear eos function of T and S: rn_alpha and rn_beta<>0, other coefficients=0
    !!              Vallis like equation: use default values of coefficients
    !!
    !! ** Action  :   compute prd , the in situ density (no units)
    !!
    !! References :   Roquet et al, Ocean Modelling, in preparation (2014)
    !!                Vallis, Atmospheric and Oceanic Fluid Dynamics, 2006
    !!                TEOS-10 Manual, 2010
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptem, psal ! temperature salinity
    INTEGER(KIND=4),                  INTENT(in) :: kpi,kpj    ! dimension of 2D arrays
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: pref       ! reference pressure (meters or db)
    REAL(KIND=8), DIMENSION(kpi,kpj) :: sigmai_dep2d           ! return value

    INTEGER(KIND=4)  :: ji, jj
    REAL(KIND=8)     :: dlt, dlh,dls,dltm, dlr, dlr0, dlr1, dlr2, dlr3, dlref
    !!----------------------------------------------------------------------
            sigmai_dep2d = 0.d0
            DO jj = 1, kpj
               DO ji = 1, kpi
                  !
                  dlh  = pref(ji,jj) * r1_Z0                                  ! depth
                  dlt  = ptem (ji,jj) * r1_T0                           ! temperature
                  dls  = SQRT( ABS( psal(ji,jj) + rdeltaS ) * r1_S0 )   ! square root salinity
                  !  Reference density at pref
                  dlref = ( ( ( ( ( R05*dlh + R04 )*dlh + R03 )*dlh + R02 )*dlh + R01 )*dlh + R00 )*dlh
                  IF ( psal(ji,jj) == 0 ) THEN
                     dltm=0.d0
                  ELSE
                     dltm=1.d0
                  ENDIF
                  !
                  dlr3 = EOS013*dlt   &
                     &   + EOS103*dls+EOS003
                     !
                  dlr2 = (EOS022*dlt   &
                     &   + EOS112*dls+EOS012)*dlt   &
                     &   + (EOS202*dls+EOS102)*dls+EOS002
                     !
                  dlr1 = (((EOS041*dlt   &
                     &   + EOS131*dls+EOS031)*dlt   &
                     &   + (EOS221*dls+EOS121)*dls+EOS021)*dlt   &
                     &   + ((EOS311*dls+EOS211)*dls+EOS111)*dls+EOS011)*dlt   &
                     &   + (((EOS401*dls+EOS301)*dls+EOS201)*dls+EOS101)*dls+EOS001
                     !
                  dlr0 = (((((EOS060*dlt   &
                     &   + EOS150*dls+EOS050)*dlt   &
                     &   + (EOS240*dls+EOS140)*dls+EOS040)*dlt   &
                     &   + ((EOS330*dls+EOS230)*dls+EOS130)*dls+EOS030)*dlt   &
                     &   + (((EOS420*dls+EOS320)*dls+EOS220)*dls+EOS120)*dls+EOS020)*dlt   &
                     &   + ((((EOS510*dls+EOS410)*dls+EOS310)*dls+EOS210)*dls+EOS110)*dls+EOS010)*dlt   &
                     &   + (((((EOS600*dls+EOS500)*dls+EOS400)*dls+EOS300)*dls+EOS200)*dls+EOS100)*dls+EOS000
                     !
                  dlr  = ( ( dlr3 * dlh + dlr2 ) * dlh + dlr1 ) * dlh + dlr0
                  !
!                 sigmai_dep2d(ji,jj) = (  dlr * r1_rau0 - 1._wp  ) * dltm  ! density anomaly (masked)
                  sigmai_dep2d(ji,jj) = (  dlr +dlref - rau0  ) * dltm  ! density anomaly (masked)
                  !
               END DO
            END DO


  END FUNCTION sigmai_dep2d

  FUNCTION sigmai_dep2d_80 ( ptem, psal, pref, kpi,kpj, ld_eos80)
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
    LOGICAL,                          INTENT(in) :: ld_eos80   ! flag used to discriminate original EOS80 for INTERFACE
    REAL(KIND=8), DIMENSION(kpi,kpj) :: sigmai_dep2d_80           ! return value

    REAL(kind=8), PARAMETER :: dpr4=4.8314d-4, dpd=-2.042967d-2 , dprau0 = 1000.d0

    INTEGER(KIND=4) :: ji, jj 
    REAL(KIND=8), DIMENSION (kpi,kpj) ::  dlrs
    REAL(KIND=8) :: dlt, dls      
    REAL(KIND=8) :: dla, dla1, dlaw, dlb, dlb1, dlbw, dlc, dle, dlk0, dlkw 
    REAL(kind=8) :: dlrhop, dlr1, dlr2, dlr3, dlref

    sigmai_dep2d_80 = 0.d0
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
          sigmai_dep2d_80(ji,jj)=dlrhop/(1.0d0-dlref/(dlk0-dlref*(dla-dlref*dlb)))&
               -dprau0

       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END FUNCTION sigmai_dep2d_80

  FUNCTION eosbn2_pol ( ptem, psal, pdep, pe3w, kpi, kpj, kup, kdown)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION eosbn2  ***
    !!
    !! ** Purpose :  Compute the local Brunt-Vaisala frequency at the time-
    !!               step of the input arguments
    !!
    !! ** Method :  The brunt-vaisala frequency is computed using the
    !!              polynomial expression of Roquet_et_al (2015) 
    !!              N^2 = grav/Rho_ref * ( alpha*dk[ t ] + beta*dk[ s ] )/e3w
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj,2), INTENT(in) :: ptem, psal ! temperaature salinity
    REAL(KIND=4)                                   :: pdep       ! reference depth
    REAL(KIND=4), DIMENSION(kpi,kpj),   INTENT(in) :: pe3w       ! e3w of the current layer
    INTEGER(KIND=4),                    INTENT(in) :: kpi, kpj   ! size of the array
    INTEGER(KIND=4),                    INTENT(in) :: kup, kdown ! index of levels up and down

    REAL(KIND=4), DIMENSION(kpi,kpj)               :: eosbn2_pol ! returned values

    INTEGER(KIND=4) :: ji, jj         ! dummy loop indices
    REAL(KIND=8)    :: dlgde3w, dlt, dls, dlh
    REAL(KIND=8)    :: dlalpha, dlbeta
    REAL(KIND=8)    :: dlgrav=9.81
    !!----------------------------------------------------------------------

    dlh = pdep
    !$OMP PARALLEL DO SCHEDULE(RUNTIME)
    DO jj = 1, kpj
       DO ji = 1, kpi
          dlgde3w = dlgrav / pe3w(ji,jj) *r1_rau0                      ! division by rho is done here
          dlt = 0.5 * ( ptem(ji,jj,kup) + ptem(ji,jj,kdown) )          ! potential temperature at w-point
          dls = 0.5 * ( psal(ji,jj,kup) + psal(ji,jj,kdown) ) - 35.0   ! salinity anomaly (s-35) at w-point
          dlalpha = alpha (dlt, dls, dlh )
          dlbeta  = beta  (dlt, dls, dlh )
          eosbn2_pol(ji,jj) = dlgde3w  *                                           &   ! N^2
               &            ( dlalpha * ( ptem(ji,jj,kup) - ptem(ji,jj,kdown) )   &
               &             -dlbeta  * ( psal(ji,jj,kup) - psal(ji,jj,kdown) ) )
       END DO
    END DO
    !$OMP END PARALLEL DO

  END FUNCTION eosbn2_pol

  FUNCTION eosbn2_80 ( ptem, psal, pdep, pe3w, kpi, kpj, kup, kdown, ld_eos80)
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
    LOGICAL,                            INTENT(in) :: ld_eos80   ! flag for differenciation eos80
    REAL(KIND=4), DIMENSION(kpi,kpj)               :: eosbn2_80  ! returned values

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

          eosbn2_80(ji,jj) = zgde3w * zbeta                                         &   ! N^2
               &          * ( zalbet * ( ptem(ji,jj,kup) - ptem(ji,jj,kdown) )   &
               &                     - ( psal(ji,jj,kup) - psal(ji,jj,kdown) ) )
       END DO
    END DO
    !$OMP END PARALLEL DO

  END FUNCTION eosbn2_80


  FUNCTION albet_2d(  ptem, psal, pdep, kpi, kpj)
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

    REAL(KIND=8), DIMENSION(kpi,kpj)             :: albet_2d   ! returned value

    INTEGER(KIND=4) :: ji, jj      ! dummy loop index
    REAL(KIND=8)    :: zt, zs, zh  ! working local variables
    !!----------------------------------------------------------------------
    zh = pdep
    !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     DO jj=1,kpj
       DO ji=1,kpi
          zt =  ptem(ji,jj)         ! potential temperature
          zs =  psal(ji,jj)- 35.0   ! salinity anomaly (s-35)

          albet_2d(ji,jj) = ( ( ( - 0.255019e-07 * zt + 0.298357e-05 ) * zt   &   ! ratio alpha/beta
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

  END FUNCTION albet_2d


  FUNCTION albet_z(  ptem, psal, pdep, kpk)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION  albet  ***
    !!
    !! ** Purpose :  Compute the ratio alpha/beta 
    !!
    !! ** Method  :  Follow Mc Dougal et al as in other functions
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpk), INTENT(in) :: ptem, psal ! temperature salinity
    REAL(KIND=4), DIMENSION(kpk), INTENT(in) :: pdep       ! refererence depth
    INTEGER(KIND=4),              INTENT(in) :: kpk        ! size of the arrays

    REAL(KIND=8), DIMENSION(kpk)             :: albet_z    ! returned value

    INTEGER(KIND=4) :: jk          ! dummy loop index
    REAL(KIND=8)    :: zt, zs, zh  ! working local variables
    !!----------------------------------------------------------------------
    !$OMP PARALLEL DO SCHEDULE(RUNTIME)
     DO jk=1,kpk
          zh = pdep(jk)
          zt = ptem(jk)         ! potential temperature
          zs = psal(jk)- 35.0   ! salinity anomaly (s-35)

          albet_z(jk) = (  ( ( - 0.255019e-07 * zt + 0.298357e-05 ) * zt   &   ! ratio alpha/beta
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
    !$OMP END PARALLEL DO

  END FUNCTION albet_z

  FUNCTION beta_2d (  ptem, psal, pdep, kpi, kpj)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION beta  ***
    !!
    !! ** Purpose :  Compute the beta 
    !!
    !! ** Method  :  Follow Roquet et Al OM 2015
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj),INTENT(in) :: ptem, psal ! temperature salinity
    REAL(KIND=4),                    INTENT(in) :: pdep       ! reference depth
    INTEGER(KIND=4),                 INTENT(in) :: kpi, kpj   ! size of the array

    REAL(KIND=8), DIMENSION(kpi,kpj)            :: beta_2d    ! returned values

    INTEGER(KIND=4) :: ji, jj      ! dummy loop index
    REAL(KIND=8) :: dln3, dln2, dln1, dln0, dln, dltm
    REAL(KIND=8) :: dls, dlt, dlh
    !!----------------------------------------------------------------------
    dlh  = pdep * r1_Z0                                ! depth
    DO jj = 1, kpj
      DO ji = 1, kpi
         dlt  = ptem (ji,jj) * r1_T0                           ! temperature
         dls  = SQRT( ABS( psal(ji,jj) + rdeltaS ) * r1_S0 )   ! square root salinity
         IF ( psal(ji,jj) == 0 ) THEN
              dltm=0.d0
         ELSE
              dltm=1.d0
         ENDIF
         dln3 = BET003
         !
         dln2 = BET012*dlt + BET102*dls+BET002
         !
         dln1 = ((BET031*dlt   &
            &   + BET121*dls+BET021)*dlt   &
            &   + (BET211*dls+BET111)*dls+BET011)*dlt   &
            &   + ((BET301*dls+BET201)*dls+BET101)*dls+BET001
            !
         dln0 = ((((BET050*dlt   &
            &   + BET140*dls+BET040)*dlt   &
            &   + (BET230*dls+BET130)*dls+BET030)*dlt   &
            &   + ((BET320*dls+BET220)*dls+BET120)*dls+BET020)*dlt   &
            &   + (((BET410*dls+BET310)*dls+BET210)*dls+BET110)*dls+BET010)*dlt   &
            &   + ((((BET500*dls+BET400)*dls+BET300)*dls+BET200)*dls+BET100)*dls+BET000
            !
         dln  = ( ( dln3 * dlh + dln2 ) * dlh + dln1 ) * dlh + dln0
         !
         beta_2d(ji,jj)  = dln / dls * dltm
       ENDDO
    ENDDO
  END FUNCTION beta_2d

  FUNCTION alpha_2d (  ptem, psal, pdep, kpi, kpj)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION beta  ***
    !!
    !! ** Purpose :  Compute alpha, the thermal expansion coeffcient
    !!
    !! ** Method  :  Follow Roquet et Al OM 2015
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj),INTENT(in) :: ptem, psal ! temperature salinity
    REAL(KIND=4),                    INTENT(in) :: pdep       ! reference depth
    INTEGER(KIND=4),                 INTENT(in) :: kpi, kpj   ! size of the array


    REAL(KIND=8), DIMENSION(kpi,kpj)            :: alpha_2d   ! returned values
    
    INTEGER(KIND=4)   :: ji, jj               ! dummy loop index
    REAL(KIND=8)      :: dlh, dlt, dls, dltm
    REAL(KIND=8)      :: dln3, dln2, dln1, dln0, dln 
    !!---------------------------------------------------------------------
    dlh  = pdep * r1_Z0                                ! depth
    DO jj = 1, kpj
      DO ji = 1, kpi
         dlt  = ptem (ji,jj) * r1_T0                           ! temperature
         dls  = SQRT( ABS( psal(ji,jj) + rdeltaS ) * r1_S0 )   ! square root salinity
         IF ( psal(ji,jj) == 0 ) THEN
              dltm=0.d0
         ELSE
              dltm=1.d0
         ENDIF
         dln3 = ALP003
         !
         dln2 = ALP012*dlt + ALP102*dls+ALP002
         !
         dln1 = ((ALP031*dlt   &
             &   + ALP121*dls+ALP021)*dlt   &
             &   + (ALP211*dls+ALP111)*dls+ALP011)*dlt   &
             &   + ((ALP301*dls+ALP201)*dls+ALP101)*dls+ALP001
         !
         dln0 = ((((ALP050*dlt   &
             &   +  ALP140*dls+ALP040)*dlt   &
             &   + (ALP230*dls+ALP130)*dls+ALP030)*dlt   &
             &   + ((ALP320*dls+ALP220)*dls+ALP120)*dls+ALP020)*dlt   &
             &   + (((ALP410*dls+ALP310)*dls+ALP210)*dls+ALP110)*dls+ALP010)*dlt   &
             &   + ((((ALP500*dls+ALP400)*dls+ALP300)*dls+ALP200)*dls+ALP100)*dls+ALP000
         !
         dln  = ( ( dln3 * dlh + dln2 ) * dlh + dln1 ) * dlh + dln0

        alpha_2d(ji,jj)=dln * dltm
      ENDDO
   ENDDO
  END FUNCTION alpha_2d

  FUNCTION alpha_0d  (  ptem, psal, pdep )
    REAL(KIND=8),                    INTENT(in) :: ptem, psal ! temperature salinity
    REAL(KIND=8),                    INTENT(in) :: pdep       ! reference depth

    REAL(KIND=8)                                :: alpha_0d   ! returned values
    
    INTEGER(KIND=4)   :: ji, jj               ! dummy loop index
    REAL(KIND=8)      :: dlh, dlt, dls, dltm
    REAL(KIND=8)      :: dln3, dln2, dln1, dln0, dln 
    !!---------------------------------------------------------------------
    dlh  = pdep * r1_Z0                                ! depth
         dlt  = ptem  * r1_T0                           ! temperature
         dls  = SQRT( ABS( psal + rdeltaS ) * r1_S0 )   ! square root salinity
         IF ( psal == 0 ) THEN
              dltm=0.d0
         ELSE
              dltm=1.d0
         ENDIF
         dln3 = ALP003
         !
         dln2 = ALP012*dlt + ALP102*dls+ALP002
         !
         dln1 = ((ALP031*dlt   &
             &   + ALP121*dls+ALP021)*dlt   &
             &   + (ALP211*dls+ALP111)*dls+ALP011)*dlt   &
             &   + ((ALP301*dls+ALP201)*dls+ALP101)*dls+ALP001
         !
         dln0 = ((((ALP050*dlt   &
             &   +  ALP140*dls+ALP040)*dlt   &
             &   + (ALP230*dls+ALP130)*dls+ALP030)*dlt   &
             &   + ((ALP320*dls+ALP220)*dls+ALP120)*dls+ALP020)*dlt   &
             &   + (((ALP410*dls+ALP310)*dls+ALP210)*dls+ALP110)*dls+ALP010)*dlt   &
             &   + ((((ALP500*dls+ALP400)*dls+ALP300)*dls+ALP200)*dls+ALP100)*dls+ALP000
         !
         dln  = ( ( dln3 * dlh + dln2 ) * dlh + dln1 ) * dlh + dln0

        alpha_0d=dln * dltm

  END FUNCTION alpha_0d

  FUNCTION alpha_z (ptem, psal, pdep, kpk )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION alpha_z  ***
    !!
    !! ** Purpose :  Compute the beta 
    !!
    !! ** Method  :  Follow Roquet et al. OM 2015
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpk),INTENT(in) :: ptem, psal ! temperature salinity
    REAL(KIND=4), DIMENSION(kpk),INTENT(in) :: pdep       ! reference depth
    INTEGER(KIND=4),             INTENT(in) :: kpk        ! size of the array

    REAL(KIND=8), DIMENSION(kpk)            :: alpha_z     ! returned values

    INTEGER(KIND=4) :: jk          ! dummy loop index
    REAL(KIND=8)    :: dlh, dlt, dls, dltm
    REAL(KIND=8)    :: dln3, dln2, dln1, dln0, dln
    !!----------------------------------------------------------------------
    DO jk = 1, kpk
         dlh  = pdep(jk) * r1_Z0                                ! depth
         dlt  = ptem (jk) * r1_T0                           ! temperature
         dls  = SQRT( ABS( psal(jk) + rdeltaS ) * r1_S0 )   ! square root salinity
         IF ( psal(jk) == 0 ) THEN
              dltm=0.d0
         ELSE
              dltm=1.d0
         ENDIF
         dln3 = ALP003
         !
         dln2 = ALP012*dlt + ALP102*dls+ALP002
         !
         dln1 = ((ALP031*dlt   &
             &   + ALP121*dls+ALP021)*dlt   &
             &   + (ALP211*dls+ALP111)*dls+ALP011)*dlt   &
             &   + ((ALP301*dls+ALP201)*dls+ALP101)*dls+ALP001
         !
         dln0 = ((((ALP050*dlt   &
             &   +  ALP140*dls+ALP040)*dlt   &
             &   + (ALP230*dls+ALP130)*dls+ALP030)*dlt   &
             &   + ((ALP320*dls+ALP220)*dls+ALP120)*dls+ALP020)*dlt   &
             &   + (((ALP410*dls+ALP310)*dls+ALP210)*dls+ALP110)*dls+ALP010)*dlt   &
             &   + ((((ALP500*dls+ALP400)*dls+ALP300)*dls+ALP200)*dls+ALP100)*dls+ALP000
         !
         dln  = ( ( dln3 * dlh + dln2 ) * dlh + dln1 ) * dlh + dln0

        alpha_z(jk)=dln * dltm
   ENDDO

    END FUNCTION alpha_z

  FUNCTION beta_0d (  ptem, psal, pdep )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION beta  ***
    !!
    !! ** Purpose :  Compute the beta 
    !!
    !! ** Method  :  Follow Roquet et Al OM 2015
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=8),                    INTENT(in) :: ptem, psal ! temperature salinity
    REAL(KIND=8),                    INTENT(in) :: pdep       ! reference depth

    REAL(KIND=8)                                :: beta_0d    ! returned values

    REAL(KIND=8) :: dln3, dln2, dln1, dln0, dln, dltm
    REAL(KIND=8) :: dls, dlt, dlh
    !!----------------------------------------------------------------------
    dlh  = pdep * r1_Z0                                ! depth
         dlt  = ptem  * r1_T0                           ! temperature
         dls  = SQRT( ABS( psal + rdeltaS ) * r1_S0 )   ! square root salinity
         IF ( psal == 0 ) THEN
              dltm=0.d0
         ELSE
              dltm=1.d0
         ENDIF
         dln3 = BET003
         !
         dln2 = BET012*dlt + BET102*dls+BET002
         !
         dln1 = ((BET031*dlt   &
            &   + BET121*dls+BET021)*dlt   &
            &   + (BET211*dls+BET111)*dls+BET011)*dlt   &
            &   + ((BET301*dls+BET201)*dls+BET101)*dls+BET001
            !
         dln0 = ((((BET050*dlt   &
            &   + BET140*dls+BET040)*dlt   &
            &   + (BET230*dls+BET130)*dls+BET030)*dlt   &
            &   + ((BET320*dls+BET220)*dls+BET120)*dls+BET020)*dlt   &
            &   + (((BET410*dls+BET310)*dls+BET210)*dls+BET110)*dls+BET010)*dlt   &
            &   + ((((BET500*dls+BET400)*dls+BET300)*dls+BET200)*dls+BET100)*dls+BET000
            !
         dln  = ( ( dln3 * dlh + dln2 ) * dlh + dln1 ) * dlh + dln0
         !
         beta_0d = dln / dls * dltm

  END FUNCTION beta_0d

  FUNCTION beta_z (  ptem, psal, pdep, kpk)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION beta  ***
    !!
    !! ** Purpose :  Compute the beta 
    !!
    !! ** Method  :  Follow Roquet et al. OM 2015
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpk),INTENT(in) :: ptem, psal ! temperature salinity
    REAL(KIND=4), DIMENSION(kpk),INTENT(in) :: pdep       ! reference depth
    INTEGER(KIND=4),             INTENT(in) :: kpk        ! size of the array

    REAL(KIND=8), DIMENSION(kpk)            :: beta_z     ! returned values

    INTEGER(KIND=4) :: jk          ! dummy loop index
    REAL(KIND=8)    :: dlh, dlt, dls, dltm
    REAL(KIND=8)    :: dln3, dln2, dln1, dln0, dln
    !!----------------------------------------------------------------------
    DO jk = 1, kpk
         dlh  = pdep(jk) * r1_Z0                                ! depth
         dlt  = ptem (jk) * r1_T0                           ! temperature
         dls  = SQRT( ABS( psal(jk) + rdeltaS ) * r1_S0 )   ! square root salinity
         IF ( psal(jk) == 0 ) THEN
              dltm=0.d0
         ELSE
              dltm=1.d0
         ENDIF
         dln3 = BET003
         !
         dln2 = BET012*dlt + BET102*dls+BET002
         !
         dln1 = ((BET031*dlt   &
            &   + BET121*dls+BET021)*dlt   &
            &   + (BET211*dls+BET111)*dls+BET011)*dlt   &
            &   + ((BET301*dls+BET201)*dls+BET101)*dls+BET001
            !
         dln0 = ((((BET050*dlt   &
            &   + BET140*dls+BET040)*dlt   &
            &   + (BET230*dls+BET130)*dls+BET030)*dlt   &
            &   + ((BET320*dls+BET220)*dls+BET120)*dls+BET020)*dlt   &
            &   + (((BET410*dls+BET310)*dls+BET210)*dls+BET110)*dls+BET010)*dlt   &
            &   + ((((BET500*dls+BET400)*dls+BET300)*dls+BET200)*dls+BET100)*dls+BET000
            !
         dln  = ( ( dln3 * dlh + dln2 ) * dlh + dln1 ) * dlh + dln0
         !
         beta_z(jk)  = dln / dls * dltm
    ENDDO
    END FUNCTION beta_z

  FUNCTION beta_2d_80 (  ptem, psal, pdep, kpi, kpj, ld_eos80)
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
    LOGICAL,                         INTENT(in) :: ld_eos80   ! flag for EOS80 discrimination

    REAL(KIND=8), DIMENSION(kpi,kpj)            :: beta_2d_80    ! returned values

    INTEGER(KIND=4) :: ji, jj      ! dummy loop index
    REAL(KIND=8)    :: zt, zs, zh  ! working variables
    !!----------------------------------------------------------------------
    zh = pdep
    !$OMP PARALLEL DO SCHEDULE(RUNTIME)
    DO jj=1,kpj
       DO ji=1,kpi
          zt =  ptem(ji,jj)         ! potential temperature
          zs =  psal(ji,jj)- 35.0   ! salinity anomaly (s-35)

          beta_2d_80(ji,jj) = ( ( -0.415613e-09 * zt + 0.555579e-07 ) * zt      &   ! beta
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

  END FUNCTION beta_2d_80

  FUNCTION beta_z_80 (  ptem, psal, pdep, kpk, ld_eos80)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION beta  ***
    !!
    !! ** Purpose :  Compute the beta 
    !!
    !! ** Method  :  Follow Mc Dougal et al as in other functions 
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpk),INTENT(in) :: ptem, psal ! temperature salinity
    REAL(KIND=4), DIMENSION(kpk),INTENT(in) :: pdep       ! reference depth
    INTEGER(KIND=4),             INTENT(in) :: kpk        ! size of the array
    LOGICAL,                     INTENT(in) :: ld_eos80   ! flag for differenciacion EOS80

    REAL(KIND=8), DIMENSION(kpk)            :: beta_z_80  ! returned values

    INTEGER(KIND=4) :: jk          ! dummy loop index
    REAL(KIND=8)    :: zt, zs, zh  ! working variables
    !!----------------------------------------------------------------------
    !$OMP PARALLEL DO SCHEDULE(RUNTIME)
    DO jk=1,kpk
          zh = pdep(jk)
          zt =  ptem(jk)         ! potential temperature
          zs =  psal(jk)- 35.0   ! salinity anomaly (s-35)

          beta_z_80(jk) = ( ( -0.415613e-09 * zt + 0.555579e-07 ) * zt        &   ! beta
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
    !$OMP END PARALLEL DO

  END FUNCTION beta_z_80

END MODULE eos
