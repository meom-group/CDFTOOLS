PROGRAM cdfmeshmask
  !!======================================================================
  !!                     ***  PROGRAM  cdfmeshmask  ***
  !!=====================================================================
  !!  ** Purpose : build mesh mask file from bathymetry
  !!
  !!  ** Method  : use nemo3.6 simplified zgr_xxx routines for zps
  !!
  !! History : 3.0  : 10/2014  : J.M. Molines 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------


  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2014
  !! $Id$
  !! Copyright (c) 2014, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
USE cdfio

   CALL zgr_z
   CALL zgr_bat
   CALL zgr_zps
   CALL zgr_bat_ctl

CONTAINS
  SUBROUTINE zgr_z
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE zgr_z  ***
    !!
    !! ** Purpose :    set the depth of model levels and the resulting 
    !!         vertical scale factors.
    !!
    !! ** Method  :   from NEMO 3.6
    !!
    !! References :  
    !!----------------------------------------------------------------------
      INTEGER  ::   jk                     ! dummy loop indices
      REAL(wp) ::   zt, zw                 ! temporary scalars
      REAL(wp) ::   zsur, za0, za1, zkth   ! Values set from parameters in
      REAL(wp) ::   zacr, zdzmin, zhmax    ! par_CONFIG_Rxx.h90
      REAL(wp) ::   zrefdep                ! depth of the reference level (~10m)
      REAL(wp) ::   za2, zkth2, zacr2      ! Values for optional double tanh function set from parameters 
   !! -----------------------------------------------------------------------------------------
      ! Set variables from parameters
      ! ------------------------------
       zkth = ppkth       ;   zacr = ppacr
       zdzmin = ppdzmin   ;   zhmax = pphmax
       zkth2 = ppkth2     ;   zacr2 = ppacr2   ! optional (ldbletanh=T) double tanh parameters

      IF(   ppa1  == pp_to_be_computed  .AND.  &
         &  ppa0  == pp_to_be_computed  .AND.  &
         &  ppsur == pp_to_be_computed           ) THEN
         !
         za1  = (  ppdzmin - pphmax / FLOAT(jpkm1)  )                                                      &
            & / ( TANH((1-ppkth)/ppacr) - ppacr/FLOAT(jpk-1) * (  LOG( COSH( (jpk - ppkth) / ppacr) )      &
            &                                                   - LOG( COSH( ( 1  - ppkth) / ppacr) )  )  )
         za0  = ppdzmin - za1 *              TANH( (1-ppkth) / ppacr )
         zsur =   - za0 - za1 * ppacr * LOG( COSH( (1-ppkth) / ppacr )  )
      ELSE
         za1 = ppa1 ;       za0 = ppa0 ;          zsur = ppsur
         za2 = ppa2                            ! optional (ldbletanh=T) double tanh parameter
      ENDIF


  END SUBROUTINE zgr_z





END PROGRAM cdfmeshmask
