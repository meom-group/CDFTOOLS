  LOGICAL               :: ll_teos10  = .FALSE.               ! teos10 flag
     PRINT *,'       [-teos10] : use TEOS10 equation of state instead of default EOS80'
     PRINT *,'                 Temperature should be conservative temperature (CT) in deg C.'
     PRINT *,'                 Salinity should be absolute salinity (SA) in g/kg.'
  CASE ( '-teos10' ) ; ll_teos10 = .TRUE. 
  CALL eos_init ( ll_teos10 )

