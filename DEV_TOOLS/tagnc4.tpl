  LOGICAL                                    :: lnc4      = .FALSE.     ! Use nc4 with chunking and deflation

     PRINT *,' usage : cdfvita U-file V_file T-file [-w W-file] [-geo ] [-cubic] [-nc4] ...'

     PRINT *,'       [ -nc4 ]     : Use netcdf4 output with chunking and deflation level 1'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'

     CASE ( '-nc4' ) ; lnc4 = .TRUE.
  
  stypvar(ivar)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
  ncout = create      (cf_out,   cf_tfil,  npiglo, npjglo, nlev     , ld_nc4=lnc4 )
  ierr  = createvar   (ncout ,   stypvar,  nvar,   ipk,    id_varout, ld_nc4=lnc4 )
