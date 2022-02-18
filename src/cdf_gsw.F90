PROGRAM cdf_gsw
  !!======================================================================
  !!                     ***  PROGRAM  cdf_gsw  ***
  !!=====================================================================
  !!  ** Purpose : A cdftool wrapper for GSW (Gibbs Sea Water eq. of state)
  !!               toolbox
  !!
  !!  ** Method  : Use libgsw from github:TEOS-10/GSW-Fortran.git
  !!               cdf_gsw takes the name of the GSW tool to use and
  !!               relevant argument, to produce the results in a netcdf
  !!               file
  !!
  !!  ** reference : http://www.teos-10.org/pubs/gsw/html/gsw_contents.html#1
  !!
  !! History :  4.0  : 10/2021  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
#if defined key_GSW
  USE cdfio
  USE modcdfnames
  USE modgsw
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2021
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER(KIND=4), PARAMETER                 :: jpvarmax=5
  INTEGER(KIND=4)                            :: narg, iargc         ! command line browsing
  INTEGER(KIND=4)                            :: ijarg               ! command line browsing
  INTEGER(KIND=4)                            :: npiglo, npjglo
  INTEGER(KIND=4)                            :: npk, npt
  INTEGER(KIND=4)                            :: ji,jj,jk,jt         ! dummy loop index
  INTEGER(KIND=4)                            :: ierr                ! error status of netcdf call
  INTEGER(KIND=4)                            :: ncout               ! ncid of output file
  INTEGER(KIND=4)                            :: nvar                ! number of output variables
  INTEGER(KIND=4), DIMENSION(jpvarmax)       :: ipk, id_varout      ! level and  varid's

  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE  :: v1_1d, v2_1d
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: v1,v2,v3, tmsk

  !  GSW functions take double precision values as input
  REAL(KIND=8)                               :: dpref               ! reference depth
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE  :: dv1_1d, dv2_1d
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE  :: dv1,dv2,dv3, dv4, dv5

  CHARACTER(LEN=256)                         :: cldum               ! dummy character variable
  CHARACTER(LEN=256)                         :: cf_out              ! output file name
  CHARACTER(LEN=256)                         :: cv_out              ! output variable name
  CHARACTER(LEN=256)                         :: cv_sal              ! Salinity variable
  CHARACTER(LEN=256)                         :: cv_tem              ! Temperature variable
  CHARACTER(LEN=256)                         :: cf_pt               ! input potential temp file
  CHARACTER(LEN=256)                         :: cf_ct               ! input conservative temp file
  CHARACTER(LEN=256)                         :: cf_it               ! input conservative temp file
  CHARACTER(LEN=256)                         :: cf_sp               ! input pratical salinity file
  CHARACTER(LEN=256)                         :: cf_sa               ! input absolute salinity file
  CHARACTER(LEN=256)                         :: cf_sk               ! input Knudsen salinity file [parts per thousand, ppt]
  CHARACTER(LEN=256)                         :: c_gsw_name          ! name of the GSW function or routine

  TYPE (variable), DIMENSION(jpvarmax)       :: stypvar            ! structure for attributes

  LOGICAL                                    :: lnc4

  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  cf_out='none'
  cv_out='none'
  cv_sal='none'
  cv_tem='none'
  c_gsw_name='none'
  dpref = -9999

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdf_gsw GSW-function  [-sa SA-file ] [-sp SP-file] [ -sk SK-file] ' 
     PRINT *,'      [-ct CT-file ] [-pt PT-file] [-t  TINSITU-file] [-pref REF-depth]' 
     PRINT *,'      [-h] [-l] [-vo VAR-name] [-vsal VAR-sal] [ -vtem VAR-temperature]' 
     PRINT *,'      [-nc4] [-o OUT-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        This tool is somehow different from other CDFTOOLS as it is just a '
     PRINT *,'       wrapper to Gibbs Sea Water toolbox (GSW). It took as a first argument,'
     PRINT *,'       the name of the GSW function or routine, and input files depends on '
     PRINT *,'       every particular function or routine. '
     PRINT *,'       AKNOWLEDGMENTS: http://www.teos-10.org/ '
     PRINT *,'       ' 
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       The command cdf_gsw -l give you the list of all GSW functions/routines.'
     PRINT *,'       while cdf_gsw <GSW FUNCTION> -h give you a short description of the GSW '
     PRINT *,'       function, (a copy of the original header of the GSW file). In particular'
     PRINT *,'       the required arguments for the gsw function are given (with their '
     PRINT *,'       expected unit.'
     PRINT *,'       The point is that the functions take salinity and/or temperature as'
     PRINT *,'       input parameters. They can be any of SA SP SK PT CT ... For instance if'
     PRINT *,'       a function requires sa, ct, then corresponding files are passed to the'
     PRINT *,'       program using -sa SA-file -ct CT-file. For other types of salinity or'
     PRINT *,'       temperature, there are -sp, -sk, -pt, -t options.  Important to note '
     PRINT *,'       that a depth is often required in the functions. It is taken from the '
     PRINT *,'       salinity or temperature file. In case of needs for a constant reference'
     PRINT *,'       deptht, it is passed to the code with -pref option.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        -sa SA-file : pass the name of a file with Absolute Salinity' 
     PRINT *,'        -sp SP-file : pass the name of a file with Practical Salinity' 
     PRINT *,'        -sk SK-file : pass the name of a file with Knudsen Salinity' 
     PRINT *,'        -vsal VAR-salinity : give variable name for input salinity.'
     PRINT *,'                 default is ',TRIM(cn_vosaline)
     PRINT *,'        -ct CT-file : pass the name of a file with Conservative Temperature' 
     PRINT *,'        -pt PT-file : pass the name of a file with Pot. Temperature (theta-0)'
     PRINT *,'        -t T-file   : pass the name of a file with in situ temperature.'
     PRINT *,'        -vtem VAR-temperature : give variable name for input temperature'
     PRINT *,'                 default is ',TRIM(cn_votemper)
     PRINT *,'        -pref REF-depth : pass a constant reference depth.'
     PRINT *,'        -o OUT-filename : give the name of the output file.'
     PRINT *,'                 default is <GSW_FUNCTiON>.nc'
     PRINT *,'        -vo VAR-name : give variable name for output file'
     PRINT *,'        -h : provide help on the function given as first argument.'
     PRINT *,'             In particular indicates the required imput files and the'
     PRINT *,'             way to specify them (which key to use). Users are encouraged to'
     PRINT *,'             visit http://www.teos-10.org/pubs/gsw/html/gsw_contents.html for'
     PRINT *,'             details and science  discussion.'
     PRINT *,'        -l : provide a list of all the GSW functions or routine.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        mask.nc'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ', TRIM(cv_out),' (    )'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      ' 
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  CALL getarg(ijarg, cldum) 
  IF ( TRIM(cldum) /= '-l') THEN
     c_gsw_name = cldum
     ijarg=ijarg +1
  ENDIF

  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-ct'  ) ; CALL getarg(ijarg, cf_ct ) ; ijarg=ijarg+1
     CASE ( '-sa'  ) ; CALL getarg(ijarg, cf_sa ) ; ijarg=ijarg+1
     CASE ( '-pt'  ) ; CALL getarg(ijarg, cf_pt ) ; ijarg=ijarg+1
     CASE ( '-sp'  ) ; CALL getarg(ijarg, cf_sp ) ; ijarg=ijarg+1
     CASE ( '-t'   ) ; CALL getarg(ijarg, cf_it ) ; ijarg=ijarg+1
        ! option
     CASE ( '-vo'  ) ; CALL getarg(ijarg, cv_out) ; ijarg=ijarg+1
     CASE ( '-vsal') ; CALL getarg(ijarg, cv_sal) ; ijarg=ijarg+1
     CASE ( '-vtem') ; CALL getarg(ijarg, cv_tem) ; ijarg=ijarg+1
     CASE ( '-pref') ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) dpref
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE ( '-h'   ) ; CALL gsw_help(c_gsw_name)  ; STOP
     CASE ( '-l'   ) ; CALL gsw_lst               ; STOP
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  IF ( cf_out == 'none' ) THEN
     cf_out=TRIM(c_gsw_name)//".nc"
  ENDIF

  SELECT CASE ( c_gsw_name ) 
  CASE ( 'gsw_add_barrier' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_add_mean' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_adiabatic_lapse_rate_from_ct' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_adiabatic_lapse_rate_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_alpha' ) 
     IF ( chkfile(cf_sa) .OR. chkfile(cf_ct) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     ALLOCATE( dv1_1d(npk) )
     IF ( cv_out == 'none') cv_out ='alpha'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=1
     stypvar(1)%cunits     = '1/K'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 30.
     stypvar(1)%clong_name = 'Thermal Expansion coefficient wrt CT and SA'

     CALL CreateOutput(cf_out,cv_out,cf_sa, nvar, pdep=dv1_1d)

     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_ct, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_alpha (dv1(ji,jj), dv2(ji,jj) , dv1_1d(jk))  * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2, v3 ,dv1_1d,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_alpha_on_beta' ) 
     IF ( chkfile(cf_sa) .OR. chkfile(cf_ct) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     ALLOCATE( dv1_1d(npk) )
     IF ( cv_out == 'none') cv_out ='alpha_on_beta'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=1
     stypvar(1)%cunits     = 'kg g^-1 K^-1'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 30.
     stypvar(1)%clong_name = 'Thermal Expansion coefficient wrt CT and SA divided by the saline contraction coefficient'

     CALL CreateOutput(cf_out,cv_out,cf_sa, nvar, pdep=dv1_1d)

     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_ct, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_alpha_on_beta (dv1(ji,jj), dv2(ji,jj) , dv1_1d(jk))  * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2, v3 ,dv1_1d,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_alpha_wrt_t_exact' ) 
     IF ( chkfile(cf_sa) .OR. chkfile(cf_ct) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     ALLOCATE( dv1_1d(npk) )
     IF ( cv_out == 'none') cv_out ='alpha_wrt_t'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=1
     stypvar(1)%cunits     = 'K^-1'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 30.
     stypvar(1)%clong_name = 'thermal expansion coefficient wrt insitu temperature'

     CALL CreateOutput(cf_out,cv_out,cf_sa, nvar, pdep=dv1_1d)

     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_it, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_alpha_wrt_t_exact (dv1(ji,jj), dv2(ji,jj) , dv1_1d(jk))  * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2, v3 ,dv1_1d,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_alpha_wrt_t_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_beta' ) 
     IF ( chkfile(cf_sa) .OR. chkfile(cf_ct) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     ALLOCATE( dv1_1d(npk) )
     IF ( cv_out == 'none') cv_out ='beta'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=1
     stypvar(1)%cunits     = 'kg/g'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 30.
     stypvar(1)%clong_name = 'Saline contraction coefficient  wrt CT and SA'

     CALL CreateOutput(cf_out,cv_out,cf_sa, nvar, pdep=dv1_1d)

     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_ct, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_beta (dv1(ji,jj), dv2(ji,jj) , dv1_1d(jk))  * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2, v3 ,dv1_1d,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_beta_const_t_exact' ) 
     IF ( chkfile(cf_sa) .OR. chkfile(cf_ct) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     ALLOCATE( dv1_1d(npk) )
     IF ( cv_out == 'none') cv_out ='beta_const_t'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=1
     stypvar(1)%cunits     = 'kg/g'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 30.
     stypvar(1)%clong_name = 'Saline contraction coefficient  wrt in situ temperature and SA'

     CALL CreateOutput(cf_out,cv_out,cf_sa, nvar, pdep=dv1_1d)

     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_it, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_beta_const_t_exact (dv1(ji,jj), dv2(ji,jj) , dv1_1d(jk))  * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2, v3 ,dv1_1d,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_c_from_sp' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_cabbeling' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_chem_potential_water_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_chem_potential_water_t_exact' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_cp_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ct_first_derivatives' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ct_first_derivatives_wrt_t_exact' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ct_freezing' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ct_freezing_exact' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ct_freezing_first_derivatives' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ct_freezing_first_derivatives_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ct_freezing_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ct_from_enthalpy' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ct_from_enthalpy_exact' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ct_from_entropy' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ct_from_pt' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ct_from_rho' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ct_from_t' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ct_maxdensity' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ct_second_derivatives' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_deltasa_atlas' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_deltasa_from_sp' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_dilution_coefficient_t_exact' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_dynamic_enthalpy' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_enthalpy' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_enthalpy_ct_exact' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_enthalpy_diff' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_enthalpy_first_derivatives' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_enthalpy_first_derivatives_ct_exact' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_enthalpy_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_enthalpy_second_derivatives' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_enthalpy_second_derivatives_ct_exact' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_enthalpy_sso_0' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_enthalpy_t_exact' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_entropy_first_derivatives' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_entropy_from_pt' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_entropy_from_t' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_entropy_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_entropy_part' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_entropy_part_zerop' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_entropy_second_derivatives' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_fdelta' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_frazil_properties' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_frazil_properties_potential' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_frazil_properties_potential_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_frazil_ratios_adiabatic' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_frazil_ratios_adiabatic_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_geo_strf_dyn_height' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_geo_strf_dyn_height_pc' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_gibbs' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_gibbs_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_gibbs_ice_part_t' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_gibbs_ice_pt0' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_gibbs_ice_pt0_pt0' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_gibbs_pt0_pt0' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_grav' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_helmholtz_energy_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_hill_ratio_at_sp2' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ice_fraction_to_freeze_seawater' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_internal_energy' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_internal_energy_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_ipv_vs_fnsquared_ratio' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_kappa' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_kappa_const_t_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_kappa_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_kappa_t_exact' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_latentheat_evap_ct' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_latentheat_evap_t' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_latentheat_melting' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_linear_interp_sa_ct' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_melting_ice_equilibrium_sa_ct_ratio' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_melting_ice_equilibrium_sa_ct_ratio_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_melting_ice_into_seawater' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_melting_ice_sa_ct_ratio' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_melting_ice_sa_ct_ratio_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_melting_seaice_equilibrium_sa_ct_ratio' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_melting_seaice_equilibrium_sa_ct_ratio_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_melting_seaice_into_seawater' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_melting_seaice_sa_ct_ratio' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_melting_seaice_sa_ct_ratio_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_mlp' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_nsquared' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_nsquared_lowerlimit' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_nsquared_min' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_nsquared_min_const_t' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_p_from_z' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pot_enthalpy_from_pt_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pot_enthalpy_from_pt_ice_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pot_enthalpy_ice_freezing' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pot_enthalpy_ice_freezing_first_derivatives' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pot_enthalpy_ice_freezing_first_derivatives_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pot_enthalpy_ice_freezing_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pot_rho_t_exact' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pressure_coefficient_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pressure_freezing_ct' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pt0_cold_ice_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pt0_from_t' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pt0_from_t_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pt_first_derivatives' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pt_from_ct' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pt_from_entropy' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pt_from_pot_enthalpy_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pt_from_pot_enthalpy_ice_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pt_from_pot_enthalpy_ice_poly_dh' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pt_from_t' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pt_from_t_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_pt_second_derivatives' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_rho' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_rho_alpha_beta' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_rho_alpha_beta_bsq' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_rho_first_derivatives' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_rho_first_derivatives_wrt_enthalpy' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_rho_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_rho_second_derivatives' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_rho_second_derivatives_wrt_enthalpy' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_rho_t_exact' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_rr68_interp_sa_ct' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sa_freezing_estimate' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sa_freezing_from_ct' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sa_freezing_from_ct_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sa_freezing_from_t' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sa_freezing_from_t_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sa_from_rho' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sa_from_sp' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sa_from_sstar' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sa_p_inrange' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_saar' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_seaice_fraction_to_freeze_seawater' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sigma0' ) 
     IF ( chkfile(cf_sa) .OR. chkfile(cf_ct) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     IF ( cv_out == 'none') cv_out ='sigma_0'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=1
     stypvar(1)%cunits     = 'kg/m3'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 30.
     stypvar(1)%clong_name = 'Potential_density:sigma-0'

     CALL CreateOutput(cf_out,cv_out,cf_sa, nvar)

     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_ct, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_sigma0 (dv1(ji,jj), dv2(ji,jj) ) * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2, v3 ,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_sigma1' ) 
     IF ( chkfile(cf_sa) .OR. chkfile(cf_ct) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     IF ( cv_out == 'none') cv_out ='sigma_1'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=1
     stypvar(1)%cunits     = 'kg/m3'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 30.
     stypvar(1)%clong_name = 'Potential_density:sigma-1'

     CALL CreateOutput(cf_out,cv_out,cf_sa, nvar)

     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_ct, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_sigma1 (dv1(ji,jj), dv2(ji,jj) ) * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2, v3 ,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_sigma2' ) 
     IF ( chkfile(cf_sa) .OR. chkfile(cf_ct) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     IF ( cv_out == 'none') cv_out ='sigma_2'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=1
     stypvar(1)%cunits     = 'kg/m3'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 50.
     stypvar(1)%clong_name = 'Potential_density:sigma-2'

     CALL CreateOutput(cf_out,cv_out,cf_sa,nvar)

     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_ct, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_sigma2 (dv1(ji,jj), dv2(ji,jj) ) * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2, v3 ,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_sigma3' ) 
     IF ( chkfile(cf_sa) .OR. chkfile(cf_ct) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     IF ( cv_out == 'none') cv_out ='sigma_3'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=1
     stypvar(1)%cunits     = 'kg/m3'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 50.
     stypvar(1)%clong_name = 'Potential_density:sigma-3'

     CALL CreateOutput(cf_out,cv_out,cf_sa,nvar)
  
     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_ct, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_sigma3 (dv1(ji,jj), dv2(ji,jj) ) * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2, v3 ,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_sigma4' ) 
     IF ( chkfile(cf_sa) .OR. chkfile(cf_ct) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     IF ( cv_out == 'none') cv_out ='sigma_4'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=1
     stypvar(1)%cunits     = 'kg/m3'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 50.
     stypvar(1)%clong_name = 'Potential_density:sigma-4'

     CALL CreateOutput(cf_out,cv_out,cf_sa,nvar)
  
     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_ct, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_sigma4 (dv1(ji,jj), dv2(ji,jj) ) * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2, v3 ,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_sound_speed' ) 
     IF ( chkfile(cf_sa) .OR. chkfile(cf_ct) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     ALLOCATE( dv1_1d(npk) )
     IF ( cv_out == 'none') cv_out ='sound_speed'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=1
     stypvar(1)%cunits     = 'm/s'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 2000.
     stypvar(1)%clong_name = 'speed_of_sound_in_seawater'

     CALL CreateOutput(cf_out,cv_out,cf_sa,nvar, pdep=dv1_1d)
    
     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_ct, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_sound_speed (dv1(ji,jj), dv2(ji,jj), dv1_1d(jk) ) * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2, v3 ,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_sound_speed_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sound_speed_t_exact' ) 
     IF ( chkfile(cf_sa) .OR. chkfile(cf_it) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     ALLOCATE( dv1_1d(npk) )
     IF ( cv_out == 'none') cv_out ='sound_speed_exact'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=1
     stypvar(1)%cunits     = 'm/s'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 2000.
     stypvar(1)%clong_name = 'sound_speed_exact'

     CALL CreateOutput(cf_out,cv_out,cf_sa,nvar, pdep=dv1_1d)

     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_it, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_sound_speed_t_exact (dv1(ji,jj), dv2(ji,jj), dv1_1d(jk) ) * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2, v3 ,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_sp_from_c' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sp_from_sa' ) 
     IF ( chkfile(cf_sa) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     ALLOCATE( dv1_1d(npk) )
     IF ( cv_out == 'none') cv_out ='vosaline'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     nvar=1
     stypvar(1)%cunits     = 'PSU'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 80.
     stypvar(1)%clong_name = 'Practical_Salinity'

     CALL CreateOutput(cf_out,cv_out,cf_sa,nvar, pdep=dv1_1d)

     dv2 = getvar(cf_sa, cn_vlon2d, 1, npiglo, npjglo)
     dv3 = getvar(cf_sa, cn_vlat2d, 1, npiglo, npjglo)

     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_sp_from_sa (dv1(ji,jj), dv1_1d(jk),dv2(ji,jj),dv3(ji,jj) ) * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2,dv3, v3 ,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_sp_from_sk' ) 
     IF ( chkfile(cf_sk) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sk)
     ALLOCATE( dv1(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     IF ( cv_out == 'none') cv_out ='vosaline'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     nvar=1
     stypvar(1)%cunits     = 'PSU'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 80.
     stypvar(1)%clong_name = 'Practical_Salinity'

     CALL CreateOutput(cf_out,cv_out,cf_sa,nvar)

     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sk, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_sp_from_sk (dv1(ji,jj) ) * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2,dv3, v3 ,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_sp_from_sr' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sp_from_sstar' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_specvol' ) 
     IF ( chkfile(cf_sa) .OR. chkfile(cf_ct) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     ALLOCATE( dv1_1d(npk) )
     IF ( cv_out == 'none') cv_out ='specvol'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=1
     stypvar(1)%cunits     = 'm3/kg'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 0.02
     stypvar(1)%clong_name = 'specific_volume'

     CALL CreateOutput(cf_out,cv_out,cf_sa,nvar, pdep=dv1_1d)

     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_ct, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_specvol (dv1(ji,jj), dv2(ji,jj), dv1_1d(jk) ) * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2, v3 ,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_specvol_alpha_beta' ) 
     IF ( chkfile(cf_sa) .OR. chkfile(cf_ct) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo), dv3(npiglo,npjglo))
     ALLOCATE( dv4(npiglo,npjglo), dv5(npiglo,npjglo))
     ALLOCATE(  v1(npiglo,npjglo),  v2(npiglo,npjglo),  v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     ALLOCATE( dv1_1d(npk) )
     IF ( cv_out == 'none') cv_out ='specvol'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=3
     stypvar(1)%cname       = 'specvol'
     stypvar(1)%cshort_name = 'specvol'
     stypvar(1)%cunits     = 'm3/kg'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 0.02
     stypvar(1)%clong_name = 'specific_volume'

     stypvar(2)%cname       = 'alpha'
     stypvar(2)%cshort_name = 'alpha'
     stypvar(2)%cunits     = '1/K'
     stypvar(2)%valid_min  =  0.
     stypvar(2)%valid_max  = 0.02
     stypvar(2)%clong_name = 'thermal_expansion_coef_CT'

     stypvar(3)%cname       = 'beta'
     stypvar(3)%cshort_name = 'beta'
     stypvar(3)%cunits     = 'kg/g'
     stypvar(3)%valid_min  =  0.
     stypvar(3)%valid_max  = 0.02
     stypvar(3)%clong_name = 'haline_contraction_coef_cst_CT'

     CALL CreateOutput(cf_out,cv_out,cf_sa, nvar, pdep=dv1_1d)

     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_ct, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 CALL gsw_specvol_alpha_beta( dv1(ji,jj), dv2(ji,jj), dv1_1d(jk), dv3(ji,jj), &
                                              dv4(ji,jj), dv5(ji,jj) ) 
                 v1(ji,jj)  = dv3(ji,jj)* tmsk(ji,jj)
                 v2(ji,jj)  = dv4(ji,jj)* tmsk(ji,jj)
                 v3(ji,jj)  = dv5(ji,jj)* tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v1, jk, npiglo, npjglo, ktime=jt )
           ierr = putvar(ncout, id_varout(2), v2, jk, npiglo, npjglo, ktime=jt )
           ierr = putvar(ncout, id_varout(3), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO 
     DEALLOCATE ( dv1, dv2, dv3 , dv4, dv5, v1, v2, v3 ,tmsk, dv1_1d)
     ierr = closeout(ncout)
  CASE ( 'gsw_specvol_anom_standard' ) 
     IF ( chkfile(cf_sa) .OR. chkfile(cf_ct) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo),v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     ALLOCATE( dv1_1d(npk) )
     IF ( cv_out == 'none') cv_out ='specvol_anom_std'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=1
     stypvar(1)%cunits     = 'm3/kg'
     stypvar(1)%valid_min  =  0.
     stypvar(1)%valid_max  = 0.02
     stypvar(1)%clong_name = 'specific_volume_anomaly_of_seawater'


     CALL CreateOutput(cf_out,cv_out,cf_sa,nvar, pdep=dv1_1d)

     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_ct, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 v3(ji,jj)  = gsw_specvol_anom_standard (dv1(ji,jj), dv2(ji,jj), dv1_1d(jk) ) * tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO
     DEALLOCATE ( dv1, dv2, v3 ,tmsk)
     ierr = closeout(ncout)
  CASE ( 'gsw_specvol_first_derivatives' ) 
!WIP
     IF ( chkfile(cf_sa) .OR. chkfile(cf_ct) .OR. chkfile( cn_fmsk) )  STOP 99
     CALL SetDims(cf_sa)
     ALLOCATE( dv1(npiglo,npjglo), dv2(npiglo,npjglo), dv3(npiglo,npjglo))
     ALLOCATE( dv4(npiglo,npjglo), dv5(npiglo,npjglo))
     ALLOCATE(  v1(npiglo,npjglo),  v2(npiglo,npjglo),  v3(npiglo,npjglo), tmsk(npiglo,npjglo) )
     ALLOCATE( dv1_1d(npk) )
     IF ( cv_out == 'none') cv_out ='specvol'
     IF ( cv_sal /= 'none') cn_vosaline = cv_sal
     IF ( cv_tem /= 'none') cn_votemper = cv_tem
     nvar=3
     stypvar(1)%cname       = 'v_SA'
     stypvar(1)%cshort_name = 'v_SA'
     stypvar(1)%cunits     = 'J/(kg (g/kg)^2)'
     stypvar(1)%valid_min  =  -10
     stypvar(1)%valid_max  = 10.
     stypvar(1)%clong_name = 'first_deriv_spec_volume_wrt_SA_at_cst_CT_P'

     stypvar(2)%cname       = 'v_CT'
     stypvar(2)%cshort_name = 'v_CT'
     stypvar(2)%cunits     = 'J/(kg K(g/kg))'
     stypvar(2)%valid_min  =  -10.
     stypvar(2)%valid_max  =  10.
     stypvar(2)%clong_name = 'first_deriv_spec_volume_wrt_CT_at_cst_SA_CT'

     stypvar(3)%cname       = 'v_P'
     stypvar(3)%cshort_name = 'v_P'
     stypvar(3)%cunits     = 'J/(kg K^2)'
     stypvar(3)%valid_min  =  0.
     stypvar(3)%valid_max  = 0.02
     stypvar(3)%clong_name = 'first_deriv_spec_volume_wrt_P_at_cst_SA_CT'

     CALL CreateOutput(cf_out,cv_out,cf_sa, nvar, pdep=dv1_1d)

     DO jt = 1, npt
        PRINT *, 'JT = ', jt,'/',npt
        DO jk = 1, npk
           tmsk(:,:)= getvar(cn_fmsk,cn_tmask, jk, npiglo, npjglo)
           dv1(:,:) = getvar(cf_sa, cn_vosaline, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           dv2(:,:) = getvar(cf_ct, cn_votemper, jk, npiglo, npjglo, ktime=jt) * tmsk(:,:)
           DO jj = 1,npjglo
              DO ji = 1,npiglo
                 CALL gsw_specvol_first_derivatives( dv1(ji,jj), dv2(ji,jj), dv1_1d(jk), dv3(ji,jj), &
                                              dv4(ji,jj), dv5(ji,jj), 14 ) 
                 v1(ji,jj)  = dv3(ji,jj)* tmsk(ji,jj)
                 v2(ji,jj)  = dv4(ji,jj)* tmsk(ji,jj)
                 v3(ji,jj)  = dv5(ji,jj)* tmsk(ji,jj)
              ENDDO
           ENDDO
           ierr = putvar(ncout, id_varout(1), v1, jk, npiglo, npjglo, ktime=jt )
           ierr = putvar(ncout, id_varout(2), v2, jk, npiglo, npjglo, ktime=jt )
           ierr = putvar(ncout, id_varout(3), v3, jk, npiglo, npjglo, ktime=jt )
        ENDDO
     ENDDO 
     DEALLOCATE ( dv1, dv2, dv3 , dv4, dv5, v1, v2, v3 ,tmsk, dv1_1d)
     ierr = closeout(ncout)

  CASE ( 'gsw_specvol_first_derivatives_wrt_enthalpy' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_specvol_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_specvol_second_derivatives' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_specvol_second_derivatives_wrt_enthalpy' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_specvol_sso_0' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_specvol_t_exact' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_spiciness0' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_spiciness1' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_spiciness2' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sr_from_sp' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sstar_from_sa' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_sstar_from_sp' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_t_deriv_chem_potential_water_t_exact' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_t_freezing' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_t_freezing_exact' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_t_freezing_first_derivatives' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_t_freezing_first_derivatives_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_t_freezing_poly' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_t_from_ct' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_t_from_pt0_ice' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_thermobaric' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_turner_rsubrho' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_util_indx' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_util_interp1q_int' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_util_xinterp1' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  CASE ( 'gsw_z_from_p' ) 
     PRINT *, "This function/subroutine is not yet ready." 
     STOP  
  END SELECT

CONTAINS

  SUBROUTINE SetDims(cd_file)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE SetDims  ***
    !!
    !! ** Purpose :   Set file dimension from file given in argument
    !!
    !! ** Method  :   use cdfio/getdim
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_file
    !!----------------------------------------------------------------------
    npiglo = getdim (cd_file, cn_x)
    npjglo = getdim (cd_file, cn_y)
    npk    = getdim (cd_file, cn_z)
    npt    = getdim (cd_file, cn_t)
  END SUBROUTINE SetDims

  SUBROUTINE CreateOutput (cd_fout, cd_vout, cd_fref, knvar, pdep)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create output netcdf file and keep variable info
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),                      INTENT(in) :: cd_fout
    CHARACTER(LEN=*),                      INTENT(in) :: cd_vout
    CHARACTER(LEN=*),                      INTENT(in) :: cd_fref
    INTEGER(KIND=4),                       INTENT(in) :: knvar
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(out) :: pdep

    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dltim
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: zdep
    INTEGER(KIND=4)                         :: jv
    !!----------------------------------------------------------------------
    ALLOCATE (dltim(npt), zdep(npk) )
    IF (knvar == 1 ) THEN
       stypvar(1)%cname             = cd_vout
       stypvar(1)%cshort_name       = cd_vout
    ENDIF
    DO jv = 1 , knvar
       ipk(jv) = npk
       stypvar(jv)%rmissing_value    = 0.
       stypvar(jv)%conline_operation = 'N/A'
       stypvar(jv)%caxis             = 'TZYX'
       stypvar(jv)%ichunk            = (/npiglo, MAX(1,npjglo/30), 1, 1 /)
    ENDDO

    ! create output fileset
    ncout = create      (cd_fout, cd_fref, npiglo, npjglo, npk,       ld_nc4=lnc4  )
    ierr  = createvar   (ncout,  stypvar, knvar,  ipk,     id_varout, ld_nc4=lnc4  )
    ierr  = putheadervar(ncout,  cd_fref, npiglo, npjglo, npk       )

    dltim = getvar1d(cd_fref, cn_vtimec, npt    )
    ierr  = putvar1d(ncout,  dltim,      npt, 'T')

    zdep = getvar1d(cd_fref, cn_vdeptht, npk     )
    ierr = putvar1d(ncout,  zdep,      npk, 'D')
    IF ( PRESENT(pdep) ) THEN 
      IF ( dpref == -9999 ) THEN
         ! return depth array from ref file
         pdep = zdep
      ELSE
         ! set pdep array to constant value dpref
         pdep = dpref
      ENDIF
    ENDIF
    DEALLOCATE( dltim, zdep )

  END SUBROUTINE CreateOutput

#else
  PRINT *,' cdf_gsw not available.'
  PRINT *,'   need to compile with GSW = -D key_GSW'
  PRINT *,'   GSW library must be compile beforehand'
#endif


END PROGRAM cdf_gsw
