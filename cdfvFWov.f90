PROGRAM cdfvFWov
  !!-------------------------------------------------------------------
  !!======================================================================
  !!                     ***  PROGRAM  cdfvFWov  ***
  !!=====================================================================
  !!  ** Purpose : from a section calculate net freshwater transport and its
  !!               overturning component
  !!               section is assumed to be 2 j lines (j and j+1)
  !!
  !!  ** Method  : compute salinity at v point
  !!               compute zonal mean of salinity-vpt and v velocity
  !!               compute total freshwater transport and overturning component
  !!
  !! History : 2.1  : 12/2011  : J. Deshayes  : Original code
  !!           3.0  : 12/2011  : J.M. Molines : Port to 3.0
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE netcdf  ! to be eliminated
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)   :: varid_time, varid_net, varid_tot, varid_ov
  INTEGER(KIND=4)   :: dimid_time 
  INTEGER(KIND=4)   :: ncin, ncout

  INTEGER(KIND=4)                            :: narg, iargc
  INTEGER(KIND=4)                            :: npt, npiglo, npjglo, npk
  INTEGER(KIND=4)                            :: ji, jk, jt                           
  INTEGER(KIND=4)                            :: istatus

  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: rmasks, rmaskn, rmaskv
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: zsaln, zsals
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: zsalv, zvitv
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: zFWv

  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE :: de3v
  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE :: de1v
  REAL(KIND=8), DIMENSION (:),   ALLOCATABLE :: dnetvFW, dtotvFW
  REAL(KIND=8), DIMENSION (:),   ALLOCATABLE :: dovFW, dtime
  REAL(KIND=8), DIMENSION (1)                :: dnetv, dnetFW
  REAL(KIND=8), DIMENSION (1)                :: darea, dareak
  REAL(KIND=8), DIMENSION (1)                :: dzonalv, dzonalFW
  REAL(KIND=8)                               :: dztrp, dcellarea
  REAL(KIND=8), PARAMETER                    :: dp_rsal=35.d0      ! reference salinity for freshwater calculation

  CHARACTER(LEN=256)                         :: cf_vfil
  CHARACTER(LEN=256)                         :: cf_sfil
  CHARACTER(LEN=256)                         :: cf_zgr
  CHARACTER(LEN=256)                         :: cf_hgr
  CHARACTER(LEN=256)                         :: cf_mask
  CHARACTER(LEN=256)                         :: cf_out    = 'vFWov.nc'
  CHARACTER(LEN=256)                         :: cv_netvFW = 'netvFW'
  CHARACTER(LEN=256)                         :: cv_totvFW = 'totvFW'
  CHARACTER(LEN=256)                         :: cv_ovFW   = 'ovFW'

  TYPE (variable), DIMENSION(4)              :: stypvar    

  LOGICAL                                    :: lchk = .FALSE.
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  !!  Read command line 
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfvFWov V-secfile S-secfile ZGR-secfile HGR-secfile MSK-secfile'
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute the fresh water transport and its overturning component through'
     PRINT *,'        a section specified by the input files (data and metrics).'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        All arguments are ''section file'', which are assumed to be files with'
     PRINT *,'        2 zonal lines of data ( j and j+1). '
     PRINT *,'        V_secfile : meridional velocity section file.'
     PRINT *,'        S_secfile : salinity section file.'
     PRINT *,'        ZGR_secfile : mesh_zgr section file '
     PRINT *,'        HGR_secfile : mesh_hgr section file '
     PRINT *,'        MSK_secfile : mask section file '
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out)
     PRINT *,'       variables : ',TRIM(cv_netvFW),', ',TRIM(cv_totvFW),', ',TRIM(cv_ovFW)
     STOP
  ENDIF

  !! get arguments
  CALL getarg (1, cf_vfil)
  CALL getarg (2, cf_sfil)
  CALL getarg (3, cf_zgr )
  CALL getarg (4, cf_hgr )
  CALL getarg (5, cf_mask)

  lchk = lchk .OR. chkfile ( cf_vfil )
  lchk = lchk .OR. chkfile ( cf_sfil )
  lchk = lchk .OR. chkfile ( cf_zgr  )
  lchk = lchk .OR. chkfile ( cf_hgr  )
  lchk = lchk .OR. chkfile ( cf_mask )

  IF ( lchk ) STOP ! missing files

  !! get dimensions of input file containing data
  npiglo = getdim(cf_vfil, cn_x)
  npjglo = getdim(cf_vfil, cn_y)
  npk    = getdim(cf_vfil, cn_z)
  npt    = getdim(cf_vfil, cn_t)  

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  IF ( npjglo /= 2 ) THEN
     PRINT *,' ERROR : This program works with section files.'
     PRINT *,'       all data should have j dimension equal to 2 '
     STOP
  ENDIF

  ALLOCATE( dnetvFW(npt), dtotvFW(npt), dovFW(npt), dtime(npt) )
  ALLOCATE( de1v(npiglo,npjglo), de3v(npiglo,npk))
  ALLOCATE( zFWv(npiglo,npk), zvitv(npiglo,npk) )
  ALLOCATE( zsals(npiglo,npk), zsaln(npiglo,npk), zsalv(npiglo,npk) )
  ALLOCATE( rmasks(npiglo,npk), rmaskn(npiglo,npk), rmaskv(npiglo,npk))

  !! load data
  dtime = getvar1d(cf_vfil, cn_vtimec, npt)

  !! define output variables
  stypvar%rmissing_value    = 0.
  stypvar%valid_min         = -1000.
  stypvar%valid_max         = 1000.
  stypvar%conline_operation = 'N/A'
  stypvar%caxis             = 'T'
  stypvar%scale_factor      = 1.
  stypvar%add_offset        = 0.
  stypvar%savelog10         = 0.

  stypvar(1)%cname          = TRIM(cv_netvFW )
  stypvar(2)%cname          = TRIM(cv_totvFW )
  stypvar(3)%cname          = TRIM(cv_ovFW   )

  stypvar(1)%cunits         ='Sv'
  stypvar(2)%cunits         ='Sv'
  stypvar(3)%cunits         ='Sv'

  stypvar(1)%clong_name     = 'Net_VFW'
  stypvar(2)%clong_name     = 'Total_VFW'
  stypvar(3)%clong_name     = 'Overturning_VFW'

  stypvar(1)%cshort_name    = TRIM(cv_netvFW )
  stypvar(2)%cshort_name    = TRIM(cv_totvFW )
  stypvar(3)%cshort_name    = TRIM(cv_ovFW   )


  !! load scale factors 
  rmasks(:,:) = getvarxz(cf_mask, 'tmask', kj=1,   kpi=npiglo, kpz=npk,    kimin=1, kkmin=1)
  rmaskn(:,:) = getvarxz(cf_mask, 'tmask', kj=2,   kpi=npiglo, kpz=npk,    kimin=1, kkmin=1)
  rmaskv(:,:) = getvarxz(cf_mask, 'vmask', kj=1,   kpi=npiglo, kpz=npk,    kimin=1, kkmin=1)
  de1v(:,:)   = getvar  (cf_hgr, cn_ve1v,  klev=1, kpi=npiglo, kpj=npjglo, kimin=1, kjmin=1)
  de3v(:,:)   = getvarxz(cf_zgr, 'e3v',    kj=1,   kpi=npiglo, kpz=npk,    kimin=1, kkmin=1)

  WHERE ( rmasks /= 0. ) rmasks = 1.
  WHERE ( rmaskn /= 0. ) rmaskn = 1.
  WHERE ( rmaskv /= 0. ) rmaskv = 1.

  !! do calculation for each time step
  DO jt=1,npt
     PRINT *,'jt =',jt
     ! reset cumulative arrays to 0
     zFWv(:,:)   = 0.e0
     dnetv       = 0.d0  ;  dnetFW = 0.d0 
     darea       = 0.d0
     dnetvFW(jt) = 0.d0  ;  dtotvFW(jt) = 0.d0  ; dovFW(jt) = 0.d0

     zvitv(:,:)= getvarxz(cf_vfil, cn_vomecrty, kj=1, kpi=npiglo, kpz=npk, kimin=1, kkmin=1, ktime=jt) 
     zsals(:,:)= getvarxz(cf_sfil, cn_vosaline, kj=1, kpi=npiglo, kpz=npk, kimin=1, kkmin=1, ktime=jt)
     zsaln(:,:)= getvarxz(cf_sfil, cn_vosaline, kj=2, kpi=npiglo, kpz=npk, kimin=1, kkmin=1, ktime=jt)

     DO jk = 1, npk
        DO ji = 1, npiglo
           IF ( rmasks(ji,jk) + rmaskn(ji,jk) /=  0 )  THEN
              zFWv(ji,jk) = ( dp_rsal - ( zsals(ji,jk) * rmasks(ji,jk) + zsaln(ji,jk) * rmaskn(ji,jk) ) &
                   &           / ( rmasks(ji,jk) + rmaskn(ji,jk) ) ) / dp_rsal ! freshwater at Vpoint
           ENDIF
           dcellarea   = de1v(ji,1) * de3v(ji,jk) * rmaskv(ji,jk)
           dztrp       = zvitv(ji,jk) * dcellarea
           dnetvFW(jt) = dnetvFW(jt) + zFWv(ji,jk) * dztrp /1.d6         ! net freshwater transport in Sv
           dnetv       = dnetv       +               dztrp
           dnetFW      = dnetFW      + zFWv(ji,jk) * dcellarea
           darea       = darea       + dcellarea
        END DO !ji
     END DO !jk

     PRINT *,'total mass transport across section =', dnetv / 1.d6,' Sv'

     dnetv  = dnetv  / darea   ! mean velocity across section
     dnetFW = dnetFW / darea   ! mean freshwater along section
     PRINT *,'mean salinity along section =', dp_rsal - dnetFW * dp_rsal,' psu'

     DO jk = 1, npk
        dzonalv  = 0.d0
        dzonalFW = 0.d0
        dareak   = 0.d0 
        DO ji = 1, npiglo
           dcellarea = de1v(ji,1) * de3v(ji,jk) * rmaskv(ji,jk)
           dzonalv   = dzonalv  + ( zvitv(ji,jk) - dnetv  ) * dcellarea
           dzonalFW  = dzonalFW + ( zFWv(ji,jk)  - dnetFW ) * dcellarea
           dareak    = dareak   +                             dcellarea

           dtotvFW(jt) = dtotvFW(jt) + ( zvitv(ji,jk) - dnetv(1) ) * zFWv(ji,jk) * dcellarea /1.d6 
        END DO !ji

        IF ( dareak(1) > 0 ) THEN
           ! overturning freshwater transport in Sv
           dovFW(jt) = dovFW(jt) + dzonalv(1) * dzonalFW(1) / dareak(1) /1.d6  
        ENDIF
     END DO !jk

     PRINT *,'netvFW = ', dnetvFW(jt), ' Sv'
     PRINT *,'totvFW = ', dtotvFW(jt), ' Sv'
     PRINT *,'ovFW   = ', dovFW(jt),   ' Sv'

  END DO  !jt

  !! create output file
  istatus = NF90_CREATE(cf_out, NF90_CLOBBER, ncout)
  istatus = NF90_DEF_DIM(ncout, cn_vtimec, NF90_UNLIMITED, dimid_time)
  istatus = NF90_OPEN(cf_vfil, NF90_NOWRITE, ncin)
  istatus = NF90_DEF_VAR(ncout, cn_vtimec, NF90_FLOAT, (/dimid_time/), varid_time)
  istatus = copyatt(cn_vtimec, varid_time, ncin, ncout)
  istatus = closeout(ncin)
  istatus = NF90_DEF_VAR(ncout, 'netvFW', NF90_FLOAT, (/dimid_time/), varid_net)
  istatus = NF90_PUT_ATT(ncout, varid_net, 'units', 'Sv')
  istatus = NF90_PUT_ATT(ncout, varid_net, 'long_name', 'Net transport of freshwater across section')
  istatus = NF90_PUT_ATT(ncout, varid_net, 'comment', 'Reference salinity 35.0, transport is positive northward')
  istatus = NF90_DEF_VAR(ncout, 'totvFW', NF90_FLOAT, (/dimid_time/), varid_tot)
  istatus = NF90_PUT_ATT(ncout, varid_tot, 'units', 'Sv')
  istatus = NF90_PUT_ATT(ncout, varid_tot, 'long_name', 'Transport of freshwater across section when net mass transport equals 0')
  istatus = NF90_PUT_ATT(ncout, varid_tot, 'comment', 'Reference salinity 35.0, transport is positive northward')
  istatus = NF90_DEF_VAR(ncout, 'ovFW', NF90_FLOAT, (/dimid_time/), varid_ov)
  istatus = NF90_PUT_ATT(ncout, varid_ov, 'units', 'Sv')
  istatus = NF90_PUT_ATT(ncout, varid_ov, 'long_name', 'Overturning component of freshwater transport across section')
  istatus = NF90_PUT_ATT(ncout, varid_ov, 'comment', 'Reference salinity 35.0, transport is positive northward')
  istatus = NF90_ENDDEF(ncout)

  istatus = NF90_PUT_VAR(ncout, varid_time, dtime)
  istatus = NF90_PUT_VAR(ncout, varid_net, dnetvFW)
  istatus = NF90_PUT_VAR(ncout, varid_tot, dtotvFW)
  istatus = NF90_PUT_VAR(ncout, varid_ov, dovFW)

  IF ( istatus /= NF90_NOERR ) PRINT *, NF90_STRERROR(istatus)
  istatus = closeout(ncout)

END PROGRAM cdfvFWov
