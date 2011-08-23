MODULE modCdfNames
  !!======================================================================
  !!                     ***  MODULE  modCdfNames  ***
  !! Declare all dimension name, variable name, attribute name as variable
  !! This will ease the generalization of CDFTOOLS
  !!=====================================================================
  !! History : 3.0  !  12/2010 ! J.M. Molines : Original code
  !! Modified: 3.0  !  08/2010 ! P.   Mathiot : Add LIM3 variables
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  PUBLIC

  ! Dimension name : cn_. [ 1 letter only ]
  CHARACTER(LEN=20) :: cn_x='x'               !: longitude, I dimension
  CHARACTER(LEN=20) :: cn_y='y'               !: latitude,  J dimension
  CHARACTER(LEN=20) :: cn_z='depth'           !: depth, z dimension
  CHARACTER(LEN=20) :: cn_t='time_counter'    !: time dimension

  ! Dimension variable
  CHARACTER(LEN=20) :: cn_vlon2d  = 'nav_lon'      !: longitude
  CHARACTER(LEN=20) :: cn_vlat2d  = 'nav_lat'      !: latitude
  CHARACTER(LEN=20) :: cn_vdeptht = 'deptht'       !: depth
  CHARACTER(LEN=20) :: cn_vdepthu = 'depthu'       !: depth
  CHARACTER(LEN=20) :: cn_vdepthv = 'depthv'       !: depth
  CHARACTER(LEN=20) :: cn_vdepthw = 'depthw'       !: depth
  CHARACTER(LEN=20) :: cn_vtimec  = 'time_counter' !: time 

  ! Attribute of a variable
  CHARACTER(LEN=20) :: cn_missing_value = 'missing_value' !: missing value (to be replaced bby _Fill_Value)

  ! Metrics
  CHARACTER(LEN=20) :: cn_ve1t='e1t', cn_ve2t='e2t'   !: e.t
  CHARACTER(LEN=20) :: cn_ve1u='e1u', cn_ve2u='e2u'   !: e.u
  CHARACTER(LEN=20) :: cn_ve1v='e1v', cn_ve2v='e2v'   !: e.v
  CHARACTER(LEN=20) :: cn_ve1f='e1f', cn_ve2f='e2f'   !: e.v
  CHARACTER(LEN=20) :: cn_ve3t='e3t', cn_ve3w='e3w'   !: e3.
  CHARACTER(LEN=20) :: cn_vff='ff'

  CHARACTER(LEN=20) :: cn_gdept='gdept', cn_gdepw='gdepw'   !: 1d dep variable
  CHARACTER(LEN=20) :: cn_hdept='hdept', cn_hdepw='hdepw'   !: 2d dep variable

  CHARACTER(LEN=20) :: cn_glamt='glamt', cn_gphit='gphit'   !:  glam gphi
  CHARACTER(LEN=20) :: cn_glamu='glamu', cn_gphiu='gphiu'   !:  glam gphi
  CHARACTER(LEN=20) :: cn_glamv='glamv', cn_gphiv='gphiv'   !:  glam gphi
  CHARACTER(LEN=20) :: cn_glamf='glamf', cn_gphif='gphif'   !:  glam gphi

  ! Generic mesh-mask file names  cn_f...
  CHARACTER(LEN=20) :: cn_fzgr='mesh_zgr.nc'
  CHARACTER(LEN=20) :: cn_fhgr='mesh_hgr.nc'
  CHARACTER(LEN=20) :: cn_fmsk='mask.nc'
  CHARACTER(LEN=20) :: cn_fcoo='coordinates.nc'
  CHARACTER(LEN=20) :: cn_fbasins='new_maskglo.nc'

  ! Variable name  : cn_v... [ starts with cn_v ]
  CHARACTER(LEN=20) :: cn_votemper='votemper' !: temperature
  CHARACTER(LEN=20) :: cn_vosaline='vosaline' !: salinity
  CHARACTER(LEN=20) :: cn_vozocrtx='vozocrtx' !: zonal velocity
  CHARACTER(LEN=20) :: cn_vomecrty='vomecrty' !: meridional velocity
  CHARACTER(LEN=20) :: cn_vomeeivv='vomeeivv' !: meridional Eddy Induced Velocity
  CHARACTER(LEN=20) :: cn_vovecrtz='vovecrtz' !: vertical velocity
  CHARACTER(LEN=20) :: cn_sossheig='sossheig' !: Sea Surface Height
  CHARACTER(LEN=20) :: cn_somxl010='somxl010' !: Mixed layer depth (density criterium)
  CHARACTER(LEN=20) :: cn_somxlt02='somxlt02' !: Mixed layer depth (temperature criterium)

  CHARACTER(LEN=20) :: cn_sohefldo='sohefldo' !: Total Heat FLux
  CHARACTER(LEN=20) :: cn_solhflup='solhflup' !: Latent Heat FLux 
  CHARACTER(LEN=20) :: cn_sosbhfup='sosbhfup' !: Sensible heat Flux
  CHARACTER(LEN=20) :: cn_solwfldo='solwfldo' !: Long Wave downward Heat Flux
  CHARACTER(LEN=20) :: cn_soshfldo='soshfldo' !: Solar Heat FLux

  CHARACTER(LEN=20) :: cn_sowaflup='sowaflup' !: Fresh Water Flux
  CHARACTER(LEN=20) :: cn_sowaflcd='sowaflcd' !: Concentration Dilution water flux
  CHARACTER(LEN=20) :: cn_sowafldp='sowafldp' !: SSS damping water Flux
  CHARACTER(LEN=20) :: cn_iowaflup='iowaflup' !: Ice Ocean Water flux ( + = freezing, - = melting)
  CHARACTER(LEN=20) :: cn_soicecov='soicecov' !: Ice cover

  ! MOC variables
  CHARACTER(LEN=20) :: cn_zomsfatl='zomsfatl' !: moc in the Atlantic
  CHARACTER(LEN=20) :: cn_zomsfglo='zomsfglo' !: moc in the Global ocean
  CHARACTER(LEN=20) :: cn_zomsfpac='zomsfpac' !: moc in the Pacific
  CHARACTER(LEN=20) :: cn_zomsfinp='zomsfinp' !: moc in the Indo-Pacific
  CHARACTER(LEN=20) :: cn_zomsfind='zomsfind' !: moc in the Indian ocean
  
  ! transport variables
  CHARACTER(LEN=20) :: cn_vozout='vozout'     !: product U x T at U point
  CHARACTER(LEN=20) :: cn_vomevt='vomevt'     !: product V x T at V point
  CHARACTER(LEN=20) :: cn_vozous='vozous'     !: product U x S at U point
  CHARACTER(LEN=20) :: cn_vomevs='vomevs'     !: product V x S at V point
  CHARACTER(LEN=20) :: cn_sozout='sozout'     !: product U x T at U point
  CHARACTER(LEN=20) :: cn_somevt='somevt'     !: product V x T at V point
  CHARACTER(LEN=20) :: cn_sozous='sozous'     !: product U x S at U point
  CHARACTER(LEN=20) :: cn_somevs='somevs'     !: product V x S at V point
  CHARACTER(LEN=20) :: cn_sozoutrp='sozoutrp' !: vertically integrated trp at U point
  CHARACTER(LEN=20) :: cn_somevtrp='somevtrp' !: vertically integrated trp at V point

  ! density, isopycnal diagnostics
  CHARACTER(LEN=20) :: cn_vosigma0='vosigma0' !: potential density refered to surface
  CHARACTER(LEN=20) :: cn_vosigmai='vosigmai' !: potential density refered to a partiular depth
  CHARACTER(LEN=20) :: cn_vodepiso='vodepiso' !: depth of isopycnal
  CHARACTER(LEN=20) :: cn_isothick='isothick' !: isopycnal tickness (from cdfsigintegr)

  ! Passive tracer variable
  CHARACTER(LEN=20) :: cn_invcfc='invcfc'     !: CFC inventory
  CHARACTER(LEN=20) :: cn_cfc11='cfc11'       !: CFC concentration
  CHARACTER(LEN=20) :: cn_pendep='pendep'     !: CFC penetration depth (from cdfpendep)
  
  ! ice variable names
  CHARACTER(LEN=20) :: cn_iicethic='iicethic' !: ice thickness
  CHARACTER(LEN=20) :: cn_ileadfra='ileadfra' !: ice concentration
  CHARACTER(LEN=20) :: cn_iicethic3='iicethic'!: ice thickness (LIM3)
  CHARACTER(LEN=20) :: cn_ileadfra3='iiceconc'!: ice concentration (LIM3)
  
  ! Bathymetry
  CHARACTER(LEN=20) :: cn_fbathymet='bathy_meter.nc' !: file Bathymetry in meters
  CHARACTER(LEN=20) :: cn_fbathylev='bathy_level.nc' !: file Bathymetry in levels

  CHARACTER(LEN=20) :: cn_bathymet='Bathymetry' !: variable Bathymetry in meters
  CHARACTER(LEN=20) :: cn_bathylev='bathy_level'!: variable Bathymetry in levels

  ! variables to be squared when performing cdfmoy
  INTEGER(KIND=4), PARAMETER :: jp_sqdvarmax=10
  INTEGER(KIND=4) :: nn_sqdvar = 4
  INTEGER(KIND=4), PRIVATE :: ji
  CHARACTER(LEN=15), DIMENSION(jp_sqdvarmax) :: cn_sqdvar = &
      & (/'vozocrtx','vomecrty','vovecrtz','sossheig',('        ', ji=jp_sqdvarmax-5,jp_sqdvarmax) /)

  ! variables eligible for 3rd moment computation when performing cdfmoy
  INTEGER(KIND=4), PARAMETER :: jp_cubvarmax=10
  INTEGER(KIND=4) :: nn_cubvar = 2
  CHARACTER(LEN=15), DIMENSION(jp_cubvarmax) :: cn_cubvar = &
      & (/'sossheig','votemper',('        ', ji=3,jp_cubvarmax) /)

! INTERFACE  
!    SUBROUTINE fdate( cldate)
!    CHARACTER(LEN=24) :: cldate
!    END SUBROUTINE fdate
! END INTERFACE

  PUBLIC :: ReadCdfNames
  PUBLIC :: PrintCdfNames

  !! NAMELIST STATEMENTS
    ! dimensions
    NAMELIST/namdim/ cn_x, cn_y, cn_z, cn_t                        ! dimensions

    ! dimension variables
    NAMELIST/namdimvar/ cn_vlon2d, cn_vlat2d, cn_vdeptht, cn_vtimec
    NAMELIST/namdimvar/ cn_vdeptht, cn_vdepthu, cn_vdepthv, cn_vdepthw
   
    ! attributes
    NAMELIST/namdimvar/ cn_missing_value

    ! metrics in coordinates, mesh_hgr
    NAMELIST/nammetrics/ cn_ve1t, cn_ve1u, cn_ve1v, cn_ve1f
    NAMELIST/nammetrics/ cn_ve2t, cn_ve2u, cn_ve2v, cn_ve2f
    NAMELIST/nammetrics/ cn_ve3t, cn_ve3w
    NAMELIST/nammetrics/ cn_vff
    NAMELIST/nammetrics/ cn_glamt, cn_glamu, cn_glamv, cn_glamf
    NAMELIST/nammetrics/ cn_gphit, cn_gphiu, cn_gphiv, cn_gphif
    !        mesh_zgr
    NAMELIST/nammetrics/ cn_gdept, cn_gdepw
    NAMELIST/nammetrics/ cn_hdept, cn_hdepw

    ! variables 
    NAMELIST/namvars/ cn_votemper, cn_vosaline
    NAMELIST/namvars/ cn_vozocrtx, cn_vomecrty, cn_vomeeivv, cn_vovecrtz
    NAMELIST/namvars/ cn_sossheig, cn_somxl010, cn_somxlt02
    NAMELIST/namvars/ cn_sohefldo, cn_solhflup, cn_sosbhfup
    NAMELIST/namvars/ cn_solwfldo, cn_soshfldo
    NAMELIST/namvars/ cn_sowaflup, cn_sowaflcd, cn_sowafldp, cn_iowaflup
    NAMELIST/namvars/ cn_zomsfatl, cn_zomsfglo, cn_zomsfpac, cn_zomsfinp, cn_zomsfind
    NAMELIST/namvars/ cn_vozout, cn_vomevt, cn_vozous, cn_vomevs
    NAMELIST/namvars/ cn_sozout, cn_somevt, cn_sozous, cn_somevs
    NAMELIST/namvars/ cn_sozoutrp, cn_somevtrp
    NAMELIST/namvars/ cn_soicecov
    NAMELIST/namvars/ cn_vosigma0, cn_vosigmai, cn_vodepiso, cn_isothick
    NAMELIST/namvars/ cn_iicethic, cn_ileadfra
    NAMELIST/namvars/ cn_invcfc,   cn_cfc11,    cn_pendep

    ! list of variable to be squared by cdfmoy
    NAMELIST/namsqdvar/ nn_sqdvar, cn_sqdvar

    ! list of variable to be cubed by cdfmoy ( option )
    NAMELIST/namcubvar/ nn_cubvar, cn_cubvar

    ! name of mesh_mask files
    NAMELIST/nammeshmask/ cn_fzgr, cn_fhgr, cn_fmsk, cn_fcoo, cn_fbasins

    ! Bathymetry
    NAMELIST/nambathy/ cn_fbathymet, cn_fbathylev, cn_bathymet, cn_bathylev
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------

CONTAINS

  SUBROUTINE ReadCdfNames ()
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ReadCdfNames  ***
    !!
    !! ** Purpose :  update the standard NetCdfName using a dedicated
    !!               namelist ( nam_cdf_names ) 
    !!
    !! ** Method  :  Look for this namelist in the following order :
    !!               1. current dir
    !!               2. HOME/CDTOOLS_cfg directory
    !!
    !!                 nam_cdf_nam can be adjusted with environment
    !!               variable NAM_CDF_NAMES
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=90) :: cl_namlist= 'nam_cdf_names'
    CHARACTER(LEN=20) :: cl_env    = 'NAM_CDF_NAMES'
    CHARACTER(LEN=90) :: cldum, cl_home
    LOGICAL           :: ll_exist
    INTEGER(KIND=4)   :: inam = 10
    !!----------------------------------------------------------------------
    CALL getenv ('HOME', cl_home) 

    ! Look for cdf namelist name
    CALL getenv (cl_env, cldum )

    IF ( cldum /= ' ' ) cl_namlist = cldum 

    ! Now look for existence of the namelist
    INQUIRE( FILE=cl_namlist, EXIST=ll_exist )

    IF ( .NOT. ll_exist ) THEN
       cldum=TRIM(cl_home)//'/CDFTOOLS_cfg/'//TRIM(cl_namlist)
       cl_namlist=cldum
       INQUIRE( FILE=cl_namlist, EXIST= ll_exist )
       IF ( .NOT. ll_exist ) THEN
          RETURN    ! assuming that there is no need to read 
                    ! a namelist for cdf names
       ENDIF
    ENDIF

    PRINT *, ' CAUTION : dim names and variable names are now set according to '
    PRINT *, ' =======   the following namelist : ', TRIM(cl_namlist) 

    OPEN(inam, FILE=cl_namlist, RECL=200)
    REWIND(inam)

    READ(inam, namdim     )
    READ(inam, namdimvar  )
    READ(inam, nammetrics )
    READ(inam, namvars    )
    READ(inam, nambathy   )
    READ(inam, namsqdvar  )
    READ(inam, nammeshmask  )
    CLOSE ( inam ) 

  END SUBROUTINE ReadCdfNames

  SUBROUTINE PrintCdfNames()
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE PrintCdfNames  ***
    !!
    !! ** Purpose :  Print a namelist like file from the actual netcdf names 
    !!
    !! ** Method  :  Use namelist facilities 
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=80) :: cl_filout='PrintCdfNames.namlist'
    CHARACTER(LEN=24) :: cl_date
    INTEGER(KIND=4)   :: iout=3
    !!----------------------------------------------------------------------
    CALL fdate(cl_date)
!   cl_date=fdate()
    OPEN(iout, FILE=cl_filout, RECL=200)
    WRITE(iout, '(a,a)' ) ' ! ', cl_date
    WRITE(iout, '(a)' ) ' ! Namelist automatically generated by PrintCdfNames '
    WRITE(iout, '(a)' ) ' ! Do not edit without changing its name ... '
    WRITE(iout, '(a)' ) ' ! ------------------------------------------'
    WRITE(iout, namdim     )
    WRITE(iout, namdimvar  )
    WRITE(iout, nammetrics )
    WRITE(iout, namvars    )
    WRITE(iout, nambathy   )
    WRITE(iout,'(a)' ) ' ! Namelist entry namsqdvar needs manual formating before'
    WRITE(iout,'(a)' ) ' ! it can be used as input : put variables names in between '' '
    WRITE(iout,'(a)' ) ' ! and separate variables by , '
    WRITE(iout, namsqdvar  )
    WRITE(iout, nammeshmask  )
    CLOSE (iout)

  END SUBROUTINE PrintCdfNames

END MODULE modCdfNames
