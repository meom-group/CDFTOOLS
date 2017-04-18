PROGRAM cdfhflx
  !!======================================================================
  !!                     ***  PROGRAM  cdfhflx  ***
  !!=====================================================================
  !!  ** Purpose : Compute the Meridional Heat Transport from the 
  !!               forcing fluxes.
  !!
  !!  ** Method  : Compute the zonaly integrated heat flux.
  !!               The program looks for the file "new_maskglo.nc". 
  !!               If it does not exist, only the calculation over all 
  !!               the whole domain is performed (this is adequate for 
  !!               a basin configuration like NATL4).
  !!               In new_maskglo.nc the masking corresponds to the global
  !!               configuration. (Global, Atlantic, Indo-Pacific, 
  !!               Indian,Pacific ocean)
  !!
  !! History : 2.1  : 07/2005  : J.M. Molines  : Original code
  !!           2.1  : 04/2006  : A.M. Treguier : adaptation to NATL4 case
  !!           2.1  : 07/2009  : R. Dussin     : Netcdf output
  !!           3.0  : 01/2011  : J.M. Molines  : Doctor norm + Lic. + generalization
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class forcing
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                             :: jbasin, ji, jj, jk, jt ! dummy loop index
  INTEGER(KIND=4)                             :: npbasins               ! number of subbasins
  INTEGER(KIND=4)                             :: ierr                   ! error status
  INTEGER(KIND=4)                             :: narg, iargc,ijarg      ! command line 
  INTEGER(KIND=4)                             :: npiglo, npjglo         ! size of the domain
  INTEGER(KIND=4)                             :: npk, npt               ! size of the domain
  INTEGER(KIND=4)                             :: ncout                  ! ncid of output file
  INTEGER(KIND=4)                             :: numout=10              ! logical unit of txt output file
  INTEGER(KIND=4)                             :: ikx=1                  ! dims of netcdf output file
  INTEGER(KIND=4), DIMENSION(:),  ALLOCATABLE :: ipk, id_varout         ! levels and varid's of output vars
  INTEGER(KIND=4), DIMENSION(2)               :: iloc                   ! used for maxloc

  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: zmask                  ! npbasins x npiglo x npjglo
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1t, e2t               ! horizontal metrics
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: gphit                  ! Latitide
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zflx                   ! fluxes read on file
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rdumlon                ! dummy longitude = 0.
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rdumlat                ! latitude for i = north pole
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: tim                    ! time counter

  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dmht                   ! cumulated heat trp
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: dhtrp                  ! MHT from fluxes

  TYPE(variable), DIMENSION(:),   ALLOCATABLE :: stypvar                ! attributes output

  CHARACTER(LEN=256)                          :: cf_tfil                ! input file
  CHARACTER(LEN=256)                          :: cf_out  ='hflx.out'    ! output txt file
  CHARACTER(LEN=256)                          :: cf_outnc='cdfhflx.nc'  ! output nc file
  CHARACTER(LEN=256)                          :: cldum                  ! working variable

  LOGICAL                                     :: lglo = .FALSE.         ! global or subbasin computation
  LOGICAL                                     :: lchk = .FALSE.         ! missing file flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfhflx  -f T-file [-o OUTNC-file ] [-ot OUTTXT-file] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the Meridional Heat Transport (MHT) from surface heat fluxes,' 
     PRINT *,'       in function of the latitude.'
     PRINT *,'       If a sub-basin file is available, MHT is computed for each sub-basin.'
     PRINT *,'       Note that the latitude is in fact a line of constant J coordinate, not'
     PRINT *,'       a true parallel, if the model grid is distorted as in the northern most'
     PRINT *,'       part of ORCA configurations.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f T-file : a file with heat fluxes (gridT). '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUTNC-file ]: specify the name of the netcdf output file, instead of'
     PRINT *,'                    ', TRIM(cf_outnc) 
     PRINT *,'       [-ot OUTTXT-file ]: specify the name of the text output file, instead of'
     PRINT *,'                    ', TRIM(cf_out) 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       Files ', TRIM(cn_fhgr),', ',TRIM(cn_fbasins),' and ',TRIM(cn_fmsk),'.' 
     PRINT *,'       If ',TRIM(cn_fbasins),' is not available, only global MHT is computed.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       ASCII file  : ', TRIM(cf_out  )
     PRINT *,'       netcdf file : ', TRIM(cf_outnc) 
     PRINT *,'         variables : hflx_glo, [hflx_atl, hflx_inp, hflx_pac, hflx_ind]'
     STOP
  ENDIF

  ijarg = 1
  DO WHILE ( ijarg <= narg )
     CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum)
     CASE ( '-f' ) ; CALL getarg (ijarg, cf_tfil ) ; ijarg=ijarg+1
        ! options
     CASE ( '-o' ) ; CALL getarg (ijarg, cf_outnc) ; ijarg=ijarg+1
     CASE ( '-ot') ; CALL getarg (ijarg, cf_out  ) ; ijarg=ijarg+1
     CASE DEFAULT  ; PRINT *,'ERROR : ',TRIM(cldum), ': unknown option.' ; STOP
     END SELECT
  ENDDO

  lchk = chkfile(cn_fhgr)
  lchk = chkfile(cn_fmsk) .OR. lchk
  lchk = chkfile(cf_tfil) .OR. lchk
  IF ( lchk ) STOP ! missing file

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  !  Detects newmaskglo file 
  lglo = .NOT. ( chkfile(cn_fbasins) )

  IF (lglo) THEN ; npbasins = 5
  ELSE           ; npbasins = 1
  ENDIF

  ! Allocate arrays
  ALLOCATE ( zmask(npbasins,npiglo,npjglo) )
  ALLOCATE ( zflx(npiglo,npjglo) )
  ALLOCATE ( e1t(npiglo,npjglo), e2t(npiglo,npjglo), gphit(npiglo,npjglo) )
  ALLOCATE ( dhtrp (npbasins,npjglo) )
  ALLOCATE ( dmht(npbasins, npjglo) )
  ALLOCATE ( rdumlon(1,npjglo), rdumlat(1,npjglo) )
  ALLOCATE ( tim(npt) )

  ALLOCATE (stypvar(npbasins), ipk(npbasins), id_varout(npbasins))

  e1t(  :,:) = getvar(cn_fhgr, cn_ve1t,  1, npiglo, npjglo) 
  e2t(  :,:) = getvar(cn_fhgr, cn_ve2t,  1, npiglo, npjglo) 
  gphit(:,:) = getvar(cn_fhgr, cn_gphit, 1, npiglo, npjglo)

  iloc         = MAXLOC(gphit)
  rdumlat(1,:) = gphit(iloc(1),:)
  rdumlon(:,:) = 0.   ! set the dummy longitude to 0

  CALL CreateOutput

  OPEN(numout, FILE=cf_out, FORM='FORMATTED', RECL=256)  ! to avoid wrapped line with ifort
  WRITE(numout,*)'! Zonal heat transport (integrated from surface fluxes) (in Pw)'

  ! reading the masks
  ! 1 : global ; 2 : Atlantic ; 3 : Indo-Pacif ; 4 : Indian ; 5 : Pacif
  zmask(1,:,:)= getvar(cn_fmsk, cn_vmask, 1, npiglo, npjglo)

  IF (lglo) THEN
     zmask(2,:,:) = getvar(cn_fbasins, cn_tmaskatl, 1, npiglo, npjglo)
     zmask(4,:,:) = getvar(cn_fbasins, cn_tmaskind, 1, npiglo, npjglo)
     zmask(5,:,:) = getvar(cn_fbasins, cn_tmaskpac, 1, npiglo, npjglo)
     zmask(3,:,:) = zmask(5,:,:) + zmask(4,:,:)
     ! ensure that there are no overlapping on the masks
     WHERE(zmask(3,:,:) > 0 ) zmask(3,:,:) = 1
     ! change global mask for GLOBAL periodic condition
     zmask(1,1,:) = 0.
     zmask(1,npiglo,:) = 0.
  ENDIF

  DO jt = 1, npt
     ! initialize dmht
     dmht(:,:)  = 0.d0
     dhtrp(:,:) = 0.d0
     WRITE(numout,*)' TIME =', jt, tim(jt)/86400.,' days'

     ! Get fluxes
     zflx(:,:)= getvar(cf_tfil, cn_sohefldo, 1, npiglo, npjglo, ktime=jt)

     ! integrates 'zonally' (along i-coordinate)
     DO ji=1,npiglo
        ! For all basins 
        DO jbasin = 1, npbasins
           dmht(jbasin,:) = dmht(jbasin,:) + e1t(ji,:)*e2t(ji,:)* zmask(jbasin,ji,:)*zflx(ji,:)*1.d0
        END DO
     END DO

     ! cumulates transport from north to south
     DO jj=npjglo-1,1,-1
        dhtrp(:,jj) = dhtrp(:,jj+1) - dmht(:,jj)
     END DO

     ! transform to peta watt
     dhtrp(:,:) = dhtrp(:,:) / 1.d15    

     IF (lglo) THEN
        WRITE(numout,*)'! J        Global          Atlantic         INDO-PACIF    INDIAN  PACIF '
        DO jj=npjglo, 1, -1
           WRITE(numout,9000) jj, &
                rdumlat(1,jj), dhtrp(1,jj),      dhtrp(2,jj),        dhtrp(3,jj), dhtrp(4,jj), dhtrp(5,jj)
        ENDDO
     ELSE
        WRITE(numout,*)'! J        Global   '
        DO jj=npjglo, 1, -1
           WRITE(numout,9000) jj, rdumlat(1,jj), dhtrp(1,jj)  
        ENDDO
     ENDIF

9000 FORMAT(I4,5(1x,f9.3,1x,f8.4))

     DO jj=1, npbasins
        ierr = putvar(ncout, id_varout(jj), REAL(dhtrp(jj,:)), ipk(jj), ikx, npjglo, ktime=jt )
     END DO
  END DO   ! time loop

  ierr = closeout(ncout)
  CLOSE(numout)
CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ! define new variables for output 
    ipk(:)                    = 1
    stypvar%cunits            = 'PW'
    stypvar%rmissing_value    = 99999.
    stypvar%valid_min         = -1000.
    stypvar%valid_max         = 1000.
    stypvar%scale_factor      = 1.
    stypvar%add_offset        = 0.
    stypvar%savelog10         = 0.
    stypvar%cunits            = 'PW' 
    stypvar%conline_operation = 'N/A'
    stypvar%caxis             = 'T'

    stypvar(1)%cname          = 'hflx_glo'
    stypvar(1)%clong_name     = 'Heat_Fluxes_Global'
    stypvar(1)%cshort_name    = 'hflx_glo'

    IF (lglo) THEN
       stypvar(2)%cname       = 'hflx_atl'             ; stypvar(3)%cname       = 'hflx_inp'
       stypvar(2)%clong_name  = 'Heat_Fluxes_Atlantic' ; stypvar(3)%clong_name  = 'Heat_Fluxes_Indo-Pacific'
       stypvar(2)%cshort_name = 'hflx_atl'             ; stypvar(3)%cshort_name = 'hflx_inp'

       stypvar(4)%cname       = 'hflx_ind'             ; stypvar(5)%cname       = 'hflx_pac'
       stypvar(4)%clong_name  = 'Heat_Fluxes_Indian'   ; stypvar(5)%clong_name  = 'Heat_Fluxes_Pacific'
       stypvar(4)%cshort_name = 'hflx_ind'             ; stypvar(5)%cshort_name = 'hflx_pac'
    ENDIF
    ! create output fileset
    ncout = create      (cf_outnc, 'none',  ikx,      npjglo, npk                                  )
    ierr  = createvar   (ncout,    stypvar, npbasins, ipk,    id_varout                            )
    ierr  = putheadervar(ncout,    cf_tfil, ikx,      npjglo, npk, pnavlon=rdumlon, pnavlat=rdumlat)

    tim  = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr = putvar1d(ncout,   tim,       npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfhflx
