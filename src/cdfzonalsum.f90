PROGRAM cdfzonalsum
  !!======================================================================
  !!                     ***  PROGRAM  cdfzonalsum  ***
  !!=====================================================================
  !!  ** Purpose : Compute the zonal sum of a file
  !!
  !!  ** Method  : In this program the 'zonal' sum is in fact a sum 
  !!               along the I coordinate. 
  !!
  !! History : 2.1  : 11/2005  : J.M. Molines : Original code
  !!                : 06/2007  : P. Mathiot   : adaptation for 2D files
  !!           3.0  : 05/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  !! * Local variables
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: ji, jj, jk ,jt      ! dummy loop index
  INTEGER(KIND=4)                               :: jbasin, jvar        ! dummy loop index
  INTEGER(KIND=4)                               :: ijvar               ! variable counter
  INTEGER(KIND=4)                               :: npbasins=1          ! number of subbasin
  INTEGER(KIND=4)                               :: ivar = 0            ! output variable counter
  INTEGER(KIND=4)                               :: narg, iargc         ! command line 
  INTEGER(KIND=4)                               :: ijarg, ireq         ! command line 
  INTEGER(KIND=4)                               :: npiglo, npjglo      ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt            ! size of the domain
  INTEGER(KIND=4)                               :: nvarin, nvar        ! number of input variables: all/valid
  INTEGER(KIND=4)                               :: nvarmx              ! number of output variables: all/valid
  INTEGER(KIND=4)                               :: nvaro=1             ! number of output variables: all/valid
  INTEGER(KIND=4)                               :: ncout               ! ncid of output file
  INTEGER(KIND=4)                               :: ierr                ! working integer
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipki, id_varin      ! jpbasin x nvar
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipko, id_varout     ! jpbasin x nvar
  INTEGER(KIND=4), DIMENSION(2)                 :: ijloc               ! working array for maxloc

  REAL(KIND=4)                                  :: ra = 6371229.       ! earth radius (m)
  REAL(KIND=4)                                  :: z2pi                ! 2 x 3.14...
  REAL(KIND=4), DIMENSION (:),      ALLOCATABLE :: tim                 ! time counter
  REAL(KIND=4), DIMENSION (:),      ALLOCATABLE :: gdep                ! gdept or gdepw
  REAL(KIND=4), DIMENSION (:),      ALLOCATABLE :: alpha               ! 
  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: e1, e2, gphi, zv    ! metrics, latitude, data value
  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: zdumlon             ! dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: zdumlat             ! latitude for i = north pole
  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: zmaskvar            ! variable mask
  REAL(KIND=4), DIMENSION (:,:,:),  ALLOCATABLE :: zmask               ! basin mask jpbasins x npiglo x npjglo

  REAL(KIND=8), DIMENSION (:,:),    ALLOCATABLE :: dzosum              ! jpbasins x npjglo x npk

  CHARACTER(LEN=2048)                            :: cf_in               ! input file name
  CHARACTER(LEN=2048)                            :: cf_out='zonalsum.nc' ! output file name
  CHARACTER(LEN=2048)                            :: cf_pdeg='zonalintdeg.nc' ! output file name
  CHARACTER(LEN=2048)                            :: cf_basins='none'    ! sub basin file name
  CHARACTER(LEN=10 )                            :: cv_e1, cv_e2        ! horizontal metrics variable names
  CHARACTER(LEN=10 )                            :: cv_phi              ! latitude variable name
  CHARACTER(LEN=10 )                            :: cv_msk              ! mask variable name
  CHARACTER(LEN=10 )                            :: cv_depi, cv_depo    ! depth variable name (input/output)
  CHARACTER(LEN=2048)                            :: cldum               ! dummy character variable
  CHARACTER(LEN=2048)                            :: ctyp                ! variable type on C-grid
  CHARACTER(LEN=2048), DIMENSION(:), ALLOCATABLE :: cv_namesi           ! input variable names
  CHARACTER(LEN=2048), DIMENSION(:), ALLOCATABLE :: cv_nameso           ! output variable names
  CHARACTER(LEN=80 ), DIMENSION(:), ALLOCATABLE :: cv_fix              ! name of the specified variable to process (-var option)
  CHARACTER(LEN=4  ), DIMENSION(5)              :: cbasin=(/'_glo','_atl','_inp','_ind','_pac'/) ! sub basin suffixes

  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvari            ! structure for input variables
  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvaro            ! structure for output variables

  LOGICAL                                       :: lpdep =.FALSE.      ! flag for depth sign (default dep < 0)
  LOGICAL                                       :: lpdeg =.FALSE.      ! flag for per degree normalization
  LOGICAL                                       :: lvar  =.FALSE.      ! flag for specified variables processing
  LOGICAL                                       :: ldebug   =.FALSE.   ! flag for activated debug print 
  LOGICAL                                       :: l2d   =.FALSE.      ! flag for 2D files
  LOGICAL                                       :: lchk  =.FALSE.      ! flag for missing files
  LOGICAL, DIMENSION(:),            ALLOCATABLE :: lbad                ! flag array for variable selection
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfzonalsum IN-file point_type [ BASIN-file] ...'
     PRINT *,'                  ... [-var var1,var2,..] [-pdep | --positive_depths]'
     PRINT *,'                  ... [-pdeg | --per_degree] [-debug]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the zonal sum of all the variables available in the' 
     PRINT *,'       input file. This program assume that all the variables are'
     PRINT *,'       located on the same C-grid point, specified on the command line.'
     PRINT *,'         Using -var option limits the variables to be processed.'
     PRINT *,'      '
     PRINT *,'       Zonal sum is in fact the integral value computed along the I coordinate.'
     PRINT *,'       The result is a vertical slice, in the meridional direction.'
     PRINT *,'      '
     PRINT *,'       REMARK : partial step are not handled properly (but probably '
     PRINT *,'                minor impact on results).'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       IN-file    : input netcdf file.' 
     PRINT *,'       point_type : indicate the location on C-grid (T|U|V|F|W)'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [BASIN-file] : netcdf file describing sub basins, similar to '
     PRINT *,'                      ', TRIM(cn_fbasins),'. If this name is not given '
     PRINT *,'                      as option, only the global zonal integral is computed.'
     PRINT *,'       [-var var1,var2,.. ] : Comma separated list of selected variables'
     PRINT *,'       [-pdep | --positive_depths ] : use positive depths in the output file.'
     PRINT *,'                      Default behaviour is to have negative depths.'
     PRINT *,'       [-pdeg | --per_degree ] : When using this option, the zonal integral'
     PRINT *,'                      is normalized per degree of latitude. This was formally'
     PRINT *,'                      done with cdfzonalintdeg program, which is now merged'
     PRINT *,'                      in this one.'
     PRINT *,'                      Default behaviour is not to normalize.'
     PRINT *,'       [-debug ] : add some print for debug'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fhgr),', ', TRIM(cn_fzgr),' and ', TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out),' or ', TRIM(cf_pdeg),' (-pdeg option)' 
     PRINT *,'         variables : output variable names are built with the following'
     PRINT *,'                     convention: zoixxxx_bas'
     PRINT *,'                      where zoi replace vo/so prefix of the input variable'
     PRINT *,'                      where bas is a suffix for each sub-basins (or glo)'
     PRINT *,'                      if a BASIN-file is used.'
     PRINT *,'            Units are modified by adding ''.m2'' at the end. Can be improved !'
     PRINT *,'            In addition, ''.degree-1'' is append to unit with -pdeg option.'
     STOP
  ENDIF

  ijarg = 1  ; ireq = 0
  DO WHILE ( ijarg <= narg ) 
    CALL getarg( ijarg, cldum ) ; ijarg=ijarg+1
    SELECT CASE (cldum)
    CASE ( '-pdep' , '--positive_depths' ) ; lpdep  =.TRUE.
    CASE ( '-pdeg' , '--per_degree'      ) ; lpdeg  =.TRUE.
    CASE ( '-debug'                      ) ; ldebug =.TRUE.
    CASE ( '-var'                        ) ; lvar   =.TRUE.
                                             CALL getarg( ijarg, cldum) ; ijarg=ijarg+1
                                             CALL ParseVars(cldum)
    CASE DEFAULT
      ireq=ireq+1
      SELECT CASE (ireq)
      CASE (1) ; cf_in      = cldum                 ! file name is the 1rst argument
      CASE (2) ; ctyp       = cldum                 ! point type is the 2nd
      CASE (3) ; cf_basins  = cldum                 ! sub basin file is the 3rd (optional)
                 npbasins   = 5
                 lchk       = chkfile (cf_basins)
      CASE DEFAULT 
        PRINT *,' Too many arguments ...' ; STOP
      END SELECT
    END SELECT
  END DO

  IF ( ldebug ) THEN  ! additional print for debuging
     PRINT *,' Option -pdep     ', lpdep
     PRINT *,' Option -pdeg     ', lpdeg
     PRINT *,' Option -debug    ', ldebug
     PRINT *,' Option -var      ', lvar
     PRINT *,' NPBASINS       = ', npbasins
  ENDIF

  ! check  files existence
  lchk = lchk .OR. chkfile (cn_fhgr)
  lchk = lchk .OR. chkfile (cn_fzgr)
  lchk = lchk .OR. chkfile (cn_fmsk)
  lchk = lchk .OR. chkfile (cf_in  )
  IF ( lchk ) STOP ! missing files

  ! set the metrics according to C grid point
  SELECT CASE (ctyp)
  CASE ('T', 't', 'S', 's')
     cv_e1   = cn_ve1t    ; cv_e2   = cn_ve2t
     cv_depi = cn_gdept   ; cv_depo = cn_vdeptht
     cv_phi  = cn_gphit   ; cv_msk  = 'tmask'
  CASE ('U', 'u')
     cv_e1   = cn_ve1u    ; cv_e2   = cn_ve2u
     cv_depi = cn_gdept   ; cv_depo = cn_vdepthu
     cv_phi  = cn_gphiu   ; cv_msk  = 'umask'
  CASE ('V', 'v')
     cv_e1   = cn_ve1v    ; cv_e2   = cn_ve2v
     cv_depi = cn_gdept   ; cv_depo = cn_vdepthv
     cv_phi  = cn_gphiv   ; cv_msk  = 'vmask'
  CASE ('F', 'f')
     cv_e1   = cn_ve1f    ; cv_e2   = cn_ve2f
     cv_depi = cn_gdept   ; cv_depo = cn_vdeptht
     cv_phi  = cn_gphif   ; cv_msk  = 'fmask'
  CASE ('W', 'w')
     cv_e1   = cn_ve1t    ; cv_e2   = cn_ve2t
     cv_depi = cn_gdepw   ; cv_depo = cn_vdepthw
     cv_phi  = cn_gphit   ; cv_msk  = 'tmask'
  CASE DEFAULT
     PRINT *, ' C grid:', TRIM(ctyp),' point not known!' ; STOP
  END SELECT

  nvarin  = getnvar(cf_in)   ! number of input variables
  IF ( .NOT. lvar ) THEN
    nvaro = nvarin
  ENDIF

  nvarmx = npbasins * nvaro
  IF ( ldebug ) PRINT *,' NVARMX (output) :', nvarmx

  ALLOCATE ( cv_namesi(nvarin), ipki(nvarin), id_varin (nvarin) )
  ALLOCATE ( cv_nameso(nvarmx), ipko(nvarmx), id_varout(nvarmx) )
  ALLOCATE ( stypvari(nvarin)                                   )
  ALLOCATE ( stypvaro(nvarmx)                                   )

  cv_namesi(1:nvarin) = getvarname(cf_in, nvarin, stypvari )
  ipki     (1:nvarin) = getipk    (cf_in, nvarin           )

  IF ( lvar )  THEN
     ALLOCATE ( lbad(nvarin) )
     lbad(:) = .true.  ! 
     ! tricks : in case of specified variables, set ipki to 0 all variables
     !          not choosen.
     DO ji = 1, nvaro
       DO jvar = 1, nvarin
         IF ( cv_namesi(jvar) == cv_fix(ji) ) lbad(jvar) = .false.
       END DO
     END DO
     WHERE ( lbad ) ipki=0
     DEALLOCATE (lbad )
  ENDIF

  ! buildt output filename
  nvar = 0  ! over all number of valid variables for zonal sum ( < nvarin)
  ivar = 0  ! over all variable counter ( nvar x basins)
  DO jvar = 1,nvarin
     ! skip variables such as nav_lon, nav_lat, time_counter deptht ...
     IF (ipki(jvar) == 0 ) THEN
        cv_namesi(jvar)='none'
     ELSE
        nvar           = nvar + 1 ! count for valid input variables
        id_varin(nvar) = jvar     ! use indirect adressing for those variables
        DO jbasin=1,npbasins
           ivar=ivar + 1      ! count for output variables
           cv_nameso(ivar)='zoi'//TRIM(cv_namesi(jvar)(3:))//TRIM(cbasin(jbasin) )
           ! intercept case of duplicate zonal name
           IF (cv_namesi(jvar) == 'iowaflup' ) cv_nameso(ivar)='zoiwaflio'  // TRIM(cbasin(jbasin) )
           IF (cv_namesi(jvar) == 'cfc11'    ) cv_nameso(ivar)='zoicfc11'   // TRIM(cbasin(jbasin) )
           IF (cv_namesi(jvar) == 'bombc14'  ) cv_nameso(ivar)='zoibc14'    // TRIM(cbasin(jbasin) )
           IF (cv_namesi(jvar) == 'invcfc'   ) cv_nameso(ivar)='zoiinvcfc'  // TRIM(cbasin(jbasin) )
           IF (cv_namesi(jvar) == 'invc14'   ) cv_nameso(ivar)='zoiinvc14'  // TRIM(cbasin(jbasin) )
           IF (cv_namesi(jvar) == 'qtrcfc'   ) cv_nameso(ivar)='zoiqtrcfc'  // TRIM(cbasin(jbasin) )
           IF (cv_namesi(jvar) == 'qtrc14'   ) cv_nameso(ivar)='zoiqtrc14'  // TRIM(cbasin(jbasin) )
           IF (cv_namesi(jvar) == 'qintcfc'  ) cv_nameso(ivar)='zoiqintcfc' // TRIM(cbasin(jbasin) )
           IF (cv_namesi(jvar) == 'qintc14'  ) cv_nameso(ivar)='zoiqintc14' // TRIM(cbasin(jbasin) )

           stypvaro(ivar)%cname             = cv_nameso(ivar)
           ! units can be build automatically: add .m2 at the end (not very nice ...)
           !  a special function to parse the unit and build the proper one is to be done
           !  this is tricky as many details are to be taken into account :
           !  eg :  mol/m2, kg.m-2, W/m2
           IF ( lpdeg ) THEN
             cf_out = cf_pdeg
             stypvaro(ivar)%cunits            = stypvari(jvar)%cunits//'.m2.degree-1'
             stypvaro(ivar)%clong_name        = 'Zonal_Integral_per_pegree_'//TRIM(stypvari(jvar)%clong_name)//TRIM(cbasin(jbasin) )
           ELSE
             stypvaro(ivar)%cunits            = stypvari(jvar)%cunits//'.m2'
             stypvaro(ivar)%clong_name        = 'Zonal_Integral_'//TRIM(stypvari(jvar)%clong_name)//TRIM(cbasin(jbasin) )
           ENDIF
           stypvaro(ivar)%rmissing_value    = stypvari(ivar)%rmissing_value
           stypvaro(ivar)%valid_min         = stypvari(jvar)%valid_min
           stypvaro(ivar)%valid_max         = stypvari(jvar)%valid_max
           stypvaro(ivar)%cshort_name       = stypvaro(ivar)%cname
           stypvaro(ivar)%conline_operation = '/N/A'

           IF (ipki(jvar) == 1 ) THEN
             stypvaro(ivar)%caxis           ='TY'
           ELSE
             stypvaro(ivar)%caxis           ='TZY'
           ENDIF

           ipko(ivar)=ipki(jvar)
        END DO
     ENDIF
  END DO
 
  npiglo = getdim (cf_in, cn_x)
  npjglo = getdim (cf_in, cn_y)
  npk    = getdim (cf_in, cn_z)
  npt    = getdim (cf_in, cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ! if 2D fields, npk=0, assume 1
  IF ( npk == 0 ) THEN
     npk = 1 
     l2d = .TRUE.
     PRINT *,' It is a 2D field, assume npk=1 and gdep=0'
  END IF
  ! Allocate arrays
  ALLOCATE ( zmask(npbasins,npiglo,npjglo)              )
  ALLOCATE ( zv(npiglo,npjglo), zmaskvar(npiglo,npjglo) )
  ALLOCATE ( e1(npiglo,npjglo), e2      (npiglo,npjglo) )
  ALLOCATE ( gphi(npiglo,npjglo), gdep(npk), tim(npt)   )
  ALLOCATE ( zdumlon(1,npjglo), zdumlat(1,npjglo)       )
  ALLOCATE ( dzosum(npjglo,npk), alpha(npjglo)          )
 

  ! get the metrics
  e1(:,:)   = getvar(cn_fhgr, cv_e1,  1, npiglo, npjglo) 
  e2(:,:)   = getvar(cn_fhgr, cv_e2,  1, npiglo, npjglo) 
  gphi(:,:) = getvar(cn_fhgr, cv_phi, 1, npiglo, npjglo)

  ! compute the size of the meridional mesh size in degree
  IF ( lpdeg ) THEN
     z2pi = 2.0 * ACOS( -1.)
     alpha(:) = e2(1,:) *360. / z2pi / ra
  ELSE
     alpha(:) = 1.e0
  ENDIF

  IF (l2d)  THEN 
     gdep(:) = 0
  ELSE
     gdep(:) = getvare3(cn_fzgr, cv_depi ,npk)
  ENDIF

  IF ( .NOT. lpdep ) gdep(:)   = -1.*  gdep(:)     ! helps for plotting the results

  ! Look for the i-index that go through the North Pole
  ijloc        = MAXLOC(gphi)
  zdumlat(1,:) = gphi(ijloc(1),:)
  zdumlon(:,:) = 0.              ! set the dummy longitude to 0

  ! create output fileset
  ncout = create      (cf_out, cf_in,    1,    npjglo, npk,  cdep=cv_depo                                )
  ierr  = createvar   (ncout,  stypvaro, ivar, ipko,   id_varout                                         )
  ierr  = putheadervar(ncout,  cf_in,    1,    npjglo, npk,  pnavlon=zdumlon, pnavlat=zdumlat, pdep=gdep )

  tim   = getvar1d(cf_in, cn_vtimec, npt     )
  ierr  = putvar1d(ncout, tim,       npt, 'T')

  ! reading the surface masks
  ! 1 : global ; 2 : Atlantic ; 3 : Indo-Pacif ; 4 : Indian ; 5 : Pacif
  zmask(1,:,:) = getvar(cn_fmsk, cv_msk, 1, npiglo, npjglo)
  IF ( cf_basins /= 'none' ) THEN
     zmask(2,:,:) = getvar(cf_basins, 'tmaskatl', 1, npiglo, npjglo )
     zmask(4,:,:) = getvar(cf_basins, 'tmaskind', 1, npiglo, npjglo )
     zmask(5,:,:) = getvar(cf_basins, 'tmaskpac', 1, npiglo, npjglo )
     zmask(3,:,:) = zmask(5,:,:) + zmask(4,:,:)
     ! ensure that there are no overlapping on the masks
     WHERE(zmask(3,:,:) > 0 ) zmask(3,:,:) = 1
  ENDIF

  ! main computing loop
  ivar = 0
  DO jvar = 1, nvar
     ijvar = id_varin(jvar)
     DO jt = 1,npt
        IF (MOD(jt,100)==0) PRINT *, jt,'/',npt
        DO jk = 1, ipki(ijvar)
           PRINT *,TRIM(cv_namesi(ijvar)), ' level ',jk
           ! Get variables and mask at level jk
           zv(:,:)       = getvar(cf_in,   cv_namesi(ijvar), jk ,npiglo, npjglo, ktime=jt)
           zmaskvar(:,:) = getvar(cn_fmsk, cv_msk,           jk ,npiglo, npjglo          )
           
           ! For all basins 
           DO jbasin = 1, npbasins
              dzosum(:,:) = 0.d0
              ! integrates 'zonally' (along i-coordinate)
              DO ji=1,npiglo
                 DO jj=1,npjglo
                    dzosum(jj,jk) = dzosum(jj,jk) + e1(ji,jj)*e2(ji,jj)* zmask(jbasin,ji,jj)*zmaskvar(ji,jj)*zv(ji,jj)*1.d0
                 END DO
              END DO
              dzosum(:,jk) = dzosum(:,jk) / alpha(:)  ! eventual normalization per degree
              ivar = (jvar-1)*npbasins + jbasin
              ierr = putvar (ncout, id_varout(ivar), REAL(dzosum(:,jk)), jk, 1, npjglo, ktime=jt)
           END DO  !next basin
        END DO  ! next k 
     END DO ! next time
  END DO ! next variable

  ierr = closeout(ncout)
CONTAINS
   SUBROUTINE ParseVars (cdum)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE ParseVars  ***
      !!
      !! ** Purpose :  Decode -var option from command line
      !!
      !! ** Method  :  look for , in the argument string and set the number of
      !!         variable (nvaro), allocate cv_fix array and fill it with the
      !!         decoded  names.
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: cdum

      CHARACTER(LEN=2048), DIMENSION(100) :: cl_dum  ! 100 is arbitrary
      INTEGER  :: ji
      INTEGER  :: inchar,  i1=1
      !!----------------------------------------------------------------------

      inchar= LEN(TRIM(cdum))
      ! scan the input string and look for ',' as separator
      DO ji=1,inchar
         IF ( cdum(ji:ji) == ',' ) THEN
            cl_dum(nvaro) = cdum(i1:ji-1)
            i1=ji+1
            nvaro=nvaro+1
         ENDIF
      ENDDO

      ! last name of the list does not have a ','
      cl_dum(nvaro) = cdum(i1:inchar)

      ALLOCATE ( cv_fix(nvaro) )
      IF ( ldebug) PRINT *,' SELECTED VARIABLES :'
      DO ji=1, nvaro
         cv_fix(ji) = cl_dum(ji)
         IF ( ldebug) PRINT *, "    ",TRIM(cv_fix(ji))
      ENDDO
   END SUBROUTINE ParseVars

END PROGRAM cdfzonalsum
