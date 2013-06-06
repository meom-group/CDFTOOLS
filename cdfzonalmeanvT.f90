PROGRAM cdfzonalmeanvT
  !!======================================================================
  !!                     ***  PROGRAM  cdfzonalmeanvT  ***
  !!=====================================================================
  !!  ** Purpose : Compute the mean product of zonal mean V by zonal mean
  !!               of tracer (T and S )
  !!
  !!  ** Method  : In this program the 'zonal' meanvT is in fact a meanvT 
  !!               along the I coordinate. 
/bin/bash: q: command not found
  !! History : 3.0  : 06/2013  : J.M. Molines : from cdfzonalmean and cdfvT
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id: cdfzonalmeanvT.f90 600 2012-05-09 14:08:31Z molines $
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
! ###########################
! WORK IN PROGRESS DO NOT USE
! ###########################
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: ji, jj, jk ,jt      ! dummy loop index
  INTEGER(KIND=4)                               :: jbasin, jvar        ! dummy loop index
  INTEGER(KIND=4)                               :: ijvar               ! variable counter
  INTEGER(KIND=4)                               :: npbasins=1          ! number of subbasin
  INTEGER(KIND=4)                               :: ivar = 0            ! output variable counter
  INTEGER(KIND=4)                               :: narg, iargc         ! command line 
  INTEGER(KIND=4)                               :: ijarg, ireq         ! command line 
  INTEGER(KIND=4)                               :: itag, ntag          ! arg index of 1rst tag, number of tags
  INTEGER(KIND=4)                               :: npiglo, npjglo      ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt            ! size of the domain
  INTEGER(KIND=4)                               :: nvarin, nvar        ! number of input variables: all/valid
  INTEGER(KIND=4)                               :: ncout               ! ncid of output file
  INTEGER(KIND=4)                               :: ierr, ik            ! working integers
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipki, id_varin      ! jpbasin x nvar
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipko, id_varout     ! jpbasin x nvar
  INTEGER(KIND=4), DIMENSION(2)                 :: ijloc               ! working array for maxloc

  REAL(KIND=4)                                  :: zspval=99999.       ! missing value 
  REAL(KIND=4), DIMENSION (:),      ALLOCATABLE :: tim                 ! time counter
  REAL(KIND=4), DIMENSION (:),      ALLOCATABLE :: gdep                ! gdept or gdepw
  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: e1, e2, gphi, zv    ! metrics, latitude, data value
  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: zdumlon             ! dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: zdumlat             ! latitude for i = north pole
  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: zmaskvar            ! variable mask
  REAL(KIND=4), DIMENSION (:,:,:),  ALLOCATABLE :: zmask               ! basin mask jpbasins x npiglo x npjglo

  REAL(KIND=8), DIMENSION (:,:),    ALLOCATABLE :: dzomean , darea     ! jpbasins x npjglo x npk

  CHARACTER(LEN=256)                            :: cf_in               ! input file name
  CHARACTER(LEN=256)                            :: cf_out='zonalmeanvt.nc' ! output file name
  CHARACTER(LEN=256)                            :: cf_basins='none'    ! sub basin file name
  CHARACTER(LEN=10 )                            :: cv_e1, cv_e2        ! horizontal metrics variable names
  CHARACTER(LEN=10 )                            :: cv_phi              ! latitude variable name
  CHARACTER(LEN=10 )                            :: cv_msk              ! mask variable name
  CHARACTER(LEN=10 )                            :: cv_depi, cv_depo    ! depth variable name (input/output)
  CHARACTER(LEN=256)                            :: cldum               ! dummy character variable
  CHARACTER(LEN=256)                            :: ctyp                ! variable type on C-grid
  CHARACTER(LEN=256)                            :: confcase            ! confcase name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_namesi           ! input variable names
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_nameso           ! output variable names
  CHARACTER(LEN=4  ), DIMENSION(5)              :: cbasin=(/'_glo','_atl','_inp','_ind','_pac'/) ! sub basin suffixes

  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvari            ! structure for input variables
  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvaro            ! structure for output variables

  LOGICAL                                       :: lpdep    =.FALSE.   ! flag for depth sign (default dep < 0)
  LOGICAL                                       :: lndep_in =.FALSE.   ! flag for depth sign (default dep < 0) in input file
  LOGICAL                                       :: ldebug   =.FALSE.   ! flag for activated debug print 
  LOGICAL                                       :: l2d      =.FALSE.   ! flag for 2D files
  LOGICAL                                       :: lchk     =.FALSE.   ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfzonalmeanvT [-b BASIN-file] [-pdep |--positive_depths] ... '
     PRINT *,'                   ...  [-ndep_in]   CONFIG-CASE  ''list_of_tags'' '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the mean product of zonal mean V by zonal mean of T and S.'
     PRINT *,'      '
     PRINT *,'       Zonal mean is in fact the mean value computed along the I coordinate.'
     PRINT *,'       The result is a vertical slice, in the meridional direction.'
     PRINT *,'      '
     PRINT *,'       REMARK : partial step are not handled properly (but probably '
     PRINT *,'                minor impact on results).'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       CONFIG-CASE is the config name of a given experiment (eg ORCA025-G70)'
     PRINT *,'            The program will look for gridT, gridU and gridV files for'
     PRINT *,'            this config ( grid_T, grid_U and grid_V are also accepted).'
     PRINT *,'            Additionaly, if gridS or grid_S file is found, it will be taken'
     PRINT *,'            in place of gridT for the salinity variable.'
     PRINT *,'       list_of_tags : a list of time tags that will be used for time'
     PRINT *,'            averaging. e.g. y2000m01d05 y2000m01d10 ...'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-b BASIN-file] : netcdf file describing sub basins, similar to '
     PRINT *,'                      ', TRIM(cn_fbasins),'. If this name is not given '
     PRINT *,'                      as option, only the global zonal mean is computed.'
     PRINT *,'       [-pdep | --positive_depths ] : use positive depths in the output file.'
     PRINT *,'                      Default behaviour is to have negative depths.'
     PRINT *,'       [-ndep_in ] : negative depths are used in the input file.'
     PRINT *,'                      Default behaviour is to have positive depths.'
     PRINT *,'       [-debug   ] : add some print for debug'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fhgr),', ', TRIM(cn_fzgr),' and ', TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : zovzotmean : mean product of zonal_mean(V) x zonal_mean(T)'
     PRINT *,'                     zovzotmean : mean product of zonal_mean(V) x zonal_mean(S)'
     PRINT *,'                       A suffix _bas is append to variable name oin order to'
     PRINT *,'                     indicate the basin (atl, inp, ind, pac) or glo for global'
     PRINT *,'         '
     STOP
  ENDIF

  ijarg = 1  ; ireq = 0
  DO WHILE ( ijarg <= narg ) 
    CALL getarg( ijarg, cldum ) ; ijarg=ijarg+1
    SELECT CASE (cldum)
    CASE ( '-pdep' , '--positive_depths' ) ; lpdep    =.TRUE.
    CASE ( '-ndep_in'                    ) ; lndep_in =.TRUE.
    CASE ( '-debug'                      ) ; ldebug   =.TRUE.
    CASE ( '-b'                          ) ; CALL getarg( ijarg,cf_basins ) ;  ijarg=ijarg+1 ;  npbasins   = 5
    CASE DEFAULT
      ireq=ireq+1
      SELECT CASE (ireq)
      CASE (1) ; confcase = cldum                 ! file name is the 1rst argument
      CASE DEFAULT 
          itag = ijarg -1
          ntag = narg - itag + 1
          EXIT             ! exit while loop after recording number of tags and arg number of 1rst tag
      END SELECT
    END SELECT
  END DO

  IF ( ldebug ) THEN
     PRINT *, ' CONFIG-CASE = ', TRIM(confcase)
     PRINT *, ' NTAGS       = ', ntag
     PRINT *, ' 1rst tag    = ', TRIM(cldum)
     PRINT *, '  arg pos    = ', itag
  ENDIF

  ! check  files existence
  lchk = lchk .OR. chkfile (cn_fhgr)
  lchk = lchk .OR. chkfile (cn_fzgr)
  lchk = lchk .OR. chkfile (cn_fmsk)
  lchk = lchk .OR. chkfile (cf_in  )
  IF ( npbasins /=1 ) THEN
     lchk = lchk .OR. chkfile (cf_basins  )
  ENDIF
  IF ( lchk ) STOP ! missing files

  cf_tfil = SetFileName( confcase, cldum, 'T')  ! look in first T file for dimensions

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)

! Allocation ...

! initialization of 2D time independant fields
  zmask(1,:,:) = getvar(cn_fmsk, 'tmask', 1, npiglo, npjglo)
  IF ( cf_basins /= 'none' ) THEN
     zmask(2,:,:) = getvar(cf_basins, 'tmaskatl', 1, npiglo, npjglo )
     zmask(4,:,:) = getvar(cf_basins, 'tmaskind', 1, npiglo, npjglo )
     zmask(5,:,:) = getvar(cf_basins, 'tmaskpac', 1, npiglo, npjglo )
     zmask(3,:,:) = zmask(5,:,:) + zmask(4,:,:)
     ! ensure that there are no overlapping on the masks
     WHERE(zmask(3,:,:) > 0 ) zmask(3,:,:) = 1
  ENDIF
  e1(:,:)   = getvar(cn_fhgr, cn_ve1v,  1, npiglo, npjglo) 
  e2(:,:)   = getvar(cn_fhgr, cn_ve2v,  1, npiglo, npjglo) 


! tag loop
 ijarg = itag
 nmoy  = 0  ! reset count index for time mean

 DO jtag = 1, ntag 
   CALL getarg ( ijarg, cldum ) ; ijarg = ijarg + 1
   cf_tfile = SetFileName( confcase, cldum, 'T',ld_stop=.TRUE.   )
   cf_sfile = SetFileName( confcase, cldum, 'S', ld_stop=.FALSE. )
   cf_vfile = SetFileName( confcase, cldum, 'V')
   IF ( chkfile (cf_sfil, ld_verbose=.FALSE.) ) cf_sfil = cf_tfil  ! do not complain if not found

   npt = getdim (cf_tfil,cn_t)   ! case of multiple time frames in a single file, assume identical of V file
   DO jt = 1, npt  
      nmoy = nmoy + 1
      DO jk = 1, npk 
         ! read variables 
         zsal(:,:) = getvar(cf_sfil,  cn_vosaline, jk, npiglo, npjglo, ktime=jt )
         ztem(:,:) = getvar(cf_tfil,  cn_votemper, jk, npiglo, npjglo, ktime=jt )
         zvel(:,:) = getvar(cf_tfil,  cn_vomecrty, jk, npiglo, npjglo, ktime=jt )
         ! do not read e3 metrics at level jk ( to do as in cdfzonal mean ... JMM : to be improved !
         ztmask(:,:) = getvar(cn_fmsk, 'tmask',    jk ,nnpiglo, npjglo          )
         zvmask(:,:) = getvar(cn_fmsk, 'vmask',    jk ,nnpiglo, npjglo          )

         ! put T and S at V points
         ztem(:,1:npjglo-1) = 0.5 * ( ztem(:,1:npjglo-1) + ztem(:,2:npjglo) ) * zvmask(:,1:npjglo-1)
         zsal(:,1:npjglo-1) = 0.5 * ( zsal(:,1:npjglo-1) + zsal(:,2:npjglo) ) * zvmask(:,1:npjglo-1)

         ! For all basins 
         DO jbasin = 1, npbasins
            dzovel(:) = 0.d0
            dzotem(:) = 0.d0
            dzosal(:) = 0.d0
            darea (:) = 0.d0
            ! integrates V 'zonally' (along i-coordinate)
            DO ji=1,npiglo
                  dzovel(:) = dzovel(:) + 1.d0*e1(ji,:)*e2(ji,:)* zmask(jbasin,ji,;)*zvel(ji,;)
                  dzotem(:) = dzotem(:) + 1.d0*e1(ji,:)*e2(ji,:)* zmask(jbasin,ji,;)*ztem(ji,;)
                  dzosal(:) = dzosal(:) + 1.d0*e1(ji,:)*e2(ji,:)* zmask(jbasin,ji,;)*zsal(ji,;)
                  darea (:) = darea (:) + 1.d0*e1(ji,:)*e2(ji,;)* zmask(jbasin,ji,;)*zvmask(ji,;)
            END DO

            ! compute the mean value if the darea is not 0, else assign spval
            WHERE (darea /= 0 )
               dzovel=dzovel/darea
               dzotem=dzotem/darea
               dzosal=dzosal/darea
            ELSEWHERE
               dzovel=zspval
               dzotem=zspval
               dzosal=zspval
            ENDWHERE
            ! cumulate in time
            dzovt(:,jk,jbasin) = dzovt(:,jk,jbasin) + dzovel(:)*dzotem(:)
            dzovs(:,jk,jbasin) = dzovs(:,jk,jbasin) + dzovel(:)*dzosal(:)

         END DO  !next basin
      END DO ! next level
    ENDDO   ! next time in file
  ENDDO  ! next file
  
  ! normalize before output
  dzovt(:,:,:) = dzovt(:,:,:) / nmoy
  dzovs(:,:,:) = dzovs(:,:,:) / nmoy

  ! output file 

  ierr = closeout(ncout)

CONTAINS 
    SUBROUTINE CreateOutput
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE CreateOutput  ***
      !!
      !! ** Purpose :  Set up all things required for the output file, create
      !!               the file and write the header part.
      !!
      !! ** Method  :  Use global module variables
      !!
      !!----------------------------------------------------------------------
      ipk(:) = npk

      ! create output fileset
      ncout = create      (cf_out, cf_tfil, 1, npjglo, npk    )
      ierr  = createvar   (ncout,  stypvar, 2, ipk, id_varout )
      ierr  = putheadervar(ncout,  cf_tfil, 1, npjglo, npk    )
    END SUBROUTINE CreateOutput

END PROGRAM cdfzonalmeanvT
