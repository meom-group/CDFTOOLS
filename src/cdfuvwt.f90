PROGRAM cdfuvwt
  !!======================================================================
  !!                     ***  PROGRAM  cdfuvwt  ***
  !!=====================================================================
  !!  ** Purpose : Compute mean values of some quantities, required for
  !!               other cdftools ( cdfbci, cdfbti and cdfnrjcomp).
  !!               At U point : ubar, u2bar
  !!               At V point : vbar, v2bar
  !!               At W point : wbar
  !!               AT T point : tbar, t2bar, uvbar, utbar, vtbar, wtbar
  !!
  !!  ** Method  : take care of double precision on product
  !!
  !! History : 2.1  : 02/2008  : A. Melet     : Original code
  !!           3.0  : 06/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils    ! SetFileName function
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class second_order_moments
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER                    :: jp_var = 11
  INTEGER(KIND=4)                               :: ji, jj, jk, jt, jtt
  INTEGER(KIND=4)                               :: ntframe
  INTEGER(KIND=4)                               :: npiglo, npjglo
  INTEGER(KIND=4)                               :: npk, npt, ntags
  INTEGER(KIND=4)                               :: iimin, iimax, ijmin, ijmax
  INTEGER(KIND=4)                               :: iup=1, idwn=2
  INTEGER(KIND=4)                               :: narg, iargc, ijarg
  INTEGER(KIND=4)                               :: ncout
  INTEGER(KIND=4)                               :: ierr
  INTEGER(KIND=4), DIMENSION(jp_var)            :: ipk, id_varout         ! 

  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: u2d, v2d, t2d
  REAL(KIND=4), DIMENSION(:,:,:),   ALLOCATABLE :: w2d

  REAL(KIND=8)                                  :: dcoef
  REAL(KIND=8)                                  :: dtotal_time
  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dtim
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dtabu, dtabv, dtabu2, dtabv2, dtabuv
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dtabw, dtabt, dtabut, dtabvt, dtabt2
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dtabwt

  CHARACTER(LEN=256)                            :: cf_ufil, cf_vfil
  CHARACTER(LEN=256)                            :: cf_wfil, cf_tfil
  CHARACTER(LEN=256)                            :: cf_out='moyuvwt.nc'
  CHARACTER(LEN=256)                            :: cldum
  CHARACTER(LEN=256)                            :: config , ctag
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: ctag_lst

  TYPE (variable), DIMENSION(jp_var)            :: stypvar         ! structure for attibutes

  LOGICAL                                       :: lnc4 = .FALSE.  ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  !!
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfuvwt -c CONFCASE -l LST-tags [-w imin imax jmin jmax] ...'
     PRINT *,'               ... [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the time-mean fields for velocity components (u,v,w) and' 
     PRINT *,'       temperature (t), as well as second order moments ( u2, v2, t2, uv, ut,'
     PRINT *,'       vt, wt).'
     PRINT *,'       These fields are required in other cdftools which computes either '
     PRINT *,'       barotropic (cdfnrj_bti) or baroclinic (cdfnrj_bci) instabilities, and a'
     PRINT *,'       global energy balance (cdfnrj_components)'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -c CONFCASE : the root name for the data files. Grid files are assumed '
     PRINT *,'                to be gridT, gridU, gridV, gridW. (grid_T, grid_U, grid_V and'
     PRINT *,'                grid_W are also supported).'
     PRINT *,'       -l LST-tags : set the list of time tags corresponding to the time serie'
     PRINT *,'                whose mean is being computed.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-w imin imax jmin jmax ] : model window limiting the area where the'
     PRINT *,'                time-means will be computed.'
     PRINT *,'       [-o OUT-file] : specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4]   : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'             This option is effective only if cdftools are compiled with'
     PRINT *,'             a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables :  There are 11 variables produced by this program.'
     PRINT *,'                 tbar, t2bar : mean t (Kelvin) and mean t^2 (K^2)   [T-point]'
     PRINT *,'                 ubar, u2bar : mean u (m/s) and mean u^2 (m2/s2)    [U-point]'
     PRINT *,'                 vbar, v2bar : mean v (m/s) and mean v^2 (m2/s2)    [V-point]'
     PRINT *,'                 wbar        : mean w (m/s)                         [W-point]'
     PRINT *,'                 uvbar       : mean product u . v (m2/s2)           [T-point]'       
     PRINT *,'                 utbar, vtbar, wtbar : mean product [uvw].t (m/s.K) [T-point]'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      cdfnrj_bti, cdfnrj_bci, cdfnrj_components and cdfnrj_transfert.' 
     PRINT *,'      '
     STOP 
  ENDIF

  iimin=0 ; ijmin=0
  iimax=0 ; ijmax=0

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum) ; ijarg = ijarg +1
     SELECT CASE (cldum)
     CASE ( '-c'  ) ; CALL getarg(ijarg, config) ; ijarg=ijarg+1
     CASE ( '-l'  ) ; CALL GetTagList
        ! options
     CASE ( '-w'  ) ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) iimin
        ;             CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) iimax
        ;             CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) ijmin
        ;             CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) ijmax
     CASE ( '-o'  ) ; CALL getarg(ijarg, cf_out) ; ijarg=ijarg+1
     CASE ( '-nc4') ; lnc4 = .TRUE.
     CASE DEFAULT   ; PRINT *,' ERROR : ',TRIM(cldum), ' : unknown option.' ; STOP 99
     END SELECT
  END DO

  cf_ufil = SetFileName( config, ctag_lst(1), 'U')

  ! assume all input files have same spatial size
  npiglo = getdim (cf_ufil, cn_x )
  npjglo = getdim (cf_ufil, cn_y )
  npk    = getdim (cf_ufil, cn_z )

  ! modify sizes with respect to working window  area
  IF (iimin /= 0 ) THEN ; npiglo=iimax -iimin + 1 ; ELSE ; iimin=1 ; iimax=npiglo ; ENDIF
  IF (ijmin /= 0 ) THEN ; npjglo=ijmax -ijmin + 1 ; ELSE ; ijmin=1 ; ijmax=npjglo ; ENDIF

  PRINT *, ' npiglo = ', npiglo
  PRINT *, ' npjglo = ', npjglo
  PRINT *, ' npk    = ', npk

  ! Allocate the memory
  ALLOCATE ( u2d(npiglo,npjglo), v2d(npiglo,npjglo)   )
  ALLOCATE ( t2d(npiglo,npjglo), w2d(npiglo,npjglo,2) )

  ALLOCATE ( dtabu(npiglo,npjglo), dtabu2(npiglo,npjglo), dtabut(npiglo,npjglo) )
  ALLOCATE ( dtabv(npiglo,npjglo), dtabv2(npiglo,npjglo), dtabvt(npiglo,npjglo) )
  ALLOCATE ( dtabt(npiglo,npjglo), dtabt2(npiglo,npjglo)                        )
  ALLOCATE ( dtabw(npiglo,npjglo),                        dtabwt(npiglo,npjglo) )
  ALLOCATE (                                              dtabuv(npiglo,npjglo) )

  CALL CreateOutput

  DO jk=1, npk-1   ! level npk is masked for T U V and is 0 for W ( bottom ) !
     PRINT *,'            level ',jk
     dtotal_time  = 0.d0 ;  ntframe=0 
     dtabu(:,:)   = 0.d0 ; dtabv(:,:)  = 0.d0 ; dtabuv(:,:) = 0.d0 
     dtabu2(:,:)  = 0.d0 ; dtabv2(:,:) = 0.d0 ; dtabt(:,:)  = 0.d0 
     dtabw(:,:)   = 0.d0 ; dtabut(:,:) = 0.d0 ; dtabvt(:,:) = 0.d0 
     dtabt2(:,:)  = 0.d0 ; dtabwt(:,:) = 0.d0

     DO jt= 1, ntags
        ctag    = ctag_lst(jt)
        cf_ufil = SetFileName( config, ctag_lst(jt), 'U')
        cf_vfil = SetFileName( config, ctag_lst(jt), 'V')
        cf_wfil = SetFileName( config, ctag_lst(jt), 'W')
        cf_tfil = SetFileName( config, ctag_lst(jt), 'T')
        
        IF ( jk == 1 ) THEN
           npt = getdim(cf_ufil, cn_t)
           ALLOCATE ( dtim(npt) )
           dtim=getvar1d(cf_ufil, cn_vtimec, npt)
           dtotal_time = dtotal_time + SUM(dtim)
           DEALLOCATE ( dtim )
        ENDIF

        DO jtt = 1, npt
           ntframe  = ntframe+1
           u2d(:,:)      = getvar(cf_ufil, cn_vozocrtx, jk,   npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=jtt )
           v2d(:,:)      = getvar(cf_vfil, cn_vomecrty, jk,   npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=jtt )
           w2d(:,:,iup)  = getvar(cf_wfil, cn_vovecrtz, jk,   npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=jtt )
           w2d(:,:,idwn) = getvar(cf_wfil, cn_vovecrtz, jk+1, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=jtt )
           t2d(:,:)      = getvar(cf_tfil, cn_votemper, jk,   npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=jtt )
           WHERE ( t2d /= 0. ) t2d = t2d + 273.15   ! from C to K

           dtabu(:,:)  = dtabu(:,:)  + u2d(:,:)
           dtabu2(:,:) = dtabu2(:,:) + u2d(:,:) * u2d(:,:) * 1.d0
           dtabv(:,:)  = dtabv(:,:)  + v2d(:,:)
           dtabv2(:,:) = dtabv2(:,:) + v2d(:,:) * v2d(:,:) * 1.d0
           dtabw(:,:)  = dtabw(:,:)  + w2d(:,:,iup)
           dtabt(:,:)  = dtabt(:,:)  + t2d(:,:)
           dtabt2(:,:) = dtabt2(:,:) + t2d(:,:) * t2d(:,:) * 1.d0

           DO jj = npjglo, 2 , -1
              DO ji = npiglo, 2 , -1
                 ! put u, v on T point ( note the loops starting from the end for using u2d and v2d as tmp array)
                 u2d(ji,jj) = 0.5 * ( u2d(ji,jj) + u2d(ji-1,jj  ) )
                 v2d(ji,jj) = 0.5 * ( v2d(ji,jj) + v2d(ji,  jj-1) )
              END DO
           END DO
           u2d(1,:) = 0. ; u2d(:,1) = 0.
           v2d(1,:) = 0. ; v2d(:,1) = 0.
           w2d(:,:,iup) = 0.5 * ( w2d(:,:,iup) + w2d(:,:,idwn) )  ! W at  T point

           dtabuv(:,:) = dtabuv(:,:) + u2d(:,:)     * v2d(:,:) * 1.d0
           dtabut(:,:) = dtabut(:,:) + u2d(:,:)     * t2d(:,:) * 1.d0
           dtabvt(:,:) = dtabvt(:,:) + v2d(:,:)     * t2d(:,:) * 1.d0
           dtabwt(:,:) = dtabwt(:,:) + w2d(:,:,iup) * t2d(:,:) * 1.d0

        END DO  ! jtt
     END DO     ! tags

     dcoef = 1.d0 / ntframe

     ! save on file
     ierr = putvar(ncout, id_varout( 1), REAL(dtabu * dcoef), jk, npiglo, npjglo )
     ierr = putvar(ncout, id_varout( 2), REAL(dtabv * dcoef), jk, npiglo, npjglo )
     ierr = putvar(ncout, id_varout( 3), REAL(dtabu2* dcoef), jk, npiglo, npjglo )
     ierr = putvar(ncout, id_varout( 4), REAL(dtabv2* dcoef), jk, npiglo, npjglo )
     ierr = putvar(ncout, id_varout( 5), REAL(dtabuv* dcoef), jk, npiglo, npjglo )
     ierr = putvar(ncout, id_varout( 6), REAL(dtabw * dcoef), jk, npiglo, npjglo )
     ierr = putvar(ncout, id_varout( 7), REAL(dtabt * dcoef), jk, npiglo, npjglo )
     ierr = putvar(ncout, id_varout( 8), REAL(dtabut* dcoef), jk, npiglo, npjglo )
     ierr = putvar(ncout, id_varout( 9), REAL(dtabvt* dcoef), jk, npiglo, npjglo )
     ierr = putvar(ncout, id_varout(10), REAL(dtabt2* dcoef), jk, npiglo, npjglo )
     ierr = putvar(ncout, id_varout(11), REAL(dtabwt* dcoef), jk, npiglo, npjglo )

  END DO   ! loop on level

  ! fill up empty last level
  dtabu = 0.d0  ! reset this dummy array to 0 for npk output
  DO jk= 1, jp_var
     ierr  = putvar(ncout, id_varout(jk), REAL(dtabu), npk, npiglo, npjglo )
  END DO

  ierr = putvar1d(ncout, (/dtotal_time*dcoef/), 1, 'T')
  ierr = closeout(ncout)

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
    INTEGER(KIND=4) :: jv
    !!----------------------------------------------------------------------
  DO jv=1,jp_var
     stypvar(jv)%ichunk   = (/npiglo,MAX(1,npjglo/30),1,1 /)
  ENDDO
  ipk(:) = npk  

  stypvar( 1)%cname       = 'ubar'
  stypvar( 1)%clong_name  = 'temporal mean of u on U point'
  stypvar( 1)%cshort_name = 'ubar'
  stypvar( 1)%cunits      = 'm/s'

  stypvar( 2)%cname       = 'vbar'
  stypvar( 2)%clong_name  = 'temporal mean of v on V point'
  stypvar( 2)%cshort_name = 'vbar'
  stypvar( 2)%cunits      = 'm/s'

  stypvar( 3)%cname       = 'u2bar'
  stypvar( 3)%clong_name  = 'temporal mean of u * u on U point'
  stypvar( 3)%cshort_name = 'u2bar'
  stypvar( 3)%cunits      = 'm2/s2'

  stypvar( 4)%cname       = 'v2bar'
  stypvar( 4)%clong_name  = 'temporal mean of v * v on V point'
  stypvar( 4)%cshort_name = 'v2bar'
  stypvar( 4)%cunits      = 'm2/s2'

  stypvar( 5)%cname       = 'uvbar'
  stypvar( 5)%clong_name  = 'temporal mean of u * v on T point'
  stypvar( 5)%cshort_name = 'uvbar'
  stypvar( 5)%cunits      = 'm2/s2'

  stypvar( 6)%cname       = 'wbar'
  stypvar( 6)%clong_name  = 'temporal mean of w on W point'
  stypvar( 6)%cshort_name = 'wbar'
  stypvar( 6)%cunits      = 'm/s'

  stypvar( 7)%cname       = 'tbar'
  stypvar( 7)%clong_name  = 'temporal mean of T on T point in K'
  stypvar( 7)%cshort_name = 'tbar'
  stypvar( 7)%cunits      = 'K'

  stypvar( 8)%cname       = 'utbar'
  stypvar( 8)%clong_name  = 'temporal mean of u * T (in K) on T point'
  stypvar( 8)%cshort_name = 'utbar'
  stypvar( 8)%cunits      = 'm/s.K'

  stypvar( 9)%cname       = 'vtbar'
  stypvar( 9)%clong_name  = 'temporal mean of v * T (in K) on T point'
  stypvar( 9)%cshort_name = 'vtbar'
  stypvar( 9)%cunits      = 'm/s.K'

  stypvar(10)%cname       = 't2bar'
  stypvar(10)%clong_name  = 'temporal mean of T * T on T point in K^2'
  stypvar(10)%cshort_name = 't2bar'
  stypvar(10)%cunits      = 'K2'

  stypvar(11)%cname       = 'wtbar'
  stypvar(11)%clong_name  = 'temporal mean of w * T (in K) on T point'
  stypvar(11)%cshort_name = 'wtbar'
  stypvar(11)%cunits      = 'm/s.K'

  stypvar%rmissing_value    = 0.
  stypvar%valid_min         = -1000.
  stypvar%valid_max         = 1000.
  stypvar%conline_operation = 'N/A'
  stypvar%caxis             = 'TYX'

  ! create output fileset
  ncout = create      (cf_out, cf_ufil, npiglo, npjglo, npk       , ld_nc4=lnc4 )
  ierr  = createvar   (ncout,  stypvar, jp_var, ipk,    id_varout , ld_nc4=lnc4 )
  ierr  = putheadervar(ncout,  cf_ufil, npiglo, npjglo, npk       )

  END SUBROUTINE CreateOutput

  SUBROUTINE GetTagList
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetTagList  ***
    !!
    !! ** Purpose :  Set up a tag list given on the command line as 
    !!               blank separated list
    !!
    !! ** Method  :  Scan the command line until a '-' is found
    !!----------------------------------------------------------------------
    INTEGER (KIND=4)  :: ji
    INTEGER (KIND=4)  :: icur
    !!----------------------------------------------------------------------
    !!
    ntags=0
    ! need to read a list of file ( number unknow ) 
    ! loop on argument till a '-' is found as first char
    icur=ijarg                          ! save current position of argument number
    DO ji = icur, narg                  ! scan arguments till - found
       CALL getarg ( ji, cldum )
       IF ( cldum(1:1) /= '-' ) THEN ; ntags = ntags+1
       ELSE                          ; EXIT
       ENDIF
    ENDDO
    ALLOCATE (ctag_lst(ntags) )
    DO ji = icur, icur + ntags -1
       CALL getarg(ji, ctag_lst( ji -icur +1 ) ) ; ijarg=ijarg+1
    END DO
  END SUBROUTINE GetTagList

END PROGRAM cdfuvwt

