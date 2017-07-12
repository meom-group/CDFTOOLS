PROGRAM cdfets
  !!======================================================================
  !!                     ***  PROGRAM  cdfets  ***
  !!=====================================================================
  !!  ** Purpose : Compute Eddy Time Scale 3D field from gridT file
  !!               and the Rosby Radius of deformation.
  !!               Store the results on a 'similar' cdf file.
  !!
  !!  ** Method  : (1) Compute the BruntVaissala frequency (N2) using eosbn2
  !!               (2) Compute the Rossby Radius as the vertical integral of N,
  !!                   scaled by |f|*pi
  !!               (3) Compute the buoyancy =-g x rho/rho0 and is horizontal 
  !!                   derivative db/dx and db/dy
  !!               (4) Compute M2 = SQRT ( (db/dx)^2 + (db/dy)^2 )
  !!               (5) Compute eddy length scale = ets = N/M2
  !!               (6) Output on netcdf file ets.nc :  
  !!                   ets = voets ;  rosby_radius = sorosrad
  !!
  !!  ** Reference : Chelton et Al.  (1998)  J. Phys.Oceanogr. 28, 433-460
  !!
  !!  ** Remarks : A special care has been taken with respect to land value 
  !!               which have been set to spval (-1000.) and not 0 as usual.
  !!               This is because a value of 0.00 has a physical meaning for N.
  !!               On the other hand, ets is N/M2. If M2 is 0, (which is likely 
  !!               not very usual), ets is set to the arbitrary value of -10.,
  !!               to flag these points.
  !!
  !! History : 2.0  : 12/2004  : J.M. Molines, J. Le Sommer  : Original code
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class energy_diagnostics
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                             :: ji, jj, jk, jt            ! dummy loop index
  INTEGER(KIND=4)                             :: it                        ! time index for vvl
  INTEGER(KIND=4)                             :: ierr                      ! working integer
  INTEGER(KIND=4)                             :: narg, iargc, ijarg        ! command line 
  INTEGER(KIND=4)                             :: npiglo, npjglo            ! size of the domain
  INTEGER(KIND=4)                             :: npk, npt                  ! size of the domain
  INTEGER(KIND=4)                             :: iup = 1, idown = 2, itmp  !
  INTEGER(KIND=4)                             :: ncout                     ! ncid of output file
  INTEGER(KIND=4), DIMENSION(2)               :: ipk, id_varout            ! 

  REAL(KIND=4)                                :: rau0  = 1000.             ! density of water
  REAL(KIND=4)                                :: grav  = 9.81              ! Gravity
  REAL(KIND=4)                                :: spval = -1000.            ! special value
  REAL(KIND=4)                                :: zpi
  REAL(KIND=4)                                :: zsps                      ! missing value for salinity
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: ztemp, zsal, zwk          ! Array to read 2 layer of data
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zn2                       ! Brunt Vaissala Frequency (N2)
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zmask, ff                 ! mask coriolis.
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1u, e2v, e3w             ! metrics
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: gdepw                     ! depth of w level

  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: dtim                      ! time counter
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dbuoy, dbu, dbv           ! Double precision
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dlda, dM2, dets           ! Double precision

  CHARACTER(LEN=256)                          :: cf_tfil                   ! out file names
  CHARACTER(LEN=256)                          :: cf_out = 'ets.nc'         ! in file names
  CHARACTER(LEN=256)                          :: cldum                     ! working variable

  TYPE (variable), DIMENSION(2)               :: stypvar                   ! structure for attribute

  LOGICAL                                     :: lchk                      ! flag for missing files
  LOGICAL                                     :: lnc4 = .FALSE.            ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfets -f T-file [-o OUT-file] [-nc4] [-vvl W-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the eddy time scale, and a proxy for rossby radius. The Rossby' 
     PRINT *,'       radius is computed as the vertical integral of N2 (Brunt Vaissala '
     PRINT *,'       frequency), scaled by |f|*pi.'
     PRINT *,'       The Eddy Time Scale is the ratio N/|grad B| where N is the square root'
     PRINT *,'       of N2 and |grad B| is the module of the horizontal buoyancy gradient.'
     PRINT *,'       B is the buoyancy computed as B=-g rho/rho0.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f T-file : netcdf input file for temperature and salinity (gridT).'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file] : specifiy the name of output file instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ] : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
     PRINT *,'       [-vvl W-file] : use time varying vertical metrics. W-file holds the '
     PRINT *,'                 time-varying e3w vertical metrics.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cn_fhgr),', ',TRIM(cn_fzgr)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out),' unless -o option is used.' 
     PRINT *,'         variables : voets (days)  and sorosrad (m)'
     STOP 
  ENDIF

  ijarg = 1
  DO WHILE ( ijarg <= narg )
     CALL getarg (ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg (ijarg, cf_tfil) ; ijarg=ijarg+1
        ! options
     CASE ( '-o'   ) ; CALL getarg (ijarg, cf_out ) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4   = .TRUE.
     CASE ( '-vvl' ) ; lg_vvl = .TRUE.
        ;              CALL getarg (ijarg, cn_fe3w) ; ijarg=ijarg+1 ! change default cn_fe3w in this case
        ;              cn_ve3w = cn_ve3wvvl                         ! change default cn_ve3w
     CASE DEFAULT    ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO
  lchk = ( chkfile (cf_tfil) .OR. chkfile( cn_fhgr ) .OR. chkfile( cn_fzgr) )
  IF ( lchk )  STOP 99 ! missing file

  ! Look for missing value for salinity
  zsps = getspval(cf_tfil, cn_vosaline)

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ! Allocate arrays
  ALLOCATE (ztemp(npiglo,npjglo,2), zsal(npiglo,npjglo,2), zwk(npiglo,npjglo,2) ,zmask(npiglo,npjglo))
  ALLOCATE (zn2(npiglo,npjglo), e1u(npiglo,npjglo), e2v(npiglo,npjglo) ,e3w(npiglo,npjglo))
  ALLOCATE (dbu(npiglo,npjglo), dbv(npiglo,npjglo),dlda(npiglo,npjglo) )
  ALLOCATE (dbuoy(npiglo,npjglo), dM2(npiglo,npjglo),dets(npiglo,npjglo) ,ff(npiglo,npjglo) )
  ALLOCATE (gdepw(npk), dtim(npt) )

  CALL CreateOutput

  zpi=ACOS(-1.)

  e1u(:,:) = getvar  (cn_fhgr, cn_ve1u,  1,  npiglo, npjglo)
  e2v(:,:) = getvar  (cn_fhgr, cn_ve2v,  1,  npiglo, npjglo)
  ff(:,:)  = getvar  (cn_fhgr, cn_vff,   1,  npiglo, npjglo)
  gdepw(:) = getvare3(cn_fzgr, cn_gdepw, npk               )

  ! eliminates zeros (which corresponds to land points where no procs were used)
  WHERE ( e1u == 0 ) 
     ff  = 1.e-6
     e1u = 1
     e2v = 1
  END WHERE

  ff(:,:) = ABS(ff(:,:))* zpi
  ! need ff at T points, zwp(:,:,iup) is used as work array here.
  DO ji = 2, npiglo
     DO jj =2, npjglo
        zwk(ji,jj,iup) = 0.25 *  ( ff(ji,jj) + ff(ji,jj-1) + ff(ji-1,jj) + ff(ji-1,jj-1) )
     END DO
  END DO
  ff(:,:) = zwk(:,:,iup)
  ff(:,1) = ff(:,2)
  ff(1,:) = ff(2,:)

  DO jt = 1, npt
     IF ( lg_vvl ) THEN ; it=jt
     ELSE               ; it=1
     ENDIF
     ! at level 1. and npk, dets is not defined
     dets(:,:) = spval
     ierr = putvar(ncout, id_varout(1) ,SNGL(dets), npk, npiglo, npjglo, ktime = jt)
     !  2 levels of T and S are required : iup,idown (with respect to W level)
     !  Compute from bottom to top (for vertical integration)
     ztemp(:,:,idown) = getvar(cf_tfil, cn_votemper,  npk-1 ,npiglo,npjglo, ktime=jt )
     zsal (:,:,idown) = getvar(cf_tfil, cn_vosaline,  npk-1 ,npiglo,npjglo, ktime=jt )
     zwk  (:,:,idown) = spval

     ! Set to 0 dlda
     dlda(:,:) = 0.d0
     DO jk = npk-1, 2, -1   ! from bottom to top 
        PRINT *,'level ',jk
        ! Get temperature and salinity at jk -1 (up )
        ztemp(:,:,iup) = getvar(cf_tfil, cn_votemper,  jk-1 ,npiglo,npjglo, ktime = jt)
        zsal (:,:,iup) = getvar(cf_tfil, cn_vosaline,  jk-1 ,npiglo,npjglo, ktime = jt)

        ! build tmask at level jk
        zmask(:,:)=1.
        WHERE(ztemp(:,:,idown) == zsps ) zmask = 0

        ! get depthw and e3w at level jk
        e3w(:,:)   = getvar(cn_fe3w, cn_ve3w, jk,npiglo,npjglo, ktime=it, ldiom=.NOT.lg_vvl )
        WHERE(e3w == 0. ) e3w = 0.1     ! avoid 0's in e3w (land points anyway)

        ! zwk will hold N2 at W level
        zwk(:,:,iup) = eosbn2 ( ztemp,zsal,gdepw(jk),e3w,npiglo,npjglo, iup, idown )   ! not masked 
        WHERE( zwk(:,:,iup) < 0 ) zwk(:,:,iup) = 0.                         ! when < 0 set N2 = 0
        WHERE( zmask == 0 ) zwk(:,:,iup) = spval                            ! set to spval on land

        ! now put zn2 at T level (k )
        WHERE ( zwk(:,:,idown) == spval )  
           zn2(:,:) =  zwk(:,:,iup)
        ELSEWHERE                         
           zn2(:,:) = 0.5 * ( zwk(:,:,iup) + zwk(:,:,idown) ) 
        END WHERE

        ! Only the square root is used in this program (work for ocean points only)
        WHERE (zmask == 1 )  
           zn2=SQRT(zn2)
        ELSEWHERE            
           zn2=spval
        END WHERE

        ! integrates vertically (ff is already ABS(ff) * pi
        dlda(:,:) = dlda(:,:) + e3w(:,:)/ff(:,:) * zn2(:,:)* zmask(:,:)

        ! Compute buoyancy at level Tk ( idown)
        dbuoy(:,:) = - grav * (sigma0 ( ztemp(:,:,idown),  zsal(:,:,idown),npiglo, npjglo) )  * zmask(:,:) / rau0

        ! Compute dB/dx (U point) and dB/dy (V point)
        DO jj =1 , npjglo -1
           DO ji= 1, npiglo -1
              dbu(ji,jj) = 1./e1u(ji,jj) *( dbuoy(ji+1,jj) - dbuoy(ji,jj) )
              dbv(ji,jj) = 1./e2v(ji,jj) *( dbuoy(ji,jj+1) - dbuoy(ji,jj) )
           END DO
        END DO

        ! dM2 at T point ( (dB/dx)^2 + (dB/dy)^2 ) ^1/2
        DO jj=2,npjglo -1
           DO ji=2,npiglo -1
              dM2(ji,jj) =  0.25*(dbu(ji,jj) + dbu(ji-1,jj)) * (dbu(ji,jj) + dbu(ji-1,jj))  &
                   + 0.25*(dbv(ji,jj) + dbv(ji,jj-1))  * (dbv(ji,jj) + dbv(ji,jj-1))
           END DO
        END DO
        dM2(:,:) = SQRT( dM2(:,:) )

        ! Eddy Time Scale = N / dM2
        dets(:,:) = spval
        WHERE (dM2 /= 0 )  
           dets =  zn2/dM2/86400.   ! in seconds
        ELSEWHERE
           dets = -10.d0            ! flag ocean points with dM2 = 0 (very few ?)
        END WHERE
        WHERE (zmask == 0 ) dets = spval

        ! write dets at level jk on the output file
        ierr = putvar(ncout, id_varout(1) ,SNGL(dets), jk, npiglo, npjglo, ktime=jt)

        ! swap up and down, next will be read in up
        itmp = idown ; idown = iup ; iup = itmp

     END DO  ! loop to next level

     ! repeat dets at the surface and level 2 (the last computed)
     ierr = putvar(ncout, id_varout(1) ,SNGL(dets), 1,npiglo, npjglo, ktime=jt)

     ! apply land mask (level 2) on dlda (level 1. and 2 have same mask, as there are  always at least 3 levels)
     WHERE (zmask == 0 ) dlda=spval
     ierr = putvar(ncout, id_varout(2) ,SNGL(dlda), 1,npiglo, npjglo, ktime=jt)

  END DO  ! time loop

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
    ! define new variables for output 
    ipk(1) = npk  ! 3D
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = 'voets'
    stypvar(1)%cunits            = 'days'
    stypvar(1)%rmissing_value    = -1000.
    stypvar(1)%valid_min         = 0
    stypvar(1)%valid_max         = 50000.
    stypvar(1)%clong_name        = 'Eddy_Time_Scale'
    stypvar(1)%cshort_name       = 'voets'
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TZYX'

    ipk(2) = 1    ! 2D
    stypvar(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname             = 'sorosrad'
    stypvar(2)%cunits            = 'm'
    stypvar(2)%rmissing_value    = -1000.
    stypvar(2)%valid_min         = 0.
    stypvar(2)%valid_max         = 50000.
    stypvar(2)%clong_name        = 'Rossby_Radius'
    stypvar(2)%cshort_name       = 'sorosrad'
    stypvar(2)%conline_operation = 'N/A'
    stypvar(2)%caxis             = 'TYX'

    ! create output fileset
    ncout = create      (cf_out, cf_tfil,  npiglo, npjglo, npk       , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar , 2,      ipk,    id_varout , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_tfil,  npiglo, npjglo, npk       )

    dtim = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr = putvar1d(ncout, dtim,       npt, 'T')

  END SUBROUTINE CreateOutput


END PROGRAM cdfets
