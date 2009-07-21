PROGRAM cdfets
  !!--------------------------------------------------------------------
  !!                      ***  PROGRAM cdfets  ***
  !!
  !!  ***  Purpose: Compute Eddy Time Scale 3D field from gridT file
  !!                and the Rosby Radius of deformation.
  !!                Store the results on a 'similar' cdf file.
  !!  
  !!  ***  Method: Try to avoid 3 d arrays.
  !!               (1) Compute the BruntVaissala frequency (N2) using eosbn2
  !!               (2) Compute the Rossby Radius as the vertical integral of N, scaled 
  !!                   by |f|*pi 
  !!               (3) Compytes the buoyancy =-g x rho/rho0 and is horizontal derivative db/dx and db/dy
  !!               (4) Computes M2 = SQRT ( (db/dx)^2 + (db/dy)^2 )
  !!               (5) Computes eddy length scale = ets = N/M2
  !!               (6) Output on netcdf file ets.nc :  ets = voets ; rosby_radius = sorosrad
  !!
  !!  *** Remarks : A special care has been taken with respect to land value which have been set to
  !!                spval (-1000.) and not 0 as usual. This is because a value of 0.00 has a physical
  !!                meaning for N. On the other hand, ets is N/M2. If M2 is 0, (which is likely not very
  !!                usual), ets is set to the arbitrary value of -10., to flag these points.
  !!
  !! history :
  !!     Original :   J.M. Molines, J. Le Sommer  (Dec. 2004 ) for ORCA025
  !!                  J.M. Molines, Apr. 2005 : use of modules
  !!--------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio
  USE eos

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk                            !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: iup = 1 , idown = 2, itmp
  INTEGER, DIMENSION(2) ::  ipk, id_varout         !

  REAL(KIND=4), DIMENSION (:,:,:), ALLOCATABLE :: ztemp, zsal,zwk    !: Array to read 2 layer of data
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zn2 , &            !:  Brunt Vaissala Frequency (N2)
       &                                   zmask, e1u, e2v, e3w, ff  !: mask, metrics, and coriolis.
  REAL(KIND=4) ,DIMENSION(1)                    ::  tim
  REAL(KIND=8), DIMENSION (:,:),   ALLOCATABLE  :: buoy, dbu,dbv, zlda, M2, ets !: Double precision
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE  :: gdepw   !: depth of w level here a 1x1 array to 
  !  be in agreement with mesh_zgr.nc

  CHARACTER(LEN=256) :: cfilet ,cfileout='ets.nc'                       !:
  CHARACTER(LEN=256) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc' !:
  TYPE (variable), DIMENSION(2) :: typvar         !: structure for attribute

  INTEGER    :: ncout
  INTEGER    :: istatus

  ! constants
  REAL(KIND=4)   ::  rau0=1000., zpi, grav= 9.81, spval=-1000.

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfets  gridT '
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc must be in te current directory'
     PRINT *,' Output on ets.nc, variables voets and sorosrad'
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')

  ! define new variables for output 
  typvar(1)%name= 'voets'
  typvar(1)%units='days'
  typvar(1)%missing_value=-1000.
  typvar(1)%valid_min= 0
  typvar(1)%valid_max= 50000.
  typvar(1)%long_name='Eddy_Time_Scale'
  typvar(1)%short_name='voets'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZYX'

  typvar(2)%name= 'sorosrad'
  typvar(2)%units='m'
  typvar(2)%missing_value=-1000.
  typvar(2)%valid_min= 0.
  typvar(2)%valid_max= 50000.
  typvar(2)%long_name='Rossby_Radius'
  typvar(2)%short_name='sorosrad'
  typvar(2)%online_operation='N/A'
  typvar(2)%axis='TYX'

  ipk(1) = npk  !  3D
  ipk(2) = 1    ! 2D

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ! Allocate arrays
  ALLOCATE (ztemp(npiglo,npjglo,2), zsal(npiglo,npjglo,2), zwk(npiglo,npjglo,2) ,zmask(npiglo,npjglo))
  ALLOCATE (zn2(npiglo,npjglo), e1u(npiglo,npjglo), e2v(npiglo,npjglo) ,e3w(npiglo,npjglo))
  ALLOCATE (dbu(npiglo,npjglo), dbv(npiglo,npjglo),zlda(npiglo,npjglo) )
  ALLOCATE (buoy(npiglo,npjglo), M2(npiglo,npjglo),ets(npiglo,npjglo) ,ff(npiglo,npjglo) )
  ALLOCATE ( gdepw(npk) )

  ! create output fileset
  ncout =create(cfileout, cfilet, npiglo,npjglo,npk)
  ierr= createvar(ncout ,typvar ,2, ipk,id_varout )
  ierr= putheadervar(ncout, cfilet, npiglo, npjglo, npk)

  zpi=ACOS(-1.)

  !  2 levels of T and S are required : iup,idown (with respect to W level)
  !  Compute from bottom to top (for vertical integration)
  ztemp(:,:,idown) = getvar(cfilet, 'votemper',  npk-1 ,npiglo,npjglo )
  zsal(:,:,idown)  = getvar(cfilet, 'vosaline',  npk-1 ,npiglo,npjglo )
  zwk(:,:,idown)   = spval

  e1u(:,:) = getvar(coordhgr, 'e1u', 1,npiglo,npjglo)
  e2v(:,:) = getvar(coordhgr, 'e2v', 1,npiglo,npjglo)
  ff(:,:)  = getvar(coordhgr, 'ff', 1,npiglo,npjglo)
  gdepw(:) = getvare3(coordzgr,'gdepw',npk)

  ! eliminates zeros (which corresponds to land points where no procs were used)
  WHERE ( e1u == 0 ) 
     ff = 1.e-6
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
  ff(:,1)=ff(:,2)
  ff(1,:)=ff(2,:)

  tim=getvar1d(cfilet,'time_counter',1)
  ierr=putvar1d(ncout,tim,1,'T')

  ! at level 1 and npk, ets is not defined
  ets(:,:) = spval
  ierr = putvar(ncout, id_varout(1) ,SNGL(ets), npk, npiglo, npjglo)

  ! Set to 0 zlda
  zlda(:,:) = 0.d0
  DO jk = npk-1, 2, -1   ! from bottom to top 
     PRINT *,'level ',jk
     ! Get temperature and salinity at jk -1 (up )
     ztemp(:,:,iup)= getvar(cfilet, 'votemper',  jk-1 ,npiglo,npjglo)
     zsal(:,:,iup) = getvar(cfilet, 'vosaline',  jk-1 ,npiglo,npjglo)

     ! build tmask at level jk
     zmask(:,:)=1.
     WHERE(ztemp(:,:,idown) == 0 ) zmask = 0

     ! get depthw and e3w at level jk
     e3w(:,:)   = getvar(coordzgr, 'e3w_ps', jk,npiglo,npjglo,ldiom=.true.)
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
        zlda(:,:) = zlda(:,:) + e3w(:,:)/ff(:,:) * zn2(:,:)* zmask(:,:)

        ! Compute buoyancy at level Tk ( idown)
        buoy(:,:) = - grav * (sigma0 ( ztemp(:,:,idown),  zsal(:,:,idown),npiglo, npjglo) )  * zmask(:,:) / rau0

        ! Compute dB/dx (U point) and dB/dy (V point)
        DO jj =1 , npjglo -1
           DO ji= 1, npiglo -1
              dbu(ji,jj) = 1./e1u(ji,jj) *( buoy(ji+1,jj) - buoy(ji,jj) )
              dbv(ji,jj) = 1./e2v(ji,jj) *( buoy(ji,jj+1) - buoy(ji,jj) )
           END DO
        END DO

        ! M2 at T point ( (dB/dx)^2 + (dB/dy)^2 ) ^1/2
        DO jj=2,npjglo -1
           DO ji=2,npiglo -1
              M2(ji,jj) =  0.25*(dbu(ji,jj) + dbu(ji-1,jj)) * (dbu(ji,jj) + dbu(ji-1,jj))  &
                   + 0.25*(dbv(ji,jj) + dbv(ji,jj-1))  * (dbv(ji,jj) + dbv(ji,jj-1))
           END DO
        END DO
        M2(:,:) = SQRT( M2(:,:) )

        ! Eddy Time Scale = N / M2
        ets(:,:) = spval
        WHERE (M2 /= 0 )  
           ets =  zn2/M2/86400.   ! in seconds
        ELSEWHERE
           ets = -10.             ! flag ocean points with M2 = 0 (very few ?)
        END WHERE
        WHERE (zmask == 0 ) ets = spval

        ! write ets at level jk on the output file
        ierr = putvar(ncout, id_varout(1) ,SNGL(ets), jk, npiglo, npjglo)

        ! swap up and down, next will be read in up
        itmp = idown ; idown = iup ; iup = itmp

     END DO  ! loop to next level

     ! repeat ets at the surface and level 2 (the last computed)
     ierr = putvar(ncout, id_varout(1) ,SNGL(ets), 1,npiglo, npjglo)
     ! apply land mask (level 2) on zlda (level 1 and 2 have same mask, as there are  always at least 3 levels)

     ! Save zlda on file
     WHERE (zmask == 0 ) zlda=spval
     ierr = putvar(ncout, id_varout(2) ,SNGL(zlda), 1,npiglo, npjglo)

     istatus = closeout(ncout)

   END PROGRAM cdfets
