PROGRAM cdfpvort
  !!-------------------------------------------------------------------
  !!                 ***  PROGRAM cdfpvort ***
  !!
  !!  **  Purpose: Compute the Ertel Potential vorticity 
  !!  
  !!  **  Method: Try to avoid 3 d arrays : work with 2 levels a a time
  !! Formula :
  !!   Qpot = drho/dz * ( f + xsi ) = Qstr + Qrel
  !!   * f is the Coriolis factor, computed from the latitudes of the T-grid :
  !!       f(i,j) = 2 * omega * sin ( phit(i,j) * pi / 180 )
  !!
  !!   * xsi is the relative vorticity (vertical component of the velocity curl),
  !!     computed from the relative vorticity of the F-points interpolated at
  !!     the T-points :
  !!     xsif(i,j) = ( ue(i,j) - ue(i,j+1) - ve(i,j) + ve(i+1,j) ) / areaf(i,j)
  !!     with : ue(i,j) = U(i,j) * e1u(i,j)
  !!            ve(i,j) = V(i,j) * e2v(i,j)
  !!            areaf(i,j) = e1f(i,j) * e2f(i,j)
  !!     xsi(i,j) = ( xsif(i-1,j-1) + xsif(i-1,j) + xsif(i,j-1) + xsif(i,j) ) / 4
  !!              = (  ue(i-1,j-1) + ue(i,j-1) - ue(i-1,j+1) - ue(i,j+1)
  !!                 - ve(i-1,j-1) - ve(i-1,j) + ve(i+1,j-1) + ve(i+1,j) )
  !!                / 4 / areat(i,j)
  !!     with : areat(i,j) = e1t(i,j) * e2t(i,j)
  !!
  !!   units : U, V in m.s-1
  !!           e1u, e2v, e1f, e2f in m
  !!           f, xsi in s-1
  !!           Qpot, Qrel, Qstr in kg.m-4.s-1
  !!
  !!
  !!           The brunt-vaisala frequency is computed using the
  !!           polynomial expression of McDougall (1987):
  !!           N^2 = grav * beta * ( alpha/beta*dk[ t ] - dk[ s ] )/e3w
  !!           N2 is then insterpolated at T levels
  !!
  !! history:
  !!     Original : A.M Treguier december 2005
  !!--------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !!
  !! * Modules used
  USE cdfio
  USE eos

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk, jj, ji,jt                          !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: 
  INTEGER   :: npiglo,npjglo, npk ,npt             !: size of the domain
  INTEGER   :: iup = 1 , idown = 2, itmp
  INTEGER, DIMENSION(3) ::  ipk, id_varout
  REAL(kind=8) , DIMENSION(:,:), ALLOCATABLE  :: e2v, e1u, e1t, e2t, gphit 
  REAL(kind=4) , DIMENSION(:,:), ALLOCATABLE  :: un, vn, rotn, zareat, z2fcor , stretch
  REAL(KIND=4) , DIMENSION (:,:,:), ALLOCATABLE :: ztemp, zsal,zwk    !: Array to read 2 layer of data
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE ::       &
       zn2 , &              !:  Brunt Vaissala Frequency (N2)
       tmask,  e3w
  REAL(KIND=4) , DIMENSION (:), ALLOCATABLE ::     gdepw  
  REAL(KIND=4),DIMENSION(:) ,ALLOCATABLE    ::  tim

  CHARACTER(LEN=80) :: cfilet , cfileu, cfilev, cfileout='vpot.nc'   !:
  CHARACTER(LEN=80) :: coordzgr='mesh_zgr.nc' !:
  CHARACTER(LEN=80) :: coord   ='mesh_hgr.nc' !:
  TYPE(variable), DIMENSION(3) :: typvar          !: structure for attribute
 
  INTEGER    :: ncout
  INTEGER    :: istatus
  REAL(KIND=4)   ::  zpi, zomega, rau0sg
  LOGICAL        ::  lprint

  rau0sg = 1020/9.81
  lprint = .false.

  !!  Read command line
  narg= iargc()
  IF ( narg /= 3 ) THEN
     PRINT *,' Usage : cdfvpot gridT gridU gridV'
     PRINT *,' Output on vpot.nc, variables vorelvor, vostrvor,vototvor'
     PRINT *,' Need mesh_zgr.nc and coordinates.nc '
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  CALL getarg (2, cfileu)
  CALL getarg (3, cfilev)

  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')
  npt   = getdim (cfilet,'time')
 
  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'npt   =', npt

  ALLOCATE ( e1u(npiglo,npjglo)  , e1t(npiglo,npjglo) )
  ALLOCATE ( e2v(npiglo,npjglo)  , e2t(npiglo,npjglo) )
  ALLOCATE ( gphit(npiglo,npjglo), z2fcor(npiglo,npjglo))
  ALLOCATE ( zareat(npiglo,npjglo), stretch(npiglo,npjglo))
  ALLOCATE ( un(npiglo,npjglo)   , vn(npiglo,npjglo) )
  ALLOCATE ( rotn(npiglo,npjglo) , tmask(npiglo,npjglo) )
  ALLOCATE (gdepw(npk),tim(npt))

  e1u=  getvar(coord, 'e1u', 1,npiglo,npjglo)
  e1t=  getvar(coord, 'e1t', 1,npiglo,npjglo)
  e2v=  getvar(coord, 'e2v', 1,npiglo,npjglo)
  e2t=  getvar(coord, 'e2t', 1,npiglo,npjglo)
  gphit=  getvar(coord, 'gphit', 1,npiglo,npjglo)
  zpi=ACOS(-1.)
  zomega = 2*zpi/(3600*24)
  z2fcor(:,:)=2.0*zomega*SIN(gphit(:,:)*zpi/180.0)
  zareat(:,:) = 4.*e1t(:,:)*e2t(:,:)  ! factor of 4 to normalize relative vorticity

  IF (lprint) print *, ' reading gdepw  in file  ', trim(coordzgr)
  gdepw(:) = getvare3(coordzgr, 'gdepw', npk)
  IF (lprint) print *, ' read gdepw  in file  ', trim(coordzgr)
  
  tim=getvar1d(cfileu,'time_counter',1)
  ierr=putvar1d(ncout,tim,1,'T')

  ALLOCATE (ztemp(npiglo,npjglo,2), zsal(npiglo,npjglo,2)) 
  ALLOCATE (zwk(npiglo,npjglo,2) )
  ALLOCATE (zn2(npiglo,npjglo) ,    e3w(npiglo,npjglo))

  ! create output fileset

  ipk(:)= npk                   ! Those three variables are  3D
  ! define variable name and attribute
  typvar(1)%name= 'vorelvor'
  typvar(2)%name= 'vostrvor'
  typvar(3)%name= 'vototvor'
  typvar%units='kg.m-4.s-1'
  typvar%missing_value=0.
  typvar%valid_min= -1000.
  typvar%valid_max= 1000.
  typvar(1)%long_name='Relative_component_of_Ertel_PV'
  typvar(2)%long_name='Stretching_component_of_Ertel_PV'
  typvar(3)%long_name='Ertel_potential_vorticity'
  typvar(1)%short_name='vorelvor'
  typvar(2)%short_name='vostrvor'
  typvar(3)%short_name='vototvor'
  typvar%online_operation='N/A'
  typvar%axis='TZYX'


  ncout =create(cfileout, cfilet, npiglo,npjglo,npk)
  ierr= createvar   (ncout ,typvar,3, ipk,id_varout )
  ierr= putheadervar(ncout, cfilet,npiglo,npjglo,npk)


  DO jt=1,npt
  !  2 levels of T and S are required : iup,idown (with respect to W level)
  !  Compute from bottom to top (for vertical integration)
     PRINT *,'time=',jt,'(days:',tim(jt)/86400.,')'
  ztemp(:,:,idown) = getvar(cfilet, 'votemper',  npk-1  ,npiglo, npjglo, ktime=jt)
  zsal( :,:,idown) = getvar(cfilet, 'vosaline',  npk-1  ,npiglo, npjglo, ktime=jt)
  IF (lprint) print *, ' read temperature and salinity at bottom  '

  tim=getvar1d(cfilet,'time_counter',npt)
  ierr=putvar1d(ncout,tim,npt,'T')

  ! -------------------------------- LOOP OVER LEVELS
  DO jk = npk-1, 1, -1 
     PRINT *,'            level ',jk
     ! ------------------------------------RELATIVE VORTICITY FIRST
     IF (lprint) print *, ' trying to read u in file:', trim(cfileu)
     un(:,:) =  getvar(cfileu, 'vozocrtx', jk ,npiglo,npjglo,ktime=jt)
     IF (lprint) print *, ' trying to read v in file:', trim(cfilev)
     vn(:,:) =  getvar(cfilev, 'vomecrty', jk ,npiglo,npjglo,ktime=jt)
     un(:,:) = un(:,:)*e1u(:,:) ; vn(:,:) = vn(:,:)*e2v(:,:) ; 
     IF (lprint) print *, ' read u and V OK'
     !     relative vorticity at T point
     rotn(:,:) = 0.
     DO jj = 2, npjglo -1 
       DO ji = 2, npiglo -1    
         rotn(ji,jj) = (   un(ji-1,jj-1)  + un(ji,jj-1)  &
                          -un(ji-1,jj+1)  - un(ji,jj+1)  &
                          -vn(ji-1,jj-1)  - vn(ji-1,jj)  &
                          +vn(ji+1,jj-1)  + vn(ji+1,jj)) &
                          /zareat(ji,jj) 
        END DO
     END DO
     IF (lprint) print *, ' curl calculated '
     !     now  tmask and Vaisala Frequency bn2
     IF ( jk > 1) then 
        tmask(:,:)=1.
        ztemp(:,:,iup)= getvar(cfilet, 'votemper',  jk-1 ,npiglo, npjglo,ktime=jt)
        WHERE(ztemp(:,:,idown) == 0 ) tmask = 0
        zsal(:,:,iup) = getvar(cfilet, 'vosaline',  jk-1 ,npiglo,npjglo,ktime=jt)
        IF (lprint) print *, ' read temperature and salinity  '
        e3w(:,:)   = getvar(coordzgr, 'e3w_ps', jk,npiglo, npjglo )
        WHERE (e3w == 0 ) e3w = 1.
        IF (lprint) print *, ' read   e3w_ps in file  ' , trim(coordzgr)


        zwk(:,:,iup) = &
   &      eosbn2 ( ztemp,zsal,gdepw(jk),e3w, npiglo,npjglo ,iup,idown)* tmask(:,:)
        IF (lprint) print *, ' bn2 calculated at w points   '
     !
     !    now put zn2 at T level (k )
        WHERE ( zwk(:,:,idown) == 0 ) 
           zn2(:,:) =  zwk(:,:,iup)
        ELSEWHERE
           zn2(:,:) = 0.5 * ( zwk(:,:,iup) + zwk(:,:,idown) ) * tmask(:,:)
        END WHERE
        IF (lprint) print *, ' bn2 put back at T points  '
     ENDIF
     !
     !   now rotn will be converted to relative vorticity and zn2 to stretching
     rotn(:,:)     = rotn(:,:)* rau0sg * zn2(:,:)
     stretch(:,:)  = zn2(:,:) * rau0sg * z2fcor(:,:)
     
     ! write the three variables on file at level k
     ierr = putvar(ncout, id_varout(1) ,rotn, jk ,npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(2) ,stretch , jk, npiglo, npjglo , ktime=jt)
     ierr = putvar(ncout, id_varout(3) ,(rotn+stretch)*1.e7  , jk, npiglo, npjglo , ktime=jt)
     IF (lprint) print *, ' three variables written   '
     itmp = idown ; idown = iup ; iup = itmp

  END DO  ! loop to next level
 
 !  set zero at bottom 
    rotn(:,:) = 0
     ierr = putvar(ncout, id_varout(1) ,rotn, npk ,npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(2) ,rotn ,npk, npiglo, npjglo , ktime=jt)
     ierr = putvar(ncout, id_varout(3) ,rotn ,npk, npiglo, npjglo , ktime=jt)
  END DO   ! loop on time

  istatus = closeout(ncout)

END PROGRAM cdfpvort
