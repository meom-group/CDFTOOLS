program cdfsections   
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          *** PROGRAM cdfsections ***
!
! ** Purpose : extract oceanic fields along a track made of several sections.
!
! ** Method : computes N sections by taking the nearest point north of 60°N 
!             and near undefined values (bottom or coasts), and interpolates 
!             between the four nearest points elsewhere.
!             
! ** Outputs : temperature, salinity, density, current (normal/tangeantial)
!              - normal current is positive northward (westward if meridional section)
!              - tangeantial current is on the right of the normal current.
!
! NB : it is recommended to put a lot of points on each section if the aim is 
!      to compute X-integrations.
!
! WARNING : 
!  - require large memory : reduce domain size with ncks if insufficient memory error.
!  - does not work if the section crosses the Greenwich line (easy to modify if needed).
!  - not yet tested north of 60°N (but should work) ...
!
! history :
! N. JOURDAIN (LEGI-MEOM), April 2009
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------

                                             
USE netcdf                                            
USE eos
                                                      
IMPLICIT NONE                                         

!--- Local variables
INTEGER   :: narg, iargc 
                                                    
!--- grid_T
INTEGER :: fidT, status, dimID_time_counter, dimID_deptht, dimID_y, dimID_x, &
& mtime_counter, mdeptht, my, mx, vosaline_ID, votemper_ID, time_counter_ID, &
& deptht_ID, nav_lat_ID, nav_lon_ID, fidM, dimID_s, X_ID, sig0_ID, sig1_ID,  &
& sig2_ID, sig4_ID 
      
!--- grid_U 
INTEGER :: fidU, dimID_depthu, mdepthu, mxu, vozocrtx_ID

!--- grid_V
INTEGER :: fidV, dimID_depthv, mdepthv, myu, vomecrty_ID
                                                      
CHARACTER(LEN=100) :: file_in_T, file_out, file_in_U, file_in_V, cdum             
           
REAL*4 :: RT, dtmp_T, dtmp_U, dtmp_V, miniT, miniU, miniV, rr, ang, pi,&
&  latinf, latsup, loninf, lonsup, a, b, c, e, missing, lonref, latref

REAL*8 :: offset

INTEGER :: N1, N2, i, j, k, s, p, iiT, jjT, Nsec, Ntot,Unorm_ID, Utang_ID, cont,&
& l, iiU, jjU, iiV, jjV, iinf, isup, jinf, jsup

INTEGER,ALLOCATABLE,DIMENSION(:) :: N

REAL*4,ALLOCATABLE,DIMENSION(:) :: lat, lon 

LOGICAL,ALLOCATABLE,DIMENSION(:) :: undefined

!---- grid_T                                           
REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:) :: vosaline, votemper, sig0, sig1, sig2, sig4
REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: somxl010, somxlt02, vosaline_sec, votemper_sec
REAL*4,ALLOCATABLE,DIMENSION(:,:) :: nav_lat_T, nav_lon_T, somxl010_sec, somxlt02_sec
REAL*4,ALLOCATABLE,DIMENSION(:) :: time_counter, deptht
                                                     
!---- grid_U
REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:) :: vozocrtx
REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: vozocrtx_sec
REAL*4,ALLOCATABLE,DIMENSION(:,:) :: nav_lat_U, nav_lon_U

!---- grid_V
REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:) :: vomecrty
REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: vomecrty_sec
REAL*4,ALLOCATABLE,DIMENSION(:,:) :: nav_lat_V, nav_lon_V

!---- grid section 
REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: Unorm, Utang, sigsec0, sigsec1, sigsec2, sigsec4
REAL*4,ALLOCATABLE,DIMENSION(:) :: lonsec, latsec
REAL*8,ALLOCATABLE,DIMENSION(:) :: d, X1

!-------------------------------------------------------------------------
! GETTING ARGUMENTS :

!-  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg.lt.10 ) THEN
     PRINT *,'Usage : '
     PRINT *,' cdfsections  Ufile Vfile Tfile larf lorf Nsec lat1 lon1 lat2 lon2 n1'
     PRINT *,'              [ lat3 lon3 n2 ] [ lat4 lon4 n3 ] ....'
     PRINT *,'   ' 
     PRINT *,' Computes temperature, salinity, sig0, sig1, sig2, sig4, Uorth, Utang '
     PRINT *,' along a section made of Nsec linear segments (see output attributes).'
     PRINT *,' Output is section.nc, var. as a function of X(km), depth(m) and time.'
     PRINT *,'   '
     PRINT *,'Arguments : '
     PRINT *,' # larf and lorf -> location of X=0 for the X-absice (may be out of section)'
     PRINT *,' # Nsec -> number of segments used to compute the whole section.'
     PRINT *,' # lat1,lat2,lat3,... -> extrema latitudes of the segments (from -90 to 90)'
     PRINT *,' # lon1,lon2,lon3,... -> extrema latitudes of the segments (from 0 to 360)'
     PRINT *,' # n1, n2, ...        -> number of output points on each segment.'
     PRINT *,'   (you have to give Nsec+1 values of lati/loni and Nsec values of ni)'
     PRINT *,' '
     PRINT *,' It is recommended to put a lot of points on each section if the aim'
     PRINT *,' is to compute X-integrations along the section (10 x the model resolution).'
     PRINT *,'NB : sections cannot cross the Greenwich line !!'
     PRINT *,'NB : Not yet tested north of 60°N.'
     PRINT *,'NB : require a large amount of memory !' 
     PRINT *,'     -> reduce domain size with  ncks -d  if insufficient memory error.'
     PRINT *,' '
     PRINT *,'Example for one linear section : '
     PRINT *,' cdfsections U.nc V.nc T.nc 48.0 305.0 1 49.0 307.0 50.5 337.5 20'
     PRINT *,'Example for a section made of 2 linear segments : '
     PRINT *,' cdfsections U.nc V.nc T.nc 48.0 305.0 2 49.0 307.0 50.5 337.5 20 40.3 305.1 50'
     STOP
  ENDIF

  CALL getarg (1, file_in_U )
  CALL getarg (2, file_in_V )
  CALL getarg (3, file_in_T )

  CALL getarg (4, cdum ); READ(cdum,*) latref
  CALL getarg (5, cdum ); READ(cdum,*) lonref
  CALL getarg (6, cdum ); READ(cdum,*) Nsec

  if ( narg.ne.8+Nsec*3) then
    PRINT *, '**!/# ERROR : wrong number of arguments in cdfsections'
    PRINT *, 'Usage : '
    PRINT *, ' cdfsections  Ufile Vfile Tfile larf lorf Nsec lat1 lon1 lat2 lon2 n1 ....'
    PRINT *, '-> please execute cdfsections without any arguments for more details.'
    STOP
  endif

  ALLOCATE( lat(Nsec+1), lon(Nsec+1) )
  ALLOCATE( N(Nsec), d(Nsec) )

   CALL getarg (7, cdum ); READ(cdum,*) lat(1)
   CALL getarg (8, cdum ); READ(cdum,*) lon(1)

  do i=1,(narg-8),3
    CALL getarg (i+8, cdum ); READ(cdum,*) lat(i/3+2)
    CALL getarg (i+9, cdum ); READ(cdum,*) lon(i/3+2)
    CALL getarg (i+10, cdum ); READ(cdum,*) N(i/3+1)
  enddo

  do i=1,Nsec+1
    if ( (lon(i).lt.0.0).or.(lonref.lt.0.0) ) then
      PRINT *, '**!/# ERROR : longitudes must be between 0 and 360'
      STOP
    endif
  enddo
   
  file_out = 'section.nc'                               

 !---- Rayon terrestre en km :
   RT = 6378

   pi = 3.1415927
   rr = pi / 180.0         
         
!---------------------------------------                   
! Read netcdf input file for grid T :                                 
         
         write(*,*) TRIM(file_in_T)
                                                  
         status = NF90_OPEN(TRIM(file_in_T),0,fidT)          
         call erreur(status,.TRUE.,"read") 
                                                           
       !Lecture des ID des dimensions qui nous interessent
         status = NF90_INQ_DIMID(fidT,"time_counter",dimID_time_counter)
         call erreur(status,.TRUE.,"inq_dimID_time_counter")
         status = NF90_INQ_DIMID(fidT,"deptht",dimID_deptht)
         call erreur(status,.TRUE.,"inq_dimID_deptht")
         status = NF90_INQ_DIMID(fidT,"y",dimID_y)
         call erreur(status,.TRUE.,"inq_dimID_y")
         status = NF90_INQ_DIMID(fidT,"x",dimID_x)
         call erreur(status,.TRUE.,"inq_dimID_x")
                                                               
       !Lecture des valeurs des dimensions qui nous interessent
         status = NF90_INQUIRE_DIMENSION(fidT,dimID_time_counter,len=mtime_counter)
         call erreur(status,.TRUE.,"inq_dim_time_counter")
         status = NF90_INQUIRE_DIMENSION(fidT,dimID_deptht,len=mdeptht)
         call erreur(status,.TRUE.,"inq_dim_deptht")
         status = NF90_INQUIRE_DIMENSION(fidT,dimID_y,len=my)
         call erreur(status,.TRUE.,"inq_dim_y")
         status = NF90_INQUIRE_DIMENSION(fidT,dimID_x,len=mx)
         call erreur(status,.TRUE.,"inq_dim_x")
         
         write(*,101) mx, my, mdeptht, mtime_counter         
101 FORMAT('   -> dimensions of arrays : (',3(i4,','),i4,')')         
             
       !Allocation of arrays : 
         ALLOCATE(  vosaline(mx,my,mdeptht,mtime_counter)  ) 
         ALLOCATE(  votemper(mx,my,mdeptht,mtime_counter)  ) 
         ALLOCATE(  time_counter(mtime_counter)  ) 
         ALLOCATE(  deptht(mdeptht)  ) 
         ALLOCATE(  nav_lat_T(mx,my)  ) 
         ALLOCATE(  nav_lon_T(mx,my)  ) 
         ALLOCATE(  undefined(mdeptht) )         
                        
       !Lecture des ID des variables qui nous interessent
         status = NF90_INQ_VARID(fidT,"vosaline",vosaline_ID)
         call erreur(status,.TRUE.,"inq_vosaline_ID")
         status = NF90_INQ_VARID(fidT,"votemper",votemper_ID)
         call erreur(status,.TRUE.,"inq_votemper_ID")
         status = NF90_INQ_VARID(fidT,"time_counter",time_counter_ID)
         call erreur(status,.TRUE.,"inq_time_counter_ID")
         status = NF90_INQ_VARID(fidT,"deptht",deptht_ID)
         call erreur(status,.TRUE.,"inq_deptht_ID")
         status = NF90_INQ_VARID(fidT,"nav_lat",nav_lat_ID)
         call erreur(status,.TRUE.,"inq_nav_lat_ID")
         status = NF90_INQ_VARID(fidT,"nav_lon",nav_lon_ID)
         call erreur(status,.TRUE.,"inq_nav_lon_ID")
                                                              
       !Lecture des valeurs des variables qui nous interessent
         status = NF90_GET_VAR(fidT,vosaline_ID,vosaline)
         call erreur(status,.TRUE.,"getvar_vosaline")
         status = NF90_GET_VAR(fidT,votemper_ID,votemper)
         call erreur(status,.TRUE.,"getvar_votemper")
         status = NF90_GET_VAR(fidT,time_counter_ID,time_counter)
         call erreur(status,.TRUE.,"getvar_time_counter")
         status = NF90_GET_VAR(fidT,deptht_ID,deptht)
         call erreur(status,.TRUE.,"getvar_deptht")
         status = NF90_GET_VAR(fidT,nav_lat_ID,nav_lat_T)
         call erreur(status,.TRUE.,"getvar_nav_lat")
         status = NF90_GET_VAR(fidT,nav_lon_ID,nav_lon_T)
         call erreur(status,.TRUE.,"getvar_nav_lon")
         
        !extract missing value for vosaline :
         status = NF90_GET_ATT(fidT,vosaline_ID,"missing_value",missing)
         call erreur(status,.TRUE.,"get_att_vosaline")
                                             
     !Fermeture du fichier lu                         
       status = NF90_CLOSE(fidT)                      
       call erreur(status,.TRUE.,"fin_lecture")     
        
!---------------------------------------
! Read netcdf input file for grid U :

         write(*,*) TRIM(file_in_U)

         status = NF90_OPEN(TRIM(file_in_U),0,fidU)          
         call erreur(status,.TRUE.,"read") 
                                                           
       !Lecture des ID des dimensions qui nous interessent
         status = NF90_INQ_DIMID(fidU,"time_counter",dimID_time_counter)
         call erreur(status,.TRUE.,"inq_dimID_time_counter")
         status = NF90_INQ_DIMID(fidU,"depthu",dimID_depthu)
         call erreur(status,.TRUE.,"inq_dimID_depthu")
         status = NF90_INQ_DIMID(fidU,"y",dimID_y)
         call erreur(status,.TRUE.,"inq_dimID_y")
         status = NF90_INQ_DIMID(fidU,"x",dimID_x)
         call erreur(status,.TRUE.,"inq_dimID_x")
                                                               
       !Lecture des valeurs des dimensions qui nous interessent
         status = NF90_INQUIRE_DIMENSION(fidU,dimID_depthu,len=mdepthu)
         call erreur(status,.TRUE.,"inq_dim_depthu")
         status = NF90_INQUIRE_DIMENSION(fidU,dimID_x,len=mxu)
         call erreur(status,.TRUE.,"inq_dim_x")

         write(*,101) mxu, my, mdepthu, mtime_counter         
                      
       !Allocation of arrays : 
         ALLOCATE(  vozocrtx(mxu,my,mdepthu,mtime_counter)  ) 
         ALLOCATE(  nav_lat_U(mxu,my)  ) 
         ALLOCATE(  nav_lon_U(mxu,my)  ) 
                                 
       !Lecture des ID des variables qui nous interessent
         status = NF90_INQ_VARID(fidU,"vozocrtx",vozocrtx_ID)
         call erreur(status,.TRUE.,"inq_vozocrtx_ID")
         status = NF90_INQ_VARID(fidU,"nav_lat",nav_lat_ID)
         call erreur(status,.TRUE.,"inq_nav_lat_ID")
         status = NF90_INQ_VARID(fidU,"nav_lon",nav_lon_ID)
         call erreur(status,.TRUE.,"inq_nav_lon_ID")
                                                              
       !Lecture des valeurs des variables qui nous interessent
         status = NF90_GET_VAR(fidU,vozocrtx_ID,vozocrtx)
         call erreur(status,.TRUE.,"getvar_vozocrtx")
         status = NF90_GET_VAR(fidU,nav_lat_ID,nav_lat_U)
         call erreur(status,.TRUE.,"getvar_nav_lat")
         status = NF90_GET_VAR(fidU,nav_lon_ID,nav_lon_U)
         call erreur(status,.TRUE.,"getvar_nav_lon")
                                                      
     !Fermeture du fichier lu                         
       status = NF90_CLOSE(fidU)                      
       call erreur(status,.TRUE.,"fin_lecture")     


!---------------------------------------
! Read netcdf input file for grid V :

         write(*,*) TRIM(file_in_V)

         status = NF90_OPEN(TRIM(file_in_V),0,fidV)          
         call erreur(status,.TRUE.,"read") 
                                                           
       !Lecture des ID des dimensions qui nous interessent
         status = NF90_INQ_DIMID(fidV,"depthv",dimID_depthv)
         call erreur(status,.TRUE.,"inq_dimID_depthv")
         status = NF90_INQ_DIMID(fidV,"y",dimID_y)
         call erreur(status,.TRUE.,"inq_dimID_y")
                                                               
       !Lecture des valeurs des dimensions qui nous interessent
         status = NF90_INQUIRE_DIMENSION(fidV,dimID_depthv,len=mdepthv)
         call erreur(status,.TRUE.,"inq_dim_depthv")
         status = NF90_INQUIRE_DIMENSION(fidV,dimID_y,len=myu)
         call erreur(status,.TRUE.,"inq_dim_y")
         write(*,101) mx, myu, mdepthv, mtime_counter
                               
       !Allocation of arrays : 
         ALLOCATE(  vomecrty(mx,myu,mdepthv,mtime_counter)  ) 
         ALLOCATE(  nav_lat_V(mx,myu)  ) 
         ALLOCATE(  nav_lon_V(mx,myu)  ) 
                                 
       !Lecture des ID des variables qui nous interessent
         status = NF90_INQ_VARID(fidV,"vomecrty",vomecrty_ID)
         call erreur(status,.TRUE.,"inq_vomecrty_ID")
         status = NF90_INQ_VARID(fidV,"nav_lat",nav_lat_ID)
         call erreur(status,.TRUE.,"inq_nav_lat_ID")
         status = NF90_INQ_VARID(fidV,"nav_lon",nav_lon_ID)
         call erreur(status,.TRUE.,"inq_nav_lon_ID")
                                                              
       !Lecture des valeurs des variables qui nous interessent
         status = NF90_GET_VAR(fidV,vomecrty_ID,vomecrty)
         call erreur(status,.TRUE.,"getvar_vomecrty")
         status = NF90_GET_VAR(fidV,nav_lat_ID,nav_lat_V)
         call erreur(status,.TRUE.,"getvar_nav_lat")
         status = NF90_GET_VAR(fidV,nav_lon_ID,nav_lon_V)
         call erreur(status,.TRUE.,"getvar_nav_lon")
                                                      
     !Fermeture du fichier lu                         
       status = NF90_CLOSE(fidV)                      
       call erreur(status,.TRUE.,"fin_lecture")     


!------------------------------------------------------------------------------- 
!-------------------------------------------------------------------------------                                                     

!-------------------------------------------------------------------------------
! Remise des longitudes de 0 à 360 (utile pour l'interpolation):

do i=1,mx
do j=1,my
  if (nav_lon_T(i,j).lt.0.0) then
    nav_lon_T(i,j) = 360.0 + nav_lon_T(i,j)
  endif
enddo
do j=1,myu
  if (nav_lon_V(i,j).lt.0.0) then
    nav_lon_V(i,j) = 360.0 + nav_lon_V(i,j)
  endif
enddo  
enddo

do i=1,mxu
do j=1,my
  if (nav_lon_U(i,j).lt.0.0) then
    nav_lon_U(i,j) = 360.0 + nav_lon_U(i,j)
  endif
enddo
enddo


!-------------------------------------------------------------------------------
! Calcul des densite avel le module eos des CDFTOOLS-2.1

 ALLOCATE( sig0(mx,my,mdeptht,mtime_counter) )
 ALLOCATE( sig1(mx,my,mdeptht,mtime_counter) )
 ALLOCATE( sig2(mx,my,mdeptht,mtime_counter) )
 ALLOCATE( sig4(mx,my,mdeptht,mtime_counter) )

 do k=1,mdeptht
 do l=1,mtime_counter
      sig0(:,:,k,l)=sigma0(votemper(:,:,k,l),vosaline(:,:,k,l),mx,my)
      sig1(:,:,k,l)=sigmai(votemper(:,:,k,l),vosaline(:,:,k,l),1000.,mx,my)
      sig2(:,:,k,l)=sigmai(votemper(:,:,k,l),vosaline(:,:,k,l),2000.,mx,my)
      sig4(:,:,k,l)=sigmai(votemper(:,:,k,l),vosaline(:,:,k,l),4000.,mx,my)
 enddo
 enddo

!-------------------------------------------------------------------------------
! Calcul de la longueur des sections et des points modeles associes                 
        
  Ntot=SUM(N(:))
  write(*,*) '********** total number of points  :', Ntot
 
  ALLOCATE( latsec(Ntot) , lonsec(Ntot) , X1(Ntot) )
  ALLOCATE( votemper_sec(Ntot,mdeptht,mtime_counter) ) 
  ALLOCATE( vosaline_sec(Ntot,mdeptht,mtime_counter) )
  ALLOCATE( vozocrtx_sec(Ntot,mdepthu,mtime_counter) )
  ALLOCATE( vomecrty_sec(Ntot,mdepthv,mtime_counter) ) 
  ALLOCATE( Unorm(Ntot,mdeptht,mtime_counter) )
  ALLOCATE( Utang(Ntot,mdeptht,mtime_counter) )
  ALLOCATE( sigsec0(Ntot,mdeptht,mtime_counter) )
  ALLOCATE( sigsec1(Ntot,mdeptht,mtime_counter) )
  ALLOCATE( sigsec2(Ntot,mdeptht,mtime_counter) )
  ALLOCATE( sigsec4(Ntot,mdeptht,mtime_counter) )

! BOUCLE SUR LE NOMBRE DE SECTIONS Nsec A ACOLLER :
N2=0
! Point de référence pour la distance de la section (exple : dans OVIDE 60N 43.25W)
offset=RT * acos(cos(latref*rr)*cos(lat(1)*rr)*cos((lonref)*rr-lon(1)*rr)+sin(latref*rr)*sin(lat(1)*rr))
cont=0
DO p=1,Nsec
 N1=N2+1
 N2=N1+N(p)-1
 !longueur de la section p en km :
  d(p) = RT * acos(cos(lat(p)*rr)*cos(lat(p+1)*rr)*cos(lon(p+1)*rr-lon(p)*rr)+sin(lat(p)*rr)*sin(lat(p+1)*rr))
  write(*,102) p,d(p)
102 FORMAT('*** Section ',i4,' = ',f8.2, 'km')
  write(*,103) lat(p), lon(p), lat(p+1), lon(p+1)
103 FORMAT('     - from (lat,lon) = (',f6.2,',',f6.2,') to (',f6.2,',',f6.2,')')
 ! "pente" de la section 1 en radians / equateur (angle algebrique)
  if (lon(p).ne.lon(p+1)) then
    ang = atan((lat(p+1)-lat(p))/(lon(p+1)-lon(p)))
  else
    ang=pi/2.
  endif
  write(*,*) '    - angle / equateur (deg) =', ang/rr
 !coordonnées de tous les points de la section p en (lon,lat) et en km :
  DO s=N1,N2
    undefined(:)=.FALSE.
    latsec(s)=(lat(p+1)-lat(p))*FLOAT(s-N1+cont)/FLOAT(N2-N1+cont) + lat(p)
    lonsec(s)=(lon(p+1)-lon(p))*FLOAT(s-N1+cont)/FLOAT(N2-N1+cont) + lon(p)
    X1(s)=d(p)*FLOAT(s-N1+cont)/FLOAT(N2-N1+cont)+offset
    miniT=1000 !km
    miniU=miniT
    miniV=miniT
   ! recherche du point le plus proche (on fait comme ça parceque la grille est bizarre vers les poles)
    do i=1,mx
    do j=1,my
      dtmp_T= RT * acos(cos(nav_lat_T(i,j)*rr)*cos(latsec(s)*rr)*cos(nav_lon_T(i,j)*rr-lonsec(s)*rr)+sin(nav_lat_T(i,j)*rr)*sin(latsec(s)*rr))
      dtmp_U= RT * acos(cos(nav_lat_U(i,j)*rr)*cos(latsec(s)*rr)*cos(nav_lon_U(i,j)*rr-lonsec(s)*rr)+sin(nav_lat_U(i,j)*rr)*sin(latsec(s)*rr))
      dtmp_V= RT * acos(cos(nav_lat_V(i,j)*rr)*cos(latsec(s)*rr)*cos(nav_lon_V(i,j)*rr-lonsec(s)*rr)+sin(nav_lat_V(i,j)*rr)*sin(latsec(s)*rr))
      if (dtmp_T.lt.miniT) then
          miniT=dtmp_T
          iiT=i
          jjT=j
      endif
      if (dtmp_U.lt.miniU) then
          miniU=dtmp_U
          iiU=i
          jjU=j
      endif
      if (dtmp_V.lt.miniV) then
          miniV=dtmp_V
          iiV=i
          jjV=j
      endif
    enddo
    enddo
   !interpolation des champs T:
    if (latsec(s).gt.60.0) then  
     !champs le plus proche de la section (U et V interpoles au point T)
      votemper_sec(s,:,:) = votemper(iiT,jjT,:,:)
      vosaline_sec(s,:,:) = vosaline(iiT,jjT,:,:)
      vozocrtx_sec(s,:,:) = vozocrtx(iiU,jjU,:,:)
      vomecrty_sec(s,:,:) = vomecrty(iiV,jjV,:,:)
     ! vitesse normale et tangeantielle a la section (section orientee vers le nord, tangeante à droite)
      Unorm(s,:,:) = vomecrty_sec(s,:,:)*cos(ang) - vozocrtx_sec(s,:,:)*sin(ang) 
      Utang(s,:,:) = vomecrty_sec(s,:,:)*sin(ang) + vozocrtx_sec(s,:,:)*cos(ang)
     ! densites :
      sigsec0(s,:,:) = MAX(sig0(iiT,jjT,:,:),10.0)
      sigsec1(s,:,:) = MAX(sig1(iiT,jjT,:,:),10.0)
      sigsec2(s,:,:) = MAX(sig2(iiT,jjT,:,:),10.0)
      sigsec4(s,:,:) = MAX(sig4(iiT,jjT,:,:),20.0)
    else
     ! Champs T interpoles
      if (lonsec(s).ge.nav_lon_T(iiT,jjT)) then
        iinf=iiT
        if ( iiT+1.le.mx ) then
          isup=iiT+1
        else
          isup=1
        endif
      else
        if ( iiT-1.ge.1 ) then
          iinf=iiT-1
        else
          iinf=mx
        endif
        isup=iiT
      endif
      if (latsec(s).ge.nav_lat_T(iiT,jjT)) then
        jinf=jjT
        jsup=jjT+1
      else
        jinf=jjT-1
        jsup=jjT
      endif
      loninf=nav_lon_T(iinf,jjT)
      lonsup=nav_lon_T(isup,jjT)
      latinf=nav_lat_T(iiT,jinf)
      latsup=nav_lat_T(iiT,jsup)
      a=(lonsec(s)-loninf)/(lonsup-loninf)
      b=(lonsup-lonsec(s))/(lonsup-loninf)
      c=(latsec(s)-latinf)/(latsup-latinf)
      e=(latsup-latsec(s))/(latsup-latinf)
      votemper_sec(s,:,:) =  c*(a*votemper(isup,jsup,:,:)+b*votemper(iinf,jsup,:,:)) &
&                           +e*(a*votemper(isup,jinf,:,:)+b*votemper(iinf,jinf,:,:))
      vosaline_sec(s,:,:) =  c*(a*vosaline(isup,jsup,:,:)+b*vosaline(iinf,jsup,:,:)) &
&                           +e*(a*vosaline(isup,jinf,:,:)+b*vosaline(iinf,jinf,:,:))
      sigsec0(s,:,:) =  c*(a*sig0(isup,jsup,:,:)+b*sig0(iinf,jsup,:,:)) &
&                      +e*(a*sig0(isup,jinf,:,:)+b*sig0(iinf,jinf,:,:))
      sigsec0(s,:,:) = MAX(sigsec0(s,:,:),10.0)
      sigsec1(s,:,:) =  c*(a*sig1(isup,jsup,:,:)+b*sig1(iinf,jsup,:,:)) &
&                      +e*(a*sig1(isup,jinf,:,:)+b*sig1(iinf,jinf,:,:))
      sigsec1(s,:,:) = MAX(sigsec1(s,:,:),10.0)
      sigsec2(s,:,:) =  c*(a*sig2(isup,jsup,:,:)+b*sig2(iinf,jsup,:,:)) &
&                      +e*(a*sig2(isup,jinf,:,:)+b*sig2(iinf,jinf,:,:))
      sigsec2(s,:,:) = MAX(sigsec2(s,:,:),10.0)
      sigsec4(s,:,:) =  c*(a*sig4(isup,jsup,:,:)+b*sig4(iinf,jsup,:,:)) &
&                      +e*(a*sig4(isup,jinf,:,:)+b*sig4(iinf,jinf,:,:))
      sigsec4(s,:,:) = MAX(sigsec4(s,:,:),20.0)
     ! test si valeurs indefinies sur un des 4 points :
      do k=1,mdeptht
        if ((vosaline(iinf,jinf,k,1).eq.missing).or.(vosaline(iinf,jsup,k,1).eq.missing).or.&
&           (vosaline(isup,jinf,k,1).eq.missing).or.(vosaline(isup,jsup,k,1).eq.missing)   ) then
          votemper_sec(s,:,:) = votemper(iiT,jjT,:,:)
          vosaline_sec(s,:,:) = vosaline(iiT,jjT,:,:)
          vozocrtx_sec(s,:,:) = vozocrtx(iiU,jjU,:,:)
          vomecrty_sec(s,:,:) = vomecrty(iiV,jjV,:,:)
          sigsec0(s,:,:) = MAX(sig0(iiT,jjT,:,:),10.0)
          sigsec1(s,:,:) = MAX(sig1(iiT,jjT,:,:),10.0)
          sigsec2(s,:,:) = MAX(sig2(iiT,jjT,:,:),10.0)
          sigsec4(s,:,:) = MAX(sig4(iiT,jjT,:,:),20.0) 
          undefined(k)=.TRUE.
        endif
      enddo
     ! Champs U interpoles
      if (lonsec(s).ge.nav_lon_U(iiU,jjU)) then
        iinf=iiU
        if ( iiU+1.le.mx ) then
          isup=iiU+1
        else
          isup=1
        endif
      else
        if ( iiU-1.ge.1 ) then
          iinf=iiU-1
        else
          iinf=mx
        endif
        isup=iiU
      endif
      if (latsec(s).ge.nav_lat_U(iiU,jjU)) then
       jinf=jjU
       jsup=jjU+1
      else
       jinf=jjU-1
       jsup=jjU
      endif
      loninf=nav_lon_U(iinf,jjU)
      lonsup=nav_lon_U(isup,jjU)
      latinf=nav_lat_U(iiU,jinf)
      latsup=nav_lat_U(iiU,jsup)
      a=(lonsec(s)-loninf)/(lonsup-loninf)
      b=(lonsup-lonsec(s))/(lonsup-loninf)
      c=(latsec(s)-latinf)/(latsup-latinf)
      e=(latsup-latsec(s))/(latsup-latinf)
      vozocrtx_sec(s,:,:) =  c*(a*vozocrtx(isup,jsup,:,:)+b*vozocrtx(iinf,jsup,:,:)) &
&                           +e*(a*vozocrtx(isup,jinf,:,:)+b*vozocrtx(iinf,jinf,:,:))
     ! Champs V interpoles
      if (lonsec(s).ge.nav_lon_U(iiU,jjU)) then
        iinf=iiU
        if ( iiU+1.le.mx ) then
          isup=iiU+1
        else
          isup=1
        endif
      else
        if ( iiU-1.ge.1 ) then
          iinf=iiU-1
        else
          iinf=mx
        endif
        isup=iiU
      endif
      if (latsec(s).ge.nav_lat_V(iiV,jjV)) then
       jinf=jjV
       jsup=jjV+1
      else
       jinf=jjV-1
       jsup=jjV
      endif
      loninf=nav_lon_V(iinf,jjV)
      lonsup=nav_lon_V(isup,jjV)
      latinf=nav_lat_V(iiV,jinf)
      latsup=nav_lat_V(iiV,jsup)
      a=(lonsec(s)-loninf)/(lonsup-loninf)
      b=(lonsup-lonsec(s))/(lonsup-loninf)
      c=(latsec(s)-latinf)/(latsup-latinf)
      e=(latsup-latsec(s))/(latsup-latinf)
      vomecrty_sec(s,:,:) =  c*(a*vomecrty(isup,jsup,:,:)+b*vomecrty(iinf,jsup,:,:)) &
&                           +e*(a*vomecrty(isup,jinf,:,:)+b*vomecrty(iinf,jinf,:,:))
     ! si l'un des 4 points de l'interpolation etait indefini :
      do k=1,mdeptht
        if (undefined(k)) then
          vozocrtx_sec(s,k,:) = vozocrtx(iiU,jjU,k,:)
          vomecrty_sec(s,k,:) = vomecrty(iiV,jjV,k,:)
        endif
      enddo
     ! vitesse normale et tangeantielle a la section (section orientee vers le nord, tangeante à droite)
      Unorm(s,:,:) = vomecrty_sec(s,:,:)*cos(ang) - vozocrtx_sec(s,:,:)*sin(ang)
      Utang(s,:,:) = vomecrty_sec(s,:,:)*sin(ang) + vozocrtx_sec(s,:,:)*cos(ang)
    endif
  ENDDO !- s nb de points sur section p
  cont=1
  offset=X1(N2)
ENDDO !- p nb de sections 
      
!----------------------------------------------------------                                                      
!----------------------------------------------------------    
! Writing new netcdf file :                                   
                                                              
        status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM)
        call erreur(status,.TRUE.,'create')                     
                                                                
        !Definition des dimensions du fichiers                  
         status = NF90_DEF_DIM(fidM,"time_counter",NF90_UNLIMITED,dimID_time_counter)
         call erreur(status,.TRUE.,"def_dimID_time_counter")
         status = NF90_DEF_DIM(fidM,"deptht",mdeptht,dimID_deptht)
         call erreur(status,.TRUE.,"def_dimID_deptht")
         status = NF90_DEF_DIM(fidM,"X",Ntot,dimID_s)
         call erreur(status,.TRUE.,"def_dimID_s")
                                                              
        !Definition des variables                             
         status = NF90_DEF_VAR(fidM,"vosaline",NF90_FLOAT,(/dimID_s,dimID_deptht,dimID_time_counter/),vosaline_ID)
         call erreur(status,.TRUE.,"def_var_vosaline_ID")
         status = NF90_DEF_VAR(fidM,"votemper",NF90_FLOAT,(/dimID_s,dimID_deptht,dimID_time_counter/),votemper_ID)
         call erreur(status,.TRUE.,"def_var_votemper_ID")
         status = NF90_DEF_VAR(fidM,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID)
         call erreur(status,.TRUE.,"def_var_time_counter_ID")
         status = NF90_DEF_VAR(fidM,"deptht",NF90_FLOAT,(/dimID_deptht/),deptht_ID)
         call erreur(status,.TRUE.,"def_var_deptht_ID")
         status = NF90_DEF_VAR(fidM,"nav_lat",NF90_FLOAT,(/dimID_s/),nav_lat_ID)
         call erreur(status,.TRUE.,"def_var_nav_lat_ID")
         status = NF90_DEF_VAR(fidM,"nav_lon",NF90_FLOAT,(/dimID_s/),nav_lon_ID)
         call erreur(status,.TRUE.,"def_var_nav_lon_ID")
         status = NF90_DEF_VAR(fidM,"X",NF90_DOUBLE,(/dimID_s/),X_ID)
         call erreur(status,.TRUE.,"def_var_X_ID")     
         status = NF90_DEF_VAR(fidM,"Uorth",NF90_FLOAT,(/dimID_s,dimID_deptht,dimID_time_counter/),Unorm_ID)
         call erreur(status,.TRUE.,"def_var_Unorm_ID")    
         status = NF90_DEF_VAR(fidM,"Utang",NF90_FLOAT,(/dimID_s,dimID_deptht,dimID_time_counter/),Utang_ID)
         call erreur(status,.TRUE.,"def_var_Utang_ID")
         status = NF90_DEF_VAR(fidM,"sig0",NF90_FLOAT,(/dimID_s,dimID_deptht,dimID_time_counter/),sig0_ID)
         call erreur(status,.TRUE.,"def_var_sig0_ID")
         status = NF90_DEF_VAR(fidM,"sig1",NF90_FLOAT,(/dimID_s,dimID_deptht,dimID_time_counter/),sig1_ID)
         call erreur(status,.TRUE.,"def_var_sig1_ID")
         status = NF90_DEF_VAR(fidM,"sig2",NF90_FLOAT,(/dimID_s,dimID_deptht,dimID_time_counter/),sig2_ID)
         call erreur(status,.TRUE.,"def_var_sig2_ID")
         status = NF90_DEF_VAR(fidM,"sig4",NF90_FLOAT,(/dimID_s,dimID_deptht,dimID_time_counter/),sig4_ID)
         call erreur(status,.TRUE.,"def_var_sig4_ID")

                          
        ! Attributs des variables :
         status = NF90_PUT_ATT(fidM,vosaline_ID,"online_operation","N/A")
         call erreur(status,.TRUE.,"put_att_vosaline_ID")
         status = NF90_PUT_ATT(fidM,vosaline_ID,"short_name","vosaline")
         call erreur(status,.TRUE.,"put_att_vosaline_ID")
         status = NF90_PUT_ATT(fidM,vosaline_ID,"long_name","sea_water_salinity")
         call erreur(status,.TRUE.,"put_att_vosaline_ID")
         status = NF90_PUT_ATT(fidM,vosaline_ID,"valid_max",45.)
         call erreur(status,.TRUE.,"put_att_vosaline_ID")
         status = NF90_PUT_ATT(fidM,vosaline_ID,"valid_min",0.)
         call erreur(status,.TRUE.,"put_att_vosaline_ID")
         status = NF90_PUT_ATT(fidM,vosaline_ID,"missing_value",0.)
         call erreur(status,.TRUE.,"put_att_vosaline_ID")
         status = NF90_PUT_ATT(fidM,vosaline_ID,"units","PSU")
         call erreur(status,.TRUE.,"put_att_vosaline_ID")

         status = NF90_PUT_ATT(fidM,votemper_ID,"online_operation","N/A")
         call erreur(status,.TRUE.,"put_att_votemper_ID")
         status = NF90_PUT_ATT(fidM,votemper_ID,"short_name","votemper")
         call erreur(status,.TRUE.,"put_att_votemper_ID")
         status = NF90_PUT_ATT(fidM,votemper_ID,"long_name","sea_water_potential_temperature")
         call erreur(status,.TRUE.,"put_att_votemper_ID")
         status = NF90_PUT_ATT(fidM,votemper_ID,"valid_max",45.)
         call erreur(status,.TRUE.,"put_att_votemper_ID")
         status = NF90_PUT_ATT(fidM,votemper_ID,"valid_min",-2.)
         call erreur(status,.TRUE.,"put_att_votemper_ID")
         status = NF90_PUT_ATT(fidM,votemper_ID,"missing_value",0.)
         call erreur(status,.TRUE.,"put_att_votemper_ID")
         status = NF90_PUT_ATT(fidM,votemper_ID,"units","C")
         call erreur(status,.TRUE.,"put_att_votemper_ID")

         status = NF90_PUT_ATT(fidM,time_counter_ID,"time_origin","2001-OCT-03 00:00:00")
         call erreur(status,.TRUE.,"put_att_time_counter_ID")
         status = NF90_PUT_ATT(fidM,time_counter_ID,"units","seconds since 2001-10-03 00:00:00")
         call erreur(status,.TRUE.,"put_att_time_counter_ID")
         status = NF90_PUT_ATT(fidM,time_counter_ID,"long_name","Time axis")
         call erreur(status,.TRUE.,"put_att_time_counter_ID")
         status = NF90_PUT_ATT(fidM,time_counter_ID,"title","Time")
         call erreur(status,.TRUE.,"put_att_time_counter_ID")
         status = NF90_PUT_ATT(fidM,time_counter_ID,"calendar","gregorian")
         call erreur(status,.TRUE.,"put_att_time_counter_ID")

         status = NF90_PUT_ATT(fidM,deptht_ID,"long_name","Vertical T levels")
         call erreur(status,.TRUE.,"put_att_deptht_ID")
         status = NF90_PUT_ATT(fidM,deptht_ID,"title","deptht")
         call erreur(status,.TRUE.,"put_att_deptht_ID")
         status = NF90_PUT_ATT(fidM,deptht_ID,"valid_max",50.)
         call erreur(status,.TRUE.,"put_att_deptht_ID")
         status = NF90_PUT_ATT(fidM,deptht_ID,"valid_min",0.)
         call erreur(status,.TRUE.,"put_att_deptht_ID")
         status = NF90_PUT_ATT(fidM,deptht_ID,"positive","unknown")
         call erreur(status,.TRUE.,"put_att_deptht_ID")
         status = NF90_PUT_ATT(fidM,deptht_ID,"units","m")
         call erreur(status,.TRUE.,"put_att_deptht_ID")

         status = NF90_PUT_ATT(fidM,nav_lat_ID,"long_name","Latitude")
         call erreur(status,.TRUE.,"put_att_nav_lat_ID")
         status = NF90_PUT_ATT(fidM,nav_lat_ID,"scale_factor",1.)
         call erreur(status,.TRUE.,"put_att_nav_lat_ID")
         status = NF90_PUT_ATT(fidM,nav_lat_ID,"add_offset",0.)
         call erreur(status,.TRUE.,"put_att_nav_lat_ID")
         status = NF90_PUT_ATT(fidM,nav_lat_ID,"valid_max",89.947868347168)
         call erreur(status,.TRUE.,"put_att_nav_lat_ID")
         status = NF90_PUT_ATT(fidM,nav_lat_ID,"valid_min",-77.0104751586914)
         call erreur(status,.TRUE.,"put_att_nav_lat_ID")
         status = NF90_PUT_ATT(fidM,nav_lat_ID,"units","degrees_north")
         call erreur(status,.TRUE.,"put_att_nav_lat_ID")

         status = NF90_PUT_ATT(fidM,nav_lon_ID,"long_name","Longitude")
         call erreur(status,.TRUE.,"put_att_nav_lon_ID")
         status = NF90_PUT_ATT(fidM,nav_lon_ID,"scale_factor",1.)
         call erreur(status,.TRUE.,"put_att_nav_lon_ID")
         status = NF90_PUT_ATT(fidM,nav_lon_ID,"add_offset",0.)
         call erreur(status,.TRUE.,"put_att_nav_lon_ID")
         status = NF90_PUT_ATT(fidM,nav_lon_ID,"valid_max",180.)
         call erreur(status,.TRUE.,"put_att_nav_lon_ID")
         status = NF90_PUT_ATT(fidM,nav_lon_ID,"valid_min",-180.)
         call erreur(status,.TRUE.,"put_att_nav_lon_ID")
         status = NF90_PUT_ATT(fidM,nav_lon_ID,"units","degrees_east")
         call erreur(status,.TRUE.,"put_att_nav_lon_ID")

         status = NF90_PUT_ATT(fidM,X_ID,"nav_model","Default grid")
         call erreur(status,.TRUE.,"put_att_X_ID")
         status = NF90_PUT_ATT(fidM,X_ID,"long_name","X")
         call erreur(status,.TRUE.,"put_att_X_ID")
         status = NF90_PUT_ATT(fidM,X_ID,"scale_factor",1.)
         call erreur(status,.TRUE.,"put_att_X_ID")
         status = NF90_PUT_ATT(fidM,X_ID,"add_offset",0.)
         call erreur(status,.TRUE.,"put_att_X_ID")
         status = NF90_PUT_ATT(fidM,X_ID,"units","km")
         call erreur(status,.TRUE.,"put_att_X_ID")

         status = NF90_PUT_ATT(fidM,Unorm_ID,"online_operation","N/A")
         call erreur(status,.TRUE.,"put_att_Unorm_ID")
         status = NF90_PUT_ATT(fidM,Unorm_ID,"short_name","Uorth")
         call erreur(status,.TRUE.,"put_att_Unorm_ID")
         status = NF90_PUT_ATT(fidM,Unorm_ID,"long_name","ocean speed orthogonal to the section oriented south-north")
         call erreur(status,.TRUE.,"put_att_Unorm_ID")
         status = NF90_PUT_ATT(fidM,Unorm_ID,"valid_max",10.)
         call erreur(status,.TRUE.,"put_att_Unorm_ID")
         status = NF90_PUT_ATT(fidM,Unorm_ID,"valid_min",-10.)
         call erreur(status,.TRUE.,"put_att_Unorm_ID")
         status = NF90_PUT_ATT(fidM,Unorm_ID,"missing_value",0.)
         call erreur(status,.TRUE.,"put_att_Unorm_ID")
         status = NF90_PUT_ATT(fidM,Unorm_ID,"units","m/s")
         call erreur(status,.TRUE.,"put_att_Unorm_ID")

         status = NF90_PUT_ATT(fidM,Utang_ID,"online_operation","N/A")
         call erreur(status,.TRUE.,"put_att_Utang_ID")
         status = NF90_PUT_ATT(fidM,Utang_ID,"short_name","Utang")
         call erreur(status,.TRUE.,"put_att_Utang_ID")
         status = NF90_PUT_ATT(fidM,Utang_ID,"long_name","ocean speed tangential to the section oriented south-north")
         call erreur(status,.TRUE.,"put_att_Utang_ID")
         status = NF90_PUT_ATT(fidM,Utang_ID,"valid_max",10.)
         call erreur(status,.TRUE.,"put_att_Utang_ID")
         status = NF90_PUT_ATT(fidM,Utang_ID,"valid_min",-10.)
         call erreur(status,.TRUE.,"put_att_Utang_ID")
         status = NF90_PUT_ATT(fidM,Utang_ID,"missing_value",0.)
         call erreur(status,.TRUE.,"put_att_Utang_ID")
         status = NF90_PUT_ATT(fidM,Utang_ID,"units","m/s")
         call erreur(status,.TRUE.,"put_att_Utang_ID")
                         
         status = NF90_PUT_ATT(fidM,sig0_ID,"short_name","sig0")
         call erreur(status,.TRUE.,"put_att_sig0_ID")
         status = NF90_PUT_ATT(fidM,sig0_ID,"long_name","Potential Density Sigma 0")
         call erreur(status,.TRUE.,"put_att_sig0_ID")
         status = NF90_PUT_ATT(fidM,sig0_ID,"valid_max",100.)
         call erreur(status,.TRUE.,"put_att_sig0_ID")
         status = NF90_PUT_ATT(fidM,sig0_ID,"valid_min",10.0)
         call erreur(status,.TRUE.,"put_att_sig0_ID")
         status = NF90_PUT_ATT(fidM,sig0_ID,"missing_value",10.0)
         call erreur(status,.TRUE.,"put_att_sig0_ID")
         status = NF90_PUT_ATT(fidM,sig0_ID,"units","kg/m3")
         call erreur(status,.TRUE.,"put_att_sig0_ID")

         status = NF90_PUT_ATT(fidM,sig1_ID,"short_name","sig1")
         call erreur(status,.TRUE.,"put_att_sig1_ID")
         status = NF90_PUT_ATT(fidM,sig1_ID,"long_name","Potential Density Sigma 0")
         call erreur(status,.TRUE.,"put_att_sig1_ID")
         status = NF90_PUT_ATT(fidM,sig1_ID,"valid_max",100.)
         call erreur(status,.TRUE.,"put_att_sig1_ID")
         status = NF90_PUT_ATT(fidM,sig1_ID,"valid_min",10.0)
         call erreur(status,.TRUE.,"put_att_sig1_ID")
         status = NF90_PUT_ATT(fidM,sig1_ID,"missing_value",10.0)
         call erreur(status,.TRUE.,"put_att_sig1_ID")
         status = NF90_PUT_ATT(fidM,sig1_ID,"units","kg/m3")
         call erreur(status,.TRUE.,"put_att_sig1_ID")

         status = NF90_PUT_ATT(fidM,sig2_ID,"short_name","sig2")
         call erreur(status,.TRUE.,"put_att_sig2_ID")
         status = NF90_PUT_ATT(fidM,sig2_ID,"long_name","Potential Density Sigma 0")
         call erreur(status,.TRUE.,"put_att_sig2_ID")
         status = NF90_PUT_ATT(fidM,sig2_ID,"valid_max",100.)
         call erreur(status,.TRUE.,"put_att_sig2_ID")
         status = NF90_PUT_ATT(fidM,sig2_ID,"valid_min",10.0)
         call erreur(status,.TRUE.,"put_att_sig2_ID")
         status = NF90_PUT_ATT(fidM,sig2_ID,"missing_value",10.0)
         call erreur(status,.TRUE.,"put_att_sig2_ID")
         status = NF90_PUT_ATT(fidM,sig2_ID,"units","kg/m3")
         call erreur(status,.TRUE.,"put_att_sig2_ID")

         status = NF90_PUT_ATT(fidM,sig4_ID,"short_name","sig4")
         call erreur(status,.TRUE.,"put_att_sig4_ID")
         status = NF90_PUT_ATT(fidM,sig4_ID,"long_name","Potential Density Sigma 0")
         call erreur(status,.TRUE.,"put_att_sig4_ID")
         status = NF90_PUT_ATT(fidM,sig4_ID,"valid_max",100.)
         call erreur(status,.TRUE.,"put_att_sig4_ID")
         status = NF90_PUT_ATT(fidM,sig4_ID,"valid_min",20.0)
         call erreur(status,.TRUE.,"put_att_sig4_ID")
         status = NF90_PUT_ATT(fidM,sig4_ID,"missing_value",20.0)
         call erreur(status,.TRUE.,"put_att_sig4_ID")
         status = NF90_PUT_ATT(fidM,sig4_ID,"units","kg/m3")
         call erreur(status,.TRUE.,"put_att_sig4_ID")

        ! Attributs globaux       :
         status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created by cdfsections (see CDFTOOLS)")
         call erreur(status,.TRUE.,"put_att_global_ID")
                                                      
        !Fin des definitions                          
         status = NF90_ENDDEF(fidM)                   
         call erreur(status,.TRUE.,"fin_definition") 
                                                      
        !Valeurs prises par les variables :           
         status = NF90_PUT_VAR(fidM,vosaline_ID,vosaline_sec)
         call erreur(status,.TRUE.,"var_vosaline_ID")
         status = NF90_PUT_VAR(fidM,votemper_ID,votemper_sec)
         call erreur(status,.TRUE.,"var_votemper_ID")
         status = NF90_PUT_VAR(fidM,time_counter_ID,time_counter)
         call erreur(status,.TRUE.,"var_time_counter_ID")
         status = NF90_PUT_VAR(fidM,deptht_ID,deptht)
         call erreur(status,.TRUE.,"var_deptht_ID")
         status = NF90_PUT_VAR(fidM,nav_lat_ID,latsec)
         call erreur(status,.TRUE.,"var_nav_lat_ID")
         status = NF90_PUT_VAR(fidM,nav_lon_ID,lonsec)
         call erreur(status,.TRUE.,"var_nav_lon_ID")
         status = NF90_PUT_VAR(fidM,X_ID,X1)
         call erreur(status,.TRUE.,"var_X_ID")
         status = NF90_PUT_VAR(fidM,Unorm_ID,Unorm)
         call erreur(status,.TRUE.,"var_Unorm_ID")
         status = NF90_PUT_VAR(fidM,Utang_ID,Utang)
         call erreur(status,.TRUE.,"var_Utang_ID") 
         status = NF90_PUT_VAR(fidM,sig0_ID,sigsec0)
         call erreur(status,.TRUE.,"var_sig0_ID")
         status = NF90_PUT_VAR(fidM,sig1_ID,sigsec1)
         call erreur(status,.TRUE.,"var_sig1_ID")
         status = NF90_PUT_VAR(fidM,sig2_ID,sigsec2)
         call erreur(status,.TRUE.,"var_sigsec2_ID") 
         status = NF90_PUT_VAR(fidM,sig4_ID,sigsec4)
         call erreur(status,.TRUE.,"var_sigsec4_ID")         
                                    
        !Fin de l'ecriture                            
         status = NF90_CLOSE(fidM)                    
         call erreur(status,.TRUE.,"final")         

end program cdfsections



SUBROUTINE erreur(iret, lstop, chaine)
  ! pour les messages d'erreur
  USE netcdf
  INTEGER, INTENT(in)                     :: iret
  LOGICAL, INTENT(in)                     :: lstop
  CHARACTER(LEN=*), INTENT(in)            :: chaine
  !
  CHARACTER(LEN=80)                       :: message
  !
  IF ( iret .NE. 0 ) THEN
    WRITE(*,*) 'ROUTINE: ', TRIM(chaine)
    WRITE(*,*) 'ERREUR: ', iret
    message=NF90_STRERROR(iret)
    WRITE(*,*) 'THIS MEANS :',TRIM(message)
    IF ( lstop ) STOP
  ENDIF
  !
END SUBROUTINE erreur
