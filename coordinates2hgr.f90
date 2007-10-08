PROGRAM coordinates2hgr
  !!-------------------------------------------------------------------------
  !!                        PROGRAM coordinates2hgr
  !!                        **********************
  !! ** Purpose:
  !!         transform a "clipper" coordinates file into an
  !!         ioipsl coordinate file
  !!
  !! History :
  !!   February 2003 : Anne de Miranda
  !!   June 2003     : J.M. Molines : modif for getarg (filename)
  !!--------------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------

  USE netcdf
  IMPLICIT NONE

  INTEGER(4) ::  jpi , jpj
  INTEGER(2) , PARAMETER :: wp = 8
  INTEGER(4) , PARAMETER :: jpk=43

  REAL(wp) :: zjpiglo, zjpjglo, znbsel, zrecl8, ztimm8
  REAL(wp) :: znt, zdim, xx1, yy1, ddx, ddy, sspval

  CHARACTER(80) :: cltextco
  INTEGER(4) :: numcoo, nummsh
  INTEGER(4) :: nrecl8
  INTEGER(4) :: ngrid
  INTEGER(4) :: ni,nj, nk
  INTEGER(4) :: jj, jk

  REAL(wp) , DIMENSION(:,:), ALLOCATABLE  ::  zlamt, zphit, zdta
  REAL(wp) , DIMENSION(1)   ::  zdep
  REAL(wp)                     zdt,zdate0 , omega, pi
  REAL(wp), DIMENSION(jpk) ::  gdept, gdepw, e3t, e3w

  CHARACTER (len=21) ::        clname
  INTEGER                      ilev, itime, iargc, narg
  LOGICAL                      clog
  CHARACTER(LEN=80)       :: cfilin, cfilout

  REAL(wp) , DIMENSION(:,:),  ALLOCATABLE ::   &
       glamt, glamu, glamv, glamf, gphit, gphiu, gphiv, gphif, &
       e1t, e1u, e1v, e1f, e2t, e2u, e2v, e2f, ff

 ! netcdf stuff
 INTEGER :: istatus, ncid
 INTEGER :: id_x, id_y, id_z, id_time, id_xa, id_ya, id_za
 INTEGER :: id_lon, id_lat, id_lev, id_tim, id_ts
 INTEGER :: id_lamt, id_lamu, id_lamv, id_lamf
 INTEGER :: id_phit, id_phiu, id_phiv, id_phif
 INTEGER :: id_e1t, id_e1u, id_e1v, id_e1f
 INTEGER :: id_e2t, id_e2u, id_e2v, id_e2f, id_ff


  numcoo = 10
  nummsh = 11
  narg=iargc()

  IF ( narg /= 1 ) THEN
   PRINT *,' >>> Usage: coordinates2hgr ''coordinates file'' '
   PRINT *,'     Output is done on mesh_hgr.nc '
   STOP 
  END IF

  CALL getarg(1,cfilin)
  cfilout='mesh_hgr.nc'
  !
  ! ... Read coordinates (only the used metric is read)
  !
  nrecl8 = 200
  OPEN( numcoo, FILE = cfilin, status='old' & 
       ,form = 'unformatted',      access ='direct',recl = nrecl8)
  READ(numcoo, rec = 1)      cltextco, zrecl8, zjpiglo, zjpjglo
  CLOSE(numcoo)
  !
  print *, cltextco, zrecl8, zjpiglo, zjpjglo
  nrecl8 = zrecl8
  print*, nrecl8

  OPEN(numcoo, FILE=cfilin, status='old', form='unformatted', &
       access='direct',recl=nrecl8)
  READ (numcoo, rec = 1) cltextco, zrecl8, zjpiglo, zjpjglo, znbsel, znt, &
       zdim, xx1, yy1, ddx, ddy, sspval &
       ,(gdept(jk),jk=1,jpk),                           &
       ztimm8,                                          &
       (gdepw(jk),jk=1,jpk),                            &
       (e3t(jk),jk=1,jpk),                              &
       (e3w(jk),jk=1,jpk)
  ni = zjpiglo
  nj = zjpjglo
  print *, 'ni=',ni,'nj=',nj
  ALLOCATE(glamt(ni,nj))
  ALLOCATE(glamu(ni,nj))
  ALLOCATE(glamv(ni,nj))
  ALLOCATE(glamf(ni,nj))
  ALLOCATE(gphit(ni,nj))
  ALLOCATE(gphiu(ni,nj))
  ALLOCATE(gphiv(ni,nj))
  ALLOCATE(gphif(ni,nj))
  ALLOCATE(e1t(ni,nj))
  ALLOCATE(e1u(ni,nj))
  ALLOCATE(e1v(ni,nj))
  ALLOCATE(e1f(ni,nj))
  ALLOCATE(e2t(ni,nj))
  ALLOCATE(e2u(ni,nj))
  ALLOCATE(e2v(ni,nj))
  ALLOCATE(e2f(ni,nj))
  ALLOCATE(ff(ni,nj))

  READ(numcoo,rec=2) glamt
  READ(numcoo,rec=3) glamu
  READ(numcoo,rec=4) glamv
  READ(numcoo,rec=5) glamf
  READ(numcoo,rec=6) gphit
  READ(numcoo,rec=7) gphiu
  READ(numcoo,rec=8) gphiv
  READ(numcoo,rec=9) gphif
  READ(numcoo,rec=10) e1t
  READ(numcoo,rec=11) e1u
  READ(numcoo,rec=12) e1v
  READ(numcoo,rec=13) e1f
  READ(numcoo,rec=14) e2t
  READ(numcoo,rec=15) e2u
  READ(numcoo,rec=16) e2v
  READ(numcoo,rec=17) e2f

  print *, 'Reading ',TRIM(cfilin),'  OK.'
  print *,'e2f:',(jj,e2f(1,jj),jj=289,nj)
  pi=acos(-1.d0)
  omega=2*pi/86400.d0
  ff=2*omega*sin(gphit*pi/180.)

  clname = cfilout
  ilev = 1
  itime = 0
  zdate0 = 0.
  zdt = 0.
  clog = .FALSE.

  ALLOCATE(zlamt(ni,nj))
  ALLOCATE(zphit(ni,nj))

  zlamt(:,:) = 0.
  zphit(:,:) = 0.
  zdep(1) = 0.

  jpi = zjpiglo
  jpj = zjpjglo
  print* , jpi, jpj

  print *,ngrid, ncid
  istatus=NF90_CREATE(cfilout,NF90_CLOBBER,ncid)
  ! define dimension x, y, z=1, time=unlimited
 istatus=NF90_DEF_DIM(ncid,'x',jpi,id_x)
 istatus=NF90_DEF_DIM(ncid,'y',jpj,id_y)
 istatus=NF90_DEF_DIM(ncid,'z',jpk  ,id_z)
 istatus=NF90_DEF_DIM(ncid,'time',NF90_UNLIMITED ,id_time)
 istatus=NF90_DEF_DIM(ncid,'x_a',1,id_xa)
 istatus=NF90_DEF_DIM(ncid,'y_a',1,id_ya)
 istatus=NF90_DEF_DIM(ncid,'z_a',1  ,id_za)

  ! define variables
  istatus=NF90_DEF_VAR(ncid,'nav_lon',NF90_FLOAT,(/id_x,id_y/),id_lon)
  istatus=NF90_DEF_VAR(ncid,'nav_lat',NF90_FLOAT,(/id_x,id_y/),id_lat)
  istatus=NF90_DEF_VAR(ncid,'nav_lev',NF90_FLOAT,(/id_z/),id_lev)
  istatus=NF90_DEF_VAR(ncid,'time'   ,NF90_FLOAT,(/id_time/),id_tim)
  istatus=NF90_DEF_VAR(ncid,'time_steps',NF90_INT,(/id_time/),id_ts)

  ! Horizontal grid-point position
  istatus=NF90_DEF_VAR(ncid,'glamt',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_lamt)
  istatus=NF90_DEF_VAR(ncid,'glamu',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_lamu)
  istatus=NF90_DEF_VAR(ncid,'glamv',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_lamv)
  istatus=NF90_DEF_VAR(ncid,'glamf',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_lamf)

  istatus=NF90_DEF_VAR(ncid,'gphit',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_phit)
  istatus=NF90_DEF_VAR(ncid,'gphiu',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_phiu)
  istatus=NF90_DEF_VAR(ncid,'gphiv',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_phiv)
  istatus=NF90_DEF_VAR(ncid,'gphif',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_phif)
 
  ! Horizontal scale factors
  istatus=NF90_DEF_VAR(ncid,'e1t',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_e1t)
  istatus=NF90_DEF_VAR(ncid,'e1u',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_e1u)
  istatus=NF90_DEF_VAR(ncid,'e1v',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_e1v)
  istatus=NF90_DEF_VAR(ncid,'e1f',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_e1f)

  istatus=NF90_DEF_VAR(ncid,'e2t',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_e2t)
  istatus=NF90_DEF_VAR(ncid,'e2u',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_e2u)
  istatus=NF90_DEF_VAR(ncid,'e2v',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_e2v)
  istatus=NF90_DEF_VAR(ncid,'e2f',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_e2f)

  istatus=NF90_DEF_VAR(ncid,'ff',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_ff)

  ! attributes
  !nav_lon:
  istatus=NF90_PUT_ATT(ncid,id_lon,'units','degrees_east')
  istatus=NF90_PUT_ATT(ncid,id_lon,'valid_min',-180.)
  istatus=NF90_PUT_ATT(ncid,id_lon,'valid_max',180.)
  istatus=NF90_PUT_ATT(ncid,id_lon,'long_name','Longitude')
  !nav_lat:
  istatus=NF90_PUT_ATT(ncid,id_lat,'units','degrees_north')
  istatus=NF90_PUT_ATT(ncid,id_lat,'valid_min',-90.)
  istatus=NF90_PUT_ATT(ncid,id_lat,'valid_max',90.)
  istatus=NF90_PUT_ATT(ncid,id_lat,'long_name','Latitude')
  !nav_lev:
  istatus=NF90_PUT_ATT(ncid,id_lev,'units','model_levels')
  istatus=NF90_PUT_ATT(ncid,id_lev,'valid_min',0.)
  istatus=NF90_PUT_ATT(ncid,id_lev,'valid_max',0.)
  istatus=NF90_PUT_ATT(ncid,id_lev,'long_name','Model levels')
  !time:
  istatus=NF90_PUT_ATT(ncid,id_tim,'units','seconds since 0000-01-01 00:00:00')
  istatus=NF90_PUT_ATT(ncid,id_tim,'calendar','gregorian')
  istatus=NF90_PUT_ATT(ncid,id_tim,'title','Time')
  istatus=NF90_PUT_ATT(ncid,id_tim,'long_name','Time axis')
  istatus=NF90_PUT_ATT(ncid,id_tim,'time_origin',' 0000-JAN-01 00:00:00')
  !time_steps:
  istatus=NF90_PUT_ATT(ncid,id_ts,'units','timesteps since 0000-01-01 00:00:00')
  istatus=NF90_PUT_ATT(ncid,id_ts,'title','Time steps')
  istatus=NF90_PUT_ATT(ncid,id_ts,'tstep_sec',0.)
  istatus=NF90_PUT_ATT(ncid,id_ts,'long_name','Time step axis')
  istatus=NF90_PUT_ATT(ncid,id_ts,'time_origin',' 0000-JAN-01 00:00:00')
 
  ! variables glamx, gphix, e?x
  istatus=NF90_PUT_ATT(ncid,id_lamt,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_lamu,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_lamv,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_lamf,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_phit,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_phiu,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_phiv,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_phif,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_e1t,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_e1u,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_e1v,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_e2f,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_e2t,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_e2u,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_e2v,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_e2f,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_ff,'missing_value',1.e+20)

  istatus=NF90_ENDDEF(ncid)

  ! Now fill the variables !
  istatus=NF90_PUT_VAR(ncid,id_lon,glamt)
  istatus=NF90_PUT_VAR(ncid,id_lat,gphit)
  istatus=NF90_PUT_VAR(ncid,id_lev,gdept)
  istatus=NF90_PUT_VAR(ncid,id_tim,0.)
  istatus=NF90_PUT_VAR(ncid,id_ts,1)

  istatus=NF90_PUT_VAR(ncid,id_lamt,glamt)
  istatus=NF90_PUT_VAR(ncid,id_lamu,glamu)
  istatus=NF90_PUT_VAR(ncid,id_lamv,glamv)
  istatus=NF90_PUT_VAR(ncid,id_lamf,glamf)
  
  istatus=NF90_PUT_VAR(ncid,id_phit,gphit)
  istatus=NF90_PUT_VAR(ncid,id_phiu,gphiu)
  istatus=NF90_PUT_VAR(ncid,id_phiv,gphiv)
  istatus=NF90_PUT_VAR(ncid,id_phif,gphif)

  istatus=NF90_PUT_VAR(ncid,id_e1t,e1t)
  istatus=NF90_PUT_VAR(ncid,id_e1u,e1u)
  istatus=NF90_PUT_VAR(ncid,id_e1v,e1v)
  istatus=NF90_PUT_VAR(ncid,id_e1f,e1f)

  istatus=NF90_PUT_VAR(ncid,id_e2t,e2t)
  istatus=NF90_PUT_VAR(ncid,id_e2u,e2u)
  istatus=NF90_PUT_VAR(ncid,id_e2v,e2v)
  istatus=NF90_PUT_VAR(ncid,id_e2f,e2f)
  istatus=NF90_PUT_VAR(ncid,id_ff,ff)

  istatus=NF90_CLOSE(ncid)

END PROGRAM coordinates2hgr
