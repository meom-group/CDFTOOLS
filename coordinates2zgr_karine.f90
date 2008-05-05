PROGRAM coordinates2zgr
  !!-------------------------------------------------------------------------
  !!                        PROGRAM coordinates2zgr
  !!                        **********************
  !! ** Purpose:
  !!         transform a "clipper" coordinates file into an
  !!         ioipsl coordinate file
  !!
  !! History :
  !!   February 2003 : Anne de Miranda
  !!   June 2003     : J.M. Molines : modif for getarg (filename)
  !!--------------------------------------------------------------------------
  !!  $Rev: 110 $
  !!  $Date: 2007-10-10 18:07:55 +0200 (mer, 10 oct 2007) $
  !!  $Id: coordinates2zgr.f90 110 2007-10-10 16:07:55Z molines $
  !!--------------------------------------------------------------

  USE netcdf
  IMPLICIT NONE

  INTEGER(4) ::  jpi , jpj
  INTEGER(2) , PARAMETER :: wp = 8
  INTEGER(4) , PARAMETER :: jpk=43

  REAL(wp) :: zjpiglo, zjpjglo, znbsel, zrecl8, ztimm8
  REAL(wp) :: znt, zdim, xx1, yy1, ddx, ddy, sspval, bidon

  CHARACTER(80) :: cltextco
  INTEGER(4) :: numcoo, nummsh, numbat
  INTEGER(4) :: nrecl8
  INTEGER(4) :: ngrid
  INTEGER(4) :: ni,nj, nk
  INTEGER(4) :: jj, jk, iim,ijm, il1, il2, ifreq, jn, ii, ji , ij
  INTEGER(KIND=4) :: ierr=0, npi, npj, npizoom=1, npjzoom=1
  INTEGER(KIND=4) , DIMENSION(:,:), ALLOCATABLE :: idata

  REAL(wp) , DIMENSION(:,:), ALLOCATABLE  ::  glamt,gphit
  REAL(wp) , DIMENSION(1)   ::  zdep
  REAL(wp)                     zdt,zdate0 , omega, pi
  REAL(wp), DIMENSION(jpk) ::  gdept, gdepw, e3t, e3w

  CHARACTER (len=21) ::        clname
  INTEGER                      ilev, itime, iargc, narg
  LOGICAL                      clog
  CHARACTER(LEN=80)       :: cfilin, cfilout, cbathy, clexp, clfmt, cdum

  REAL(wp) , DIMENSION(:,:),  ALLOCATABLE ::  bathy  
  LOGICAL :: ln_glo

 ! netcdf stuff
 INTEGER :: istatus, ncid
 INTEGER :: id_x, id_y, id_z, id_time, id_xa, id_ya, id_za
 INTEGER :: id_lon, id_lat, id_lev, id_tim, id_ts
 INTEGER :: id_bat, id_dept, id_depw, id_e3t, id_e3w


  numcoo = 10
  nummsh = 11
  numbat = 12
  narg=iargc()

  IF ( narg < 2 ) THEN
   PRINT *,' >>> Usage: coordinates2zgr ''coordinates file'' '' ascii bathy'' [ jpizoom  jpjzoom]'
   PRINT *,'     Output is done on mesh_zgr.nc '
   PRINT *,'  If optional arguments jpizoom and jpjzoom are given, bathy is extracted with regard to these values'
   PRINT *,'  the global domain size is then read from the header of bathy file '
   STOP 
  END IF

  CALL getarg(1,cfilin)
  CALL getarg(2,cbathy)
  IF ( narg > 2 ) THEN
   CALL getarg(3, cdum) ; READ(cdum,*) npizoom
   CALL getarg(4, cdum) ; READ(cdum,*) npjzoom
  ENDIF
  cfilout='mesh_zgr.nc'

  ! Read ASCII BATHY_LEVEL
         OPEN( UNIT=numbat, FILE=cbathy, FORM='FORMATTED', &
                 ACCESS='SEQUENTIAL', STATUS='OLD' )
         ! read bathymetry file
         REWIND numbat
         READ(numbat,9001) clexp, iim, ijm
         nrecl8 = 200
         ! open and read header of coordinates files ( 
         OPEN( numcoo, FILE = cfilin, status='old' & 
         ,form = 'unformatted',      access ='direct',recl = nrecl8)
         READ(numcoo, rec = 1)      cltextco, zrecl8, zjpiglo, zjpjglo
         CLOSE(numcoo)
         print *, iim,ijm, INT(zjpiglo), INT(zjpjglo)
       
         ! BATHY ASCII is alway on the full domain  
         ! coordinates.diags may be on the full domain or on the zoomed domain: check dimensions !
          npi=INT(zjpiglo) ; npj=INT(zjpjglo) ; ierr = 0
         IF ( iim == npi  .AND. ijm == npj ) THEN
          ln_glo=.true.
         ELSE
          ln_glo=.false.
         ENDIF
       
         ! ALLOCATE bathy array
         ALLOCATE(idata(iim,ijm), bathy(npi,npj), glamt(npi,npj),gphit(npi,npj) )
         READ(numbat,'(/)')
         clfmt =  '(i3,41i3)'
         IF ( ijm >= 1000 ) clfmt = '(i4,41i3)'
         ifreq=40
         il1=1
         DO jn=1,iim/ifreq+1
            READ(numbat,'(/)')
            il2 = MIN( iim, il1+ifreq-1 )
            READ(numbat,9002) ( ii, ji = il1, il2, 5 )
            READ(numbat,'(/)')
            DO jj = ijm, 1, -1
               READ(numbat,clfmt) ij, ( idata(ji,jj), ji = il1, il2 )
            END DO
            il1 = il1 + ifreq
         END DO
         CLOSE(numbat)
         bathy(:,:)=idata(npizoom:npizoom+npi-1, npjzoom:npjzoom+npj-1)

9001     FORMAT(1x,a15,2i8)
9002     FORMAT(3x,13(i3,12x))
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
  READ(numcoo,rec=1) &
       cltextco,zrecl8,zjpiglo,zjpjglo,znbsel,znt,zdim, &
       xx1,yy1,ddx,ddy,sspval,                          &
       (bidon,jk=1,INT(znbsel)),                          &
       ztimm8,                                          &
       (gdept(jk),jk=1,jpk),                           &
       (bidon,jk=1,jpk),                          &
       (gdepw(jk),jk=1,jpk),                            &
       (bidon,jk=1,jpk),                          &
       (e3t(jk),jk=1,jpk),                              &
       (bidon,jk=1,jpk),                          &
       (e3w(jk),jk=1,jpk)
  ni = zjpiglo
  nj = zjpjglo
  print *, 'ni=',ni,'nj=',nj
  READ(numcoo,rec=3) glamt
  READ(numcoo,rec=7) gphit

  clname = cfilout
  ilev = 1
  itime = 0
  zdate0 = 0.
  zdt = 0.
  clog = .FALSE.

  jpi = zjpiglo
  jpj = zjpjglo
  print* , jpi, jpj

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
  istatus=NF90_DEF_VAR(ncid,'mbathy',NF90_DOUBLE,(/id_x,id_y,id_za,id_time/),id_bat)
  istatus=NF90_DEF_VAR(ncid,'gdept',NF90_DOUBLE,(/id_xa,id_ya,id_z,id_time/),id_dept)
  istatus=NF90_DEF_VAR(ncid,'gdepw',NF90_DOUBLE,(/id_xa,id_ya,id_z,id_time/),id_depw)
  istatus=NF90_DEF_VAR(ncid,'e3t',NF90_DOUBLE,(/id_xa,id_ya,id_z,id_time/),id_e3t)
  istatus=NF90_DEF_VAR(ncid,'e3w',NF90_DOUBLE,(/id_xa,id_ya,id_z,id_time/),id_e3w)

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
  istatus=NF90_PUT_ATT(ncid,id_bat,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_dept,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_depw,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_e3t,'missing_value',1.e+20)
  istatus=NF90_PUT_ATT(ncid,id_e3w,'missing_value',1.e+20)

  istatus=NF90_ENDDEF(ncid)

  ! Now fill the variables !
  istatus=NF90_PUT_VAR(ncid,id_lon,glamt)
  istatus=NF90_PUT_VAR(ncid,id_lat,gphit)
  istatus=NF90_PUT_VAR(ncid,id_lev,gdept)
  istatus=NF90_PUT_VAR(ncid,id_tim,0.)
  istatus=NF90_PUT_VAR(ncid,id_ts,1)

  istatus=NF90_PUT_VAR(ncid,id_bat,bathy )
  istatus=NF90_PUT_VAR(ncid,id_dept,gdept, start=(/1,1,1,1/),count=(/1,1,jpk,1/))
  istatus=NF90_PUT_VAR(ncid,id_depw,gdepw , start=(/1,1,1,1/),count=(/1,1,jpk,1/))
  istatus=NF90_PUT_VAR(ncid,id_e3t,e3t,  start=(/1,1,1,1/),count=(/1,1,jpk,1/))
  istatus=NF90_PUT_VAR(ncid,id_e3w,e3w, start=(/1,1,1,1/),count=(/1,1,jpk,1/))

  istatus=NF90_CLOSE(ncid)

END PROGRAM coordinates2zgr
