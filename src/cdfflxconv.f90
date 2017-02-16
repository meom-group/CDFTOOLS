PROGRAM cdfflxconv
  !!-------------------------------------------------------------------
  !!              PROGRAM CDFFLXCONV
  !!              ******************
  !!
  !!  **  Purpose: Convert a set of fluxes dimgfile (Clipper like)
  !!               to a set of CDF files (Drakkar like )
  !!  
  !!  **  Method: takes the current year as input, and config name
  !!              automatically read 
  !!                  ECMWF.Y${year}.M??.FLUX.${config}.dimg (daily, 1 file per month)
  !!                  ECMWF.Y${year}.M??.STRESS.${config}.dimg (daily, 1 file per month)
  !!                  REYNOLDS.Y${year}.SST.${config}.dimg ( weekly, 1 file per year ) ! Danger !
  !!              creates 6 netcdf daily files :
  !!                  ECMWF_emp_1d_${year}.${config}.nc                
  !!                  ECMWF_qnet_1d_${year}.${config}.nc                
  !!                  ECMWF_qsr_1d_${year}.${config}.nc                
  !!                  ECMWF_sst_1d_${year}.${config}.nc                
  !!                  ECMWF_taux_1d_${year}.${config}.nc                
  !!                  ECMWF_tauy_1d_${year}.${config}.nc                
  !!              Requires  coordinates.diags file (to be input consistent)
  !!
  !! history:
  !!    Original:  J.M. Molines (Feb. 2007 )
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !!
  !! * Modules used
  USE cdfio 

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk, jvar, jmonth, jdim, jday, jt
  INTEGER   :: narg, iargc, nvar
  INTEGER   :: npiglo,npjglo, npk                           !: size of the domain
  INTEGER   :: iyear, icurrday, jul, jul1, jul2
  INTEGER   :: id1, id2, ii1, ii2, ntime, ntp, ntn, itt
  INTEGER   :: january1, december31
  INTEGER, DIMENSION(:), ALLOCATABLE :: itime

  REAL(KIND=4) , DIMENSION (:,:,:), ALLOCATABLE :: v2d
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE ::  glam, gphi, z2d, v2daily
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE ::  glamu, gphiu
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE ::  glamv, gphiv
  REAL(KIND=4) , DIMENSION (:), ALLOCATABLE :: dep, timetab
  REAL(KIND=8) , DIMENSION (:), ALLOCATABLE ::  timetag, timetagp,timetagn
  REAL(KIND=4) ,DIMENSION(1)                  :: timean 

  CHARACTER(LEN=2048) :: ctag, confcase

  ! Dimg stuff
  INTEGER   :: irecl, ii, nt, ndim, irec
  INTEGER   :: numflx=10, numcoo=11, numtau=12, numsst=14, numsstp=15, numsstn=16
  CHARACTER(LEN=2048) :: cflux, ctau, csstr,csstrp, csstrn
  CHARACTER(LEN=2048) :: coord='coordinates.diags'
  CHARACTER(LEN=2048) :: cheader, cdum, config
  CHARACTER(LEN=4) :: cver
  REAL(KIND=4) :: x1,y1, dx,dy, spval
  ! coordinates.diags
  INTEGER :: nrecl8
  REAL(KIND=8) :: zrecl8, zpiglo,zpjglo
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE ::  dzvar
  CHARACTER(LEN=2048) :: cltextco
  LOGICAL :: lexist

  ! Netcdf Stuff
  CHARACTER(LEN=2048) :: cemp, cqnet, cqsr, ctaux, ctauy, csst
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvaremp,typvarqnet,typvarqsr
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvartaux,typvartauy,typvarsst
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipkemp, ipkqnet, ipkqsr, id_varoutemp,id_varoutqnet, id_varoutqsr
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipktaux, ipktauy, ipksst, id_varouttaux,id_varouttauy, id_varoutsst
  INTEGER    :: ncoutemp, ncoutqnet, ncoutqsr, ncouttaux, ncouttauy, ncoutsst
  INTEGER    :: istatus

  !!  Read command line
  narg= iargc()
  IF ( narg /= 2 ) THEN
     PRINT *,' Usage : cdfflxconv YEAR config '
     PRINT *,'    Output 6 cdf files : for emp, qnet, qsr, sst, taux, tauy with standard var name :'
     PRINT *,'        sowaflup, sohefldo, soshfldo, sst, sozotaux, sometauy '
     PRINT *,'    coordinates.diags ( clipper like) is required in current dir '
     STOP
  ENDIF
  !!
  CALL getarg (1, cdum)
  READ(cdum,*) iyear
  CALL getarg (2, config)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!     .....    FLUXES    FLUXES   FLUXES .....                                                     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PRINT *, 'Doing fluxes ... '

  !! read glam gphi in the coordinates file for T point (fluxes)
  nrecl8=200
  OPEN(numcoo,FILE=coord,status='old' ,form='unformatted', access='direct',recl=nrecl8)
  READ(numcoo,rec=1) cltextco,zrecl8,zpiglo,zpjglo
  CLOSE(numcoo)
  nrecl8=zrecl8 ;  npiglo=zpiglo ; npjglo=zpjglo
  ALLOCATE ( glam(npiglo,npjglo), gphi(npiglo,npjglo) ,dzvar(npiglo,npjglo) )
  OPEN(numcoo,FILE=coord,status='old' ,form='unformatted', access='direct',recl=nrecl8)
  READ(numcoo,REC=2)((dzvar(ji,jj),ji=1,npiglo),jj=1,npjglo)  ; glam(:,:) = dzvar(:,:)
  READ(numcoo,REC=6)((dzvar(ji,jj),ji=1,npiglo),jj=1,npjglo)  ; gphi(:,:) = dzvar(:,:)
  DEALLOCATE ( dzvar )
  CLOSE(numcoo)

  !! build nc output files
  WRITE(cemp,'(a,I4.4,a)') 'ECMWF_emp_1d_',iyear,'.'//TRIM(config)//'.nc'
  WRITE(cqnet,'(a,I4.4,a)') 'ECMWF_qnet_1d_',iyear,'.'//TRIM(config)//'.nc'
  WRITE(cqsr,'(a,I4.4,a)') 'ECMWF_qsr_1d_',iyear,'.'//TRIM(config)//'.nc'

  jmonth=1
  !! Build dimg file names
  WRITE(cflux,'(a,I4.4,a,I2.2,a)') 'ECMWF.Y',iyear,'.M',jmonth,'.FLUX.'//TRIM(config)//'.dimg'
  ! WRITE(ctau ,'(a,I4.4,a,I2.2,a)') 'ECMWF.Y',iyear,'.M',imonth,'.STRESS.'//TRIM(config)//'.dimg'
  ! WRITE(csst ,'(a,I4.4,a,I2.2,a)') 'REYNOLDS.Y',iyear,'.SST.'//TRIM(config)//'.dimg'

  ! open (and check ?? if they exists )
  irecl=isdirect(cflux)  ; OPEN( numflx,FILE=cflux, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )

  READ(numflx,REC=1) cver, cheader, ii, npiglo, npjglo, npk

  ALLOCATE (v2d(npiglo, npjglo,4), dep(npk) )
  ALLOCATE (z2d(npiglo, npjglo) )
  READ(numflx,REC=1) cver, cheader, ii, npiglo, npjglo, npk, nt, ndim, &
       x1,y1,dx,dy,spval,   &
       (dep(jk),jk=1,npk), &
       timean(1)
  CLOSE(numflx)

  ! Build  cdf files output
  nvar = 1 ! 1 var but many files ... (OK  ... 3 actually )
  ALLOCATE ( typvaremp(nvar), ipkemp(nvar), id_varoutemp(nvar) )
  ALLOCATE ( typvarqnet(nvar), ipkqnet(nvar), id_varoutqnet(nvar) )
  ALLOCATE ( typvarqsr(nvar), ipkqsr(nvar), id_varoutqsr(nvar) )
  jvar=1
  ipkemp(jvar)      = 1
  typvaremp(jvar)%cname='sowaflup'    ! E - P = dim 3 - dim 4 dimgfile
  typvaremp(jvar)%cunits='kg/m2/s'
  typvaremp(jvar)%rmissing_value=0.
  typvaremp(jvar)%valid_min= -0.002
  typvaremp(jvar)%valid_max= 0.002
  typvaremp(jvar)%clong_name='E-P Upward water flux'
  typvaremp(jvar)%cshort_name='sowaflup'
  typvaremp(jvar)%conline_operation='N/A'
  typvaremp(jvar)%caxis='TYX'

  ipkqnet(jvar)      = 1
  typvarqnet(jvar)%cname='sohefldo'   !  QNET = dim 1 dimgfile
  typvarqnet(jvar)%cunits='W/m2'
  typvarqnet(jvar)%rmissing_value=0.
  typvarqnet(jvar)%valid_min= -1000.
  typvarqnet(jvar)%valid_max= 1000.
  typvarqnet(jvar)%clong_name='Net_Downward_Heat_Flux'
  typvarqnet(jvar)%cshort_name='sohefldo'
  typvarqnet(jvar)%conline_operation='N/A'
  typvarqnet(jvar)%caxis='TYX'

  ipkqsr(jvar)      = 1
  typvarqsr(jvar)%cname='soshfldo'   ! QSR = dim 2 dimgfile
  typvarqsr(jvar)%cunits='W/m2'
  typvarqsr(jvar)%rmissing_value=0.
  typvarqsr(jvar)%valid_min= -1000.
  typvarqsr(jvar)%valid_max= 1000.
  typvarqsr(jvar)%clong_name='Short_Wave_Radiation'
  typvarqsr(jvar)%cshort_name='soshfldo'
  typvarqsr(jvar)%conline_operation='N/A'
  typvarqsr(jvar)%caxis='TYX'

  ncoutemp =create(cemp, 'none',npiglo,npjglo,npk,cdep='deptht' )
  istatus= createvar(ncoutemp ,typvaremp,nvar, ipkemp,id_varoutemp )
  istatus= putheadervar(ncoutemp, 'none', npiglo, npjglo,npk, pnavlon=glam,pnavlat=gphi,pdep=dep )
  ncoutqnet =create(cqnet, 'none',npiglo,npjglo,npk,cdep='deptht' )
  istatus= createvar(ncoutqnet ,typvarqnet,nvar, ipkqnet,id_varoutqnet )
  istatus= putheadervar(ncoutqnet, 'none', npiglo, npjglo,npk, pnavlon=glam,pnavlat=gphi,pdep=dep )
  ncoutqsr =create(cqsr, 'none',npiglo,npjglo,npk,cdep='deptht' )
  istatus= createvar(ncoutqsr ,typvarqsr,nvar, ipkqsr,id_varoutqsr )
  istatus= putheadervar(ncoutqsr, 'none', npiglo, npjglo,npk, pnavlon=glam,pnavlat=gphi,pdep=dep )

  ! Ready for time loop on month
  icurrday=0
  DO jmonth = 1, 12
     WRITE(cflux,'(a,I4.4,a,I2.2,a)') 'ECMWF.Y',iyear,'.M',jmonth,'.FLUX.'//TRIM(config)//'.dimg'
     irecl=isdirect(cflux)  ; OPEN( numflx,FILE=cflux, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
     READ(numflx,REC=1) cver, cheader, ii, npiglo, npjglo, npk, nt, ndim
     ! loop for days in files
     DO jday=1,nt
        icurrday=icurrday +1
        DO jdim=1,ndim
           irec=1+(jday-1)*ndim +jdim
           READ(numflx,REC=irec) (( v2d(ji,jj,jdim),ji=1,npiglo),jj=1,npjglo) 
        END DO
        ! emp
        z2d=(v2d(:,:,3) - v2d(:,:,4) )/ 86400.   ! scaling from mm/d to kg/m2/s
        istatus = putvar(ncoutemp,id_varoutemp(1),z2d,icurrday,npiglo,npjglo)
        ! qnet
        istatus = putvar(ncoutqnet,id_varoutqnet(1),v2d(:,:,1),icurrday,npiglo,npjglo)
        ! qsr
        istatus = putvar(ncoutqsr,id_varoutqsr(1),v2d(:,:,2),icurrday,npiglo,npjglo)
     END DO ! loop on days
     CLOSE(numflx)
  END DO ! loop on month

  ! update time_counter
  ALLOCATE( timetab (icurrday) )
  timetab=(/(jt,jt=1,icurrday)/)
  istatus=putvar1d(ncoutemp,timetab,icurrday,'T')
  istatus=putvar1d(ncoutqnet,timetab,icurrday,'T')
  istatus=putvar1d(ncoutqsr,timetab,icurrday,'T')
  ! close fluxes files
  istatus=closeout(ncoutemp)
  istatus=closeout(ncoutqnet)
  istatus=closeout(ncoutqsr)
  DEALLOCATE (v2d , dep, z2d , timetab )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!     .....    STRESSES    STRESSES    STRESSES   ......                                       !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PRINT *,' Doing Stresses ...'

  !! read glam gphi in the coordinates file for U point (fluxes)
  nrecl8=200
  OPEN(numcoo,FILE=coord,status='old' ,form='unformatted', access='direct',recl=nrecl8)
  READ(numcoo,rec=1) cltextco,zrecl8,zpiglo,zpjglo
  CLOSE(numcoo)
  nrecl8=zrecl8 ;  npiglo=zpiglo ; npjglo=zpjglo
  ALLOCATE ( glamu(npiglo,npjglo), gphiu(npiglo,npjglo) ,dzvar(npiglo,npjglo) )
  ALLOCATE ( glamv(npiglo,npjglo), gphiv(npiglo,npjglo) )
  OPEN(numcoo,FILE=coord,status='old' ,form='unformatted', access='direct',recl=nrecl8)
  READ(numcoo,REC=3)((dzvar(ji,jj),ji=1,npiglo),jj=1,npjglo)  ; glamu(:,:) = dzvar(:,:)
  READ(numcoo,REC=7)((dzvar(ji,jj),ji=1,npiglo),jj=1,npjglo)  ; gphiu(:,:) = dzvar(:,:)
  READ(numcoo,REC=4)((dzvar(ji,jj),ji=1,npiglo),jj=1,npjglo)  ; glamv(:,:) = dzvar(:,:)
  READ(numcoo,REC=8)((dzvar(ji,jj),ji=1,npiglo),jj=1,npjglo)  ; gphiv(:,:) = dzvar(:,:)
  DEALLOCATE ( dzvar )
  CLOSE(numcoo)

  !! build nc output files
  WRITE(ctaux,'(a,I4.4,a)') 'ECMWF_taux_1d_',iyear,'.'//TRIM(config)//'.nc'
  WRITE(ctauy,'(a,I4.4,a)') 'ECMWF_tauy_1d_',iyear,'.'//TRIM(config)//'.nc'

  jmonth=1
  !! Build dimg file names
  WRITE(ctau ,'(a,I4.4,a,I2.2,a)') 'ECMWF.Y',iyear,'.M',jmonth,'.STRESS.'//TRIM(config)//'.dimg'
  ! WRITE(csst ,'(a,I4.4,a,I2.2,a)') 'REYNOLDS.Y',iyear,'.SST.'//TRIM(config)//'.dimg'

  ! open (and check ?? if they exists )
  irecl=isdirect(ctau)  ; OPEN( numtau,FILE=ctau, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )

  READ(numtau,REC=1) cver, cheader, ii, npiglo, npjglo, npk

  ALLOCATE (v2d(npiglo, npjglo,2), dep(npk) )
  ALLOCATE (z2d(npiglo, npjglo) )
  READ(numtau,REC=1) cver, cheader, ii, npiglo, npjglo, npk, nt, ndim, &
       x1,y1,dx,dy,spval,   &
       (dep(jk),jk=1,npk), &
       timean(1)
  CLOSE(numtau)

  ! Build  cdf files output
  nvar = 1 ! 1 var but many files ... (OK  ... 3 actually )
  ALLOCATE ( typvartaux(nvar), ipktaux(nvar), id_varouttaux(nvar) )
  ALLOCATE ( typvartauy(nvar), ipktauy(nvar), id_varouttauy(nvar) )
  jvar=1
  ipktaux(jvar)      = 1
  typvartaux(jvar)%cname='sozotaux'    ! taux dim 1 of dimgfile
  typvartaux(jvar)%cunits='N/m2'
  typvartaux(jvar)%rmissing_value=0.
  typvartaux(jvar)%valid_min= -0.1
  typvartaux(jvar)%valid_max= 0.1
  typvartaux(jvar)%clong_name='Zonal Wind Stress'
  typvartaux(jvar)%cshort_name='sozotaux'
  typvartaux(jvar)%conline_operation='N/A'
  typvartaux(jvar)%caxis='TYX'

  ipktauy(jvar)      = 1
  typvartauy(jvar)%cname='sometauy'    ! tauy dim 2 of dimgfile
  typvartauy(jvar)%cunits='N/m2'
  typvartauy(jvar)%rmissing_value=0.
  typvartauy(jvar)%valid_min= -0.1
  typvartauy(jvar)%valid_max= 0.1
  typvartauy(jvar)%clong_name='Meridional Wind Stress'
  typvartauy(jvar)%cshort_name='sometauy'
  typvartauy(jvar)%conline_operation='N/A'
  typvartauy(jvar)%caxis='TYX'

  ncouttaux =create(ctaux, 'none',npiglo,npjglo,npk,cdep='deptht' )
  istatus= createvar(ncouttaux ,typvartaux,nvar, ipktaux,id_varouttaux )
  istatus= putheadervar(ncouttaux, 'none', npiglo, npjglo,npk, pnavlon=glam,pnavlat=gphi,pdep=dep )

  ncouttauy =create(ctauy, 'none',npiglo,npjglo,npk,cdep='deptht' )
  istatus= createvar(ncouttauy ,typvartauy,nvar, ipktauy,id_varouttauy )
  istatus= putheadervar(ncouttauy, 'none', npiglo, npjglo,npk, pnavlon=glam,pnavlat=gphi,pdep=dep )

  ! Ready for time loop on month
  icurrday=0
  DO jmonth = 1, 12
     WRITE(ctau,'(a,I4.4,a,I2.2,a)') 'ECMWF.Y',iyear,'.M',jmonth,'.STRESS.'//TRIM(config)//'.dimg'
     irecl=isdirect(ctau)  ; OPEN( numtau,FILE=ctau, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
     READ(numtau,REC=1) cver, cheader, ii, npiglo, npjglo, npk, nt, ndim
     ! loop for days in files
     DO jday=1,nt
        icurrday=icurrday +1
        DO jdim=1,ndim
           irec=1+(jday-1)*ndim +jdim
           READ(numtau,REC=irec) (( v2d(ji,jj,jdim),ji=1,npiglo),jj=1,npjglo) 
        END DO
        ! taux
        istatus = putvar(ncouttaux,id_varouttaux(1),v2d(:,:,1),icurrday,npiglo,npjglo)
        ! tauy
        istatus = putvar(ncouttauy,id_varouttauy(1),v2d(:,:,2),icurrday,npiglo,npjglo)
     END DO ! loop on days
     CLOSE(numtau)
  END DO ! loop on month

  ! update time_counter
  ALLOCATE( timetab (icurrday) )
  timetab=(/(jt,jt=1,icurrday)/)
  istatus=putvar1d(ncouttaux,timetab,icurrday,'T')
  istatus=putvar1d(ncouttauy,timetab,icurrday,'T')
  ! close fluxes files
  istatus=closeout(ncouttaux)
  istatus=closeout(ncouttauy)
  DEALLOCATE (v2d , dep, z2d , timetab)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!     .....    SST    SST     SST     .....
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PRINT *,' Doing SST ...'

  !!  glam gphi are already read ( T point)

  !! build nc output files
  WRITE(csst,'(a,I4.4,a)') 'REYNOLDS_sst_1d_',iyear,'.'//TRIM(config)//'.nc'

  !! Build dimg file names
  WRITE(csstr ,'(a,I4.4,a,I2.2,a)') 'REYNOLDS.Y',iyear,'.SST.'//TRIM(config)//'.dimg'

  ! open (and check ?? if they exists )
  irecl=isdirect(csstr)  ; OPEN( numsst,FILE=csstr, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )

  READ(numsst,REC=1) cver, cheader, ii, npiglo, npjglo, npk, nt

  ALLOCATE (v2d(npiglo, npjglo,nt+2),itime(nt+2), dep(npk) ,timetab(nt), timetag(nt) )
  ALLOCATE (z2d(npiglo, npjglo) ,v2daily(npiglo,npjglo)  )
  READ(numsst,REC=1) cver, cheader, ii, npiglo, npjglo, npk, nt, ndim, &
       x1,y1,dx,dy,spval,   &
       (dep(jk),jk=1,npk), &
       (timetab(jt), jt=1,nt)
       timetag=timetab ! convert to dble precision
       DEALLOCATE(timetab)

  ! Build  cdf files output
  nvar = 1 ! 1 var but many files ... (OK  ... 3 actually )
  ALLOCATE ( typvarsst(nvar), ipksst(nvar), id_varoutsst(nvar) )
  jvar=1
  ipksst(jvar)      = 1
  typvarsst(jvar)%cname='sst'    ! sst dim 1 of dimgfile
  typvarsst(jvar)%cunits='C'
  typvarsst(jvar)%rmissing_value=0.
  typvarsst(jvar)%valid_min= -10.
  typvarsst(jvar)%valid_max= 50.
  typvarsst(jvar)%clong_name='Reynolds SST'
  typvarsst(jvar)%cshort_name='SST'
  typvarsst(jvar)%conline_operation='N/A'
  typvarsst(jvar)%caxis='TYX'

  ncoutsst =create(csst, 'none',npiglo,npjglo,npk,cdep='deptht' )
  istatus= createvar(ncoutsst ,typvarsst,nvar, ipksst,id_varoutsst )
  istatus= putheadervar(ncoutsst, 'none', npiglo, npjglo,npk, pnavlon=glam,pnavlat=gphi,pdep=dep )

  ! We want to interpolate the data for every day. (weekly in the file)
  ! if first day of the file is not 01/01, needs to read previous year
  !  Clipper SST files are not y2k compliant ...
  IF (timetag (1) < 10000 ) THEN 
     timetag(:)=timetag(:)+20000000.
  ELSE
     timetag(:)=timetag(:)+19000000.
  ENDIF
  january1=iyear*10000+01*100+01
  december31=iyear*10000+12*100+31
  jul1=julday(january1)
  jul2=julday(december31)

  itt=0
  IF (jul1 < julday(INT(timetag(1)))  ) THEN
     ! need to read previous year 
     WRITE(csstrp ,'(a,I4.4,a,I2.2,a)') 'REYNOLDS.Y',iyear-1,'.SST.'//TRIM(config)//'.dimg'
     irecl=isdirect(csstrp)  ; OPEN( numsstp,FILE=csstrp, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
     READ(numsstp,REC=1) cver, cheader, ii, npiglo, npjglo, npk, ntp
     ALLOCATE (timetagp (ntp) ,timetab(ntp))
     READ(numsstp,REC=1) cver, cheader, ii, npiglo, npjglo, npk, ntp, ndim, &
          x1,y1,dx,dy,spval,   &
          (dep(jk),jk=1,npk), &
          (timetab(jt), jt=1,ntp)
         timetagp=timetab 
         DEALLOCATE(timetab)
     IF (timetagp (1) < 10000 ) THEN 
        timetagp(:)=timetagp(:)+20000000.
     ELSE
        timetagp(:)=timetagp(:)+19000000.
     ENDIF
     IF ( julday(INT(timetagp (ntp)))  <= jul1 ) THEN
        !read ntp sst as 1 data
        itt = itt +1
        READ(numsstp,REC=ntp+1) (( v2d(ji,jj,itt),ji=1,npiglo),jj=1,npjglo)
        itime(itt)=julday( INT(timetagp(ntp)) )
     ELSE IF  ( julday(INT(timetagp (ntp-1))  ) <= jul1 ) THEN
        itt = itt +1
        READ(numsstp,REC=ntp) (( v2d(ji,jj,itt),ji=1,npiglo),jj=1,npjglo)
        itime(itt)=julday(INT(timetagp(ntp-1)) )
     ELSE IF  ( julday(INT(timetagp (ntp-2) )) <= jul1 ) THEN
        itt = itt +1
        READ(numsstp,REC=ntp-1) (( v2d(ji,jj,itt),ji=1,npiglo),jj=1,npjglo)
        itime(itt)=julday(INT(timetagp(ntp-2)) )
     ELSE
        PRINT *,' Something is wrong in previous file SST ' ; STOP
     ENDIF
  ENDIF
  DO jt=1,nt
     itt = itt +1
     READ(numsst,REC=jt+1) (( v2d(ji,jj,itt),ji=1,npiglo),jj=1,npjglo)
     itime(itt)=julday(INT(timetag(jt)) )
  END DO

  IF ( jul2 > julday(INT(timetag(nt))) ) THEN
     ! need to read next year
     IF ( iyear == 2000 ) THEN   ! persistance ...
        itt=itt+1 ;  v2d(:,:,itt)= v2d(:,:,itt-1)  ; itime(itt)=jul2
     ELSE
        WRITE(csstrn ,'(a,I4.4,a,I2.2,a)') 'REYNOLDS.Y',iyear+1,'.SST.'//TRIM(config)//'.dimg'
        irecl=isdirect(csstrn)  ; OPEN( numsstn,FILE=csstrn, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=irecl )
        READ(numsstn,REC=1) cver, cheader, ii, npiglo, npjglo, npk, ntn
        ALLOCATE (timetagn (ntn) ,timetab(ntn))
        READ(numsstn,REC=1) cver, cheader, ii, npiglo, npjglo, npk, ntn, ndim, &
             x1,y1,dx,dy,spval,   &
             (dep(jk),jk=1,npk), &
             (timetab(jt), jt=1,ntn)
        timetagn=timetab
        DEALLOCATE( timetab)
        IF (timetagn (1) < 10000 ) THEN 
           timetagn(:)=INT(timetagn(:))+20000000
        ELSE
           timetagn(:)=INT(timetagn(:))+19000000
        ENDIF

        IF ( julday(INT(timetagn (1) )) >= jul2 ) THEN
           !read 1 sst as 1 data
           itt = itt +1
           READ(numsstn,REC=2) (( v2d(ji,jj,itt),ji=1,npiglo),jj=1,npjglo)
           itime(itt)=julday(INT(timetagn(1)) )
        ELSE IF  ( julday(INT(timetagn (2)) ) >= jul2 ) THEN
           itt = itt +1
           READ(numsstn,REC=3) (( v2d(ji,jj,itt),ji=1,npiglo),jj=1,npjglo)
           itime(itt)=julday(INT( timetagn(2)) )
        ELSE IF  ( julday(INT(timetagn (3))) >= jul2 ) THEN
           itt = itt +1
           READ(numsstn,REC=4) (( v2d(ji,jj,itt),ji=1,npiglo),jj=1,npjglo)
           itime(itt)=julday(INT(timetagn(3)) )
        ELSE
           PRINT *,' Something is wrong in next file SST ' ; STOP
        ENDIF
     ENDIF
  ENDIF
     ntime=itt

     icurrday=0
     ii1=1 ; ii2 = 2 ; id1=itime(ii1) ; id2=itime(ii2) 
     DO jul = jul1, jul2
        icurrday=icurrday + 1
        IF ( jul > id2 ) THEN
           ii1=ii1+1 ; ii2=ii2+1 ;  id1=itime(ii1) ; id2=itime(ii2)
        ENDIF
        v2daily(:,:)=FLOAT((jul - id1 ))/(FLOAT(id2-id1))*(v2d(:,:,ii2) - v2d(:,:,ii1) ) +  v2d(:,:,ii1)
        istatus = putvar(ncoutsst,id_varoutsst(1),v2daily(:,:),icurrday,npiglo,npjglo)
     END DO

     ! update time_counter
     ALLOCATE( timetab (icurrday) )
     timetab=(/(jt,jt=1,icurrday)/)
     istatus=putvar1d(ncoutsst,timetab,icurrday,'T')
     ! close fluxes files
     istatus=closeout(ncoutsst)
     istatus=closeout(ncoutsst)
     DEALLOCATE (v2d , dep, z2d )

   CONTAINS
     INTEGER FUNCTION isdirect(clname)
!!!                     FUNCTION ISDIRECT
!!!                     *****************
!!!
!!!    PURPOSE : This integer function returns the record length if clname
!!!              is a valid dimg file, it returns 0 either.
!!!
!!!    METHOD : Open the file and look for the key characters (@!01) for
!!!             identification.
!!!
!!!    AUTHOR : Jean-Marc Molines (Apr. 1998)
!!! -------------------------------------------------------------------------
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(in) ::  clname
       CHARACTER(LEN=4)  ::  cver
       CHARACTER(LEN=2048) ::  clheader
       !
       INTEGER :: irecl

       !
       OPEN(100,FILE=clname, FORM   ='UNFORMATTED', ACCESS ='DIRECT', RECL   =88)
       READ(100,REC=1) cver ,clheader,irecl
       CLOSE(100)
       !
       IF (cver ==  '@!01' ) THEN
          isdirect=irecl
       ELSE
          isdirect=0
       END IF
       !
     END FUNCTION isdirect

     FUNCTION julday(kdastp)
       !! ------------------------------------------------------------------
       !!          ***        FUNCTION JULDAY    ***
       !!
       !!   Purpose:   This routine returns the julian day number which begins at noon
       !!         of the calendar date specified by month kmm, day kid, and year kiyyy.
       !!         positive year signifies a.d.; negative, b.c.  (remember that the
       !!         year after 1 b.c. was 1 a.d.)
       !!         routine handles changeover to gregorian calendar on oct. 15, 1582.
       !!
       !!   Method:  This routine comes directly from the Numerical Recipe Book,
       !!           press et al., numerical recipes, cambridge univ. press, 1986.
       !!
       !!   Arguments:
       !!     kdastp  : OPA date yyyymmdd (instead of kmm kid kiyyy)
       !!     kmm     : input, corresponding month
       !!     kid     : input, corresponding day
       !!     kiyyy   : input, corresponding year, positive IF a.d, negative b.c.
       !!      
       !!     
       !!   history
       !!     1998: J.M. Molines for the Doctor form. 
       !!     2007 : J.M. Molines in F90
       !! -----------------------------------------------------------------
       !  *  Declarations
       !
       INTEGER :: julday, kiyyy,kid,kmm
       INTEGER, INTENT(in)  ::kdastp
       !  * Local 
       INTEGER, PARAMETER ::jpgreg=15+31*(10+12*1582)
       INTEGER  :: iy, im, ia
       ! ... Year 0 never existed ...
       kiyyy=kdastp/10000
       kmm=(kdastp - kiyyy*10000)/100
       kid= kdastp - kiyyy*10000 - kmm*100
       IF (kiyyy == 0) STOP 101
       !
       IF (kiyyy < 0) kiyyy = kiyyy + 1
       IF (kmm > 2) THEN
          iy = kiyyy
          im = kmm + 1
       ELSE
          iy = kiyyy - 1
          im = kmm + 13
       END IF
       !
       julday = INT(365.25*iy) + INT(30.6001*im) + kid + 1720995 
       IF (kid+31*(kmm+12*kiyyy).GE.jpgreg) THEN
          ia = INT(0.01*iy)
          julday = julday + 2 - ia + INT(0.25*ia) 
       END IF
     END FUNCTION JULDAY
   END PROGRAM cdfflxconv
