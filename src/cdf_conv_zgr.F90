PROGRAM cdf_conv_zgr
  !!======================================================================
  !!                     ***  PROGRAM  cdf_conv_zgr  ***
  !!=====================================================================
  !!  ** Purpose : Convert mesh_zgr file from 3.0 Mercator style to 
  !!               3.6 Drakkar style.
  !!
  !!  ** Method  : Mercator style use 2D e3t_ps and e3w_ps + 1d e3t_0, e3w_0
  !!               Drakkar style use 3D e3t_0, e3u_0 e3v_0 e3w_0 + all
  !!               standard features of standard 3.6 ( It corresponds to
  !!               mesh creation with nn_msh=6 in NEMO namelist).
  !!
  !!  ** Remark  : This code is not using cdfio module and has netcdf 
  !!               primitive directly on it...
  !!
  !! History :  4.0  : 10/2018  : J.M. Molines : 
  !!----------------------------------------------------------------------

  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2018
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class file_operations
  !!----------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE

  INTEGER(KIND=4) ::  ji,jj,jk
  INTEGER(KIND=4) :: narg, ijarg, icky
  INTEGER(KIND=4) :: ncid, id, ierr, ik
  INTEGER(KIND=4) :: ncido
  INTEGER(KIND=4) :: id_lon, id_lat, id_lev, id_tim
  INTEGER(KIND=4) :: id_mbt, id_hdt, id_hdw, id_gdt, id_gdw
  INTEGER(KIND=4) :: id_e3t1d, id_e3w1d
  INTEGER(KIND=4) :: idx, idy, idz, idt
  INTEGER(KIND=4) :: id_e3t, id_e3u, id_e3v, id_e3w
  INTEGER(KIND=4) :: npiglo, npjglo, npk, npt
  INTEGER(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: mbathy

  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: e3t, e3x,e3w
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: e3t_ps, e3w_ps, v2d
  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: e3t_1d, e3w_1d, v1d, v1dt

  CHARACTER(LEN=80) :: cf_zgr_in
  CHARACTER(LEN=80) :: cf_zgr_out
  CHARACTER(LEN=80) :: cldum

  LOGICAL           :: lnc4=.FALSE.
  !!----------------------------------------------------------------------

  narg=iargc() 
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdf_conv_zgr -i <IN-file> -o <OUT-file> [-nc4] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Convert mesh_zgr files whith 2D +1D variables, into a full 3D '
     PRINT *,'       mesh_zgr file compliant with NEMO_3.6 format (names)' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        -i IN-file : name of the input mesh_zgr file to convert.'
     PRINT *,'        -o OUT-file : name of the resulting converted mesh_zgr.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        -nc4 : use netcdf4 with chunking and compression'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'         none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'         Output file name is specified on the command line.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      ' 
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-i'   ) ; CALL getarg(ijarg, cf_zgr_in  ) ; ijarg=ijarg+1
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_zgr_out ) ; ijarg=ijarg+1
        ! option
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  ! open file
  ierr=NF90_OPEN(cf_zgr_in, NF90_NOWRITE, ncid)
  ! read dimensions
  ierr=NF90_INQ_DIMID(ncid,"x",id) ; ierr=NF90_INQUIRE_DIMENSION(ncid,id,len=npiglo)
  ierr=NF90_INQ_DIMID(ncid,"y",id) ; ierr=NF90_INQUIRE_DIMENSION(ncid,id,len=npjglo)
  ierr=NF90_INQ_DIMID(ncid,"z",id) ; ierr=NF90_INQUIRE_DIMENSION(ncid,id,len=npk)
  ierr=NF90_INQ_DIMID(ncid,"t",id) ; ierr=NF90_INQUIRE_DIMENSION(ncid,id,len=npt)

  ! allocate space
  ALLOCATE( e3t(npiglo, npjglo, npk), e3x(npiglo, npjglo, npk) )
  ALLOCATE( e3w(npiglo, npjglo, npk) )
  ALLOCATE( e3t_ps(npiglo, npjglo) , mbathy(npiglo, npjglo),e3w_ps(npiglo, npjglo), v2d(npiglo, npjglo))
  ALLOCATE( e3t_1d(npk), e3w_1d(npk), v1d(npk), v1dt(npt) )

  ! read e3t_ps
  ierr=NF90_INQ_VARID(ncid,"e3t_ps",id) ; ierr=NF90_GET_VAR(ncid, id, e3t_ps,start=(/1,1,1/), count=(/npiglo,npjglo,1/) )
  PRINT * ,NF90_STRERROR(ierr)
  ! read e3w_ps
  ierr=NF90_INQ_VARID(ncid,"e3w_ps",id) ; ierr=NF90_GET_VAR(ncid, id, e3w_ps,start=(/1,1,1/), count=(/npiglo,npjglo,1/) )
  PRINT * ,NF90_STRERROR(ierr)
  ! read e3t_1d
  ierr=NF90_INQ_VARID(ncid,"e3t_0",id) ; ierr=NF90_GET_VAR(ncid, id, e3t_1d,start=(/1,1/), count=(/npk,1/) )
  PRINT * ,NF90_STRERROR(ierr), e3t_1d
  ! read e3w_1d
  ierr=NF90_INQ_VARID(ncid,"e3w_0",id) ; ierr=NF90_GET_VAR(ncid, id, e3w_1d,start=(/1,1/), count=(/npk,1/) )
  PRINT * ,NF90_STRERROR(ierr), e3w_1d
  ! read mbathy : give the k-index for the last point in the sea
  ierr=NF90_INQ_VARID(ncid,"mbathy",id) ; ierr=NF90_GET_VAR(ncid, id, mbathy,start=(/1,1,1/), count=(/npiglo,npjglo,1/) )
  PRINT * ,NF90_STRERROR(ierr), MAXVAL(mbathy)

  ! open the output file
#if defined key_netcdf4
  IF ( lnc4 ) THEN
     ierr=NF90_CREATE(cf_zgr_out,cmode=or(NF90_CLOBBER,NF90_NETCDF4     ), ncid=ncido)
  ELSE
     ierr=NF90_CREATE(cf_zgr_out,cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncido)
  ENDIF
#else
  ierr=NF90_CREATE(cf_zgr_out,cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncido)
#endif
  ierr=NF90_DEF_DIM(ncido,"x",npiglo, idx)
  ierr=NF90_DEF_DIM(ncido,"y",npjglo, idy)
  ierr=NF90_DEF_DIM(ncido,"z",npk, idz)
  ierr=NF90_DEF_DIM(ncido,"t",NF90_UNLIMITED, idt)

#if defined key_netcdf4
  IF ( lnc4 ) THEN
     ierr=NF90_DEF_VAR(ncido,"nav_lon",NF90_FLOAT,(/idx,idy/), id_lon,chunksizes=(/npiglo,1/),deflate_level=1 )
     ierr=NF90_DEF_VAR(ncido,"nav_lat",NF90_FLOAT,(/idx,idy/), id_lat,chunksizes=(/npiglo,1/),deflate_level=1 )
     ierr=NF90_DEF_VAR(ncido,"nav_lev",NF90_FLOAT,(/idz/)    , id_lev)
     ierr=NF90_DEF_VAR(ncido,"time_counter",NF90_DOUBLE,(/idt/), id_tim)
     !
     icky=MAX(npiglo/30,10)
     ierr=NF90_DEF_VAR(ncido,"e3t_0",NF90_DOUBLE,(/idx,idy,idz,idt/), id_e3t,chunksizes=(/npiglo,icky,1,1/),deflate_level=1 )
     ierr=NF90_DEF_VAR(ncido,"e3u_0",NF90_DOUBLE,(/idx,idy,idz,idt/), id_e3u,chunksizes=(/npiglo,icky,1,1/),deflate_level=1 )
     ierr=NF90_DEF_VAR(ncido,"e3v_0",NF90_DOUBLE,(/idx,idy,idz,idt/), id_e3v,chunksizes=(/npiglo,icky,1,1/),deflate_level=1 )
     ierr=NF90_DEF_VAR(ncido,"e3w_0",NF90_DOUBLE,(/idx,idy,idz,idt/), id_e3w,chunksizes=(/npiglo,icky,1,1/),deflate_level=1 )
     ierr=NF90_DEF_VAR(ncido,"mbathy",NF90_SHORT,(/idx,idy,idt/)    , id_mbt,chunksizes=(/npiglo,icky,1/)  ,deflate_level=1 )
     ierr=NF90_DEF_VAR(ncido,"hdept",NF90_FLOAT,(/idx,idy,idt/)     , id_hdt,chunksizes=(/npiglo,icky,1/)  ,deflate_level=1 )
     ierr=NF90_DEF_VAR(ncido,"hdepw",NF90_FLOAT,(/idx,idy,idt/)     , id_hdw,chunksizes=(/npiglo,icky,1/)  ,deflate_level=1 )
     ierr=NF90_DEF_VAR(ncido,"gdept_1d",NF90_DOUBLE,(/idz,idt/)     , id_gdt)
     ierr=NF90_DEF_VAR(ncido,"gdepw_1d",NF90_DOUBLE,(/idz,idt/)     , id_gdw)
     ierr=NF90_DEF_VAR(ncido,"e3t_1d",NF90_DOUBLE,(/idz,idt/)       , id_e3t1d)
     ierr=NF90_DEF_VAR(ncido,"e3w_1d",NF90_DOUBLE,(/idz,idt/)       , id_e3w1d)
  ELSE
     ierr=NF90_DEF_VAR(ncido,"nav_lon",NF90_FLOAT,(/idx,idy/), id_lon)
     ierr=NF90_DEF_VAR(ncido,"nav_lat",NF90_FLOAT,(/idx,idy/), id_lat)
     ierr=NF90_DEF_VAR(ncido,"nav_lev",NF90_FLOAT,(/idz/)    , id_lev)
     ierr=NF90_DEF_VAR(ncido,"time_counter",NF90_DOUBLE,(/idt/), id_tim)
     !
     ierr=NF90_DEF_VAR(ncido,"e3t_0",NF90_DOUBLE,(/idx,idy,idz,idt/), id_e3t)
     ierr=NF90_DEF_VAR(ncido,"e3u_0",NF90_DOUBLE,(/idx,idy,idz,idt/), id_e3u)
     ierr=NF90_DEF_VAR(ncido,"e3v_0",NF90_DOUBLE,(/idx,idy,idz,idt/), id_e3v)
     ierr=NF90_DEF_VAR(ncido,"e3w_0",NF90_DOUBLE,(/idx,idy,idz,idt/), id_e3w)
     ierr=NF90_DEF_VAR(ncido,"mbathy",NF90_SHORT,(/idx,idy,idt/), id_mbt)
     ierr=NF90_DEF_VAR(ncido,"hdept",NF90_FLOAT,(/idx,idy,idt/), id_hdt)
     ierr=NF90_DEF_VAR(ncido,"hdepw",NF90_FLOAT,(/idx,idy,idt/), id_hdw)
     ierr=NF90_DEF_VAR(ncido,"gdept_1d",NF90_DOUBLE,(/idz,idt/), id_gdt)
     ierr=NF90_DEF_VAR(ncido,"gdepw_1d",NF90_DOUBLE,(/idz,idt/), id_gdw)
     ierr=NF90_DEF_VAR(ncido,"e3t_1d",NF90_DOUBLE,(/idz,idt/), id_e3t1d)
     ierr=NF90_DEF_VAR(ncido,"e3w_1d",NF90_DOUBLE,(/idz,idt/), id_e3w1d)
  ENDIF
#else
  ierr=NF90_DEF_VAR(ncido,"nav_lon",NF90_FLOAT,(/idx,idy/), id_lon)
  ierr=NF90_DEF_VAR(ncido,"nav_lat",NF90_FLOAT,(/idx,idy/), id_lat)
  ierr=NF90_DEF_VAR(ncido,"nav_lev",NF90_FLOAT,(/idz/)    , id_lev)
  ierr=NF90_DEF_VAR(ncido,"time_counter",NF90_DOUBLE,(/idt/), id_tim)
  !
  ierr=NF90_DEF_VAR(ncido,"e3t_0",NF90_DOUBLE,(/idx,idy,idz,idt/), id_e3t)
  ierr=NF90_DEF_VAR(ncido,"e3u_0",NF90_DOUBLE,(/idx,idy,idz,idt/), id_e3u)
  ierr=NF90_DEF_VAR(ncido,"e3v_0",NF90_DOUBLE,(/idx,idy,idz,idt/), id_e3v)
  ierr=NF90_DEF_VAR(ncido,"e3w_0",NF90_DOUBLE,(/idx,idy,idz,idt/), id_e3w)
  ierr=NF90_DEF_VAR(ncido,"mbathy",NF90_SHORT,(/idx,idy,idt/), id_mbt)
  ierr=NF90_DEF_VAR(ncido,"hdept",NF90_FLOAT,(/idx,idy,idt/), id_hdt)
  ierr=NF90_DEF_VAR(ncido,"hdepw",NF90_FLOAT,(/idx,idy,idt/), id_hdw)
  ierr=NF90_DEF_VAR(ncido,"gdept_1d",NF90_DOUBLE,(/idz,idt/), id_gdt)
  ierr=NF90_DEF_VAR(ncido,"gdepw_1d",NF90_DOUBLE,(/idz,idt/), id_gdw)
  ierr=NF90_DEF_VAR(ncido,"e3t_1d",NF90_DOUBLE,(/idz,idt/), id_e3t1d)
  ierr=NF90_DEF_VAR(ncido,"e3w_1d",NF90_DOUBLE,(/idz,idt/), id_e3w1d)
#endif
  ierr=NF90_ENDDEF(ncido)

  ! fill the dimension variables :

  ierr=NF90_INQ_VARID(ncid,"nav_lon",id) ; ierr=NF90_GET_VAR(ncid,id,v2d) ; ierr=NF90_PUT_VAR(ncido,id_lon,v2d)
  ierr=NF90_INQ_VARID(ncid,"nav_lat",id) ; ierr=NF90_GET_VAR(ncid,id,v2d) ; ierr=NF90_PUT_VAR(ncido,id_lat,v2d)
  ierr=NF90_INQ_VARID(ncid,"nav_lev",id) ; ierr=NF90_GET_VAR(ncid,id,v1d) ; ierr=NF90_PUT_VAR(ncido,id_lev,v1d)
  ierr=NF90_INQ_VARID(ncid,"time_counter",id) ; ierr=NF90_GET_VAR(ncid,id,v1dt) ; ierr=NF90_PUT_VAR(ncido,id_tim,v1dt)

  ! now depth and metrics
  ierr=NF90_INQ_VARID(ncid,"hdept",id) ; ierr=NF90_GET_VAR(ncid,id,     v2d,start=(/1,1,1/), count=(/npiglo,npjglo,1/) ) 
  ;                                    ; ierr=NF90_PUT_VAR(ncido,id_hdt,v2d,start=(/1,1,1/), count=(/npiglo,npjglo,1/) )
  ierr=NF90_INQ_VARID(ncid,"hdepw",id) ; ierr=NF90_GET_VAR(ncid,id,     v2d,start=(/1,1,1/), count=(/npiglo,npjglo,1/) ) 
  ;                                    ; ierr=NF90_PUT_VAR(ncido,id_hdw,v2d,start=(/1,1,1/), count=(/npiglo,npjglo,1/) )

  ierr=NF90_INQ_VARID(ncid,"gdept_0",id) ; ierr=NF90_GET_VAR(ncid,id,     v1d,start=(/1,1/), count=(/npk,1/) ) 
  ;                                      ; ierr=NF90_PUT_VAR(ncido,id_gdt,v1d,start=(/1,1/), count=(/npk,1/) )
  ierr=NF90_INQ_VARID(ncid,"gdepw_0",id) ; ierr=NF90_GET_VAR(ncid,id,     v1d,start=(/1,1/), count=(/npk,1/) ) 
  ;                                      ; ierr=NF90_PUT_VAR(ncido,id_gdw,v1d,start=(/1,1/), count=(/npk,1/) )

  ierr=NF90_INQ_VARID(ncid,"e3t_0",id) ; ierr=NF90_GET_VAR(ncid,id,       v1d,start=(/1,1/), count=(/npk,1/) ) 
  ;                                    ; ierr=NF90_PUT_VAR(ncido,id_e3t1d,v1d,start=(/1,1/), count=(/npk,1/) )
  ierr=NF90_INQ_VARID(ncid,"e3w_0",id) ; ierr=NF90_GET_VAR(ncid,id,       v1d,start=(/1,1/), count=(/npk,1/) ) 
  ;                                    ; ierr=NF90_PUT_VAR(ncido,id_e3w1d,v1d,start=(/1,1/), count=(/npk,1/) )

  ierr=NF90_PUT_VAR(ncido,id_mbt, mbathy, start=(/1,1,1/), count=(/npiglo,npjglo,1/) )

  ! build the 3d e3t
  ! 1 init
  DO jk=1,npk
     e3t(:,:,jk) = e3t_1d(jk)
  ENDDO
  ! 2 ps
  DO jj=1,npjglo
     DO ji=1,npiglo
        ik=MAX(mbathy(ji,jj),1)
        e3t(ji,jj,ik)=MAX(e3t_ps(ji,jj),e3t_1d(1))
     ENDDO
  ENDDO
  ierr=NF90_PUT_VAR(ncido,id_e3t,e3t,start=(/1,1,1,1/), count=(/npiglo,npjglo, npk,1/) )
  ! Now build e3u 
  e3x=100.d0
  DO jk=1, npk
     DO jj=1,npjglo
        DO ji=1, npiglo-1
           e3x(ji,jj,jk) = MIN(e3t(ji,jj,jk), e3t(ji+1,jj,jk))
        ENDDO
     ENDDO
  ENDDO
  ! putvar ...
  ierr=NF90_PUT_VAR(ncido,id_e3u,e3x,start=(/1,1,1,1/), count=(/npiglo,npjglo, npk,1/) )

  ! Now build e3v 
  e3x=100.d0
  DO jk=1, npk
     DO jj=1,npjglo-1
        DO ji=1, npiglo
           e3x(ji,jj,jk) = MIN(e3t(ji,jj,jk), e3t(ji,jj+1,jk))
        ENDDO
     ENDDO
  ENDDO

  ! putvar ...
  ierr=NF90_PUT_VAR(ncido,id_e3v,e3x,start=(/1,1,1,1/), count=(/npiglo,npjglo, npk,1/) )

  ! build the 3d e3w
  ! 1 init
  DO jk=1,npk
     e3w(:,:,jk) = e3w_1d(jk)
  ENDDO
  ! 2 ps
  DO jj=1,npjglo
     DO ji=1,npiglo
        ik=MAX(mbathy(ji,jj),1)
        e3w(ji,jj,ik)=MAX(e3w_ps(ji,jj),e3w_1d(1))
     ENDDO
  ENDDO
  ierr=NF90_PUT_VAR(ncido,id_e3w,e3w,start=(/1,1,1,1/), count=(/npiglo,npjglo, npk,1/) )

  ierr=NF90_CLOSE(ncido)
  ierr=NF90_CLOSE(ncid)


END PROGRAM cdf_conv_zgr


