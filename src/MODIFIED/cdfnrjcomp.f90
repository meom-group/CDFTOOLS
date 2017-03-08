PROGRAM cdfnrjcomp
  !!======================================================================
  !!                     ***  PROGRAM  cdfnrjcomp  ***
  !!=====================================================================
  !!  ** Purpose : Compute the terms for energy components
  !!               (Mean Kinetic Energy, Eddy Kinetic Energy,
  !!                Mean Potential Energy, Eddy Potential Energy )
  !!               compute : tbar, ubar, vbar, anotsqrt, anousqrt, anovsqrt
  !!
  !! History : 2.1  : 02/2008  : A. Melet     : Original code
  !!         :  4.0  : 03/2017  : J.M. Molines  
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class energy_diagnostics
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: ji, jj, jk           ! dummy loop index
  INTEGER(KIND=4)                            :: npiglo, npjglo       ! domain size
  INTEGER(KIND=4)                            :: npk, npt             ! domain size
  INTEGER(KIND=4)                            :: narg, iargc          ! browse line
  INTEGER(KIND=4)                            :: ierr                 ! error status
  INTEGER(KIND=4)                            :: ncout                ! ncid of output file
  INTEGER(KIND=4), DIMENSION(6)              :: ipk, id_varout       ! level and varid's of output var

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: un, vn, u2n, v2n
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: tn, t2n, anotsqrt
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: umask, vmask
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: anousqrt, anovsqrt 
  REAL(KIND=4), DIMENSION(1)                 :: tim                   ! time counter

  CHARACTER(LEN=256)                         :: cf_in                 ! input filename
  CHARACTER(LEN=256)                         :: cf_out='nrjcomp.nc'   ! output file name
  TYPE (variable), DIMENSION(6)              :: stypvar               ! structure for attibutes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  !!
  narg = iargc()
  IF ( narg /= 1 ) THEN
     PRINT *,' usage : cdfnrjcomp IN-file'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute contributing terms of the energy equation at T-points.'
     PRINT *,'       Input file contains mean values processed by cdfmoyuvwt.' 
     PRINT *,'       The means must have been computed on long enough period'
     PRINT *,'       for the statistics to be meaningful'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       IN-file   : netcdf file produced by cdfmoyuvwt.' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         all variables are located at T point.'
     PRINT *,'         variables : tbar : mean temperature '
     PRINT *,'                     ubar : mean zonal velocity'
     PRINT *,'                     vbar : mean meridional velocity'
     PRINT *,'                     anotsqrt : mean squared temperature anomaly'
     PRINT *,'                     anousqrt : mean squared zonal velocity anomaly'
     PRINT *,'                     anovsqrt : mean squared meridional velocity anomaly'
     STOP
  ENDIF

  CALL getarg(1, cf_in)

  IF ( chkfile(cf_in) ) STOP ! missing file

  npiglo = getdim(cf_in,cn_x)
  npjglo = getdim(cf_in,cn_y)
  npk    = getdim(cf_in,cn_z)
  npt    = getdim(cf_in,cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt
 
  ! define new variables for output 
  ipk(:)                 = npk  
  stypvar(1)%cname       = 'tbar'
  stypvar(1)%clong_name  = 'temporal mean of the temperature on T point'
  stypvar(1)%cshort_name = 'tbar'

  stypvar(2)%cname       = 'ubar'
  stypvar(2)%clong_name  = 'temporal mean of the zonal velocity on T point'
  stypvar(2)%cshort_name = 'ubar'
 
  stypvar(3)%cname       = 'vbar'
  stypvar(3)%clong_name  = 'temporal mean of the meridional velocity on T point'
  stypvar(3)%cshort_name = 'vbar'
  
  stypvar(4)%cname       = 'anotsqrt'
  stypvar(4)%clong_name  = 'temporal mean of the square of the temperature anomaly on T point (*1000)'
  stypvar(4)%cshort_name = 'anotsqrt'

  stypvar(5)%cname       = 'anousqrt'
  stypvar(5)%clong_name  = 'temporal mean of the square of the zonal speed anomaly on T point (*1000)'
  stypvar(5)%cshort_name = 'anousqrt'

  stypvar(6)%cname       = 'anovsqrt'
  stypvar(6)%clong_name  = 'temporal mean of the square of the meridional speed anomaly on T point (*1000)'
  stypvar(6)%cshort_name = 'anovsqrt'

  stypvar%cunits            = ' '
  stypvar%rmissing_value    = 0.
  stypvar%valid_min         = -1000.
  stypvar%valid_max         = 1000.
  stypvar%conline_operation = 'N/A'
  stypvar%caxis             = 'TZYX'
  
  ! create output fileset
  ncout = create      (cf_out, cf_in,   npiglo, npjglo, npk       )
  ierr  = createvar   (ncout,  stypvar, 6,      ipk,    id_varout )
  ierr  = putheadervar(ncout,  cf_in,   npiglo, npjglo, npk       )

  ! Allocate the memory
  ALLOCATE ( un(npiglo,npjglo), vn(npiglo,npjglo)  )
  ALLOCATE ( umask(npiglo,npjglo), vmask(npiglo,npjglo) )
  ALLOCATE ( u2n(npiglo,npjglo), v2n(npiglo,npjglo)  )
  ALLOCATE ( anousqrt(npiglo,npjglo), anovsqrt(npiglo,npjglo)  )
  ALLOCATE ( tn(npiglo,npjglo), t2n(npiglo,npjglo)  )
  ALLOCATE ( anotsqrt(npiglo,npjglo) )

  tim  = getvar1d(cf_in,cn_vtimec, npt     )
  ierr = putvar1d(ncout, tim,      npt, 'T')
     
  DO jk=1, npk
     PRINT *,'            level ',jk
     anousqrt(:,:) = 0.0
     anovsqrt(:,:) = 0.0      
     anotsqrt(:,:) = 0.0
     
     un(:,:)  = getvar(cf_in, 'ubar',  jk, npiglo, npjglo, ktime=1)
     vn(:,:)  = getvar(cf_in, 'vbar',  jk, npiglo, npjglo, ktime=1)
     u2n(:,:) = getvar(cf_in, 'u2bar', jk, npiglo, npjglo, ktime=1)
     v2n(:,:) = getvar(cf_in, 'v2bar', jk, npiglo, npjglo, ktime=1)
     tn(:,:)  = getvar(cf_in, 'tbar',  jk, npiglo, npjglo, ktime=1)
     t2n(:,:) = getvar(cf_in, 't2bar', jk, npiglo, npjglo, ktime=1)

     ! compute the mask
     DO jj = 2, npjglo
        DO ji = 2, npiglo
           umask(ji,jj) = 0.
           vmask(ji,jj) = 0.
           umask(ji,jj) = un(ji,jj)*un(ji-1,jj) 
           vmask(ji,jj) = vn(ji,jj)*vn(ji,jj-1)
           IF (umask(ji,jj) /= 0.) umask(ji,jj)=1.
           IF (vmask(ji,jj) /= 0.) vmask(ji,jj)=1.
        ENDDO
     ENDDO

     DO jj = 2, npjglo  
        DO ji = 2, npiglo    ! vector opt.
           anotsqrt(ji,jj) = 1000. * ( t2n(ji,jj) - tn(ji,jj) * tn(ji,jj) )
           anousqrt(ji,jj) = 1000./2. * umask(ji,jj)*( ( u2n(ji,jj) - un(ji,jj)*un(ji,jj) ) &
                &                     + ( u2n(ji-1,jj) - un(ji-1,jj)*un(ji-1,jj) ) )       
                
           anovsqrt(ji,jj) = 1000./2. * vmask(ji,jj)*( ( v2n(ji,jj) - vn(ji,jj)*vn(ji,jj) ) &
                &                     + ( v2n(ji,jj-1) - vn(ji,jj)*vn(ji,jj-1) ) ) 
        END DO
     END DO
     ! 
     ierr = putvar(ncout, id_varout(1), tn,       jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(2), un,       jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(3), vn,       jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(4), anotsqrt, jk, npiglo, npjglo, ktime=1)
     ierr = putvar(ncout, id_varout(5), anousqrt, jk, npiglo, npjglo, ktime=1) 
     ierr = putvar(ncout, id_varout(6), anovsqrt, jk, npiglo, npjglo, ktime=1)
  END DO

  ierr = closeout(ncout)

END PROGRAM cdfnrjcomp

