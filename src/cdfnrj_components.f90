PROGRAM cdfnrj_components
  !!======================================================================
  !!                     ***  PROGRAM  cdfnrj_components  ***
  !!=====================================================================
  !!  ** Purpose : Compute the terms for energy components
  !!               (Mean Kinetic Energy, Eddy Kinetic Energy,
  !!                Mean Potential Energy, Eddy Potential Energy )
  !!               compute : tbar, ubar, vbar, anotsqrt, anousqrt, anovsqrt
  !!
  !! History : 2.1  : 02/2008  : A. Melet     : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic. 
  !!         : 4.0  : 03/2017  : J.M. Molines  
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

  INTEGER(KIND=4)                            :: ji, jj, jk ,jt       ! dummy loop index
  INTEGER(KIND=4)                            :: npiglo, npjglo       ! domain size
  INTEGER(KIND=4)                            :: npk, npt             ! domain size
  INTEGER(KIND=4)                            :: narg, iargc, ijarg   ! browse line
  INTEGER(KIND=4)                            :: ierr                 ! error status
  INTEGER(KIND=4)                            :: ncout                ! ncid of output file
  INTEGER(KIND=4), DIMENSION(6)              :: ipk, id_varout       ! level and varid's of output var

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: un, vn, u2n, v2n
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: tn, t2n, anotsqrt
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: umask, vmask
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: anousqrt, anovsqrt 

  REAL(KIND=8), DIMENSION(1)                 :: dtim                  ! time counter

  CHARACTER(LEN=256)                         :: cf_in                 ! input filename
  CHARACTER(LEN=256)                         :: cf_out='nrjcomp.nc'   ! output file name
  CHARACTER(LEN=256)                         :: cldum                 ! working char variable

  TYPE (variable), DIMENSION(6)              :: stypvar               ! structure for attibutes

  LOGICAL                                    :: lnc4 = .FALSE.        ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  !!
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfnrj_components -f UVWT-file [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute contributing terms of the energy equation at T-points.'
     PRINT *,'       Input file contains mean values processed by cdfuvwt.' 
     PRINT *,'       The means must have been computed on long enough period for the'
     PRINT *,'       statistics to be meaningful.'
     PRINT *,'       Note : this program was formerly named cdfnrjcomp. It needs some'
     PRINT *,'       additional cleaning and revision. Use with care!'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f UVWT-file: netcdf file produced by cdfuvwt.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file]: specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ] : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'             This option is effective only if cdftools are compiled with'
     PRINT *,'             a netcdf library supporting chunking and deflation.'
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
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfuvwt, cdfnrj_bti, cdfnrj_bci, cdfnrj_transfert but also'
     PRINT *,'       cdfeke, cdfmoy, cdfstdevt etc...'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_in ) ; ijarg=ijarg+1
        ! option
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 99
     END SELECT
  ENDDO

  IF ( chkfile(cf_in) ) STOP 99 ! missing file

  npiglo = getdim(cf_in,cn_x)
  npjglo = getdim(cf_in,cn_y)
  npk    = getdim(cf_in,cn_z)
  npt    = getdim(cf_in,cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ! Allocate the memory
  ALLOCATE ( un(npiglo,npjglo), vn(npiglo,npjglo)  )
  ALLOCATE ( umask(npiglo,npjglo), vmask(npiglo,npjglo) )
  ALLOCATE ( u2n(npiglo,npjglo), v2n(npiglo,npjglo)  )
  ALLOCATE ( anousqrt(npiglo,npjglo), anovsqrt(npiglo,npjglo)  )
  ALLOCATE ( tn(npiglo,npjglo), t2n(npiglo,npjglo)  )
  ALLOCATE ( anotsqrt(npiglo,npjglo) )

  CALL CreateOutput

  DO jt =1, npt
     DO jk=1, npk
        PRINT *,'            level ',jk
        umask   (:,:) = 0.0
        vmask   (:,:) = 0.0
        anousqrt(:,:) = 0.0
        anovsqrt(:,:) = 0.0      
        anotsqrt(:,:) = 0.0

        ! JMM COMMENT : this program read data from one file and copy it to
        !               another one. Only the anomalies are computed.
        un(:,:)  = getvar(cf_in, 'ubar',  jk, npiglo, npjglo, ktime=jt)  ! U-point
        vn(:,:)  = getvar(cf_in, 'vbar',  jk, npiglo, npjglo, ktime=jt)  ! V-point
        tn(:,:)  = getvar(cf_in, 'tbar',  jk, npiglo, npjglo, ktime=jt)  ! T-point
        u2n(:,:) = getvar(cf_in, 'u2bar', jk, npiglo, npjglo, ktime=jt)  ! U-point
        v2n(:,:) = getvar(cf_in, 'v2bar', jk, npiglo, npjglo, ktime=jt)  ! V-point
        t2n(:,:) = getvar(cf_in, 't2bar', jk, npiglo, npjglo, ktime=jt)  ! T-point

        ! compute velocity  mask at T-point (??)
        !  JMM COMMENT : this will mask value in the first ocean point from the coast.
        DO jj = 2, npjglo
           DO ji = 2, npiglo
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
        ierr = putvar(ncout, id_varout(1), tn,       jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(2), un,       jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(3), vn,       jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(4), anotsqrt, jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncout, id_varout(5), anousqrt, jk, npiglo, npjglo, ktime=jt) 
        ierr = putvar(ncout, id_varout(6), anovsqrt, jk, npiglo, npjglo, ktime=jt)
     END DO
  ENDDO

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
    ipk(:)                 = npk  
    stypvar(1)%ichunk      = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname       = 'tbar'
    stypvar(1)%clong_name  = 'temporal mean of the temperature on T point'
    stypvar(1)%cshort_name = 'tbar'

    stypvar(2)%ichunk      = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname       = 'ubar'
    stypvar(2)%clong_name  = 'temporal mean of the zonal velocity on T point'
    stypvar(2)%cshort_name = 'ubar'

    stypvar(3)%ichunk      = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(3)%cname       = 'vbar'
    stypvar(3)%clong_name  = 'temporal mean of the meridional velocity on T point'
    stypvar(3)%cshort_name = 'vbar'

    stypvar(4)%ichunk      = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(4)%cname       = 'anotsqrt'
    stypvar(4)%clong_name  = 'temporal mean of the square of the temperature anomaly on T point (*1000)'
    stypvar(4)%cshort_name = 'anotsqrt'

    stypvar(5)%ichunk      = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(5)%cname       = 'anousqrt'
    stypvar(5)%clong_name  = 'temporal mean of the square of the zonal speed anomaly on T point (*1000)'
    stypvar(5)%cshort_name = 'anousqrt'

    stypvar(6)%ichunk      = (/npiglo,MAX(1,npjglo/30),1,1 /)
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
    ncout = create      (cf_out, cf_in,   npiglo, npjglo, npk      , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar, 6,      ipk,    id_varout, ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_in,   npiglo, npjglo, npk       )

    dtim = getvar1d(cf_in,cn_vtimec, npt     )
    ierr = putvar1d(ncout, dtim,     npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfnrj_components

