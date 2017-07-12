PROGRAM cdfnrj_transfert
  !!======================================================================
  !!                     ***  PROGRAM  cdfnrj_transfert  ***
  !!=====================================================================
  !!  ** Purpose : Compute the term of energetic transfert from mean kinetic
  !!               energy to mean potential energy (T1) and from eddy 
  !!               potential energy to eddy kinetic energy (T3)
  !!
  !! History : 2.1  : 03/2008  : A. Melet     : Original code
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

  INTEGER(KIND=4)                             :: jj, jk             ! dummy loop index
  INTEGER(KIND=4)                             :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                             :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                             :: narg, iargc, ijarg ! browse line
  INTEGER(KIND=4)                             :: ncout, ierr        ! ncid of outputfile, error status
  INTEGER(KIND=4), DIMENSION(2)               :: ipk, id_varout     ! levels and varid's of output vars

  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: wbartbar
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: anowt 
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: t1mask
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: w1mask
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: txz
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: wxz 
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: wtxz
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: wbartbarxz
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: anowtxz

  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: dtim               ! time counter (dummy)

  CHARACTER(LEN=256)                          :: cf_uvwtfil         ! input file
  CHARACTER(LEN=256)                          :: cf_out='trf_t1t3.nc' ! output file
  CHARACTER(LEN=256)                          :: cldum              ! working char variable

  TYPE  (variable), DIMENSION(2)              :: stypvar            ! structure for attibutes

  LOGICAL                                     :: lnc4 = .FALSE.     ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfnrj_transfert -f UVWT-file [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       This program computes the energy transfert term from previously'
     PRINT *,'       computed high order moments (cdfuvwt). High order moments must'
     PRINT *,'       have been evaluated on a long enough period, in order to get '
     PRINT *,'       meaningfull statistics.'
     PRINT *,'       Note : this program was formerly named cdfkempemekeepe. (no idea of the'
     PRINT *,'       pronunciation :) ).'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f UVWT-file : Input file is the output of cdfuvwt, holding the required'
     PRINT *,'             high order moments.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file] : Specify the output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ] : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'             This option is effective only if cdftools are compiled with'
     PRINT *,'             a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless option -o is used.'
     PRINT *,'         variables : WT    : temporal mean of Wbar x Tbar at T point.'
     PRINT *,'                     anoWT : temporal mean of W''xT'' at T points.'
     PRINT *,'             units : 1000 x Celsius x m/s'
     PRINT *,'     SEE ALSO :'
     PRINT *,'      cdfuvwt, cdfnrj_bci, cdfnrj_bti, cdfnrj_components' 
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_uvwtfil ) ; ijarg=ijarg+1
        ! option
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out     ) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 99
     END SELECT
  ENDDO

  IF (chkfile(cf_uvwtfil) ) STOP 99 ! missing file
  npiglo = getdim(cf_uvwtfil, cn_x)
  npjglo = getdim(cf_uvwtfil, cn_y)
  npk    = getdim(cf_uvwtfil, cn_z)
  npt    = getdim(cf_uvwtfil, cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ! Allocate the memory
  ALLOCATE ( wbartbar(  npiglo, npjglo, npk) )  ! 3D can be huge !
  ALLOCATE ( anowt(     npiglo, npjglo, npk) )  ! 3D can be huge
  ALLOCATE ( t1mask(    npiglo,npk) )
  ALLOCATE ( w1mask(    npiglo,npk) )
  ALLOCATE ( txz(       npiglo,npk) )
  ALLOCATE ( wxz(       npiglo,npk) )
  ALLOCATE ( wtxz(      npiglo,npk) )
  ALLOCATE ( anowtxz(   npiglo,npk) )
  ALLOCATE ( wbartbarxz(npiglo,npk) )  
  ALLOCATE ( dtim(npt)  )  

  CALL CreateOutput

  DO jj = 1, npjglo
     wbartbarxz(:,:) = 0.0
     anowtxz(:,:)    = 0.0
     wtxz(:,:)       = 0.0
     wxz(:,:)        = 0.0
     txz(:,:)        = 0.0

     txz( :,:) = getvarxz(cf_uvwtfil, 'tbar',  jj, npiglo, npk, kimin=1, kkmin=1, ktime=1)
     wxz( :,:) = getvarxz(cf_uvwtfil, 'wbar',  jj, npiglo, npk, kimin=1, kkmin=1, ktime=1)
     wtxz(:,:) = getvarxz(cf_uvwtfil, 'wtbar', jj, npiglo, npk, kimin=1, kkmin=1, ktime=1)

     DO jk=1, npk-1
        w1mask(:,jk) = wxz(:,jk) * wxz(:,jk+1)
        t1mask(:,jk) = txz(:,jk)
        WHERE ( w1mask(:,jk) /= 0.) w1mask(:,jk)=1.
        WHERE ( t1mask(:,jk) /= 0.) t1mask(:,jk)=1. 
        wbartbarxz(:,jk) = 1000. * t1mask(:,jk) * txz(:,jk) * 0.5 * w1mask(:,jk) * ( wxz(:,jk) + wxz(:,jk+1) )
        anowtxz(   :,jk) = 1000. * ( wtxz(:,jk) - wbartbarxz(:,jk)*0.001 )      
     END DO
     wbartbar(:,jj,:) = wbartbarxz(:,:)
     anowt(   :,jj,:) = anowtxz(   :,:)
  END DO

  DO jk=1,npk
     ierr = putvar(ncout, id_varout(1), wbartbar(:,:,jk), jk, npiglo, npjglo )
     ierr = putvar(ncout, id_varout(2), anowt(   :,:,jk), jk, npiglo, npjglo )
  END DO

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

    ipk(:)                    = npk  
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname          = 'wT'
    stypvar(1)%clong_name     = 'temporal mean of w times temporal mean of T on T point (*1000)'
    stypvar(1)%cshort_name    = 'wT'

    stypvar(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname          = 'anowT'
    stypvar(2)%clong_name     = 'temporal mean of anomaly of w times ano of T on T point (*1000)'
    stypvar(2)%cshort_name    = 'anowT'

    stypvar%cunits            = '1000 m/s.K'
    stypvar%rmissing_value    = 0.
    stypvar%valid_min         = -1000.
    stypvar%valid_max         = 1000.
    stypvar%conline_operation = 'N/A'
    stypvar%caxis             = 'TYX'

    ! create output fileset
    ncout = create      (cf_out, cf_uvwtfil,   npiglo, npjglo, npk   , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,    stypvar, 2,      ipk,    id_varout, ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,    cf_uvwtfil,   npiglo, npjglo, npk       )

    dtim = getvar1d(cf_uvwtfil, cn_vtimec, npt )
    ierr = putvar1d(ncout, dtim,      npt, 'T' )

  END SUBROUTINE CreateOutput

END PROGRAM cdfnrj_transfert

