PROGRAM cdfkempemekeepe
  !!======================================================================
  !!                     ***  PROGRAM  cdfkempemekeepe  ***
  !!=====================================================================
  !!  ** Purpose : Compute the term of energetic transfert from mean kinetic
  !!               energy to mean potential energy (T1) and from eddy 
  !!               potential energy to eddy kinetic energy (T3)
  !!
  !! History : 2.1  : 03/2008  : A. Melet     : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                             :: jj, jk         ! dummy loop index
  INTEGER(KIND=4)                             :: npiglo, npjglo ! size of the domain
  INTEGER(KIND=4)                             :: npk, npt       ! size of the domain
  INTEGER(KIND=4)                             :: narg, iargc    ! browse line
  INTEGER(KIND=4)                             :: ncout, ierr    ! ncid of outputfile, error status
  INTEGER(KIND=4), DIMENSION(2)               :: ipk, id_varout ! levels and varid's of output vars

  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: wbartbar
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: anowt 
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: t1mask
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: w1mask
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: txz
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: wxz 
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: wtxz
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: wbartbarxz
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: anowtxz
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: tim            ! time counter (dummy)

  CHARACTER(LEN=256)                          :: cf_uvwtfil          ! input file
  CHARACTER(LEN=256)                          :: cf_out='transfertst1t3.nc'

  TYPE  (variable), DIMENSION(2)              :: stypvar        ! structure for attibutes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg /= 1 ) THEN
     PRINT *,'usage : cdfkempemekeepe file'
     PRINT *,'     Produce a cdf file transfertst1t3.nc with wT and anowT variables'
     PRINT *,'     file is from cdfmoyuvwt'
     PRINT *,'     the mean must have been computed on a period long enough'
     PRINT *,'     for the statistics to be meaningful'
     PRINT *,'                         '
     STOP
  ENDIF

  CALL getarg(1, cf_uvwtfil)

  IF (chkfile(cf_uvwtfil) ) STOP ! missing file
  npiglo = getdim(cf_uvwtfil, cn_x)
  npjglo = getdim(cf_uvwtfil, cn_y)
  npk    = getdim(cf_uvwtfil, cn_z)
  npt    = getdim(cf_uvwtfil, cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ! define new variables for output ( must update att.txt)
  ipk(:)                    = npk  
  stypvar(1)%cname          = 'wT'
  stypvar(1)%clong_name     = 'temporal mean of w times temporal mean of T on T point (*1000)'
  stypvar(1)%cshort_name    = 'wT'

  stypvar(2)%cname          = 'anowT'
  stypvar(2)%clong_name     = 'temporal mean of anomaly of w times ano of T on T point (*1000)'
  stypvar(2)%cshort_name    = 'anowT'
 
  stypvar%cunits            = '1000 m.K'
  stypvar%rmissing_value    = 0.
  stypvar%valid_min         = -1000.
  stypvar%valid_max         = 1000.
  stypvar%conline_operation = 'N/A'
  stypvar%caxis             = 'TYX'
  
  ! create output fileset
  ncout = create      (cf_out, cf_uvwtfil,   npiglo, npjglo, npk       )
  ierr  = createvar   (ncout,    stypvar, 2,      ipk,    id_varout )
  ierr  = putheadervar(ncout,    cf_uvwtfil,   npiglo, npjglo, npk       )

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
  ALLOCATE ( tim(npt)  )  
  
  tim  = getvar1d(cf_uvwtfil, cn_vtimec, npt     )
  ierr = putvar1d(ncout, tim,       npt, 'T')

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

END PROGRAM cdfkempemekeepe

