PROGRAM cdfmkmask
  !!======================================================================
  !!                     ***  PROGRAM  cdfmkmask  ***
  !!=====================================================================
  !!  ** Purpose : Build mask file from a salinity output
  !!
  !!  ** Method  : Read vosaline and set tmask to 1 where sal is not 0
  !!               then umask, vmask and fmask are deduced from tmask
  !!               REM: the result may be locally different for fmask than
  !!                   fmask produced online as there are computed on line
  !!               merged with cdfmkmask-zone by adding a zoom option. When
  !!               used with -zoom option, the mask is 0 outside the zoom
  !!               area.
  !!
  !! History : 2.1  : 11/2005  : J.M. Molines : Original code
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

  INTEGER(KIND=4)                           :: ji, jj, jk               ! dummy loop index
  INTEGER(KIND=4)                           :: ierr                     ! working integer
  INTEGER(KIND=4)                           :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                           :: npiglo, npjglo, npk      ! size of the domain
  INTEGER(KIND=4)                           :: ncout                    ! ncid of output file
  INTEGER(KIND=4), DIMENSION(4)             :: ipk, id_varout           ! outptut variables : number of levels,

  
  REAL(KIND=4)                              :: rlonmin, rlonmax         ! limit in longitude
  REAL(KIND=4)                              :: rlatmin, rlatmax         ! limit in latitude
  REAL(KIND=4), DIMENSION(1)                :: tim                      ! time counter
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask, zmask             ! 2D mask at current level
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlon, rlat               ! latitude and longitude

  CHARACTER(LEN=256)                        :: cf_tfil                  ! file name
  CHARACTER(LEN=256)                        :: cf_out = 'mask_sal.nc'   ! output file
  CHARACTER(LEN=256)                        :: cldum                    ! dummy string

  TYPE (variable), DIMENSION(4)             :: stypvar                  ! output attribute
 
  LOGICAL                                   :: lzoom = .false.          ! zoom flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmkmask T-file [-zoom lonmin lonmax latmin latmax] ...'
     PRINT *,'                   ... [-o OUT-file ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Build a mask file from vosaline array read from the input file.' 
     PRINT *,'       It assumes that land salinity values are set to 0.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file : netcdf file with salinity.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-zoom lonmin lonmax latmin latmax] : geographical windows used to'
     PRINT *,'                        limit the area where the mask is builded. Outside'
     PRINT *,'                        this area, the mask is set to 0.'
     PRINT *,'       [-o OUT-file ] : output file name to be used in place of standard'
     PRINT *,'                        name [ ',TRIM(cf_out),' ]'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out), ' or OUT-file.'
     PRINT *,'         variables : tmask, umask, vmask, fmask'
     PRINT *,'                fmask can differ from standard fmask because it does not'
     PRINT *,'                reflect the slip/noslip lateral condition.'
     STOP
  ENDIF

  ijarg = 1
  CALL getarg (ijarg, cf_tfil) ; ijarg = ijarg + 1

  DO WHILE ( ijarg <= narg ) 
    CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
    SELECT CASE ( cldum )
    !
    CASE ( '-zoom' )  ! read a zoom area
       lzoom = .true.
       CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rlonmin
       CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rlonmax
       CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rlatmin
       CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rlatmax
    !
    CASE ( '-o'    )  ! change output file name
       CALL getarg (ijarg, cf_out) ; ijarg = ijarg + 1
    !
    CASE DEFAULT
        PRINT *, 'ERROR : unknown option :', TRIM(cldum)
        STOP
    END SELECT
  ENDDO

  IF ( chkfile(cf_tfil) ) STOP ! missing file

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)

  PRINT *,' npiglo = ', npiglo
  PRINT *,' npjglo = ', npjglo
  PRINT *,' npk    = ', npk


  ipk(1:4)                       = npk
  stypvar(1)%cname               = 'tmask'
  stypvar(2)%cname               = 'umask'
  stypvar(3)%cname               = 'vmask'
  stypvar(4)%cname               = 'fmask'

  stypvar(1:4)%cunits            = '1/0'
  stypvar(1:4)%rmissing_value    = 9999.
  stypvar(1:4)%valid_min         = 0.
  stypvar(1:4)%valid_max         = 1.

  stypvar(1)%clong_name          = 'tmask'
  stypvar(2)%clong_name          = 'umask'
  stypvar(3)%clong_name          = 'vmask'
  stypvar(4)%clong_name          = 'fmask'

  stypvar(1)%cshort_name         = 'tmask'
  stypvar(2)%cshort_name         = 'umask'
  stypvar(3)%cshort_name         = 'vmask'
  stypvar(4)%cshort_name         = 'fmask'

  stypvar(1:4)%conline_operation = 'N/A'
  stypvar(1:4)%caxis             = 'TZYX'
  stypvar(1:4)%cprecision        = 'i2'

  ncout = create      (cf_out, cf_tfil,  npiglo, npjglo, npk)
  ierr  = createvar   (ncout,    stypvar, 4,      ipk,    id_varout )
  ierr  = putheadervar(ncout,    cf_tfil,  npiglo, npjglo, npk)

  ALLOCATE (tmask(npiglo,npjglo), zmask(npiglo,npjglo))

  IF ( lzoom ) THEN
    ALLOCATE (rlon(npiglo,npjglo), rlat(npiglo,npjglo))
    rlon(:,:) = getvar(cf_tfil, cn_vlon2d, 1, npiglo, npjglo)
    rlat(:,:) = getvar(cf_tfil, cn_vlat2d, 1, npiglo, npjglo)
  ENDIF

  DO jk=1, npk
     ! tmask
     tmask(:,:) = getvar(cf_tfil, 'vosaline',  jk, npiglo, npjglo)
     WHERE (tmask > 0 ) tmask = 1
     WHERE (tmask <=0 ) tmask = 0

     IF ( lzoom ) THEN
        IF (rlonmax > rlonmin) THEN
           WHERE (rlon > rlonmax ) tmask = 0
           WHERE (rlon < rlonmin ) tmask = 0
        ELSE
           WHERE (rlon < rlonmin .AND. rlon > rlonmax ) tmask = 0
        END IF

        WHERE (rlat > rlatmax ) tmask = 0
        WHERE (rlat < rlatmin ) tmask = 0
     ENDIF
     ierr       = putvar(ncout, id_varout(1), tmask, jk ,npiglo, npjglo)
     ! umask
     zmask = 0.
     DO ji=1,npiglo-1
       DO jj=1,npjglo
        zmask(ji,jj) = tmask(ji,jj)*tmask(ji+1,jj)
       END DO
     END DO
     ierr       = putvar(ncout, id_varout(2), zmask, jk ,npiglo, npjglo)
    ! vmask
     zmask=0.
     DO ji=1,npiglo
       DO jj=1,npjglo-1
        zmask(ji,jj) = tmask(ji,jj)*tmask(ji,jj+1)
       END DO
     END DO
     ierr       = putvar(ncout, id_varout(3), zmask, jk, npiglo, npjglo)
     !fmask
     zmask=0.
     DO ji=1,npiglo-1
       DO jj=1,npjglo-1
        zmask(ji,jj) = tmask(ji,jj)*tmask(ji,jj+1)*tmask(ji+1,jj)*tmask(ji+1,jj+1)
       END DO
     END DO
     ierr       = putvar(ncout, id_varout(4), zmask, jk, npiglo, npjglo)
  END DO  ! loop to next level

  tim(:) = 0.
  ierr   = putvar1d(ncout, tim, 1,'T')
  ierr   = closeout(ncout              )

END PROGRAM cdfmkmask
