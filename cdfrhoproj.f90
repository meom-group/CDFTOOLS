PROGRAM cdfrhoproj
  !!======================================================================
  !!                     ***  PROGRAM  cdfrhoproj  ***
  !!=====================================================================
  !!  ** Purpose : This program is used to project any scalar on the A grid
  !!               onto given isopycnic surfaces.
  !!
  !!  ** Method  : Linear interpolation is used on the vertical to define
  !!               the depth of the given isopycn and linear interpolation
  !!               is also performed on the scalar to determine its value at
  !!               this depth.
  !!
  !! History :      :  1996    : J.M. Molines for SPEM in Dynamo
  !!                :  1999    : J.O. Beismann for OPA
  !!                :  2000    : J.M. Molines for normalization
  !!           2.1  : 11/2005  : J.M. Molines : netcdf 
  !!           3.0  : 05/2011  : J.M. Molines : Doctor norm + Lic.
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

  INTEGER(KIND=4)                               :: ji,jj,jk,jsig,jfich, jvar
  INTEGER(KIND=4)                               :: npiglo, npjglo
  INTEGER(KIND=4)                               :: npk, npsig, npt
  INTEGER(KIND=4)                               :: nvars, nvout=2
  INTEGER(KIND=4)                               :: narg, iargc
  INTEGER(KIND=4)                               :: ijarg, ireq
  INTEGER(KIND=4)                               :: ik0, ijk
  INTEGER(KIND=4)                               :: istartarg = 1
  INTEGER(KIND=4)                               :: nfilin
  INTEGER(KIND=4)                               :: numlev=10
  INTEGER(KIND=4)                               :: ncout, ierr
  INTEGER(KIND=4), DIMENSION(2)                 :: ipk, id_varout ! for output variables
  !
  REAL(KIND=4), DIMENSION(:,:,:),   ALLOCATABLE :: zsig, alpha
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2dint, zint, v2d
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: zi, tim, h1d
  REAL(KIND=4)                                  :: P1, P2
  REAL(KIND=4)                                  :: zalpha
  REAL(KIND=4)                                  :: zspvalo=999999.
  REAL(KIND=4)                                  :: zspvali=0.

  CHARACTER(LEN=256)                            :: cf_rholev='rho_lev'
  CHARACTER(LEN=256)                            :: cf_dta
  CHARACTER(LEN=256)                            :: cf_rhofil
  CHARACTER(LEN=256)                            :: cf_out
  CHARACTER(LEN=256)                            :: cv_in
  CHARACTER(LEN=256)                            :: cv_sig
  CHARACTER(LEN=256)                            :: ctype='T'
  CHARACTER(LEN=256)                            :: cldum
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names     ! temporary arry for variable name in file
 
  TYPE(variable), DIMENSION(2)                  :: stypvar      ! structure for attributes
  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypzvar      ! structure for attributes
  !
  LOGICAL                                       :: lsingle =.FALSE.
  LOGICAL                                       :: lchk    =.FALSE.
  LOGICAL                                       :: lisodep =.FALSE.
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  cv_sig = cn_vosigma0

  narg=iargc()
  IF ( narg < 3 ) THEN
     PRINT *,' usage : cdfrhoproj IN-var RHO-file List_of_IN-files [VAR-type] ... '
     PRINT *,'                ... [-s0 sig0 ] [-sig sigma_name] [-isodep ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Project IN-var on isopycnal surfaces defined either by sig0 given' 
     PRINT *,'       as argument or on all sigma surfaces defined in ',TRIM(cf_rholev),' ascii file.'
     PRINT *,'       IN-var will be interpolated on the T point of the C-grid, previous'
     PRINT *,'       to projection on isopycnal.'
     PRINT *,'       This cdftool is one of the few using 3D arrays. Further development is '
     PRINT *,'       required to work with vertical slabs instead.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       IN-var   : name of the input variable to be projected' 
     PRINT *,'       RHO-file : netcdf file with potential density field. If not a sigma0'
     PRINT *,'                  file, use -sig option to indicate the name of the density'
     PRINT *,'                  variable.'
     PRINT *,'       List_of_IN-file  : netcdf files with IN-var '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-s0 sigma ] : define a single sigma surface on the command line' 
     PRINT *,'                    instead of reading rho_lev ascii file.'
     PRINT *,'       [VAR-type] : position of IN-var on the C-grid ( either T U V F W )'
     PRINT *,'                    default is ''T''.'
     PRINT *,'       [-sig sigma_name] : name of the density variable in RHO_file.'
     PRINT *,'                    default is ', TRIM(cv_sig)
     PRINT *,'       [-isodep ] : only compute the isopycnal depth. In this case you must'
     PRINT *,'                    still specify a IN-var variable (in fact a dummy name).'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       no metrics, information is taken from depth variable in input files.'
     PRINT *,'       ', TRIM(cf_rholev),' if not using -s0 option.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       There are as many output files as input files.'
     PRINT *,'       netcdf file : IN-file.interp'
     PRINT *,'         variables : VAR-in (unit is the same as input var)'
     PRINT *,'                     ', TRIM(cn_vodepiso),' (m) : depth of isopycnal.'
     PRINT *,'      '
     PRINT *,'       If option -isodep is used, only isopycnal depth is output :'
     PRINT *,'       netcdf file : isopycdep.nc'
     PRINT *,'         variables : ',TRIM(cn_vodepiso),' (m) '
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       replace cdfisopycdep when using -isodep option.'
     PRINT *,'       '
     STOP
  ENDIF

  ijarg = 1 ; ireq=0 ; nfilin=0

  DO WHILE ( ijarg <= narg )
    CALL getarg( ijarg, cldum) ; ijarg=ijarg+1
    SELECT CASE ( cldum )
    CASE ('-s0') 
       npsig = 1 ; lsingle=.TRUE. ; ALLOCATE (zi(npsig) )
       CALL getarg( ijarg, cldum) ; ijarg=ijarg+1 ; READ(cldum,*) zi(1)
    CASE ( 'T','t','U','u','V','v','W','w','F','f' )
       ctype=cldum
    CASE ('-sig') 
       CALL getarg( ijarg, cv_sig) ; ijarg=ijarg+1 
    CASE ('-isodep')  ; lisodep = .TRUE. ; nvout=1 ; cf_out='isopycdep.nc'
    CASE DEFAULT 
       ireq=ireq+1
       SELECT CASE (ireq )
       CASE ( 1 ) ; cv_in     = cldum
       CASE ( 2 ) ; cf_rhofil = cldum
       CASE DEFAULT 
         ! count the input files
         nfilin=nfilin+1
         IF ( nfilin == 1 ) istartarg=ijarg-1
       END SELECT
    END SELECT
  END DO

  lchk = chkfile(cf_rhofil)
  IF ( .NOT. lsingle ) lchk = lchk .OR. chkfile(cf_rholev)
  IF ( lchk ) STOP ! missing file

  IF ( .NOT.  lsingle ) THEN
     OPEN(numlev,FILE=cf_rholev)
     READ(numlev,*) npsig
     ALLOCATE ( zi(npsig) )
     DO jsig=1,npsig
        READ(numlev,*) zi(jsig)
        PRINT *,zi(jsig)
     END DO
     CLOSE(numlev)
  ENDIF

  ! Read Rho file
  npiglo = getdim(cf_rhofil,cn_x)
  npjglo = getdim(cf_rhofil,cn_y)
  npk    = getdim(cf_rhofil,cn_z)
  npt    = getdim(cf_rhofil,cn_t)

  CALL getarg(istartarg, cf_dta)
  nvars=getnvar(cf_dta)
  ALLOCATE(cv_names(nvars), stypzvar(nvars))

  cv_names(:)=getvarname(cf_dta, nvars, stypzvar)

  ALLOCATE( zsig(npiglo,npjglo,npk), alpha(npiglo, npjglo, npsig)            )
  ALLOCATE( v2dint(npiglo, npjglo), v2d(npiglo,npjglo), zint(npiglo,npjglo)  )
  ALLOCATE( tim(npt), h1d(npk)                                               )

  tim(:)=getvar1d(cf_rhofil, cn_vtimec,  npt)
  h1d(:)=getvar1d(cf_rhofil, cn_vdeptht, npk)
   
  DO jk=1,npk
     zsig(:,:,jk) = getvar(cf_rhofil, cv_sig, jk, npiglo, npjglo)
  END DO

  !! ** Compute interpolation coefficients as well as the level used
  !!    to interpolate between
  DO ji=1,npiglo
     DO jj = 1, npjglo
        ijk = 1
        DO jsig=1,npsig
        !  Assume that rho (z) is increasing downward (no inversion)
        !     Caution with sigma0 at great depth !
           DO WHILE (zi(jsig) >=  zsig(ji,jj,ijk) .AND. ijk <= npk &
                  &   .AND. zsig(ji,jj,ijk) /=  zspvali )
              ijk=ijk+1
           END DO
           ijk=ijk-1
           ik0=ijk
           IF (ijk == 0) THEN
              ijk=1
              alpha(ji,jj,jsig) = 0.
           ELSE IF (zsig(ji,jj,ijk+1) == zspvali ) THEN
              ik0=0
              alpha(ji,jj,jsig) = 0.
           ELSE 
           ! ... alpha is always in [0,1]. Adding ik0 ( >=1 ) for saving space for ik0
              alpha(ji,jj,jsig)= &
                  &  (zi(jsig)-zsig(ji,jj,ijk))/(zsig(ji,jj,ijk+1)-zsig(ji,jj,ijk)) + ik0
           ENDIF
        END DO
     END DO
  END DO

  IF ( lisodep ) THEN
     ipk(1)                       = npsig
     stypvar(1)%cname             = cn_vodepiso
     stypvar(1)%cunits            = 'm'
     stypvar(1)%rmissing_value    = 999999.
     stypvar(1)%valid_min         = 0.
     stypvar(1)%valid_max         = 7000.
     stypvar(1)%clong_name        = 'Depth_of_Isopycnals'
     stypvar(1)%cshort_name       = cn_vodepiso
     stypvar(1)%conline_operation = 'N/A'
     stypvar(1)%caxis             = 'TRYX'

     ncout = create      (cf_out, cf_rhofil, npiglo, npjglo, npsig          )
     ierr  = createvar   (ncout,  stypvar,   nvout,  ipk,    id_varout      )
     ierr  = putheadervar(ncout , cf_rhofil, npiglo, npjglo, npsig, pdep=zi )

     DO jsig=1,npsig
        DO ji=1,npiglo
           DO jj=1,npjglo
             ! ik0 is retrieved from alpha, taking the integer part.
             ! The remnant is alpha.
             ik0 = INT(alpha(ji,jj,jsig))
             zalpha =  alpha(ji,jj,jsig) - ik0
             IF (ik0 /= 0) THEN
              P1 = zsig(ji,jj,ik0  )
              P2 = zsig(ji,jj,ik0+1)
                IF (P1 /= zspvali .AND. P2 /= zspvali) THEN
                    zint (ji,jj) = zalpha *h1d(ik0+1) &
                     &         +(1-zalpha)*h1d(ik0  )
                ELSE
                   zint  (ji,jj)=zspvalo
               ENDIF
             ELSE
               zint  (ji,jj)=zspvalo
             ENDIF
           END DO
        END DO
        ierr = putvar(ncout, id_varout(1), zint , jsig, npiglo, npjglo)
     END DO
     ierr = closeout(ncout    )
  ENDIF

  !! ** Loop on the scalar files to project on choosen isopycnics surfaces
  DO jfich= 1, nfilin
     ijarg = istartarg + jfich - 1

     CALL getarg(ijarg, cf_dta)
     PRINT *,'working with ', TRIM(cf_dta)

     npt    = getdim(cf_dta, cn_t)
     IF (npt /= 1 ) THEN
        PRINT *,' This program has to be modified for multiple'
        PRINT *,' time frames.'
        STOP ' Error : npt # 1'
     ENDIF
     tim(:)=getvar1d(cf_dta, cn_vtimec, 1)
 
     DO jk=1,npk
        v2d(:,:) = getvar(cf_dta,cv_in,jk,npiglo,npjglo)
        SELECT CASE ( ctype )
        CASE ('T', 't' )
           zsig(:,:,jk) = v2d(:,:)
        CASE ('U','u' )
           DO ji=2,npiglo
              DO jj=1, npjglo
                 zsig(ji,jj,jk)=0.5*( v2d(ji,jj) + v2d(ji-1,jj) )  ! put variable on T point
              END DO
           END DO
        CASE ('V','v' )
           DO jj=2,npjglo
              DO ji=1, npiglo
                 zsig(ji,jj,jk)=0.5*( v2d(ji,jj) + v2d(ji,jj-1) )  ! put variable on T point
              END DO
           END DO
         CASE('W','w' )
           zint(:,:) = getvar(cf_dta, cv_in, jk+1, npiglo, npjglo)
           DO jj=1,npjglo
              DO ji=1, npiglo
                 zsig(ji,jj,jk)=0.5*( v2d(ji,jj) + zint(ji,jj) )  ! put variable on T point
              END DO
           END DO
         CASE('F','f' )
           DO jj=2,npjglo
              DO ji=2, npiglo
                 zsig(ji,jj,jk)=0.25*( v2d(ji,jj) + v2d(ji,jj-1) + v2d(ji-1,jj) + v2d(ji-1,jj-1 )) ! put variable on T point
              END DO
           END DO
         END SELECT
     END DO

     ! ... open output file and write header
     ipk(:)=npsig
     DO jvar=1,nvars
       IF ( cv_in == stypzvar(jvar)%cname ) THEN 
          stypvar(2)=stypzvar(jvar)
          EXIT
       ENDIF
     END DO
     stypvar(2)%clong_name        = TRIM(stypvar(2)%clong_name)//' on iso sigma'
     stypvar(2)%caxis             = 'TRYX'

     stypvar(1)%cname             = cn_vodepiso
     stypvar(1)%cunits            = 'm'
     stypvar(1)%rmissing_value    = 999999.
     stypvar(1)%valid_min         = 0.
     stypvar(1)%valid_max         = 7000.
     stypvar(1)%clong_name        = 'Depth_of_Isopycnals'
     stypvar(1)%cshort_name       = cn_vodepiso
     stypvar(1)%conline_operation = 'N/A'
     stypvar(1)%caxis             = 'TRYX'

     cf_out=TRIM(cf_dta)//'.interp'
 
     ncout = create      (cf_out, cf_rhofil, npiglo, npjglo, npsig          )
     ierr  = createvar   (ncout,  stypvar,   nvout,  ipk,    id_varout      )
     ierr  = putheadervar(ncout , cf_rhofil, npiglo, npjglo, npsig, pdep=zi )
 
     DO jsig=1,npsig
        DO ji=1,npiglo
           DO jj=1,npjglo
             ! ik0 is retrieved from alpha, taking the integer part.
             ! The remnant is alpha. 
             ik0    = INT(alpha(ji,jj,jsig))
             zalpha =     alpha(ji,jj,jsig) - ik0
             IF (ik0 /= 0) THEN
              P1 = zsig(ji,jj,ik0  )
              P2 = zsig(ji,jj,ik0+1)
                IF (P1 /= zspvali .AND. P2 /= zspvali) THEN
                   v2dint(ji,jj) = zalpha *P2  &
                     &         +(1-zalpha)*P1
                    zint (ji,jj) = zalpha *h1d(ik0+1) &
                     &         +(1-zalpha)*h1d(ik0  )
                ELSE 
                   v2dint(ji,jj)=zspvalo
                   zint  (ji,jj)=zspvalo
               ENDIF
             ELSE 
               v2dint(ji,jj)=zspvalo
               zint  (ji,jj)=zspvalo
             ENDIF
           END DO
        END DO
        ierr = putvar(ncout, id_varout(1), zint  , jsig, npiglo, npjglo)
        ierr = putvar(ncout, id_varout(2), v2dint, jsig, npiglo, npjglo)
     END DO
     ierr = putvar1d(ncout, tim, 1, 'T')
     ierr = closeout(ncout             )
  END DO  ! loop on scalar files
        PRINT *,'Projection on isopycns completed successfully'
END  PROGRAM cdfrhoproj
