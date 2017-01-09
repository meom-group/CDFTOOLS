PROGRAM cdfcensus
  !!======================================================================
  !!                     ***  PROGRAM  cdfcensus  ***
  !!=====================================================================
  !!  ** Purpose : Build an array giving the volume of water in a TS cell.
  !!
  !!  ** Method  : T-file and S-file are scanned for a given region 
  !!              (eventually limited in depth) and the volume of water in 
  !!              a (T,S) cell such that T < Tmodele < T+dt and 
  !!              S < Smodele < S+ds.
  !!              If Smodel or T model are out of the bound they are 
  !!              cumulated in the nearest (T,S) cell.
  !!                  The output is done on a bimg file where S is given as
  !!              the x-direction and T the y-direction, the field value 
  !!              being the volume of water. Due to a very large range in 
  !!              the water volume over the TS field the field is indeed 
  !!              the LOG (1 + VOLUME), and even, the scale can be made
  !!              more non-linear by repeating the LOG operation, ie, for 
  !!              example, field=LOG(1 + LOG (1 + VOLUME)). The parameter 
  !!              nlog, passed as command argument can be used to fix the
  !!              number of LOG. If nlog = 0, the true volume is saved.
  !!                 Depending on the user purpose, limiting values tmin,
  !!              tmax, and smin,smax as well as the increments dt, ds can 
  !!              be adjusted.
  !!              output is STILL a dimg file
  !!
  !! History : --   : 02/1997  : J.M. Molines as bimgtools in DYNAMO
  !!           --   : 09/1999  : A. de Miranda for OPA
  !!                : 01/2002  : J.M. Molines : DOctor norm
  !!                : 01/2006  : C. Langlais  : CDF I and partial cell
  !!           2.0  : 03/2006  : J.M. Molines : integration in CDFTOOLS
  !!           2.1  : 12/2006  : J.M. Molines : add sigma-2 and sigma-4 O
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jk, jt, jlog
  INTEGER(KIND=4)                           :: npiglo, npjglo            ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt                  ! size of the domain
  INTEGER(KIND=4)                           :: nlog
  INTEGER(KIND=4)                           :: narg, iargc, ijarg
  INTEGER(KIND=4)                           :: it, is
  INTEGER(KIND=4)                           :: ii1, ii2
  INTEGER(KIND=4)                           :: ij1, ij2
  INTEGER(KIND=4)                           :: ik1, ik2
  INTEGER(KIND=4)                           :: nt, ns
  INTEGER(KIND=4)                           :: ncout, ierr
  INTEGER(KIND=4), DIMENSION(2)             :: ijloc
  INTEGER(KIND=4), DIMENSION(4)             :: ipk, id_varout

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zt, zs, rsigma0, rsigma2, rsigma4
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1t, e2t
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e3t
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsx, zty
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: rdumdep, tim
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: e31d
  REAL(KIND=4)                              :: ztmin, ztmax, zdt, ztm
  REAL(KIND=4)                              :: zsmin, zsmax, zds, zsm
  REAL(KIND=4)                              :: ztpoint, zspoint, rcmax

  REAL(KIND=8)                              :: dvoltotal, dvolpoint
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcensus, ddump

  CHARACTER(LEN=256)                        :: cf_tfil
  CHARACTER(LEN=256)                        :: cf_bimg='censusopa.bimg'
  CHARACTER(LEN=256)                        :: cf_out='census.nc'
  CHARACTER(LEN=256)                        :: cglobal
  CHARACTER(LEN=256)                        :: cline1, cline2, cline3, cline4
  CHARACTER(LEN=256)                        :: cldum

  TYPE(variable), DIMENSION(4)              :: stypvar

  LOGICAL                                   :: lcdf=.TRUE. , lbimg=.FALSE.
  LOGICAL                                   :: lchk
  LOGICAL                                   :: lfull = .FALSE.   ! flag for full step

  ! Initialisations
  DATA ztmin, ztmax, zdt /-2.0, 38.0, 0.05/
  DATA zsmin, zsmax, zds /25.0, 40.0, 0.02/
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()


  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfcensus T-file nlog [-zoom imin imax jmin jmax] ...'
     PRINT *,'                ... [-klim kmin kmax]  [-full] [-bimg] ... '
     PRINT *,'                ... [-srange smin smax ds ] ...'
     PRINT *,'                ... [-trange tmin tmax dt ] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute the volumetric water mass census: the ocean is divided in'
     PRINT *,'        T,S bins; the program gives the volume of water for each bin.'
     PRINT *,'        A sub-area can be specified, both horizontaly and vertically.'
     PRINT *,'        Temperature and salinity ranges can be also adapted, as well as the'
     PRINT *,'        width of the bins. Default values are provided. In order to attenuate'
     PRINT *,'        the huge maximum values, a log10 operator can be applied many times,'
     PRINT *,'        the number of filter passes being set on the command line.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file  : netcdf file name for temperature and salinity' 
     PRINT *,'       nlog    : number of log10 filter to perform. Can be 0.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-zoom imin imax jmin jmax] : define a model sub-area, in model '
     PRINT *,'                                     coordinates' 
     PRINT *,'       [-klim ik1 ik2            ] : set limits on the vertical.'
     PRINT *,'       [-srange smin smax ds     ] : define the size of the salinity bin'
     PRINT '(a,2f5.1,x,f6.3)','                         defaut is :', zsmin, zsmax, zds
     PRINT *,'       [-trange tmin tmax dt     ] : define the size of the temperatude bin'
     PRINT '(a,2f5.1,x,f6.3)','                         defaut is :', ztmin, ztmax, zdt
     PRINT *,'       [-full                    ] : use for full step computation'
     PRINT *,'       [-bimg                    ] : output on bimg files (to be deprecated).'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fhgr),'  and ',TRIM(cn_fzgr) 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       - netcdf file : ', TRIM(cf_out) 
     PRINT *,'           variables : volcensus  (10^15 m3 )'
     PRINT *,'                       sigma0  (kg/m3 -1000 )'
     PRINT *,'                       sigma2  (kg/m3 -1000 )'
     PRINT *,'                       sigma3  (kg/m3 -1000 )'
     PRINT *,'       - bimg file   : According to options.'
     STOP
  ENDIF

  ijarg = 1
  CALL getarg(ijarg, cf_tfil) ; ijarg = ijarg + 1
  CALL getarg(ijarg, cldum  ) ; ijarg = ijarg + 1  ; READ(cldum,*) nlog
  cglobal = 'Census computed from '//TRIM(cf_tfil)

  lchk =           chkfile ( cn_fzgr )
  lchk = lchk .OR. chkfile ( cn_fhgr )
  lchk = lchk .OR. chkfile ( cf_tfil  )
  IF ( lchk ) STOP  ! some compulsory files are missing

  PRINT *,' TS_FILE = ',TRIM(cf_tfil)
  PRINT *,' NLOG    = ', nlog

  ! set domain size from TS file
  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  ! Allocate memory
  ALLOCATE (zt(npiglo,npjglo),zs(npiglo,npjglo))
  ALLOCATE (e1t(npiglo,npjglo),e2t(npiglo,npjglo),e31d(npk),e3t(npiglo,npjglo))

  ! Read metrics
  e1t(:,:) = getvar  (cn_fhgr, cn_ve1t, 1, npiglo, npjglo)
  e2t(:,:) = getvar  (cn_fhgr, cn_ve2t, 1, npiglo, npjglo)
  e31d(:)  = getvare3(cn_fzgr, cn_ve3t, npk              )  ! used in full step case

  ! default is full domain, full depth
  ii1 = 1 ; ii2 = npiglo
  ij1 = 1 ; ij2 = npjglo
  ik1 = 1 ; ik2 = npk

  ! Read additional optional argument (zoom)
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg,cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum)
     CASE ( '-zoom' )
        CALL getarg(ijarg,cldum) ; READ(cldum,*) ii1  ; ijarg = ijarg+1
        CALL getarg(ijarg,cldum) ; READ(cldum,*) ii2  ; ijarg = ijarg+1
        CALL getarg(ijarg,cldum) ; READ(cldum,*) ij1  ; ijarg = ijarg+1
        CALL getarg(ijarg,cldum) ; READ(cldum,*) ij2  ; ijarg = ijarg+1
     CASE ( '-klim' )
        CALL getarg(ijarg,cldum) ; READ(cldum,*) ik1  ; ijarg = ijarg+1
        CALL getarg(ijarg,cldum) ; READ(cldum,*) ik2  ; ijarg = ijarg+1
     CASE ( '-bimg' )
        lbimg = .TRUE.
        lcdf  = .FALSE.
     CASE ( '-srange' )
        CALL getarg(ijarg,cldum) ; READ(cldum,*) zsmin ; ijarg = ijarg+1
        CALL getarg(ijarg,cldum) ; READ(cldum,*) zsmax ; ijarg = ijarg+1
        CALL getarg(ijarg,cldum) ; READ(cldum,*) zds   ; ijarg = ijarg+1
     CASE ( '-trange' )
        CALL getarg(ijarg,cldum) ; READ(cldum,*) ztmin ; ijarg = ijarg+1
        CALL getarg(ijarg,cldum) ; READ(cldum,*) ztmax ; ijarg = ijarg+1
        CALL getarg(ijarg,cldum) ; READ(cldum,*) zdt   ; ijarg = ijarg+1
     CASE ( '-full' )
        lfull = .TRUE.
     CASE DEFAULT
        PRINT *,' Unknown option :',TRIM(cldum)
        STOP
     END SELECT
  END DO

  ! Extra checking for over bound
  ii1 = MAX(ii1,1) ; ii2 = MIN(ii2,npiglo)
  ij1 = MAX(ij1,1) ; ij2 = MIN(ij2,npjglo)
  ik1 = MAX(ik1,1) ; ik2 = MIN(ik2,npk   )

  PRINT '(a,6i5)','indices:',ii1, ii2, ij1, ij2, ik1, ik2

  ! Compute the census on the requested domain
  PRINT *,' Water mass census on the file '
  PRINT *, TRIM(cf_tfil)
  PRINT *, ' running .........'
  nt = NINT( (ztmax - ztmin )/zdt + 1 )
  ns = NINT( (zsmax - zsmin )/zds + 1 )

  ! Allocate arrays
  ALLOCATE ( dcensus (ns,nt), ddump(ns,nt) )
  ALLOCATE ( rsigma0(ns,nt), rsigma2(ns,nt), rsigma4(ns,nt) )
  ALLOCATE ( zsx (ns,nt), zty(ns,nt), rdumdep(1), tim(npt))

  ! fill up rsigma0 array with theoretical density
  DO ji=1,ns
     DO jj=1,nt
        zsx(ji,jj) = zsmin + (ji-1)*zds
        zty(ji,jj) = ztmin + (jj-1)*zdt
     END DO
  END DO
  
  rsigma0 = sigma0(zty, zsx,        ns, nt)
  rsigma2 = sigmai(zty, zsx, 2000., ns, nt)
  rsigma4 = sigmai(zty, zsx, 4000., ns, nt)
  rdumdep(1) = 0.

  IF ( lcdf ) THEN
     ! create output fileset
     ipk(:)= 1                  
     stypvar%rmissing_value    = -100.
     stypvar%valid_min         = 0.
     stypvar%valid_max         = 1.e20
     stypvar%conline_operation = 'N/A'
     stypvar%caxis             = 'TYX'

     stypvar(1)%cname          = 'volcensus'
     stypvar(2)%cname          = 'sigma0'
     stypvar(3)%cname          = 'sigma2'
     stypvar(4)%cname          = 'sigma4'

     stypvar(1)%cunits         = 'm3'
     stypvar(2:4)%cunits       = 'kg/m3'

     stypvar(1)%clong_name     = 'Volume_Census_TS'
     stypvar(2)%clong_name     = 'Sigma0_TS'
     stypvar(3)%clong_name     = 'Sigma2_TS'
     stypvar(4)%clong_name     = 'Sigma4_TS'

     stypvar(1)%cshort_name    = 'volcensus'
     stypvar(2)%cshort_name    = 'sigma0'
     stypvar(3)%cshort_name    = 'sigma2'
     stypvar(4)%cshort_name    = 'sigma4'

     ncout = create      (cf_out, cf_tfil,  ns, nt,  1                                           )
     ierr  = createvar   (ncout,  stypvar,  4,  ipk, id_varout, cdglobal=cglobal                 )
     ierr  = putheadervar(ncout,  cf_tfil,  ns, nt,  1,   pnavlon=zsx, pnavlat=zty, pdep=rdumdep )
  ENDIF

  DO jt = 1, npt
     ! reset cumulating variables to 0
     dcensus(:,:) = 0.d0
     dvoltotal    = 0.d0
     ! Enter main loop
     DO jk=ik1,ik2
        zt(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime = jt)
        zs(:,:) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime = jt)

        IF ( lfull ) THEN
           e3t(:,:) = e31d(jk)
        ELSE
           e3t(:,:) = getvar(cn_fzgr, 'e3t_ps', jk, npiglo, npjglo, ldiom=.TRUE.)
        ENDIF

        DO ji=ii1,ii2
           DO jj=ij1,ij2
              ztpoint   = zt(ji,jj)
              zspoint   = zs(ji,jj)
              dvolpoint = e1t(ji,jj)*e2t(ji,jj)*e3t(ji,jj)*1.d0

              ! salinity = 0 on masked points ( OPA !!! )
              IF (zspoint /= 0) THEN
                 it=NINT( (ztpoint-ztmin)/zdt) + 1
                 is=NINT( (zspoint-zsmin)/zds) + 1
                 ! check for out of bound values
                 it = MIN ( MAX(it,1), nt )
                 is = MIN ( MAX(is,1), ns )

                 dcensus(is,it) = dcensus(is,it) + dvolpoint*1.d-15
                 dvoltotal      = dvoltotal      + dvolpoint*1.d-15
              END IF
           END DO
        END DO

     END DO  ! Main loop

     ! Computes some statistics
     rcmax = MAXVAL ( dcensus )
     ijloc = MAXLOC ( dcensus )
     zsm    = zsmin + (ijloc(1) -1 ) * zds
     ztm    = ztmin + (ijloc(2) -1 ) * zdt

     PRINT *,'  Total Volume of the domain in  10^15 m3:', REAL(dvoltotal)
     PRINT *,'  Volume of the most represented water mass :', rcmax 
     PRINT '(a,f6.2,a)' ,'       this is about ', rcmax/dvoltotal *100,' % of the total'
     PRINT *,'     Salinity   = ', zsm
     PRINT *,'     Temperature= ', ztm

     ! use a distorsion function ( n x log ) to reduce extrema in the output file.
     ddump(:,:) = dcensus(:,:)
     DO jlog = 1, nlog
        ddump(:,:) = LOG10 (1.d0 + ddump(:,:) )
     ENDDO

     IF ( lcdf ) THEN
        ! Output on census.nc file
        ierr = putvar(ncout, id_varout(1), REAL(ddump), 1, ns, nt, ktime=jt)
        ierr = putvar(ncout, id_varout(2), rsigma0,     1, ns, nt, ktime=jt)
        ierr = putvar(ncout, id_varout(3), rsigma2,     1, ns, nt, ktime=jt)
        ierr = putvar(ncout, id_varout(4), rsigma4,     1, ns, nt, ktime=jt)
     ENDIF
  ENDDO  ! time loop

  IF ( lcdf ) THEN
     tim  = getvar1d(cf_tfil, cn_vtimec, npt     )
     ierr = putvar1d(ncout,   tim,       npt, 'T')
     ierr = closeout(ncout)
  ENDIF

  IF (lbimg ) THEN
     ! Output on bimg file
     OPEN (10,file=cf_bimg,form='UNFORMATTED')

     WRITE(cline1,942)' Water Masses Census [10-15 m3] on',ii1,ii2,ij1,ij2
942  FORMAT(a,4i5)
     cline2 = '   computed from the following T-S files:'
     cline3 = cf_tfil
     cline4 = ''
     !
     WRITE(10) cline1
     WRITE(10) cline2
     WRITE(10) cline3
     WRITE(10) cline4
     WRITE(10) ns, nt, 1, 1, 4, nlog
     WRITE(10) zsmin, ztmin, zds, zdt, 0.
     WRITE(10) 0.
     WRITE(10) 0.
     WRITE(10) ((REAL(ddump(ji,jj)),ji=1,ns),jj=1,nt)
     WRITE(10) ((rsigma0(ji,jj),    ji=1,ns),jj=1,nt)
     WRITE(10) ((rsigma2(ji,jj),    ji=1,ns),jj=1,nt)
     WRITE(10) ((rsigma4(ji,jj),    ji=1,ns),jj=1,nt)
     CLOSE(10)
  ENDIF

  PRINT *,' Done.'
END PROGRAM cdfcensus
