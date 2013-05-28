PROGRAM cdf_xtract_brokenline
   !!======================================================================
   !!                     ***  PROGRAM  cdf_xtract_brokenline  ***
   !!=====================================================================
   !!  ** Purpose : Extract temperature, Salinity and velocity components
   !!               along a broken line formed by various legs
   !!
   !!  ** Method  : A broken line is defined by various segments or leg.
   !!               Each leg is defined by its starting and endig point
   !!               given as geographical coordinates on standard input.
   !!
   !!
   !! History : 2.1  : 12/2009  : R. Dussin    : Original code
   !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
   !!           3.0  : 05/2013  : T. Penduff & R. Dussin  : Saving new variables
   !!           3.0  : 05/2013  : J.M. Molines : Code review, generalization 
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   routines      : description
   !!----------------------------------------------------------------------
============================
 WORK IN PROGRESS DO NOT USE 
============================
   USE cdfio
   USE cdftools
   USE modcdfnames
   !!----------------------------------------------------------------------
   !! CDFTOOLS_3.0 , MEOM 2011
   !! $Id: cdfovide.f90 539 2011-07-11 10:33:35Z molines $
   !! Copyright (c) 2011, J.-M. Molines
   !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
   !!----------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER(KIND=4) :: jleg,  jk,  jipt
   INTEGER(KIND=4) :: narg, iargc         ! command line
   INTEGER(KIND=4) :: npiglo, npjglo, npk ! size of the domain
   INTEGER(KIND=4) :: niter               ! number of iteration of the search procedure
   INTEGER(KIND=4) :: iimin, iimax, ijmin, ijmax
   INTEGER(KIND=4) :: ii, ij , ipoint
   INTEGER(KIND=4) :: iloop

   INTEGER(KIND=4) :: nsec=0 ! total number of points on the broken line
   INTEGER(KIND=4), DIMENSION (:), ALLOCATABLE :: isec, jsec ! indices des points a recuperer

   ! broken line definition
   INTEGER(KIND=4), PARAMETER :: nsta=5    ! number of intermediate points

   INTEGER(KIND=4), DIMENSION(nsta  ) :: iista, ijsta
   INTEGER(KIND=4), DIMENSION(nsta-1) :: ikeepn

   ! broken line stuff
   INTEGER(KIND=4), PARAMETER :: jpseg=10000
   INTEGER(KIND=4) :: i0,j0,i1,j1, i, j
   INTEGER(KIND=4) :: n,nn, jseg, kk
   INTEGER(KIND=4) :: norm_u, norm_v, ist, jst
   INTEGER(KIND=4) :: nxtarg
   INTEGER(KIND=4), DIMENSION(nsta-1,jpseg) :: xlegs=0, ylegs=0

   REAL(KIND=4), DIMENSION(nsta)             :: rlonsta, rlatsta
   REAL(KIND=4)                              :: xmin, xmax, ymin, ymax, rdis
   REAL(KIND=4)                              :: glamfound, glamin, glamax
   REAL(KIND=8)                              :: glam0, emax

   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: glam, gphi, e1, e2
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1t, e2t, e1u, e2v, e3t

   REAL(KIND=4) :: rxi0, ryj0, rxi1, ryj1
   REAL(KIND=4) :: ai, bi, aj,bj,d
   REAL(KIND=4) :: rxx(jpseg), ryy(jpseg)
   REAL(KIND=4), DIMENSION(jpseg) :: gla !, gphi

   COMPLEX yypt(jpseg), yypti

   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1v, e3v ,gphiv, zv, zvt, zvs !: mask, metrics
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e2u, e3u ,gphiu, zu, zut, zus !: mask, metrics
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: temper, saline, uzonal, vmerid, rlon, rlat
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: ovidetemper, ovidesaline, ovidezonalu, ovidemeridv
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rlonsec, rlatsec, risec, rjsec
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1tsec, e1usec, e1vsec, e2tsec, e2usec, e2vsec
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e3tsec, e3usec, e3vsec
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e2join, e3join, rvmask
   REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: glamu, glamv
   REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: gdepw

   REAL(KIND=8)                              :: dtmp,dbarot

   CHARACTER(LEN=80) :: cf_tfil , cf_ufil, cf_vfil, csection 

   LOGICAL  :: lagain, lbord
   LOGICAL  :: ltest=.FALSE.
   LOGICAL  :: lchk

   ! cdf output stuff
   CHARACTER(LEN=80) :: cf_out='ovide.nc'
   TYPE (variable), DIMENSION(:), ALLOCATABLE   :: stypvar
   INTEGER(KIND=4)    :: ierr, ncout
   REAL(KIND=4), DIMENSION(1)                   ::  tim
   INTEGER(KIND=4)    :: nfield=14
   INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout
   !!----------------------------------------------------------------------
   CALL ReadCdfNames()

   narg = iargc()

   IF ( narg < 3 ) THEN
      PRINT *,' usage :  cdfovide T-file U-file V-file '
      PRINT *,'      '
      PRINT *,'     PURPOSE :'
      PRINT *,'       ' 
      PRINT *,'      '
      PRINT *,'     ARGUMENTS :'
      PRINT *,'       T-file :  model gridT file '
      PRINT *,'       U-file :  model gridU file '
      PRINT *,'       V-file :  model gridV file '
      PRINT *,'      ' 
      PRINT *,'     OPTIONS :'
      PRINT *,'       none ' 
      PRINT *,'      '
      PRINT *,'     REQUIRED FILES :'
      PRINT *,'       ', TRIM(cn_fhgr),' and ',TRIM(cn_fzgr),' must be in the current directory ' 
      PRINT *,'      '
      PRINT *,'     OUTPUT : '
      PRINT *,'       netcdf file : ', TRIM(cf_out) 
!     PRINT *,'         variables : ', TRIM(cv_out),' (    )'
      PRINT *,'      '
      PRINT *,'     SEE ALSO :'
      PRINT *,'      ' 
      PRINT *,'      '
      STOP
   ENDIF

   CALL getarg (1, cf_tfil)
   CALL getarg (2, cf_ufil)
   CALL getarg (3, cf_vfil)

   ! check file existence
   lchk = chkfile(cn_fhgr )
   lchk = chkfile(cn_fzgr ) .OR. lchk
   lchk = chkfile(cf_tfil ) .OR. lchk
   lchk = chkfile(cf_ufil ) .OR. lchk
   lchk = chkfile(cf_vfil ) .OR. lchk
   IF ( lchk ) STOP ! missing files

   ! R. Dussin : Location of leg points that define the 3 legs of OVIDE section
   !rlonsta(1) = -43.00 ; rlatsta(1) = 60.60    ! Greenland
   !rlonsta(2) = -31.30 ; rlatsta(2) = 58.90    ! Reykjanes Ridge
   !rlonsta(3) = -12.65 ; rlatsta(3) = 40.33    ! Off Portugal
   !rlonsta(4) =  -8.70 ; rlatsta(4) = 40.33    ! Lisboa

   ! D. Desbruyeres : Location of leg points that define the 4 legs of the OVIDE section
   rlonsta(1) = -43.70 ; rlatsta(1) = 59.90    ! 
   rlonsta(2) = -30.30 ; rlatsta(2) = 58.90    ! 
   rlonsta(3) = -19.40 ; rlatsta(3) = 44.90    ! 
   rlonsta(4) = -12.65 ; rlatsta(4) = 40.33    ! 
   rlonsta(5) = -08.70 ; rlatsta(5) = 40.33    ! 

   !!---------------------------------------------------------------------
   !!  Find the indexes of the legs (from cdffindij) 
   !!---------------------------------------------------------------------

   npiglo = getdim (cn_fhgr, cn_x)
   npjglo = getdim (cn_fhgr, cn_y)

   !! loop on the legs
   DO jleg = 1, nsta-1

      xmin = rlonsta(jleg)
      ymin = rlatsta(jleg)
      xmax = rlonsta(jleg+1)
      ymax = rlatsta(jleg+1)

      CALL cdf_findij ( xmin, xmax, ymin, ymax, iimin, iimax, ijmin, ijmax, &
           &            cd_coord=cn_fhgr, cd_point='F', cd_verbose='y')

      ! save leg information
      iista(jleg  ) = iimin
      ijsta(jleg  ) = ijmin
      iista(jleg+1) = iimax
      ijsta(jleg+1) = ijmax

      !! Find the broken line between P1 (iimin,ijmin) and P2 (iimax, ijmax)
      !! ---------------------------------------------------------------
      CALL broken_line( iimin, iimax, ijmin, ijmax, rxx, ryy, nn, npiglo, npjglo )

      IF (rxx(1) < rxx(nn) ) THEN ! leg is oriented eastward
         xlegs(jleg,:)=rxx
         ylegs(jleg,:)=ryy
      ELSE                        ! leg is oriented westward
         xlegs(jleg,:)=rxx(nn:1:-1)
         ylegs(jleg,:)=ryy(nn:1:-1)
      END IF

      ! compute the number of total points
      ikeepn(jleg) = nn  ! point on leg jleg
      nsec = nsec + nn
   END DO !! loop on the legs

   ! fancy control print
   WRITE(*,*) '------------------------------------------------------------'
   WRITE(*,9100) 'leg 1 start at ', rlonsta(1) ,'N ', rlatsta(1), 'W and ends at ', rlonsta(2) ,'N ', rlatsta(2), 'W'
   WRITE(*,9101) 'corresponding to F-gridpoints(', iista(1),',',ijsta(1),') and (', iista(2),',',ijsta(2),')' 
   WRITE(*,*) '------------------------------------------------------------'
   WRITE(*,*) '------------------------------------------------------------'
   WRITE(*,9100) 'leg 2 start at ', rlonsta(2) ,'N ', rlatsta(2), 'W and ends at ', rlonsta(3) ,'N ', rlatsta(3), 'W'
   WRITE(*,9101) 'corresponding to F-gridpoints(', iista(2),',',ijsta(2),') and (', iista(3),',',ijsta(3),')' 
   WRITE(*,*) '------------------------------------------------------------'
   WRITE(*,*) '------------------------------------------------------------'
   WRITE(*,9100) 'leg 3 start at ', rlonsta(3) ,'N ', rlatsta(3), 'W and ends at ', rlonsta(4) ,'N ', rlatsta(4), 'W'
   WRITE(*,9101) 'corresponding to F-gridpoints(', iista(3),',',ijsta(3),') and (', iista(4),',',ijsta(4),')' 
   WRITE(*,*) '------------------------------------------------------------'
   WRITE(*,*) '------------------------------------------------------------'
   WRITE(*,9100) 'leg 4 start at ', rlonsta(4) ,'N ', rlatsta(4), 'W and ends at ', rlonsta(5) ,'N ', rlatsta(5), 'W'
   WRITE(*,9101) 'corresponding to F-gridpoints(', iista(4),',',ijsta(4),') and (', iista(5),',',ijsta(5),')' 
   WRITE(*,*) '------------------------------------------------------------'

9100 FORMAT(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a)
9101 FORMAT(a,i4,a,i4,a,i4,a,i4,a)

   ALLOCATE (isec(nsec), jsec(nsec)) 

   ipoint = 0
   DO jleg=1, nsta-1  ! loop on legs 
      DO jipt=1, ikeepn(jleg) 
         ipoint = ipoint + 1
         isec(ipoint)=xlegs(jleg,jipt)  ! i-index
         jsec(ipoint)=ylegs(jleg,jipt)  ! j-index
      END DO
   END DO

   npk = getdim (cf_tfil, cn_z)

   ! input fields
   ALLOCATE(rlon(npiglo,npjglo), rlat(npiglo,npjglo))
   ALLOCATE(temper(npiglo,npjglo), saline(npiglo,npjglo))
   ALLOCATE(uzonal(npiglo,npjglo), vmerid(npiglo,npjglo))
   ALLOCATE(e1v(npiglo,npjglo))
   ALLOCATE(e2u(npiglo,npjglo))
   ALLOCATE(e3u(npiglo,npjglo), e3v(npiglo,npjglo))

   ! output fields
   ALLOCATE(rlonsec(1,nsec), rlatsec(1,nsec) )
   ALLOCATE(risec(1,nsec), rjsec(1,nsec) )
   ALLOCATE(e2usec(1,nsec-1), e3usec(nsec-1,npk) )
   ALLOCATE(e1vsec(1,nsec-1), e3vsec(nsec-1,npk) )
   ALLOCATE(ovidetemper(nsec-1,npk), ovidesaline(nsec-1,npk) )
   ALLOCATE(ovidezonalu(nsec-1,npk), ovidemeridv(nsec-1,npk) )
   ALLOCATE(rvmask(nsec-1,npk))

   e1vsec = -9999.
   e2usec = -9999.

   risec(:,:) = 0.
   rjsec(:,:) = 0.

   rlon(:,:) = getvar(cn_fhgr, cn_glamt, 1, npiglo, npjglo)
   rlat(:,:) = getvar(cn_fhgr, cn_gphit, 1, npiglo, npjglo)
   e1v(:,:)  = getvar(cn_fhgr, cn_ve1v,  1, npiglo, npjglo)
   e2u(:,:)  = getvar(cn_fhgr, cn_ve2u,  1, npiglo, npjglo)

   ! il faut faire un test sur la continuité des segments
   ! on va prendre T et S comme etant la moyenne du point
   ! en dessous et au-dessus du segment pour pouvoir calculer
   ! les fluxs de maniere optimales...

   ! loop on 2d arrays
   DO jipt = 1,nsec
      ii = isec(jipt)
      ij = jsec(jipt)

      risec  (1,jipt) = ii
      rjsec  (1,jipt) = ij
      rlonsec(1,jipt) = rlon(ii, ij)
      rlatsec(1,jipt) = rlat(ii, ij)
   END DO

   jk=1

   DO iloop=1,nsec-1
      !PRINT*, 'iloop=', iloop

      IF ( jsec(iloop+1) == jsec(iloop) ) THEN ! horizontal segment
         IF ( isec(iloop+1) > isec(iloop) ) THEN ! eastward

            e2usec(jk,iloop) = 0.
            e1vsec(jk,iloop) = e1v(isec(iloop)+1,jsec(iloop))

         ELSE

            e2usec(jk,iloop) = 0.
            e1vsec(jk,iloop) = e1v(isec(iloop),jsec(iloop))

         ENDIF
      ELSEIF ( isec(iloop+1) == isec(iloop) ) THEN ! vertical segment
         IF ( jsec(iloop+1) < jsec(iloop) ) THEN ! southward

            e2usec(jk,iloop) = e2u(isec(iloop),jsec(iloop))
            e1vsec(jk,iloop) = 0.

         ELSE

            e2usec(jk,iloop) = e2u(isec(iloop),jsec(iloop)+1)
            e1vsec(jk,iloop) = 0.

         ENDIF
      ELSE
         PRINT *, 'problem'
         exit 
      ENDIF
   END DO

   !	PRINT*,nsec
   !	PRINT*, MINVAL(e1v),MAXVAL(e1v),MINVAL(e2u),MAXVAL(e2u)  
   !	PRINT*, MINVAL(e1vsec),MAXVAL(e1vsec),MINVAL(e2usec),MAXVAL(e2usec)  
   !	PAUSE		
   ! loop on 3d arrays
   DO jk=1,npk
      temper(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo)
      saline(:,:) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo)
      uzonal(:,:) = getvar(cf_ufil, cn_vozocrtx, jk, npiglo, npjglo)
      vmerid(:,:) = getvar(cf_vfil, cn_vomecrty, jk, npiglo, npjglo)
      e3u(:,:)    = getvar(cn_fzgr, 'e3u_ps',    jk, npiglo, npjglo, ldiom=.true.)
      e3v(:,:)    = getvar(cn_fzgr, 'e3v_ps',    jk, npiglo, npjglo, ldiom=.true.)

      DO iloop=1,nsec-1
         IF ( jsec(iloop+1) == jsec(iloop) ) THEN ! horizontal segment
            IF ( isec(iloop+1) > isec(iloop) ) THEN ! eastward

               IF ( min( temper(isec(iloop)+1,jsec(iloop)) , temper(isec(iloop)+1,jsec(iloop)+1) ) == 0. ) THEN
                  ovidetemper(iloop,jk) = 0. ; ovidesaline(iloop,jk) = 0.
               ELSE
                  ovidetemper(iloop,jk) = 0.5 * ( temper(isec(iloop)+1,jsec(iloop)) + temper(isec(iloop)+1,jsec(iloop)+1) )
                  ovidesaline(iloop,jk) = 0.5 * ( saline(isec(iloop)+1,jsec(iloop)) + saline(isec(iloop)+1,jsec(iloop)+1) )
               ENDIF
               ovidezonalu(iloop,jk) = 0.
               ovidemeridv(iloop,jk) = vmerid(isec(iloop)+1,jsec(iloop))
               e3usec(iloop,jk) = 0.
               e3vsec(iloop,jk) = e3v(isec(iloop)+1,jsec(iloop))

            ELSE ! westward

               IF ( min( temper(isec(iloop),jsec(iloop)) , temper(isec(iloop),jsec(iloop)+1) ) == 0. ) THEN
                  ovidetemper(iloop,jk) = 0. ; ovidesaline(iloop,jk) = 0.
               ELSE
                  ovidetemper(iloop,jk) = 0.5 * ( temper(isec(iloop),jsec(iloop)) + temper(isec(iloop),jsec(iloop)+1) )
                  ovidesaline(iloop,jk) = 0.5 * ( saline(isec(iloop),jsec(iloop)) + saline(isec(iloop),jsec(iloop)+1) )
               ENDIF
               ovidezonalu(iloop,jk) = 0.
               ovidemeridv(iloop,jk) = vmerid(isec(iloop),jsec(iloop))
               e3usec(iloop,jk) = 0.
               e3vsec(iloop,jk) = e3v(isec(iloop),jsec(iloop))

            ENDIF
         ELSEIF ( isec(iloop+1) == isec(iloop) ) THEN ! vertical segment
            IF ( jsec(iloop+1) < jsec(iloop) ) THEN ! southward

               IF ( min( temper(isec(iloop),jsec(iloop)) , temper(isec(iloop)+1,jsec(iloop)) ) == 0. ) THEN
                  ovidetemper(iloop,jk) = 0. ; ovidesaline(iloop,jk) = 0.
               ELSE
                  ovidetemper(iloop,jk) = 0.5 * ( temper(isec(iloop),jsec(iloop)) + temper(isec(iloop)+1,jsec(iloop)) )
                  ovidesaline(iloop,jk) = 0.5 * ( saline(isec(iloop),jsec(iloop)) + saline(isec(iloop)+1,jsec(iloop)) )
               ENDIF
               ovidezonalu(iloop,jk) = uzonal(isec(iloop),jsec(iloop))
               ovidemeridv(iloop,jk) = 0.
               e3usec(iloop,jk) = e3u(isec(iloop),jsec(iloop))
               e3vsec(iloop,jk) = 0.

            ELSE ! northward

               IF ( min( temper(isec(iloop),jsec(iloop)+1) , temper(isec(iloop)+1,jsec(iloop)+1) ) == 0. ) THEN
                  ovidetemper(iloop,jk) = 0. ; ovidesaline(iloop,jk) = 0.
               ELSE
                  ovidetemper(iloop,jk) = 0.5 * ( temper(isec(iloop),jsec(iloop)+1) + temper(isec(iloop)+1,jsec(iloop)+1) )
                  ovidesaline(iloop,jk) = 0.5 * ( saline(isec(iloop),jsec(iloop)+1) + saline(isec(iloop)+1,jsec(iloop)+1) )
               ENDIF
               ovidezonalu(iloop,jk) = uzonal(isec(iloop),jsec(iloop)+1)
               ovidemeridv(iloop,jk) = 0.
               e3usec(iloop,jk) = e3u(isec(iloop),jsec(iloop)+1)
               e3vsec(iloop,jk) = 0.

            ENDIF

         ELSE
            PRINT *, 'problem'
            exit 
         ENDIF

      END DO
   END DO


   ALLOCATE ( stypvar(nfield), ipk(nfield), id_varout(nfield) )

   DO iloop=1,nfield
      ipk(iloop) = npk
   END DO

   ! define new variables for output 
   stypvar(1)%cname= 'votemper'
   stypvar(1)%cunits='deg C'
   stypvar%rmissing_value=0.
   stypvar(1)%valid_min= -2.
   stypvar(1)%valid_max= 40.
   stypvar%scale_factor= 1.
   stypvar%add_offset= 0.
   stypvar%savelog10= 0.
   stypvar(1)%clong_name='Temperature along OVIDE section'
   stypvar(1)%cshort_name='votemper'
   stypvar%conline_operation='N/A'
   stypvar%caxis='TXZ'

   stypvar(2)%cname= 'vosaline'
   stypvar(2)%cunits='PSU'
   stypvar(2)%valid_min= 0.
   stypvar(2)%valid_max= 50.
   stypvar(2)%clong_name='Salinity along OVIDE section'
   stypvar(2)%cshort_name='vosaline'

   stypvar(3)%cname= 'vozocrtx_native'
   stypvar(3)%cunits='m.s-1'
   stypvar(3)%valid_min= -20.
   stypvar(3)%valid_max= 20.
   stypvar(3)%clong_name='Zonal velocity along OVIDE section'
   stypvar(3)%cshort_name='vozocrtx'

   stypvar(4)%cname= 'vomecrty_native'
   stypvar(4)%cunits='m.s-1'
   stypvar(4)%valid_min= -20.
   stypvar(4)%valid_max= 20.
   stypvar(4)%clong_name='Meridionnal velocity along OVIDE section'
   stypvar(4)%cshort_name='vomecrty'

   stypvar(5)%cname= 'isec'
   stypvar(5)%valid_min= 0.
   stypvar(5)%valid_max= npiglo 

   stypvar(6)%cname= 'jsec'
   stypvar(6)%valid_min= 0.
   stypvar(6)%valid_max= npjglo 

   stypvar(7)%cname= 'e2u_native'
   stypvar(7)%valid_min= MINVAL(e2usec(1,:))
   stypvar(7)%valid_max= MAXVAL(e2usec(1,:)) 

   stypvar(8)%cname= 'e1v_native'
   stypvar(8)%valid_min= MINVAL(e1vsec(1,:))
   stypvar(8)%valid_max= MAXVAL(e1vsec(1,:))

   stypvar(9)%cname= 'e3u_native'
   stypvar(9)%valid_min= MINVAL(e3usec(:,:))
   stypvar(9)%valid_max= MAXVAL(e3usec(:,:)) 

   stypvar(10)%cname= 'e3v_native'
   stypvar(10)%valid_min= MINVAL(e3vsec(:,:))
   stypvar(10)%valid_max= MAXVAL(e3vsec(:,:)) 

   stypvar(11)%cname= 'vomecrty'
   stypvar(11)%cunits='m.s-1'
   stypvar(11)%valid_min= -20.
   stypvar(11)%valid_max= 20.
   stypvar(11)%clong_name='Normal velocity along OVIDE section'
   stypvar(11)%cshort_name='vomecrty'

   stypvar(12)%cname= 'e1v'
   stypvar(12)%cunits='m'
   stypvar(12)%valid_min= 0.
   stypvar(12)%valid_max= 1000000.
   stypvar(12)%clong_name='Local horiz. resolution along OVIDE section'
   stypvar(12)%cshort_name='e1v'

   stypvar(13)%cname= 'e3v_ps'
   stypvar(13)%cunits='m'
   stypvar(13)%valid_min= 0.
   stypvar(13)%valid_max= 100000000.
   stypvar(13)%clong_name='Local vert. resolution along OVIDE section'
   stypvar(13)%cshort_name='e3v_ps'

   stypvar(14)%cname= 'vmask'
   stypvar(12)%cunits=''
   stypvar(12)%valid_min= 0.
   stypvar(12)%valid_max= 1.
   stypvar(12)%clong_name='Mask along OVIDE section'
   stypvar(12)%cshort_name='vmask'

   ! create output fileset
   ncout =create(cf_out, 'none', nsec,1, npk,cdep='depthw')
   ierr= createvar(ncout ,stypvar,nfield, ipk,id_varout )
   ierr= putheadervar(ncout, cf_tfil,nsec,1, npk,pnavlon=rlonsec,pnavlat=rlatsec,pdep=gdepw)
   tim=getvar1d(cf_tfil,'time_counter',1)
   ierr=putvar1d(ncout,tim,1,'T')


   rvmask(:,:) = 1.
   WHERE( ovidesaline(:,:) == 0. ) rvmask(:,:) = 0.


   !PRINT*, MINVAL(e1v),MAXVAL(e1v),MINVAL(e2u),MAXVAL(e2u)  
   !PRINT*, MINVAL(e1vsec),MAXVAL(e1vsec),MINVAL(e2usec),MAXVAL(e2usec)  
   !PAUSE

   !------------------- BAROTROPIC TRANSPORT
   dbarot = 0.d0

   DO iloop=1,nsec-1
      DO jk=1,npk
         dtmp=1.d0* (ovidezonalu(iloop,jk) + ovidemeridv(iloop,jk))*&
              (e2usec(1,iloop )+ e1vsec(1,iloop ))*                 &
              (e3usec(iloop,jk)+ e3vsec(iloop,jk))*                 &
              rvmask(iloop,jk)
         dbarot=dbarot+dtmp
      ENDDO
      !jk=1
      !PRINT*,iloop,(ovidezonalu(iloop,jk)+ovidemeridv(iloop,jk)),(e2usec(1,iloop)+e1vsec(1,iloop)),&
      !(e3usec(iloop,jk)+e3vsec(iloop,jk)),rvmask(iloop,jk),dbarot
   ENDDO
   PRINT*, 'BAROTROPIC TRANSPORT = ', dbarot/1.e6, ' Sv.'
   !--------------------------------------------


   ! netcdf output 
   DO jk =1, npk
      ierr = putvar (ncout, id_varout(1), ovidetemper(:,jk), jk, nsec-1, 1 )
      ierr = putvar (ncout, id_varout(2), ovidesaline(:,jk), jk, nsec-1, 1 )
      ierr = putvar (ncout, id_varout(3), ovidezonalu(:,jk), jk, nsec-1, 1 )
      ierr = putvar (ncout, id_varout(4), ovidemeridv(:,jk), jk, nsec-1, 1 )
      ierr = putvar (ncout, id_varout(5), risec(1,:), jk,1,nsec)
      ierr = putvar (ncout, id_varout(6), rjsec(1,:), jk,1,nsec)
      ierr = putvar (ncout, id_varout(7), e2usec(1,:), jk,nsec-1,1)
      ierr = putvar (ncout, id_varout(8), e1vsec(1,:), jk,nsec-1,1)
      ierr = putvar (ncout, id_varout(9), e3usec(:,jk), jk,nsec-1,1)
      ierr = putvar (ncout, id_varout(10),e3vsec(:,jk), jk,nsec-1,1)
      ! along-track normal velocity, horiz. and vert. resolution, and mask
      ierr = putvar (ncout, id_varout(11),ovidezonalu(:,jk) + ovidemeridv(:,jk), jk,nsec-1,1)
      ierr = putvar (ncout, id_varout(12),e2usec(1,:) + e1vsec(1,:), jk,nsec-1,1)
      ierr = putvar (ncout, id_varout(13),e3usec(:,jk) + e3vsec(:,jk), jk,nsec-1,1)
      ierr = putvar (ncout, id_varout(14),rvmask(:,jk), jk,nsec-1,1)
   END DO

   ierr = closeout(ncout)

   !!--------------------------------------------------------------------
   !!
   !!   SUBROUTINES USED
   !!
   !!--------------------------------------------------------------------
CONTAINS

   SUBROUTINE interm_pt (ydpt,k,pai,pbi,paj,pbj,ydpti)
      !! -----------------------------------------------------
      !!           SUBROUTINE INTERM_PT
      !!           ********************
      !!
      !!   PURPOSE:
      !!   --------
      !!     Find the best intermediate points on a pathway.
      !!
      !!    ARGUMENTS:
      !!    ----------
      !!      ydpt : complex vector of the positions of the nearest points
      !!         k : current working index
      !!       pai ,pbi : slope and original ordinate of x(y)
      !!       paj ,pbj : slope and original ordinate of y(x)
      !!      ydpti : Complex holding the position of intermediate point
      !!
      !!    AUTHOR:
      !!    -------
      !!      19/07/1999 : Jean-Marc MOLINES
      !!      14/01/2005 : J M M in F90
      !!
      !!--------------------------------------------------------------
      !!
      !! 0. Declarations:
      !! ----------------
      IMPLICIT NONE
      COMPLEX, INTENT(in) :: ydpt(*)
      COMPLEX, INTENT(out) :: ydpti
      REAL(KIND=4), INTENT(IN) ::  pai,pbi,paj,pbj
      INTEGER(KIND=4) ,INTENT(in) :: k
      ! ... local
      COMPLEX      :: ylptmp1, ylptmp2
      REAL(KIND=4) ::  za0,zb0,za1,zb1,zd1,zd2
      REAL(KIND=4) ::  zxm,zym
      REAL(KIND=4) ::  zxp,zyp
      !!
      !! 1. Compute intermediate points
      !! ------------------------------
      !
      ! ... Determines whether we use y(x) or x(y):
      IF (ABS(paj) <=  1) THEN
         ! ..... y(x)
         ! ... possible intermediate points:
         ylptmp1=ydpt(k-1)+(1.,0.)
         ylptmp2=ydpt(k-1)+CMPLX(0.,SIGN(1.,paj))
         !
         ! ... M is the candidate point:
         zxm=REAL(ylptmp1)
         zym=IMAG(ylptmp1)
         za0=paj
         zb0=pbj
         !
         za1=-1./za0
         zb1=zym - za1*zxm
         ! ... P is the projection of M on the strait line
         zxp=-(zb1-zb0)/(za1-za0)
         zyp=za0*zxp + zb0
         ! ... zd1 is the distance MP
         zd1=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
         !
         ! ... M is the candidate point:
         zxm=REAL(ylptmp2)
         zym=IMAG(ylptmp2)
         za1=-1./za0
         zb1=zym - za1*zxm
         ! ... P is the projection of M on the strait line
         zxp=-(zb1-zb0)/(za1-za0)
         zyp=za0*zxp + zb0
         ! ... zd2 is the distance MP
         zd2=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
         ! ... chose the smallest (zd1,zd2)
         IF (zd2 <=  zd1) THEN
            ydpti=ylptmp2
         ELSE
            ydpti=ylptmp1
         END IF
         !
      ELSE
         !
         ! ... x(y)
         ylptmp1=ydpt(k-1)+CMPLX(SIGN(1.,pai),0.)
         ylptmp2=ydpt(k-1)+(0.,1.)
         zxm=REAL(ylptmp1)
         zym=IMAG(ylptmp1)
         za0=pai
         zb0=pbi
         !
         za1=-1./za0
         zb1=zxm - za1*zym
         zyp=-(zb1-zb0)/(za1-za0)
         zxp=za0*zyp + zb0
         zd1=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
         !
         zxm=REAL(ylptmp2)
         zym=IMAG(ylptmp2)
         za1=-1./za0
         zb1=zxm - za1*zym
         zyp=-(zb1-zb0)/(za1-za0)
         zxp=za0*zyp + zb0
         zd2=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
         IF (zd2 <=  zd1) THEN
            ydpti=ylptmp2
         ELSE
            ydpti=ylptmp1
         END IF
      END IF
   END SUBROUTINE interm_pt


   SUBROUTINE  broken_line( kimin, kimax, kjmin, kjmax, pxx, pyy, knn , kpi, kpj)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE broken_line  ***
      !!
      !! ** Purpose :  determine the broken line between 2 points 
      !!
      !! ** Method  :   
      !!
      !!----------------------------------------------------------------------
      INTEGER(KIND=4),            INTENT(in ) :: kimin, kimax, kjmin, kjmax
      INTEGER(KIND=4),            INTENT(out) :: knn
      INTEGER(KIND=4),            INTENT(in ) :: kpi, kpj
      REAL(KIND=4), DIMENSION(kpi+kpj), INTENT(out) :: pxx, pyy

      INTEGER(KIND=4) :: ji, jj   ! dummy loop index
      INTEGER(KIND=4) :: ii0, ij0, ii1, ij1, iitmp, ijtmp
      INTEGER(KIND=4) :: iist, ijst
      INTEGER(KIND=4) :: ii, ij
      INTEGER(KIND=4) :: idirx, idiry
      INTEGER(KIND=4) :: ipts

      REAL(KIND=4)    :: zxi0, zyj0, zxi1, zyj1
      REAL(KIND=4)    :: zai, zbi
      REAL(KIND=4)    :: zaj, zbj
      REAL(KIND=4)    :: zd

      COMPLEX, DIMENSION(kpi+kpj) :: ylpt           ! array of points coordinates in a section
      COMPLEX                     :: ylpti          ! working point
      !!----------------------------------------------------------------------

      !! Find the broken line between P1 (kimin,kjmin) and P2 (kimax, kjmax)
      ! ... Initialization
      ii0  = kimin ; ij0  = kjmin ; ii1  = kimax ;  ij1 = kjmax
      zxi0 = ii0   ; zyj0 = ij0   ; zxi1 = ii1   ; zyj1 = ij1

      ! compute direction of integrations and signs
      !The transport across the section is the dot product of
      !integral(line){(Mx,My)*dS}
      !Mx=integral(u*dz)  My=integral(v*dz)) and dS=(dy,-dx)}

      !By defining the direction of the integration as
      idirx = SIGN(1,ii1-ii0) !positive to the east or if ii1=ii0
      idiry = SIGN(1,ij1-ij0) !positive to the north or if ij1=ij0

      !Then dS=(e2u*idiry,-e1v*idirx)
      !This will produce the following sign convention:
      !    West-to-est line (dx>0, dy=0)=> -My*dx (-ve for a northward flow)
      !    South-to-north   (dy>0, dx=0)=>  Mx*dy (+ve for an eastward flow)
      norm_u =  idiry
      norm_v = -idirx

      ! .. Compute equation:  ryj = zaj rxi + zbj [valid in the (i,j) plane]
      IF ( (zxi1 -zxi0) /=  0 ) THEN
         zaj = (zyj1 - zyj0 ) / (zxi1 -zxi0)
         zbj = zyj0 - zaj * zxi0
      ELSE
         zaj = 10000.  ! flag value
         zbj = 0.
      END IF
      ! .. Compute equation:  rxi = zai ryj + zbi [valid in the (i,j) plane]
      IF ( (zyj1 -zyj0) /=  0 ) THEN
         zai = (zxi1 - zxi0 ) / ( zyj1 -zyj0 )
         zbi = zxi0 - zai * zyj0
      ELSE
         zai = 10000. ! flag value
         zbi = 0.
      END IF

      ! ..  Compute the integer pathway: a succession of F points
      ipts=0
      ! .. Chose the strait line with the smallest slope
      IF (ABS(zaj) <=  1 ) THEN
         ! ... Here, the best line is y(x)
         ! ... If ii1 < ii0 swap points [ always describe section from left to right ]
         IF (ii1 <  ii0 ) THEN
            iitmp = ii0   ; ijtmp = ij0
            ii0   = ii1   ; ij0   = ij1
            ii1   = iitmp ; ij1   = ijtmp
         END IF

         ! iist,ijst is the grid offset to pass from F point to either U/V point
         IF ( ij1 >= ij0 ) THEN     ! line heading NE
            iist = 1 ; ijst = 1
         ELSE                       ! line heading SE
            iist = 1 ; ijst = 0
         END IF

         ! ... compute the nearest ji point on the line crossing at ji
         DO ji=ii0, ii1
            ipts=ipts+1
            IF (ipts > kpi+kpj) STOP 'ipts > kpi+kpj !'
            ij=NINT(zaj*ji + zbj )
            ylpt(ipts) = CMPLX(ji,ij)
         END DO
      ELSE
         ! ... Here, the best line is x(y)
         ! ... If ij1 < ij0 swap points [ always describe section from bottom to top ]
         IF (ij1 <  ij0 ) THEN
            iitmp = ii0   ; ijtmp = ij0
            ii0   = ii1   ; ij0   = ij1
            ii1   = iitmp ; ij1   = ijtmp
         END IF

         ! iist,ijst is the grid offset to pass from F point to either U/V point
         IF ( ii1 >=  ii0 ) THEN
            iist = 1 ; ijst = 1
         ELSE
            iist = 0 ; ijst = 1
         END IF

         ! ... compute the nearest ji point on the line crossing at jj
         DO jj=ij0,ij1
            ipts=ipts+1
            IF (ipts > kpi+kpj) STOP 'ipts > kpi+kpj !'
            ii=NINT(zai*jj + zbi)
            ylpt(ipts) = CMPLX(ii,jj)
         END DO
      END IF
      !!
      !! Look for intermediate points to be added.
      !  ..  The final positions are saved in pxx,pyy
      pxx(1) = REAL(ylpt(1))
      pyy(1) = IMAG(ylpt(1))
      knn     = 1

      DO jk=2,ipts
         ! .. distance between 2 neighbour points
         zd=ABS(ylpt(jk)-ylpt(jk-1))
         ! .. intermediate points required if zd > 1
         IF ( zd > 1 ) THEN
            CALL interm_pt(ylpt, jk, zai, zbi, zaj, zbj, ylpti)
            knn=knn+1
            IF (knn > kpi+kpj) STOP 'knn>kpi+kpj !'
            pxx(knn) = REAL(ylpti)
            pyy(knn) = IMAG(ylpti)
         END IF
         knn=knn+1
         IF (knn > kpi+kpj) STOP 'knn>kpi+kpj !'
         pxx(knn) = REAL(ylpt(jk))
         pyy(knn) = IMAG(ylpt(jk))
      END DO

   END SUBROUTINE broken_line


END PROGRAM cdf_xtract_brokenline
