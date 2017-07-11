MODULE modutils
  !!======================================================================
  !!                     ***  MODULE  modutils  ***
  !! Hold functions and subroutine dedicated to common utility task
  !!=====================================================================
  !! History : 3.0  : 04/2011  : J.M. Molines : Original code
  !!                : 10/2012  : N. Ferry, E. Durand, F. Hernandez : add shapiro
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   SetGlobalAtt  : Set Global Attribute to the command line
  !!   SetFilename   : Build standard name from confname
  !!   shapiro_fill_smooth : shapiro smoother or filler
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class system
  !!----------------------------------------------------------------------
  USE cdfio

  IMPLICIT NONE

  PRIVATE
  PUBLIC SetGlobalAtt
  PUBLIC SetFileName
  PUBLIC GetList
  PUBLIC shapiro_fill_smooth
  PUBLIC FillPool2D
  PUBLIC FillPool3D
  PUBLIC heading         ! compute true heading between point A and B

CONTAINS
  SUBROUTINE SetGlobalAtt(cdglobal, cd_append)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE SetGlobalAtt  ***
    !!
    !! ** Purpose : Append command line to the string given as argument.
    !!              This is basically used for setting a global attribute 
    !!              in the output files 
    !!
    !! ** Method  : Decrypt line command with getarg  
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),           INTENT(inout) :: cdglobal
    CHARACTER(LEN=1), OPTIONAL, INTENT(in   ) :: cd_append

    INTEGER(KIND=4)    :: iargc, inarg
    INTEGER(KIND=4)    :: jarg
    CHARACTER(LEN=100) :: cl_arg
    CHARACTER(LEN=1  ) :: cl_app
    !!----------------------------------------------------------------------
    cl_app = 'N'
    IF ( PRESENT( cd_append ) ) THEN 
       cl_app = 'A'
    ENDIF

    CALL getarg(0, cl_arg)
    SELECT CASE ( cl_app)
    CASE ('A') 
       cdglobal = TRIM(cdglobal)//' ; '//TRIM(cl_arg) 
    CASE ('N') 
       cdglobal = TRIM(cl_arg) 
    END SELECT

    inarg = iargc()
    DO jarg=1, inarg
       CALL getarg(jarg,cl_arg) 
       cdglobal = TRIM(cdglobal)//' '//TRIM(cl_arg) 
    END DO

  END SUBROUTINE SetGlobalAtt

  CHARACTER(LEN=256) FUNCTION SetFileName(cdconf, cdtag, cdgrid ,ld_stop )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION SetFileName  ***
    !!
    !! ** Purpose :  Build filename from cdconf, tag and grid
    !!
    !! ** Method  :  Check 2 forms of file names and return
    !!               error is file is missing
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),  INTENT(in) :: cdconf, cdtag, cdgrid
    LOGICAL, OPTIONAL, INTENT(in) :: ld_stop

    LOGICAL :: ll_stop
    !!----------------------------------------------------------------------
    IF ( PRESENT(ld_stop) ) THEN
       ll_stop = ld_stop
    ELSE
       ll_stop = .TRUE.
    ENDIF

    WRITE( SetFileName,'(a,"_",a,"_grid",a,".nc")') TRIM(cdconf), TRIM(cdtag), TRIM(cdgrid)
    IF ( chkfile(SetFileName ,ld_verbose=.FALSE.) ) THEN ! look for another name
       WRITE(SetFileName,'(a,"_",a,"_grid_",a,".nc")') TRIM(cdconf), TRIM(cdtag), TRIM(cdgrid)
       IF ( chkfile( SetFileName, ld_verbose=.FALSE.) .AND. ll_stop  ) THEN
          PRINT *,' ERROR : missing grid',TRIM(cdgrid),'or even grid_',TRIM(cdgrid),' file '
          STOP 97
       ENDIF
    ENDIF
  END FUNCTION SetFileName

  SUBROUTINE GetList ( cd_list, klist, ksiz )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE getlist  ***
    !!
    !! ** Purpose :   Expand list described with input string like
    !!                k1,k2,k3  or k1-k2,k3 or any valid combination
    !!
    !! ** Method  :   Look for ',' and '-' in cd_list and interprets
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),                           INTENT(in ) :: cd_list ! list to decipher
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE, INTENT(out) :: klist   ! deciphered integer list
    INTEGER(KIND=4),                            INTENT(out) :: ksiz    ! size of the list

    INTEGER(KIND=4),             PARAMETER :: jp_maxlist=500
    INTEGER(KIND=4)                        :: ji, jk
    INTEGER(KIND=4)                        :: inlev, ipos, iposm,  ik1, ik2
    INTEGER(KIND=4), DIMENSION(jp_maxlist) :: itmp
    CHARACTER(LEN=80)                      :: cldum
    CHARACTER(LEN=80)                      :: cldum2
    !----------------------------------------------------------------------------
    cldum=cd_list
    ipos=1
    ksiz=0
    DO WHILE (ipos /= 0 )
       ipos=INDEX(cldum,',')
       IF (ipos == 0 ) THEN
          cldum2=cldum
       ELSE
          cldum2=cldum(1:ipos-1)
       ENDIF

       iposm=INDEX(cldum2,'-')
       IF ( iposm == 0 ) THEN
          ksiz=ksiz+1 
          IF (ksiz > jp_maxlist) THEN ; PRINT *'jp_maxlist too small in getlist ' ; STOP 97
          ENDIF
          READ(cldum2,* ) itmp(ksiz)
       ELSE
          READ(cldum2(1:iposm-1),*) ik1
          READ(cldum2(iposm+1:),* ) ik2
          DO jk = ik1,ik2
             ksiz=ksiz+1 
             IF (ksiz > jp_maxlist) THEN ; PRINT *'jp_maxlist too small in getlist ' ; STOP 97
             ENDIF
             itmp(ksiz)=jk
          ENDDO
       ENDIF
       cldum=cldum(ipos+1:)
    ENDDO
    ALLOCATE (klist(ksiz) )
    klist(:)=itmp(1:ksiz)

  END SUBROUTINE GetList

  SUBROUTINE shapiro_fill_smooth ( psig, kpi, kpj, kpass, cdfs, pbad, klmasktrue, psigf )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE shapiro_fill_smooth  ***
    !!
    !! ** Purpose : Shapiro smoother or filler
    !!
    !! ** Method  : Shapiro algorithm 
    !!           psig    : variable to be filtered 2D
    !!           kpi,kpj : dimension of psig
    !!           kpass   : number of passes of the filter
    !!           cdfs    : 'smooth' or 'fill' according to choice
    !!           pbad    : psig Fill_Value
    !!           klmasktrue : mask flag for continent.
    !!                If land extrapolation is desired, set klmasktrue=1 everywhere
    !!
    !!           psigf   : filtered/filled variable (output)
    !!
    !!  code history:
    !!      original  : 05-11 (N. Ferry)
    !!      additions : 05-12 (E. Durand)
    !!      correction: 07-12 (F. Hernandez)
    !!      cdftools norm : 11-12 (J.M. Molines)
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),                     INTENT(in ) :: kpi, kpj, kpass
    INTEGER(KIND=4), DIMENSION(kpi,kpj), INTENT(in ) :: klmasktrue

    REAL(KIND=4),                        INTENT(in ) :: pbad
    REAL(KIND=4), DIMENSION(kpi,kpj),    INTENT(in ) :: psig
    REAL(KIND=4), DIMENSION(kpi,kpj),    INTENT(out) :: psigf

    CHARACTER(LEN=6),                    INTENT(in ) :: cdfs

    INTEGER(KIND=4)                               :: ji, jj, jp    ! dummy loop index
    INTEGER(KIND=4), DIMENSION(0:kpi+1,kpj)       :: ilmask_e     ! extra i-point for E-W periodicity
    INTEGER(KIND=4), DIMENSION(0:kpi+1,kpj)       :: ilmask0_e    ! extra i-point for E-W periodicity
    INTEGER(KIND=4), DIMENSION(0:kpi+1,kpj)       :: ilmasktrue_e ! extra i-point for E-W periodicity

    REAL(KIND=4), DIMENSION(0:kpi+1,kpj)          :: zsigf_e      ! extra i-point for E-W periodicity
    REAL(KIND=4), DIMENSION(0:kpi+1,kpj)          :: zsig_e       ! extra i-point for E-W periodicity
    REAL(KIND=4)                                  :: znum, zden, zsum

    !!----------------------------------------------------------------------
    ! ... Initialization : 
    zsig_e      (1:kpi,:) = psig      (:,:)
    ilmasktrue_e(1:kpi,:) = klmasktrue(:,:)
    !  E-W periodic
    zsig_e      (0,:)     = zsig_e      (kpi,:)      
    ilmasktrue_e(0,:)     = ilmasktrue_e(kpi,:)      
    zsig_e      (kpi+1,:) = zsig_e      (1,:)      
    ilmasktrue_e(kpi+1,:) = ilmasktrue_e(1,:)      

    ! check cdfs compliance
    IF ( cdfs(1:4)  /= 'fill' .AND. cdfs(1:6) /= 'smooth' ) THEN
       PRINT*, 'cdfs = ',cdfs ,' <> fill or smooth'
       STOP 97
    ENDIF
    !
    ! ... Shapiro filter : 
    !
    DO jp = 1, kpass          ! number of passes for the filter
       !
       ! in both cases 'smooth' and ' fill' we check points w/o values
       ilmask_e(:,:) = 0 ; ilmask0_e(:,:) = 0
       WHERE ( zsig_e(:,:) /= pbad )
          !   set ilmask_e = 1 when field is already filled
          ilmask_e (:,:) = 1 
          ilmask0_e(:,:) = 1 
       ENDWHERE

       ! case 'fill'
       IF ( cdfs(1:4) == 'fill' ) THEN
          ilmask0_e(:,:) = 0
          DO ji=1,kpi
             DO jj=2,kpj-1
                zsum = ilmask_e(ji+1,jj) + ilmask_e(ji-1,jj) + ilmask_e(ji,jj+1) + ilmask_e(ji,jj-1)
                ! set ilmask0_e = 1 if it is possible to do a 4-point interpolation (N-S-E-W)
                ! not on  land
                IF ( ( zsum                >= 1 ) .AND. &
                     ( ilmask_e    (ji,jj) == 0 ) .AND. &
                     ( ilmasktrue_e(ji,jj) == 1 ) )  THEN
                          ilmask0_e(ji,jj) = 1
                ENDIF
             ENDDO
             ! for the northernmost line
             zsum = ilmask_e(ji+1,kpj) + ilmask_e(ji-1,kpj) + ilmask_e(ji,kpj-1)
             IF ( ( zsum                 >= 1 ) .AND. &
                  ( ilmask_e    (ji,kpj) == 0 ) .AND. &
                  ( ilmasktrue_e(ji,kpj) == 1 ) )  THEN 
                       ilmask0_e(ji,kpj) = 1
             ENDIF
          ENDDO
       ENDIF
       !
       ! loop on data points for both cases
       DO ji = 1, kpi
          DO jj = 2, kpj-1
             IF ( ilmask0_e(ji,jj) == 1. )  THEN
                znum =  zsig_e(ji-1,jj  )*ilmask_e(ji-1,jj  ) &
                      + zsig_e(ji+1,jj  )*ilmask_e(ji+1,jj  ) &
                      + zsig_e(ji  ,jj-1)*ilmask_e(ji  ,jj-1) &
                      + zsig_e(ji  ,jj+1)*ilmask_e(ji  ,jj+1)  
                zden =  ilmask_e(ji-1,jj  ) &
                      + ilmask_e(ji+1,jj  ) &
                      + ilmask_e(ji  ,jj-1) &
                      + ilmask_e(ji  ,jj+1) 
                zsigf_e(ji,jj) = znum/zden
             ELSE
                zsigf_e(ji,jj) = zsig_e(ji,jj)
             ENDIF
          ENDDO
          ! for the northernmost line, we do not take kpj+1 into account
          IF ( ilmask0_e(ji,kpj) == 1. )  THEN
             znum =  zsig_e(ji-1,kpj  )*ilmask_e(ji-1,kpj  ) &
                   + zsig_e(ji+1,kpj  )*ilmask_e(ji+1,kpj  ) &
                   + zsig_e(ji  ,kpj-1)*ilmask_e(ji  ,kpj-1) 
             zden =  ilmask_e(ji-1,kpj  ) &
                   + ilmask_e(ji+1,kpj  ) &
                   + ilmask_e(ji  ,kpj-1) 
             zsigf_e(ji,kpj) = znum/zden
          ELSE
             zsigf_e(ji,kpj) = zsig_e(ji,kpj)
          ENDIF
       ENDDO
       !
       !    fill or smooth ?
       !
       IF ( cdfs(1:6) == 'smooth' ) THEN
          WHERE ( ilmasktrue_e(:,:) == 1 )
             zsig_e(:,:) = zsigf_e(:,:)
          END WHERE
       ENDIF
       !
       IF ( cdfs(1:4) == 'fill' ) THEN
          WHERE ( ilmask0_e(:,:) == 1 )
             zsig_e(:,:) = zsigf_e(:,:)
          END WHERE
       ENDIF
       ! Boundary condition  : E-W  (simplifie)
       zsig_e(0,:) = zsig_e(kpi,:)
       zsig_e(kpi+1,:) = zsig_e(1,:)

       !
    ENDDO                     ! jp

    psigf(:,:) = zsig_e(1:kpi,:)

  END SUBROUTINE shapiro_fill_smooth

  SUBROUTINE FillPool2D(kiseed, kjseed, kdta, kifill)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE FillPool2D  ***
    !!  
    !! ** Purpose :  Replace all area surrounding by mask value by mask value
    !!  
    !! ** Method  :  flood fill algorithm
    !!  
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),                 INTENT(in)    :: kiseed, kjseed
    INTEGER(KIND=4),                 INTENT(in)    :: kifill         ! pool value
    INTEGER(KIND=2), DIMENSION(:,:), INTENT(inout) :: kdta           ! mask

    INTEGER :: ik                       ! number of point change
    INTEGER :: ip                       ! size of the pile
    INTEGER :: ji, jj                   ! loop index
    INTEGER :: iip1, iim1, ii, ij       ! working integer
    INTEGER :: ipiglo, ipjglo           ! size of the domain, infered from kdta size

    INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: ipile    ! pile variable
    INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: idata   ! new bathymetry
    !!----------------------------------------------------------------------
    ! infer domain size from input array
    ipiglo = SIZE(kdta,1)
    ipjglo = SIZE(kdta,2)

    ! allocate variable
    ALLOCATE(ipile(2*ipiglo*ipjglo,2))
    ALLOCATE(idata(ipiglo,ipjglo))

    ! initialise variables
    idata=kdta
    ipile(:,:)=0
    ipile(1,:)=[kiseed,kjseed]
    ip=1; ik=0

    ! loop until the pile size is 0 or if the pool is larger than the critical size
    DO WHILE ( ip /= 0 ) ! .AND. ik < 600000);
       ik=ik+1
       ii=ipile(ip,1); ij=ipile(ip,2)
       IF ( MOD(ik, 10000) == 0 ) PRINT *, 'IP =', ip, ik, ii,ij

       ! update bathy and update pile size
       idata(ii,ij) =kifill
       ipile(ip,:)  =[0,0]; ip=ip-1

       ! check neighbour cells and update pile ( assume E-W periodicity )
       iip1=ii+1; IF ( iip1 == ipiglo+1 ) iip1=2
       iim1=ii-1; IF ( iim1 == 0        ) iim1=ipiglo-1

       IF (idata(ii, ij+1) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[ii  ,ij+1]
       END IF
       IF (idata(ii, ij-1) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[ii  ,ij-1]
       END IF
       IF (idata(iip1, ij) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[iip1,ij  ]
       END IF
       IF (idata(iim1, ij) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[iim1,ij  ]
       END IF

    END DO
    kdta=idata;

    DEALLOCATE(ipile); DEALLOCATE(idata)

  END SUBROUTINE FillPool2D


  SUBROUTINE FillPool3D(kiseed, kjseed, kkseed, kdta, kifill)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE FillPool3D  ***
    !!
    !! ** Purpose :  Replace all area surrounding by mask value by mask value
    !!
    !! ** Method  :  flood fill algorithm
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),                   INTENT(in)    :: kiseed, kjseed, kkseed
    INTEGER(KIND=4),                   INTENT(in)    :: kifill   ! new bathymetry
    INTEGER(KIND=2), DIMENSION(:,:,:), INTENT(inout) :: kdta     ! new bathymetry

    INTEGER :: ik,iik                   ! number of point change
    INTEGER :: ip                       ! size of the pile
    INTEGER :: ji, jj                   ! loop index
    INTEGER :: ipiglo, ipjglo, ipk      ! size of the domain inferred from kdta
    INTEGER :: iip1, iim1, ii, ij       ! working integer
    INTEGER :: ijp1, ijm1, ikp1, ikm1
    INTEGER(KIND=2), DIMENSION(:,:),   ALLOCATABLE :: ipile    ! pile variable
    INTEGER(KIND=2), DIMENSION(:,:,:), ALLOCATABLE :: idata   ! new bathymetry
    !!----------------------------------------------------------------------
    ! infer domain size from input array
    ipiglo = SIZE(kdta,1)
    ipjglo = SIZE(kdta,2)
    ipk    = SIZE(kdta,3)

    ! allocate variable
    ALLOCATE(ipile(2*ipiglo*ipjglo*ipk,3))
    ALLOCATE(idata(ipiglo,ipjglo,ipk))

    ! initialise variables
    idata=kdta
    ipile(:,:)=0
    ipile(1,:)=[kiseed,kjseed,kkseed]
    ip=1; iik=0

    ! loop until the pile size is 0 or if the pool is larger than the critical size
    DO WHILE ( ip /= 0 ) !.AND. iik < 600000);
       iik=iik+1
       ii=ipile(ip,1); ij=ipile(ip,2) ; ik=ipile(ip,3)
       IF ( MOD( ip, 1000000) == 0 ) PRINT *, 'IP =', ip, iik, ii,ij, ik

       ! update bathy and update pile size
       idata(ii,ij,ik) = kifill
       ipile(ip,:)  =[0,0,0]; ip=ip-1

       ! check neighbour cells and update pile ( assume E-W periodicity )
       iip1=ii+1; IF ( iip1 == ipiglo+1 ) iip1=2
       iim1=ii-1; IF ( iim1 == 0        ) iim1=ipiglo-1
       ijp1=ij+1 ; ijm1 =ij-1
       ikp1=ik+1 ; IF (ikp1 == ipk+1 ) ikp1=ik
       ikm1=ik-1 ; IF (ikm1 == 0     ) ikm1=ik


       IF (idata(ii, ijp1,ik) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[ii  ,ijp1,ik]
       END IF
       IF (idata(ii, ijm1,ik) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[ii  ,ijm1,ik]
       END IF

       IF (idata(iip1, ij,ik) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[iip1,ij  ,ik]
       END IF
       IF (idata(iim1, ij,ik) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[iim1,ij  ,ik]
       END IF
       IF (idata(ii, ij,ikp1) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[ii  ,ij,ikp1]
       END IF

       IF (idata(ii, ij,ikm1) > 0 ) THEN
          ip=ip+1; ipile(ip,:)=[ii  ,ij,ikm1]
       END IF

    END DO

    kdta=idata;

    DEALLOCATE(ipile); DEALLOCATE(idata)

  END SUBROUTINE FillPool3D

  FUNCTION heading(dplona, dplonb, dplata, dplatb)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION heading  ***
    !!
    !! ** Purpose : Compute true heading between point a and b
    !!
    !! ** Method  : Suppose that the 2 points are not too far away 
    !!              from each other so that heading can be computed 
    !!              with loxodromy.
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=8), INTENT(in) :: dplata, dplona ! lat lon of point a
    REAL(KIND=8), INTENT(in) :: dplatb, dplonb ! lat lon of point b
    REAL(KIND=8)             :: heading        ! return value in degree

    REAL(KIND=8)             :: dlpi, dlconv   ! pi and conversion factor
    REAL(KIND=8)             :: dlxa,dlya      ! working variable
    REAL(KIND=8)             :: dlxb,dlyb      ! working variable
    REAL(KIND=8)             :: dlxb_xa        !  ""        ""
    !!----------------------------------------------------------------------

    dlpi   = ACOS(-1.d0)
    dlconv = dlpi/180.d0  ! for degree to radian conversion

    ! there is a problem if the Greenwich meridian pass between a and b
    dlxa = dplona*dlconv
    dlxb = dplonb*dlconv

    dlya = -LOG(TAN(dlpi/4.-dlconv*dplata/2.d0))
    dlyb = -LOG(TAN(dlpi/4.-dlconv*dplatb/2.d0))

    dlxb_xa = MOD((dlxb-dlxa),2*dlpi)

    IF ( dlxb_xa >=  dlpi ) dlxb_xa = dlxb_xa -2*dlpi
    IF ( dlxb_xa <= -dlpi ) dlxb_xa = dlxb_xa +2*dlpi

    heading=ATAN2(dlxb_xa,(dlyb-dlya))*180.d0/dlpi

    IF (heading < 0) heading=heading+360.d0
  END FUNCTION heading


END MODULE modutils
