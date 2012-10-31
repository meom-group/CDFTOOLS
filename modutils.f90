MODULE modutils
  !!======================================================================
  !!                     ***  MODULE  modutils  ***
  !! Hold functions and subroutine dedicated to common utility task
  !!=====================================================================
  !! History : 3.0  : 04/2011  : J.M. Molines : Original code
  !!                : 10/2012  : N. Ferry, E. Durand, F. Hernandez : add shapiro
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   SetGlobalAtt  : Set Global Attribute to the command line
  !!   SetFilename   : Build standard name from confname
  !!   shapiro_fill_smooth : shapiro smoother or filler
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  USE cdfio

  IMPLICIT NONE

  PRIVATE
  PUBLIC SetGlobalAtt
  PUBLIC SetFileName
  PUBLIC shapiro_fill_smooth

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

  CHARACTER(LEN=256) FUNCTION SetFileName(cdconf, cdtag, cdgrid )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION SetFileName  ***
    !!
    !! ** Purpose :  Build filename from cdconf, tag and grid
    !!
    !! ** Method  :  Check 2 forms of file names and return
    !!               error is file is missing
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdconf, cdtag, cdgrid
    !!----------------------------------------------------------------------
    WRITE( SetFileName,'(a,"_",a,"_grid",a,".nc")') TRIM(cdconf), TRIM(cdtag), TRIM(cdgrid)
    IF ( chkfile(SetFileName ) ) THEN ! look for another name
       WRITE(SetFileName,'(a,"_",a,"_grid_",a,".nc")') TRIM(cdconf), TRIM(cdtag), TRIM(cdgrid)
       IF ( chkfile( SetFileName)  ) THEN
          PRINT *,' ERROR : missing grid',TRIM(cdgrid),'or even grid_',TRIM(cdgrid),' file '
          STOP
       ENDIF
    ENDIF
  END FUNCTION SetFileName

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
    IF ( cdfs(1:4)  .NE. 'fill' .AND. cdfs(1:6) .NE. 'smooth' ) THEN
       PRINT*, 'cdfs = ',cdfs ,' <> fill or smooth'
       STOP
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

END MODULE modutils
