PROGRAM cdfpolymask
  !!-------------------------------------------------------------------
  !!              ***  PROGRAM CDFPOLYMASK   ***
  !!
  !!  **  Purpose: Create a nc file with 1 into subareas defined as a polygone
  !!  
  !!  **  Method:  Use polylib routine (from finite element mesh generator Trigrid)
  !!               Read vertices of polygone in an ascii file an produce a resulting
  !!               file the same shape as file givent in argumment (used only for size and
  !!               header )
  !!
  !! history:
  !!    Original:  J.M. Molines (July 2007 )
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !!
  !! * Modules used
  USE cdfio 

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk
  INTEGER   :: narg, iargc                                  !: 
  INTEGER   :: npiglo,npjglo, npk                                !: size of the domain
  INTEGER, DIMENSION(1) ::  ipk, id_varout
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: rpmask
  REAL(KIND=4) ,DIMENSION(1)                  :: timean

  CHARACTER(LEN=80) :: cfile, cpoly, cfileout='polymask.nc'            !: file name
  TYPE(variable), DIMENSION(1) :: typvar

  INTEGER    :: ncout
  INTEGER    :: istatus, ierr

  !!  Read command line
  narg= iargc()
  IF ( narg /= 2 ) THEN
     PRINT *,' Usage : cdfpolymask ''polygon file'' ''reference ncfile'' '
     PRINT *,'   polygons are defined on the I,J grid'
     PRINT *,'   Output on polymask.nc ,variable polymask'
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cpoly)
  CALL getarg (2, cfile)

  npiglo = getdim (cfile,'x')
  npjglo = getdim (cfile,'y')
  npk = 1

  ipk(1)      = 1
  typvar(1)%name='polymask'
  typvar(1)%units='1/0'
  typvar(1)%missing_value=999.
  typvar(1)%valid_min= 0.
  typvar(1)%valid_max= 1.
  typvar(1)%long_name='Polymask'
  typvar(1)%short_name='polymask'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo

  ALLOCATE( rpmask(npiglo,npjglo) )

  ncout =create(cfileout, cfile,npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfile, npiglo, npjglo,npk)

  CALL polymask(cpoly, rpmask) 

  ierr=putvar(ncout,id_varout(1), rpmask, 1 ,npiglo, npjglo)
  timean=getvar1d(cfile,'time_counter',1)
  ierr=putvar1d(ncout,timean,1,'T')
  istatus = closeout(ncout)

CONTAINS
  SUBROUTINE polymask( cdpoly, pmask)
    USE modpoly
    REAL(KIND=4), DIMENSION(:,:), INTENT(out) :: pmask
    CHARACTER(LEN=*), INTENT(in) :: cdpoly 

    ! *Local variables
    INTEGER :: ji,jj, nfront, jjpoly
    LOGICAL :: l_in
    CHARACTER(LEN=80), DIMENSION(jpolys) :: carea
    pmask=0.
    CALL ReadPoly(cdpoly,nfront, carea)
    DO jjpoly=1, nfront
       CALL PrepPoly(jjpoly)
       DO jj=npjglo, 1, -1
          DO ji=1,npiglo
             CALL InPoly(jjpoly,float(ji), float(jj), l_in)
             IF (l_in) pmask(ji,jj)=1.
          ENDDO
!      IF ( jj < 405 .AND. jj > 335 ) THEN
!         print '(i4,100i2)', jj, NINT(pmask(170:260,jj))
!      ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE polymask

END PROGRAM cdfpolymask
