MODULE modutils
  !!======================================================================
  !!                     ***  MODULE  modutils  ***
  !! Hold functions and subroutine dedicated to common utility task
  !!=====================================================================
  !! History : 3.0  : 04/2011  : J.M. Molines : Original code
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   SetGlobalAtt  : Set Global Attribute to the command line
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


END MODULE modutils
