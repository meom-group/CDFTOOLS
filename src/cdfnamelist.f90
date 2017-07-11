PROGRAM cdfnamelist
  !!======================================================================
  !!                     ***  PROGRAM  cdfnamelist  ***
  !!=====================================================================
  !!  ** Purpose : Give informations on the namelist mechanism implemented
  !!               in CDFTOOLS_3. 
  !! CDFTOOLS_4.0 , MEOM 2017 
  !!               to change default file names, variable  or dimension
  !!               names.
  !!
  !!  ** Method  : 
  !!
  !! History :  3.0  : 01/2011  : J.M. Molines : Original code
  !!         :  4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class file_informations
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)   :: narg, iargc, ijarg

  CHARACTER(LEN=80) :: cldum
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfnamelist [-i] [-p]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Give information [-i option] on the namelist mechanism implemented ' 
     PRINT *,'       since CDFTOOLS v3. Write a namelist template [-p option] to initialize'
     PRINT *,'       the mechanism.'
     PRINT *,'      '
     PRINT *,'       Note also than since CDFTOOLS v4, the name of the mesh/mask files as'
     PRINT *,'       well as the name of the coordinates file can be modified from the '
     PRINT *,'       default, setting corresponding environment variable :'
     PRINT *,'           CDFT_MESH_HGR   to change the default ',TRIM(cn_fhgr)
     PRINT *,'           CDFT_MESH_ZGR   to change the default ',TRIM(cn_fzgr)
     PRINT *,'           CDFT_MASK       to change the default ',TRIM(cn_fmsk)
     PRINT *,'           CDFT_BASINS     to change the default ',TRIM(cn_fbasins)
     PRINT *,'           CDFT_COORD      to change the default ',TRIM(cn_fcoo)
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [ -i ] : print informations ' 
     PRINT *,'       [ -p ] : write a template namelist.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       with option -p, print a template namelist : PrintCdfNames.namlist'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1
  DO WHILE ( ijarg <= narg ) 
     CALL getarg( ijarg, cldum ) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-i' ) 
        CALL InfoUseNamelist()
     CASE ( '-p' )
        CALL PrintCdfNames()
     CASE DEFAULT
        PRINT *, TRIM(cldum),' : unknown option in cdfnamelist '
     END SELECT
  END DO

CONTAINS

  SUBROUTINE InfoUseNamelist()
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE InfoUseNamelist  ***
    !!
    !! ** Purpose :  Print detailed info on the use of namelist in
  !! CDFTOOLS_4.0 , MEOM 2017 
    !!
    !!----------------------------------------------------------------------
    PRINT *,'   In CDFTOOLS_3 the variable names, dimension names, mesh_mask'
    PRINT *,' file names can be customized via a system of namelist.'
    PRINT *,'   A call to ReadCdfNames at the begining of the program allows'
    PRINT *,' the update of the names used in the program.'
    PRINT *,'   If there is no need for changing names, then it is not necessary'
    PRINT *,' to give a namelist, the default values are OK.'
    PRINT *,'   '
    PRINT *,'   If you need to change any of the default values, then you can'
    PRINT *,' use the namelist system to make this change effective. Doing do'
    PRINT *,' some rules are to be followed for proper use.'
    PRINT *,'   '
    PRINT *,'NAMELIST EDITING'
    PRINT *,'   To have a template of a CDFTOOLS_3 namelist, use the statement'
    PRINT *,'     cdfnamelist -p '
    PRINT *,' This will give you a template namelist (PrintCdfNames.namlist)'
    PRINT *,' that you have to customized for your application.'
    PRINT *,'   Some comments are made within this namelist for particular blocks.'
    PRINT *,' '
    PRINT *,'NAME AND LOCATION OF THE NAMELIST'
    PRINT *,'   The default name of the namelist read by ReadCdfNames is '
    PRINT *,'     nam_cdf_names'
    PRINT *,'   ReadCdfNames look for the namelist in the current directory (./)'
    PRINT *,' and, if not found there, in the $HOME/CDFTOOLS_cfg/ directory'
    PRINT *,'   The name of the namelist can be changed setting the environment'
    PRINT *,' variable NAM_CDF_NAMES to the desired value.' 
    PRINT *,'    '

  END SUBROUTINE InfoUseNamelist

END PROGRAM cdfnamelist
