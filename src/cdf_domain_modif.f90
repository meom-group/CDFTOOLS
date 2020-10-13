PROGRAM cdf_domain_modif
  !!======================================================================
  !!                     ***  PROGRAM  cdf_domain_modif  ***
  !!=====================================================================
  !!  ** Purpose : adjust scalar variable in (jpiglo, jpjglo, ... jperio)
  !!               in an existing domain_cfg file, for instance obtained
  !!               when extracting a sub area with ncks
  !!
  !!  ** Method  :  READ and Write ... 
  !!
  !! History :  4.0  : 03/2017  : J.M. Molines : 
  !!----------------------------------------------------------------------
  USE netcdf

  IMPLICIT NONE
  INTEGER(KIND=4) :: jperio, jpiglo, jpjglo, jpkglo
  INTEGER(KIND=4) :: jpisfcav, jpsco, jpzps, jpzco
  INTEGER(KIND=4) :: narg, ijarg, ncid, id, ierr

  LOGICAL         :: ll_jperio = .FALSE.
  LOGICAL         :: ll_jpiglo = .FALSE.
  LOGICAL         :: ll_jpjglo = .FALSE.
  LOGICAL         :: ll_jpkglo = .FALSE.
  LOGICAL         :: ll_isfcav = .FALSE.
  LOGICAL         :: ll_sco    = .FALSE.
  LOGICAL         :: ll_zps    = .FALSE.
  LOGICAL         :: ll_zco    = .FALSE.


  CHARACTER(LEN=80) :: cf_domain , cldum

  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2020
  !! $Id$
  !! Copyright (c) 2020, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class domain_file
  !!----------------------------------------------------------------------
  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdf_domain_modif -d DOMAIN_cfg [-jperio JPERIO] [-jpiglo JPIGLO] ...'
     PRINT *,'                          [-jpjglo JPJGLO] [-jpkglo JPKGLO   ...]'
     PRINT *,'                          [-ln_isfcav LN_ISFCAV] [-ln_sco LN_SCO] ...'
     PRINT *,'                          [-ln_zps LN_ZPS ] [-ln_zco LN_ZCO]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Adjust scalar variable in (jpiglo, jpjglo, jpkglo, jperio) and boolean'
     PRINT *,'       variables (ln_ifscav, ln_sco, ln_zps, ln_zco) in an existing domain_cfg'
     PRINT *,'       file, for instance obtained  when extracting a sub area with ncks.' 

     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -d DOMAIN_cfg : specify the domain_cfg file to work with' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -jperio JPERIO : pass the new value for jperio. If option not used,'
     PRINT *,'                        no changes.' 
     PRINT *,'       -jpiglo JPIGLO : pass the new value for jpiglo. If option not used,'
     PRINT *,'                        no changes.' 
     PRINT *,'       -jpjglo JPJGLO : pass the new value for jpjglo. If option not used,'
     PRINT *,'                        no changes.' 
     PRINT *,'       -jpkglo JPKGLO : pass the new value for jpkglo. If option not used,'
     PRINT *,'                        no changes.' 
     PRINT *,'       -ln_isfcav LN_ISFCAV : pass the new value for ln_isfcav ( 1 or 0 ).'
     PRINT *,'                 If option not used, no changes.'
     PRINT *,'       -ln_sco LN_SCO : pass the new value for ln_sco ( 1 or 0 ).'
     PRINT *,'                 If option not used, no changes.'
     PRINT *,'       -ln_zps LN_ZPS : pass the new value for ln_zps ( 1 or 0 ).'
     PRINT *,'                 If option not used, no changes.'
     PRINT *,'       -ln_zco LN_ZCO : pass the new value for ln_zco ( 1 or 0 ).'
     PRINT *,'                 If option not used, no changes.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : Input file is modified on the output.'
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-d'   ) ; CALL getarg(ijarg, cf_domain ) ; ijarg=ijarg+1
        ! option
     CASE ( '-jperio'    ) ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) jperio ; ll_jperio=.TRUE.
     CASE ( '-jpiglo'    ) ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) jpiglo ; ll_jpiglo=.TRUE.
     CASE ( '-jpjglo'    ) ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) jpjglo ; ll_jpjglo=.TRUE.
     CASE ( '-jpkglo'    ) ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) jpkglo ; ll_jpkglo=.TRUE.

     CASE ( '-ln_isfcav' ) ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) jpisfcav ; ll_isfcav=.TRUE. 
     CASE ( '-ln_sco'    ) ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) jpsco    ; ll_sco=.TRUE.
     CASE ( '-ln_zps'    ) ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) jpzps    ; ll_zps=.TRUE.
     CASE ( '-ln_zco'    ) ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) jpzco    ; ll_zco=.TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  ierr = NF90_OPEN(cf_domain, NF90_WRITE,ncid)

  IF ( ll_jperio) THEN
     ierr = NF90_INQ_VARID(ncid,'jperio',id )
     ierr = NF90_PUT_VAR(ncid,id,jperio)
  ENDIF
  IF ( ll_jpiglo) THEN
     ierr = NF90_INQ_VARID(ncid,'jpiglo',id )
     ierr = NF90_PUT_VAR(ncid,id,jpiglo)
  ENDIF
  IF ( ll_jpjglo) THEN
     ierr = NF90_INQ_VARID(ncid,'jpjglo',id )
     ierr = NF90_PUT_VAR(ncid,id,jpjglo)
  ENDIF
  IF ( ll_jpkglo) THEN
     ierr = NF90_INQ_VARID(ncid,'jpkglo',id )
     ierr = NF90_PUT_VAR(ncid,id,jpkglo)
  ENDIF
  IF ( ll_isfcav) THEN
     ierr = NF90_INQ_VARID(ncid,'ln_isfcav',id )
     ierr = NF90_PUT_VAR(ncid,id,jpisfcav)
  ENDIF
  IF ( ll_sco) THEN
     ierr = NF90_INQ_VARID(ncid,'ln_sco',id )
     ierr = NF90_PUT_VAR(ncid,id,jpsco)
  ENDIF
  IF ( ll_zps) THEN
     ierr = NF90_INQ_VARID(ncid,'ln_zps',id )
     ierr = NF90_PUT_VAR(ncid,id,jpzps)
  ENDIF
  IF ( ll_zco) THEN
     ierr = NF90_INQ_VARID(ncid,'ln_zco',id )
     ierr = NF90_PUT_VAR(ncid,id,jpzco)
  ENDIF

  ierr = NF90_CLOSE(ncid) 


END PROGRAM cdf_domain_modif
