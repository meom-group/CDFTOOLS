  IF ( narg == 0 ) THEN
     PRINT *,' usage : 
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ', TRIM(cv_out),' (    )'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_uvwtfil ) ; ijarg=ijarg+1
        ! option
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out     ) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

