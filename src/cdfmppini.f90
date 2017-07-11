PROGRAM cdfmppini 
  !!======================================================================
  !!                     ***  PROGRAM  cdfmppini  ***
  !!=====================================================================
  !!  ** Purpose : off line domain decomposition using mesh_hgr
  !!
  !!  ** Method  : just an incapsulation of mpp_ini from NEMO
  !!
  !! History : 2.1  : 05/2010  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   mpp_init2     Nemo routine for mpp initialisation
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class preprocessing
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  ! REM : some of the doctor rules are not followed because we want to use
  !       the mpp_init2 routine almost out of the box, thus we need to define
  !       variables which are parameters in NEMO, with the same name (jp..)

  INTEGER(KIND=4), PARAMETER                   :: wp=8   ! working precision
  INTEGER(KIND=4)                              :: jpni, jpnj, jpnij
  INTEGER(KIND=4)                              :: jpreci=1 , jprecj=1
  INTEGER(KIND=4)                              :: jpi, jpj, jpiglo, jpjglo
  INTEGER(KIND=4)                              :: jperio=6, jv
  INTEGER(KIND=4)                              :: narg, iargc, numout=6
  INTEGER(KIND=4)                              :: ijarg
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: imask
  INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nimppt, njmppt, nlcit, nlcjt
  INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nldit, nldjt, nleit, nlejt
  INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: nbondi, nbondj, icount


  CHARACTER(LEN=80)                            :: cf_msk='m'
  CHARACTER(LEN=80)                            :: cf_out='mppini.txt'
  CHARACTER(LEN=80)                            :: cv_in
  CHARACTER(LEN=80)                            :: cldum

  LOGICAL                                      :: lwp=.TRUE.
  !----------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmppini -i jpni -j jpnj [-m/-b/-z] [-jperio jperio]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Performs the mpp initialisation with NEMO routine mpp_init2 and gives' 
     PRINT *,'       some statistics about the domains. Save the layout on a text file.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -i jpni : number of domains in the i direction.' 
     PRINT *,'       -j jpnj : number of domains in the j direction.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-m/-b/-z] : use one of these option to choose the land/sea mask.'
     PRINT *,'               m  : take mask from ',TRIM(cn_fmsk),' (tmask) [ default ]'
     PRINT *,'               b  : take mask from ',TRIM(cn_fbathymet),' (Bathymetry)'
     PRINT *,'               z  : take mask from ',TRIM(cn_fzgr),' (mbathy)'
     PRINT *,'                   Default is ',TRIM(cf_msk)
     PRINT *,'       [-jperio jperio ] : specify jperio. '
     PRINT '(a,i2)','                          default value is ', jperio
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       one of ',TRIM(cn_fmsk),', ',TRIM(cn_fbathymet),' or ',TRIM(cn_fzgr),' according to option' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       - Standard output'
     PRINT *,'       - ASCII file ', TRIM(cf_out)
     STOP 
  ENDIF

  cf_msk = cn_fmsk ; cv_in=cn_tmask
  ijarg=1 
  DO WHILE ( ijarg <= narg )
     CALL getarg( ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-i'     ) ; CALL getarg( ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) jpni
     CASE ( '-j'     ) ; CALL getarg( ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) jpnj
        ! options
     CASE( '-m'      ) ; cf_msk = cn_fmsk       ; cv_in = cn_tmask
     CASE( '-b'      ) ; cf_msk = cn_fbathymet  ; cv_in = cn_bathymet
     CASE( '-z'      ) ; cf_msk = cn_fzgr       ; cv_in = cn_mbathy
     CASE( '-jperio' ) ; CALL getarg( ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) jperio
     CASE DEFAULT      ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  IF ( chkfile (cf_msk ))  STOP 99 ! missing file

  jpiglo = getdim (cf_msk,cn_x)
  jpjglo = getdim (cf_msk,cn_y)

  jpi = ( jpiglo-2*jpreci + (jpni-1) ) / jpni + 2*jpreci 
  jpj = ( jpjglo-2*jprecj + (jpnj-1) ) / jpnj + 2*jprecj 

  ALLOCATE ( imask(jpiglo,jpjglo) )
  imask(:,:) = getvar(cf_msk, cv_in, 1, jpiglo, jpjglo)
  WHERE (imask <= 0 ) imask = 0
  WHERE (imask > 0  ) imask = 1
  CALL mpp_init2
  PRINT *, 'JPIGLO= ', jpiglo
  PRINT *, 'JPJGLO= ', jpjglo
  PRINT *, 'JPI   = ', jpi
  PRINT *, 'JPJ   = ', jpj
  PRINT *, 'JPNI  = ', jpni
  PRINT *, 'JPNJ  = ', jpnj
  PRINT *, 'JPNIJ = ', jpnij

  PRINT *, 'NBONDI between : ',MINVAL(nbondi),' AND ', MAXVAL(nbondi)
  PRINT *, 'NBONDJ between : ',MINVAL(nbondj),' AND ', MAXVAL(nbondj)
  PRINT *,' Accounting ...'
  ALLOCATE (icount(jpnij))
  DO jv=-1,2
     icount=0
     WHERE(nbondi == jv ) icount=1
     PRINT *,' NBONDI = ', jv,' : ', SUM(icount)
  ENDDO

  DO jv=-1,2
     icount=0
     WHERE(nbondj == jv ) icount=1
     PRINT *,' NBONDJ = ', jv,' : ', SUM(icount)
  ENDDO


CONTAINS

  SUBROUTINE mpp_init2
    !!----------------------------------------------------------------------
    !!                  ***  ROUTINE mpp_init2  ***
    !!
    !! * Purpose :   Lay out the global domain over processors.
    !!     FOR USING THIS VERSION, A PREPROCESSING TRAITMENT IS RECOMMENDED
    !!     FOR DEFINING BETTER CUTTING OUT.
    !!       This routine is used with a the bathymetry file.
    !!       In this version, the land processors are avoided and the adress
    !!     processor (nproc, narea,noea, ...) are calculated again.
    !!     The jpnij parameter can be lesser than jpni x jpnj
    !!     and this jpnij parameter must be calculated before with an
    !!     algoritmic preprocessing program.
    !!
    !! ** Method  :   Global domain is distributed in smaller local domains.
    !!      Periodic condition is a function of the local domain position
    !!      (global boundary or neighbouring domain) and of the global
    !!      periodic
    !!      Type :         jperio global periodic condition
    !!                     nperio local  periodic condition
    !!
    !! ** Action :        nimpp     : longitudinal index 
    !!                    njmpp     : latitudinal  index
    !!                    nperio    : lateral condition type 
    !!                    narea     : number for local area
    !!                    nlci      : first dimension
    !!                    nlcj      : second dimension
    !!                    nproc     : number for local processor
    !!                    noea      : number for local neighboring processor
    !!                    nowe      : number for local neighboring processor
    !!                    noso      : number for local neighboring processor
    !!                    nono      : number for local neighboring processor
    !!
    !! History :
    !!         :  4.0  : 03/2017  : J.M. Molines  
    !!        !  94-11  (M. Guyon)  Original code
    !!        !  95-04  (J. Escobar, M. Imbard)
    !!        !  98-02  (M. Guyon)  FETI method
    !!        !  98-05  (M. Imbard, J. Escobar, L. Colombet )  SHMEM and MPI versions
    !!   9.0  !  04-01  (G. Madec, J.M Molines)  F90 : free form , north fold jpni > 1
    !!----------------------------------------------------------------------
    !! 
    INTEGER(KIND=4) :: ji, jj, jn, jproc, jarea     ! dummy loop indices
    INTEGER(KIND=4) ::  inum = 99                   ! temporary logical unit
    INTEGER(KIND=4) ::   &
         ii, ij, ifreq, il1, il2,          &  ! temporary integers
         icont, ili, ilj,                  &  !    "          "
         isurf, ijm1, imil,                &  !    "          "
         iino, ijno, iiso, ijso,           &  !    "          " 
         iiea, ijea, iiwe, ijwe,           &  !    "          "
         iresti, irestj, iproc                !    "          "
    INTEGER(KIND=4) :: nreci, nrecj,  nperio
    INTEGER(KIND=4), DIMENSION(10000)          ::    iint, ijnt          
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE ::    iin, ijn          
    INTEGER(KIND=4), DIMENSION(jpni,jpnj) ::   &
         iimppt, ijmppt, ilci  , ilcj  ,   &  ! temporary workspace
         ipproc, ibondj, ibondi,           &  !    "           "
         ilei  , ilej  , ildi  , ildj  ,   &  !    "           "
         ioea  , iowe  , ioso  , iono         !    "           "
    REAL(wp) ::   zidom , zjdom          ! temporary scalars

    INTEGER(KIND=4) :: nono, noso, noea, nowe
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ii_nono, ii_noso, ii_noea, ii_nowe

    ! 0. initialisation
    ! -----------------

    !  1. Dimension arrays for subdomains
    ! -----------------------------------

    !  Computation of local domain sizes ilci() ilcj()
    !  These dimensions depend on global sizes jpni,jpnj and jpiglo,jpjglo
    !  The subdomains are squares leeser than or equal to the global
    !  dimensions divided by the number of processors minus the overlap
    !  array.

    nreci=2*jpreci
    nrecj=2*jprecj
    iresti = 1 + MOD( jpiglo - nreci -1 , jpni )
    irestj = 1 + MOD( jpjglo - nrecj -1 , jpnj )

    ilci(1:iresti      ,:) = jpi
    ilci(iresti+1:jpni ,:) = jpi-1

    ilcj(:,      1:irestj) = jpj
    ilcj(:, irestj+1:jpnj) = jpj-1

    IF(lwp) WRITE(numout,*)
    IF(lwp) WRITE(numout,*) ' mpp_init2: defines mpp subdomains'
    IF(lwp) WRITE(numout,*) ' ~~~~~~  ----------------------'
    IF(lwp) WRITE(numout,*)
    IF(lwp) WRITE(numout,*) 'iresti=',iresti,' irestj=',irestj
    IF(lwp) WRITE(numout,*)
    IF(lwp) WRITE(numout,*) 'jpni=',jpni,' jpnj=',jpnj

    zidom = nreci + SUM(ilci(:,1) - nreci ) 
    IF(lwp) WRITE(numout,*)
    IF(lwp) WRITE(numout,*)' sum ilci(i,1)=',zidom,' jpiglo=',jpiglo

    zjdom = nrecj + SUM(ilcj(1,:) - nrecj ) 
    IF(lwp) WRITE(numout,*) ' sum ilcj(1,j)=',zjdom,' jpjglo=',jpjglo
    IF(lwp) WRITE(numout,*)


    !  2. Index arrays for subdomains
    ! -------------------------------

    iimppt(:,:) = 1
    ijmppt(:,:) = 1
    ipproc(:,:) = -1

    IF( jpni > 1 )THEN
       DO jj = 1, jpnj
          DO ji = 2, jpni
             iimppt(ji,jj) = iimppt(ji-1,jj) + ilci(ji-1,jj) - nreci
          END DO
       END DO
    ENDIF

    IF( jpnj > 1 )THEN
       DO jj = 2, jpnj
          DO ji = 1, jpni
             ijmppt(ji,jj) = ijmppt(ji,jj-1) + ilcj(ji,jj-1) - nrecj
          END DO
       END DO
    ENDIF


    ! 3. Subdomain description in the Regular Case
    ! --------------------------------------------

    nperio = 0
    icont = -1
    DO jarea = 1, jpni*jpnj
       ii = 1 + MOD(jarea-1,jpni)
       ij = 1 +    (jarea-1)/jpni
       ili = ilci(ii,ij)
       ilj = ilcj(ii,ij)

       ibondj(ii,ij) = -1
       IF( jarea >  jpni          )   ibondj(ii,ij) = 0
       IF( jarea >  (jpnj-1)*jpni )   ibondj(ii,ij) = 1
       IF( jpnj  == 1             )   ibondj(ii,ij) = 2

       ibondi(ii,ij) = 0
       IF( MOD(jarea,jpni) == 1 )   ibondi(ii,ij) = -1
       IF( MOD(jarea,jpni) == 0 )   ibondi(ii,ij) =  1
       IF( jpni            == 1 )   ibondi(ii,ij) =  2

       ! 2.4 Subdomain neighbors

       iproc = jarea - 1
       ioso(ii,ij) = iproc - jpni
       iowe(ii,ij) = iproc - 1
       ioea(ii,ij) = iproc + 1
       iono(ii,ij) = iproc + jpni

       ildi(ii,ij) = 1 + jpreci
       ilei(ii,ij) = ili -jpreci

       IF( ibondi(ii,ij) == -1 .OR. ibondi(ii,ij) == 2 ) ildi(ii,ij) = 1
       IF( ibondi(ii,ij) ==  1 .OR. ibondi(ii,ij) == 2 ) ilei(ii,ij) = ili

       ildj(ii,ij) =  1  + jprecj
       ilej(ii,ij) = ilj - jprecj
       IF( ibondj(ii,ij) == -1 .OR. ibondj(ii,ij) == 2 ) ildj(ii,ij) = 1
       IF( ibondj(ii,ij) ==  1 .OR. ibondj(ii,ij) == 2 ) ilej(ii,ij) = ilj

       ! warning ii*ij (zone) /= nproc (processors)!

       IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 ) THEN
          IF( jpni == 1 )THEN
             ibondi(ii,ij) = 2
             nperio = 1
          ELSE
             ibondi(ii,ij) = 0
          ENDIF
          IF( MOD(jarea,jpni) == 0 ) THEN
             ioea(ii,ij) = iproc - (jpni-1)
          ENDIF
          IF( MOD(jarea,jpni) == 1 ) THEN
             iowe(ii,ij) = iproc + jpni - 1
          ENDIF
       ENDIF

       isurf = 0
       DO jj = 1+jprecj, ilj-jprecj
          DO  ji = 1+jpreci, ili-jpreci
             IF( imask(ji+iimppt(ii,ij)-1, jj+ijmppt(ii,ij)-1) == 1) isurf = isurf+1
          END DO
       END DO
       IF(isurf /= 0) THEN
          icont = icont + 1
          ipproc(ii,ij) = icont
          iint(icont+1) = ii
          ijnt(icont+1) = ij
       ENDIF
    END DO
    jpnij=icont+1
    ALLOCATE(iin(jpnij),ijn(jpnij),nimppt(jpnij), njmppt(jpnij), nlcit(jpnij), nlcjt(jpnij)  )
    ALLOCATE(nldit(jpnij), nldjt(jpnij)  )
    ALLOCATE(nleit(jpnij), nlejt(jpnij)  )
    ALLOCATE(nbondi(jpnij), nbondj(jpnij)  )
    ALLOCATE(ii_nono(jpnij), ii_noso(jpnij), ii_noea(jpnij) , ii_nowe(jpnij) )

    iin(:)=iint(1:jpnij)
    ijn(:)=ijnt(1:jpnij)

    ! Control
    ! 4. Subdomain print
    ! ------------------

    IF(lwp) THEN
       ifreq = 4
       il1 = 1
       DO jn = 1,(jpni-1)/ifreq+1
          il2 = MIN(jpni,il1+ifreq-1)
          WRITE(numout,*)
          WRITE(numout,9400) ('***',ji=il1,il2-1)
          DO jj = jpnj, 1, -1
             WRITE(numout,9403) ('   ',ji=il1,il2-1)
             WRITE(numout,9402) jj, (ilci(ji,jj),ilcj(ji,jj),ji=il1,il2)
             WRITE(numout,9404) (ipproc(ji,jj),ji=il1,il2)
             WRITE(numout,9403) ('   ',ji=il1,il2-1)
             WRITE(numout,9400) ('***',ji=il1,il2-1)
          END DO
          WRITE(numout,9401) (ji,ji=il1,il2)
          il1 = il1+ifreq
       END DO
9400   FORMAT('     ***',20('*************',a3))
9403   FORMAT('     *     ',20('         *   ',a3))
9401   FORMAT('        ',20('   ',i3,'          '))
9402   FORMAT(' ',i3,' *  ',20(i3,'  x',i3,'   *   '))
9404   FORMAT('     *  ',20('     ',i4,'   *   '))
    ENDIF


    ! 5. neighbour treatment
    ! ----------------------

    DO jarea = 1, jpni*jpnj
       iproc = jarea-1
       ii = 1 + MOD(jarea-1,jpni)
       ij = 1 +    (jarea-1)/jpni
       IF( ipproc(ii,ij) == -1 .AND. iono(ii,ij) >= 0   &
            .AND. iono(ii,ij) <= jpni*jpnj-1 ) THEN
          iino = 1 + MOD(iono(ii,ij),jpni)
          ijno = 1 +    (iono(ii,ij))/jpni
          IF( ibondj(iino,ijno) == 1 ) ibondj(iino,ijno)=2
          IF( ibondj(iino,ijno) == 0 ) ibondj(iino,ijno) = -1
       ENDIF
       IF( ipproc(ii,ij) == -1 .AND. ioso(ii,ij) >= 0   &
            .AND. ioso(ii,ij) <= jpni*jpnj-1 ) THEN
          iiso = 1 + MOD(ioso(ii,ij),jpni)
          ijso = 1 +    (ioso(ii,ij))/jpni
          IF( ibondj(iiso,ijso) == -1 ) ibondj(iiso,ijso) = 2
          IF( ibondj(iiso,ijso) ==  0 ) ibondj(iiso,ijso) = 1
       ENDIF
       IF( ipproc(ii,ij) == -1 .AND. ioea(ii,ij) >= 0   &
            .AND. ioea(ii,ij) <= jpni*jpnj-1) THEN
          iiea = 1 + MOD(ioea(ii,ij),jpni)
          ijea = 1 +    (ioea(ii,ij))/jpni
          IF( ibondi(iiea,ijea) == 1 ) ibondi(iiea,ijea) = 2
          IF( ibondi(iiea,ijea) == 0 ) ibondi(iiea,ijea) = -1
       ENDIF
       IF( ipproc(ii,ij) == -1 .AND. iowe(ii,ij) >= 0   &
            .AND. iowe(ii,ij) <= jpni*jpnj-1) THEN
          iiwe = 1 + MOD(iowe(ii,ij),jpni)
          ijwe = 1 +    (iowe(ii,ij))/jpni
          IF( ibondi(iiwe,ijwe) == -1 ) ibondi(iiwe,ijwe) = 2
          IF( ibondi(iiwe,ijwe) ==  0 ) ibondi(iiwe,ijwe) = 1
       ENDIF
    END DO


    ! just to save nono etc for all proc
    DO jarea = 1, jpnij
       ii = iin(jarea)
       ij = ijn(jarea)
       IF( ioso(ii,ij) >= 0 .AND. ioso(ii,ij) <= (jpni*jpnj-1) ) THEN
          iiso = 1 + MOD(ioso(ii,ij),jpni)
          ijso = 1 +    (ioso(ii,ij))/jpni
          noso = ipproc(iiso,ijso)
          ii_noso(jarea)= noso
       ENDIF
       IF( iowe(ii,ij) >= 0 .AND. iowe(ii,ij) <= (jpni*jpnj-1) ) THEN
          iiwe = 1 + MOD(iowe(ii,ij),jpni)
          ijwe = 1 +    (iowe(ii,ij))/jpni
          nowe = ipproc(iiwe,ijwe)
          ii_nowe(jarea)= nowe
       ENDIF
       IF( ioea(ii,ij) >= 0 .AND. ioea(ii,ij) <= (jpni*jpnj-1) ) THEN
          iiea = 1 + MOD(ioea(ii,ij),jpni)
          ijea = 1 +    (ioea(ii,ij))/jpni
          noea = ipproc(iiea,ijea)
          ii_noea(jarea)= noea
       ENDIF
       IF( iono(ii,ij) >= 0 .AND. iono(ii,ij) <= (jpni*jpnj-1) ) THEN
          iino = 1 + MOD(iono(ii,ij),jpni)
          ijno = 1 +    (iono(ii,ij))/jpni
          nono = ipproc(iino,ijno)
          ii_nono(jarea)= nono
       ENDIF
    END DO
    ! 6. Change processor name
    ! ------------------------

    DO jproc = 1, jpnij
       ii = iin(jproc)
       ij = ijn(jproc)

       nimppt(jproc) = iimppt(ii,ij)  
       njmppt(jproc) = ijmppt(ii,ij)  

       nlcit(jproc) = ilci(ii,ij)
       nlcjt(jproc) = ilcj(ii,ij)

       nldit(jproc) = ildi(ii,ij)
       nldjt(jproc) = ildj(ii,ij)

       nleit(jproc) = ilei(ii,ij)
       nlejt(jproc) = ilej(ii,ij)
    END DO

    ! Save processor layout in ascii file
    IF (lwp) THEN
       OPEN (inum, FILE=cf_out, FORM='FORMATTED', RECL=255)
       WRITE(inum,'(6i8)') jpnij,jpi,jpj,jpiglo,jpjglo
       WRITE(inum,'(a)') 'NAREA   nlci nlcj nldi nldj nlei nlej nimpp njmpp nono noso nowe noea nbondi nbondj '

       DO  jproc = 1, jpnij
          ii = iin(jproc)
          ij = ijn(jproc)
          nbondi(jproc) = ibondi(ii,ij)
          nbondj(jproc) = ibondj(ii,ij)


          WRITE(inum,'(15i5)') jproc  , nlcit(jproc), nlcjt(jproc), &
               nldit(jproc), nldjt(jproc), &
               nleit(jproc), nlejt(jproc), &
               nimppt(jproc), njmppt(jproc),& 
               ii_nono(jproc), ii_noso(jproc), ii_nowe(jproc), ii_noea(jproc) ,&
               nbondi(jproc),  nbondj(jproc) 
       END DO
       CLOSE(inum)   
    END IF


  END SUBROUTINE mpp_init2

END PROGRAM cdfmppini
