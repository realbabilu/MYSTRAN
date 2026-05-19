! BD_RBE2 field contract:
!   first line:        fields 5-9  -> dependent grids
!   continuation line: fields 2-9  -> dependent grids
!
!   FIELD   ITEM
!   -----   ------------
!    2      RELID, Rigid Elem ID
!    3      IGID, Independent Grid ID
!    4      DDOF, Dependent DOF's for the Grids following
!   5-9     DGID, Dependent Grid ID's
!
! on optional continuation cards:
!   2-9     DGID, Dependent grid ID's

      DO J=5,9
         IF (JCARD(J)(1:) == ' ') THEN
            CYCLE
         ELSE
            CALL I4FLD ( JCARD(J), JF(J), DGID )
            IF ((IERRFL(J) == 'N') .AND. (JERR == 0)) THEN
               WRITE(L1F) RTYPE
               WRITE(L1F) RELID,DGID,DDOF,IGID
            ENDIF
         ENDIF
      ENDDO

      DO
         CALL NEXTC  ( CARD, ICONT, IERR )
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN
            DO J=2,9
               IF (JCARD(J)(1:) == ' ') THEN
                  CYCLE
               ELSE
                  CALL I4FLD ( JCARD(J), JF(J), DGID )
                  IF ((IERRFL(J) == 'N') .AND. (JERR == 0)) THEN
                     WRITE(L1F) RTYPE
                     WRITE(L1F) RELID,DGID,DDOF,IGID
                  ENDIF
               ENDIF
            ENDDO
         ELSE
            EXIT
         ENDIF
      ENDDO
