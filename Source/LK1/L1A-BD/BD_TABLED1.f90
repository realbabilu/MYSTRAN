      SUBROUTINE BD_TABLED1 ( CARD )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, JCARD_LEN, JF, IERRFL
      USE RESPONSE_SPECTRA_STUF, ONLY :  RS_SET_TABLED1_ID, RS_ADD_TABLED1_POINT
      USE MKJCARD_Interface
      USE I4FLD_Interface
      USE R8FLD_Interface
      USE CHAR_FLD_Interface
      USE LEFT_ADJ_BDFLD_Interface
      USE NEXTC_Interface

      IMPLICIT NONE

      CHARACTER(LEN=BD_ENTRY_LEN), INTENT(INOUT) :: CARD
      CHARACTER(LEN=JCARD_LEN)                    :: JCARD(10)
      CHARACTER(LEN=JCARD_LEN)                    :: TOKEN, TOKEN2
      INTEGER(LONG)                               :: TID, ICONT, IERR, IFLD
      REAL(DOUBLE)                                :: F1, A1, F2, A2

      CALL MKJCARD ( 'BD_TABLED1', CARD, JCARD )
      CALL I4FLD ( JCARD(2), JF(2), TID )
      IF (IERRFL(2) == 'N') CALL RS_SET_TABLED1_ID ( TID )

      DO IFLD=3,9,2
         TOKEN  = JCARD(IFLD)
         TOKEN2 = JCARD(IFLD+1)
         CALL LEFT_ADJ_BDFLD ( TOKEN  )
         CALL LEFT_ADJ_BDFLD ( TOKEN2 )
         IF ((TOKEN(1:4) == 'ENDT') .OR. (TOKEN2(1:4) == 'ENDT')) EXIT
         CALL R8FLD ( JCARD(IFLD  ), JF(IFLD  ), F1 )
         CALL R8FLD ( JCARD(IFLD+1), JF(IFLD+1), A1 )
         IF ((JCARD(IFLD)(1:) /= ' ') .AND. (JCARD(IFLD+1)(1:) /= ' ') .AND. (IERRFL(IFLD) == 'N') .AND.                      &
             (IERRFL(IFLD+1) == 'N')) THEN
            CALL RS_ADD_TABLED1_POINT ( TID, F1, A1 )
         ENDIF
      ENDDO

      DO
         CALL NEXTC ( CARD, ICONT, IERR )
         IF (ICONT /= 1) EXIT
         JCARD = ' '
         CALL MKJCARD ( 'BD_TABLED1', CARD, JCARD )

         TOKEN = ' '
         CALL CHAR_FLD ( JCARD(2), JF(2), TOKEN )
         CALL LEFT_ADJ_BDFLD ( TOKEN )
         IF (TOKEN(1:4) == 'ENDT') EXIT

         DO IFLD=2,8,2
            TOKEN  = JCARD(IFLD)
            TOKEN2 = JCARD(IFLD+1)
            CALL LEFT_ADJ_BDFLD ( TOKEN  )
            CALL LEFT_ADJ_BDFLD ( TOKEN2 )
            IF ((TOKEN(1:4) == 'ENDT') .OR. (TOKEN2(1:4) == 'ENDT')) EXIT
            CALL R8FLD ( JCARD(IFLD  ), JF(IFLD  ), F1 )
            CALL R8FLD ( JCARD(IFLD+1), JF(IFLD+1), A1 )
            IF ((JCARD(IFLD)(1:) /= ' ') .AND. (JCARD(IFLD+1)(1:) /= ' ') .AND. (IERRFL(IFLD) == 'N') .AND.                   &
                (IERRFL(IFLD+1) == 'N')) THEN
               CALL RS_ADD_TABLED1_POINT ( TID, F1, A1 )
            ENDIF
         ENDDO
      ENDDO

      END SUBROUTINE BD_TABLED1
