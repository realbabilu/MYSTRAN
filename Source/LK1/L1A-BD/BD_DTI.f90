      SUBROUTINE BD_DTI ( CARD )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, JCARD_LEN, JF, IERRFL
      USE RESPONSE_SPECTRA_STUF, ONLY :  RS_ADD_SPECSEL_PAIR
      USE MKJCARD_Interface
      USE I4FLD_Interface
      USE R8FLD_Interface
      USE CHAR_FLD_Interface
      USE LEFT_ADJ_BDFLD_Interface
      USE NEXTC_Interface

      IMPLICIT NONE

      CHARACTER(LEN=BD_ENTRY_LEN), INTENT(INOUT) :: CARD
      CHARACTER(LEN=JCARD_LEN)                   :: JCARD(10)
      CHARACTER(16)                              :: DTI_NAME, KIND
      INTEGER(LONG)                              :: LINE_ID, TID, ICONT, IERR
      REAL(DOUBLE)                               :: DAMP

      CALL MKJCARD ( 'BD_DTI', CARD, JCARD )
      DTI_NAME = ' '
      CALL CHAR_FLD ( JCARD(2), JF(2), DTI_NAME )
      CALL LEFT_ADJ_BDFLD ( DTI_NAME )

      IF (DTI_NAME(1:7) /= 'SPECSEL') RETURN

      CALL I4FLD ( JCARD(3), JF(3), LINE_ID )
      IF ((IERRFL(3) /= 'N') .OR. (LINE_ID <= 0)) RETURN

      KIND = 'A       '
      IF (JCARD(5)(1:) /= ' ') THEN
         CALL CHAR_FLD ( JCARD(5), JF(5), KIND )
         CALL LEFT_ADJ_BDFLD ( KIND )
      ENDIF

      TID  = 0
      DAMP = 0.0D0
      CALL I4FLD ( JCARD(6), JF(6), TID )
      CALL R8FLD ( JCARD(7), JF(7), DAMP )
      IF ((IERRFL(6) == 'N') .AND. (IERRFL(7) == 'N')) CALL RS_ADD_SPECSEL_PAIR ( LINE_ID, KIND, TID, DAMP )

      TID  = 0
      DAMP = 0.0D0
      CALL I4FLD ( JCARD(8), JF(8), TID )
      CALL R8FLD ( JCARD(9), JF(9), DAMP )
      IF ((IERRFL(8) == 'N') .AND. (IERRFL(9) == 'N')) CALL RS_ADD_SPECSEL_PAIR ( LINE_ID, KIND, TID, DAMP )

      DO
         CALL NEXTC ( CARD, ICONT, IERR )
         IF (ICONT /= 1) EXIT
         JCARD = ' '
         CALL MKJCARD ( 'BD_DTI', CARD, JCARD )

         TID  = 0
         DAMP = 0.0D0
         CALL I4FLD ( JCARD(2), JF(2), TID )
         CALL R8FLD ( JCARD(3), JF(3), DAMP )
         IF ((IERRFL(2) == 'N') .AND. (IERRFL(3) == 'N')) CALL RS_ADD_SPECSEL_PAIR ( LINE_ID, KIND, TID, DAMP )

         TID  = 0
         DAMP = 0.0D0
         CALL I4FLD ( JCARD(4), JF(4), TID )
         CALL R8FLD ( JCARD(5), JF(5), DAMP )
         IF ((IERRFL(4) == 'N') .AND. (IERRFL(5) == 'N')) CALL RS_ADD_SPECSEL_PAIR ( LINE_ID, KIND, TID, DAMP )
      ENDDO

      END SUBROUTINE BD_DTI
