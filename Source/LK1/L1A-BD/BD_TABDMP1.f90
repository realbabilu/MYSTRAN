      SUBROUTINE BD_TABDMP1 ( CARD )

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, JCARD_LEN, JF, IERRFL
      USE RESPONSE_SPECTRA_STUF, ONLY :  RS_SET_TABDMP1_DAMP
      USE MKJCARD_Interface
      USE I4FLD_Interface
      USE R8FLD_Interface
      USE NEXTC_Interface

      IMPLICIT NONE

      CHARACTER(LEN=BD_ENTRY_LEN), INTENT(INOUT) :: CARD
      CHARACTER(LEN=JCARD_LEN)                   :: JCARD(10)
      CHARACTER(LEN=JCARD_LEN)                   :: FIELD_FREQ, FIELD_DAMP
      INTEGER(LONG)                              :: SID, ICONT, IERR, J
      REAL(DOUBLE)                               :: FREQ, DAMP
      LOGICAL                                    :: HAVE_DAMP

      CALL MKJCARD ( 'BD_TABDMP1', CARD, JCARD )
      CALL I4FLD ( JCARD(2), JF(2), SID )
      IF ((IERRFL(2) /= 'N') .OR. (SID <= 0)) RETURN

      HAVE_DAMP = .FALSE.

      DO J=4,8,2
         FREQ = 0.0D0
         DAMP = -1.0D0
         FIELD_FREQ = JCARD(J)
         FIELD_DAMP = JCARD(J+1)
         IF (INDEX(FIELD_FREQ,'ENDT') > 0) FIELD_FREQ = FIELD_FREQ(1:INDEX(FIELD_FREQ,'ENDT')-1)
         IF (INDEX(FIELD_DAMP,'ENDT') > 0) FIELD_DAMP = FIELD_DAMP(1:INDEX(FIELD_DAMP,'ENDT')-1)
         IF ((LEN_TRIM(FIELD_FREQ) == 0) .OR. (LEN_TRIM(FIELD_DAMP) == 0)) CYCLE
         CALL R8FLD ( FIELD_FREQ, JF(J),   FREQ )
         CALL R8FLD ( FIELD_DAMP, JF(J+1), DAMP )
         IF ((IERRFL(J) == 'N') .AND. (IERRFL(J+1) == 'N')) THEN
            CALL RS_SET_TABDMP1_DAMP ( SID, FREQ, DAMP )
            HAVE_DAMP = .TRUE.
         ENDIF
      ENDDO

      DO
         CALL NEXTC ( CARD, ICONT, IERR )
         IF (ICONT /= 1) EXIT
         JCARD = ' '
         CALL MKJCARD ( 'BD_TABDMP1', CARD, JCARD )
         DO J=2,8,2
            FREQ = 0.0D0
            DAMP = -1.0D0
            FIELD_FREQ = JCARD(J)
            FIELD_DAMP = JCARD(J+1)
            IF (INDEX(FIELD_FREQ,'ENDT') > 0) FIELD_FREQ = FIELD_FREQ(1:INDEX(FIELD_FREQ,'ENDT')-1)
            IF (INDEX(FIELD_DAMP,'ENDT') > 0) FIELD_DAMP = FIELD_DAMP(1:INDEX(FIELD_DAMP,'ENDT')-1)
            IF ((LEN_TRIM(FIELD_FREQ) == 0) .OR. (LEN_TRIM(FIELD_DAMP) == 0)) CYCLE
            CALL R8FLD ( FIELD_FREQ, JF(J),   FREQ )
            CALL R8FLD ( FIELD_DAMP, JF(J+1), DAMP )
            IF ((IERRFL(J) == 'N') .AND. (IERRFL(J+1) == 'N')) THEN
               CALL RS_SET_TABDMP1_DAMP ( SID, FREQ, DAMP )
               HAVE_DAMP = .TRUE.
            ENDIF
         ENDDO
      ENDDO

      END SUBROUTINE BD_TABDMP1
