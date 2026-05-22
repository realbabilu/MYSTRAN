! --- rsa_nastran begin --- !
! ##################################################################################################################################
      SUBROUTINE BD_TABDMP1 ( CARD )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, JCARD_LEN, JF, IERRFL
      USE RESPONSE_SPECTRA_STUF, ONLY :  RS_SET_TABDMP1_ID, RS_NUM_DAMP, RS_DAMP_FREQ, RS_DAMP_VAL, MAX_RS_DAMP
      USE MKJCARD_Interface
      USE I4FLD_Interface
      USE NEXTC_Interface
      USE R8FLD_Interface

      IMPLICIT NONE

      CHARACTER(LEN=BD_ENTRY_LEN), INTENT(INOUT) :: CARD
      CHARACTER(LEN=JCARD_LEN)                   :: JCARD(10)
      CHARACTER(LEN=JCARD_LEN)                   :: WORK
      INTEGER(LONG)                              :: TID
      INTEGER(LONG)                              :: ICONT, IERR
      REAL(DOUBLE)                               :: F1, A1, F2, A2

! **********************************************************************************************************************************
      CALL MKJCARD ( 'BD_TABDMP1', CARD, JCARD )
      CALL I4FLD ( JCARD(2), JF(2), TID )
      IF (IERRFL(2) == 'N') CALL RS_SET_TABDMP1_ID ( TID )

      CALL CLEAN_ENDT_FIELD ( JCARD(4), WORK )
      CALL R8FLD ( WORK, JF(4), F1 )
      CALL CLEAN_ENDT_FIELD ( JCARD(5), WORK )
      CALL R8FLD ( WORK, JF(5), A1 )
      IF ((JCARD(4)(1:) /= ' ') .AND. (JCARD(5)(1:) /= ' ') .AND. (IERRFL(4) == 'N') .AND. (IERRFL(5) == 'N')) THEN
         CALL STORE_DAMP_POINT ( F1, A1 )
      ENDIF

      CALL CLEAN_ENDT_FIELD ( JCARD(6), WORK )
      CALL R8FLD ( WORK, JF(6), F2 )
      CALL CLEAN_ENDT_FIELD ( JCARD(7), WORK )
      CALL R8FLD ( WORK, JF(7), A2 )
      IF ((JCARD(6)(1:) /= ' ') .AND. (JCARD(7)(1:) /= ' ') .AND. (IERRFL(6) == 'N') .AND. (IERRFL(7) == 'N')) THEN
         CALL STORE_DAMP_POINT ( F2, A2 )
      ENDIF

      DO
         CALL NEXTC ( CARD, ICONT, IERR )
         IF (ICONT /= 1) EXIT
         JCARD = ' '
         CALL MKJCARD ( 'BD_TABDMP1', CARD, JCARD )

         CALL CLEAN_ENDT_FIELD ( JCARD(2), WORK )
         CALL R8FLD ( WORK, JF(2), F1 )
         CALL CLEAN_ENDT_FIELD ( JCARD(3), WORK )
         CALL R8FLD ( WORK, JF(3), A1 )
         IF ((JCARD(2)(1:) /= ' ') .AND. (JCARD(3)(1:) /= ' ') .AND. (IERRFL(2) == 'N') .AND. (IERRFL(3) == 'N')) THEN
            CALL STORE_DAMP_POINT ( F1, A1 )
         ENDIF

         CALL CLEAN_ENDT_FIELD ( JCARD(4), WORK )
         CALL R8FLD ( WORK, JF(4), F2 )
         CALL CLEAN_ENDT_FIELD ( JCARD(5), WORK )
         CALL R8FLD ( WORK, JF(5), A2 )
         IF ((JCARD(4)(1:) /= ' ') .AND. (JCARD(5)(1:) /= ' ') .AND. (IERRFL(4) == 'N') .AND. (IERRFL(5) == 'N')) THEN
            CALL STORE_DAMP_POINT ( F2, A2 )
         ENDIF
      ENDDO

      RETURN

      CONTAINS

      SUBROUTINE STORE_DAMP_POINT ( FREQ, DAMP )
      REAL(DOUBLE), INTENT(IN) :: FREQ, DAMP
      IF (RS_NUM_DAMP < MAX_RS_DAMP) THEN
         RS_NUM_DAMP = RS_NUM_DAMP + 1
         RS_DAMP_FREQ(RS_NUM_DAMP) = FREQ
         RS_DAMP_VAL (RS_NUM_DAMP) = DAMP
      ENDIF
      END SUBROUTINE STORE_DAMP_POINT

      SUBROUTINE CLEAN_ENDT_FIELD ( FIELD_IN, FIELD_OUT )
      CHARACTER(LEN=*), INTENT(IN)  :: FIELD_IN
      CHARACTER(LEN=*), INTENT(OUT) :: FIELD_OUT
      INTEGER(LONG)                 :: IPOS

      FIELD_OUT = FIELD_IN
      IPOS = INDEX(FIELD_OUT,'ENDT')
      IF (IPOS > 0) THEN
         FIELD_OUT(IPOS:) = ' '
      ENDIF
      END SUBROUTINE CLEAN_ENDT_FIELD

      END SUBROUTINE BD_TABDMP1
! --- rsa_nastran end --- !
