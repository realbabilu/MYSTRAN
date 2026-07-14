      SUBROUTINE BD_DLOAD ( CARD )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, JCARD_LEN, JF, IERRFL
      USE RESPONSE_SPECTRA_STUF, ONLY :  RS_DLOAD_SID, RS_ADD_DLOAD_TERM
      USE MKJCARD_Interface
      USE I4FLD_Interface
      USE R8FLD_Interface
      USE NEXTC_Interface

      IMPLICIT NONE

      CHARACTER(LEN=BD_ENTRY_LEN), INTENT(INOUT) :: CARD
      CHARACTER(LEN=JCARD_LEN)                :: JCARD(10)
      INTEGER(LONG)                           :: SID, RLOAD1_SID, IFIELD, ICONT, IERR
      REAL(DOUBLE)                            :: S0, SI

      CALL MKJCARD ( 'BD_DLOAD', CARD, JCARD )
      CALL I4FLD ( JCARD(2), JF(2), SID )
      IF (IERRFL(2) == 'N') THEN
         RS_DLOAD_SID = SID
         S0 = 1.0D0
         CALL R8FLD ( JCARD(3), JF(3), S0 )
         DO IFIELD=4,8,2
            SI = 1.0D0
            RLOAD1_SID = 0
            CALL R8FLD ( JCARD(IFIELD),   JF(IFIELD),   SI )
            CALL I4FLD ( JCARD(IFIELD+1), JF(IFIELD+1), RLOAD1_SID )
            IF (RLOAD1_SID > 0) CALL RS_ADD_DLOAD_TERM ( SID, RLOAD1_SID, S0*SI )
         ENDDO

         DO
            CALL NEXTC ( CARD, ICONT, IERR )
            IF (ICONT /= 1) EXIT
            JCARD = ' '
            CALL MKJCARD ( 'BD_DLOAD', CARD, JCARD )
            DO IFIELD=2,8,2
               SI = 1.0D0
               RLOAD1_SID = 0
               CALL R8FLD ( JCARD(IFIELD),   JF(IFIELD),   SI )
               CALL I4FLD ( JCARD(IFIELD+1), JF(IFIELD+1), RLOAD1_SID )
               IF (RLOAD1_SID > 0) CALL RS_ADD_DLOAD_TERM ( SID, RLOAD1_SID, S0*SI )
            ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE BD_DLOAD
