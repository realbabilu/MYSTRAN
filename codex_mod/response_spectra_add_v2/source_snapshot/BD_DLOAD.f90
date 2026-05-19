! --- response_spectra_add begin --- !
! ##################################################################################################################################
      SUBROUTINE BD_DLOAD ( CARD )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, JCARD_LEN, JF, IERRFL
      USE RESPONSE_SPECTRA_STUF, ONLY :  RS_DLOAD_SID, RS_ADD_DLOAD_TERM
      USE MKJCARD_Interface
      USE I4FLD_Interface

      IMPLICIT NONE

      CHARACTER(LEN=BD_ENTRY_LEN), INTENT(IN) :: CARD
      CHARACTER(LEN=JCARD_LEN)                :: JCARD(10)
      INTEGER(LONG)                           :: SID, RLOAD1_SID, IFIELD

! **********************************************************************************************************************************
      CALL MKJCARD ( 'BD_DLOAD', CARD, JCARD )
      CALL I4FLD ( JCARD(2), JF(2), SID )
      IF (IERRFL(2) == 'N') THEN
         RS_DLOAD_SID = SID
         DO IFIELD=5,9,2
            RLOAD1_SID = 0
            CALL I4FLD ( JCARD(IFIELD), JF(IFIELD), RLOAD1_SID )
            IF (RLOAD1_SID > 0) THEN
               CALL RS_ADD_DLOAD_TERM ( SID, RLOAD1_SID )
            ENDIF
         ENDDO
      ENDIF

      RETURN

      END SUBROUTINE BD_DLOAD
! --- response_spectra_add end --- !
