! --- response_spectra_add begin --- !
! ##################################################################################################################################
      SUBROUTINE BD_RLOAD1 ( CARD )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, JCARD_LEN, JF, IERRFL
      USE RESPONSE_SPECTRA_STUF, ONLY :  RS_SET_RLOAD1
      USE MKJCARD_Interface
      USE I4FLD_Interface

      IMPLICIT NONE

      CHARACTER(LEN=BD_ENTRY_LEN), INTENT(IN) :: CARD
      CHARACTER(LEN=JCARD_LEN)                :: JCARD(10)
      INTEGER(LONG)                           :: SID, EXCITE_SID, TABLED1_ID

! **********************************************************************************************************************************
      CALL MKJCARD ( 'BD_RLOAD1', CARD, JCARD )
      CALL I4FLD ( JCARD(2), JF(2), SID )
      EXCITE_SID = 0
      TABLED1_ID = 0
      CALL I4FLD ( JCARD(3), JF(3), EXCITE_SID )
      CALL I4FLD ( JCARD(6), JF(6), TABLED1_ID )
      IF (IERRFL(2) == 'N') THEN
         CALL RS_SET_RLOAD1 ( SID, EXCITE_SID, TABLED1_ID )
      ENDIF

      RETURN

      END SUBROUTINE BD_RLOAD1
! --- response_spectra_add end --- !
