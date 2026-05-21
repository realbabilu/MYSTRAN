! --- response_spectra_add begin --- !
! ##################################################################################################################################
      SUBROUTINE BD_DAREA ( CARD )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, JCARD_LEN, JF, IERRFL
      USE RESPONSE_SPECTRA_STUF, ONLY :  RS_ACCUM_DAREA
      USE MKJCARD_Interface
      USE I4FLD_Interface
      USE R8FLD_Interface

      IMPLICIT NONE

      CHARACTER(LEN=BD_ENTRY_LEN), INTENT(IN) :: CARD
      CHARACTER(LEN=JCARD_LEN)                :: JCARD(10)
      INTEGER(LONG)                           :: SID, GRID_ID, COMP, IFIELD
      REAL(DOUBLE)                            :: SCALE

! **********************************************************************************************************************************
      CALL MKJCARD ( 'BD_DAREA', CARD, JCARD )
      SID = 0
      CALL I4FLD ( JCARD(2), JF(2), SID )
      IF (IERRFL(2) /= 'N') RETURN

      DO IFIELD=3,6,3
         GRID_ID = 0
         COMP    = 0
         SCALE   = 0.0D0
         CALL I4FLD ( JCARD(IFIELD  ), JF(IFIELD  ), GRID_ID )
         CALL I4FLD ( JCARD(IFIELD+1), JF(IFIELD+1), COMP    )
         CALL R8FLD ( JCARD(IFIELD+2), JF(IFIELD+2), SCALE   )
         IF ((IERRFL(IFIELD) == 'N') .AND. (IERRFL(IFIELD+1) == 'N') .AND. (IERRFL(IFIELD+2) == 'N')) THEN
            CALL RS_ACCUM_DAREA ( SID, COMP, SCALE )
         ENDIF
      ENDDO

      RETURN

      END SUBROUTINE BD_DAREA
! --- response_spectra_add end --- !
