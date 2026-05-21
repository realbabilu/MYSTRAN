! --- response_spectra_add begin --- !
! ##################################################################################################################################
      SUBROUTINE BD_FREQ1 ( CARD )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, JCARD_LEN, JF, IERRFL
      USE RESPONSE_SPECTRA_STUF, ONLY :  RS_FREQ1_SID, RS_NUM_FREQ, RS_FREQS, MAX_RS_FREQ
      USE MKJCARD_Interface
      USE I4FLD_Interface
      USE R8FLD_Interface

      IMPLICIT NONE

      CHARACTER(LEN=BD_ENTRY_LEN), INTENT(IN) :: CARD
      CHARACTER(LEN=JCARD_LEN)                :: JCARD(10)
      INTEGER(LONG)                           :: SID, NPTS, I
      REAL(DOUBLE)                            :: F1, DF

! **********************************************************************************************************************************
      CALL MKJCARD ( 'BD_FREQ1', CARD, JCARD )
      CALL I4FLD ( JCARD(2), JF(2), SID )
      CALL R8FLD ( JCARD(3), JF(3), F1 )
      CALL R8FLD ( JCARD(4), JF(4), DF )
      CALL I4FLD ( JCARD(5), JF(5), NPTS )

      IF ((IERRFL(2) == 'N') .AND. (IERRFL(3) == 'N') .AND. (IERRFL(4) == 'N') .AND. (IERRFL(5) == 'N')) THEN
         RS_FREQ1_SID = SID
         DO I=1,NPTS
            IF (RS_NUM_FREQ < MAX_RS_FREQ) THEN
               RS_NUM_FREQ = RS_NUM_FREQ + 1
               RS_FREQS(RS_NUM_FREQ) = F1 + DF*DBLE(I-1)
            ENDIF
         ENDDO
      ENDIF

      RETURN

      END SUBROUTINE BD_FREQ1
! --- response_spectra_add end --- !
