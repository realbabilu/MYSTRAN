! --- response_spectra_add begin --- !
! ##################################################################################################################################
      SUBROUTINE BD_TABLED1 ( CARD )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, JCARD_LEN, JF, IERRFL
      USE RESPONSE_SPECTRA_STUF, ONLY :  RS_SET_TABLED1_ID, RS_APPEND_TABLED1_POINT
      USE MKJCARD_Interface
      USE I4FLD_Interface
      USE R8FLD_Interface
      USE NEXTC_Interface

      IMPLICIT NONE

      CHARACTER(LEN=BD_ENTRY_LEN), INTENT(INOUT) :: CARD
      CHARACTER(LEN=JCARD_LEN)                :: JCARD(10)
      INTEGER(LONG)                           :: TID
      INTEGER(LONG)                           :: ICONT, IERR
      REAL(DOUBLE)                            :: F1, A1, F2, A2

! **********************************************************************************************************************************
      CALL MKJCARD ( 'BD_TABLED1', CARD, JCARD )
      CALL I4FLD ( JCARD(2), JF(2), TID )
      IF (IERRFL(2) == 'N') CALL RS_SET_TABLED1_ID ( TID )

      CALL R8FLD ( JCARD(3), JF(3), F1 )
      CALL R8FLD ( JCARD(4), JF(4), A1 )
      IF ((JCARD(3)(1:) /= ' ') .AND. (JCARD(4)(1:) /= ' ') .AND. (IERRFL(3) == 'N') .AND. (IERRFL(4) == 'N')) THEN
         CALL RS_APPEND_TABLED1_POINT ( F1, A1 )
      ENDIF

      CALL R8FLD ( JCARD(5), JF(5), F2 )
      CALL R8FLD ( JCARD(6), JF(6), A2 )
      IF ((JCARD(5)(1:) /= ' ') .AND. (JCARD(6)(1:) /= ' ') .AND. (IERRFL(5) == 'N') .AND. (IERRFL(6) == 'N')) THEN
         CALL RS_APPEND_TABLED1_POINT ( F2, A2 )
      ENDIF

      DO
         CALL NEXTC ( CARD, ICONT, IERR )
         IF (ICONT /= 1) EXIT
         JCARD = ' '
         CALL MKJCARD ( 'BD_TABLED1', CARD, JCARD )

         CALL R8FLD ( JCARD(2), JF(2), F1 )
         CALL R8FLD ( JCARD(3), JF(3), A1 )
         IF ((JCARD(2)(1:) /= ' ') .AND. (JCARD(3)(1:) /= ' ') .AND. (IERRFL(2) == 'N') .AND. (IERRFL(3) == 'N')) THEN
            CALL RS_APPEND_TABLED1_POINT ( F1, A1 )
         ENDIF

         CALL R8FLD ( JCARD(4), JF(4), F2 )
         CALL R8FLD ( JCARD(5), JF(5), A2 )
         IF ((JCARD(4)(1:) /= ' ') .AND. (JCARD(5)(1:) /= ' ') .AND. (IERRFL(4) == 'N') .AND. (IERRFL(5) == 'N')) THEN
            CALL RS_APPEND_TABLED1_POINT ( F2, A2 )
         ENDIF
      ENDDO

      RETURN

      END SUBROUTINE BD_TABLED1
! --- response_spectra_add end --- !
