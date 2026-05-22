! --- rsa_nastran begin --- !
      SUBROUTINE BD_DTI_SPECSEL ( CARD )

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, JCARD_LEN, JF, IERRFL, FATAL_ERR, RSA_NX_SEMODES
      USE RESPONSE_SPECTRA_STUF, ONLY :  RS_ADD_SPECSEL_ENTRY
      USE MKJCARD_Interface
      USE I4FLD_Interface
      USE R8FLD_Interface

      IMPLICIT NONE

      CHARACTER(LEN=BD_ENTRY_LEN), INTENT(IN) :: CARD
      CHARACTER(LEN=JCARD_LEN)                :: JCARD(10)
      INTEGER(LONG)                           :: TABLED1_ID
      REAL(DOUBLE)                            :: DAMP

! **********************************************************************************************************************************
      JCARD = ' '
      CALL MKJCARD ( 'BD_DTI_SPECSEL', CARD, JCARD )

      IF (JCARD(2)(1:7) /= 'SPECSEL') THEN
         IF (RSA_NX_SEMODES == 'Y') THEN
            WRITE(ERR,1902) JCARD(2)
            WRITE(F06,1902) JCARD(2)
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
         RETURN
      ENDIF

!     Header line like "DTI,SPECSEL,0" is allowed and intentionally ignored.
      IF (JCARD(3)(1:) == '0') RETURN

!     Minimal NX-style support: read (table_id, damping) pairs from fields 6/7 and 8/9.
      TABLED1_ID = 0
      CALL I4FLD ( JCARD(6), JF(6), TABLED1_ID )
      CALL R8FLD ( JCARD(7), JF(7), DAMP )
      IF ((JCARD(6)(1:) /= ' ') .AND. (JCARD(7)(1:) /= ' ') .AND. (IERRFL(6) == 'N') .AND. (IERRFL(7) == 'N')) THEN
         CALL RS_ADD_SPECSEL_ENTRY ( TABLED1_ID, DAMP )
      ENDIF

      TABLED1_ID = 0
      CALL I4FLD ( JCARD(8), JF(8), TABLED1_ID )
      CALL R8FLD ( JCARD(9), JF(9), DAMP )
      IF ((JCARD(8)(1:) /= ' ') .AND. (JCARD(9)(1:) /= ' ') .AND. (IERRFL(8) == 'N') .AND. (IERRFL(9) == 'N')) THEN
         CALL RS_ADD_SPECSEL_ENTRY ( TABLED1_ID, DAMP )
      ENDIF

      RETURN

 1902 FORMAT(' *ERROR  1902: NX-style SOL SEMODES DTI option "',A8,'" is not supported yet. Only DTI,SPECSEL is implemented.')

      END SUBROUTINE BD_DTI_SPECSEL
! --- rsa_nastran end --- !
