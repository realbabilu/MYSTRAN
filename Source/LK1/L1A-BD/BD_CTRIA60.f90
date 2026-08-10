! ##################################################################################################################################
      SUBROUTINE BD_CTRIA60 ( CARD, LARGE_FLD_INP )

! Sizing pass for CTRIA6 Bulk Data cards.

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, LMATANGLE, LPLATEOFF, LPLATETHICK

      USE MKJCARD_Interface
      USE NEXTC0_Interface
      USE NEXTC20_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CTRIA60'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP
      CHARACTER(LEN(CARD))            :: CHILD
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)

      INTEGER(LONG)                   :: ICONT = 0
      INTEGER(LONG)                   :: IERR  = 0

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC0  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF

      IF (ICONT == 1) THEN
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (JCARD(2)(1:) /= ' ') LMATANGLE = LMATANGLE + 1
         IF (JCARD(3)(1:) /= ' ') LPLATEOFF = LPLATEOFF + 1
         IF ((JCARD(4)(1:) /= ' ') .OR. (JCARD(5)(1:) /= ' ') .OR. (JCARD(6)(1:) /= ' ')) THEN
            LPLATETHICK = LPLATETHICK + 3
         ENDIF
      ENDIF

      RETURN

      END SUBROUTINE BD_CTRIA60
