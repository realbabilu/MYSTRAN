      SUBROUTINE BD_PYRAM0 ( CARD, LARGE_FLD_INP, DELTA_LEDAT )

! Counts EDAT storage for PYRAM5/PYRAM14. No continuation means PYRAM5; any
! nonblank continuation promotes the entry to the 14-node pyramid.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, MEDAT_PYRAM5, MEDAT_PYRAM14
      USE TIMDAT, ONLY                :  TSEC

      USE BD_PYRAM0_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PYRAM0'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP
      CHARACTER(LEN(CARD))            :: CHILD

      INTEGER(LONG)                   :: ICONT = 0
      INTEGER(LONG)                   :: IERR  = 0
      INTEGER(LONG), INTENT(OUT)      :: DELTA_LEDAT

      SUBR_NAME = SUBR_NAME

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC0  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF

      IF (ICONT == 1) THEN
         DELTA_LEDAT = MEDAT_PYRAM14
      ELSE
         DELTA_LEDAT = MEDAT_PYRAM5
      ENDIF

      RETURN

      END SUBROUTINE BD_PYRAM0
