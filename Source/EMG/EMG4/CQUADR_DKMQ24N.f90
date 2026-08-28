! #################################################################################################################################
! Backup compatibility wrapper for the historical DKMQ24N kernel name.

      SUBROUTINE CQUADR_DKMQ24N ( OPT, INT_ELEM_ID )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE CQUADR_DKM24EA_Interface

      IMPLICIT NONE

      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      CALL CQUADR_DKM24EA ( OPT, INT_ELEM_ID )

      RETURN

      END SUBROUTINE CQUADR_DKMQ24N
