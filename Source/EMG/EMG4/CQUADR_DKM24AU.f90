! #################################################################################################################################
! CQUADR adapter for PARAM,QUADRTYP,DKM24AU.

      SUBROUTINE CQUADR_DKM24AU ( OPT, INT_ELEM_ID )

! Legacy AU DKMQ24 branch. Kept as a separate entry point so existing decks
! using QUADRTYP=DKM24AU still map to the historical CQUADR_DKMQ24R kernel
! while QUADRTYP=DKMQ24 can point at the plain CQUADR_DKMQ24 implementation.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE CQUADR_DKMQ24_Interface
      USE CQUADR_DKMQ24R_Interface

      IMPLICIT NONE

      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      CALL CQUADR_DKMQ24R ( OPT, INT_ELEM_ID )

      RETURN

      END SUBROUTINE CQUADR_DKM24AU
