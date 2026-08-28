! #################################################################################################################################
! CQUADR adapter for PARAM,QUADRTYP,DKM24AU.

      SUBROUTINE CQUADR_DKM24AU ( OPT, INT_ELEM_ID )

! Legacy AU DKMQ24 branch. Kept as a separate entry point so the QUADRTYP
! dispatcher maps one public element name to one replaceable source file.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE CQUADR_DKMQ24_Interface
      USE CQUADR_DKMQ24R_Interface

      IMPLICIT NONE

      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      CHARACTER(1*BYTE)               :: REC_OPT(6)

      CALL CQUADR_DKMQ24 ( OPT, INT_ELEM_ID )

      IF (OPT(3) == 'Y') THEN
         REC_OPT = 'N'
         REC_OPT(3) = 'Y'
         CALL CQUADR_DKMQ24R ( REC_OPT, INT_ELEM_ID )
      ENDIF

      RETURN

      END SUBROUTINE CQUADR_DKM24AU
