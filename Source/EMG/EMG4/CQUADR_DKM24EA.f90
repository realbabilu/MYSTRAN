! #################################################################################################################################
! CQUADR adapter for PARAM,QUADRTYP,DKM24EA.

      SUBROUTINE CQUADR_DKM24EA ( OPT, INT_ELEM_ID )

! DKMQ24 EAS branch based on DKMQ24_EAS4_ShellElement_RHR.py.  The existing
! kernel file still carries the historical DKMQ24N name; this adapter isolates
! that implementation detail from the public PARAM name.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE CQUADR_DKMQ24N_Interface

      IMPLICIT NONE

      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      CALL CQUADR_DKMQ24N ( OPT, INT_ELEM_ID )

      RETURN

      END SUBROUTINE CQUADR_DKM24EA
