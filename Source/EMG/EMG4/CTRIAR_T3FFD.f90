! ##################################################################################################################################
! CTRIAR T3FFD selector for PARAM,TRIARTYP,T3FFD.

      SUBROUTINE CTRIAR_T3FFD ( OPT, INT_ELEM_ID )

! Port target: D:\18a\bending_only\Shell\gemini2\shit\T3FFA_ShellElement_Std.py
! T3FFA Std default stiffness is the same SNORM-aware three-stage T3FF core:
! K_E -> K_A, fictitious drilling in A-basis, then K_G.  The shared core detects
! TRIARTYP='T3FFD' for the T3FFA-specific shear recovery convention.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE CTRIA3_T3FF_Interface

      IMPLICIT NONE

      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      CALL CTRIA3_T3FF ( OPT, INT_ELEM_ID )

      RETURN

      END SUBROUTINE CTRIAR_T3FFD
