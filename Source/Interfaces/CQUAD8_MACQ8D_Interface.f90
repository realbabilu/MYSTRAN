! #################################################################################################################################
! Interface for CQUAD8 MacNeal Q8 + drilling penalty branch.

   MODULE CQUAD8_MACQ8D_Interface

   INTERFACE

      SUBROUTINE CQUAD8_MACQ8D ( OPT, INT_ELEM_ID )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CQUAD8_MACQ8D'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      END SUBROUTINE CQUAD8_MACQ8D

   END INTERFACE

   END MODULE CQUAD8_MACQ8D_Interface
