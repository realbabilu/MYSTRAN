! #################################################################################################################################
! Interface for experimental CQUAD8 SIMOEAS1 branch.

   MODULE CQUAD8_SIMOEAS1_Interface

   INTERFACE

      SUBROUTINE CQUAD8_SIMOEAS1 ( OPT, INT_ELEM_ID )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CQUAD8_SIMOEAS1'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      END SUBROUTINE CQUAD8_SIMOEAS1

   END INTERFACE

   END MODULE CQUAD8_SIMOEAS1_Interface
