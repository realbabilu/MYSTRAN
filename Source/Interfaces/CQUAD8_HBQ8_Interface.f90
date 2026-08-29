! #################################################################################################################################
! Interface for CQUAD8 HBQ8 branch.

   MODULE CQUAD8_HBQ8_Interface

   INTERFACE

      SUBROUTINE CQUAD8_HBQ8 ( OPT, INT_ELEM_ID )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CQUAD8_HBQ8'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      END SUBROUTINE CQUAD8_HBQ8

   END INTERFACE

   END MODULE CQUAD8_HBQ8_Interface
