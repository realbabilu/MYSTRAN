! #################################################################################################################################
   MODULE CTRIA6_REZAIEE_Interface

   INTERFACE

      SUBROUTINE CTRIA6_REZAIEE ( OPT, INT_ELEM_ID )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CTRIA6_REZAIEE'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      END SUBROUTINE CTRIA6_REZAIEE

   END INTERFACE

   END MODULE CTRIA6_REZAIEE_Interface
