! ##################################################################################################################################
   MODULE BD_CTRIA60_Interface

   INTERFACE

      SUBROUTINE BD_CTRIA60 ( CARD, LARGE_FLD_INP )

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, LMATANGLE, LPLATEOFF, LPLATETHICK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CTRIA60'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP
      CHARACTER(LEN(CARD))            :: CHILD
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)

      INTEGER(LONG)                   :: ICONT
      INTEGER(LONG)                   :: IERR

      END SUBROUTINE BD_CTRIA60

   END INTERFACE

   END MODULE BD_CTRIA60_Interface
