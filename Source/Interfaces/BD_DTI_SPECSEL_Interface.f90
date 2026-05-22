! --- rsa_nastran begin --- !
   MODULE BD_DTI_SPECSEL_Interface

      INTERFACE

      SUBROUTINE BD_DTI_SPECSEL ( CARD )

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE SCONTR, ONLY                :  BD_ENTRY_LEN

      IMPLICIT NONE

      CHARACTER(LEN=BD_ENTRY_LEN), INTENT(IN) :: CARD

      END SUBROUTINE BD_DTI_SPECSEL

      END INTERFACE

   END MODULE BD_DTI_SPECSEL_Interface
! --- rsa_nastran end --- !
