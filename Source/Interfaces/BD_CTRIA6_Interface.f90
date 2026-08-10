! ##################################################################################################################################
   MODULE BD_CTRIA6_Interface

   INTERFACE

      SUBROUTINE BD_CTRIA6 ( CARD, LARGE_FLD_INP, NUM_GRD )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, IERRFL, FATAL_ERR, JCARD_LEN, JF, LMATANGLE, LPLATEOFF, LPLATETHICK,       &
                                         MEDAT_CTRIA6, NCTRIA6, NEDAT, NELE, NMATANGLE, NPLATEOFF, NPLATETHICK
      USE MODEL_STUF, ONLY            :  EDAT, ETYPE, MATANGLE, PLATEOFF, PLATETHICK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CTRIA6'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP
      INTEGER(LONG), INTENT(OUT)      :: NUM_GRD

      END SUBROUTINE BD_CTRIA6

   END INTERFACE

   END MODULE BD_CTRIA6_Interface
