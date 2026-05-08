! ###############################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________

! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)

! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:

! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________

! End MIT license text.

   MODULE TPLT_MITC3P_Interface

   INTERFACE

      SUBROUTINE TPLT_MITC3P ( OPT, AREA, X2E, X3E, Y3E, CALC_EMATS, IERROR, KV, PTV, PPV, B2V, B3V, S2V, S3V, BIG_BB,            &
                               MN4T_QD, TRIA_NUM, PSI )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NSUB, NTSUB
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  KE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TPLT_MITC3P'
      CHARACTER(1*BYTE), INTENT(IN)   :: CALC_EMATS
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      CHARACTER(LEN=*) , INTENT(IN)   :: MN4T_QD

      INTEGER(LONG), INTENT(OUT)      :: IERROR
      INTEGER(LONG), INTENT(IN)       :: TRIA_NUM

      REAL(DOUBLE) , INTENT(IN)       :: AREA
      REAL(DOUBLE) , INTENT(IN)       :: PSI
      REAL(DOUBLE) , INTENT(IN)       :: X2E
      REAL(DOUBLE) , INTENT(IN)       :: X3E
      REAL(DOUBLE) , INTENT(IN)       :: Y3E
      REAL(DOUBLE) , INTENT(OUT)      :: BIG_BB(3,18,1)
      REAL(DOUBLE) , INTENT(OUT)      :: B2V(3,9)
      REAL(DOUBLE) , INTENT(OUT)      :: B3V(3,9)
      REAL(DOUBLE) , INTENT(OUT)      :: KV(9,9)
      REAL(DOUBLE) , INTENT(OUT)      :: PPV(9,NSUB)
      REAL(DOUBLE) , INTENT(OUT)      :: PTV(9,NTSUB)
      REAL(DOUBLE) , INTENT(OUT)      :: S2V(3,9)
      REAL(DOUBLE) , INTENT(OUT)      :: S3V(3,9)

      END SUBROUTINE TPLT_MITC3P

   END INTERFACE

   END MODULE TPLT_MITC3P_Interface
