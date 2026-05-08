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

      SUBROUTINE TPLT_MITC3P ( OPT, AREA, X2E, X3E, Y3E, BIG_BB )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE MODEL_STUF, ONLY            :  ELDOF

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TPLT_MITC3P'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)

      REAL(DOUBLE), INTENT(IN)        :: AREA
      REAL(DOUBLE), INTENT(IN)        :: X2E
      REAL(DOUBLE), INTENT(IN)        :: X3E
      REAL(DOUBLE), INTENT(IN)        :: Y3E
      REAL(DOUBLE), INTENT(OUT)       :: BIG_BB(3,ELDOF,1)

      END SUBROUTINE TPLT_MITC3P

   END INTERFACE

   END MODULE TPLT_MITC3P_Interface
