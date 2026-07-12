! ##################################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________
!
! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
! LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
! EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
! USE OR OTHER DEALINGS IN THE SOFTWARE.
! _______________________________________________________________________________________________________
!
! End MIT license text.

      SUBROUTINE MATMULX_FFF_T ( A, B, NROWA, NCOLA, NCOLB, C )

! BLAS-first dense matrix multiply for hotspot callers.
! Multiplies A' x B and returns C. Caller is responsible for conformity.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE

      USE MATMULX_FFF_T_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATMULX_FFF_T'
      CHARACTER( 1*BYTE), PARAMETER   :: TRANSA = 'T'
      CHARACTER( 1*BYTE), PARAMETER   :: TRANSB = 'N'
      CHARACTER( 1*BYTE), PARAMETER   :: TRANS  = 'T'

      INTEGER(LONG), INTENT(IN)       :: NROWA
      INTEGER(LONG), INTENT(IN)       :: NCOLA
      INTEGER(LONG), INTENT(IN)       :: NCOLB

      REAL(DOUBLE), INTENT(IN)        :: A(NROWA,NCOLA)
      REAL(DOUBLE), INTENT(IN)        :: B(NROWA,NCOLB)
      REAL(DOUBLE), INTENT(OUT)       :: C(NCOLA,NCOLB)
      REAL(DOUBLE), PARAMETER         :: ALPHA = ONE
      REAL(DOUBLE), PARAMETER         :: BETA  = ZERO

      EXTERNAL                        :: DGEMM, DGEMV

! **********************************************************************************************************************************

      IF (NCOLB == 1) THEN
         CALL DGEMV ( TRANSA, NROWA, NCOLA, ALPHA, A, NROWA, B(1,1), 1, BETA, C(1,1), 1 )

      ELSE IF (NCOLA == 1) THEN
         CALL DGEMV ( TRANS, NROWA, NCOLB, ALPHA, B, NROWA, A(1,1), 1, BETA, C(1,1), 1 )

      ELSE
         CALL DGEMM ( TRANSA, TRANSB, NCOLA, NCOLB, NROWA, ALPHA, A, NROWA, B, NROWA, BETA, C, NCOLA )

      ENDIF

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE MATMULX_FFF_T
