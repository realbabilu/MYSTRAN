! ##################################################################################################################################
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

      SUBROUTINE MATMULT_FFF_T ( A, B, NROWA, NCOLA, NCOLB, C )

! Multiplies two matrices: A' x B (A is transposed). Returns result, matrix C. All matrices are in full format
! NOTE: User is responsible for making sure that A(t) and B are conformable

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE

      USE MATMULT_FFF_T_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATMULT_FFF_T'
      CHARACTER( 1*BYTE), PARAMETER   :: TRANSA = 'T'
      CHARACTER( 1*BYTE), PARAMETER   :: TRANSB = 'N'
      CHARACTER( 1*BYTE), PARAMETER   :: TRANS  = 'T'

      INTEGER(LONG), INTENT(IN)       :: NROWA             ! No. rows in input matrix A (NOT A')
      INTEGER(LONG), INTENT(IN)       :: NCOLA             ! No. cols in input matrix A (NOT A')
      INTEGER(LONG), INTENT(IN)       :: NCOLB             ! No. cols in input matrix B
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: MIN_DIM
      INTEGER(LONG)                   :: NROWB             ! No. rows in input matrix B
      INTEGER(LONG)                   :: NROWA_T           ! No. rows in A' (ranspose)
      INTEGER(LONG)                   :: NCOLA_T           ! No. cols in A' (ranspose)
      REAL(DOUBLE), PARAMETER         :: DGEMM_MIN_WORK = 32768.D0
      INTEGER(LONG), PARAMETER        :: DGEMM_MIN_DIM  = 24
      INTEGER(LONG), PARAMETER        :: DGEMV_MIN_LEN  = 64


      REAL(DOUBLE) , INTENT(IN)       :: A(NROWA,NCOLA)    ! Input  matrix A
      REAL(DOUBLE) , INTENT(IN)       :: B(NROWA,NCOLB)    ! Input  matrix B
      REAL(DOUBLE) , INTENT(OUT)      :: C(NCOLA,NCOLB)    ! Output matrix C
      REAL(DOUBLE) , PARAMETER        :: ALPHA = ONE
      REAL(DOUBLE) , PARAMETER        :: BETA  = ZERO
      REAL(DOUBLE)                    :: WORK_EST

      EXTERNAL                        :: DGEMM, DGEMV



! **********************************************************************************************************************************
      NROWA_T = NCOLA
      NCOLA_T = NROWA
      NROWB   = NCOLA_T
      MIN_DIM = MIN(NROWA_T,MIN(NROWB,NCOLB))
      WORK_EST = DBLE(NROWA_T)*DBLE(NROWB)*DBLE(NCOLB)

! Use BLAS for medium/large dense multiplies. Small multiplies stay in the explicit loops to avoid DGEMM call overhead.

      IF ((NCOLB == 1) .AND. (NROWB >= DGEMV_MIN_LEN)) THEN
         CALL DGEMV ( TRANSA, NROWA, NCOLA, ALPHA, A, NROWA, B(1,1), 1, BETA, C(1,1), 1 )
         RETURN
      ELSE IF ((NCOLA == 1) .AND. (NROWB >= DGEMV_MIN_LEN)) THEN
         CALL DGEMV ( TRANS, NROWB, NCOLB, ALPHA, B, NROWB, A(1,1), 1, BETA, C(1,1), 1 )
         RETURN
      ENDIF

      IF ((MIN_DIM >= DGEMM_MIN_DIM) .AND. (WORK_EST >= DGEMM_MIN_WORK)) THEN
         CALL DGEMM ( TRANSA, TRANSB, NROWA_T, NCOLB, NROWB, ALPHA, A, NROWA, B, NROWA, BETA, C, NROWA_T )
         RETURN
      ENDIF

! Initialize outputs

! Multiply A' x B

      DO I =1,NROWA_T
         DO J = 1,NCOLB
            C(I,J) = ZERO
            DO K = 1,NROWB
               C(I,J) = C(I,J) + A(K,I)*B(K,J)
            ENDDO
         ENDDO
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE MATMULT_FFF_T
