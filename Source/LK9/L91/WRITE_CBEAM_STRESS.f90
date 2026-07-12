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
! LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
! _______________________________________________________________________________________________________
!
! End MIT license text.

      SUBROUTINE WRITE_CBEAM_STRESS (NUM, WRITE_F06)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE LINK9_STUFF, ONLY           :  CBEAM_XL_OUT, EID_OUT_ARRAY, GID_OUT_ARRAY, OGEL
      USE FAST_OUTPUT_FORMATTERS, ONLY:  FAST_FMT_F06_E14_6, FAST_FMT_F9_3, FAST_FMT_I8_RJ

      USE WRITE_CBEAM_STRESS_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_CBEAM_STRESS'

      INTEGER(LONG), INTENT(IN)       :: NUM
      LOGICAL,       INTENT(IN)       :: WRITE_F06
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: IBEG
      INTEGER(LONG)                   :: IEND
      INTEGER(LONG)                   :: ISTA
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: NSTA
      INTEGER(LONG)                   :: ELEMENT_ID
      INTEGER(LONG)                   :: GRID_ID

! --- cbeam_stations begin --- !
      IF (.NOT. WRITE_F06) THEN
         RETURN
      ENDIF

      I = 1
      DO WHILE (I <= NUM)
         IBEG = I
         ELEMENT_ID = EID_OUT_ARRAY(I,1)
         DO WHILE ((I <= NUM) .AND. (EID_OUT_ARRAY(I,1) == ELEMENT_ID))
            I = I + 1
         ENDDO
         IEND = I - 1
         NSTA = IEND - IBEG + 1

         IF (NSTA >= 1) THEN
            WRITE(F06,*)
            WRITE(F06,'(I1,8X,I8)') 0, ELEMENT_ID
            DO ISTA = 1, NSTA
               IF (ISTA == 1) THEN
                   GRID_ID = GID_OUT_ARRAY(IBEG,2)
               ELSE IF (ISTA == NSTA) THEN
                  GRID_ID = GID_OUT_ARRAY(IBEG,3)
               ELSE
                  GRID_ID = 0
               ENDIF

               ! Beam stress/strain results are stored in OGEL as a pair of
               ! rows per output station. The first row contains the SXC..M.S.-T
               ! values and the second paired row carries the compressive margin.
               K = 2*(IBEG + ISTA - 2) + 1

               CALL WRITE_CBEAM_STATION_LINE ( GRID_ID, CBEAM_XL_OUT(IBEG + ISTA - 1),                         &
                                               OGEL(K,1), OGEL(K,2), OGEL(K,3), OGEL(K,4),                     &
                                               OGEL(K,6), OGEL(K,7), OGEL(K,8), OGEL(K + 1,8) )
            ENDDO
         ENDIF
      ENDDO
! --- cbeam_stations end --- !

      RETURN

      CONTAINS

      SUBROUTINE WRITE_CBEAM_STATION_LINE ( GRID_ID, XL, V1, V2, V3, V4, V5, V6, V7, V8 )

      INTEGER(LONG), INTENT(IN)       :: GRID_ID
      REAL(DOUBLE), INTENT(IN)        :: XL, V1, V2, V3, V4, V5, V6, V7, V8

      CHARACTER(128*BYTE)             :: LINE_BUF
      CHARACTER(8*BYTE)               :: I8_TEXT
      CHARACTER(9*BYTE)               :: F9_TEXT
      CHARACTER(14*BYTE)              :: E14_TEXT
      INTEGER(LONG)                   :: POS

      LINE_BUF = ' '
      CALL FAST_FMT_I8_RJ ( GRID_ID, I8_TEXT )
      LINE_BUF(2:9) = I8_TEXT
      LINE_BUF(12:20) = ' '
      CALL FAST_FMT_F9_3 ( XL, F9_TEXT )
      LINE_BUF(12:20) = F9_TEXT

      POS = 21
      CALL FAST_FMT_F06_E14_6 ( V1, E14_TEXT ); LINE_BUF(POS:POS+13) = E14_TEXT; POS = POS + 14
      CALL FAST_FMT_F06_E14_6 ( V2, E14_TEXT ); LINE_BUF(POS:POS+13) = E14_TEXT; POS = POS + 14
      CALL FAST_FMT_F06_E14_6 ( V3, E14_TEXT ); LINE_BUF(POS:POS+13) = E14_TEXT; POS = POS + 14
      CALL FAST_FMT_F06_E14_6 ( V4, E14_TEXT ); LINE_BUF(POS:POS+13) = E14_TEXT; POS = POS + 14
      CALL FAST_FMT_F06_E14_6 ( V5, E14_TEXT ); LINE_BUF(POS:POS+13) = E14_TEXT; POS = POS + 14
      CALL FAST_FMT_F06_E14_6 ( V6, E14_TEXT ); LINE_BUF(POS:POS+13) = E14_TEXT; POS = POS + 14
      CALL FAST_FMT_F06_E14_6 ( V7, E14_TEXT ); LINE_BUF(POS:POS+13) = E14_TEXT; POS = POS + 14
      CALL FAST_FMT_F06_E14_6 ( V8, E14_TEXT ); LINE_BUF(POS:POS+13) = E14_TEXT; POS = POS + 14

      WRITE(F06,'(A)') LINE_BUF(1:POS-1)

      END SUBROUTINE WRITE_CBEAM_STATION_LINE

      END SUBROUTINE WRITE_CBEAM_STRESS
