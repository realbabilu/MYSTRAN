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
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
! THE USE OR OTHER DEALINGS IN THE SOFTWARE.
! _______________________________________________________________________________________________________
!
! End MIT license text.

      MODULE FEMAP_NEU_WRITE_HELPERS

! Lightweight string-assembly helpers for FEMAP-neutral output. These routines
! build each NEU record into a character buffer first, then emit it with a
! simple WRITE(NEU,'(A)'). This follows the faster "assemble string then write"
! style that benchmarks better than repeated formatted file writes.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE FAST_OUTPUT_FORMATTERS, ONLY:  FAST_FMT_E17_6, FAST_FMT_I8_RJ

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: NEU_WRITE_TEXT
      PUBLIC :: NEU_WRITE_BLOCK_END
      PUBLIC :: NEU_WRITE_SET_ID
      PUBLIC :: NEU_WRITE_SET_VEC_HEADER
      PUBLIC :: NEU_WRITE_ANALYSIS_IDS
      PUBLIC :: NEU_WRITE_ZERO_REAL
      PUBLIC :: NEU_WRITE_ZERO_INT
      PUBLIC :: NEU_WRITE_TITLES
      PUBLIC :: NEU_WRITE_TRIPLE_REAL
      PUBLIC :: NEU_WRITE_TEN_IDS
      PUBLIC :: NEU_WRITE_GRID_RANGE
      PUBLIC :: NEU_WRITE_GRID_VALUE
      PUBLIC :: NEU_WRITE_VECTOR_END
      PUBLIC :: NEU_WRITE_ELEM_RANGE
      PUBLIC :: NEU_WRITE_ELEM_FLAGS
      PUBLIC :: NEU_WRITE_ELEM_VECTOR

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_TEXT ( TEXT )

      USE IOUNT1, ONLY                :  NEU

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: TEXT

      WRITE(NEU,'(A)') TEXT

      END SUBROUTINE NEU_WRITE_TEXT

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_BLOCK_END

      CALL NEU_WRITE_TEXT('   -1')

      END SUBROUTINE NEU_WRITE_BLOCK_END

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_SET_ID ( SET_ID )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: SET_ID

      CHARACTER(9*BYTE)               :: LINE
      CHARACTER(8*BYTE)               :: SET_ID_TEXT

      CALL FAST_FMT_I8_RJ(SET_ID, SET_ID_TEXT)
      LINE = SET_ID_TEXT // ','
      CALL NEU_WRITE_TEXT(LINE)

      END SUBROUTINE NEU_WRITE_SET_ID

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_SET_VEC_HEADER ( SET_ID, VEC_ID )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: SET_ID, VEC_ID

      CHARACTER(27*BYTE)              :: LINE
      CHARACTER(8*BYTE)               :: SET_ID_TEXT, VEC_ID_TEXT

      CALL FAST_FMT_I8_RJ(SET_ID, SET_ID_TEXT)
      CALL FAST_FMT_I8_RJ(VEC_ID, VEC_ID_TEXT)
      LINE = SET_ID_TEXT // ',' // VEC_ID_TEXT // ',       1,'
      CALL NEU_WRITE_TEXT(LINE)

      END SUBROUTINE NEU_WRITE_SET_VEC_HEADER

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_ANALYSIS_IDS ( FROM_PROG, ANAL_TYPE )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: FROM_PROG, ANAL_TYPE

      CHARACTER(18*BYTE)              :: LINE
      CHARACTER(8*BYTE)               :: FROM_PROG_TEXT, ANAL_TYPE_TEXT

      CALL FAST_FMT_I8_RJ(FROM_PROG, FROM_PROG_TEXT)
      CALL FAST_FMT_I8_RJ(ANAL_TYPE, ANAL_TYPE_TEXT)
      LINE = FROM_PROG_TEXT // ',' // ANAL_TYPE_TEXT // ','
      CALL NEU_WRITE_TEXT(LINE)

      END SUBROUTINE NEU_WRITE_ANALYSIS_IDS

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_ZERO_REAL

      CALL NEU_WRITE_TEXT('0.,')

      END SUBROUTINE NEU_WRITE_ZERO_REAL

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_ZERO_INT

      CALL NEU_WRITE_TEXT('0,')

      END SUBROUTINE NEU_WRITE_ZERO_INT

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_TITLES ( TITLE1, TITLE2 )

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: TITLE1, TITLE2

      CALL NEU_WRITE_TEXT(TRIM(TITLE1)//TITLE2)

      END SUBROUTINE NEU_WRITE_TITLES

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_TRIPLE_REAL ( V1, V2, V3 )

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: V1, V2, V3

      CHARACTER(54*BYTE)              :: LINE
      CHARACTER(17*BYTE)              :: V1_TEXT, V2_TEXT, V3_TEXT

      CALL FAST_FMT_E17_6(V1, V1_TEXT)
      CALL FAST_FMT_E17_6(V2, V2_TEXT)
      CALL FAST_FMT_E17_6(V3, V3_TEXT)
      LINE = V1_TEXT // ',' // V2_TEXT // ',' // V3_TEXT // ','
      CALL NEU_WRITE_TEXT(LINE)

      END SUBROUTINE NEU_WRITE_TRIPLE_REAL

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_TEN_IDS ( IDS )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: IDS(10)

      CHARACTER(90*BYTE)              :: LINE
      CHARACTER(8*BYTE)               :: ID_TEXT
      INTEGER(LONG)                   :: I, POS

      LINE = ' '
      POS  = 1
      DO I=1,10
         CALL FAST_FMT_I8_RJ(IDS(I), ID_TEXT)
         LINE(POS:POS+7) = ID_TEXT
         LINE(POS+8:POS+8) = ','
         POS = POS + 9
      ENDDO
      CALL NEU_WRITE_TEXT(LINE)

      END SUBROUTINE NEU_WRITE_TEN_IDS

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_GRID_RANGE ( GRID_MIN, GRID_MAX )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: GRID_MIN, GRID_MAX

      CHARACTER(36*BYTE)              :: LINE
      CHARACTER(8*BYTE)               :: GRID_MIN_TEXT, GRID_MAX_TEXT

      CALL FAST_FMT_I8_RJ(GRID_MIN, GRID_MIN_TEXT)
      CALL FAST_FMT_I8_RJ(GRID_MAX, GRID_MAX_TEXT)
      LINE = GRID_MIN_TEXT // ',' // GRID_MAX_TEXT // ',       1,       7,'
      CALL NEU_WRITE_TEXT(LINE)
      CALL NEU_WRITE_TEXT('       1,       1,       1')

      END SUBROUTINE NEU_WRITE_GRID_RANGE

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_GRID_VALUE ( GRID_ID, VALUE )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: GRID_ID
      REAL(DOUBLE), INTENT(IN)        :: VALUE

      CHARACTER(27*BYTE)              :: LINE
      CHARACTER(8*BYTE)               :: GRID_ID_TEXT
      CHARACTER(17*BYTE)              :: VALUE_TEXT

      CALL FAST_FMT_I8_RJ(GRID_ID, GRID_ID_TEXT)
      CALL FAST_FMT_E17_6(VALUE, VALUE_TEXT)
      LINE = GRID_ID_TEXT // ',' // VALUE_TEXT // ','
      CALL NEU_WRITE_TEXT(LINE)

      END SUBROUTINE NEU_WRITE_GRID_VALUE

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_ELEM_RANGE ( ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: ELEM_MIN, ELEM_MAX
      CHARACTER(LEN=*), INTENT(IN)    :: OUT_TYPE, ENT_TYPE

      CHARACTER(96*BYTE)              :: LINE
      CHARACTER(8*BYTE)               :: ELEM_MIN_TEXT, ELEM_MAX_TEXT
      INTEGER(LONG)                   :: POS

      LINE = ' '
      CALL FAST_FMT_I8_RJ(ELEM_MIN, ELEM_MIN_TEXT)
      CALL FAST_FMT_I8_RJ(ELEM_MAX, ELEM_MAX_TEXT)
      LINE(1:8)   = ELEM_MIN_TEXT
      LINE(9:9)   = ','
      LINE(10:17) = ELEM_MAX_TEXT
      LINE(18:18) = ','
      POS = 25
      LINE(POS:POS+LEN(OUT_TYPE)-1) = OUT_TYPE
      LINE(POS+LEN(OUT_TYPE):POS+LEN(OUT_TYPE)) = ','
      POS = POS + 8 + LEN(OUT_TYPE)
      LINE(POS:POS+LEN(ENT_TYPE)-1) = ENT_TYPE
      LINE(POS+LEN(ENT_TYPE):POS+LEN(ENT_TYPE)) = ','
      CALL NEU_WRITE_TEXT(TRIM(LINE))

      END SUBROUTINE NEU_WRITE_ELEM_RANGE

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_ELEM_FLAGS ( CALC_WARN, COMP_DIR, CENT_TOTAL )

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: CALC_WARN, COMP_DIR, CENT_TOTAL

      CHARACTER(64*BYTE)              :: LINE

      WRITE(LINE,"(3(7X,A,','))") CALC_WARN, COMP_DIR, CENT_TOTAL
      CALL NEU_WRITE_TEXT(TRIM(LINE))

      END SUBROUTINE NEU_WRITE_ELEM_FLAGS

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_ELEM_VECTOR ( SET_ID, VEC_ID, TITLE_A, TITLE_B, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX,        &
                                         OUT_TYPE, ENT_TYPE, CALC_WARN, COMP_DIR, CENT_TOTAL, ID, NROWS, ELEM_NUMS, ELEM_VEC )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: SET_ID, VEC_ID, ELEM_MIN, ELEM_MAX, NROWS
      INTEGER(LONG), INTENT(IN)       :: ID(20), ELEM_NUMS(NROWS)
      REAL(DOUBLE), INTENT(IN)        :: ELEM_VEC(NROWS)
      REAL(DOUBLE), INTENT(IN)        :: VEC_MIN, VEC_MAX, VEC_ABS

      CHARACTER(LEN=*), INTENT(IN)    :: TITLE_A, TITLE_B, OUT_TYPE, ENT_TYPE, CALC_WARN, COMP_DIR, CENT_TOTAL

      INTEGER(LONG)                   :: I

      CALL NEU_WRITE_SET_VEC_HEADER(SET_ID, VEC_ID)
      CALL NEU_WRITE_TITLES(TITLE_A, TITLE_B)
      CALL NEU_WRITE_TRIPLE_REAL(VEC_MIN, VEC_MAX, VEC_ABS)
      CALL NEU_WRITE_TEN_IDS(ID(1:10))
      CALL NEU_WRITE_TEN_IDS(ID(11:20))
      CALL NEU_WRITE_ELEM_RANGE(ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE)
      CALL NEU_WRITE_ELEM_FLAGS(CALC_WARN, COMP_DIR, CENT_TOTAL)
      DO I=1,NROWS
         CALL NEU_WRITE_GRID_VALUE(ELEM_NUMS(I), ELEM_VEC(I))
      ENDDO
      CALL NEU_WRITE_VECTOR_END

      END SUBROUTINE NEU_WRITE_ELEM_VECTOR

! ##################################################################################################################################

      SUBROUTINE NEU_WRITE_VECTOR_END

      CALL NEU_WRITE_TEXT('      -1,     0.          ,')

      END SUBROUTINE NEU_WRITE_VECTOR_END

      END MODULE FEMAP_NEU_WRITE_HELPERS
