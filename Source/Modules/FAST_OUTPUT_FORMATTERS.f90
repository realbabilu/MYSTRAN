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

      MODULE FAST_OUTPUT_FORMATTERS

      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_CHAR, C_DOUBLE, C_INT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: FAST_FMT_F06_E14_6
      PUBLIC :: FAST_FMT_E17_6
      PUBLIC :: FAST_FMT_ES13_5
      PUBLIC :: FAST_FMT_ES14_5
      PUBLIC :: FAST_FMT_ES11_3
      PUBLIC :: FAST_FMT_ES10_2
      PUBLIC :: FAST_FMT_F8_2
      PUBLIC :: FAST_FMT_F9_3
      PUBLIC :: FAST_FMT_E9_1
      PUBLIC :: FAST_FMT_I8_RJ
      PUBLIC :: FAST_BUILD_GRID_F06_LINE
      PUBLIC :: FAST_BUILD_QUAD_1403_LINE
      PUBLIC :: FAST_BUILD_QUAD_1404_LINE
      PUBLIC :: FAST_BUILD_QUAD_1405_LINE
      PUBLIC :: FAST_BUILD_QUAD_1406_LINE
      PUBLIC :: FAST_BUILD_TRIA_1703_LINE
      PUBLIC :: FAST_BUILD_TRIA_1704_LINE
      PUBLIC :: FAST_BUILD_TRIA_1706_LINE
      PUBLIC :: FAST_BUILD_PLY_FIRST_LINE
      PUBLIC :: FAST_BUILD_PLY_CONT_LINE

      INTERFACE
         SUBROUTINE C_MYSTRAN_FMT_E14_6_F06 ( VALUE, TEXT ) BIND(C, NAME='mystran_fmt_e14_6_f06')
            IMPORT :: C_CHAR, C_DOUBLE
            REAL(C_DOUBLE), VALUE       :: VALUE
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_FMT_E14_6_F06

         SUBROUTINE C_MYSTRAN_FMT_E17_6 ( VALUE, TEXT ) BIND(C, NAME='mystran_fmt_e17_6')
            IMPORT :: C_CHAR, C_DOUBLE
            REAL(C_DOUBLE), VALUE       :: VALUE
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_FMT_E17_6

         SUBROUTINE C_MYSTRAN_FMT_ES13_5 ( VALUE, TEXT ) BIND(C, NAME='mystran_fmt_es13_5')
            IMPORT :: C_CHAR, C_DOUBLE
            REAL(C_DOUBLE), VALUE       :: VALUE
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_FMT_ES13_5

         SUBROUTINE C_MYSTRAN_FMT_ES14_5 ( VALUE, TEXT ) BIND(C, NAME='mystran_fmt_es14_5')
            IMPORT :: C_CHAR, C_DOUBLE
            REAL(C_DOUBLE), VALUE       :: VALUE
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_FMT_ES14_5

         SUBROUTINE C_MYSTRAN_FMT_ES11_3 ( VALUE, TEXT ) BIND(C, NAME='mystran_fmt_es11_3')
            IMPORT :: C_CHAR, C_DOUBLE
            REAL(C_DOUBLE), VALUE       :: VALUE
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_FMT_ES11_3

         SUBROUTINE C_MYSTRAN_FMT_ES10_2 ( VALUE, TEXT ) BIND(C, NAME='mystran_fmt_es10_2')
            IMPORT :: C_CHAR, C_DOUBLE
            REAL(C_DOUBLE), VALUE       :: VALUE
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_FMT_ES10_2

         SUBROUTINE C_MYSTRAN_FMT_F8_2 ( VALUE, TEXT ) BIND(C, NAME='mystran_fmt_f8_2')
            IMPORT :: C_CHAR, C_DOUBLE
            REAL(C_DOUBLE), VALUE       :: VALUE
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_FMT_F8_2

         SUBROUTINE C_MYSTRAN_FMT_F9_3 ( VALUE, TEXT ) BIND(C, NAME='mystran_fmt_f9_3')
            IMPORT :: C_CHAR, C_DOUBLE
            REAL(C_DOUBLE), VALUE       :: VALUE
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_FMT_F9_3

         SUBROUTINE C_MYSTRAN_FMT_E9_1 ( VALUE, TEXT ) BIND(C, NAME='mystran_fmt_e9_1')
            IMPORT :: C_CHAR, C_DOUBLE
            REAL(C_DOUBLE), VALUE       :: VALUE
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_FMT_E9_1

         SUBROUTINE C_MYSTRAN_FMT_I8_RJ ( VALUE, TEXT ) BIND(C, NAME='mystran_fmt_i8_rj')
            IMPORT :: C_CHAR, C_INT
            INTEGER(C_INT), VALUE       :: VALUE
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_FMT_I8_RJ

         SUBROUTINE C_MYSTRAN_BUILD_GRID_F06_LINE ( GID, CID, VALUES, TEXT ) BIND(C, NAME='mystran_build_grid_f06_line')
            IMPORT :: C_CHAR, C_DOUBLE, C_INT
            INTEGER(C_INT), VALUE       :: GID
            INTEGER(C_INT), VALUE       :: CID
            REAL(C_DOUBLE)              :: VALUES(6)
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_BUILD_GRID_F06_LINE

         SUBROUTINE C_MYSTRAN_BUILD_QUAD_1403_LINE ( EID, VALUES, TEXT ) BIND(C, NAME='mystran_build_quad_1403_line')
            IMPORT :: C_CHAR, C_DOUBLE, C_INT
            INTEGER(C_INT), VALUE       :: EID
            REAL(C_DOUBLE)              :: VALUES(10)
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_BUILD_QUAD_1403_LINE

         SUBROUTINE C_MYSTRAN_BUILD_QUAD_1404_LINE ( VALUES, TEXT ) BIND(C, NAME='mystran_build_quad_1404_line')
            IMPORT :: C_CHAR, C_DOUBLE
            REAL(C_DOUBLE)              :: VALUES(8)
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_BUILD_QUAD_1404_LINE

         SUBROUTINE C_MYSTRAN_BUILD_QUAD_1405_LINE ( GID, VALUES, POLY_ERR, POLY_IDX, TEXT ) BIND(C, NAME='mystran_build_quad_1405_line')
            IMPORT :: C_CHAR, C_DOUBLE, C_INT
            INTEGER(C_INT), VALUE       :: GID
            REAL(C_DOUBLE)              :: VALUES(10)
            REAL(C_DOUBLE), VALUE       :: POLY_ERR
            INTEGER(C_INT), VALUE       :: POLY_IDX
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_BUILD_QUAD_1405_LINE

         SUBROUTINE C_MYSTRAN_BUILD_QUAD_1406_LINE ( GID, VALUES, POLY_ERR, TEXT ) BIND(C, NAME='mystran_build_quad_1406_line')
            IMPORT :: C_CHAR, C_DOUBLE, C_INT
            INTEGER(C_INT), VALUE       :: GID
            REAL(C_DOUBLE)              :: VALUES(10)
            REAL(C_DOUBLE), VALUE       :: POLY_ERR
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_BUILD_QUAD_1406_LINE

         SUBROUTINE C_MYSTRAN_BUILD_TRIA_1703_LINE ( EID, VALUES, TEXT ) BIND(C, NAME='mystran_build_tria_1703_line')
            IMPORT :: C_CHAR, C_DOUBLE, C_INT
            INTEGER(C_INT), VALUE       :: EID
            REAL(C_DOUBLE)              :: VALUES(10)
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_BUILD_TRIA_1703_LINE

         SUBROUTINE C_MYSTRAN_BUILD_TRIA_1704_LINE ( VALUES, TEXT ) BIND(C, NAME='mystran_build_tria_1704_line')
            IMPORT :: C_CHAR, C_DOUBLE
            REAL(C_DOUBLE)              :: VALUES(10)
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_BUILD_TRIA_1704_LINE

         SUBROUTINE C_MYSTRAN_BUILD_TRIA_1706_LINE ( EID, GID, VALUES, TEXT ) BIND(C, NAME='mystran_build_tria_1706_line')
            IMPORT :: C_CHAR, C_DOUBLE, C_INT
            INTEGER(C_INT), VALUE       :: EID
            INTEGER(C_INT), VALUE       :: GID
            REAL(C_DOUBLE)              :: VALUES(10)
            CHARACTER(C_CHAR)           :: TEXT(*)
         END SUBROUTINE C_MYSTRAN_BUILD_TRIA_1706_LINE
      END INTERFACE

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE FAST_FMT_F06_E14_6 ( VALUE, TEXT )

      REAL(DOUBLE), INTENT(IN)         :: VALUE
      CHARACTER(14*BYTE), INTENT(OUT)  :: TEXT
      CHARACTER(C_CHAR)                :: C_TEXT(14)

      CALL C_MYSTRAN_FMT_E14_6_F06(REAL(VALUE,C_DOUBLE), C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 14_LONG)

      END SUBROUTINE FAST_FMT_F06_E14_6

! ##################################################################################################################################

      SUBROUTINE FAST_FMT_E17_6 ( VALUE, TEXT )

      REAL(DOUBLE), INTENT(IN)         :: VALUE
      CHARACTER(17*BYTE), INTENT(OUT)  :: TEXT
      CHARACTER(C_CHAR)                :: C_TEXT(17)

      CALL C_MYSTRAN_FMT_E17_6(REAL(VALUE,C_DOUBLE), C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 17_LONG)

      END SUBROUTINE FAST_FMT_E17_6

! ##################################################################################################################################

      SUBROUTINE FAST_FMT_ES13_5 ( VALUE, TEXT )

      REAL(DOUBLE), INTENT(IN)         :: VALUE
      CHARACTER(13*BYTE), INTENT(OUT)  :: TEXT
      CHARACTER(C_CHAR)                :: C_TEXT(13)

      CALL C_MYSTRAN_FMT_ES13_5(REAL(VALUE,C_DOUBLE), C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 13_LONG)

      END SUBROUTINE FAST_FMT_ES13_5

! ##################################################################################################################################

      SUBROUTINE FAST_FMT_ES14_5 ( VALUE, TEXT )

      REAL(DOUBLE), INTENT(IN)         :: VALUE
      CHARACTER(14*BYTE), INTENT(OUT)  :: TEXT

      CHARACTER(C_CHAR)                :: C_TEXT(14)

      CALL C_MYSTRAN_FMT_ES14_5(REAL(VALUE,C_DOUBLE), C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 14_LONG)

      END SUBROUTINE FAST_FMT_ES14_5

! ##################################################################################################################################

      SUBROUTINE FAST_FMT_ES11_3 ( VALUE, TEXT )

      REAL(DOUBLE), INTENT(IN)         :: VALUE
      CHARACTER(11*BYTE), INTENT(OUT)  :: TEXT
      CHARACTER(C_CHAR)                :: C_TEXT(11)

      CALL C_MYSTRAN_FMT_ES11_3(REAL(VALUE,C_DOUBLE), C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 11_LONG)

      END SUBROUTINE FAST_FMT_ES11_3

! ##################################################################################################################################

      SUBROUTINE FAST_FMT_ES10_2 ( VALUE, TEXT )

      REAL(DOUBLE), INTENT(IN)         :: VALUE
      CHARACTER(10*BYTE), INTENT(OUT)  :: TEXT

      CHARACTER(C_CHAR)                :: C_TEXT(10)

      CALL C_MYSTRAN_FMT_ES10_2(REAL(VALUE,C_DOUBLE), C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 10_LONG)

      END SUBROUTINE FAST_FMT_ES10_2

! ##################################################################################################################################

      SUBROUTINE FAST_FMT_F8_2 ( VALUE, TEXT )

      REAL(DOUBLE), INTENT(IN)         :: VALUE
      CHARACTER(8*BYTE), INTENT(OUT)   :: TEXT
      CHARACTER(C_CHAR)                :: C_TEXT(8)

      CALL C_MYSTRAN_FMT_F8_2(REAL(VALUE,C_DOUBLE), C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 8_LONG)

      END SUBROUTINE FAST_FMT_F8_2

! ##################################################################################################################################

      SUBROUTINE FAST_FMT_F9_3 ( VALUE, TEXT )

      REAL(DOUBLE), INTENT(IN)         :: VALUE
      CHARACTER(9*BYTE), INTENT(OUT)   :: TEXT

      CHARACTER(C_CHAR)                :: C_TEXT(9)

      CALL C_MYSTRAN_FMT_F9_3(REAL(VALUE,C_DOUBLE), C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 9_LONG)

      END SUBROUTINE FAST_FMT_F9_3

! ##################################################################################################################################

      SUBROUTINE FAST_FMT_E9_1 ( VALUE, TEXT )

      REAL(DOUBLE), INTENT(IN)         :: VALUE
      CHARACTER(9*BYTE), INTENT(OUT)   :: TEXT
      CHARACTER(C_CHAR)                :: C_TEXT(9)

      CALL C_MYSTRAN_FMT_E9_1(REAL(VALUE,C_DOUBLE), C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 9_LONG)

      END SUBROUTINE FAST_FMT_E9_1

! ##################################################################################################################################

      SUBROUTINE FAST_FMT_I8_RJ ( VALUE, TEXT )

      INTEGER(LONG), INTENT(IN)        :: VALUE
      CHARACTER(8*BYTE), INTENT(OUT)   :: TEXT
      CHARACTER(C_CHAR)                :: C_TEXT(8)

      CALL C_MYSTRAN_FMT_I8_RJ(INT(VALUE,C_INT), C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 8_LONG)

      END SUBROUTINE FAST_FMT_I8_RJ

! ##################################################################################################################################

      SUBROUTINE FAST_BUILD_GRID_F06_LINE ( GID, CID, VALUES, TEXT )

      INTEGER(LONG), INTENT(IN)        :: GID, CID
      REAL(DOUBLE), INTENT(IN)         :: VALUES(6)
      CHARACTER(108*BYTE), INTENT(OUT) :: TEXT

      CHARACTER(C_CHAR)                :: C_TEXT(108)
      REAL(C_DOUBLE)                   :: C_VALUES(6)

      CALL COPY_REAL_VECTOR_6(VALUES, C_VALUES)
      CALL C_MYSTRAN_BUILD_GRID_F06_LINE(INT(GID,C_INT), INT(CID,C_INT), C_VALUES, C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 108_LONG)

      END SUBROUTINE FAST_BUILD_GRID_F06_LINE

! ##################################################################################################################################

      SUBROUTINE FAST_BUILD_QUAD_1403_LINE ( EID, VALUES, TEXT )

      INTEGER(LONG), INTENT(IN)        :: EID
      REAL(DOUBLE), INTENT(IN)         :: VALUES(10)
      CHARACTER(145*BYTE), INTENT(OUT) :: TEXT

      CHARACTER(C_CHAR)                :: C_TEXT(145)
      REAL(C_DOUBLE)                   :: C_VALUES(10)

      CALL COPY_REAL_VECTOR_10(VALUES, C_VALUES)
      CALL C_MYSTRAN_BUILD_QUAD_1403_LINE(INT(EID,C_INT), C_VALUES, C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 145_LONG)

      END SUBROUTINE FAST_BUILD_QUAD_1403_LINE

! ##################################################################################################################################

      SUBROUTINE FAST_BUILD_QUAD_1404_LINE ( VALUES, TEXT )

      REAL(DOUBLE), INTENT(IN)         :: VALUES(8)
      CHARACTER(119*BYTE), INTENT(OUT) :: TEXT

      CHARACTER(C_CHAR)                :: C_TEXT(119)
      REAL(C_DOUBLE)                   :: C_VALUES(8)

      CALL COPY_REAL_VECTOR_8(VALUES, C_VALUES)
      CALL C_MYSTRAN_BUILD_QUAD_1404_LINE(C_VALUES, C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 119_LONG)

      END SUBROUTINE FAST_BUILD_QUAD_1404_LINE

! ##################################################################################################################################

      SUBROUTINE FAST_BUILD_QUAD_1405_LINE ( GID, VALUES, POLY_ERR, POLY_IDX, TEXT )

      INTEGER(LONG), INTENT(IN)        :: GID, POLY_IDX
      REAL(DOUBLE), INTENT(IN)         :: VALUES(10), POLY_ERR
      CHARACTER(157*BYTE), INTENT(OUT) :: TEXT

      CHARACTER(C_CHAR)                :: C_TEXT(157)
      REAL(C_DOUBLE)                   :: C_VALUES(10)

      CALL COPY_REAL_VECTOR_10(VALUES, C_VALUES)
      CALL C_MYSTRAN_BUILD_QUAD_1405_LINE(INT(GID,C_INT), C_VALUES, REAL(POLY_ERR,C_DOUBLE), INT(POLY_IDX,C_INT), C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 157_LONG)

      END SUBROUTINE FAST_BUILD_QUAD_1405_LINE

! ##################################################################################################################################

      SUBROUTINE FAST_BUILD_QUAD_1406_LINE ( GID, VALUES, POLY_ERR, TEXT )

      INTEGER(LONG), INTENT(IN)        :: GID
      REAL(DOUBLE), INTENT(IN)         :: VALUES(10), POLY_ERR
      CHARACTER(154*BYTE), INTENT(OUT) :: TEXT

      CHARACTER(C_CHAR)                :: C_TEXT(154)
      REAL(C_DOUBLE)                   :: C_VALUES(10)

      CALL COPY_REAL_VECTOR_10(VALUES, C_VALUES)
      CALL C_MYSTRAN_BUILD_QUAD_1406_LINE(INT(GID,C_INT), C_VALUES, REAL(POLY_ERR,C_DOUBLE), C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 154_LONG)

      END SUBROUTINE FAST_BUILD_QUAD_1406_LINE

! ##################################################################################################################################

      SUBROUTINE FAST_BUILD_TRIA_1703_LINE ( EID, VALUES, TEXT )

      INTEGER(LONG), INTENT(IN)        :: EID
      REAL(DOUBLE), INTENT(IN)         :: VALUES(10)
      CHARACTER(159*BYTE), INTENT(OUT) :: TEXT

      CHARACTER(C_CHAR)                :: C_TEXT(159)
      REAL(C_DOUBLE)                   :: C_VALUES(10)

      CALL COPY_REAL_VECTOR_10(VALUES, C_VALUES)
      CALL C_MYSTRAN_BUILD_TRIA_1703_LINE(INT(EID,C_INT), C_VALUES, C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 159_LONG)

      END SUBROUTINE FAST_BUILD_TRIA_1703_LINE

! ##################################################################################################################################

      SUBROUTINE FAST_BUILD_TRIA_1704_LINE ( VALUES, TEXT )

      REAL(DOUBLE), INTENT(IN)         :: VALUES(10)
      CHARACTER(159*BYTE), INTENT(OUT) :: TEXT

      CHARACTER(C_CHAR)                :: C_TEXT(159)
      REAL(C_DOUBLE)                   :: C_VALUES(10)

      CALL COPY_REAL_VECTOR_10(VALUES, C_VALUES)
      CALL C_MYSTRAN_BUILD_TRIA_1704_LINE(C_VALUES, C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 159_LONG)

      END SUBROUTINE FAST_BUILD_TRIA_1704_LINE

! ##################################################################################################################################

      SUBROUTINE FAST_BUILD_TRIA_1706_LINE ( EID, GID, VALUES, TEXT )

      INTEGER(LONG), INTENT(IN)        :: EID, GID
      REAL(DOUBLE), INTENT(IN)         :: VALUES(10)
      CHARACTER(159*BYTE), INTENT(OUT) :: TEXT

      CHARACTER(C_CHAR)                :: C_TEXT(159)
      REAL(C_DOUBLE)                   :: C_VALUES(10)

      CALL COPY_REAL_VECTOR_10(VALUES, C_VALUES)
      CALL C_MYSTRAN_BUILD_TRIA_1706_LINE(INT(EID,C_INT), INT(GID,C_INT), C_VALUES, C_TEXT)
      CALL COPY_C_TEXT(C_TEXT, TEXT, 159_LONG)

      END SUBROUTINE FAST_BUILD_TRIA_1706_LINE

! ##################################################################################################################################

      SUBROUTINE FAST_BUILD_PLY_FIRST_LINE ( EID, PLY_NUM, VALUES, NVALS, FTNAME, TEXT )

      INTEGER(LONG), INTENT(IN)        :: EID, PLY_NUM, NVALS
      REAL(DOUBLE), INTENT(IN)         :: VALUES(NVALS)
      CHARACTER(LEN=*), INTENT(IN)     :: FTNAME
      CHARACTER(LEN=*), INTENT(OUT)    :: TEXT

      CHARACTER(8*BYTE)                :: EID_TEXT
      CHARACTER(8*BYTE)                :: PLY_TEXT
      CHARACTER(13*BYTE)               :: ES13_TEXT
      CHARACTER(14*BYTE)               :: ES14_TEXT
      CHARACTER(9*BYTE)                :: F9_TEXT
      CHARACTER(10*BYTE)               :: ES10_TEXT
      INTEGER(LONG)                    :: POS

      TEXT = ' '
      CALL FAST_FMT_I8_RJ(EID, EID_TEXT)
      CALL FAST_FMT_I8_RJ(PLY_NUM, PLY_TEXT)

      TEXT(2:9)   = EID_TEXT
      TEXT(12:17) = PLY_TEXT(3:8)

      POS = 18
      CALL FAST_FMT_ES13_5(VALUES(1), ES13_TEXT); TEXT(POS:POS+12) = ES13_TEXT; POS = POS + 13
      CALL FAST_FMT_ES13_5(VALUES(2), ES13_TEXT); TEXT(POS:POS+12) = ES13_TEXT; POS = POS + 13
      CALL FAST_FMT_ES13_5(VALUES(3), ES13_TEXT); TEXT(POS:POS+12) = ES13_TEXT; POS = POS + 13
      TEXT(POS:POS+1) = '  '; POS = POS + 2
      CALL FAST_FMT_ES14_5(VALUES(4), ES14_TEXT); TEXT(POS:POS+13) = ES14_TEXT; POS = POS + 14
      CALL FAST_FMT_ES14_5(VALUES(5), ES14_TEXT); TEXT(POS:POS+13) = ES14_TEXT; POS = POS + 14
      CALL FAST_FMT_F9_3 (VALUES(6), F9_TEXT ); TEXT(POS:POS+8 ) = F9_TEXT ; POS = POS + 9
      CALL FAST_FMT_ES13_5(VALUES(7), ES13_TEXT); TEXT(POS:POS+12) = ES13_TEXT; POS = POS + 13
      CALL FAST_FMT_ES13_5(VALUES(8), ES13_TEXT); TEXT(POS:POS+12) = ES13_TEXT; POS = POS + 13
      CALL FAST_FMT_ES13_5(VALUES(9), ES13_TEXT); TEXT(POS:POS+12) = ES13_TEXT; POS = POS + 13

      IF (NVALS >= 10) THEN
         CALL FAST_FMT_ES10_2(VALUES(10), ES10_TEXT)
         TEXT(POS:POS+9) = ES10_TEXT
         POS = POS + 10
      ENDIF

      IF (LEN_TRIM(FTNAME) > 0) THEN
         TEXT(POS:POS+1) = '  '
         POS = POS + 2
         TEXT(POS:POS+LEN_TRIM(FTNAME)-1) = FTNAME(1:LEN_TRIM(FTNAME))
      ENDIF

      END SUBROUTINE FAST_BUILD_PLY_FIRST_LINE

! ##################################################################################################################################

      SUBROUTINE FAST_BUILD_PLY_CONT_LINE ( PLY_NUM, VALUES, NVALS, FTNAME, TEXT )

      INTEGER(LONG), INTENT(IN)        :: PLY_NUM, NVALS
      REAL(DOUBLE), INTENT(IN)         :: VALUES(NVALS)
      CHARACTER(LEN=*), INTENT(IN)     :: FTNAME
      CHARACTER(LEN=*), INTENT(OUT)    :: TEXT

      CHARACTER(8*BYTE)                :: PLY_TEXT
      CHARACTER(13*BYTE)               :: ES13_TEXT
      CHARACTER(14*BYTE)               :: ES14_TEXT
      CHARACTER(9*BYTE)                :: F9_TEXT
      CHARACTER(10*BYTE)               :: ES10_TEXT
      INTEGER(LONG)                    :: POS

      TEXT = ' '
      CALL FAST_FMT_I8_RJ(PLY_NUM, PLY_TEXT)
      TEXT(18:23) = PLY_TEXT(3:8)

      POS = 24
      CALL FAST_FMT_ES13_5(VALUES(1), ES13_TEXT); TEXT(POS:POS+12) = ES13_TEXT; POS = POS + 13
      CALL FAST_FMT_ES13_5(VALUES(2), ES13_TEXT); TEXT(POS:POS+12) = ES13_TEXT; POS = POS + 13
      CALL FAST_FMT_ES13_5(VALUES(3), ES13_TEXT); TEXT(POS:POS+12) = ES13_TEXT; POS = POS + 13
      TEXT(POS:POS+1) = '  '; POS = POS + 2
      CALL FAST_FMT_ES14_5(VALUES(4), ES14_TEXT); TEXT(POS:POS+13) = ES14_TEXT; POS = POS + 14
      CALL FAST_FMT_ES14_5(VALUES(5), ES14_TEXT); TEXT(POS:POS+13) = ES14_TEXT; POS = POS + 14
      CALL FAST_FMT_F9_3 (VALUES(6), F9_TEXT ); TEXT(POS:POS+8 ) = F9_TEXT ; POS = POS + 9
      CALL FAST_FMT_ES13_5(VALUES(7), ES13_TEXT); TEXT(POS:POS+12) = ES13_TEXT; POS = POS + 13
      CALL FAST_FMT_ES13_5(VALUES(8), ES13_TEXT); TEXT(POS:POS+12) = ES13_TEXT; POS = POS + 13
      CALL FAST_FMT_ES13_5(VALUES(9), ES13_TEXT); TEXT(POS:POS+12) = ES13_TEXT; POS = POS + 13

      IF (NVALS >= 10) THEN
         CALL FAST_FMT_ES10_2(VALUES(10), ES10_TEXT)
         TEXT(POS:POS+9) = ES10_TEXT
         POS = POS + 10
      ENDIF

      IF (LEN_TRIM(FTNAME) > 0) THEN
         TEXT(POS:POS+1) = '  '
         POS = POS + 2
         TEXT(POS:POS+LEN_TRIM(FTNAME)-1) = FTNAME(1:LEN_TRIM(FTNAME))
      ENDIF

      END SUBROUTINE FAST_BUILD_PLY_CONT_LINE

! ##################################################################################################################################

      SUBROUTINE COPY_C_TEXT ( SOURCE, DEST, NCHARS )

      CHARACTER(C_CHAR), INTENT(IN)    :: SOURCE(*)
      CHARACTER(*), INTENT(OUT)        :: DEST
      INTEGER(LONG), INTENT(IN)        :: NCHARS

      INTEGER(LONG)                    :: I

      DO I=1,NCHARS
         DEST(I:I) = TRANSFER(SOURCE(I),' ')
      ENDDO

      END SUBROUTINE COPY_C_TEXT

! ##################################################################################################################################

      SUBROUTINE COPY_REAL_VECTOR_6 ( SOURCE, DEST )

      REAL(DOUBLE), INTENT(IN)         :: SOURCE(6)
      REAL(C_DOUBLE), INTENT(OUT)      :: DEST(6)
      INTEGER(LONG)                    :: I

      DO I=1,6
         DEST(I) = REAL(SOURCE(I),C_DOUBLE)
      ENDDO

      END SUBROUTINE COPY_REAL_VECTOR_6

! ##################################################################################################################################

      SUBROUTINE COPY_REAL_VECTOR_8 ( SOURCE, DEST )

      REAL(DOUBLE), INTENT(IN)         :: SOURCE(8)
      REAL(C_DOUBLE), INTENT(OUT)      :: DEST(8)
      INTEGER(LONG)                    :: I

      DO I=1,8
         DEST(I) = REAL(SOURCE(I),C_DOUBLE)
      ENDDO

      END SUBROUTINE COPY_REAL_VECTOR_8

! ##################################################################################################################################

      SUBROUTINE COPY_REAL_VECTOR_10 ( SOURCE, DEST )

      REAL(DOUBLE), INTENT(IN)         :: SOURCE(10)
      REAL(C_DOUBLE), INTENT(OUT)      :: DEST(10)
      INTEGER(LONG)                    :: I

      DO I=1,10
         DEST(I) = REAL(SOURCE(I),C_DOUBLE)
      ENDDO

      END SUBROUTINE COPY_REAL_VECTOR_10

      END MODULE FAST_OUTPUT_FORMATTERS
