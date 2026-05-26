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

      SUBROUTINE WRITE_FEMAP_STRE_VECS ( ELEM_TYP, IS_PCOMP, NUM_FEMAP_ROWS, FEMAP_SET_ID )

! Writes elem stress to FEMAP neutral file for ELAS, ROD, BAR, BEAM, TRIA3, QUAD4, SHEAR, HEXA, PENTA, TETRA4

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, NEU
      USE PARAMS, ONLY                :  SUPWARN
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NGRID, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRE_OPT
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_NUMS, FEMAP_EL_VECS

      USE WRITE_FEMAP_STRE_VECS_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_FEMAP_STRE_VECS'
      CHARACTER(LEN=*), INTENT(IN)    :: ELEM_TYP               ! Element type
      CHARACTER(LEN=*), INTENT(IN)    :: IS_PCOMP               ! 'Y'/'N' for whether elements are PCOMP
      CHARACTER( 1*BYTE)              :: CALC_WARN(22)          ! FEMAP value for record 7
      CHARACTER( 1*BYTE)              :: CENT_TOTAL(22)         ! FEMAP value for record 7
      CHARACTER( 1*BYTE)              :: COMP_DIR(22)           ! FEMAP value for record 7
      CHARACTER( 1*BYTE)              :: ENT_TYPE = '8'         ! FEMAP value for record 6
      CHARACTER(LEN=LEN(ELEM_TYP))    :: ELEM_NAME              ! ELEM_TYP with trailing blanks stripped
      CHARACTER( 1*BYTE)              :: OUT_TYPE = '4'         ! FEMAP value for record 6
      CHARACTER(25*BYTE)              :: TITLE_E(22)            ! Titles for vectors written to NEU

      INTEGER(LONG), INTENT(IN)       :: NUM_FEMAP_ROWS         ! Number of rows of FEMAP data to write
      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID           ! FEMAP set ID to write out
      INTEGER(LONG)                   :: ELEM_MAX               ! Elem ID where vector is max
      INTEGER(LONG)                   :: ELEM_MIN               ! Elem ID where vector is min

                                                                ! Col from FEMAP_EL_NUMS (elem ID's)
      INTEGER(LONG), ALLOCATABLE      :: ELEM_NUMS(:)

      INTEGER(LONG)                   :: ELEM_NAME_LEN          ! Length of ELEM_TYP without trailing blanks
      INTEGER(LONG)                   :: I,J                    ! DO loop indices
      INTEGER(LONG)                   :: ID(22)                 ! Vector ID's for FEMAP output
      INTEGER(LONG)                   :: VEC_ID_OFFSET          ! Offset in determining output vector ID
      INTEGER(LONG)                   :: VEC_ID                 ! Vector ID for FEMAP output


                                                                ! One column from FEMAP_EL_VECS
      REAL(DOUBLE), ALLOCATABLE       :: ELEM_VEC(:)

      REAL(DOUBLE)                    :: VEC_ABS                ! Abs value in vector
      REAL(DOUBLE)                    :: VEC_MAX                ! Max value in vector
      REAL(DOUBLE)                    :: VEC_MIN                ! Min value in vector



! **********************************************************************************************************************************
      ELEM_NAME_LEN = LEN(ELEM_TYP)
      ELEM_NAME(1:) = ELEM_TYP(1:)
      DO I=LEN(ELEM_TYP),1,-1
         IF (ELEM_TYP(I:I) == ' ') THEN
            CYCLE
         ELSE
            ELEM_NAME_LEN = I
            EXIT
         ENDIF
      ENDDO

      IF ((ELEM_TYP(1:5) == 'TRIA3') .OR. (ELEM_TYP(1:5) == 'QUAD4') .OR. (ELEM_TYP == 'QUADR   ') .OR.                         &
          (ELEM_TYP == 'SHEAR   ') .OR. (ELEM_TYP == 'HEXA8   ') .OR. (ELEM_TYP == 'HEXA20  ') .OR.                             &
          (ELEM_TYP == 'PENTA6  ') .OR. (ELEM_TYP == 'PENTA15 ') .OR. (ELEM_TYP == 'TETRA4  ') .OR.                             &
          (ELEM_TYP == 'TETRA10 ') .OR. (ELEM_TYP == 'CPYRAM5 ') .OR. (ELEM_TYP == 'CPYRAM14')) THEN
         ELEM_NAME(1:1) = ' '
         ELEM_NAME_LEN  = 1
      ENDIF

      ALLOCATE ( ELEM_NUMS(NUM_FEMAP_ROWS), ELEM_VEC(NUM_FEMAP_ROWS) )

      IF     (ELEM_TYP(1:4) == 'ELAS') THEN
         VEC_ID_OFFSET = 60100
      ELSE IF (ELEM_TYP == 'ROD     ') THEN
         VEC_ID_OFFSET = 60200
      ELSE IF (ELEM_TYP == 'BAR     ') THEN
         VEC_ID_OFFSET = 60300
      ELSE IF (ELEM_TYP == 'BEAM    ') THEN
         VEC_ID_OFFSET = 3138
      ELSE IF (ELEM_TYP == 'TRIA3K  ') THEN
         VEC_ID_OFFSET = 60400
      ELSE IF (ELEM_TYP == 'TRIA3   ') THEN
         VEC_ID_OFFSET = 60500
      ELSE IF (ELEM_TYP == 'QUAD4K  ') THEN
         VEC_ID_OFFSET = 60600
      ELSE IF ((ELEM_TYP == 'QUAD4   ') .OR. (ELEM_TYP == 'QUADR   ')) THEN
         VEC_ID_OFFSET = 60700
      ELSE IF (ELEM_TYP == 'HEXA8   ') THEN
         VEC_ID_OFFSET = 60800
      ELSE IF (ELEM_TYP == 'HEXA20  ') THEN
         VEC_ID_OFFSET = 60900
      ELSE IF (ELEM_TYP == 'PENTA6  ') THEN
         VEC_ID_OFFSET = 61000
      ELSE IF (ELEM_TYP == 'PENTA15 ') THEN
         VEC_ID_OFFSET = 61100
      ELSE IF (ELEM_TYP == 'TETRA4  ') THEN
         VEC_ID_OFFSET = 61200
      ELSE IF (ELEM_TYP == 'TETRA10 ') THEN
         VEC_ID_OFFSET = 61300
      ELSE IF (ELEM_TYP == 'CPYRAM5 ') THEN
         VEC_ID_OFFSET = 61500
      ELSE IF (ELEM_TYP == 'CPYRAM14') THEN
         VEC_ID_OFFSET = 61600
      ELSE IF (ELEM_TYP == 'SHEAR   ') THEN
         VEC_ID_OFFSET = 61400
      ELSE
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,943) TRIM(ELEM_TYP), 'STRESS', TRIM(SUBR_NAME)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,943) TRIM(ELEM_TYP), 'STRESS', TRIM(SUBR_NAME)
         ENDIF
      ENDIF

! Process elements

      IF (ELEM_TYP(1:4) == 'ELAS')  THEN

         TITLE_E(1) = 'EndA Stress';   CALC_WARN(1) = '0';   COMP_DIR(1) = '0';   CENT_TOTAL(1) = '1'
         TITLE_E(2) = 'EndB Stress';   CALC_WARN(2) = '0';   COMP_DIR(2) = '0';   CENT_TOTAL(2) = '1'

         DO J=1,2,2

            VEC_ID = VEC_ID_OFFSET + J
            WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
            WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J)
            DO I=1,NUM_FEMAP_ROWS
               ELEM_VEC(I)  = FEMAP_EL_VECS(I,J)
               ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
            ENDDO
            CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
            WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
            ID(1) = VEC_ID
            ID(2) = VEC_ID + 1
            DO I=3,20
               ID(I) = 0
            ENDDO
            WRITE(NEU,1004) (ID(I),I= 1,10)
            WRITE(NEU,1004) (ID(I),I=11,20)
            WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
            WRITE(NEU,1006) CALC_WARN(J), COMP_DIR(J), CENT_TOTAL(J)
            DO I=1,NUM_FEMAP_ROWS
               WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
            ENDDO
            WRITE(NEU,1008)

            VEC_ID = VEC_ID_OFFSET + J
            WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID+1
            WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J+1)
            DO I=1,NUM_FEMAP_ROWS
               ELEM_VEC(I)  = FEMAP_EL_VECS(I,J+1)
               ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
            ENDDO
            CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
            WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
            ID(1) = VEC_ID
            ID(2) = VEC_ID + 1
            DO I=3,20
               ID(I) = 0
            ENDDO
            WRITE(NEU,1004) (ID(I),I= 1,10)
            WRITE(NEU,1004) (ID(I),I=11,20)
            WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
            WRITE(NEU,1006) CALC_WARN(J+1), COMP_DIR(J+1), CENT_TOTAL(J+1)
            DO I=1,NUM_FEMAP_ROWS
               WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
            ENDDO
            WRITE(NEU,1008)

         ENDDO


      ELSE IF (ELEM_TYP == 'ROD     ') THEN

         TITLE_E(1) = 'EndA Axial Stress'    ;   CALC_WARN(1) = '0';   COMP_DIR(1) = '0';   CENT_TOTAL(1) = '1'
         TITLE_E(2) = 'EndB Axial Stress'    ;   CALC_WARN(2) = '0';   COMP_DIR(2) = '0';   CENT_TOTAL(2) = '1'
         TITLE_E(3) = 'EndA Torsional Stress';   CALC_WARN(3) = '0';   COMP_DIR(3) = '0';   CENT_TOTAL(3) = '1'
         TITLE_E(4) = 'EndB Torsional Stress';   CALC_WARN(4) = '0';   COMP_DIR(4) = '0';   CENT_TOTAL(4) = '1'

         DO J=1,4,2

            VEC_ID = VEC_ID_OFFSET + J
            WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
            WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J)
            DO I=1,NUM_FEMAP_ROWS
               ELEM_VEC(I)  = FEMAP_EL_VECS(I,J)
               ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
            ENDDO
            CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
            WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
            ID(1) = VEC_ID
            ID(2) = VEC_ID + 1
            DO I=3,20
               ID(I) = 0
            ENDDO
            WRITE(NEU,1004) (ID(I),I= 1,10)
            WRITE(NEU,1004) (ID(I),I=11,20)
            WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
            WRITE(NEU,1006) CALC_WARN(J), COMP_DIR(J), CENT_TOTAL(J)
            DO I=1,NUM_FEMAP_ROWS
               WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
            ENDDO
            WRITE(NEU,1008)

            VEC_ID = VEC_ID_OFFSET + J
            WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID+1
            WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J+1)
            DO I=1,NUM_FEMAP_ROWS
               ELEM_VEC(I)  = FEMAP_EL_VECS(I,J+1)
               ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
            ENDDO
            CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
            WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
            ID(1) = VEC_ID
            ID(2) = VEC_ID + 1
            DO I=3,20
               ID(I) = 0
            ENDDO
            WRITE(NEU,1004) (ID(I),I= 1,10)
            WRITE(NEU,1004) (ID(I),I=11,20)
            WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
            WRITE(NEU,1006) CALC_WARN(J+1), COMP_DIR(J+1), CENT_TOTAL(J+1)
            DO I=1,NUM_FEMAP_ROWS
               WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
            ENDDO
            WRITE(NEU,1008)

         ENDDO

      ELSE IF (ELEM_TYP == 'BAR     ') THEN

         TITLE_E( 1) = 'EndA Pt1 Comb Stress';   CALC_WARN( 1) = '0';   COMP_DIR( 1) = '3';   CENT_TOTAL( 1) = '1'
         TITLE_E( 2) = 'EndB Pt1 Comb Stress';   CALC_WARN( 2) = '0';   COMP_DIR( 2) = '3';   CENT_TOTAL( 2) = '1'
         TITLE_E( 3) = 'EndA Pt2 Comb Stress';   CALC_WARN( 3) = '0';   COMP_DIR( 3) = '3';   CENT_TOTAL( 3) = '1'
         TITLE_E( 4) = 'EndB Pt2 Comb Stress';   CALC_WARN( 4) = '0';   COMP_DIR( 4) = '3';   CENT_TOTAL( 4) = '1'
         TITLE_E( 5) = 'EndA Pt3 Comb Stress';   CALC_WARN( 5) = '0';   COMP_DIR( 5) = '3';   CENT_TOTAL( 5) = '1'
         TITLE_E( 6) = 'EndB Pt3 Comb Stress';   CALC_WARN( 6) = '0';   COMP_DIR( 6) = '3';   CENT_TOTAL( 6) = '1'
         TITLE_E( 7) = 'EndA Pt4 Comb Stress';   CALC_WARN( 7) = '0';   COMP_DIR( 7) = '3';   CENT_TOTAL( 7) = '1'
         TITLE_E( 8) = 'EndB Pt4 Comb Stress';   CALC_WARN( 8) = '0';   COMP_DIR( 8) = '3';   CENT_TOTAL( 8) = '1'
         TITLE_E( 9) = 'EndA Max Stress'     ;   CALC_WARN( 9) = '1';   COMP_DIR( 9) = '3';   CENT_TOTAL( 9) = '1'
         TITLE_E(10) = 'EndB Max Stress'     ;   CALC_WARN(10) = '1';   COMP_DIR(10) = '3';   CENT_TOTAL(10) = '1'
         TITLE_E(11) = 'EndA Min Stress'     ;   CALC_WARN(11) = '1';   COMP_DIR(11) = '3';   CENT_TOTAL(11) = '1'
         TITLE_E(12) = 'EndB Min Stress'     ;   CALC_WARN(12) = '1';   COMP_DIR(12) = '3';   CENT_TOTAL(12) = '1'

         DO J=1,12,2

            VEC_ID = VEC_ID_OFFSET + J
            WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
            WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J)
            DO I=1,NUM_FEMAP_ROWS
               ELEM_VEC(I)  = FEMAP_EL_VECS(I,J)
               ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
            ENDDO
            CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
            WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
            ID(1) = VEC_ID
            ID(2) = VEC_ID + 1
            DO I=3,20
               ID(I) = 0
            ENDDO
            WRITE(NEU,1004) (ID(I),I= 1,10)
            WRITE(NEU,1004) (ID(I),I=11,20)
            WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
            WRITE(NEU,1006) CALC_WARN(J), COMP_DIR(J), CENT_TOTAL(J)
            DO I=1,NUM_FEMAP_ROWS
               WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
            ENDDO
            WRITE(NEU,1008)

            VEC_ID = VEC_ID_OFFSET + J
            WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID+1
            WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J+1)
            DO I=1,NUM_FEMAP_ROWS
               ELEM_VEC(I)  = FEMAP_EL_VECS(I,J+1)
               ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
            ENDDO
            CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
            WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
            ID(1) = VEC_ID
            ID(2) = VEC_ID + 1
            DO I=3,20
               ID(I) = 0
            ENDDO
            WRITE(NEU,1004) (ID(I),I= 1,10)
            WRITE(NEU,1004) (ID(I),I=11,20)
            WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
            WRITE(NEU,1006) CALC_WARN(J+1), COMP_DIR(J+1), CENT_TOTAL(J+1)
            DO I=1,NUM_FEMAP_ROWS
               WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
            ENDDO
            WRITE(NEU,1008)

         ENDDO

      ELSE IF (ELEM_TYP == 'BEAM    ') THEN
! --- neu_upgrade begin --- !
! --- CBEAM_standard begin --- !
         TITLE_E( 1) = 'EndA A (+y +z) Stress'; CALC_WARN( 1) = '0';   COMP_DIR( 1) = '3';   CENT_TOTAL( 1) = '1'
         TITLE_E( 2) = 'EndB A (+y +z) Stress'; CALC_WARN( 2) = '0';   COMP_DIR( 2) = '3';   CENT_TOTAL( 2) = '1'
         TITLE_E( 3) = 'EndA B (-y +z) Stress'; CALC_WARN( 3) = '0';   COMP_DIR( 3) = '3';   CENT_TOTAL( 3) = '1'
         TITLE_E( 4) = 'EndB B (-y +z) Stress'; CALC_WARN( 4) = '0';   COMP_DIR( 4) = '3';   CENT_TOTAL( 4) = '1'
         TITLE_E( 5) = 'EndA C (-y -z) Stress'; CALC_WARN( 5) = '0';   COMP_DIR( 5) = '3';   CENT_TOTAL( 5) = '1'
         TITLE_E( 6) = 'EndB C (-y -z) Stress'; CALC_WARN( 6) = '0';   COMP_DIR( 6) = '3';   CENT_TOTAL( 6) = '1'
         TITLE_E( 7) = 'EndA D (+y -z) Stress'; CALC_WARN( 7) = '0';   COMP_DIR( 7) = '3';   CENT_TOTAL( 7) = '1'
         TITLE_E( 8) = 'EndB D (+y -z) Stress'; CALC_WARN( 8) = '0';   COMP_DIR( 8) = '3';   CENT_TOTAL( 8) = '1'
         TITLE_E( 9) = 'EndA Max Comb Stress';   CALC_WARN( 9) = '1';   COMP_DIR( 9) = '3';   CENT_TOTAL( 9) = '1'
         TITLE_E(10) = 'EndB Max Comb Stress';   CALC_WARN(10) = '1';   COMP_DIR(10) = '3';   CENT_TOTAL(10) = '1'
         TITLE_E(11) = 'EndA Min Comb Stress';   CALC_WARN(11) = '1';   COMP_DIR(11) = '3';   CENT_TOTAL(11) = '1'
         TITLE_E(12) = 'EndB Min Comb Stress';   CALC_WARN(12) = '1';   COMP_DIR(12) = '3';   CENT_TOTAL(12) = '1'
         TITLE_E(13) = 'Beam Tension M.S.'   ;   CALC_WARN(13) = '1';   COMP_DIR(13) = '3';   CENT_TOTAL(13) = '1'
         TITLE_E(14) = 'Beam Compression M.S.';  CALC_WARN(14) = '1';   COMP_DIR(14) = '3';   CENT_TOTAL(14) = '1'
! --- CBEAM_standard end --- !

         DO J=1,14
            DO I=1,20
               ID(I) = 0
            ENDDO
            IF (J == 1) THEN
               VEC_ID = 3139
               ID(1) = 3139
               ID(2) = 3151
            ELSE IF (J == 2) THEN
               VEC_ID = 3151
               ID(1) = 3139
               ID(2) = 3151
            ELSE IF (J == 3) THEN
               VEC_ID = 3140
               ID(1) = 3140
               ID(2) = 3152
            ELSE IF (J == 4) THEN
               VEC_ID = 3152
               ID(1) = 3140
               ID(2) = 3152
            ELSE IF (J == 5) THEN
               VEC_ID = 3141
               ID(1) = 3141
               ID(2) = 3153
            ELSE IF (J == 6) THEN
               VEC_ID = 3153
               ID(1) = 3141
               ID(2) = 3153
            ELSE IF (J == 7) THEN
               VEC_ID = 3142
               ID(1) = 3142
               ID(2) = 3154
            ELSE IF (J == 8) THEN
               VEC_ID = 3154
               ID(1) = 3142
               ID(2) = 3154
            ELSE IF (J == 9) THEN
               VEC_ID = 3143
               ID(1) = 3143
               ID(2) = 3145
            ELSE IF (J == 10) THEN
               VEC_ID = 3145
               ID(1) = 3143
               ID(2) = 3145
            ELSE IF (J == 11) THEN
               VEC_ID = 3144
               ID(1) = 3144
               ID(2) = 3146
            ELSE IF (J == 12) THEN
               VEC_ID = 3146
               ID(1) = 3144
               ID(2) = 3146
            ELSE IF (J == 13) THEN
               VEC_ID = 3155
               ID(1) = 3155
               ID(2) = 0
            ELSE IF (J == 14) THEN
               VEC_ID = 3156
               ID(1) = 3156
               ID(2) = 0
            ENDIF
            WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
            WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J)
            DO I=1,NUM_FEMAP_ROWS
               IF (J <= 12) THEN
                  ELEM_VEC(I) = FEMAP_EL_VECS(I,J)
               ELSE
                  ELEM_VEC(I) = ZERO
               ENDIF
               ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
            ENDDO
            CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
            WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
            WRITE(NEU,1004) (ID(I),I= 1,10)
            WRITE(NEU,1004) (ID(I),I=11,20)
            WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
            WRITE(NEU,1006) CALC_WARN(J), COMP_DIR(J), CENT_TOTAL(J)
            DO I=1,NUM_FEMAP_ROWS
               WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
            ENDDO
            WRITE(NEU,1008)
         ENDDO
! --- neu_upgrade end --- !

      ELSE IF ((ELEM_TYP(1:5) == 'TRIA3') .OR. (ELEM_TYP(1:5) == 'QUAD4') .OR. (ELEM_TYP == 'QUADR   ')) THEN

         IF (IS_PCOMP == 'N') THEN

            TITLE_E( 1) = 'Plate Top Fiber'           ;   CALC_WARN( 1) = '0';   COMP_DIR( 1) = '0';   CENT_TOTAL( 1) = '1'
            TITLE_E( 2) = 'Plate Bottom Fiber'        ;   CALC_WARN( 2) = '0';   COMP_DIR( 2) = '0';   CENT_TOTAL( 2) = '1'
            TITLE_E( 3) = 'Plate Top X Normal Stress' ;   CALC_WARN( 3) = '0';   COMP_DIR( 3) = '0';   CENT_TOTAL( 3) = '1'
            TITLE_E( 4) = 'Plate Top Y Normal Stress' ;   CALC_WARN( 4) = '0';   COMP_DIR( 4) = '0';   CENT_TOTAL( 4) = '1'
            TITLE_E( 5) = 'Plate Top XY Shear Stress' ;   CALC_WARN( 5) = '0';   COMP_DIR( 5) = '0';   CENT_TOTAL( 5) = '1'
            TITLE_E( 6) = 'Plate Top MajorPrn Stress' ;   CALC_WARN( 6) = '1';   COMP_DIR( 6) = '0';   CENT_TOTAL( 6) = '1'
            TITLE_E( 7) = 'Plate Top MinorPrn Stress' ;   CALC_WARN( 7) = '1';   COMP_DIR( 7) = '0';   CENT_TOTAL( 7) = '1'
            TITLE_E( 8) = 'Plate Top PrnStress Angle' ;   CALC_WARN( 8) = '1';   COMP_DIR( 8) = '0';   CENT_TOTAL( 8) = '1'
            TITLE_E( 9) = 'Plate Top Mean Stress'     ;   CALC_WARN( 9) = '1';   COMP_DIR( 9) = '0';   CENT_TOTAL( 9) = '1'
            TITLE_E(10) = 'Plate Top MaxShear Stress' ;   CALC_WARN(10) = '1';   COMP_DIR(10) = '0';   CENT_TOTAL(10) = '1'
            TITLE_E(11) = 'Plate Top VonMises Stress' ;   CALC_WARN(11) = '1';   COMP_DIR(11) = '0';   CENT_TOTAL(11) = '1'
            TITLE_E(12) = 'Plate Bot X Normal Stress' ;   CALC_WARN(12) = '0';   COMP_DIR(12) = '0';   CENT_TOTAL(12) = '1'
            TITLE_E(13) = 'Plate Bot Y Normal Stress' ;   CALC_WARN(13) = '0';   COMP_DIR(13) = '0';   CENT_TOTAL(13) = '1'
            TITLE_E(14) = 'Plate Bot XY Shear Stress' ;   CALC_WARN(14) = '0';   COMP_DIR(14) = '0';   CENT_TOTAL(14) = '1'
            TITLE_E(15) = 'Plate Bot MajorPrn Stress' ;   CALC_WARN(15) = '1';   COMP_DIR(15) = '0';   CENT_TOTAL(15) = '1'
            TITLE_E(16) = 'Plate Bot MinorPrn Stress' ;   CALC_WARN(16) = '1';   COMP_DIR(16) = '0';   CENT_TOTAL(16) = '1'
            TITLE_E(17) = 'Plate Bot PrnStress Angle' ;   CALC_WARN(17) = '1';   COMP_DIR(17) = '0';   CENT_TOTAL(17) = '1'
            TITLE_E(18) = 'Plate Bot Mean Stress'     ;   CALC_WARN(18) = '1';   COMP_DIR(18) = '0';   CENT_TOTAL(18) = '1'
            TITLE_E(19) = 'Plate Bot MaxShear Stress' ;   CALC_WARN(19) = '1';   COMP_DIR(19) = '0';   CENT_TOTAL(19) = '1'
            TITLE_E(20) = 'Plate Bot VonMises Stress' ;   CALC_WARN(20) = '1';   COMP_DIR(20) = '0';   CENT_TOTAL(20) = '1'

            DO J=1,20
               VEC_ID = VEC_ID_OFFSET + J
               WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
               WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J)
               DO I=1,NUM_FEMAP_ROWS
                  IF (J == 1) THEN
                     ELEM_VEC(I) = FEMAP_EL_VECS(I, 1)
                  ELSE IF (J == 2) THEN
                     ELEM_VEC(I) = FEMAP_EL_VECS(I,13)
                  ELSE IF ((J >= 3) .AND. (J <= 11)) THEN
                     ELEM_VEC(I) = FEMAP_EL_VECS(I,J-1)
                  ELSE
                     ELEM_VEC(I) = FEMAP_EL_VECS(I,J+2)
                  ENDIF
                  ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
               ENDDO
               CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
               WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
               DO I=1,20
                  ID(I) = 0
               ENDDO
               WRITE(NEU,1004) (ID(I),I= 1,10)
               WRITE(NEU,1004) (ID(I),I=11,20)
               WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
               WRITE(NEU,1006) CALC_WARN(J), COMP_DIR(J), CENT_TOTAL(J)
               DO I=1,NUM_FEMAP_ROWS
                  WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
               ENDDO
               WRITE(NEU,1008)
            ENDDO

         ELSE

            WRITE(ERR,*) ' *WARNING    : CODE NOT WRITTEN FOR FEMAP PROCESSING OF STRESSES FOR PCOMP TYPE ELEMWMTS'
            WRITE(F06,*) ' *WARNING    : CODE NOT WRITTEN FOR FEMAP PROCESSING OF STRESSES FOR PCOMP TYPE ELEMWMTS'

         ENDIF

      ELSE IF (ELEM_TYP(1:5) == 'SHEAR') THEN

         IF (IS_PCOMP == 'N') THEN

            TITLE_E( 1) = 'X Direct Stress' ;   CALC_WARN( 1) = '0';   COMP_DIR( 1) = '0';   CENT_TOTAL( 1) = '1'
            TITLE_E( 2) = 'Y Direct Stress' ;   CALC_WARN( 2) = '0';   COMP_DIR( 2) = '0';   CENT_TOTAL( 2) = '1'
            TITLE_E( 3) = 'XY Shear Stress' ;   CALC_WARN( 3) = '0';   COMP_DIR( 3) = '0';   CENT_TOTAL( 3) = '1'
            TITLE_E( 4) = 'Maj Prn Stress'  ;   CALC_WARN( 4) = '1';   COMP_DIR( 4) = '0';   CENT_TOTAL( 4) = '1'
            TITLE_E( 5) = 'Min Prn Stress'  ;   CALC_WARN( 5) = '1';   COMP_DIR( 5) = '0';   CENT_TOTAL( 5) = '1'
            TITLE_E( 6) = 'Prn Str Angle'   ;   CALC_WARN( 6) = '1';   COMP_DIR( 6) = '0';   CENT_TOTAL( 6) = '1'
            TITLE_E( 7) = 'Mean Stress'     ;   CALC_WARN( 7) = '1';   COMP_DIR( 7) = '0';   CENT_TOTAL( 7) = '1'
            TITLE_E( 8) = 'Max Shear Stress';   CALC_WARN( 8) = '1';   COMP_DIR( 8) = '0';   CENT_TOTAL( 8) = '1'
            TITLE_E( 9) = 'Von Mises Stress';   CALC_WARN( 9) = '1';   COMP_DIR( 9) = '0';   CENT_TOTAL( 9) = '1'
            TITLE_E(10) = 'XZ Shear Stress' ;   CALC_WARN(10) = '0';   COMP_DIR(10) = '0';   CENT_TOTAL(10) = '1'
            TITLE_E(11) = 'YZ Shear Stress' ;   CALC_WARN(11) = '0';   COMP_DIR(11) = '0';   CENT_TOTAL(11) = '1'

            DO J=1,11
               VEC_ID = VEC_ID_OFFSET + J
               WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
               WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J)
               DO I=1,NUM_FEMAP_ROWS
                  ELEM_VEC(I)  = FEMAP_EL_VECS(I,J)
                  ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
               ENDDO

               CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
               WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
               DO I=1,20
                  ID(I) = 0
               ENDDO
               WRITE(NEU,1004) (ID(I),I= 1,10)
               WRITE(NEU,1004) (ID(I),I=11,20)
               WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
               WRITE(NEU,1006) CALC_WARN(J), COMP_DIR(J), CENT_TOTAL(J)
               DO I=1,NUM_FEMAP_ROWS
                  WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
               ENDDO
               WRITE(NEU,1008)
            ENDDO

         ELSE

            WRITE(ERR,*) ' *WARNING    : CODE NOT WRITTEN FOR FEMAP PROCESSING OF STRESSES FOR PCOMP TYPE ELEMWMTS'
            WRITE(F06,*) ' *WARNING    : CODE NOT WRITTEN FOR FEMAP PROCESSING OF STRESSES FOR PCOMP TYPE ELEMWMTS'

         ENDIF

      ELSE IF ((ELEM_TYP == 'HEXA8   ') .OR. (ELEM_TYP == 'PENTA6  ') .OR. (ELEM_TYP == 'TETRA4  ') .OR.                           &
               (ELEM_TYP == 'HEXA20  ') .OR. (ELEM_TYP == 'PENTA15 ') .OR. (ELEM_TYP == 'TETRA10 ') .OR.                           &
               (ELEM_TYP == 'CPYRAM5 ') .OR. (ELEM_TYP == 'CPYRAM14')) THEN

         TITLE_E( 1) = 'Solid X Normal Stress' ; CALC_WARN( 1) = '0'; COMP_DIR( 1) = '0'; CENT_TOTAL( 1) = '1'
         TITLE_E( 2) = 'Solid Y Normal Stress' ; CALC_WARN( 2) = '0'; COMP_DIR( 2) = '0'; CENT_TOTAL( 2) = '1'
         TITLE_E( 3) = 'Solid Z Normal Stress' ; CALC_WARN( 3) = '0'; COMP_DIR( 3) = '0'; CENT_TOTAL( 3) = '1'
         TITLE_E( 4) = 'Solid XY Shear Stress' ; CALC_WARN( 4) = '0'; COMP_DIR( 4) = '0'; CENT_TOTAL( 4) = '1'
         TITLE_E( 5) = 'Solid YZ Shear Stress' ; CALC_WARN( 5) = '0'; COMP_DIR( 5) = '0'; CENT_TOTAL( 5) = '1'
         TITLE_E( 6) = 'Solid ZX Shear Stress' ; CALC_WARN( 6) = '0'; COMP_DIR( 6) = '0'; CENT_TOTAL( 6) = '1'

         TITLE_E( 7) = 'Solid MajorPrn Stress' ; CALC_WARN( 7) = '0'; COMP_DIR( 7) = '0'; CENT_TOTAL( 7) = '1'
         TITLE_E( 8) = 'Solid MidPrn Stress'   ; CALC_WARN( 8) = '0'; COMP_DIR( 8) = '0'; CENT_TOTAL( 8) = '1'
         TITLE_E( 9) = 'Solid MinorPrn Stress' ; CALC_WARN( 9) = '0'; COMP_DIR( 9) = '0'; CENT_TOTAL( 9) = '1'
         TITLE_E(10) = 'Solid Mean Stress'     ; CALC_WARN(10) = '0'; COMP_DIR(10) = '0'; CENT_TOTAL(10) = '1'

         IF (STRE_OPT == 'VONMISES') THEN
            TITLE_E(11) = 'Solid VonMises Stress'; CALC_WARN(11) = '1'; COMP_DIR(11) = '0'; CENT_TOTAL(11) = '1'
            TITLE_E(12) = '  (null field)  ';      CALC_WARN(12) = '1'; COMP_DIR(12) = '0'; CENT_TOTAL(12) = '1'
         ELSE
            TITLE_E(11) = 'Solid OctDir Stress'  ; CALC_WARN(11) = '1'; COMP_DIR(11) = '0'; CENT_TOTAL(11) = '1'
            TITLE_E(12) = 'Solid OctShear Stress'; CALC_WARN(12) = '1'; COMP_DIR(12) = '0'; CENT_TOTAL(12) = '1'
         ENDIF

         DO J=1,12
            VEC_ID = VEC_ID_OFFSET + J
            WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
            WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J)
            DO I=1,NUM_FEMAP_ROWS
               ELEM_VEC(I)  = FEMAP_EL_VECS(I,J)
               ELEM_NUMS(I) = FEMAP_EL_NUMS(I,1)
            ENDDO
            CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
            WRITE(NEU,1003) VEC_MIN, VEC_MAX, VEC_ABS
            DO I=1,20
               ID(I) = 0
            ENDDO
            WRITE(NEU,1004) (ID(I),I= 1,10)
            WRITE(NEU,1004) (ID(I),I=11,20)
            WRITE(NEU,1005) ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE
            WRITE(NEU,1006) CALC_WARN(J), COMP_DIR(J), CENT_TOTAL(J)
            DO I=1,NUM_FEMAP_ROWS
               WRITE(NEU,1007) FEMAP_EL_NUMS(I,1), ELEM_VEC(I)
            ENDDO
            WRITE(NEU,1008)
         ENDDO

      ENDIF




      IF (ALLOCATED(ELEM_NUMS)) DEALLOCATE(ELEM_NUMS)
      IF (ALLOCATED(ELEM_VEC )) DEALLOCATE(ELEM_VEC)

      RETURN

! **********************************************************************************************************************************
  943 FORMAT(' *WARNING    : ELEMENT TYPE = "',A,'" FOR FEMAP ',A,' OUTPUT IN SUBROUTINE ',A,' HAS NOT BEEN PROGRAMMED')

 1001 FORMAT(2(I8,','),'       1,')

 1002 FORMAT(A,1X,A)

 1003 FORMAT(3(1ES17.6,','))

 1004 FORMAT(10(I8,','))

 1005 FORMAT(2(I8,','),2(7X,A,','))

 1006 FORMAT(3(7X,A,','))

 1007 FORMAT(I8,',',1ES17.6,',')

 1008 FORMAT('      -1,     0.          ,')

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_FEMAP_STRE_VECS
