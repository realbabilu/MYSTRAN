!--- cbeam add --- begin!
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

      SUBROUTINE WRITE_FEMAP_ELFO_VECS ( ELEM_TYP, NUM_FEMAP_ROWS, FEMAP_SET_ID )

! Writes element engineering forces to FEMAP neutral file for ROD, BAR, TRIA3, QUAD4, SHEAR

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_LOG, ERR, F04, F06
      USE PARAMS, ONLY                :  SUPWARN
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_NUMS, FEMAP_EL_VECS
      USE FEMAP_NEU_WRITE_HELPERS, ONLY : NEU_WRITE_ELEM_VECTOR
      USE SUBR_BEGEND_LEVELS, ONLY    :  WRITE_FEMAP_ELFO_VECS_BEGEND

      USE WRITE_FEMAP_ELFO_VECS_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_FEMAP_ELFO_VECS'
      CHARACTER( 1*BYTE)              :: CALC_WARN
      CHARACTER( 1*BYTE)              :: CENT_TOTAL
      CHARACTER( 1*BYTE)              :: COMP_DIR
      CHARACTER( 1*BYTE)              :: ENT_TYPE = '8'
      CHARACTER( 1*BYTE)              :: OUT_TYPE = '3'
      CHARACTER(LEN=*), INTENT(IN)    :: ELEM_TYP
      CHARACTER(LEN=LEN(ELEM_TYP))    :: ELEM_NAME
      CHARACTER(25*BYTE)              :: TITLE_E(12)

      INTEGER(LONG), INTENT(IN)       :: NUM_FEMAP_ROWS
      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID
      INTEGER(LONG)                   :: ELEM_MAX
      INTEGER(LONG)                   :: ELEM_MIN
      INTEGER(LONG), ALLOCATABLE      :: ELEM_NUMS(:)
      INTEGER(LONG)                   :: ELEM_NAME_LEN
      INTEGER(LONG)                   :: I, J
      INTEGER(LONG)                   :: ID(20)
      INTEGER(LONG)                   :: VEC_ID_OFFSET
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = WRITE_FEMAP_ELFO_VECS_BEGEND

      REAL(DOUBLE), ALLOCATABLE       :: ELEM_VECS(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: ELEM_VEC(:)
      REAL(DOUBLE)                    :: VEC_ABS
      REAL(DOUBLE)                    :: VEC_MAX
      REAL(DOUBLE)                    :: VEC_MIN

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************
      ELEM_NAME_LEN = LEN(ELEM_TYP)
      ELEM_NAME(1:) = ELEM_TYP(1:)
      DO I=LEN(ELEM_TYP),1,-1
         IF (ELEM_TYP(I:I) == ' ') CYCLE
         ELEM_NAME_LEN = I
         EXIT
      ENDDO

      ALLOCATE ( ELEM_NUMS(NUM_FEMAP_ROWS) )
      ALLOCATE ( ELEM_VECS(NUM_FEMAP_ROWS,12) )
      ALLOCATE ( ELEM_VEC(NUM_FEMAP_ROWS) )

      IF      (ELEM_TYP == 'ROD     ') THEN
         VEC_ID_OFFSET = 50100
      ELSE IF (ELEM_TYP == 'BAR     ') THEN
         VEC_ID_OFFSET = 50200
      ELSE IF (ELEM_TYP == 'TRIA3K  ') THEN
         VEC_ID_OFFSET = 50300
      ELSE IF (ELEM_TYP == 'TRIA3   ') THEN
         VEC_ID_OFFSET = 50400
      ELSE IF (ELEM_TYP == 'QUAD4K  ') THEN
         VEC_ID_OFFSET = 50500
      ELSE IF (ELEM_TYP == 'QUAD4   ') THEN
         VEC_ID_OFFSET = 50600
      ELSE IF (ELEM_TYP == 'SHEAR   ') THEN
         VEC_ID_OFFSET = 50700
      ELSE IF (ELEM_TYP == 'ELAS1   ') THEN
         VEC_ID_OFFSET = 50800
      ELSE IF (ELEM_TYP == 'ELAS2   ') THEN
         VEC_ID_OFFSET = 50900
      ELSE IF (ELEM_TYP == 'ELAS3   ') THEN
         VEC_ID_OFFSET = 51000
      ELSE IF (ELEM_TYP == 'ELAS4   ') THEN
         VEC_ID_OFFSET = 51100
      ELSE IF (ELEM_TYP == 'BUSH    ') THEN
         VEC_ID_OFFSET = 51200
      ELSE IF (ELEM_TYP == 'BEAM    ') THEN
         VEC_ID_OFFSET = 51300
      ELSE
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,943) TRIM(ELEM_TYP), 'ELEM FORCE', TRIM(SUBR_NAME)
         IF (SUPWARN == 'N') WRITE(F06,943) TRIM(ELEM_TYP), 'ELEM FORCE', TRIM(SUBR_NAME)
      ENDIF

      IF (ELEM_TYP == 'BEAM    ') THEN

         TITLE_E( 1) = 'EndA Plane1 Moment'
         TITLE_E( 2) = 'EndA Plane2 Moment'
         TITLE_E( 3) = 'EndB Plane1 Moment'
         TITLE_E( 4) = 'EndB Plane2 Moment'
         TITLE_E( 5) = 'EndA Pl1 Shear Force'
         TITLE_E( 6) = 'EndA Pl2 Shear Force'
         TITLE_E( 7) = 'EndB Pl1 Shear Force'
         TITLE_E( 8) = 'EndB Pl2 Shear Force'
         TITLE_E( 9) = 'EndA Axial Force'
         TITLE_E(10) = 'EndB Axial Force'
         TITLE_E(11) = 'EndA Torque'
         TITLE_E(12) = 'EndB Torque'

         CALC_WARN  = '0'
         COMP_DIR   = '3'
         CENT_TOTAL = '1'
         DO J=1,12
            ID = 0
            CALL WRITE_ELFO_COLUMN ( J, VEC_ID_OFFSET + J, ID )
         ENDDO

      ELSE IF ((ELEM_TYP == 'BAR     ') .OR. (ELEM_TYP == 'ROD     ')) THEN

         TITLE_E( 1) = 'EndA Plane1 Moment'
         TITLE_E( 2) = 'EndB Plane1 Moment'
         TITLE_E( 3) = 'EndA Plane2 Moment'
         TITLE_E( 4) = 'EndB Plane2 Moment'
         TITLE_E( 5) = 'EndA Pl1 Shear Force'
         TITLE_E( 6) = 'EndB Pl1 Shear Force'
         TITLE_E( 7) = 'EndA Pl2 Shear Force'
         TITLE_E( 8) = 'EndB Pl2 Shear Force'
         TITLE_E( 9) = 'EndA Axial Force'
         TITLE_E(10) = 'EndB Axial Force'
         TITLE_E(11) = 'EndA Torque'
         TITLE_E(12) = 'EndB Torque'

         DO I=1,NUM_FEMAP_ROWS
            ELEM_VECS(I, 1) = FEMAP_EL_VECS(I,1)
            ELEM_VECS(I, 2) = FEMAP_EL_VECS(I,2)
            ELEM_VECS(I, 3) = FEMAP_EL_VECS(I,3)
            ELEM_VECS(I, 4) = FEMAP_EL_VECS(I,4)
            ELEM_VECS(I, 5) = FEMAP_EL_VECS(I,5)
            ELEM_VECS(I, 6) = FEMAP_EL_VECS(I,5)
            ELEM_VECS(I, 7) = FEMAP_EL_VECS(I,6)
            ELEM_VECS(I, 8) = FEMAP_EL_VECS(I,6)
            ELEM_VECS(I, 9) = FEMAP_EL_VECS(I,7)
            ELEM_VECS(I,10) = FEMAP_EL_VECS(I,7)
            ELEM_VECS(I,11) = FEMAP_EL_VECS(I,8)
            ELEM_VECS(I,12) = FEMAP_EL_VECS(I,8)
         ENDDO

         IF (ELEM_TYP == 'ROD     ') THEN
            CALC_WARN  = '0'
            COMP_DIR   = '0'
            CENT_TOTAL = '1'
            DO J=9,12,2
               ID = 0
               ID(1) = VEC_ID_OFFSET + J
               ID(2) = VEC_ID_OFFSET + J + 1
               CALL WRITE_ELFO_COPY_COLUMN ( J  , VEC_ID_OFFSET + J    , ID )
               CALL WRITE_ELFO_COPY_COLUMN ( J+1, VEC_ID_OFFSET + J + 1, ID )
            ENDDO
         ELSE IF (ELEM_TYP == 'BAR     ') THEN
            CALC_WARN  = '0'
            COMP_DIR   = '3'
            CENT_TOTAL = '1'
            DO J=1,12,2
               ID = 0
               ID(1) = VEC_ID_OFFSET + J
               ID(2) = VEC_ID_OFFSET + J + 1
               CALL WRITE_ELFO_COPY_COLUMN ( J  , VEC_ID_OFFSET + J    , ID )
               CALL WRITE_ELFO_COPY_COLUMN ( J+1, VEC_ID_OFFSET + J + 1, ID )
            ENDDO
         ENDIF

      ELSE IF ((ELEM_TYP == 'TRIA3K  ') .OR. (ELEM_TYP == 'TRIA3   ') .OR.                                                         &
               (ELEM_TYP == 'QUAD4K  ') .OR. (ELEM_TYP == 'QUAD4   ') .OR. (ELEM_TYP == 'SHEAR   ')) THEN

         TITLE_E( 1) = 'X  Membrane Force'
         TITLE_E( 2) = 'Y  Membrane Force'
         TITLE_E( 3) = 'XY Membrane Force'
         TITLE_E( 4) = 'X  Moment'
         TITLE_E( 5) = 'Y  Moment'
         TITLE_E( 6) = 'XY Moment'
         TITLE_E( 7) = 'X  TransShear Force'
         TITLE_E( 8) = 'Y  TransShear Force'

         CALC_WARN  = '0'
         COMP_DIR   = '0'
         CENT_TOTAL = '1'
         DO J=1,8
            ID = 0
            CALL WRITE_ELFO_COLUMN ( J, VEC_ID_OFFSET + J, ID )
         ENDDO

      ELSE IF ((ELEM_TYP == 'ELAS1   ') .OR. (ELEM_TYP == 'ELAS2   ') .OR.                                                         &
               (ELEM_TYP == 'ELAS3   ') .OR. (ELEM_TYP == 'ELAS4   ')) THEN

         TITLE_E( 1) = 'Spring Force'

         CALC_WARN  = '0'
         COMP_DIR   = '0'
         CENT_TOTAL = '1'
         ID = 0
         CALL WRITE_ELFO_COLUMN ( 1, VEC_ID_OFFSET + 1, ID )

      ELSE IF (ELEM_TYP == 'BUSH    ') THEN

         TITLE_E( 1) = 'Force XE'
         TITLE_E( 2) = 'Force YE'
         TITLE_E( 3) = 'Force ZE'
         TITLE_E( 4) = 'Moment XE'
         TITLE_E( 5) = 'Moment YE'
         TITLE_E( 6) = 'Moment ZE'

         CALC_WARN  = '0'
         COMP_DIR   = '0'
         CENT_TOTAL = '1'
         DO J=1,6
            ID = 0
            CALL WRITE_ELFO_COLUMN ( J, VEC_ID_OFFSET + J, ID )
         ENDDO
      ENDIF

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      IF (ALLOCATED(ELEM_NUMS)) DEALLOCATE ( ELEM_NUMS )
      IF (ALLOCATED(ELEM_VECS)) DEALLOCATE ( ELEM_VECS )
      IF (ALLOCATED(ELEM_VEC )) DEALLOCATE ( ELEM_VEC  )

      RETURN

! **********************************************************************************************************************************
  943 FORMAT(' *WARNING    : ELEMENT TYPE = "',A,'" FOR FEMAP ',A,' OUTPUT IN SUBROUTINE ',A,' HAS NOT BEEN PROGRAMMED')

      CONTAINS

      SUBROUTINE WRITE_ELFO_COLUMN ( COL_NUM, CUR_VEC_ID, CUR_ID )

      INTEGER(LONG), INTENT(IN)       :: COL_NUM, CUR_VEC_ID
      INTEGER(LONG), INTENT(IN)       :: CUR_ID(20)
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: LOCAL_ID(20)

      LOCAL_ID = CUR_ID
      DO K=1,NUM_FEMAP_ROWS
         ELEM_VEC(K)  = FEMAP_EL_VECS(K,COL_NUM)
         ELEM_NUMS(K) = FEMAP_EL_NUMS(K,1)
      ENDDO
      CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
      CALL NEU_WRITE_ELEM_VECTOR ( FEMAP_SET_ID, CUR_VEC_ID, ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(COL_NUM),                        &
                                   VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE,                             &
                                   CALC_WARN, COMP_DIR, CENT_TOTAL, LOCAL_ID, NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC )

      END SUBROUTINE WRITE_ELFO_COLUMN

      SUBROUTINE WRITE_ELFO_COPY_COLUMN ( COL_NUM, CUR_VEC_ID, CUR_ID )

      INTEGER(LONG), INTENT(IN)       :: COL_NUM, CUR_VEC_ID
      INTEGER(LONG), INTENT(IN)       :: CUR_ID(20)
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: LOCAL_ID(20)

      LOCAL_ID = CUR_ID
      DO K=1,NUM_FEMAP_ROWS
         ELEM_VEC(K)  = ELEM_VECS(K,COL_NUM)
         ELEM_NUMS(K) = FEMAP_EL_NUMS(K,1)
      ENDDO
      CALL GET_VEC_MIN_MAX_ABS ( NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC, VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX )
      CALL NEU_WRITE_ELEM_VECTOR ( FEMAP_SET_ID, CUR_VEC_ID, ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(COL_NUM),                        &
                                   VEC_MIN, VEC_MAX, VEC_ABS, ELEM_MIN, ELEM_MAX, OUT_TYPE, ENT_TYPE,                             &
                                   CALC_WARN, COMP_DIR, CENT_TOTAL, LOCAL_ID, NUM_FEMAP_ROWS, ELEM_NUMS, ELEM_VEC )

      END SUBROUTINE WRITE_ELFO_COPY_COLUMN

      END SUBROUTINE WRITE_FEMAP_ELFO_VECS

!---  cbeam add --- end!
