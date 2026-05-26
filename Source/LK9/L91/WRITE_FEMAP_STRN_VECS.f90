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

      SUBROUTINE WRITE_FEMAP_STRN_VECS ( ELEM_TYP, IS_PCOMP, NUM_FEMAP_ROWS, FEMAP_SET_ID )

! Writes elem strain to FEMAP neutral file for TRIA3, QUAD4, SHEAR, HEXA, PENTA, TETRA4

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, NEU
      USE PARAMS, ONLY                :  SUPWARN
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NGRID, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRN_OPT
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_NUMS, FEMAP_EL_VECS

      USE WRITE_FEMAP_STRN_VECS_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_FEMAP_STRN_VECS'
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
      INTEGER(LONG)                   :: ELEM_NAME_LEN          ! Length of ELEM_TYP without trailing blanks

                                                                ! Col from FEMAP_EL_NUMS (elem ID's)
      INTEGER(LONG), ALLOCATABLE      :: ELEM_NUMS(:)

      INTEGER(LONG)                   :: I,J                    ! DO loop indices
      INTEGER(LONG)                   :: ID(20)                 ! Vector ID's for FEMAP output
      INTEGER(LONG)                   :: VEC_ID_OFFSET          ! Offset in determining output vector ID
      INTEGER(LONG)                   :: VEC_ID                 ! Vector ID for FEMAP output


      REAL(DOUBLE), ALLOCATABLE       :: ELEM_VEC(:)                 ! One column from FEMAP_EL_VECS
      REAL(DOUBLE)                    :: CURV_X
      REAL(DOUBLE)                    :: CURV_XY
      REAL(DOUBLE)                    :: CURV_Y
      REAL(DOUBLE)                    :: DZ
      REAL(DOUBLE)                    :: MID_ANGLE
      REAL(DOUBLE)                    :: MID_MEAN
      REAL(DOUBLE)                    :: MID_SMAJ
      REAL(DOUBLE)                    :: MID_SMIN
      REAL(DOUBLE)                    :: MID_SX
      REAL(DOUBLE)                    :: MID_SXY
      REAL(DOUBLE)                    :: MID_SXYMAX
      REAL(DOUBLE)                    :: MID_SY
      REAL(DOUBLE)                    :: MID_VONMISES
      REAL(DOUBLE)                    :: VEC_ABS                ! Abs value in vector
      REAL(DOUBLE)                    :: VEC_MAX                ! Max value in vector
      REAL(DOUBLE)                    :: VEC_MIN                ! Min value in vector
      REAL(DOUBLE)                    :: Z_BOT
      REAL(DOUBLE)                    :: Z_TOP


! **********************************************************************************************************************************
! --- warning_reduce-v2 begin --- !
      VEC_ID_OFFSET = 0
! --- warning_reduce-v2 end --- !
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

      IF      (ELEM_TYP == 'TRIA3K  ') THEN
         VEC_ID_OFFSET = 70400
      ELSE IF (ELEM_TYP == 'TRIA3   ') THEN
         VEC_ID_OFFSET = 70500
      ELSE IF (ELEM_TYP == 'QUAD4K  ') THEN
         VEC_ID_OFFSET = 70600
      ELSE IF ((ELEM_TYP == 'QUAD4   ') .OR. (ELEM_TYP == 'QUADR   ')) THEN
         VEC_ID_OFFSET = 70700
      ELSE IF (ELEM_TYP == 'HEXA8   ') THEN
         VEC_ID_OFFSET = 70800
      ELSE IF (ELEM_TYP == 'HEXA20  ') THEN
         VEC_ID_OFFSET = 70900
      ELSE IF (ELEM_TYP == 'PENTA6  ') THEN
         VEC_ID_OFFSET = 71000
      ELSE IF (ELEM_TYP == 'PENTA15 ') THEN
         VEC_ID_OFFSET = 71100
      ELSE IF (ELEM_TYP == 'TETRA4  ') THEN
         VEC_ID_OFFSET = 71200
      ELSE IF (ELEM_TYP == 'TETRA10 ') THEN
         VEC_ID_OFFSET = 71300
      ELSE IF (ELEM_TYP == 'CPYRAM5 ') THEN
         VEC_ID_OFFSET = 71500
      ELSE IF (ELEM_TYP == 'CPYRAM14') THEN
         VEC_ID_OFFSET = 71600
      ELSE IF (ELEM_TYP == 'SHEAR   ') THEN
         VEC_ID_OFFSET = 71400
      ELSE
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,943) TRIM(ELEM_TYP), 'STRAIN', TRIM(SUBR_NAME)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,943) TRIM(ELEM_TYP), 'STRAIN', TRIM(SUBR_NAME)
         ENDIF
      ENDIF

! Process elements

      IF ((ELEM_TYP(1:5) == 'TRIA3') .OR. (ELEM_TYP(1:5) == 'QUAD4') .OR. (ELEM_TYP == 'QUADR   ')) THEN

         IF (IS_PCOMP == 'N') THEN

            TITLE_E( 1) = 'Plate Mid X Normal Strain' ;   CALC_WARN( 1) = '0';   COMP_DIR( 1) = '0';   CENT_TOTAL( 1) = '1'
            TITLE_E( 2) = 'Plate Mid Y Normal Strain' ;   CALC_WARN( 2) = '0';   COMP_DIR( 2) = '0';   CENT_TOTAL( 2) = '1'
            TITLE_E( 3) = 'Plate Mid XY Shear Strain' ;   CALC_WARN( 3) = '0';   COMP_DIR( 3) = '0';   CENT_TOTAL( 3) = '1'
            TITLE_E( 4) = 'Plate Mid MajorPrn Strain' ;   CALC_WARN( 4) = '1';   COMP_DIR( 4) = '0';   CENT_TOTAL( 4) = '1'
            TITLE_E( 5) = 'Plate Mid MinorPrn Strain' ;   CALC_WARN( 5) = '1';   COMP_DIR( 5) = '0';   CENT_TOTAL( 5) = '1'
            TITLE_E( 6) = 'Plate Mid PrnStrain Angle' ;   CALC_WARN( 6) = '1';   COMP_DIR( 6) = '0';   CENT_TOTAL( 6) = '1'
            TITLE_E( 7) = 'Plate Mid Mean Strain'     ;   CALC_WARN( 7) = '1';   COMP_DIR( 7) = '0';   CENT_TOTAL( 7) = '1'
            TITLE_E( 8) = 'Plate Mid MaxShear Strain' ;   CALC_WARN( 8) = '1';   COMP_DIR( 8) = '0';   CENT_TOTAL( 8) = '1'
            TITLE_E( 9) = 'Plate Mid VonMises Strain' ;   CALC_WARN( 9) = '1';   COMP_DIR( 9) = '0';   CENT_TOTAL( 9) = '1'
            TITLE_E(10) = 'Plate X Normal Curvature'  ;   CALC_WARN(10) = '0';   COMP_DIR(10) = '0';   CENT_TOTAL(10) = '1'
            TITLE_E(11) = 'Plate Y Normal Curvature'  ;   CALC_WARN(11) = '0';   COMP_DIR(11) = '0';   CENT_TOTAL(11) = '1'
            TITLE_E(12) = 'Plate XY Shear Curvature'  ;   CALC_WARN(12) = '0';   COMP_DIR(12) = '0';   CENT_TOTAL(12) = '1'

            DO J=1,12
               VEC_ID = VEC_ID_OFFSET + J
               WRITE(NEU,1001) FEMAP_SET_ID, VEC_ID
               WRITE(NEU,1002) ELEM_NAME(1:ELEM_NAME_LEN), TITLE_E(J)
               DO I=1,NUM_FEMAP_ROWS
                  Z_TOP = FEMAP_EL_VECS(I, 1)
                  Z_BOT = FEMAP_EL_VECS(I,13)
                  DZ    = Z_BOT - Z_TOP
                  IF (DABS(DZ) > 1.0D-12) THEN
                     MID_SX  = (FEMAP_EL_VECS(I, 2)*Z_BOT - FEMAP_EL_VECS(I,14)*Z_TOP) / DZ
                     MID_SY  = (FEMAP_EL_VECS(I, 3)*Z_BOT - FEMAP_EL_VECS(I,15)*Z_TOP) / DZ
                     MID_SXY = (FEMAP_EL_VECS(I, 4)*Z_BOT - FEMAP_EL_VECS(I,16)*Z_TOP) / DZ
                     CURV_X  = (FEMAP_EL_VECS(I,14) - FEMAP_EL_VECS(I, 2)) / DZ
                     CURV_Y  = (FEMAP_EL_VECS(I,15) - FEMAP_EL_VECS(I, 3)) / DZ
                     CURV_XY = (FEMAP_EL_VECS(I,16) - FEMAP_EL_VECS(I, 4)) / DZ
                  ELSE
                     MID_SX  = 0.5D0*(FEMAP_EL_VECS(I, 2) + FEMAP_EL_VECS(I,14))
                     MID_SY  = 0.5D0*(FEMAP_EL_VECS(I, 3) + FEMAP_EL_VECS(I,15))
                     MID_SXY = 0.5D0*(FEMAP_EL_VECS(I, 4) + FEMAP_EL_VECS(I,16))
                     CURV_X  = ZERO
                     CURV_Y  = ZERO
                     CURV_XY = ZERO
                  ENDIF
                  CALL PRINCIPAL_2D ( MID_SX, MID_SY, MID_SXY, MID_ANGLE, MID_SMAJ, MID_SMIN, MID_SXYMAX, MID_MEAN, MID_VONMISES )
                  SELECT CASE (J)
                  CASE (1);  ELEM_VEC(I) = MID_SX
                  CASE (2);  ELEM_VEC(I) = MID_SY
                  CASE (3);  ELEM_VEC(I) = MID_SXY
                  CASE (4);  ELEM_VEC(I) = MID_SMAJ
                  CASE (5);  ELEM_VEC(I) = MID_SMIN
                  CASE (6);  ELEM_VEC(I) = MID_ANGLE
                  CASE (7);  ELEM_VEC(I) = MID_MEAN
                  CASE (8);  ELEM_VEC(I) = MID_SXYMAX
                  CASE (9);  ELEM_VEC(I) = MID_VONMISES
                  CASE (10); ELEM_VEC(I) = CURV_X
                  CASE (11); ELEM_VEC(I) = CURV_Y
                  CASE (12); ELEM_VEC(I) = CURV_XY
                  END SELECT
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

            WRITE(ERR,*) ' *WARNING    : CODE NOT WRITTEN FOR FEMAP PROCESSING OF STRAINS FOR PCOMP TYPE ELEMWMTS'
            WRITE(F06,*) ' *WARNING    : CODE NOT WRITTEN FOR FEMAP PROCESSING OF STRAINS FOR PCOMP TYPE ELEMWMTS'

         ENDIF

      ELSE IF (ELEM_TYP(1:5) == 'SHEAR') THEN

         IF (IS_PCOMP == 'N') THEN

            TITLE_E( 1) = 'Top X Direct Strain' ;   CALC_WARN( 1) = '0';   COMP_DIR( 1) = '0';   CENT_TOTAL( 1) = '1'
            TITLE_E( 2) = 'Top Y Direct Strain' ;   CALC_WARN( 2) = '0';   COMP_DIR( 2) = '0';   CENT_TOTAL( 2) = '1'
            TITLE_E( 3) = 'Top XY Shear Strain' ;   CALC_WARN( 3) = '0';   COMP_DIR( 3) = '0';   CENT_TOTAL( 3) = '1'
            TITLE_E( 4) = 'Top Maj Prn Strain'  ;   CALC_WARN( 4) = '1';   COMP_DIR( 4) = '0';   CENT_TOTAL( 4) = '1'
            TITLE_E( 5) = 'Top Min Prn Strain'  ;   CALC_WARN( 5) = '1';   COMP_DIR( 5) = '0';   CENT_TOTAL( 5) = '1'
            TITLE_E( 6) = 'Top Prn Str Angle'   ;   CALC_WARN( 6) = '1';   COMP_DIR( 6) = '0';   CENT_TOTAL( 6) = '1'
            TITLE_E( 7) = 'Top Mean Strain'     ;   CALC_WARN( 7) = '1';   COMP_DIR( 7) = '0';   CENT_TOTAL( 7) = '1'
            TITLE_E( 8) = 'Top Max Shear Strain';   CALC_WARN( 8) = '1';   COMP_DIR( 8) = '0';   CENT_TOTAL( 8) = '1'
            TITLE_E( 9) = 'Top Von Mises Strain';   CALC_WARN( 9) = '1';   COMP_DIR( 9) = '0';   CENT_TOTAL( 9) = '1'
            TITLE_E(10) = 'Top XZ Shear Strain' ;   CALC_WARN(10) = '0';   COMP_DIR(10) = '0';   CENT_TOTAL(10) = '1'
            TITLE_E(11) = 'Top YZ Shear Strain' ;   CALC_WARN(11) = '0';   COMP_DIR(11) = '0';   CENT_TOTAL(11) = '1'

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

            WRITE(ERR,*) ' *WARNING    : CODE NOT WRITTEN FOR FEMAP PROCESSING OF STRAINS FOR PCOMP TYPE ELEMWMTS'
            WRITE(F06,*) ' *WARNING    : CODE NOT WRITTEN FOR FEMAP PROCESSING OF STRAINS FOR PCOMP TYPE ELEMWMTS'

         ENDIF

      ELSE IF ((ELEM_TYP == 'HEXA8   ') .OR. (ELEM_TYP == 'PENTA6  ') .OR. (ELEM_TYP == 'TETRA4  ') .OR.                           &
               (ELEM_TYP == 'HEXA20  ') .OR. (ELEM_TYP == 'PENTA15 ') .OR. (ELEM_TYP == 'TETRA10 ') .OR.                           &
               (ELEM_TYP == 'CPYRAM5 ') .OR. (ELEM_TYP == 'CPYRAM14')) THEN

         TITLE_E( 1) = 'Solid X Normal Strain' ; CALC_WARN( 1) = '0'; COMP_DIR( 1) = '0'; CENT_TOTAL( 1) = '1'
         TITLE_E( 2) = 'Solid Y Normal Strain' ; CALC_WARN( 2) = '0'; COMP_DIR( 2) = '0'; CENT_TOTAL( 2) = '1'
         TITLE_E( 3) = 'Solid Z Normal Strain' ; CALC_WARN( 3) = '0'; COMP_DIR( 3) = '0'; CENT_TOTAL( 3) = '1'
         TITLE_E( 4) = 'Solid XY Shear Strain' ; CALC_WARN( 4) = '0'; COMP_DIR( 4) = '0'; CENT_TOTAL( 4) = '1'
         TITLE_E( 5) = 'Solid YZ Shear Strain' ; CALC_WARN( 5) = '0'; COMP_DIR( 5) = '0'; CENT_TOTAL( 5) = '1'
         TITLE_E( 6) = 'Solid ZX Shear Strain' ; CALC_WARN( 6) = '0'; COMP_DIR( 6) = '0'; CENT_TOTAL( 6) = '1'

         TITLE_E( 7) = 'Solid MajorPrn Strain' ; CALC_WARN( 7) = '0'; COMP_DIR( 7) = '0'; CENT_TOTAL( 7) = '1'
         TITLE_E( 8) = 'Solid MidPrn Strain'   ; CALC_WARN( 8) = '0'; COMP_DIR( 8) = '0'; CENT_TOTAL( 8) = '1'
         TITLE_E( 9) = 'Solid MinorPrn Strain' ; CALC_WARN( 9) = '0'; COMP_DIR( 9) = '0'; CENT_TOTAL( 9) = '1'
         TITLE_E(10) = 'Solid Mean Strain'     ; CALC_WARN(10) = '0'; COMP_DIR(10) = '0'; CENT_TOTAL(10) = '1'

         IF (STRN_OPT == 'VONMISES') THEN
            TITLE_E(11) = 'Solid VonMises Strain'; CALC_WARN(11) = '1'; COMP_DIR(11) = '0'; CENT_TOTAL(11) = '1'
            TITLE_E(12) = '  (null field)  ';      CALC_WARN(12) = '1'; COMP_DIR(12) = '0'; CENT_TOTAL(12) = '1'
         ELSE
            TITLE_E(11) = 'Solid OctDir Strain'  ; CALC_WARN(11) = '1'; COMP_DIR(11) = '0'; CENT_TOTAL(11) = '1'
            TITLE_E(12) = 'Solid OctShear Strain'; CALC_WARN(12) = '1'; COMP_DIR(12) = '0'; CENT_TOTAL(12) = '1'
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

      END SUBROUTINE WRITE_FEMAP_STRN_VECS
