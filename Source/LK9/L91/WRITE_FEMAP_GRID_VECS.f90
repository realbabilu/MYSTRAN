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

      SUBROUTINE WRITE_FEMAP_GRID_VECS ( GRID_VEC, FEMAP_SET_ID, WHAT )

! Writes grid related vectors to FEMAP neutral file (displ, applied load, SPC and MPC forces)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NCORD, NDOFG, NGRID
      USE CONSTANTS_1, ONLY           :  ZERO
      USE FEMAP_NEU_WRITE_HELPERS, ONLY : NEU_WRITE_SET_VEC_HEADER, NEU_WRITE_TITLES, NEU_WRITE_TRIPLE_REAL,                     &
                                           NEU_WRITE_TEN_IDS, NEU_WRITE_GRID_RANGE, NEU_WRITE_GRID_VALUE, NEU_WRITE_VECTOR_END
      USE MODEL_STUF, ONLY            :  CORD, GRID, GRID_ID, INV_GRID_SEQ

      USE WRITE_FEMAP_GRID_VECS_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_FEMAP_GRID_VECS'
      CHARACTER(LEN=*), INTENT(IN)    :: WHAT              ! Indicator if GRID_VEC is DISP, VELO, ACCE, OLOA, SPCF or MPCF
      CHARACTER(LEN= 5*BYTE)          :: TITLE1(4,2)       ! Titles for vectors written to NEU
      CHARACTER(LEN=20*BYTE)          :: TITLE2(2)         ! Titles for vectors written to NEU

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID      ! FEMAP set ID to write out
      INTEGER(LONG)                   :: ACID_G            ! Actual coordinate system ID for a grid
      INTEGER(LONG)                   :: GRID_NUMS(NGRID)  ! Grid IDs in global order
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: ICID              ! Internal coord sys no. corresponding to an actual coord sys no.
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID for a grid in array GRID_NUMS
      INTEGER(LONG)                   :: ID(20)            ! Vector ID's for FEMAP output
      INTEGER(LONG)                   :: IDOFG             ! A G-set DOF number
      INTEGER(LONG)                   :: J                 ! Counter
      INTEGER(LONG)                   :: NUM_COMPS         ! 6 if GRID_NUM is a physical grid, 1 if an SPOINT
      INTEGER(LONG)                   :: VEC_ID            ! Vector ID for FEMAP output
      INTEGER(LONG)                   :: VEC_ID_OFFSET     ! Offset in determining output vector ID

      REAL(DOUBLE) , INTENT(IN)       :: GRID_VEC(NDOFG)   ! G-set Vector to process
      REAL(DOUBLE)                    :: DIS(3)            ! Array of 3 translation components
      REAL(DOUBLE)                    :: PHID, THETAD      ! Outputs from subr GEN_T0L
      REAL(DOUBLE)                    :: ROT(3)            ! Array of 3 rotation components
      REAL(DOUBLE)                    :: R1_VEC(NGRID)     ! R1 rotation    component from GRID_VEC
      REAL(DOUBLE)                    :: R2_VEC(NGRID)     ! R2 rotation    component from GRID_VEC
      REAL(DOUBLE)                    :: R3_VEC(NGRID)     ! R3 rotation    component from GRID_VEC
      REAL(DOUBLE)                    :: T0G(3,3)          ! Matrix to transform offsets from global to basic coords
      REAL(DOUBLE)                    :: T1_VEC(NGRID)     ! T1 translation component from GRID_VEC
      REAL(DOUBLE)                    :: T2_VEC(NGRID)     ! T2 translation component from GRID_VEC
      REAL(DOUBLE)                    :: T3_VEC(NGRID)     ! T3 translation component from GRID_VEC
      REAL(DOUBLE)                    :: TOTR_VEC(NGRID)   ! RSS of 3 rotation    components in GRID_VEC
      REAL(DOUBLE)                    :: TOTT_VEC(NGRID)   ! RSS of 3 translation components in GRID_VEC

! **********************************************************************************************************************************
      TITLE1(1,1) = 'RSS'
      TITLE1(2,1) = 'T1'
      TITLE1(3,1) = 'T2'
      TITLE1(4,1) = 'T3'
      TITLE1(1,2) = 'RSS'
      TITLE1(2,2) = 'R1'
      TITLE1(3,2) = 'R2'
      TITLE1(4,2) = 'R3'

      IF      (WHAT == 'DISP') THEN
         VEC_ID_OFFSET = 0
         TITLE1(1,1) = 'Total'
         TITLE1(1,2) = 'Total'
         TITLE2(1) = ' Translation'
         TITLE2(2) = ' Rotation'
      ELSE IF (WHAT == 'VELO') THEN
         VEC_ID_OFFSET = 50000
         TITLE1(1,1) = 'Total'
         TITLE1(1,2) = 'Total'
         TITLE2(1) = ' Velocity'
         TITLE2(2) = ' RotVelocity'
      ELSE IF (WHAT == 'ACCE') THEN
         VEC_ID_OFFSET = 60000
         TITLE1(1,1) = 'Total'
         TITLE1(1,2) = 'Total'
         TITLE2(1) = ' Accel'
         TITLE2(2) = ' RotAccel'
      ELSE IF (WHAT == 'OLOA') THEN
         VEC_ID_OFFSET = 20000
         TITLE1(1,1) = 'Total'
         TITLE1(1,2) = 'Total'
         TITLE2(1) = ' Applied Force'
         TITLE2(2) = ' Applied Moment'
      ELSE IF (WHAT == 'SPCF') THEN
         VEC_ID_OFFSET = 30000
         TITLE1(1,1) = 'Total'
         TITLE1(1,2) = 'Total'
         TITLE2(1) = ' Constraint Force'
         TITLE2(2) = ' Constraint Moment'
      ELSE IF (WHAT == 'MPCF') THEN
         VEC_ID_OFFSET = 40000
         TITLE2(1) = ' MPC force'
         TITLE2(2) = ' MPC moment'
      ELSE
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,939) SUBR_NAME, WHAT
         WRITE(F06,939) SUBR_NAME, WHAT
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IDOFG = 0
      DO I=1,NGRID

         T1_VEC(I) = ZERO
         T2_VEC(I) = ZERO
         T3_VEC(I) = ZERO
         R1_VEC(I) = ZERO
         R2_VEC(I) = ZERO
         R3_VEC(I) = ZERO

         GRID_NUMS(I) = GRID_ID(INV_GRID_SEQ(I))
         CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
         IF (NUM_COMPS == 6) THEN
            IDOFG = IDOFG + 1 ; T1_VEC(I) = GRID_VEC(IDOFG)
            IDOFG = IDOFG + 1 ; T2_VEC(I) = GRID_VEC(IDOFG)
            IDOFG = IDOFG + 1 ; T3_VEC(I) = GRID_VEC(IDOFG)
            IDOFG = IDOFG + 1 ; R1_VEC(I) = GRID_VEC(IDOFG)
            IDOFG = IDOFG + 1 ; R2_VEC(I) = GRID_VEC(IDOFG)
            IDOFG = IDOFG + 1 ; R3_VEC(I) = GRID_VEC(IDOFG)
         ELSE
            IDOFG = IDOFG + 1 ; T1_VEC(I) = GRID_VEC(IDOFG)
         ENDIF

         TOTR_VEC(I) = DSQRT( R1_VEC(I)*R1_VEC(I) + R2_VEC(I)*R2_VEC(I) + R3_VEC(I)*R3_VEC(I) )
         TOTT_VEC(I) = DSQRT( T1_VEC(I)*T1_VEC(I) + T2_VEC(I)*T2_VEC(I) + T3_VEC(I)*T3_VEC(I) )

         IF ((WHAT == 'DISP') .OR. (WHAT == 'VELO') .OR. (WHAT == 'ACCE')) THEN
            DIS(1) = T1_VEC(I)
            DIS(2) = T2_VEC(I)
            DIS(3) = T3_VEC(I)
            ROT(1) = R1_VEC(I)
            ROT(2) = R2_VEC(I)
            ROT(3) = R3_VEC(I)

            IGRID  = INV_GRID_SEQ(I)
            ACID_G = GRID(IGRID,3)
            IF (ACID_G /= 0) THEN
               ICID = 0
               DO J=1,NCORD
                  IF (ACID_G == CORD(J,2)) THEN
                     ICID = J
                     EXIT
                  ENDIF
               ENDDO
               CALL GEN_T0L ( IGRID, ICID, THETAD, PHID, T0G )
               T1_VEC(I) = T0G(1,1)*DIS(1) + T0G(1,2)*DIS(2) + T0G(1,3)*DIS(3)
               T2_VEC(I) = T0G(2,1)*DIS(1) + T0G(2,2)*DIS(2) + T0G(2,3)*DIS(3)
               T3_VEC(I) = T0G(3,1)*DIS(1) + T0G(3,2)*DIS(2) + T0G(3,3)*DIS(3)
               R1_VEC(I) = T0G(1,1)*ROT(1) + T0G(1,2)*ROT(2) + T0G(1,3)*ROT(3)
               R2_VEC(I) = T0G(2,1)*ROT(1) + T0G(2,2)*ROT(2) + T0G(2,3)*ROT(3)
               R3_VEC(I) = T0G(3,1)*ROT(1) + T0G(3,2)*ROT(2) + T0G(3,3)*ROT(3)
            ELSE
               T1_VEC(I) = DIS(1)
               T2_VEC(I) = DIS(2)
               T3_VEC(I) = DIS(3)
               R1_VEC(I) = ROT(1)
               R2_VEC(I) = ROT(2)
               R3_VEC(I) = ROT(3)
            ENDIF
         ENDIF
      ENDDO

      VEC_ID = VEC_ID_OFFSET + 1
      ID = 0
      ID(1) = VEC_ID + 1
      ID(2) = VEC_ID + 2
      ID(3) = VEC_ID + 3
      CALL WRITE_ONE_FEMAP_GRID_VEC ( SUBR_NAME, FEMAP_SET_ID, VEC_ID, TITLE1(1,1), TITLE2(1), NGRID, GRID_NUMS, TOTT_VEC, ID )

      VEC_ID = VEC_ID_OFFSET + 2
      ID = 0
      ID(1) = VEC_ID
      CALL WRITE_ONE_FEMAP_GRID_VEC ( SUBR_NAME, FEMAP_SET_ID, VEC_ID, TITLE1(2,1), TITLE2(1), NGRID, GRID_NUMS, T1_VEC, ID )

      VEC_ID = VEC_ID_OFFSET + 3
      ID = 0
      ID(2) = VEC_ID
      CALL WRITE_ONE_FEMAP_GRID_VEC ( SUBR_NAME, FEMAP_SET_ID, VEC_ID, TITLE1(3,1), TITLE2(1), NGRID, GRID_NUMS, T2_VEC, ID )

      VEC_ID = VEC_ID_OFFSET + 4
      ID = 0
      ID(3) = VEC_ID
      CALL WRITE_ONE_FEMAP_GRID_VEC ( SUBR_NAME, FEMAP_SET_ID, VEC_ID, TITLE1(4,1), TITLE2(1), NGRID, GRID_NUMS, T3_VEC, ID )

      VEC_ID = VEC_ID_OFFSET + 5
      ID = 0
      ID(1) = VEC_ID + 1
      ID(2) = VEC_ID + 2
      ID(3) = VEC_ID + 3
      CALL WRITE_ONE_FEMAP_GRID_VEC ( SUBR_NAME, FEMAP_SET_ID, VEC_ID, TITLE1(1,2), TITLE2(2), NGRID, GRID_NUMS, TOTR_VEC, ID )

      VEC_ID = VEC_ID_OFFSET + 6
      ID = 0
      ID(1) = VEC_ID
      CALL WRITE_ONE_FEMAP_GRID_VEC ( SUBR_NAME, FEMAP_SET_ID, VEC_ID, TITLE1(2,2), TITLE2(2), NGRID, GRID_NUMS, R1_VEC, ID )

      VEC_ID = VEC_ID_OFFSET + 7
      ID = 0
      ID(2) = VEC_ID
      CALL WRITE_ONE_FEMAP_GRID_VEC ( SUBR_NAME, FEMAP_SET_ID, VEC_ID, TITLE1(3,2), TITLE2(2), NGRID, GRID_NUMS, R2_VEC, ID )

      VEC_ID = VEC_ID_OFFSET + 8
      ID = 0
      ID(3) = VEC_ID
      CALL WRITE_ONE_FEMAP_GRID_VEC ( SUBR_NAME, FEMAP_SET_ID, VEC_ID, TITLE1(4,2), TITLE2(2), NGRID, GRID_NUMS, R3_VEC, ID )

      RETURN

! **********************************************************************************************************************************
  939 FORMAT(' *ERROR   939: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' WRONG VALUE = ',A,' FOR ARGUMENT "WHAT"')

      END SUBROUTINE WRITE_FEMAP_GRID_VECS

! ##################################################################################################################################

      SUBROUTINE WRITE_ONE_FEMAP_GRID_VEC ( SUBR_NAME, FEMAP_SET_ID, VEC_ID, TITLE_A, TITLE_B, NGRID, GRID_NUMS, DATA_VEC, ID )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE FEMAP_NEU_WRITE_HELPERS, ONLY : NEU_WRITE_SET_VEC_HEADER, NEU_WRITE_TITLES, NEU_WRITE_TRIPLE_REAL,                     &
                                           NEU_WRITE_TEN_IDS, NEU_WRITE_GRID_RANGE, NEU_WRITE_GRID_VALUE, NEU_WRITE_VECTOR_END

      USE WRITE_FEMAP_GRID_VECS_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: SUBR_NAME
      CHARACTER(LEN=*), INTENT(IN)    :: TITLE_A
      CHARACTER(LEN=*), INTENT(IN)    :: TITLE_B

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID, VEC_ID, NGRID
      INTEGER(LONG), INTENT(IN)       :: GRID_NUMS(NGRID), ID(20)
      INTEGER(LONG)                   :: GRID_MIN, GRID_MAX, I
      INTEGER(LONG)                   :: IARRAY(NGRID)

      REAL(DOUBLE), INTENT(IN)        :: DATA_VEC(NGRID)
      REAL(DOUBLE)                    :: SORTED_VEC(NGRID)
      REAL(DOUBLE)                    :: VEC_ABS, VEC_MAX, VEC_MIN

      CALL NEU_WRITE_SET_VEC_HEADER(FEMAP_SET_ID, VEC_ID)
      CALL NEU_WRITE_TITLES(TITLE_A, TITLE_B)
      CALL GET_VEC_MIN_MAX_ABS ( NGRID, GRID_NUMS, DATA_VEC, VEC_MIN, VEC_MAX, VEC_ABS, GRID_MIN, GRID_MAX )
      CALL NEU_WRITE_TRIPLE_REAL(VEC_MIN, VEC_MAX, VEC_ABS)
      CALL NEU_WRITE_TEN_IDS(ID(1:10))
      CALL NEU_WRITE_TEN_IDS(ID(11:20))
      CALL NEU_WRITE_GRID_RANGE(GRID_MIN, GRID_MAX)

      DO I=1,NGRID
         IARRAY(I)     = GRID_NUMS(I)
         SORTED_VEC(I) = DATA_VEC(I)
      ENDDO

      CALL SORT_INT1_REAL1 ( SUBR_NAME, 'FEMAP ARRAYS: GRID_NUMS, DATA_VEC', NGRID, IARRAY, SORTED_VEC )
      DO I=1,NGRID
         CALL NEU_WRITE_GRID_VALUE(IARRAY(I), SORTED_VEC(I))
      ENDDO
      CALL NEU_WRITE_VECTOR_END

      END SUBROUTINE WRITE_ONE_FEMAP_GRID_VEC
