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

      SUBROUTINE WRITE_ELEM_STRESSES ( JSUB, NUM, IHEADER, NUM_PTS, ITABLE )

      ! Writes blocks of element stresses for one subcase and one element type for elements that do not have PCOMP properties, including
      ! all 1-D, 2-D, 3-D elements.
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, OP2
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, BARTOR, INT_SC_NUM, MAX_NUM_STR, NDOFR, NUM_CB_DOFS,             &
                                         NVEC, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  STR_CID
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL
      USE LINK9_STUFF, ONLY           :  CBEAM_XL_OUT, EID_OUT_ARRAY, GID_OUT_ARRAY, OGEL, SHELL_OUT_TE, POLY_FIT_ERR,          &
                                         POLY_FIT_ERR_INDEX
      USE MODEL_STUF, ONLY            :  ELEM_ONAME, ELMTYP, LABEL, SCNUM, STITLE, TITLE, TYPE
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRE_LOC, STRE_OPT, STRE_OUT, GPSTRESS_OUT, GPSTRESS_REQ, STRESS_USER_REQ
      USE FAST_OUTPUT_FORMATTERS, ONLY:  FAST_FMT_F06_E14_6, FAST_FMT_I8_RJ,                                            &
                                         FAST_BUILD_QUAD_1403_LINE, FAST_BUILD_QUAD_1404_LINE,                            &
                                         FAST_BUILD_QUAD_1405_LINE, FAST_BUILD_QUAD_1406_LINE,                            &
                                         FAST_BUILD_TRIA_1703_LINE, FAST_BUILD_TRIA_1704_LINE,                            &
                                         FAST_BUILD_TRIA_1706_LINE

      USE WRITE_ELEM_STRESSES_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_ELEM_STRESSES'
      CHARACTER(LEN=*), INTENT(IN)    :: IHEADER           ! Indicator of whether to write an output header

                                                           ! Array of different notes to write regarding poly fit errors
      CHARACTER( 50*BYTE)             :: ERR_INDEX_NOTE(MAX_NUM_STR)
      CHARACTER(128*BYTE)             :: FILL              ! Padding for output format
      CHARACTER(160*BYTE)             :: LINE_BUF
      CHARACTER(145*BYTE)             :: QUAD_CENTER_LINE
      CHARACTER(119*BYTE)             :: QUAD_LOWER_LINE
      CHARACTER(157*BYTE)             :: QUAD_GRID_NOTE_LINE
      CHARACTER(154*BYTE)             :: QUAD_GRID_LINE
      CHARACTER(LEN=LEN(ELEM_ONAME))  :: ONAME             ! Element name to write out in F06 file
      CHARACTER( 1*BYTE)              :: WRITE_NOTES = 'N' ! Indicator of whether to write any WRT_ERR_INDEX_NOTE(i)

                                                           ! Indicators of whether to write note on indices of POLY_FIT_ERR
      CHARACTER( 1*BYTE)              :: WRT_ERR_INDEX_NOTE(MAX_NUM_STR)

      INTEGER(LONG), INTENT(IN)       :: JSUB              ! Solution vector number
      INTEGER(LONG), INTENT(IN)       :: NUM               ! The number of rows of OGEL to write out
      INTEGER(LONG), INTENT(IN)       :: NUM_PTS           ! Num diff stress points for one element (3rd dim in arrays SEi, STEi)
      INTEGER(LONG), INTENT(INOUT)    :: ITABLE            ! the current op2 subtable, should be -3, -5, ...
      INTEGER(LONG)                   :: BDY_COMP          ! Component (1-6) for a boundary DOF in CB analyses
      INTEGER(LONG)                   :: BDY_GRID          ! Grid for a boundary DOF in CB analyses
      INTEGER(LONG)                   :: BDY_DOF_NUM       ! DOF number for BDY_GRID/BDY_COMP
      INTEGER(LONG)                   :: I,J,L             ! DO loop indices
      INTEGER(LONG)                   :: K                 ! Counter
      INTEGER(LONG)                   :: IBEG, IEND, IELEM, ISTA, NSTA_ELEM, NROW_ELEM
      INTEGER(LONG)                   :: NCOLS             ! Num of cols to write out


      REAL(DOUBLE)                    :: ABS_ANS(11)       ! Max ABS for all element output
      REAL(DOUBLE)                    :: ANGLE
      REAL(DOUBLE)                    :: MAX_ANS(11)       ! Max for all element output
      REAL(DOUBLE)                    :: MEAN
      REAL(DOUBLE)                    :: MIN_ANS(11)       ! Min for all element output
      REAL(DOUBLE)                    :: QUAD_VALUES_10(10)! Local contiguous copy to avoid strided slice temporaries
      REAL(DOUBLE)                    :: QUAD_VALUES_8(8)  ! Local contiguous copy to avoid strided slice temporaries
      REAL(DOUBLE)                    :: ROW_CURV(10)
      REAL(DOUBLE)                    :: ROW_MEM(10)
      REAL(DOUBLE)                    :: SMAJ
      REAL(DOUBLE)                    :: SMIN
      REAL(DOUBLE)                    :: SXYMAX
      REAL(DOUBLE)                    :: TINT, XI_STD, XI0, XI1
      REAL(DOUBLE)                    :: VONMISES
      REAL(DOUBLE)                    :: XI_RAW(11), SXC_RAW(11), SXD_RAW(11), SXE_RAW(11), SXF_RAW(11), SMAX_RAW(11),            &
                                         SMIN_RAW(11), MST_RAW(11), MSC_RAW(11)
      REAL(DOUBLE)                    :: Z1, Z2, Z_DEN
      REAL(DOUBLE), ALLOCATABLE       :: BEAM_XI(:,:), BEAM_SXC(:,:), BEAM_SXD(:,:), BEAM_SXE(:,:), BEAM_SXF(:,:),               &
                                         BEAM_SMAX(:,:), BEAM_SMIN(:,:), BEAM_MST(:,:), BEAM_MSC(:,:)
      INTEGER(LONG), ALLOCATABLE      :: BEAM_EID(:), BEAM_GRID(:,:)

      ! op2 info
      CHARACTER( 8*BYTE)              :: TABLE_NAME             ! the name of the op2 table
      INTEGER(LONG)                   :: NNODES                 ! number of nodes for the element

      ! table -3 info
      INTEGER(LONG)                   :: ANALYSIS_CODE          ! static/modal/time/etc. flag
      INTEGER(LONG)                   :: ELEMENT_TYPE           ! the OP2 flag for the element
      LOGICAL                         :: FIELD_5_INT_FLAG       ! flag to trigger FIELD5_INT_MODE vs. FIELD5_FLOAT_TIME_FREQ
      LOGICAL                         :: WRITE_F06, WRITE_OP2, OP2_OPENED   ! flag
      LOGICAL                         :: WRITE_GPSTRESS_F06
      INTEGER(LONG)                   :: FIELD5_INT_MODE        ! int value for field 5
      REAL(DOUBLE)                    :: FIELD5_FLOAT_TIME_FREQ ! float value for field 5
      REAL(DOUBLE)                    :: FIELD6_EIGENVALUE      ! float value for field 6
      CHARACTER(LEN=128)              :: TITLEI                 ! the model TITLE
      CHARACTER(LEN=128)              :: STITLEI                ! the subcase SUBTITLE
      CHARACTER(LEN=128)              :: LABELI                 ! the subcase LABEL
      INTEGER(LONG)                   :: STRESS_CODE            ! flag for type of stress; see GET_STRESS_CODE

!     op2 specific flags
      INTEGER(LONG)                   :: DEVICE_CODE  ! PLOT, PRINT, PUNCH flag
      INTEGER(LONG)                   :: NUM_WIDE         ! the number of "words" for an element
      INTEGER(LONG)                   :: NVALUES          ! the number of "words" for all the elments
      INTEGER(LONG)                   :: NTOTAL           ! the number of bytes for all NVALUES
      INTEGER(LONG)                   :: ISUBCASE         ! the subcase ID
      INTEGER(LONG)                   :: NELEMENTS
      INTEGER(LONG)                   :: ISUBCASE_INDEX   ! the index into SCNUM
      INTEGER(LONG)                   :: CID              ! coordinate system
      CHARACTER(4*BYTE)               :: CEN_WORD         ! the word "CEN/" (we need to cast the length)



! **********************************************************************************************************************************
      ! Initialize
      DEVICE_CODE = 1  ! PLOT
      STRESS_CODE = 0
      FILL(1:) = ' '

      DO I=1,MAX_NUM_STR
         WRT_ERR_INDEX_NOTE(I) = 'N'
      ENDDO

      ERR_INDEX_NOTE(1) = ' (1): Polynomial fit error is in X  normal  stress'
      ERR_INDEX_NOTE(2) = ' (2): Polynomial fit error is in Y  normal  stress'
      ERR_INDEX_NOTE(3) = ' (3): Polynomial fit error is in XY shear   stress'
      ERR_INDEX_NOTE(4) = ' (4): Polynomial fit error is in X  bending stress'
      ERR_INDEX_NOTE(5) = ' (5): Polynomial fit error is in Y  bending stress'
      ERR_INDEX_NOTE(6) = ' (6): Polynomial fit error is in XY twist   stress'
      ERR_INDEX_NOTE(7) = ' (7): Polynomial fit error is in XZ shear   stress '
      ERR_INDEX_NOTE(8) = ' (8): Polynomial fit error is in YZ shear   stress'
      ERR_INDEX_NOTE(9) = ''

      FILL(1:) = ' '

      ! Get element output name
      ONAME(1:) = ' '
      CALL GET_ELEM_ONAME ( ONAME )

      ! Write output headers if this is not the first use of this subr.
      ANALYSIS_CODE = -1
      FIELD_5_INT_FLAG = .TRUE.
      FIELD5_INT_MODE = 0
      !FIELD5_FLOAT_TIME_FREQ = 0.0
      FIELD6_EIGENVALUE = 0.0
      WRITE_F06 = STRESS_USER_REQ .AND. (STRE_OUT(1:1) == 'Y')
      WRITE_GPSTRESS_F06 = GPSTRESS_REQ .AND. (GPSTRESS_OUT(1:1) == 'Y')
      INQUIRE ( UNIT=OP2, OPENED=OP2_OPENED )
      WRITE_OP2 = STRESS_USER_REQ .AND. OP2_OPENED

      IF (IHEADER == 'Y') THEN
         IF (WRITE_F06) THEN
             WRITE(F06,*)
             WRITE(F06,*)
         ENDIF
         ! -- F06 header: OUTPUT FOR SUBCASE, EIGENVECTOR or CRAIG-BAMPTON DOF
         ISUBCASE_INDEX = 0
         IF    (SOL_NAME(1:7) == 'STATICS') THEN
            ISUBCASE_INDEX = JSUB
            ANALYSIS_CODE = 1
            FIELD5_INT_MODE = 1  ! temp
            FIELD5_INT_MODE = SCNUM(JSUB)
            IF (WRITE_F06) WRITE(F06,101) SCNUM(JSUB)
         ELSE IF (SOL_NAME(1:8) == 'NLSTATIC') THEN
            ISUBCASE_INDEX = 1  ! statics
            ANALYSIS_CODE = 10
            FIELD5_INT_MODE = SCNUM(JSUB)
            IF (WRITE_F06) WRITE(F06,101) SCNUM(JSUB)

         ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1)) THEN
            ISUBCASE_INDEX = 1  ! statics
            ANALYSIS_CODE = 1
            FIELD5_INT_MODE = SCNUM(JSUB)
            IF (WRITE_F06) WRITE(F06,101) SCNUM(JSUB)

         ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
            ISUBCASE_INDEX = 2  ! modes
            ANALYSIS_CODE = 7
            FIELD5_INT_MODE = JSUB
            FIELD6_EIGENVALUE = EIGEN_VAL(JSUB)
            IF (WRITE_F06) WRITE(F06,102) JSUB

         ELSE IF (SOL_NAME(1:5) == 'MODES') THEN
            ISUBCASE_INDEX = 1  ! modes
            ANALYSIS_CODE = 2
            FIELD5_INT_MODE = JSUB
            FIELD6_EIGENVALUE = EIGEN_VAL(JSUB)
            IF (WRITE_F06) WRITE(F06,102) JSUB

         ELSE IF (SOL_NAME(1:8) == 'MFREQ') THEN
            ISUBCASE_INDEX = INT_SC_NUM
            ANALYSIS_CODE = 5
            FIELD5_INT_MODE = JSUB
            FIELD6_EIGENVALUE = EIGEN_VAL(JSUB)
            IF (WRITE_F06) WRITE(F06,102) JSUB

         ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
            ISUBCASE_INDEX = 1  ! modes
            IF ((JSUB <= NDOFR) .OR. (JSUB >= NDOFR+NVEC)) THEN
               IF (JSUB <= NDOFR) THEN
                  BDY_DOF_NUM = JSUB
               ELSE
                  BDY_DOF_NUM = JSUB-(NDOFR+NVEC)
               ENDIF
               CALL GET_GRID_AND_COMP ( 'R ', BDY_DOF_NUM, BDY_GRID, BDY_COMP  )
            ENDIF

            IF (WRITE_F06) THEN
               IF (JSUB <= NDOFR) THEN
                   WRITE(F06,103) JSUB, NUM_CB_DOFS, 'acceleration', BDY_GRID, BDY_COMP
               ELSE IF ((JSUB > NDOFR) .AND. (JSUB <= NDOFR+NVEC)) THEN
                   WRITE(F06,104) JSUB, NUM_CB_DOFS, JSUB-NDOFR
               ELSE
                   WRITE(F06,103) JSUB, NUM_CB_DOFS, 'displacement', BDY_GRID, BDY_COMP
               ENDIF
            ENDIF  ! write f06

         ENDIF
         ISUBCASE = SCNUM(ISUBCASE_INDEX)

         ! -- F06 header for TITLE, SUBTITLE, LABEL (but only to F06)
         TITLEI = TITLE(INT_SC_NUM)
         STITLEI = STITLE(INT_SC_NUM)
         LABELI = LABEL(INT_SC_NUM)

         IF (WRITE_F06) THEN
            IF (TITLE(INT_SC_NUM)(1:)  /= ' ') THEN
               WRITE(F06,201) TITLE(INT_SC_NUM)
            ENDIF

            IF (STITLE(INT_SC_NUM)(1:) /= ' ') THEN
               WRITE(F06,201) STITLE(INT_SC_NUM)
            ENDIF

            IF (LABEL(INT_SC_NUM)(1:)  /= ' ') THEN
               WRITE(F06,201) LABEL(INT_SC_NUM)
            ENDIF
            WRITE(F06,*)

           ! -- F06 1st 2 header lines for stress output description
            IF     ((TYPE(1:3) == 'BAR') .OR. (TYPE(1:4) == 'BEAM')) THEN
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  WRITE(F06,302) FILL(1: 20)
               ELSE
                  WRITE(F06,301) FILL(1: 13)
               ENDIF
               WRITE(F06,401) FILL(1: 42), ONAME

            ELSE IF (TYPE(1:4) == 'BUSH') THEN
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  WRITE(F06,302) FILL(1: 20)
               ELSE
                  WRITE(F06,301) FILL(1: 11)
               ENDIF
               WRITE(F06,401) FILL(1: 40), ONAME

            ELSE IF (TYPE(1:4) == 'ELAS') THEN
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  WRITE(F06,302) FILL(1: 20)
               ELSE
                  WRITE(F06,301) FILL(1: 11)
               ENDIF
               WRITE(F06,401) FILL(1: 40), ONAME

            ELSE IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PYRAM') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
               IF (STRE_OPT == 'VONMISES') THEN
                  IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                     IF(STR_CID == -2) THEN
                        WRITE(F06,312) FILL(1: 20)
                     ELSE
                        WRITE(F06,302) FILL(1: 15)
                     ENDIF
                  ELSE
                     IF(STR_CID == -2) THEN
                        WRITE(F06,311) FILL(1: 32)
                     ELSE
                        WRITE(F06,301) FILL(1: 27)
                     ENDIF
                  ENDIF
                  WRITE(F06,401) FILL(1: 55), ONAME
               ELSE
                  IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                     IF(STR_CID == -2) THEN
                        WRITE(F06,312) FILL(1: 27)
                     ELSE
                        WRITE(F06,302) FILL(1: 22)
                     ENDIF
                  ELSE
                     IF(STR_CID == -2) THEN
                        WRITE(F06,311) FILL(1: 38)
                     ELSE
                        WRITE(F06,301) FILL(1: 33)
                     ENDIF
                  ENDIF
                  WRITE(F06,401) FILL(1: 61), ONAME
               ENDIF

! --- CQUADR_DKMQ24 begin --- !
            ELSE IF (((TYPE(1:5) == 'QUAD4') .OR. (TYPE == 'QUADR   ')) .OR. (TYPE(1:5) == 'QUAD8')) THEN
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  WRITE(F06,302) FILL(1: 20)
               ELSE
                  IF (STR_CID == 0) THEN
                     WRITE(F06,321) FILL(1: 45)
                  ELSE IF (STR_CID > 0) THEN
                     WRITE(F06,331) FILL(1: 41), STR_CID
                  ELSE
                     WRITE(F06,301) FILL(1: 42)
                  ENDIF
               ENDIF
               WRITE(F06,401) FILL(1: 71), ONAME

            ELSE IF (TYPE(1:3) == 'ROD') THEN
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  WRITE(F06,302) FILL(1: 20)
               ELSE
                  WRITE(F06,301) FILL(1: 13)
               ENDIF
               WRITE(F06,401) FILL(1: 42), ONAME

            ELSE IF (TYPE(1:5) == 'SHEAR') THEN
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  WRITE(F06,302) FILL(1: 20)
               ELSE
                  WRITE(F06,301) FILL(1: 13)
               ENDIF
               WRITE(F06,401) FILL(1: 42), ONAME

            ELSE IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'TRIA6')) THEN
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  WRITE(F06,302) FILL(1: 20)
               ELSE
                  IF (STR_CID == 0) THEN
                     WRITE(F06,321) FILL(1: 39)
                  ELSE IF (STR_CID > 0) THEN
                     WRITE(F06,331) FILL(1: 35), STR_CID
                  ELSE
                     WRITE(F06,301) FILL(1: 36)
                  ENDIF
               ENDIF
               WRITE(F06,401) FILL(1: 65), ONAME
            ENDIF

            ! -- F06 header lines describing stress columns
            IF      (TYPE == 'BAR     ') THEN
               IF (BARTOR == 'Y') THEN
                  WRITE(F06,1101) FILL(1:1), FILL(1:1)
               ELSE
                  WRITE(F06,1102) FILL(1:1), FILL(1:1)
               ENDIF

            ELSE IF (TYPE == 'BEAM    ') THEN
               IF (BARTOR == 'Y') THEN
                  WRITE(F06,1104) FILL(1:1), FILL(1:1)
               ELSE
                  WRITE(F06,1105) FILL(1:1), FILL(1:1)
               ENDIF

            ELSE IF (TYPE(1:4) == 'ELAS') THEN
               WRITE(F06,1201) FILL(1:1), FILL(1:1)

            ELSE IF((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PYRAM') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
               IF (STRE_OPT == 'VONMISES') THEN
                  WRITE(F06,1301) FILL(1: 1), FILL(1: 1)
               ELSE
                  WRITE(F06,1302) FILL(1: 1), FILL(1: 1)
               ENDIF

            ELSE IF (((TYPE(1:5) == 'QUAD4') .OR. (TYPE == 'QUADR   ')) .OR. (TYPE(1:5) == 'QUAD8')) THEN
               IF (STRE_OPT == 'VONMISES') THEN
                  WRITE(F06,1401) FILL(1: 1), FILL(1: 1), FILL(1: 1)
               ELSE
                  WRITE(F06,1402) FILL(1: 1), FILL(1: 1)
               ENDIF

            ELSE IF  (TYPE == 'ROD     ') THEN
               WRITE(F06,1501) FILL(1: 1), FILL(1: 1)

            ELSE IF (TYPE(1:5) == 'SHEAR') THEN
               WRITE(F06,1601) FILL(1: 1), FILL(1: 1)

            ELSE IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'TRIA6')) THEN
               IF (STRE_OPT == 'VONMISES') THEN
                  WRITE(F06,1701) FILL(1: 1), FILL(1: 1), FILL(1: 1)
               ELSE
                  WRITE(F06,1702) FILL(1: 1), FILL(1: 1)
               ENDIF

            ELSE IF  (TYPE == 'BUSH    ') THEN
               WRITE(F06,1801) FILL(1: 1), FILL(1: 1)

            ELSE IF  (TYPE == 'USERIN  ') THEN
               WRITE(F06,1901) FILL(1: 1), FILL(1: 1)

            ENDIF


         ENDIF  ! write f06






      ENDIF

      ! Write the element stress output
! --- cbeam_stations begin --- !
      IF      (TYPE == 'BAR     ') THEN
         CALL WRITE_BAR(NUM, FILL(1:1), ISUBCASE, ITABLE, TITLEI, STITLEI, LABELI, &
                        FIELD5_INT_MODE, FIELD6_EIGENVALUE, WRITE_F06)

      ELSE IF (TYPE == 'BEAM    ') THEN
! --- CBEAM_standard begin --- !
         IF (WRITE_OP2) THEN
            NELEMENTS = 0
            I = 1
            DO WHILE (I <= NUM)
               NELEMENTS = NELEMENTS + 1
               J = EID_OUT_ARRAY(I,1)
               DO WHILE ((I <= NUM) .AND. (EID_OUT_ARRAY(I,1) == J))
                  I = I + 1
               ENDDO
            ENDDO

            ALLOCATE ( BEAM_EID(NELEMENTS), BEAM_GRID(NELEMENTS,11), BEAM_XI(NELEMENTS,11), BEAM_SXC(NELEMENTS,11),               &
                       BEAM_SXD(NELEMENTS,11), BEAM_SXE(NELEMENTS,11), BEAM_SXF(NELEMENTS,11), BEAM_SMAX(NELEMENTS,11),            &
                       BEAM_SMIN(NELEMENTS,11), BEAM_MST(NELEMENTS,11), BEAM_MSC(NELEMENTS,11) )

            BEAM_GRID(:,:) = 0
            BEAM_XI(:,:)   = 0.0D0
            BEAM_SXC(:,:)  = 0.0D0
            BEAM_SXD(:,:)  = 0.0D0
            BEAM_SXE(:,:)  = 0.0D0
            BEAM_SXF(:,:)  = 0.0D0
            BEAM_SMAX(:,:) = 0.0D0
            BEAM_SMIN(:,:) = 0.0D0
            BEAM_MST(:,:)  = 0.0D0
            BEAM_MSC(:,:)  = 0.0D0

            I = 1
            IELEM = 0
            DO WHILE (I <= NUM)
               IELEM = IELEM + 1
               BEAM_EID(IELEM) = EID_OUT_ARRAY(I,1)
               IBEG = I
               DO WHILE ((I <= NUM) .AND. (EID_OUT_ARRAY(I,1) == BEAM_EID(IELEM)))
                  I = I + 1
               ENDDO
               IEND = I - 1
               NROW_ELEM = IEND - IBEG + 1
               NSTA_ELEM = NROW_ELEM
               IF (NSTA_ELEM > 11) NSTA_ELEM = 11

! --- CBEAM_standard begin --- !
               ! Populate station grid ids for OP2 CBEAM stress records.
               ! Keep interior stations at zero and identify only the beam
               ! endpoint stations with the element end grids.
               BEAM_GRID(IELEM,:) = 0
               BEAM_GRID(IELEM,1)  = GID_OUT_ARRAY(IBEG,2)
               BEAM_GRID(IELEM,11) = GID_OUT_ARRAY(IBEG,3)
! --- CBEAM_standard end --- !

               DO ISTA=1,NSTA_ELEM
! --- CBEAM_standard begin --- !
                  ! For beam stress, EID_OUT_ARRAY/CBEAM_XL_OUT are stored once per
                  ! station, but OGEL stores two rows per station (top/bottom style
                  ! beam section stress output). Convert the station index into the
                  ! corresponding OGEL row pair explicitly.
                  K = 2*(IBEG + ISTA - 2) + 1
                  XI_RAW  (ISTA) = CBEAM_XL_OUT(IBEG + ISTA - 1)
                  SXC_RAW (ISTA) = OGEL(K    ,1)
                  SXD_RAW (ISTA) = OGEL(K    ,2)
                  SXE_RAW (ISTA) = OGEL(K    ,3)
                  SXF_RAW (ISTA) = OGEL(K    ,4)
                  SMAX_RAW(ISTA) = OGEL(K    ,6)
                  SMIN_RAW(ISTA) = OGEL(K    ,7)
                  MST_RAW (ISTA) = OGEL(K    ,8)
                  MSC_RAW (ISTA) = OGEL(K + 1,8)
                  IF (MST_RAW(ISTA) <= -0.999D0) MST_RAW(ISTA) = 0.0D0
                  IF (MSC_RAW(ISTA) <= -0.999D0) MSC_RAW(ISTA) = 0.0D0
! --- CBEAM_standard end --- !
               ENDDO

               DO ISTA=1,11
                  XI_STD = DBLE(ISTA - 1)/10.0D0
                  BEAM_XI(IELEM,ISTA) = XI_STD
                  IF (NSTA_ELEM <= 1) THEN
                     BEAM_SXC (IELEM,ISTA) = SXC_RAW (1)
                     BEAM_SXD (IELEM,ISTA) = SXD_RAW (1)
                     BEAM_SXE (IELEM,ISTA) = SXE_RAW (1)
                     BEAM_SXF (IELEM,ISTA) = SXF_RAW (1)
                     BEAM_SMAX(IELEM,ISTA) = SMAX_RAW(1)
                     BEAM_SMIN(IELEM,ISTA) = SMIN_RAW(1)
                     BEAM_MST (IELEM,ISTA) = MST_RAW (1)
                     BEAM_MSC (IELEM,ISTA) = MSC_RAW (1)
                  ELSE IF (XI_STD <= XI_RAW(1)) THEN
                     BEAM_SXC (IELEM,ISTA) = SXC_RAW (1)
                     BEAM_SXD (IELEM,ISTA) = SXD_RAW (1)
                     BEAM_SXE (IELEM,ISTA) = SXE_RAW (1)
                     BEAM_SXF (IELEM,ISTA) = SXF_RAW (1)
                     BEAM_SMAX(IELEM,ISTA) = SMAX_RAW(1)
                     BEAM_SMIN(IELEM,ISTA) = SMIN_RAW(1)
                     BEAM_MST (IELEM,ISTA) = MST_RAW (1)
                     BEAM_MSC (IELEM,ISTA) = MSC_RAW (1)
                  ELSE IF (XI_STD >= XI_RAW(NSTA_ELEM)) THEN
                     BEAM_SXC (IELEM,ISTA) = SXC_RAW (NSTA_ELEM)
                     BEAM_SXD (IELEM,ISTA) = SXD_RAW (NSTA_ELEM)
                     BEAM_SXE (IELEM,ISTA) = SXE_RAW (NSTA_ELEM)
                     BEAM_SXF (IELEM,ISTA) = SXF_RAW (NSTA_ELEM)
                     BEAM_SMAX(IELEM,ISTA) = SMAX_RAW(NSTA_ELEM)
                     BEAM_SMIN(IELEM,ISTA) = SMIN_RAW(NSTA_ELEM)
                     BEAM_MST (IELEM,ISTA) = MST_RAW (NSTA_ELEM)
                     BEAM_MSC (IELEM,ISTA) = MSC_RAW (NSTA_ELEM)
                  ELSE
                     DO K=1,NSTA_ELEM-1
                        XI0 = XI_RAW(K)
                        XI1 = XI_RAW(K+1)
                        IF ((XI_STD >= XI0) .AND. (XI_STD <= XI1)) THEN
                           TINT = (XI_STD - XI0)/(XI1 - XI0)
                           BEAM_SXC (IELEM,ISTA) = (1.0D0 - TINT)*SXC_RAW (K) + TINT*SXC_RAW (K+1)
                           BEAM_SXD (IELEM,ISTA) = (1.0D0 - TINT)*SXD_RAW (K) + TINT*SXD_RAW (K+1)
                           BEAM_SXE (IELEM,ISTA) = (1.0D0 - TINT)*SXE_RAW (K) + TINT*SXE_RAW (K+1)
                           BEAM_SXF (IELEM,ISTA) = (1.0D0 - TINT)*SXF_RAW (K) + TINT*SXF_RAW (K+1)
                           BEAM_SMAX(IELEM,ISTA) = (1.0D0 - TINT)*SMAX_RAW(K) + TINT*SMAX_RAW(K+1)
                           BEAM_SMIN(IELEM,ISTA) = (1.0D0 - TINT)*SMIN_RAW(K) + TINT*SMIN_RAW(K+1)
                           BEAM_MST (IELEM,ISTA) = (1.0D0 - TINT)*MST_RAW (K) + TINT*MST_RAW (K+1)
                           BEAM_MSC (IELEM,ISTA) = (1.0D0 - TINT)*MSC_RAW (K) + TINT*MSC_RAW (K+1)
                           EXIT
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO

            ELEMENT_TYPE = 2
            NUM_WIDE = 111
            NVALUES = NELEMENTS*NUM_WIDE
            CALL GET_STRESS_CODE( STRESS_CODE, 1, 0, 0 )
            CALL WRITE_OES3_STATIC_AC(ITABLE, ISUBCASE, DEVICE_CODE, ANALYSIS_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                    TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
            WRITE(OP2) NVALUES
            WRITE(OP2) (BEAM_EID(IELEM)*10+DEVICE_CODE,                                                                           &
                        (BEAM_GRID(IELEM,ISTA), REAL(BEAM_XI(IELEM,ISTA),4), REAL(BEAM_SXC(IELEM,ISTA),4),                      &
                          REAL(BEAM_SXD(IELEM,ISTA),4), REAL(BEAM_SXE(IELEM,ISTA),4), REAL(BEAM_SXF(IELEM,ISTA),4),            &
                          REAL(BEAM_SMAX(IELEM,ISTA),4), REAL(BEAM_SMIN(IELEM,ISTA),4),                                          &
                          REAL(BEAM_MST(IELEM,ISTA),4), REAL(BEAM_MSC(IELEM,ISTA),4), ISTA=1,11),                               &
                        IELEM=1,NELEMENTS)
            DEALLOCATE ( BEAM_EID, BEAM_GRID, BEAM_XI, BEAM_SXC, BEAM_SXD, BEAM_SXE, BEAM_SXF, BEAM_SMAX, BEAM_SMIN,              &
                         BEAM_MST, BEAM_MSC )
         ENDIF
! --- CBEAM_standard end --- !
         CALL WRITE_CBEAM_STRESS(NUM, WRITE_F06)
! --- cbeam_stations end --- !

      ELSE IF (TYPE(1:4) == 'ELAS') THEN
         IF (WRITE_OP2) THEN
             CALL GET_SPRING_OP2_ELEMENT_TYPE(ELEMENT_TYPE)

             NUM_WIDE = 2 ! eid, spring_stress
             NVALUES = NUM_WIDE * NUM

             DEVICE_CODE = 1   ! PLOT

             !CALL GET_STRESS_CODE(STRESS_CODE, IS_VON_MISES, IS_STRAIN, IS_FIBER_DISTANCE)
             CALL GET_STRESS_CODE( STRESS_CODE, 1,            0,         1)
             CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                    TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)

             WRITE(OP2) NVALUES
             WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, REAL(OGEL(I,1), 4), I=1,NUM)
         ENDIF   ! end of op2

         IF (WRITE_F06) THEN
            DO I=1,NUM,5
               CALL WRITE_STRESS_ELAS_GROUP_LINE ( I, MIN(I+4,NUM) )
            ENDDO
         ENDIF

      ELSE IF((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PYRAM') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
         !       12345
         ! 39  : CTETRA
         ! 67  : CHEXA
         ! 68  : CPENTA
         ! 255 : PYRAM
         IF (TYPE(1:4) == "HEXA") THEN
             ELEMENT_TYPE = 67
             NNODES = 9
         ELSE IF (TYPE(1:5) == "TETRA") THEN
             ELEMENT_TYPE = 39
             NNODES = 5
         ELSE IF (TYPE(1:5) == "PENTA") THEN
             ELEMENT_TYPE = 68
             NNODES = 7
         ELSE IF (TYPE(1:5) == "PYRAM") THEN
             ELEMENT_TYPE = 255
             NNODES = 6
         ENDIF
         NUM_WIDE = 4 + 21*NNODES
         NELEMENTS = NUM / NUM_PTS
         NVALUES = NUM_WIDE * NELEMENTS

         IF (WRITE_OP2) THEN
           !CALL GET_STRESS_CODE(STRESS_CODE, IS_VON_MISES, IS_STRAIN, IS_FIBER_DISTANCE)
           CALL GET_STRESS_CODE( STRESS_CODE, 1,            0,         0)
           CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                  TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
           WRITE(OP2) NVALUES
           CEN_WORD = "CEN/"

          ! See the CHEXA, CPENTA, or CTETRA entry for the definition of the element coordinate systems.
          ! The material coordinate system (CORDM) may be the basic system (0 or blank), any defined system
          ! (Integer > 0), or the standard internal coordinate system of the element designated as:
          ! -1: element coordinate system (-1)
          ! -2: element system based on eigenvalue techniques to insure non bias in the element formulation.

          ! TODO hardcoded
           CID = -1

          ! setting:
          !  - CTETRA: [element_device, cid, 'CEN/', 4]
          !  - PYRAM: [element_device, cid, 'CEN/', 5]
          !  - CPENTA: [element_device, cid, 'CEN/', 6]
          !  - CHEXA:  [element_device, cid, 'CEN/', 8]

           !                 1             2             3            4            5               6             7
           !  Element    Sigma-xx      Sigma-yy      Sigma-zz       Tau-xy        Tau-yz        Tau-zx      von Mises
           !     ID

           WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, CID, CEN_WORD, NNODES-1,                                       &
                       (GID_OUT_ARRAY(I,J),                                                                              &
                        REAL(OGEL(I+J-1,1),4), REAL(OGEL(I+J-1,4),4), REAL(OGEL(I+J-1,9), 4), 0., 0., 0.,              &
                        REAL(OGEL(I+J-1,12),4), REAL(OGEL(I+J-1,7),4),                                                    &
                        REAL(OGEL(I+J-1,2),4), REAL(OGEL(I+J-1,5),4), REAL(OGEL(I+J-1,10),4), 0., 0., 0.,              &
                        REAL(OGEL(I+J-1,3),4), REAL(OGEL(I+J-1,6),4), REAL(OGEL(I+J-1,11),4), 0., 0., 0.,              &
                        J=1,NNODES), I=1,NUM,NUM_PTS)
         ENDIF  ! end of op2

         IF (STRE_OPT == 'VONMISES') THEN
            NCOLS = 7
         ELSE
            NCOLS = 8
         ENDIF

         IF (WRITE_F06) THEN
            K = 0
            DO I=1,NUM,NUM_PTS
               K = K + 1
               ! Center
               CALL WRITE_STRESS_SOLID_CENTER_LINE ( EID_OUT_ARRAY(I,1), OGEL(K,1:NCOLS), NCOLS )
               ! Corner
               DO L=1,NUM_PTS-1
                  K = K + 1
                  CALL WRITE_STRESS_SOLID_GRID_LINE ( GID_OUT_ARRAY(I,L+1), OGEL(K,1:NCOLS), NCOLS )
               ENDDO
            ENDDO
         ENDIF

         CALL GET_MAX_MIN_ABS_STR ( NUM, NCOLS, 'N', MAX_ANS, MIN_ANS, ABS_ANS )

         IF (WRITE_F06) THEN
            IF (STRE_OPT == 'VONMISES') THEN
               WRITE(F06,1304) (MAX_ANS(J),J=1,7), (MIN_ANS(J),J=1,7), (ABS_ANS(J),J=1,7)
            ELSE
               WRITE(F06,1305) (MAX_ANS(J),J=1,8), (MIN_ANS(J),J=1,8), (ABS_ANS(J),J=1,8)
            ENDIF
         ENDIF

      ELSE IF (((TYPE(1:5) == 'QUAD4') .OR. (TYPE == 'QUADR   ')) .OR. (TYPE(1:5) == 'QUAD8')) THEN
         !CALL WRITE_OES_CQUAD4 ( NUM, FILL, ISUBCASE, ITABLE, TITLEI, STITLEI, LABELI )

         IF (WRITE_OP2) THEN
           IF ((STRE_LOC /= 'CORNER  ') .AND. (TYPE(1:5) /= 'QUAD8')) THEN
              CALL GET_STRESS_CODE( STRESS_CODE, 1,            0,         0)
              ! CQUAD4-33
               !(eid_device,
              ! fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
              ! fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out; n=17
              NUM_WIDE = 17
              ELEMENT_TYPE = 33
              NELEMENTS = NUM / MAX(1_LONG,NUM_PTS)
              NVALUES = NUM_WIDE * NELEMENTS
              CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                     TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
              !NUM_PTS = 1
              ! just a copy of the CTRIA3 code
              ! op2 version of the upper & lower layers all in one call, but without the transverse shear
              WRITE(OP2) NVALUES
              WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, (REAL(OGEL(2*I-1,J),4), J=1,8), (REAL(OGEL(2*I,J),4), J=1,8), &
                          I=1,NUM,MAX(1_LONG,NUM_PTS))
           ELSE
              CALL GET_STRESS_CODE( STRESS_CODE, 1,            0,         1)
              ! CQUAD4-144 / CQUAD8-64
              IF (TYPE(1:5) == 'QUAD8') THEN
                 ELEMENT_TYPE = 64
              ELSE
                 ELEMENT_TYPE = 144
              ENDIF
              NUM_WIDE = 87 ! 2 + 17 * (4+1)  ! 4 nodes + 1 centroid

              ! TODO: probably wrong...divide NUM by NUM_PTS?
              NELEMENTS = NUM / NUM_PTS
              NVALUES = NUM_WIDE * NELEMENTS
              ! NUM=  10 NUM_PTS=   5
              !(eid_device, "CEN/", 4, # "CEN/4"
              ! fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
              ! fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = n = 17+2
              !
              ! (grid,
              !  fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
              !  fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)*4 = n = 17*4
              CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                     TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
              WRITE(OP2) NVALUES
              ! see the CQUAD4-33 stress/strain (the IF part of this IF-ELSE block)
              ! writing before trying to understand this...
              !
              ! basically a one-liner version of the F06 writing
              ! we broke out the L=1,NUM_PTS-1 loop to 4 lines (the GID_OUT_ARRAY lines)
              ! to avoid an additional hard to write loop
              WRITE(OP2) (EID_OUT_ARRAY(5*I+1,1)*10+DEVICE_CODE, "CEN/", 4,                                           &
                                                  (REAL(OGEL(10*I+1,J),4), J=1,8), (REAL(OGEL(10*I+2,  J),4), J=1,8), &
                          GID_OUT_ARRAY(5*I+1,2), (REAL(OGEL(10*I+3,J),4), J=1,8), (REAL(OGEL(10*I+4,  J),4), J=1,8), &
                          GID_OUT_ARRAY(5*I+1,3), (REAL(OGEL(10*I+5,J),4), J=1,8), (REAL(OGEL(10*I+6,  J),4), J=1,8), &
                          GID_OUT_ARRAY(5*I+1,4), (REAL(OGEL(10*I+7,J),4), J=1,8), (REAL(OGEL(10*I+8,  J),4), J=1,8), &
                          GID_OUT_ARRAY(5*I+1,5), (REAL(OGEL(10*I+9,J),4), J=1,8), (REAL(OGEL(10*(I+1),J),4), J=1,8), &
                          I=0,NELEMENTS-1)
           ENDIF
         ENDIF  ! end of op2

         K = 0
         DO I=1,NUM,NUM_PTS
            K = K + 1
            IF (WRITE_F06) WRITE(F06,*)
            IF (WRITE_F06) THEN
               QUAD_VALUES_10(1:10) = OGEL(K,1:10)
               CALL TRANSFORM_SHELL_OUTPUT_ROW_TO_BASIC ( SHELL_OUT_TE(1:3,1:3,K), QUAD_VALUES_10 )
               CALL FAST_BUILD_QUAD_1403_LINE ( EID_OUT_ARRAY(I,1), QUAD_VALUES_10, QUAD_CENTER_LINE )
               WRITE(F06,'(A)') QUAD_CENTER_LINE
            ENDIF
            K = K + 1
            IF (WRITE_F06) THEN
               QUAD_VALUES_10(1:10) = OGEL(K,1:10)
               CALL TRANSFORM_SHELL_OUTPUT_ROW_TO_BASIC ( SHELL_OUT_TE(1:3,1:3,K), QUAD_VALUES_10 )
               QUAD_VALUES_8(1:8) = QUAD_VALUES_10(1:8)
               CALL FAST_BUILD_QUAD_1404_LINE ( QUAD_VALUES_8, QUAD_LOWER_LINE )
               WRITE(F06,'(A)') QUAD_LOWER_LINE
            ENDIF

            IF ((STRE_LOC == 'CORNER  ') .OR. (TYPE(1:5) == 'QUAD8')) THEN
               DO L=1,NUM_PTS-1
                  K = K + 1
                  IF (WRITE_F06) WRITE(F06,*)
                  IF (DABS(POLY_FIT_ERR(I+L)) >= 0.01D0) THEN
                     IF (WRITE_F06) THEN
                        QUAD_VALUES_10(1:10) = OGEL(K,1:10)
                        CALL TRANSFORM_SHELL_OUTPUT_ROW_TO_BASIC ( SHELL_OUT_TE(1:3,1:3,K), QUAD_VALUES_10 )
                        CALL FAST_BUILD_QUAD_1405_LINE ( GID_OUT_ARRAY(I,L+1), QUAD_VALUES_10, POLY_FIT_ERR(I+L),       &
                                                         POLY_FIT_ERR_INDEX(I+L), QUAD_GRID_NOTE_LINE )
                        WRITE(F06,'(A)') QUAD_GRID_NOTE_LINE
                     ENDIF
                     WRT_ERR_INDEX_NOTE(POLY_FIT_ERR_INDEX(I+L)) = 'Y'
                  ELSE
                     IF (WRITE_F06) THEN
                        QUAD_VALUES_10(1:10) = OGEL(K,1:10)
                        CALL TRANSFORM_SHELL_OUTPUT_ROW_TO_BASIC ( SHELL_OUT_TE(1:3,1:3,K), QUAD_VALUES_10 )
                        CALL FAST_BUILD_QUAD_1406_LINE ( GID_OUT_ARRAY(I,L+1), QUAD_VALUES_10, POLY_FIT_ERR(I+L),       &
                                                         QUAD_GRID_LINE )
                        WRITE(F06,'(A)') QUAD_GRID_LINE
                     ENDIF
                  ENDIF

                  K = K + 1
                  IF (WRITE_F06) THEN
                     QUAD_VALUES_10(1:10) = OGEL(K,1:10)
                     CALL TRANSFORM_SHELL_OUTPUT_ROW_TO_BASIC ( SHELL_OUT_TE(1:3,1:3,K), QUAD_VALUES_10 )
                     QUAD_VALUES_8(1:8) = QUAD_VALUES_10(1:8)
                     CALL FAST_BUILD_QUAD_1404_LINE ( QUAD_VALUES_8, QUAD_LOWER_LINE )
                     WRITE(F06,'(A)') QUAD_LOWER_LINE
                  ENDIF

               ENDDO
            ELSE
               K = K + 2*(NUM_PTS - 1)
            ENDIF
         ENDDO

         MAX_ANS(1:10) = -HUGE(1.0D0)
         MIN_ANS(1:10) =  HUGE(1.0D0)
         ABS_ANS(1:10) = ZERO
         IF ((STRE_LOC == 'CORNER  ') .OR. (TYPE(1:5) == 'QUAD8')) THEN
            DO K=1,2*NUM
               QUAD_VALUES_10(1:10) = OGEL(K,1:10)
               CALL TRANSFORM_SHELL_OUTPUT_ROW_TO_BASIC ( SHELL_OUT_TE(1:3,1:3,K), QUAD_VALUES_10 )
               DO J=2,10
                  IF (QUAD_VALUES_10(J) > MAX_ANS(J)) MAX_ANS(J) = QUAD_VALUES_10(J)
                  IF (QUAD_VALUES_10(J) < MIN_ANS(J)) MIN_ANS(J) = QUAD_VALUES_10(J)
                  ABS_ANS(J) = MAX( ABS_ANS(J), DABS(QUAD_VALUES_10(J)) )
               ENDDO
            ENDDO
         ELSE
            DO I=1,NUM,MAX(1_LONG,NUM_PTS)
               DO L=0,1
                  K = 2*I + L - 1
                  QUAD_VALUES_10(1:10) = OGEL(K,1:10)
                  CALL TRANSFORM_SHELL_OUTPUT_ROW_TO_BASIC ( SHELL_OUT_TE(1:3,1:3,K), QUAD_VALUES_10 )
                  DO J=2,10
                     IF (QUAD_VALUES_10(J) > MAX_ANS(J)) MAX_ANS(J) = QUAD_VALUES_10(J)
                     IF (QUAD_VALUES_10(J) < MIN_ANS(J)) MIN_ANS(J) = QUAD_VALUES_10(J)
                     ABS_ANS(J) = MAX( ABS_ANS(J), DABS(QUAD_VALUES_10(J)) )
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         MAX_ANS(1) = ZERO
         MIN_ANS(1) = ZERO
         ABS_ANS(1) = ZERO

             ! Get max POLY_FIT_ERR
         MAX_ANS(11) = ZERO
         K = 0
         DO I=1,NUM
            K = K + 1
            IF (POLY_FIT_ERR(I) > MAX_ANS(11)) THEN
               MAX_ANS(11) = POLY_FIT_ERR(I)
            ENDIF
            K = K + 1
         ENDDO

         MIN_ANS(11) = MAX_ANS(11)

             ! Get min POLY_FIT_ERR
         K = 0
         DO I=1,NUM
            K = K + 1
            IF (POLY_FIT_ERR(I) < MIN_ANS(11)) THEN
               MIN_ANS(11) = POLY_FIT_ERR(I)
            ENDIF
            K = K + 1
         ENDDO
             ! Get abs POLY_FIT_ERR
         ABS_ANS(11) = MAX( DABS(MAX_ANS(11)), DABS(MIN_ANS(11)) )

         IF (WRITE_F06) THEN
            IF ((STRE_LOC == 'CORNER  ') .OR. (TYPE(1:5) == 'QUAD8')) THEN
               WRITE(F06,1408) FILL(1: 0),                                                                                      &
                               FILL(1: 0),       MAX_ANS(2),MAX_ANS(3),MAX_ANS(4),MAX_ANS(6),MAX_ANS(7),MAX_ANS(8),MAX_ANS(9),  &
                                                 MAX_ANS(10),MAX_ANS(11),                                                       &
                               FILL(1: 0),       MIN_ANS(2),MIN_ANS(3),MIN_ANS(4),MIN_ANS(6),MIN_ANS(7),MIN_ANS(8),MIN_ANS(9),  &
                                                 MIN_ANS(10),MIN_ANS(11),                                                       &
                               FILL(1: 0),       ABS_ANS(2),ABS_ANS(3),ABS_ANS(4),ABS_ANS(6),ABS_ANS(7),ABS_ANS(8),ABS_ANS(9),  &
                                                 ABS_ANS(10),ABS_ANS(11), FILL(1: 0)
            ELSE
               WRITE(F06,1408) FILL(1: 0),                                                                                      &
                               FILL(1: 0),       MAX_ANS(2),MAX_ANS(3),MAX_ANS(4),MAX_ANS(6),MAX_ANS(7),MAX_ANS(8),MAX_ANS(9),  &
                                                 MAX_ANS(10),MAX_ANS(11),                                                       &
                               FILL(1: 0),       MIN_ANS(2),MIN_ANS(3),MIN_ANS(4),MIN_ANS(6),MIN_ANS(7),MIN_ANS(8),MIN_ANS(9),  &
                                                 MIN_ANS(10),MIN_ANS(11),                                                       &
                               FILL(1: 0),       ABS_ANS(2),ABS_ANS(3),ABS_ANS(4),ABS_ANS(6),ABS_ANS(7),ABS_ANS(8),ABS_ANS(9),  &
                                                 ABS_ANS(10),ABS_ANS(11), FILL(1: 0)
            ENDIF
         ENDIF

         WRITE_NOTES = 'N'
         DO I=1,MAX_NUM_STR
            IF (WRT_ERR_INDEX_NOTE(I) == 'Y') THEN
               WRITE_NOTES = 'Y'
            ENDIF
         ENDDO

         IF ((WRITE_NOTES == 'Y') .AND. (WRITE_F06)) THEN
            WRITE(F06,1498)
            DO I=1,MAX_NUM_STR
               IF (WRT_ERR_INDEX_NOTE(I) == 'Y') THEN
                  WRITE(F06,1499) ERR_INDEX_NOTE(I)
               ENDIF
            ENDDO
         ENDIF

         IF ((OP2_OPENED .OR. WRITE_GPSTRESS_F06) .AND. GPSTRESS_REQ .AND. (NUM_PTS > 1)) THEN
            CALL WRITE_OGS1_SURFACE_STRESS ( ITABLE, ISUBCASE, NUM, NUM_PTS, DEVICE_CODE, ANALYSIS_CODE, FIELD5_INT_MODE,       &
                                             FIELD6_EIGENVALUE, TITLEI, STITLEI, LABELI, 'QUAD', WRITE_GPSTRESS_F06, OP2_OPENED )
         ENDIF

      ELSE IF (TYPE == 'ROD     ') THEN
         CALL WRITE_ROD (ISUBCASE, NUM, FILL(1:1), ITABLE, TITLEI, STITLEI, LABELI,  &
                         FIELD5_INT_MODE, FIELD6_EIGENVALUE, WRITE_OP2 )

      ELSE IF (TYPE(1:5) == 'SHEAR') THEN
         CALL WRITE_OES_CSHEAR(NUM, FILL, ISUBCASE, ITABLE, TITLEI, STITLEI, LABELI, &
                               FIELD5_INT_MODE, FIELD6_EIGENVALUE,                   &
                               WRITE_F06, WRITE_OP2)

      ELSE IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'TRIA6')) THEN
         CALL WRITE_OES_CTRIA3(NUM, NUM_PTS, FILL, ISUBCASE, ITABLE, TITLEI, STITLEI, LABELI, &
                               FIELD5_INT_MODE, FIELD6_EIGENVALUE,                   &
                               WRITE_F06, WRITE_OP2, OP2_OPENED, WRITE_GPSTRESS_F06)

      ELSE IF (TYPE == 'BUSH    ') THEN
         IF (WRITE_OP2) THEN
             ELEMENT_TYPE = 102 ! CBUSH
             NUM_WIDE = 7       ! eid, tx, ty, tz, rx, ry, rz
             STRESS_CODE = 1    ! dunno
             NVALUES = NUM * NUM_WIDE

             CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                    TITLEI, STITLEI, LABELI, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
             WRITE(OP2) NVALUES
             WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE,(REAL(OGEL(I,J),4),J=1,6), I=1,NUM)
         ENDIF

         IF (WRITE_F06) THEN
            DO I=1,NUM
               CALL WRITE_STRESS_I8_PLUS_R14_LINE ( 1_LONG, EID_OUT_ARRAY(I,1), OGEL(I,1:6), 6_LONG )
            ENDDO
         ENDIF

      ELSE IF (TYPE == 'USERIN  ') THEN
         IF (WRITE_F06) THEN
            DO I=1,NUM
               CALL WRITE_STRESS_I8_PLUS_R14_LINE ( 1_LONG, EID_OUT_ARRAY(I,1), OGEL(I,1:6), 6_LONG )
            ENDDO
         ENDIF

      ELSE
         WRITE(ERR,9300) SUBR_NAME,TYPE
         WRITE(F06,9300) SUBR_NAME,TYPE
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error (elem type not valid) , so quit
! --- CQUADR_DKMQ24 end --- !
      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' OUTPUT FOR SUBCASE ',I8)

  102 FORMAT(' OUTPUT FOR EIGENVECTOR ',I8)

  103 FORMAT(' OUTPUT FOR CRAIG-BAMPTON DOF ',I8,' OF ',I8,' (boundary ',A,' for grid',I8,' component',I2,')')

  104 FORMAT(' OUTPUT FOR CRAIG-BAMPTON DOF ',I8,' OF ',I8,' (modal acceleration for mode ',I8,')')

  201 FORMAT(1X,A)

  301 FORMAT(A,'E L E M E N T   S T R E S S E S   I N   L O C A L   E L E M E N T   C O O R D I N A T E   S Y S T E M')

  302 FORMAT(A,'C B   E L E M E N T   S T R E S S E S   O T M   I N   L O C A L   E L E M E N T   C O O R D I N A T E',            &
  '   S Y S T E M')

  311 FORMAT(A,'E L E M E N T   S T R E S S E S   I N   M A T E R I A L   C O O R D I N A T E   S Y S T E M')

  312 FORMAT(A,'C B   E L E M E N T   S T R E S S E S   O T M   I N   M A T E R I A L   C O O R D I N A T E   S Y S T E M')

  321 FORMAT(A,'E L E M E N T   S T R E S S E S   I N   B A S I C   C O O R D I N A T E   S Y S T E M')

  331 FORMAT(A,'E L E M E N T   S T R E S S E S   I N   C O O R D I N A T E   S Y S T E M ',I8)

  401 FORMAT(A,'F O R   E L E M E N T   T Y P E   ',A11)



! BAR >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1101 FORMAT(                                                                                                                      &
          1X,A,'Element      SA1           SA2           SA3           SA4           Axial        SA-Max        SA-Min      M.S.-T'&
          ,'     Torsional'                                                                                                        &
       ,/,1X,A,'   ID        SB1           SB2           SB3           SB4          Stress        SB-Max        SB-Min      M.S.-C'&
          ,'   Stress/Margin')

 1102 FORMAT(  &
         1X,A,'Element      SA1           SA2           SA3           SA4          Axial         SA-Max        SA-Min      M.S.-T' &
      ,/,1X,A,'   ID        SB1           SB2           SB3           SB4          Stress        SB-Max        SB-Min      M.S.-C')

 1104 FORMAT(                                                                                                                      &
         1X,A,'                       S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )'                   &
      ,/,1X,A,'GRID   ELEMENT-ID       SXC           SXD           SXE           SXF         S-MAX          S-MIN         M.S.-T'      &
         ,'      M.S.-C'                                                                                                            &
      ,/,10X,'x/L')

 1105 FORMAT(                                                                                                                      &
         1X,A,'                       S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )'                   &
      ,/,1X,A,'GRID   ELEMENT-ID       SXC           SXD           SXE           SXF         S-MAX          S-MIN         M.S.-T'      &
         ,'      M.S.-C'                                                                                                            &
      ,/,10X,'x/L')

! ELAS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1201 FORMAT(1X,A,'Element     Stress     Element     Stress     Element     Stress     Element     Stress     Element     Stress' &
          ,/,1X,A,'   ID                     ID                     ID                     ID                     ID')

 1103 FORMAT(5(A,I8,1ES14.6))

! 3D Elems >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 1301 FORMAT(1X,A,'  Elem  Location            Sigma-xx      Sigma-yy      Sigma-zz       Tau-xy        Tau-yz        Tau-zx    ', &
             '   von Mises'                                                                                                        &
          ,/,1X,A,'   ID')

 1302 FORMAT(1X,A,'  Elem  Location            Sigma-xx      Sigma-yy      Sigma-zz       Tau-xy        Tau-yz        Tau-zx    ', &
             '      Octahedral Stress'                                                                                             &
          ,/,1X,A,'   ID',109X,'Direct        Shear')

 1303 FORMAT(1X,I8,2X,'CENTER  ',8X,8(1ES14.6))

 1304 FORMAT(28X,'------------- ------------- ------------- ------------- ------------- ------------- -------------',/,            &
             16X,'MAX* :     ',7(ES14.6),/,                                                                                        &
             16X,'MIN* :     ',7(ES14.6),//,                                                                                       &
             16X,'ABS* :     ',7(ES14.6),/                                                                                         &
             16X,'* for output set')

 1305 FORMAT(27X,' ------------- ------------- ------------- ------------- ------------- ------------- -------------',             &
                 ' -------------',/,                                                                                               &
             16X,'MAX* :     ',8(ES14.6),/,                                                                                        &
             16X,'MIN* :     ',8(ES14.6),//,                                                                                       &
             16X,'ABS* :     ',8(ES14.6),/                                                                                         &
             16X,'* for output set')

 1306 FORMAT(1X,A,10X,'GRD',I8,5X,8(1ES14.6))

! QUAD4 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1401 FORMAT(1X,A,'  Elem  Location       Stress      Stresses In Reported Coord System    Principal Stresses (Zero Shear)',       &
  '               Transverse   Transverse   % Poly',/,1X,A,                                                                        &
  '   ID                Distance     Normal-X     Normal-Y     Shear-XY     Angle     Major        Minor      von Mises',          &
  '    Shear-XZ     Shear-YZ    Fit Err',/,1X,A,119X,'(max through thickness)')

 1402 FORMAT(1X,A,'Elem  Location         Stress      Stresses In Reported Coord System    Principal Stresses (Zero Shear)',       &
  '       Max     Transverse   Transverse   % Poly',/,1X,A,                                                                        &
  ' ID                  Distance     Normal-X     Normal-Y     Shear-XY     Angle     Major        Minor      Shear-XY',           &
  '     Shear-XZ     Shear-YZ    Fit Err',/,1X,A,119X,'(max through thickness)')

 1403 FORMAT(1X,A,I8,2X,'CENTER  ',3X,1ES11.3,3(1ES13.5),0PF8.2,5(1ES13.5))

 1404 FORMAT(1X,A,21X,1ES11.3,3(1ES13.5),0PF8.2,3(1ES13.5))

 1405 FORMAT(1X,A,10X,'GRD',I8,1ES11.3,3(1ES13.5),0PF8.2,5(1ES13.5),E9.1,'(',I1,')')

 1406 FORMAT(1X,A,10X,'GRD',I8,1ES11.3,3(1ES13.5),0PF8.2,5(1ES13.5),E9.1)

 1407 FORMAT(1X,A,21X,1ES11.3,3(1ES13.5),0PF8.2,3(1ES13.5))

 1408 FORMAT(1X,A,32X,' ------------ ------------ ------------         ------------ ------------ ------------ ------------',       &
                 ' ------------ --------',/,                                                                                       &
             1X,A,'MAX* : ',25x,3(ES13.5),8X,5(ES13.5),E9.1,/,                                                                     &
             1X,A,'MIN* : ',25x,3(ES13.5),8X,5(ES13.5),E9.1,//,                                                                    &
             1X,A,'ABS* : ',25x,3(ES13.5),8X,5(ES13.5),E9.1,/,                                                                     &
             1X,A,'*for output set')

 1498 FORMAT(' NOTE: Explanation of errors in the polynomial fit to extrapolate element corner point stresses from values at the', &
                   ' Gauss points:')

 1499 FORMAT(6X,A)

! ROD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1501 FORMAT(  &
          1X,A,'Element     Axial       Safety     Torsional     Safety    Element    Axial       Safety     Torsional     Safety' &
       ,/,1X,A,'   ID       Stress      Margin       Stress      Margin       ID      Stress      Margin       Stress      Margin')

! SHEAR ----------------------------------------------------------------------------------------------------------------------------
 1601 FORMAT(1X,A,'Element              S t r e s s e s                           Element              S t r e s s e s'   &
          ,/,1X,A,'   ID      Normal-X      Normal-Y      Shear-XY                   ID      Normal-X      Normal-Y'               &
                 ,'      Shear-XY')


! TRIA3 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1701 FORMAT(1X,A,'Element    Location      Stress       Stresses In Reported Coord System      Principal Stresses (Zero Shear)',  &
                '               Transverse   Transverse'                                                                           &
          ,/,1X,A,'   ID                  Distance      Normal-X     Normal-Y      Shear-XY     Angle     Major        Minor'      &
          ,'      von Mises     Shear-XZ     Shear-YZ'                                                                             &
          ,/,1X,A,123X,'(max through thickness)')

 1702 FORMAT(1X,A,'Element    Location      Stress       Stresses In Reported Coord System      Principal Stresses (Zero Shear)',  &
  '      Max      Transverse   Transverse'                                                                                         &
          ,/,1X,A,'   ID                  Distance      Normal-X     Normal-Y      Shear-XY     Angle     Major        Minor',     &
          '      Shear-XY     Shear-XZ     Shear-YZ',/,1X,123X,'(max through thickness)')

 1703 FORMAT(1X,I8,4X,'CENTER  ',4X,4(1ES13.5),0PF9.3,5(1ES13.5))

 1704 FORMAT(23X,4(1ES13.5),0PF9.3,5(1ES13.5))

 1705 FORMAT(37X,'------------ ------------ ------------          ------------ ------------ ------------ ------------',            &
                 ' ------------',/,                                                                                                &
             1X,'MAX* : ',28x,3(ES13.5),9X,5(ES13.5),/,                                                                            &
             1X,'MIN* : ',28x,3(ES13.5),9X,5(ES13.5),//,                                                                           &
             1X,'ABS* : ',28x,3(ES13.5),9X,5(ES13.5),/,                                                                            &
             1X,'*for output set')

! BUSH >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1801 FORMAT(20X,A,'Element   Stress-1      Stress-2      Stress-3      Stress-4      Stress-5      Stress-6'                      &
          ,/,20X,A,'   ID')

 1802 FORMAT(19X,I8,8(1ES14.6))

! USERIN >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 1901 FORMAT(20X,A,'Element   Stress-1      Stress-2      Stress-3      Stress-4      Stress-5      Stress-6'                      &
          ,/,20X,A,'   ID')

 1902 FORMAT(19X,I8,8(1ES14.6))

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 9300 FORMAT(' *ERROR  9300: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NO OUTPUT FORMAT AVAILABLE FOR ELEMENT TYPE = ',A)
! **********************************************************************************************************************************
      END SUBROUTINE WRITE_ELEM_STRESSES
!==============================================================================

      SUBROUTINE WRITE_OES_CSHEAR(NUM, FILL, ISUBCASE, ITABLE, TITLE, SUBTITLE, LABEL, &
                                  FIELD5_INT_MODE, FIELD6_EIGENVALUE,                  &
                                  WRITE_F06, WRITE_OP2)
!     TODO: calculate margin
!
      USE PENTIUM_II_KIND, ONLY     :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY :  F06, OP2, ERR
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, OGEL
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_Value, IEEE_QUIET_NAN
      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: REAL32
      IMPLICIT NONE
      !
      INTEGER(LONG), INTENT(IN)       :: NUM               ! the number of elements
      INTEGER(LONG), INTENT(IN)       :: ISUBCASE          ! the current subcase
      CHARACTER(LEN=128), INTENT(IN)  :: TITLE             ! the model TITLE
      CHARACTER(LEN=128), INTENT(IN)  :: SUBTITLE          ! the subcase SUBTITLE
      CHARACTER(LEN=128), INTENT(IN)  :: LABEL             ! the subcase LABEL
      LOGICAL, INTENT(IN)             :: WRITE_F06, WRITE_OP2
      CHARACTER(128*BYTE)             :: FILL              ! Padding for output format

      INTEGER(LONG), INTENT(INOUT) :: ITABLE          ! the current subtable number
      !LOGICAL                         :: FIELD_5_INT_FLAG       ! flag to trigger FIELD5_INT_MODE vs. FIELD5_FLOAT_TIME_FREQ
      INTEGER(LONG)                   :: FIELD5_INT_MODE        ! int value for field 5
      !REAL(DOUBLE)                    :: FIELD5_FLOAT_TIME_FREQ ! float value for field 5
      REAL(DOUBLE)                    :: FIELD6_EIGENVALUE      ! float value for field 6

      INTEGER(LONG)               :: DEVICE_CODE  ! PLOT, PRINT, PUNCH flag
      INTEGER(LONG)               :: ANALYSIS_CODE = 1 ! static fallback for OP2 table-3 header
      INTEGER(LONG)               :: NUM_WIDE = 4     ! the number of "words" for an element
      INTEGER(LONG)               :: NVALUES          ! the number of "words" for all the elments
      INTEGER(LONG)               :: NTOTAL           ! the number of bytes for all NVALUES
      INTEGER(LONG)               :: ELEMENT_TYPE = 4 ! the OP2 flag for the element
      INTEGER(LONG)               :: STRESS_CODE      ! the OP2 flag for the stress
      REAL(DOUBLE)                :: ABS_ANS(3)       ! Max ABS for output
      REAL(DOUBLE)                :: MAX_ANS(3)       ! Max for output
      REAL(DOUBLE)                :: MIN_ANS(3)       ! Min for output
      INTEGER(LONG)               :: I, J             ! DO loop indices
      REAL(REAL32)  :: NAN
      NAN = IEEE_VALUE(NAN, IEEE_QUIET_NAN)

      IF (WRITE_OP2) THEN
          DEVICE_CODE = 1   ! PLOT
          NVALUES = NUM * NUM_WIDE
          NTOTAL = NVALUES * 4

          ! the below call is a placeholder to prevent a memory bug
          CALL GET_STRESS_CODE( STRESS_CODE, 1,            0,         1)
          ! eid, max_shear, avg_shear, margin
          CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                 TITLE, SUBTITLE, LABEL, FIELD5_INT_MODE, FIELD6_EIGENVALUE)

 100      FORMAT("*DEBUG: WRITE_CSHEAR    ITABLE=",I8, "; NUM=",I8,"; NVALUES=",I8,"; NTOTAL=",I8)
 101      FORMAT("*DEBUG: WRITE_CSHEAR    ITABLE=",I8," (should be -5, -7,...)")
          NVALUES = NUM * NUM_WIDE
          NTOTAL = NVALUES * 4
          WRITE(ERR,100) ITABLE,NUM,NVALUES,NTOTAL
          WRITE(OP2) NVALUES

          ! Nastran OP2 requires this write call be a one liner...so it's a little weird...
          ! translating:
          !    DO I=1,NUM
          !        WRITE(OP2) EID_OUT_ARRAY(I,1)*10+DEVICE_CODE  ! Nastran is weird and requires scaling the ELEMENT_ID
          !
          !        convert from float64 (double precision) to float32 (single precision)
          !        RE1 = REAL(OGEL(I,1), 4)
          !        RE2 = REAL(OGEL(I,2), 4)
          !        RE3 = REAL(OGEL(I,3), 4)
          !
          !        write the max_shear, avg_shear,
          !        WRITE(OP2) RE1, RE2, RE3
          !    ENDDO
          !
          ! write the CSHEAR stress/strain data
          !Normal-X      Normal-Y      Shear-XY -> max_shear, avg_shear, margin
          WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, REAL(OGEL(I,3), 4), REAL(OGEL(I,3), 4), &
                                                     NAN, I=1,NUM)
          WRITE(ERR,100) ITABLE
      ENDIF  ! write op2

      DO I=1,NUM,2
         IF (I+1 <= NUM) THEN
            CALL WRITE_STRESS_CSHEAR_PAIR_LINE ( EID_OUT_ARRAY(I,1), OGEL(I,1:3), EID_OUT_ARRAY(I+1,1), OGEL(I+1,1:3), .TRUE. )
         ELSE
            CALL WRITE_STRESS_CSHEAR_PAIR_LINE ( EID_OUT_ARRAY(I,1), OGEL(I,1:3), 0_LONG, OGEL(I,1:3), .FALSE. )
         ENDIF
      ENDDO

      MAX_ANS(1:3) = -HUGE(1.0D0)
      MIN_ANS(1:3) =  HUGE(1.0D0)
      ABS_ANS(1:3) =  0.0D0
      DO I=1,NUM
         DO J=1,3
            IF (OGEL(I,J) > MAX_ANS(J)) MAX_ANS(J) = OGEL(I,J)
            IF (OGEL(I,J) < MIN_ANS(J)) MIN_ANS(J) = OGEL(I,J)
            ABS_ANS(J) = MAX(ABS_ANS(J), DABS(OGEL(I,J)))
         ENDDO
      ENDDO

      WRITE(F06,1604) FILL(1: 0), FILL(1: 0), MAX_ANS(1),MAX_ANS(2),MAX_ANS(3),                                                 &
                      FILL(1: 0),             MIN_ANS(1),MIN_ANS(2),MIN_ANS(3),                                                 &
                      FILL(1: 0),             ABS_ANS(1),ABS_ANS(2),ABS_ANS(3)


 1603 FORMAT(1X,A,I8,3(1ES14.6),13X,I8,3(1ES14.6))
 1604 FORMAT(1X,A,'         ------------- ------------- ------------- ',20X,' ------------- ------------- ------------- ',/,       &
             1X,A,'MAX* : ',1X,3(ES14.6),/,                                                                                        &
             1X,A,'MIN* : ',1X,3(ES14.6),//,                                                                                       &
             1X,A,'ABS* : ',1X,3(ES14.6),/,                                                                                        &
             1X,A,'*for output set')
      END SUBROUTINE WRITE_OES_CSHEAR

!==============================================================================
      SUBROUTINE WRITE_OES_CTRIA3 ( NUM, NUM_PTS, FILL, ISUBCASE, ITABLE, TITLE, SUBTITLE, LABEL, &
                                    FIELD5_INT_MODE, FIELD6_EIGENVALUE ,                 &
                                    WRITE_F06, WRITE_OP2, OP2_OPENED, WRITE_GPSTRESS_F06)
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, OP2
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, GID_OUT_ARRAY, OGEL, SHELL_OUT_TE
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRE_LOC, GPSTRESS_REQ
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE GET_MAX_MIN_ABS_STR_Interface
      USE FAST_OUTPUT_FORMATTERS, ONLY:  FAST_BUILD_TRIA_1703_LINE, FAST_BUILD_TRIA_1704_LINE, FAST_BUILD_TRIA_1706_LINE
      IMPLICIT NONE
      !
      INTEGER(LONG), INTENT(IN)       :: NUM               ! the number of elements
      INTEGER(LONG), INTENT(IN)       :: NUM_PTS           ! number of recovered stress points per element
      INTEGER(LONG), INTENT(IN)       :: ISUBCASE          ! the current subcase
      CHARACTER(LEN=128), INTENT(IN)  :: TITLE             ! the model TITLE
      CHARACTER(LEN=128), INTENT(IN)  :: SUBTITLE          ! the subcase SUBTITLE
      CHARACTER(LEN=128), INTENT(IN)  :: LABEL             ! the subcase LABEL
      LOGICAL, INTENT(IN)             :: WRITE_F06, WRITE_OP2, OP2_OPENED, WRITE_GPSTRESS_F06  ! flags
      CHARACTER(128*BYTE)             :: FILL              ! Padding for output format

      INTEGER(LONG), INTENT(INOUT) :: ITABLE           ! the current subtable number
      !LOGICAL                         :: FIELD_5_INT_FLAG       ! flag to trigger FIELD5_INT_MODE vs. FIELD5_FLOAT_TIME_FREQ
      INTEGER(LONG)                   :: FIELD5_INT_MODE        ! int value for field 5
      !REAL(DOUBLE)                    :: FIELD5_FLOAT_TIME_FREQ ! float value for field 5
      REAL(DOUBLE)                    :: FIELD6_EIGENVALUE      ! float value for field 6

      INTEGER(LONG)               :: DEVICE_CODE = 1   ! PLOT, PRINT, PUNCH flag; set as PLOT
      INTEGER(LONG)               :: ANALYSIS_CODE = 1 ! static fallback for OP2 table-3 header
      INTEGER(LONG)               :: NUM_WIDE          ! the number of "words" for an element
      INTEGER(LONG)               :: NVALUES           ! the number of "words" for all the elments
      INTEGER(LONG)               :: NTOTAL            ! the number of bytes for all NVALUES
      INTEGER(LONG)               :: ELEMENT_TYPE      ! the OP2 flag for the element
      INTEGER(LONG)               :: STRESS_CODE = 1   ! the OP2 flag for the stress; preallocate
      REAL(DOUBLE)                :: ABS_ANS(11)       ! Max ABS for output
      REAL(DOUBLE)                :: ANGLE
      REAL(DOUBLE)                :: MAX_ANS(11)       ! Max for output
      REAL(DOUBLE)                :: MEAN
      REAL(DOUBLE)                :: MIN_ANS(11)       ! Min for output
      REAL(DOUBLE)                :: ROW_CURV(10)
      REAL(DOUBLE)                :: ROW_MEM(10)
      REAL(DOUBLE)                :: ROW_TMP(10)
      REAL(DOUBLE)                :: SMAJ
      REAL(DOUBLE)                :: SMIN
      REAL(DOUBLE)                :: SXYMAX
      REAL(DOUBLE)                :: VONMISES
      REAL(DOUBLE)                :: Z1, Z2, Z_DEN
      INTEGER(LONG)               :: I, J, K, L, IP       ! DO loop indices
      CHARACTER(159*BYTE)         :: TRIA_CENTER_LINE
      CHARACTER(159*BYTE)         :: TRIA_LOWER_LINE
      CHARACTER(159*BYTE)         :: TRIA_GRID_LINE

      ! [eid, fiber_dist/curvature, oxx, oyy, txy, angle, omax, omin, ovm/max_shear,   ! upper
      !       fiber_dist/curvature, oxx, oyy, txy, angle, omax, omin, ovm/max_shear,   ! lower
      !]
      K = 0

!100  FORMAT("*DEBUG: WRITE_CTRIA3    ITABLE=",I8, "; NUM=",I8,"; NVALUES=",I8,"; NTOTAL=",I8)
!101  FORMAT("*DEBUG: WRITE_CTRIA3    ITABLE=",I8," (should be -5, -7,...)")
      NUM_WIDE = 17
      ELEMENT_TYPE = 74
      NVALUES = NUM * NUM_WIDE
      NTOTAL = NVALUES * 4
!      WRITE(ERR,100) ITABLE,NUM,NVALUES,NTOTAL

      IF (WRITE_OP2) THEN
          !CALL GET_STRESS_CODE(STRESS_CODE, IS_VON_MISES, IS_STRAIN, IS_FIBER_DISTANCE)
          CALL GET_STRESS_CODE( STRESS_CODE, 1,            0,         1)
          CALL WRITE_OES3_STATIC(ITABLE, ISUBCASE, DEVICE_CODE, ELEMENT_TYPE, NUM_WIDE, STRESS_CODE, &
                                 TITLE, SUBTITLE, LABEL, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
          WRITE(OP2) NVALUES
          WRITE(OP2) (EID_OUT_ARRAY(I,1)*10+DEVICE_CODE, (REAL(OGEL(2*I-1,J),4), J=1,8),                                 &
                     (REAL(OGEL(2*I,J),4), J=1,8), I=1,NUM)
      ENDIF  ! write op2

 1703 FORMAT(1X,I8,4X,'CENTER  ',4X,4(1ES13.5),0PF9.3,5(1ES13.5))

 1704 FORMAT(23X,4(1ES13.5),0PF9.3,5(1ES13.5))
 1706 FORMAT(1X,A,I8,4X,4(1ES13.5),0PF9.3,5(1ES13.5))

 1705 FORMAT(37X,'------------ ------------ ------------          ------------ ------------ ------------ ------------',            &
                 ' ------------',/,                                                                                                &
             1X,'MAX* : ',28x,3(ES13.5),9X,5(ES13.5),/,                                                                            &
             1X,'MIN* : ',28x,3(ES13.5),9X,5(ES13.5),//,                                                                           &
             1X,'ABS* : ',28x,3(ES13.5),9X,5(ES13.5),/,                                                                            &
             1X,'*for output set')

      IF (STRE_LOC == 'CENTER  ') THEN
         DO I=1,NUM,MAX(1_LONG,NUM_PTS)
            K = 2*I - 1
            IF (WRITE_F06) WRITE(F06,*)
            DO J=1,10
               ROW_TMP(J) = OGEL(K,J)
            ENDDO
            IF (WRITE_F06) THEN
               CALL FAST_BUILD_TRIA_1703_LINE ( EID_OUT_ARRAY(I,1), ROW_TMP, TRIA_CENTER_LINE )
               WRITE(F06,'(A)') TRIA_CENTER_LINE(1:159)
            ENDIF
            DO J=1,10
               ROW_TMP(J) = OGEL(K+1,J)
            ENDDO
            IF (WRITE_F06) THEN
               CALL FAST_BUILD_TRIA_1704_LINE ( ROW_TMP, TRIA_LOWER_LINE )
               WRITE(F06,'(A)') TRIA_LOWER_LINE(1:159)
            ENDIF
         ENDDO
      ELSE
         DO I=1,NUM,NUM_PTS
            DO L=1,NUM_PTS-1
               IP = I + L
               K = 2*IP - 1
               IF (WRITE_F06) WRITE(F06,*)
               DO J=1,10
                  ROW_TMP(J) = OGEL(K,J)
               ENDDO
               IF (WRITE_F06) THEN
                  CALL FAST_BUILD_TRIA_1706_LINE ( EID_OUT_ARRAY(I,1), GID_OUT_ARRAY(I,L+1), ROW_TMP, TRIA_GRID_LINE )
                  WRITE(F06,'(A)') TRIA_GRID_LINE(1:159)
               ENDIF
               DO J=1,10
                  ROW_TMP(J) = OGEL(K+1,J)
               ENDDO
               IF (WRITE_F06) THEN
                  CALL FAST_BUILD_TRIA_1704_LINE ( ROW_TMP, TRIA_LOWER_LINE )
                  WRITE(F06,'(A)') TRIA_LOWER_LINE(1:159)
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      CALL GET_MAX_MIN_ABS_STR ( NUM, 10, 'Y', MAX_ANS, MIN_ANS, ABS_ANS )

      IF (WRITE_F06) THEN
         WRITE(F06,1705) MAX_ANS(2),MAX_ANS(3),MAX_ANS(4),MAX_ANS(6),MAX_ANS(7),MAX_ANS(8),MAX_ANS(9),MAX_ANS(10),              &
                         MIN_ANS(2),MIN_ANS(3),MIN_ANS(4),MIN_ANS(6),MIN_ANS(7),MIN_ANS(8),MIN_ANS(9),MIN_ANS(10),              &
                         ABS_ANS(2),ABS_ANS(3),ABS_ANS(4),ABS_ANS(6),ABS_ANS(7),ABS_ANS(8),ABS_ANS(9),ABS_ANS(10)
      ENDIF

      IF ((OP2_OPENED .OR. WRITE_GPSTRESS_F06) .AND. GPSTRESS_REQ .AND. (NUM_PTS > 1)) THEN
         CALL WRITE_OGS1_SURFACE_STRESS ( ITABLE, ISUBCASE, NUM, NUM_PTS, DEVICE_CODE, ANALYSIS_CODE, FIELD5_INT_MODE,          &
                                          FIELD6_EIGENVALUE, TITLE, SUBTITLE, LABEL, 'TRIA', WRITE_GPSTRESS_F06, OP2_OPENED )
      ENDIF

      END SUBROUTINE WRITE_OES_CTRIA3

!==============================================================================
      SUBROUTINE WRITE_OGS1_SURFACE_STRESS ( ITABLE, ISUBCASE, NUM, NUM_PTS, DEVICE_CODE, ANALYSIS_CODE, FIELD5_INT_MODE,        &
                                             FIELD6_EIGENVALUE, TITLE, SUBTITLE, LABEL, FAMILY, WRITE_F06, WRITE_OP2 )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  OP2, F06
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, GID_OUT_ARRAY, OGEL, SHELL_OUT_TE
      USE MODEL_STUF, ONLY            :  GRID_ID, RGRID
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  GPSTRESS_REQ, NUM_GP_SURFACE, MAX_GP_SURFACES, GP_SURFACE_IDS,                 &
                                         GP_SURFACE_NORMAL_MODE
      USE GPSTRESS_SURFACE_UTILS, ONLY:  GPSTRESS_COLLECT_SURFACE_PATCH
      USE FAST_OUTPUT_FORMATTERS, ONLY:  FAST_FMT_I8_RJ, FAST_FMT_F06_E14_6

      USE GET_ARRAY_ROW_NUM_Interface

      IMPLICIT NONE

      INTEGER(LONG), INTENT(INOUT)    :: ITABLE
      INTEGER(LONG), INTENT(IN)       :: ISUBCASE
      INTEGER(LONG), INTENT(IN)       :: NUM
      INTEGER(LONG), INTENT(IN)       :: NUM_PTS
      INTEGER(LONG), INTENT(IN)       :: DEVICE_CODE
      INTEGER(LONG), INTENT(IN)       :: ANALYSIS_CODE
      INTEGER(LONG), INTENT(IN)       :: FIELD5_INT_MODE
      REAL(DOUBLE), INTENT(IN)        :: FIELD6_EIGENVALUE
      CHARACTER(LEN=128), INTENT(IN)  :: TITLE
      CHARACTER(LEN=128), INTENT(IN)  :: SUBTITLE
      CHARACTER(LEN=128), INTENT(IN)  :: LABEL
      CHARACTER(LEN=*), INTENT(IN)    :: FAMILY
      LOGICAL, INTENT(IN)             :: WRITE_F06
      LOGICAL, INTENT(IN)             :: WRITE_OP2

      CHARACTER(8*BYTE)               :: TABLE_NAME
      CHARACTER(LEN=128)              :: TITLE2
      CHARACTER(LEN=128)              :: SUBTITLE2
      CHARACTER(LEN=128)              :: LABEL2
      INTEGER(LONG)                   :: APPROACH_CODE
      INTEGER(LONG)                   :: AXIS
      INTEGER(LONG)                   :: FORMAT_CODE
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: IERR
      INTEGER(LONG)                   :: L
      INTEGER(LONG)                   :: NOUT
      INTEGER(LONG)                   :: NROWS
      INTEGER(LONG)                   :: NUM_WIDE
      INTEGER(LONG)                   :: NVALUES
      INTEGER(LONG)                   :: OCOORD
      INTEGER(LONG)                   :: OGS_ID
      INTEGER(LONG)                   :: OGS_ITABLE
      INTEGER(LONG)                   :: REFID
      INTEGER(LONG)                   :: S_CODE
      INTEGER(LONG)                   :: SURF
      INTEGER(LONG)                   :: TABLE_CODE
      INTEGER(LONG)                   :: THERMAL
      INTEGER(LONG), ALLOCATABLE      :: PATCH_ELEMS(:)
      INTEGER(LONG), ALLOCATABLE      :: PATCH_GRIDS(:)
      INTEGER(LONG), ALLOCATABLE      :: OUT_EIDS(:)
      INTEGER(LONG), ALLOCATABLE      :: OUT_GRIDS(:)
      INTEGER(LONG)                   :: SURF_START(MAX_GP_SURFACES)
      INTEGER(LONG)                   :: SURF_END(MAX_GP_SURFACES)
      REAL(DOUBLE), ALLOCATABLE       :: OUT_Z1(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: OUT_Z2(:,:)
      REAL(DOUBLE)                    :: FIELD7
      REAL(DOUBLE)                    :: MID_ROW(8)
      REAL(DOUBLE)                    :: MID_STRESS3(3)
      LOGICAL                         :: USED_RECOVERY

      IF (NUM <= 0) RETURN

      ! Keep OGS1 as a separate OP2 result table; do not interleave it inside an open OES table.
      OGS_ITABLE = -1
      IF (WRITE_OP2) THEN
         IF (ITABLE < -1) THEN
            CALL END_OP2_TABLE(ITABLE)
            ITABLE = -1
         ENDIF

         TABLE_NAME = 'OGS1    '
         CALL WRITE_TABLE_HEADER(TABLE_NAME)
         OGS_ITABLE = -3
         CALL WRITE_ITABLE(OGS_ITABLE)
      ENDIF

      IF ((ANALYSIS_CODE == 1) .OR. (ANALYSIS_CODE == 10)) THEN
         FIELD7 = ZERO
      ELSE
         FIELD7 = SQRT(ABS(FIELD6_EIGENVALUE))
      ENDIF

      APPROACH_CODE = ANALYSIS_CODE * 10 + DEVICE_CODE
      TABLE_CODE = 26
      OGS_ID = 100
      REFID = 0
      FORMAT_CODE = 1
      NUM_WIDE = 11
      S_CODE = 0
      OCOORD = 2
      AXIS = 0
      THERMAL = 0
      TITLE2 = TITLE(1:100)
      SUBTITLE2 = SUBTITLE(1:67)
      LABEL2 = LABEL(1:100)

      IF (WRITE_OP2) THEN
         WRITE(OP2) 146
         WRITE(OP2) APPROACH_CODE, TABLE_CODE, OGS_ID, ISUBCASE, FIELD5_INT_MODE,                                        &
               REAL(FIELD6_EIGENVALUE, 4), REAL(FIELD7, 4), REFID, FORMAT_CODE, NUM_WIDE, S_CODE, OCOORD, AXIS, 0, 0,    &
               0, 0, 0, 0, 0,                                                                                            &
               0, 0, THERMAL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                                                              &
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                                                                    &
               0, 0, 0, 0,                                                                                               &
               TITLE2, SUBTITLE2, LABEL2

         OGS_ITABLE = OGS_ITABLE - 1
         CALL WRITE_ITABLE(OGS_ITABLE)
         OGS_ITABLE = OGS_ITABLE - 1
      ENDIF

      USED_RECOVERY = .FALSE.
      NOUT = 0
      SURF_START = 0
      SURF_END = 0
      IF (GPSTRESS_REQ .AND. (NUM_GP_SURFACE > 0) .AND. ((FAMILY(1:4) == 'TRIA') .OR. (FAMILY(1:4) == 'QUAD'))) THEN
         ALLOCATE(PATCH_ELEMS(MAX(1_LONG,NUM)))
         ALLOCATE(PATCH_GRIDS(MAX(1_LONG,4*NUM)))
         ALLOCATE(OUT_EIDS(MAX(1_LONG,4*NUM)))
         ALLOCATE(OUT_GRIDS(MAX(1_LONG,4*NUM)))
         ALLOCATE(OUT_Z1(8,MAX(1_LONG,4*NUM)))
         ALLOCATE(OUT_Z2(8,MAX(1_LONG,4*NUM)))
         OUT_EIDS = 0
         OUT_GRIDS = 0
         OUT_Z1 = ZERO
         OUT_Z2 = ZERO
         DO SURF=1,NUM_GP_SURFACE
            SURF_START(SURF) = NOUT + 1
            CALL BUILD_SURFACE_PATCH_OUTPUT ( SURF, PATCH_ELEMS, PATCH_GRIDS, OUT_EIDS, OUT_GRIDS, OUT_Z1, OUT_Z2, NOUT, IERR )
            IF (IERR == 0) USED_RECOVERY = .TRUE.
            SURF_END(SURF) = NOUT
         ENDDO
      ENDIF

      IF (USED_RECOVERY .AND. (NOUT > 0)) THEN
         NROWS = 2 * NOUT
      ELSE
         IF (FAMILY(1:4) == 'TRIA') THEN
            NROWS = 6 * NUM
         ELSE
            NROWS = 2 * (NUM_PTS - 1) * (NUM / NUM_PTS)
         ENDIF
      ENDIF
      NVALUES = NUM_WIDE * NROWS
      IF (WRITE_OP2) WRITE(OP2) NVALUES

      IF (WRITE_OP2 .AND. USED_RECOVERY .AND. (NOUT > 0)) THEN
         WRITE(OP2) (OUT_GRIDS(I)*10+DEVICE_CODE, OUT_EIDS(I), 'Z1  ',                                                   &
                     REAL(OUT_Z1(1,I),4), REAL(OUT_Z1(2,I),4), REAL(OUT_Z1(3,I),4), REAL(OUT_Z1(4,I),4),                &
                     REAL(OUT_Z1(5,I),4), REAL(OUT_Z1(6,I),4), REAL(OUT_Z1(7,I),4), REAL(OUT_Z1(8,I),4),                &
                     OUT_GRIDS(I)*10+DEVICE_CODE, OUT_EIDS(I), 'Z2  ',                                                   &
                     REAL(OUT_Z2(1,I),4), REAL(OUT_Z2(2,I),4), REAL(OUT_Z2(3,I),4), REAL(OUT_Z2(4,I),4),                &
                     REAL(OUT_Z2(5,I),4), REAL(OUT_Z2(6,I),4), REAL(OUT_Z2(7,I),4), REAL(OUT_Z2(8,I),4), I=1,NOUT)
      ELSE IF (WRITE_OP2 .AND. (FAMILY(1:4) == 'TRIA')) THEN
         WRITE(OP2) ((GID_OUT_ARRAY(I,L+1)*10+DEVICE_CODE, EID_OUT_ARRAY(I,1), 'Z1  ',                                  &
                      REAL(OGEL(2*I-1,2),4), REAL(OGEL(2*I-1,3),4), REAL(OGEL(2*I-1,4),4),                              &
                      REAL(OGEL(2*I-1,6),4), REAL(OGEL(2*I-1,7),4), REAL(OGEL(2*I-1,8),4),                              &
                      REAL(0.5D0*ABS(OGEL(2*I-1,7)-OGEL(2*I-1,8)),4), REAL(OGEL(2*I-1,9),4),                            &
                      GID_OUT_ARRAY(I,L+1)*10+DEVICE_CODE, EID_OUT_ARRAY(I,1), 'Z2  ',                                  &
                      REAL(OGEL(2*I,2),4), REAL(OGEL(2*I,3),4), REAL(OGEL(2*I,4),4), REAL(OGEL(2*I,6),4),               &
                      REAL(OGEL(2*I,7),4), REAL(OGEL(2*I,8),4), REAL(0.5D0*ABS(OGEL(2*I,7)-OGEL(2*I,8)),4),             &
                      REAL(OGEL(2*I,9),4), L=1,3), I=1,NUM)
      ELSE IF (WRITE_OP2) THEN
         WRITE(OP2) ((GID_OUT_ARRAY(I,L+1)*10+DEVICE_CODE, EID_OUT_ARRAY(I,1), 'Z1  ',                                 &
                      REAL(OGEL(I+2*L,2),4), REAL(OGEL(I+2*L,3),4), REAL(OGEL(I+2*L,4),4),                              &
                      REAL(OGEL(I+2*L,6),4), REAL(OGEL(I+2*L,7),4), REAL(OGEL(I+2*L,8),4),                              &
                      REAL(0.5D0*ABS(OGEL(I+2*L,7)-OGEL(I+2*L,8)),4), REAL(OGEL(I+2*L,9),4),                            &
                      GID_OUT_ARRAY(I,L+1)*10+DEVICE_CODE, EID_OUT_ARRAY(I,1), 'Z2  ',                                  &
                      REAL(OGEL(I+2*L+1,2),4), REAL(OGEL(I+2*L+1,3),4), REAL(OGEL(I+2*L+1,4),4),                        &
                      REAL(OGEL(I+2*L+1,6),4), REAL(OGEL(I+2*L+1,7),4), REAL(OGEL(I+2*L+1,8),4),                        &
                      REAL(0.5D0*ABS(OGEL(I+2*L+1,7)-OGEL(I+2*L+1,8)),4), REAL(OGEL(I+2*L+1,9),4),                      &
                      L=1,NUM_PTS-1), I=1,NUM,NUM_PTS)
      ENDIF

      IF (WRITE_F06 .AND. USED_RECOVERY .AND. (NOUT > 0)) THEN
         DO SURF=1,NUM_GP_SURFACE
            IF ((SURF_START(SURF) <= 0) .OR. (SURF_END(SURF) < SURF_START(SURF))) CYCLE
            WRITE(F06,*)
            WRITE(F06,'(1X,A)') TITLE
            WRITE(F06,*)
            WRITE(F06,'(''0     '',A,101X,''SUBCASE '',I8)') TRIM(SUBTITLE), ISUBCASE
            WRITE(F06,'(34X,''S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E'',I8)') GP_SURFACE_IDS(SURF)
            WRITE(F06,'(''0'',23X,''SURFACE X-AXIS X  NORMAL(Z-AXIS)  Z         REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID'',I9)') 0
            WRITE(F06,'(5X,''GRID'',6X,''ELEMENT'',12X,''STRESSES IN SURFACE SYSTEM'',11X,''PRINCIPAL STRESSES'',12X,''MAX'',/,&
     &                  5X,''ID'',10X,''ID'',4X,''FIBER'',3X,''NORMAL-X'',3X,''NORMAL-Y'',3X,''SHEAR-XY'',5X,''ANGLE'',6X,''MAJOR'',6X,''MINOR'',6X,''SHEAR'',5X,''VON MISES'')')
            DO I=SURF_START(SURF),SURF_END(SURF)
               MID_STRESS3(1:3) = 0.5D0 * (OUT_Z1(1:3,I) + OUT_Z2(1:3,I))
               CALL BUILD_SURFACE_RESULT_ROW ( MID_STRESS3, MID_ROW )
               CALL WRITE_OGS1_F06_ROW ( OUT_GRIDS(I), OUT_EIDS(I), 'Z1 ', OUT_Z1(:,I) )
               CALL WRITE_OGS1_F06_ROW ( 0_LONG,       OUT_EIDS(I), 'Z2 ', OUT_Z2(:,I) )
               CALL WRITE_OGS1_F06_ROW ( 0_LONG,       OUT_EIDS(I), 'MID', MID_ROW )
            ENDDO
         ENDDO
      ENDIF

      IF (WRITE_OP2) THEN
         CALL END_OP2_TABLE(OGS_ITABLE)
         ITABLE = 0
      ENDIF

      IF (ALLOCATED(PATCH_ELEMS)) DEALLOCATE(PATCH_ELEMS)
      IF (ALLOCATED(PATCH_GRIDS)) DEALLOCATE(PATCH_GRIDS)
      IF (ALLOCATED(OUT_EIDS)) DEALLOCATE(OUT_EIDS)
      IF (ALLOCATED(OUT_GRIDS)) DEALLOCATE(OUT_GRIDS)
      IF (ALLOCATED(OUT_Z1)) DEALLOCATE(OUT_Z1)
      IF (ALLOCATED(OUT_Z2)) DEALLOCATE(OUT_Z2)

      CONTAINS

      SUBROUTINE BUILD_SURFACE_PATCH_OUTPUT ( SURF_INDEX, PATCH_ELEMS, PATCH_GRIDS, OUT_EIDS, OUT_GRIDS, OUT_Z1, OUT_Z2, NOUT, IERR )

      INTEGER(LONG), INTENT(IN)       :: SURF_INDEX
      INTEGER(LONG), INTENT(INOUT)    :: PATCH_ELEMS(:)
      INTEGER(LONG), INTENT(INOUT)    :: PATCH_GRIDS(:)
      INTEGER(LONG), INTENT(INOUT)    :: OUT_EIDS(:)
      INTEGER(LONG), INTENT(INOUT)    :: OUT_GRIDS(:)
      REAL(DOUBLE), INTENT(INOUT)     :: OUT_Z1(:,:)
      REAL(DOUBLE), INTENT(INOUT)     :: OUT_Z2(:,:)
      INTEGER(LONG), INTENT(INOUT)    :: NOUT
      INTEGER(LONG), INTENT(OUT)      :: IERR

      INTEGER(LONG)                   :: EID
      INTEGER(LONG)                   :: GID
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: IELEM
      INTEGER(LONG)                   :: IPATCH
      INTEGER(LONG)                   :: J
      INTEGER(LONG)                   :: NPATCH_ELEMS
      INTEGER(LONG)                   :: NPATCH_GRIDS
      INTEGER(LONG)                   :: ROW1
      INTEGER(LONG)                   :: ROW2
      REAL(DOUBLE)                    :: AREA
      REAL(DOUBLE)                    :: GXYZ(3)
      REAL(DOUBLE)                    :: SUMW1, SUMW2
      REAL(DOUBLE)                    :: VALS1(3), VALS2(3), OUTVAL(8)

      IERR = 0
      CALL GPSTRESS_COLLECT_SURFACE_PATCH ( SURF_INDEX, PATCH_ELEMS, NPATCH_ELEMS, PATCH_GRIDS, NPATCH_GRIDS, IERR )
      IF (IERR /= 0) RETURN

      DO IPATCH=1,NPATCH_GRIDS
         GID = PATCH_GRIDS(IPATCH)
         IF (FIND_INT(GID, OUT_GRIDS, NOUT) > 0) CYCLE

         SUMW1 = ZERO
         SUMW2 = ZERO
         VALS1 = ZERO
         VALS2 = ZERO
         EID = 0

         CALL GET_GRID_BASIC_COORDS ( GID, GXYZ )

         IF (FAMILY(1:4) == 'TRIA') THEN
            DO I=1,NUM
               IELEM = EID_OUT_ARRAY(I,1)
               IF (FIND_INT(IELEM, PATCH_ELEMS, NPATCH_ELEMS) == 0) CYCLE
               AREA = GET_TRIA_AREA(I)
               IF (AREA <= ZERO) CYCLE
               DO J=1,3
                  IF (GID_OUT_ARRAY(I,J+1) /= GID) CYCLE
                  ROW1 = 2*I - 1
                  ROW2 = 2*I
                  CALL GET_SURFACE_STRESS3 ( SURF_INDEX, I, (/ OGEL(ROW1,2), OGEL(ROW1,3), OGEL(ROW1,4) /), OUTVAL(1:3) )
                  VALS1 = VALS1 + (AREA/3.0D0) * OUTVAL(1:3)
                  CALL GET_SURFACE_STRESS3 ( SURF_INDEX, I, (/ OGEL(ROW2,2), OGEL(ROW2,3), OGEL(ROW2,4) /), OUTVAL(1:3) )
                  VALS2 = VALS2 + (AREA/3.0D0) * OUTVAL(1:3)
                  SUMW1 = SUMW1 + AREA/3.0D0
                  SUMW2 = SUMW2 + AREA/3.0D0
                  IF (EID == 0) EID = IELEM
               ENDDO
            ENDDO
         ELSE
            DO I=1,NUM,NUM_PTS
               IELEM = EID_OUT_ARRAY(I,1)
               IF (FIND_INT(IELEM, PATCH_ELEMS, NPATCH_ELEMS) == 0) CYCLE
               AREA = GET_QUAD_AREA(I)
               IF (AREA <= ZERO) CYCLE
               DO J=1,NUM_PTS-1
                  IF (GID_OUT_ARRAY(I,J+1) /= GID) CYCLE
                  ROW1 = 2*I + 2*J - 1
                  ROW2 = 2*I + 2*J
                  CALL GET_SURFACE_STRESS3 ( SURF_INDEX, I+J, (/ OGEL(ROW1,2), OGEL(ROW1,3), OGEL(ROW1,4) /), OUTVAL(1:3) )
                  VALS1 = VALS1 + (AREA/4.0D0) * OUTVAL(1:3)
                  CALL GET_SURFACE_STRESS3 ( SURF_INDEX, I+J, (/ OGEL(ROW2,2), OGEL(ROW2,3), OGEL(ROW2,4) /), OUTVAL(1:3) )
                  VALS2 = VALS2 + (AREA/4.0D0) * OUTVAL(1:3)
                  SUMW1 = SUMW1 + AREA/4.0D0
                  SUMW2 = SUMW2 + AREA/4.0D0
                  IF (EID == 0) EID = IELEM
               ENDDO
            ENDDO
         ENDIF

         IF ((SUMW1 > ZERO) .AND. (SUMW2 > ZERO)) THEN
            IF (NOUT >= SIZE(OUT_GRIDS)) EXIT
            NOUT = NOUT + 1
            OUT_GRIDS(NOUT) = GID
            OUT_EIDS(NOUT) = 0_LONG
            CALL BUILD_SURFACE_RESULT_ROW ( VALS1 / SUMW1, OUTVAL )
            OUT_Z1(:,NOUT) = OUTVAL
            CALL BUILD_SURFACE_RESULT_ROW ( VALS2 / SUMW2, OUTVAL )
            OUT_Z2(:,NOUT) = OUTVAL
         ENDIF
      ENDDO

      END SUBROUTINE BUILD_SURFACE_PATCH_OUTPUT

!----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG) FUNCTION GET_FIRST_PATCH_ELEM_FOR_GRID ( GID, PATCH_ELEMS, NPATCH_ELEMS )

      INTEGER(LONG), INTENT(IN)       :: GID
      INTEGER(LONG), INTENT(IN)       :: PATCH_ELEMS(:)
      INTEGER(LONG), INTENT(IN)       :: NPATCH_ELEMS

      INTEGER(LONG)                   :: I, J, IELEM

      GET_FIRST_PATCH_ELEM_FOR_GRID = 0
      IF (FAMILY(1:4) == 'TRIA') THEN
         DO I=1,NUM
            IELEM = EID_OUT_ARRAY(I,1)
            IF (FIND_INT(IELEM, PATCH_ELEMS, NPATCH_ELEMS) == 0) CYCLE
            DO J=1,3
               IF (GID_OUT_ARRAY(I,J+1) == GID) THEN
                  GET_FIRST_PATCH_ELEM_FOR_GRID = IELEM
                  RETURN
               ENDIF
            ENDDO
         ENDDO
      ELSE
         DO I=1,NUM,NUM_PTS
            IELEM = EID_OUT_ARRAY(I,1)
            IF (FIND_INT(IELEM, PATCH_ELEMS, NPATCH_ELEMS) == 0) CYCLE
            DO J=1,NUM_PTS-1
               IF (GID_OUT_ARRAY(I,J+1) == GID) THEN
                  GET_FIRST_PATCH_ELEM_FOR_GRID = IELEM
                  RETURN
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      END FUNCTION GET_FIRST_PATCH_ELEM_FOR_GRID

!----------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE BUILD_SURFACE_RESULT_ROW ( STRESS3, RESULT8 )

      REAL(DOUBLE), INTENT(IN)        :: STRESS3(3)
      REAL(DOUBLE), INTENT(OUT)       :: RESULT8(8)

      REAL(DOUBLE)                    :: ANGLE, AVG, DIFF, RAD, SMAJ, SMIN, SXYMAX, VM
      REAL(DOUBLE), PARAMETER         :: RAD2DEG = 180.0D0 / 3.1415926535897932384626433832795D0

      AVG = 0.5D0 * (STRESS3(1) + STRESS3(2))
      DIFF = 0.5D0 * (STRESS3(1) - STRESS3(2))
      RAD = DSQRT(DIFF*DIFF + STRESS3(3)*STRESS3(3))
      SMAJ = AVG + RAD
      SMIN = AVG - RAD
      ANGLE = 0.5D0 * DATAN2(2.0D0*STRESS3(3), STRESS3(1) - STRESS3(2)) * RAD2DEG
      SXYMAX = 0.5D0 * DABS(SMAJ - SMIN)
      VM = DSQRT(STRESS3(1)*STRESS3(1) - STRESS3(1)*STRESS3(2) + STRESS3(2)*STRESS3(2) + 3.0D0*STRESS3(3)*STRESS3(3))

      RESULT8(1) = STRESS3(1)
      RESULT8(2) = STRESS3(2)
      RESULT8(3) = STRESS3(3)
      RESULT8(4) = ANGLE
      RESULT8(5) = SMAJ
      RESULT8(6) = SMIN
      RESULT8(7) = SXYMAX
      RESULT8(8) = VM

      END SUBROUTINE BUILD_SURFACE_RESULT_ROW

!----------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE GET_SURFACE_STRESS3 ( SURF_INDEX, POINT_INDEX, LOCAL_STRESS3, SURF_STRESS3 )

      INTEGER(LONG), INTENT(IN)       :: SURF_INDEX
      INTEGER(LONG), INTENT(IN)       :: POINT_INDEX
      REAL(DOUBLE), INTENT(IN)        :: LOCAL_STRESS3(3)
      REAL(DOUBLE), INTENT(OUT)       :: SURF_STRESS3(3)

      REAL(DOUBLE)                    :: SURF_BASIS(3,3)
      REAL(DOUBLE)                    :: LOCAL_TENSOR(3,3)
      REAL(DOUBLE)                    :: SURF_TENSOR(3,3)

      SURF_STRESS3 = LOCAL_STRESS3
      IF (POINT_INDEX < 1) RETURN
      CALL GET_SURFACE_BASIS ( SURF_INDEX, SURF_BASIS )

      LOCAL_TENSOR = ZERO
      LOCAL_TENSOR(1,1) = LOCAL_STRESS3(1)
      LOCAL_TENSOR(2,2) = LOCAL_STRESS3(2)
      LOCAL_TENSOR(1,2) = LOCAL_STRESS3(3)
      LOCAL_TENSOR(2,1) = LOCAL_STRESS3(3)

      SURF_TENSOR = MATMUL(SURF_BASIS, MATMUL(LOCAL_TENSOR, TRANSPOSE(SURF_BASIS)))
      SURF_STRESS3(1) = SURF_TENSOR(1,1)
      SURF_STRESS3(2) = SURF_TENSOR(2,2)
      SURF_STRESS3(3) = SURF_TENSOR(1,2)

      END SUBROUTINE GET_SURFACE_STRESS3

!----------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE GET_SURFACE_BASIS ( SURF_INDEX, SURF_BASIS )

      INTEGER(LONG), INTENT(IN)       :: SURF_INDEX
      REAL(DOUBLE), INTENT(OUT)       :: SURF_BASIS(3,3)

      CHARACTER(8*BYTE)               :: NMODE

      SURF_BASIS = ZERO
      NMODE = GP_SURFACE_NORMAL_MODE(SURF_INDEX)

      IF (NMODE(1:1) == 'X') THEN
         SURF_BASIS(1,2) = ONE
         SURF_BASIS(2,3) = ONE
         SURF_BASIS(3,1) = ONE
      ELSE IF (NMODE(1:1) == 'Y') THEN
         SURF_BASIS(1,3) = ONE
         SURF_BASIS(2,1) = ONE
         SURF_BASIS(3,2) = ONE
      ELSE
         SURF_BASIS(1,1) = ONE
         SURF_BASIS(2,2) = ONE
         SURF_BASIS(3,3) = ONE
      ENDIF

      END SUBROUTINE GET_SURFACE_BASIS

!----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE) FUNCTION GET_TRIA_AREA ( ISTART )

      INTEGER(LONG), INTENT(IN)       :: ISTART
      REAL(DOUBLE)                    :: X1(3), X2(3), X3(3), C(3)

      CALL GET_GRID_BASIC_COORDS ( GID_OUT_ARRAY(ISTART,2), X1 )
      CALL GET_GRID_BASIC_COORDS ( GID_OUT_ARRAY(ISTART,3), X2 )
      CALL GET_GRID_BASIC_COORDS ( GID_OUT_ARRAY(ISTART,4), X3 )
      C(1) = (X2(2)-X1(2))*(X3(3)-X1(3)) - (X2(3)-X1(3))*(X3(2)-X1(2))
      C(2) = (X2(3)-X1(3))*(X3(1)-X1(1)) - (X2(1)-X1(1))*(X3(3)-X1(3))
      C(3) = (X2(1)-X1(1))*(X3(2)-X1(2)) - (X2(2)-X1(2))*(X3(1)-X1(1))
      GET_TRIA_AREA = 0.5D0 * DSQRT(C(1)*C(1) + C(2)*C(2) + C(3)*C(3))

      END FUNCTION GET_TRIA_AREA

!----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE) FUNCTION GET_QUAD_AREA ( ISTART )

      INTEGER(LONG), INTENT(IN)       :: ISTART
      REAL(DOUBLE)                    :: X1(3), X2(3), X3(3), X4(3)

      CALL GET_GRID_BASIC_COORDS ( GID_OUT_ARRAY(ISTART,2), X1 )
      CALL GET_GRID_BASIC_COORDS ( GID_OUT_ARRAY(ISTART,3), X2 )
      CALL GET_GRID_BASIC_COORDS ( GID_OUT_ARRAY(ISTART,4), X3 )
      CALL GET_GRID_BASIC_COORDS ( GID_OUT_ARRAY(ISTART,5), X4 )
      GET_QUAD_AREA = TRI_AREA_FROM_XYZ ( X1, X2, X3 ) + TRI_AREA_FROM_XYZ ( X1, X3, X4 )

      END FUNCTION GET_QUAD_AREA

!----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE) FUNCTION TRI_AREA_FROM_XYZ ( X1, X2, X3 )

      REAL(DOUBLE), INTENT(IN)        :: X1(3), X2(3), X3(3)
      REAL(DOUBLE)                    :: C(3)

      C(1) = (X2(2)-X1(2))*(X3(3)-X1(3)) - (X2(3)-X1(3))*(X3(2)-X1(2))
      C(2) = (X2(3)-X1(3))*(X3(1)-X1(1)) - (X2(1)-X1(1))*(X3(3)-X1(3))
      C(3) = (X2(1)-X1(1))*(X3(2)-X1(2)) - (X2(2)-X1(2))*(X3(1)-X1(1))
      TRI_AREA_FROM_XYZ = 0.5D0 * DSQRT(C(1)*C(1) + C(2)*C(2) + C(3)*C(3))

      END FUNCTION TRI_AREA_FROM_XYZ

!----------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE GET_GRID_BASIC_COORDS ( GRID_NUM, XYZ )

      INTEGER(LONG), INTENT(IN)       :: GRID_NUM
      REAL(DOUBLE), INTENT(OUT)       :: XYZ(3)

      CHARACTER(32*BYTE)              :: SUBR_NAME = 'WRITE_OGS1_SURFACE_STRESS'
      INTEGER(LONG)                   :: IGRID

      CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, SIZE(GRID_ID), GRID_ID, GRID_NUM, IGRID )
      IF (IGRID > 0) THEN
         XYZ(1) = RGRID(IGRID,1)
         XYZ(2) = RGRID(IGRID,2)
         XYZ(3) = RGRID(IGRID,3)
      ELSE
         XYZ = ZERO
      ENDIF

      END SUBROUTINE GET_GRID_BASIC_COORDS

!----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG) FUNCTION FIND_INT ( VALUE, ARRAY, NUSED )

      INTEGER(LONG), INTENT(IN)       :: VALUE
      INTEGER(LONG), INTENT(IN)       :: ARRAY(:)
      INTEGER(LONG), INTENT(IN)       :: NUSED

      INTEGER(LONG)                   :: I

      FIND_INT = 0
      DO I=1,NUSED
         IF (ARRAY(I) == VALUE) THEN
            FIND_INT = I
            EXIT
         ENDIF
      ENDDO

      END FUNCTION FIND_INT

!----------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE WRITE_OGS1_F06_ROW ( GRID_NUM, ELEM_NUM, FIBER, VALUES )

      INTEGER(LONG), INTENT(IN)       :: GRID_NUM
      INTEGER(LONG), INTENT(IN)       :: ELEM_NUM
      CHARACTER(LEN=*), INTENT(IN)    :: FIBER
      REAL(DOUBLE), INTENT(IN)        :: VALUES(8)

      CHARACTER(160*BYTE)             :: LINE_BUF
      INTEGER(LONG)                   :: POS, J

      LINE_BUF = ' '
      IF (GRID_NUM > 0) THEN
         CALL FAST_FMT_I8_RJ ( GRID_NUM, LINE_BUF(2:9) )
         CALL FAST_FMT_I8_RJ ( ELEM_NUM, LINE_BUF(11:18) )
      ENDIF
      LINE_BUF(24:26) = FIBER(1:3)
      POS = 31
      DO J=1,3
         CALL FAST_FMT_F06_E14_6 ( VALUES(J), LINE_BUF(POS:POS+13) )
         POS = POS + 14
      ENDDO
      WRITE(LINE_BUF(73:81),'(F9.4)') VALUES(4)
      POS = 83
      DO J=5,8
         CALL FAST_FMT_F06_E14_6 ( VALUES(J), LINE_BUF(POS:POS+13) )
         POS = POS + 14
      ENDDO
      WRITE(F06,'(A)') TRIM(LINE_BUF)

      END SUBROUTINE WRITE_OGS1_F06_ROW

      END SUBROUTINE WRITE_OGS1_SURFACE_STRESS

!==============================================================================
      SUBROUTINE WRITE_STRESS_I8_PLUS_R14_LINE ( NLEAD, IDVAL, VALUES, NVALS )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE FAST_OUTPUT_FORMATTERS, ONLY:  FAST_FMT_F06_E14_6, FAST_FMT_I8_RJ

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: NLEAD, IDVAL, NVALS
      REAL(DOUBLE), INTENT(IN)        :: VALUES(NVALS)

      CHARACTER(160*BYTE)             :: LINE_BUF
      CHARACTER(8*BYTE)               :: ID_TEXT
      CHARACTER(14*BYTE)              :: VAL_TEXT
      INTEGER(LONG)                   :: I, POS

      LINE_BUF = ' '
      CALL FAST_FMT_I8_RJ ( IDVAL, ID_TEXT )
      LINE_BUF(NLEAD+1:NLEAD+8) = ID_TEXT
      POS = NLEAD + 9
      DO I=1,NVALS
         CALL FAST_FMT_F06_E14_6 ( VALUES(I), VAL_TEXT )
         LINE_BUF(POS:POS+13) = VAL_TEXT
         POS = POS + 14
      ENDDO
      WRITE(F06,'(A)') LINE_BUF(1:POS-1)

      END SUBROUTINE WRITE_STRESS_I8_PLUS_R14_LINE

! ##################################################################################################################################

      SUBROUTINE WRITE_STRESS_SOLID_CENTER_LINE ( EID, VALUES, NVALS )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE FAST_OUTPUT_FORMATTERS, ONLY:  FAST_FMT_F06_E14_6, FAST_FMT_I8_RJ

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: EID, NVALS
      REAL(DOUBLE), INTENT(IN)        :: VALUES(NVALS)

      CHARACTER(139*BYTE)             :: LINE_BUF
      CHARACTER(8*BYTE)               :: ID_TEXT
      CHARACTER(14*BYTE)              :: VAL_TEXT
      INTEGER(LONG)                   :: I, POS

      LINE_BUF = ' '
      CALL FAST_FMT_I8_RJ ( EID, ID_TEXT )
      LINE_BUF(2:9) = ID_TEXT
      LINE_BUF(12:19) = 'CENTER  '
      POS = 28
      DO I=1,NVALS
         CALL FAST_FMT_F06_E14_6 ( VALUES(I), VAL_TEXT )
         LINE_BUF(POS:POS+13) = VAL_TEXT
         POS = POS + 14
      ENDDO
      WRITE(F06,'(A)') LINE_BUF(1:POS-1)

      END SUBROUTINE WRITE_STRESS_SOLID_CENTER_LINE

! ##################################################################################################################################

      SUBROUTINE WRITE_STRESS_SOLID_GRID_LINE ( GRID_ID, VALUES, NVALS )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE FAST_OUTPUT_FORMATTERS, ONLY:  FAST_FMT_F06_E14_6, FAST_FMT_I8_RJ

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: GRID_ID, NVALS
      REAL(DOUBLE), INTENT(IN)        :: VALUES(NVALS)

      CHARACTER(139*BYTE)             :: LINE_BUF
      CHARACTER(8*BYTE)               :: ID_TEXT
      CHARACTER(14*BYTE)              :: VAL_TEXT
      INTEGER(LONG)                   :: I, POS

      LINE_BUF = ' '
      LINE_BUF(12:14) = 'GRD'
      CALL FAST_FMT_I8_RJ ( GRID_ID, ID_TEXT )
      LINE_BUF(15:22) = ID_TEXT
      POS = 28
      DO I=1,NVALS
         CALL FAST_FMT_F06_E14_6 ( VALUES(I), VAL_TEXT )
         LINE_BUF(POS:POS+13) = VAL_TEXT
         POS = POS + 14
      ENDDO
      WRITE(F06,'(A)') LINE_BUF(1:POS-1)

      END SUBROUTINE WRITE_STRESS_SOLID_GRID_LINE

! ##################################################################################################################################

      SUBROUTINE WRITE_STRESS_CSHEAR_PAIR_LINE ( EID1, VALUES1, EID2, VALUES2, HAS_SECOND )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE FAST_OUTPUT_FORMATTERS, ONLY:  FAST_FMT_F06_E14_6, FAST_FMT_I8_RJ

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: EID1, EID2
      REAL(DOUBLE), INTENT(IN)        :: VALUES1(3), VALUES2(3)
      LOGICAL, INTENT(IN)             :: HAS_SECOND

      CHARACTER(128*BYTE)             :: LINE_BUF
      CHARACTER(8*BYTE)               :: ID_TEXT
      CHARACTER(14*BYTE)              :: VAL_TEXT
      INTEGER(LONG)                   :: IVAL, POS

      LINE_BUF = ' '
      POS = 2
      CALL FAST_FMT_I8_RJ ( EID1, ID_TEXT )
      LINE_BUF(POS:POS+7) = ID_TEXT
      POS = POS + 8
      DO IVAL=1,3
         CALL FAST_FMT_F06_E14_6 ( VALUES1(IVAL), VAL_TEXT )
         LINE_BUF(POS:POS+13) = VAL_TEXT
         POS = POS + 14
      ENDDO

      IF (HAS_SECOND) THEN
         LINE_BUF(POS:POS+13) = '              '
         POS = POS + 14
         CALL FAST_FMT_I8_RJ ( EID2, ID_TEXT )
         LINE_BUF(POS:POS+7) = ID_TEXT
         POS = POS + 8
         DO IVAL=1,3
            CALL FAST_FMT_F06_E14_6 ( VALUES2(IVAL), VAL_TEXT )
            LINE_BUF(POS:POS+13) = VAL_TEXT
            POS = POS + 14
         ENDDO
      ENDIF

      WRITE(F06,'(A)') LINE_BUF(1:POS-1)

      END SUBROUTINE WRITE_STRESS_CSHEAR_PAIR_LINE

! ##################################################################################################################################

      SUBROUTINE WRITE_STRESS_ELAS_GROUP_LINE ( IBEG, IEND )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, OGEL
      USE FAST_OUTPUT_FORMATTERS, ONLY:  FAST_FMT_F06_E14_6, FAST_FMT_I8_RJ

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: IBEG, IEND

      CHARACTER(160*BYTE)             :: LINE_BUF
      CHARACTER(8*BYTE)               :: ID_TEXT
      CHARACTER(14*BYTE)              :: VAL_TEXT
      INTEGER(LONG)                   :: IROW, POS

      LINE_BUF = ' '
      POS = 1
      DO IROW=IBEG,IEND
         CALL FAST_FMT_I8_RJ ( EID_OUT_ARRAY(IROW,1), ID_TEXT )
         CALL FAST_FMT_F06_E14_6 ( OGEL(IROW,1), VAL_TEXT )
         LINE_BUF(POS:POS)       = ' '
         LINE_BUF(POS+1:POS+8)   = ID_TEXT
         LINE_BUF(POS+9:POS+22)  = VAL_TEXT
         POS = POS + 23
      ENDDO

      WRITE(F06,'(A)') LINE_BUF(1:POS-1)

      END SUBROUTINE WRITE_STRESS_ELAS_GROUP_LINE

!==============================================================================
      SUBROUTINE GET_SPRING_OP2_ELEMENT_TYPE(ELEMENT_TYPE)
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY      :  ERR
      USE MODEL_STUF, ONLY  :  TYPE
      IMPLICIT NONE
      INTEGER(LONG)    :: ELEMENT_TYPE ! the OP2 flag for the element
      !               12345
      ! 11 : CELAS1 - ELAS1
      ! 12 : CELAS2
      ! 13 : CELAS3
      ! 14 : CELAS4
      IF (TYPE(5:5) == "1") THEN
          ELEMENT_TYPE = 11
      ELSE IF (TYPE(5:5) == "2") THEN
          ELEMENT_TYPE = 12
      ELSE IF (TYPE(5:5) == "3") THEN
          ELEMENT_TYPE = 13
      ELSE IF (TYPE(5:5) == "4") THEN
          ELEMENT_TYPE = 14
      ELSE
 42       FORMAT("TYPE(4:4)=",A," TYPE(5:5)=",A," TYPE(6:6)=",A)
          WRITE(ERR,42) TYPE(4:4),TYPE(5:5),TYPE(6:6)
          ELEMENT_TYPE = -1
      ENDIF
      END SUBROUTINE GET_SPRING_OP2_ELEMENT_TYPE

!==============================================================================
      SUBROUTINE TRANSFORM_SHELL_OUTPUT_ROW_TO_BASIC ( TE_LOCAL, VALUES )

      USE PENTIUM_II_KIND, ONLY       :  DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: TE_LOCAL(3,3)
      REAL(DOUBLE), INTENT(INOUT)     :: VALUES(10)

      REAL(DOUBLE)                    :: TBL(3,3)
      REAL(DOUBLE)                    :: SIG_LOCAL(3,3)
      REAL(DOUBLE)                    :: SIG_BASIC(3,3)
      REAL(DOUBLE)                    :: AVG, DIFF, RAD, SMAJ, SMIN, SXYMAX, VM
      REAL(DOUBLE), PARAMETER         :: RAD2DEG = 180.0D0 / 3.1415926535897932384626433832795D0

      IF (MAXVAL(ABS(TE_LOCAL)) <= ZERO) RETURN

      TBL(1:3,1:3) = TRANSPOSE(TE_LOCAL(1:3,1:3))

      SIG_LOCAL = ZERO
      SIG_LOCAL(1,1) = VALUES(2)
      SIG_LOCAL(2,2) = VALUES(3)
      SIG_LOCAL(1,2) = VALUES(4)
      SIG_LOCAL(2,1) = VALUES(4)
      SIG_LOCAL(1,3) = VALUES(9)
      SIG_LOCAL(3,1) = VALUES(9)
      SIG_LOCAL(2,3) = VALUES(10)
      SIG_LOCAL(3,2) = VALUES(10)

      SIG_BASIC = MATMUL( TBL, MATMUL(SIG_LOCAL, TRANSPOSE(TBL)) )

      VALUES(2)  = SIG_BASIC(1,1)
      VALUES(3)  = SIG_BASIC(2,2)
      VALUES(4)  = SIG_BASIC(1,2)
      VALUES(9)  = SIG_BASIC(1,3)
      VALUES(10) = SIG_BASIC(2,3)

      AVG = 0.5D0 * (VALUES(2) + VALUES(3))
      DIFF = 0.5D0 * (VALUES(2) - VALUES(3))
      RAD = DSQRT(DIFF*DIFF + VALUES(4)*VALUES(4))
      SMAJ = AVG + RAD
      SMIN = AVG - RAD
      SXYMAX = 0.5D0 * DABS(SMAJ - SMIN)
      VM = DSQRT(VALUES(2)*VALUES(2) - VALUES(2)*VALUES(3) + VALUES(3)*VALUES(3) + 3.0D0*VALUES(4)*VALUES(4))

      VALUES(5) = 0.5D0 * DATAN2(2.0D0*VALUES(4), VALUES(2) - VALUES(3)) * RAD2DEG
      VALUES(6) = SMAJ
      VALUES(7) = SMIN
      VALUES(8) = VM

      END SUBROUTINE TRANSFORM_SHELL_OUTPUT_ROW_TO_BASIC
