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

      MODULE LINK9_STUFF

! Grid point and element solution variables for data recovery LINK9

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      IMPLICIT NONE

      SAVE

      CHARACTER( 1*BYTE), ALLOCATABLE :: MSPRNT(:,:)           ! Flags for whether to print margins of safety for ROD, BAR
      CHARACTER( 4*BYTE), ALLOCATABLE :: FTNAME(:)             ! Stress failure index name output with stresses/strains
      LOGICAL                         :: SMART_OUTPUT_MODE = .FALSE.
      LOGICAL                         :: WRITE_NEU_GEOM    = .FALSE.
      LOGICAL                         :: WRITE_NEU_DISP    = .FALSE.
      LOGICAL                         :: WRITE_NEU_OLOA    = .FALSE.
      LOGICAL                         :: WRITE_NEU_SPCF    = .FALSE.
      LOGICAL                         :: WRITE_NEU_MPCF    = .FALSE.
      LOGICAL                         :: WRITE_NEU_ELFO    = .FALSE.
      LOGICAL                         :: WRITE_NEU_STRE    = .FALSE.
      LOGICAL                         :: WRITE_NEU_STRN    = .FALSE.
      LOGICAL                         :: RSA_ELFE_CAPTURE  = .FALSE.

      INTEGER(LONG)                   :: MAXREQ                ! Max number of rows needed for array OGEL
      INTEGER(LONG)                   :: RSA_ELFE_NUM_ROWS = 0

      INTEGER(LONG)     , ALLOCATABLE :: GID_OUT_ARRAY(:,:)    ! Array of integer grid no's for some output in LINK9

      INTEGER(LONG)     , ALLOCATABLE :: EID_OUT_ARRAY(:,:)    ! Array of elem no's (col 1) and num of plies (col 2) for that elem
!                                                                that are printed with certain outputs IN LINK9

      INTEGER(LONG)     , ALLOCATABLE :: POLY_FIT_ERR_INDEX(:)! Index num for POLY_FIT_ERR (i.e. which of the 1 through 9 stress
!                                                                or strain values has the largest error in polynomial fit


      REAL(DOUBLE)      , ALLOCATABLE :: OGEL(:,:)             ! Master array for holding outputs in LINK9 until they are printed
! --- cbeam_stations begin --- !
      REAL(DOUBLE)      , ALLOCATABLE :: CBEAM_XL_OUT(:)       ! x/L station metadata for beam-style output rows in LINK9
! --- cbeam_stations end --- !
      REAL(DOUBLE)      , ALLOCATABLE :: SHELL_OUT_TE(:,:,:)   ! Shell output basis per stored stress/strain point row.
!                                                                Dimensions are (3,3,MAXREQ) and TE maps basic -> shell local.
      LOGICAL           , ALLOCATABLE :: SHELL_STRESS_IN_LOCAL(:) ! True when shell stress rows need local -> surface projection.

      REAL(DOUBLE)      , ALLOCATABLE :: POLY_FIT_ERR(:)       ! Array of polynom fit errors for elems that extrapolate stress or
!                                                                strain values from one set of output points to another
      REAL(DOUBLE)                    :: RSA_MODE_SCALE = 0.0D0 ! Modal RSA scale applied to non-grid result recovery in LINK9
      REAL(DOUBLE)      , ALLOCATABLE :: RSA_ELFE_SUMSQ(:)     ! Row-wise SRSS accumulator for RSA ELFORCE(ENGR) summaries
      REAL(DOUBLE)      , ALLOCATABLE :: RSA_ELFE_SUMABS(:)    ! Row-wise ABS accumulator for RSA ELFORCE(ENGR) summaries
      CHARACTER(132*BYTE), ALLOCATABLE:: RSA_ELFE_DESC(:)      ! Row descriptors for RSA ELFORCE(ENGR) summaries

      CONTAINS

      PURE SUBROUTINE SHELL_POINT_TENSOR_TO_FRAME ( POINT_BASIS, OUTPUT_BASIS, POINT_VALUES, OUTPUT_TENSOR )
! Basis vectors are rows (basic -> local). POINT_VALUES must be recovered at
! the SAME point as POINT_BASIS. An identity OUTPUT_BASIS returns the full
! physical basic-coordinate tensor; the OP2 convention uses the center basis.
! This output-only operation does not change shell operators or stiffness.
      REAL(DOUBLE), INTENT(IN)  :: POINT_BASIS(3,3), OUTPUT_BASIS(3,3), POINT_VALUES(3)
      REAL(DOUBLE), INTENT(OUT) :: OUTPUT_TENSOR(3,3)
      REAL(DOUBLE) :: LOCAL_TENSOR(3,3), ROTATION(3,3)

      LOCAL_TENSOR = 0.0_DOUBLE
      LOCAL_TENSOR(1,1) = POINT_VALUES(1)
      LOCAL_TENSOR(2,2) = POINT_VALUES(2)
      LOCAL_TENSOR(1,2) = POINT_VALUES(3)
      LOCAL_TENSOR(2,1) = POINT_VALUES(3)
      ROTATION = MATMUL(OUTPUT_BASIS, TRANSPOSE(POINT_BASIS))
      OUTPUT_TENSOR = MATMUL(ROTATION, MATMUL(LOCAL_TENSOR, TRANSPOSE(ROTATION)))
      END SUBROUTINE SHELL_POINT_TENSOR_TO_FRAME

      END MODULE LINK9_STUFF
