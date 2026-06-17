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

      SUBROUTINE BD_PBEAM ( CARD, LARGE_FLD_INP )

! Processes PBEAM Bulk Data Cards.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE PARAMS, ONLY                :  EPSIL
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, BEAMTOR, FATAL_ERR, IERRFL, JCARD_LEN, JF, LPBEAM, NPBEAM, MPBEAM_STATIONS, WARN_ERR
      USE CONSTANTS_1, ONLY           :  ZERO
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  PBEAM, PBEAM_NSTATIONS, PBEAM_XL, PBEAM_RPROPS, RPBEAM
      USE PARAMS, ONLY                :  SUPINFO

      USE BD_PBEAM_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME   =   'BD_PBEAM'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD               ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)          ! The 10 fields of 8 characters making up CARD
      CHARACTER(LEN(JCARD))           :: CHRINP             ! Character field read from CARD
      CHARACTER(LEN(CARD))            :: RAW_CONT           ! Raw continuation text for NX free-field station detection
      CHARACTER(LEN(JCARD))           :: ID                 ! Property ID for this PBEAM
      CHARACTER(LEN(JCARD))           :: NAME               ! Char name for output error purposes

      INTEGER(LONG)                   :: ICONT       = 0    ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR        = 0    ! Error indicator
      INTEGER(LONG)                   :: J                  ! DO loop index
      INTEGER(LONG)                   :: MATERIAL_ID = 0    ! Material ID (field 3 of this property card)
      INTEGER(LONG)                   :: PROPERTY_ID = 0    ! Property ID (field 2 of this property card)
! --- CBEAM_standard begin --- !
      INTEGER(LONG)                   :: ISTATION     = 0    ! Count of stored NX-style continuation stations for this PBEAM
      INTEGER(LONG)                   :: SO_FIELD     = 0
      INTEGER(LONG)                   :: XL_FIELD     = 0
      INTEGER(LONG)                   :: PROP_FIELD0  = 0
! --- CBEAM_standard end --- !


      REAL(DOUBLE)                    :: AREA_A      = ZERO ! Cross sectional area at end A
      REAL(DOUBLE)                    :: I1_A        = ZERO ! Moment of inertia, plane 1 at end A
      REAL(DOUBLE)                    :: I2_A        = ZERO ! Moment of inertia, plane 2 at end A
      REAL(DOUBLE)                    :: I12_A       = ZERO ! Product of inertia at end A
      REAL(DOUBLE)                    :: JTOR_A      = ZERO ! Torsional constantr at end A
      REAL(DOUBLE)                    :: NSM_A       = ZERO ! Nonstructural mass at end A

      REAL(DOUBLE)                    :: AREA        = ZERO ! Cross sectional area at any location along beam
      REAL(DOUBLE)                    :: I1          = ZERO ! Moment of inertia, plane 1 at any location along beam
      REAL(DOUBLE)                    :: I2          = ZERO ! Moment of inertia, plane 2 at any location along beam
      REAL(DOUBLE)                    :: I12         = ZERO ! Product of inertia at any location along beam
      REAL(DOUBLE)                    :: JTOR        = ZERO ! Torsional constantr at any location along beam
      REAL(DOUBLE)                    :: NSM         = ZERO ! Nonstructural mass at any location along beam
! --- CBEAM_standard begin --- !
      REAL(DOUBLE)                    :: STATION_XL  = ZERO ! Current station x/L value read from continuation chain
      LOGICAL                         :: CONT_IS_STATION = .FALSE.
      LOGICAL                         :: STATION_WARNED  = .FALSE.
! --- CBEAM_standard end --- !



! **********************************************************************************************************************************
! PBEAM Bulk Data Card routine

!  FIELD          ITEM                                                ARRAY ELEMENT
!  -----          ----                                             ------------------
!    2      Property ID                                 PID          PBEAM(npbeam, 1)
!    3      Material ID                                 MID          PBEAM(npbeam, 2)
!    4      Cross sectional area          : end A       A(A)        RPBEAM(npbeam, 1)
!    5      Area moment of inertia 1      :   "         I1(A)       RPBEAM(npbeam, 2)
!    6      Area moment of inertia 2      :   "         I2(A)       RPBEAM(npbeam, 3)
!    7      Area product of inertia 12    :   "         I12(A)      RPBEAM(npbeam, 4)
!    8      Torsion constant              :   "         J(A)        RPBEAM(npbeam, 5)
!    9      Non structural mass           :   "         NSM(A)      RPBEAM(npbeam, 6)

! on mandatory 2nd card:
!    2      Stress coefficient            :   "         C1(A)       RPBEAM(npbeam, 7)
!    3      Stress coefficient            :   "         C2(A)       RPBEAM(npbeam, 8)
!    4      Stress coefficient            :   "         D1(A)       RPBEAM(npbeam, 9)
!    5      Stress coefficient            :   "         D2(A)       RPBEAM(npbeam,10)
!    6      Stress coefficient            :   "         E1(A)       RPBEAM(npbeam,11)
!    7      Stress coefficient            :   "         E2(A)       RPBEAM(npbeam,12)
!    8      Stress coefficient            :   "         F1(A)       RPBEAM(npbeam,13)
!    9      Stress coefficient            :   "         F2(A)       RPBEAM(npbeam,14)

! on mandatory 3rd card:
!    2      Stress output request option                S0           PBEAM(npbeam, 3)
!    3      Loc of next set of data (X/XB)              X/XB        RPBEAM(NPBEAM,15)
!    4      Cross sectional area          : end B       A(B)        RPBEAM(npbeam,16)
!    5      Area moment of inertia 1      :   "         I1(B)       RPBEAM(npbeam,17)
!    6      Area moment of inertia 2      :   "         I2(B)       RPBEAM(npbeam,18)
!    7      Area produc of inertia 12     :   "         I12(B)      RPBEAM(npbeam,19)
!    8      Torsion constant              :   "         J(B)        RPBEAM(npbeam,20)
!    9      Non structural mass           :   "         NSM(B)      RPBEAM(npbeam,21)

! on optional  4th card:
!    2      Stress coefficient            :   "         C1(B)       RPBEAM(npbeam,22)
!    3      Stress coefficient            :   "         C2(B)       RPBEAM(npbeam,23)
!    4      Stress coefficient            :   "         D1(B)       RPBEAM(npbeam,24)
!    5      Stress coefficient            :   "         D2(B)       RPBEAM(npbeam,25)
!    6      Stress coefficient            :   "         E1(B)       RPBEAM(npbeam,26)
!    7      Stress coefficient            :   "         E2(B)       RPBEAM(npbeam,27)
!    8      Stress coefficient            :   "         F1(B)       RPBEAM(npbeam,28)
!    9      Stress coefficient            :   "         F2(B)       RPBEAM(npbeam,29)

! on optional  5th card:
!    2      Shear factor for plane 1                    K1          RPBEAM(npbeam,30)
!    3      Shear factor for plane 2                    K2          RPBEAM(npbeam,31)
!    4      Shear relief coeff due to taper for plane 1 S1          RPBEAM(npbeam,32)
!    5      Shear relief coeff due to taper for plane 2 S2          RPBEAM(npbeam,33)
!    6      NSM MOI/length about NSM C.G. at end A      NSI(A)      RPBEAM(npbeam,34)
!    7      NSM MOI/length about NSM C.G. at end B      NSI(B)      RPBEAM(npbeam,35)
!    8      Warping coefficient for end A               CW(A)       RPBEAM(npbeam,36)
!    9      Warping coefficient for end B               CW(B)       RPBEAM(npbeam,37)

! on optional  6th card:
!    2      y coord of C.G. of NSM at end A             M1(A)       RPBEAM(npbeam,38)
!    3      z coord of C.G. of NSM at end A             M2(A)       RPBEAM(npbeam,39)
!    4      y coord of C.G. of NSM at end B             M1(B)       RPBEAM(npbeam,40)
!    5      z coord of C.G. of NSM at end B             M2(B)       RPBEAM(npbeam,41)
!    6      y coord of neutral axis for end A           N1(A)       RPBEAM(npbeam,42)
!    7      z coord of neutral axis for end A           N2(A)       RPBEAM(npbeam,43)
!    8      y coord of neutral axis for end B           N1(B)       RPBEAM(npbeam,44)
!    9      z coord of neutral axis for end B           N2(B)       RPBEAM(npbeam,45)


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Increment NPBEAM

      NPBEAM = NPBEAM + 1
! --- cbeam_tapered_add begin --- !
      PBEAM_NSTATIONS(NPBEAM) = 1
      PBEAM_XL(NPBEAM,1) = ZERO
! --- cbeam_tapered_add end --- !

! Read and check data on parent card

      NAME = JCARD(1)
      ID   = JCARD(2)
      CALL I4FLD ( JCARD(2), JF(2), PROPERTY_ID )          ! Read property ID and enter into array PBEAM
      IF (IERRFL(2) == 'N') THEN
         DO J=1,NPBEAM-1
            IF (PROPERTY_ID == PBEAM(J,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),PROPERTY_ID
               WRITE(F06,1145) JCARD(1),PROPERTY_ID
               EXIT
             ENDIF
         ENDDO
         PBEAM(NPBEAM,1) = PROPERTY_ID
      ENDIF

      CALL I4FLD ( JCARD(3), JF(3), MATERIAL_ID )          ! Read material ID and enter into array PBEAM
      IF (IERRFL(3) == 'N') THEN
         IF (MATERIAL_ID <= 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', MATERIAL_ID
            WRITE(F06,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', MATERIAL_ID
         ELSE
            PBEAM(NPBEAM,2) = MATERIAL_ID
         ENDIF
      ENDIF

      DO J = 1,6                                           ! Read real property values in fields 4-8
         CALL R8FLD ( JCARD(J+3), JF(J+3), RPBEAM(NPBEAM,J) )
      ENDDO
      IF (IERRFL(4)  == 'N') THEN
         AREA_A = RPBEAM(NPBEAM,1)
      ENDIF
      IF (IERRFL(5)  == 'N') THEN
         I1_A   = RPBEAM(NPBEAM,2)
      ENDIF
      IF (IERRFL(6)  == 'N') THEN
         I2_A   = RPBEAM(NPBEAM,3)
      ENDIF
      IF (IERRFL(7)  == 'N') THEN
         I12_A  = RPBEAM(NPBEAM,4)
      ENDIF
      IF (IERRFL(8)  == 'N') THEN
         JTOR_A = RPBEAM(NPBEAM,5)
      ENDIF
      IF (IERRFL(9)  == 'N') THEN
         NSM_A  = RPBEAM(NPBEAM,6)
      ENDIF

! Call subr to check sensibility of I1, I2, I12 combinations

      CALL CHECK_BAR_MOIs ( 'PBEAM', ID, I1_A, I2_A, I12_A, IERR )
      RPBEAM(NPBEAM,2) = I1_A
      RPBEAM(NPBEAM,3) = I2_A
      RPBEAM(NPBEAM,4) = I12_A
      I1 = I1_A
      I2 = I2_A
      I12 = I12_A
      PBEAM_RPROPS(NPBEAM,1,1) = AREA_A
      PBEAM_RPROPS(NPBEAM,1,2) = I1_A
      PBEAM_RPROPS(NPBEAM,1,3) = I2_A
      PBEAM_RPROPS(NPBEAM,1,4) = I12_A
      PBEAM_RPROPS(NPBEAM,1,5) = JTOR_A
      PBEAM_RPROPS(NPBEAM,1,6) = NSM_A
      IF (IERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
      ENDIF

! Read and check data on mandatory 2nd card (needed even if all fields are blank since 3rd cont is mandatory):

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN
         DO J = 7,14                                       ! Read real property values in fields 2-9 of 2nd card
            CALL R8FLD ( JCARD(J-5), JF(J-5), RPBEAM(NPBEAM,J) )
         ENDDO
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-9
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields
      ELSE
         DO J = 7,14
            RPBEAM(NPBEAM,J) = ZERO
         ENDDO
      ENDIF

! --- cbeam_tapered_add begin --- !
! Read station continuation chain. Phase-1 beam redevelopment stores every x/L station explicitly for NX-oriented CBEAM work,
! while still preserving the legacy RPBEAM end-B snapshot in cols 15:29 using the final station that is read.
! --- cbeam_tapered_add end --- !
! Read and check data on mandatory 3rd card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN

         AREA = AREA_A
         I1   = I1_A
         I2   = I2_A
         I12  = I12_A
         JTOR = JTOR_A
         NSM  = NSM_A
         ISTATION = 1

station_loop: DO
            RAW_CONT = CARD
            CALL LEFT_ADJ_BDFLD ( RAW_CONT )
            CALL CHAR_FLD ( JCARD(2), JF(2), CHRINP )
            CALL LEFT_ADJ_BDFLD ( CHRINP )
            SO_FIELD    = 2
            XL_FIELD    = 3
            PROP_FIELD0 = 4
            CONT_IS_STATION = ((CHRINP(1:3) == 'YES') .OR. (CHRINP(1:2) == 'NO'))
            IF (.NOT. CONT_IS_STATION) THEN
               CALL CHAR_FLD ( JCARD(3), JF(3), CHRINP )
               CALL LEFT_ADJ_BDFLD ( CHRINP )
               SO_FIELD    = 3
               XL_FIELD    = 4
               PROP_FIELD0 = 5
               CONT_IS_STATION = ((CHRINP(1:3) == 'YES') .OR. (CHRINP(1:2) == 'NO'))
            ENDIF
            IF (.NOT. CONT_IS_STATION) THEN
               CONT_IS_STATION = ((RAW_CONT(1:5) == '+,YES') .OR. (RAW_CONT(1:4) == '+,NO'))
               IF (CONT_IS_STATION) THEN
                  SO_FIELD    = 3
                  XL_FIELD    = 4
                  PROP_FIELD0 = 5
               ENDIF
            ENDIF
            IF (.NOT. CONT_IS_STATION) EXIT station_loop

            PBEAM(NPBEAM,3) = 1

            CALL R8FLD ( JCARD(XL_FIELD), JF(XL_FIELD), STATION_XL )
            IF (IERRFL(XL_FIELD) == 'N') THEN
               IF (ISTATION < MPBEAM_STATIONS) THEN
                  ISTATION = ISTATION + 1
                  PBEAM_XL(NPBEAM,ISTATION) = STATION_XL
                  PBEAM_NSTATIONS(NPBEAM) = ISTATION
               ELSE IF (.NOT. STATION_WARNED) THEN
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,1197) NAME, ID, MPBEAM_STATIONS
                  WRITE(F06,1197) NAME, ID, MPBEAM_STATIONS
                  STATION_WARNED = .TRUE.
               ENDIF
               RPBEAM(NPBEAM,15) = STATION_XL
            ENDIF

            IF (JCARD(PROP_FIELD0)(1:) /= ' ') THEN
               CALL R8FLD ( JCARD(PROP_FIELD0), JF(PROP_FIELD0), RPBEAM(NPBEAM,16) )
               IF (IERRFL(PROP_FIELD0) == 'N') AREA = RPBEAM(NPBEAM,16)
            ELSE
               RPBEAM(NPBEAM,16) = AREA
            ENDIF

            IF (JCARD(PROP_FIELD0+1)(1:) /= ' ') THEN
               CALL R8FLD ( JCARD(PROP_FIELD0+1), JF(PROP_FIELD0+1), RPBEAM(NPBEAM,17) )
               IF (IERRFL(PROP_FIELD0+1) == 'N') I1 = RPBEAM(NPBEAM,17)
            ELSE
               RPBEAM(NPBEAM,17) = I1
            ENDIF

            IF (JCARD(PROP_FIELD0+2)(1:) /= ' ') THEN
               CALL R8FLD ( JCARD(PROP_FIELD0+2), JF(PROP_FIELD0+2), RPBEAM(NPBEAM,18) )
               IF (IERRFL(PROP_FIELD0+2) == 'N') I2 = RPBEAM(NPBEAM,18)
            ELSE
               RPBEAM(NPBEAM,18) = I2
            ENDIF

            IF (JCARD(PROP_FIELD0+3)(1:) /= ' ') THEN
               CALL R8FLD ( JCARD(PROP_FIELD0+3), JF(PROP_FIELD0+3), RPBEAM(NPBEAM,19) )
               IF (IERRFL(PROP_FIELD0+3) == 'N') I12 = RPBEAM(NPBEAM,19)
            ELSE
               RPBEAM(NPBEAM,19) = I12
            ENDIF

            IF (JCARD(PROP_FIELD0+4)(1:) /= ' ') THEN
               CALL R8FLD ( JCARD(PROP_FIELD0+4), JF(PROP_FIELD0+4), RPBEAM(NPBEAM,20) )
               IF (IERRFL(PROP_FIELD0+4) == 'N') JTOR = RPBEAM(NPBEAM,20)
            ELSE
               RPBEAM(NPBEAM,20) = JTOR
            ENDIF

            IF (JCARD(PROP_FIELD0+5)(1:) /= ' ') THEN
               CALL R8FLD ( JCARD(PROP_FIELD0+5), JF(PROP_FIELD0+5), RPBEAM(NPBEAM,21) )
               IF (IERRFL(PROP_FIELD0+5) == 'N') NSM = RPBEAM(NPBEAM,21)
            ELSE
               RPBEAM(NPBEAM,21) = NSM
            ENDIF

            IF (ISTATION <= MPBEAM_STATIONS) THEN
               PBEAM_RPROPS(NPBEAM,ISTATION,1) = AREA
               PBEAM_RPROPS(NPBEAM,ISTATION,2) = I1
               PBEAM_RPROPS(NPBEAM,ISTATION,3) = I2
               PBEAM_RPROPS(NPBEAM,ISTATION,4) = I12
               PBEAM_RPROPS(NPBEAM,ISTATION,5) = JTOR
               PBEAM_RPROPS(NPBEAM,ISTATION,6) = NSM
            ENDIF

            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )
            CALL CRDERR ( CARD )

            IF (LARGE_FLD_INP == 'N') THEN
               CALL NEXTC  ( CARD, ICONT, IERR )
            ELSE
               CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
               CARD = CHILD
            ENDIF
            CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
            IF (ICONT /= 1) THEN
               DO J = 22,29
                  RPBEAM(NPBEAM,J) = ZERO
               ENDDO
               EXIT station_loop
            ENDIF

            RAW_CONT = CARD
            CALL LEFT_ADJ_BDFLD ( RAW_CONT )
            CALL CHAR_FLD ( JCARD(2), JF(2), CHRINP )
            CALL LEFT_ADJ_BDFLD ( CHRINP )
            CONT_IS_STATION = ((CHRINP(1:3) == 'YES') .OR. (CHRINP(1:2) == 'NO'))
            IF (.NOT. CONT_IS_STATION) THEN
               CALL CHAR_FLD ( JCARD(3), JF(3), CHRINP )
               CALL LEFT_ADJ_BDFLD ( CHRINP )
               CONT_IS_STATION = ((CHRINP(1:3) == 'YES') .OR. (CHRINP(1:2) == 'NO'))
            ENDIF
            IF (.NOT. CONT_IS_STATION) THEN
               CONT_IS_STATION = ((RAW_CONT(1:5) == '+,YES') .OR. (RAW_CONT(1:4) == '+,NO'))
            ENDIF
            IF (CONT_IS_STATION) THEN
               DO J = 22,29
                  RPBEAM(NPBEAM,J) = ZERO
               ENDDO
               CYCLE station_loop
            ENDIF

            DO J = 22,29
               CALL R8FLD ( JCARD(J-20), JF(J-20), RPBEAM(NPBEAM,J) )
            ENDDO
            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )
            CALL CRDERR ( CARD )

            IF (LARGE_FLD_INP == 'N') THEN
               CALL NEXTC  ( CARD, ICONT, IERR )
            ELSE
               CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
               CARD = CHILD
            ENDIF
            IF (ICONT /= 1) EXIT station_loop
            CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         ENDDO station_loop

      ELSE

         PBEAM(NPBEAM,3)  = 0
! --- cbeam_tapered_add begin --- !
! Standard NX-style single-line PBEAM cards can legally stop after the end-A
! property line. When that happens, MYSTRAN must still zero the unused stress,
! shear, and end-B tail storage explicitly so later CBEAM property handoff does
! not see uninitialized values.
         DO J = 7,14
            RPBEAM(NPBEAM,J) = ZERO
         ENDDO
         RPBEAM(NPBEAM,15) = 1.0D0
         RPBEAM(NPBEAM,16) = AREA_A
         RPBEAM(NPBEAM,17) = I1_A
         RPBEAM(NPBEAM,18) = I2_A
         RPBEAM(NPBEAM,19) = I12_A
         RPBEAM(NPBEAM,20) = JTOR_A
         RPBEAM(NPBEAM,21) = NSM_A
         DO J = 22,45
            RPBEAM(NPBEAM,J) = ZERO
         ENDDO
! --- cbeam_tapered_add end --- !

      ENDIF

! Read and check data on optional shear/neutral-axis tail cards, if present after the NX-style station chain.

      IF (ICONT == 1) THEN
         DO J = 30,37
            CALL R8FLD ( JCARD(J-28), JF(J-28), RPBEAM(NPBEAM,J) )
         ENDDO
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )
         CALL CRDERR ( CARD )

         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN
            DO J = 38,45
               CALL R8FLD ( JCARD(J-36), JF(J-36), RPBEAM(NPBEAM,J) )
            ENDDO
            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )
            CALL CRDERR ( CARD )
         ENDIF
      ENDIF

      RETURN

! **********************************************************************************************************************************
 1136 FORMAT(' *ERROR  1136: REQUIRED CONTINUATION FOR ',A,' ID = ',A,' MISSING')

 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

1167 FORMAT(' *ERROR  1167: INCORRECT VALUE IN FIELD ',I2,' OF ',A,A,' = ',A,' MUST BE "YES", "YESA" OR "NO"')

1197 FORMAT(' *WARNING 1197: ',A,' ID = ',A,' HAS MORE THAN ',I8,' STORED STATIONS. EXTRA x/L VALUES WILL BE IGNORED')

1192 FORMAT(' *ERROR  1192: ID IN FIELD ',I3,' OF ',A,A,' MUST BE ',A,' BUT IS = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE BD_PBEAM
