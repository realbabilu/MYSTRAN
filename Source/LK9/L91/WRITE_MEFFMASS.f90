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

      SUBROUTINE WRITE_MEFFMASS

      ! Writes output for modal effective mass
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NVEC, NSUB, MODE_SUBCASE
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, ONE_HUNDRED, PI
      USE PARAMS, ONLY                :  PRTF06, PRTOP2
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL, MEFFMASS
      USE MODEL_STUF, ONLY            :  MEFM_RB_MASS, LABEL, STITLE, TITLE, SCNUM, MEFFMASS_REQ_SUMMARY,                        &
                                         MEFFMASS_REQ_MEFFM, MEFFMASS_REQ_MEFFW, MEFFMASS_REQ_FRACSUM,                            &
                                         MEFFMASS_REQ_SUMMARY_SUB, MEFFMASS_REQ_MEFFM_SUB, MEFFMASS_REQ_MEFFW_SUB,                &
                                         MEFFMASS_REQ_FRACSUM_SUB, MEFMLOC_SUB, MEFMGRID_SUB, MCG, MODEL_MASS,                    &
                                         MODEL_XCG, MODEL_YCG, MODEL_ZCG, GRID_ID, RGRID
      USE PARAMS, ONLY                :  EPSIL, GRDPNT, MEFMCORD, MEFMGRID, MEFMLOC, SUPINFO, WTMASS

      USE WRITE_MEFFMASS_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_MEFFMASS'
      CHARACTER(14*BYTE)              :: CHAR_PCT(6)       ! Character representation of MEFFMASS sum percents of total model mass
      CHARACTER(1*BYTE)               :: IHDR   = 'Y'      ! Indicator of whether to write an output header

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: ISUB
      INTEGER(LONG)                   :: MODE_NUM_OUT
      INTEGER(LONG)                   :: NUM_SUB_LOOPS


      REAL(DOUBLE)                    :: CYCLES            ! Circular frequency of a mode
      REAL(DOUBLE)                    :: EPS1              ! Small number to compare against zero
      REAL(DOUBLE)                    :: ICG(3,3)
      REAL(DOUBLE)                    :: IREF(3,3)
      REAL(DOUBLE)                    :: LOCAL_MEFM_RB_MASS(6,6)
      REAL(DOUBLE)                    :: MEFM_TOTALS(6)    ! Totals for the 6 modal effective masses over all modes
      REAL(DOUBLE)                    :: MODES_PCT(6)      ! Modal mass as % of total mass
      REAL(DOUBLE)                    :: R2
      REAL(DOUBLE)                    :: XREF(3)
      REAL(DOUBLE)                    :: XB(3)
      !LOGICAL                        :: WRITE_F06  ! flag
      !LOGICAL                        :: WRITE_OP2  ! flag
      LOGICAL                         :: IS_LOW_PRECISION  ! Print MPFACTOR, MEFFMASS values with 2 decimal places of accuracy rather than 6
      LOGICAL                         :: WRITE_SUMMARY
      LOGICAL                         :: WRITE_MEFFM
      LOGICAL                         :: WRITE_MEFFW
      LOGICAL                         :: WRITE_FRACSUM
      LOGICAL                         :: HAS_SUBCASE_MAP
      CHARACTER(6*BYTE)               :: LOCAL_MEFMLOC
      INTEGER(LONG)                   :: LOCAL_MEFMGRID
      INTEGER(LONG)                   :: XREF_GRID_ROW




! **********************************************************************************************************************************
      IS_LOW_PRECISION = (DEBUG(174) == 0)
      EPS1 = EPSIL(1)
      HAS_SUBCASE_MAP = ALLOCATED(MODE_SUBCASE)
      NUM_SUB_LOOPS = MAX(1,NSUB)

      DO ISUB=1,NUM_SUB_LOOPS

         IF (ALLOCATED(MEFFMASS_REQ_SUMMARY_SUB)) THEN
            WRITE_SUMMARY = (MEFFMASS_REQ_SUMMARY_SUB(ISUB) == 'Y')
            WRITE_MEFFM   = (MEFFMASS_REQ_MEFFM_SUB(ISUB)   == 'Y')
            WRITE_MEFFW   = (MEFFMASS_REQ_MEFFW_SUB(ISUB)   == 'Y')
            WRITE_FRACSUM = (MEFFMASS_REQ_FRACSUM_SUB(ISUB) == 'Y')
         ELSE
            IF (ISUB > 1) CYCLE
            WRITE_SUMMARY = (MEFFMASS_REQ_SUMMARY == 'Y')
            WRITE_MEFFM   = (MEFFMASS_REQ_MEFFM   == 'Y')
            WRITE_MEFFW   = (MEFFMASS_REQ_MEFFW   == 'Y')
            WRITE_FRACSUM = (MEFFMASS_REQ_FRACSUM == 'Y')
         ENDIF

         IF (.NOT. (WRITE_SUMMARY .OR. WRITE_MEFFM .OR. WRITE_MEFFW .OR. WRITE_FRACSUM)) CYCLE

         LOCAL_MEFMLOC = MEFMLOC
         LOCAL_MEFMGRID = MEFMGRID
         IF (ALLOCATED(MEFMLOC_SUB)) THEN
            IF (MEFMLOC_SUB(ISUB) /= '      ') LOCAL_MEFMLOC = MEFMLOC_SUB(ISUB)
         ENDIF
         IF (ALLOCATED(MEFMGRID_SUB)) LOCAL_MEFMGRID = MEFMGRID_SUB(ISUB)

         XREF(1) = ZERO
         XREF(2) = ZERO
         XREF(3) = ZERO
         IF (LOCAL_MEFMLOC == 'CG    ') THEN
            XREF(1) = MODEL_XCG
            XREF(2) = MODEL_YCG
            XREF(3) = MODEL_ZCG
         ELSE
            IF (LOCAL_MEFMLOC == 'GRDPNT') LOCAL_MEFMGRID = GRDPNT
            IF (LOCAL_MEFMGRID > 0) THEN
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, SIZE(GRID_ID), GRID_ID, LOCAL_MEFMGRID, XREF_GRID_ROW )
               IF (XREF_GRID_ROW > 0) THEN
                  XREF(1) = RGRID(XREF_GRID_ROW,1)
                  XREF(2) = RGRID(XREF_GRID_ROW,2)
                  XREF(3) = RGRID(XREF_GRID_ROW,3)
               ENDIF
            ENDIF
         ENDIF

         XB(1) = MODEL_XCG - XREF(1)
         XB(2) = MODEL_YCG - XREF(2)
         XB(3) = MODEL_ZCG - XREF(3)
         ICG = ZERO
         IREF = ZERO
         LOCAL_MEFM_RB_MASS = ZERO
         DO I=1,3
            DO J=1,3
               ICG(I,J) = MCG(I+3,J+3)
            ENDDO
         ENDDO
         R2 = XB(1)*XB(1) + XB(2)*XB(2) + XB(3)*XB(3)
         DO I=1,3
            DO J=1,3
               IREF(I,J) = ICG(I,J) - MODEL_MASS*XB(I)*XB(J)
            ENDDO
            IREF(I,I) = IREF(I,I) + MODEL_MASS*R2
         ENDDO
         DO I=1,3
            LOCAL_MEFM_RB_MASS(I,I) = MODEL_MASS
         ENDDO
         LOCAL_MEFM_RB_MASS(1,5) =  MODEL_MASS*XB(3)
         LOCAL_MEFM_RB_MASS(1,6) = -MODEL_MASS*XB(2)
         LOCAL_MEFM_RB_MASS(2,4) = -MODEL_MASS*XB(3)
         LOCAL_MEFM_RB_MASS(2,6) =  MODEL_MASS*XB(1)
         LOCAL_MEFM_RB_MASS(3,4) =  MODEL_MASS*XB(2)
         LOCAL_MEFM_RB_MASS(3,5) = -MODEL_MASS*XB(1)
         LOCAL_MEFM_RB_MASS(4,1) =  LOCAL_MEFM_RB_MASS(1,4)
         LOCAL_MEFM_RB_MASS(5,1) =  LOCAL_MEFM_RB_MASS(1,5)
         LOCAL_MEFM_RB_MASS(6,1) =  LOCAL_MEFM_RB_MASS(1,6)
         LOCAL_MEFM_RB_MASS(4,2) =  LOCAL_MEFM_RB_MASS(2,4)
         LOCAL_MEFM_RB_MASS(5,2) =  LOCAL_MEFM_RB_MASS(2,5)
         LOCAL_MEFM_RB_MASS(6,2) =  LOCAL_MEFM_RB_MASS(2,6)
         LOCAL_MEFM_RB_MASS(4,3) =  LOCAL_MEFM_RB_MASS(3,4)
         LOCAL_MEFM_RB_MASS(5,3) =  LOCAL_MEFM_RB_MASS(3,5)
         LOCAL_MEFM_RB_MASS(6,3) =  LOCAL_MEFM_RB_MASS(3,6)
         DO I=1,3
            DO J=1,3
               LOCAL_MEFM_RB_MASS(I+3,J+3) = IREF(I,J)
            ENDDO
         ENDDO

         IF (IHDR == 'Y') THEN
            WRITE(F06,900)
            IF (NSUB > 1) WRITE(F06,9101) SCNUM(ISUB)
            WRITE(F06,909) TITLE(ISUB)
            WRITE(F06,909) STITLE(ISUB)
            WRITE(F06,909) LABEL(ISUB)
            WRITE(F06,*)
         ENDIF

         WRITE(F06,9102) MEFMCORD

         IF      (LOCAL_MEFMLOC == 'GRDPNT') THEN
            IF (LOCAL_MEFMGRID == 0) THEN
               WRITE(F06,9103)
            ELSE
               WRITE(F06,9104) GRDPNT
            ENDIF
         ELSE IF (LOCAL_MEFMLOC == 'CG    ') THEN
            WRITE(F06,9105)
         ELSE IF (LOCAL_MEFMLOC == 'GRID  ') THEN
            WRITE(F06,9106) LOCAL_MEFMGRID
         ENDIF

         DO J=1,6
            MEFM_TOTALS(J) = ZERO
         ENDDO

         IF (WRITE_MEFFM) THEN
            IF (IS_LOW_PRECISION) THEN
               WRITE(F06,9107)
            ELSE
               WRITE(F06,9108)
            ENDIF
         ENDIF

         MODE_NUM_OUT = 0
         DO I=1,NVEC
            IF (HAS_SUBCASE_MAP) THEN
               IF (I > SIZE(MODE_SUBCASE)) EXIT
               IF (MODE_SUBCASE(I) /= ISUB) CYCLE
            ELSE
               IF (ISUB > 1) CYCLE
            ENDIF
            MODE_NUM_OUT = MODE_NUM_OUT + 1
            CYCLES = DSQRT(DABS(EIGEN_VAL(I)))/(TWO*PI)

            IF (WRITE_MEFFM) THEN
               IF (IS_LOW_PRECISION) THEN
                  WRITE(F06,9110) MODE_NUM_OUT, CYCLES, (MEFFMASS(I,J)/WTMASS,J=1,6)
               ELSE
                  WRITE(F06,9111) MODE_NUM_OUT, CYCLES, (MEFFMASS(I,J)/WTMASS,J=1,6)
               ENDIF
            ENDIF

            DO J=1,6
               MEFM_TOTALS(J) = MEFM_TOTALS(J) + MEFFMASS(I,J)/WTMASS
            ENDDO
         ENDDO

         IF (WRITE_MEFFW) THEN
            WRITE(F06,*)
            WRITE(F06,9120) MEFMCORD

            IF      (LOCAL_MEFMLOC == 'GRDPNT') THEN
               IF (LOCAL_MEFMGRID == 0) THEN
                  WRITE(F06,9103)
               ELSE
                  WRITE(F06,9104) GRDPNT
               ENDIF
            ELSE IF (LOCAL_MEFMLOC == 'CG    ') THEN
               WRITE(F06,9105)
            ELSE IF (LOCAL_MEFMLOC == 'GRID  ') THEN
               WRITE(F06,9106) LOCAL_MEFMGRID
            ENDIF

            IF (IS_LOW_PRECISION) THEN
               WRITE(F06,9107)
            ELSE
               WRITE(F06,9108)
            ENDIF

            MODE_NUM_OUT = 0
            DO I=1,NVEC
               IF (HAS_SUBCASE_MAP) THEN
                  IF (I > SIZE(MODE_SUBCASE)) EXIT
                  IF (MODE_SUBCASE(I) /= ISUB) CYCLE
               ELSE
                  IF (ISUB > 1) CYCLE
               ENDIF
               MODE_NUM_OUT = MODE_NUM_OUT + 1
               CYCLES = DSQRT(DABS(EIGEN_VAL(I)))/(TWO*PI)
               IF (IS_LOW_PRECISION) THEN
                  WRITE(F06,9110) MODE_NUM_OUT, CYCLES, (MEFFMASS(I,J),J=1,6)
               ELSE
                  WRITE(F06,9111) MODE_NUM_OUT, CYCLES, (MEFFMASS(I,J),J=1,6)
               ENDIF
            ENDDO

            IF (WRITE_SUMMARY) THEN
               IF (IS_LOW_PRECISION) THEN
               WRITE(F06,9112) (WTMASS*MEFM_TOTALS(J),J=1,6)
                  WRITE(F06,9122) (WTMASS*LOCAL_MEFM_RB_MASS(I,I),I=1,6)
               ELSE
                  WRITE(F06,9113) (WTMASS*MEFM_TOTALS(J),J=1,6)
                  WRITE(F06,9123) (WTMASS*LOCAL_MEFM_RB_MASS(I,I),I=1,6)
               ENDIF
            ENDIF
         ENDIF

         IF (WRITE_SUMMARY) THEN
            IF (IS_LOW_PRECISION) THEN
               WRITE(F06,9112) (MEFM_TOTALS(J),J=1,6)
            ELSE
               WRITE(F06,9113) (MEFM_TOTALS(J),J=1,6)
            ENDIF
            IF (IS_LOW_PRECISION) THEN
               WRITE(F06,9116) (LOCAL_MEFM_RB_MASS(I,I),I=1,6)
            ELSE
               WRITE(F06,9117) (LOCAL_MEFM_RB_MASS(I,I),I=1,6)
            ENDIF
         ENDIF

         IF (DABS(LOCAL_MEFM_RB_MASS(1,1)) > EPS1) THEN
            MODES_PCT(1) = ONE_HUNDRED*MEFM_TOTALS(1)/LOCAL_MEFM_RB_MASS(1,1)
            WRITE(CHAR_PCT(1),1001) MODES_PCT(1), '%'
         ELSE
            CHAR_PCT(1)(1:) = ' '
            IF (MEFM_TOTALS(1) > EPS1) THEN
               WRITE(ERR,2001)
               IF (SUPINFO == 'N') WRITE(F06,2001)
            ENDIF
         ENDIF

         IF (DABS(LOCAL_MEFM_RB_MASS(2,2)) > EPS1) THEN
            MODES_PCT(2) = ONE_HUNDRED*MEFM_TOTALS(2)/LOCAL_MEFM_RB_MASS(2,2)
            WRITE(CHAR_PCT(2),1001) MODES_PCT(2), '%'
         ELSE
            CHAR_PCT(2)(1:) = ' '
            IF (MEFM_TOTALS(2) > EPS1) THEN
               WRITE(ERR,2002)
               IF (SUPINFO == 'N') WRITE(F06,2002)
            ENDIF
         ENDIF

         IF (DABS(LOCAL_MEFM_RB_MASS(3,3)) > EPS1) THEN
            MODES_PCT(3) = ONE_HUNDRED*MEFM_TOTALS(3)/LOCAL_MEFM_RB_MASS(3,3)
            WRITE(CHAR_PCT(3),1001) MODES_PCT(3), '%'
         ELSE
            CHAR_PCT(3)(1:) = ' '
            IF (MEFM_TOTALS(3) > EPS1) THEN
               WRITE(ERR,2003)
               IF (SUPINFO == 'N') WRITE(F06,2003)
            ENDIF
         ENDIF

         IF (DABS(LOCAL_MEFM_RB_MASS(4,4)) > EPS1) THEN
            MODES_PCT(4) = ONE_HUNDRED*MEFM_TOTALS(4)/LOCAL_MEFM_RB_MASS(4,4)
            WRITE(CHAR_PCT(4),1001) MODES_PCT(4), '%'
         ELSE
            CHAR_PCT(4)(1:) = ' '
            IF (MEFM_TOTALS(4) > EPS1) THEN
               WRITE(ERR,2004)
               IF (SUPINFO == 'N') WRITE(F06,2004)
            ENDIF
         ENDIF

         IF (DABS(LOCAL_MEFM_RB_MASS(5,5)) > EPS1) THEN
            MODES_PCT(5) = ONE_HUNDRED*MEFM_TOTALS(5)/LOCAL_MEFM_RB_MASS(5,5)
            WRITE(CHAR_PCT(5),1001) MODES_PCT(5), '%'
         ELSE
            CHAR_PCT(5)(1:) = ' '
            IF (MEFM_TOTALS(5) > EPS1) THEN
               WRITE(ERR,2005)
               IF (SUPINFO == 'N') WRITE(F06,2005)
            ENDIF
         ENDIF

         IF (DABS(LOCAL_MEFM_RB_MASS(6,6)) > EPS1) THEN
            MODES_PCT(6) = ONE_HUNDRED*MEFM_TOTALS(6)/LOCAL_MEFM_RB_MASS(6,6)
            WRITE(CHAR_PCT(6),1001) MODES_PCT(6), '%'
         ELSE
            CHAR_PCT(6)(1:) = ' '
            IF (MEFM_TOTALS(6) > EPS1) THEN
               WRITE(ERR,2006)
               IF (SUPINFO == 'N') WRITE(F06,2006)
            ENDIF
         ENDIF

         IF (WRITE_FRACSUM) THEN
            WRITE(F06,9118) (CHAR_PCT(I),I=1,6)
         ENDIF

      ENDDO



      RETURN

! **********************************************************************************************************************************
  900 FORMAT('--------------------------------------------------------------------------------------------------------------------'&
            ,'----------------')

  909 FORMAT(1X,A)

 1001 FORMAT(F13.2,A)

 2001 FORMAT(' *INFORMATION: CANNOT CALCULATE T1 MODAL EFFECTIVE MASS PERCENT OF MODEL MASS SO IT IS LEFT BLANK')

 2002 FORMAT(' *INFORMATION: CANNOT CALCULATE T2 MODAL EFFECTIVE MASS PERCENT OF MODEL MASS SO IT IS LEFT BLANK')

 2003 FORMAT(' *INFORMATION: CANNOT CALCULATE T3 MODAL EFFECTIVE MASS PERCENT OF MODEL MASS SO IT IS LEFT BLANK')

 2004 FORMAT(' *INFORMATION: CANNOT CALCULATE R1 MODAL EFFECTIVE MASS PERCENT OF MODEL IXX  SO IT IS LEFT BLANK')

  2005 FORMAT(' *INFORMATION: CANNOT CALCULATE R2 MODAL EFFECTIVE MASS PERCENT OF MODEL IYY  SO IT IS LEFT BLANK')

  2006 FORMAT(' *INFORMATION: CANNOT CALCULATE R3 MODAL EFFECTIVE MASS PERCENT OF MODEL IZZ  SO IT IS LEFT BLANK')

 9101 FORMAT(1X,'OUTPUT FOR SUBCASE ',I8)

 9102 FORMAT(14X,'                    E F F E C T I V E   M O D A L   M A S S E S   O R   W E I G H T S',/,                        &
             14X,'                                     (in coordinate system ',I8,')',/,                                           &
             14X,'                       Units are same as units for mass input in the Bulk Data Deck')

 9103 FORMAT(14X,'                          Reference point is the basic coordinate system origin',/)

 9104 FORMAT(14X,'                            Reference point is the PARAM GRDPNT grid: ',I8,/)

 9105 FORMAT(14X,'                              Reference point is the model center of gravity',/)

 9106 FORMAT(14X,'                                    Reference point is grid ',I8,/)

 9107 FORMAT(13X,'MODE     CYCLES          T1            T2            T3            R1            R2            R3',/,            &
             13X,' NUM')

 9108 FORMAT(13X,'MODE       CYCLES          T1            T2            T3            R1            R2            R3',/,          &
             13X,' NUM')

 9110 FORMAT(9X,I8,7(1ES14.6))

 9111 FORMAT(9X,I8,7(1ES14.2))

 9112 FORMAT(32X,' ------------  ------------  ------------  ------------  ------------  ------------',/,                          &
             17X,'Sum all modes:',6(1ES14.6))

 9113 FORMAT(32X,'     --------      --------      --------      --------      --------      --------',/,                          &
             17X,'Sum all modes:',6(1ES14.2))

 9116 FORMAT(14X,'Total model mass:',6(1ES14.6))

 9117 FORMAT(14X,'Total model mass:',6(1ES14.2))

 9118 FORMAT(8X,'Modes % of total mass*:',6(A14),//,' *If all modes are calculated the % of total mass should be 100% of the '     &
                ,'free mass (i.e. not counting mass at constrained DOF''s).',/,                                                     &
                '  Percentages are only printed for components that have finite model mass.',/,                                     &
                '                                                               -----')

 9120 FORMAT(14X,'                   E F F E C T I V E   M O D A L   W E I G H T S',/,                                              &
             14X,'                                     (in coordinate system ',I8,')',/,                                           &
             14X,'                     Units are same as the weighted mass matrix after PARAM WTMASS is applied')

 9122 FORMAT(14X,'Total model weight:',6(1ES14.6))

 9123 FORMAT(14X,'Total model weight:',6(1ES14.2))

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_MEFFMASS
