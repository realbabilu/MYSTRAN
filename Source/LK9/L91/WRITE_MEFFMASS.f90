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
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NVEC
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, ONE_HUNDRED, PI
      USE PARAMS, ONLY                :  PRTF06, PRTOP2
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL, MEFFMASS, MPFACTOR_N6
      USE MODEL_STUF, ONLY            :  MEFM_RB_MASS, LABEL, STITLE, TITLE
      USE PARAMS, ONLY                :  EPSIL, GRDPNT, MEFMCORD, MEFMGRID, MEFMLOC, SUPINFO, WTMASS

      USE WRITE_MEFFMASS_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_MEFFMASS'
      CHARACTER(14*BYTE)              :: CHAR_PCT(6)       ! Character representation of MEFFMASS sum percents of total model mass
      CHARACTER(1*BYTE)               :: IHDR   = 'Y'      ! Indicator of whether to write an output header
      CHARACTER(18*BYTE)              :: DIR_NAME(6)
      CHARACTER(8*BYTE)               :: DIR_CODE(6)

      INTEGER(LONG)                   :: I,J               ! DO loop indices


      REAL(DOUBLE)                    :: CYCLES            ! Circular frequency of a mode
      REAL(DOUBLE)                    :: PERIOD            ! Period of a mode
      REAL(DOUBLE)                    :: EPS1              ! Small number to compare against zero
      REAL(DOUBLE)                    :: MEFM_TOTALS(6)    ! Totals for the 6 modal effective masses over all modes
      REAL(DOUBLE)                    :: MODES_PCT(6)      ! Modal mass as % of total mass
      REAL(DOUBLE)                    :: MASS_RATIO        ! Effective mass ratio to total model mass in a direction
      REAL(DOUBLE)                    :: CUM_MASS_RATIO(6) ! Cumulative effective mass ratio
      !LOGICAL                        :: WRITE_F06  ! flag
      !LOGICAL                        :: WRITE_OP2  ! flag
      LOGICAL                         :: IS_LOW_PRECISION  ! Print MPFACTOR, MEFFMASS values with 2 decimal places of accuracy rather than 6




! **********************************************************************************************************************************
      IS_LOW_PRECISION = (DEBUG(174) == 0)
      !--------------------------------------------------

      EPS1 = EPSIL(1)
      DIR_NAME(1) = 'X TRANSLATION'
      DIR_NAME(2) = 'Y TRANSLATION'
      DIR_NAME(3) = 'Z TRANSLATION'
      DIR_NAME(4) = 'X ROTATION'
      DIR_NAME(5) = 'Y ROTATION'
      DIR_NAME(6) = 'Z ROTATION'
      DIR_CODE(1) = 'T1'
      DIR_CODE(2) = 'T2'
      DIR_CODE(3) = 'T3'
      DIR_CODE(4) = 'R1'
      DIR_CODE(5) = 'R2'
      DIR_CODE(6) = 'R3'

      ! Write output headers.
      IF (IHDR == 'Y') THEN
         WRITE(F06,900)
         ! There is always a TITLE(1), etc (even if they are blank)
         WRITE(F06,909) TITLE(1)
         WRITE(F06,909) STITLE(1)
         WRITE(F06,909) LABEL(1)
         WRITE(F06,*)
      ENDIF

      ! Write modal effective masses
      WRITE(F06,9102) MEFMCORD

      IF      (MEFMLOC == 'GRDPNT') THEN
         IF (MEFMGRID == 0) THEN
            WRITE(F06,9103)
         ELSE
            WRITE(F06,9104) GRDPNT
         ENDIF
      ELSE IF (MEFMLOC == 'CG    ') THEN
         WRITE(F06,9105)
      ELSE IF (MEFMLOC == 'GRID  ') THEN
         WRITE(F06,9106) MEFMGRID
      ENDIF

      IF (IS_LOW_PRECISION) THEN
         WRITE(F06,9107)
      ELSE
         WRITE(F06,9108)
      ENDIF

      DO J=1,6
         MEFM_TOTALS(J) = ZERO
      ENDDO

      DO I=1,NVEC
         CYCLES = DSQRT(DABS(EIGEN_VAL(I)))/(TWO*PI)
         IF (CYCLES > ZERO) THEN
            PERIOD = ONE/CYCLES
         ELSE
            PERIOD = ZERO
         ENDIF

         IF (IS_LOW_PRECISION) THEN ! 6 digits
            WRITE(F06,9110) I, CYCLES, PERIOD, (MEFFMASS(I,J)/WTMASS,J=1,6)
         ELSE ! low precision (2 digits)
            WRITE(F06,9111) I, CYCLES, PERIOD, (MEFFMASS(I,J)/WTMASS,J=1,6)
         ENDIF

         DO J=1,6
            MEFM_TOTALS(J) = MEFM_TOTALS(J) + MEFFMASS(I,J)/WTMASS
         ENDDO

      ENDDO

      IF (IS_LOW_PRECISION) THEN
         WRITE(F06,9112) (MEFM_TOTALS(J),J=1,6)
      ELSE
         WRITE(F06,9113) (MEFM_TOTALS(J),J=1,6)
      ENDIF
                                                           ! MEFM_RB_MASS is in the same units as in the DAT file
      IF (IS_LOW_PRECISION) THEN
         WRITE(F06,9116) (MEFM_RB_MASS(I,I),I=1,6)
      ELSE
         WRITE(F06,9117) (MEFM_RB_MASS(I,I),I=1,6)
      ENDIF

      ! For each of the 6 modal masses, calc % of total mass.
      ! A character variable is used to store the % so that blank percentages can
      ! be printed if zero modal mass exists for a component (T1 - R3) or
      ! if a denominator in the % expression is zero
      IF (DABS(MEFM_RB_MASS(1,1)) > EPS1) THEN
         MODES_PCT(1) = ONE_HUNDRED*MEFM_TOTALS(1)/MEFM_RB_MASS(1,1)
         WRITE(CHAR_PCT(1),1001) MODES_PCT(1), '%'
      ELSE
         ! Denominator is zero so leave % blank
         CHAR_PCT(1)(1:) = ' '
         IF (MEFM_TOTALS(1) > EPS1) THEN
            ! Error: denominator is zero but numerator is not, so write message
            WRITE(ERR,2001)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,2001)
            ENDIF
         ENDIF
      ENDIF

      IF (DABS(MEFM_RB_MASS(2,2)) > EPS1) THEN
         MODES_PCT(2) = ONE_HUNDRED*MEFM_TOTALS(2)/MEFM_RB_MASS(2,2)
         WRITE(CHAR_PCT(2),1001) MODES_PCT(2), '%'
      ELSE
         ! Denominator is zero so leave % blank
         CHAR_PCT(2)(1:) = ' '
         IF (MEFM_TOTALS(2) > EPS1) THEN
            ! Error: denominator is zero but numerator is not, so write message
            WRITE(ERR,2002)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,2002)
            ENDIF
         ENDIF
      ENDIF

      IF (DABS(MEFM_RB_MASS(3,3)) > EPS1) THEN
         MODES_PCT(3) = ONE_HUNDRED*MEFM_TOTALS(3)/MEFM_RB_MASS(3,3)
         WRITE(CHAR_PCT(3),1001) MODES_PCT(3), '%'
      ELSE
         ! Denominator is zero so leave % blank
         CHAR_PCT(3)(1:) = ' '
         IF (MEFM_TOTALS(3) > EPS1) THEN
            ! Error: denominator is zero but numerator is not, so write message
            WRITE(ERR,2003)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,2003)
            ENDIF
         ENDIF
      ENDIF

      IF (DABS(MEFM_RB_MASS(4,4)) > EPS1) THEN
         MODES_PCT(4) = ONE_HUNDRED*MEFM_TOTALS(4)/MEFM_RB_MASS(4,4)
         WRITE(CHAR_PCT(4),1001) MODES_PCT(4), '%'
      ELSE
         ! Denominator is zero so leave % blank
         CHAR_PCT(4)(1:) = ' '
         IF (MEFM_TOTALS(4) > EPS1) THEN
            ! Error: denominator is zero but numerator is not, so write message
            WRITE(ERR,2004)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,2004)
            ENDIF
         ENDIF
      ENDIF

      IF (DABS(MEFM_RB_MASS(5,5)) > EPS1) THEN
         MODES_PCT(5) = ONE_HUNDRED*MEFM_TOTALS(5)/MEFM_RB_MASS(5,5)
         WRITE(CHAR_PCT(5),1001) MODES_PCT(5), '%'
      ELSE
         ! Denominator is zero so leave % blank
         CHAR_PCT(5)(1:) = ' '
         IF (MEFM_TOTALS(5) > EPS1) THEN
            ! Error: denominator is zero but numerator is not, so write message
            WRITE(ERR,2005)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,2005)
            ENDIF
         ENDIF
      ENDIF

      IF (DABS(MEFM_RB_MASS(6,6)) > EPS1) THEN
         MODES_PCT(6) = ONE_HUNDRED*MEFM_TOTALS(6)/MEFM_RB_MASS(6,6)
         WRITE(CHAR_PCT(6),1001) MODES_PCT(6), '%'
      ELSE
         ! Denominator is zero so leave % blank
         CHAR_PCT(6)(1:) = ' '
         IF (MEFM_TOTALS(6) > EPS1) THEN
            ! Error: denominator is zero but numerator is not, so write message
            WRITE(ERR,2006)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,2006)
            ENDIF
         ENDIF
      ENDIF

      WRITE(F06,9118) (CHAR_PCT(I),I=1,6)
      WRITE(F06,*)
      WRITE(F06,9120)
      DO J=1,6
         CUM_MASS_RATIO(J) = ZERO
         WRITE(F06,9121) DIR_NAME(J), DIR_CODE(J)
         WRITE(F06,9122)
         DO I=1,NVEC
            CYCLES = DSQRT(DABS(EIGEN_VAL(I)))/(TWO*PI)
            IF (CYCLES > ZERO) THEN
               PERIOD = ONE/CYCLES
            ELSE
               PERIOD = ZERO
            ENDIF
            MASS_RATIO = ZERO
            IF (DABS(MEFM_RB_MASS(J,J)) > EPS1) THEN
               MASS_RATIO = (MEFFMASS(I,J)/WTMASS)/MEFM_RB_MASS(J,J)
            ENDIF
            CUM_MASS_RATIO(J) = CUM_MASS_RATIO(J) + MASS_RATIO
            IF (IS_LOW_PRECISION) THEN
               WRITE(F06,9123) I, CYCLES, PERIOD, MPFACTOR_N6(I,J), MEFFMASS(I,J)/WTMASS, MASS_RATIO, CUM_MASS_RATIO(J)
            ELSE
               WRITE(F06,9124) I, CYCLES, PERIOD, MPFACTOR_N6(I,J), MEFFMASS(I,J)/WTMASS, MASS_RATIO, CUM_MASS_RATIO(J)
            ENDIF
         ENDDO
         WRITE(F06,*)
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

 9102 FORMAT(14X,'                    E F F E C T I V E   M O D A L   M A S S E S   O R   W E I G H T S',/,                        &
             14X,'                                     (in coordinate system ',I8,')',/,                                           &
             14X,'                       Units are same as units for mass input in the Bulk Data Deck')

 9103 FORMAT(14X,'                          Reference point is the basic coordinate system origin',/)

 9104 FORMAT(14X,'                            Reference point is the PARAM GRDPNT grid: ',I8,/)

 9105 FORMAT(14X,'                              Reference point is the model center of gravity',/)

 9106 FORMAT(14X,'                                    Reference point is grid ',I8,/)

 9107 FORMAT(13X,'MODE    FREQ(Hz)     PERIOD(s)        T1            T2            T3            R1            R2            R3',/,            &
             13X,' NUM')

 9108 FORMAT(13X,'MODE      FREQ(Hz)     PERIOD(s)        T1            T2            T3            R1            R2            R3',/,          &
             13X,' NUM')

 9110 FORMAT(9X,I8,8(1ES14.6))

 9111 FORMAT(9X,I8,8(1ES14.2))

 9112 FORMAT(45X,' ------------  ------------  ------------  ------------  ------------  ------------',/,&
             31X,'Sum all modes:',6(1ES14.6))

 9113 FORMAT(45X,'     --------      --------      --------      --------      --------      --------',/,&
             31X,'Sum all modes:',6(1ES14.2))

 9116 FORMAT(28X,'Total model mass:',6(1ES14.6))

 9117 FORMAT(28X,'Total model mass:',6(1ES14.2))

  9118 FORMAT(22X,'Modes % of total mass*:',6(A14),//,' *If all modes are calculated the % of total mass should be 100% of the '     &
               ,'free mass (i.e. not counting mass at constrained DOF''s).',/,                                                     &
               '  Percentages are only printed for components that have finite model mass.',/,                                     &
               '                                                               -----')

  9120 FORMAT(14X,'D E T A I L E D   M O D A L   M A S S   P A R T I C I P A T I O N   T A B L E S',/)

  9121 FORMAT(14X,A,'  (',A,')')

  9122 FORMAT(13X,'MODE    FREQ(Hz)     PERIOD(s)    PARTIC.FACTOR    EFFECTIVE MASS     MASS FRACTION      CUMULATIVE',/,        &
             13X,' NUM')

  9123 FORMAT(9X,I8,6(1ES16.6))

  9124 FORMAT(9X,I8,6(1ES16.8))

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_MEFFMASS

