!--- cbeam_add begin ---!
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

      SUBROUTINE BEAM ( OPT, L, AREA, I1, I2, JTOR, SCOEFF, K1, K2, I12, E, G, ALPHA, TREF )

! Calculates, for 1-D general beam element using a DSB/Timoshenko bending formulation:
!
!  1) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  2) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  3) KE        = element linear stiffness matrix      , if OPT(6) = 'N' (i.e. always calc KE linear unless OPT(6) = 'Y')
!  4) KED       = element differen stiff matrix        , if OPT(6) = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  NSUB, NTSUB, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX, TEN, TWELVE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  EPSIL, ART_KED, ART_ROT_KED, ART_TRAN_KED
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  DOFPIN, DT, EID, KE, KED, PEL, PPE, PRESS, PTE, SE1, SE2, STE1, STE2

      USE BEAM_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BEAM'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)

      INTEGER(LONG)                   :: I,J
      INTEGER(LONG)                   :: NUM_PFLAG_DOFS

      REAL(DOUBLE), INTENT(IN)        :: ALPHA
      REAL(DOUBLE), INTENT(IN)        :: AREA
      REAL(DOUBLE), INTENT(IN)        :: E
      REAL(DOUBLE), INTENT(IN)        :: G
      REAL(DOUBLE), INTENT(IN)        :: I1
      REAL(DOUBLE), INTENT(IN)        :: I12
      REAL(DOUBLE), INTENT(IN)        :: I2
      REAL(DOUBLE), INTENT(IN)        :: JTOR
      REAL(DOUBLE), INTENT(IN)        :: K1
      REAL(DOUBLE), INTENT(IN)        :: K2
      REAL(DOUBLE), INTENT(IN)        :: L
      REAL(DOUBLE), INTENT(IN)        :: SCOEFF
      REAL(DOUBLE), INTENT(IN)        :: TREF
      REAL(DOUBLE)                    :: ABAR(6,5)
      REAL(DOUBLE)                    :: B1(3,6)
      REAL(DOUBLE)                    :: B2(3,6)
      REAL(DOUBLE)                    :: BT1(3,5)
      REAL(DOUBLE)                    :: BT2(3,5)
      REAL(DOUBLE)                    :: BTA(6,5)
      REAL(DOUBLE)                    :: BTB(6,5)
      REAL(DOUBLE)                    :: C01
      REAL(DOUBLE)                    :: DELTA1
      REAL(DOUBLE)                    :: DELTA2
      REAL(DOUBLE)                    :: DELTA12
      REAL(DOUBLE)                    :: DEN
      REAL(DOUBLE)                    :: DUM1(3,NTSUB)
      REAL(DOUBLE)                    :: DUM2(3,NTSUB)
      REAL(DOUBLE)                    :: EPS1
      REAL(DOUBLE)                    :: FAC1
      REAL(DOUBLE)                    :: FAC2
      REAL(DOUBLE)                    :: FX
      REAL(DOUBLE)                    :: KAA(6,6)
      REAL(DOUBLE)                    :: KAB(6,6)
      REAL(DOUBLE)                    :: KBA(6,6)
      REAL(DOUBLE)                    :: M1A
      REAL(DOUBLE)                    :: M1B
      REAL(DOUBLE)                    :: M2A
      REAL(DOUBLE)                    :: M2B
      REAL(DOUBLE)                    :: N1
      REAL(DOUBLE)                    :: N2
      REAL(DOUBLE)                    :: N3
      REAL(DOUBLE)                    :: N4
      REAL(DOUBLE)                    :: P1
      REAL(DOUBLE)                    :: P2
      REAL(DOUBLE)                    :: PC
      REAL(DOUBLE)                    :: PHI1
      REAL(DOUBLE)                    :: PHI2
      REAL(DOUBLE)                    :: PTA(6,NTSUB)
      REAL(DOUBLE)                    :: PTB(6,NTSUB)
      REAL(DOUBLE)                    :: QT
      REAL(DOUBLE)                    :: RG
      REAL(DOUBLE)                    :: S11(3,6)
      REAL(DOUBLE)                    :: S12(3,6)
      REAL(DOUBLE)                    :: S21(3,6)
      REAL(DOUBLE)                    :: S22(3,6)
      REAL(DOUBLE)                    :: TBAR
      REAL(DOUBLE)                    :: TPRIME(5,NTSUB)
      REAL(DOUBLE)                    :: V1
      REAL(DOUBLE)                    :: V2
      REAL(DOUBLE)                    :: WGT
      REAL(DOUBLE)                    :: X1L
      REAL(DOUBLE)                    :: X2L
      REAL(DOUBLE)                    :: XI
      REAL(DOUBLE)                    :: XI_GAUSS(3)
      REAL(DOUBLE)                    :: XI_SCALE
      REAL(DOUBLE)                    :: XI_WGT(3)

! MYSTRAN's 1D force/recovery conventions tie:
!   plane 1 -> DOFs (UY,RZ) with I1 and K2
!   plane 2 -> DOFs (UZ,RY) with I2 and K1
! The DSB shear-influence factors are therefore derived using those same pairings.

      INTRINSIC DABS

! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

      DO I=1,12
         DO J=1,12
            KE(I,J)  = ZERO
            KED(I,J) = ZERO
         ENDDO
      ENDDO

      FAC1 = E*I1/(L*L*L)
      FAC2 = E*I2/(L*L*L)
      RG   = G*JTOR/L

      IF (DEBUG(203) > 0) CALL DEBUG_BEAM ( 1 )

      PHI1 = ZERO
      PHI2 = ZERO
      IF ((DABS(K1) <= EPS1) .AND. (DABS(K2) <= EPS1)) THEN
         PHI1 = ZERO
         PHI2 = ZERO
      ELSE
         IF (DABS(K2*G*AREA*L*L) > EPS1) THEN
            PHI1 = TWELVE*E*I1/(K2*G*AREA*L*L)
         ENDIF
         IF (DABS(K1*G*AREA*L*L) > EPS1) THEN
            PHI2 = TWELVE*E*I2/(K1*G*AREA*L*L)
         ENDIF
      ENDIF

      FAC1 = FAC1/(ONE + PHI1)
      FAC2 = FAC2/(ONE + PHI2)

      DEN     = I1*I2 - I12*I12
      DELTA1  = ZERO
      DELTA2  = ZERO
      DELTA12 = ZERO
      IF (DABS(DEN) > EPS1) THEN
         DELTA1  = I2/DEN
         DELTA2  = I1/DEN
         DELTA12 = I12/DEN
      ENDIF

      IF (DEBUG(203) > 0) CALL DEBUG_BEAM ( 2 )

! **********************************************************************************************************************************
! Stiffness matrix

      KE( 1, 1) = AREA*E/L
      KE( 1, 7) =-KE(1,1)
      KE( 7, 7) = KE(1,1)

      KE( 4, 4) = RG
      KE( 4,10) =-RG
      KE(10,10) = RG

! Plane 1 DSB bending block on DOFs (UY, RZ, UY, RZ)
      KE( 2, 2) =  TWELVE*FAC1
      KE( 2, 6) =  SIX*L*FAC1
      KE( 2, 8) = -TWELVE*FAC1
      KE( 2,12) =  SIX*L*FAC1

      KE( 6, 6) = (FOUR + PHI1)*E*I1/(L*(ONE + PHI1))
      KE( 6, 8) = -SIX*L*FAC1
      KE( 6,12) = (TWO - PHI1)*E*I1/(L*(ONE + PHI1))

      KE( 8, 8) =  TWELVE*FAC1
      KE( 8,12) = -SIX*L*FAC1

      KE(12,12) = (FOUR + PHI1)*E*I1/(L*(ONE + PHI1))

! Plane 2 DSB bending block on DOFs (UZ, RY, UZ, RY)
      KE( 3, 3) =  TWELVE*FAC2
      KE( 3, 5) = -SIX*L*FAC2
      KE( 3, 9) = -TWELVE*FAC2
      KE( 3,11) = -SIX*L*FAC2

      KE( 5, 5) = (FOUR + PHI2)*E*I2/(L*(ONE + PHI2))
      KE( 5, 9) =  SIX*L*FAC2
      KE( 5,11) = (TWO - PHI2)*E*I2/(L*(ONE + PHI2))

      KE( 9, 9) =  TWELVE*FAC2
      KE( 9,11) =  SIX*L*FAC2

      KE(11,11) = (FOUR + PHI2)*E*I2/(L*(ONE + PHI2))

      DO I=2,12
         DO J=1,I-1
            KE(I,J) = KE(J,I)
         ENDDO
      ENDDO

! **********************************************************************************************************************************
! Pin flags

      NUM_PFLAG_DOFS = 0
      DO I=1,12
         IF (DOFPIN(I) > 0) THEN
            NUM_PFLAG_DOFS = NUM_PFLAG_DOFS + 1
         ENDIF
      ENDDO
      IF (NUM_PFLAG_DOFS /= 0) THEN
         CALL PINFLG ( NUM_PFLAG_DOFS )
      ENDIF

      DO I=1,6
         DO J=1,6
            KAA(I,J) = KE(I,J)
            KAB(I,J) = KE(I,J+6)
            KBA(I,J) = KE(I+6,J)
         ENDDO
      ENDDO

! **********************************************************************************************************************************
! Temperatures

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         IF (NTSUB > 0) THEN
            DO J=1,NTSUB
               TBAR        = (DT(1,J) + DT(2,J))/TWO
               TPRIME(1,J) = TBAR - TREF
               TPRIME(2,J) = DT(3,J)
               TPRIME(3,J) = DT(4,J)
               TPRIME(4,J) = DT(5,J)
               TPRIME(5,J) = DT(6,J)
            ENDDO
         ENDIF
      ENDIF

! **********************************************************************************************************************************
! Thermal loads and thermal stress recovery use the existing 1D beam thermal convention.

      IF (NTSUB > 0) THEN

         DO I=1,6
            DO J=1,5
               ABAR(I,J) = ZERO
            ENDDO
         ENDDO

         ABAR(1,1) =  ONE
         ABAR(2,2) =  DELTA1*I1*L/SIX
         ABAR(2,3) =  DELTA1*I1*L/THREE
         ABAR(2,4) = -DELTA12*I2*L/SIX
         ABAR(2,5) = -DELTA12*I2*L/THREE
         ABAR(3,2) = -DELTA12*I1*L/SIX
         ABAR(3,3) = -DELTA12*I1*L/THREE
         ABAR(3,4) =  DELTA2*I2*L/SIX
         ABAR(3,5) =  DELTA2*I2*L/THREE
         ABAR(5,2) = -DELTA12*I1/TWO
         ABAR(5,3) = -DELTA12*I1/TWO
         ABAR(5,4) =  DELTA2*I2/TWO
         ABAR(5,5) =  DELTA2*I2/TWO
         ABAR(6,2) = -DELTA1*I1/TWO
         ABAR(6,3) = -DELTA1*I1/TWO
         ABAR(6,4) =  DELTA12*I2/TWO
         ABAR(6,5) =  DELTA12*I2/TWO

         DO I=1,6
            DO J=1,5
               ABAR(I,J) = -ALPHA*L*ABAR(I,J)
            ENDDO
         ENDDO

         CALL MATMULT_FFF ( KAA, ABAR, 6, 6, 5, BTA )
         CALL MATMULT_FFF ( KBA, ABAR, 6, 6, 5, BTB )

         CALL MATMULT_FFF ( BTA, TPRIME, 6, 5, NTSUB, PTA )
         CALL MATMULT_FFF ( BTB, TPRIME, 6, 5, NTSUB, PTB )
         DO I=1,6
            DO J=1,NTSUB
               PTE(I,J)   = PTA(I,J)
               PTE(I+6,J) = PTB(I,J)
            ENDDO
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Determine element load vector PPE from local beam-axis PLOAD1 data:
! each component uses [P1,P2,X1,X2], where X values are fractional [0,1].
! X1 = X2 is treated as concentrated-in-element.

      IF (OPT(5) == 'Y') THEN
         XI_GAUSS(1) = -0.774596669241483D0
         XI_GAUSS(2) =  ZERO
         XI_GAUSS(3) =  0.774596669241483D0
         XI_WGT(1)   =  0.555555555555556D0
         XI_WGT(2)   =  0.888888888888889D0
         XI_WGT(3)   =  0.555555555555556D0

         DO J=1,NSUB
            P1  = PRESS(1,J)
            P2  = PRESS(2,J)
            X1L = PRESS(3,J)
            X2L = PRESS(4,J)
            IF (X1L >= ZERO) THEN
               IF (DABS(X2L - X1L) <= EPS1) THEN
                  PC = P1
                  N1 = ONE - THREE*X1L*X1L + TWO*X1L*X1L*X1L
                  N2 = L*(X1L - TWO*X1L*X1L + X1L*X1L*X1L)
                  N3 = THREE*X1L*X1L - TWO*X1L*X1L*X1L
                  N4 = L*(-X1L*X1L + X1L*X1L*X1L)
                  PPE( 2,J) = PPE( 2,J) + PC*N1
                  PPE( 6,J) = PPE( 6,J) + PC*N2
                  PPE( 8,J) = PPE( 8,J) + PC*N3
                  PPE(12,J) = PPE(12,J) + PC*N4
               ELSE
                  XI_SCALE = (X2L - X1L)/TWO
                  DO I=1,3
                     XI  = XI_SCALE*XI_GAUSS(I) + (X2L + X1L)/TWO
                     WGT = XI_WGT(I)
                     QT  = P1 + (P2-P1)*(XI-X1L)/(X2L-X1L)
                     N1 = ONE - THREE*XI*XI + TWO*XI*XI*XI
                     N2 = L*(XI - TWO*XI*XI + XI*XI*XI)
                     N3 = THREE*XI*XI - TWO*XI*XI*XI
                     N4 = L*(-XI*XI + XI*XI*XI)
                     PPE( 2,J) = PPE( 2,J) + QT*L*WGT*XI_SCALE*N1
                     PPE( 6,J) = PPE( 6,J) + QT*L*WGT*XI_SCALE*N2
                     PPE( 8,J) = PPE( 8,J) + QT*L*WGT*XI_SCALE*N3
                     PPE(12,J) = PPE(12,J) + QT*L*WGT*XI_SCALE*N4
                  ENDDO
               ENDIF
            ENDIF

            P1  = PRESS(5,J)
            P2  = PRESS(6,J)
            X1L = PRESS(7,J)
            X2L = PRESS(8,J)
            IF (X1L >= ZERO) THEN
               IF (DABS(X2L - X1L) <= EPS1) THEN
                  PC = P1
                  N1 = ONE - THREE*X1L*X1L + TWO*X1L*X1L*X1L
                  N2 = L*(X1L - TWO*X1L*X1L + X1L*X1L*X1L)
                  N3 = THREE*X1L*X1L - TWO*X1L*X1L*X1L
                  N4 = L*(-X1L*X1L + X1L*X1L*X1L)
                  PPE( 3,J) = PPE( 3,J) + PC*N1
                  PPE( 5,J) = PPE( 5,J) - PC*N2
                  PPE( 9,J) = PPE( 9,J) + PC*N3
                  PPE(11,J) = PPE(11,J) - PC*N4
               ELSE
                  XI_SCALE = (X2L - X1L)/TWO
                  DO I=1,3
                     XI  = XI_SCALE*XI_GAUSS(I) + (X2L + X1L)/TWO
                     WGT = XI_WGT(I)
                     QT  = P1 + (P2-P1)*(XI-X1L)/(X2L-X1L)
                     N1 = ONE - THREE*XI*XI + TWO*XI*XI*XI
                     N2 = L*(XI - TWO*XI*XI + XI*XI*XI)
                     N3 = THREE*XI*XI - TWO*XI*XI*XI
                     N4 = L*(-XI*XI + XI*XI*XI)
                     PPE( 3,J) = PPE( 3,J) + QT*L*WGT*XI_SCALE*N1
                     PPE( 5,J) = PPE( 5,J) - QT*L*WGT*XI_SCALE*N2
                     PPE( 9,J) = PPE( 9,J) + QT*L*WGT*XI_SCALE*N3
                     PPE(11,J) = PPE(11,J) - QT*L*WGT*XI_SCALE*N4
                  ENDDO
               ENDIF
            ENDIF

            P1  = PRESS(9 ,J)
            P2  = PRESS(10,J)
            X1L = PRESS(11,J)
            X2L = PRESS(12,J)
            IF (X1L >= ZERO) THEN
               IF (DABS(X2L - X1L) <= EPS1) THEN
                  PC = P1
                  PPE( 1,J) = PPE( 1,J) + PC*(ONE - X1L)
                  PPE( 7,J) = PPE( 7,J) + PC*X1L
               ELSE
                  XI_SCALE = (X2L - X1L)/TWO
                  DO I=1,3
                     XI  = XI_SCALE*XI_GAUSS(I) + (X2L + X1L)/TWO
                     WGT = XI_WGT(I)
                     QT  = P1 + (P2-P1)*(XI-X1L)/(X2L-X1L)
                     PPE( 1,J) = PPE( 1,J) + QT*L*WGT*XI_SCALE*(ONE - XI)
                     PPE( 7,J) = PPE( 7,J) + QT*L*WGT*XI_SCALE*XI
                  ENDDO
               ENDIF
            ENDIF

            P1  = PRESS(13,J)
            P2  = PRESS(14,J)
            X1L = PRESS(15,J)
            X2L = PRESS(16,J)
            IF (X1L >= ZERO) THEN
               IF (DABS(X2L - X1L) <= EPS1) THEN
                  PC = P1
                  PPE( 4,J) = PPE( 4,J) + PC*(ONE - X1L)
                  PPE(10,J) = PPE(10,J) + PC*X1L
               ELSE
                  XI_SCALE = (X2L - X1L)/TWO
                  DO I=1,3
                     XI  = XI_SCALE*XI_GAUSS(I) + (X2L + X1L)/TWO
                     WGT = XI_WGT(I)
                     QT  = P1 + (P2-P1)*(XI-X1L)/(X2L-X1L)
                     PPE( 4,J) = PPE( 4,J) + QT*L*WGT*XI_SCALE*(ONE - XI)
                     PPE(10,J) = PPE(10,J) + QT*L*WGT*XI_SCALE*XI
                  ENDDO
               ENDIF
            ENDIF

            P1  = PRESS(17,J)
            P2  = PRESS(18,J)
            X1L = PRESS(19,J)
            X2L = PRESS(20,J)
            IF (X1L >= ZERO) THEN
               IF (DABS(X2L - X1L) <= EPS1) THEN
                  PC = P1
                  PPE( 5,J) = PPE( 5,J) + PC*(ONE - X1L)
                  PPE(11,J) = PPE(11,J) + PC*X1L
               ELSE
                  XI_SCALE = (X2L - X1L)/TWO
                  DO I=1,3
                     XI  = XI_SCALE*XI_GAUSS(I) + (X2L + X1L)/TWO
                     WGT = XI_WGT(I)
                     QT  = P1 + (P2-P1)*(XI-X1L)/(X2L-X1L)
                     PPE( 5,J) = PPE( 5,J) + QT*L*WGT*XI_SCALE*(ONE - XI)
                     PPE(11,J) = PPE(11,J) + QT*L*WGT*XI_SCALE*XI
                  ENDDO
               ENDIF
            ENDIF

            P1  = PRESS(21,J)
            P2  = PRESS(22,J)
            X1L = PRESS(23,J)
            X2L = PRESS(24,J)
            IF (X1L >= ZERO) THEN
               IF (DABS(X2L - X1L) <= EPS1) THEN
                  PC = P1
                  PPE( 6,J) = PPE( 6,J) + PC*(ONE - X1L)
                  PPE(12,J) = PPE(12,J) + PC*X1L
               ELSE
                  XI_SCALE = (X2L - X1L)/TWO
                  DO I=1,3
                     XI  = XI_SCALE*XI_GAUSS(I) + (X2L + X1L)/TWO
                     WGT = XI_WGT(I)
                     QT  = P1 + (P2-P1)*(XI-X1L)/(X2L-X1L)
                     PPE( 6,J) = PPE( 6,J) + QT*L*WGT*XI_SCALE*(ONE - XI)
                     PPE(12,J) = PPE(12,J) + QT*L*WGT*XI_SCALE*XI
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDIF

! **********************************************************************************************************************************
! Stress recovery matrices

      DO I=1,3
         DO J=1,6
            B1(I,J) = ZERO
            B2(I,J) = ZERO
         ENDDO
      ENDDO

      IF (DABS(AREA) > EPS1) THEN
         B1(1,1) = -ONE/AREA
      ENDIF

      B1(2,5) = -DELTA12
      B1(2,6) = -DELTA1

      B1(3,5) =  DELTA2
      B1(3,6) =  DELTA12

      B2(1,2) =  DELTA1*L
      B2(1,3) = -DELTA12*L
      B2(1,5) = -DELTA12
      B2(1,6) = -DELTA1

      B2(2,2) = -DELTA12*L
      B2(2,3) =  DELTA2*L
      B2(2,5) =  DELTA2
      B2(2,6) =  DELTA12

      IF (DABS(JTOR) > EPS1) THEN
         B2(3,4) = -SCOEFF/JTOR
      ENDIF

      CALL MATMULT_FFF ( B1, KAA, 3, 6, 6, S11 )
      CALL MATMULT_FFF ( B1, KAB, 3, 6, 6, S12 )
      CALL MATMULT_FFF ( B2, KAA, 3, 6, 6, S21 )
      CALL MATMULT_FFF ( B2, KAB, 3, 6, 6, S22 )

      DO I=1,3
         DO J=1,6
            SE1(I,J,1) = S11(I,J)
            SE2(I,J,1) = S21(I,J)
         ENDDO
         DO J=7,12
            SE1(I,J,1) = S12(I,J-6)
            SE2(I,J,1) = S22(I,J-6)
         ENDDO
      ENDDO

      IF (NTSUB > 0) THEN
         DO I=1,3
            DO J=1,5
               BT1(I,J) = ZERO
               BT2(I,J) = ZERO
            ENDDO
         ENDDO

         CALL MATMULT_FFF ( S11, ABAR, 3, 6, 5, BT1 )
         CALL MATMULT_FFF ( S21, ABAR, 3, 6, 5, BT2 )

         CALL MATMULT_FFF ( BT1, TPRIME, 3, 5, NTSUB, DUM1 )
         CALL MATMULT_FFF ( BT2, TPRIME, 3, 5, NTSUB, DUM2 )
         DO I=1,3
            DO J=1,NTSUB
               STE1(I,J,1) = DUM1(I,J)
               STE2(I,J,1) = DUM2(I,J)
            ENDDO
         ENDDO
      ENDIF

! **********************************************************************************************************************************
! Reuse the existing 1D geometric stiffness pattern so nonlinear branches do not fail.

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN

         CALL ELMDIS
         CALL CALC_ELEM_NODE_FORCES

         M1A = -PEL(6)
         M2A =  PEL(5)
         M1B = -PEL(6) + PEL(2)*L
         M2B =  PEL(5) + PEL(3)*L
         V1  = -PEL(2)
         V2  = -PEL(3)
         FX  = -PEL(1)

         IF (ART_KED == 'Y') THEN
            KED( 1, 1) = ART_TRAN_KED
            KED( 4, 4) = ART_ROT_KED
            KED( 7, 7) = ART_TRAN_KED
            KED(10,10) = ART_ROT_KED
         ENDIF

         C01 = FX/L

         KED( 2, 2) =  (SIX/FIVE)*FX/L
         KED( 2, 4) =  M2B/L
         KED( 2, 6) =  FX/TEN
         KED( 2, 8) = -KED( 2, 2)
         KED( 2,10) =  M2A/L
         KED( 2,12) =  KED( 2, 6)

         KED( 3, 3) =  KED( 2, 2)
         KED( 3, 4) =  M1B/L
         KED( 3, 5) = -KED( 2, 6)
         KED( 3, 9) = -KED( 2, 2)
         KED( 3,10) =  M1A/L
         KED( 3,11) = -KED( 2, 6)

         KED( 4, 4) =  JTOR*FX/(L*L*AREA)
         KED( 4, 5) = -V1*L/SIX
         KED( 4, 6) = -V2*L/SIX
         KED( 4, 8) = -M2B/L
         KED( 4, 9) = -M1B/L
         KED( 4,10) = -KED( 4, 4)
         KED( 4,11) = -KED( 4, 5)
         KED( 4,12) = -KED( 4, 6)

         KED( 5, 5) =  FOUR*FX*L/(THREE*TEN)
         KED( 5, 9) =  KED( 2, 6)
         KED( 5,10) = -KED( 4, 5)
         KED( 5,11) = -FX*L/(THREE*TEN)

         KED( 6, 6) =  KED( 5, 5)
         KED( 6, 8) = -KED( 2, 6)
         KED( 6,10) = -KED( 4, 6)
         KED( 6,12) =  KED( 5,11)

         KED( 8, 8) =  KED( 2, 2)
         KED( 8,10) = -KED( 2,10)
         KED( 8,12) = -KED( 2, 6)

         KED( 9, 9) =  KED( 2, 2)
         KED( 9,10) = -KED( 3,10)
         KED( 9,11) =  KED( 2, 6)

         KED(10,10) =  KED( 4, 4)
         KED(10,11) =  KED( 4, 5)
         KED(10,12) =  KED( 4, 6)

         KED(11,11) =  KED( 5, 5)
         KED(12,12) =  KED( 6, 6)

         DO I=2,12
            DO J=1,I-1
               KED(I,J) = KED(J,I)
            ENDDO
         ENDDO

      ENDIF

      RETURN

! **********************************************************************************************************************************

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DEBUG_BEAM ( WHAT )

      USE PENTIUM_II_KIND

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: WHAT

      IF (WHAT == 1) THEN
         WRITE(F06,*)
         WRITE(F06,1997)
         WRITE(F06,'(A,I8)') 'In subr BEAM with BEAM element ',EID
         WRITE(F06,*) '--------------------------------------'
         WRITE(F06,1998) 'L       = ',L
         WRITE(F06,1998) 'AREA    = ',AREA
         WRITE(F06,1998) 'I1      = ',I1
         WRITE(F06,1998) 'I2      = ',I2
         WRITE(F06,1998) 'I12     = ',I12
         WRITE(F06,1998) 'JTOR    = ',JTOR
         WRITE(F06,1998) 'K1      = ',K1
         WRITE(F06,1998) 'K2      = ',K2
         WRITE(F06,1998) 'E       = ',E
         WRITE(F06,1998) 'G       = ',G
         WRITE(F06,1998) 'SCOEFF  = ',SCOEFF
         WRITE(F06,1998) 'ALPHA   = ',ALPHA
         WRITE(F06,1998) 'TREF    = ',TREF
      ELSE IF (WHAT == 2) THEN
         WRITE(F06,1999) 'plane1 (UY,RZ) uses I1/K2; phi1 = ',PHI1
         WRITE(F06,1999) 'plane2 (UZ,RY) uses I2/K1; phi2 = ',PHI2
         WRITE(F06,1998) 'FAC1    = ',FAC1
         WRITE(F06,1998) 'FAC2    = ',FAC2
         WRITE(F06,1998) 'DEN     = ',DEN
         WRITE(F06,1998) 'DELTA1  = ',DELTA1
         WRITE(F06,1998) 'DELTA2  = ',DELTA2
         WRITE(F06,1998) 'DELTA12 = ',DELTA12
         WRITE(F06,*)
         WRITE(F06,1997)
         WRITE(F06,*)
      ENDIF

 1997 FORMAT('************************************************************')
 1998 FORMAT(A, 1ES14.6)
 1999 FORMAT(A, 1ES14.6)

      END SUBROUTINE DEBUG_BEAM

      END SUBROUTINE BEAM

!--- cbeam_add end ---!
