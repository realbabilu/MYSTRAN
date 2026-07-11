! #################################################################################################################################
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
      SUBROUTINE CALC_K6ROT ()

! Builds the stiffness matrix for a spring connecting the drilling DOF to the translational DOFs of the adjacent nodes
! of the element. Adds that to the existing stiffness matrix KE in the element coordinate system.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  TYPE, ELGP, INTL_MID, XEL, SHELL_A, KE
      USE PARAMS, ONLY                :  K6ROT
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE SCONTR, ONLY                :  MAX_ORDER_GAUSS
      USE CROSS_Interface
      USE MATMULT_FFF_T_Interface

      IMPLICIT NONE

      REAL(DOUBLE)                    :: KROT(6*ELGP,6*ELGP)
      REAL(DOUBLE)                    :: N(3)
      REAL(DOUBLE)                    :: X_PREV(3)
      REAL(DOUBLE)                    :: X_NEXT(3)
      REAL(DOUBLE)                    :: TERM_PREV(3)
      REAL(DOUBLE)                    :: TERM_NEXT(3)
      REAL(DOUBLE)                    :: B(6*ELGP)
      REAL(DOUBLE)                    :: STIFFNESS
      REAL(DOUBLE)                    :: AREA
      REAL(DOUBLE)                    :: DETJ
      INTEGER(LONG)                   :: I,J
      INTEGER(LONG)                   :: GP
      INTEGER(LONG)                   :: GP_PREV
      INTEGER(LONG)                   :: GP_NEXT
      REAL(DOUBLE)                    :: JAC(2,2)
      REAL(DOUBLE)                    :: JACI(2,2)
      REAL(DOUBLE)                    :: X2E
      REAL(DOUBLE)                    :: Y3E
      REAL(DOUBLE)                    :: XSD(4)
      REAL(DOUBLE)                    :: YSD(4)
      REAL(DOUBLE)                    :: HHH(MAX_ORDER_GAUSS)
      REAL(DOUBLE)                    :: SSS(MAX_ORDER_GAUSS)

! **********************************************************************************************************************************

      IF (INTL_MID(2) > 0) THEN

         AREA = ZERO

         IF ((TYPE(1:5) == "QUAD4")) THEN

            XSD(1) = XEL(1,1) - XEL(2,1)
            XSD(2) = XEL(2,1) - XEL(3,1)
            XSD(3) = XEL(3,1) - XEL(4,1)
            XSD(4) = XEL(4,1) - XEL(1,1)

            YSD(1) = XEL(1,2) - XEL(2,2)
            YSD(2) = XEL(2,2) - XEL(3,2)
            YSD(3) = XEL(3,2) - XEL(4,2)
            YSD(4) = XEL(4,2) - XEL(1,2)

            CALL ORDER_GAUSS ( 2, SSS, HHH )
            DO I=1,2
               DO J=1,2
                  CALL JAC2D ( SSS(I), SSS(J), XSD, YSD, 'N', JAC, JACI, DETJ )
                  AREA = AREA + HHH(I)*HHH(J)*DETJ
               ENDDO
            ENDDO

         ELSEIF (TYPE(1:5) == "TRIA3") THEN

            X2E  = XEL(2,1)
            Y3E  = XEL(3,2)
            AREA = X2E*Y3E

         ENDIF

         STIFFNESS = 10.0D0**(-6.0D0) * K6ROT * SHELL_A(3,3) * ABS(AREA)

         N = [ZERO, ZERO, ONE]

         DO GP=1,ELGP

            B = ZERO

            GP_PREV = GP - 1
            IF (GP_PREV < 1) THEN
               GP_PREV = ELGP
            ENDIF

            GP_NEXT = GP + 1
            IF (GP_NEXT > ELGP) THEN
               GP_NEXT = 1
            ENDIF

            X_PREV = XEL(GP_PREV,1:3) - XEL(GP,1:3)
            X_NEXT = XEL(GP_NEXT,1:3) - XEL(GP,1:3)

            CALL CROSS(N, X_PREV, TERM_PREV)
            TERM_PREV = TERM_PREV * 1 / (2 * (X_PREV(1)**2 + X_PREV(2)**2 + X_PREV(3)**2) )

            B((GP_PREV - 1) * 6 + 1) = -TERM_PREV(1)
            B((GP_PREV - 1) * 6 + 2) = -TERM_PREV(2)
            B((GP_PREV - 1) * 6 + 3) = -TERM_PREV(3)

            CALL CROSS(N, X_NEXT, TERM_NEXT)
            TERM_NEXT = TERM_NEXT * 1 / (2 * (X_NEXT(1)**2 + X_NEXT(2)**2 + X_NEXT(3)**2) )

            B((GP_NEXT - 1) * 6 + 1) = -TERM_NEXT(1)
            B((GP_NEXT - 1) * 6 + 2) = -TERM_NEXT(2)
            B((GP_NEXT - 1) * 6 + 3) = -TERM_NEXT(3)

            B((GP - 1) * 6 + 1) = TERM_PREV(1) + TERM_NEXT(1)
            B((GP - 1) * 6 + 2) = TERM_PREV(2) + TERM_NEXT(2)
            B((GP - 1) * 6 + 3) = TERM_PREV(3) + TERM_NEXT(3)
            B((GP - 1) * 6 + 4) = N(1)
            B((GP - 1) * 6 + 5) = N(2)
            B((GP - 1) * 6 + 6) = N(3)

            CALL MATMULT_FFF_T(B, B, 1, 6*ELGP, 6*ELGP, KROT)
            KROT = KROT * STIFFNESS

            KE(1:6*ELGP, 1:6*ELGP) = KE(1:6*ELGP, 1:6*ELGP) + KROT(1:6*ELGP, 1:6*ELGP)

         ENDDO

      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE CALC_K6ROT
