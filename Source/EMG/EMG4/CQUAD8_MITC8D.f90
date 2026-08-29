! #################################################################################################################################
! CQUAD8 MITC8D shell for PARAM,QUAD8TYP,MITC8D.

      SUBROUTINE CQUAD8_MITC8D ( OPT, INT_ELEM_ID )

! Reuse the base MITC8 path, then replace the non-D drilling penalty with the
! pure Wilson theta_z penalty used by the Python MITC8D variant.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  MAX_ORDER_GAUSS
      USE MODEL_STUF, ONLY            :  ELGP, KE, XEL
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE MITC8_Interface
      USE ORDER_GAUSS_Interface
      USE MITC_SHAPE_FUNCTIONS_Interface
      USE MITC8_CARTESIAN_LOCAL_BASIS_Interface
      USE MITC_DETJ_Interface
      USE MITC_ELASTICITY_Interface
      USE MATMULT_FFF_T_Interface

      IMPLICIT NONE

      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      INTEGER(LONG), PARAMETER        :: IORD_IJ = 3
      INTEGER(LONG), PARAMETER        :: IORD_K  = 2
      REAL(DOUBLE), PARAMETER         :: BETA_DRILL = 1.0D-6

      INTEGER(LONG)                   :: I, J, K
      REAL(DOUBLE)                    :: E(6,6), GDRILL
      REAL(DOUBLE)                    :: HH_IJ(MAX_ORDER_GAUSS), SS_IJ(MAX_ORDER_GAUSS)
      REAL(DOUBLE)                    :: HH_K(MAX_ORDER_GAUSS),  SS_K(MAX_ORDER_GAUSS)
      REAL(DOUBLE)                    :: R, S, T, DETJ, INTFAC
      REAL(DOUBLE)                    :: BD_STD(1,6*ELGP), BD_PURE(1,6*ELGP), KD(6*ELGP,6*ELGP)

      CALL MITC8 ( OPT, INT_ELEM_ID )

      IF (OPT(4) /= 'Y') RETURN

      E = MITC_ELASTICITY()
      GDRILL = E(4,4)

      CALL ORDER_GAUSS ( IORD_IJ, SS_IJ, HH_IJ )
      CALL ORDER_GAUSS ( IORD_K , SS_K , HH_K  )

      DO I=1,IORD_IJ
         DO J=1,IORD_IJ
            R = SS_IJ(I)
            S = SS_IJ(J)
            CALL MITC8_STD_DRILL_B  ( R, S, BD_STD  )
            CALL MITC8_PURE_DRILL_B ( R, S, BD_PURE )
            CALL MATMULT_FFF_T ( BD_PURE, BD_PURE, 1, 6*ELGP, 6*ELGP, KD )
            CALL MATMULT_FFF_T ( BD_STD , BD_STD , 1, 6*ELGP, 6*ELGP, KE )
            KD = KD - KE
            DO K=1,IORD_K
               T = SS_K(K)
               DETJ = MITC_DETJ ( R, S, T )
               INTFAC = DETJ*HH_IJ(I)*HH_IJ(J)*HH_K(K)
               KE(1:6*ELGP,1:6*ELGP) = KE(1:6*ELGP,1:6*ELGP) + BETA_DRILL*GDRILL*KD*INTFAC
            ENDDO
         ENDDO
      ENDDO

      RETURN

      CONTAINS

      SUBROUTINE MITC8_STD_DRILL_B ( R, S, BDOUT )

      REAL(DOUBLE), INTENT(IN)        :: R, S
      REAL(DOUBLE), INTENT(OUT)       :: BDOUT(1,6*ELGP)

      INTEGER(LONG)                   :: II
      REAL(DOUBLE)                    :: PSH(ELGP), DPSHG(2,ELGP), CLB(3,3)
      REAL(DOUBLE)                    :: XI_LOC(ELGP), ETA_LOC(ELGP)
      REAL(DOUBLE)                    :: J11, J12, J21, J22, DET2
      REAL(DOUBLE)                    :: DNDX, DNDY

      BDOUT = ZERO
      CALL MITC_SHAPE_FUNCTIONS ( R, S, PSH, DPSHG )
      CLB = MITC8_CARTESIAN_LOCAL_BASIS ( R, S )

      DO II=1,ELGP
         XI_LOC(II)  = DOT_PRODUCT( XEL(II,:), CLB(:,1) )
         ETA_LOC(II) = DOT_PRODUCT( XEL(II,:), CLB(:,2) )
      ENDDO

      J11 = DOT_PRODUCT( DPSHG(1,:), XI_LOC  )
      J12 = DOT_PRODUCT( DPSHG(1,:), ETA_LOC )
      J21 = DOT_PRODUCT( DPSHG(2,:), XI_LOC  )
      J22 = DOT_PRODUCT( DPSHG(2,:), ETA_LOC )
      DET2 = J11*J22 - J12*J21
      IF (DABS(DET2) < 1.0D-14) RETURN

      DO II=1,ELGP
         DNDX = ( J22*DPSHG(1,II) - J12*DPSHG(2,II) ) / DET2
         DNDY = (-J21*DPSHG(1,II) + J11*DPSHG(2,II) ) / DET2
         BDOUT(1,6*(II-1)+1) = -0.5D0 * DNDY
         BDOUT(1,6*(II-1)+2) =  0.5D0 * DNDX
         BDOUT(1,6*(II-1)+6) =  PSH(II)
      ENDDO

      END SUBROUTINE MITC8_STD_DRILL_B

      SUBROUTINE MITC8_PURE_DRILL_B ( R, S, BDOUT )

      REAL(DOUBLE), INTENT(IN)        :: R, S
      REAL(DOUBLE), INTENT(OUT)       :: BDOUT(1,6*ELGP)

      INTEGER(LONG)                   :: II
      REAL(DOUBLE)                    :: PSH(ELGP), DPSHG(2,ELGP)

      BDOUT = ZERO
      CALL MITC_SHAPE_FUNCTIONS ( R, S, PSH, DPSHG )
      DO II=1,ELGP
         BDOUT(1,6*(II-1)+6) = PSH(II)
      ENDDO

      END SUBROUTINE MITC8_PURE_DRILL_B

      END SUBROUTINE CQUAD8_MITC8D
