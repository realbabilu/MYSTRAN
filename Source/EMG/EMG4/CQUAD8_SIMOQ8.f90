! #################################################################################################################################
! CQUAD8 Python-aligned Simo1993 Q8 shell for PARAM,QUAD8TYP,SIMOQ8.

      SUBROUTINE CQUAD8_SIMOQ8 ( OPT, INT_ELEM_ID )

! Ported from:
!   D:\18a\python\quadratic\Simo1993_Q8_ShellElement_v1p8_standalone.py
!   D:\18a\python\quadratic\Simo1993_Q8_thermal_buckling.py
!
! Static stiffness path:
!   Q8 serendipity Simo/Fox director kinematics + Hughes-Brezzi drilling,
!   3x3 membrane/bending/drilling/shear, and one EAS shear bubble
!   phi=(1-r^2)(1-s^2) condensed at element level.
!
! This branch is the active Simo Q8 path and follows the Python benchmark
! conventions more closely, including the center-constant geometric stiffness
! used by the thermal-buckling helper.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_STRESS_POINTS, SOL_NAME
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  ALPVEC, BGRID, DT, EID, ELGP, GRID_SNORM, KE, KED, ME, BE1, BE2, BE3, EPROP, FCONV,    &
                                         MASS_PER_UNIT_AREA, NUM_EMG_FATAL_ERRS, PCOMP_PROPS, PPE, PRESS, PTE, SHELL_A,        &
                                         SHELL_D, SHELL_T, TREF, UEL, XEB
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO
      USE PARAMS, ONLY                :  COUPMASS
      USE ELMDIS_Interface
      USE OUTA_HERE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CQUAD8_SIMOQ8'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      INTEGER(LONG)                   :: I, J, GP, IA, IB, K, L, JSUB
      REAL(DOUBLE)                    :: XYZ(8,3), NORMALS(8,3)
      REAL(DOUBLE)                    :: KDD(48,48), KOUT(48,48), KDA(48), KAA
      REAL(DOUBLE)                    :: BM(3,48), BB(3,48), BS(2,48), BD(1,48), BSE(2)
      REAL(DOUBLE)                    :: GP3(3), W3(3), R, S, WT, JAC, CDRILL, FAC
      REAL(DOUBLE)                    :: M1(8,8), N8(8), DN8(2,8), MASS_ELEM, MASS_NODE
      REAL(DOUBLE)                    :: UNIT_PPE(48), UNIT_PTE(48), DXDR(3), DXDS(3), SURF_VEC(3), TBAR
      REAL(DOUBLE)                    :: CTE(3), THERMAL_RESULTANT(3)
      REAL(DOUBLE)                    :: DLOC(2,8), SIG0(2,2), DUM28(2,8), KG8(8,8), STRAIN0(3), N0(3)
      REAL(DOUBLE)                    :: G1(3), G2(3), A11, A22, A12, DET, AI11, AI22, AI12
      INTEGER(LONG)                   :: KI, KJ

      IF (ELGP /= 8) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9001) SUBR_NAME, EID, ELGP
         WRITE(F06,9001) SUBR_NAME, EID, ELGP
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IF (PCOMP_PROPS == 'Y') THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,*) ' *ERROR: Code not written for composite material with PARAM,QUAD8TYP,SIMOQ8'
         WRITE(F06,*) ' *ERROR: Code not written for composite material with PARAM,QUAD8TYP,SIMOQ8'
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL LOAD_BASIC_COORDS_Q8 ( XYZ )
      CALL CALC_NODAL_NORMALS_Q8 ( XYZ, NORMALS )

      GP3 = (/-DSQRT(3.0D0/5.0D0), ZERO, DSQRT(3.0D0/5.0D0)/)
      W3  = (/5.0D0/9.0D0, 8.0D0/9.0D0, 5.0D0/9.0D0/)

      IF (OPT(1) == 'Y') THEN
         M1 = ZERO
         MASS_ELEM = ZERO
         DO I=1,3
            DO J=1,3
               R = GP3(I)
               S = GP3(J)
               WT = W3(I)*W3(J)
               CALL SHAPE_Q8 ( R, S, N8, DN8 )
               CALL BM_Q8_AT ( XYZ, R, S, BM, JAC )
               MASS_ELEM = MASS_ELEM + MASS_PER_UNIT_AREA*WT*JAC
               DO K=1,8
                  DO L=1,8
                     M1(K,L) = M1(K,L) + N8(K)*N8(L)*MASS_PER_UNIT_AREA*WT*JAC
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         ME = ZERO
         IF ((SOL_NAME(1:5) == 'MODES') .AND. (COUPMASS > 0)) THEN
            DO K=1,8
               DO L=1,8
                  DO IA=1,3
                     ME(6*(K-1)+IA,6*(L-1)+IA) = M1(K,L)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            MASS_NODE = MASS_ELEM/8.0D0
            DO K=1,8
               DO IA=1,3
                  ME(6*(K-1)+IA,6*(K-1)+IA) = MASS_NODE
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      IF (OPT(2) == 'Y') THEN
         UNIT_PTE = ZERO
         DO I=1,3
            DO J=1,3
               R = GP3(I)
               S = GP3(J)
               WT = W3(I)*W3(J)
               CALL BM_Q8_AT ( XYZ, R, S, BM, JAC )
               CTE(1) = ALPVEC(1,1)
               CTE(2) = ALPVEC(2,1)
               CTE(3) = ALPVEC(4,1)
               THERMAL_RESULTANT = MATMUL(SHELL_A, CTE)
               UNIT_PTE = UNIT_PTE + MATMUL( TRANSPOSE(BM), THERMAL_RESULTANT ) * WT * JAC
            ENDDO
         ENDDO
         DO JSUB=1,SIZE(PTE,2)
            TBAR = ZERO
            DO J=1,8
               TBAR = TBAR + DT(J,JSUB)
            ENDDO
            TBAR = TBAR / 8.0D0 - TREF(1)
            PTE(1:48,JSUB) = UNIT_PTE(1:48) * TBAR
         ENDDO
      ENDIF

      IF (OPT(3) == 'Y') THEN
         GP = 1
         DO I=1,2
            DO J=1,2
               GP = GP + 1
               R = 0.577350269189626D0
               S = 0.577350269189626D0
               IF (I == 1) R = -R
               IF (J == 1) S = -S
               CALL BM_Q8_AT ( XYZ, R, S, BM, JAC )
               CALL BB_Q8_AT ( XYZ, NORMALS, R, S, BB, JAC )
               CALL BS_Q8_AT ( XYZ, NORMALS, R, S, BS, JAC )
               IF (GP <= MAX_STRESS_POINTS) THEN
                  BE1(1:3,1:48,GP) = BM
                  BE2(1:3,1:48,GP) = BB
                  BE3(1:2,1:48,GP) = BS
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (OPT(4) == 'Y') THEN
         KDD = ZERO
         KDA = ZERO
         KAA = ZERO

         CDRILL = 2.0D-2*SHELL_T(1,1)/(5.0D0/6.0D0)

         DO I=1,3
            DO J=1,3
               R = GP3(I)
               S = GP3(J)
               WT = W3(I)*W3(J)

               CALL BM_Q8_AT ( XYZ, R, S, BM, JAC )
               CALL BB_Q8_AT ( XYZ, NORMALS, R, S, BB, JAC )
               CALL BS_Q8_AT ( XYZ, NORMALS, R, S, BS, JAC )
               CALL BDRILL_Q8_AT ( XYZ, NORMALS, R, S, BD, JAC )
               CALL EAS1_SHEAR_Q8_AT ( XYZ, R, S, BSE )

               KDD = KDD + WT*JAC*MATMUL(TRANSPOSE(BM), MATMUL(SHELL_A, BM))
               KDD = KDD + WT*JAC*MATMUL(TRANSPOSE(BB), MATMUL(SHELL_D, BB))
               KDD = KDD + WT*JAC*MATMUL(TRANSPOSE(BS), MATMUL(SHELL_T, BS))
               KDD = KDD + WT*JAC*CDRILL*MATMUL(TRANSPOSE(BD), BD)

               KDA = KDA + WT*JAC*MATMUL(TRANSPOSE(BS), MATMUL(SHELL_T, BSE))
               KAA = KAA + WT*JAC*DOT_PRODUCT(BSE, MATMUL(SHELL_T, BSE))
            ENDDO
         ENDDO

         KOUT = KDD
         IF (DABS(KAA) > 1.0D-14) THEN
            DO IA=1,48
               DO IB=1,48
                  KOUT(IA,IB) = KOUT(IA,IB) - KDA(IA)*KDA(IB)/KAA
               ENDDO
            ENDDO
         ENDIF

         DO IA=1,48
            DO IB=1,48
               FAC = 0.5D0*(KOUT(IA,IB) + KOUT(IB,IA))
               KE(IA,IB) = FAC
            ENDDO
         ENDDO
      ENDIF

      IF (OPT(5) == 'Y') THEN
         UNIT_PPE = ZERO
         DO I=1,3
            DO J=1,3
               R = GP3(I)
               S = GP3(J)
               WT = W3(I)*W3(J)
               CALL SHAPE_Q8 ( R, S, N8, DN8 )
               DXDR = ZERO
               DXDS = ZERO
               DO K=1,8
                  DXDR(:) = DXDR(:) + DN8(1,K)*XYZ(K,:)
                  DXDS(:) = DXDS(:) + DN8(2,K)*XYZ(K,:)
               ENDDO
               CALL CROSS3 ( DXDR, DXDS, SURF_VEC )
               DO K=1,8
                  UNIT_PPE(6*(K-1)+1) = UNIT_PPE(6*(K-1)+1) + N8(K)*SURF_VEC(1)*WT
                  UNIT_PPE(6*(K-1)+2) = UNIT_PPE(6*(K-1)+2) + N8(K)*SURF_VEC(2)*WT
                  UNIT_PPE(6*(K-1)+3) = UNIT_PPE(6*(K-1)+3) + N8(K)*SURF_VEC(3)*WT
               ENDDO
            ENDDO
         ENDDO
         DO J=1,SIZE(PPE,2)
            PPE(1:48,J) = PPE(1:48,J) + UNIT_PPE(1:48)*PRESS(3,J)
         ENDDO
      ENDIF

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN
         CALL ELMDIS
         CALL BM_Q8_AT ( XYZ, ZERO, ZERO, BM, JAC )
         STRAIN0 = MATMUL(BM, UEL(1:48))
         N0 = MATMUL(SHELL_A, STRAIN0)

         SIG0 = ZERO
         SIG0(1,1) = N0(1)
         SIG0(2,2) = N0(2)
         SIG0(1,2) = N0(3)
         SIG0(2,1) = N0(3)

         KG8 = ZERO
         DO I=1,3
            DO J=1,3
               R = GP3(I)
               S = GP3(J)
               WT = W3(I)*W3(J)

               CALL SHAPE_Q8 ( R, S, N8, DN8 )
               G1 = MATMUL(DN8(1,:), XYZ)
               G2 = MATMUL(DN8(2,:), XYZ)
               CALL CROSS3 ( G1, G2, SURF_VEC )
               JAC = VNORM(SURF_VEC)

               A11 = DOT_PRODUCT(G1,G1)
               A22 = DOT_PRODUCT(G2,G2)
               A12 = DOT_PRODUCT(G1,G2)
               DET = A11*A22 - A12*A12
               IF (DABS(DET) <= 1.0D-30) CYCLE

               AI11 =  A22/DET
               AI22 =  A11/DET
               AI12 = -A12/DET
               DLOC(1,:) = AI11*DN8(1,:) + AI12*DN8(2,:)
               DLOC(2,:) = AI12*DN8(1,:) + AI22*DN8(2,:)

               DUM28 = MATMUL(SIG0, DLOC)
               KG8 = KG8 + MATMUL(TRANSPOSE(DLOC), DUM28) * WT * JAC
            ENDDO
         ENDDO

         KED(1:48,1:48) = ZERO
         DO I=1,8
            DO J=1,8
               KI = 6*(I-1)
               KJ = 6*(J-1)
               KED(KI+1,KJ+1) = KG8(I,J)
               KED(KI+2,KJ+2) = KG8(I,J)
               KED(KI+3,KJ+3) = KG8(I,J)
            ENDDO
         ENDDO
      ENDIF

      RETURN

 9001 FORMAT(' *ERROR: ',A,' expects ELGP=8 for element ',I8,' but got ',I8)

      CONTAINS

      SUBROUTINE LOAD_BASIC_COORDS_Q8 ( XYZOUT )
      REAL(DOUBLE), INTENT(OUT) :: XYZOUT(8,3)
      INTEGER(LONG) :: II, JJ
      DO II=1,8
         DO JJ=1,3
            XYZOUT(II,JJ) = XEB(II,JJ)
         ENDDO
      ENDDO
      END SUBROUTINE LOAD_BASIC_COORDS_Q8

      SUBROUTINE SHAPE_Q8 ( XI, ETA, NVAL, DN )
      REAL(DOUBLE), INTENT(IN)  :: XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: NVAL(8), DN(2,8)
      NVAL(1) = 0.25D0*(ONE-XI)*(ONE-ETA)*(-XI-ETA-ONE)
      NVAL(2) = 0.25D0*(ONE+XI)*(ONE-ETA)*( XI-ETA-ONE)
      NVAL(3) = 0.25D0*(ONE+XI)*(ONE+ETA)*( XI+ETA-ONE)
      NVAL(4) = 0.25D0*(ONE-XI)*(ONE+ETA)*(-XI+ETA-ONE)
      NVAL(5) = 0.5D0*(ONE-XI*XI)*(ONE-ETA)
      NVAL(6) = 0.5D0*(ONE+XI)*(ONE-ETA*ETA)
      NVAL(7) = 0.5D0*(ONE-XI*XI)*(ONE+ETA)
      NVAL(8) = 0.5D0*(ONE-XI)*(ONE-ETA*ETA)

      DN(1,1) = 0.25D0*(ONE-ETA)*(TWO*XI+ETA)
      DN(1,2) = 0.25D0*(ONE-ETA)*(TWO*XI-ETA)
      DN(1,3) = 0.25D0*(ONE+ETA)*(TWO*XI+ETA)
      DN(1,4) = 0.25D0*(ONE+ETA)*(TWO*XI-ETA)
      DN(1,5) = -XI*(ONE-ETA)
      DN(1,6) = 0.5D0*(ONE-ETA*ETA)
      DN(1,7) = -XI*(ONE+ETA)
      DN(1,8) = -0.5D0*(ONE-ETA*ETA)

      DN(2,1) = 0.25D0*(ONE-XI)*(TWO*ETA+XI)
      DN(2,2) = 0.25D0*(ONE+XI)*(TWO*ETA-XI)
      DN(2,3) = 0.25D0*(ONE+XI)*(TWO*ETA+XI)
      DN(2,4) = 0.25D0*(ONE-XI)*(TWO*ETA-XI)
      DN(2,5) = -0.5D0*(ONE-XI*XI)
      DN(2,6) = -(ONE+XI)*ETA
      DN(2,7) = 0.5D0*(ONE-XI*XI)
      DN(2,8) = -(ONE-XI)*ETA
      END SUBROUTINE SHAPE_Q8

      SUBROUTINE CALC_NODAL_NORMALS_Q8 ( XYZN, NORMS )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3)
      REAL(DOUBLE), INTENT(OUT) :: NORMS(8,3)
      REAL(DOUBLE) :: RS(8,2), NVAL(8), DN(2,8), G1(3), G2(3), N(3), NM, SN(3)
      INTEGER(LONG) :: II, BIDX
      RS(1,:) = (/-ONE, -ONE/)
      RS(2,:) = (/ ONE, -ONE/)
      RS(3,:) = (/ ONE,  ONE/)
      RS(4,:) = (/-ONE,  ONE/)
      RS(5,:) = (/ZERO, -ONE/)
      RS(6,:) = (/ ONE, ZERO/)
      RS(7,:) = (/ZERO,  ONE/)
      RS(8,:) = (/-ONE, ZERO/)
      DO II=1,8
         CALL SHAPE_Q8(RS(II,1), RS(II,2), NVAL, DN)
         G1 = MATMUL(DN(1,:), XYZN)
         G2 = MATMUL(DN(2,:), XYZN)
         CALL CROSS3(G1, G2, N)
         NM = VNORM(N)
         IF (NM <= 1.0D-12) THEN
            CALL SHAPE_Q8(ZERO, ZERO, NVAL, DN)
            G1 = MATMUL(DN(1,:), XYZN)
            G2 = MATMUL(DN(2,:), XYZN)
            CALL CROSS3(G1, G2, N)
            NM = VNORM(N)
         ENDIF
         IF (NM > 1.0D-15) THEN
            NORMS(II,:) = N/NM
         ELSE
            NORMS(II,:) = (/ZERO, ZERO, ONE/)
         ENDIF
         IF (ALLOCATED(GRID_SNORM)) THEN
            BIDX = 0
            IF (II <= SIZE(BGRID)) BIDX = BGRID(II)
            IF ((BIDX > 0) .AND. (BIDX <= SIZE(GRID_SNORM,1))) THEN
               SN = GRID_SNORM(BIDX,:)
               NM = VNORM(SN)
               IF (NM > 1.0D-15) THEN
                  SN = SN/NM
                  IF (DOT_PRODUCT(SN, NORMS(II,:)) < ZERO) SN = -SN
                  NORMS(II,:) = SN
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      END SUBROUTINE CALC_NODAL_NORMALS_Q8

      SUBROUTINE SURFACE_BASIS_Q8 ( XYZN, XI, ETA, G1, G2, E1, E2, E3, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: G1(3), G2(3), E1(3), E2(3), E3(3), JAC
      REAL(DOUBLE) :: NVAL(8), DN(2,8), G3(3), TMP(3), NM
      CALL SHAPE_Q8(XI, ETA, NVAL, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      CALL CROSS3(G1, G2, G3)
      JAC = VNORM(G3)
      IF (JAC > 1.0D-15) THEN
         E3 = G3/JAC
      ELSE
         E3 = (/ZERO, ZERO, ONE/)
      ENDIF
      CALL CROSS3(G2, E3, TMP)
      NM = VNORM(TMP)
      IF (NM < 1.0D-12) THEN
         CALL CROSS3(G1, E3, TMP)
         NM = VNORM(TMP)
      ENDIF
      IF (NM > 1.0D-15) THEN
         E1 = TMP/NM
      ELSE
         E1 = (/ONE, ZERO, ZERO/)
      ENDIF
      CALL CROSS3(E3, E1, E2)
      NM = VNORM(E2)
      IF (NM > 1.0D-15) E2 = E2/NM
      END SUBROUTINE SURFACE_BASIS_Q8

      SUBROUTINE FIXED_FRAME_Q8 ( XYZN, E1F, E2F, E3F )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3)
      REAL(DOUBLE), INTENT(OUT) :: E1F(3), E2F(3), E3F(3)
      REAL(DOUBLE) :: G1(3), G2(3), JAC, ZSGN
      CALL SURFACE_BASIS_Q8(XYZN, ZERO, ZERO, G1, G2, E1F, E2F, E3F, JAC)
      IF ((DABS(E3F(1)) + DABS(E3F(2))) <= 1.0D-12) THEN
         ZSGN = ONE
         IF (E3F(3) < ZERO) ZSGN = -ONE
         E1F = (/ONE, ZERO, ZERO/)
         E2F = (/ZERO, ZSGN, ZERO/)
         E3F = (/ZERO, ZERO, ZSGN/)
      ENDIF
      END SUBROUTINE FIXED_FRAME_Q8

      SUBROUTINE COV_MAP_Q8 ( XYZN, XI, ETA, E1F, E2F, G1, G2, C, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), XI, ETA, E1F(3), E2F(3)
      REAL(DOUBLE), INTENT(OUT) :: G1(3), G2(3), C(4), JAC
      REAL(DOUBLE) :: E1(3), E2(3), E3(3), A11, A22, A12, DET, AI11, AI22, AI12, GC1(3), GC2(3)
      CALL SURFACE_BASIS_Q8(XYZN, XI, ETA, G1, G2, E1, E2, E3, JAC)
      A11 = DOT_PRODUCT(G1,G1)
      A22 = DOT_PRODUCT(G2,G2)
      A12 = DOT_PRODUCT(G1,G2)
      DET = A11*A22 - A12*A12
      IF (DABS(DET) <= 1.0D-30) THEN
         C = ZERO
         RETURN
      ENDIF
      AI11 =  A22/DET
      AI22 =  A11/DET
      AI12 = -A12/DET
      GC1 = AI11*G1 + AI12*G2
      GC2 = AI12*G1 + AI22*G2
      C(1) = DOT_PRODUCT(E1F, GC1)
      C(2) = DOT_PRODUCT(E1F, GC2)
      C(3) = DOT_PRODUCT(E2F, GC1)
      C(4) = DOT_PRODUCT(E2F, GC2)
      END SUBROUTINE COV_MAP_Q8

      SUBROUTINE TENSOR_PHYS_Q8 ( V11, V22, V12, C, VOUT )
      REAL(DOUBLE), INTENT(IN)  :: V11(3), V22(3), V12(3), C(4)
      REAL(DOUBLE), INTENT(OUT) :: VOUT(3,3)
      VOUT(1,:) = C(1)*C(1)*V11 + C(2)*C(2)*V22 + TWO*C(1)*C(2)*V12
      VOUT(2,:) = C(3)*C(3)*V11 + C(4)*C(4)*V22 + TWO*C(3)*C(4)*V12
      VOUT(3,:) = TWO*(C(1)*C(3)*V11 + C(2)*C(4)*V22 + (C(1)*C(4)+C(2)*C(3))*V12)
      END SUBROUTINE TENSOR_PHYS_Q8

      SUBROUTINE BM_Q8_AT ( XYZN, XI, ETA, BMOUT, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BMOUT(3,48), JAC
      REAL(DOUBLE) :: NVAL(8), DN(2,8), E1F(3), E2F(3), E3F(3), G1(3), G2(3), C(4), VP(3,3)
      INTEGER(LONG) :: II, COL
      CALL SHAPE_Q8(XI, ETA, NVAL, DN)
      CALL FIXED_FRAME_Q8(XYZN, E1F, E2F, E3F)
      CALL COV_MAP_Q8(XYZN, XI, ETA, E1F, E2F, G1, G2, C, JAC)
      BMOUT = ZERO
      DO II=1,8
         COL = (II-1)*6
         CALL TENSOR_PHYS_Q8(DN(1,II)*G1, DN(2,II)*G2, 0.5D0*(DN(1,II)*G2 + DN(2,II)*G1), C, VP)
         BMOUT(1:3,COL+1:COL+3) = VP
      ENDDO
      END SUBROUTINE BM_Q8_AT

      SUBROUTINE BB_Q8_AT ( XYZN, NORMS, XI, ETA, BBOUT, JAC, NORMS_EXT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), NORMS(8,3), XI, ETA
      REAL(DOUBLE), INTENT(IN), OPTIONAL :: NORMS_EXT(8,3)
      REAL(DOUBLE), INTENT(OUT) :: BBOUT(3,48), JAC
      REAL(DOUBLE) :: NVAL(8), DN(2,8), E1F(3), E2F(3), E3F(3), G1(3), G2(3), C(4), VP(3,3)
      REAL(DOUBLE) :: T1(3), T2(3), T0(3), CG1(3), CG2(3)
      REAL(DOUBLE) :: NORMS_LOC(8,3)
      INTEGER(LONG) :: II, COL
      NORMS_LOC = NORMS
      IF (PRESENT(NORMS_EXT)) NORMS_LOC = NORMS_EXT
      CALL SHAPE_Q8(XI, ETA, NVAL, DN)
      CALL FIXED_FRAME_Q8(XYZN, E1F, E2F, E3F)
      CALL COV_MAP_Q8(XYZN, XI, ETA, E1F, E2F, G1, G2, C, JAC)
      T1 = MATMUL(DN(1,:), NORMS_LOC)
      T2 = MATMUL(DN(2,:), NORMS_LOC)
      BBOUT = ZERO
      DO II=1,8
         COL = (II-1)*6
         T0 = NORMS_LOC(II,:)
         CALL TENSOR_PHYS_Q8(DN(1,II)*T1, DN(2,II)*T2, 0.5D0*(DN(1,II)*T2 + DN(2,II)*T1), C, VP)
         BBOUT(1:3,COL+1:COL+3) = VP
         CALL CROSS3(T0, G1, CG1)
         CALL CROSS3(T0, G2, CG2)
         CALL TENSOR_PHYS_Q8(DN(1,II)*CG1, DN(2,II)*CG2, 0.5D0*(DN(1,II)*CG2 + DN(2,II)*CG1), C, VP)
         BBOUT(1:3,COL+4:COL+6) = VP
      ENDDO
      END SUBROUTINE BB_Q8_AT

      SUBROUTINE BS_Q8_AT ( XYZN, NORMS, XI, ETA, BSOUT, JAC, NORMS_EXT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), NORMS(8,3), XI, ETA
      REAL(DOUBLE), INTENT(IN), OPTIONAL :: NORMS_EXT(8,3)
      REAL(DOUBLE), INTENT(OUT) :: BSOUT(2,48), JAC
      REAL(DOUBLE) :: NVAL(8), DN(2,8), E1F(3), E2F(3), E3F(3), G1(3), G2(3), C(4), BSN(2,48), T0(3), T0I(3), C1(3), C2(3), NM
      REAL(DOUBLE) :: NORMS_LOC(8,3)
      INTEGER(LONG) :: II, COL
      NORMS_LOC = NORMS
      IF (PRESENT(NORMS_EXT)) NORMS_LOC = NORMS_EXT
      CALL SHAPE_Q8(XI, ETA, NVAL, DN)
      CALL FIXED_FRAME_Q8(XYZN, E1F, E2F, E3F)
      CALL COV_MAP_Q8(XYZN, XI, ETA, E1F, E2F, G1, G2, C, JAC)
      T0 = MATMUL(NVAL, NORMS_LOC)
      NM = VNORM(T0)
      IF (NM > 1.0D-15) T0 = T0/NM
      BSN = ZERO
      DO II=1,8
         COL = (II-1)*6
         T0I = NORMS_LOC(II,:)
         BSN(1,COL+1:COL+3) = BSN(1,COL+1:COL+3) + DN(1,II)*T0
         BSN(2,COL+1:COL+3) = BSN(2,COL+1:COL+3) + DN(2,II)*T0
         CALL CROSS3(T0I, G1, C1)
         CALL CROSS3(T0I, G2, C2)
         BSN(1,COL+4:COL+6) = BSN(1,COL+4:COL+6) + NVAL(II)*C1
         BSN(2,COL+4:COL+6) = BSN(2,COL+4:COL+6) + NVAL(II)*C2
      ENDDO
      BSOUT(1,:) = C(1)*BSN(1,:) + C(2)*BSN(2,:)
      BSOUT(2,:) = C(3)*BSN(1,:) + C(4)*BSN(2,:)
      END SUBROUTINE BS_Q8_AT

      SUBROUTINE BDRILL_Q8_AT ( XYZN, NORMS, XI, ETA, BDOUT, JAC, NORMS_EXT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), NORMS(8,3), XI, ETA
      REAL(DOUBLE), INTENT(IN), OPTIONAL :: NORMS_EXT(8,3)
      REAL(DOUBLE), INTENT(OUT) :: BDOUT(1,48), JAC
      REAL(DOUBLE) :: NVAL(8), DN(2,8), G1(3), G2(3), E1(3), E2(3), E3(3), A(2,2), AINV(2,2), DLOC(2,8)
      REAL(DOUBLE) :: NORMS_LOC(8,3)
      INTEGER(LONG) :: II, COL
      NORMS_LOC = NORMS
      IF (PRESENT(NORMS_EXT)) NORMS_LOC = NORMS_EXT
      CALL SHAPE_Q8(XI, ETA, NVAL, DN)
      CALL SURFACE_BASIS_Q8(XYZN, XI, ETA, G1, G2, E1, E2, E3, JAC)
      A(1,1) = DOT_PRODUCT(G1,E1)
      A(1,2) = DOT_PRODUCT(G1,E2)
      A(2,1) = DOT_PRODUCT(G2,E1)
      A(2,2) = DOT_PRODUCT(G2,E2)
      CALL INV2(A, AINV)
      DLOC = MATMUL(AINV, DN)
      BDOUT = ZERO
      DO II=1,8
         COL = (II-1)*6
         BDOUT(1,COL+1:COL+3) = 0.5D0*(DLOC(1,II)*E2 - DLOC(2,II)*E1)
         BDOUT(1,COL+4:COL+6) = BDOUT(1,COL+4:COL+6) - NVAL(II)*NORMS_LOC(II,:)
      ENDDO
      END SUBROUTINE BDRILL_Q8_AT

      SUBROUTINE EAS1_SHEAR_Q8_AT ( XYZN, XI, ETA, BSE )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BSE(2)
      REAL(DOUBLE) :: E1F(3), E2F(3), E3F(3), G1(3), G2(3), C(4), JAC, DPHIR, DPHIS
      CALL FIXED_FRAME_Q8(XYZN, E1F, E2F, E3F)
      CALL COV_MAP_Q8(XYZN, XI, ETA, E1F, E2F, G1, G2, C, JAC)
      DPHIR = -TWO*XI*(ONE-ETA*ETA)
      DPHIS = -TWO*ETA*(ONE-XI*XI)
      BSE(1) = DPHIR*C(1) + DPHIS*C(2)
      BSE(2) = DPHIR*C(3) + DPHIS*C(4)
      END SUBROUTINE EAS1_SHEAR_Q8_AT

      SUBROUTINE CROSS3 ( A, B, C )
      REAL(DOUBLE), INTENT(IN)  :: A(3), B(3)
      REAL(DOUBLE), INTENT(OUT) :: C(3)
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)
      END SUBROUTINE CROSS3

      FUNCTION VNORM ( V ) RESULT(NM)
      REAL(DOUBLE), INTENT(IN) :: V(3)
      REAL(DOUBLE) :: NM
      NM = DSQRT(MAX(ZERO, DOT_PRODUCT(V,V)))
      END FUNCTION VNORM

      SUBROUTINE INV2 ( A, AINV )
      REAL(DOUBLE), INTENT(IN)  :: A(2,2)
      REAL(DOUBLE), INTENT(OUT) :: AINV(2,2)
      REAL(DOUBLE) :: DET
      DET = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      IF (DABS(DET) < 1.0D-20) THEN
         AINV = ZERO
         AINV(1,1) = ONE
         AINV(2,2) = ONE
      ELSE
         AINV(1,1) =  A(2,2)/DET
         AINV(1,2) = -A(1,2)/DET
         AINV(2,1) = -A(2,1)/DET
         AINV(2,2) =  A(1,1)/DET
      ENDIF
      END SUBROUTINE INV2

      END SUBROUTINE CQUAD8_SIMOQ8
