! #################################################################################################################################
! CTRIA6 Rezaiee2017 quadratic triangular shell.

      SUBROUTINE CTRIA6_REZAIEE ( OPT, INT_ELEM_ID )

! Ported from:
!   D:\18a\bending_only\Shell\gemini2\shit\validation\q8\Rezaiee2017_Tri6_Final.py
!
! Formulation notes:
!   6-node quadratic triangle with the same director, bending, drilling,
!   mass, pressure and thermal framework as CTRIA6_SIMO1993.
!   Membrane uses the same assumed covariant reconstruction as MITC6.
!   The Rezaiee branch differs mainly in the transverse-shear common point,
!   where rc = sc = 1/sqrt(3) is used in the tied shear interpolation.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_STRESS_POINTS, SOL_NAME
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  ALPVEC, BGRID, DT, EID, ELGP, GRID_SNORM, KE, ME, BE1, BE2, BE3, MASS_PER_UNIT_AREA,    &
                                         NUM_EMG_FATAL_ERRS, PCOMP_PROPS, PPE, PRESS, PTE, RGRID, SHELL_A, SHELL_D, SHELL_T,     &
                                         TREF, XEB
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO
      USE PARAMS, ONLY                :  COUPMASS
      USE OUTA_HERE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CTRIA6_REZAIEE'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      INTEGER(LONG)                   :: I, J, IA, IB, K, L, JSUB
      REAL(DOUBLE)                    :: XYZ(6,3), NORMALS(6,3)
      REAL(DOUBLE)                    :: KOUT(36,36)
      REAL(DOUBLE)                    :: BM(3,36), BB(3,36), BS(2,36), BD(1,36)
      REAL(DOUBLE)                    :: R3(3), S3(3), W3(3), R6(6), S6(6), W6(6)
      REAL(DOUBLE)                    :: R, S, WT, JAC, CDRILL, FAC, SQRT2, SQRT3, R1MITC, R2MITC, RCMITC, RCSHEAR
      REAL(DOUBLE)                    :: M1(6,6), N6(6), DN6(2,6), MASS_ELEM, MASS_NODE
      REAL(DOUBLE)                    :: UNIT_PPE(36), UNIT_PTE(36), DXDR(3), DXDS(3), SURF_VEC(3), TBAR
      REAL(DOUBLE)                    :: CTE(3), THERMAL_RESULTANT(3)

      IF (ELGP /= 6) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9001) SUBR_NAME, EID, ELGP
         WRITE(F06,9001) SUBR_NAME, EID, ELGP
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IF (PCOMP_PROPS == 'Y') THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,*) ' *ERROR: Code not written for composite material with CTRIA6 REZAIEE'
         WRITE(F06,*) ' *ERROR: Code not written for composite material with CTRIA6 REZAIEE'
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL LOAD_BASIC_COORDS_T6 ( XYZ )
      CALL CALC_NODAL_NORMALS_T6 ( XYZ, NORMALS )
      SQRT2 = DSQRT(TWO)
      SQRT3 = DSQRT(3.0D0)
      R1MITC = 0.5D0 - 0.5D0/SQRT3
      R2MITC = 0.5D0 + 0.5D0/SQRT3
      RCMITC = ONE/3.0D0
      RCSHEAR = ONE/SQRT3

      R3 = (/ONE/6.0D0, TWO/3.0D0, ONE/6.0D0/)
      S3 = (/ONE/6.0D0, ONE/6.0D0, TWO/3.0D0/)
      W3 = (/ONE/6.0D0, ONE/6.0D0, ONE/6.0D0/)

      R6 = (/0.445948490144588D0, 0.108103018168070D0, 0.445948490144588D0,                                      &
             0.091576213509771D0, 0.816847572980459D0, 0.091576213509771D0/)
      S6 = (/0.445948490144588D0, 0.445948490144588D0, 0.108103018168070D0,                                      &
             0.091576213509771D0, 0.091576213509771D0, 0.816847572980459D0/)
      W6 = (/0.111690794839050D0, 0.111690794839050D0, 0.111690794839050D0,                                      &
             0.054975871827661D0, 0.054975871827661D0, 0.054975871827661D0/)

      IF (OPT(1) == 'Y') THEN
         M1 = ZERO
         MASS_ELEM = ZERO
         DO I=1,6
            R = R6(I)
            S = S6(I)
            WT = W6(I)
            CALL SHAPE_T6 ( R, S, N6, DN6 )
            CALL BM_REZAIEE_AT ( XYZ, R, S, BM, JAC, R1MITC, R2MITC, RCMITC )
            MASS_ELEM = MASS_ELEM + MASS_PER_UNIT_AREA*WT*JAC
            DO K=1,6
               DO L=1,6
                  M1(K,L) = M1(K,L) + N6(K)*N6(L)*MASS_PER_UNIT_AREA*WT*JAC
               ENDDO
            ENDDO
         ENDDO

         ME = ZERO
         IF ((SOL_NAME(1:5) == 'MODES') .AND. (COUPMASS > 0)) THEN
            DO K=1,6
               DO L=1,6
                  DO IA=1,3
                     ME(6*(K-1)+IA,6*(L-1)+IA) = M1(K,L)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            MASS_NODE = MASS_ELEM/6.0D0
            DO K=1,6
               DO IA=1,3
                  ME(6*(K-1)+IA,6*(K-1)+IA) = MASS_NODE
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      IF (OPT(2) == 'Y') THEN
         UNIT_PTE = ZERO
         DO I=1,6
            R = R6(I)
            S = S6(I)
            WT = W6(I)
            CALL BM_REZAIEE_AT ( XYZ, R, S, BM, JAC, R1MITC, R2MITC, RCMITC )
            CTE(1) = ALPVEC(1,1)
            CTE(2) = ALPVEC(2,1)
            CTE(3) = ALPVEC(4,1)
            THERMAL_RESULTANT = MATMUL(SHELL_A, CTE)
            UNIT_PTE = UNIT_PTE + MATMUL( TRANSPOSE(BM), THERMAL_RESULTANT ) * WT * JAC
         ENDDO
         DO JSUB=1,SIZE(PTE,2)
            TBAR = ZERO
            DO J=1,6
               TBAR = TBAR + DT(J,JSUB)
            ENDDO
            TBAR = TBAR / 6.0D0 - TREF(1)
            PTE(1:36,JSUB) = UNIT_PTE(1:36) * TBAR
         ENDDO
      ENDIF

      IF (OPT(3) == 'Y') THEN
         DO I=1,3
            CALL BM_REZAIEE_AT ( XYZ, R3(I), S3(I), BM, JAC, R1MITC, R2MITC, RCMITC )
            CALL BB_T6_AT ( XYZ, NORMALS, R3(I), S3(I), BB, JAC )
            CALL BS_REZAIEE_AT ( XYZ, NORMALS, R3(I), S3(I), BS, JAC, R1MITC, R2MITC, RCSHEAR, SQRT2 )
            IF (I <= MAX_STRESS_POINTS) THEN
               BE1(1:3,1:36,I) = BM
               BE2(1:3,1:36,I) = BB
               BE3(1:2,1:36,I) = BS
            ENDIF
         ENDDO
      ENDIF

      IF (OPT(4) == 'Y') THEN
         KOUT = ZERO
         CDRILL = 2.0D-2*SHELL_T(1,1)/(5.0D0/6.0D0)

         DO I=1,3
            R = R3(I)
            S = S3(I)
            WT = W3(I)
            CALL BM_REZAIEE_AT ( XYZ, R, S, BM, JAC, R1MITC, R2MITC, RCMITC )
            CALL BB_T6_AT ( XYZ, NORMALS, R, S, BB, JAC )
            CALL BS_REZAIEE_AT ( XYZ, NORMALS, R, S, BS, JAC, R1MITC, R2MITC, RCSHEAR, SQRT2 )
            KOUT = KOUT + WT*JAC*MATMUL(TRANSPOSE(BM), MATMUL(SHELL_A, BM))
            KOUT = KOUT + WT*JAC*MATMUL(TRANSPOSE(BB), MATMUL(SHELL_D, BB))
            KOUT = KOUT + WT*JAC*MATMUL(TRANSPOSE(BS), MATMUL(SHELL_T, BS))
         ENDDO

         DO I=1,6
            R = R6(I)
            S = S6(I)
            WT = W6(I)
            CALL BDRILL_T6_AT ( XYZ, NORMALS, R, S, BD, JAC )
            KOUT = KOUT + WT*JAC*CDRILL*MATMUL(TRANSPOSE(BD), BD)
         ENDDO

         DO IA=1,36
            DO IB=1,36
               FAC = 0.5D0*(KOUT(IA,IB) + KOUT(IB,IA))
               KE(IA,IB) = FAC
            ENDDO
         ENDDO
      ENDIF

      IF (OPT(5) == 'Y') THEN
         UNIT_PPE = ZERO
         DO I=1,6
            R = R6(I)
            S = S6(I)
            WT = W6(I)
            CALL SHAPE_T6 ( R, S, N6, DN6 )
            DXDR = ZERO
            DXDS = ZERO
            DO K=1,6
               DXDR(:) = DXDR(:) + DN6(1,K)*XYZ(K,:)
               DXDS(:) = DXDS(:) + DN6(2,K)*XYZ(K,:)
            ENDDO
            CALL CROSS3 ( DXDR, DXDS, SURF_VEC )
            DO K=1,6
               UNIT_PPE(6*(K-1)+1) = UNIT_PPE(6*(K-1)+1) + N6(K)*SURF_VEC(1)*WT
               UNIT_PPE(6*(K-1)+2) = UNIT_PPE(6*(K-1)+2) + N6(K)*SURF_VEC(2)*WT
               UNIT_PPE(6*(K-1)+3) = UNIT_PPE(6*(K-1)+3) + N6(K)*SURF_VEC(3)*WT
            ENDDO
         ENDDO
         DO J=1,SIZE(PPE,2)
            PPE(1:36,J) = PPE(1:36,J) + UNIT_PPE(1:36)*PRESS(3,J)
         ENDDO
      ENDIF

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN
         WRITE(ERR,*) ' *ERROR: Code not written for CTRIA6 REZAIEE differential stiffness matrix'
         WRITE(F06,*) ' *ERROR: Code not written for CTRIA6 REZAIEE differential stiffness matrix'
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      RETURN

 9001 FORMAT(' *ERROR: ',A,' expects ELGP=6 for element ',I8,' but got ',I8)

      CONTAINS

      SUBROUTINE LOAD_BASIC_COORDS_T6 ( XYZOUT )
      REAL(DOUBLE), INTENT(OUT) :: XYZOUT(6,3)
      INTEGER(LONG) :: II, JJ
      DO II=1,6
         IF ((II <= SIZE(BGRID)) .AND. (BGRID(II) > 0)) THEN
            DO JJ=1,3
               XYZOUT(II,JJ) = RGRID(BGRID(II),JJ)
            ENDDO
         ELSE
            DO JJ=1,3
               XYZOUT(II,JJ) = XEB(II,JJ)
            ENDDO
         ENDIF
      ENDDO
      END SUBROUTINE LOAD_BASIC_COORDS_T6

      SUBROUTINE SHAPE_T6 ( R, S, NVAL, DN )
      REAL(DOUBLE), INTENT(IN)  :: R, S
      REAL(DOUBLE), INTENT(OUT) :: NVAL(6), DN(2,6)
      REAL(DOUBLE) :: L1, L2, L3
      L1 = ONE - R - S
      L2 = R
      L3 = S
      NVAL(1) = L1*(TWO*L1 - ONE)
      NVAL(2) = L2*(TWO*L2 - ONE)
      NVAL(3) = L3*(TWO*L3 - ONE)
      NVAL(4) = 4.0D0*L1*L2
      NVAL(5) = 4.0D0*L2*L3
      NVAL(6) = 4.0D0*L3*L1

      DN(1,1) = 4.0D0*R + 4.0D0*S - 3.0D0
      DN(1,2) = 4.0D0*R - ONE
      DN(1,3) = ZERO
      DN(1,4) = 4.0D0 - 8.0D0*R - 4.0D0*S
      DN(1,5) = 4.0D0*S
      DN(1,6) = -4.0D0*S

      DN(2,1) = 4.0D0*R + 4.0D0*S - 3.0D0
      DN(2,2) = ZERO
      DN(2,3) = 4.0D0*S - ONE
      DN(2,4) = -4.0D0*R
      DN(2,5) = 4.0D0*R
      DN(2,6) = 4.0D0 - 4.0D0*R - 8.0D0*S
      END SUBROUTINE SHAPE_T6

      SUBROUTINE CALC_NODAL_NORMALS_T6 ( XYZN, NORMS )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3)
      REAL(DOUBLE), INTENT(OUT) :: NORMS(6,3)
      REAL(DOUBLE) :: RS(6,2), NVAL(6), DN(2,6), G1(3), G2(3), N(3), NM, SN(3), SDOT
      INTEGER(LONG) :: II, BIDX
      RS(1,:) = (/ZERO, ZERO/)
      RS(2,:) = (/ONE, ZERO/)
      RS(3,:) = (/ZERO, ONE/)
      RS(4,:) = (/0.5D0, ZERO/)
      RS(5,:) = (/0.5D0, 0.5D0/)
      RS(6,:) = (/ZERO, 0.5D0/)
      DO II=1,6
         CALL SHAPE_T6(RS(II,1), RS(II,2), NVAL, DN)
         G1 = MATMUL(DN(1,:), XYZN)
         G2 = MATMUL(DN(2,:), XYZN)
         CALL CROSS3(G1, G2, N)
         NM = VNORM(N)
         IF (NM <= 1.0D-12) THEN
            CALL SHAPE_T6(ONE/3.0D0, ONE/3.0D0, NVAL, DN)
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
                  SDOT = DOT_PRODUCT(SN, NORMS(II,:))
                  IF (SDOT < ZERO) THEN
                     SN = -SN
                     SDOT = -SDOT
                  ENDIF
                  IF (SDOT >= 0.85D0) THEN
                     NORMS(II,:) = SN
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      END SUBROUTINE CALC_NODAL_NORMALS_T6

      SUBROUTINE SURFACE_BASIS_T6 ( XYZN, R, S, G1, G2, E1, E2, E3, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3), R, S
      REAL(DOUBLE), INTENT(OUT) :: G1(3), G2(3), E1(3), E2(3), E3(3), JAC
      REAL(DOUBLE) :: NVAL(6), DN(2,6), G3(3), TMP(3), NM
      CALL SHAPE_T6(R, S, NVAL, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      CALL CROSS3(G1, G2, G3)
      JAC = VNORM(G3)
      IF (JAC > 1.0D-15) THEN
         E3 = G3/JAC
      ELSE
         E3 = (/ZERO, ZERO, ONE/)
      ENDIF
      NM = VNORM(G1)
      IF (NM > 1.0D-15) THEN
         E1 = G1/NM
      ELSE
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
      ENDIF
      CALL CROSS3(E3, E1, E2)
      NM = VNORM(E2)
      IF (NM > 1.0D-15) E2 = E2/NM
      END SUBROUTINE SURFACE_BASIS_T6

      SUBROUTINE FIXED_FRAME_T6 ( XYZN, E1F, E2F, E3F )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3)
      REAL(DOUBLE), INTENT(OUT) :: E1F(3), E2F(3), E3F(3)
      REAL(DOUBLE) :: G1(3), G2(3), JAC, ZSGN
      CALL SURFACE_BASIS_T6(XYZN, ONE/3.0D0, ONE/3.0D0, G1, G2, E1F, E2F, E3F, JAC)
      IF ((DABS(E3F(1)) + DABS(E3F(2))) <= 1.0D-12) THEN
         ZSGN = ONE
         IF (E3F(3) < ZERO) ZSGN = -ONE
         E1F = (/ONE, ZERO, ZERO/)
         E2F = (/ZERO, ZSGN, ZERO/)
         E3F = (/ZERO, ZERO, ZSGN/)
      ENDIF
      END SUBROUTINE FIXED_FRAME_T6

      SUBROUTINE COV_MAP_T6 ( XYZN, R, S, E1F, E2F, G1, G2, C, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3), R, S, E1F(3), E2F(3)
      REAL(DOUBLE), INTENT(OUT) :: G1(3), G2(3), C(4), JAC
      REAL(DOUBLE) :: E1(3), E2(3), E3(3), A11, A22, A12, DET, AI11, AI22, AI12, GC1(3), GC2(3)
      CALL SURFACE_BASIS_T6(XYZN, R, S, G1, G2, E1, E2, E3, JAC)
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
      C(1) = DOT_PRODUCT(E1, GC1)
      C(2) = DOT_PRODUCT(E1, GC2)
      C(3) = DOT_PRODUCT(E2, GC1)
      C(4) = DOT_PRODUCT(E2, GC2)
      END SUBROUTINE COV_MAP_T6

      SUBROUTINE TENSOR_PHYS_T6 ( V11, V22, V12, C, VOUT )
      REAL(DOUBLE), INTENT(IN)  :: V11(3), V22(3), V12(3), C(4)
      REAL(DOUBLE), INTENT(OUT) :: VOUT(3,3)
      VOUT(1,:) = C(1)*C(1)*V11 + C(2)*C(2)*V22 + TWO*C(1)*C(2)*V12
      VOUT(2,:) = C(3)*C(3)*V11 + C(4)*C(4)*V22 + TWO*C(3)*C(4)*V12
      VOUT(3,:) = TWO*(C(1)*C(3)*V11 + C(2)*C(4)*V22 + (C(1)*C(4)+C(2)*C(3))*V12)
      END SUBROUTINE TENSOR_PHYS_T6

      SUBROUTINE BM_REZAIEE_AT ( XYZN, R, S, BMOUT, JAC, R1MITC, R2MITC, RCMITC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3), R, S, R1MITC, R2MITC, RCMITC
      REAL(DOUBLE), INTENT(OUT) :: BMOUT(3,36), JAC
      REAL(DOUBLE) :: RR11(36), SS11(36), RS11(36), RR21(36), SS21(36), RS21(36), RR12(36), SS12(36), RS12(36)
      REAL(DOUBLE) :: RRC(36), SSC(36), RSC(36), QQ21(36), QQ12(36), QQC(36), A1(36), B1(36), C1V(36)
      REAL(DOUBLE) :: A2(36), B2(36), C2V(36), A3(36), B3(36), C3V(36), BERR(36), BESS(36), BEQQ(36), BERS(36)
      REAL(DOUBLE) :: E1F(3), E2F(3), E3F(3), G1(3), G2(3), C(4), VP(3,3)
      REAL(DOUBLE) :: VALS3(3,36), PTS3(3,2)
      CALL FIXED_FRAME_T6(XYZN, E1F, E2F, E3F)
      CALL NAT_MEM_ROWS_T6(XYZN, R1MITC, R1MITC, RR11, SS11, RS11)
      CALL NAT_MEM_ROWS_T6(XYZN, R2MITC, R1MITC, RR21, SS21, RS21)
      CALL NAT_MEM_ROWS_T6(XYZN, R1MITC, R2MITC, RR12, SS12, RS12)
      CALL NAT_MEM_ROWS_T6(XYZN, RCMITC, RCMITC, RRC, SSC, RSC)
      QQ21 = 0.5D0*(RR21 + SS21) - RS21
      QQ12 = 0.5D0*(RR12 + SS12) - RS12
      QQC  = 0.5D0*(RRC  + SSC ) - RSC

      VALS3(1,:) = RR11
      VALS3(2,:) = RR21
      VALS3(3,:) = RRC
      PTS3(1,:) = (/R1MITC, R1MITC/)
      PTS3(2,:) = (/R2MITC, R1MITC/)
      PTS3(3,:) = (/RCMITC, RCMITC/)
      CALL FIT_AFFINE36_T6(VALS3, PTS3, 1, A1, B1, C1V)

      VALS3(1,:) = SS11
      VALS3(2,:) = SS12
      VALS3(3,:) = SSC
      PTS3(1,:) = (/R1MITC, R1MITC/)
      PTS3(2,:) = (/R1MITC, R2MITC/)
      PTS3(3,:) = (/RCMITC, RCMITC/)
      CALL FIT_AFFINE36_T6(VALS3, PTS3, 1, A2, B2, C2V)

      VALS3(1,:) = QQ21
      VALS3(2,:) = QQ12
      VALS3(3,:) = QQC
      PTS3(1,:) = (/R2MITC, R1MITC/)
      PTS3(2,:) = (/R1MITC, R2MITC/)
      PTS3(3,:) = (/RCMITC, RCMITC/)
      CALL FIT_AFFINE36_T6(VALS3, PTS3, 2, A3, B3, C3V)

      BERR = A1 + B1*R + C1V*S
      BESS = A2 + B2*R + C2V*S
      BEQQ = A3 + B3*R + C3V*(ONE - R - S)
      BERS = 0.5D0*(BERR + BESS) - BEQQ

      CALL COV_MAP_T6(XYZN, R, S, E1F, E2F, G1, G2, C, JAC)
      CALL TENSOR_PHYS_T6(BERR(1:3), BESS(1:3), BERS(1:3), C, VP)
      BMOUT = ZERO
      BMOUT(:,1:3) = VP
      CALL TENSOR_PHYS_T6(BERR(4:6), BESS(4:6), BERS(4:6), C, VP)
      BMOUT(:,4:6) = VP
      CALL TENSOR_PHYS_T6(BERR(7:9), BESS(7:9), BERS(7:9), C, VP)
      BMOUT(:,7:9) = VP
      CALL TENSOR_PHYS_T6(BERR(10:12), BESS(10:12), BERS(10:12), C, VP)
      BMOUT(:,10:12) = VP
      CALL TENSOR_PHYS_T6(BERR(13:15), BESS(13:15), BERS(13:15), C, VP)
      BMOUT(:,13:15) = VP
      CALL TENSOR_PHYS_T6(BERR(16:18), BESS(16:18), BERS(16:18), C, VP)
      BMOUT(:,16:18) = VP
      CALL TENSOR_PHYS_T6(BERR(19:21), BESS(19:21), BERS(19:21), C, VP)
      BMOUT(:,19:21) = VP
      CALL TENSOR_PHYS_T6(BERR(22:24), BESS(22:24), BERS(22:24), C, VP)
      BMOUT(:,22:24) = VP
      CALL TENSOR_PHYS_T6(BERR(25:27), BESS(25:27), BERS(25:27), C, VP)
      BMOUT(:,25:27) = VP
      CALL TENSOR_PHYS_T6(BERR(28:30), BESS(28:30), BERS(28:30), C, VP)
      BMOUT(:,28:30) = VP
      CALL TENSOR_PHYS_T6(BERR(31:33), BESS(31:33), BERS(31:33), C, VP)
      BMOUT(:,31:33) = VP
      CALL TENSOR_PHYS_T6(BERR(34:36), BESS(34:36), BERS(34:36), C, VP)
      BMOUT(:,34:36) = VP
      END SUBROUTINE BM_REZAIEE_AT

      SUBROUTINE BB_T6_AT ( XYZN, NORMS, R, S, BBOUT, JAC, NORMS_EXT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3), NORMS(6,3), R, S
      REAL(DOUBLE), INTENT(IN), OPTIONAL :: NORMS_EXT(6,3)
      REAL(DOUBLE), INTENT(OUT) :: BBOUT(3,36), JAC
      REAL(DOUBLE) :: NVAL(6), DN(2,6), E1F(3), E2F(3), E3F(3), G1(3), G2(3), C(4), VP(3,3)
      REAL(DOUBLE) :: T1(3), T2(3), T0(3), CG1(3), CG2(3)
      REAL(DOUBLE) :: NORMS_LOC(6,3)
      INTEGER(LONG) :: II, COL
      NORMS_LOC = NORMS
      IF (PRESENT(NORMS_EXT)) NORMS_LOC = NORMS_EXT
      CALL SHAPE_T6(R, S, NVAL, DN)
      CALL FIXED_FRAME_T6(XYZN, E1F, E2F, E3F)
      CALL COV_MAP_T6(XYZN, R, S, E1F, E2F, G1, G2, C, JAC)
      T1 = MATMUL(DN(1,:), NORMS_LOC)
      T2 = MATMUL(DN(2,:), NORMS_LOC)
      BBOUT = ZERO
      DO II=1,6
         COL = (II-1)*6
         T0 = NORMS_LOC(II,:)
         CALL TENSOR_PHYS_T6(DN(1,II)*T1, DN(2,II)*T2, 0.5D0*(DN(1,II)*T2 + DN(2,II)*T1), C, VP)
         BBOUT(1:3,COL+1:COL+3) = VP
         CALL CROSS3(T0, G1, CG1)
         CALL CROSS3(T0, G2, CG2)
         CALL TENSOR_PHYS_T6(DN(1,II)*CG1, DN(2,II)*CG2, 0.5D0*(DN(1,II)*CG2 + DN(2,II)*CG1), C, VP)
         BBOUT(1:3,COL+4:COL+6) = VP
      ENDDO
      END SUBROUTINE BB_T6_AT

      SUBROUTINE BS_REZAIEE_AT ( XYZN, NORMS, R, S, BSOUT, JAC, R1MITC, R2MITC, RCSHEAR, SQRT2 )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3), NORMS(6,3), R, S, R1MITC, R2MITC, RCSHEAR, SQRT2
      REAL(DOUBLE), INTENT(OUT) :: BSOUT(2,36), JAC
      REAL(DOUBLE) :: RHS(8,36), MAT8(8,8), MAT8INV(8,8), COEFF(8,36), BRT(36), BST(36)
      REAL(DOUBLE) :: A1(36), B1(36), C1V(36), A2(36), B2(36), C2V(36), D(36), E(36)
      REAL(DOUBLE) :: E1F(3), E2F(3), E3F(3), G1(3), G2(3), C(4), BRTPHYS(36), BSTPHYS(36)
      INTEGER(LONG) :: INFO
      CALL FIXED_FRAME_T6(XYZN, E1F, E2F, E3F)

      CALL NAT_SHEAR_ROWS_T6(XYZN, NORMS, R1MITC, ZERO,   BRT, BST, JAC)
      RHS(1,:) = BRT
      CALL NAT_SHEAR_ROWS_T6(XYZN, NORMS, R2MITC, ZERO,   BRT, BST, JAC)
      RHS(2,:) = BRT
      CALL NAT_SHEAR_ROWS_T6(XYZN, NORMS, R1MITC, R1MITC, BRT, BST, JAC)
      RHS(3,:) = BST
      CALL NAT_SHEAR_ROWS_T6(XYZN, NORMS, R1MITC, R2MITC, BRT, BST, JAC)
      RHS(4,:) = BST
      CALL NAT_SHEAR_ROWS_T6(XYZN, NORMS, RCSHEAR, RCSHEAR, BRT, BST, JAC)
      RHS(5,:) = BRT
      RHS(6,:) = BST
      CALL NAT_SHEAR_ROWS_T6(XYZN, NORMS, R2MITC, R1MITC, BRT, BST, JAC)
      RHS(7,:) = (BST - BRT)/SQRT2
      CALL NAT_SHEAR_ROWS_T6(XYZN, NORMS, R1MITC, R2MITC, BRT, BST, JAC)
      RHS(8,:) = (BST - BRT)/SQRT2

      MAT8 = ZERO
      MAT8(1,:) = (/ONE, R1MITC, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/)
      MAT8(2,:) = (/ONE, R2MITC, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/)
      MAT8(3,:) = (/ZERO, ZERO, ZERO, ONE, R1MITC, R1MITC, -(R1MITC*R1MITC), -(R1MITC*R1MITC)/)
      MAT8(4,:) = (/ZERO, ZERO, ZERO, ONE, R1MITC, R2MITC, -(R1MITC*R1MITC), -(R1MITC*R2MITC)/)
      MAT8(5,:) = (/ONE, RCSHEAR, RCSHEAR, ZERO, ZERO, ZERO, RCSHEAR*RCSHEAR, RCSHEAR*RCSHEAR/)
      MAT8(6,:) = (/ZERO, ZERO, ZERO, ONE, RCSHEAR, RCSHEAR, -(RCSHEAR*RCSHEAR), -(RCSHEAR*RCSHEAR)/)
      MAT8(7,:) = (/ -ONE/SQRT2, -R2MITC/SQRT2, -R1MITC/SQRT2, ONE/SQRT2, R2MITC/SQRT2, R1MITC/SQRT2, &
                      (-(R2MITC*R2MITC) - R2MITC*R1MITC)/SQRT2, (-(R2MITC*R1MITC) - R1MITC*R1MITC)/SQRT2 /)
      MAT8(8,:) = (/ -ONE/SQRT2, -R1MITC/SQRT2, -R2MITC/SQRT2, ONE/SQRT2, R1MITC/SQRT2, R2MITC/SQRT2, &
                      (-(R1MITC*R1MITC) - R1MITC*R2MITC)/SQRT2, (-(R1MITC*R2MITC) - R2MITC*R2MITC)/SQRT2 /)

      CALL INV_SMALL_T6(MAT8, 8, MAT8INV, INFO)
      COEFF = MATMUL(MAT8INV, RHS)
      A1  = COEFF(1,:)
      B1  = COEFF(2,:)
      C1V = COEFF(3,:)
      A2  = COEFF(4,:)
      B2  = COEFF(5,:)
      C2V = COEFF(6,:)
      D   = COEFF(7,:)
      E   = COEFF(8,:)

      BRT = A1 + B1*R + C1V*S + D*R*S + E*S*S
      BST = A2 + B2*R + C2V*S - D*R*R - E*R*S
      CALL COV_MAP_T6(XYZN, R, S, E1F, E2F, G1, G2, C, JAC)
      BRTPHYS = C(1)*BRT + C(2)*BST
      BSTPHYS = C(3)*BRT + C(4)*BST
      BSOUT(1,:) = BRTPHYS
      BSOUT(2,:) = BSTPHYS
      END SUBROUTINE BS_REZAIEE_AT

      SUBROUTINE NAT_MEM_ROWS_T6 ( XYZN, R, S, BRR, BSS, BRS )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3), R, S
      REAL(DOUBLE), INTENT(OUT) :: BRR(36), BSS(36), BRS(36)
      REAL(DOUBLE) :: NVAL(6), DN(2,6), G1(3), G2(3)
      INTEGER(LONG) :: II, COL
      CALL SHAPE_T6(R, S, NVAL, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      BRR = ZERO
      BSS = ZERO
      BRS = ZERO
      DO II=1,6
         COL = (II-1)*6
         BRR(COL+1:COL+3) = DN(1,II)*G1
         BSS(COL+1:COL+3) = DN(2,II)*G2
         BRS(COL+1:COL+3) = 0.5D0*(DN(1,II)*G2 + DN(2,II)*G1)
      ENDDO
      END SUBROUTINE NAT_MEM_ROWS_T6

      SUBROUTINE NAT_SHEAR_ROWS_T6 ( XYZN, NORMS, R, S, BRT, BST, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3), NORMS(6,3), R, S
      REAL(DOUBLE), INTENT(OUT) :: BRT(36), BST(36), JAC
      REAL(DOUBLE) :: NVAL(6), DN(2,6), G1(3), G2(3), E1(3), E2(3), E3(3), T0(3), T0I(3), C1(3), C2(3), NM
      INTEGER(LONG) :: II, COL
      CALL SHAPE_T6(R, S, NVAL, DN)
      CALL SURFACE_BASIS_T6(XYZN, R, S, G1, G2, E1, E2, E3, JAC)
      T0 = MATMUL(NVAL, NORMS)
      NM = VNORM(T0)
      IF (NM > 1.0D-15) T0 = T0/NM
      BRT = ZERO
      BST = ZERO
      DO II=1,6
         COL = (II-1)*6
         T0I = NORMS(II,:)
         BRT(COL+1:COL+3) = BRT(COL+1:COL+3) + DN(1,II)*T0
         BST(COL+1:COL+3) = BST(COL+1:COL+3) + DN(2,II)*T0
         CALL CROSS3(T0I, G1, C1)
         CALL CROSS3(T0I, G2, C2)
         BRT(COL+4:COL+6) = BRT(COL+4:COL+6) + NVAL(II)*C1
         BST(COL+4:COL+6) = BST(COL+4:COL+6) + NVAL(II)*C2
      ENDDO
      END SUBROUTINE NAT_SHEAR_ROWS_T6

      SUBROUTINE FIT_AFFINE36_T6 ( VALS, PTS, BASIS_MODE, AVEC, BVEC, CVEC )
      REAL(DOUBLE), INTENT(IN)  :: VALS(3,36), PTS(3,2)
      INTEGER(LONG), INTENT(IN) :: BASIS_MODE
      REAL(DOUBLE), INTENT(OUT) :: AVEC(36), BVEC(36), CVEC(36)
      REAL(DOUBLE) :: A3(8,8), AINV(8,8)
      INTEGER(LONG) :: IROW, INFO
      A3 = ZERO
      DO IROW=1,3
         A3(IROW,1) = ONE
         A3(IROW,2) = PTS(IROW,1)
         IF (BASIS_MODE == 1) THEN
            A3(IROW,3) = PTS(IROW,2)
         ELSE
            A3(IROW,3) = ONE - PTS(IROW,1) - PTS(IROW,2)
         ENDIF
      ENDDO
      CALL INV_SMALL_T6(A3, 3, AINV, INFO)
      AVEC = AINV(1,1)*VALS(1,:) + AINV(1,2)*VALS(2,:) + AINV(1,3)*VALS(3,:)
      BVEC = AINV(2,1)*VALS(1,:) + AINV(2,2)*VALS(2,:) + AINV(2,3)*VALS(3,:)
      CVEC = AINV(3,1)*VALS(1,:) + AINV(3,2)*VALS(2,:) + AINV(3,3)*VALS(3,:)
      END SUBROUTINE FIT_AFFINE36_T6

      SUBROUTINE INV_SMALL_T6 ( AIN, N, AINV, INFO )
      REAL(DOUBLE), INTENT(IN)  :: AIN(8,8)
      INTEGER(LONG), INTENT(IN) :: N
      REAL(DOUBLE), INTENT(OUT) :: AINV(8,8)
      INTEGER(LONG), INTENT(OUT):: INFO
      REAL(DOUBLE) :: AW(8,8), IW(8,8), PIV, FACT, ABSMAX, TMPROW(8)
      INTEGER(LONG) :: I, J, K, IPIV
      AW = ZERO
      IW = ZERO
      AW(1:N,1:N) = AIN(1:N,1:N)
      DO I=1,N
         IW(I,I) = ONE
      ENDDO
      INFO = 0
      DO I=1,N
         ABSMAX = ZERO
         IPIV = I
         DO K=I,N
            IF (DABS(AW(K,I)) > ABSMAX) THEN
               ABSMAX = DABS(AW(K,I))
               IPIV = K
            ENDIF
         ENDDO
         IF (ABSMAX <= 1.0D-20) THEN
            INFO = 1
            AINV = ZERO
            DO K=1,N
               AINV(K,K) = ONE
            ENDDO
            RETURN
         ENDIF
         IF (IPIV /= I) THEN
            TMPROW(1:N) = AW(I,1:N)
            AW(I,1:N) = AW(IPIV,1:N)
            AW(IPIV,1:N) = TMPROW(1:N)
            TMPROW(1:N) = IW(I,1:N)
            IW(I,1:N) = IW(IPIV,1:N)
            IW(IPIV,1:N) = TMPROW(1:N)
         ENDIF
         PIV = AW(I,I)
         AW(I,1:N) = AW(I,1:N)/PIV
         IW(I,1:N) = IW(I,1:N)/PIV
         DO J=1,N
            IF (J /= I) THEN
               FACT = AW(J,I)
               AW(J,1:N) = AW(J,1:N) - FACT*AW(I,1:N)
               IW(J,1:N) = IW(J,1:N) - FACT*IW(I,1:N)
            ENDIF
         ENDDO
      ENDDO
      AINV = ZERO
      AINV(1:N,1:N) = IW(1:N,1:N)
      END SUBROUTINE INV_SMALL_T6

      SUBROUTINE BDRILL_T6_AT ( XYZN, NORMS, R, S, BDOUT, JAC, NORMS_EXT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3), NORMS(6,3), R, S
      REAL(DOUBLE), INTENT(IN), OPTIONAL :: NORMS_EXT(6,3)
      REAL(DOUBLE), INTENT(OUT) :: BDOUT(1,36), JAC
      REAL(DOUBLE) :: NVAL(6), DN(2,6), G1(3), G2(3), E1(3), E2(3), E3(3), A(2,2), AINV(2,2), DLOC(2,6)
      REAL(DOUBLE) :: NORMS_LOC(6,3)
      INTEGER(LONG) :: II, COL
      NORMS_LOC = NORMS
      IF (PRESENT(NORMS_EXT)) NORMS_LOC = NORMS_EXT
      CALL SHAPE_T6(R, S, NVAL, DN)
      CALL SURFACE_BASIS_T6(XYZN, R, S, G1, G2, E1, E2, E3, JAC)
      A(1,1) = DOT_PRODUCT(G1,E1)
      A(1,2) = DOT_PRODUCT(G1,E2)
      A(2,1) = DOT_PRODUCT(G2,E1)
      A(2,2) = DOT_PRODUCT(G2,E2)
      CALL INV2(A, AINV)
      DLOC = MATMUL(AINV, DN)
      BDOUT = ZERO
      DO II=1,6
         COL = (II-1)*6
         BDOUT(1,COL+1:COL+3) = 0.5D0*(DLOC(1,II)*E2 - DLOC(2,II)*E1)
         BDOUT(1,COL+4:COL+6) = BDOUT(1,COL+4:COL+6) - NVAL(II)*NORMS_LOC(II,:)
      ENDDO
      END SUBROUTINE BDRILL_T6_AT

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

      END SUBROUTINE CTRIA6_REZAIEE
