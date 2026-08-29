! #################################################################################################################################
! CQUAD8 ANS8BDG6 shell for PARAM,QUAD8TYP,ANS8BDG6.

      SUBROUTINE CQUAD8_ANS8BDG6 ( OPT, INT_ELEM_ID )

! Ported from:
!   D:\18a\python\quadratic\ANS8_BDG6_ShellElement.py

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_STRESS_POINTS, SOL_NAME
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  ALPVEC, DT, EID, ELGP, KE, KED, ME, BE1, BE2, BE3, EPROP, MASS_PER_UNIT_AREA, PPE,   &
                                         PRESS, PTE, SHELL_A, SHELL_D, SHELL_T, TREF, UEL, XEB, NUM_EMG_FATAL_ERRS,           &
                                         PCOMP_PROPS
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  COUPMASS
      USE ELMDIS_Interface
      USE OUTA_HERE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CQUAD8_ANS8BDG6'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      INTEGER(LONG), PARAMETER        :: NNODE = 8
      INTEGER(LONG), PARAMETER        :: NDOF  = 48
      REAL(DOUBLE), PARAMETER         :: KT_DRILL = 1.0D-1
      INTEGER(LONG)                   :: I, J, K, L, IA, JSUB, GP
      REAL(DOUBLE)                    :: XYZ(8,3), XY8(8,2), E1F(3), E2F(3), E3F(3)
      REAL(DOUBLE)                    :: GP3(3), W3(3), R, S, WT, DETJ, THICK, TBAR, GVAL
      REAL(DOUBLE)                    :: BM(3,NDOF), BB(3,NDOF), BS(2,NDOF), BD(1,NDOF)
      REAL(DOUBLE)                    :: M1(8,8), N8(8), DNDX(8), DNDY(8), MASS_ELEM, MASS_NODE
      REAL(DOUBLE)                    :: UNIT_PPE(NDOF), UNIT_PTE(NDOF), DXDR(3), DXDS(3), SURF_VEC(3)
      REAL(DOUBLE)                    :: CTE(3), THERMAL_RESULTANT(3), CDRILL
      REAL(DOUBLE)                    :: SIG0(2,2), KG8(8,8), STRAIN0(3), N0V(3)
      INTEGER(LONG)                   :: KI, KJ

      IF (ELGP /= NNODE) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9001) SUBR_NAME, EID, ELGP
         WRITE(F06,9001) SUBR_NAME, EID, ELGP
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IF (PCOMP_PROPS == 'Y') THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,*) ' *ERROR: Code not written for composite material with PARAM,QUAD8TYP,ANS8BDG6'
         WRITE(F06,*) ' *ERROR: Code not written for composite material with PARAM,QUAD8TYP,ANS8BDG6'
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL LOAD_BASIC_COORDS_Q8(XYZ)
      CALL FIXED_FRAME_Q8(XYZ, E1F, E2F, E3F)
      CALL BUILD_LOCAL_XY_Q8(XYZ, E1F, E2F, XY8)

      GP3 = (/-DSQRT(3.0D0/5.0D0), ZERO, DSQRT(3.0D0/5.0D0)/)
      W3  = (/5.0D0/9.0D0, 8.0D0/9.0D0, 5.0D0/9.0D0/)
      THICK = EPROP(1)
      GVAL = SHELL_T(1,1)
      CDRILL = KT_DRILL*GVAL

      IF (OPT(1) == 'Y') THEN
         M1 = ZERO
         MASS_ELEM = ZERO
         DO I=1,3
            DO J=1,3
               R = GP3(I)
               S = GP3(J)
               WT = W3(I)*W3(J)
               CALL Q8_DXY(XY8, R, S, N8, DNDX, DNDY, DETJ)
               MASS_ELEM = MASS_ELEM + MASS_PER_UNIT_AREA*WT*DABS(DETJ)
               DO K=1,8
                  DO L=1,8
                     M1(K,L) = M1(K,L) + N8(K)*N8(L)*MASS_PER_UNIT_AREA*WT*DABS(DETJ)
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
               CALL BM_ANS_Q8_AT(XY8, R, S, BM, DETJ)
               CTE(1) = ALPVEC(1,1)
               CTE(2) = ALPVEC(2,1)
               CTE(3) = ALPVEC(4,1)
               THERMAL_RESULTANT = MATMUL(SHELL_A, CTE)
               UNIT_PTE = UNIT_PTE + MATMUL(TRANSPOSE(BM), THERMAL_RESULTANT)*WT*DABS(DETJ)
            ENDDO
         ENDDO
         DO JSUB=1,SIZE(PTE,2)
            TBAR = ZERO
            DO J=1,8
               TBAR = TBAR + DT(J,JSUB)
            ENDDO
            TBAR = TBAR/8.0D0 - TREF(1)
            PTE(1:NDOF,JSUB) = UNIT_PTE(1:NDOF)*TBAR
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
               CALL BM_ANS_Q8_AT(XY8, R, S, BM, DETJ)
               CALL BB_Q8_AT(XY8, R, S, BB, DETJ)
               CALL BS_ANS_Q8_AT(XY8, R, S, BS, DETJ)
               IF (GP <= MAX_STRESS_POINTS) THEN
                  BE1(1:3,1:NDOF,GP) = BM
                  BE2(1:3,1:NDOF,GP) = BB
                  BE3(1:2,1:NDOF,GP) = BS
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (OPT(4) == 'Y') THEN
         KE = ZERO
         DO I=1,3
            DO J=1,3
               R = GP3(I)
               S = GP3(J)
               WT = W3(I)*W3(J)
               CALL BM_ANS_Q8_AT(XY8, R, S, BM, DETJ)
               CALL BB_Q8_AT(XY8, R, S, BB, DETJ)
               CALL BS_ANS_Q8_AT(XY8, R, S, BS, DETJ)
               CALL BDRILL_Q8_AT(XY8, R, S, BD, DETJ)
               KE = KE + WT*DABS(DETJ)*MATMUL(TRANSPOSE(BM), MATMUL(SHELL_A, BM))
               KE = KE + WT*DABS(DETJ)*MATMUL(TRANSPOSE(BB), MATMUL(SHELL_D, BB))
               KE = KE + WT*DABS(DETJ)*MATMUL(TRANSPOSE(BS), MATMUL(SHELL_T, BS))
               KE = KE + WT*DABS(DETJ)*THICK*CDRILL*MATMUL(TRANSPOSE(BD), BD)
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
               CALL SHAPE_Q8_STD(R, S, N8, DXDR, DXDS, XYZ, DETJ)
               CALL CROSS3(DXDR, DXDS, SURF_VEC)
               DO K=1,8
                  UNIT_PPE(6*(K-1)+1) = UNIT_PPE(6*(K-1)+1) + N8(K)*SURF_VEC(1)*WT
                  UNIT_PPE(6*(K-1)+2) = UNIT_PPE(6*(K-1)+2) + N8(K)*SURF_VEC(2)*WT
                  UNIT_PPE(6*(K-1)+3) = UNIT_PPE(6*(K-1)+3) + N8(K)*SURF_VEC(3)*WT
               ENDDO
            ENDDO
         ENDDO
         DO J=1,SIZE(PPE,2)
            PPE(1:NDOF,J) = PPE(1:NDOF,J) + UNIT_PPE(1:NDOF)*PRESS(3,J)
         ENDDO
      ENDIF

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN
         CALL ELMDIS
         CALL BM_ANS_Q8_AT(XY8, ZERO, ZERO, BM, DETJ)
         STRAIN0 = MATMUL(BM, UEL(1:NDOF))
         N0V = MATMUL(SHELL_A, STRAIN0)

         SIG0 = ZERO
         SIG0(1,1) = N0V(1)
         SIG0(2,2) = N0V(2)
         SIG0(1,2) = N0V(3)
         SIG0(2,1) = N0V(3)

         KG8 = ZERO
         DO I=1,3
            DO J=1,3
               R = GP3(I)
               S = GP3(J)
               WT = W3(I)*W3(J)
               CALL Q8_DXY(XY8, R, S, N8, DNDX, DNDY, DETJ)
               DO K=1,8
                  DO L=1,8
                     KG8(K,L) = KG8(K,L) + ( DNDX(K)*(SIG0(1,1)*DNDX(L) + SIG0(1,2)*DNDY(L)) +                     &
                                             DNDY(K)*(SIG0(2,1)*DNDX(L) + SIG0(2,2)*DNDY(L)) ) * WT * DABS(DETJ)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         KED(1:NDOF,1:NDOF) = ZERO
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

      SUBROUTINE SHAPE_Q8_STD_DERIVS ( XI, ETA, NVAL, DN )
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

      DN(1,1) = 0.25D0*(ONE-ETA)*(2.0D0*XI+ETA)
      DN(1,2) = 0.25D0*(ONE-ETA)*(2.0D0*XI-ETA)
      DN(1,3) = 0.25D0*(ONE+ETA)*(2.0D0*XI+ETA)
      DN(1,4) = 0.25D0*(ONE+ETA)*(2.0D0*XI-ETA)
      DN(1,5) = -XI*(ONE-ETA)
      DN(1,6) = 0.5D0*(ONE-ETA*ETA)
      DN(1,7) = -XI*(ONE+ETA)
      DN(1,8) = -0.5D0*(ONE-ETA*ETA)

      DN(2,1) = 0.25D0*(ONE-XI)*(2.0D0*ETA+XI)
      DN(2,2) = 0.25D0*(ONE+XI)*(2.0D0*ETA-XI)
      DN(2,3) = 0.25D0*(ONE+XI)*(2.0D0*ETA+XI)
      DN(2,4) = 0.25D0*(ONE-XI)*(2.0D0*ETA-XI)
      DN(2,5) = -0.5D0*(ONE-XI*XI)
      DN(2,6) = -(ONE+XI)*ETA
      DN(2,7) = 0.5D0*(ONE-XI*XI)
      DN(2,8) = -(ONE-XI)*ETA
      END SUBROUTINE SHAPE_Q8_STD_DERIVS

      SUBROUTINE SHAPE_Q8_STD ( XI, ETA, NVAL, G1, G2, XYZN, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XI, ETA, XYZN(8,3)
      REAL(DOUBLE), INTENT(OUT) :: NVAL(8), G1(3), G2(3), DETJ
      REAL(DOUBLE) :: DN(2,8), GV(3)
      CALL SHAPE_Q8_STD_DERIVS(XI, ETA, NVAL, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      CALL CROSS3(G1, G2, GV)
      DETJ = VNORM(GV)
      END SUBROUTINE SHAPE_Q8_STD

      SUBROUTINE FIXED_FRAME_Q8 ( XYZN, E1OUT, E2OUT, E3OUT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3)
      REAL(DOUBLE), INTENT(OUT) :: E1OUT(3), E2OUT(3), E3OUT(3)
      REAL(DOUBLE) :: C(4,3), TMP(3), NM
      C = XYZN(1:4,:)
      E1OUT = 0.5D0*((C(2,:)-C(1,:)) + (C(3,:)-C(4,:)))
      NM = VNORM(E1OUT)
      IF (NM > 1.0D-15) THEN
         E1OUT = E1OUT/NM
      ELSE
         E1OUT = (/ONE, ZERO, ZERO/)
      ENDIF
      CALL CROSS3(C(3,:)-C(1,:), C(4,:)-C(2,:), E3OUT)
      NM = VNORM(E3OUT)
      IF (NM > 1.0D-15) THEN
         E3OUT = E3OUT/NM
      ELSE
         E3OUT = (/ZERO, ZERO, ONE/)
      ENDIF
      CALL CROSS3(E3OUT, E1OUT, E2OUT)
      NM = VNORM(E2OUT)
      IF (NM > 1.0D-15) E2OUT = E2OUT/NM
      CALL CROSS3(E2OUT, E3OUT, TMP)
      NM = VNORM(TMP)
      IF (NM > 1.0D-15) E1OUT = TMP/NM
      END SUBROUTINE FIXED_FRAME_Q8

      SUBROUTINE BUILD_LOCAL_XY_Q8 ( XYZN, E1IN, E2IN, XYOUT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), E1IN(3), E2IN(3)
      REAL(DOUBLE), INTENT(OUT) :: XYOUT(8,2)
      REAL(DOUBLE) :: ORG(3)
      INTEGER(LONG) :: II
      ORG = 0.25D0*(XYZN(1,:) + XYZN(2,:) + XYZN(3,:) + XYZN(4,:))
      DO II=1,8
         XYOUT(II,1) = DOT_PRODUCT(XYZN(II,:) - ORG, E1IN)
         XYOUT(II,2) = DOT_PRODUCT(XYZN(II,:) - ORG, E2IN)
      ENDDO
      END SUBROUTINE BUILD_LOCAL_XY_Q8

      SUBROUTINE Q8_DXY ( XYLOC, XI, ETA, NVAL, DNDX, DNDY, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: NVAL(8), DNDX(8), DNDY(8), DETJ
      REAL(DOUBLE) :: DN(2,8), JACM(2,2), JINV(2,2), DET
      INTEGER(LONG) :: II
      CALL SHAPE_Q8_STD_DERIVS(XI, ETA, NVAL, DN)
      JACM = MATMUL(DN, XYLOC)
      DET = JACM(1,1)*JACM(2,2) - JACM(1,2)*JACM(2,1)
      DETJ = DET
      IF (DABS(DET) <= 1.0D-20) THEN
         DNDX = ZERO
         DNDY = ZERO
      ELSE
         JINV(1,1) =  JACM(2,2)/DET
         JINV(1,2) = -JACM(1,2)/DET
         JINV(2,1) = -JACM(2,1)/DET
         JINV(2,2) =  JACM(1,1)/DET
         DO II=1,8
            DNDX(II) = JINV(1,1)*DN(1,II) + JINV(1,2)*DN(2,II)
            DNDY(II) = JINV(2,1)*DN(1,II) + JINV(2,2)*DN(2,II)
         ENDDO
      ENDIF
      END SUBROUTINE Q8_DXY

      SUBROUTINE BM_STD_Q8_AT ( XYLOC, XI, ETA, BMOUT, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BMOUT(3,48), DETJ
      REAL(DOUBLE) :: NVAL(8), DNDX(8), DNDY(8)
      INTEGER(LONG) :: II, COL
      CALL Q8_DXY(XYLOC, XI, ETA, NVAL, DNDX, DNDY, DETJ)
      BMOUT = ZERO
      DO II=1,8
         COL = (II-1)*6
         BMOUT(1,COL+1) = DNDX(II)
         BMOUT(2,COL+2) = DNDY(II)
         BMOUT(3,COL+1) = DNDY(II)
         BMOUT(3,COL+2) = DNDX(II)
      ENDDO
      END SUBROUTINE BM_STD_Q8_AT

      SUBROUTINE BB_Q8_AT ( XYLOC, XI, ETA, BBOUT, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BBOUT(3,48), DETJ
      REAL(DOUBLE) :: NVAL(8), DNDX(8), DNDY(8)
      INTEGER(LONG) :: II, COL
      CALL Q8_DXY(XYLOC, XI, ETA, NVAL, DNDX, DNDY, DETJ)
      BBOUT = ZERO
      DO II=1,8
         COL = (II-1)*6
         BBOUT(1,COL+5) =  DNDX(II)
         BBOUT(2,COL+4) = -DNDY(II)
         BBOUT(3,COL+5) =  DNDY(II)
         BBOUT(3,COL+4) = -DNDX(II)
      ENDDO
      END SUBROUTINE BB_Q8_AT

      SUBROUTINE BS_STD_Q8_AT ( XYLOC, XI, ETA, BSOUT, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BSOUT(2,48), DETJ
      REAL(DOUBLE) :: NVAL(8), DNDX(8), DNDY(8)
      INTEGER(LONG) :: II, COL
      CALL Q8_DXY(XYLOC, XI, ETA, NVAL, DNDX, DNDY, DETJ)
      BSOUT = ZERO
      DO II=1,8
         COL = (II-1)*6
         BSOUT(1,COL+3) = DNDX(II)
         BSOUT(1,COL+5) = NVAL(II)
         BSOUT(2,COL+3) = DNDY(II)
         BSOUT(2,COL+4) = -NVAL(II)
      ENDDO
      END SUBROUTINE BS_STD_Q8_AT

      SUBROUTINE BDRILL_Q8_AT ( XYLOC, XI, ETA, BDOUT, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BDOUT(1,48), DETJ
      REAL(DOUBLE) :: NVAL(8), DNDX(8), DNDY(8)
      INTEGER(LONG) :: II, COL
      CALL Q8_DXY(XYLOC, XI, ETA, NVAL, DNDX, DNDY, DETJ)
      BDOUT = ZERO
      DO II=1,8
         COL = (II-1)*6
         BDOUT(1,COL+1) = -0.5D0*DNDY(II)
         BDOUT(1,COL+2) =  0.5D0*DNDX(II)
         BDOUT(1,COL+6) =  NVAL(II)
      ENDDO
      END SUBROUTINE BDRILL_Q8_AT

      SUBROUTINE BM_ANS_Q8_AT ( XYLOC, XI, ETA, BMOUT, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BMOUT(3,48), DETJ
      REAL(DOUBLE) :: B1(1,48), B2(1,48), B3(1,48)
      CALL INTERP_EXX_Q8(XYLOC, XI, ETA, B1, DETJ)
      CALL INTERP_EYY_Q8(XYLOC, XI, ETA, B2, DETJ)
      CALL INTERP_EXY_Q8(XYLOC, XI, ETA, B3, DETJ)
      BMOUT = ZERO
      BMOUT(1,:) = B1(1,:)
      BMOUT(2,:) = B2(1,:)
      BMOUT(3,:) = B3(1,:)
      END SUBROUTINE BM_ANS_Q8_AT

      SUBROUTINE BS_ANS_Q8_AT ( XYLOC, XI, ETA, BSOUT, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BSOUT(2,48), DETJ
      REAL(DOUBLE) :: B1(1,48), B2(1,48)
      CALL INTERP_GXZ_Q8(XYLOC, XI, ETA, B1, DETJ)
      CALL INTERP_GYZ_Q8(XYLOC, XI, ETA, B2, DETJ)
      BSOUT = ZERO
      BSOUT(1,:) = B1(1,:)
      BSOUT(2,:) = B2(1,:)
      END SUBROUTINE BS_ANS_Q8_AT

      SUBROUTINE INTERP_EXX_Q8 ( XYLOC, XI, ETA, BOUT, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BOUT(1,48), DETJ
      REAL(DOUBLE) :: BM(3,48), ETA_N(3), B, HL, HR
      INTEGER(LONG) :: K1
      ETA_N = (/-DSQRT(3.0D0/5.0D0), ZERO, DSQRT(3.0D0/5.0D0)/)
      HL = 0.5D0*(ONE-XI)
      HR = 0.5D0*(ONE+XI)
      BOUT = ZERO
      DETJ = ZERO
      DO K1=1,3
         CALL BM_STD_Q8_AT(XYLOC, -ONE, ETA_N(K1), BM, DETJ)
         B = L3(ETA, ETA_N, K1)
         BOUT = BOUT + HL*B*BM(1:1,:)
         CALL BM_STD_Q8_AT(XYLOC,  ONE, ETA_N(K1), BM, DETJ)
         BOUT = BOUT + HR*B*BM(1:1,:)
      ENDDO
      END SUBROUTINE INTERP_EXX_Q8

      SUBROUTINE INTERP_EYY_Q8 ( XYLOC, XI, ETA, BOUT, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BOUT(1,48), DETJ
      REAL(DOUBLE) :: BM(3,48), XI_N(3), B, HB, HT
      INTEGER(LONG) :: K1
      XI_N = (/-DSQRT(3.0D0/5.0D0), ZERO, DSQRT(3.0D0/5.0D0)/)
      HB = 0.5D0*(ONE-ETA)
      HT = 0.5D0*(ONE+ETA)
      BOUT = ZERO
      DETJ = ZERO
      DO K1=1,3
         CALL BM_STD_Q8_AT(XYLOC, XI_N(K1), -ONE, BM, DETJ)
         B = L3(XI, XI_N, K1)
         BOUT = BOUT + HB*B*BM(2:2,:)
         CALL BM_STD_Q8_AT(XYLOC, XI_N(K1),  ONE, BM, DETJ)
         BOUT = BOUT + HT*B*BM(2:2,:)
      ENDDO
      END SUBROUTINE INTERP_EYY_Q8

      SUBROUTINE INTERP_EXY_Q8 ( XYLOC, XI, ETA, BOUT, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BOUT(1,48), DETJ
      REAL(DOUBLE) :: BM(3,48), A, H(4)
      A = ONE/DSQRT(3.0D0)
      H(1) = 0.25D0*(ONE-XI/A)*(ONE-ETA/A)
      H(2) = 0.25D0*(ONE+XI/A)*(ONE-ETA/A)
      H(3) = 0.25D0*(ONE+XI/A)*(ONE+ETA/A)
      H(4) = 0.25D0*(ONE-XI/A)*(ONE+ETA/A)
      BOUT = ZERO
      CALL BM_STD_Q8_AT(XYLOC, -A, -A, BM, DETJ); BOUT = BOUT + H(1)*BM(3:3,:)
      CALL BM_STD_Q8_AT(XYLOC,  A, -A, BM, DETJ); BOUT = BOUT + H(2)*BM(3:3,:)
      CALL BM_STD_Q8_AT(XYLOC,  A,  A, BM, DETJ); BOUT = BOUT + H(3)*BM(3:3,:)
      CALL BM_STD_Q8_AT(XYLOC, -A,  A, BM, DETJ); BOUT = BOUT + H(4)*BM(3:3,:)
      END SUBROUTINE INTERP_EXY_Q8

      SUBROUTINE INTERP_GXZ_Q8 ( XYLOC, XI, ETA, BOUT, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BOUT(1,48), DETJ
      REAL(DOUBLE) :: BS(2,48), A, H(4)
      A = ONE/DSQRT(3.0D0)
      H(1) = 0.25D0*(ONE-XI/A)*(ONE-ETA)
      H(2) = 0.25D0*(ONE+XI/A)*(ONE-ETA)
      H(3) = 0.25D0*(ONE+XI/A)*(ONE+ETA)
      H(4) = 0.25D0*(ONE-XI/A)*(ONE+ETA)
      BOUT = ZERO
      CALL BS_STD_Q8_AT(XYLOC, -A, -ONE, BS, DETJ); BOUT = BOUT + H(1)*BS(1:1,:)
      CALL BS_STD_Q8_AT(XYLOC,  A, -ONE, BS, DETJ); BOUT = BOUT + H(2)*BS(1:1,:)
      CALL BS_STD_Q8_AT(XYLOC,  A,  ONE, BS, DETJ); BOUT = BOUT + H(3)*BS(1:1,:)
      CALL BS_STD_Q8_AT(XYLOC, -A,  ONE, BS, DETJ); BOUT = BOUT + H(4)*BS(1:1,:)
      END SUBROUTINE INTERP_GXZ_Q8

      SUBROUTINE INTERP_GYZ_Q8 ( XYLOC, XI, ETA, BOUT, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BOUT(1,48), DETJ
      REAL(DOUBLE) :: BS(2,48), A, H(4)
      A = ONE/DSQRT(3.0D0)
      H(1) = 0.25D0*(ONE-XI)*(ONE-ETA/A)
      H(2) = 0.25D0*(ONE+XI)*(ONE-ETA/A)
      H(3) = 0.25D0*(ONE+XI)*(ONE+ETA/A)
      H(4) = 0.25D0*(ONE-XI)*(ONE+ETA/A)
      BOUT = ZERO
      CALL BS_STD_Q8_AT(XYLOC, -ONE, -A, BS, DETJ); BOUT = BOUT + H(1)*BS(2:2,:)
      CALL BS_STD_Q8_AT(XYLOC,  ONE, -A, BS, DETJ); BOUT = BOUT + H(2)*BS(2:2,:)
      CALL BS_STD_Q8_AT(XYLOC,  ONE,  A, BS, DETJ); BOUT = BOUT + H(3)*BS(2:2,:)
      CALL BS_STD_Q8_AT(XYLOC, -ONE,  A, BS, DETJ); BOUT = BOUT + H(4)*BS(2:2,:)
      END SUBROUTINE INTERP_GYZ_Q8

      FUNCTION L3 ( X, XN, K1 ) RESULT(VAL)
      REAL(DOUBLE), INTENT(IN) :: X, XN(3)
      INTEGER(LONG), INTENT(IN):: K1
      REAL(DOUBLE) :: VAL
      INTEGER(LONG) :: J1
      REAL(DOUBLE) :: NUM, DEN
      NUM = ONE
      DEN = ONE
      DO J1=1,3
         IF (J1 /= K1) THEN
            NUM = NUM*(X - XN(J1))
            DEN = DEN*(XN(K1) - XN(J1))
         ENDIF
      ENDDO
      IF (DABS(DEN) > 1.0D-14) THEN
         VAL = NUM/DEN
      ELSE
         VAL = ZERO
      ENDIF
      END FUNCTION L3

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

      END SUBROUTINE CQUAD8_ANS8BDG6
