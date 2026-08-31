! #################################################################################################################################
! CQUAD8 MacNeal Q8 + drilling penalty shell for PARAM,QUAD8TYP,MACQ8D.

      SUBROUTINE CQUAD8_MACQ8D ( OPT, INT_ELEM_ID )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_STRESS_POINTS, SOL_NAME
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  ALPVEC, DT, EID, ELGP, KE, KED, ME, BE1, BE2, BE3, EPROP, FCONV, MASS_PER_UNIT_AREA,  &
                                         PPE, PRESS, PTE, SHELL_A, SHELL_D, SHELL_T, TREF, UEL, XEB, NUM_EMG_FATAL_ERRS,        &
                                         PCOMP_PROPS
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO
      USE PARAMS, ONLY                :  COUPMASS
      USE ELMDIS_Interface
      USE OUTA_HERE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CQUAD8_MACQ8D'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      INTEGER(LONG), PARAMETER        :: NNODE = 8
      INTEGER(LONG), PARAMETER        :: NDOF = 48
      INTEGER(LONG)                   :: I, J, K, L, IA, IB, JSUB, GP
      REAL(DOUBLE), PARAMETER         :: BETA_DRILL = 9.0D-6
      REAL(DOUBLE)                    :: XYZ(8,3), XY8(8,2), E1F(3), E2F(3), E3F(3)
      REAL(DOUBLE)                    :: BM(3,NDOF), BB(3,NDOF), BS(2,NDOF), BD(1,NDOF)
      REAL(DOUBLE)                    :: GP3(3), W3(3), R, S, WT, DETJ, THICK, GVAL, CDRILL, TBAR
      REAL(DOUBLE)                    :: M1(8,8), N8(8), DNDX(8), DNDY(8), UNIT_PPE(NDOF), UNIT_PTE(NDOF)
      REAL(DOUBLE)                    :: DXDR(3), DXDS(3), SURF_VEC(3), MASS_ELEM, MASS_NODE
      REAL(DOUBLE)                    :: CTE(3), THERMAL_RESULTANT(3)
      REAL(DOUBLE)                    :: SIG0(2,2), KG8(8,8), STRAIN0(3), N0V(3)
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
         WRITE(ERR,*) ' *ERROR: Code not written for composite material with PARAM,QUAD8TYP,MACQ8D'
         WRITE(F06,*) ' *ERROR: Code not written for composite material with PARAM,QUAD8TYP,MACQ8D'
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL LOAD_BASIC_COORDS_Q8(XYZ)
      CALL FIXED_FRAME_Q8(XYZ, E1F, E2F, E3F)
      CALL BUILD_LOCAL_XY_Q8(XYZ, E1F, E2F, XY8)

      GP3 = (/-DSQRT(3.0D0/5.0D0), ZERO, DSQRT(3.0D0/5.0D0)/)
      W3  = (/5.0D0/9.0D0, 8.0D0/9.0D0, 5.0D0/9.0D0/)
      THICK = EPROP(1)

      IF (OPT(1) == 'Y') THEN
         M1 = ZERO
         MASS_ELEM = ZERO
         DO I=1,3
            DO J=1,3
               R = GP3(I)
               S = GP3(J)
               WT = W3(I)*W3(J)
               CALL MACQ8_DXY(XY8, R, S, N8, DNDX, DNDY, DETJ)
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
               CALL BM_Q8_AT(XY8, E1F, E2F, R, S, BM, DETJ)
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
            TBAR = TBAR / 8.0D0 - TREF(1)
            PTE(1:NDOF,JSUB) = UNIT_PTE(1:NDOF) * TBAR
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
               CALL BM_Q8_AT(XY8, E1F, E2F, R, S, BM, DETJ)
               CALL BB_Q8_AT(XY8, E1F, E2F, R, S, BB, DETJ)
               CALL BS_Q8_AT(XY8, E1F, E2F, E3F, R, S, BS, DETJ)
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
         GVAL = ZERO
         IF (DABS(THICK) > 1.0D-20) GVAL = SHELL_T(1,1) / ((5.0D0/6.0D0)*THICK)
         CDRILL = BETA_DRILL * GVAL * THICK
         DO I=1,3
            DO J=1,3
               R = GP3(I)
               S = GP3(J)
               WT = W3(I)*W3(J)
               CALL BM_Q8_AT(XY8, E1F, E2F, R, S, BM, DETJ)
               CALL BB_Q8_AT(XY8, E1F, E2F, R, S, BB, DETJ)
               CALL BS_Q8_AT(XY8, E1F, E2F, E3F, R, S, BS, DETJ)
               CALL BDRILL_Q8_AT(XY8, R, S, BD, DETJ)
               KE = KE + WT*DABS(DETJ)*MATMUL(TRANSPOSE(BM), MATMUL(SHELL_A, BM))
               KE = KE + WT*DABS(DETJ)*MATMUL(TRANSPOSE(BB), MATMUL(SHELL_D, BB))
               KE = KE + WT*DABS(DETJ)*MATMUL(TRANSPOSE(BS), MATMUL(SHELL_T, BS))
               KE = KE + WT*DABS(DETJ)*CDRILL*MATMUL(TRANSPOSE(BD), BD)
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
         CALL BM_Q8_AT(XY8, E1F, E2F, ZERO, ZERO, BM, DETJ)
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
               CALL MACQ8_DXY(XY8, R, S, N8, DNDX, DNDY, DETJ)
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
      END SUBROUTINE SHAPE_Q8_STD_DERIVS

      SUBROUTINE FIXED_FRAME_Q8 ( XYZN, E1OUT, E2OUT, E3OUT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3)
      REAL(DOUBLE), INTENT(OUT) :: E1OUT(3), E2OUT(3), E3OUT(3)
      REAL(DOUBLE) :: NVAL(8), DN(2,8), G1(3), G2(3), TMP(3), NM
      CALL SHAPE_Q8_STD_DERIVS(ZERO, ZERO, NVAL, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      CALL CROSS3(G1, G2, E3OUT)
      NM = VNORM(E3OUT)
      IF (NM > 1.0D-15) THEN
         E3OUT = E3OUT / NM
      ELSE
         E3OUT = (/ZERO, ZERO, ONE/)
      ENDIF
      IF (E3OUT(3) < ZERO) E3OUT = -E3OUT
      E1OUT = G1 - DOT_PRODUCT(G1, E3OUT)*E3OUT
      NM = VNORM(E1OUT)
      IF (NM > 1.0D-15) THEN
         E1OUT = E1OUT / NM
      ELSE
         E1OUT = (/ONE, ZERO, ZERO/)
      ENDIF
      CALL CROSS3(E3OUT, E1OUT, E2OUT)
      NM = VNORM(E2OUT)
      IF (NM > 1.0D-15) E2OUT = E2OUT / NM
      END SUBROUTINE FIXED_FRAME_Q8

      SUBROUTINE BUILD_LOCAL_XY_Q8 ( XYZN, E1IN, E2IN, XYOUT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), E1IN(3), E2IN(3)
      REAL(DOUBLE), INTENT(OUT) :: XYOUT(8,2)
      REAL(DOUBLE) :: CENT(3), REL(3)
      INTEGER(LONG) :: II
      CENT = ZERO
      DO II=1,4
         CENT = CENT + XYZN(II,:)
      ENDDO
      CENT = CENT / 4.0D0
      DO II=1,8
         REL = XYZN(II,:) - CENT
         XYOUT(II,1) = DOT_PRODUCT(REL, E1IN)
         XYOUT(II,2) = DOT_PRODUCT(REL, E2IN)
      ENDDO
      END SUBROUTINE BUILD_LOCAL_XY_Q8

      SUBROUTINE MACQ8_DXY ( XYLOC, XI, ETA, NVAL, DNDX, DNDY, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: NVAL(8), DNDX(8), DNDY(8), DETJ
      REAL(DOUBLE) :: NSTD(8), DNSTD(2,8), DNMOD(2,8), N0(8), CORR(8), N9, DN9R, DN9S
      REAL(DOUBLE) :: X9Y9(2), XYC(8,2), XIM(8,8), AINV(8,8), J(2,2), JINV(2,2), DET2
      REAL(DOUBLE) :: PHI(8,2)
      INTEGER(LONG) :: II

      CALL SHAPE_Q8_STD_DERIVS(XI, ETA, NSTD, DNSTD)
      CALL SHAPE_Q8_STD_DERIVS(ZERO, ZERO, N0, DNMOD)

      X9Y9 = ZERO
      DO II=1,8
         X9Y9(1) = X9Y9(1) + N0(II)*XYLOC(II,1)
         X9Y9(2) = X9Y9(2) + N0(II)*XYLOC(II,2)
      ENDDO
      DO II=1,8
         XYC(II,1) = XYLOC(II,1) - X9Y9(1)
         XYC(II,2) = XYLOC(II,2) - X9Y9(2)
      ENDDO

      DO II=1,8
         PHI(II,1) = PARAM_R(II)
         PHI(II,2) = PARAM_S(II)
         XIM(II,1) = ONE
         XIM(II,2) = XYC(II,1)
         XIM(II,3) = XYC(II,2)
         XIM(II,4) = XYC(II,1)*XYC(II,1)
         XIM(II,5) = XYC(II,1)*XYC(II,2)
         XIM(II,6) = XYC(II,2)*XYC(II,2)
         XIM(II,7) = PHI(II,1)*PHI(II,1)*PHI(II,2)
         XIM(II,8) = PHI(II,1)*PHI(II,2)*PHI(II,2)
      ENDDO
      CALL INV8(XIM, AINV)

      DO II=1,8
         CORR(II) = AINV(1,II) - N0(II)
      ENDDO

      N9 = (ONE-XI*XI)*(ONE-ETA*ETA)
      DN9R = -TWO*XI*(ONE-ETA*ETA)
      DN9S = -TWO*ETA*(ONE-XI*XI)
      DO II=1,8
         NVAL(II) = NSTD(II) + N9*CORR(II)
         DNMOD(1,II) = DNSTD(1,II) + DN9R*CORR(II)
         DNMOD(2,II) = DNSTD(2,II) + DN9S*CORR(II)
      ENDDO

      J(1,1) = DOT_PRODUCT(DNSTD(1,:), XYLOC(:,1))
      J(1,2) = DOT_PRODUCT(DNSTD(1,:), XYLOC(:,2))
      J(2,1) = DOT_PRODUCT(DNSTD(2,:), XYLOC(:,1))
      J(2,2) = DOT_PRODUCT(DNSTD(2,:), XYLOC(:,2))
      CALL INV2(J, JINV, DET2)
      DETJ = DET2
      DNDX = JINV(1,1)*DNMOD(1,:) + JINV(1,2)*DNMOD(2,:)
      DNDY = JINV(2,1)*DNMOD(1,:) + JINV(2,2)*DNMOD(2,:)
      END SUBROUTINE MACQ8_DXY

      SUBROUTINE BM_Q8_AT ( XYLOC, E1IN, E2IN, XI, ETA, BMOUT, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), E1IN(3), E2IN(3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BMOUT(3,NDOF), DETJ
      REAL(DOUBLE) :: NVAL(8), DNDX(8), DNDY(8)
      INTEGER(LONG) :: II, COL
      CALL MACQ8_DXY(XYLOC, XI, ETA, NVAL, DNDX, DNDY, DETJ)
      BMOUT = ZERO
      DO II=1,8
         COL = 6*(II-1)
         BMOUT(1,COL+1:COL+3) = DNDX(II)*E1IN
         BMOUT(2,COL+1:COL+3) = DNDY(II)*E2IN
         BMOUT(3,COL+1:COL+3) = DNDY(II)*E1IN + DNDX(II)*E2IN
      ENDDO
      END SUBROUTINE BM_Q8_AT

      SUBROUTINE BB_Q8_AT ( XYLOC, E1IN, E2IN, XI, ETA, BBOUT, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), E1IN(3), E2IN(3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BBOUT(3,NDOF), DETJ
      REAL(DOUBLE) :: NVAL(8), DNDX(8), DNDY(8)
      INTEGER(LONG) :: II, COL
      CALL MACQ8_DXY(XYLOC, XI, ETA, NVAL, DNDX, DNDY, DETJ)
      BBOUT = ZERO
      DO II=1,8
         COL = 6*(II-1)
         BBOUT(1,COL+4:COL+6) = -DNDX(II)*E2IN
         BBOUT(2,COL+4:COL+6) =  DNDY(II)*E1IN
         BBOUT(3,COL+4:COL+6) =  DNDX(II)*E1IN - DNDY(II)*E2IN
      ENDDO
      END SUBROUTINE BB_Q8_AT

      SUBROUTINE BS_Q8_AT ( XYLOC, E1IN, E2IN, E3IN, XI, ETA, BSOUT, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), E1IN(3), E2IN(3), E3IN(3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BSOUT(2,NDOF), DETJ
      REAL(DOUBLE) :: NVAL(8), DNDX(8), DNDY(8)
      INTEGER(LONG) :: II, COL
      CALL MACQ8_DXY(XYLOC, XI, ETA, NVAL, DNDX, DNDY, DETJ)
      BSOUT = ZERO
      DO II=1,8
         COL = 6*(II-1)
         BSOUT(1,COL+1:COL+3) = DNDX(II)*E3IN
         BSOUT(1,COL+4:COL+6) = NVAL(II)*E2IN
         BSOUT(2,COL+1:COL+3) = DNDY(II)*E3IN
         BSOUT(2,COL+4:COL+6) = -NVAL(II)*E1IN
      ENDDO
      END SUBROUTINE BS_Q8_AT

      SUBROUTINE BDRILL_Q8_AT ( XYLOC, XI, ETA, BDOUT, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XYLOC(8,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BDOUT(1,NDOF), DETJ
      REAL(DOUBLE) :: NVAL(8), DNDX(8), DNDY(8)
      INTEGER(LONG) :: II, COL
      CALL MACQ8_DXY(XYLOC, XI, ETA, NVAL, DNDX, DNDY, DETJ)
      BDOUT = ZERO
      DO II=1,8
         COL = 6*(II-1)
         BDOUT(1,COL+4:COL+6) = NVAL(II)*E3F
      ENDDO
      END SUBROUTINE BDRILL_Q8_AT

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

      SUBROUTINE INV2 ( A, AINV, DET )
      REAL(DOUBLE), INTENT(IN)  :: A(2,2)
      REAL(DOUBLE), INTENT(OUT) :: AINV(2,2), DET
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

      SUBROUTINE INV8 ( A, AINV )
      REAL(DOUBLE), INTENT(IN)  :: A(8,8)
      REAL(DOUBLE), INTENT(OUT) :: AINV(8,8)
      REAL(DOUBLE) :: AUG(8,16), FAC, PIV, ROWTMP(16)
      INTEGER(LONG) :: II, JJ, PP
      AUG = ZERO
      AUG(:,1:8) = A
      DO II=1,8
         AUG(II,8+II) = ONE
      ENDDO
      DO II=1,8
         PP = II
         DO JJ=II+1,8
            IF (DABS(AUG(JJ,II)) > DABS(AUG(PP,II))) PP = JJ
         ENDDO
         IF (DABS(AUG(PP,II)) < 1.0D-20) THEN
            AINV = ZERO
            RETURN
         ENDIF
         IF (PP /= II) THEN
            ROWTMP = AUG(II,:)
            AUG(II,:) = AUG(PP,:)
            AUG(PP,:) = ROWTMP
         ENDIF
         PIV = AUG(II,II)
         AUG(II,:) = AUG(II,:) / PIV
         DO JJ=1,8
            IF (JJ == II) CYCLE
            FAC = AUG(JJ,II)
            AUG(JJ,:) = AUG(JJ,:) - FAC*AUG(II,:)
         ENDDO
      ENDDO
      AINV = AUG(:,9:16)
      END SUBROUTINE INV8

      REAL(DOUBLE) FUNCTION PARAM_R ( I )
      INTEGER(LONG), INTENT(IN) :: I
      REAL(DOUBLE), PARAMETER :: RV(8) = (/-ONE, ONE, ONE, -ONE, ZERO, ONE, ZERO, -ONE/)
      PARAM_R = RV(I)
      END FUNCTION PARAM_R

      REAL(DOUBLE) FUNCTION PARAM_S ( I )
      INTEGER(LONG), INTENT(IN) :: I
      REAL(DOUBLE), PARAMETER :: SV(8) = (/-ONE, -ONE, ONE, ONE, -ONE, ZERO, ONE, ZERO/)
      PARAM_S = SV(I)
      END FUNCTION PARAM_S

      END SUBROUTINE CQUAD8_MACQ8D



