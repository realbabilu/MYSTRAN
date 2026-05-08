! --- MITC3+_add begin --- !
      SUBROUTINE TPLT_MITC3P ( OPT, AREA, X2E, X3E, Y3E, BIG_BB )

! Experimental MITC3+ triangular shell bending/shear kernel for PARAM,TRIA3TYP,MITC3+.
! This is a first MYSTRAN port of the uploaded Python implementation in:
!   D:\mystran2\mitc3+\mitc3plus.py
!
! Scope in this first port:
!   - 18 shell DOF: u, v, w, rx, ry, rz at each of 3 CTRIA3 grids
!   - adds only the plate/shell bending + transverse shear + light drilling penalty block
!   - uses 2 internal bubble rotational DOF condensed at element level
!   - membrane remains supplied by TMEM1 in TREL1 for legacy CTRIA3
!   - pressure/thermal/stress recovery are intentionally left to legacy paths for now

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG, BUGOUT, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, THREE
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF, ONLY            :  EID, ELDOF, KE, SHELL_D, SHELL_T, TYPE
      USE OUTA_HERE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TPLT_MITC3P'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      REAL(DOUBLE), INTENT(IN)        :: AREA
      REAL(DOUBLE), INTENT(IN)        :: X2E
      REAL(DOUBLE), INTENT(IN)        :: X3E
      REAL(DOUBLE), INTENT(IN)        :: Y3E
      REAL(DOUBLE), INTENT(OUT)       :: BIG_BB(3,ELDOF,1)

      INTEGER(LONG), PARAMETER        :: IDX_B(8) = (/ 4, 5, 10, 11, 16, 17, 19, 20 /)
      INTEGER(LONG), PARAMETER        :: IDX_S(11)= (/ 3, 4, 5,  9, 10, 11, 15, 16, 17, 19, 20 /)
      INTEGER(LONG), PARAMETER        :: IDX_RZ(3)= (/ 6, 12, 18 /)

      INTEGER(LONG)                   :: I, J, K, L, GP
      REAL(DOUBLE)                    :: XY(3,2), JMAT(2,2), JINV(2,2), DETJ
      REAL(DOUBLE)                    :: COV_S(2,2), DRILL_PEN
      REAL(DOUBLE)                    :: KFULL(20,20), KAA(18,18), KAB(18,2), KBA(2,18), KBB(2,2), KBB_INV(2,2), KCOND(18,18)
      REAL(DOUBLE)                    :: BB(3,8), BS(2,11), KEI, FACTOR
      REAL(DOUBLE)                    :: GAUSS_R(7), GAUSS_S(7), GAUSS_W(7)
      REAL(DOUBLE)                    :: DABS_SHELL
      REAL(DOUBLE)                    :: EPS1

      BIG_BB = ZERO
      IF (OPT(4) /= 'Y') RETURN

      EPS1 = EPSIL(1)
      IF (AREA <= EPS1) THEN
         WRITE(ERR,*) ' *ERROR: ', TYPE, ' ELEMENT ', EID, ' HAS NONPOSITIVE AREA IN TPLT_MITC3P'
         WRITE(F06,*) ' *ERROR: ', TYPE, ' ELEMENT ', EID, ' HAS NONPOSITIVE AREA IN TPLT_MITC3P'
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      XY = ZERO
      XY(2,1) = X2E
      XY(3,1) = X3E
      XY(3,2) = Y3E

! --- MITC3+_add begin --- !
! Match the uploaded Python implementation's Jacobian layout:
!   J = d(x,y)/d(r,s) with rows = parametric directions (r,s)
! so:
!   row 1 = [dx/dr, dy/dr] = [X2E, 0]
!   row 2 = [dx/ds, dy/ds] = [X3E, Y3E]
! This is the transpose of the legacy "component-first" layout and is required
! for the covariant MITC3+ shear operator to match the Python reference.
      JMAT(1,1) = X2E
      JMAT(1,2) = ZERO
      JMAT(2,1) = X3E
      JMAT(2,2) = Y3E
! --- MITC3+_add end --- !
      DETJ      = JMAT(1,1)*JMAT(2,2) - JMAT(1,2)*JMAT(2,1)
      IF (DABS(DETJ) <= EPS1) THEN
         WRITE(ERR,*) ' *ERROR: ', TYPE, ' ELEMENT ', EID, ' HAS SINGULAR JACOBIAN IN TPLT_MITC3P'
         WRITE(F06,*) ' *ERROR: ', TYPE, ' ELEMENT ', EID, ' HAS SINGULAR JACOBIAN IN TPLT_MITC3P'
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      JINV(1,1) =  JMAT(2,2)/DETJ
      JINV(1,2) = -JMAT(1,2)/DETJ
      JINV(2,1) = -JMAT(2,1)/DETJ
      JINV(2,2) =  JMAT(1,1)/DETJ

      COV_S = MATMUL(TRANSPOSE(JINV), MATMUL(SHELL_T, JINV))
      DABS_SHELL = MAX(DABS(SHELL_T(1,1)), DABS(SHELL_T(2,2)))
      DRILL_PEN  = 1.0D-05*DABS_SHELL

      CALL MITC3P_GAUSS_7PT(GAUSS_R, GAUSS_S, GAUSS_W)

      KFULL = ZERO
      DO GP=1,7
         CALL MITC3P_BB_AT(GAUSS_R(GP), GAUSS_S(GP), JINV, BB)
         CALL MITC3P_BS_AT(GAUSS_R(GP), GAUSS_S(GP), XY, JINV, JMAT, BS)
         FACTOR = GAUSS_W(GP)*DETJ

         DO I=1,8
            DO J=1,8
               KEI = ZERO
               DO K=1,3
                  DO L=1,3
                     KEI = KEI + BB(K,I)*SHELL_D(K,L)*BB(L,J)
                  ENDDO
               ENDDO
               KFULL(IDX_B(I),IDX_B(J)) = KFULL(IDX_B(I),IDX_B(J)) + FACTOR*KEI
            ENDDO
         ENDDO

         DO I=1,11
            DO J=1,11
               KEI = ZERO
               DO K=1,2
                  DO L=1,2
                     KEI = KEI + BS(K,I)*COV_S(K,L)*BS(L,J)
                  ENDDO
               ENDDO
               KFULL(IDX_S(I),IDX_S(J)) = KFULL(IDX_S(I),IDX_S(J)) + FACTOR*KEI
            ENDDO
         ENDDO
      ENDDO

      DO I=1,3
         KFULL(IDX_RZ(I),IDX_RZ(I)) = KFULL(IDX_RZ(I),IDX_RZ(I)) + DRILL_PEN
      ENDDO

      KAA = KFULL(1:18,1:18)
      KAB = KFULL(1:18,19:20)
      KBA = KFULL(19:20,1:18)
      KBB = KFULL(19:20,19:20)
      CALL MITC3P_INV2(KBB, KBB_INV)
      KCOND = KAA - MATMUL(KAB, MATMUL(KBB_INV, KBA))

! --- MITC3+_add begin --- !
      BUGOUT = 'Y'
      WRITE(ERR,*) 'MITC3P_DEBUG_BEGIN EID=', EID
      CALL MITC3P_WRITE_MAT(ERR, 'SHELL_D', SHELL_D, 3, 3)
      CALL MITC3P_WRITE_MAT(ERR, 'SHELL_T', SHELL_T, 2, 2)
      CALL MITC3P_WRITE_MAT(ERR, 'COV_S',   COV_S  , 2, 2)
      WRITE(ERR,'(A,1P,E16.8)') 'MITC3P_DRILL_PEN ', DRILL_PEN
      CALL MITC3P_WRITE_MAT(ERR, 'KAA',   KAA,   18, 18)
      CALL MITC3P_WRITE_MAT(ERR, 'KAB',   KAB,   18,  2)
      CALL MITC3P_WRITE_MAT(ERR, 'KBB',   KBB,    2,  2)
      CALL MITC3P_WRITE_MAT(ERR, 'KCOND', KCOND, 18, 18)
      WRITE(ERR,*) 'MITC3P_DEBUG_END EID=', EID
! --- MITC3+_add end --- !

      DO I=1,18
         DO J=1,18
            KE(I,J) = KE(I,J) + KCOND(I,J)
         ENDDO
      ENDDO

      CALL MITC3P_BB_AT(ONE/THREE, ONE/THREE, JINV, BB)
      DO I=1,3
         BIG_BB(I, 4,1)  = BB(I,1)
         BIG_BB(I, 5,1)  = BB(I,2)
         BIG_BB(I,10,1)  = BB(I,3)
         BIG_BB(I,11,1)  = BB(I,4)
         BIG_BB(I,16,1)  = BB(I,5)
         BIG_BB(I,17,1)  = BB(I,6)
      ENDDO

      CONTAINS

      SUBROUTINE MITC3P_GAUSS_7PT(R, S, W)
         REAL(DOUBLE), INTENT(OUT) :: R(7), S(7), W(7)
         REAL(DOUBLE) :: A1, B1, A2, B2, C, W1, W2, W3
         A1 = 0.1012865073235D0
         B1 = 0.7974269853531D0
         A2 = 0.4701420641051D0
         B2 = 0.0597158717898D0
         C  = ONE/THREE
         W1 = 0.1259391805448D0
         W2 = 0.1323941527885D0
         W3 = 0.2250000000000D0
         R = (/ A1, B1, A1, A2, B2, A2, C /)
         S = (/ A1, A1, B1, A2, A2, B2, C /)
         W = 0.5D0*(/ W1, W1, W1, W2, W2, W2, W3 /)
      END SUBROUTINE MITC3P_GAUSS_7PT

      SUBROUTINE MITC3P_INV2(A, AINV)
         REAL(DOUBLE), INTENT(IN)  :: A(2,2)
         REAL(DOUBLE), INTENT(OUT) :: AINV(2,2)
         REAL(DOUBLE) :: DET2
         DET2 = A(1,1)*A(2,2) - A(1,2)*A(2,1)
         IF (DABS(DET2) <= EPS1) THEN
            WRITE(ERR,*) ' *ERROR: ', TYPE, ' ELEMENT ', EID, ' HAS SINGULAR BUBBLE BLOCK IN TPLT_MITC3P'
            WRITE(F06,*) ' *ERROR: ', TYPE, ' ELEMENT ', EID, ' HAS SINGULAR BUBBLE BLOCK IN TPLT_MITC3P'
            CALL OUTA_HERE ( 'Y' )
         ENDIF
         AINV(1,1) =  A(2,2)/DET2
         AINV(1,2) = -A(1,2)/DET2
         AINV(2,1) = -A(2,1)/DET2
         AINV(2,2) =  A(1,1)/DET2
      END SUBROUTINE MITC3P_INV2

      SUBROUTINE MITC3P_DH(DH)
         REAL(DOUBLE), INTENT(OUT) :: DH(3,2)
         DH(1,1) = -ONE ; DH(1,2) = -ONE
         DH(2,1) =  ONE ; DH(2,2) =  ZERO
         DH(3,1) =  ZERO; DH(3,2) =  ONE
      END SUBROUTINE MITC3P_DH

      SUBROUTINE MITC3P_FI(R, S, FIALL)
         REAL(DOUBLE), INTENT(IN)  :: R, S
         REAL(DOUBLE), INTENT(OUT) :: FIALL(4)
         REAL(DOUBLE) :: H(3), F4V
         H(1) = ONE - R - S
         H(2) = R
         H(3) = S
         F4V  = 27.0D0*R*S*(ONE - R - S)
         FIALL(1) = H(1) - F4V/THREE
         FIALL(2) = H(2) - F4V/THREE
         FIALL(3) = H(3) - F4V/THREE
         FIALL(4) = F4V
      END SUBROUTINE MITC3P_FI

      SUBROUTINE MITC3P_DFI(R, S, DFI)
         REAL(DOUBLE), INTENT(IN)  :: R, S
         REAL(DOUBLE), INTENT(OUT) :: DFI(4,2)
         REAL(DOUBLE) :: DH(3,2), DF4(2), T
         CALL MITC3P_DH(DH)
         T = ONE - R - S
         DF4(1) = 27.0D0*S*(T - R)
         DF4(2) = 27.0D0*R*(T - S)
         DFI(1,1) = DH(1,1) - DF4(1)/THREE
         DFI(1,2) = DH(1,2) - DF4(2)/THREE
         DFI(2,1) = DH(2,1) - DF4(1)/THREE
         DFI(2,2) = DH(2,2) - DF4(2)/THREE
         DFI(3,1) = DH(3,1) - DF4(1)/THREE
         DFI(3,2) = DH(3,2) - DF4(2)/THREE
         DFI(4,1) = DF4(1)
         DFI(4,2) = DF4(2)
      END SUBROUTINE MITC3P_DFI

      SUBROUTINE MITC3P_BB_AT(R, S, JI, BBOUT)
         REAL(DOUBLE), INTENT(IN)  :: R, S, JI(2,2)
         REAL(DOUBLE), INTENT(OUT) :: BBOUT(3,8)
         REAL(DOUBLE) :: DFI(4,2), DFI_XY(4,2)
         INTEGER(LONG) :: II
         CALL MITC3P_DFI(R, S, DFI)
         DFI_XY = MATMUL(DFI, TRANSPOSE(JI))
         BBOUT = ZERO
         DO II=1,4
            BBOUT(1,2*II  ) =  DFI_XY(II,1)
            BBOUT(2,2*II-1) = -DFI_XY(II,2)
            BBOUT(3,2*II  ) =  DFI_XY(II,2)
            BBOUT(3,2*II-1) = -DFI_XY(II,1)
         ENDDO
      END SUBROUTINE MITC3P_BB_AT

      SUBROUTINE MITC3P_ERT_EST_ROWS(R, S, XYL, JI, JM, ERT, EST)
         REAL(DOUBLE), INTENT(IN)  :: R, S, XYL(3,2), JI(2,2), JM(2,2)
         REAL(DOUBLE), INTENT(OUT) :: ERT(11), EST(11)
         REAL(DOUBLE) :: DH(3,2), DH_XY(3,2), FIALL(4)
         INTEGER(LONG) :: II
         CALL MITC3P_DH(DH)
         CALL MITC3P_FI(R, S, FIALL)
         DH_XY = MATMUL(DH, TRANSPOSE(JI))
         ERT = ZERO
         EST = ZERO
         DO II=1,3
            ERT(3*II-2) = JM(1,1)*DH_XY(II,1) + JM(1,2)*DH_XY(II,2)
            ERT(3*II-1) = JM(1,2)*FIALL(II)
            ERT(3*II  ) = -JM(1,1)*FIALL(II)
            EST(3*II-2) = JM(2,1)*DH_XY(II,1) + JM(2,2)*DH_XY(II,2)
            EST(3*II-1) = JM(2,2)*FIALL(II)
            EST(3*II  ) = -JM(2,1)*FIALL(II)
         ENDDO
         ERT(10) =  JM(1,2)*FIALL(4)
         ERT(11) = -JM(1,1)*FIALL(4)
         EST(10) =  JM(2,2)*FIALL(4)
         EST(11) = -JM(2,1)*FIALL(4)
      END SUBROUTINE MITC3P_ERT_EST_ROWS

      SUBROUTINE MITC3P_BS_AT(R, S, XYL, JI, JM, BSOUT)
         REAL(DOUBLE), INTENT(IN)  :: R, S, XYL(3,2), JI(2,2), JM(2,2)
         REAL(DOUBLE), INTENT(OUT) :: BSOUT(2,11)
         REAL(DOUBLE) :: APT(2), BPT(2), CPT(2), DPT(2), EPT(2), FPT(2)
         REAL(DOUBLE) :: ERTA(11), ESTA(11), ERTB(11), ESTB(11), ERTC(11), ESTC(11)
         REAL(DOUBLE) :: ERTD(11), ESTD(11), ERTE(11), ESTE(11), ERTF(11), ESTF(11)
         REAL(DOUBLE) :: CONST_ERT(11), CONST_EST(11), CHAT(11), FAC_RT, FAC_ST, DPAR

         APT = (/ ONE/6.0D0, TWO/THREE /)
         BPT = (/ TWO/THREE, ONE/6.0D0 /)
         CPT = (/ ONE/6.0D0, ONE/6.0D0 /)
         DPAR = 1.0D0/10000.0D0
         DPT = (/ ONE/THREE + DPAR   , ONE/THREE - TWO*DPAR /)
         EPT = (/ ONE/THREE - TWO*DPAR, ONE/THREE + DPAR    /)
         FPT = (/ ONE/THREE + DPAR   , ONE/THREE + DPAR    /)

         CALL MITC3P_ERT_EST_ROWS(APT(1), APT(2), XYL, JI, JM, ERTA, ESTA)
         CALL MITC3P_ERT_EST_ROWS(BPT(1), BPT(2), XYL, JI, JM, ERTB, ESTB)
         CALL MITC3P_ERT_EST_ROWS(CPT(1), CPT(2), XYL, JI, JM, ERTC, ESTC)
         CALL MITC3P_ERT_EST_ROWS(DPT(1), DPT(2), XYL, JI, JM, ERTD, ESTD)
         CALL MITC3P_ERT_EST_ROWS(EPT(1), EPT(2), XYL, JI, JM, ERTE, ESTE)
         CALL MITC3P_ERT_EST_ROWS(FPT(1), FPT(2), XYL, JI, JM, ERTF, ESTF)

         CONST_ERT = (TWO/THREE)*(ERTB - 0.5D0*ESTB) + (ONE/THREE)*(ERTC + ESTC)
         CONST_EST = (TWO/THREE)*(ESTA - 0.5D0*ERTA) + (ONE/THREE)*(ERTC + ESTC)
         CHAT      = (ERTF - ERTD) - (ESTF - ESTE)

         FAC_RT = (THREE*S - ONE)/THREE
         FAC_ST = (ONE - THREE*R)/THREE
         BSOUT(1,:) = CONST_ERT + FAC_RT*CHAT
         BSOUT(2,:) = CONST_EST + FAC_ST*CHAT
      END SUBROUTINE MITC3P_BS_AT

! --- MITC3+_add begin --- !
      SUBROUTINE MITC3P_WRITE_MAT(LU, NAME, MAT, NROW, NCOL)
         INTEGER(LONG), INTENT(IN)         :: LU, NROW, NCOL
         CHARACTER(*), INTENT(IN)          :: NAME
         REAL(DOUBLE), INTENT(IN)          :: MAT(NROW,NCOL)
         INTEGER(LONG)                     :: IR, IC
         WRITE(LU,*) 'MITC3P_MATRIX_BEGIN ', TRIM(NAME), NROW, NCOL
         DO IR=1,NROW
            WRITE(LU,'(A,I0,A,1P,100E16.8)') 'ROW ', IR, ' :', (MAT(IR,IC), IC=1,NCOL)
         ENDDO
         WRITE(LU,*) 'MITC3P_MATRIX_END ', TRIM(NAME)
      END SUBROUTINE MITC3P_WRITE_MAT
! --- MITC3+_add end --- !

      END SUBROUTINE TPLT_MITC3P
! --- MITC3+_add end --- !
