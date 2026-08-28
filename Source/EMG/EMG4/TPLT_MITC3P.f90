! --- shell_renovation begin --- !
      SUBROUTINE TPLT_MITC3P ( OPT, AREA, X2E, X3E, Y3E, BIG_BB )

! Experimental MITC3+ triangular shell bending/shear kernel for PARAM,TRIA3TYP,MITC3+.
! This is a first MYSTRAN port of the uploaded Python implementation in:
!   D:\mystran2\mitc3+\mitc3plus.py
!
! Scope in this first port:
!   - 18 shell DOF: u, v, w, rx, ry, rz at each of 3 CTRIA3 grids
!   - adds only the plate/shell bending + transverse shear
!   - uses 2 internal bubble rotational DOF condensed at element level
!   - membrane remains supplied by TMEM1 in TREL1 for legacy CTRIA3
!   - pressure/thermal are intentionally left to legacy paths for now

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG, BUGOUT, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, THREE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF, ONLY            :  BE2, BE3, EID, ELDOF, KE, PHI_SQ, SE2, SE3, SHELL_A, SHELL_D, SHELL_T, TYPE
      USE MITC_STUF, ONLY             :  DIRECTOR
      USE OUTA_HERE_Interface
      USE CROSS_Interface
      USE MITC_COVARIANT_BASIS_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TPLT_MITC3P'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      REAL(DOUBLE), INTENT(IN)        :: AREA
      REAL(DOUBLE), INTENT(IN)        :: X2E
      REAL(DOUBLE), INTENT(IN)        :: X3E
      REAL(DOUBLE), INTENT(IN)        :: Y3E
      REAL(DOUBLE), INTENT(OUT)       :: BIG_BB(3,ELDOF,1)

      INTEGER(LONG), PARAMETER        :: IDX_M(6) = (/ 1, 2, 7, 8, 13, 14 /)

      INTEGER(LONG)                   :: I, J, K, L, GP
      REAL(DOUBLE)                    :: XY(3,2), JMAT(2,2), JINV(2,2), DETJ
      REAL(DOUBLE)                    :: COV_S(2,2)
      REAL(DOUBLE)                    :: KFULL(20,20), KAA(18,18), KAB(18,2), KBA(2,18), KBB(2,2), KBB_INV(2,2), KCOND(18,18)
      REAL(DOUBLE)                    :: KPHYS(18,18), KA(18,18), KOUT(18,18), TAE(18,18), TEA(18,18)
      REAL(DOUBLE)                    :: BB(3,8), BS(2,11), BM(3,6), KEI, FACTOR
      REAL(DOUBLE)                    :: BB_REC(3,18), BS_REC(2,18), DUM318(3,18), DUM218(2,18)
      REAL(DOUBLE)                    :: GAUSS_R(7), GAUSS_S(7), GAUSS_W(7)
      REAL(DOUBLE)                    :: EPS1

      BIG_BB = ZERO
      IF ((OPT(3) /= 'Y') .AND. (OPT(4) /= 'Y') .AND. (OPT(6) /= 'Y')) RETURN

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

! --- shell_renovation begin --- !
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
! --- shell_renovation end --- !
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
      PHI_SQ = ONE

      CALL MITC3P_GAUSS_7PT(GAUSS_R, GAUSS_S, GAUSS_W)

      KFULL = ZERO
      DO GP=1,7
         CALL MITC3P_BM_AT(JINV, BM)
         CALL MITC3P_BB_AT(GAUSS_R(GP), GAUSS_S(GP), JINV, BB)
         CALL MITC3P_BS_AT(GAUSS_R(GP), GAUSS_S(GP), XY, JINV, JMAT, BS)
         FACTOR = GAUSS_W(GP)*DETJ

         DO I=1,6
            DO J=1,6
               KEI = ZERO
               DO K=1,3
                  DO L=1,3
                     KEI = KEI + BM(K,I)*SHELL_A(K,L)*BM(L,J)
                  ENDDO
               ENDDO
               KFULL(IDX_M(I),IDX_M(J)) = KFULL(IDX_M(I),IDX_M(J)) + FACTOR*KEI
            ENDDO
         ENDDO

         DO I=1,20
            DO J=1,20
               KEI = ZERO
               DO K=1,3
                  DO L=1,3
                     KEI = KEI + MITC3P_BB_TERM(BB,K,I)*SHELL_D(K,L)*MITC3P_BB_TERM(BB,L,J)
                  ENDDO
               ENDDO
               KFULL(I,J) = KFULL(I,J) + FACTOR*KEI
            ENDDO
         ENDDO

         DO I=1,20
            DO J=1,20
               KEI = ZERO
               DO K=1,2
                  DO L=1,2
                     KEI = KEI + MITC3P_BS_TERM(BS,K,I)*COV_S(K,L)*MITC3P_BS_TERM(BS,L,J)
                  ENDDO
               ENDDO
               KFULL(I,J) = KFULL(I,J) + FACTOR*KEI
            ENDDO
         ENDDO
      ENDDO

      KAA = KFULL(1:18,1:18)
      KAB = KFULL(1:18,19:20)
      KBA = KFULL(19:20,1:18)
      KBB = KFULL(19:20,19:20)
      CALL MITC3P_INV2(KBB, KBB_INV)
      KCOND = KAA - MATMUL(KAB, MATMUL(KBB_INV, KBA))

      CALL MITC3P_BB_AT(ONE/THREE, ONE/THREE, JINV, BB)
      CALL MITC3P_BS_AT(ONE/THREE, ONE/THREE, XY, JINV, JMAT, BS)
      CALL MITC3P_RECOVERY_MATS(BB, BS, KBB_INV, KBA, BB_REC, BS_REC)

      KPHYS = KCOND
      DO I=1,3
         KPHYS(6*(I-1)+4,:) = -KPHYS(6*(I-1)+4,:)
         KPHYS(:,6*(I-1)+4) = -KPHYS(:,6*(I-1)+4)
         KPHYS(6*(I-1)+5,:) = -KPHYS(6*(I-1)+5,:)
         KPHYS(:,6*(I-1)+5) = -KPHYS(:,6*(I-1)+5)
      ENDDO

      IF (OPT(4) == 'Y') THEN
         CALL MITC3P_TRANSFORMS(JINV, TAE, TEA)
         KA = MATMUL(TRANSPOSE(TAE), MATMUL(KPHYS, TAE))
         KOUT = MATMUL(TRANSPOSE(TEA), MATMUL(KA, TEA))

         DO I=1,18
            DO J=1,18
               KE(I,J) = KE(I,J) + 0.5D0*(KOUT(I,J) + KOUT(J,I))
            ENDDO
         ENDDO
      ENDIF

      IF ((OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         CALL MITC3P_BUILD_RECOVERY ( BB_REC, BS_REC, BIG_BB )
      ENDIF

      CONTAINS

      SUBROUTINE MITC3P_BUILD_RECOVERY ( BBIN, BSIN, BIG_BB_OUT )
         REAL(DOUBLE), INTENT(IN)     :: BBIN(3,18), BSIN(2,18)
         REAL(DOUBLE), INTENT(INOUT)  :: BIG_BB_OUT(3,ELDOF,1)

         BE2(1:3,1:18,1) = BBIN
         BE3(1:2,1:18,1) = BSIN
         DUM318 = MATMUL(SHELL_D, BBIN)
         DUM218 = MATMUL(SHELL_T, BSIN)
         SE2(1:3,1:18,1) = DUM318
         SE3(1:2,1:18,1) = DUM218
         BIG_BB_OUT(:,:,1) = BBIN

         IF (DEBUG(233) > 0) THEN
            WRITE(F06,'(A,I8)') 'TPLT_MITC3P RECOVERY EID=', EID
            CALL MITC3P_WRITE_MAT(F06, 'BB_REC', BBIN, 3, 18)
            CALL MITC3P_WRITE_MAT(F06, 'BS_REC', BSIN, 2, 18)
         ENDIF
      END SUBROUTINE MITC3P_BUILD_RECOVERY

      SUBROUTINE MITC3P_RECOVERY_MATS(BBIN, BSIN, KBBI, KBAI, BBOUT, BSOUT)
         REAL(DOUBLE), INTENT(IN)  :: BBIN(3,8), BSIN(2,11), KBBI(2,2), KBAI(2,18)
         REAL(DOUBLE), INTENT(OUT) :: BBOUT(3,18), BSOUT(2,18)
         REAL(DOUBLE)              :: BUB_MAP(2,18)
         INTEGER(LONG)             :: IC, IR
         BUB_MAP = -MATMUL(KBBI, KBAI)
         BBOUT = ZERO
         BSOUT = ZERO
         DO IC=1,18
            DO IR=1,3
               BBOUT(IR,IC) = MITC3P_BB_TERM(BBIN,IR,IC) + MITC3P_BB_TERM(BBIN,IR,19)*BUB_MAP(1,IC) + &
                              MITC3P_BB_TERM(BBIN,IR,20)*BUB_MAP(2,IC)
            ENDDO
            DO IR=1,2
               BSOUT(IR,IC) = MITC3P_BS_TERM(BSIN,IR,IC) + MITC3P_BS_TERM(BSIN,IR,19)*BUB_MAP(1,IC) + &
                              MITC3P_BS_TERM(BSIN,IR,20)*BUB_MAP(2,IC)
            ENDDO
         ENDDO
      END SUBROUTINE MITC3P_RECOVERY_MATS

      SUBROUTINE MITC3P_BM_AT(JI, BMOUT)
         REAL(DOUBLE), INTENT(IN)  :: JI(2,2)
         REAL(DOUBLE), INTENT(OUT) :: BMOUT(3,6)
         REAL(DOUBLE) :: DH(3,2), DH_XY(3,2)
         INTEGER(LONG) :: II
         CALL MITC3P_DH(DH)
         DH_XY = MATMUL(DH, TRANSPOSE(JI))
         BMOUT = ZERO
         DO II=1,3
            BMOUT(1,2*II-1) = DH_XY(II,1)
            BMOUT(2,2*II  ) = DH_XY(II,2)
            BMOUT(3,2*II-1) = DH_XY(II,2)
            BMOUT(3,2*II  ) = DH_XY(II,1)
         ENDDO
      END SUBROUTINE MITC3P_BM_AT

      FUNCTION MITC3P_BB_TERM(BBIN, IR, IC) RESULT(VAL)
         REAL(DOUBLE), INTENT(IN) :: BBIN(3,8)
         INTEGER(LONG), INTENT(IN):: IR, IC
         REAL(DOUBLE)             :: VAL
         INTEGER(LONG)            :: INODE, LD
         VAL = ZERO
         IF (IC <= 18) THEN
            INODE = (IC-1)/6 + 1
            LD = IC - 6*(INODE-1)
            IF (LD == 4) VAL = BBIN(IR,2*INODE-1)
            IF (LD == 5) VAL = BBIN(IR,2*INODE  )
         ELSEIF (IC == 19) THEN
            VAL = BBIN(IR,7)
         ELSEIF (IC == 20) THEN
            VAL = BBIN(IR,8)
         ENDIF
      END FUNCTION MITC3P_BB_TERM

      FUNCTION MITC3P_BS_TERM(BSIN, IR, IC) RESULT(VAL)
         REAL(DOUBLE), INTENT(IN) :: BSIN(2,11)
         INTEGER(LONG), INTENT(IN):: IR, IC
         REAL(DOUBLE)             :: VAL
         INTEGER(LONG)            :: INODE, LD
         VAL = ZERO
         IF (IC <= 18) THEN
            INODE = (IC-1)/6 + 1
            LD = IC - 6*(INODE-1)
            IF (LD == 3) VAL = BSIN(IR,3*INODE-2)
            IF (LD == 4) VAL = BSIN(IR,3*INODE-1)
            IF (LD == 5) VAL = BSIN(IR,3*INODE  )
         ELSEIF (IC == 19) THEN
            VAL = BSIN(IR,10)
         ELSEIF (IC == 20) THEN
            VAL = BSIN(IR,11)
         ENDIF
      END FUNCTION MITC3P_BS_TERM

      SUBROUTINE MITC3P_TRANSFORMS(JI, TAE, TEA)
         REAL(DOUBLE), INTENT(IN)  :: JI(2,2)
         REAL(DOUBLE), INTENT(OUT) :: TAE(18,18), TEA(18,18)
         REAL(DOUBLE)              :: DH(3,2), GRAD(3,2), A(3,3,3), A33, M1, M2, A3
         INTEGER(LONG)             :: INODE, JNODE, K3, RO, CO, CL, RW

         CALL MITC3P_DH(DH)
         GRAD = MATMUL(DH, TRANSPOSE(JI))

         TAE = ZERO
         TEA = ZERO

         DO INODE=1,3
            CALL MITC3P_ROTATION_FROM_E3(INODE, A(:,:,INODE))
         ENDDO

         DO INODE=1,3
            RO = 6*(INODE-1)
            DO RW=1,3
               DO CL=1,3
                  TAE(RO+RW,RO+CL) = A(RW,CL,INODE)
                  TEA(RO+RW,RO+CL) = A(CL,RW,INODE)
                  TEA(RO+3+RW,RO+3+CL) = A(CL,RW,INODE)
               ENDDO
            ENDDO

            A33 = A(3,3,INODE)
            IF (DABS(A33) < 1.0D-12) THEN
               IF (A33 < ZERO) THEN
                  A33 = -1.0D-12
               ELSE
                  A33 =  1.0D-12
               ENDIF
            ENDIF

            DO CL=1,2
               DO RW=1,2
                  TAE(RO+3+RW,RO+3+CL) = A(RW,CL,INODE) - (A(RW,3,INODE)*A(CL,3,INODE))/A33
               ENDDO
            ENDDO

            M1 = A(1,3,INODE)/A33
            M2 = A(2,3,INODE)/A33
            DO JNODE=1,3
               CO = 6*(JNODE-1)
               DO K3=1,3
                  A3 = 0.5D0*(A(2,K3,INODE)*GRAD(JNODE,1) - A(1,K3,INODE)*GRAD(JNODE,2))
                  TAE(RO+4,CO+K3) = TAE(RO+4,CO+K3) + M1*A3
                  TAE(RO+5,CO+K3) = TAE(RO+5,CO+K3) + M2*A3
               ENDDO
            ENDDO
         ENDDO
      END SUBROUTINE MITC3P_TRANSFORMS

      SUBROUTINE MITC3P_ROTATION_FROM_E3(INODE, A)
         INTEGER(LONG), INTENT(IN)  :: INODE
         REAL(DOUBLE), INTENT(OUT)  :: A(3,3)
         REAL(DOUBLE)               :: N(3), V(3), VX(3,3), VX2(3,3), C, NN
         REAL(DOUBLE), PARAMETER    :: COS_LIMIT = 0.8191520442889918D0

         N = DIRECTOR(:,INODE)
         NN = DSQRT(DOT_PRODUCT(N,N))
         IF (NN <= EPS1) THEN
            N = (/ ZERO, ZERO, ONE /)
         ELSE
            N = N/NN
         ENDIF

         C = N(3)
         IF (C < COS_LIMIT) THEN
            A = ZERO
            A(1,1) = ONE
            A(2,2) = ONE
            A(3,3) = ONE
            RETURN
         ENDIF

         IF (C > ONE - 1.0D-12) THEN
            A = ZERO
            A(1,1) = ONE
            A(2,2) = ONE
            A(3,3) = ONE
            RETURN
         ENDIF

         IF (C < -ONE + 1.0D-12) THEN
            A = ZERO
            A(1,1) =  ONE
            A(2,2) = -ONE
            A(3,3) = -ONE
            RETURN
         ENDIF

         V = (/ -N(2), N(1), ZERO /)
         VX = ZERO
         VX(1,2) = -V(3)
         VX(1,3) =  V(2)
         VX(2,1) =  V(3)
         VX(2,3) = -V(1)
         VX(3,1) = -V(2)
         VX(3,2) =  V(1)
         VX2 = MATMUL(VX,VX)

         A = ZERO
         A(1,1) = ONE
         A(2,2) = ONE
         A(3,3) = ONE
         A = A + VX + VX2/(ONE + C)
      END SUBROUTINE MITC3P_ROTATION_FROM_E3

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

! --- shell_renovation begin --- !
      SUBROUTINE MITC3P_NODE_BASIS(INODE, REF1, REF2, V1, V2)
         INTEGER(LONG), INTENT(IN)  :: INODE
         REAL(DOUBLE), INTENT(IN)   :: REF1(3), REF2(3)
         REAL(DOUBLE), INTENT(OUT)  :: V1(3), V2(3)
         REAL(DOUBLE)               :: VN(3), REF(3), NORMV

         VN = DIRECTOR(:,INODE)
         NORMV = DSQRT(DOT_PRODUCT(VN,VN))
         IF (NORMV <= EPS1) THEN
            VN = (/ ZERO, ZERO, ONE /)
         ELSE
            VN = VN/NORMV
         ENDIF

         REF = REF1
         V1 = REF - VN*DOT_PRODUCT(REF,VN)
         NORMV = DSQRT(DOT_PRODUCT(V1,V1))
         IF (NORMV <= 1.0D-12) THEN
            REF = REF2
            V1 = REF - VN*DOT_PRODUCT(REF,VN)
            NORMV = DSQRT(DOT_PRODUCT(V1,V1))
         ENDIF
         IF (NORMV <= 1.0D-12) THEN
            REF = (/ ONE, ZERO, ZERO /)
            V1 = REF - VN*DOT_PRODUCT(REF,VN)
            NORMV = DSQRT(DOT_PRODUCT(V1,V1))
         ENDIF
         IF (NORMV <= 1.0D-12) THEN
            REF = (/ ZERO, ONE, ZERO /)
            V1 = REF - VN*DOT_PRODUCT(REF,VN)
            NORMV = DSQRT(DOT_PRODUCT(V1,V1))
         ENDIF
         V1 = V1/NORMV

         CALL CROSS(VN, V1, V2)
         NORMV = DSQRT(DOT_PRODUCT(V2,V2))
         IF (NORMV > EPS1) V2 = V2/NORMV
      END SUBROUTINE MITC3P_NODE_BASIS

      SUBROUTINE MITC3P_EXPAND_BB_DIRECTOR(BBIN, REF1, REF2, BBOUT)
         REAL(DOUBLE), INTENT(IN)  :: BBIN(3,8)
         REAL(DOUBLE), INTENT(IN)  :: REF1(3), REF2(3)
         REAL(DOUBLE), INTENT(OUT) :: BBOUT(3,20)
         REAL(DOUBLE)              :: V1(3), V2(3)
         INTEGER(LONG)             :: INODE, IR, IC

         BBOUT = ZERO
         DO INODE=1,3
            CALL MITC3P_NODE_BASIS(INODE, REF1, REF2, V1, V2)
            DO IR=1,3
               DO IC=1,3
                  BBOUT(IR,6*(INODE-1)+3+IC) = BBOUT(IR,6*(INODE-1)+3+IC) + BBIN(IR,2*INODE-1)*V1(IC) &
                                                                                 + BBIN(IR,2*INODE  )*V2(IC)
               ENDDO
            ENDDO
         ENDDO
         BBOUT(:,19) = BBIN(:,7)
         BBOUT(:,20) = BBIN(:,8)
      END SUBROUTINE MITC3P_EXPAND_BB_DIRECTOR

      SUBROUTINE MITC3P_EXPAND_BS_DIRECTOR(BSIN, REF1, REF2, BSOUT)
         REAL(DOUBLE), INTENT(IN)  :: BSIN(2,11)
         REAL(DOUBLE), INTENT(IN)  :: REF1(3), REF2(3)
         REAL(DOUBLE), INTENT(OUT) :: BSOUT(2,20)
         REAL(DOUBLE)              :: V1(3), V2(3)
         INTEGER(LONG)             :: INODE, IR, IC

         BSOUT = ZERO
         DO INODE=1,3
            CALL MITC3P_NODE_BASIS(INODE, REF1, REF2, V1, V2)
            DO IR=1,2
               BSOUT(IR,6*(INODE-1)+3) = BSOUT(IR,6*(INODE-1)+3) + BSIN(IR,3*INODE-2)
               DO IC=1,3
                  BSOUT(IR,6*(INODE-1)+3+IC) = BSOUT(IR,6*(INODE-1)+3+IC) + BSIN(IR,3*INODE-1)*V1(IC) &
                                                                                 + BSIN(IR,3*INODE  )*V2(IC)
               ENDDO
            ENDDO
         ENDDO
         BSOUT(:,19) = BSIN(:,10)
         BSOUT(:,20) = BSIN(:,11)
      END SUBROUTINE MITC3P_EXPAND_BS_DIRECTOR
! --- shell_renovation end --- !

! --- shell_renovation begin --- !
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
! --- shell_renovation end --- !

      END SUBROUTINE TPLT_MITC3P
! --- shell_renovation end --- !
