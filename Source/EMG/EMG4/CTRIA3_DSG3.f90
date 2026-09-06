! ##################################################################################################################################
! Standalone CTRIA3 DSG3 branch for PARAM,TRIA3TYP,DSG3.
! Direct port of D:\18a_Sept\python\linear\DSG3_Pure5DOF_ShellElement.py.

      SUBROUTINE CTRIA3_DSG3 ( OPT, INT_ELEM_ID )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NSUB, SOL_NAME, NTSUB
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, THREE, FIVE, SIX
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  COUPMASS
      USE MODEL_STUF, ONLY            :  EID, ELGP, KE, KED, ME, BE1, BE2, BE3, SE1, SE2, SE3, STE1, EPROP,                     &
                                         MASS_PER_UNIT_AREA, PRESS, PPE, PTE, TE, NUM_EMG_FATAL_ERRS, SHELL_A, SHELL_D, SHELL_T, &
                                         XEB, UEL, ALPVEC, DT, TREF
      USE ELMDIS_Interface
      USE OUTA_HERE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CTRIA3_DSG3'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      INTEGER(LONG), PARAMETER        :: NNODE = 3
      INTEGER(LONG), PARAMETER        :: NDOFN = 6
      INTEGER(LONG), PARAMETER        :: NDOF  = NNODE*NDOFN

      INTEGER(LONG)                   :: I, J, K, JSUB, IA, IB
      REAL(DOUBLE)                    :: XYZ(NNODE,3), EG(3,3), ECOORDS(NNODE,2), GRADN(NNODE,2), AREA, AE, THICK
      REAL(DOUBLE)                    :: T18(NDOF,NDOF), BM(3,NDOF), BB(3,NDOF), BS(2,NDOF), BSREC(2,NDOF)
      REAL(DOUBLE)                    :: KLOCAL(NDOF,NDOF), KGLOBAL(NDOF,NDOF), MLOCAL(NDOF,NDOF), MGLOBAL(NDOF,NDOF), DRILL
      REAL(DOUBLE)                    :: MASS_NODE, UNIT_PPE(NDOF), UNIT_PTE(NDOF), CTE3(3), NTH(3), TBAR
      REAL(DOUBLE)                    :: EPSM(3), NRES(3), KGVAL

      IF (INT_ELEM_ID < 0) THEN
         CONTINUE
      ENDIF

      IF (ELGP /= 3) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9001) SUBR_NAME, EID, ELGP
         WRITE(F06,9001) SUBR_NAME, EID, ELGP
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      THICK = EPROP(1)
      IF (THICK <= ZERO) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9002) SUBR_NAME, EID, THICK
         WRITE(F06,9002) SUBR_NAME, EID, THICK
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL LOAD_BASIC_COORDS ( XYZ )
      CALL DSG3_GEOMETRY ( XYZ, EG, ECOORDS, GRADN, AE, AREA )
      CALL DSG3_T_MATRIX ( EG, T18 )
      CALL DSG3_BM_BB ( GRADN, BM, BB )
      CALL DSG3_BS_RECOVERY ( ECOORDS, AE, BSREC )

!     DSG3 recovery operators explicitly include the basic-to-local rotation
!     (BM*T18, BB*T18, BSREC*T18). Keep TE as identity so ELMDIS returns basic
!     element displacements and the Python-style recovery remains unchanged.
      TE = ZERO
      TE(1,1) = ONE
      TE(2,2) = ONE
      TE(3,3) = ONE

      IF (OPT(4) == 'Y') THEN
         KLOCAL = ZERO
         KLOCAL = KLOCAL + MATMUL(TRANSPOSE(BM), MATMUL(SHELL_A, BM)) * AREA
         KLOCAL = KLOCAL + MATMUL(TRANSPOSE(BB), MATMUL(SHELL_D, BB)) * AREA

         CALL DSG3_BS_ORDER ( ECOORDS, AE, 1, 2, 3, BS )
         KLOCAL = KLOCAL + MATMUL(TRANSPOSE(BS), MATMUL(SHELL_T, BS)) * (AREA/THREE)
         CALL DSG3_BS_ORDER ( ECOORDS, AE, 2, 3, 1, BS )
         KLOCAL = KLOCAL + MATMUL(TRANSPOSE(BS), MATMUL(SHELL_T, BS)) * (AREA/THREE)
         CALL DSG3_BS_ORDER ( ECOORDS, AE, 3, 1, 2, BS )
         KLOCAL = KLOCAL + MATMUL(TRANSPOSE(BS), MATMUL(SHELL_T, BS)) * (AREA/THREE)

         DRILL = 1.0D-3 * (SIX/FIVE) * MAX(SHELL_T(1,1), SHELL_T(2,2)) * AREA
         DO I=1,NNODE
            KLOCAL(6*(I-1)+6,6*(I-1)+6) = KLOCAL(6*(I-1)+6,6*(I-1)+6) + DRILL
         ENDDO

         KLOCAL  = 0.5D0*(KLOCAL + TRANSPOSE(KLOCAL))
         KGLOBAL = MATMUL(TRANSPOSE(T18), MATMUL(KLOCAL, T18))
         KE(1:NDOF,1:NDOF) = 0.5D0*(KGLOBAL + TRANSPOSE(KGLOBAL))

         IF (DEBUG(233) > 0) THEN
            WRITE(F06,'(A,I8,A,ES15.7)') 'CTRIA3_DSG3 EID=', EID, ' KE_NORM=', DSQRT(SUM(KE(1:NDOF,1:NDOF)*KE(1:NDOF,1:NDOF)))
         ENDIF
      ENDIF

      IF ((OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         DO K=1,ELGP
            BE1(1:3,1:NDOF,K) = MATMUL(BM, T18)
            BE2(1:3,1:NDOF,K) = MATMUL(BB, T18)
            BE3(1:2,1:NDOF,K) = MATMUL(BSREC, T18)
            SE1(1:3,1:NDOF,K) = MATMUL(SHELL_A, BE1(1:3,1:NDOF,K))
            SE2(1:3,1:NDOF,K) = MATMUL(SHELL_D, BE2(1:3,1:NDOF,K))
            SE3(1:2,1:NDOF,K) = MATMUL(SHELL_T, BE3(1:2,1:NDOF,K))
         ENDDO

         CTE3(1) = ALPVEC(1,1)
         CTE3(2) = ALPVEC(2,1)
         CTE3(3) = ALPVEC(4,1)
         NTH = MATMUL(SHELL_A, CTE3)
         DO JSUB=1,NTSUB
            TBAR = (DT(1,JSUB) + DT(2,JSUB) + DT(3,JSUB))/THREE - TREF(1)
            STE1(1:3,JSUB,1:ELGP) = SPREAD(NTH*TBAR,2,ELGP)
         ENDDO
      ENDIF

      IF (OPT(1) == 'Y') THEN
         MLOCAL = ZERO
         IF ((SOL_NAME(1:5) == 'MODES') .AND. (COUPMASS > 0)) THEN
            DO I=1,NNODE
               DO J=1,NNODE
                  IF (I == J) THEN
                     MASS_NODE = MASS_PER_UNIT_AREA * AREA / 6.0D0
                  ELSE
                     MASS_NODE = MASS_PER_UNIT_AREA * AREA / 12.0D0
                  ENDIF
                  DO K=1,3
                     MLOCAL(6*(I-1)+K,6*(J-1)+K) = MASS_NODE
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            MASS_NODE = MASS_PER_UNIT_AREA * AREA / THREE
            DO I=1,NNODE
               DO K=1,3
                  MLOCAL(6*(I-1)+K,6*(I-1)+K) = MASS_NODE
               ENDDO
            ENDDO
         ENDIF
         MGLOBAL = MATMUL(TRANSPOSE(T18), MATMUL(MLOCAL, T18))
         ME(1:NDOF,1:NDOF) = 0.5D0*(MGLOBAL + TRANSPOSE(MGLOBAL))
      ENDIF

      IF (OPT(5) == 'Y') THEN
         UNIT_PPE = ZERO
         DO I=1,NNODE
            UNIT_PPE(6*(I-1)+1) = AREA*EG(1,3)/THREE
            UNIT_PPE(6*(I-1)+2) = AREA*EG(2,3)/THREE
            UNIT_PPE(6*(I-1)+3) = AREA*EG(3,3)/THREE
         ENDDO
         DO JSUB=1,NSUB
            PPE(1:NDOF,JSUB) = PPE(1:NDOF,JSUB) + UNIT_PPE(1:NDOF) * PRESS(3,JSUB)
         ENDDO
      ENDIF

      IF (OPT(2) == 'Y') THEN
         CTE3(1) = ALPVEC(1,1)
         CTE3(2) = ALPVEC(2,1)
         CTE3(3) = ALPVEC(4,1)
         NTH = MATMUL(SHELL_A, CTE3)
         UNIT_PTE = MATMUL(TRANSPOSE(MATMUL(BM,T18)), NTH) * AREA
         DO JSUB=1,NSUB
            TBAR = (DT(1,JSUB) + DT(2,JSUB) + DT(3,JSUB))/THREE - TREF(1)
            PTE(1:NDOF,JSUB) = UNIT_PTE(1:NDOF) * TBAR
         ENDDO
      ENDIF

      IF (OPT(6) == 'Y') THEN
         CALL ELMDIS
         EPSM = MATMUL(MATMUL(BM,T18), UEL(1:NDOF))
         NRES = MATMUL(SHELL_A, EPSM)
         KED(1:NDOF,1:NDOF) = ZERO
         DO IA=1,NNODE
            DO IB=1,NNODE
               KGVAL = AREA*( NRES(1)*GRADN(IA,1)*GRADN(IB,1) + NRES(2)*GRADN(IA,2)*GRADN(IB,2) +                    &
                              NRES(3)*(GRADN(IA,1)*GRADN(IB,2) + GRADN(IA,2)*GRADN(IB,1)) )
               KLOCAL = ZERO
               KLOCAL(6*(IA-1)+1,6*(IB-1)+1) = KGVAL
               KLOCAL(6*(IA-1)+2,6*(IB-1)+2) = KGVAL
               KGLOBAL = MATMUL(TRANSPOSE(T18), MATMUL(KLOCAL, T18))
               KED(1:NDOF,1:NDOF) = KED(1:NDOF,1:NDOF) + KGLOBAL
            ENDDO
         ENDDO
      ENDIF

      RETURN

 9001 FORMAT(' *ERROR: ',A,' expects ELGP=3 for element ',I8,' but got ',I8)
 9002 FORMAT(' *ERROR: ',A,' element ',I8,' has nonpositive thickness ',ES15.7)

      CONTAINS

      SUBROUTINE LOAD_BASIC_COORDS ( XYZOUT )
      REAL(DOUBLE), INTENT(OUT) :: XYZOUT(NNODE,3)
      INTEGER(LONG) :: II, JJ
      DO II=1,NNODE
         DO JJ=1,3
            XYZOUT(II,JJ) = XEB(II,JJ)
         ENDDO
      ENDDO
      END SUBROUTINE LOAD_BASIC_COORDS

      SUBROUTINE DSG3_GEOMETRY ( XYZN, EG, ECOORDS, GRADN, AE, AREA )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(NNODE,3)
      REAL(DOUBLE), INTENT(OUT) :: EG(3,3), ECOORDS(NNODE,2), GRADN(NNODE,2), AE, AREA
      REAL(DOUBLE) :: V12(3), V13(3), VPERP(3), NM, X2, X3, Y3
      V12 = XYZN(2,:) - XYZN(1,:)
      V13 = XYZN(3,:) - XYZN(1,:)
      X2 = VNORM(V12)
      IF (X2 <= 1.0D-14) THEN
         WRITE(ERR,*) ' *ERROR: CTRIA3_DSG3 degenerate edge on element ', EID
         WRITE(F06,*) ' *ERROR: CTRIA3_DSG3 degenerate edge on element ', EID
         CALL OUTA_HERE ( 'Y' )
      ENDIF
      EG(:,1) = V12/X2
      X3 = DOT_PRODUCT(V13, EG(:,1))
      VPERP = V13 - X3*EG(:,1)
      Y3 = VNORM(VPERP)
      IF (Y3 <= 1.0D-14) THEN
         WRITE(ERR,*) ' *ERROR: CTRIA3_DSG3 degenerate area on element ', EID
         WRITE(F06,*) ' *ERROR: CTRIA3_DSG3 degenerate area on element ', EID
         CALL OUTA_HERE ( 'Y' )
      ENDIF
      EG(:,2) = VPERP/Y3
      CALL CROSSV(EG(:,1), EG(:,2), EG(:,3))

      ECOORDS = ZERO
      ECOORDS(2,1) = X2
      ECOORDS(3,1) = X3
      ECOORDS(3,2) = Y3
      AE   = 0.5D0*X2*Y3
      AREA = AE

      GRADN(1,1) = -ONE/X2
      GRADN(2,1) =  ONE/X2
      GRADN(3,1) =  ZERO
      GRADN(1,2) = (X3 - X2)/(X2*Y3)
      GRADN(2,2) = -X3/(X2*Y3)
      GRADN(3,2) =  ONE/Y3
      END SUBROUTINE DSG3_GEOMETRY

      SUBROUTINE DSG3_T_MATRIX ( EG, T18 )
      REAL(DOUBLE), INTENT(IN)  :: EG(3,3)
      REAL(DOUBLE), INTENT(OUT) :: T18(NDOF,NDOF)
      INTEGER(LONG) :: INODE, ROFF
      T18 = ZERO
      DO INODE=1,NNODE
         ROFF = 6*(INODE-1)
         T18(ROFF+1:ROFF+3,ROFF+1:ROFF+3) = TRANSPOSE(EG)
         T18(ROFF+4:ROFF+6,ROFF+4:ROFF+6) = TRANSPOSE(EG)
      ENDDO
      END SUBROUTINE DSG3_T_MATRIX

      SUBROUTINE DSG3_BM_BB ( GRADN, BM, BB )
      REAL(DOUBLE), INTENT(IN)  :: GRADN(NNODE,2)
      REAL(DOUBLE), INTENT(OUT) :: BM(3,NDOF), BB(3,NDOF)
      INTEGER(LONG) :: II, OFF
      BM = ZERO
      BB = ZERO
      DO II=1,NNODE
         OFF = 6*(II-1)
         BM(1,OFF+1) = GRADN(II,1)
         BM(2,OFF+2) = GRADN(II,2)
         BM(3,OFF+1) = GRADN(II,2)
         BM(3,OFF+2) = GRADN(II,1)
         BB(1,OFF+5) =  GRADN(II,1)
         BB(2,OFF+4) = -GRADN(II,2)
         BB(3,OFF+4) = -GRADN(II,1)
         BB(3,OFF+5) =  GRADN(II,2)
      ENDDO
      END SUBROUTINE DSG3_BM_BB

      SUBROUTINE DSG3_BS_RECOVERY ( ECOORDS, AE, BSREC )
      REAL(DOUBLE), INTENT(IN)  :: ECOORDS(NNODE,2), AE
      REAL(DOUBLE), INTENT(OUT) :: BSREC(2,NDOF)
      REAL(DOUBLE) :: BSTMP(2,NDOF)
      BSREC = ZERO
      CALL DSG3_BS_ORDER ( ECOORDS, AE, 1, 2, 3, BSTMP )
      BSREC = BSREC + BSTMP/THREE
      CALL DSG3_BS_ORDER ( ECOORDS, AE, 2, 3, 1, BSTMP )
      BSREC = BSREC + BSTMP/THREE
      CALL DSG3_BS_ORDER ( ECOORDS, AE, 3, 1, 2, BSTMP )
      BSREC = BSREC + BSTMP/THREE
      END SUBROUTINE DSG3_BS_RECOVERY

      SUBROUTINE DSG3_BS_ORDER ( ECOORDS, AE, IS, IP, IQ, BS )
      REAL(DOUBLE), INTENT(IN)  :: ECOORDS(NNODE,2), AE
      INTEGER(LONG), INTENT(IN) :: IS, IP, IQ
      REAL(DOUBLE), INTENT(OUT) :: BS(2,NDOF)
      REAL(DOUBLE) :: XS, YS, XP, YP, XQ, YQ, A, B, C, D, M
      BS = ZERO
      XS = ECOORDS(IS,1); YS = ECOORDS(IS,2)
      XP = ECOORDS(IP,1); YP = ECOORDS(IP,2)
      XQ = ECOORDS(IQ,1); YQ = ECOORDS(IQ,2)
      A = XP - XS
      B = YP - YS
      C = XQ - XS
      D = YQ - YS
      M = ONE/(TWO*AE)
      CALL ADD_BS_NODE_S ( BS, IS, M, A, B, C, D, AE )
      CALL ADD_BS_NODE_P ( BS, IP, M, A, B, C, D )
      CALL ADD_BS_NODE_Q ( BS, IQ, M, A, B, C, D )
      END SUBROUTINE DSG3_BS_ORDER

      SUBROUTINE ADD_BS_NODE_S ( BS, INODE, M, A, B, C, D, AE )
      REAL(DOUBLE), INTENT(INOUT) :: BS(2,NDOF)
      INTEGER(LONG), INTENT(IN) :: INODE
      REAL(DOUBLE), INTENT(IN) :: M, A, B, C, D, AE
      INTEGER(LONG) :: OFF
      OFF = 6*(INODE-1)
      BS(1,OFF+3) = BS(1,OFF+3) + M*(B - D)
      BS(1,OFF+5) = BS(1,OFF+5) + M*AE
      BS(2,OFF+3) = BS(2,OFF+3) + M*(C - A)
      BS(2,OFF+4) = BS(2,OFF+4) - M*AE
      END SUBROUTINE ADD_BS_NODE_S

      SUBROUTINE ADD_BS_NODE_P ( BS, INODE, M, A, B, C, D )
      REAL(DOUBLE), INTENT(INOUT) :: BS(2,NDOF)
      INTEGER(LONG), INTENT(IN) :: INODE
      REAL(DOUBLE), INTENT(IN) :: M, A, B, C, D
      INTEGER(LONG) :: OFF
      OFF = 6*(INODE-1)
      BS(1,OFF+3) = BS(1,OFF+3) + M*D
      BS(1,OFF+4) = BS(1,OFF+4) - M*B*D/TWO
      BS(1,OFF+5) = BS(1,OFF+5) + M*A*D/TWO
      BS(2,OFF+3) = BS(2,OFF+3) - M*C
      BS(2,OFF+4) = BS(2,OFF+4) + M*B*C/TWO
      BS(2,OFF+5) = BS(2,OFF+5) - M*A*C/TWO
      END SUBROUTINE ADD_BS_NODE_P

      SUBROUTINE ADD_BS_NODE_Q ( BS, INODE, M, A, B, C, D )
      REAL(DOUBLE), INTENT(INOUT) :: BS(2,NDOF)
      INTEGER(LONG), INTENT(IN) :: INODE
      REAL(DOUBLE), INTENT(IN) :: M, A, B, C, D
      INTEGER(LONG) :: OFF
      OFF = 6*(INODE-1)
      BS(1,OFF+3) = BS(1,OFF+3) - M*B
      BS(1,OFF+4) = BS(1,OFF+4) + M*B*D/TWO
      BS(1,OFF+5) = BS(1,OFF+5) - M*B*C/TWO
      BS(2,OFF+3) = BS(2,OFF+3) + M*A
      BS(2,OFF+4) = BS(2,OFF+4) - M*A*D/TWO
      BS(2,OFF+5) = BS(2,OFF+5) + M*A*C/TWO
      END SUBROUTINE ADD_BS_NODE_Q

      SUBROUTINE CROSSV(A,B,C)
      REAL(DOUBLE), INTENT(IN) :: A(3), B(3)
      REAL(DOUBLE), INTENT(OUT) :: C(3)
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      END SUBROUTINE CROSSV

      FUNCTION VNORM ( V ) RESULT(NM)
      REAL(DOUBLE), INTENT(IN) :: V(3)
      REAL(DOUBLE) :: NM
      NM = DSQRT(DOT_PRODUCT(V,V))
      END FUNCTION VNORM

      END SUBROUTINE CTRIA3_DSG3
