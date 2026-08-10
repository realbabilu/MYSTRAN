! ##################################################################################################################################
! CTRIA3 T3FF selector for PARAM,TRIA3TYP,T3FF.

      SUBROUTINE CTRIA3_T3FF ( OPT, INT_ELEM_ID )

! Port target: D:\18a\bending_only\Shell\gemini2\shit\T3FF_ShellElement_Std.py
!              D:\18a\bending_only\Shell\gemini2\shit\T3FF_ShellElement_fixed.py
! Generated SNORM preprocessing is handled in LINK0 before this routine is reached.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NSUB, SOL_NAME
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, THREE, TWELVE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  COUPMASS, TRIARTYP
      USE MODEL_STUF, ONLY            :  EID, ELGP, KE, KED, ME, BE1, BE2, BE3, EPROP, MASS_PER_UNIT_AREA, PRESS, PPE,             &
                                         TE, NUM_EMG_FATAL_ERRS, SHELL_A, SHELL_D, SHELL_T, BGRID, GRID_SNORM, XEB, UEL
      USE ELMDIS_Interface
      USE OUTA_HERE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CTRIA3_T3FF'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      INTEGER(LONG), PARAMETER        :: NNODE = 3
      INTEGER(LONG), PARAMETER        :: NDOFN = 6
      INTEGER(LONG), PARAMETER        :: NDOF  = NNODE*NDOFN

      INTEGER(LONG)                   :: I, J, JSUB, IA, IB, RR
      REAL(DOUBLE)                    :: XYZ(NNODE,3), NORMALS(NNODE,3), EG(3,3)
      REAL(DOUBLE)                    :: ECOORDS(NNODE,2), GRADN(NNODE,2), AREA, AE, HCHAR, THICK, STAB
      REAL(DOUBLE)                    :: TGA(NDOF,NDOF), TAE(NDOF,NDOF), TALL(NDOF,NDOF)
      REAL(DOUBLE)                    :: BME(3,NDOF), BBE(3,NDOF), BSE(2,NDOF), BSG(2,NDOF)
      REAL(DOUBLE)                    :: KE_E(NDOF,NDOF), KE_A(NDOF,NDOF), KGLOBAL(NDOF,NDOF)
      REAL(DOUBLE)                    :: KAVG, MASS_NODE, UNIT_PPE(NDOF), NVEC(3)
      REAL(DOUBLE)                    :: BMG(3,NDOF), EPSM(3), NRES(3), KGVAL
      LOGICAL                         :: T3FFD_RECOVERY

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
      CALL T3FF_GEOMETRY ( XYZ, EG, ECOORDS, GRADN, AE, AREA, HCHAR )
      CALL T3FF_NODAL_NORMALS ( EG, NORMALS )
      CALL T3FF_T_MATRICES ( NORMALS, GRADN, EG, TGA, TAE )
      TALL = MATMUL(TAE, TGA)
      T3FFD_RECOVERY = (TRIARTYP == 'T3FFD   ')

      TE = ZERO
      TE(1,1) = ONE
      TE(2,2) = ONE
      TE(3,3) = ONE

      STAB = THICK*THICK / (THICK*THICK + (5.0D0/18.0D0)*HCHAR*HCHAR)

      IF ((OPT(3) == 'Y') .OR. (OPT(4) == 'Y')) THEN
         CALL T3FF_BM_BB_E ( GRADN, BME, BBE )

         KE_E = MATMUL(TRANSPOSE(BME), MATMUL(SHELL_A, BME)) * AREA
         KE_E = KE_E + MATMUL(TRANSPOSE(BBE), MATMUL(SHELL_D, BBE)) * AREA

         CALL T3FF_BS_ORDER ( ECOORDS, AE, 1, 2, 3, BSE )
         KE_E = KE_E + MATMUL(TRANSPOSE(BSE), MATMUL(SHELL_T, BSE)) * (STAB*AREA/THREE)
         CALL T3FF_BS_ORDER ( ECOORDS, AE, 2, 3, 1, BSE )
         KE_E = KE_E + MATMUL(TRANSPOSE(BSE), MATMUL(SHELL_T, BSE)) * (STAB*AREA/THREE)
         CALL T3FF_BS_ORDER ( ECOORDS, AE, 3, 1, 2, BSE )
         KE_E = KE_E + MATMUL(TRANSPOSE(BSE), MATMUL(SHELL_T, BSE)) * (STAB*AREA/THREE)

         KE_E = 0.5D0*(KE_E + TRANSPOSE(KE_E))
         KE_A = MATMUL(TRANSPOSE(TAE), MATMUL(KE_E, TAE))

         KAVG = (KE_A(4,4) + KE_A(10,10) + KE_A(16,16) + KE_A(5,5) + KE_A(11,11) + KE_A(17,17))/6.0D0
         IF ((KAVG > ZERO) .AND. (KAVG == KAVG)) THEN
            KE_A(6,6)   = KE_A(6,6)   + KAVG
            KE_A(12,12) = KE_A(12,12) + KAVG
            KE_A(18,18) = KE_A(18,18) + KAVG
         ENDIF

         KGLOBAL = MATMUL(TRANSPOSE(TGA), MATMUL(KE_A, TGA))
         KGLOBAL = 0.5D0*(KGLOBAL + TRANSPOSE(KGLOBAL))

         IF (OPT(4) == 'Y') THEN
            KE(1:NDOF,1:NDOF) = KGLOBAL
         ENDIF

         IF (OPT(3) == 'Y') THEN
            BE1(1:3,1:NDOF,1) = MATMUL(BME, TALL)
            BE2(1:3,1:NDOF,1) = MATMUL(BBE, TALL)
            IF (T3FFD_RECOVERY) THEN
               CALL T3FF_BS_ORDER ( ECOORDS, AE, 1, 2, 3, BSE )
               BE3(1:2,1:NDOF,1) = MATMUL(BSE, TALL)
            ELSE
               BSG = ZERO
               CALL T3FF_BS_ORDER ( ECOORDS, AE, 1, 2, 3, BSE )
               BSG = BSG + MATMUL(BSE, TALL)
               CALL T3FF_BS_ORDER ( ECOORDS, AE, 2, 3, 1, BSE )
               BSG = BSG + MATMUL(BSE, TALL)
               CALL T3FF_BS_ORDER ( ECOORDS, AE, 3, 1, 2, BSE )
               BSG = BSG + MATMUL(BSE, TALL)
               BE3(1:2,1:NDOF,1) = BSG / THREE
            ENDIF
         ENDIF

         IF ((DEBUG(233) > 0) .AND. (OPT(4) == 'Y')) THEN
            WRITE(F06,'(A,I8,A,ES15.7)') 'CTRIA3_T3FF EID=', EID, ' KE_NORM=', DSQRT(SUM(KE(1:NDOF,1:NDOF)*KE(1:NDOF,1:NDOF)))
         ENDIF
      ENDIF

      IF (OPT(1) == 'Y') THEN
         ME(1:NDOF,1:NDOF) = ZERO
         IF ((SOL_NAME(1:5) == 'MODES') .AND. (COUPMASS > 0)) THEN
            DO I=1,NNODE
               DO J=1,NNODE
                  IF (I == J) THEN
                     MASS_NODE = MASS_PER_UNIT_AREA * AREA / 6.0D0
                  ELSE
                     MASS_NODE = MASS_PER_UNIT_AREA * AREA / 12.0D0
                  ENDIF
                  DO RR=1,3
                     ME(6*(I-1)+RR,6*(J-1)+RR) = MASS_NODE
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            MASS_NODE = MASS_PER_UNIT_AREA * AREA / THREE
            DO I=1,NNODE
               DO J=1,3
                  ME(6*(I-1)+J,6*(I-1)+J) = MASS_NODE
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      IF (OPT(5) == 'Y') THEN
         NVEC = EG(:,3)
         UNIT_PPE = ZERO
         DO I=1,NNODE
            UNIT_PPE(6*(I-1)+1) = AREA*NVEC(1)/THREE
            UNIT_PPE(6*(I-1)+2) = AREA*NVEC(2)/THREE
            UNIT_PPE(6*(I-1)+3) = AREA*NVEC(3)/THREE
         ENDDO
         DO JSUB=1,NSUB
            PPE(1:NDOF,JSUB) = PPE(1:NDOF,JSUB) + UNIT_PPE(1:NDOF) * PRESS(3,JSUB)
         ENDDO
      ENDIF

      IF (OPT(6) == 'Y') THEN
!        Buckling path from D:\18a\bending_only\Shell\gemini2\shit\buckling\T3FF_ShellElement_buckling.py.
!        One-point triangular geometric stiffness based on the current membrane resultant.
         CALL T3FF_BM_BB_E ( GRADN, BME, BBE )
         BMG = MATMUL(BME, TALL)
         CALL ELMDIS
         EPSM = MATMUL(BMG, UEL(1:NDOF))
         NRES = MATMUL(SHELL_A, EPSM)
         KED(1:NDOF,1:NDOF) = ZERO
         DO IA=1,NNODE
            DO IB=1,NNODE
               KGVAL = AREA*( NRES(1)*GRADN(IA,1)*GRADN(IB,1) + NRES(2)*GRADN(IA,2)*GRADN(IB,2) +                    &
                              NRES(3)*(GRADN(IA,1)*GRADN(IB,2) + GRADN(IA,2)*GRADN(IB,1)) )
               DO RR=1,3
                  KED(6*(IA-1)+RR,6*(IB-1)+RR) = KED(6*(IA-1)+RR,6*(IB-1)+RR) + KGVAL
               ENDDO
            ENDDO
         ENDDO
         IF ((DEBUG(233) > 0) .AND. (EID <= 8)) THEN
            WRITE(F06,'(A,I8,A,3(1X,ES15.7))') 'CTRIA3_T3FF KGGD EID=', EID, ' NRES=', NRES(1), NRES(2), NRES(3)
            WRITE(F06,'(A,I8,A,ES15.7)') 'CTRIA3_T3FF KGGD EID=', EID, ' KED_NORM=', DSQRT(SUM(KED(1:NDOF,1:NDOF)*KED(1:NDOF,1:NDOF)))
         ENDIF
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

      SUBROUTINE T3FF_GEOMETRY ( XYZN, EG, ECOORDS, GRADN, AE, AREA, HCHAR )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(NNODE,3)
      REAL(DOUBLE), INTENT(OUT) :: EG(3,3), ECOORDS(NNODE,2), GRADN(NNODE,2), AE, AREA, HCHAR
      REAL(DOUBLE) :: G1(3), G2(3), E3(3), NM, A, B, C, D, DETJ
      G1 = XYZN(2,:) - XYZN(1,:)
      G2 = XYZN(3,:) - XYZN(1,:)
      NM = VNORM(G1)
      IF (NM <= 1.0D-14) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,*) ' *ERROR: CTRIA3_T3FF degenerate edge on element ', EID
         WRITE(F06,*) ' *ERROR: CTRIA3_T3FF degenerate edge on element ', EID
         CALL OUTA_HERE ( 'Y' )
      ENDIF
      EG(:,1) = G1/NM
      CALL CROSSV(EG(:,1), G2, E3)
      NM = VNORM(E3)
      IF (NM <= 1.0D-14) THEN
         CALL CROSSV(G1, G2, E3)
         NM = VNORM(E3)
      ENDIF
      IF (NM <= 1.0D-14) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,*) ' *ERROR: CTRIA3_T3FF degenerate area on element ', EID
         WRITE(F06,*) ' *ERROR: CTRIA3_T3FF degenerate area on element ', EID
         CALL OUTA_HERE ( 'Y' )
      ENDIF
      EG(:,3) = E3/NM
      CALL CROSSV(EG(:,3), EG(:,1), EG(:,2))

      ECOORDS = ZERO
      ECOORDS(2,1) = DOT_PRODUCT(G1, EG(:,1))
      ECOORDS(2,2) = DOT_PRODUCT(G1, EG(:,2))
      ECOORDS(3,1) = DOT_PRODUCT(G2, EG(:,1))
      ECOORDS(3,2) = DOT_PRODUCT(G2, EG(:,2))

      A = ECOORDS(2,1) - ECOORDS(1,1)
      B = ECOORDS(2,2) - ECOORDS(1,2)
      C = ECOORDS(3,1) - ECOORDS(1,1)
      D = ECOORDS(3,2) - ECOORDS(1,2)
      DETJ = A*D - B*C
      IF (DABS(DETJ) <= 1.0D-14) DETJ = SIGN(1.0D-14, DETJ + 1.0D-30)
      GRADN(1,1) = (B - D)/DETJ
      GRADN(2,1) = D/DETJ
      GRADN(3,1) = -B/DETJ
      GRADN(1,2) = (C - A)/DETJ
      GRADN(2,2) = -C/DETJ
      GRADN(3,2) = A/DETJ
      AE = 0.5D0*DETJ
      AREA = DABS(AE)
      HCHAR = DSQRT(TWO*AREA)
      END SUBROUTINE T3FF_GEOMETRY

      SUBROUTINE T3FF_NODAL_NORMALS ( EG, NORMALS )
      REAL(DOUBLE), INTENT(IN)  :: EG(3,3)
      REAL(DOUBLE), INTENT(OUT) :: NORMALS(NNODE,3)
      REAL(DOUBLE) :: SN(3), NM
      INTEGER(LONG) :: II
      DO II=1,NNODE
         NORMALS(II,:) = EG(:,3)
      ENDDO
      IF (ALLOCATED(GRID_SNORM)) THEN
         DO II=1,NNODE
            IF ((BGRID(II) > 0) .AND. (BGRID(II) <= SIZE(GRID_SNORM,1))) THEN
               SN = GRID_SNORM(BGRID(II),:)
               NM = VNORM(SN)
               IF (NM > 1.0D-15) NORMALS(II,:) = SN/NM
            ENDIF
         ENDDO
      ENDIF
      END SUBROUTINE T3FF_NODAL_NORMALS

      SUBROUTINE T3FF_BM_BB_E ( GRADN, BM, BB )
      REAL(DOUBLE), INTENT(IN)  :: GRADN(NNODE,2)
      REAL(DOUBLE), INTENT(OUT) :: BM(3,NDOF), BB(3,NDOF)
      INTEGER(LONG) :: II, OFF
      BM = ZERO
      BB = ZERO
      DO II=1,NNODE
         OFF = (II-1)*6
         BM(1,OFF+1) = GRADN(II,1)
         BM(2,OFF+2) = GRADN(II,2)
         BM(3,OFF+1) = GRADN(II,2)
         BM(3,OFF+2) = GRADN(II,1)
         BB(1,OFF+5) =  GRADN(II,1)
         BB(2,OFF+4) = -GRADN(II,2)
         BB(3,OFF+4) = -GRADN(II,1)
         BB(3,OFF+5) =  GRADN(II,2)
      ENDDO
      END SUBROUTINE T3FF_BM_BB_E

      SUBROUTINE T3FF_BS_ORDER ( ECOORDS, AE, IS, IP, IQ, BS )
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
      IF (DABS(AE) <= 1.0D-14) THEN
         M = ZERO
      ELSE
         M = ONE/(TWO*AE)
      ENDIF
      CALL ADD_BS_NODE_S ( BS, IS, M, A, B, C, D, AE )
      CALL ADD_BS_NODE_P ( BS, IP, M, A, B, C, D )
      CALL ADD_BS_NODE_Q ( BS, IQ, M, A, B, C, D )
      END SUBROUTINE T3FF_BS_ORDER

      SUBROUTINE ADD_BS_NODE_S ( BS, INODE, M, A, B, C, D, AE )
      REAL(DOUBLE), INTENT(INOUT) :: BS(2,NDOF)
      INTEGER(LONG), INTENT(IN) :: INODE
      REAL(DOUBLE), INTENT(IN) :: M, A, B, C, D, AE
      INTEGER(LONG) :: OFF
      OFF = (INODE-1)*6
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
      OFF = (INODE-1)*6
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
      OFF = (INODE-1)*6
      BS(1,OFF+3) = BS(1,OFF+3) - M*B
      BS(1,OFF+4) = BS(1,OFF+4) + M*B*D/TWO
      BS(1,OFF+5) = BS(1,OFF+5) - M*B*C/TWO
      BS(2,OFF+3) = BS(2,OFF+3) + M*A
      BS(2,OFF+4) = BS(2,OFF+4) - M*A*D/TWO
      BS(2,OFF+5) = BS(2,OFF+5) + M*A*C/TWO
      END SUBROUTINE ADD_BS_NODE_Q

      SUBROUTINE T3FF_T_MATRICES ( NORMS, GRADN, EG, TGA, TAE )
      REAL(DOUBLE), INTENT(IN)  :: NORMS(NNODE,3), GRADN(NNODE,2), EG(3,3)
      REAL(DOUBLE), INTENT(OUT) :: TGA(NDOF,NDOF), TAE(NDOF,NDOF)
      REAL(DOUBLE) :: AROT(3,3), NK(3), TBLOCK(3,3), A33, M1, M2, A3
      INTEGER(LONG) :: INODE, JNODE, R, C, ROFF, COFF
      TGA = ZERO
      TAE = ZERO
      DO INODE=1,NNODE
         NK = MATMUL(TRANSPOSE(EG), NORMS(INODE,:))
         IF (VNORM(NK) <= 1.0D-12) NK = (/ZERO,ZERO,ONE/)
         NK = NK / VNORM(NK)
         CALL ROT_F3_TO_N(NK, AROT)
         TBLOCK = MATMUL(TRANSPOSE(AROT), TRANSPOSE(EG))
         ROFF = (INODE-1)*6
         TGA(ROFF+1:ROFF+3,ROFF+1:ROFF+3) = TBLOCK
         TGA(ROFF+4:ROFF+6,ROFF+4:ROFF+6) = TBLOCK
         TAE(ROFF+1:ROFF+3,ROFF+1:ROFF+3) = AROT
         A33 = AROT(3,3)
         DO R=1,2
            DO C=1,2
               IF (DABS(A33) > 1.0D-12) THEN
                  TAE(ROFF+3+R,ROFF+3+C) = AROT(R,C) - (AROT(R,3)*AROT(C,3))/A33
               ELSE
                  TAE(ROFF+3+R,ROFF+3+C) = AROT(R,C)
               ENDIF
            ENDDO
         ENDDO
         IF (DABS(A33) > 1.0D-12) THEN
            M1 = AROT(1,3)/A33
            M2 = AROT(2,3)/A33
            DO JNODE=1,NNODE
               COFF = (JNODE-1)*6
               DO C=1,3
                  A3 = 0.5D0*(AROT(2,C)*GRADN(JNODE,1) - AROT(1,C)*GRADN(JNODE,2))
                  TAE(ROFF+4,COFF+C) = TAE(ROFF+4,COFF+C) + M1*A3
                  TAE(ROFF+5,COFF+C) = TAE(ROFF+5,COFF+C) + M2*A3
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      END SUBROUTINE T3FF_T_MATRICES

      SUBROUTINE ROT_F3_TO_N ( N, R )
      REAL(DOUBLE), INTENT(IN) :: N(3)
      REAL(DOUBLE), INTENT(OUT) :: R(3,3)
      REAL(DOUBLE) :: A(3), V(3), KX(3,3), C, NM, AXIS(3)
      A = (/ZERO,ZERO,ONE/)
      CALL CROSSV(A,N,V)
      C = DOT_PRODUCT(A,N)
      NM = VNORM(V)
      IF (NM <= 1.0D-12) THEN
         IF (C > ZERO) THEN
            R = ID3()
         ELSE
            AXIS = (/ONE,ZERO,ZERO/)
            CALL SKEW(AXIS,KX)
            R = ID3() + TWO*MATMUL(KX,KX)
         ENDIF
         RETURN
      ENDIF
      CALL SKEW(V,KX)
      R = ID3() + KX + MATMUL(KX,KX)/(ONE + C)
      END SUBROUTINE ROT_F3_TO_N

      FUNCTION ID3() RESULT(A)
      REAL(DOUBLE) :: A(3,3)
      A=ZERO; A(1,1)=ONE; A(2,2)=ONE; A(3,3)=ONE
      END FUNCTION ID3

      SUBROUTINE SKEW(V,A)
      REAL(DOUBLE), INTENT(IN) :: V(3)
      REAL(DOUBLE), INTENT(OUT) :: A(3,3)
      A=ZERO
      A(1,2)=-V(3); A(1,3)= V(2)
      A(2,1)= V(3); A(2,3)=-V(1)
      A(3,1)=-V(2); A(3,2)= V(1)
      END SUBROUTINE SKEW

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

      END SUBROUTINE CTRIA3_T3FF
