! #################################################################################################################################
! CQUADR Q4RS shell for PARAM,QUADRTYP,Q4RS.

      SUBROUTINE CQUADR_Q4RS ( OPT, INT_ELEM_ID )

! Port seed from D:\18a\bending_only\Shell\gemini2\shit\Q4RS_ShellElement.py and
! D:\18a\bending_only\Shell\gemini2\shit\buckling\Q4RS_ShellElement_buckling.py.
! This file is intentionally independent from the Simo/DKMQ kernels so Q4RS can
! be swapped and validated without editing the CQUADR dispatcher.  OPT(6) follows
! the Python k_geom_local() path directly from UEL instead of stress recovery.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_GAUSS, NSUB, SOL_NAME
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, FOUR
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  COUPMASS
      USE MODEL_STUF, ONLY            :  EID, ELGP, KE, ME, BE1, BE2, BE3, EPROP, MASS_PER_UNIT_AREA, PRESS, PPE, PTE,            &
                                         TE, NUM_EMG_FATAL_ERRS, SHELL_A, SHELL_D, SHELL_T, BGRID, GRID_SNORM, XEB,                &
                                         KED, UEL, ALPVEC, DT, TREF
      USE CQUADR_DKMQ24R_Interface
      USE ELMDIS_Interface
      USE ORDER_GAUSS_Interface
      USE OUTA_HERE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CQUADR_Q4RS'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID
      CHARACTER(1*BYTE)                :: REC_OPT(6)

      INTEGER(LONG)                   :: I,J,K,GP,JSUB,IA,IB,RR
      REAL(DOUBLE)                    :: XYZ(4,3), NORMALS(4,3), T24(24,24), T24T(24,24)
      REAL(DOUBLE)                    :: SS(MAX_ORDER_GAUSS), HH(MAX_ORDER_GAUSS), XI, ETA, WT
      REAL(DOUBLE)                    :: REC_XI(5), REC_ETA(5)
      REAL(DOUBLE)                    :: BMB(3,24), BBB(3,24), BSB(2,24), BML(3,24), BBL(3,24), BSL(2,24)
      REAL(DOUBLE)                    :: GBE1(3,24,4), GBE2(3,24,4), GBE3(2,24,4)
      REAL(DOUBLE)                    :: KBASIC(24,24), KLOCAL(24,24), KDRILL(24,24), MLOCAL(24,24), MBASIC(24,24)
      REAL(DOUBLE)                    :: KGLOCAL(24,24), KGVAL, EPSM(3), NRES(3), UE_BASIC(24)
      REAL(DOUBLE)                    :: JAC, STAB, THICK, DIAM, NVG(4), MASS_DIAG(4), M1(4,4)
      REAL(DOUBLE)                    :: GRADN(4,2), ECOORDS(4,2)
      REAL(DOUBLE)                    :: EG(3,3), UNIT_PPE_B(24), UNIT_PPE_L(24), UNIT_PTE_B(24)
      REAL(DOUBLE)                    :: CTE3(3), NTH(3), TBAR

      IF (ELGP /= 4) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9001) SUBR_NAME, EID, ELGP
         WRITE(F06,9001) SUBR_NAME, EID, ELGP
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL LOAD_BASIC_COORDS ( XYZ )
      CALL CALC_Q4RS_NODAL_NORMALS ( XYZ, NORMALS )
      CALL BUILD_T24 ( TE, T24 )
      T24T = TRANSPOSE(T24)

      THICK = EPROP(1)
      IF (THICK <= ZERO) THICK = ONE
      DIAM = QUAD_DIAMETER(XYZ)
      STAB = (THICK*THICK)/(THICK*THICK + 0.2D0*DIAM*DIAM)

      IF ((OPT(3) == 'Y') .OR. (OPT(4) == 'Y') .OR. (OPT(1) == 'Y') .OR. (OPT(2) == 'Y') .OR. (OPT(5) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         CALL ORDER_GAUSS(2, SS, HH)
      ENDIF

      IF ((OPT(3) == 'Y') .OR. (OPT(4) == 'Y')) THEN
         KBASIC = ZERO

         DO I=1,2
            DO J=1,2
               XI = SS(I)
               ETA = SS(J)
               WT = HH(I)*HH(J)

               CALL Q4RS_B_MATRICES ( XYZ, NORMALS, XI, ETA, BMB, BBB, BSB, EG, JAC )
               KBASIC = KBASIC + WT*JAC*MATMUL(TRANSPOSE(BMB), MATMUL(SHELL_A, BMB))
               KBASIC = KBASIC + WT*JAC*MATMUL(TRANSPOSE(BBB), MATMUL(SHELL_D, BBB))
               KBASIC = KBASIC + WT*JAC*STAB*MATMUL(TRANSPOSE(BSB), MATMUL(SHELL_T, BSB))

               GP = GP_INDEX(I,J)
               BML = MATMUL(BMB, T24T)
               BBL = MATMUL(BBB, T24T)
               BSL = MATMUL(BSB, T24T)
               GBE1(1:3,1:24,GP) = BML
               GBE2(1:3,1:24,GP) = BBL
               GBE3(1:2,1:24,GP) = BSL
            ENDDO
         ENDDO

         KDRILL = Q4RS_DRILL_STIFFNESS(KBASIC, NORMALS)
         KBASIC = 0.5D0*(KBASIC + TRANSPOSE(KBASIC)) + KDRILL
         KLOCAL = MATMUL(T24, MATMUL(KBASIC, T24T))

         IF (OPT(4) == 'Y') THEN
            KE(1:24,1:24) = KLOCAL
         ENDIF

         IF ((DEBUG(190) > 0) .AND. (OPT(4) == 'Y')) THEN
            WRITE(F06,'(A,I8,A,ES15.7)') 'CQUADR_Q4RS EID=', EID, ' KE_NORM=', DSQRT(SUM(KE(1:24,1:24)*KE(1:24,1:24)))
         ENDIF
      ENDIF

      IF (OPT(1) == 'Y') THEN
         MBASIC = ZERO
         M1 = ZERO
         MASS_DIAG = ZERO
         DO I=1,2
            DO J=1,2
               XI = SS(I)
               ETA = SS(J)
               WT = HH(I)*HH(J)
               CALL Q4RS_GEOM ( XYZ, XI, ETA, EG, JAC )
               CALL SHAPE_N ( XI, ETA, NVG )
               DO GP=1,4
                  MASS_DIAG(GP) = MASS_DIAG(GP) + NVG(GP)*MASS_PER_UNIT_AREA*WT*JAC
                  DO K=1,4
                     M1(GP,K) = M1(GP,K) + NVG(GP)*NVG(K)*MASS_PER_UNIT_AREA*WT*JAC
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         IF ((SOL_NAME(1:5) == 'MODES') .AND. (COUPMASS > 0)) THEN
            DO I=1,4
               DO J=1,4
                  DO K=1,3
                     MBASIC((I-1)*6+K,(J-1)*6+K) = M1(I,J)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO I=1,4
               DO K=1,3
                  MBASIC((I-1)*6+K,(I-1)*6+K) = MASS_DIAG(I)
               ENDDO
            ENDDO
         ENDIF
         MLOCAL = MATMUL(T24, MATMUL(MBASIC, T24T))
         ME(1:24,1:24) = MLOCAL
      ENDIF

      IF (OPT(5) == 'Y') THEN
         UNIT_PPE_B = ZERO
         DO I=1,2
            DO J=1,2
               XI = SS(I)
               ETA = SS(J)
               WT = HH(I)*HH(J)
               CALL Q4RS_GEOM ( XYZ, XI, ETA, EG, JAC )
               CALL SHAPE_N ( XI, ETA, NVG )
               DO GP=1,4
                  UNIT_PPE_B((GP-1)*6+1:(GP-1)*6+3) = UNIT_PPE_B((GP-1)*6+1:(GP-1)*6+3) + EG(:,3)*NVG(GP)*WT*JAC
               ENDDO
            ENDDO
         ENDDO
         UNIT_PPE_L = MATMUL(T24, UNIT_PPE_B)
         DO JSUB=1,NSUB
            PPE(1:24,JSUB) = UNIT_PPE_L(1:24) * PRESS(3,JSUB)
         ENDDO
      ENDIF

      IF (OPT(2) == 'Y') THEN
         UNIT_PTE_B = ZERO
         CTE3(1) = ALPVEC(1,1)
         CTE3(2) = ALPVEC(2,1)
         CTE3(3) = ALPVEC(4,1)
         NTH = MATMUL(SHELL_A, CTE3)
         DO I=1,2
            DO J=1,2
               XI = SS(I)
               ETA = SS(J)
               WT = HH(I)*HH(J)
               CALL Q4RS_B_MATRICES ( XYZ, NORMALS, XI, ETA, BMB, BBB, BSB, EG, JAC )
               UNIT_PTE_B = UNIT_PTE_B + WT*JAC*MATMUL(TRANSPOSE(BMB), NTH)
            ENDDO
         ENDDO
         DO JSUB=1,NSUB
            TBAR = (DT(1,JSUB) + DT(2,JSUB) + DT(3,JSUB) + DT(4,JSUB))/FOUR - TREF(1)
            PTE(1:24,JSUB) = MATMUL(T24, UNIT_PTE_B) * TBAR
         ENDDO
      ENDIF

      IF (OPT(6) == 'Y') THEN
         CALL ELMDIS
         UE_BASIC = MATMUL(T24T, UEL(1:24))

         KGLOCAL = ZERO
         DO I=1,2
            DO J=1,2
               XI = SS(I)
               ETA = SS(J)
               WT = HH(I)*HH(J)

               CALL Q4RS_B_MATRICES ( XYZ, NORMALS, XI, ETA, BMB, BBB, BSB, EG, JAC )
               EPSM = MATMUL(BMB, UE_BASIC)
               NRES = MATMUL(SHELL_A, EPSM)
               CALL Q4RS_GRAD_ECOORDS ( XYZ, XI, ETA, GRADN, ECOORDS, EG, JAC )
               DO IA=1,4
                  DO IB=1,4
                     KGVAL = WT*JAC*( NRES(1)*GRADN(IA,1)*GRADN(IB,1) + NRES(2)*GRADN(IA,2)*GRADN(IB,2) +                    &
                                      NRES(3)*(GRADN(IA,1)*GRADN(IB,2) + GRADN(IA,2)*GRADN(IB,1)) )
                     DO RR=1,3
                        KGLOCAL(6*(IA-1)+RR,6*(IB-1)+RR) = KGLOCAL(6*(IA-1)+RR,6*(IB-1)+RR) + KGVAL
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         KED(1:24,1:24) = MATMUL(T24, MATMUL(KGLOCAL, T24T))
         IF ((DEBUG(233) > 0) .AND. (EID <= 8)) THEN
            WRITE(F06,'(A,I8,A,ES15.7)') 'CQUADR_Q4RS KGGD EID=', EID, ' KED_NORM=', DSQRT(SUM(KED(1:24,1:24)*KED(1:24,1:24)))
         ENDIF
      ENDIF

      IF (OPT(3) == 'Y') THEN
         REC_OPT = 'N'
         REC_OPT(3) = 'Y'
         CALL CQUADR_DKMQ24R ( REC_OPT, INT_ELEM_ID )
      ENDIF

      RETURN

 9001 FORMAT(' *ERROR: ',A,' expects ELGP=4 for element ',I8,' but got ',I8)

      CONTAINS

      SUBROUTINE LOAD_BASIC_COORDS ( XYZOUT )
      REAL(DOUBLE), INTENT(OUT) :: XYZOUT(4,3)
      INTEGER(LONG) :: II, JJ
      DO II=1,4
         DO JJ=1,3
            XYZOUT(II,JJ) = XEB(II,JJ)
         ENDDO
      ENDDO
      END SUBROUTINE LOAD_BASIC_COORDS

      SUBROUTINE SHAPE_N ( XI, ETA, NVAL )
      REAL(DOUBLE), INTENT(IN) :: XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: NVAL(4)
      NVAL(1) = 0.25D0*(ONE - XI)*(ONE - ETA)
      NVAL(2) = 0.25D0*(ONE + XI)*(ONE - ETA)
      NVAL(3) = 0.25D0*(ONE + XI)*(ONE + ETA)
      NVAL(4) = 0.25D0*(ONE - XI)*(ONE + ETA)
      END SUBROUTINE SHAPE_N

      SUBROUTINE SHAPE_DN ( XI, ETA, DN )
      REAL(DOUBLE), INTENT(IN) :: XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: DN(2,4)
      DN(1,1) = -0.25D0*(ONE - ETA)
      DN(1,2) =  0.25D0*(ONE - ETA)
      DN(1,3) =  0.25D0*(ONE + ETA)
      DN(1,4) = -0.25D0*(ONE + ETA)
      DN(2,1) = -0.25D0*(ONE - XI)
      DN(2,2) = -0.25D0*(ONE + XI)
      DN(2,3) =  0.25D0*(ONE + XI)
      DN(2,4) =  0.25D0*(ONE - XI)
      END SUBROUTINE SHAPE_DN

      SUBROUTINE Q4RS_GEOM ( XYZN, XI, ETA, EG, JAC )
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: EG(3,3), JAC
      REAL(DOUBLE) :: DN(2,4), G1(3), G2(3), E3(3), NM
      CALL SHAPE_DN(XI, ETA, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      NM = VNORM(G1)
      IF (NM <= 1.0D-14) THEN
         EG(:,1) = (/ONE,ZERO,ZERO/)
      ELSE
         EG(:,1) = G1/NM
      ENDIF
      CALL CROSSV(EG(:,1), G2, E3)
      NM = VNORM(E3)
      IF (NM <= 1.0D-14) THEN
         CALL CROSSV(G1, G2, E3)
         NM = VNORM(E3)
      ENDIF
      IF (NM <= 1.0D-14) THEN
         E3 = (/ZERO,ZERO,ONE/)
         NM = ONE
      ENDIF
      EG(:,3) = E3/NM
      CALL CROSSV(EG(:,3), EG(:,1), EG(:,2))
      JAC = VNORM(CROSS_PRODUCT(G1,G2))
      END SUBROUTINE Q4RS_GEOM

      SUBROUTINE Q4RS_GRAD_ECOORDS ( XYZN, XI, ETA, GRADN, ECOORDS, EG, JAC )
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: GRADN(4,2), ECOORDS(4,2), EG(3,3), JAC
      REAL(DOUBLE) :: DN(2,4), G1(3), G2(3), GM(2,2), GMI(2,2), TMP2(2), GV(3), CENT(3), DETG
      INTEGER(LONG) :: A
      CALL SHAPE_DN(XI, ETA, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      CALL Q4RS_GEOM(XYZN, XI, ETA, EG, JAC)
      GM(1,1) = DOT_PRODUCT(G1,G1)
      GM(1,2) = DOT_PRODUCT(G1,G2)
      GM(2,1) = GM(1,2)
      GM(2,2) = DOT_PRODUCT(G2,G2)
      DETG = GM(1,1)*GM(2,2) - GM(1,2)*GM(2,1)
      IF (DABS(DETG) <= 1.0D-14) DETG = SIGN(1.0D-14, DETG + 1.0D-30)
      GMI(1,1) =  GM(2,2)/DETG
      GMI(1,2) = -GM(1,2)/DETG
      GMI(2,1) = -GM(2,1)/DETG
      GMI(2,2) =  GM(1,1)/DETG
      DO A=1,4
         TMP2(1) = GMI(1,1)*DN(1,A) + GMI(1,2)*DN(2,A)
         TMP2(2) = GMI(2,1)*DN(1,A) + GMI(2,2)*DN(2,A)
         GV = G1*TMP2(1) + G2*TMP2(2)
         GRADN(A,1) = DOT_PRODUCT(EG(:,1), GV)
         GRADN(A,2) = DOT_PRODUCT(EG(:,2), GV)
      ENDDO
      CENT = (XYZN(1,:) + XYZN(2,:) + XYZN(3,:) + XYZN(4,:))/FOUR
      DO A=1,4
         GV = XYZN(A,:) - CENT
         ECOORDS(A,1) = DOT_PRODUCT(GV, EG(:,1))
         ECOORDS(A,2) = DOT_PRODUCT(GV, EG(:,2))
      ENDDO
      END SUBROUTINE Q4RS_GRAD_ECOORDS

      SUBROUTINE Q4RS_B_MATRICES ( XYZN, NORMS, XI, ETA, BM, BB, BS, EG, JAC )
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), NORMS(4,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BM(3,24), BB(3,24), BS(2,24), EG(3,3), JAC
      REAL(DOUBLE) :: GRADN(4,2), ECOORDS(4,2), T(24,24), TEMP(3,24), TEMPS(2,24)
      CALL Q4RS_GRAD_ECOORDS(XYZN, XI, ETA, GRADN, ECOORDS, EG, JAC)
      CALL Q4RS_T_MATRIX(NORMS, GRADN, EG, T)

      TEMP = ZERO
      DO K=1,4
         TEMP(1,(K-1)*6+1) = GRADN(K,1)
         TEMP(2,(K-1)*6+2) = GRADN(K,2)
         TEMP(3,(K-1)*6+1) = GRADN(K,2)
         TEMP(3,(K-1)*6+2) = GRADN(K,1)
      ENDDO
      BM = MATMUL(TEMP, T)

      TEMP = ZERO
      DO K=1,4
         TEMP(1,(K-1)*6+5) =  GRADN(K,1)
         TEMP(2,(K-1)*6+4) = -GRADN(K,2)
         TEMP(3,(K-1)*6+4) = -GRADN(K,1)
         TEMP(3,(K-1)*6+5) =  GRADN(K,2)
      ENDDO
      BB = MATMUL(TEMP, T)

      CALL Q4RS_BS_TEMP(XI, ETA, ECOORDS, TEMPS)
      BS = MATMUL(TEMPS, T)
      END SUBROUTINE Q4RS_B_MATRICES

      SUBROUTINE Q4RS_T_MATRIX ( NORMS, GRADN, EG, T )
      REAL(DOUBLE), INTENT(IN) :: NORMS(4,3), GRADN(4,2), EG(3,3)
      REAL(DOUBLE), INTENT(OUT) :: T(24,24)
      REAL(DOUBLE) :: TGA(24,24), TAE(24,24), A(3,3), NK(3), TBLOCK(3,3), A33, M1, M2, A3
      INTEGER(LONG) :: INODE, JNODE, R, C, ROFF, COFF
      TGA = ZERO
      TAE = ZERO
      DO INODE=1,4
         NK = MATMUL(TRANSPOSE(EG), NORMS(INODE,:))
         IF (VNORM(NK) <= 1.0D-12) NK = (/ZERO,ZERO,ONE/)
         NK = NK / VNORM(NK)
         CALL ROT_F3_TO_N(NK, A)
         TBLOCK = MATMUL(TRANSPOSE(A), TRANSPOSE(EG))
         ROFF = (INODE-1)*6
         TGA(ROFF+1:ROFF+3,ROFF+1:ROFF+3) = TBLOCK
         TGA(ROFF+4:ROFF+6,ROFF+4:ROFF+6) = TBLOCK
         TAE(ROFF+1:ROFF+3,ROFF+1:ROFF+3) = A
         A33 = A(3,3)
         DO R=1,2
            DO C=1,2
               IF (DABS(A33) > 1.0D-12) THEN
                  TAE(ROFF+3+R,ROFF+3+C) = A(R,C) - (A(R,3)*A(C,3))/A33
               ELSE
                  TAE(ROFF+3+R,ROFF+3+C) = A(R,C)
               ENDIF
            ENDDO
         ENDDO
         IF (DABS(A33) > 1.0D-12) THEN
            M1 = A(1,3)/A33
            M2 = A(2,3)/A33
            DO JNODE=1,4
               COFF = (JNODE-1)*6
               DO K=1,3
                  A3 = 0.5D0*(A(2,K)*GRADN(JNODE,1) - A(1,K)*GRADN(JNODE,2))
                  TAE(ROFF+4,COFF+K) = TAE(ROFF+4,COFF+K) + M1*A3
                  TAE(ROFF+5,COFF+K) = TAE(ROFF+5,COFF+K) + M2*A3
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      T = MATMUL(TAE, TGA)
      END SUBROUTINE Q4RS_T_MATRIX

      SUBROUTINE Q4RS_BS_TEMP ( XI, ETA, ECOORDS, BS )
      REAL(DOUBLE), INTENT(IN) :: XI, ETA, ECOORDS(4,2)
      REAL(DOUBLE), INTENT(OUT) :: BS(2,24)
      REAL(DOUBLE) :: X1,Y1,X2,Y2,X3,Y3,X4,Y4,J11,J21,J12,J22,ALEN,BLEN,CA,SA,CB,SB,DETJ
      REAL(DOUBLE) :: AX,AY,BX,BY,CX,CY,SSS,SR,R,S
      BS = ZERO
      R = XI
      S = ETA
      X1=ECOORDS(1,1); Y1=ECOORDS(1,2); X2=ECOORDS(2,1); Y2=ECOORDS(2,2)
      X3=ECOORDS(3,1); Y3=ECOORDS(3,2); X4=ECOORDS(4,1); Y4=ECOORDS(4,2)
      J11 = (X1*(S-ONE)/FOUR - X2*(S-ONE)/FOUR + X3*(S+ONE)/FOUR - X4*(S+ONE)/FOUR)
      J21 = (Y1*(S-ONE)/FOUR - Y2*(S-ONE)/FOUR + Y3*(S+ONE)/FOUR - Y4*(S+ONE)/FOUR)
      J12 = (X1*(R-ONE)/FOUR - X2*(R+ONE)/FOUR + X3*(R+ONE)/FOUR - X4*(R-ONE)/FOUR)
      J22 = (Y1*(R-ONE)/FOUR - Y2*(R+ONE)/FOUR + Y3*(R+ONE)/FOUR - Y4*(R-ONE)/FOUR)
      ALEN = DSQRT(J11*J11 + J21*J21)
      BLEN = DSQRT(J12*J12 + J22*J22)
      IF ((ALEN <= 1.0D-12) .OR. (BLEN <= 1.0D-12)) RETURN
      CA = J11/ALEN; SA = J21/ALEN; CB = J12/BLEN; SB = J22/BLEN
      DETJ = J11*J22 - J12*J21
      IF (DABS(DETJ) <= 1.0D-14) DETJ = SIGN(1.0D-14, DETJ + 1.0D-30)
      AX=X1-X2-X3+X4; AY=Y1-Y2-Y3+Y4
      BX=X1-X2+X3-X4; BY=Y1-Y2+Y3-Y4
      CX=X1+X2-X3-X4; CY=Y1+Y2-Y3-Y4
      SSS = DSQRT((AX+BX*S)**2 + (AY+BY*S)**2)
      SR  = DSQRT((BX*R+CX)**2 + (BY*R+CY)**2)

      CALL ADD_BS_ROW(BS, 1, 1,  SA*(R+ONE)*SSS - SB*(S+ONE)*SR, -SA*(Y1-Y4)*(R+ONE)*SSS + SB*(Y1-Y2)*(S+ONE)*SR, &
                               SA*(X1-X4)*(R+ONE)*SSS - SB*(X1-X2)*(S+ONE)*SR, DETJ)
      CALL ADD_BS_ROW(BS, 1, 2, -SA*(R-ONE)*SSS + SB*(S+ONE)*SR,  SA*(Y2-Y3)*(R-ONE)*SSS + SB*(Y1-Y2)*(S+ONE)*SR, &
                              -SA*(X2-X3)*(R-ONE)*SSS - SB*(X1-X2)*(S+ONE)*SR, DETJ)
      CALL ADD_BS_ROW(BS, 1, 3,  SA*(R-ONE)*SSS - SB*(S-ONE)*SR,  SA*(Y2-Y3)*(R-ONE)*SSS + SB*(Y3-Y4)*(S-ONE)*SR, &
                              -SA*(X2-X3)*(R-ONE)*SSS - SB*(X3-X4)*(S-ONE)*SR, DETJ)
      CALL ADD_BS_ROW(BS, 1, 4, -SA*(R+ONE)*SSS + SB*(S-ONE)*SR, -SA*(Y1-Y4)*(R+ONE)*SSS + SB*(Y3-Y4)*(S-ONE)*SR, &
                               SA*(X1-X4)*(R+ONE)*SSS - SB*(X3-X4)*(S-ONE)*SR, DETJ)

      CALL ADD_BS_ROW(BS, 2, 1, -CA*(R+ONE)*SSS + CB*(S+ONE)*SR,  CA*(Y1-Y4)*(R+ONE)*SSS - CB*(Y1-Y2)*(S+ONE)*SR, &
                              -CA*(X1-X4)*(R+ONE)*SSS + CB*(X1-X2)*(S+ONE)*SR, DETJ)
      CALL ADD_BS_ROW(BS, 2, 2,  CA*(R-ONE)*SSS - CB*(S+ONE)*SR, -CA*(Y2-Y3)*(R-ONE)*SSS - CB*(Y1-Y2)*(S+ONE)*SR, &
                               CA*(X2-X3)*(R-ONE)*SSS + CB*(X1-X2)*(S+ONE)*SR, DETJ)
      CALL ADD_BS_ROW(BS, 2, 3, -CA*(R-ONE)*SSS + CB*(S-ONE)*SR, -CA*(Y2-Y3)*(R-ONE)*SSS - CB*(Y3-Y4)*(S-ONE)*SR, &
                               CA*(X2-X3)*(R-ONE)*SSS + CB*(X3-X4)*(S-ONE)*SR, DETJ)
      CALL ADD_BS_ROW(BS, 2, 4,  CA*(R+ONE)*SSS - CB*(S-ONE)*SR,  CA*(Y1-Y4)*(R+ONE)*SSS - CB*(Y3-Y4)*(S-ONE)*SR, &
                              -CA*(X1-X4)*(R+ONE)*SSS + CB*(X3-X4)*(S-ONE)*SR, DETJ)
      END SUBROUTINE Q4RS_BS_TEMP

      SUBROUTINE ADD_BS_ROW(BS, IROW, INODE, WCOEF, RXCOEF, RYCOEF, DETJ)
      REAL(DOUBLE), INTENT(INOUT) :: BS(2,24)
      INTEGER(LONG), INTENT(IN) :: IROW, INODE
      REAL(DOUBLE), INTENT(IN) :: WCOEF, RXCOEF, RYCOEF, DETJ
      INTEGER(LONG) :: OFF
      OFF = (INODE-1)*6
      BS(IROW,OFF+3) = BS(IROW,OFF+3) + WCOEF/(16.0D0*DETJ)
      BS(IROW,OFF+4) = BS(IROW,OFF+4) + RXCOEF/(32.0D0*DETJ)
      BS(IROW,OFF+5) = BS(IROW,OFF+5) + RYCOEF/(32.0D0*DETJ)
      END SUBROUTINE ADD_BS_ROW

      SUBROUTINE CALC_Q4RS_NODAL_NORMALS ( XYZN, NORMS )
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3)
      REAL(DOUBLE), INTENT(OUT) :: NORMS(4,3)
      REAL(DOUBLE) :: EG(3,3), JAC, SN(3), NM
      INTEGER(LONG) :: II
      DO II=1,4
         CALL Q4RS_GEOM(XYZN, NODE_XI(II), NODE_ETA(II), EG, JAC)
         NORMS(II,:) = EG(:,3)
      ENDDO
      IF (ALLOCATED(GRID_SNORM)) THEN
         DO II=1,4
            IF ((BGRID(II) > 0) .AND. (BGRID(II) <= SIZE(GRID_SNORM,1))) THEN
               SN = GRID_SNORM(BGRID(II),:)
               NM = VNORM(SN)
               IF (NM > 1.0D-15) NORMS(II,:) = SN/NM
            ENDIF
         ENDDO
      ENDIF
      END SUBROUTINE CALC_Q4RS_NODAL_NORMALS

      FUNCTION Q4RS_DRILL_STIFFNESS ( KIN, NORMS ) RESULT(KD)
      REAL(DOUBLE), INTENT(IN) :: KIN(24,24), NORMS(4,3)
      REAL(DOUBLE) :: KD(24,24), N(3), P(3,3), KRR(3,3), KT(3,3), KAVG, VALS(4), NN
      INTEGER(LONG) :: II, BASE, NVAL
      KD = ZERO
      VALS = ZERO
      NVAL = 0
      DO II=1,4
         BASE = (II-1)*6 + 4
         N = NORMS(II,:)
         NN = VNORM(N)
         IF (NN <= 1.0D-12) CYCLE
         N = N/NN
         P = -OUTER3(N,N)
         P(1,1)=P(1,1)+ONE; P(2,2)=P(2,2)+ONE; P(3,3)=P(3,3)+ONE
         KRR = KIN(BASE:BASE+2,BASE:BASE+2)
         KT = MATMUL(P, MATMUL(KRR, P))
         NVAL = NVAL + 1
         VALS(NVAL) = MAX(ZERO, (KT(1,1)+KT(2,2)+KT(3,3))/TWO)
      ENDDO
      IF (NVAL > 0) THEN
         KAVG = SUM(VALS(1:NVAL))/REAL(NVAL,DOUBLE)
         DO II=1,4
            BASE = (II-1)*6 + 4
            N = NORMS(II,:)
            NN = VNORM(N)
            IF (NN <= 1.0D-12) CYCLE
            N = N/NN
            KD(BASE:BASE+2,BASE:BASE+2) = KD(BASE:BASE+2,BASE:BASE+2) + KAVG*OUTER3(N,N)
         ENDDO
      ENDIF
      END FUNCTION Q4RS_DRILL_STIFFNESS

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
            R = ZERO; R(1,1)=ONE; R(2,2)=ONE; R(3,3)=ONE
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

      FUNCTION OUTER3(A,B) RESULT(C)
      REAL(DOUBLE), INTENT(IN) :: A(3), B(3)
      REAL(DOUBLE) :: C(3,3)
      INTEGER(LONG) :: R,COL
      DO R=1,3
         DO COL=1,3
            C(R,COL)=A(R)*B(COL)
         ENDDO
      ENDDO
      END FUNCTION OUTER3

      SUBROUTINE CROSSV(A,B,C)
      REAL(DOUBLE), INTENT(IN) :: A(3), B(3)
      REAL(DOUBLE), INTENT(OUT) :: C(3)
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      END SUBROUTINE CROSSV

      FUNCTION CROSS_PRODUCT(A,B) RESULT(C)
      REAL(DOUBLE), INTENT(IN) :: A(3), B(3)
      REAL(DOUBLE) :: C(3)
      CALL CROSSV(A,B,C)
      END FUNCTION CROSS_PRODUCT

      FUNCTION VNORM ( V ) RESULT(NM)
      REAL(DOUBLE), INTENT(IN) :: V(3)
      REAL(DOUBLE) :: NM
      NM = DSQRT(DOT_PRODUCT(V,V))
      END FUNCTION VNORM

      FUNCTION QUAD_DIAMETER ( XYZN ) RESULT(D)
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3)
      REAL(DOUBLE) :: D, DD
      INTEGER(LONG) :: II
      D = ZERO
      DO II=2,4
         DD = VNORM(XYZN(II,:) - XYZN(1,:))
         IF (DD > D) D = DD
      ENDDO
      IF (D <= ZERO) D = ONE
      END FUNCTION QUAD_DIAMETER

      FUNCTION NODE_XI ( INODE ) RESULT(VAL)
      INTEGER(LONG), INTENT(IN) :: INODE
      REAL(DOUBLE) :: VAL
      IF ((INODE == 1) .OR. (INODE == 4)) THEN
         VAL = -ONE
      ELSE
         VAL = ONE
      ENDIF
      END FUNCTION NODE_XI

      FUNCTION NODE_ETA ( INODE ) RESULT(VAL)
      INTEGER(LONG), INTENT(IN) :: INODE
      REAL(DOUBLE) :: VAL
      IF ((INODE == 1) .OR. (INODE == 2)) THEN
         VAL = -ONE
      ELSE
         VAL = ONE
      ENDIF
      END FUNCTION NODE_ETA

      FUNCTION GP_INDEX ( IGP, JGP ) RESULT(IDX)
      INTEGER(LONG), INTENT(IN) :: IGP, JGP
      INTEGER(LONG) :: IDX
      IDX = 2*(IGP-1) + JGP
      END FUNCTION GP_INDEX

      SUBROUTINE BUILD_T24 ( T3, TOUT )
      REAL(DOUBLE), INTENT(IN) :: T3(3,3)
      REAL(DOUBLE), INTENT(OUT) :: TOUT(24,24)
      INTEGER(LONG) :: BI, R, C, OFF
      TOUT = ZERO
      DO BI=1,4
         OFF = 6*(BI-1)
         DO R=1,3
            DO C=1,3
               TOUT(OFF+R,OFF+C) = T3(R,C)
               TOUT(OFF+3+R,OFF+3+C) = T3(R,C)
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE BUILD_T24

      END SUBROUTINE CQUADR_Q4RS
