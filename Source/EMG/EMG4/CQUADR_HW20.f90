! #################################################################################################################################
! CQUADR HW20 shell for PARAM,QUADRTYP,HW20.

      SUBROUTINE CQUADR_HW20 ( OPT, INT_ELEM_ID )

! Port seed from D:\18a\bending_only\Shell\gemini2\shit\buckling\HW_20ShellElement_buckling.py.
! This kernel keeps the HW/Wagner-Gruttmann kinematics in a dedicated file and
! follows the Python reference mixed Hu-Washizu condensation (Nsigma/Neps/G/F/H)
! plus k_geom_local-style geometric stiffness for OPT(6).

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_GAUSS, NSUB
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, FOUR
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE MODEL_STUF, ONLY            :  EID, ELGP, KE, ME, BE1, BE2, BE3, EPROP, MASS_PER_UNIT_AREA, PRESS, PPE,                 &
                                         TE, NUM_EMG_FATAL_ERRS, SHELL_A, SHELL_D, SHELL_T, XEB, KED, UEL
      USE ELMDIS_Interface
      USE ELEM_STRE_STRN_ARRAYS_Interface
      USE ORDER_GAUSS_Interface
      USE OUTA_HERE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CQUADR_HW20'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      INTEGER(LONG), PARAMETER        :: NSTRESS = 14
      INTEGER(LONG), PARAMETER        :: NENH = 11
      INTEGER(LONG), PARAMETER        :: NSTRAIN = 25
      INTEGER(LONG), PARAMETER        :: NMIX = 39

      INTEGER(LONG)                   :: I,J,K,GP,JSUB,IA,IB,RR
      REAL(DOUBLE)                    :: XYZ(4,3), T1(3), T2(3), T3(3), T24(24,24), T24T(24,24)
      REAL(DOUBLE)                    :: SS(MAX_ORDER_GAUSS), HH(MAX_ORDER_GAUSS), XI, ETA, WT
      REAL(DOUBLE)                    :: BM(3,24), BB(3,24), BS(2,24), BML(3,24), BBL(3,24), BSL(2,24)
      REAL(DOUBLE)                    :: KBASIC(24,24), KLOCAL(24,24), KDRILL(24,24), MBASIC(24,24), MLOCAL(24,24)
      REAL(DOUBLE)                    :: KGLOCAL(24,24), KGVAL, EPSM(3), NRES(3), UE_BASIC(24)
      REAL(DOUBLE)                    :: JAC, NVG(4), MASS_DIAG(4), UNIT_PPE_B(24), UNIT_PPE_L(24), GRADN(4,2)
      REAL(DOUBLE)                    :: TSIG0(3,3), TEPS0(3,3), TTILDE0(2,2), J0DET, CSHAPE
      REAL(DOUBLE)                    :: B8(8,24), NSIG(8,NSTRESS), NEPS(8,NSTRAIN)
      REAL(DOUBLE)                    :: GMIX(NSTRESS,24), FMIX(NSTRAIN,NSTRESS), HMIX(NSTRAIN,NSTRAIN)
      REAL(DOUBLE)                    :: AMIX(NMIX,NMIX), AINV(NMIX,NMIX), BMIX(NMIX,24)

      IF (ELGP /= 4) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9001) SUBR_NAME, EID, ELGP
         WRITE(F06,9001) SUBR_NAME, EID, ELGP
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL LOAD_BASIC_COORDS(XYZ)
      CALL HW20_LOCAL_BASIS(XYZ, T1, T2, T3)
      CALL BUILD_T24(TE, T24)
      T24T = TRANSPOSE(T24)

      IF ((OPT(3) == 'Y') .OR. (OPT(4) == 'Y') .OR. (OPT(1) == 'Y') .OR. (OPT(5) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         CALL ORDER_GAUSS(2, SS, HH)
      ENDIF

      IF ((OPT(3) == 'Y') .OR. (OPT(4) == 'Y')) THEN
         CALL HW20_METRIC_AND_C(XYZ, T1, T2, TSIG0, TEPS0, TTILDE0, J0DET, CSHAPE)
         GMIX = ZERO
         FMIX = ZERO
         HMIX = ZERO

         DO I=1,3
            DO J=1,3
               XI = GAUSS3_X(I)
               ETA = GAUSS3_X(J)
               WT = GAUSS3_W(I)*GAUSS3_W(J)
               CALL HW20_B_MATRICES(XYZ, T1, T2, T3, XI, ETA, BM, BB, BS, JAC)
               CALL HW20_BUILD_NSIGMA(XI, ETA, TSIG0, TTILDE0, NSIG)
               CALL HW20_BUILD_NEPS(XYZ, XI, ETA, TEPS0, TTILDE0, J0DET, CSHAPE, NEPS)
               B8(1:3,:) = BM
               B8(4:6,:) = BB
               B8(7:8,:) = BS
               GMIX = GMIX + WT*JAC*MATMUL(TRANSPOSE(NSIG), B8)
               FMIX = FMIX - WT*JAC*MATMUL(TRANSPOSE(NEPS), NSIG)
               HMIX = HMIX + WT*JAC*HW20_HMAT(NEPS)

               IF ((I <= 2) .AND. (J <= 2)) THEN
                  GP = GP_INDEX(I,J)
                  BML = MATMUL(BM, T24T)
                  BBL = MATMUL(BB, T24T)
                  BSL = MATMUL(BS, T24T)
                  BE1(1:3,1:24,GP) = BML
                  BE2(1:3,1:24,GP) = BBL
                  BE3(1:2,1:24,GP) = BSL
               ENDIF
            ENDDO
         ENDDO

         FMIX(15:NSTRAIN,:) = ZERO
         AMIX = ZERO
         AMIX(1:NSTRESS,NSTRESS+1:NMIX) = TRANSPOSE(FMIX)
         AMIX(NSTRESS+1:NMIX,1:NSTRESS) = FMIX
         AMIX(NSTRESS+1:NMIX,NSTRESS+1:NMIX) = HMIX
         DO I=1,NMIX
            AMIX(I,I) = AMIX(I,I) + 1.0D-12
         ENDDO
         BMIX = ZERO
         BMIX(1:NSTRESS,:) = GMIX
         CALL INVN(AMIX, AINV, NMIX)
         KBASIC = -MATMUL(TRANSPOSE(BMIX), MATMUL(AINV, BMIX))

         KDRILL = HW20_DRILL_STIFFNESS_3GP(XYZ, T3)
         KBASIC = 0.5D0*(KBASIC + TRANSPOSE(KBASIC)) + KDRILL
         KLOCAL = MATMUL(T24, MATMUL(KBASIC, T24T))
         IF (OPT(4) == 'Y') KE(1:24,1:24) = KLOCAL

         IF (OPT(3) == 'Y') THEN
            BE1(:,:,1) = (BE1(:,:,1) + BE1(:,:,2) + BE1(:,:,3) + BE1(:,:,4))/FOUR
            BE2(:,:,1) = (BE2(:,:,1) + BE2(:,:,2) + BE2(:,:,3) + BE2(:,:,4))/FOUR
            BE3(1:2,:,1) = (BE3(1:2,:,1) + BE3(1:2,:,2) + BE3(1:2,:,3) + BE3(1:2,:,4))/FOUR
         ENDIF
      ENDIF

      IF (OPT(1) == 'Y') THEN
         MBASIC = ZERO
         MASS_DIAG = ZERO
         DO I=1,2
            DO J=1,2
               XI = SS(I)
               ETA = SS(J)
               WT = HH(I)*HH(J)
               CALL SHAPE_N(XI, ETA, NVG)
               JAC = HW20_DETJ(XYZ, XI, ETA)
               DO GP=1,4
                  MASS_DIAG(GP) = MASS_DIAG(GP) + NVG(GP)*MASS_PER_UNIT_AREA*WT*JAC
               ENDDO
            ENDDO
         ENDDO
         DO I=1,4
            DO K=1,3
               MBASIC((I-1)*6+K,(I-1)*6+K) = MASS_DIAG(I)
            ENDDO
         ENDDO
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
               CALL SHAPE_N(XI, ETA, NVG)
               JAC = HW20_DETJ(XYZ, XI, ETA)
               DO GP=1,4
                  UNIT_PPE_B((GP-1)*6+1:(GP-1)*6+3) = UNIT_PPE_B((GP-1)*6+1:(GP-1)*6+3) + T3*NVG(GP)*WT*JAC
               ENDDO
            ENDDO
         ENDDO
         UNIT_PPE_L = MATMUL(T24, UNIT_PPE_B)
         DO JSUB=1,NSUB
            PPE(1:24,JSUB) = UNIT_PPE_L(1:24)*PRESS(3,JSUB)
         ENDDO
      ENDIF

      IF (OPT(6) == 'Y') THEN
         CALL ELMDIS
         UE_BASIC = MATMUL(T24T, UEL(1:24))
         KGLOCAL = ZERO
         DO I=1,3
            DO J=1,3
               XI = GAUSS3_X(I)
               ETA = GAUSS3_X(J)
               WT = GAUSS3_W(I)*GAUSS3_W(J)
               CALL HW20_B_MATRICES(XYZ, T1, T2, T3, XI, ETA, BM, BB, BS, JAC)
               EPSM = MATMUL(BM, UE_BASIC)
               NRES = MATMUL(SHELL_A, EPSM)
               CALL HW20_GRADN_LOCAL(XYZ, T1, T2, XI, ETA, GRADN)
               DO IA=1,4
                  DO IB=1,4
                     KGVAL = WT*JAC*(NRES(1)*GRADN(IA,1)*GRADN(IB,1) + NRES(2)*GRADN(IA,2)*GRADN(IB,2) +                     &
                                     NRES(3)*(GRADN(IA,1)*GRADN(IB,2) + GRADN(IA,2)*GRADN(IB,1)))
                     DO RR=1,3
                        KGLOCAL(6*(IA-1)+RR,6*(IB-1)+RR) = KGLOCAL(6*(IA-1)+RR,6*(IB-1)+RR) + KGVAL
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         KED(1:24,1:24) = MATMUL(T24, MATMUL(KGLOCAL, T24T))
      ENDIF

      RETURN

 9001 FORMAT(' *ERROR: ',A,' expects ELGP=4 for element ',I8,' but got ',I8)

      CONTAINS

      SUBROUTINE LOAD_BASIC_COORDS(XYZOUT)
      REAL(DOUBLE), INTENT(OUT) :: XYZOUT(4,3)
      INTEGER(LONG) :: II, JJ
      DO II=1,4
         DO JJ=1,3
            XYZOUT(II,JJ) = XEB(II,JJ)
         ENDDO
      ENDDO
      END SUBROUTINE LOAD_BASIC_COORDS

      SUBROUTINE SHAPE_N(XI, ETA, NVAL)
      REAL(DOUBLE), INTENT(IN) :: XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: NVAL(4)
      NVAL(1)=0.25D0*(ONE-XI)*(ONE-ETA)
      NVAL(2)=0.25D0*(ONE+XI)*(ONE-ETA)
      NVAL(3)=0.25D0*(ONE+XI)*(ONE+ETA)
      NVAL(4)=0.25D0*(ONE-XI)*(ONE+ETA)
      END SUBROUTINE SHAPE_N

      SUBROUTINE SHAPE_DN(XI, ETA, DN)
      REAL(DOUBLE), INTENT(IN) :: XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: DN(2,4)
      DN(1,1)=-0.25D0*(ONE-ETA); DN(1,2)=0.25D0*(ONE-ETA); DN(1,3)=0.25D0*(ONE+ETA); DN(1,4)=-0.25D0*(ONE+ETA)
      DN(2,1)=-0.25D0*(ONE-XI);  DN(2,2)=-0.25D0*(ONE+XI); DN(2,3)=0.25D0*(ONE+XI); DN(2,4)=0.25D0*(ONE-XI)
      END SUBROUTINE SHAPE_DN

      SUBROUTINE HW20_LOCAL_BASIS(XYZ, T1, T2, T3)
      REAL(DOUBLE), INTENT(IN) :: XYZ(4,3)
      REAL(DOUBLE), INTENT(OUT) :: T1(3), T2(3), T3(3)
      REAL(DOUBLE) :: D1(3), D2(3), E1(3), E2(3), NM
      D1 = XYZ(3,:) - XYZ(1,:)
      D2 = XYZ(2,:) - XYZ(4,:)
      E1 = D1/MAX(VNORM(D1),1.0D-15)
      E2 = D2/MAX(VNORM(D2),1.0D-15)
      T1 = E1 + E2
      NM = VNORM(T1)
      IF (NM <= 1.0D-15) T1 = E1
      T1 = T1/MAX(VNORM(T1),1.0D-15)
      T2 = E1 - E2
      NM = VNORM(T2)
      IF (NM <= 1.0D-15) THEN
         CALL CROSSV((/ZERO,ZERO,ONE/), T1, T2)
      ENDIF
      T2 = T2/MAX(VNORM(T2),1.0D-15)
      CALL CROSSV(T1, T2, T3)
      T3 = T3/MAX(VNORM(T3),1.0D-15)
      CALL CROSSV(T3, T1, T2)
      T2 = T2/MAX(VNORM(T2),1.0D-15)
      END SUBROUTINE HW20_LOCAL_BASIS

      FUNCTION HW20_DETJ(XYZ, XI, ETA) RESULT(JAC)
      REAL(DOUBLE), INTENT(IN) :: XYZ(4,3), XI, ETA
      REAL(DOUBLE) :: JAC, DN(2,4), G1(3), G2(3), C(3)
      CALL SHAPE_DN(XI, ETA, DN)
      G1 = MATMUL(DN(1,:), XYZ)
      G2 = MATMUL(DN(2,:), XYZ)
      CALL CROSSV(G1, G2, C)
      JAC = VNORM(C)
      END FUNCTION HW20_DETJ

      SUBROUTINE HW20_GRADN_LOCAL(XYZ, T1, T2, XI, ETA, GRADN)
      REAL(DOUBLE), INTENT(IN) :: XYZ(4,3), T1(3), T2(3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: GRADN(4,2)
      REAL(DOUBLE) :: DN(2,4), G1(3), G2(3), A11,A22,A12,DET,GC1(3),GC2(3)
      INTEGER(LONG) :: A
      CALL SHAPE_DN(XI, ETA, DN)
      G1 = MATMUL(DN(1,:), XYZ)
      G2 = MATMUL(DN(2,:), XYZ)
      A11=DOT_PRODUCT(G1,G1); A22=DOT_PRODUCT(G2,G2); A12=DOT_PRODUCT(G1,G2)
      DET=A11*A22-A12*A12
      IF (DABS(DET) <= 1.0D-12) DET = SIGN(1.0D-12, DET + 1.0D-30)
      GC1=(A22*G1-A12*G2)/DET
      GC2=(-A12*G1+A11*G2)/DET
      DO A=1,4
         GRADN(A,1) = DN(1,A)*DOT_PRODUCT(T1,GC1) + DN(2,A)*DOT_PRODUCT(T1,GC2)
         GRADN(A,2) = DN(1,A)*DOT_PRODUCT(T2,GC1) + DN(2,A)*DOT_PRODUCT(T2,GC2)
      ENDDO
      END SUBROUTINE HW20_GRADN_LOCAL

      SUBROUTINE HW20_B_MATRICES(XYZ, T1, T2, T3, XI, ETA, BM, BB, BS, JAC)
      REAL(DOUBLE), INTENT(IN) :: XYZ(4,3), T1(3), T2(3), T3(3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BM(3,24), BB(3,24), BS(2,24), JAC
      REAL(DOUBLE) :: DN(2,4), G1(3), G2(3), A11,A22,A12,DET,GC1(3),GC2(3)
      REAL(DOUBLE) :: C11,C12,C21,C22, A1,A2, DE11(3),DE22(3),DE12(3), XX(3),YY(3),XY(3)
      REAL(DOUBLE) :: CG1(3),CG2(3), DR11(3),DR22(3),DR12(3), T0XI1(3),T0XI2(3)
      INTEGER(LONG) :: A, OFF
      CALL SHAPE_DN(XI, ETA, DN)
      G1 = MATMUL(DN(1,:), XYZ)
      G2 = MATMUL(DN(2,:), XYZ)
      A11=DOT_PRODUCT(G1,G1); A22=DOT_PRODUCT(G2,G2); A12=DOT_PRODUCT(G1,G2)
      DET=A11*A22-A12*A12
      IF (DABS(DET) <= 1.0D-12) DET = SIGN(1.0D-12, DET + 1.0D-30)
      GC1=(A22*G1-A12*G2)/DET
      GC2=(-A12*G1+A11*G2)/DET
      C11=DOT_PRODUCT(T1,GC1); C12=DOT_PRODUCT(T1,GC2)
      C21=DOT_PRODUCT(T2,GC1); C22=DOT_PRODUCT(T2,GC2)
      JAC = HW20_DETJ(XYZ, XI, ETA)
      BM = ZERO; BB = ZERO
      T0XI1 = ZERO; T0XI2 = ZERO
      DO A=1,4
         T0XI1 = T0XI1 + DN(1,A)*T3
         T0XI2 = T0XI2 + DN(2,A)*T3
      ENDDO
      DO A=1,4
         OFF=(A-1)*6
         A1=DN(1,A); A2=DN(2,A)
         DE11=A1*G1; DE22=A2*G2; DE12=0.5D0*(A1*G2+A2*G1)
         XX=C11*C11*DE11 + C12*C12*DE22 + TWO*C11*C12*DE12
         YY=C21*C21*DE11 + C22*C22*DE22 + TWO*C21*C22*DE12
         XY=TWO*(C11*C21*DE11 + C12*C22*DE22 + (C11*C22+C12*C21)*DE12)
         BM(1,OFF+1:OFF+3)=XX; BM(2,OFF+1:OFF+3)=YY; BM(3,OFF+1:OFF+3)=XY
         CALL CROSSV(T3,G1,CG1); CALL CROSSV(T3,G2,CG2)
         DR11=A1*CG1; DR22=A2*CG2; DR12=0.5D0*(A1*CG2+A2*CG1)
         XX=C11*C11*DR11 + C12*C12*DR22 + TWO*C11*C12*DR12
         YY=C21*C21*DR11 + C22*C22*DR22 + TWO*C21*C22*DR12
         XY=TWO*(C11*C21*DR11 + C12*C22*DR22 + (C11*C22+C12*C21)*DR12)
         BB(1,OFF+4:OFF+6)=XX; BB(2,OFF+4:OFF+6)=YY; BB(3,OFF+4:OFF+6)=XY
      ENDDO
      CALL HW20_BS_DIRECT(XYZ, T1, T2, T3, XI, ETA, BS)
      END SUBROUTINE HW20_B_MATRICES

      SUBROUTINE HW20_BS_DIRECT(XYZ, T1, T2, T3, XI, ETA, BS)
      REAL(DOUBLE), INTENT(IN) :: XYZ(4,3), T1(3), T2(3), T3(3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BS(2,24)
      REAL(DOUBLE) :: BSA(2,24), BSB(2,24), BSC(2,24), BSD(2,24), BSN(2,24)
      REAL(DOUBLE) :: DN(2,4), N(4), G1(3), G2(3), A11,A22,A12,DET,GC1(3),GC2(3),C11,C12,C21,C22
      CALL HW20_BS_NAT(XYZ,T3,ZERO,-ONE,BSA)
      CALL HW20_BS_NAT(XYZ,T3,ZERO, ONE,BSB)
      CALL HW20_BS_NAT(XYZ,T3,-ONE,ZERO,BSC)
      CALL HW20_BS_NAT(XYZ,T3, ONE,ZERO,BSD)
      BSN = ZERO
      BSN(1,:) = 0.5D0*(ONE-ETA)*BSA(1,:) + 0.5D0*(ONE+ETA)*BSB(1,:)
      BSN(2,:) = 0.5D0*(ONE-XI )*BSC(2,:) + 0.5D0*(ONE+XI )*BSD(2,:)
      CALL SHAPE_DN(XI, ETA, DN)
      CALL SHAPE_N(XI, ETA, N)
      G1=MATMUL(DN(1,:),XYZ); G2=MATMUL(DN(2,:),XYZ)
      A11=DOT_PRODUCT(G1,G1); A22=DOT_PRODUCT(G2,G2); A12=DOT_PRODUCT(G1,G2)
      DET=A11*A22-A12*A12
      IF (DABS(DET) <= 1.0D-12) DET = SIGN(1.0D-12, DET + 1.0D-30)
      GC1=(A22*G1-A12*G2)/DET; GC2=(-A12*G1+A11*G2)/DET
      C11=DOT_PRODUCT(T1,GC1); C12=DOT_PRODUCT(T1,GC2)
      C21=DOT_PRODUCT(T2,GC1); C22=DOT_PRODUCT(T2,GC2)
      BS(1,:) = C11*BSN(1,:) + C12*BSN(2,:)
      BS(2,:) = C21*BSN(1,:) + C22*BSN(2,:)
      END SUBROUTINE HW20_BS_DIRECT

      SUBROUTINE HW20_BS_NAT(XYZ,T3,XI,ETA,BS)
      REAL(DOUBLE), INTENT(IN) :: XYZ(4,3), T3(3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BS(2,24)
      REAL(DOUBLE) :: DN(2,4), N(4), G1(3), G2(3), C1(3), C2(3)
      INTEGER(LONG) :: A, OFF
      CALL SHAPE_DN(XI,ETA,DN); CALL SHAPE_N(XI,ETA,N)
      G1=MATMUL(DN(1,:),XYZ); G2=MATMUL(DN(2,:),XYZ)
      BS=ZERO
      DO A=1,4
         OFF=(A-1)*6
         BS(1,OFF+1:OFF+3)=BS(1,OFF+1:OFF+3)+DN(1,A)*T3
         BS(2,OFF+1:OFF+3)=BS(2,OFF+1:OFF+3)+DN(2,A)*T3
         CALL CROSSV(T3,G1,C1); CALL CROSSV(T3,G2,C2)
         BS(1,OFF+4:OFF+6)=BS(1,OFF+4:OFF+6)+N(A)*C1
         BS(2,OFF+4:OFF+6)=BS(2,OFF+4:OFF+6)+N(A)*C2
      ENDDO
      END SUBROUTINE HW20_BS_NAT

      FUNCTION HW20_DRILL_STIFFNESS_3GP(XYZ,T3) RESULT(KD)
      REAL(DOUBLE), INTENT(IN) :: XYZ(4,3), T3(3)
      REAL(DOUBLE) :: KD(24,24), AREA, HAVG, DK, JAC
      INTEGER(LONG) :: A, B, OFF
      AREA=ZERO
      DO A=1,3
         DO B=1,3
            JAC=HW20_DETJ(XYZ, GAUSS3_X(A), GAUSS3_X(B))
            AREA=AREA+GAUSS3_W(A)*GAUSS3_W(B)*JAC
         ENDDO
      ENDDO
      HAVG=EPROP(1)
      DK=0.02D0*(SHELL_T(1,1)+SHELL_T(2,2))/TWO*AREA/FOUR*0.001D0
      KD=ZERO
      DO A=1,4
         OFF=(A-1)*6+4
         KD(OFF:OFF+2,OFF:OFF+2)=KD(OFF:OFF+2,OFF:OFF+2)+DK*OUTER3(T3,T3)
      ENDDO
      END FUNCTION HW20_DRILL_STIFFNESS_3GP

      FUNCTION GAUSS3_X(I) RESULT(X)
      INTEGER(LONG), INTENT(IN) :: I
      REAL(DOUBLE) :: X
      IF (I == 1) THEN
         X = -DSQRT(3.0D0/5.0D0)
      ELSE IF (I == 2) THEN
         X = ZERO
      ELSE
         X = DSQRT(3.0D0/5.0D0)
      ENDIF
      END FUNCTION GAUSS3_X

      FUNCTION GAUSS3_W(I) RESULT(W)
      INTEGER(LONG), INTENT(IN) :: I
      REAL(DOUBLE) :: W
      IF (I == 2) THEN
         W = 8.0D0/9.0D0
      ELSE
         W = 5.0D0/9.0D0
      ENDIF
      END FUNCTION GAUSS3_W

      SUBROUTINE HW20_METRIC_AND_C(XYZ, T1, T2, TSIG, TEPS, TTILDE, J0, CSH)
      REAL(DOUBLE), INTENT(IN) :: XYZ(4,3), T1(3), T2(3)
      REAL(DOUBLE), INTENT(OUT) :: TSIG(3,3), TEPS(3,3), TTILDE(2,2), J0, CSH
      REAL(DOUBLE) :: DN(2,4), G1(3), G2(3), CR(3), G11, G22, G12, TRC, DET, DISC, L1, L2
      REAL(DOUBLE) :: J11,J12,J21,J22
      CALL SHAPE_DN(ZERO, ZERO, DN)
      G1 = MATMUL(DN(1,:), XYZ)
      G2 = MATMUL(DN(2,:), XYZ)
      G11 = DOT_PRODUCT(G1,G1)
      G22 = DOT_PRODUCT(G2,G2)
      G12 = DOT_PRODUCT(G1,G2)
      TRC = G11 + G22
      DET = G11*G22 - G12*G12
      DISC = MAX(TRC*TRC - FOUR*DET, ZERO)
      L1 = 0.5D0*(TRC + DSQRT(DISC))
      L2 = 0.5D0*(TRC - DSQRT(DISC))
      IF (L2 > 1.0D-12) THEN
         CSH = L1/L2
      ELSE
         CSH = ONE
      ENDIF
      CSH = MIN(MAX(CSH, 0.25D0), FOUR)
      J11 = DOT_PRODUCT(G1,T1)
      J12 = DOT_PRODUCT(G1,T2)
      J21 = DOT_PRODUCT(G2,T1)
      J22 = DOT_PRODUCT(G2,T2)
      TTILDE(1,1)=J11; TTILDE(1,2)=J12
      TTILDE(2,1)=J21; TTILDE(2,2)=J22
      CALL CROSSV(G1, G2, CR)
      J0 = VNORM(CR)
      TSIG = ZERO
      TSIG(1,1)=J11*J11;       TSIG(1,2)=J12*J12;       TSIG(1,3)=J11*J12
      TSIG(2,1)=J21*J21;       TSIG(2,2)=J22*J22;       TSIG(2,3)=J21*J22
      TSIG(3,1)=TWO*J11*J21;   TSIG(3,2)=TWO*J12*J22;   TSIG(3,3)=J11*J22+J12*J21
      TEPS = ZERO
      TEPS(1,1)=J11*J11;       TEPS(1,2)=J12*J12;       TEPS(1,3)=TWO*J11*J12
      TEPS(2,1)=J21*J21;       TEPS(2,2)=J22*J22;       TEPS(2,3)=TWO*J21*J22
      TEPS(3,1)=J11*J21;       TEPS(3,2)=J12*J22;       TEPS(3,3)=J11*J22+J12*J21
      END SUBROUTINE HW20_METRIC_AND_C

      SUBROUTINE HW20_BUILD_NSIGMA(XI, ETA, TSIG, TTILDE, NSIG)
      REAL(DOUBLE), INTENT(IN) :: XI, ETA, TSIG(3,3), TTILDE(2,2)
      REAL(DOUBLE), INTENT(OUT) :: NSIG(8,NSTRESS)
      REAL(DOUBLE) :: MM(3,2), MS(2,2), VARM(3,2), VARS(2,2)
      NSIG = ZERO
      NSIG(1,1)=ONE; NSIG(2,2)=ONE; NSIG(3,3)=ONE
      NSIG(4,4)=ONE; NSIG(5,5)=ONE; NSIG(6,6)=ONE
      NSIG(7,7)=ONE; NSIG(8,8)=ONE
      MM = ZERO
      MM(1,1)=ETA
      MM(2,2)=XI
      MS = ZERO
      MS(1,1)=ETA
      MS(2,2)=XI
      VARM = MATMUL(TSIG, MM)
      VARS = MATMUL(TTILDE, MS)
      NSIG(1:3,9:10) = VARM
      NSIG(4:6,11:12) = VARM
      NSIG(7:8,13:14) = VARS
      END SUBROUTINE HW20_BUILD_NSIGMA

      SUBROUTINE HW20_BUILD_NEPS(XYZ, XI, ETA, TEPS, TTILDE, J0, CSH, NEPS)
      REAL(DOUBLE), INTENT(IN) :: XYZ(4,3), XI, ETA, TEPS(3,3), TTILDE(2,2), J0, CSH
      REAL(DOUBLE), INTENT(OUT) :: NEPS(8,NSTRAIN)
      REAL(DOUBLE) :: MM(3,2), MS(2,2), VARM(3,2), VARS(2,2), MN(3,NENH), NENH_MAT(8,NENH)
      REAL(DOUBLE) :: DN(2,4), G1(3), G2(3), CR(3), JAC, FAC, XI2, ETA2
      NEPS = ZERO
      NEPS(1,1)=ONE; NEPS(2,2)=ONE; NEPS(3,3)=ONE
      NEPS(4,4)=ONE; NEPS(5,5)=ONE; NEPS(6,6)=ONE
      NEPS(7,7)=ONE; NEPS(8,8)=ONE
      MM = ZERO
      MM(1,1)=ETA
      MM(2,2)=XI
      MS = ZERO
      MS(1,1)=ETA
      MS(2,2)=XI
      VARM = MATMUL(TEPS, MM)
      VARS = MATMUL(TTILDE, MS)
      NEPS(1:3,9:10) = VARM
      NEPS(4:6,11:12) = VARM
      NEPS(7:8,13:14) = VARS
      XI2 = XI*XI
      ETA2 = ETA*ETA
      MN = ZERO
      MN(1,1)=XI
      MN(2,2)=ETA
      MN(3,3)=XI
      MN(3,4)=ETA
      MN(1,5)=XI*ETA
      MN(2,6)=XI*ETA
      MN(3,7)=XI*ETA
      MN(1,8)=(XI2-CSH)*ETA
      MN(2,9)=(ETA2-CSH)*XI
      MN(1,10)=ETA2*XI
      MN(2,11)=XI2*ETA
      CALL SHAPE_DN(XI, ETA, DN)
      G1 = MATMUL(DN(1,:), XYZ)
      G2 = MATMUL(DN(2,:), XYZ)
      CALL CROSSV(G1, G2, CR)
      JAC = VNORM(CR)
      IF (JAC > 1.0D-12) THEN
         FAC = J0/JAC
      ELSE
         FAC = ONE
      ENDIF
      NENH_MAT = ZERO
      NENH_MAT(1:3,:) = FAC*MATMUL(TEPS, MN)
      NEPS(1:8,15:25) = NENH_MAT
      END SUBROUTINE HW20_BUILD_NEPS

      FUNCTION HW20_HMAT(NEPS) RESULT(HOUT)
      REAL(DOUBLE), INTENT(IN) :: NEPS(8,NSTRAIN)
      REAL(DOUBLE) :: HOUT(NSTRAIN,NSTRAIN)
      HOUT = ZERO
      HOUT = HOUT + MATMUL(TRANSPOSE(NEPS(1:3,:)), MATMUL(SHELL_A, NEPS(1:3,:)))
      HOUT = HOUT + MATMUL(TRANSPOSE(NEPS(4:6,:)), MATMUL(SHELL_D, NEPS(4:6,:)))
      HOUT = HOUT + MATMUL(TRANSPOSE(NEPS(7:8,:)), MATMUL(SHELL_T, NEPS(7:8,:)))
      END FUNCTION HW20_HMAT

      SUBROUTINE INVN(A, AINV, N)
      INTEGER(LONG), INTENT(IN) :: N
      REAL(DOUBLE), INTENT(IN) :: A(N,N)
      REAL(DOUBLE), INTENT(OUT) :: AINV(N,N)
      REAL(DOUBLE) :: AUG(N,2*N), PIV, FACT, ROWTMP(2*N)
      INTEGER(LONG) :: I, J, PIVROW
      AUG = ZERO
      AUG(:,1:N) = A
      DO I=1,N
         AUG(I,N+I) = ONE
      ENDDO
      DO I=1,N
         PIVROW = I
         DO J=I+1,N
            IF (DABS(AUG(J,I)) > DABS(AUG(PIVROW,I))) PIVROW = J
         ENDDO
         IF (PIVROW /= I) THEN
            ROWTMP = AUG(I,:)
            AUG(I,:) = AUG(PIVROW,:)
            AUG(PIVROW,:) = ROWTMP
         ENDIF
         PIV = AUG(I,I)
         IF (DABS(PIV) < 1.0D-20) PIV = SIGN(1.0D-20, PIV + 1.0D-30)
         AUG(I,:) = AUG(I,:)/PIV
         DO J=1,N
            IF (J /= I) THEN
               FACT = AUG(J,I)
               AUG(J,:) = AUG(J,:) - FACT*AUG(I,:)
            ENDIF
         ENDDO
      ENDDO
      AINV = AUG(:,N+1:2*N)
      END SUBROUTINE INVN

      SUBROUTINE CROSSV(A,B,C)
      REAL(DOUBLE), INTENT(IN) :: A(3), B(3)
      REAL(DOUBLE), INTENT(OUT) :: C(3)
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      END SUBROUTINE CROSSV

      FUNCTION VNORM(V) RESULT(NM)
      REAL(DOUBLE), INTENT(IN) :: V(3)
      REAL(DOUBLE) :: NM
      NM=DSQRT(DOT_PRODUCT(V,V))
      END FUNCTION VNORM

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

      FUNCTION GP_INDEX(IGP,JGP) RESULT(IDX)
      INTEGER(LONG), INTENT(IN) :: IGP,JGP
      INTEGER(LONG) :: IDX
      IDX=2*(IGP-1)+JGP
      END FUNCTION GP_INDEX

      SUBROUTINE BUILD_T24(T3IN,TOUT)
      REAL(DOUBLE), INTENT(IN) :: T3IN(3,3)
      REAL(DOUBLE), INTENT(OUT) :: TOUT(24,24)
      INTEGER(LONG) :: BI,R,C,OFF
      TOUT=ZERO
      DO BI=1,4
         OFF=6*(BI-1)
         DO R=1,3
            DO C=1,3
               TOUT(OFF+R,OFF+C)=T3IN(R,C)
               TOUT(OFF+3+R,OFF+3+C)=T3IN(R,C)
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE BUILD_T24

      END SUBROUTINE CQUADR_HW20
