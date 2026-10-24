! Experimental Nguyen2024 composite curvature reconstruction. Explicit opt-in only.
      SUBROUTINE CTRIAR_KMITC3P(OPT,INT_ELEM_ID)
      USE PENTIUM_II_KIND, ONLY: BYTE,LONG,DOUBLE
      USE CONSTANTS_1, ONLY: ZERO,ONE,TWO,THREE
      USE IOUNT1, ONLY: ERR,F06
      USE SCONTR, ONLY: BLNK_SUB_NAM,NSUB,NTSUB,SOL_NAME
      USE PARAMS, ONLY: COUPMASS
      USE MODEL_STUF, ONLY: EID,TYPE,ELGP,ELDOF,EPROP,XEB,TE,KE,KED,ME,BE1,BE2,BE3,SE1,SE2,SE3, &
       STE1,STE2,SHELL_A,SHELL_D,SHELL_T,SHELL_DALP,EM,EB,ET,MASS_PER_UNIT_AREA,PTE,PPE,PRESS,DT,ALPVEC,TREF,UEL
      USE ELMDIS_Interface
      USE OUTA_HERE_Interface
      IMPLICIT NONE
      CHARACTER(1*BYTE),INTENT(IN):: OPT(6)
      INTEGER(LONG),INTENT(IN):: INT_ELEM_ID
      INTEGER(LONG):: I,J,K,G,IA,IB,GP
      REAL(DOUBLE):: XYZ(3,3),EX(3),EY(3),EZ(3),VP(3),XY(3,2),JI(2,2),DETJ,AREA,THICK,GVAL,EPS1
      REAL(DOUBLE):: GR(3,2),DH(3,2),RGP(3),SGP(3),SUBR(3),SUBS(3),PW(3),BN(3)
      REAL(DOUBLE):: KF(20,20),KBBI(2,2),KBC(2,18),BMAP(2,18),KC(18,18)
      REAL(DOUBLE):: BM(3,20),BB(3,20),BS(2,20),BBS(3,20,3),BC(3,8),SC(2,11),DR(20)
      REAL(DOUBLE):: BT(3,18),ST(2,18),UNIT_TM(20),UNIT_TG(20),FTG(18),CTE(3),EPSM(3),NRES(3),TBAR,WT,KG
      LOGICAL,SAVE:: ANNOUNCED=.FALSE.
      EPS1=1.0D-30
      THICK=EPROP(1)
      IF(ELGP/=3.OR.THICK<=ZERO) THEN
       WRITE(F06,*) ' *ERROR: KMITC3+ expects three nodes and positive thickness, EID=',EID
       CALL OUTA_HERE('Y')
      ENDIF
      XYZ=XEB(1:3,1:3)
      EX=XYZ(2,:)-XYZ(1,:); WT=SQRT(SUM(EX*EX))
      IF(WT<=1.0D-14) THEN
       WRITE(F06,*) ' *ERROR: KMITC3+ degenerate first edge, EID=',EID
       CALL OUTA_HERE('Y')
      ENDIF
      EX=EX/WT
      VP=XYZ(3,:)-XYZ(1,:); EY=VP-DOT_PRODUCT(VP,EX)*EX
      IF(SUM(EY*EY)<1.0D-28) CALL OUTA_HERE('Y')
      EY=EY/SQRT(SUM(EY*EY))
      EZ=(/ EX(2)*EY(3)-EX(3)*EY(2),EX(3)*EY(1)-EX(1)*EY(3),EX(1)*EY(2)-EX(2)*EY(1) /)
      TE(1,:)=EX; TE(2,:)=EY; TE(3,:)=EZ
      XY=ZERO
      DO I=2,3
       XY(I,1)=DOT_PRODUCT(XYZ(I,:)-XYZ(1,:),EX)
       XY(I,2)=DOT_PRODUCT(XYZ(I,:)-XYZ(1,:),EY)
      ENDDO
      DETJ=XY(2,1)*XY(3,2); AREA=DETJ/TWO
      JI=ZERO; JI(1,1)=ONE/XY(2,1); JI(2,1)=-XY(3,1)/DETJ; JI(2,2)=ONE/XY(3,2)
      CALL MITC3P_DH(DH); GR=MATMUL(DH,TRANSPOSE(JI))
      RGP=(/ONE/6.0D0,TWO/THREE,ONE/6.0D0/);SGP=(/ONE/6.0D0,ONE/6.0D0,TWO/THREE/)
! Three-point L2 projection is interpolation through the three Gauss samples.
! This is exactly f^T F^-1 Q with F and Q evaluated by the specified three-point rule.
      DO G=1,3
       CALL MITC3P_BB_AT(RGP(G),SGP(G),JI,BC)
       DO I=1,3
        DO J=1,20
         BBS(I,J,G)=MITC3P_BB_TERM(BC,I,J)
        ENDDO
       ENDDO
      ENDDO
      KF=ZERO;UNIT_TM=ZERO;UNIT_TG=ZERO
      CTE=(/ALPVEC(1,1),ALPVEC(2,1),ALPVEC(4,1)/)
      DO GP=1,3
       CALL ALLMAN(RGP(GP),SGP(GP),BM)
       KF=KF+AREA/THREE*MATMUL(TRANSPOSE(BM),MATMUL(SHELL_A,BM))
       UNIT_TM=UNIT_TM+AREA/THREE*MATMUL(TRANSPOSE(BM),MATMUL(SHELL_A,CTE))
       CALL MITC3P_BS_LEE(RGP(GP),SGP(GP),XY,JI,SC)
       DO I=1,2
        DO J=1,20
         BS(I,J)=MITC3P_BS_TERM(SC,I,J)
        ENDDO
       ENDDO
       KF=KF+AREA/THREE*MATMUL(TRANSPOSE(BS),MATMUL(SHELL_T,BS))
      ENDDO
! Nguyen independent true-drill stabilizer K_Add, coefficient G*h*A.
      DR=ZERO
      DO I=1,3
       DR(6*I-5)=0.5D0*GR(I,2); DR(6*I-4)=-0.5D0*GR(I,1);DR(6*I)=ONE/THREE
      ENDDO
      GVAL=SHELL_T(1,1)/((5.0D0/6.0D0)*THICK)
      DO I=1,20
       DO J=1,20
        KF(I,J)=KF(I,J)+GVAL*THICK*AREA*DR(I)*DR(J)
       ENDDO
      ENDDO
      SUBR=(/4.0D0/9.0D0,4.0D0/9.0D0,ONE/9.0D0/)
      SUBS=(/ONE/9.0D0,4.0D0/9.0D0,4.0D0/9.0D0/)
      DO GP=1,3
       CALL PROJECT_BB(SUBR(GP),SUBS(GP),BB)
       KF=KF+AREA/THREE*MATMUL(TRANSPOSE(BB),MATMUL(SHELL_D,BB))
       UNIT_TG=UNIT_TG+AREA/THREE*MATMUL(TRANSPOSE(BB),SHELL_DALP)
      ENDDO
      KBC=KF(19:20,1:18);CALL MITC3P_INV2(KF(19:20,19:20),KBBI)
      BMAP=-MATMUL(KBBI,KBC)
      KC=KF(1:18,1:18)+MATMUL(KF(1:18,19:20),BMAP)
      IF(OPT(4)=='Y') KE(1:18,1:18)=0.5D0*(KC+TRANSPOSE(KC))
      IF(OPT(3)=='Y'.OR.OPT(6)=='Y') THEN
       CALL ALLMAN(ONE/THREE,ONE/THREE,BM)
       CALL PROJECT_BB(ONE/THREE,ONE/THREE,BB)
       CALL MITC3P_BS_LEE(ONE/THREE,ONE/THREE,XY,JI,SC)
       DO I=1,2
        DO J=1,20
         BS(I,J)=MITC3P_BS_TERM(SC,I,J)
        ENDDO
       ENDDO
       BT=BB(:,1:18)+MATMUL(BB(:,19:20),BMAP)
       ST=BS(:,1:18)+MATMUL(BS(:,19:20),BMAP)
       DO GP=1,ELGP
        BE1(1:3,1:18,GP)=BM(:,1:18);BE2(1:3,1:18,GP)=BT;BE3(1:2,1:18,GP)=ST
        SE1(1:3,1:18,GP)=MATMUL(EM,BM(:,1:18))
        SE2(1:3,1:18,GP)=MATMUL(EB,BT);SE3(1:2,1:18,GP)=MATMUL(ET(1:2,1:2),ST)
       ENDDO
       DO J=1,NTSUB
        TBAR=SUM(DT(1:3,J))/THREE-TREF(1)
        DO GP=1,ELGP
         STE1(1:3,J,GP)=MATMUL(EM,CTE)*TBAR
         STE2(1:3,J,GP)=MATMUL(EB,ALPVEC(1:3,2))*DT(4,J)
        ENDDO
       ENDDO
      ENDIF
      IF(OPT(2)=='Y') THEN
       FTG=UNIT_TG(1:18)+MATMUL(TRANSPOSE(BMAP),UNIT_TG(19:20))
       DO J=1,NTSUB
        TBAR=SUM(DT(1:3,J))/THREE-TREF(1)
        PTE(1:18,J)=UNIT_TM(1:18)*TBAR+FTG*DT(4,J)
       ENDDO
      ENDIF
      IF(OPT(5)=='Y') THEN
       DO J=1,NSUB
        DO I=1,3
         PPE(6*I-3,J)=AREA/THREE*PRESS(3,J)
        ENDDO
       ENDDO
      ENDIF
      IF(OPT(1)=='Y') THEN
       ME(1:18,1:18)=ZERO
       DO I=1,3
        DO J=1,3
         WT=ZERO
         IF(SOL_NAME(1:5)=='MODES'.AND.COUPMASS>0) THEN
          WT=MASS_PER_UNIT_AREA*AREA/12.0D0
          IF(I==J) WT=TWO*WT
         ELSE IF(I==J) THEN
          WT=MASS_PER_UNIT_AREA*AREA/THREE
         ENDIF
         DO K=1,3
          ME(6*(I-1)+K,6*(J-1)+K)=WT
         ENDDO
        ENDDO
       ENDDO
      ENDIF
      IF(OPT(6)=='Y') THEN
       CALL ELMDIS;KED(1:18,1:18)=ZERO
       DO GP=1,3
        CALL ALLMAN(RGP(GP),SGP(GP),BM)
        EPSM=MATMUL(BM(:,1:18),UEL(1:18));NRES=MATMUL(SHELL_A,EPSM)
        DO IA=1,3
         DO IB=1,3
          KG=AREA/THREE*(NRES(1)*GR(IA,1)*GR(IB,1)+NRES(2)*GR(IA,2)*GR(IB,2)+ &
           NRES(3)*(GR(IA,1)*GR(IB,2)+GR(IA,2)*GR(IB,1)))
          DO K=1,3
           KED(6*(IA-1)+K,6*(IB-1)+K)=KED(6*(IA-1)+K,6*(IB-1)+K)+KG
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDIF
      IF(.NOT.ANNOUNCED) THEN
       WRITE(F06,*) ' *INFORMATION: CTRIAR TRIARTYP=KMITC3+; Nguyen2024 research reconstruction EXPERIMENTAL/HOLD.'
       ANNOUNCED=.TRUE.
      ENDIF
      CONTAINS
      SUBROUTINE ALLMAN(R,S,B)
       REAL(DOUBLE),INTENT(IN)::R,S
       REAL(DOUBLE),INTENT(OUT)::B(3,20)
       REAL(DOUBLE):: N(3),H(3),BI,CI
       INTEGER::II,JJ,KK
       N=(/ONE-R-S,R,S/);B=ZERO
       DO II=1,3
        B(1,6*II-5)=GR(II,1);B(2,6*II-4)=GR(II,2)
        B(3,6*II-5)=GR(II,2);B(3,6*II-4)=GR(II,1)
       ENDDO
       DO II=1,3
        JJ=MOD(II,3)+1;KK=MOD(II+1,3)+1
        BI=XY(JJ,2)-XY(KK,2);CI=XY(KK,1)-XY(JJ,1)
        H(1)=BI*(GR(JJ,1)*N(KK)+GR(KK,1)*N(JJ))
        H(2)=CI*(GR(JJ,2)*N(KK)+GR(KK,2)*N(JJ))
        H(3)=(BI*GR(KK,2)+CI*GR(KK,1))*N(JJ)+(BI*GR(JJ,2)+CI*GR(JJ,1))*N(KK)
        B(:,6*JJ)=B(:,6*JJ)+(TWO/THREE)*H
        B(:,6*KK)=B(:,6*KK)-(TWO/THREE)*H
       ENDDO
      END SUBROUTINE
      SUBROUTINE PROJECT_BB(R,S,B)
       REAL(DOUBLE),INTENT(IN):: R,S
       REAL(DOUBLE),INTENT(OUT):: B(3,20)
       REAL(DOUBLE):: W(3)
       INTEGER:: GG
       W=(/5.0D0/THREE-TWO*R-TWO*S,TWO*R-ONE/THREE,TWO*S-ONE/THREE/)
       B=ZERO
       DO GG=1,3
        B=B+W(GG)*BBS(:,:,GG)
       ENDDO
      END SUBROUTINE
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

      SUBROUTINE MITC3P_BB_AT(R, S, JI, BBOUT)
         REAL(DOUBLE), INTENT(IN)  :: R, S, JI(2,2)
         REAL(DOUBLE), INTENT(OUT) :: BBOUT(3,8)
         REAL(DOUBLE) :: DH(3,2), DH_XY(3,2), DF4(2), DF4_XY(2), T
         INTEGER(LONG) :: II
         CALL MITC3P_DH(DH)
         DH_XY = MATMUL(DH, TRANSPOSE(JI))
         T = ONE - R - S
         DF4(1) = 27.0D0*S*(T - R)
         DF4(2) = 27.0D0*R*(T - S)
! Lee candidate uses the Python physical bubble gradient, not legacy transpose.
         DF4_XY = MATMUL(JI, DF4)
         BBOUT = ZERO
         DO II=1,3
            BBOUT(1,2*II  ) =  DH_XY(II,1)
            BBOUT(2,2*II-1) = -DH_XY(II,2)
            BBOUT(3,2*II  ) =  DH_XY(II,2)
            BBOUT(3,2*II-1) = -DH_XY(II,1)
         ENDDO
         BBOUT(1,8) =  DF4_XY(1)
         BBOUT(2,7) = -DF4_XY(2)
         BBOUT(3,8) =  DF4_XY(2)
         BBOUT(3,7) = -DF4_XY(1)
      END SUBROUTINE MITC3P_BB_AT

      SUBROUTINE MITC3P_COV_LEE(R,S,XYL,JI,ER,ES)
         REAL(DOUBLE), INTENT(IN) :: R,S,XYL(3,2),JI(2,2)
         REAL(DOUBLE), INTENT(OUT) :: ER(11),ES(11)
         REAL(DOUBLE) :: DH(3,2),GRAD(3,2),N(3),GX(11),GY(11),BUB
         INTEGER(LONG) :: II
         CALL MITC3P_DH(DH)
         GRAD=MATMUL(DH,TRANSPOSE(JI))
         N=(/ ONE-R-S,R,S /)
         GX=ZERO
         GY=ZERO
         DO II=1,3
            GX(3*II-2)=GRAD(II,1)
            GY(3*II-2)=GRAD(II,2)
            GX(3*II)=N(II)
            GY(3*II-1)=-N(II)
         ENDDO
         BUB=27.0D0*R*S*(ONE-R-S)
         GX(11)=BUB
         GY(10)=-BUB
         ER=XYL(2,1)*GX+XYL(2,2)*GY
         ES=XYL(3,1)*GX+XYL(3,2)*GY
      END SUBROUTINE MITC3P_COV_LEE

      SUBROUTINE MITC3P_BS_LEE(R,S,XYL,JI,BSOUT)
         REAL(DOUBLE), INTENT(IN) :: R,S,XYL(3,2),JI(2,2)
         REAL(DOUBLE), INTENT(OUT) :: BSOUT(2,11)
         REAL(DOUBLE), PARAMETER :: DIST=1.0D-4
         REAL(DOUBLE) :: ERA(11),ESA(11),ERB(11),ESB(11),ERC(11),ESC(11)
         REAL(DOUBLE) :: ERD(11),ESD(11),ERE(11),ESE(11),ERF(11),ESF(11),C(11),TMP(2,11)
         CALL MITC3P_COV_LEE(ONE/6.0D0,TWO/THREE,XYL,JI,ERA,ESA)
         CALL MITC3P_COV_LEE(TWO/THREE,ONE/6.0D0,XYL,JI,ERB,ESB)
         CALL MITC3P_COV_LEE(ONE/6.0D0,ONE/6.0D0,XYL,JI,ERC,ESC)
         CALL MITC3P_COV_LEE(ONE/THREE+DIST,ONE/THREE-TWO*DIST,XYL,JI,ERD,ESD)
         CALL MITC3P_COV_LEE(ONE/THREE-TWO*DIST,ONE/THREE+DIST,XYL,JI,ERE,ESE)
         CALL MITC3P_COV_LEE(ONE/THREE+DIST,ONE/THREE+DIST,XYL,JI,ERF,ESF)
         C=(ERF-ERD)-(ESF-ESE)
         TMP(1,:)=(TWO*ERB-ESB+ERC+ESC)/THREE+C*(THREE*S-ONE)/THREE
         TMP(2,:)=(TWO*ESA-ERA+ERC+ESC)/THREE+C*(ONE-THREE*R)/THREE
         BSOUT=MATMUL(JI,TMP)
      END SUBROUTINE MITC3P_BS_LEE

      END SUBROUTINE CTRIAR_KMITC3P
