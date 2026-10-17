      SUBROUTINE PYRAM ( OPT, INT_ELEM_ID, RED_INT_SHEAR, WRITE_WARN )

! PYRAM solid element kernel.
!
! PYRAM5 : standard linear five-node pyramid in both SOLIDTYP modes.
! PYRAM13: Bedrosian rational interpolation, eight-point stiffness.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NTSUB, MAX_ORDER_GAUSS, SOL_NAME, INT_SC_NUM, JTSUB
      USE CONSTANTS_1, ONLY           :  HALF, ZERO, ONE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  ALPVEC, BE1, BE2, DT, EID, ELGP, NUM_EMG_FATAL_ERRS, ES, KE, KED, ME, PTE, RHO,         &
                                         SE1, SE2, STE1, STE2, STRESS, TREF, TYPE, XEL, UEL, ISOLID, SUBLOD
      USE PARAMS, ONLY                :  EPSIL, COUPMASS
      USE PYRAM13_BEDROSIAN, ONLY     :  BEDROSIAN13_SHAPE, BEDROSIAN13_QUADRATURE, &
                                         BEDROSIAN13_RECOVER, BEDROSIAN13_KG

      USE PYRAM_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'PYRAM'
      CHARACTER( 1*BYTE), INTENT(IN)  :: RED_INT_SHEAR
      CHARACTER( 1*BYTE), INTENT(IN)  :: OPT(6)
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID
      INTEGER(LONG)                   :: I,J,K,L
      INTEGER(LONG)                   :: ID(3*ELGP)
      INTEGER(LONG)                   :: II,JJ
      INTEGER(LONG)                   :: NGP, KGORDER, KGERR, KGNP
      INTEGER(LONG)                   :: OUT_PT

      REAL(DOUBLE)                    :: GPSTRESS(6,27), GPSTRAIN(6,27), UTEMP(39), TNODES(13)
      REAL(DOUBLE)                    :: ALP(6)
      REAL(DOUBLE)                    :: B(6,3*ELGP)
      REAL(DOUBLE)                    :: CBAR(3,3*ELGP)
      REAL(DOUBLE)                    :: DETJ
      REAL(DOUBLE)                    :: DNX(3,ELGP)
      REAL(DOUBLE)                    :: DUM0(3*ELGP)
      REAL(DOUBLE)                    :: DUM3(3*ELGP,3*ELGP)
      REAL(DOUBLE)                    :: DUM4(6,3*ELGP)
      REAL(DOUBLE)                    :: DUM5(3*ELGP,3*ELGP)
      REAL(DOUBLE)                    :: DUM6(3,3*ELGP)
      REAL(DOUBLE)                    :: EALP(6)
      REAL(DOUBLE)                    :: EPS1
      REAL(DOUBLE)                    :: INT_R(64)
      REAL(DOUBLE)                    :: INT_S(64)
      REAL(DOUBLE)                    :: INT_T(64)
      REAL(DOUBLE)                    :: INT_W(64)
      REAL(DOUBLE)                    :: KWW(3,3)
      REAL(DOUBLE)                    :: M_1DOF(ELGP,ELGP)
      REAL(DOUBLE)                    :: N(ELGP)
      REAL(DOUBLE)                    :: RHO1
      REAL(DOUBLE)                    :: SIGxx
      REAL(DOUBLE)                    :: SIGxy
      REAL(DOUBLE)                    :: SIGyy
      REAL(DOUBLE)                    :: SIGyz
      REAL(DOUBLE)                    :: SIGzx
      REAL(DOUBLE)                    :: SIGzz
      REAL(DOUBLE)                    :: TBAR(NTSUB)
      REAL(DOUBLE)                    :: TEMP
      REAL(DOUBLE)                    :: TREF1
      REAL(DOUBLE)                    :: VOLUME
      REAL(DOUBLE)                    :: OUT_POINTS(3,6), SAMPLE_B(6,3*ELGP,8), SHAPE_DN(3,ELGP)
      REAL(DOUBLE)                    :: G, ZLO, ZHI, WEIGHT, ZWEIGHT

      IF ((RED_INT_SHEAR == 'Y') .AND. (WRITE_WARN == 'Y') .AND. (INT_ELEM_ID < 0)) THEN
         RETURN
      ENDIF

      EPS1 = EPSIL(1)
      TREF1 = TREF(1)

      DO I=1,6
         ALP(I) = ALPVEC(I,1)
      ENDDO
      CALL MATMULT_FFF ( ES, ALP, 6, 6, 1, EALP )

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         DO J=1,NTSUB
            TBAR(J) = ZERO
            DO K=1,ELGP
               TBAR(J) = TBAR(J) + DT(K,J)
            ENDDO
            TBAR(J) = TBAR(J)/ELGP - TREF1
         ENDDO
      ENDIF

      JJ = 0
      DO I=1,ELGP
         II = 6*(I - 1)
         DO J=1,3
            II = II + 1
            JJ = JJ + 1
            ID(JJ) = II
         ENDDO
      ENDDO

      IF (TYPE == 'PYRAM5  ') THEN
         CALL PYRAM5_POINTS ( NGP, INT_R, INT_S, INT_T, INT_W )
      ELSE IF (TYPE == 'PYRAM13 ') THEN
         CALL PYRAM13_POINTS ( NGP, INT_R, INT_S, INT_T, INT_W )
      ELSE
         RETURN
      ENDIF

      VOLUME = ZERO
      DO I=1,NGP
         CALL PYRAM_BMAT ( INT_R(I), INT_S(I), INT_T(I), DETJ, B, N )
         VOLUME = VOLUME + INT_W(I)*DETJ
      ENDDO

      IF (VOLUME < EPS1) THEN
         WRITE(ERR,1925) EID, TYPE, 'VOLUME', VOLUME
         WRITE(F06,1925) EID, TYPE, 'VOLUME', VOLUME
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

      IF (OPT(1) == 'Y') THEN
         RHO1 = RHO(1)
         IF (RHO1 /= ZERO) THEN
            IF ((TYPE == 'PYRAM13 ') .AND. (SOL_NAME(1:5) == 'MODES') .AND. (COUPMASS <= 0)) THEN
               WRITE(ERR,'(A)') ' *ERROR: PYRAM13 MODAL MASS REQUIRES PARAM,COUPMASS,1; ROW-SUM MASS IS INDEFINITE'
               WRITE(F06,'(A)') ' *ERROR: PYRAM13 MODAL MASS REQUIRES PARAM,COUPMASS,1; ROW-SUM MASS IS INDEFINITE'
               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               FATAL_ERR = FATAL_ERR + 1
               RETURN
            ENDIF
! Stiffness underintegration must not be reused for quadratic mass: eight
! samples give only rank 24/39. Integrate mass independently, then restore
! the validated eight-point stiffness/thermal/recovery quadrature.
            IF (TYPE == 'PYRAM13 ') CALL PYRAM13_MASS_POINTS ( NGP, INT_R, INT_S, INT_T, INT_W )
            M_1DOF(:,:) = ZERO
            DO I=1,NGP
               CALL PYRAM_BMAT ( INT_R(I), INT_S(I), INT_T(I), DETJ, B, N )
               DO J=1,ELGP
                  DO K=1,ELGP
                     M_1DOF(J,K) = M_1DOF(J,K) + RHO1*INT_W(I)*DETJ*N(J)*N(K)
                  ENDDO
               ENDDO
            ENDDO
            CALL EXPAND_MASS_DOFS ( M_1DOF )
            IF (TYPE == 'PYRAM13 ') CALL PYRAM13_POINTS ( NGP, INT_R, INT_S, INT_T, INT_W )
         ENDIF
      ENDIF

      IF (OPT(2) == 'Y') THEN
         DO L=1,NTSUB
            DO J=1,3*ELGP
               DUM0(J) = ZERO
            ENDDO
            DO I=1,NGP
               CALL PYRAM_BMAT ( INT_R(I), INT_S(I), INT_T(I), DETJ, B, N )
               CALL MATMULT_FFF_T ( B, EALP, 6, 3*ELGP, 1, DUM3 )
               TEMP = ZERO
               DO J=1,ELGP
                  TEMP = TEMP + N(J)*DT(J,L)
               ENDDO
               TEMP = TEMP - TREF1
               DO J=1,3*ELGP
                  DUM0(J) = DUM0(J) + DUM3(J,1)*TEMP*INT_W(I)*DETJ
               ENDDO
            ENDDO
            DO J=1,3*ELGP
               PTE(ID(J),L) = DUM0(J)
            ENDDO
         ENDDO
      ENDIF

      IF (OPT(4) == 'Y') THEN
! Standard PYRAM5 is deliberately retained under NEWSOLID. Experimental
! EAS54 is not selected. PYRAM13 has no rank-one stabilization.
            DUM3(:,:) = ZERO
            DO I=1,NGP
               CALL PYRAM_BMAT ( INT_R(I), INT_S(I), INT_T(I), DETJ, B, N )
               CALL MATMULT_FFF ( ES, B, 6, 6, 3*ELGP, DUM4 )
               CALL MATMULT_FFF_T ( B, DUM4, 6, 3*ELGP, 3*ELGP, DUM5 )
               DO J=1,3*ELGP
                  DO K=1,3*ELGP
                     DUM3(J,K) = DUM3(J,K) + DUM5(J,K)*INT_W(I)*DETJ
                  ENDDO
               ENDDO
            ENDDO
            DO J=1,3*ELGP
               DO K=1,3*ELGP
                  KE(ID(J),ID(K)) = DUM3(J,K)
               ENDDO
            ENDDO
         DO I=2,6*ELGP
            DO J=1,I-1
               KE(I,J) = KE(J,I)
            ENDDO
         ENDDO
      ENDIF

      IF ((OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN
! Center plus five corner tensors. Extrapolate the valid eight-point field
! in collapsed coordinates; the apex itself has no invertible Jacobian.
         OUT_POINTS=RESHAPE([ZERO,ZERO,0.25D0, -ONE,-ONE,ZERO, ONE,-ONE,ZERO, &
                             ONE,ONE,ZERO, -ONE,ONE,ZERO, ZERO,ZERO,ONE],[3,6])
         DO I=1,NGP
            CALL PYRAM_BMAT(INT_R(I),INT_S(I),INT_T(I),DETJ,SAMPLE_B(:,:,I),N)
         ENDDO
         IF (NUM_EMG_FATAL_ERRS > 0) RETURN
         G=ONE/DSQRT(3.0D0); ZLO=MINVAL(INT_T(1:NGP)); ZHI=MAXVAL(INT_T(1:NGP))
         DO OUT_PT=1,6
            IF (OUT_PT==1) THEN
               CALL PYRAM_BMAT(ZERO,ZERO,0.25D0,DETJ,B,N)
            ELSE
               B=ZERO
               DO I=1,NGP
                  IF (DABS(INT_T(I)-ZLO)<EPS1) THEN
                     ZWEIGHT=(ZHI-OUT_POINTS(3,OUT_PT))/(ZHI-ZLO)
                  ELSE
                     ZWEIGHT=(OUT_POINTS(3,OUT_PT)-ZLO)/(ZHI-ZLO)
                  ENDIF
                  WEIGHT=(ONE+OUT_POINTS(1,OUT_PT)*INT_R(I)/(G*G))* &
                         (ONE+OUT_POINTS(2,OUT_PT)*INT_S(I)/(G*G))*ZWEIGHT/4.0D0
                  B=B+WEIGHT*SAMPLE_B(:,:,I)
               ENDDO
               CALL PYRAM_SHAPE(OUT_POINTS(1,OUT_PT),OUT_POINTS(2,OUT_PT),OUT_POINTS(3,OUT_PT),N,SHAPE_DN)
            ENDIF
            CALL MATMULT_FFF ( ES, B, 6, 6, 3*ELGP, DUM4 )
            DO I=1,3
               DO J=1,3*ELGP
                  BE1(I,ID(J),OUT_PT) = B(I  ,J)
                  BE2(I,ID(J),OUT_PT) = B(I+3,J)
                  SE1(I,ID(J),OUT_PT) = DUM4(I  ,J)
                  SE2(I,ID(J),OUT_PT) = DUM4(I+3,J)
               ENDDO
            ENDDO
            DO L=1,NTSUB
               TEMP=DOT_PRODUCT(N,DT(1:ELGP,L))-TREF1
               DO I=1,3
                  STE1(I,L,OUT_PT) = EALP(I)*TEMP
                  STE2(I,L,OUT_PT) = EALP(I+3)*TEMP
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN
         CALL ELMDIS
         IF (TYPE == 'PYRAM13 ') THEN
! PSOLID IN=2/3 selects only initial-stress quadrature here; K remains K8.
            KGORDER = 2
            IF (ISOLID(4) == 3) KGORDER = 3
            DO I=1,39
               UTEMP(I) = UEL(ID(I))
            ENDDO
            TNODES = TREF1
            IF (SUBLOD(INT_SC_NUM,2) > 0) THEN
               IF ((JTSUB < 1) .OR. (JTSUB > NTSUB)) THEN
                  WRITE(ERR,'(A,I8)') ' *ERROR: PYRAM13 INVALID THERMAL SUBCASE FOR KG, ELEMENT ', EID
                  WRITE(F06,'(A,I8)') ' *ERROR: PYRAM13 INVALID THERMAL SUBCASE FOR KG, ELEMENT ', EID
                  FATAL_ERR = FATAL_ERR + 1
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                  RETURN
               ENDIF
               TNODES = DT(1:13,JTSUB)
            ENDIF
            CALL BEDROSIAN13_RECOVER(XEL(1:13,1:3),ES,UTEMP,ALP,TNODES,TREF1,KGORDER,GPSTRESS,GPSTRAIN,KGNP,KGERR)
            IF (KGERR == 0) CALL BEDROSIAN13_KG(XEL(1:13,1:3),GPSTRESS,KGORDER,DUM3,KGERR)
            IF (KGERR /= 0) THEN
               WRITE(ERR,'(A,I8)') ' *ERROR: PYRAM13 INVALID INTEGRATION POINT FOR KG, ELEMENT ', EID
               WRITE(F06,'(A,I8)') ' *ERROR: PYRAM13 INVALID INTEGRATION POINT FOR KG, ELEMENT ', EID
               FATAL_ERR = FATAL_ERR + 1
               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               RETURN
            ENDIF
         ELSE
         CALL ELEM_STRE_STRN_ARRAYS ( 1 )

         SIGxx = STRESS(1)
         SIGyy = STRESS(2)
         SIGzz = STRESS(3)
         SIGxy = STRESS(4)
         SIGyz = STRESS(5)
         SIGzx = STRESS(6)

         KWW(1,1) = SIGyy + SIGzz
         KWW(1,2) = -SIGxy
         KWW(1,3) = -SIGzx
         KWW(2,2) = SIGxx + SIGzz
         KWW(2,3) = -SIGyz
         KWW(3,3) = SIGxx + SIGyy
         KWW(2,1) = KWW(1,2)
         KWW(3,1) = KWW(1,3)
         KWW(3,2) = KWW(2,3)

         DUM3(:,:) = ZERO
         DO I=1,NGP
            CALL PYRAM_BMAT ( INT_R(I), INT_S(I), INT_T(I), DETJ, B, N, DNX )
            CBAR(:,:) = ZERO
            DO L=1,ELGP
               CBAR(1,3*(L-1)+1) =  ZERO
               CBAR(1,3*(L-1)+2) = -HALF*DNX(3,L)
               CBAR(1,3*(L-1)+3) =  HALF*DNX(2,L)
               CBAR(2,3*(L-1)+1) =  HALF*DNX(3,L)
               CBAR(2,3*(L-1)+2) =  ZERO
               CBAR(2,3*(L-1)+3) = -HALF*DNX(1,L)
               CBAR(3,3*(L-1)+1) = -HALF*DNX(2,L)
               CBAR(3,3*(L-1)+2) =  HALF*DNX(1,L)
               CBAR(3,3*(L-1)+3) =  ZERO
            ENDDO
            CALL MATMULT_FFF ( KWW, CBAR, 3, 3, 3*ELGP, DUM6 )
            CALL MATMULT_FFF_T ( CBAR, DUM6, 3, 3*ELGP, 3*ELGP, DUM5 )
            DO J=1,3*ELGP
               DO K=1,3*ELGP
                  DUM3(J,K) = DUM3(J,K) + DUM5(J,K)*INT_W(I)*DETJ
               ENDDO
            ENDDO
         ENDDO

         ENDIF

         DO I=1,3*ELGP
            DO J=1,3*ELGP
               KED(ID(I),ID(J)) = DUM3(I,J)
            ENDDO
         ENDDO

         DO I=2,6*ELGP
            DO J=1,I-1
               KED(I,J) = KED(J,I)
            ENDDO
         ENDDO

         IF (DEBUG(178) >= 1) THEN
            WRITE(F06,2001) TRIM(SUBR_NAME), OPT(6), LOAD_ISTEP, TYPE, EID
         ENDIF
      ENDIF

      RETURN

 1925 FORMAT(' *ERROR  1925: ELEMENT ',I8,', TYPE ',A,' HAS NONPOSITIVE ',A,' = ',1ES14.6)
 2001 FORMAT(' In ', A, ' with OPT(6) = ', A, ', and LOAD_ISTEP = ', I8, ': KED Differential Stiffness Matrix for ', A,         &
                    ' element ', I8)

      CONTAINS

      SUBROUTINE PYRAM5_POINTS ( NPT, RR, SS, TT, WW )

      INTEGER(LONG), INTENT(OUT)      :: NPT
      REAL(DOUBLE),  INTENT(OUT)      :: RR(:),SS(:),TT(:),WW(:)

      INTEGER(LONG)                   :: A,B,C
      REAL(DOUBLE)                    :: G(2),GT(2),W(2),WT(2)

      G(1)  = -ONE/DSQRT(3.0D0)
      G(2)  =  ONE/DSQRT(3.0D0)
      W(1)  =  ONE
      W(2)  =  ONE
      GT(1) =  0.5D0*(ONE - ONE/DSQRT(3.0D0))
      GT(2) =  0.5D0*(ONE + ONE/DSQRT(3.0D0))
      WT(1) =  0.5D0
      WT(2) =  0.5D0

      NPT = 0
      DO A=1,2
         DO B=1,2
            DO C=1,2
               NPT = NPT + 1
               RR(NPT) = G(A)
               SS(NPT) = G(B)
               TT(NPT) = GT(C)
               WW(NPT) = W(A)*W(B)*WT(C)
            ENDDO
         ENDDO
      ENDDO

      END SUBROUTINE PYRAM5_POINTS

      SUBROUTINE PYRAM13_POINTS ( NPT, RR, SS, TT, WW )
! 2x2 Legendre x 2 Jacobi(alpha=2,beta=0), as test_v7.assemble.
! Jacobi weights already include (1-z)^2; divide it out because
! PYRAM_BMAT uses the collapsed-coordinate Jacobian containing that factor.
      INTEGER(LONG), INTENT(OUT) :: NPT
      REAL(DOUBLE), INTENT(OUT) :: RR(:),SS(:),TT(:),WW(:)
      REAL(DOUBLE) :: G(2), Z(2), WZ(2), SQ10
      INTEGER(LONG) :: A,B,C
      G = (/-ONE,ONE/)/DSQRT(3.0D0)
      SQ10 = DSQRT(10.0D0)
      Z = (/5.0D0-SQ10,5.0D0+SQ10/)/15.0D0
      WZ = (/ONE/6.0D0+SQ10/48.0D0,ONE/6.0D0-SQ10/48.0D0/)
      NPT = 0
      DO A=1,2
         DO B=1,2
            DO C=1,2
               NPT = NPT + 1
               RR(NPT) = G(A)
               SS(NPT) = G(B)
               TT(NPT) = Z(C)
               WW(NPT) = WZ(C)/(ONE-Z(C))**2
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE PYRAM13_POINTS

      SUBROUTINE PYRAM13_MASS_POINTS ( NPT, RR, SS, TT, WW )
! Independent 3x3x3 Jacobi quadrature, v3 Python mass_order=3.
      INTEGER(LONG), INTENT(OUT) :: NPT
      REAL(DOUBLE), INTENT(OUT) :: RR(:),SS(:),TT(:),WW(:)
      REAL(DOUBLE) :: P(3,27), W(27)
      INTEGER(LONG) :: IERR
      CALL BEDROSIAN13_QUADRATURE(3,P,W,NPT,IERR)
      RR(1:NPT)=P(1,1:NPT); SS(1:NPT)=P(2,1:NPT); TT(1:NPT)=P(3,1:NPT); WW(1:NPT)=W(1:NPT)
      END SUBROUTINE PYRAM13_MASS_POINTS

      SUBROUTINE PYRAM_SHAPE ( R, S, T, N, DN )

      REAL(DOUBLE), INTENT(IN)        :: R,S,T
      REAL(DOUBLE), INTENT(OUT)       :: N(ELGP), DN(3,ELGP)

      IF (TYPE == 'PYRAM5  ') THEN
         CALL PYRAM5_SHAPE ( R, S, T, N, DN )
      ELSE IF (TYPE == 'PYRAM13 ') THEN
         CALL BEDROSIAN13_SHAPE ( R, S, T, N, DN )
      ENDIF

      END SUBROUTINE PYRAM_SHAPE

      SUBROUTINE PYRAM5_SHAPE ( R, S, T, N, DN )

      REAL(DOUBLE), INTENT(IN)        :: R,S,T
      REAL(DOUBLE), INTENT(OUT)       :: N(ELGP), DN(3,ELGP)
      REAL(DOUBLE)                    :: Q

      Q = ONE - T
      N(1) = 0.25D0*(ONE - R)*(ONE - S)*Q
      N(2) = 0.25D0*(ONE + R)*(ONE - S)*Q
      N(3) = 0.25D0*(ONE + R)*(ONE + S)*Q
      N(4) = 0.25D0*(ONE - R)*(ONE + S)*Q
      N(5) = T

      DN(:,:) = ZERO
      DN(1,1) = -0.25D0*(ONE - S)*Q
      DN(1,2) =  0.25D0*(ONE - S)*Q
      DN(1,3) =  0.25D0*(ONE + S)*Q
      DN(1,4) = -0.25D0*(ONE + S)*Q
      DN(2,1) = -0.25D0*(ONE - R)*Q
      DN(2,2) = -0.25D0*(ONE + R)*Q
      DN(2,3) =  0.25D0*(ONE + R)*Q
      DN(2,4) =  0.25D0*(ONE - R)*Q
      DN(3,1) = -0.25D0*(ONE - R)*(ONE - S)
      DN(3,2) = -0.25D0*(ONE + R)*(ONE - S)
      DN(3,3) = -0.25D0*(ONE + R)*(ONE + S)
      DN(3,4) = -0.25D0*(ONE - R)*(ONE + S)
      DN(3,5) =  ONE

      END SUBROUTINE PYRAM5_SHAPE

      SUBROUTINE PYRAM_BMAT ( R, S, T, DETJ, B, N, DNX_OUT )

      REAL(DOUBLE), INTENT(IN)            :: R,S,T
      REAL(DOUBLE), INTENT(OUT)           :: DETJ
      REAL(DOUBLE), INTENT(OUT)           :: B(6,3*ELGP)
      REAL(DOUBLE), INTENT(OUT)           :: N(ELGP)
      REAL(DOUBLE), INTENT(OUT), OPTIONAL :: DNX_OUT(3,ELGP)

      INTEGER(LONG)                   :: A,C,P,Q
      REAL(DOUBLE)                    :: DN(3,ELGP)
      REAL(DOUBLE)                    :: DNX(3,ELGP)
      REAL(DOUBLE)                    :: INVJ(3,3)
      REAL(DOUBLE)                    :: JAC(3,3)

      CALL PYRAM_SHAPE ( R, S, T, N, DN )

      JAC(:,:) = ZERO
      DO A=1,ELGP
         DO P=1,3
            DO Q=1,3
               JAC(P,Q) = JAC(P,Q) + DN(P,A)*XEL(A,Q)
            ENDDO
         ENDDO
      ENDDO

      CALL INV3 ( JAC, DETJ, INVJ )
      IF (DETJ <= EPS1) THEN
         WRITE(ERR,2925) EID, TYPE, 'JACOBIAN', DETJ
         WRITE(F06,2925) EID, TYPE, 'JACOBIAN', DETJ
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         B(:,:) = ZERO
         RETURN
      ENDIF

      DO A=1,ELGP
         DO P=1,3
            DNX(P,A) = ZERO
            DO Q=1,3
               DNX(P,A) = DNX(P,A) + INVJ(P,Q)*DN(Q,A)
            ENDDO
         ENDDO
      ENDDO

      IF (PRESENT(DNX_OUT)) DNX_OUT(:,:) = DNX(:,:)

      B(:,:) = ZERO
      DO A=1,ELGP
         C = 3*(A-1)
         B(1,C+1) = DNX(1,A)
         B(2,C+2) = DNX(2,A)
         B(3,C+3) = DNX(3,A)
         B(4,C+1) = DNX(2,A)
         B(4,C+2) = DNX(1,A)
         B(5,C+2) = DNX(3,A)
         B(5,C+3) = DNX(2,A)
         B(6,C+1) = DNX(3,A)
         B(6,C+3) = DNX(1,A)
      ENDDO

 2925 FORMAT(' *ERROR  1925: ELEMENT ',I8,', TYPE ',A,' HAS NONPOSITIVE ',A,' = ',1ES14.6)

      END SUBROUTINE PYRAM_BMAT

      SUBROUTINE INV3 ( A, DET, AI )

      REAL(DOUBLE), INTENT(IN)        :: A(3,3)
      REAL(DOUBLE), INTENT(OUT)       :: DET,AI(3,3)

      DET = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) +                     &
            A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
      IF (DABS(DET) <= EPS1) THEN
         AI(:,:) = ZERO
         RETURN
      ENDIF
      AI(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))/DET
      AI(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))/DET
      AI(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2))/DET
      AI(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))/DET
      AI(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))/DET
      AI(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))/DET
      AI(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1))/DET
      AI(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))/DET
      AI(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))/DET

      END SUBROUTINE INV3

      END SUBROUTINE PYRAM
