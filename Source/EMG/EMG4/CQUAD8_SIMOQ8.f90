! #################################################################################################################################
! CQUAD8 Python-aligned Simo1993 Q8 shell for PARAM,QUAD8TYP,SIMOQ8.

      SUBROUTINE CQUAD8_SIMOQ8 ( OPT, INT_ELEM_ID )

! Ported from:
!   C:/PROJECTAI/18a/python/Simo1993_Q8_ShellElement_v12_standalone.py
!   D:\18a\python\quadratic\Simo1993_Q8_thermal_buckling.py
!
! Static stiffness path:
!   Q8 serendipity Simo/Fox director kinematics + Hughes-Brezzi drilling,
!   3x3 membrane/bending/drilling/shear. v12 auto selects covariant ANS
!   for curved/distorted elements; otherwise one EAS shear bubble
!   phi=(1-r^2)(1-s^2) is condensed at element level.
!
! This branch is the active Simo Q8 path and follows the Python benchmark
! conventions more closely. Flat membrane geometric stiffness recovers
! active prestress at each integration point in basic coordinates.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_STRESS_POINTS, SOL_NAME, NSNORM
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  ALPVEC, BGRID, DT, EID, ELGP, GRID_SNORM, GRID_ID, SNORM, KE, KED, ME, BE1, BE2, BE3, EPROP, FCONV,    &
                                         Q8_POINT_BASIS, MASS_PER_UNIT_AREA, NUM_EMG_FATAL_ERRS, PCOMP_PROPS, PPE, PRESS, PTE, SHELL_A,        &
                                         SHELL_D, SHELL_T, TREF, UEL, XEB
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO
      USE PARAMS, ONLY                :  COUPMASS
      USE ELMDIS_Interface
      USE OUTA_HERE_Interface

      USE QUADRATIC_SURFACE_MASS_Interface
      USE QUADRATIC_SURFACE_PRESSURE_Interface
      USE MODEL_STUF, ONLY : UEB, KED
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CQUAD8_SIMOQ8'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      INTEGER(LONG)                   :: I, J, GP, IA, IB, K, L, JSUB
      REAL(DOUBLE)                    :: XYZ(8,3), NORMALS(8,3)
      REAL(DOUBLE)                    :: KDD(48,48), KOUT(48,48), KDA(48), KAA
      REAL(DOUBLE)                    :: BM(3,48), BB(3,48), BS(2,48), BD(1,48), BSE(2)
      REAL(DOUBLE)                    :: GP3(3), W3(3), R, S, WT, JAC, CDRILL, FAC
      REAL(DOUBLE)                    :: M1(8,8), N8(8), DN8(2,8), MASS_ELEM, MASS_NODE
      REAL(DOUBLE)                    :: UNIT_PPE(48), UNIT_PTE(48), DXDR(3), DXDS(3), SURF_VEC(3), TBAR
      REAL(DOUBLE)                    :: CTE(3), THERMAL_RESULTANT(3)
      REAL(DOUBLE)                    :: DLOC(2,8), SIG0(2,2), DUM28(2,8), KG8(8,8), STRAIN0(3), N0(3)
      REAL(DOUBLE)                    :: G1(3), G2(3), A11, A22, A12, DET, AI11, AI22, AI12
      REAL(DOUBLE) :: RR(9),SS(9),E1OUT(3),E2OUT(3),E3OUT(3)
      INTEGER(LONG)                   :: KI, KJ
      REAL(DOUBLE) :: NORMAL_SIGN
      LOGICAL :: USE_ANS

      REAL(DOUBLE) :: UNIT_PTG(48), TEMP_GRAD_SIGN, TEMP_NORMAL(3), TEMP_G1(3), TEMP_G2(3)

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
         WRITE(ERR,*) ' *ERROR: Code not written for composite material with PARAM,QUAD8TYP,SIMOQ8'
         WRITE(F06,*) ' *ERROR: Code not written for composite material with PARAM,QUAD8TYP,SIMOQ8'
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL LOAD_BASIC_COORDS_Q8 ( XYZ )
      CALL SET_NORMAL_SIGN_Q8(XYZ)
      CALL CALC_NODAL_NORMALS_Q8 ( XYZ, NORMALS )
      CALL SELECT_ANS_Q8(XYZ, NORMALS)

      GP3 = (/-DSQRT(3.0D0/5.0D0), ZERO, DSQRT(3.0D0/5.0D0)/)
      W3  = (/5.0D0/9.0D0, 8.0D0/9.0D0, 5.0D0/9.0D0/)

      IF (OPT(1) == 'Y') THEN
         CALL QUADRATIC_SURFACE_MASS(8,XYZ,(/(ZERO,J=1,8)/))
      ENDIF

      IF (OPT(2) == 'Y') THEN
! TEMPP1 gradient is along connectivity normal; directors use the normalized frame.
         TEMP_G1=XYZ(2,:)-XYZ(1,:)
         TEMP_G2=XYZ(3,:)-XYZ(1,:)
         TEMP_NORMAL(1)=TEMP_G1(2)*TEMP_G2(3)-TEMP_G1(3)*TEMP_G2(2)
         TEMP_NORMAL(2)=TEMP_G1(3)*TEMP_G2(1)-TEMP_G1(1)*TEMP_G2(3)
         TEMP_NORMAL(3)=TEMP_G1(1)*TEMP_G2(2)-TEMP_G1(2)*TEMP_G2(1)
         TEMP_GRAD_SIGN=SIGN(ONE,DOT_PRODUCT(TEMP_NORMAL,NORMALS(1,:)))
         CTE=(/ALPVEC(1,1),ALPVEC(2,1),ALPVEC(4,1)/)
         UNIT_PTE=ZERO
         UNIT_PTG=ZERO
         DO I=1,3
            DO J=1,3
               R=GP3(I)
               S=GP3(J)
               WT=W3(I)*W3(J)
               CALL BM_Q8_AT(XYZ,R,S,BM,JAC)
               CALL BB_Q8_AT(XYZ,NORMALS,R,S,BB,JAC)
               UNIT_PTE=UNIT_PTE+MATMUL(TRANSPOSE(BM),MATMUL(SHELL_A,CTE))*WT*JAC
! eps(z)=Bm*u-z*Bb*u: bending thermal force has the physical minus sign.
               UNIT_PTG=UNIT_PTG-MATMUL(TRANSPOSE(BB),MATMUL(SHELL_D,CTE))*WT*JAC
            ENDDO
         ENDDO
         DO JSUB=1,SIZE(PTE,2)
            TBAR=SUM(DT(1:ELGP,JSUB))/REAL(ELGP,DOUBLE)-TREF(1)
            PTE(1:48,JSUB)=UNIT_PTE*TBAR+UNIT_PTG*TEMP_GRAD_SIGN*DT(ELGP+1,JSUB)
         ENDDO
      ENDIF

      IF (OPT(3) == 'Y') THEN
         RR=(/ZERO,-ONE,ONE,ONE,-ONE,ZERO,ONE,ZERO,-ONE/)
         SS=(/ZERO,-ONE,-ONE,ONE,ONE,-ONE,ZERO,ONE,ZERO/)
         DO GP=1,9
            R=RR(GP)
            S=SS(GP)
            CALL BM_Q8_AT(XYZ,R,S,BM,JAC)
            CALL BB_Q8_AT(XYZ,NORMALS,R,S,BB,JAC)
            CALL BS_Q8_AT(XYZ,NORMALS,R,S,BS,JAC)
            BE1(1:3,1:48,GP)=BM
! Physical fiber strain is membrane-z*Bb*u; Bb is +Hessian(w).
            BE2(1:3,1:48,GP)=BB
            BE3(1:2,1:48,GP)=BS
            CALL LOCAL_BASIS_AT_Q8(XYZ,R,S,E1OUT,E2OUT,E3OUT,JAC)
            Q8_POINT_BASIS(1,:,GP)=E1OUT
            Q8_POINT_BASIS(2,:,GP)=E2OUT
            Q8_POINT_BASIS(3,:,GP)=E3OUT
         ENDDO
      ENDIF

      IF (OPT(4) == 'Y') THEN
         KDD = ZERO
         KDA = ZERO
         KAA = ZERO


         CDRILL = 2.0D-2*SHELL_T(1,1)/(5.0D0/6.0D0)

         DO I=1,3
            DO J=1,3
               R = GP3(I)
               S = GP3(J)
               WT = W3(I)*W3(J)

               CALL BM_Q8_AT ( XYZ, R, S, BM, JAC )
               CALL BB_Q8_AT ( XYZ, NORMALS, R, S, BB, JAC )
               CALL BS_Q8_AT ( XYZ, NORMALS, R, S, BS, JAC )
               CALL BDRILL_Q8_AT ( XYZ, NORMALS, R, S, BD, JAC )
               CALL EAS1_SHEAR_Q8_AT ( XYZ, R, S, BSE )


               KDD = KDD + WT*JAC*MATMUL(TRANSPOSE(BM), MATMUL(SHELL_A, BM))
               KDD = KDD + WT*JAC*MATMUL(TRANSPOSE(BB), MATMUL(SHELL_D, BB))
               KDD = KDD + WT*JAC*MATMUL(TRANSPOSE(BS), MATMUL(SHELL_T, BS))
               KDD = KDD + WT*JAC*CDRILL*MATMUL(TRANSPOSE(BD), BD)

               KDA = KDA + WT*JAC*MATMUL(TRANSPOSE(BS), MATMUL(SHELL_T, BSE))
               KAA = KAA + WT*JAC*DOT_PRODUCT(BSE, MATMUL(SHELL_T, BSE))
            ENDDO
         ENDDO

         KOUT = KDD
         IF (.NOT.USE_ANS .AND. DABS(KAA) > 1.0D-14) THEN
            DO IA=1,48
               DO IB=1,48
                  KOUT(IA,IB) = KOUT(IA,IB) - KDA(IA)*KDA(IB)/KAA
               ENDDO
            ENDDO
         ENDIF


         DO IA=1,48
            DO IB=1,48
               FAC = 0.5D0*(KOUT(IA,IB) + KOUT(IB,IA))
               KE(IA,IB) = FAC
            ENDDO
         ENDDO
      ENDIF

      IF (OPT(5) == 'Y') THEN
         CALL QUADRATIC_SURFACE_PRESSURE(INT_ELEM_ID,8,XYZ,(/(ZERO,J=1,8)/))
      ENDIF

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN
         CALL NATIVE_MEMBRANE_KG
      ENDIF

      RETURN

 9001 FORMAT(' *ERROR: ',A,' expects ELGP=8 for element ',I8,' but got ',I8)

      CONTAINS

! Flat-shell membrane initial-stress stiffness, tension-positive resultants.
! Four-point Gauss/Duffy; standard geometry area, active translation field.
! Mechanical linear reference state only; no director or follower tangent.
      SUBROUTINE NATIVE_MEMBRANE_KG
      REAL(DOUBLE) :: GX(4),GW(4),RG,SG,WG,NVAL(8),DG(2,8),DF(2,8)
      REAL(DOUBLE) :: TG(2,3),TF(2,3),CV(3),AREA,AJ,MT(2,2),INV(2,2),DETMT
      REAL(DOUBLE) :: E1(3),E2(3),E3(3),GRAD(2,8),BMG(3,48),NV(3),SIG(2,2),BLOCK(8,8)
      REAL(DOUBLE) :: SCALE_GEOM,NORMAL(3)
      INTEGER(LONG) :: IG,JG,IN,JN,ID
      GX=(/-0.8611363115940526D0,-0.3399810435848563D0, &
            0.3399810435848563D0,0.8611363115940526D0/)
      GW=(/0.3478548451374538D0,0.6521451548625461D0, &
           0.6521451548625461D0,0.3478548451374538D0/)
      TG(1,:)=XYZ(2,:)-XYZ(1,:)
      TG(2,:)=XYZ(3,:)-XYZ(1,:)
      NORMAL=(/TG(1,2)*TG(2,3)-TG(1,3)*TG(2,2), &
                TG(1,3)*TG(2,1)-TG(1,1)*TG(2,3), &
                TG(1,1)*TG(2,2)-TG(1,2)*TG(2,1)/)
      AREA=SQRT(SUM(NORMAL*NORMAL))
      SCALE_GEOM=MAXVAL(ABS(XYZ-SPREAD(XYZ(1,:),1,8)))
      IF (AREA <= 1.0D-14) THEN
         CALL KG_GEOMETRY_ERROR
      ENDIF
      NORMAL=NORMAL/AREA
      DO IN=1,8
         IF (ABS(DOT_PRODUCT(XYZ(IN,:)-XYZ(1,:),NORMAL)) > 1.0D-9*MAX(SCALE_GEOM,ONE)) THEN
            CALL KG_GEOMETRY_ERROR
         ENDIF
      ENDDO
      CALL ELMDIS
      BLOCK=ZERO
      DO IG=1,4
         DO JG=1,4
               RG=GX(IG)
               SG=GX(JG)
               WG=GW(IG)*GW(JG)
               CALL SHAPE_Q8(RG,SG,NVAL,DG)
               TG=MATMUL(DG,XYZ)
               CV=(/TG(1,2)*TG(2,3)-TG(1,3)*TG(2,2), &
                     TG(1,3)*TG(2,1)-TG(1,1)*TG(2,3), &
                     TG(1,1)*TG(2,2)-TG(1,2)*TG(2,1)/)
               AREA=SQRT(SUM(CV*CV))
               DF=DG
               TF=MATMUL(DF,XYZ)
               MT=MATMUL(TF,TRANSPOSE(TF))
               DETMT=MT(1,1)*MT(2,2)-MT(1,2)*MT(2,1)
               IF (AREA <= 1.0D-14 .OR. DETMT <= 1.0D-30) CALL KG_GEOMETRY_ERROR
               INV(1,1)=MT(2,2)/DETMT
               INV(2,2)=MT(1,1)/DETMT
               INV(1,2)=-MT(1,2)/DETMT
               INV(2,1)=INV(1,2)
               CALL LOCAL_BASIS_AT_Q8(XYZ,RG,SG,E1,E2,E3,AJ)
               GRAD(1,:)=MATMUL(MATMUL(INV,MATMUL(TF,E1)),DF)
               GRAD(2,:)=MATMUL(MATMUL(INV,MATMUL(TF,E2)),DF)
               CALL BM_Q8_AT(XYZ,RG,SG,BMG,AJ)
               NV=MATMUL(SHELL_A,MATMUL(BMG,UEB(1:48)))
               SIG(1,:)=(/NV(1),NV(3)/)
               SIG(2,:)=(/NV(3),NV(2)/)
               BLOCK=BLOCK+MATMUL(TRANSPOSE(GRAD),MATMUL(SIG,GRAD))*AREA*WG
         ENDDO
      ENDDO
      BLOCK=(BLOCK+TRANSPOSE(BLOCK))/TWO
      KED(1:48,1:48)=ZERO
      DO IN=1,8
         DO JN=1,8
            DO ID=1,3
               KED(6*(IN-1)+ID,6*(JN-1)+ID)=BLOCK(IN,JN)
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE NATIVE_MEMBRANE_KG

      SUBROUTINE KG_GEOMETRY_ERROR
      WRITE(ERR,*) ' *ERROR: native Q8/T6 membrane buckling requires nondegenerate flat geometry. EID=',EID
      WRITE(F06,*) ' *ERROR: native Q8/T6 membrane buckling requires nondegenerate flat geometry. EID=',EID
      FATAL_ERR=FATAL_ERR+1
      NUM_EMG_FATAL_ERRS=NUM_EMG_FATAL_ERRS+1
      CALL OUTA_HERE('Y')
      END SUBROUTINE KG_GEOMETRY_ERROR



      SUBROUTINE LOAD_BASIC_COORDS_Q8 ( XYZOUT )
      REAL(DOUBLE), INTENT(OUT) :: XYZOUT(8,3)
      INTEGER(LONG) :: II, JJ
      DO II=1,8
         DO JJ=1,3
            XYZOUT(II,JJ) = XEB(II,JJ)
         ENDDO
      ENDDO
      END SUBROUTINE LOAD_BASIC_COORDS_Q8

      SUBROUTINE SHAPE_Q8 ( XI, ETA, NVAL, DN )
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
      END SUBROUTINE SHAPE_Q8

      SUBROUTINE CALC_NODAL_NORMALS_Q8 ( XYZN, NORMS )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3)
      REAL(DOUBLE), INTENT(OUT) :: NORMS(8,3)
      REAL(DOUBLE) :: RS(8,2), NVAL(8), DN(2,8), G1(3), G2(3), N(3), NM, SN(3)
      INTEGER(LONG) :: II, BIDX, ISN
      RS(1,:) = (/-ONE, -ONE/)
      RS(2,:) = (/ ONE, -ONE/)
      RS(3,:) = (/ ONE,  ONE/)
      RS(4,:) = (/-ONE,  ONE/)
      RS(5,:) = (/ZERO, -ONE/)
      RS(6,:) = (/ ONE, ZERO/)
      RS(7,:) = (/ZERO,  ONE/)
      RS(8,:) = (/-ONE, ZERO/)
      DO II=1,8
         CALL SHAPE_Q8(RS(II,1), RS(II,2), NVAL, DN)
         G1 = MATMUL(DN(1,:), XYZN)
         G2 = MATMUL(DN(2,:), XYZN)
         CALL CROSS3(G1, G2, N)
         NM = VNORM(N)
         IF (NM <= 1.0D-12) THEN
            CALL SHAPE_Q8(ZERO, ZERO, NVAL, DN)
            G1 = MATMUL(DN(1,:), XYZN)
            G2 = MATMUL(DN(2,:), XYZN)
            CALL CROSS3(G1, G2, N)
            NM = VNORM(N)
         ENDIF
         IF (NM > 1.0D-15) THEN
            NORMS(II,:) = NORMAL_SIGN*N/NM
         ELSE
            NORMS(II,:) = (/ZERO, ZERO, ONE/)
         ENDIF
         IF (NSNORM > 0 .AND. ALLOCATED(GRID_SNORM) .AND. ALLOCATED(SNORM)) THEN
            DO ISN=1,NSNORM
            IF (SNORM(ISN,1) /= GRID_ID(BGRID(II))) CYCLE
            BIDX = 0
            IF (II <= SIZE(BGRID)) BIDX = BGRID(II)
            IF ((BIDX > 0) .AND. (BIDX <= SIZE(GRID_SNORM,1))) THEN
               SN = GRID_SNORM(BIDX,:)
               NM = VNORM(SN)
               IF (NM > 1.0D-15) THEN
                  SN = SN/NM
                  IF (DOT_PRODUCT(SN, NORMS(II,:)) < ZERO) SN = -SN
                  NORMS(II,:) = SN
               ENDIF
            ENDIF
            ENDDO
         ENDIF
      ENDDO
      END SUBROUTINE CALC_NODAL_NORMALS_Q8

      SUBROUTINE SURFACE_BASIS_Q8 ( XYZN, XI, ETA, G1, G2, E1, E2, E3, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: G1(3), G2(3), E1(3), E2(3), E3(3), JAC
      REAL(DOUBLE) :: NVAL(8), DN(2,8), G3(3), TMP(3), NM
      CALL SHAPE_Q8(XI, ETA, NVAL, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      CALL CROSS3(G1, G2, G3)
      JAC = VNORM(G3)
      IF (JAC > 1.0D-15) THEN
         E3 = NORMAL_SIGN*G3/JAC
      ELSE
         E3 = (/ZERO, ZERO, ONE/)
      ENDIF
      CALL CROSS3(G2, E3, TMP)
      NM = VNORM(TMP)
      IF (NM < 1.0D-12) THEN
         CALL CROSS3(G1, E3, TMP)
         NM = VNORM(TMP)
      ENDIF
      IF (NM > 1.0D-15) THEN
         E1 = TMP/NM
      ELSE
         E1 = (/ONE, ZERO, ZERO/)
      ENDIF
      CALL CROSS3(E3, E1, E2)
      NM = VNORM(E2)
      IF (NM > 1.0D-15) E2 = E2/NM
      END SUBROUTINE SURFACE_BASIS_Q8

      SUBROUTINE FIXED_FRAME_Q8 ( XYZN, E1F, E2F, E3F )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3)
      REAL(DOUBLE), INTENT(OUT) :: E1F(3), E2F(3), E3F(3)
      REAL(DOUBLE) :: G1(3), G2(3), JAC, ZSGN
      CALL SURFACE_BASIS_Q8(XYZN, ZERO, ZERO, G1, G2, E1F, E2F, E3F, JAC)
      IF ((DABS(E3F(1)) + DABS(E3F(2))) <= 1.0D-12) THEN
         ZSGN = ONE
         IF (E3F(3) < ZERO) ZSGN = -ONE
         E1F = (/ONE, ZERO, ZERO/)
         E2F = (/ZERO, ZSGN, ZERO/)
         E3F = (/ZERO, ZERO, ZSGN/)
      ENDIF
      END SUBROUTINE FIXED_FRAME_Q8

      SUBROUTINE LOCAL_BASIS_AT_Q8 ( XYZN, XI, ETA, E1, E2, E3, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: E1(3), E2(3), E3(3), JAC
      REAL(DOUBLE) :: NVAL(8), DN(2,8), G1(3), G2(3), G3(3), NM
      CALL SHAPE_Q8(XI, ETA, NVAL, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      CALL CROSS3(G1, G2, G3)
      JAC = VNORM(G3)
      IF (JAC > 1.0D-15) THEN
         E3 = NORMAL_SIGN*G3/JAC
      ELSE
         E3 = (/ZERO, ZERO, ONE/)
      ENDIF
      NM = VNORM(G1)
      IF (NM > 1.0D-15) THEN
         E1 = G1/NM
      ELSE
         E1 = (/ONE, ZERO, ZERO/)
      ENDIF
      CALL CROSS3(E3, E1, E2)
      NM = VNORM(E2)
      IF (NM > 1.0D-15) THEN
         E2 = E2/NM
      ELSE
         E2 = (/ZERO, ONE, ZERO/)
      ENDIF
      END SUBROUTINE LOCAL_BASIS_AT_Q8

      SUBROUTINE COV_MAP_Q8 ( XYZN, XI, ETA, G1, G2, C, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: G1(3), G2(3), C(4), JAC
      REAL(DOUBLE) :: E1(3), E2(3), E3(3), NVAL(8), DN(2,8)
      REAL(DOUBLE) :: A11, A22, A12, DET, AI11, AI22, AI12, GC1(3), GC2(3)
      CALL SHAPE_Q8(XI, ETA, NVAL, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      CALL LOCAL_BASIS_AT_Q8(XYZN, XI, ETA, E1, E2, E3, JAC)
      A11 = DOT_PRODUCT(G1,G1)
      A22 = DOT_PRODUCT(G2,G2)
      A12 = DOT_PRODUCT(G1,G2)
      DET = A11*A22 - A12*A12
      IF (DABS(DET) <= 1.0D-30) THEN
         C = ZERO
         RETURN
      ENDIF
      AI11 =  A22/DET
      AI22 =  A11/DET
      AI12 = -A12/DET
      GC1 = AI11*G1 + AI12*G2
      GC2 = AI12*G1 + AI22*G2
      ! Python _covariant_maps order:
      !   C = (c1_1, c1_2, c2_1, c2_2)
      ! where c1_* maps the first physical in-plane basis direction and
      ! c2_* maps the second physical in-plane basis direction.
      C(1) = DOT_PRODUCT(E1, GC1)
      C(2) = DOT_PRODUCT(E2, GC1)
      C(3) = DOT_PRODUCT(E1, GC2)
      C(4) = DOT_PRODUCT(E2, GC2)
      END SUBROUTINE COV_MAP_Q8

      SUBROUTINE COV_MAP_CENTER_Q8 ( XYZN, XI, ETA, E1F, E2F, G1, G2, C, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), XI, ETA, E1F(3), E2F(3)
      REAL(DOUBLE), INTENT(OUT) :: G1(3), G2(3), C(4), JAC
      REAL(DOUBLE) :: E1(3), E2(3), E3(3)
      REAL(DOUBLE) :: A11, A22, A12, DET, AI11, AI22, AI12, GC1(3), GC2(3)
      CALL SURFACE_BASIS_Q8(XYZN, XI, ETA, G1, G2, E1, E2, E3, JAC)
      A11 = DOT_PRODUCT(G1,G1)
      A22 = DOT_PRODUCT(G2,G2)
      A12 = DOT_PRODUCT(G1,G2)
      DET = A11*A22 - A12*A12
      IF (DABS(DET) <= 1.0D-30) THEN
         C = ZERO
         RETURN
      ENDIF
      AI11 =  A22/DET
      AI22 =  A11/DET
      AI12 = -A12/DET
      GC1 = AI11*G1 + AI12*G2
      GC2 = AI12*G1 + AI22*G2
      ! Same C order as COV_MAP_Q8 / Python _covariant_maps:
      !   C = (c1_1, c1_2, c2_1, c2_2)
      C(1) = DOT_PRODUCT(E1F, GC1)
      C(2) = DOT_PRODUCT(E2F, GC1)
      C(3) = DOT_PRODUCT(E1F, GC2)
      C(4) = DOT_PRODUCT(E2F, GC2)
      END SUBROUTINE COV_MAP_CENTER_Q8

      SUBROUTINE TENSOR_PHYS_Q8 ( V11, V22, V12, C, VOUT )
      REAL(DOUBLE), INTENT(IN)  :: V11(3), V22(3), V12(3), C(4)
      REAL(DOUBLE), INTENT(OUT) :: VOUT(3,3)
      VOUT(1,:) = C(1)*C(1)*V11 + C(3)*C(3)*V22 + TWO*C(1)*C(3)*V12
      VOUT(2,:) = C(2)*C(2)*V11 + C(4)*C(4)*V22 + TWO*C(2)*C(4)*V12
      VOUT(3,:) = TWO*(C(1)*C(2)*V11 + C(3)*C(4)*V22 + (C(1)*C(4)+C(2)*C(3))*V12)
      END SUBROUTINE TENSOR_PHYS_Q8

      SUBROUTINE BM_Q8_AT ( XYZN, XI, ETA, BMOUT, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BMOUT(3,48), JAC
      REAL(DOUBLE) :: NVAL(8), DN(2,8), G1(3), G2(3), C(4), VP(3,3)
      INTEGER(LONG) :: II, COL
      CALL SHAPE_Q8(XI, ETA, NVAL, DN)
      CALL COV_MAP_Q8(XYZN, XI, ETA, G1, G2, C, JAC)
      IF (USE_ANS) THEN
         CALL ANS_BM_Q8(XYZN,XI,ETA,C,BMOUT)
         RETURN
      ENDIF
      BMOUT = ZERO
      DO II=1,8
         COL = (II-1)*6
         CALL TENSOR_PHYS_Q8(DN(1,II)*G1, DN(2,II)*G2, 0.5D0*(DN(1,II)*G2 + DN(2,II)*G1), C, VP)
         BMOUT(1:3,COL+1:COL+3) = VP
      ENDDO
      END SUBROUTINE BM_Q8_AT

      SUBROUTINE BB_Q8_AT ( XYZN, NORMS, XI, ETA, BBOUT, JAC, NORMS_EXT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), NORMS(8,3), XI, ETA
      REAL(DOUBLE), INTENT(IN), OPTIONAL :: NORMS_EXT(8,3)
      REAL(DOUBLE), INTENT(OUT) :: BBOUT(3,48), JAC
      REAL(DOUBLE) :: NVAL(8), DN(2,8), G1(3), G2(3), C(4), VP(3,3)
      REAL(DOUBLE) :: T1(3), T2(3), T0(3), CG1(3), CG2(3)
      REAL(DOUBLE) :: NORMS_LOC(8,3)
      INTEGER(LONG) :: II, COL
      NORMS_LOC = NORMS
      IF (PRESENT(NORMS_EXT)) NORMS_LOC = NORMS_EXT
      CALL SHAPE_Q8(XI, ETA, NVAL, DN)
      CALL COV_MAP_Q8(XYZN, XI, ETA, G1, G2, C, JAC)
      T1 = MATMUL(DN(1,:), NORMS_LOC)
      T2 = MATMUL(DN(2,:), NORMS_LOC)
      BBOUT = ZERO
      DO II=1,8
         COL = (II-1)*6
         T0 = NORMS_LOC(II,:)
         CALL TENSOR_PHYS_Q8(DN(1,II)*T1, DN(2,II)*T2, 0.5D0*(DN(1,II)*T2 + DN(2,II)*T1), C, VP)
         BBOUT(1:3,COL+1:COL+3) = -VP
         CALL CROSS3(G1, T0, CG1)
         CALL CROSS3(G2, T0, CG2)
         CALL TENSOR_PHYS_Q8(DN(1,II)*CG1, DN(2,II)*CG2, 0.5D0*(DN(1,II)*CG2 + DN(2,II)*CG1), C, VP)
         BBOUT(1:3,COL+4:COL+6) = VP
      ENDDO
      END SUBROUTINE BB_Q8_AT

      SUBROUTINE BS_Q8_AT ( XYZN, NORMS, XI, ETA, BSOUT, JAC, NORMS_EXT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), NORMS(8,3), XI, ETA
      REAL(DOUBLE), INTENT(IN), OPTIONAL :: NORMS_EXT(8,3)
      REAL(DOUBLE), INTENT(OUT) :: BSOUT(2,48), JAC
      REAL(DOUBLE) :: NVAL(8), DN(2,8), G1(3), G2(3), C(4), BSN(2,48), T0(3), T0I(3), C1(3), C2(3), NM
      REAL(DOUBLE) :: NORMS_LOC(8,3)
      INTEGER(LONG) :: II, COL
      NORMS_LOC = NORMS
      IF (PRESENT(NORMS_EXT)) NORMS_LOC = NORMS_EXT
      CALL SHAPE_Q8(XI, ETA, NVAL, DN)
      CALL COV_MAP_Q8(XYZN, XI, ETA, G1, G2, C, JAC)
      T0 = MATMUL(NVAL, NORMS_LOC)
      IF (USE_ANS) THEN
         CALL ANS_BS_Q8(XYZN,NORMS_LOC,XI,ETA,C,BSOUT)
         RETURN
      ENDIF
      BSN = ZERO
      DO II=1,8
         COL = (II-1)*6
         T0I = NORMS_LOC(II,:)
         BSN(1,COL+1:COL+3) = BSN(1,COL+1:COL+3) + DN(1,II)*T0
         BSN(2,COL+1:COL+3) = BSN(2,COL+1:COL+3) + DN(2,II)*T0
         CALL CROSS3(T0I, G1, C1)
         CALL CROSS3(T0I, G2, C2)
         BSN(1,COL+4:COL+6) = BSN(1,COL+4:COL+6) + NVAL(II)*C1
         BSN(2,COL+4:COL+6) = BSN(2,COL+4:COL+6) + NVAL(II)*C2
      ENDDO
      BSOUT(1,:) = C(1)*BSN(1,:) + C(3)*BSN(2,:)
      BSOUT(2,:) = C(2)*BSN(1,:) + C(4)*BSN(2,:)
      END SUBROUTINE BS_Q8_AT

      SUBROUTINE BDRILL_Q8_AT ( XYZN, NORMS, XI, ETA, BDOUT, JAC, NORMS_EXT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), NORMS(8,3), XI, ETA
      REAL(DOUBLE), INTENT(IN), OPTIONAL :: NORMS_EXT(8,3)
      REAL(DOUBLE), INTENT(OUT) :: BDOUT(1,48), JAC
      REAL(DOUBLE) :: NVAL(8), DN(2,8), G1(3), G2(3), E1(3), E2(3), E3(3)
      REAL(DOUBLE) :: A(2,2), AINV(2,2), DLOC(2,8), DX(8), DY(8)
      REAL(DOUBLE) :: NORMS_LOC(8,3)
      INTEGER(LONG) :: II, COL
      NORMS_LOC = NORMS
      IF (PRESENT(NORMS_EXT)) NORMS_LOC = NORMS_EXT
      CALL SHAPE_Q8(XI, ETA, NVAL, DN)
      CALL LOCAL_BASIS_AT_Q8(XYZN, XI, ETA, E1, E2, E3, JAC)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      A(1,1) = DOT_PRODUCT(G1,G1)
      A(1,2) = DOT_PRODUCT(G1,G2)
      A(2,1) = DOT_PRODUCT(G1,G2)
      A(2,2) = DOT_PRODUCT(G2,G2)
      CALL INV2(A, AINV)
      DLOC = MATMUL(AINV, DN)
      DX=DLOC(1,:)*DOT_PRODUCT(G1,E1)+DLOC(2,:)*DOT_PRODUCT(G2,E1)
      DY=DLOC(1,:)*DOT_PRODUCT(G1,E2)+DLOC(2,:)*DOT_PRODUCT(G2,E2)
      BDOUT = ZERO
      DO II=1,8
         COL = (II-1)*6
         BDOUT(1,COL+1:COL+3) = 0.5D0*(DX(II)*E2 - DY(II)*E1)
         BDOUT(1,COL+4:COL+6) = BDOUT(1,COL+4:COL+6) - NVAL(II)*E3
      ENDDO
      END SUBROUTINE BDRILL_Q8_AT

      SUBROUTINE EAS1_SHEAR_Q8_AT ( XYZN, XI, ETA, BSE )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BSE(2)
      REAL(DOUBLE) :: E1F(3), E2F(3), E3F(3), G1(3), G2(3), C(4), JAC, DPHIR, DPHIS
      CALL COV_MAP_Q8(XYZN,XI,ETA,G1,G2,C,JAC)
      DPHIR = -TWO*XI*(ONE-ETA*ETA)
      DPHIS = -TWO*ETA*(ONE-XI*XI)
      BSE(1) = DPHIR*C(1) + DPHIS*C(3)
      BSE(2) = DPHIR*C(2) + DPHIS*C(4)
      END SUBROUTINE EAS1_SHEAR_Q8_AT

      SUBROUTINE SET_NORMAL_SIGN_Q8(XYZN)
      REAL(DOUBLE), INTENT(IN) :: XYZN(8,3)
      REAL(DOUBLE) :: N(8),DN(2,8),NC(3),G1(3),G2(3)
      CALL SHAPE_Q8(ZERO,ZERO,N,DN)
      G1=MATMUL(DN(1,:),XYZN)
      G2=MATMUL(DN(2,:),XYZN)
      CALL CROSS3(G1,G2,NC)
      NORMAL_SIGN=ONE
      IF (NC(3) < -1.0D-6*VNORM(NC)) NORMAL_SIGN=-ONE
      END SUBROUTINE

      SUBROUTINE SELECT_ANS_Q8(XYZN,NORMS)
      REAL(DOUBLE), INTENT(IN) :: XYZN(8,3),NORMS(8,3)
      REAL(DOUBLE) :: E1(3),E2(3),E3(3),JAC,ANG,DEV,LENGTH
      INTEGER(LONG) :: II,JJ
      CALL LOCAL_BASIS_AT_Q8(XYZN,ZERO,ZERO,E1,E2,E3,JAC)
      ANG=ZERO
      DEV=ZERO
      DO II=1,8
         ANG=MAX(ANG,DACOS(MAX(-ONE,MIN(ONE,DOT_PRODUCT(NORMS(II,:),E3)))))
      ENDDO
      DO II=1,4
         JJ=MOD(II,4)+1
         LENGTH=VNORM(XYZN(JJ,:)-XYZN(II,:))
         DEV=MAX(DEV,VNORM(XYZN(II+4,:)-0.5D0*(XYZN(II,:)+XYZN(JJ,:)))/MAX(LENGTH,1.0D-30))
      ENDDO
      USE_ANS=(ANG > 0.01D0 .OR. DEV > 0.001D0)
      END SUBROUTINE

      SUBROUTINE ANS_ROWS_Q8(XYZN,NORMS,R,S,EM,ES)
      REAL(DOUBLE), INTENT(IN) :: XYZN(8,3),NORMS(8,3),R,S
      REAL(DOUBLE), INTENT(OUT) :: EM(3,48),ES(2,48)
      REAL(DOUBLE) :: N(8),DN(2,8),G1(3),G2(3),T0(3),CV(3)
      INTEGER(LONG) :: II,COL
      CALL SHAPE_Q8(R,S,N,DN)
      G1=MATMUL(DN(1,:),XYZN)
      G2=MATMUL(DN(2,:),XYZN)
      T0=MATMUL(N,NORMS)
      EM=ZERO
      ES=ZERO
      DO II=1,8
         COL=6*(II-1)
         EM(1,COL+1:COL+3)=DN(1,II)*G1
         EM(2,COL+1:COL+3)=DN(2,II)*G2
         EM(3,COL+1:COL+3)=0.5D0*(DN(1,II)*G2+DN(2,II)*G1)
         ES(1,COL+1:COL+3)=DN(1,II)*T0
         ES(2,COL+1:COL+3)=DN(2,II)*T0
         CALL CROSS3(NORMS(II,:),G1,CV)
         ES(1,COL+4:COL+6)=N(II)*CV
         CALL CROSS3(NORMS(II,:),G2,CV)
         ES(2,COL+4:COL+6)=N(II)*CV
      ENDDO
      END SUBROUTINE

      SUBROUTINE ANS_INTERP_Q8(XYZN,NORMS,R,S,COMP,ISMEM,ROW)
      REAL(DOUBLE), INTENT(IN) :: XYZN(8,3),NORMS(8,3),R,S
      INTEGER(LONG), INTENT(IN) :: COMP
      LOGICAL, INTENT(IN) :: ISMEM
      REAL(DOUBLE), INTENT(OUT) :: ROW(48)
      REAL(DOUBLE) :: A,B,PA(2),PB(3),LR(2),LS(2),QR(3),QS(3),EM(3,48),ES(2,48),P,Q,W
      INTEGER(LONG) :: II,JJ,NJ
      A=ONE/DSQRT(3.0D0)
      B=DSQRT(0.6D0)
      PA=(/-A,A/)
      PB=(/-B,ZERO,B/)
      LR=(/0.5D0*(ONE-R/A),0.5D0*(ONE+R/A)/)
      LS=(/0.5D0*(ONE-S/A),0.5D0*(ONE+S/A)/)
      QR=(/R*(R-B)/(TWO*B*B),ONE-R*R/(B*B),R*(R+B)/(TWO*B*B)/)
      QS=(/S*(S-B)/(TWO*B*B),ONE-S*S/(B*B),S*(S+B)/(TWO*B*B)/)
      ROW=ZERO
      NJ=3
      IF (COMP == 3) NJ=2
      DO II=1,2
         DO JJ=1,NJ
            IF (COMP == 1) THEN
               P=PA(II); Q=PB(JJ); W=LR(II)*QS(JJ)
            ELSE IF (COMP == 2) THEN
               P=PB(JJ); Q=PA(II); W=LS(II)*QR(JJ)
            ELSE
               P=PA(II); Q=PA(JJ); W=LR(II)*LS(JJ)
            ENDIF
            CALL ANS_ROWS_Q8(XYZN,NORMS,P,Q,EM,ES)
            IF (ISMEM) THEN
               ROW=ROW+W*EM(COMP,:)
            ELSE
               ROW=ROW+W*ES(COMP,:)
            ENDIF
         ENDDO
      ENDDO
      END SUBROUTINE

      SUBROUTINE ANS_BM_Q8(XYZN,R,S,C,BMOUT)
      REAL(DOUBLE), INTENT(IN) :: XYZN(8,3),R,S,C(4)
      REAL(DOUBLE), INTENT(OUT) :: BMOUT(3,48)
      REAL(DOUBLE) :: ROWS(3,48),NORMS(8,3)
      INTEGER(LONG) :: II
      NORMS=ZERO
      DO II=1,3
         CALL ANS_INTERP_Q8(XYZN,NORMS,R,S,II,.TRUE.,ROWS(II,:))
      ENDDO
      BMOUT(1,:)=C(1)**2*ROWS(1,:)+C(3)**2*ROWS(2,:)+TWO*C(1)*C(3)*ROWS(3,:)
      BMOUT(2,:)=C(2)**2*ROWS(1,:)+C(4)**2*ROWS(2,:)+TWO*C(2)*C(4)*ROWS(3,:)
      BMOUT(3,:)=TWO*(C(1)*C(2)*ROWS(1,:)+C(3)*C(4)*ROWS(2,:)+(C(1)*C(4)+C(3)*C(2))*ROWS(3,:))
      END SUBROUTINE

      SUBROUTINE ANS_BS_Q8(XYZN,NORMS,R,S,C,BSOUT)
      REAL(DOUBLE), INTENT(IN) :: XYZN(8,3),NORMS(8,3),R,S,C(4)
      REAL(DOUBLE), INTENT(OUT) :: BSOUT(2,48)
      REAL(DOUBLE) :: ROWS(2,48)
      INTEGER(LONG) :: II
      DO II=1,2
         CALL ANS_INTERP_Q8(XYZN,NORMS,R,S,II,.FALSE.,ROWS(II,:))
      ENDDO
      BSOUT(1,:)=C(1)*ROWS(1,:)+C(3)*ROWS(2,:)
      BSOUT(2,:)=C(2)*ROWS(1,:)+C(4)*ROWS(2,:)
      END SUBROUTINE

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

      SUBROUTINE INV2 ( A, AINV )
      REAL(DOUBLE), INTENT(IN)  :: A(2,2)
      REAL(DOUBLE), INTENT(OUT) :: AINV(2,2)
      REAL(DOUBLE) :: DET
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

      END SUBROUTINE CQUAD8_SIMOQ8
