! #################################################################################################################################
! CTRIA6 Simo1993 quadratic triangular shell.

      SUBROUTINE CTRIA6_SIMO1993 ( OPT, INT_ELEM_ID )

! Ported from:
!   D:\18a\bending_only\Shell\gemini2\shit\validation\q8\Simo1993_Tri6_ShellElement_v1p8.py
!
! Static stiffness path:
!   6-node quadratic triangle, Simo/Fox director kinematics,
!   Hughes-Brezzi drilling penalty, no EAS condensation.
!   Membrane/bending/shear use the 3-point degree-2 triangle rule.
!   Drilling uses the 6-point degree-4 triangle rule from the Python source
!   to avoid rank-deficient theta-z drilling modes.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_STRESS_POINTS, SOL_NAME, NSNORM
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  ALPVEC, BGRID, DT, EID, ELGP, EPROP, GRID_SNORM, GRID_ID, SNORM, KE, KED, ME, BE1, BE2, BE3, &
                                         MASS_PER_UNIT_AREA, NUM_EMG_FATAL_ERRS, PCOMP_PROPS, PPE, PRESS, PTE, RGRID, SHELL_A,   &
                                         SHELL_D, SHELL_T, TREF, UEL, XEB
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO
      USE PARAMS, ONLY                :  COUPMASS
      USE ELMDIS_Interface
      USE OUTA_HERE_Interface

      USE QUADRATIC_SURFACE_MASS_Interface
      USE QUADRATIC_SURFACE_PRESSURE_Interface
      USE MODEL_STUF, ONLY : UEB, KED
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CTRIA6_SIMO1993'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      INTEGER(LONG)                   :: I, J, IA, IB, K, L, JSUB
      REAL(DOUBLE)                    :: XYZ(6,3), NORMALS(6,3)
      REAL(DOUBLE)                    :: KOUT(36,36)
      REAL(DOUBLE)                    :: BM(3,36), BB(3,36), BS(2,36), BD(1,36)
      REAL(DOUBLE)                    :: R3(3), S3(3), W3(3), R6(6), S6(6), W6(6), RN7(7), SN7(7)
      REAL(DOUBLE)                    :: R, S, WT, JAC, CDRILL, FAC
      REAL(DOUBLE)                    :: M1(6,6), N6(6), DN6(2,6), MASS_ELEM, MASS_NODE
      REAL(DOUBLE)                    :: UNIT_PPE(36), UNIT_PTE(36), DXDR(3), DXDS(3), SURF_VEC(3), TBAR
      REAL(DOUBLE)                    :: CTE(3), THERMAL_RESULTANT(3)
      REAL(DOUBLE)                    :: BDNORM
      REAL(DOUBLE)                    :: DLOC(2,6), SIG0(2,2), KG6(6,6), STRAIN0(3), N0(3)
      REAL(DOUBLE)                    :: G1K(3), G2K(3), A11, A22, A12, DET, AI11, AI22, AI12
      INTEGER(LONG)                   :: KI, KJ

      REAL(DOUBLE) :: UNIT_PTG(36), TEMP_GRAD_SIGN, TEMP_NORMAL(3), TEMP_G1(3), TEMP_G2(3)

      IF (ELGP /= 6) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9001) SUBR_NAME, EID, ELGP
         WRITE(F06,9001) SUBR_NAME, EID, ELGP
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IF (PCOMP_PROPS == 'Y') THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,*) ' *ERROR: Code not written for composite material with CTRIA6 Simo1993'
         WRITE(F06,*) ' *ERROR: Code not written for composite material with CTRIA6 Simo1993'
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL LOAD_BASIC_COORDS_T6 ( XYZ )
      CALL CALC_NODAL_NORMALS_T6 ( XYZ, NORMALS )

      R3 = (/ONE/6.0D0, TWO/3.0D0, ONE/6.0D0/)
      S3 = (/ONE/6.0D0, ONE/6.0D0, TWO/3.0D0/)
      W3 = (/ONE/6.0D0, ONE/6.0D0, ONE/6.0D0/)

      ! Nodal coordinates for T6: center, 3 corners, 2 midside

      R6 = (/0.445948490144588D0, 0.108103018168070D0, 0.445948490144588D0,                                      &
             0.091576213509771D0, 0.816847572980459D0, 0.091576213509771D0/)
      S6 = (/0.445948490144588D0, 0.445948490144588D0, 0.108103018168070D0,                                      &
             0.091576213509771D0, 0.091576213509771D0, 0.816847572980459D0/)
      W6 = (/0.111690794839050D0, 0.111690794839050D0, 0.111690794839050D0,                                      &
             0.054975871827661D0, 0.054975871827661D0, 0.054975871827661D0/)

      RN7 = (/ONE/3.0D0, 0.0D0, 1.0D0, 0.0D0, 0.5D0, 0.5D0, 0.0D0/)
      SN7 = (/ONE/3.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.5D0, 0.5D0/)

      IF (OPT(1) == 'Y') THEN
         CALL QUADRATIC_SURFACE_MASS(6,XYZ,(/(ZERO,J=1,8)/))
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
            R=R3(I)
            S=S3(I)
            WT=W3(I)
               CALL BM_T6_AT ( XYZ, R, S, BM, JAC )
               CALL BB_T6_AT ( XYZ, NORMALS, R, S, BB, JAC )
               UNIT_PTE=UNIT_PTE+MATMUL(TRANSPOSE(BM),MATMUL(SHELL_A,CTE))*WT*JAC
! eps(z)=Bm*u-z*Bb*u: bending thermal force has the physical minus sign.
               UNIT_PTG=UNIT_PTG-MATMUL(TRANSPOSE(BB),MATMUL(SHELL_D,CTE))*WT*JAC
         ENDDO
         DO JSUB=1,SIZE(PTE,2)
            TBAR=SUM(DT(1:ELGP,JSUB))/REAL(ELGP,DOUBLE)-TREF(1)
            PTE(1:36,JSUB)=UNIT_PTE*TBAR+UNIT_PTG*TEMP_GRAD_SIGN*DT(ELGP+1,JSUB)
         ENDDO
      ENDIF

      IF (OPT(3) == 'Y') THEN
         DO I=1,7
            CALL BM_T6_AT ( XYZ, RN7(I), SN7(I), BM, JAC )
            CALL BB_T6_AT ( XYZ, NORMALS, RN7(I), SN7(I), BB, JAC )
            CALL BS_T6_AT ( XYZ, NORMALS, RN7(I), SN7(I), BS, JAC )
            IF (I <= MAX_STRESS_POINTS+1) THEN
               BE1(1:3,1:36,I) = BM
! MYSTRAN forms fiber stress as membrane - z*BE2*u and uses a negative
! bending force conversion: physical moment is minus the Bb-conjugate moment.
               BE2(1:3,1:36,I) =BB
               BE3(1:2,1:36,I) = BS
            ENDIF
         ENDDO
      ENDIF


      IF (OPT(4) == 'Y') THEN
         KOUT = ZERO
         CDRILL = 2.0D-2*SHELL_T(1,1)/(5.0D0/6.0D0)

         DO I=1,3
            R = R3(I)
            S = S3(I)
            WT = W3(I)
            CALL BM_T6_AT ( XYZ, R, S, BM, JAC )
            CALL BB_T6_AT ( XYZ, NORMALS, R, S, BB, JAC )
            CALL BS_T6_AT ( XYZ, NORMALS, R, S, BS, JAC )
            IF (EID == 1) THEN
               CALL BDRILL_T6_AT ( XYZ, NORMALS, R, S, BD, BDNORM )
            ENDIF
            KOUT = KOUT + WT*JAC*MATMUL(TRANSPOSE(BM), MATMUL(SHELL_A, BM))
            KOUT = KOUT + WT*JAC*MATMUL(TRANSPOSE(BB), MATMUL(SHELL_D, BB))
            KOUT = KOUT + WT*JAC*MATMUL(TRANSPOSE(BS), MATMUL(SHELL_T, BS))
         ENDDO

         DO I=1,6
            R = R6(I)
            S = S6(I)
            WT = W6(I)
            CALL BDRILL_T6_AT ( XYZ, NORMALS, R, S, BD, JAC )
            KOUT = KOUT + WT*JAC*CDRILL*MATMUL(TRANSPOSE(BD), BD)
         ENDDO

         DO IA=1,36
            DO IB=1,36
               FAC = 0.5D0*(KOUT(IA,IB) + KOUT(IB,IA))
               KE(IA,IB) = FAC
            ENDDO
         ENDDO
      ENDIF

      IF (OPT(5) == 'Y') THEN
         CALL QUADRATIC_SURFACE_PRESSURE(INT_ELEM_ID,6,XYZ,(/(ZERO,J=1,8)/))
      ENDIF

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN
         CALL NATIVE_MEMBRANE_KG
      ENDIF

      RETURN

 9001 FORMAT(' *ERROR: ',A,' expects ELGP=6 for element ',I8,' but got ',I8)

      CONTAINS

! Flat-shell membrane initial-stress stiffness, tension-positive resultants.
! Four-point Gauss/Duffy; standard geometry area, active translation field.
! Mechanical linear reference state only; no director or follower tangent.
      SUBROUTINE NATIVE_MEMBRANE_KG
      REAL(DOUBLE) :: GX(4),GW(4),RG,SG,WG,NVAL(6),DG(2,6),DF(2,6)
      REAL(DOUBLE) :: TG(2,3),TF(2,3),CV(3),AREA,AJ,MT(2,2),INV(2,2),DETMT
      REAL(DOUBLE) :: E1(3),E2(3),E3(3),GRAD(2,6),BMG(3,36),NV(3),SIG(2,2),BLOCK(6,6)
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
      SCALE_GEOM=MAXVAL(ABS(XYZ-SPREAD(XYZ(1,:),1,6)))
      IF (AREA <= 1.0D-14) THEN
         CALL KG_GEOMETRY_ERROR
      ENDIF
      NORMAL=NORMAL/AREA
      DO IN=1,6
         IF (ABS(DOT_PRODUCT(XYZ(IN,:)-XYZ(1,:),NORMAL)) > 1.0D-9*MAX(SCALE_GEOM,ONE)) THEN
            CALL KG_GEOMETRY_ERROR
         ENDIF
      ENDDO
      CALL ELMDIS
      BLOCK=ZERO
      DO IG=1,4
         DO JG=1,4
               RG=(GX(IG)+ONE)/TWO
               SG=(ONE-RG)*(GX(JG)+ONE)/TWO
               WG=GW(IG)*GW(JG)*(ONE-RG)/4.0D0
               CALL SHAPE_T6(RG,SG,NVAL,DG)
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
               CALL SURFACE_BASIS_T6(XYZ,RG,SG,TG(1,:),TG(2,:),E1,E2,E3,AJ)
               GRAD(1,:)=MATMUL(MATMUL(INV,MATMUL(TF,E1)),DF)
               GRAD(2,:)=MATMUL(MATMUL(INV,MATMUL(TF,E2)),DF)
               CALL BM_T6_AT(XYZ,RG,SG,BMG,AJ)
               NV=MATMUL(SHELL_A,MATMUL(BMG,UEB(1:36)))
               SIG(1,:)=(/NV(1),NV(3)/)
               SIG(2,:)=(/NV(3),NV(2)/)
               BLOCK=BLOCK+MATMUL(TRANSPOSE(GRAD),MATMUL(SIG,GRAD))*AREA*WG
         ENDDO
      ENDDO
      BLOCK=(BLOCK+TRANSPOSE(BLOCK))/TWO
      KED(1:36,1:36)=ZERO
      DO IN=1,6
         DO JN=1,6
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



      SUBROUTINE LOAD_BASIC_COORDS_T6 ( XYZOUT )
      REAL(DOUBLE), INTENT(OUT) :: XYZOUT(6,3)
      INTEGER(LONG) :: II, JJ
      DO II=1,6
         IF ((II <= SIZE(BGRID)) .AND. (BGRID(II) > 0)) THEN
            DO JJ=1,3
               XYZOUT(II,JJ) = RGRID(BGRID(II),JJ)
            ENDDO
         ELSE
            DO JJ=1,3
               XYZOUT(II,JJ) = XEB(II,JJ)
            ENDDO
         ENDIF
      ENDDO
      END SUBROUTINE LOAD_BASIC_COORDS_T6

      SUBROUTINE SHAPE_T6 ( R, S, NVAL, DN )
      REAL(DOUBLE), INTENT(IN)  :: R, S
      REAL(DOUBLE), INTENT(OUT) :: NVAL(6), DN(2,6)
      REAL(DOUBLE) :: L1, L2, L3
      L1 = ONE - R - S
      L2 = R
      L3 = S
      NVAL(1) = L1*(TWO*L1 - ONE)
      NVAL(2) = L2*(TWO*L2 - ONE)
      NVAL(3) = L3*(TWO*L3 - ONE)
      NVAL(4) = 4.0D0*L1*L2
      NVAL(5) = 4.0D0*L2*L3
      NVAL(6) = 4.0D0*L3*L1

      DN(1,1) = 4.0D0*R + 4.0D0*S - 3.0D0
      DN(1,2) = 4.0D0*R - ONE
      DN(1,3) = ZERO
      DN(1,4) = 4.0D0 - 8.0D0*R - 4.0D0*S
      DN(1,5) = 4.0D0*S
      DN(1,6) = -4.0D0*S

      DN(2,1) = 4.0D0*R + 4.0D0*S - 3.0D0
      DN(2,2) = ZERO
      DN(2,3) = 4.0D0*S - ONE
      DN(2,4) = -4.0D0*R
      DN(2,5) = 4.0D0*R
      DN(2,6) = 4.0D0 - 4.0D0*R - 8.0D0*S
      END SUBROUTINE SHAPE_T6

      SUBROUTINE CALC_NODAL_NORMALS_T6 ( XYZN, NORMS )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3)
      REAL(DOUBLE), INTENT(OUT) :: NORMS(6,3)
      REAL(DOUBLE) :: RS(6,2), NVAL(6), DN(2,6), G1(3), G2(3), N(3), NC(3), NM, SN(3), SDOT, NORMAL_SIGN
      INTEGER(LONG) :: II, BIDX, ISN
      RS(1,:) = (/ZERO, ZERO/)
      RS(2,:) = (/ONE, ZERO/)
      RS(3,:) = (/ZERO, ONE/)
      RS(4,:) = (/0.5D0, ZERO/)
      RS(5,:) = (/0.5D0, 0.5D0/)
      RS(6,:) = (/ZERO, 0.5D0/)
!     v2 (Simo1993_Tri6_ShellElement_v2.py): coherent director sign.
!     Compute the centroid normal once; apply its Z sign to every nodal
!     normal so the directors agree across the element. This replaces the
!     v1p8 "tie-break +Z per node" which forced every warped midside
!     director to +Z and broke downstream cross(g, t) vs cross(t, g) signs.
      CALL SHAPE_T6(ONE/3.0D0, ONE/3.0D0, NVAL, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      CALL CROSS3(G1, G2, NC)
      NM = VNORM(NC)
      NORMAL_SIGN = -ONE
      IF (NM > 1.0D-15) THEN
         IF (NC(3) >= -1.0D-06*NM) NORMAL_SIGN = ONE
      ENDIF
      DO II=1,6
         CALL SHAPE_T6(RS(II,1), RS(II,2), NVAL, DN)
         G1 = MATMUL(DN(1,:), XYZN)
         G2 = MATMUL(DN(2,:), XYZN)
         CALL CROSS3(G1, G2, N)
         NM = VNORM(N)
         IF (NM <= 1.0D-12) THEN
            N = NC
            NM = VNORM(NC)
         ENDIF
         IF (NM > 1.0D-15) THEN
            NORMS(II,:) = NORMAL_SIGN*N/NM
         ELSE
            NORMS(II,:) = (/ZERO, ZERO, NORMAL_SIGN/)
         ENDIF
! Explicit SNORM overrides geometric directors; automatically averaged
! GRID_SNORM values do not replace the final Python geometric defaults.
         IF (NSNORM > 0 .AND. ALLOCATED(GRID_SNORM) .AND. ALLOCATED(SNORM)) THEN
            BIDX = BGRID(II)
            IF (BIDX > 0) THEN
               DO ISN=1,NSNORM
                  IF (SNORM(ISN,1) /= GRID_ID(BIDX)) CYCLE
                  SN = GRID_SNORM(BIDX,:)
                  NM = VNORM(SN)
                  IF (NM > 1.0D-15) THEN
                     SN = SN/NM
                     IF (DOT_PRODUCT(SN,NORMS(II,:)) < ZERO) SN = -SN
                     NORMS(II,:) = SN
                  ENDIF
                  EXIT
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      END SUBROUTINE CALC_NODAL_NORMALS_T6

      SUBROUTINE SURFACE_BASIS_T6 ( XYZN, R, S, G1, G2, E1, E2, E3, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3), R, S
      REAL(DOUBLE), INTENT(OUT) :: G1(3), G2(3), E1(3), E2(3), E3(3), JAC
      REAL(DOUBLE) :: NVAL(6), DN(2,6), G3(3), TMP(3), NC(3), NM, SGN
      CALL SHAPE_T6(R, S, NVAL, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      CALL CROSS3(G1, G2, G3)
      JAC = VNORM(G3)
      IF (JAC > 1.0D-15) THEN
         E3 = G3/JAC
      ELSE
         E3 = (/ZERO, ZERO, ONE/)
      ENDIF
      CALL SHAPE_T6(ONE/3.0D0, ONE/3.0D0, NVAL, DN)
      NC = MATMUL(DN(1,:), XYZN)
      TMP = MATMUL(DN(2,:), XYZN)
      CALL CROSS3(NC, TMP, G3)
      SGN = ONE
      IF (G3(3) < -1.0D-6*VNORM(G3)) SGN = -ONE
      E3 = SGN*E3
      NM = VNORM(G1)
      IF (NM > 1.0D-15) THEN
         E1 = G1/NM
      ELSE
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
      ENDIF
      CALL CROSS3(E3, E1, E2)
      NM = VNORM(E2)
      IF (NM > 1.0D-15) E2 = E2/NM
      END SUBROUTINE SURFACE_BASIS_T6

      SUBROUTINE FIXED_FRAME_T6 ( XYZN, E1F, E2F, E3F )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3)
      REAL(DOUBLE), INTENT(OUT) :: E1F(3), E2F(3), E3F(3)
      REAL(DOUBLE) :: G1(3), G2(3), JAC
      CALL SURFACE_BASIS_T6(XYZN, ONE/3.0D0, ONE/3.0D0, G1, G2, E1F, E2F, E3F, JAC)
!     Use the local basis from SURFACE_BASIS_T6 directly (e1 = g1/|g1|, e2 = e3 x e1).
!     Do NOT override with a fixed global frame for planar elements - the local basis
!     must follow the element geometry to match the Python reference implementation.
      END SUBROUTINE FIXED_FRAME_T6

      SUBROUTINE COV_MAP_T6 ( XYZN, R, S, E1F, E2F, G1, G2, C, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3), R, S, E1F(3), E2F(3)
      REAL(DOUBLE), INTENT(OUT) :: G1(3), G2(3), C(4), JAC
      REAL(DOUBLE) :: E1(3), E2(3), E3(3), A11, A22, A12, DET, AI11, AI22, AI12, GC1(3), GC2(3)
      CALL SURFACE_BASIS_T6(XYZN, R, S, G1, G2, E1, E2, E3, JAC)
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
      C(1) = DOT_PRODUCT(E1, GC1)
      C(2) = DOT_PRODUCT(E1, GC2)
      C(3) = DOT_PRODUCT(E2, GC1)
      C(4) = DOT_PRODUCT(E2, GC2)
      END SUBROUTINE COV_MAP_T6

      SUBROUTINE TENSOR_PHYS_T6 ( V11, V22, V12, C, VOUT )
      REAL(DOUBLE), INTENT(IN)  :: V11(3), V22(3), V12(3), C(4)
      REAL(DOUBLE), INTENT(OUT) :: VOUT(3,3)
      VOUT(1,:) = C(1)*C(1)*V11 + C(2)*C(2)*V22 + TWO*C(1)*C(2)*V12
      VOUT(2,:) = C(3)*C(3)*V11 + C(4)*C(4)*V22 + TWO*C(3)*C(4)*V12
      VOUT(3,:) = TWO*(C(1)*C(3)*V11 + C(2)*C(4)*V22 + (C(1)*C(4)+C(2)*C(3))*V12)
      END SUBROUTINE TENSOR_PHYS_T6

      SUBROUTINE BM_T6_AT ( XYZN, R, S, BMOUT, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3), R, S
      REAL(DOUBLE), INTENT(OUT) :: BMOUT(3,36), JAC
      REAL(DOUBLE) :: NVAL(6), DN(2,6), E1F(3), E2F(3), E3F(3), G1(3), G2(3), C(4), VP(3,3)
      INTEGER(LONG) :: II, COL
      CALL SHAPE_T6(R, S, NVAL, DN)
      CALL FIXED_FRAME_T6(XYZN, E1F, E2F, E3F)
      CALL COV_MAP_T6(XYZN, R, S, E1F, E2F, G1, G2, C, JAC)
      BMOUT = ZERO
      DO II=1,6
         COL = (II-1)*6
         CALL TENSOR_PHYS_T6(DN(1,II)*G1, DN(2,II)*G2, 0.5D0*(DN(1,II)*G2 + DN(2,II)*G1), C, VP)
         BMOUT(1:3,COL+1:COL+3) = VP
      ENDDO
      END SUBROUTINE BM_T6_AT

      SUBROUTINE BB_T6_AT ( XYZN, NORMS, R, S, BBOUT, JAC, NORMS_EXT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3), NORMS(6,3), R, S
      REAL(DOUBLE), INTENT(IN), OPTIONAL :: NORMS_EXT(6,3)
      REAL(DOUBLE), INTENT(OUT) :: BBOUT(3,36), JAC
      REAL(DOUBLE) :: NVAL(6), DN(2,6), E1F(3), E2F(3), E3F(3), G1(3), G2(3), C(4), VP(3,3)
      REAL(DOUBLE) :: T1(3), T2(3), T0(3), CG1(3), CG2(3), DRHO11(3), DRHO22(3), DRHO12(3)
      REAL(DOUBLE) :: NORMS_LOC(6,3)
      INTEGER(LONG) :: II, COL
      NORMS_LOC = NORMS
      IF (PRESENT(NORMS_EXT)) NORMS_LOC = NORMS_EXT
      CALL SHAPE_T6(R, S, NVAL, DN)
      CALL FIXED_FRAME_T6(XYZN, E1F, E2F, E3F)
      CALL COV_MAP_T6(XYZN, R, S, E1F, E2F, G1, G2, C, JAC)
      T1 = MATMUL(DN(1,:), NORMS_LOC)
      T2 = MATMUL(DN(2,:), NORMS_LOC)
      BBOUT = ZERO
      DO II=1,6
         COL = (II-1)*6
         T0 = NORMS_LOC(II,:)
         CALL CROSS3(G1, T0, CG1)
         CALL CROSS3(G2, T0, CG2)
         ! Director bending (director curvature) — Python: cross(g, t0_I)
         CALL TENSOR_PHYS_T6(DN(1,II)*CG1, DN(2,II)*CG2, 0.5D0*(DN(1,II)*CG2 + DN(2,II)*CG1), C, VP)
         BBOUT(1:3,COL+4:COL+6) = VP
         ! Membrane coupling — v2 Simo1993_Tri6_ShellElement_v2.py uses -[dNr[k], dNs[k]]
         ! (sign-flipped translation block paired with the coherent director sign).
         ! Without this flip, the Bb translation block has the wrong sign vs. the Python ref.
         DRHO11 = -(/ DN(1,II)*T1(1), DN(1,II)*T1(2), DN(1,II)*T1(3) /)
         DRHO22 = -(/ DN(2,II)*T2(1), DN(2,II)*T2(2), DN(2,II)*T2(3) /)
         DRHO12 = -(/ 0.5D0*(DN(1,II)*T2(1) + DN(2,II)*T1(1)), &
                     0.5D0*(DN(1,II)*T2(2) + DN(2,II)*T1(2)), &
                     0.5D0*(DN(1,II)*T2(3) + DN(2,II)*T1(3)) /)
         CALL TENSOR_PHYS_T6(DRHO11, DRHO22, DRHO12, C, VP)
         BBOUT(1:3,COL+1:COL+3) = BBOUT(1:3,COL+1:COL+3) + VP
      ENDDO
      END SUBROUTINE BB_T6_AT

      SUBROUTINE BS_T6_AT ( XYZN, NORMS, R, S, BSOUT, JAC, NORMS_EXT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3), NORMS(6,3), R, S
      REAL(DOUBLE), INTENT(IN), OPTIONAL :: NORMS_EXT(6,3)
      REAL(DOUBLE), INTENT(OUT) :: BSOUT(2,36), JAC
      REAL(DOUBLE) :: NVAL(6), DN(2,6), E1F(3), E2F(3), E3F(3), G1(3), G2(3), C(4), BSN(2,36)
      REAL(DOUBLE) :: T0(3), T0I(3), C1(3), C2(3)
      REAL(DOUBLE) :: NORMS_LOC(6,3)
      INTEGER(LONG) :: II, COL
      NORMS_LOC = NORMS
      IF (PRESENT(NORMS_EXT)) NORMS_LOC = NORMS_EXT
      CALL SHAPE_T6(R, S, NVAL, DN)
      CALL FIXED_FRAME_T6(XYZN, E1F, E2F, E3F)
      CALL COV_MAP_T6(XYZN, R, S, E1F, E2F, G1, G2, C, JAC)
      ! v2 Simo1993_Tri6_ShellElement_v2.py uses the UN-normalized interpolated
! director (matches v2 Python `_compute_Bs` which omits `t0 /= np.linalg.norm(t0)`).
! v1p8 normalized T0 to a unit vector; v2 keeps the magnitude so that the
! director-rotation shear strain γ = (∇u + ∇uᵀ)·t + (∇t)·u scales correctly.
      T0 = MATMUL(NVAL, NORMS_LOC)
      BSN = ZERO
      DO II=1,6
         COL = (II-1)*6
         T0I = NORMS_LOC(II,:)
         BSN(1,COL+1:COL+3) = BSN(1,COL+1:COL+3) + DN(1,II)*T0
         BSN(2,COL+1:COL+3) = BSN(2,COL+1:COL+3) + DN(2,II)*T0
         CALL CROSS3(T0I, G1, C1)
         CALL CROSS3(T0I, G2, C2)
         BSN(1,COL+4:COL+6) = BSN(1,COL+4:COL+6) + NVAL(II)*C1
         BSN(2,COL+4:COL+6) = BSN(2,COL+4:COL+6) + NVAL(II)*C2
      ENDDO
      BSOUT(1,:) = C(1)*BSN(1,:) + C(2)*BSN(2,:)
      BSOUT(2,:) = C(3)*BSN(1,:) + C(4)*BSN(2,:)
      END SUBROUTINE BS_T6_AT

      SUBROUTINE BDRILL_T6_AT ( XYZN, NORMS, R, S, BDOUT, JAC, NORMS_EXT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(6,3), NORMS(6,3), R, S
      REAL(DOUBLE), INTENT(IN), OPTIONAL :: NORMS_EXT(6,3)
      REAL(DOUBLE), INTENT(OUT) :: BDOUT(1,36), JAC
      REAL(DOUBLE) :: NVAL(6), DN(2,6), G1(3), G2(3), E1(3), E2(3), E3(3)
      REAL(DOUBLE) :: DR(6), DS(6), A11, A22, A12, DET, AI11, AI22, AI12
      REAL(DOUBLE) :: DXII, DYII
      REAL(DOUBLE) :: NORMS_LOC(6,3)
      INTEGER(LONG) :: II, COL
      NORMS_LOC = NORMS
      IF (PRESENT(NORMS_EXT)) NORMS_LOC = NORMS_EXT
      CALL SHAPE_T6(R, S, NVAL, DN)
      CALL SURFACE_BASIS_T6(XYZN, R, S, G1, G2, E1, E2, E3, JAC)
!     Simo1993_Tri6_ShellElement_v2.py _compute_Bdrill (pointwise drilling):
!       dr, ds = inv(metric [[g·g]]) @ [dNr; dNs]
!       dx,dy  = dr*(g·e_axis) + ds*(g·e_axis)
!       Bd[trans] = 0.5*(dx*E2 - dy*E1)
!       Bd[rot]   = -N*nor_I
      A11 = DOT_PRODUCT(G1,G1); A22 = DOT_PRODUCT(G2,G2)
      A12 = DOT_PRODUCT(G1,G2); DET = A11*A22 - A12*A12
      IF (DABS(DET) <= 1.0D-30) THEN
         DR = ZERO; DS = ZERO
      ELSE
         AI11 =  A22/DET; AI22 =  A11/DET; AI12 = -A12/DET
         DO II=1,6
            DR(II) = AI11*DN(1,II) + AI12*DN(2,II)
            DS(II) = AI12*DN(1,II) + AI22*DN(2,II)
         ENDDO
      ENDIF
      BDOUT = ZERO
      DO II=1,6
         DXII = DR(II)*DOT_PRODUCT(G1,E1) + DS(II)*DOT_PRODUCT(G2,E1)
         DYII = DR(II)*DOT_PRODUCT(G1,E2) + DS(II)*DOT_PRODUCT(G2,E2)
         COL = (II-1)*6
         BDOUT(1,COL+1) = 0.5D0*(DXII*E2(1) - DYII*E1(1))
         BDOUT(1,COL+2) = 0.5D0*(DXII*E2(2) - DYII*E1(2))
         BDOUT(1,COL+3) = 0.5D0*(DXII*E2(3) - DYII*E1(3))
         BDOUT(1,COL+4) = -NVAL(II)*E3(1)
         BDOUT(1,COL+5) = -NVAL(II)*E3(2)
         BDOUT(1,COL+6) = -NVAL(II)*E3(3)
      ENDDO
      END SUBROUTINE BDRILL_T6_AT

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

      END SUBROUTINE CTRIA6_SIMO1993
