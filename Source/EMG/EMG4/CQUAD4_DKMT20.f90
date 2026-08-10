! #################################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________
!
! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________
!
! End MIT license text.

      SUBROUTINE CQUAD4_DKMT20 ( OPT, INT_ELEM_ID )

! DKMT20 shell element.
! This branch follows the DKMQ24-style EAS formulation, but stays compatible
! with the existing MYSTRAN CQUAD4/K6ROT workflow:
!   - drilling DOF is left soft/zero in the kernel so K6ROT can stabilize it
!   - CALC_K6ROT remains an optional additional stabilizer at EMG level
!   - the formulation is kept separate from legacy DKMQ20 so we can tune it
!     independently while preserving the old path

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_GAUSS
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, FOUR
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE MODEL_STUF, ONLY            :  BGRID, EID, ELGP, GRID_SNORM, KE, KED, ME, BE1, BE2, BE3, EM, EB, ET, EPROP,       &
                                         MASS_PER_UNIT_AREA, PRESS, PPE, TE, NUM_EMG_FATAL_ERRS, SHELL_A, SHELL_D, SHELL_T, XEB

      USE ELMDIS_Interface
      USE ELEM_STRE_STRN_ARRAYS_Interface
      USE ORDER_GAUSS_Interface
      USE OUTA_HERE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CQUAD4_DKMT20'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      INTEGER(LONG), PARAMETER        :: NNODE = 4
      INTEGER(LONG)                   :: I, J, K, GP, JSUB
      REAL(DOUBLE)                    :: XYZ(4,3), NORMALS(4,3), T24(24,24), T24T(24,24)
      REAL(DOUBLE)                    :: SS(MAX_ORDER_GAUSS), HH(MAX_ORDER_GAUSS)
      REAL(DOUBLE)                    :: XI, ETA, WT
      REAL(DOUBLE)                    :: TV1(3), TV2(3), NVEC(3), JDET, CO(2,2), BCMAT(2,2)
      REAL(DOUBLE)                    :: RT(3)
      REAL(DOUBLE)                    :: AU(4,24), ADELTA(4,4), AINV_AU(4,24)
      REAL(DOUBLE)                    :: BMB(3,24), BBB(3,24), BSB(2,24)
      REAL(DOUBLE)                    :: BML(3,24), BBL(3,24), BSL(2,24)
      REAL(DOUBLE)                    :: KLOCAL(24,24), KBASIC(24,24), KMEM(24,24), KBEND(24,24), KSHEAR(24,24), KDRILL(24,24)
      REAL(DOUBLE)                    :: KUA(24,4), KAA(4,4), KAAINV(4,4), MEAS(3,4), MEAS_MEAN(3,4)
      REAL(DOUBLE)                    :: MEAS_RAW(3,4,2,2), MEAS_DA(2,2), DA_SUM
      REAL(DOUBLE)                    :: M1(4,4), MBASIC(24,24), MLOCAL(24,24), NVG(4), MDIAG(4)
      REAL(DOUBLE)                    :: UNIT_PPE_B(24), UNIT_PPE_L(24)
      REAL(DOUBLE)                    :: GBE1(3,24,4), GBE2(3,24,4), GBE3(2,24,4)

! **********************************************************************************************************************************

      IF (ELGP /= 4) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9001) SUBR_NAME, EID, ELGP
         WRITE(F06,9001) SUBR_NAME, EID, ELGP
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      XYZ = ZERO
      CALL LOAD_BASIC_COORDS ( XYZ )
      CALL CALC_NODAL_NORMALS ( XYZ, NORMALS )
      CALL BUILD_T24 ( TE, T24 )
      T24T = TRANSPOSE(T24)

      AU = BUILD_AU(XYZ, NORMALS)
      ADELTA = BUILD_ADELTA(XYZ, EPROP(1))
      AINV_AU = ZERO
      DO I=1,4
         IF (DABS(ADELTA(I,I)) > 1.0D-14) THEN
            AINV_AU(I,1:24) = AU(I,1:24) / ADELTA(I,I)
         ENDIF
      ENDDO

      KBASIC = ZERO
      KMEM = ZERO
      KBEND = ZERO
      KSHEAR = ZERO
      KLOCAL = ZERO
      GBE1 = ZERO
      GBE2 = ZERO
      GBE3 = ZERO

      IF ((OPT(3) == 'Y') .OR. (OPT(4) == 'Y') .OR. (OPT(5) == 'Y') .OR. (OPT(1) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         CALL ORDER_GAUSS(2, SS, HH)
      ENDIF

      IF ((OPT(3) == 'Y') .OR. (OPT(4) == 'Y')) THEN
         KUA = ZERO
         KAA = ZERO
         MEAS_MEAN = ZERO
         MEAS_RAW = ZERO
         MEAS_DA = ZERO
         DA_SUM = ZERO

         DO I=1,2
            DO J=1,2
               XI  = SS(I)
               ETA = SS(J)
               WT  = HH(I)*HH(J)
               CALL GEOMETRY_AT(XYZ, NORMALS, XI, ETA, TV1, TV2, NVEC, JDET, CO, BCMAT)
               MEAS_RAW(:,:,I,J) = EAS_AT(XYZ, NORMALS, XI, ETA)
               MEAS_DA(I,J) = WT*JDET
               MEAS_MEAN = MEAS_MEAN + MEAS_RAW(:,:,I,J)*MEAS_DA(I,J)
               DA_SUM = DA_SUM + MEAS_DA(I,J)
            ENDDO
         ENDDO
         IF (DA_SUM > 1.0D-30) THEN
            MEAS_MEAN = MEAS_MEAN / DA_SUM
         ENDIF

         DO I=1,2
            DO J=1,2
               XI  = SS(I)
               ETA = SS(J)
               WT  = HH(I)*HH(J)
               CALL GEOMETRY_AT(XYZ, NORMALS, XI, ETA, TV1, TV2, NVEC, JDET, CO, BCMAT)
               BMB = BM_AT(XI, ETA, TV1, TV2, CO)
               BBB = BB_AT(XYZ, NORMALS, XI, ETA, TV1, TV2, CO, BCMAT, AINV_AU)
               BSB = BS_AT(XYZ, XI, ETA, CO, AINV_AU, EPROP(1))

               IF ((DEBUG(190) > 0) .AND. (I == 1) .AND. (J == 1)) THEN
                  CALL DEBUG_PRINT_MATRIX('CQUAD4_DKMT20 GP11 BMB', BMB)
                  CALL DEBUG_PRINT_MATRIX('CQUAD4_DKMT20 GP11 BBB', BBB)
                  CALL DEBUG_PRINT_MATRIX('CQUAD4_DKMT20 GP11 BSB', BSB)
               ENDIF

               BML = MATMUL(BMB, T24T)
               BBL = MATMUL(BBB, T24T)
               BSL = MATMUL(BSB, T24T)

               KMEM   = KMEM   + WT*JDET*MATMUL(TRANSPOSE(BMB), MATMUL(SHELL_A, BMB))
               MEAS = MEAS_RAW(:,:,I,J) - MEAS_MEAN
               KUA = KUA + WT*JDET*MATMUL(TRANSPOSE(BMB), MATMUL(SHELL_A, MEAS))
               KAA = KAA + WT*JDET*MATMUL(TRANSPOSE(MEAS), MATMUL(SHELL_A, MEAS))
               KBEND  = KBEND  + WT*JDET*MATMUL(TRANSPOSE(BBB), MATMUL(SHELL_D, BBB))
               KSHEAR = KSHEAR + WT*JDET*MATMUL(TRANSPOSE(BSB), MATMUL(SHELL_T, BSB))
               KBASIC = KMEM + KBEND + KSHEAR

               GP = (I - 1)*2 + J
               GBE1(1:3,1:24,GP) = BML
               GBE2(1:3,1:24,GP) = BBL
               GBE3(1:2,1:24,GP) = BSL
            ENDDO
         ENDDO

         KDRILL = ZERO
         CALL INV4(KAA, KAAINV)
         IF (MAXVAL(ABS(KAAINV)) > 0.0D0) THEN
            KMEM = KMEM - MATMUL(KUA, MATMUL(KAAINV, TRANSPOSE(KUA)))
         ENDIF
         KBASIC = KMEM + KBEND + KSHEAR
         KBASIC = KBASIC + KDRILL
         KLOCAL = MATMUL(T24, MATMUL(KBASIC, T24T))

         IF (OPT(4) == 'Y') THEN
            KE(1:24,1:24) = KLOCAL
         ENDIF

         IF (OPT(3) == 'Y') THEN
            BE1(:,:,1) = (GBE1(:,:,1) + GBE1(:,:,2) + GBE1(:,:,3) + GBE1(:,:,4)) / FOUR
            BE2(:,:,1) = (GBE2(:,:,1) + GBE2(:,:,2) + GBE2(:,:,3) + GBE2(:,:,4)) / FOUR
            BE3(1:2,:,1) = (GBE3(1:2,:,1) + GBE3(1:2,:,2) + GBE3(1:2,:,3) + GBE3(1:2,:,4)) / FOUR
            BE1(:,:,2) = GBE1(:,:,4); BE1(:,:,3) = GBE1(:,:,3); BE1(:,:,4) = GBE1(:,:,2); BE1(:,:,5) = GBE1(:,:,1)
            BE2(:,:,2) = GBE2(:,:,4); BE2(:,:,3) = GBE2(:,:,3); BE2(:,:,4) = GBE2(:,:,2); BE2(:,:,5) = GBE2(:,:,1)
            BE3(1:2,:,2) = GBE3(1:2,:,4); BE3(1:2,:,3) = GBE3(1:2,:,3)
            BE3(1:2,:,4) = GBE3(1:2,:,2); BE3(1:2,:,5) = GBE3(1:2,:,1)
         ENDIF
      ENDIF

      IF (OPT(1) == 'Y') THEN
         M1 = ZERO
         DO I=1,2
            DO J=1,2
               XI  = SS(I)
               ETA = SS(J)
               WT  = HH(I)*HH(J)
               CALL GEOMETRY_AT(XYZ, NORMALS, XI, ETA, TV1, TV2, NVEC, JDET, CO, BCMAT)
               CALL SHAPE_N(XI, ETA, NVG)
               DO GP=1,4
                  DO K=1,4
                     M1(GP,K) = M1(GP,K) + NVG(GP)*NVG(K)*MASS_PER_UNIT_AREA*WT*JDET
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         MBASIC = ZERO
         MDIAG = ZERO
         DO I=1,4
            MDIAG(I) = SUM(M1(I,1:4))
         ENDDO
         DO I=1,4
            DO K=1,3
               MBASIC((I-1)*6+K,(I-1)*6+K) = MDIAG(I)
            ENDDO
         ENDDO
         MLOCAL = MATMUL(T24, MATMUL(MBASIC, T24T))
         ME(1:24,1:24) = MLOCAL
      ENDIF

      IF (OPT(5) == 'Y') THEN
         UNIT_PPE_B = ZERO
         DO I=1,2
            DO J=1,2
               XI  = SS(I)
               ETA = SS(J)
               WT  = HH(I)*HH(J)
               CALL GEOMETRY_AT(XYZ, NORMALS, XI, ETA, TV1, TV2, NVEC, JDET, CO, BCMAT)
               CALL SHAPE_N(XI, ETA, NVG)
               DO GP=1,4
                  UNIT_PPE_B((GP-1)*6+1:(GP-1)*6+3) = UNIT_PPE_B((GP-1)*6+1:(GP-1)*6+3) + NVEC*NVG(GP)*WT*JDET
               ENDDO
            ENDDO
         ENDDO
         UNIT_PPE_L = MATMUL(T24, UNIT_PPE_B)
         DO JSUB=1,SIZE(PPE,2)
            PPE(1:24,JSUB) = UNIT_PPE_L(1:24) * PRESS(3,JSUB)
         ENDDO
      ENDIF

      IF (OPT(6) == 'Y') THEN
         CALL ELMDIS
         KED(1:24,1:24) = ZERO
         DO I=1,2
            DO J=1,2
               XI  = SS(I)
               ETA = SS(J)
               WT  = HH(I)*HH(J)
               CALL GEOMETRY_AT(XYZ, NORMALS, XI, ETA, TV1, TV2, NVEC, JDET, CO, BCMAT)
               BMB = BM_AT(XI, ETA, TV1, TV2, CO)
               KED(1:24,1:24) = KED(1:24,1:24) + WT*JDET*MATMUL(TRANSPOSE(BMB), MATMUL(SHELL_A, BMB))
            ENDDO
         ENDDO
      ENDIF

      RETURN

 9001 FORMAT(' *ERROR: ',A,' expects ELGP=4 for element ',I8,' but got ',I8)

! **********************************************************************************************************************************

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

      SUBROUTINE SHAPE_DP ( XI, ETA, DP )
      REAL(DOUBLE), INTENT(IN) :: XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: DP(2,4)
      DP(1,1) = -XI*(ONE - ETA)
      DP(1,2) =  0.5D0*(ONE - ETA*ETA)
      DP(1,3) = -XI*(ONE + ETA)
      DP(1,4) = -0.5D0*(ONE - ETA*ETA)
      DP(2,1) = -0.5D0*(ONE - XI*XI)
      DP(2,2) = -ETA*(ONE + XI)
      DP(2,3) =  0.5D0*(ONE - XI*XI)
      DP(2,4) = -ETA*(ONE - XI)
      END SUBROUTINE SHAPE_DP

      SUBROUTINE CROSS3 ( A, B, C )
      REAL(DOUBLE), INTENT(IN) :: A(3), B(3)
      REAL(DOUBLE), INTENT(OUT) :: C(3)
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)
      END SUBROUTINE CROSS3

      SUBROUTINE INV2 ( A, AINV )
      REAL(DOUBLE), INTENT(IN)  :: A(2,2)
      REAL(DOUBLE), INTENT(OUT) :: AINV(2,2)
      REAL(DOUBLE) :: DET
      DET = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      IF (DABS(DET) < 1.0D-20) THEN
         AINV = ZERO
      ELSE
         AINV(1,1) =  A(2,2)/DET
         AINV(1,2) = -A(1,2)/DET
         AINV(2,1) = -A(2,1)/DET
         AINV(2,2) =  A(1,1)/DET
      ENDIF
      END SUBROUTINE INV2

      SUBROUTINE CALC_NODAL_NORMALS ( XYZN, NORMS )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(4,3)
      REAL(DOUBLE), INTENT(OUT) :: NORMS(4,3)
      INTEGER(LONG) :: II
      REAL(DOUBLE) :: A(3), B(3), N(3), NM, SN(3), SDOT
      A = XYZN(2,:) - XYZN(1,:)
      B = XYZN(4,:) - XYZN(1,:)
      CALL CROSS3(A, B, N)
      NM = VNORM(N)
      IF (NM <= 1.0D-15) THEN
         N = (/ZERO, ZERO, ONE/)
         NM = ONE
      ENDIF
      NORMS(1,:) = N / NM
      A = XYZN(3,:) - XYZN(2,:)
      B = XYZN(1,:) - XYZN(2,:)
      CALL CROSS3(A, B, N)
      NM = VNORM(N)
      IF (NM <= 1.0D-15) THEN
         N = (/ZERO, ZERO, ONE/)
         NM = ONE
      ENDIF
      NORMS(2,:) = N / NM
      A = XYZN(4,:) - XYZN(3,:)
      B = XYZN(2,:) - XYZN(3,:)
      CALL CROSS3(A, B, N)
      NM = VNORM(N)
      IF (NM <= 1.0D-15) THEN
         N = (/ZERO, ZERO, ONE/)
         NM = ONE
      ENDIF
      NORMS(3,:) = N / NM
      A = XYZN(1,:) - XYZN(4,:)
      B = XYZN(3,:) - XYZN(4,:)
      CALL CROSS3(A, B, N)
      NM = VNORM(N)
      IF (NM <= 1.0D-15) THEN
         N = (/ZERO, ZERO, ONE/)
         NM = ONE
      ENDIF
      NORMS(4,:) = N / NM
      IF (ALLOCATED(GRID_SNORM)) THEN
         DO II=1,4
            IF ((BGRID(II) > 0) .AND. (BGRID(II) <= SIZE(GRID_SNORM,1))) THEN
               SN = GRID_SNORM(BGRID(II),:)
               NM = VNORM(SN)
               IF (NM > 1.0D-15) THEN
                  SN = SN / NM
                  SDOT = DOT_PRODUCT(SN, NORMS(II,:))
                  IF (SDOT < 1.0D-2) THEN
                     NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,'(A,A,A,I8,A,I2,A,ES14.6)') ' *ERROR: ', TRIM(SUBR_NAME), ' EID=', EID,                       &
                        ' SNORM AT NODE ', II, ' IS TOO FAR FROM CQUAD4_DKMT20 MIDSURFACE NORMAL. DOT=', SDOT
                     WRITE(F06,'(A,A,A,I8,A,I2,A,ES14.6)') ' *ERROR: ', TRIM(SUBR_NAME), ' EID=', EID,                       &
                        ' SNORM AT NODE ', II, ' IS TOO FAR FROM CQUAD4_DKMT20 MIDSURFACE NORMAL. DOT=', SDOT
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF
                  NORMS(II,:) = SN
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      END SUBROUTINE CALC_NODAL_NORMALS

      SUBROUTINE GEOMETRY_AT ( XYZN, NORMS, XI, ETA, T1, T2, NORMV, JAC, CO, BCM )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(4,3), NORMS(4,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: T1(3), T2(3), NORMV(3), JAC, CO(2,2), BCM(2,2)
      REAL(DOUBLE) :: DN(2,4), A1(3), A2(3), AXB(3), AMAT(2,2), INVA(2,2), A1C(3), A2C(3)
      REAL(DOUBLE) :: NHATXI(3), NHATETA(3), BNHAT(2,2), REF(3), TMP(3), NM
      CALL SHAPE_DN(XI, ETA, DN)
      A1 = MATMUL(DN(1,:), XYZN)
      A2 = MATMUL(DN(2,:), XYZN)
      CALL CROSS3(A1, A2, AXB)
      JAC = VNORM(AXB)
      IF (JAC > 1.0D-15) THEN
         NORMV = AXB / JAC
      ELSE
         NORMV = (/ZERO, ZERO, ONE/)
      ENDIF
      AMAT(1,1) = DOT_PRODUCT(A1, A1)
      AMAT(1,2) = DOT_PRODUCT(A1, A2)
      AMAT(2,1) = AMAT(1,2)
      AMAT(2,2) = DOT_PRODUCT(A2, A2)
      CALL INV2(AMAT, INVA)
      A1C = INVA(1,1)*A1 + INVA(1,2)*A2
      A2C = INVA(2,1)*A1 + INVA(2,2)*A2
      REF = (/ZERO, ZERO, ONE/)
      CALL CROSS3(NORMV, REF, TMP)
      NM = VNORM(TMP)
      IF (NM < 1.0D-10) THEN
         REF = (/ZERO, ONE, ZERO/)
         CALL CROSS3(NORMV, REF, TMP)
         NM = VNORM(TMP)
      ENDIF
      T1 = TMP / NM
      CALL CROSS3(NORMV, T1, T2)
      NM = VNORM(T2)
      IF (NM > 1.0D-15) T2 = T2 / NM
      CO(1,1) = DOT_PRODUCT(A1C, T1)
      CO(1,2) = DOT_PRODUCT(A1C, T2)
      CO(2,1) = DOT_PRODUCT(A2C, T1)
      CO(2,2) = DOT_PRODUCT(A2C, T2)
      NHATXI = MATMUL(DN(1,:), NORMS)
      NHATETA = MATMUL(DN(2,:), NORMS)
      BNHAT(1,1) = DOT_PRODUCT(A2C, NHATETA)
      BNHAT(1,2) = -DOT_PRODUCT(A1C, NHATETA)
      BNHAT(2,1) = -DOT_PRODUCT(A2C, NHATXI)
      BNHAT(2,2) = DOT_PRODUCT(A1C, NHATXI)
      BCM = MATMUL(BNHAT, CO)
      END SUBROUTINE GEOMETRY_AT

      FUNCTION BUILD_AU ( XYZN, NORMS ) RESULT(AUOUT)
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), NORMS(4,3)
      REAL(DOUBLE) :: AUOUT(4,24)
      REAL(DOUBLE) :: XJI(3), LK, TSK(3), NK(3), NORMK, RNI(3,3), RNJ(3,3)
      INTEGER(LONG) :: KK, I1, J1, CI, CJ
      INTEGER(LONG), PARAMETER :: SIDE_I(4) = (/1,2,3,4/)
      INTEGER(LONG), PARAMETER :: SIDE_J(4) = (/2,3,4,1/)
      AUOUT = ZERO
      DO KK=1,4
         I1 = SIDE_I(KK)
         J1 = SIDE_J(KK)
         XJI = XYZN(J1,:) - XYZN(I1,:)
         LK = VNORM(XJI)
         IF (LK < 1.0D-15) CYCLE
         TSK = XJI / LK
         NK = 0.5D0*(NORMS(I1,:) + NORMS(J1,:))
         NORMK = VNORM(NK)
         IF (NORMK > 1.0D-15) NK = NK / NORMK
         CALL RNMAT(NORMS(I1,1), NORMS(I1,2), NORMS(I1,3), RNI)
         CALL RNMAT(NORMS(J1,1), NORMS(J1,2), NORMS(J1,3), RNJ)
         CI = (I1-1)*6
         CJ = (J1-1)*6
         AUOUT(KK,CI+1:CI+3) = -NK / LK
         AUOUT(KK,CJ+1:CJ+3) =  NK / LK
         RT = 0.5D0*MATMUL(TRANSPOSE(RNI), TSK)
         AUOUT(KK,CI+4:CI+5) = RT(1:2)
         RT = 0.5D0*MATMUL(TRANSPOSE(RNJ), TSK)
         AUOUT(KK,CJ+4:CJ+5) = RT(1:2)
         AUOUT(KK,CI+6) = ZERO
         AUOUT(KK,CJ+6) = ZERO
      ENDDO
      END FUNCTION BUILD_AU

      FUNCTION BUILD_ADELTA ( XYZN, THICK ) RESULT(ADELOUT)
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), THICK
      REAL(DOUBLE) :: ADELOUT(4,4)
      INTEGER(LONG), PARAMETER :: SIDE_I(4) = (/1,2,3,4/)
      INTEGER(LONG), PARAMETER :: SIDE_J(4) = (/2,3,4,1/)
      INTEGER(LONG) :: KK, I1, J1
      REAL(DOUBLE) :: LK, PHI
      ADELOUT = ZERO
      DO KK=1,4
         I1 = SIDE_I(KK)
         J1 = SIDE_J(KK)
         LK = VNORM(XYZN(J1,:) - XYZN(I1,:))
         PHI = PHI_SIDE(LK, THICK)
         ADELOUT(KK,KK) = (TWO/3.0D0)*(ONE + PHI)
      ENDDO
      END FUNCTION BUILD_ADELTA

      FUNCTION PHI_SIDE ( LK, THICK ) RESULT(PHI)
      REAL(DOUBLE), INTENT(IN) :: LK, THICK
      REAL(DOUBLE) :: PHI
      IF (LK < 1.0D-15) THEN
         PHI = ZERO
      ELSE
         PHI = 2.0D0/( (5.0D0/6.0D0)*(ONE - EM(1,2)/MAX(EM(1,1),1.0D-30)) ) * (THICK*THICK/(LK*LK))
      ENDIF
      END FUNCTION PHI_SIDE

      FUNCTION BM_AT ( XI, ETA, T1, T2, CO ) RESULT(BMOUT)
      REAL(DOUBLE), INTENT(IN) :: XI, ETA, T1(3), T2(3), CO(2,2)
      REAL(DOUBLE) :: BMOUT(3,24), DN(2,4), NIX(4), NIY(4)
      INTEGER(LONG) :: II, DD, COL
      CALL SHAPE_DN(XI, ETA, DN)
      NIX = DN(1,:)*CO(1,1) + DN(2,:)*CO(2,1)
      NIY = DN(1,:)*CO(1,2) + DN(2,:)*CO(2,2)
      BMOUT = ZERO
      DO II=1,4
         COL = (II-1)*6
         DO DD=1,3
            BMOUT(1,COL+DD) = T1(DD)*NIX(II)
            BMOUT(2,COL+DD) = T2(DD)*NIY(II)
            BMOUT(3,COL+DD) = T1(DD)*NIY(II) + T2(DD)*NIX(II)
         ENDDO
      ENDDO
      END FUNCTION BM_AT

      FUNCTION BB_AT ( XYZN, NORMS, XI, ETA, T1, T2, CO, BCM, AIAU ) RESULT(BBOUT)
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), NORMS(4,3), XI, ETA, T1(3), T2(3), CO(2,2), BCM(2,2), AIAU(4,24)
      REAL(DOUBLE) :: BBOUT(3,24), DN(2,4), DP(2,4), NIX(4), NIY(4), PKX(4), PKY(4), NBC1(4), NBC2(4)
      REAL(DOUBLE) :: BBETA(3,24), BDEL(3,4), RNI(3,3), V1(3), V2(3), XJI(3), LK, TSK(3)
      INTEGER(LONG) :: II, DD, COL
      INTEGER(LONG), PARAMETER :: SIDE_I(4) = (/1,2,3,4/)
      INTEGER(LONG), PARAMETER :: SIDE_J(4) = (/2,3,4,1/)
      CALL SHAPE_DN(XI, ETA, DN)
      CALL SHAPE_DP(XI, ETA, DP)
      NIX = DN(1,:)*CO(1,1) + DN(2,:)*CO(2,1)
      NIY = DN(1,:)*CO(1,2) + DN(2,:)*CO(2,2)
      PKX = DP(1,:)*CO(1,1) + DP(2,:)*CO(2,1)
      PKY = DP(1,:)*CO(1,2) + DP(2,:)*CO(2,2)
      NBC1 = DN(1,:)*BCM(1,1) + DN(2,:)*BCM(2,1)
      NBC2 = DN(1,:)*BCM(1,2) + DN(2,:)*BCM(2,2)
      BBETA = ZERO
      DO II=1,4
         COL = (II-1)*6
         CALL RNMAT(NORMS(II,1), NORMS(II,2), NORMS(II,3), RNI)
         V1 = MATMUL(TRANSPOSE(RNI), T1)
         V2 = MATMUL(TRANSPOSE(RNI), T2)
         DO DD=1,3
            BBETA(1,COL+DD) = T1(DD)*NBC1(II)
            BBETA(2,COL+DD) = T2(DD)*NBC2(II)
            BBETA(3,COL+DD) = T1(DD)*NBC2(II) + T2(DD)*NBC1(II)
         ENDDO
         BBETA(1,COL+4:COL+5) = V1(1:2)*NIX(II)
         BBETA(2,COL+4:COL+5) = V2(1:2)*NIY(II)
         BBETA(3,COL+4:COL+5) = V1(1:2)*NIY(II) + V2(1:2)*NIX(II)
         BBETA(:,COL+6) = ZERO
      ENDDO
      BDEL = ZERO
      DO II=1,4
         XJI = XYZN(SIDE_J(II),:) - XYZN(SIDE_I(II),:)
         LK = VNORM(XJI)
         IF (LK < 1.0D-15) CYCLE
         TSK = XJI / LK
         BDEL(1,II) = DOT_PRODUCT(T1, TSK)*PKX(II)
         BDEL(2,II) = DOT_PRODUCT(T2, TSK)*PKY(II)
         BDEL(3,II) = DOT_PRODUCT(T1, TSK)*PKY(II) + DOT_PRODUCT(T2, TSK)*PKX(II)
      ENDDO
      BBOUT = BBETA + MATMUL(BDEL, AIAU)
      END FUNCTION BB_AT

      FUNCTION BS_AT ( XYZN, XI, ETA, CO, AIAU, THICK ) RESULT(BSOUT)
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), XI, ETA, CO(2,2), AIAU(4,24), THICK
      REAL(DOUBLE) :: BSOUT(2,24), NGAM(2,4), AG(4,4), APHI(4,4), BSG(2,4), LK, PHI, BGDEL(2,4), BSDEL(2,4)
      INTEGER(LONG) :: II
      INTEGER(LONG), PARAMETER :: SIDE_I(4) = (/1,2,3,4/)
      INTEGER(LONG), PARAMETER :: SIDE_J(4) = (/2,3,4,1/)
      REAL(DOUBLE), PARAMETER :: SGN(4) = (/ONE, ONE, -ONE, -ONE/)
      BGDEL = ZERO
      DO II=1,4
         LK = VNORM(XYZN(SIDE_J(II),:) - XYZN(SIDE_I(II),:))
         PHI = PHI_SIDE(LK, THICK)
         IF (II == 1) BGDEL(1,II) = -((ONE-ETA)*LK*PHI) / 6.0D0
         IF (II == 3) BGDEL(1,II) =  ((ONE+ETA)*LK*PHI) / 6.0D0
         IF (II == 2) BGDEL(2,II) = -((ONE+XI )*LK*PHI) / 6.0D0
         IF (II == 4) BGDEL(2,II) =  ((ONE-XI )*LK*PHI) / 6.0D0
      ENDDO
      BSDEL(1,:) = CO(1,1)*BGDEL(1,:) + CO(2,1)*BGDEL(2,:)
      BSDEL(2,:) = CO(1,2)*BGDEL(1,:) + CO(2,2)*BGDEL(2,:)
      BSOUT = MATMUL(BSDEL, AIAU)
      END FUNCTION BS_AT

      FUNCTION EAS_AT ( XYZN, NORMS, XI, ETA ) RESULT(MOUT)
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), NORMS(4,3), XI, ETA
      REAL(DOUBLE) :: MOUT(3,4), MHAT(3,4), T0EAS(3,3)
      REAL(DOUBLE) :: T1C(3), T2C(3), NVC(3), JJC, COC(2,2), BCMC(2,2)
      REAL(DOUBLE) :: I11, I12, I21, I22, DETCO, DETT0
      MHAT = ZERO
      MHAT(1,1) = XI
      MHAT(2,2) = ETA
      MHAT(3,3) = XI
      MHAT(3,4) = ETA
      CALL GEOMETRY_AT(XYZN, NORMS, ZERO, ZERO, T1C, T2C, NVC, JJC, COC, BCMC)
      DETCO = COC(1,1)*COC(2,2) - COC(1,2)*COC(2,1)
      IF (DABS(DETCO) < 1.0D-30) THEN
         MOUT = ZERO
         RETURN
      ENDIF
      DETT0 = ONE / DETCO
      I11 = COC(1,1)
      I12 = COC(1,2)
      I21 = COC(2,1)
      I22 = COC(2,2)
      T0EAS(1,1) = I11*I11 / DETT0
      T0EAS(1,2) = I21*I21 / DETT0
      T0EAS(1,3) = I11*I21 / DETT0
      T0EAS(2,1) = I12*I12 / DETT0
      T0EAS(2,2) = I22*I22 / DETT0
      T0EAS(2,3) = I12*I22 / DETT0
      T0EAS(3,1) = TWO*I11*I12 / DETT0
      T0EAS(3,2) = TWO*I21*I22 / DETT0
      T0EAS(3,3) = (I11*I22 + I12*I21) / DETT0
      MOUT = MATMUL(T0EAS, MHAT)
      END FUNCTION EAS_AT

      SUBROUTINE INV4 ( A, AINV )
      REAL(DOUBLE), INTENT(IN)  :: A(4,4)
      REAL(DOUBLE), INTENT(OUT) :: AINV(4,4)
      REAL(DOUBLE) :: AUG(4,8), PIV, FAC, ROWTMP(8)
      INTEGER(LONG) :: II, JJ, PP
      AUG = ZERO
      AUG(:,1:4) = A
      DO II=1,4
         AUG(II,4+II) = ONE
      ENDDO
      DO II=1,4
         PP = II
         DO JJ=II+1,4
            IF (DABS(AUG(JJ,II)) > DABS(AUG(PP,II))) PP = JJ
         ENDDO
         IF (DABS(AUG(PP,II)) < 1.0D-20) THEN
            AINV = ZERO
            RETURN
         ENDIF
         IF (PP /= II) THEN
            ROWTMP = AUG(II,:)
            AUG(II,:) = AUG(PP,:)
            AUG(PP,:) = ROWTMP
         ENDIF
         PIV = AUG(II,II)
         AUG(II,:) = AUG(II,:) / PIV
         DO JJ=1,4
            IF (JJ == II) CYCLE
            FAC = AUG(JJ,II)
            AUG(JJ,:) = AUG(JJ,:) - FAC*AUG(II,:)
         ENDDO
      ENDDO
      AINV = AUG(:,5:8)
      END SUBROUTINE INV4

      FUNCTION VNORM ( V ) RESULT(NM)
      REAL(DOUBLE), INTENT(IN) :: V(3)
      REAL(DOUBLE) :: NM
      NM = DSQRT(MAX(ZERO, DOT_PRODUCT(V,V)))
      END FUNCTION VNORM

      SUBROUTINE BUILD_T24 ( T3, TOUT )
      REAL(DOUBLE), INTENT(IN) :: T3(3,3)
      REAL(DOUBLE), INTENT(OUT) :: TOUT(24,24)
      INTEGER(LONG) :: II, RR, CC, BASE
      TOUT = ZERO
      DO II=1,4
         BASE = (II-1)*6
         DO RR=1,3
            DO CC=1,3
               TOUT(BASE+RR,   BASE+CC  ) = T3(RR,CC)
               TOUT(BASE+RR+3, BASE+CC+3) = T3(RR,CC)
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE BUILD_T24

      SUBROUTINE RNMAT ( NX, NY, NZ, RN )
      REAL(DOUBLE), INTENT(IN) :: NX, NY, NZ
      REAL(DOUBLE), INTENT(OUT) :: RN(3,3)
      RN = ZERO
      RN(1,2) =  NZ
      RN(1,3) = -NY
      RN(2,1) = -NZ
      RN(2,3) =  NX
      RN(3,1) =  NY
      RN(3,2) = -NX
      END SUBROUTINE RNMAT

      SUBROUTINE DEBUG_PRINT_MATRIX ( TITLE, MAT )
      CHARACTER(LEN=*), INTENT(IN) :: TITLE
      REAL(DOUBLE), INTENT(IN) :: MAT(:,:)
      INTEGER(LONG) :: II
      WRITE(F06,'(A)') ' '
      WRITE(F06,'(A)') TRIM(TITLE)
      DO II=1,SIZE(MAT,1)
         WRITE(F06,'(100(1X,ES15.7))') MAT(II,:)
      ENDDO
      END SUBROUTINE DEBUG_PRINT_MATRIX

      END SUBROUTINE CQUAD4_DKMT20
