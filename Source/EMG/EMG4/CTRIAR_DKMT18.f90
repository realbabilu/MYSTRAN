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
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
! LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
! NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
! _______________________________________________________________________________________________________
!
! End MIT license text.

      SUBROUTINE CTRIAR_DKMT18 ( OPT, INT_ELEM_ID )

! --- shell_renovation begin --- !
! DKMT18 triangular shell now tracks the special SNORM-aware reference path.
! The CTRIAR DKMT18 kernel keeps the existing flat-shell algebra, but the
! local element frame is built from nodal SNORM data when that data exists.
! --- shell_renovation end --- !
! --- cquadr_ctriar_composite begin --- !
! Composite routing contract:
!   - when PCOMP_PROPS = 'Y', SHELL_ABD_MATRICES has already populated the
!     laminate-driven shell matrices.
!   - CTRIAR should keep consuming that laminate basis here through DKMT18, so
!     composite tri shells do not fall back to the generic legacy path.
! --- cquadr_ctriar_composite end --- !

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, THREE, FOUR, SIX, TWELVE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE MODEL_STUF, ONLY            :  EID, ELGP, KE, KED, ME, BE1, BE2, BE3, EPROP, MASS_PER_UNIT_AREA,                       &
                                         NUM_EMG_FATAL_ERRS, SHELL_A, TE, XEB, FCONV, STRESS, BGRID, GRID_SNORM

      USE ELMDIS_Interface
      USE ELEM_STRE_STRN_ARRAYS_Interface
      USE OUTA_HERE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CTRIAR_DKMT18'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      INTEGER(LONG), PARAMETER        :: NNODE = 3
      INTEGER(LONG), PARAMETER        :: NDOFN = 6
      INTEGER(LONG), PARAMETER        :: NDOF = NNODE*NDOFN
      REAL(DOUBLE), PARAMETER         :: HALF = 5.0D-1
      REAL(DOUBLE), PARAMETER         :: KAPPA = 5.0D0/6.0D0

      INTEGER(LONG)                   :: I, J, K, IA, IB, RR, CC
      REAL(DOUBLE)                    :: H, E, NU, AREA, XI0, ETA0
      REAL(DOUBLE)                    :: KFAC, ROT_KG_FAC
      REAL(DOUBLE)                    :: K18(18,18), KMEM(18,18), KBEND(18,18), KSHEAR(18,18), KFICT(18,18)
      REAL(DOUBLE)                    :: SIG0(2,2), DNDX3(3), DNDY3(3), DNP3(2,3), KG18(18,18), KGVAL
      REAL(DOUBLE)                    :: BM(3,18), BB(3,18), BS(2,18), KTMP(18,18)
      REAL(DOUBLE)                    :: T1(3), T2(3), NVEC(3), A1(3), A2(3), A1C(3), A2C(3), NORMALS(3,3)
      REAL(DOUBLE)                    :: JAC, CO(2,2), BC(2,2), ADINV_AU(3,18)
      REAL(DOUBLE)                    :: Hm(3,3), Hb(3,3), Hs(2,2)
      REAL(DOUBLE)                    :: MASS_NODE

      IF (ELGP /= 3) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9001) SUBR_NAME, EID, ELGP
         WRITE(F06,9001) SUBR_NAME, EID, ELGP
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      H = EPROP(1)
      IF (H <= 0.0D0) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9002) SUBR_NAME, EID, H
         WRITE(F06,9002) SUBR_NAME, EID, H
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL BUILD_GEOMETRY_AT ( T1, T2, NVEC, JAC, CO, BC, NORMALS, AREA )
      IF (AREA <= 1.0D-14) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9003) SUBR_NAME, EID, AREA
         WRITE(F06,9003) SUBR_NAME, EID, AREA
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IF (DABS(SHELL_A(1,1)) <= 1.0D-14) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9004) SUBR_NAME, EID
         WRITE(F06,9004) SUBR_NAME, EID
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      NU = SHELL_A(1,2) / SHELL_A(1,1)
      E  = SHELL_A(1,1) * (ONE - NU*NU) / H

      CALL BUILD_CONSTITUTIVE ( E, NU, H, Hb, Hm, Hs )
      CALL BUILD_ADELTA_INV_AU ( H, NU, NORMALS, ADINV_AU )
      CALL BUILD_MAKNUN_K_STAGES ( H, E, NU, T1, T2, NVEC, JAC, CO, BC, NORMALS, ADINV_AU, Hm, Hb, Hs,                   &
                                   KMEM, KBEND, KSHEAR, KFICT, K18 )

      TE = ZERO
      TE(1,1) = ONE
      TE(2,2) = ONE
      TE(3,3) = ONE
      KE(1:18,1:18) = K18

      IF (OPT(3) == 'Y') THEN
         XI0  = ONE/THREE
         ETA0 = ONE/THREE
         CALL BUILD_STAGE_B_MAKNUN ( T1, T2, CO, BM )
         CALL BUILD_STAGE_BB_MAKNUN ( XI0, ETA0, T1, T2, NVEC, CO, BC, NORMALS, ADINV_AU, BB )
         CALL BUILD_STAGE_BS_MAKNUN ( T1, T2, ADINV_AU, BS )
         BE1(1:3,1:18,1) = BM
         BE2(1:3,1:18,1) = BB
         BE3(1:2,1:18,1) = BS
      ENDIF

      IF (OPT(1) == 'Y') THEN
         ME(1:18,1:18) = ZERO
         MASS_NODE = MASS_PER_UNIT_AREA * AREA / THREE
         DO I=0,2
            DO J=1,3
               ME(6*I+J,6*I+J) = MASS_NODE
            ENDDO
         ENDDO
      ENDIF

! --- shell_renovation begin --- !
! Differential stiffness for linear buckling. Harahap 2021 Eq. (15)
! includes gradients of u, v, w and rotations beta_x, beta_y.  Since SIG0 is
! handled here as a membrane resultant, the rotational terms carry h^2/12.
! --- shell_renovation end --- !
      IF (OPT(6) == 'Y') THEN
         CALL BUILD_STAGE_B_MAKNUN ( T1, T2, CO, BM )
         BE1(1:3,1:18,1) = BM
         CALL ELMDIS
         CALL ELEM_STRE_STRN_ARRAYS ( 1 )

         SIG0(1,1) = FCONV(1)*STRESS(1)
         SIG0(2,2) = FCONV(1)*STRESS(2)
         SIG0(1,2) = FCONV(1)*STRESS(3)
         SIG0(2,1) = SIG0(1,2)

         IF ((DEBUG(233) > 0) .AND. (EID <= 8)) THEN
            WRITE(F06,'(A,I8,A,3(1X,ES15.7))') 'CTRIAR KGGD EID=', EID, ' SIG0=', SIG0(1,1), SIG0(2,2), SIG0(1,2)
            WRITE(F06,'(A,I8,A,3(1X,ES15.7))') 'CTRIAR KGGD EID=', EID, ' NORMAL=', NVEC(1), NVEC(2), NVEC(3)
         ENDIF

         DNP3 = ZERO
         DNP3(1,1) = -ONE
         DNP3(1,2) =  ONE
         DNP3(2,1) = -ONE
         DNP3(2,3) =  ONE
         DO I=1,3
            DNDX3(I) = DNP3(1,I)*CO(1,1) + DNP3(2,I)*CO(2,1)
            DNDY3(I) = DNP3(1,I)*CO(1,2) + DNP3(2,I)*CO(2,2)
         ENDDO

         KG18 = ZERO
         ROT_KG_FAC = H*H/TWELVE
         DO IA=1,3
            DO IB=1,3
               KGVAL = AREA*( DNDX3(IA)*(SIG0(1,1)*DNDX3(IB) + SIG0(1,2)*DNDY3(IB)) +                         &
                              DNDY3(IA)*(SIG0(2,1)*DNDX3(IB) + SIG0(2,2)*DNDY3(IB)) )
               DO RR=1,3
                  KG18(6*(IA-1)+RR,6*(IB-1)+RR) = KG18(6*(IA-1)+RR,6*(IB-1)+RR) + KGVAL
               ENDDO
               KG18(6*(IA-1)+4,6*(IB-1)+4) = KG18(6*(IA-1)+4,6*(IB-1)+4) + ROT_KG_FAC*KGVAL
               KG18(6*(IA-1)+5,6*(IB-1)+5) = KG18(6*(IA-1)+5,6*(IB-1)+5) + ROT_KG_FAC*KGVAL
            ENDDO
         ENDDO
         KED(1:18,1:18) = KG18
         IF ((DEBUG(233) > 0) .AND. (EID <= 8)) THEN
            WRITE(F06,'(A,I8,A,ES15.7)') 'CTRIAR KGGD EID=', EID, ' KED_NORM=', DSQRT(SUM(KED(1:18,1:18)*KED(1:18,1:18)))
         ENDIF
      ENDIF

      IF (DEBUG(233) > 0) THEN
         XI0  = ONE/THREE
         ETA0 = ONE/THREE
         CALL BUILD_STAGE_B_MAKNUN ( T1, T2, CO, BM )
         CALL DEBUG_PRINT_MATRIX('CTRIAR BM PRE', BM)
         CALL DEBUG_PRINT_MATRIX('CTRIAR Hm', Hm)
         CALL DEBUG_PRINT_MATRIX('CTRIAR Hb', Hb)
         CALL DEBUG_PRINT_MATRIX('CTRIAR Hs', Hs)
         CALL BUILD_STAGE_BB_MAKNUN ( XI0, ETA0, T1, T2, NVEC, CO, BC, NORMALS, ADINV_AU, BB )
         CALL BUILD_STAGE_BS_MAKNUN ( T1, T2, ADINV_AU, BS )
         CALL DEBUG_PRINT_MATRIX('CTRIAR BM POST', BM)
         CALL DEBUG_PRINT_MATRIX('CTRIAR BB', BB)
         CALL DEBUG_PRINT_MATRIX('CTRIAR BS', BS)
         CALL DEBUG_PRINT_MATRIX('CTRIAR TRIAD', RESHAPE((/ T1, T2, NVEC /), (/3,3/) ))
         CALL DEBUG_PRINT_MATRIX('CTRIAR K18 BASIC', K18)
         CALL DEBUG_PRINT_MATRIX('CTRIAR KE LOCAL', KE(1:18,1:18))
         CALL DEBUG_PRINT_MATRIX('CTRIAR CO', CO)
         CALL DEBUG_PRINT_MATRIX('CTRIAR ADINV_AU', ADINV_AU)
         CALL DEBUG_PRINT_MATRIX('CTRIAR BM', BM)
         CALL DEBUG_PRINT_MATRIX('CTRIAR KMEM', KMEM)
         CALL DEBUG_PRINT_MATRIX('CTRIAR KBEND', KBEND)
         CALL DEBUG_PRINT_MATRIX('CTRIAR KSHEAR', KSHEAR)
         CALL DEBUG_PRINT_MATRIX('CTRIAR KFICT', KFICT)
      ENDIF

      RETURN

 9001 FORMAT(' *ERROR: ',A,' expects ELGP=3 for element ',I8,' but got ',I8)
 9002 FORMAT(' *ERROR: ',A,' got nonpositive thickness for element ',I8,': ',1ES14.6)
 9003 FORMAT(' *ERROR: ',A,' got degenerate area for element ',I8,': ',1ES14.6)
 9004 FORMAT(' *ERROR: ',A,' could not derive isotropic shell constants from SHELL_A for element ',I8)

      CONTAINS

      SUBROUTINE BUILD_GEOMETRY_AT ( T1O, T2O, NVO, JACO, COO, BCO, NORMSO, AREAO )
      REAL(DOUBLE), INTENT(OUT) :: T1O(3), T2O(3), NVO(3), JACO, COO(2,2), BCO(2,2), NORMSO(3,3), AREAO
      REAL(DOUBLE)              :: DNP(2,3), XYZ(3,3), A1L(3), A2L(3), AXB(3), A11, A12, A22, DET
      REAL(DOUBLE)              :: INVA(2,2)
      REAL(DOUBLE)              :: A1CLOC(3), A2CLOC(3), KHAT(3), YTMP(3), NXI(3), NETA(3), BNHAT(2,2)
      REAL(DOUBLE)              :: BN11, BN12, BN21, BN22
      INTEGER(LONG)             :: II

      DNP = ZERO
      DNP(1,1) = -ONE
      DNP(1,2) =  ONE
      DNP(2,1) = -ONE
      DNP(2,3) =  ONE

      DO II=1,3
         XYZ(II,1:3) = XEB(II,1:3)
      ENDDO

      A1L = MATMUL(DNP(1,1:3), XYZ)
      A2L = MATMUL(DNP(2,1:3), XYZ)
      CALL CROSS3 ( A1L, A2L, AXB )
      JACO = DSQRT(DOT_PRODUCT(AXB,AXB))

      A11 = DOT_PRODUCT(A1L,A1L)
      A12 = DOT_PRODUCT(A1L,A2L)
      A22 = DOT_PRODUCT(A2L,A2L)
      DET = A11*A22 - A12*A12
      INVA(1,1) =  A22 / DET
      INVA(1,2) = -A12 / DET
      INVA(2,1) = -A12 / DET
      INVA(2,2) =  A11 / DET
      A1CLOC = INVA(1,1)*A1L + INVA(1,2)*A2L
      A2CLOC = INVA(2,1)*A1L + INVA(2,2)*A2L

      NVO = AXB / JACO
      DO II=1,3
         NORMSO(II,1:3) = GRID_SNORM(BGRID(II),1:3)
         IF (DSQRT(DOT_PRODUCT(NORMSO(II,1:3), NORMSO(II,1:3))) <= 1.0D-14) THEN
            NORMSO(II,1:3) = NVO
         ELSE
            NORMSO(II,1:3) = NORMSO(II,1:3) / DSQRT(DOT_PRODUCT(NORMSO(II,1:3), NORMSO(II,1:3)))
            IF (DOT_PRODUCT(NORMSO(II,1:3), NVO) < ZERO) NORMSO(II,1:3) = -NORMSO(II,1:3)
         ENDIF
      ENDDO
      KHAT = ZERO
      KHAT(3) = ONE
      T1O(1) = NVO(2)*KHAT(3) - NVO(3)*KHAT(2)
      T1O(2) = NVO(3)*KHAT(1) - NVO(1)*KHAT(3)
      T1O(3) = NVO(1)*KHAT(2) - NVO(2)*KHAT(1)
      IF (DSQRT(DOT_PRODUCT(T1O,T1O)) <= 1.0D-10) THEN
         YTMP = ZERO
         YTMP(2) = ONE
         T1O(1) = NVO(2)*YTMP(3) - NVO(3)*YTMP(2)
         T1O(2) = NVO(3)*YTMP(1) - NVO(1)*YTMP(3)
         T1O(3) = NVO(1)*YTMP(2) - NVO(2)*YTMP(1)
      ENDIF
      T1O = T1O / DSQRT(DOT_PRODUCT(T1O,T1O))
      T2O(1) = NVO(2)*T1O(3) - NVO(3)*T1O(2)
      T2O(2) = NVO(3)*T1O(1) - NVO(1)*T1O(3)
      T2O(3) = NVO(1)*T1O(2) - NVO(2)*T1O(1)
      T2O = T2O / DSQRT(DOT_PRODUCT(T2O,T2O))

      COO(1,1) = DOT_PRODUCT(A1CLOC, T1O)
      COO(1,2) = DOT_PRODUCT(A1CLOC, T2O)
      COO(2,1) = DOT_PRODUCT(A2CLOC, T1O)
      COO(2,2) = DOT_PRODUCT(A2CLOC, T2O)

      NXI  = MATMUL(DNP(1,1:3), NORMSO)
      NETA = MATMUL(DNP(2,1:3), NORMSO)
      BN11 = DOT_PRODUCT(A1CLOC, NXI)
      BN12 = DOT_PRODUCT(A1CLOC, NETA)
      BN21 = DOT_PRODUCT(A2CLOC, NXI)
      BN22 = DOT_PRODUCT(A2CLOC, NETA)
      BNHAT(1,1) = BN22
      BNHAT(1,2) = -BN12
      BNHAT(2,1) = -BN21
      BNHAT(2,2) = BN11
      BCO = MATMUL(BNHAT, COO)

      AREAO = HALF * JACO
      END SUBROUTINE BUILD_GEOMETRY_AT

      SUBROUTINE CROSS3 ( A, B, C )
      REAL(DOUBLE), INTENT(IN)  :: A(3), B(3)
      REAL(DOUBLE), INTENT(OUT) :: C(3)
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)
      END SUBROUTINE CROSS3

      SUBROUTINE BUILD_ADELTA_INV_AU ( THICK, POIS, NORMSO, AOUT )
      REAL(DOUBLE), INTENT(IN)  :: THICK, POIS, NORMSO(3,3)
      REAL(DOUBLE), INTENT(OUT) :: AOUT(3,18)
      REAL(DOUBLE)              :: LOUT(3), PHIOUT(3), DIAG(3), AU(3,18), DX(3), TSK(3), NK(3), EN(3), RNI(3,3), RNJ(3,3), NVI(3), NVJ(3)
      INTEGER(LONG)             :: II, I, J

      CALL BUILD_SIDE_LENGTHS ( LOUT )
      CALL BUILD_SIDE_PHI ( THICK, POIS, LOUT, PHIOUT )

      EN = ZERO
      DX = XEB(2,1:3) - XEB(1,1:3)
      CALL CROSS3 ( XEB(2,1:3) - XEB(1,1:3), XEB(3,1:3) - XEB(1,1:3), EN )
      EN = EN / DSQRT(DOT_PRODUCT(EN,EN))

      AU = ZERO
      DO II=1,3
         I = II
         IF (II == 1) THEN
            I = 1
            J = 2
         ELSE IF (II == 2) THEN
            I = 2
            J = 3
         ELSE
            I = 3
            J = 1
         ENDIF
         DX = XEB(J,1:3) - XEB(I,1:3)
         TSK = DX / DSQRT(DOT_PRODUCT(DX,DX))
         NK = 0.5D0 * (NORMSO(I,1:3) + NORMSO(J,1:3))
         IF (DSQRT(DOT_PRODUCT(NK,NK)) <= 1.0D-14) NK = EN
         NK = NK / DSQRT(DOT_PRODUCT(NK,NK))
         NVI = NORMSO(I,1:3)
         NVJ = NORMSO(J,1:3)
         CALL BUILD_RN_MATRIX ( NVI, RNI )
         CALL BUILD_RN_MATRIX ( NVJ, RNJ )
         AU(II,6*(I-1)+1:6*(I-1)+3) = AU(II,6*(I-1)+1:6*(I-1)+3) - NK / LOUT(II)
         AU(II,6*(J-1)+1:6*(J-1)+3) = AU(II,6*(J-1)+1:6*(J-1)+3) + NK / LOUT(II)
         AU(II,6*(I-1)+4:6*(I-1)+6) = AU(II,6*(I-1)+4:6*(I-1)+6) + 0.5D0 * MATMUL(TRANSPOSE(RNI), TSK)
         AU(II,6*(J-1)+4:6*(J-1)+6) = AU(II,6*(J-1)+4:6*(J-1)+6) + 0.5D0 * MATMUL(TRANSPOSE(RNJ), TSK)
      ENDDO

      DO II=1,3
         DIAG(II) = -(TWO/THREE) * (ONE + PHIOUT(II))
         AOUT(II,1:18) = AU(II,1:18) / DIAG(II)
      ENDDO
      END SUBROUTINE BUILD_ADELTA_INV_AU

      SUBROUTINE BUILD_SIDE_LENGTHS ( LOUT )
      REAL(DOUBLE), INTENT(OUT) :: LOUT(3)
      LOUT(1) = DSQRT(DOT_PRODUCT(XEB(2,1:3) - XEB(1,1:3), XEB(2,1:3) - XEB(1,1:3)))
      LOUT(2) = DSQRT(DOT_PRODUCT(XEB(3,1:3) - XEB(2,1:3), XEB(3,1:3) - XEB(2,1:3)))
      LOUT(3) = DSQRT(DOT_PRODUCT(XEB(1,1:3) - XEB(3,1:3), XEB(1,1:3) - XEB(3,1:3)))
      END SUBROUTINE BUILD_SIDE_LENGTHS

      SUBROUTINE BUILD_SIDE_PHI ( THICK, POIS, LOUT, PHIOUT )
      REAL(DOUBLE), INTENT(IN)  :: THICK, POIS, LOUT(3)
      REAL(DOUBLE), INTENT(OUT) :: PHIOUT(3)
      INTEGER(LONG)             :: II
      DO II=1,3
         IF (LOUT(II) > 1.0D-14) THEN
            PHIOUT(II) = TWO/(KAPPA*(ONE-POIS)) * (THICK/LOUT(II))**2
         ELSE
            PHIOUT(II) = ZERO
         ENDIF
      ENDDO
      END SUBROUTINE BUILD_SIDE_PHI

      SUBROUTINE BUILD_RN_MATRIX ( NV, RNM )
      REAL(DOUBLE), INTENT(IN)  :: NV(3)
      REAL(DOUBLE), INTENT(OUT) :: RNM(3,3)
      RNM(1,1) = ZERO
      RNM(1,2) =  NV(3)
      RNM(1,3) = -NV(2)
      RNM(2,1) = -NV(3)
      RNM(2,2) = ZERO
      RNM(2,3) =  NV(1)
      RNM(3,1) =  NV(2)
      RNM(3,2) = -NV(1)
      RNM(3,3) = ZERO
      END SUBROUTINE BUILD_RN_MATRIX

      SUBROUTINE BUILD_MAKNUN_K_STAGES ( THICK, YOUNG, POIS, T1, T2, NVO, JACO, COO, BCO, NORMSO, ADINV_AU, HMO, HBO, HSO,   &
                                         KMEMO, KBENDO, KSHEARO, KFICTO, K18O )
      REAL(DOUBLE), INTENT(IN)  :: THICK, YOUNG, POIS, T1(3), T2(3), NVO(3), JACO, COO(2,2), BCO(2,2), NORMSO(3,3), ADINV_AU(3,18)
      REAL(DOUBLE), INTENT(IN)  :: HMO(3,3), HBO(3,3), HSO(2,2)
      REAL(DOUBLE), INTENT(OUT) :: KMEMO(18,18), KBENDO(18,18), KSHEARO(18,18), KFICTO(18,18), K18O(18,18)
      REAL(DOUBLE)              :: XI, ETA, WT, BM(3,18), BB(3,18), BS(2,18), FAC, KTMP(18,18)

      KMEMO = ZERO
      KBENDO = ZERO
      KSHEARO = ZERO
      KFICTO = ZERO
      K18O = ZERO

      CALL BUILD_FICTITIOUS_K_MAKNUN ( THICK, YOUNG, POIS, NORMSO, HMO, T1, T2, JACO, COO, BCO, KFICTO )

      CALL PLATE_GP ( 1, XI, ETA, WT )
         CALL BUILD_STAGE_B_MAKNUN ( T1, T2, COO, BM )
         CALL BUILD_STAGE_BB_MAKNUN ( XI, ETA, T1, T2, NVO, COO, BCO, NORMSO, ADINV_AU, BB )
         CALL BUILD_STAGE_BS_MAKNUN ( T1, T2, ADINV_AU, BS )
         FAC = WT * JACO
         KMEMO = KMEMO + FAC * MATMUL(TRANSPOSE(BM), MATMUL(HMO, BM))
         KTMP(1:3,1:18) = MATMUL(HBO, BB)
         KBENDO = KBENDO + FAC * MATMUL(TRANSPOSE(BB), KTMP(1:3,1:18))
         KTMP(1:2,1:18) = MATMUL(HSO, BS)
         KSHEARO = KSHEARO + FAC * MATMUL(TRANSPOSE(BS), KTMP(1:2,1:18))
         IF (DEBUG(233) > 0) THEN
            CALL DEBUG_PRINT_MATRIX('CTRIAR GP1 KMEM', FAC * MATMUL(TRANSPOSE(BM), MATMUL(HMO, BM)))
            CALL DEBUG_PRINT_MATRIX('CTRIAR GP1 KBEND', FAC * MATMUL(TRANSPOSE(BB), MATMUL(HBO, BB)))
            CALL DEBUG_PRINT_MATRIX('CTRIAR GP1 KSHEAR', FAC * MATMUL(TRANSPOSE(BS), MATMUL(HSO, BS)))
         ENDIF

      CALL PLATE_GP ( 2, XI, ETA, WT )
      CALL BUILD_STAGE_B_MAKNUN ( T1, T2, COO, BM )
      CALL BUILD_STAGE_BB_MAKNUN ( XI, ETA, T1, T2, NVO, COO, BCO, NORMSO, ADINV_AU, BB )
      CALL BUILD_STAGE_BS_MAKNUN ( T1, T2, ADINV_AU, BS )
      FAC = WT * JACO
      KMEMO = KMEMO + FAC * MATMUL(TRANSPOSE(BM), MATMUL(HMO, BM))
         KTMP(1:3,1:18) = MATMUL(HBO, BB)
         KBENDO = KBENDO + FAC * MATMUL(TRANSPOSE(BB), KTMP(1:3,1:18))
         KTMP(1:2,1:18) = MATMUL(HSO, BS)
         KSHEARO = KSHEARO + FAC * MATMUL(TRANSPOSE(BS), KTMP(1:2,1:18))
      IF (DEBUG(233) > 0) THEN
         CALL DEBUG_PRINT_MATRIX('CTRIAR GP2 KMEM', FAC * MATMUL(TRANSPOSE(BM), MATMUL(HMO, BM)))
         CALL DEBUG_PRINT_MATRIX('CTRIAR GP2 KBEND', FAC * MATMUL(TRANSPOSE(BB), MATMUL(HBO, BB)))
         CALL DEBUG_PRINT_MATRIX('CTRIAR GP2 KSHEAR', FAC * MATMUL(TRANSPOSE(BS), MATMUL(HSO, BS)))
      ENDIF

      CALL PLATE_GP ( 3, XI, ETA, WT )
      CALL BUILD_STAGE_B_MAKNUN ( T1, T2, COO, BM )
      CALL BUILD_STAGE_BB_MAKNUN ( XI, ETA, T1, T2, NVO, COO, BCO, NORMSO, ADINV_AU, BB )
      CALL BUILD_STAGE_BS_MAKNUN ( T1, T2, ADINV_AU, BS )
      FAC = WT * JACO
      KMEMO = KMEMO + FAC * MATMUL(TRANSPOSE(BM), MATMUL(HMO, BM))
         KTMP(1:3,1:18) = MATMUL(HBO, BB)
         KBENDO = KBENDO + FAC * MATMUL(TRANSPOSE(BB), KTMP(1:3,1:18))
         KTMP(1:2,1:18) = MATMUL(HSO, BS)
         KSHEARO = KSHEARO + FAC * MATMUL(TRANSPOSE(BS), KTMP(1:2,1:18))
      IF (DEBUG(233) > 0) THEN
         CALL DEBUG_PRINT_MATRIX('CTRIAR GP3 KMEM', FAC * MATMUL(TRANSPOSE(BM), MATMUL(HMO, BM)))
         CALL DEBUG_PRINT_MATRIX('CTRIAR GP3 KBEND', FAC * MATMUL(TRANSPOSE(BB), MATMUL(HBO, BB)))
         CALL DEBUG_PRINT_MATRIX('CTRIAR GP3 KSHEAR', FAC * MATMUL(TRANSPOSE(BS), MATMUL(HSO, BS)))
      ENDIF

      K18O = KMEMO + KBENDO + KSHEARO + KFICTO
      K18O = HALF * (K18O + TRANSPOSE(K18O))
      END SUBROUTINE BUILD_MAKNUN_K_STAGES

      SUBROUTINE BUILD_STAGE_B_MAKNUN ( T1, T2, CO, BMO )
      REAL(DOUBLE), INTENT(IN)  :: T1(3), T2(3), CO(2,2)
      REAL(DOUBLE), INTENT(OUT) :: BMO(3,18)
      REAL(DOUBLE)              :: DNP(2,3), NIX(3), NIY(3)
      INTEGER(LONG)             :: I

      DNP = ZERO
      DNP(1,1) = -ONE
      DNP(1,2) =  ONE
      DNP(2,1) = -ONE
      DNP(2,3) =  ONE

      NIX = DNP(1,1:3)*CO(1,1) + DNP(2,1:3)*CO(2,1)
      NIY = DNP(1,1:3)*CO(1,2) + DNP(2,1:3)*CO(2,2)

      BMO = ZERO
      DO I=1,3
         BMO(1,6*(I-1)+1:6*(I-1)+3) = T1(1:3) * NIX(I)
         BMO(2,6*(I-1)+1:6*(I-1)+3) = T2(1:3) * NIY(I)
         BMO(3,6*(I-1)+1:6*(I-1)+3) = T1(1:3) * NIY(I) + T2(1:3) * NIX(I)
      ENDDO
      END SUBROUTINE BUILD_STAGE_B_MAKNUN

      SUBROUTINE BUILD_STAGE_BB_MAKNUN ( XI, ETA, T1, T2, NVO, CO, BCO, NORMSO, ADINV_AU, BBO )
      REAL(DOUBLE), INTENT(IN)  :: XI, ETA, T1(3), T2(3), NVO(3), CO(2,2), BCO(2,2), NORMSO(3,3), ADINV_AU(3,18)
      REAL(DOUBLE), INTENT(OUT) :: BBO(3,18)
      REAL(DOUBLE)              :: DNP(2,3), DPP(2,3), NIX(3), NIY(3), PKX(3), PKY(3), NBC1(3), NBC2(3)
      REAL(DOUBLE)              :: BBETA(3,18), BBD(3,3), RNI(3,3), V1(3), V2(3), TSK(3), NVI(3), T1S, T2S
      INTEGER(LONG)             :: I, K

      DNP = ZERO
      DNP(1,1) = -ONE
      DNP(1,2) =  ONE
      DNP(2,1) = -ONE
      DNP(2,3) =  ONE
      CALL SHAPE_DP ( XI, ETA, DPP )

      NIX = DNP(1,1:3)*CO(1,1) + DNP(2,1:3)*CO(2,1)
      NIY = DNP(1,1:3)*CO(1,2) + DNP(2,1:3)*CO(2,2)
      PKX = DPP(1,1:3)*CO(1,1) + DPP(2,1:3)*CO(2,1)
      PKY = DPP(1,1:3)*CO(1,2) + DPP(2,1:3)*CO(2,2)
      NBC1 = DNP(1,1:3)*BCO(1,1) + DNP(2,1:3)*BCO(2,1)
      NBC2 = DNP(1,1:3)*BCO(1,2) + DNP(2,1:3)*BCO(2,2)

      BBETA = ZERO
      DO I=1,3
         K = 6*(I-1)
         NVI = NORMSO(I,1:3)
         CALL BUILD_RN_MATRIX ( NVI, RNI )
         V1 = MATMUL(TRANSPOSE(RNI), T1)
         V2 = MATMUL(TRANSPOSE(RNI), T2)
         BBETA(1,K+1:K+3) = T1(1:3) * NBC1(I)
         BBETA(2,K+1:K+3) = T2(1:3) * NBC2(I)
         BBETA(3,K+1:K+3) = T1(1:3) * NBC2(I) + T2(1:3) * NBC1(I)
         BBETA(1,K+4:K+6) = V1 * NIX(I)
         BBETA(2,K+4:K+6) = V2 * NIY(I)
         BBETA(3,K+4:K+6) = V1 * NIY(I) + V2 * NIX(I)
      ENDDO

      BBD = ZERO
      DO I=1,3
         IF (I == 1) THEN
            TSK = XEB(2,1:3) - XEB(1,1:3)
         ELSE IF (I == 2) THEN
            TSK = XEB(3,1:3) - XEB(2,1:3)
         ELSE
            TSK = XEB(1,1:3) - XEB(3,1:3)
         ENDIF
         TSK = TSK / DSQRT(DOT_PRODUCT(TSK,TSK))
         T1S = DOT_PRODUCT(T1, TSK)
         T2S = DOT_PRODUCT(T2, TSK)
         BBD(1,I) = T1S * PKX(I)
         BBD(2,I) = T2S * PKY(I)
         BBD(3,I) = T1S * PKY(I) + T2S * PKX(I)
      ENDDO
      BBO = BBETA + MATMUL(BBD, ADINV_AU)
      END SUBROUTINE BUILD_STAGE_BB_MAKNUN

      SUBROUTINE BUILD_STAGE_BS_MAKNUN ( T1, T2, ADINV_AU, BSO )
      REAL(DOUBLE), INTENT(IN)  :: T1(3), T2(3), ADINV_AU(3,18)
      REAL(DOUBLE), INTENT(OUT) :: BSO(2,18)
      REAL(DOUBLE)              :: COUT(3), SOUT(3), BSG(2,3), LOUT(3), PHIOUT(3), TSK(3), TMP(2,3)
      INTEGER(LONG)             :: I

      DO I=1,3
         IF (I == 1) THEN
            TSK = XEB(2,1:3) - XEB(1,1:3)
         ELSE IF (I == 2) THEN
            TSK = XEB(3,1:3) - XEB(2,1:3)
         ELSE
            TSK = XEB(1,1:3) - XEB(3,1:3)
         ENDIF
         TSK = TSK / DSQRT(DOT_PRODUCT(TSK,TSK))
         COUT(I) = DOT_PRODUCT(T1, TSK)
         SOUT(I) = DOT_PRODUCT(T2, TSK)
      ENDDO
      CALL BUILD_SIDE_LENGTHS ( LOUT )
      CALL BUILD_SIDE_PHI ( EPROP(1), SHELL_A(1,2)/SHELL_A(1,1), LOUT, PHIOUT )
      CALL BUILD_BSG_FROM_CS ( COUT, SOUT, BSG )
      TMP = BSG
      DO I=1,3
         TMP(1,I) = TMP(1,I) * PHIOUT(I)
         TMP(2,I) = TMP(2,I) * PHIOUT(I)
      ENDDO
      BSO = MATMUL(TMP, ADINV_AU)
      END SUBROUTINE BUILD_STAGE_BS_MAKNUN

      SUBROUTINE BUILD_BSG_FROM_CS ( COUT, SOUT, BSGO )
      REAL(DOUBLE), INTENT(IN)  :: COUT(3), SOUT(3)
      REAL(DOUBLE), INTENT(OUT) :: BSGO(2,3)
      REAL(DOUBLE)              :: A1, A2, A3
      A1 = COUT(1)*SOUT(3) - COUT(3)*SOUT(1)
      A2 = COUT(2)*SOUT(1) - COUT(1)*SOUT(2)
      A3 = COUT(3)*SOUT(2) - COUT(2)*SOUT(3)
      BSGO(1,1) = (SOUT(2)/A2 - SOUT(3)/A1)
      BSGO(1,2) = (SOUT(3)/A3 - SOUT(1)/A2)
      BSGO(1,3) = (SOUT(1)/A1 - SOUT(2)/A3)
      BSGO(2,1) = (COUT(3)/A1 - COUT(2)/A2)
      BSGO(2,2) = (COUT(1)/A2 - COUT(3)/A3)
      BSGO(2,3) = (COUT(2)/A3 - COUT(1)/A1)
      BSGO = (TWO/THREE) * BSGO
      END SUBROUTINE BUILD_BSG_FROM_CS

      SUBROUTINE BUILD_FICTITIOUS_K_MAKNUN ( THICK, YOUNG, POIS, NORMSO, HMO, T1, T2, JACO, COO, BCO, KFO )
      REAL(DOUBLE), INTENT(IN)  :: THICK, YOUNG, POIS, NORMSO(3,3), HMO(3,3), T1(3), T2(3), JACO, COO(2,2), BCO(2,2)
      REAL(DOUBLE), INTENT(OUT) :: KFO(18,18)
      REAL(DOUBLE)              :: ALPHAG, ALPHAM, XI, ETA, WT, N(3), DNP(2,3), NIX(3), NIY(3), GTH(2,18), HTHV(18)
      REAL(DOUBLE)              :: FAC, KTMP(18,18)
      INTEGER(LONG)             :: I

      KFO = ZERO
      ALPHAG = 1.0D-3 * YOUNG * THICK**3 / TWELVE
      ALPHAM = 1.0D-3 * YOUNG * THICK / (TWO*(ONE+POIS))
      DNP = ZERO
      DNP(1,1) = -ONE
      DNP(1,2) =  ONE
      DNP(2,1) = -ONE
      DNP(2,3) =  ONE

      DO I=1,3
         CALL PLATE_GP ( I, XI, ETA, WT )
         N(1) = ONE - XI - ETA
         N(2) = XI
         N(3) = ETA
         NIX = DNP(1,1:3)*COO(1,1) + DNP(2,1:3)*COO(2,1)
         NIY = DNP(1,1:3)*COO(1,2) + DNP(2,1:3)*COO(2,2)
         GTH = ZERO
         HTHV = ZERO
         DO J=1,3
            GTH(1,6*(J-1)+4:6*(J-1)+6) = NIX(J) * NORMSO(J,1:3)
            GTH(2,6*(J-1)+4:6*(J-1)+6) = NIY(J) * NORMSO(J,1:3)
            HTHV(6*(J-1)+4:6*(J-1)+6) = N(J) * NORMSO(J,1:3)
         ENDDO
         FAC = WT * JACO
         KTMP = MATMUL(TRANSPOSE(GTH), GTH)
         KFO = KFO + FAC * (ALPHAG * KTMP + ALPHAM * OUTER18_V ( HTHV ))
      ENDDO
      END SUBROUTINE BUILD_FICTITIOUS_K_MAKNUN

      FUNCTION OUTER18_V ( V ) RESULT (M)
      REAL(DOUBLE), INTENT(IN) :: V(18)
      REAL(DOUBLE)             :: M(18,18)
      INTEGER(LONG)            :: I, J
      DO I=1,18
         DO J=1,18
            M(I,J) = V(I) * V(J)
         ENDDO
      ENDDO
      END FUNCTION OUTER18_V

      SUBROUTINE BUILD_CONSTITUTIVE ( YOUNG, POIS, THICK, DBO, DMO, DSO )
      REAL(DOUBLE), INTENT(IN)  :: YOUNG, POIS, THICK
      REAL(DOUBLE), INTENT(OUT) :: DBO(3,3), DMO(3,3), DSO(2,2)
      REAL(DOUBLE) :: DBF, DMF, GSH

      DBO = ZERO
      DMO = ZERO
      DSO = ZERO

      DBF = YOUNG*THICK**3 / (TWELVE*(ONE-POIS*POIS))
      DBO(1,1) = DBF
      DBO(1,2) = DBF*POIS
      DBO(2,1) = DBO(1,2)
      DBO(2,2) = DBF
      DBO(3,3) = DBF*(ONE-POIS)/TWO

      DMF = YOUNG*THICK / (ONE-POIS*POIS)
      DMO(1,1) = DMF
      DMO(1,2) = DMF*POIS
      DMO(2,1) = DMO(1,2)
      DMO(2,2) = DMF
      DMO(3,3) = DMF*(ONE-POIS)/TWO

      GSH = KAPPA*YOUNG*THICK / (TWO*(ONE+POIS))
      DSO(1,1) = GSH
      DSO(2,2) = GSH
      END SUBROUTINE BUILD_CONSTITUTIVE

      SUBROUTINE BUILD_MEMBRANE_K ( AR, XLOC, YLOC, DMO, KMO )
      REAL(DOUBLE), INTENT(IN)  :: AR, XLOC(3), YLOC(3), DMO(3,3)
      REAL(DOUBLE), INTENT(OUT) :: KMO(6,6)
      REAL(DOUBLE) :: B(3,6), TMP(3,6), BV(3), CV(3), F

      BV = (/ YLOC(2)-YLOC(3), YLOC(3)-YLOC(1), YLOC(1)-YLOC(2) /)
      CV = (/ XLOC(3)-XLOC(2), XLOC(1)-XLOC(3), XLOC(2)-XLOC(1) /)
      F = ONE/(TWO*AR)

      B = ZERO
      B(1,1) = F*BV(1)
      B(1,3) = F*BV(2)
      B(1,5) = F*BV(3)
      B(2,2) = F*CV(1)
      B(2,4) = F*CV(2)
      B(2,6) = F*CV(3)
      B(3,1) = F*CV(1)
      B(3,2) = F*BV(1)
      B(3,3) = F*CV(2)
      B(3,4) = F*BV(2)
      B(3,5) = F*CV(3)
      B(3,6) = F*BV(3)

      TMP = MATMUL(DMO, B)
      KMO = AR*MATMUL(TRANSPOSE(B), TMP)
      END SUBROUTINE BUILD_MEMBRANE_K

      SUBROUTINE BUILD_PLATE_K ( AR, XLOC, YLOC, LOUT, COUT, SOUT, PHIOUT, DBO, DSO, KPO )
      REAL(DOUBLE), INTENT(IN)  :: AR, XLOC(3), YLOC(3), LOUT(3), COUT(3), SOUT(3), PHIOUT(3), DBO(3,3), DSO(2,2)
      REAL(DOUBLE), INTENT(OUT) :: KPO(9,9)
      REAL(DOUBLE) :: XI, ETA, WT, AW(3,9), AN(3,9), BSG(2,3), BB0(3,9), BBD(3,3), BB(3,9), BS(2,9)
      REAL(DOUBLE) :: TMPB(3,9), TMPS(2,9)

      KPO = ZERO
      CALL BUILD_AW_AN ( LOUT, COUT, SOUT, PHIOUT, AW, AN )
      CALL BUILD_BSG ( COUT, SOUT, BSG )
      CALL BUILD_BS ( AN, PHIOUT, BSG, BS )

      CALL PLATE_GP(1, XI, ETA, WT)
      CALL BUILD_BB ( XI, ETA, AR, XLOC, YLOC, COUT, SOUT, AN, BB0, BBD, BB )
      TMPB = MATMUL(DBO, BB)
      TMPS = MATMUL(DSO, BS)
      KPO  = KPO + WT*AR*MATMUL(TRANSPOSE(BB), TMPB) + WT*AR*MATMUL(TRANSPOSE(BS), TMPS)

      CALL PLATE_GP(2, XI, ETA, WT)
      CALL BUILD_BB ( XI, ETA, AR, XLOC, YLOC, COUT, SOUT, AN, BB0, BBD, BB )
      TMPB = MATMUL(DBO, BB)
      KPO  = KPO + WT*AR*MATMUL(TRANSPOSE(BB), TMPB) + WT*AR*MATMUL(TRANSPOSE(BS), TMPS)

      CALL PLATE_GP(3, XI, ETA, WT)
      CALL BUILD_BB ( XI, ETA, AR, XLOC, YLOC, COUT, SOUT, AN, BB0, BBD, BB )
      TMPB = MATMUL(DBO, BB)
      KPO  = KPO + WT*AR*MATMUL(TRANSPOSE(BB), TMPB) + WT*AR*MATMUL(TRANSPOSE(BS), TMPS)
      END SUBROUTINE BUILD_PLATE_K

      SUBROUTINE BUILD_DRILLING_K ( D66, KZO )
      REAL(DOUBLE), INTENT(IN)  :: D66
      REAL(DOUBLE), INTENT(OUT) :: KZO(3,3)
      REAL(DOUBLE) :: KFAC

      KFAC = FOUR*D66 / DSQRT(THREE)
      KZO = ZERO
      KZO(1,1) =  KFAC
      KZO(1,2) = -HALF*KFAC
      KZO(1,3) = -HALF*KFAC
      KZO(2,1) = -HALF*KFAC
      KZO(2,2) =  KFAC
      KZO(2,3) = -HALF*KFAC
      KZO(3,1) = -HALF*KFAC
      KZO(3,2) = -HALF*KFAC
      KZO(3,3) =  KFAC
      END SUBROUTINE BUILD_DRILLING_K

      SUBROUTINE BUILD_AW_AN ( LOUT, COUT, SOUT, PHIOUT, AWO, ANO )
      REAL(DOUBLE), INTENT(IN)  :: LOUT(3), COUT(3), SOUT(3), PHIOUT(3)
      REAL(DOUBLE), INTENT(OUT) :: AWO(3,9), ANO(3,9)
      INTEGER(LONG) :: II

      AWO = ZERO
      AWO(1,1) = -ONE/LOUT(1)
      AWO(1,2) =  HALF*COUT(1)
      AWO(1,3) =  HALF*SOUT(1)
      AWO(1,4) =  ONE/LOUT(1)
      AWO(1,5) =  HALF*COUT(1)
      AWO(1,6) =  HALF*SOUT(1)

      AWO(2,4) = -ONE/LOUT(2)
      AWO(2,5) =  HALF*COUT(2)
      AWO(2,6) =  HALF*SOUT(2)
      AWO(2,7) =  ONE/LOUT(2)
      AWO(2,8) =  HALF*COUT(2)
      AWO(2,9) =  HALF*SOUT(2)

      AWO(3,7) = -ONE/LOUT(3)
      AWO(3,8) =  HALF*COUT(3)
      AWO(3,9) =  HALF*SOUT(3)
      AWO(3,1) =  ONE/LOUT(3)
      AWO(3,2) =  HALF*COUT(3)
      AWO(3,3) =  HALF*SOUT(3)

      ANO = ZERO
      DO II=1,3
         ANO(II,1:9) = AWO(II,1:9) / ((TWO/THREE)*LOUT(II)*(ONE + PHIOUT(II)))
      ENDDO
      END SUBROUTINE BUILD_AW_AN

      SUBROUTINE BUILD_BSG ( COUT, SOUT, BSGO )
      REAL(DOUBLE), INTENT(IN)  :: COUT(3), SOUT(3)
      REAL(DOUBLE), INTENT(OUT) :: BSGO(2,3)
      REAL(DOUBLE) :: A1, A2, A3

      A1 = COUT(1)*SOUT(3) - COUT(3)*SOUT(1)
      A2 = COUT(2)*SOUT(1) - COUT(1)*SOUT(2)
      A3 = COUT(3)*SOUT(2) - COUT(2)*SOUT(3)

      BSGO(1,1) = (SOUT(2)/A2 - SOUT(3)/A1)
      BSGO(1,2) = (SOUT(3)/A3 - SOUT(1)/A2)
      BSGO(1,3) = (SOUT(1)/A1 - SOUT(2)/A3)
      BSGO(2,1) = (COUT(3)/A1 - COUT(2)/A2)
      BSGO(2,2) = (COUT(1)/A2 - COUT(3)/A3)
      BSGO(2,3) = (COUT(2)/A3 - COUT(1)/A1)
      BSGO = (TWO/THREE) * BSGO
      END SUBROUTINE BUILD_BSG

      SUBROUTINE BUILD_BS ( ANO, PHIOUT, BSGO, BSO )
      REAL(DOUBLE), INTENT(IN)  :: ANO(3,9), PHIOUT(3), BSGO(2,3)
      REAL(DOUBLE), INTENT(OUT) :: BSO(2,9)
      REAL(DOUBLE) :: TMP(2,3)
      INTEGER(LONG) :: II

      TMP = BSGO
      DO II=1,3
         TMP(1,II) = TMP(1,II) * PHIOUT(II)
         TMP(2,II) = TMP(2,II) * PHIOUT(II)
      ENDDO
      BSO = MATMUL(TMP, ANO)
      END SUBROUTINE BUILD_BS

      SUBROUTINE BUILD_BB ( XI, ETA, AR, XLOC, YLOC, COUT, SOUT, ANO, BB0O, BBDO, BBO )
      REAL(DOUBLE), INTENT(IN)  :: XI, ETA, AR, XLOC(3), YLOC(3), COUT(3), SOUT(3), ANO(3,9)
      REAL(DOUBLE), INTENT(OUT) :: BB0O(3,9), BBDO(3,3), BBO(3,9)
      REAL(DOUBLE) :: DP(2,3), DPX(3), DPY(3), BV(3), CV(3)
      INTEGER(LONG) :: II

      BV = (/ YLOC(2)-YLOC(3), YLOC(3)-YLOC(1), YLOC(1)-YLOC(2) /)
      CV = (/ XLOC(3)-XLOC(2), XLOC(1)-XLOC(3), XLOC(2)-XLOC(1) /)

      BB0O = ZERO
      DO II=1,3
         BB0O(1,3*II-1) = BV(II)/(TWO*AR)
         BB0O(2,3*II  ) = CV(II)/(TWO*AR)
         BB0O(3,3*II-1) = CV(II)/(TWO*AR)
         BB0O(3,3*II  ) = BV(II)/(TWO*AR)
      ENDDO

      CALL SHAPE_DP ( XI, ETA, DP )
      DPX = (DP(1,:)*(YLOC(3)-YLOC(1)) + DP(2,:)*(YLOC(1)-YLOC(2))) / (TWO*AR)
      DPY = (-DP(1,:)*(XLOC(3)-XLOC(1)) - DP(2,:)*(XLOC(1)-XLOC(2))) / (TWO*AR)

      BBDO = ZERO
      DO II=1,3
         BBDO(1,II) = COUT(II)*DPX(II)
         BBDO(2,II) = SOUT(II)*DPY(II)
         BBDO(3,II) = COUT(II)*DPY(II) + SOUT(II)*DPX(II)
      ENDDO

      BBO = BB0O + MATMUL(BBDO, ANO)
      END SUBROUTINE BUILD_BB

      SUBROUTINE BUILD_CENTER_B ( XI, ETA, AR, XLOC, YLOC, LOUT, COUT, SOUT, PHIOUT, BMO, BBO, BSO )
      REAL(DOUBLE), INTENT(IN)  :: XI, ETA, AR, XLOC(3), YLOC(3), LOUT(3), COUT(3), SOUT(3), PHIOUT(3)
      REAL(DOUBLE), INTENT(OUT) :: BMO(3,6), BBO(3,9), BSO(2,9)
      REAL(DOUBLE) :: ANO(3,9), AWO(3,9), BSGO(2,3), BB0O(3,9), BBDO(3,3), BV(3), CV(3), F

      CALL BUILD_AW_AN ( LOUT, COUT, SOUT, PHIOUT, AWO, ANO )
      CALL BUILD_BSG ( COUT, SOUT, BSGO )
      CALL BUILD_BS ( ANO, PHIOUT, BSGO, BSO )
      CALL BUILD_BB ( XI, ETA, AR, XLOC, YLOC, COUT, SOUT, ANO, BB0O, BBDO, BBO )

      BV = (/ YLOC(2)-YLOC(3), YLOC(3)-YLOC(1), YLOC(1)-YLOC(2) /)
      CV = (/ XLOC(3)-XLOC(2), XLOC(1)-XLOC(3), XLOC(2)-XLOC(1) /)
      F = ONE/(TWO*AR)
      BMO = ZERO
      BMO(1,1) = F*BV(1)
      BMO(1,3) = F*BV(2)
      BMO(1,5) = F*BV(3)
      BMO(2,2) = F*CV(1)
      BMO(2,4) = F*CV(2)
      BMO(2,6) = F*CV(3)
      BMO(3,1) = F*CV(1)
      BMO(3,2) = F*BV(1)
      BMO(3,3) = F*CV(2)
      BMO(3,4) = F*BV(2)
      BMO(3,5) = F*CV(3)
      BMO(3,6) = F*BV(3)
      END SUBROUTINE BUILD_CENTER_B

      SUBROUTINE SHAPE_DP ( XI, ETA, DP )
      REAL(DOUBLE), INTENT(IN)  :: XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: DP(2,3)
      REAL(DOUBLE) :: LAM
      LAM = ONE - XI - ETA
      DP(1,1) =  FOUR*(LAM - XI)
      DP(1,2) =  FOUR*ETA
      DP(1,3) = -FOUR*ETA
      DP(2,1) = -FOUR*XI
      DP(2,2) =  FOUR*XI
      DP(2,3) =  FOUR*(LAM - ETA)
      END SUBROUTINE SHAPE_DP

      SUBROUTINE PLATE_GP ( IGP, XI, ETA, WT )
      INTEGER(LONG), INTENT(IN) :: IGP
      REAL(DOUBLE), INTENT(OUT) :: XI, ETA, WT
      IF (IGP == 1) THEN
         XI = HALF
         ETA = ZERO
      ELSE IF (IGP == 2) THEN
         XI = ZERO
         ETA = HALF
      ELSE
         XI = HALF
         ETA = HALF
      ENDIF
      WT = ONE/SIX
      END SUBROUTINE PLATE_GP

      SUBROUTINE BUILD_T18 ( TE3, TT )
      REAL(DOUBLE), INTENT(IN)  :: TE3(3,3)
      REAL(DOUBLE), INTENT(OUT) :: TT(18,18)
      INTEGER(LONG) :: II, JJ, NN, R0, C0

      TT = ZERO
      DO NN=0,2
         R0 = 6*NN
         C0 = 6*NN
         DO II=1,3
            DO JJ=1,3
               TT(R0+II    ,C0+JJ    ) = TE3(II,JJ)
               TT(R0+II+3  ,C0+JJ+3  ) = TE3(II,JJ)
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE BUILD_T18

      SUBROUTINE DEBUG_PRINT_MATRIX ( TITLE, MAT )
      CHARACTER(LEN=*), INTENT(IN) :: TITLE
      REAL(DOUBLE), INTENT(IN)     :: MAT(:,:)
      INTEGER(LONG) :: II
      WRITE(F06,'(A)') TRIM(TITLE)
      DO II=1,SIZE(MAT,1)
         WRITE(F06,'(1X,100(1X,ES14.6))') MAT(II,:)
      ENDDO
      WRITE(F06,*)
      END SUBROUTINE DEBUG_PRINT_MATRIX

      END SUBROUTINE CTRIAR_DKMT18
