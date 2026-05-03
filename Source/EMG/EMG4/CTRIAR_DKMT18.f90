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

! --- CQUAD4R_CTRIAR_add begin --- !
! DKMT18 triangular shell based on the 2026-05-03 Python reference:
!   codex_mod\dkmt18_ctriar\dkmt18 files 3-5-2026\dkmt18_element.py
!
! This phase-1 MYSTRAN port keeps the earlier fast integration strategy:
!   - user-facing card name: CTRIAR
!   - internal ETYPE remains TRIA3
!   - EDAT(T3 thickness key) = -18 selects this kernel
!
! Current scope:
!   - isotropic flat-shell stiffness
!   - center-point shell recovery matrices
!   - translational lumped mass
! --- CQUAD4R_CTRIAR_add end --- !

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, THREE, FOUR, SIX, TWELVE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE MODEL_STUF, ONLY            :  EID, ELGP, KE, ME, BE1, BE2, BE3, EPROP, MASS_PER_UNIT_AREA,                           &
                                         NUM_EMG_FATAL_ERRS, SHELL_A, TE, XEB

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

      INTEGER(LONG)                   :: I, J, K
      INTEGER(LONG)                   :: MAP_MEM(6), MAP_PLT(9), MAP_TZ(3)
      REAL(DOUBLE)                    :: H, E, NU, AREA, XI0, ETA0
      REAL(DOUBLE)                    :: XL(3), YL(3), LVAL(3), CVAL(3), SVAL(3), PHI(3), A_DIAG(3)
      REAL(DOUBLE)                    :: DB(3,3), DM(3,3), DS(2,2)
      REAL(DOUBLE)                    :: KM(6,6), KP(9,9), KZ(3,3), K15(15,15), K18(18,18)
      REAL(DOUBLE)                    :: T18(18,18), T18T(18,18)
      REAL(DOUBLE)                    :: BMC(3,6), BBC(3,9), BSC(2,9)
      REAL(DOUBLE)                    :: BEM(3,18), BEB(3,18), BES(2,18), BEML(3,18), BEBL(3,18), BESL(2,18)
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

      CALL BUILD_LOCAL_XY ( XL, YL, AREA )
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

      CALL BUILD_T18 ( TE, T18 )
      T18T = TRANSPOSE(T18)

      CALL BUILD_SIDE_GEOM ( XL, YL, LVAL, CVAL, SVAL )
      CALL BUILD_PHI ( H, NU, LVAL, PHI )
      CALL BUILD_CONSTITUTIVE ( E, NU, H, DB, DM, DS )
      CALL BUILD_MEMBRANE_K ( AREA, XL, YL, DM, KM )
      CALL BUILD_PLATE_K ( AREA, XL, YL, LVAL, CVAL, SVAL, PHI, DB, DS, KP )
      CALL BUILD_DRILLING_K ( DB(3,3), KZ )

      K15 = ZERO
      K15(1:6,1:6)     = KM
      K15(7:15,7:15)   = KP

      MAP_MEM = (/ 1, 2, 7, 8, 13, 14 /)
      MAP_PLT = (/ 3, 4, 5, 9, 10, 11, 15, 16, 17 /)
      MAP_TZ  = (/ 6, 12, 18 /)

      K18 = ZERO
      DO I=1,6
         DO J=1,6
            K18(MAP_MEM(I), MAP_MEM(J)) = KM(I,J)
         ENDDO
      ENDDO
      DO I=1,9
         DO J=1,9
            K18(MAP_PLT(I), MAP_PLT(J)) = KP(I,J)
         ENDDO
      ENDDO
      DO I=1,3
         DO J=1,3
            K18(MAP_TZ(I), MAP_TZ(J)) = K18(MAP_TZ(I), MAP_TZ(J)) + KZ(I,J)
         ENDDO
      ENDDO

      IF (OPT(4) == 'Y') THEN
         KE(1:18,1:18) = K18
      ENDIF

      IF (OPT(3) == 'Y') THEN
         XI0  = ONE/THREE
         ETA0 = ONE/THREE
         CALL BUILD_CENTER_B ( XI0, ETA0, AREA, XL, YL, LVAL, CVAL, SVAL, PHI, BMC, BBC, BSC )

         BEM = ZERO
         BEB = ZERO
         BES = ZERO
         DO I=1,3
            BEM(:,MAP_MEM(2*I-1)) = BMC(:,2*I-1)
            BEM(:,MAP_MEM(2*I  )) = BMC(:,2*I  )
            BEB(:,MAP_PLT(3*I-2)) = BBC(:,3*I-2)
            BEB(:,MAP_PLT(3*I-1)) = BBC(:,3*I-1)
            BEB(:,MAP_PLT(3*I  )) = BBC(:,3*I  )
            BES(:,MAP_PLT(3*I-2)) = BSC(:,3*I-2)
            BES(:,MAP_PLT(3*I-1)) = BSC(:,3*I-1)
            BES(:,MAP_PLT(3*I  )) = BSC(:,3*I  )
         ENDDO

         BEML = BEM
         BEBL = BEB
         BESL = BES

         BE1(1:3,1:18,1) = BEML
         BE2(1:3,1:18,1) = BEBL
         BE3(1:2,1:18,1) = BESL
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

      IF (DEBUG(233) > 0) THEN
         CALL DEBUG_PRINT_MATRIX('CTRIAR K18 BASIC', K18)
         CALL DEBUG_PRINT_MATRIX('CTRIAR KE LOCAL', KE(1:18,1:18))
         CALL DEBUG_PRINT_MATRIX('CTRIAR KMEM', KM)
         CALL DEBUG_PRINT_MATRIX('CTRIAR KPLATE', KP)
         CALL DEBUG_PRINT_MATRIX('CTRIAR KDRILL', KZ)
      ENDIF

      RETURN

 9001 FORMAT(' *ERROR: ',A,' expects ELGP=3 for element ',I8,' but got ',I8)
 9002 FORMAT(' *ERROR: ',A,' got nonpositive thickness for element ',I8,': ',1ES14.6)
 9003 FORMAT(' *ERROR: ',A,' got degenerate area for element ',I8,': ',1ES14.6)
 9004 FORMAT(' *ERROR: ',A,' could not derive isotropic shell constants from SHELL_A for element ',I8)

      CONTAINS

      SUBROUTINE BUILD_LOCAL_XY ( XLOC, YLOC, AR )
      REAL(DOUBLE), INTENT(OUT) :: XLOC(3), YLOC(3), AR
      REAL(DOUBLE) :: DX(3)
      INTEGER(LONG) :: II
      XLOC = ZERO
      YLOC = ZERO
      DO II=1,3
         DX = XEB(II,1:3) - XEB(1,1:3)
         XLOC(II) = DOT_PRODUCT(TE(1,1:3), DX)
         YLOC(II) = DOT_PRODUCT(TE(2,1:3), DX)
      ENDDO
      AR = HALF*DABS((XLOC(2)-XLOC(1))*(YLOC(3)-YLOC(1)) - (XLOC(3)-XLOC(1))*(YLOC(2)-YLOC(1)))
      END SUBROUTINE BUILD_LOCAL_XY

      SUBROUTINE BUILD_SIDE_GEOM ( XLOC, YLOC, LOUT, COUT, SOUT )
      REAL(DOUBLE), INTENT(IN)  :: XLOC(3), YLOC(3)
      REAL(DOUBLE), INTENT(OUT) :: LOUT(3), COUT(3), SOUT(3)
      REAL(DOUBLE) :: DX, DY

      DX = XLOC(2) - XLOC(1)
      DY = YLOC(2) - YLOC(1)
      LOUT(1) = DSQRT(DX*DX + DY*DY)
      COUT(1) = DX/LOUT(1)
      SOUT(1) = DY/LOUT(1)

      DX = XLOC(3) - XLOC(2)
      DY = YLOC(3) - YLOC(2)
      LOUT(2) = DSQRT(DX*DX + DY*DY)
      COUT(2) = DX/LOUT(2)
      SOUT(2) = DY/LOUT(2)

      DX = XLOC(1) - XLOC(3)
      DY = YLOC(1) - YLOC(3)
      LOUT(3) = DSQRT(DX*DX + DY*DY)
      COUT(3) = DX/LOUT(3)
      SOUT(3) = DY/LOUT(3)
      END SUBROUTINE BUILD_SIDE_GEOM

      SUBROUTINE BUILD_PHI ( THICK, POIS, LOUT, PHIOUT )
      REAL(DOUBLE), INTENT(IN)  :: THICK, POIS, LOUT(3)
      REAL(DOUBLE), INTENT(OUT) :: PHIOUT(3)
      INTEGER(LONG) :: II
      DO II=1,3
         PHIOUT(II) = TWO/(KAPPA*(ONE-POIS)) * (THICK/LOUT(II))**2
      ENDDO
      END SUBROUTINE BUILD_PHI

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
      WT = ONE/THREE
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
