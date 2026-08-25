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

      SUBROUTINE CQUAD4_DSQK_RHR ( OPT, INT_ELEM_ID )

! Standalone DSQK CQUAD4 branch ported from:
!   D:\18a\bending_only\Shell\gemini2\shit\claudeBuckling\got\DSQK_ShellElement_RHR.py
!   D:\18a\bending_only\Shell\gemini2\shit\claudeBuckling\got\DSQK_ShellElement_RHR_6DOF.py
!
! This file intentionally keeps its own helpers local so the DSQK family can be
! swapped in/out without entangling other CQUAD4 branches.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_GAUSS, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, FOUR
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  COUPMASS
      USE MODEL_STUF, ONLY            :  ALPVEC, BGRID, DT, EID, ELGP, EPROP, FCONV, GRID_SNORM, KE, KED, ME, BE1, BE2, BE3,    &
                                         MASS_PER_UNIT_AREA, NUM_EMG_FATAL_ERRS, PCOMP_PROPS, PPE, PRESS, PTE, SHELL_A,         &
                                         SHELL_D, SHELL_T, STRESS, TE, TREF, XEB

      USE ELMDIS_Interface
      USE ELEM_STRE_STRN_ARRAYS_Interface
      USE ORDER_GAUSS_Interface
      USE OUTA_HERE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CQUAD4_DSQK_RHR'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      INTEGER(LONG), PARAMETER        :: NNODE = 4
      INTEGER(LONG), PARAMETER        :: NDOF6 = 24
      INTEGER(LONG), PARAMETER        :: NDOF5 = 20
      REAL(DOUBLE), PARAMETER         :: SHEAR_KAPPA = 5.0D0/6.0D0
      REAL(DOUBLE), PARAMETER         :: DRILL_SCALE = 3.0D-2

      INTEGER(LONG)                   :: I, J, K, GP, JSUB, IA, IB, RR
      REAL(DOUBLE)                    :: XYZ(4,3), NORMS(4,3), XLOCAL(4,2), T24(24,24), T24T(24,24), T56(20,24)
      REAL(DOUBLE)                    :: T1(3), T2(3), T3(3), TV1(3), TV2(3), NVEC(3), SURF(3)
      REAL(DOUBLE)                    :: SS(MAX_ORDER_GAUSS), HH(MAX_ORDER_GAUSS), XI, ETA, WT
      REAL(DOUBLE)                    :: JDET2, JDET3, AVERAGE_THICK, G12, TBAR
      REAL(DOUBLE)                    :: SIDE_C(4), SIDE_S(4), SIDE_L(4)
      REAL(DOUBLE)                    :: BM5(3,20), BBU5(3,20), BB5(3,20), BS5(2,20), BBD(3,4), AN4(4,20)
      REAL(DOUBLE)                    :: BM24(3,24), BB24(3,24), BS24(2,24), BML(3,24), BBL(3,24), BSL(2,24)
      REAL(DOUBLE)                    :: KM5(20,20), KBU5(20,20), KBD4(4,4), KB5(20,20), KS5(20,20), KTOT5(20,20)
      REAL(DOUBLE)                    :: KBASIC(24,24), KLOCAL(24,24), KDRILL(24,24), KG20(20,20), KGBASIC(24,24), KGLOCAL(24,24)
      REAL(DOUBLE)                    :: M5(20,20), MBASIC(24,24), MLOCAL(24,24), UNIT_PPE_B(24), UNIT_PPE_L(24), UNIT_PTE_B(24)
      REAL(DOUBLE)                    :: CTE3(3), NTH(3), SIG0(2,2), GUX(2,20), GVY(2,20), GWW(2,20)
      REAL(DOUBLE)                    :: GU24(2,24), GV24(2,24), GW24(2,24), S2(2,2)
      REAL(DOUBLE)                    :: GBE1(3,24,4), GBE2(3,24,4), GBE3(2,24,4), REC_XI(5), REC_ETA(5)

! **********************************************************************************************************************************

      IF (ELGP /= 4) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9001) SUBR_NAME, EID, ELGP
         WRITE(F06,9001) SUBR_NAME, EID, ELGP
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IF (PCOMP_PROPS == 'Y') THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9002) TRIM(SUBR_NAME), EID
         WRITE(F06,9002) TRIM(SUBR_NAME), EID
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL LOAD_BASIC_COORDS ( XYZ )
      CALL CALC_NODAL_NORMALS ( XYZ, NORMS )
      CALL BUILD_DSQK_BASIS ( XYZ, NORMS, T1, T2, T3, XLOCAL, SIDE_C, SIDE_S, SIDE_L )
      CALL BUILD_T56 ( T1, T2, T3, T56 )
      CALL BUILD_T24 ( TE, T24 )
      T24T = TRANSPOSE(T24)

      AVERAGE_THICK = MAX(EPROP(1), 1.0D-12)
      G12 = SHELL_T(1,1) / MAX(SHEAR_KAPPA*AVERAGE_THICK, 1.0D-30)

      KM5 = ZERO
      KBU5 = ZERO
      KBD4 = ZERO
      KB5 = ZERO
      KS5 = ZERO
      KTOT5 = ZERO
      KBASIC = ZERO
      KLOCAL = ZERO
      KDRILL = ZERO
      KG20 = ZERO
      KGBASIC = ZERO
      KGLOCAL = ZERO
      M5 = ZERO
      MBASIC = ZERO
      MLOCAL = ZERO
      UNIT_PPE_B = ZERO
      UNIT_PPE_L = ZERO
      UNIT_PTE_B = ZERO
      GBE1 = ZERO
      GBE2 = ZERO
      GBE3 = ZERO

      IF ((OPT(1) == 'Y') .OR. (OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(4) == 'Y') .OR. (OPT(5) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         CALL ORDER_GAUSS ( 2, SS, HH )
      ENDIF

      AN4 = AN_DSQK ( SIDE_C, SIDE_S, SIDE_L )

      IF ((OPT(3) == 'Y') .OR. (OPT(4) == 'Y')) THEN
         DO I=1,2
            DO J=1,2
               XI  = SS(I)
               ETA = SS(J)
               WT  = HH(I)*HH(J)

               JDET2 = DETJ_DSQK ( XLOCAL, XI, ETA )
               BM5 = BM_DSQK ( XLOCAL, XI, ETA )
               BBU5 = BBU_DSQK ( XLOCAL, XI, ETA )
               BBD = BBDELTA_DSQK ( XLOCAL, SIDE_C, SIDE_S, XI, ETA )
               BB5 = BBU5 + MATMUL(BBD, AN4)
               BS5 = BS_DSQK ( XLOCAL, XI, ETA )

               BM24 = MATMUL(BM5, T56)
               BB24 = MATMUL(BB5, T56)
               BS24 = MATMUL(BS5, T56)
               BML  = MATMUL(BM24, T24T)
               BBL  = MATMUL(BB24, T24T)
               BSL  = MATMUL(BS24, T24T)

               IF ((DEBUG(233) > 0) .AND. (EID <= 8) .AND. (I == 1) .AND. (J == 1)) THEN
                  WRITE(F06,'(A,I8)') 'CQUAD4_DSQK BASIS EID=', EID
                  WRITE(F06,'(A,3(1X,ES15.7))') '  T1', T1
                  WRITE(F06,'(A,3(1X,ES15.7))') '  T2', T2
                  WRITE(F06,'(A,3(1X,ES15.7))') '  T3', T3
                  WRITE(F06,'(A,8(1X,ES15.7))') '  XLOCAL', XLOCAL(1,1), XLOCAL(1,2), XLOCAL(2,1), XLOCAL(2,2),                       &
                                                 XLOCAL(3,1), XLOCAL(3,2), XLOCAL(4,1), XLOCAL(4,2)
                  CALL DEBUG_PRINT_MATRIX('CQUAD4_DSQK TE', TE)
                  CALL DEBUG_PRINT_MATRIX('CQUAD4_DSQK GP11 AN4', AN4)
                  CALL DEBUG_PRINT_MATRIX('CQUAD4_DSQK GP11 T56', T56)
                  CALL DEBUG_PRINT_MATRIX('CQUAD4_DSQK GP11 BM5', BM5)
                  CALL DEBUG_PRINT_MATRIX('CQUAD4_DSQK GP11 BBU5', BBU5)
                  CALL DEBUG_PRINT_MATRIX('CQUAD4_DSQK GP11 BBD', BBD)
                  CALL DEBUG_PRINT_MATRIX('CQUAD4_DSQK GP11 BB5', BB5)
                  CALL DEBUG_PRINT_MATRIX('CQUAD4_DSQK GP11 BS5', BS5)
                  CALL DEBUG_PRINT_MATRIX('CQUAD4_DSQK GP11 BM24', BM24)
                  CALL DEBUG_PRINT_MATRIX('CQUAD4_DSQK GP11 BB24', BB24)
               ENDIF

               KM5  = KM5  + WT*JDET2*MATMUL(TRANSPOSE(BM5 ), MATMUL(SHELL_A, BM5 ))
               KBU5 = KBU5 + WT*JDET2*MATMUL(TRANSPOSE(BBU5), MATMUL(SHELL_D, BBU5))
               KBD4 = KBD4 + WT*JDET2*MATMUL(TRANSPOSE(BBD ), MATMUL(SHELL_D, BBD ))
               KS5  = KS5  + WT*JDET2*MATMUL(TRANSPOSE(BS5 ), MATMUL(SHELL_T, BS5 ))

               GP = GP_INDEX(I, J)
               GBE1(1:3,1:24,GP) = BML
               GBE2(1:3,1:24,GP) = BBL
               GBE3(1:2,1:24,GP) = BSL
            ENDDO
         ENDDO

         KB5 = KBU5 + MATMUL(TRANSPOSE(AN4), MATMUL(KBD4, AN4))
         KTOT5 = KM5 + KB5 + KS5
         KBASIC = MATMUL(TRANSPOSE(T56), MATMUL(KTOT5, T56))
         KDRILL = ZERO
! DEBUG(250) isolates the DSQK reduced-kernel port from the added drilling
! stabilizer. This helps separate a bad 5DOF->24DOF lift from drill-only
! contamination in patch tests and simple bending checks.
         IF (DEBUG(250) <= 0) THEN
            CALL BUILD_DSQK_DRILL ( T3, G12, AVERAGE_THICK, XLOCAL, KDRILL )
            KBASIC = KBASIC + KDRILL
         ENDIF
         KLOCAL = MATMUL(T24, MATMUL(KBASIC, T24T))

         IF (OPT(4) == 'Y') THEN
            KE(1:24,1:24) = KLOCAL
         ENDIF

         IF (OPT(3) == 'Y') THEN
            BE1(:,:,1) = (GBE1(:,:,1) + GBE1(:,:,2) + GBE1(:,:,3) + GBE1(:,:,4)) / FOUR
            BE2(:,:,1) = (GBE2(:,:,1) + GBE2(:,:,2) + GBE2(:,:,3) + GBE2(:,:,4)) / FOUR
            BE3(1:2,:,1) = (GBE3(1:2,:,1) + GBE3(1:2,:,2) + GBE3(1:2,:,3) + GBE3(1:2,:,4)) / FOUR

            BE1(:,:,2) = GBE1(:,:,4)
            BE1(:,:,3) = GBE1(:,:,3)
            BE1(:,:,4) = GBE1(:,:,2)
            BE1(:,:,5) = GBE1(:,:,1)

            BE2(:,:,2) = GBE2(:,:,4)
            BE2(:,:,3) = GBE2(:,:,3)
            BE2(:,:,4) = GBE2(:,:,2)
            BE2(:,:,5) = GBE2(:,:,1)
            DO GP=2,5
               BE2(:,:,GP) = BE2(:,:,1)
            ENDDO

            BE3(1:2,:,2) = GBE3(1:2,:,4)
            BE3(1:2,:,3) = GBE3(1:2,:,3)
            BE3(1:2,:,4) = GBE3(1:2,:,2)
            BE3(1:2,:,5) = GBE3(1:2,:,1)

         ENDIF
      ENDIF

      IF (OPT(2) == 'Y') THEN
         CTE3(1) = ALPVEC(1,1)
         CTE3(2) = ALPVEC(2,1)
         CTE3(3) = ALPVEC(4,1)
         NTH = MATMUL(SHELL_A, CTE3)

         DO I=1,2
            DO J=1,2
               XI  = SS(I)
               ETA = SS(J)
               WT  = HH(I)*HH(J)
               JDET2 = DETJ_DSQK ( XLOCAL, XI, ETA )
               BM24 = MATMUL( BM_DSQK(XLOCAL, XI, ETA), T56 )
               UNIT_PTE_B = UNIT_PTE_B + WT*JDET2*MATMUL( TRANSPOSE(BM24), NTH )
            ENDDO
         ENDDO

         DO JSUB=1,SIZE(PTE,2)
            TBAR = (DT(1,JSUB) + DT(2,JSUB) + DT(3,JSUB) + DT(4,JSUB))/FOUR - TREF(1)
            PTE(1:24,JSUB) = MATMUL(T24, UNIT_PTE_B) * TBAR
         ENDDO
      ENDIF

      IF (OPT(1) == 'Y') THEN
         DO I=1,2
            DO J=1,2
               XI  = SS(I)
               ETA = SS(J)
               WT  = HH(I)*HH(J)
               JDET2 = DETJ_DSQK ( XLOCAL, XI, ETA )
               CALL MASS_DSQK_AT ( XI, ETA, WT*JDET2, MASS_PER_UNIT_AREA, AVERAGE_THICK, M5 )
            ENDDO
         ENDDO
         MBASIC = MATMUL( TRANSPOSE(T56), MATMUL(M5, T56) )
         MLOCAL = MATMUL(T24, MATMUL(MBASIC, T24T))
         ME(1:24,1:24) = MLOCAL
      ENDIF

      IF (OPT(5) == 'Y') THEN
         DO I=1,2
            DO J=1,2
               XI  = SS(I)
               ETA = SS(J)
               WT  = HH(I)*HH(J)
               CALL GEOM_SURFACE_AT ( XYZ, XI, ETA, SURF, JDET3, NVEC )
               CALL PRESSURE_VECTOR_AT ( XI, ETA, WT*JDET3, NVEC, UNIT_PPE_B )
            ENDDO
         ENDDO
         UNIT_PPE_L = MATMUL(T24, UNIT_PPE_B)
         DO JSUB=1,SIZE(PPE,2)
            PPE(1:24,JSUB) = UNIT_PPE_L(1:24) * PRESS(3,JSUB)
         ENDDO
      ENDIF

      IF (OPT(6) == 'Y') THEN
         CALL ELMDIS

         S2 = ZERO
         DO I=1,2
            DO J=1,2
               XI  = SS(I)
               ETA = SS(J)
               WT  = HH(I)*HH(J)
               JDET2 = DETJ_DSQK ( XLOCAL, XI, ETA )

               BM24 = MATMUL( BM_DSQK(XLOCAL, XI, ETA), T56 )
               BML = MATMUL(BM24, T24T)
               BE1(1:3,1:24,1) = BML
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )

               SIG0(1,1) = FCONV(1)*STRESS(1)
               SIG0(2,2) = FCONV(1)*STRESS(2)
               SIG0(1,2) = FCONV(1)*STRESS(3)
               SIG0(2,1) = SIG0(1,2)
               S2 = SIG0

               CALL BUILD_G_DSQK ( XLOCAL, XI, ETA, GUX, GVY, GWW )
               GU24 = MATMUL(GUX, T56)
               GV24 = MATMUL(GVY, T56)
               GW24 = MATMUL(GWW, T56)
               KGBASIC = KGBASIC + WT*JDET2*( MATMUL( TRANSPOSE(GU24), MATMUL(S2, GU24) ) +                                      &
                                             MATMUL( TRANSPOSE(GV24), MATMUL(S2, GV24) ) +                                      &
                                             MATMUL( TRANSPOSE(GW24), MATMUL(S2, GW24) ) )
            ENDDO
         ENDDO

         KGLOCAL = MATMUL(T24, MATMUL(KGBASIC, T24T))
         KED(1:24,1:24) = KGLOCAL
      ENDIF

      IF (DEBUG(233) > 0) THEN
         WRITE(F06,'(A,I8,A,ES15.7)') 'CQUAD4_DSQK EID=', EID, ' KM5_NORM=', DSQRT(SUM(KM5*KM5))
         WRITE(F06,'(A,I8,A,ES15.7)') 'CQUAD4_DSQK EID=', EID, ' KB5_NORM=', DSQRT(SUM(KB5*KB5))
         WRITE(F06,'(A,I8,A,ES15.7)') 'CQUAD4_DSQK EID=', EID, ' KS5_NORM=', DSQRT(SUM(KS5*KS5))
         WRITE(F06,'(A,I8,A,ES15.7)') 'CQUAD4_DSQK EID=', EID, ' KTOT5_NORM=', DSQRT(SUM(KTOT5*KTOT5))
         WRITE(F06,'(A,I8,A,ES15.7)') 'CQUAD4_DSQK EID=', EID, ' KDRILL_NORM=', DSQRT(SUM(KDRILL*KDRILL))
         CALL DEBUG_PRINT_MATRIX('CQUAD4_DSQK_RHR KBASIC', KBASIC)
         CALL DEBUG_PRINT_MATRIX('CQUAD4_DSQK_RHR KLOCAL', KLOCAL)
         CALL DEBUG_PRINT_MATRIX('CQUAD4_DSQK_RHR KDRILL', KDRILL)
      ENDIF

      RETURN

 9001 FORMAT(' *ERROR: ',A,' expects ELGP=4 for element ',I8,' but got ',I8)
 9002 FORMAT(' *ERROR: ',A,' does not yet support PCOMP on element ',I8)

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

      SUBROUTINE CALC_NODAL_NORMALS ( XYZN, NORMSOUT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(4,3)
      REAL(DOUBLE), INTENT(OUT) :: NORMSOUT(4,3)
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
      NORMSOUT(1,:) = N / NM

      A = XYZN(3,:) - XYZN(2,:)
      B = XYZN(1,:) - XYZN(2,:)
      CALL CROSS3(A, B, N)
      NM = VNORM(N)
      IF (NM <= 1.0D-15) THEN
         N = (/ZERO, ZERO, ONE/)
         NM = ONE
      ENDIF
      NORMSOUT(2,:) = N / NM

      A = XYZN(4,:) - XYZN(3,:)
      B = XYZN(2,:) - XYZN(3,:)
      CALL CROSS3(A, B, N)
      NM = VNORM(N)
      IF (NM <= 1.0D-15) THEN
         N = (/ZERO, ZERO, ONE/)
         NM = ONE
      ENDIF
      NORMSOUT(3,:) = N / NM

      A = XYZN(1,:) - XYZN(4,:)
      B = XYZN(3,:) - XYZN(4,:)
      CALL CROSS3(A, B, N)
      NM = VNORM(N)
      IF (NM <= 1.0D-15) THEN
         N = (/ZERO, ZERO, ONE/)
         NM = ONE
      ENDIF
      NORMSOUT(4,:) = N / NM

      IF (ALLOCATED(GRID_SNORM)) THEN
         DO II=1,4
            IF ((BGRID(II) > 0) .AND. (BGRID(II) <= SIZE(GRID_SNORM,1))) THEN
               SN = GRID_SNORM(BGRID(II),:)
               NM = VNORM(SN)
               IF (NM > 1.0D-15) THEN
                  SN = SN / NM
                  SDOT = DOT_PRODUCT(SN, NORMSOUT(II,:))
                  IF (SDOT < 1.0D-2) THEN
                     NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,'(A,A,A,I8,A,I2,A,ES14.6)') ' *ERROR: ', TRIM(SUBR_NAME), ' EID=', EID,                            &
                        ' SNORM AT NODE ', II, ' IS TOO FAR FROM DSQK MIDSURFACE NORMAL. DOT=', SDOT
                     WRITE(F06,'(A,A,A,I8,A,I2,A,ES14.6)') ' *ERROR: ', TRIM(SUBR_NAME), ' EID=', EID,                            &
                        ' SNORM AT NODE ', II, ' IS TOO FAR FROM DSQK MIDSURFACE NORMAL. DOT=', SDOT
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF
                  NORMSOUT(II,:) = SN
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      END SUBROUTINE CALC_NODAL_NORMALS

      SUBROUTINE BUILD_DSQK_BASIS ( XYZN, NORMSIN, E1, E2, E3, XLOC, CK, SK, LK )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(4,3), NORMSIN(4,3)
      REAL(DOUBLE), INTENT(OUT) :: E1(3), E2(3), E3(3), XLOC(4,2), CK(4), SK(4), LK(4)
      INTEGER(LONG) :: II, JJ
      REAL(DOUBLE) :: D1(3), D2(3), U1(3), U2(3), E1F(3), E2F(3), E3F(3), TMP(3), AVG(3), NM, DX, DY
      INTEGER(LONG), PARAMETER :: EDGE_I(4) = (/1,2,3,4/)
      INTEGER(LONG), PARAMETER :: EDGE_J(4) = (/2,3,4,1/)

      D1 = XYZN(3,:) - XYZN(1,:)
      D2 = XYZN(2,:) - XYZN(4,:)
      U1 = SAFE_UNIT(D1, (/ONE, ZERO, ZERO/))
      U2 = SAFE_UNIT(D2, (/ZERO, ONE, ZERO/))
      E1F = SAFE_UNIT(U1 + U2, U1)
      E2F = SAFE_UNIT(U1 - U2, U2)
      CALL CROSS3(E1F, E2F, E3F)
      E3F = SAFE_UNIT(E3F, (/ZERO, ZERO, ONE/))

      AVG = ZERO
      DO II=1,4
         AVG = AVG + NORMSIN(II,:)
      ENDDO
      NM = VNORM(AVG)
      IF (NM > 1.0D-12) THEN
         E3 = AVG / NM
      ELSE
         E3 = E3F
      ENDIF
      IF (DOT_PRODUCT(E3, E3F) < ZERO) E3 = -E3

      TMP = E1F - DOT_PRODUCT(E1F, E3)*E3
      NM = VNORM(TMP)
      IF (NM < 1.0D-8) THEN
         TMP = E2F - DOT_PRODUCT(E2F, E3)*E3
         NM = VNORM(TMP)
      ENDIF
      E1 = SAFE_UNIT(TMP, E1F)
      CALL CROSS3(E3, E1, E2)
      E2 = SAFE_UNIT(E2, E2F)

      DO II=1,4
         XLOC(II,1) = DOT_PRODUCT(XYZN(II,:), E1)
         XLOC(II,2) = DOT_PRODUCT(XYZN(II,:), E2)
      ENDDO

      DO II=1,4
         JJ = EDGE_J(II)
         DX = XLOC(JJ,1) - XLOC(EDGE_I(II),1)
         DY = XLOC(JJ,2) - XLOC(EDGE_I(II),2)
         LK(II) = DSQRT(MAX(DX*DX + DY*DY, 1.0D-30))
         CK(II) = DX / LK(II)
         SK(II) = DY / LK(II)
      ENDDO
      END SUBROUTINE BUILD_DSQK_BASIS

      SUBROUTINE BUILD_T56 ( E1, E2, E3, TOUT )
      REAL(DOUBLE), INTENT(IN)  :: E1(3), E2(3), E3(3)
      REAL(DOUBLE), INTENT(OUT) :: TOUT(20,24)
      INTEGER(LONG) :: II, BASE5, BASE6
      TOUT = ZERO
      DO II=1,4
         BASE5 = (II-1)*5
         BASE6 = (II-1)*6
         TOUT(BASE5+1, BASE6+1:BASE6+3) = E1
         TOUT(BASE5+2, BASE6+1:BASE6+3) = E2
         TOUT(BASE5+3, BASE6+1:BASE6+3) = E3
         TOUT(BASE5+4, BASE6+4:BASE6+6) = E1
         TOUT(BASE5+5, BASE6+4:BASE6+6) = E2
      ENDDO
      END SUBROUTINE BUILD_T56

      SUBROUTINE SHAPE_N ( XI, ETA, NVAL )
      REAL(DOUBLE), INTENT(IN)  :: XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: NVAL(4)
      NVAL(1) = 0.25D0*(ONE - XI)*(ONE - ETA)
      NVAL(2) = 0.25D0*(ONE + XI)*(ONE - ETA)
      NVAL(3) = 0.25D0*(ONE + XI)*(ONE + ETA)
      NVAL(4) = 0.25D0*(ONE - XI)*(ONE + ETA)
      END SUBROUTINE SHAPE_N

      SUBROUTINE SHAPE_DN ( XI, ETA, DN )
      REAL(DOUBLE), INTENT(IN)  :: XI, ETA
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

      FUNCTION DETJ_DSQK ( XLOC, XI, ETA ) RESULT(DETJ)
      REAL(DOUBLE), INTENT(IN) :: XLOC(4,2), XI, ETA
      REAL(DOUBLE) :: DETJ, J11, J12, J21, J22
      CALL JAC_DSQK ( XLOC, XI, ETA, J11, J12, J21, J22, DETJ )
      END FUNCTION DETJ_DSQK

      SUBROUTINE JAC_DSQK ( XLOC, XI, ETA, J11, J12, J21, J22, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XLOC(4,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: J11, J12, J21, J22, DETJ
      REAL(DOUBLE) :: DN(2,4)
      CALL SHAPE_DN ( XI, ETA, DN )
      J11 = DOT_PRODUCT(DN(1,:), XLOC(:,1))
      J12 = DOT_PRODUCT(DN(1,:), XLOC(:,2))
      J21 = DOT_PRODUCT(DN(2,:), XLOC(:,1))
      J22 = DOT_PRODUCT(DN(2,:), XLOC(:,2))
      DETJ = J11*J22 - J12*J21
      END SUBROUTINE JAC_DSQK

      SUBROUTINE DNXY_DSQK ( XLOC, XI, ETA, DNDX, DNDY )
      REAL(DOUBLE), INTENT(IN)  :: XLOC(4,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: DNDX(4), DNDY(4)
      REAL(DOUBLE) :: DN(2,4), J11, J12, J21, J22, DETJ, I11, I12, I21, I22
      INTEGER(LONG) :: II
      CALL SHAPE_DN ( XI, ETA, DN )
      CALL JAC_DSQK ( XLOC, XI, ETA, J11, J12, J21, J22, DETJ )
      I11 =  J22 / DETJ
      I12 = -J12 / DETJ
      I21 = -J21 / DETJ
      I22 =  J11 / DETJ
      DO II=1,4
         DNDX(II) = I11*DN(1,II) + I12*DN(2,II)
         DNDY(II) = I21*DN(1,II) + I22*DN(2,II)
      ENDDO
      END SUBROUTINE DNXY_DSQK

      FUNCTION BM_DSQK ( XLOC, XI, ETA ) RESULT(BMOUT)
      REAL(DOUBLE), INTENT(IN) :: XLOC(4,2), XI, ETA
      REAL(DOUBLE) :: BMOUT(3,20), DNDX(4), DNDY(4)
      INTEGER(LONG) :: II, COL
      CALL DNXY_DSQK ( XLOC, XI, ETA, DNDX, DNDY )
      BMOUT = ZERO
      DO II=1,4
         COL = (II-1)*5
         BMOUT(1,COL+1) = DNDX(II)
         BMOUT(2,COL+2) = DNDY(II)
         BMOUT(3,COL+1) = DNDY(II)
         BMOUT(3,COL+2) = DNDX(II)
      ENDDO
      END FUNCTION BM_DSQK

      FUNCTION BBU_DSQK ( XLOC, XI, ETA ) RESULT(BBOUT)
      REAL(DOUBLE), INTENT(IN) :: XLOC(4,2), XI, ETA
      REAL(DOUBLE) :: BBOUT(3,20), DNDX(4), DNDY(4)
      INTEGER(LONG) :: II, COL
      CALL DNXY_DSQK ( XLOC, XI, ETA, DNDX, DNDY )
      BBOUT = ZERO
      DO II=1,4
         COL = (II-1)*5
         BBOUT(1,COL+5) = DNDX(II)
         BBOUT(2,COL+4) = -DNDY(II)
         BBOUT(3,COL+5) = DNDY(II)
         BBOUT(3,COL+4) = -DNDX(II)
      ENDDO
      END FUNCTION BBU_DSQK

      FUNCTION BBDELTA_DSQK ( XLOC, CK, SK, XI, ETA ) RESULT(BDOUT)
      REAL(DOUBLE), INTENT(IN) :: XLOC(4,2), CK(4), SK(4), XI, ETA
      REAL(DOUBLE) :: BDOUT(3,4), DPXI(4), DPETA(4), DNDX(4), DNDY(4)
      REAL(DOUBLE) :: J11, J12, J21, J22, DETJ, I11, I12, I21, I22
      INTEGER(LONG) :: II
      CALL PDERIV_DSQK ( XI, ETA, DPXI, DPETA )
      CALL JAC_DSQK ( XLOC, XI, ETA, J11, J12, J21, J22, DETJ )
      I11 =  J22 / DETJ
      I12 = -J12 / DETJ
      I21 = -J21 / DETJ
      I22 =  J11 / DETJ
      BDOUT = ZERO
      DO II=1,4
         DNDX(II) = I11*DPXI(II) + I12*DPETA(II)
         DNDY(II) = I21*DPXI(II) + I22*DPETA(II)
         BDOUT(1,II) = CK(II)*DNDX(II)
         BDOUT(2,II) = SK(II)*DNDY(II)
         BDOUT(3,II) = CK(II)*DNDY(II) + SK(II)*DNDX(II)
      ENDDO
      END FUNCTION BBDELTA_DSQK

      SUBROUTINE PDERIV_DSQK ( XI, ETA, DPXI, DPETA )
      REAL(DOUBLE), INTENT(IN)  :: XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: DPXI(4), DPETA(4)
      DPXI(1)  = -XI*(ONE - ETA)
      DPXI(2)  =  0.5D0*(ONE - ETA*ETA)
      DPXI(3)  = -XI*(ONE + ETA)
      DPXI(4)  = -0.5D0*(ONE - ETA*ETA)
      DPETA(1) = -0.5D0*(ONE - XI*XI)
      DPETA(2) = -ETA*(ONE + XI)
      DPETA(3) =  0.5D0*(ONE - XI*XI)
      DPETA(4) = -ETA*(ONE - XI)
      END SUBROUTINE PDERIV_DSQK

      FUNCTION AN_DSQK ( CK, SK, LK ) RESULT(ANOUT)
      REAL(DOUBLE), INTENT(IN) :: CK(4), SK(4), LK(4)
      REAL(DOUBLE) :: ANOUT(4,20)
      INTEGER(LONG) :: K, I1, J1, CI, CJ
      INTEGER(LONG), PARAMETER :: EDGE_I(4) = (/1,2,3,4/)
      INTEGER(LONG), PARAMETER :: EDGE_J(4) = (/2,3,4,1/)
      ANOUT = ZERO
      DO K=1,4
         I1 = EDGE_I(K)
         J1 = EDGE_J(K)
         CI = (I1-1)*5
         CJ = (J1-1)*5
         ANOUT(K,CI+3) = ANOUT(K,CI+3) + 1.5D0/LK(K)
         ANOUT(K,CJ+3) = ANOUT(K,CJ+3) - 1.5D0/LK(K)
         ANOUT(K,CI+5) = ANOUT(K,CI+5) - 0.75D0*CK(K)
         ANOUT(K,CJ+5) = ANOUT(K,CJ+5) - 0.75D0*CK(K)
         ANOUT(K,CI+4) = ANOUT(K,CI+4) + 0.75D0*SK(K)
         ANOUT(K,CJ+4) = ANOUT(K,CJ+4) + 0.75D0*SK(K)
      ENDDO
      END FUNCTION AN_DSQK

      FUNCTION BS_TYING_DSQK ( XLOC, XI, ETA ) RESULT(BSOUT)
      REAL(DOUBLE), INTENT(IN) :: XLOC(4,2), XI, ETA
      REAL(DOUBLE) :: BSOUT(2,20), NVAL(4), DNDX(4), DNDY(4)
      INTEGER(LONG) :: II, COL
      CALL SHAPE_N ( XI, ETA, NVAL )
      CALL DNXY_DSQK ( XLOC, XI, ETA, DNDX, DNDY )
      BSOUT = ZERO
      DO II=1,4
         COL = (II-1)*5
         BSOUT(1,COL+3) = DNDX(II)
         BSOUT(1,COL+5) = NVAL(II)
         BSOUT(2,COL+3) = DNDY(II)
         BSOUT(2,COL+4) = -NVAL(II)
      ENDDO
      END FUNCTION BS_TYING_DSQK

      FUNCTION BS_DSQK ( XLOC, XI, ETA ) RESULT(BSOUT)
      REAL(DOUBLE), INTENT(IN) :: XLOC(4,2), XI, ETA
      REAL(DOUBLE) :: BSOUT(2,20), BSA(2,20), BSB(2,20), BSC(2,20), BSD(2,20)
      BSA = BS_TYING_DSQK ( XLOC, ZERO, -ONE )
      BSB = BS_TYING_DSQK ( XLOC, ZERO,  ONE )
      BSC = BS_TYING_DSQK ( XLOC, -ONE, ZERO )
      BSD = BS_TYING_DSQK ( XLOC,  ONE, ZERO )
      BSOUT = ZERO
      BSOUT(1,:) = 0.5D0*(ONE - ETA)*BSA(1,:) + 0.5D0*(ONE + ETA)*BSB(1,:)
      BSOUT(2,:) = 0.5D0*(ONE - XI )*BSC(2,:) + 0.5D0*(ONE + XI )*BSD(2,:)
      END FUNCTION BS_DSQK

      SUBROUTINE MASS_DSQK_AT ( XI, ETA, FACTOR, MPA, THK, MOUT )
      REAL(DOUBLE), INTENT(IN)    :: XI, ETA, FACTOR, MPA, THK
      REAL(DOUBLE), INTENT(INOUT) :: MOUT(20,20)
      REAL(DOUBLE) :: NVAL(4), RHOI
      INTEGER(LONG) :: II, COL
      CALL SHAPE_N ( XI, ETA, NVAL )
      RHOI = MPA*THK*THK / 12.0D0
      DO II=1,4
         COL = (II-1)*5
         MOUT(COL+1,COL+1) = MOUT(COL+1,COL+1) + MPA*NVAL(II)*FACTOR
         MOUT(COL+2,COL+2) = MOUT(COL+2,COL+2) + MPA*NVAL(II)*FACTOR
         MOUT(COL+3,COL+3) = MOUT(COL+3,COL+3) + MPA*NVAL(II)*FACTOR
         MOUT(COL+4,COL+4) = MOUT(COL+4,COL+4) + RHOI*NVAL(II)*FACTOR
         MOUT(COL+5,COL+5) = MOUT(COL+5,COL+5) + RHOI*NVAL(II)*FACTOR
      ENDDO
      END SUBROUTINE MASS_DSQK_AT

      SUBROUTINE GEOM_SURFACE_AT ( XYZN, XI, ETA, SVEC, JDET, NORMV )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(4,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: SVEC(3), JDET, NORMV(3)
      REAL(DOUBLE) :: DN(2,4), A1(3), A2(3)
      CALL SHAPE_DN ( XI, ETA, DN )
      A1 = MATMUL(DN(1,:), XYZN)
      A2 = MATMUL(DN(2,:), XYZN)
      CALL CROSS3(A1, A2, SVEC)
      JDET = VNORM(SVEC)
      IF (JDET > 1.0D-15) THEN
         NORMV = SVEC / JDET
      ELSE
         NORMV = (/ZERO, ZERO, ONE/)
      ENDIF
      END SUBROUTINE GEOM_SURFACE_AT

      SUBROUTINE PRESSURE_VECTOR_AT ( XI, ETA, FACTOR, NORMV, POUT )
      REAL(DOUBLE), INTENT(IN)    :: XI, ETA, FACTOR, NORMV(3)
      REAL(DOUBLE), INTENT(INOUT) :: POUT(24)
      REAL(DOUBLE) :: NVAL(4)
      INTEGER(LONG) :: II, BASE
      CALL SHAPE_N ( XI, ETA, NVAL )
      DO II=1,4
         BASE = (II-1)*6
         POUT(BASE+1:BASE+3) = POUT(BASE+1:BASE+3) + NVAL(II)*FACTOR*NORMV
      ENDDO
      END SUBROUTINE PRESSURE_VECTOR_AT

      SUBROUTINE BUILD_G_DSQK ( XLOC, XI, ETA, GU, GV, GW )
      REAL(DOUBLE), INTENT(IN)  :: XLOC(4,2), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: GU(2,20), GV(2,20), GW(2,20)
      REAL(DOUBLE) :: DNDX(4), DNDY(4)
      INTEGER(LONG) :: II, COL
      CALL DNXY_DSQK ( XLOC, XI, ETA, DNDX, DNDY )
      GU = ZERO
      GV = ZERO
      GW = ZERO
      DO II=1,4
         COL = (II-1)*5
         GU(1,COL+1) = DNDX(II)
         GU(2,COL+1) = DNDY(II)
         GV(1,COL+2) = DNDX(II)
         GV(2,COL+2) = DNDY(II)
         GW(1,COL+3) = DNDX(II)
         GW(2,COL+3) = DNDY(II)
      ENDDO
      END SUBROUTINE BUILD_G_DSQK

      SUBROUTINE BUILD_DSQK_DRILL ( E3, GVAL, THK, XLOC, KOUT )
      REAL(DOUBLE), INTENT(IN)  :: E3(3), GVAL, THK, XLOC(4,2)
      REAL(DOUBLE), INTENT(OUT) :: KOUT(24,24)
      REAL(DOUBLE) :: AREA, DRK
      INTEGER(LONG) :: II, BASE, R1, C1
      KOUT = ZERO
      AREA = AREA_DSQK ( XLOC )
      DRK = DRILL_SCALE * GVAL * THK * AREA
      DO II=1,4
         BASE = (II-1)*6
         DO R1=1,3
            DO C1=1,3
               KOUT(BASE+3+R1, BASE+3+C1) = KOUT(BASE+3+R1, BASE+3+C1) + DRK*E3(R1)*E3(C1)
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE BUILD_DSQK_DRILL

      FUNCTION AREA_DSQK ( XLOC ) RESULT(AREA)
      REAL(DOUBLE), INTENT(IN) :: XLOC(4,2)
      REAL(DOUBLE) :: AREA, SS2(MAX_ORDER_GAUSS), HH2(MAX_ORDER_GAUSS)
      INTEGER(LONG) :: I1, J1
      CALL ORDER_GAUSS ( 2, SS2, HH2 )
      AREA = ZERO
      DO I1=1,2
         DO J1=1,2
            AREA = AREA + HH2(I1)*HH2(J1)*DETJ_DSQK(XLOC, SS2(I1), SS2(J1))
         ENDDO
      ENDDO
      END FUNCTION AREA_DSQK

      FUNCTION GP_INDEX ( IGP, JGP ) RESULT(IDX)
      INTEGER(LONG), INTENT(IN) :: IGP, JGP
      INTEGER(LONG) :: IDX
      IF ((IGP == 1) .AND. (JGP == 1)) THEN
         IDX = 1
      ELSE IF ((IGP == 1) .AND. (JGP == 2)) THEN
         IDX = 2
      ELSE IF ((IGP == 2) .AND. (JGP == 1)) THEN
         IDX = 3
      ELSE
         IDX = 4
      ENDIF
      END FUNCTION GP_INDEX

      SUBROUTINE BUILD_T24 ( T3IN, TOUT )
      REAL(DOUBLE), INTENT(IN)  :: T3IN(3,3)
      REAL(DOUBLE), INTENT(OUT) :: TOUT(24,24)
      INTEGER(LONG) :: II, RR, CC, BASE
      TOUT = ZERO
      DO II=1,4
         BASE = (II-1)*6
         DO RR=1,3
            DO CC=1,3
               TOUT(BASE+RR,   BASE+CC  ) = T3IN(RR,CC)
               TOUT(BASE+RR+3, BASE+CC+3) = T3IN(RR,CC)
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE BUILD_T24

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

      FUNCTION SAFE_UNIT ( V, VFALL ) RESULT(UV)
      REAL(DOUBLE), INTENT(IN) :: V(3), VFALL(3)
      REAL(DOUBLE) :: UV(3), NM
      NM = VNORM(V)
      IF (NM > 1.0D-15) THEN
         UV = V / NM
      ELSE
         NM = VNORM(VFALL)
         IF (NM > 1.0D-15) THEN
            UV = VFALL / NM
         ELSE
            UV = (/ONE, ZERO, ZERO/)
         ENDIF
      ENDIF
      END FUNCTION SAFE_UNIT

      SUBROUTINE DEBUG_PRINT_MATRIX ( TITLE, MAT )
      CHARACTER(LEN=*), INTENT(IN) :: TITLE
      REAL(DOUBLE), INTENT(IN)     :: MAT(:,:)
      INTEGER(LONG) :: II
      WRITE(F06,'(A)') ' '
      WRITE(F06,'(A)') TRIM(TITLE)
      DO II=1,SIZE(MAT,1)
         WRITE(F06,'(100(1X,ES15.7))') MAT(II,:)
      ENDDO
      END SUBROUTINE DEBUG_PRINT_MATRIX

      END SUBROUTINE CQUAD4_DSQK_RHR
