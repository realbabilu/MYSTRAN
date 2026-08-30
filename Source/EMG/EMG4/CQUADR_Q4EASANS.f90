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

      SUBROUTINE CQUADR_Q4EASANS ( OPT, INT_ELEM_ID )

! --- shell_renovation begin --- !
! Q4 EAS + ANS Hughes-Brezzi shell element port based on the Python reference:
!   D:\18a\bending_only\Shell\gemini2\shit\Q4EAS_ANS_HB_ShellElement.py
!
! Phase-1 scope:
!   - linear stiffness
!   - shell stress recovery matrices
!   - translational lumped mass for selfweight/static use
!   - face pressure vector
!   - debug dumps for parity work
!
! Notes on mass/load policy:
!   - This branch keeps MYSTRAN's common CQUADR mass and pressure handling.
!   - The stiffness path for PARAM,QUADRTYP,Q4EASANS uses the Q4EAS/ANS/HB B/EAS
!     routines below.
! --- shell_renovation end --- !

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_GAUSS, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, FOUR
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  COUPMASS, QUADRTYP
      USE MODEL_STUF, ONLY            :  EID, ELGP, KE, KED, ME, BE1, BE2, BE3, EM, EB, ET, EPROP, MASS_PER_UNIT_AREA, PRESS, PPE,&
                                         PTE, TE, TYPE, NUM_EMG_FATAL_ERRS, SHELL_A, SHELL_D, SHELL_T, FCONV, STRESS, ALPVEC, DT, &
                                         TREF

      USE ELMDIS_Interface
      USE ELEM_STRE_STRN_ARRAYS_Interface
      USE ORDER_GAUSS_Interface
      USE OUTA_HERE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CQUADR_Q4EASANS'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      INTEGER(LONG), PARAMETER        :: NNODE = 4
      INTEGER(LONG), PARAMETER        :: NDOFN = 6
      INTEGER(LONG), PARAMETER        :: NDOF = NNODE*NDOFN
      INTEGER(LONG), PARAMETER        :: NSTRA = 3
      INTEGER(LONG), PARAMETER        :: NSHEAR = 2
      REAL(DOUBLE), PARAMETER         :: CQUADR_DRILL_SCALE = 1.0D0

      INTEGER(LONG)                   :: I, J, K, GP, JSUB, STR_PT_NUM, IA, IB, RR
      REAL(DOUBLE)                    :: XYZ(4,3), NORMALS(4,3), T24(24,24), T24T(24,24)
      REAL(DOUBLE)                    :: SS(MAX_ORDER_GAUSS), HH(MAX_ORDER_GAUSS)
      REAL(DOUBLE)                    :: XI, ETA, WT
      REAL(DOUBLE)                    :: TV1(3), TV2(3), NVEC(3), JDET, CO(2,2), BCMAT(2,2)
      REAL(DOUBLE)                    :: DN_G(2,4), DNDX(4), DNDY(4), SIG0(2,2), KGVAL
      REAL(DOUBLE)                    :: BMB(3,24), BBB(3,24), BSB(2,24), BMEFF(3,24)
      REAL(DOUBLE)                    :: BML(3,24), BBL(3,24), BSL(2,24)
      REAL(DOUBLE)                    :: KLOCAL(24,24), KBASIC(24,24), KMEM(24,24), KBEND(24,24), KSHEAR(24,24), KDRILL(24,24)
      REAL(DOUBLE)                    :: KUA(24,4), KAA(4,4), KAAINV(4,4), MEAS(3,4)
      REAL(DOUBLE)                    :: KGBASIC(24,24), KGLOCAL(24,24)
      REAL(DOUBLE)                    :: M1(4,4), MBASIC(24,24), MLOCAL(24,24), NVG(4), MDIAG(4), DENS_A
      REAL(DOUBLE)                    :: MASS_AREA_INT, MASS_ELEM_SUM
      REAL(DOUBLE)                    :: UNIT_PPE_B(24), UNIT_PPE_L(24), UNIT_PTE_B(24)
      REAL(DOUBLE)                    :: CTE3(3), NTH(3), TBAR
      REAL(DOUBLE)                    :: GBE1(3,24,4), GBE2(3,24,4), GBE3(2,24,4)
      LOGICAL                         :: Q4EASANS_MODE

! **********************************************************************************************************************************

      Q4EASANS_MODE = ((TYPE == 'QUADR   ') .AND. (QUADRTYP == 'Q4EASANS'))

      IF (ELGP /= 4) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9001) SUBR_NAME, EID, ELGP
         WRITE(F06,9001) SUBR_NAME, EID, ELGP
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      XYZ = ZERO

! The Q4EASANS reference B-operators below are evaluated from basic/global
! coordinates and nodal directors.  KBASIC is intentionally comparable to
! Q4EAS_ANS_HB_ShellElement.py for parity/debugging, then transformed into
! MYSTRAN's element-local assembly basis with T24.

      CALL LOAD_BASIC_COORDS ( XYZ )
      CALL CALC_NODAL_NORMALS ( XYZ, NORMALS )
      CALL BUILD_T24 ( TE, T24 )
      T24T = TRANSPOSE(T24)

      KBASIC = ZERO
      KMEM = ZERO
      KBEND = ZERO
      KSHEAR = ZERO
      KLOCAL = ZERO
      GBE1 = ZERO
      GBE2 = ZERO
      GBE3 = ZERO

      IF ((OPT(3) == 'Y') .OR. (OPT(4) == 'Y') .OR. (OPT(5) == 'Y') .OR. (OPT(2) == 'Y') .OR. (OPT(1) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         CALL ORDER_GAUSS(2, SS, HH)
      ENDIF

      IF ((OPT(3) == 'Y') .OR. (OPT(4) == 'Y')) THEN
         KUA = ZERO
         KAA = ZERO

         DO I=1,2
            DO J=1,2
               XI  = SS(I)
               ETA = SS(J)
               WT  = HH(I)*HH(J)

               CALL GEOMETRY_AT(XYZ, NORMALS, XI, ETA, TV1, TV2, NVEC, JDET, CO, BCMAT)
                BMB = BM_Q4EASANS_AT(XYZ, XI, ETA)
                BBB = BB_Q4EASANS_AT(XYZ, NORMALS, XI, ETA)
                BSB = BS_Q4EASANS_ANS_AT(XYZ, NORMALS, XI, ETA)

                IF ((DEBUG(190) > 0) .AND. (I == 1) .AND. (J == 1)) THEN
                   CALL DEBUG_PRINT_MATRIX('CQUADR GP11 BMB', BMB)
                   CALL DEBUG_PRINT_MATRIX('CQUADR GP11 BBB', BBB)
                   CALL DEBUG_PRINT_MATRIX('CQUADR GP11 BSB', BSB)
                ENDIF

                BML = MATMUL(BMB, T24T)
                BBL = MATMUL(BBB, T24T)
                BSL = MATMUL(BSB, T24T)

               KMEM   = KMEM   + WT*JDET*MATMUL(TRANSPOSE(BMB), MATMUL(SHELL_A, BMB))
               MEAS = EAS_Q4EASANS_AT(XYZ, XI, ETA)
               KUA = KUA + WT*JDET*MATMUL(TRANSPOSE(BMB), MATMUL(SHELL_A, MEAS))
               KAA = KAA + WT*JDET*MATMUL(TRANSPOSE(MEAS), MATMUL(SHELL_A, MEAS))
               KBEND  = KBEND  + WT*JDET*MATMUL(TRANSPOSE(BBB), MATMUL(SHELL_D, BBB))
               KSHEAR = KSHEAR + WT*JDET*MATMUL(TRANSPOSE(BSB), MATMUL(SHELL_T, BSB))
               KBASIC = KMEM + KBEND + KSHEAR

               GP = GP_INDEX(I, J)
               GBE1(1:3,1:24,GP) = BML
               GBE2(1:3,1:24,GP) = BBL
               GBE3(1:2,1:24,GP) = BSL
            ENDDO
         ENDDO

         KDRILL = DRILL_STIFFNESS_Q4EASANS(XYZ, NORMALS)
         CALL INV4(KAA, KAAINV)
         KMEM = KMEM - MATMUL(KUA, MATMUL(KAAINV, TRANSPOSE(KUA)))
         KBASIC = KMEM + KBEND + KSHEAR
         KBASIC = KBASIC + KDRILL
        KLOCAL = MATMUL(T24, MATMUL(KBASIC, T24T))

         IF (OPT(3) == 'Y') THEN
            DO I=1,2
               DO J=1,2
                  XI  = SS(I)
                  ETA = SS(J)
                  BMB = BM_Q4EASANS_AT(XYZ, XI, ETA)
                  MEAS = EAS_Q4EASANS_AT(XYZ, XI, ETA)
                  BMEFF = BMB - MATMUL(MEAS, MATMUL(KAAINV, TRANSPOSE(KUA)))
                  GP = GP_INDEX(I, J)
                  GBE1(1:3,1:24,GP) = MATMUL(BMEFF, T24T)
               ENDDO
            ENDDO
         ENDIF

         IF (OPT(4) == 'Y') THEN
            KE(1:24,1:24) = KLOCAL
         ENDIF

         IF (OPT(3) == 'Y') THEN
            BE1(:,:,1) = (GBE1(:,:,1) + GBE1(:,:,2) + GBE1(:,:,3) + GBE1(:,:,4)) / FOUR
            BE2(:,:,1) = (GBE2(:,:,1) + GBE2(:,:,2) + GBE2(:,:,3) + GBE2(:,:,4)) / FOUR
            BE3(1:2,:,1) = (GBE3(1:2,:,1) + GBE3(1:2,:,2) + GBE3(1:2,:,3) + GBE3(1:2,:,4)) / FOUR

            ! CQUAD4 shell output ordering: center, then (+,+), (+,-), (-,+), (-,-)
            BE1(:,:,2) = GBE1(:,:,4)
            BE1(:,:,3) = GBE1(:,:,3)
            BE1(:,:,4) = GBE1(:,:,2)
            BE1(:,:,5) = GBE1(:,:,1)

            BE2(:,:,2) = GBE2(:,:,4)
            BE2(:,:,3) = GBE2(:,:,3)
            BE2(:,:,4) = GBE2(:,:,2)
            BE2(:,:,5) = GBE2(:,:,1)

            BE3(1:2,:,2) = GBE3(1:2,:,4)
            BE3(1:2,:,3) = GBE3(1:2,:,3)
            BE3(1:2,:,4) = GBE3(1:2,:,2)
            BE3(1:2,:,5) = GBE3(1:2,:,1)
         ENDIF
      ENDIF

      IF (OPT(1) == 'Y') THEN
         M1 = ZERO
         MASS_AREA_INT = ZERO
         DO I=1,2
            DO J=1,2
               XI  = SS(I)
               ETA = SS(J)
               WT  = HH(I)*HH(J)
               CALL GEOMETRY_AT(XYZ, NORMALS, XI, ETA, TV1, TV2, NVEC, JDET, CO, BCMAT)
               CALL SHAPE_N(XI, ETA, NVG)
               MASS_AREA_INT = MASS_AREA_INT + WT*JDET
               DO GP=1,4
                  DO K=1,4
                     M1(GP,K) = M1(GP,K) + NVG(GP)*NVG(K)*MASS_PER_UNIT_AREA*WT*JDET
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         MBASIC = ZERO
         IF ((SOL_NAME(1:5) == 'MODES') .AND. (COUPMASS > 0)) THEN
            DO I=1,4
               DO J=1,4
                  DO K=1,3
                     MBASIC((I-1)*6+K,(J-1)*6+K) = M1(I,J)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            MDIAG = ZERO
            DO I=1,4
               MDIAG(I) = SUM(M1(I,1:4))
            ENDDO
            DO I=1,4
               DO K=1,3
                  MBASIC((I-1)*6+K,(I-1)*6+K) = MDIAG(I)
               ENDDO
            ENDDO
         ENDIF
         MASS_ELEM_SUM = SUM(M1)
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

      IF (OPT(2) == 'Y') THEN
         UNIT_PTE_B = ZERO
         CTE3(1) = ALPVEC(1,1)
         CTE3(2) = ALPVEC(2,1)
         CTE3(3) = ALPVEC(4,1)
         NTH = MATMUL(SHELL_A, CTE3)
         DO I=1,2
            DO J=1,2
               XI  = SS(I)
               ETA = SS(J)
               WT  = HH(I)*HH(J)
               CALL GEOMETRY_AT(XYZ, NORMALS, XI, ETA, TV1, TV2, NVEC, JDET, CO, BCMAT)
               BMB = BM_Q4EASANS_AT(XYZ, XI, ETA)
               UNIT_PTE_B = UNIT_PTE_B + WT*JDET*MATMUL(TRANSPOSE(BMB), NTH)
            ENDDO
         ENDDO
         DO JSUB=1,SIZE(PTE,2)
            TBAR = (DT(1,JSUB) + DT(2,JSUB) + DT(3,JSUB) + DT(4,JSUB))/FOUR - TREF(1)
            PTE(1:24,JSUB) = MATMUL(T24, UNIT_PTE_B) * TBAR
         ENDDO
      ENDIF

      IF (OPT(6) == 'Y') THEN
         CALL ELMDIS

         KGLOCAL = ZERO
         DO I=1,2
            DO J=1,2
               XI  = SS(I)
               ETA = SS(J)
               WT  = HH(I)*HH(J)
               CALL SHAPE_DN(XI, ETA, DN_G)
               CALL GEOMETRY_AT(XYZ, NORMALS, XI, ETA, TV1, TV2, NVEC, JDET, CO, BCMAT)
               BMB = BM_Q4EASANS_AT(XYZ, XI, ETA)
               BML = MATMUL(BMB, T24T)
               BE1(1:3,1:24,1) = BML
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )

               SIG0(1,1) = FCONV(1)*STRESS(1)
               SIG0(2,2) = FCONV(1)*STRESS(2)
               SIG0(1,2) = FCONV(1)*STRESS(3)
               SIG0(2,1) = SIG0(1,2)
! --- shell_renovation begin --- !
! DEBUG(235) probes whether ELEM_STRE_STRN_ARRAYS is returning membrane stress
! instead of membrane resultant for this path, matching the MITC4 KGGD
! TPLY*STRESS pattern.  It is intentionally debug-gated until validated.
               IF (DEBUG(235) > 0) THEN
                  SIG0(1,1) = EPROP(1)*FCONV(1)*STRESS(1)
                  SIG0(2,2) = EPROP(1)*FCONV(1)*STRESS(2)
                  SIG0(1,2) = EPROP(1)*FCONV(1)*STRESS(3)
                  SIG0(2,1) = SIG0(1,2)
               ENDIF

! DEBUG(234) isolates this shell buckling operator from MYSTRAN prebuckling
! stress recovery by imposing the same constant pure-Nx resultant used in the
! Python SAC068 comparison decks.
               IF (DEBUG(234) > 0) THEN
                  SIG0(1,1) = -ONE
                  SIG0(2,2) = ZERO
                  SIG0(1,2) = ZERO
                  SIG0(2,1) = ZERO
               ENDIF
! --- shell_renovation end --- !

               IF ((DEBUG(233) > 0) .AND. (EID <= 8)) THEN
                  WRITE(F06,'(A,I8,A,I2,A,I2,A,3(1X,ES15.7))') 'CQUADR KGGD EID=', EID, ' I=', I, ' J=', J,                      &
                                                               ' SIG0=', SIG0(1,1), SIG0(2,2), SIG0(1,2)
                  WRITE(F06,'(A,I8,A,3(1X,ES15.7))') 'CQUADR KGGD EID=', EID, ' NORMAL=', NVEC(1), NVEC(2), NVEC(3)
               ENDIF

               CALL BUILD_SPECIAL_KGLOCAL_Q4EASANS(DN_G, CO, SIG0, WT, JDET, KGLOCAL)
            ENDDO
         ENDDO

         KED(1:24,1:24) = KGLOCAL
         IF ((DEBUG(233) > 0) .AND. (EID <= 8)) THEN
            WRITE(F06,'(A,I8,A,ES15.7)') 'CQUADR KGGD EID=', EID, ' KED_NORM=', DSQRT(SUM(KED(1:24,1:24)*KED(1:24,1:24)))
            CALL DEBUG_PRINT_MATRIX('CQUADR KGGD KED', KED(1:24,1:24))
         ENDIF
      ENDIF

      IF (DEBUG(233) > 0) THEN
         IF (EID <= 8) THEN
            WRITE(F06,'(A,I8,A,ES15.7)') 'CQUADR KE EID=', EID, ' KBASIC_NORM=', DSQRT(SUM(KBASIC*KBASIC))
            WRITE(F06,'(A,I8,A,ES15.7)') 'CQUADR KE EID=', EID, ' KMEM_NORM=', DSQRT(SUM(KMEM*KMEM))
            WRITE(F06,'(A,I8,A,ES15.7)') 'CQUADR KE EID=', EID, ' KBEND_NORM=', DSQRT(SUM(KBEND*KBEND))
            WRITE(F06,'(A,I8,A,ES15.7)') 'CQUADR KE EID=', EID, ' KSHEAR_NORM=', DSQRT(SUM(KSHEAR*KSHEAR))
         ENDIF
         CALL DEBUG_PRINT_MATRIX('CQUADR KBASIC', KBASIC)
         CALL DEBUG_PRINT_MATRIX('CQUADR KLOCAL', KLOCAL)
         CALL DEBUG_PRINT_MATRIX('CQUADR KMEM', KMEM)
         CALL DEBUG_PRINT_MATRIX('CQUADR KBEND', KBEND)
         CALL DEBUG_PRINT_MATRIX('CQUADR KSHEAR', KSHEAR)
         CALL DEBUG_PRINT_MATRIX('CQUADR KDRILL', KDRILL)
          WRITE(F06,'(A,6(1X,ES14.6))') 'CQUADR EPROP1-6', (EPROP(I), I=1,6)
          CALL DEBUG_PRINT_MATRIX('CQUADR EM', EM)
          CALL DEBUG_PRINT_MATRIX('CQUADR EB', EB)
         CALL DEBUG_PRINT_MATRIX('CQUADR ET', ET)
         CALL DEBUG_PRINT_MATRIX('CQUADR SHELL_A', SHELL_A)
         CALL DEBUG_PRINT_MATRIX('CQUADR SHELL_D', SHELL_D)
         CALL DEBUG_PRINT_MATRIX('CQUADR SHELL_T', SHELL_T)
         IF (OPT(1) == 'Y') THEN
            WRITE(F06,'(A,1X,ES15.7)') 'CQUADR MASS_PER_UNIT_AREA', MASS_PER_UNIT_AREA
            WRITE(F06,'(A,1X,ES15.7)') 'CQUADR MASS_AREA_INT', MASS_AREA_INT
            WRITE(F06,'(A,1X,ES15.7)') 'CQUADR MASS_ELEM_SUM', MASS_ELEM_SUM
            WRITE(F06,'(A,1X,ES15.7)') 'CQUADR MBASIC_NORM', DSQRT(SUM(MBASIC*MBASIC))
            WRITE(F06,'(A,1X,ES15.7)') 'CQUADR MLOCAL_NORM', DSQRT(SUM(MLOCAL*MLOCAL))
         ENDIF
      ENDIF

      RETURN

 9001 FORMAT(' *ERROR: ',A,' expects ELGP=4 for element ',I8,' but got ',I8)

      CONTAINS

      SUBROUTINE LOAD_BASIC_COORDS ( XYZOUT )
      USE MODEL_STUF, ONLY : XEB
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

      SUBROUTINE CALC_NODAL_NORMALS ( XYZN, NORMS )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(4,3)
      REAL(DOUBLE), INTENT(OUT) :: NORMS(4,3)
      INTEGER(LONG) :: II
      REAL(DOUBLE) :: A(3), B(3), N(3), NM

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

! Match Q4EAS_ANS_HB_ShellElement.py: directors are generated from this
! element's own midsurface geometry, not from user/generated SNORM.
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
      IF (NM > 1.0D-15) THEN
         T2 = T2 / NM
      ENDIF

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

      FUNCTION BM_Q4EASANS_AT ( XYZN, XI, ETA ) RESULT(BMOUT)
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), XI, ETA
      REAL(DOUBLE) :: BMOUT(3,24), BMCOV(3,24), DN(2,4), G1(3), G2(3), E1C(3), E2C(3), E3C(3), TSTR(3,3)
      INTEGER(LONG) :: II, COL

      CALL SHAPE_DN(XI, ETA, DN)
      CALL SURFACE_BASIS(XYZN, XI, ETA, G1, G2, E1C, E2C, E3C)
      CALL SIMO_FIXED_FRAME(XYZN, E1C, E2C, E3C)
      CALL T_STRAIN_FIXED_AT(XYZN, XI, ETA, E1C, E2C, TSTR)
      CALL SURFACE_BASIS(XYZN, XI, ETA, G1, G2, E1C, E2C, E3C)

      BMCOV = ZERO
      DO II=1,4
         COL = (II-1)*6
         BMCOV(1,COL+1:COL+3) = DN(1,II)*G1
         BMCOV(2,COL+1:COL+3) = DN(2,II)*G2
         BMCOV(3,COL+1:COL+3) = DN(1,II)*G2 + DN(2,II)*G1
      ENDDO

      BMOUT = MATMUL(TSTR, BMCOV)
      END FUNCTION BM_Q4EASANS_AT

      FUNCTION BB_Q4EASANS_AT ( XYZN, NORMS, XI, ETA ) RESULT(BBOUT)
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), NORMS(4,3), XI, ETA
      REAL(DOUBLE) :: BBOUT(3,24), BBCOV(3,24), DN(2,4), G1(3), G2(3), E1C(3), E2C(3), E3C(3), TSTR(3,3)
      REAL(DOUBLE) :: T_1(3), T_2(3), T0(3), CROSS_G1(3), CROSS_G2(3)
      INTEGER(LONG) :: II, COL

      CALL SHAPE_DN(XI, ETA, DN)
      CALL SURFACE_BASIS(XYZN, XI, ETA, G1, G2, E1C, E2C, E3C)
      CALL SIMO_FIXED_FRAME(XYZN, E1C, E2C, E3C)
      CALL T_STRAIN_FIXED_AT(XYZN, XI, ETA, E1C, E2C, TSTR)

      T_1 = MATMUL(DN(1,:), NORMS)
      T_2 = MATMUL(DN(2,:), NORMS)

      BBCOV = ZERO
      DO II=1,4
         COL = (II-1)*6
         T0 = NORMS(II,:)
         BBCOV(1,COL+1:COL+3) = DN(1,II)*T_1
         BBCOV(2,COL+1:COL+3) = DN(2,II)*T_2
         BBCOV(3,COL+1:COL+3) = DN(1,II)*T_2 + DN(2,II)*T_1

         CALL CROSS3(T0, G1, CROSS_G1)
         CALL CROSS3(T0, G2, CROSS_G2)
         BBCOV(1,COL+4:COL+6) = DN(1,II)*CROSS_G1
         BBCOV(2,COL+4:COL+6) = DN(2,II)*CROSS_G2
         BBCOV(3,COL+4:COL+6) = DN(2,II)*CROSS_G1 + DN(1,II)*CROSS_G2
      ENDDO

      BBOUT = MATMUL(TSTR, BBCOV)
      END FUNCTION BB_Q4EASANS_AT

      FUNCTION BS_Q4EASANS_ANS_AT ( XYZN, NORMS, XI, ETA ) RESULT(BSOUT)
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), NORMS(4,3), XI, ETA
      REAL(DOUBLE) :: BSOUT(2,24), BSNAT(2,24), BSA(2,24), BSB(2,24), BSC(2,24), BSD(2,24)
      REAL(DOUBLE) :: G1C(3), G2C(3), E1C(3), E2C(3), E3C(3), TSH(2,2)

      BSA = BS_Q4EASANS_LOC_AT(XYZN, NORMS, ZERO, -ONE)
      BSB = BS_Q4EASANS_LOC_AT(XYZN, NORMS,  ONE, ZERO)
      BSC = BS_Q4EASANS_LOC_AT(XYZN, NORMS, ZERO,  ONE)
      BSD = BS_Q4EASANS_LOC_AT(XYZN, NORMS, -ONE, ZERO)

      BSNAT = ZERO
      BSNAT(1,:) = 0.5D0*(ONE - ETA)*BSA(1,:) + 0.5D0*(ONE + ETA)*BSC(1,:)
      BSNAT(2,:) = 0.5D0*(ONE - XI )*BSD(2,:) + 0.5D0*(ONE + XI )*BSB(2,:)

      CALL SIMO_FIXED_FRAME(XYZN, E1C, E2C, E3C)
      CALL T_SHEAR_FIXED_AT(XYZN, XI, ETA, E1C, E2C, TSH)
      BSOUT = MATMUL(TSH, BSNAT)
      END FUNCTION BS_Q4EASANS_ANS_AT

      FUNCTION BS_Q4EASANS_LOC_AT ( XYZN, NORMS, XI, ETA ) RESULT(BSOUT)
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), NORMS(4,3), XI, ETA
      REAL(DOUBLE) :: BSOUT(2,24), DN(2,4), NVAL(4), G1(3), G2(3), E1(3), E2(3), E3(3)
      REAL(DOUBLE) :: T0PT(3), T0I(3), C1(3), C2(3), NM
      INTEGER(LONG) :: II, COL

      CALL SHAPE_N(XI, ETA, NVAL)
      CALL SHAPE_DN(XI, ETA, DN)
      CALL SURFACE_BASIS(XYZN, XI, ETA, G1, G2, E1, E2, E3)

      T0PT = MATMUL(NVAL, NORMS)
      NM = VNORM(T0PT)
      IF (NM > 1.0D-15) THEN
         T0PT = T0PT / NM
      ELSE
         T0PT = E3
      ENDIF

      BSOUT = ZERO
      DO II=1,4
         COL = (II-1)*6
         T0I = NORMS(II,:)
         BSOUT(1,COL+1:COL+3) = BSOUT(1,COL+1:COL+3) + DN(1,II)*T0PT
         BSOUT(2,COL+1:COL+3) = BSOUT(2,COL+1:COL+3) + DN(2,II)*T0PT
         CALL CROSS3(T0I, G1, C1)
         CALL CROSS3(T0I, G2, C2)
         BSOUT(1,COL+4:COL+6) = BSOUT(1,COL+4:COL+6) + NVAL(II)*C1
         BSOUT(2,COL+4:COL+6) = BSOUT(2,COL+4:COL+6) + NVAL(II)*C2
      ENDDO
      END FUNCTION BS_Q4EASANS_LOC_AT

      FUNCTION EAS_Q4EASANS_AT ( XYZN, XI, ETA ) RESULT(MOUT)
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), XI, ETA
      REAL(DOUBLE) :: MOUT(3,4), MHAT(3,4), TSTR0(3,3)
      REAL(DOUBLE) :: G1(3), G2(3), E1(3), E2(3), E3(3), G10(3), G20(3), E10(3), E20(3), E30(3)
      REAL(DOUBLE) :: DETJ0, DETJ

      MHAT = ZERO
      MHAT(1,1) = XI
      MHAT(2,2) = ETA
      MHAT(3,3) = XI
      MHAT(3,4) = ETA

      CALL SIMO_FIXED_FRAME(XYZN, E10, E20, E30)
      CALL T_STRAIN_FIXED_AT(XYZN, ZERO, ZERO, E10, E20, TSTR0)
      CALL DET_FIXED_AT(XYZN, ZERO, ZERO, E10, E20, DETJ0)

      CALL SURFACE_BASIS(XYZN, XI, ETA, G1, G2, E1, E2, E3)
      CALL DET_FIXED_AT(XYZN, XI, ETA, E10, E20, DETJ)

      IF (DETJ > 1.0D-30) THEN
         MOUT = (DETJ0/DETJ)*MATMUL(TSTR0, MHAT)
      ELSE
         MOUT = ZERO
      ENDIF
      END FUNCTION EAS_Q4EASANS_AT

      FUNCTION DRILL_STIFFNESS_Q4EASANS ( XYZN, NORMS ) RESULT(KD)
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), NORMS(4,3)
      REAL(DOUBLE) :: KD(24,24), BD(1,24), XI, ETA, WT, G1(3), G2(3), E1(3), E2(3), E3(3), C(3), JAC
      REAL(DOUBLE) :: CDRILL
      INTEGER(LONG) :: I1, J1

      KD = ZERO
      CDRILL = 1.0D-4 * SHELL_T(1,1) / (5.0D0/6.0D0)
      DO I1=1,2
         DO J1=1,2
            XI = SS(I1)
            ETA = SS(J1)
            WT = HH(I1)*HH(J1)
            CALL SURFACE_BASIS(XYZN, XI, ETA, G1, G2, E1, E2, E3)
            CALL CROSS3(G1, G2, C)
            JAC = VNORM(C)
            BD = BDRILL_Q4EASANS_AT(XYZN, NORMS, XI, ETA)
            KD = KD + WT*JAC*CDRILL*MATMUL(TRANSPOSE(BD), BD)
         ENDDO
      ENDDO
      END FUNCTION DRILL_STIFFNESS_Q4EASANS

      FUNCTION BDRILL_Q4EASANS_AT ( XYZN, NORMS, XI, ETA ) RESULT(BD)
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), NORMS(4,3), XI, ETA
      REAL(DOUBLE) :: BD(1,24), DN(2,4), NVAL(4), DLOC(2,4), G1(3), G2(3), E1(3), E2(3), E3(3)
      INTEGER(LONG) :: II, COL

      CALL SHAPE_N(XI, ETA, NVAL)
      CALL SHAPE_DN(XI, ETA, DN)
      CALL SURFACE_BASIS(XYZN, XI, ETA, G1, G2, E1, E2, E3)
      CALL LOCAL_DERIVS(DN, G1, G2, E1, E2, DLOC)

      BD = ZERO
      DO II=1,4
         COL = (II-1)*6
         BD(1,COL+1:COL+3) = 0.5D0*(DLOC(1,II)*E2 - DLOC(2,II)*E1)
         BD(1,COL+4:COL+6) = BD(1,COL+4:COL+6) - NVAL(II)*NORMS(II,:)
      ENDDO
      END FUNCTION BDRILL_Q4EASANS_AT

      SUBROUTINE SURFACE_BASIS ( XYZN, XI, ETA, G1, G2, E1, E2, E3 )
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: G1(3), G2(3), E1(3), E2(3), E3(3)
      REAL(DOUBLE) :: DN(2,4), G3(3), TMP(3), NM

      CALL SHAPE_DN(XI, ETA, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      CALL CROSS3(G1, G2, G3)
      NM = VNORM(G3)
      IF (NM > 1.0D-15) THEN
         E3 = G3 / NM
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
         E1 = TMP / NM
      ELSE
         E1 = (/ONE, ZERO, ZERO/)
      ENDIF
      CALL CROSS3(E3, E1, E2)
      NM = VNORM(E2)
      IF (NM > 1.0D-15) E2 = E2 / NM
      END SUBROUTINE SURFACE_BASIS

      SUBROUTINE SIMO_FIXED_FRAME ( XYZN, E1F, E2F, E3F )
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3)
      REAL(DOUBLE), INTENT(OUT) :: E1F(3), E2F(3), E3F(3)
      REAL(DOUBLE) :: G1C(3), G2C(3), ZSGN

      CALL SURFACE_BASIS(XYZN, ZERO, ZERO, G1C, G2C, E1F, E2F, E3F)

! For flat XY patch-test decks, report/recover in the global XY frame so
! CQUADR/SIMO F06 engineering forces are directly comparable with DKMQ24 and
! the SAP/MSC reference values. For general shell geometry, keep the v1p6
! center-fixed frame.
      IF ((DABS(E3F(1)) + DABS(E3F(2))) <= 1.0D-12) THEN
         ZSGN = ONE
         IF (E3F(3) < ZERO) ZSGN = -ONE
         E1F = (/ONE, ZERO, ZERO/)
         E2F = (/ZERO, ZSGN, ZERO/)
         E3F = (/ZERO, ZERO, ZSGN/)
      ENDIF
      END SUBROUTINE SIMO_FIXED_FRAME

      SUBROUTINE T_STRAIN_FIXED_AT ( XYZN, XI, ETA, E1F, E2F, TSTR )
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), XI, ETA, E1F(3), E2F(3)
      REAL(DOUBLE), INTENT(OUT) :: TSTR(3,3)
      REAL(DOUBLE) :: G1(3), G2(3), E1TMP(3), E2TMP(3), E3TMP(3)
      REAL(DOUBLE) :: A11, A22, A12, DET, AINV11, AINV22, AINV12
      REAL(DOUBLE) :: GC1(3), GC2(3), C11, C12, C21, C22

      CALL SURFACE_BASIS(XYZN, XI, ETA, G1, G2, E1TMP, E2TMP, E3TMP)
      A11 = DOT_PRODUCT(G1, G1)
      A22 = DOT_PRODUCT(G2, G2)
      A12 = DOT_PRODUCT(G1, G2)
      DET = A11*A22 - A12*A12
      IF (DABS(DET) <= 1.0D-30) THEN
         TSTR = ZERO
         RETURN
      ENDIF

      AINV11 =  A22 / DET
      AINV22 =  A11 / DET
      AINV12 = -A12 / DET
      GC1 = AINV11*G1 + AINV12*G2
      GC2 = AINV12*G1 + AINV22*G2

      C11 = DOT_PRODUCT(E1F, GC1)
      C12 = DOT_PRODUCT(E1F, GC2)
      C21 = DOT_PRODUCT(E2F, GC1)
      C22 = DOT_PRODUCT(E2F, GC2)

      TSTR(1,1) = C11*C11
      TSTR(1,2) = C12*C12
      TSTR(1,3) = C11*C12
      TSTR(2,1) = C21*C21
      TSTR(2,2) = C22*C22
      TSTR(2,3) = C21*C22
      TSTR(3,1) = TWO*C11*C21
      TSTR(3,2) = TWO*C12*C22
      TSTR(3,3) = C11*C22 + C12*C21
      END SUBROUTINE T_STRAIN_FIXED_AT

      SUBROUTINE T_SHEAR_FIXED_AT ( XYZN, XI, ETA, E1F, E2F, TSH )
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), XI, ETA, E1F(3), E2F(3)
      REAL(DOUBLE), INTENT(OUT) :: TSH(2,2)
      REAL(DOUBLE) :: G1(3), G2(3), E1TMP(3), E2TMP(3), E3TMP(3)
      REAL(DOUBLE) :: A11, A22, A12, DET, AINV11, AINV22, AINV12, GC1(3), GC2(3)

      CALL SURFACE_BASIS(XYZN, XI, ETA, G1, G2, E1TMP, E2TMP, E3TMP)
      A11 = DOT_PRODUCT(G1, G1)
      A22 = DOT_PRODUCT(G2, G2)
      A12 = DOT_PRODUCT(G1, G2)
      DET = A11*A22 - A12*A12
      IF (DABS(DET) <= 1.0D-30) THEN
         TSH = ZERO
         RETURN
      ENDIF

      AINV11 =  A22 / DET
      AINV22 =  A11 / DET
      AINV12 = -A12 / DET
      GC1 = AINV11*G1 + AINV12*G2
      GC2 = AINV12*G1 + AINV22*G2

      TSH(1,1) = DOT_PRODUCT(E1F, GC1)
      TSH(1,2) = DOT_PRODUCT(E1F, GC2)
      TSH(2,1) = DOT_PRODUCT(E2F, GC1)
      TSH(2,2) = DOT_PRODUCT(E2F, GC2)
      END SUBROUTINE T_SHEAR_FIXED_AT

      SUBROUTINE DET_FIXED_AT ( XYZN, XI, ETA, E1F, E2F, DETJ )
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), XI, ETA, E1F(3), E2F(3)
      REAL(DOUBLE), INTENT(OUT) :: DETJ
      REAL(DOUBLE) :: G1(3), G2(3), E1TMP(3), E2TMP(3), E3TMP(3)

      CALL SURFACE_BASIS(XYZN, XI, ETA, G1, G2, E1TMP, E2TMP, E3TMP)
      DETJ = DOT_PRODUCT(G1,E1F)*DOT_PRODUCT(G2,E2F) - DOT_PRODUCT(G1,E2F)*DOT_PRODUCT(G2,E1F)
      END SUBROUTINE DET_FIXED_AT

      SUBROUTINE LOCAL_DERIVS ( DN, G1, G2, E1, E2, DLOC )
      REAL(DOUBLE), INTENT(IN) :: DN(2,4), G1(3), G2(3), E1(3), E2(3)
      REAL(DOUBLE), INTENT(OUT) :: DLOC(2,4)
      REAL(DOUBLE) :: A(2,2), AINV(2,2)

      A(1,1) = DOT_PRODUCT(G1, E1)
      A(1,2) = DOT_PRODUCT(G1, E2)
      A(2,1) = DOT_PRODUCT(G2, E1)
      A(2,2) = DOT_PRODUCT(G2, E2)
      CALL INV2(A, AINV)
      DLOC = MATMUL(AINV, DN)
      END SUBROUTINE LOCAL_DERIVS

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

      SUBROUTINE CROSS3 ( A, B, C )
      REAL(DOUBLE), INTENT(IN) :: A(3), B(3)
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
      REAL(DOUBLE), INTENT(IN) :: A(2,2)
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

      SUBROUTINE BUILD_SPECIAL_KGLOCAL_Q4EASANS ( DN_G, CO, SIG0, WT, JDET, KGOUT )
      REAL(DOUBLE), INTENT(IN) :: DN_G(2,4), CO(2,2), SIG0(2,2), WT, JDET
      REAL(DOUBLE), INTENT(INOUT) :: KGOUT(24,24)
      REAL(DOUBLE) :: DNDX(4), DNDY(4), KGVAL
      INTEGER(LONG) :: IA, IB, RR

      DO IA=1,4
         DNDX(IA) = DN_G(1,IA)*CO(1,1) + DN_G(2,IA)*CO(2,1)
         DNDY(IA) = DN_G(1,IA)*CO(1,2) + DN_G(2,IA)*CO(2,2)
      ENDDO

      DO IA=1,4
         DO IB=1,4
            KGVAL = WT*JDET*( DNDX(IA)*(SIG0(1,1)*DNDX(IB) + SIG0(1,2)*DNDY(IB)) + &
                             DNDY(IA)*(SIG0(2,1)*DNDX(IB) + SIG0(2,2)*DNDY(IB)) )
            IF (DEBUG(237) > 0) THEN
               KGVAL = KGVAL * 7.147442330726D-01
            ENDIF
            DO RR=1,3
               KGOUT(6*(IA-1)+RR,6*(IB-1)+RR) = KGOUT(6*(IA-1)+RR,6*(IB-1)+RR) + KGVAL
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE BUILD_SPECIAL_KGLOCAL_Q4EASANS

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

      END SUBROUTINE CQUADR_Q4EASANS
