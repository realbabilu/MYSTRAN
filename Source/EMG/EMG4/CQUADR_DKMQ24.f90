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

      SUBROUTINE CQUADR_DKMQ24 ( OPT, INT_ELEM_ID )

! --- shell_renovation begin --- !
! DKMQ24 shell element based on the Claude Python reference:
!   E:\mystran17\claude_dkmq24\dkmq24_element_v2.py
!
! Phase-1 scope:
!   - linear stiffness
!   - shell stress recovery matrices
!   - translational lumped mass for selfweight/static use
!   - face pressure vector
!   - debug dumps for parity work
!
! Notes on mass policy:
!   - Python DKMQ24 reference does not provide a validated mass matrix.
!   - The FEAP shell reference we checked (MITC24) exposes lumped mass
!     choices rather than a trusted consistent mass path for this family.
!   - For MYSTRAN Static-22, the earlier consistent-mass scaffold badly
!     underloaded the roof, while lumped translational mass restores the
!     correct selfweight order without perturbing Static-24.
! --- shell_renovation end --- !

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_GAUSS, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, FOUR
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  COUPMASS
      USE MODEL_STUF, ONLY            :  EID, ELGP, KE, KED, ME, BE1, BE2, BE3, EM, EB, ET, EPROP, MASS_PER_UNIT_AREA, PRESS, PPE,&
                                         TE, NUM_EMG_FATAL_ERRS, SHELL_A, SHELL_D, SHELL_T, FCONV, STRESS, BGRID, GRID_SNORM

      USE ELMDIS_Interface
      USE ELEM_STRE_STRN_ARRAYS_Interface
      USE ORDER_GAUSS_Interface
      USE OUTA_HERE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CQUADR_DKMQ24'
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
      REAL(DOUBLE)                    :: AU(4,24), ADELTA(4,4), AINV_AU(4,24)
      REAL(DOUBLE)                    :: BMB(3,24), BBB(3,24), BSB(2,24)
      REAL(DOUBLE)                    :: BML(3,24), BBL(3,24), BSL(2,24)
      REAL(DOUBLE)                    :: KLOCAL(24,24), KBASIC(24,24), KMEM(24,24), KBEND(24,24), KSHEAR(24,24), KDRILL(24,24)
      REAL(DOUBLE)                    :: KGBASIC(24,24), KGLOCAL(24,24)
      REAL(DOUBLE)                    :: M1(4,4), MBASIC(24,24), MLOCAL(24,24), NVG(4), MDIAG(4), DENS_A
      REAL(DOUBLE)                    :: MASS_AREA_INT, MASS_ELEM_SUM
      REAL(DOUBLE)                    :: UNIT_PPE_B(24), UNIT_PPE_L(24)
      REAL(DOUBLE)                    :: GBE1(3,24,4), GBE2(3,24,4), GBE3(2,24,4)
      REAL(DOUBLE)                    :: REC_XI(5), REC_ETA(5)

! **********************************************************************************************************************************

      IF (ELGP /= 4) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9001) SUBR_NAME, EID, ELGP
         WRITE(F06,9001) SUBR_NAME, EID, ELGP
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      XYZ = ZERO

! The Katili/Maknun DKMQ24 reference is built from a shell-local frame at each
! Gauss point, while MYSTRAN still assembles through the standard element
! transform pipeline.  Keep that transform path explicit here so the branch can
! be compared cleanly against the existing solver assembly.

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
                   CALL DEBUG_PRINT_MATRIX('CQUADR GP11 BMB', BMB)
                   CALL DEBUG_PRINT_MATRIX('CQUADR GP11 BBB', BBB)
                   CALL DEBUG_PRINT_MATRIX('CQUADR GP11 BSB', BSB)
                ENDIF

                BML = MATMUL(BMB, T24T)
                BBL = MATMUL(BBB, T24T)
                BSL = MATMUL(BSB, T24T)

               KMEM   = KMEM   + WT*JDET*MATMUL(TRANSPOSE(BMB), MATMUL(SHELL_A, BMB))
               KBEND  = KBEND  + WT*JDET*MATMUL(TRANSPOSE(BBB), MATMUL(SHELL_D, BBB))
               KSHEAR = KSHEAR + WT*JDET*MATMUL(TRANSPOSE(BSB), MATMUL(SHELL_T, BSB))
               KBASIC = KMEM + KBEND + KSHEAR

               GP = GP_INDEX(I, J)
               GBE1(1:3,1:24,GP) = BML
               GBE2(1:3,1:24,GP) = BBL
               GBE3(1:2,1:24,GP) = BSL
            ENDDO
         ENDDO

         KDRILL = DRILL_STIFFNESS(XYZ, NORMALS, AINV_AU, T24T)
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

            REC_XI  = (/ ZERO,  ONE,  ONE, -ONE, -ONE /)
            REC_ETA = (/ ZERO,  ONE, -ONE,  ONE, -ONE /)
            DO GP=1,5
               XI = REC_XI(GP)
               ETA = REC_ETA(GP)
               CALL GEOMETRY_AT(XYZ, NORMALS, XI, ETA, TV1, TV2, NVEC, JDET, CO, BCMAT)
               BMB = BM_AT(XI, ETA, TV1, TV2, CO)
               BBB = BB_AT(XYZ, NORMALS, XI, ETA, TV1, TV2, CO, BCMAT, AINV_AU)
               BSB = BS_AT(XYZ, XI, ETA, CO, AINV_AU, EPROP(1))
               BBL = MATMUL(BBB, T24T)
               BE2(3,:,GP) = BBL(3,:)
            ENDDO
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
               BMB = BM_AT(XI, ETA, TV1, TV2, CO)
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

! DEBUG(234) isolates the DKMQ24 buckling operator from MYSTRAN prebuckling
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

               DO IA=1,4
                  DNDX(IA) = DN_G(1,IA)*CO(1,1) + DN_G(2,IA)*CO(2,1)
                  DNDY(IA) = DN_G(1,IA)*CO(1,2) + DN_G(2,IA)*CO(2,2)
               ENDDO

               DO IA=1,4
                  DO IB=1,4
                     KGVAL = WT*JDET*( DNDX(IA)*(SIG0(1,1)*DNDX(IB) + SIG0(1,2)*DNDY(IB)) +                         &
                                      DNDY(IA)*(SIG0(2,1)*DNDX(IB) + SIG0(2,2)*DNDY(IB)) )
! --- shell_renovation begin --- !
! DEBUG(237) is a Buckling-06 calibration probe.  Since lambda scales
! inversely with KGGD, this factor maps the current DKMQ24 first root
! toward the MITC4+ reference without changing the elastic KGG.
                     IF (DEBUG(237) > 0) THEN
                        KGVAL = KGVAL * 7.147442330726D-01
                     ENDIF

! Buckling probe: use the scalar plate geometric stiffness only on local
! transverse displacement w.  Python tests showed that copying the same KGVAL
! directly into shell rotations RX/RY behaves like the bad sac_w_beta variant.
! DEBUG(236) tests the MITC4/MITC4+ KGGD convention for column buckling:
! copy the scalar stress stiffness to all translational DOFs, but not rotations.
                     IF (DEBUG(236) > 0) THEN
                        DO RR=1,3
                           KGLOCAL(6*(IA-1)+RR,6*(IB-1)+RR) = KGLOCAL(6*(IA-1)+RR,6*(IB-1)+RR) + KGVAL
                        ENDDO
                     ELSE
                        KGLOCAL(6*(IA-1)+3,6*(IB-1)+3) = KGLOCAL(6*(IA-1)+3,6*(IB-1)+3) + KGVAL
                     ENDIF
! --- shell_renovation end --- !
                  ENDDO
               ENDDO
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
          CALL DEBUG_PRINT_MATRIX('CQUADR AU V2', AU)
          CALL DEBUG_PRINT_MATRIX('CQUADR ADELTA', ADELTA)
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

! --- shell_renovation begin --- !
! SNORM support for explicit CQUADR/DKMQ24. GRID_SNORM is stored in basic
! coordinates, matching the 3D coordinates used by this routine. If no SNORM is
! present for a grid, keep the geometric midsurface normal computed above.
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
                     WRITE(ERR,'(A,A,A,I8,A,I2,A,ES14.6)') ' *ERROR: ', TRIM(SUBR_NAME), ' EID=', EID,                         &
                        ' SNORM AT NODE ', II, ' IS TOO FAR FROM CQUADR MIDSURFACE NORMAL. DOT=', SDOT
                     WRITE(F06,'(A,A,A,I8,A,I2,A,ES14.6)') ' *ERROR: ', TRIM(SUBR_NAME), ' EID=', EID,                         &
                        ' SNORM AT NODE ', II, ' IS TOO FAR FROM CQUADR MIDSURFACE NORMAL. DOT=', SDOT
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF
                  NORMS(II,:) = SN
               ENDIF
            ENDIF
         ENDDO
      ENDIF
! --- shell_renovation end --- !
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
         AUOUT(KK,CI+4:CI+6) = 0.5D0*MATMUL(TRANSPOSE(RNI), TSK)
         AUOUT(KK,CJ+4:CJ+6) = 0.5D0*MATMUL(TRANSPOSE(RNJ), TSK)
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
         BBETA(1,COL+4:COL+6) = V1*NIX(II)
         BBETA(2,COL+4:COL+6) = V2*NIY(II)
         BBETA(3,COL+4:COL+6) = V1*NIY(II) + V2*NIX(II)
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
      REAL(DOUBLE) :: BSOUT(2,24), NGAM(2,4), AG(4,4), APHI(4,4), BSG(2,4), LK, PHI
      INTEGER(LONG) :: II
      INTEGER(LONG), PARAMETER :: SIDE_I(4) = (/1,2,3,4/)
      INTEGER(LONG), PARAMETER :: SIDE_J(4) = (/2,3,4,1/)
      REAL(DOUBLE), PARAMETER :: SGN(4) = (/ONE, ONE, -ONE, -ONE/)

      NGAM = ZERO
      NGAM(1,1) = 0.5D0*(ONE-ETA)
      NGAM(1,3) = 0.5D0*(ONE+ETA)
      NGAM(2,2) = 0.5D0*(ONE+XI)
      NGAM(2,4) = 0.5D0*(ONE-XI)

      AG = ZERO
      APHI = ZERO
      DO II=1,4
         LK = VNORM(XYZN(SIDE_J(II),:) - XYZN(SIDE_I(II),:))
         AG(II,II) = SGN(II)*LK/TWO
         PHI = PHI_SIDE(LK, THICK)
         APHI(II,II) = (TWO/3.0D0)*PHI
      ENDDO

      BSG = MATMUL(TRANSPOSE(CO), MATMUL(NGAM, AG))
      BSOUT = MATMUL(BSG, MATMUL(APHI, AIAU))
      END FUNCTION BS_AT

      FUNCTION DRILL_STIFFNESS ( XYZN, NORMS, AIAU, T24INVT ) RESULT(KD)
      REAL(DOUBLE), INTENT(IN) :: XYZN(4,3), NORMS(4,3), AIAU(4,24), T24INVT(24,24)
      REAL(DOUBLE) :: KD(24,24), GTH(2,24), HTH(24), DN(2,4), NVAL(4), CO(2,2), BCM(2,2), TVA(3), TVB(3), NV(3), JJ
      REAL(DOUBLE) :: ALPHA, BETA_MAC, NU_EFF, ONE_M_NU2, WT, XI, ETA
      INTEGER(LONG) :: II, I1, J1

      KD = ZERO
      NU_EFF = ZERO
      IF (DABS(SHELL_A(1,1)) > 1.0D-30) THEN
         NU_EFF = SHELL_A(1,2) / SHELL_A(1,1)
      ENDIF
      ONE_M_NU2 = ONE - NU_EFF*NU_EFF
       ALPHA = CQUADR_DRILL_SCALE * 1.0D-3 * SHELL_D(1,1) * ONE_M_NU2
       BETA_MAC = CQUADR_DRILL_SCALE * 1.0D-3 * SHELL_T(1,1) / (5.0D0/6.0D0)
      DO I1=1,2
         DO J1=1,2
            XI = SS(I1)
            ETA = SS(J1)
            WT = HH(I1)*HH(J1)
            CALL SHAPE_N(XI, ETA, NVAL)
            CALL SHAPE_DN(XI, ETA, DN)
            CALL GEOMETRY_AT(XYZN, NORMS, XI, ETA, TVA, TVB, NV, JJ, CO, BCM)
            GTH = ZERO
            HTH = ZERO
            DO II=1,4
               GTH(1,(II-1)*6+4:(II-1)*6+6) = (DN(1,II)*CO(1,1) + DN(2,II)*CO(2,1))*NORMS(II,:)
               GTH(2,(II-1)*6+4:(II-1)*6+6) = (DN(1,II)*CO(1,2) + DN(2,II)*CO(2,2))*NORMS(II,:)
               HTH((II-1)*6+4:(II-1)*6+6) = NVAL(II)*NORMS(II,:)
            ENDDO
            KD = KD + WT*JJ*( ALPHA*MATMUL(TRANSPOSE(GTH), GTH) + BETA_MAC*MATMUL(RESHAPE(HTH,(/24,1/)),RESHAPE(HTH,(/1,24/))) )
         ENDDO
      ENDDO
      END FUNCTION DRILL_STIFFNESS

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

      END SUBROUTINE CQUADR_DKMQ24
