! ##################################################################################################################################
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

      SUBROUTINE EIG_LANCZOS_FEAST

! FEAST backend for Lanczos entry point.
! Native path currently targets generalized symmetric real sparse problem:
!     KLL * x = lambda * MLL * x
! in modal solutions. Fallback remains MGIV for unsupported configurations.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFL, NUM_EIGENS, NVEC, NTERM_KLL, NTERM_MLL, WARN_ERR, SOL_NAME
      USE CONSTANTS_1, ONLY           :  ZERO, PI
      USE PARAMS, ONLY                :  EPSIL, SUPINFO
      USE MODEL_STUF, ONLY            :  EIG_FRQ1, EIG_FRQ2, EIG_METH, EIG_MSGLVL, EIG_N2
      USE SPARSE_MATRICES, ONLY       :  I_KLL, J_KLL, KLL, I_MLL, J_MLL, MLL
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL, EIGEN_VEC, MODE_NUM

      USE EIG_LANCZOS_FEAST_USE_IFs
      USE LINK_MESSAGE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EIG_LANCZOS_FEAST'
      CHARACTER(8*BYTE)               :: EIG_METH_SAVE = ' '
      CHARACTER(1*BYTE)               :: UPLO

      INTEGER(LONG), PARAMETER        :: MAX_FEAST_TRY = 6
      INTEGER(LONG), PARAMETER        :: FPM_SIZE = 64
      INTEGER(LONG)                   :: I, J, K, IT
      INTEGER(LONG)                   :: M0
      INTEGER(LONG)                   :: MODE
      INTEGER(LONG)                   :: INFO
      INTEGER(LONG)                   :: LOOP
      INTEGER(LONG)                   :: NEV_TARGET
      INTEGER(LONG)                   :: NKEEP
      INTEGER(LONG)                   :: IMIN
      INTEGER(LONG)                   :: TMP_IDX
      INTEGER(LONG)                   :: MIN_NEEDED
      INTEGER(LONG)                   :: NMASS_DIAG_POS
      INTEGER(LONG)                   :: FPM(FPM_SIZE)
      INTEGER(LONG), ALLOCATABLE      :: PERM(:)

      REAL(DOUBLE)                    :: EPS1
      REAL(DOUBLE)                    :: EPSOUT
      REAL(DOUBLE)                    :: EMIN
      REAL(DOUBLE)                    :: EMAX
      REAL(DOUBLE)                    :: MAX_RATIO
      REAL(DOUBLE)                    :: KDIAG
      REAL(DOUBLE)                    :: MDIAG
      REAL(DOUBLE), ALLOCATABLE       :: LAMBDA(:)
      REAL(DOUBLE), ALLOCATABLE       :: RES(:)
      REAL(DOUBLE), ALLOCATABLE       :: Q(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: DIAG_K(:)
      REAL(DOUBLE), ALLOCATABLE       :: DIAG_M(:)

#ifdef MYSTRAN_HAVE_EXTERNAL_FEAST
      INTERFACE
         SUBROUTINE FEASTINIT(FPM)
            USE PENTIUM_II_KIND, ONLY : LONG
            INTEGER(LONG) :: FPM(*)
         END SUBROUTINE FEASTINIT
         SUBROUTINE DFEAST_SCSRGV(UPLO,N,SA,ISA,JSA,SB,ISB,JSB,FPM,EPSOUT,LOOP,EMIN,EMAX,M0,LAMBDA,Q,MODE,RES,INFO)
            USE PENTIUM_II_KIND, ONLY : BYTE, LONG, DOUBLE
            CHARACTER(1*BYTE) :: UPLO
            INTEGER(LONG) :: N, ISA(*), JSA(*), ISB(*), JSB(*), FPM(*), LOOP, M0, MODE, INFO
            REAL(DOUBLE) :: SA(*), SB(*), EPSOUT, EMIN, EMAX, LAMBDA(*), Q(N,*), RES(*)
         END SUBROUTINE DFEAST_SCSRGV
      END INTERFACE
#endif

      EPS1 = EPSIL(1)

#ifdef MYSTRAN_HAVE_EXTERNAL_FEAST

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - FEAST FALLBACK (BUCKLING NOT YET WIRED)')
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,4917)
         IF (SUPINFO == 'N') WRITE(F06,4917)
         GOTO 900
      ENDIF

      IF (EIG_FRQ2 <= EPS1) THEN
         CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - FEAST FALLBACK (NO EIGENVALUE INTERVAL)')
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,4915)
         IF (SUPINFO == 'N') WRITE(F06,4915)
         GOTO 900
      ENDIF

      CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - FEAST NATIVE (GENERALIZED SYMMETRIC SPARSE)')

      WARN_ERR = WARN_ERR + 1
      WRITE(ERR,4918)
      IF (SUPINFO == 'N') WRITE(F06,4918)

      NEV_TARGET = EIG_N2
      IF (NEV_TARGET < 1) NEV_TARGET = 1
      IF (NEV_TARGET > NDOFL) NEV_TARGET = NDOFL

      EMIN = ZERO
      IF (EIG_FRQ1 > EPS1) EMIN = (2.0D0*PI*EIG_FRQ1)**2
      EMAX = (2.0D0*PI*EIG_FRQ2)**2

      ALLOCATE(DIAG_K(NDOFL), DIAG_M(NDOFL))
      DO I=1,NDOFL
         DIAG_K(I) = ZERO
         DIAG_M(I) = ZERO
      ENDDO

      DO I=1,NDOFL
         DO K=I_KLL(I),I_KLL(I+1)-1
            IF (J_KLL(K) == I) THEN
               DIAG_K(I) = KLL(K)
               EXIT
            ENDIF
         ENDDO
         DO K=I_MLL(I),I_MLL(I+1)-1
            IF (J_MLL(K) == I) THEN
               DIAG_M(I) = MLL(K)
               EXIT
            ENDIF
         ENDDO
      ENDDO

      MAX_RATIO = ZERO
      NMASS_DIAG_POS = 0
      DO I=1,NDOFL
         KDIAG = DIAG_K(I)
         MDIAG = DIAG_M(I)
         IF (MDIAG > EPS1) THEN
            NMASS_DIAG_POS = NMASS_DIAG_POS + 1
            IF (KDIAG/MDIAG > MAX_RATIO) MAX_RATIO = KDIAG/MDIAG
         ENDIF
      ENDDO

      IF (NMASS_DIAG_POS < NDOFL) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,4919) NMASS_DIAG_POS, NDOFL
         IF (SUPINFO == 'N') WRITE(F06,4919) NMASS_DIAG_POS, NDOFL
         DEALLOCATE(DIAG_K, DIAG_M)
         GOTO 900
      ENDIF

      IF (MAX_RATIO > ZERO) THEN
         IF (EMAX < 1.5D0*MAX_RATIO) EMAX = 1.5D0*MAX_RATIO
      ENDIF

      DEALLOCATE(DIAG_K, DIAG_M)

      M0 = MAX(2*NEV_TARGET, NEV_TARGET + 8)
      IF (M0 > NDOFL) M0 = NDOFL
      IF (M0 < 8) M0 = MIN(8, NDOFL)

      INFO = -999
      MODE = 0

      DO IT=1,MAX_FEAST_TRY

         IF (ALLOCATED(LAMBDA)) DEALLOCATE(LAMBDA)
         IF (ALLOCATED(RES))    DEALLOCATE(RES)
         IF (ALLOCATED(Q))      DEALLOCATE(Q)
         ALLOCATE(LAMBDA(M0), RES(M0), Q(NDOFL,M0))

         CALL FEASTINIT(FPM)

         IF (EIG_MSGLVL <= 0) THEN
            FPM(1) = 0
         ELSE
            FPM(1) = 1
         ENDIF
         FPM(3) = 12
         FPM(7) = 7

         UPLO = 'U'
         CALL DFEAST_SCSRGV(UPLO, NDOFL, KLL, I_KLL, J_KLL, &
                            MLL, I_MLL, J_MLL, FPM, EPSOUT, LOOP, &
                            EMIN, EMAX, M0, LAMBDA, Q, MODE, RES, INFO)

         MIN_NEEDED = MIN(NEV_TARGET, M0)
         IF ((INFO == 0) .AND. (MODE >= MIN_NEEDED)) THEN
            EXIT
         ENDIF

         IF ((MODE >= M0) .OR. (INFO == 3)) THEN
            M0 = MIN(NDOFL, MAX(2*M0, M0 + 8))
         ELSE
            EMAX = 4.0D0*EMAX
         ENDIF
      ENDDO

      IF ((INFO /= 0) .OR. (MODE <= 0)) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,4916) INFO
         IF (SUPINFO == 'N') WRITE(F06,4916) INFO
         IF (ALLOCATED(LAMBDA)) DEALLOCATE(LAMBDA)
         IF (ALLOCATED(RES))    DEALLOCATE(RES)
         IF (ALLOCATED(Q))      DEALLOCATE(Q)
         GOTO 900
      ENDIF

      NKEEP = MIN(MODE, NEV_TARGET)
      IF (NKEEP < 1) NKEEP = 1

      CALL ALLOCATE_EIGEN1_MAT('EIGEN_VEC', NDOFL, NKEEP, SUBR_NAME)
      CALL ALLOCATE_EIGEN1_MAT('MODE_NUM' , NDOFL, 1, SUBR_NAME)
      CALL ALLOCATE_EIGEN1_MAT('EIGEN_VAL', NDOFL, 1, SUBR_NAME)

      ALLOCATE(PERM(MODE))
      DO I=1,MODE
         PERM(I) = I
      ENDDO

      DO I=1,MODE-1
         IMIN = I
         DO J=I+1,MODE
            IF (LAMBDA(PERM(J)) < LAMBDA(PERM(IMIN))) IMIN = J
         ENDDO
         IF (IMIN /= I) THEN
            TMP_IDX = PERM(I)
            PERM(I) = PERM(IMIN)
            PERM(IMIN) = TMP_IDX
         ENDIF
      ENDDO

      DO I=1,NKEEP
         EIGEN_VAL(I) = LAMBDA(PERM(I))
         MODE_NUM(I)  = I
         DO J=1,NDOFL
            EIGEN_VEC(J,I) = Q(J,PERM(I))
         ENDDO
      ENDDO

      NUM_EIGENS = NKEEP
      NVEC       = NKEEP

      DEALLOCATE(PERM)
      DEALLOCATE(LAMBDA, RES, Q)
      RETURN

#else
      CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - FEAST SURROGATE (MGIV CORE)')

      WARN_ERR = WARN_ERR + 1
      WRITE(ERR,4911)
      IF (SUPINFO == 'N') THEN
         WRITE(F06,4911)
      ENDIF
#endif

 900  CONTINUE
      EIG_METH_SAVE = EIG_METH
      EIG_METH      = 'MGIV    '
      CALL EIG_GIV_MGIV
      EIG_METH = EIG_METH_SAVE
      RETURN

 4911 FORMAT(' *WARNING 4911: FEAST EXTERNAL BACKEND NOT LINKED. USING MGIV SURROGATE (NON-ARPACK).')
 4915 FORMAT(' *WARNING 4915: FEAST NATIVE REQUIRES EIGRL FREQUENCY INTERVAL (V1/V2). USING MGIV SURROGATE.')
 4916 FORMAT(' *WARNING 4916: FEAST NATIVE FAILED (INFO=',I8,'). USING MGIV SURROGATE.')
 4917 FORMAT(' *WARNING 4917: FEAST NATIVE FOR BUCKLING (KLL,KLLD) NOT WIRED YET. USING MGIV SURROGATE.')
 4918 FORMAT(' *WARNING 4918: FEAST EXTERNAL BACKEND LINKED. RUNNING NATIVE FEAST PATH FOR MODES.')
 4919 FORMAT(' *WARNING 4919: FEAST NATIVE REQUIRES FULL-RANK POSITIVE MASS DIAGONAL. FOUND ',I8,' OF ',I8, &
             ' POSITIVE MASS DIAGONALS IN MLL. USING MGIV SURROGATE.')

      END SUBROUTINE EIG_LANCZOS_FEAST
