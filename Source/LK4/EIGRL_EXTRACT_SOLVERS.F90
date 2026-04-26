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

! --- chase_feast_add --- begin !
      MODULE EIGRL_EXTRACT_SOLVERS

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFL, NUM_EIGENS, NVEC, SOL_NAME, WARN_ERR
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL, SPARSTOR, SUPINFO
      USE MODEL_STUF, ONLY            :  EIG_CHASE_DEG, EIG_CHASE_MAX_ITER, EIG_CHASE_NEX, EIG_CHASE_TOL,                  &
                                         EIG_DENSE_NEX, EIG_EXTRACT_METHOD, EIG_FEAST_M0, EIG_FEAST_MAX_LOOP,                &
                                         EIG_FEAST_N_CONTOUR, EIG_FEAST_SEARCH_SCALE, EIG_FEAST_TOL_DIGITS, EIG_FRQ1,       &
                                         EIG_FRQ2, EIG_N2, EIG_SUBSPACE_MAX_ITER, EIG_SUBSPACE_NSUB, EIG_SUBSPACE_TOL
      USE SPARSE_MATRICES, ONLY       :  I_KLL, I_MLL, J_KLL, J_MLL, KLL, MLL
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL, EIGEN_VEC, MODE_NUM

      USE ALLOCATE_EIGEN1_MAT_Interface
      USE EIG_LANCZOS_ARPACK_Interface
      USE LINK_MESSAGE_Interface
#ifdef MYSTRAN_HAVE_EXTERNAL_CHASE
      USE chase_diag
#endif

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: EIG_LANCZOS_CHASE, EIG_LANCZOS_FEAST, EIG_LANCZOS_SUBSPACE, EIG_LANCZOS_DENSE

      INTERFACE
         SUBROUTINE DPOTRF(UPLO, N, A, LDA, INFO)
            CHARACTER(1), INTENT(IN) :: UPLO
            INTEGER, INTENT(IN) :: N, LDA
            DOUBLE PRECISION, INTENT(INOUT) :: A(LDA,*)
            INTEGER, INTENT(OUT) :: INFO
         END SUBROUTINE DPOTRF

         SUBROUTINE DPOTRS(UPLO, N, NRHS, A, LDA, B, LDB, INFO)
            CHARACTER(1), INTENT(IN) :: UPLO
            INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
            DOUBLE PRECISION, INTENT(IN) :: A(LDA,*)
            DOUBLE PRECISION, INTENT(INOUT) :: B(LDB,*)
            INTEGER, INTENT(OUT) :: INFO
         END SUBROUTINE DPOTRS

         SUBROUTINE DSYEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
            CHARACTER(1), INTENT(IN) :: JOBZ, UPLO
            INTEGER, INTENT(IN) :: N, LDA, LWORK
            DOUBLE PRECISION, INTENT(INOUT) :: A(LDA,*), WORK(*)
            DOUBLE PRECISION, INTENT(OUT) :: W(*)
            INTEGER, INTENT(OUT) :: INFO
         END SUBROUTINE DSYEV
      END INTERFACE

      CONTAINS

!***********************************************************************************************************************************
      SUBROUTINE EIG_LANCZOS_DENSE

      CALL SOLVE_CONDENSED_MODAL('DENSE')

      END SUBROUTINE EIG_LANCZOS_DENSE

!***********************************************************************************************************************************
      SUBROUTINE EIG_LANCZOS_SUBSPACE

      CALL SOLVE_CONDENSED_MODAL('SUBSP')

      END SUBROUTINE EIG_LANCZOS_SUBSPACE

!***********************************************************************************************************************************
      SUBROUTINE EIG_LANCZOS_CHASE

      CALL SOLVE_CONDENSED_MODAL('CHASE')

      END SUBROUTINE EIG_LANCZOS_CHASE

!***********************************************************************************************************************************
      SUBROUTINE EIG_LANCZOS_FEAST

      CALL SOLVE_CONDENSED_MODAL('FEAST')

      END SUBROUTINE EIG_LANCZOS_FEAST

!***********************************************************************************************************************************
      SUBROUTINE SOLVE_CONDENSED_MODAL ( METHOD )

      CHARACTER(LEN=*), INTENT(IN)    :: METHOD

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: INFO
      INTEGER(LONG)                   :: KEEP_COUNT
      INTEGER(LONG)                   :: NACTIVE
      INTEGER(LONG)                   :: NZERO

      INTEGER(LONG), ALLOCATABLE      :: ACTIVE_MAP(:)
      INTEGER(LONG), ALLOCATABLE      :: ACTIVE_POS(:)
      INTEGER(LONG), ALLOCATABLE      :: KEEP_IDX(:)
      INTEGER(LONG), ALLOCATABLE      :: ZERO_MAP(:)
      INTEGER(LONG), ALLOCATABLE      :: ZERO_POS(:)
      INTEGER(LONG), ALLOCATABLE      :: KZA_COLPTR(:)
      INTEGER(LONG), ALLOCATABLE      :: KZA_ROW(:)

      REAL(DOUBLE)                    :: EMAX
      REAL(DOUBLE)                    :: EMIN
      REAL(DOUBLE), ALLOCATABLE       :: ASTD(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: EVAL_ALL(:)
      REAL(DOUBLE), ALLOCATABLE       :: EVEC_ACTIVE(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: EVEC_FULL(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: EVEC_STD(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: KAA(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: KCOND(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: KZZ(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: MASS_A(:)
      REAL(DOUBLE), ALLOCATABLE       :: MDIAG(:)
      REAL(DOUBLE), ALLOCATABLE       :: KZA_VAL(:)

      SUBR_NAME = 'EIGRL_EXTRACT_SOLVERS'

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CALL WRITE_FALLBACK_WARNING(METHOD, 4941, 'CONDENSED EXTRACT METHODS ARE CURRENTLY MODES-ONLY. USING ARPACK LANCZOS.')
         CALL EIG_LANCZOS_ARPACK
         RETURN
      ENDIF

      CALL BUILD_CONDENSED_MODAL_CORE(NACTIVE, NZERO, ACTIVE_POS, ZERO_POS, ACTIVE_MAP, ZERO_MAP, MASS_A, KAA, KZZ, KZA_COLPTR,   &
                                      KZA_ROW, KZA_VAL, MDIAG, INFO)
      IF (INFO /= 0) THEN
         CALL WRITE_FALLBACK_WARNING(METHOD, 4942, 'FAILED TO BUILD CONDENSED MODAL CORE. USING ARPACK LANCZOS.')
         CALL EIG_LANCZOS_ARPACK
         RETURN
      ENDIF

      IF (NACTIVE <= 0) THEN
         CALL WRITE_FALLBACK_WARNING(METHOD, 4943, 'NO POSITIVE-MASS ACTIVE DOF FOUND. USING ARPACK LANCZOS.')
         CALL EIG_LANCZOS_ARPACK
         RETURN
      ENDIF

      CALL FORM_CONDENSED_STIFFNESS(NACTIVE, NZERO, KAA, KZZ, KZA_COLPTR, KZA_ROW, KZA_VAL, KCOND, INFO)
      IF (INFO /= 0) THEN
         CALL WRITE_FALLBACK_WARNING(METHOD, 4944, 'FAILED TO CONDENSE ZERO-MASS DOF. USING ARPACK LANCZOS.')
         CALL EIG_LANCZOS_ARPACK
         RETURN
      ENDIF

      CALL BUILD_STANDARDIZED_OPERATOR(KCOND, MASS_A, ASTD)

      EMIN = ZERO
      EMAX = -ONE()
      IF (EIG_FRQ1 > EPSIL(1)) EMIN = FREQ_TO_LAMBDA(EIG_FRQ1)
      IF (EIG_FRQ2 > EPSIL(1)) EMAX = FREQ_TO_LAMBDA(EIG_FRQ2)

      SELECT CASE (METHOD(1:MIN(LEN(METHOD),5)))
      CASE ('DENSE')
         CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - DENSE REFERENCE (CONDENSED DSYEV)')
         CALL RUN_DENSE_BACKEND(ASTD, EVAL_ALL, EVEC_STD, INFO)
         IF (INFO /= 0) THEN
            CALL WRITE_FALLBACK_WARNING(METHOD, 4945, 'DENSE DSYEV FAILED. USING ARPACK LANCZOS.')
            CALL EIG_LANCZOS_ARPACK
            RETURN
         ENDIF
         CALL STANDARDIZED_TO_ACTIVE(MASS_A, EVEC_STD, EVEC_ACTIVE)

      CASE ('SUBSP')
         IF (EIG_FRQ1 > EPSIL(1)) THEN
            CALL WRITE_FALLBACK_WARNING(METHOD, 4946, 'SUBSPACE V1/V2 INTERVAL DISPATCH FALLS BACK TO DENSE IN V1.')
            CALL RUN_DENSE_BACKEND(ASTD, EVAL_ALL, EVEC_STD, INFO)
            IF (INFO /= 0) THEN
               CALL EIG_LANCZOS_ARPACK
               RETURN
            ENDIF
            CALL STANDARDIZED_TO_ACTIVE(MASS_A, EVEC_STD, EVEC_ACTIVE)
         ELSE
            CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - SUBSPACE (CONDENSED STANDARDIZED)')
            CALL RUN_SUBSPACE_BACKEND(ASTD, MAX(1,EIG_N2), MAX(EIG_SUBSPACE_NSUB, EIG_N2), EIG_SUBSPACE_TOL,                      &
                                      EIG_SUBSPACE_MAX_ITER, EVAL_ALL, EVEC_STD, INFO)
            IF (INFO /= 0) THEN
               CALL WRITE_FALLBACK_WARNING(METHOD, 4947, 'SUBSPACE ITERATION FAILED. USING DENSE REFERENCE FALLBACK.')
               CALL RUN_DENSE_BACKEND(ASTD, EVAL_ALL, EVEC_STD, INFO)
               IF (INFO /= 0) THEN
                  CALL EIG_LANCZOS_ARPACK
                  RETURN
               ENDIF
            ENDIF
            CALL STANDARDIZED_TO_ACTIVE(MASS_A, EVEC_STD, EVEC_ACTIVE)
         ENDIF

      CASE ('CHASE')
#ifdef MYSTRAN_HAVE_EXTERNAL_CHASE
         IF (EIG_FRQ1 > EPSIL(1)) THEN
            CALL WRITE_FALLBACK_WARNING(METHOD, 4948, 'CHASE INTERVAL SEARCH FALLS BACK TO DENSE IN V1.')
            CALL RUN_DENSE_BACKEND(ASTD, EVAL_ALL, EVEC_STD, INFO)
            IF (INFO /= 0) THEN
               CALL EIG_LANCZOS_ARPACK
               RETURN
            ENDIF
         ELSE
            CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - CHASE NATIVE (CONDENSED STANDARDIZED)')
            CALL RUN_CHASE_BACKEND(ASTD, MAX(1,EIG_N2), MAX(1,EIG_CHASE_NEX), EIG_CHASE_TOL, EIG_CHASE_MAX_ITER,                  &
                                   EVAL_ALL, EVEC_STD, INFO)
            IF (INFO /= 0) THEN
               CALL WRITE_FALLBACK_WARNING(METHOD, 4949, 'CHASE NATIVE FAILED. USING DENSE REFERENCE FALLBACK.')
               CALL RUN_DENSE_BACKEND(ASTD, EVAL_ALL, EVEC_STD, INFO)
               IF (INFO /= 0) THEN
                  CALL EIG_LANCZOS_ARPACK
                  RETURN
               ENDIF
            ENDIF
         ENDIF
#else
         CALL WRITE_FALLBACK_WARNING(METHOD, 4950, 'CHASE EXTERNAL BACKEND NOT LINKED. USING DENSE REFERENCE FALLBACK.')
         CALL RUN_DENSE_BACKEND(ASTD, EVAL_ALL, EVEC_STD, INFO)
         IF (INFO /= 0) THEN
            CALL EIG_LANCZOS_ARPACK
            RETURN
         ENDIF
#endif
         CALL STANDARDIZED_TO_ACTIVE(MASS_A, EVEC_STD, EVEC_ACTIVE)

      CASE ('FEAST')
#ifdef MYSTRAN_HAVE_EXTERNAL_FEAST
         IF (EIG_FRQ2 <= EPSIL(1)) THEN
            CALL WRITE_FALLBACK_WARNING(METHOD, 4951, 'FEAST REQUIRES V2 INTERVAL UPPER BOUND. USING DENSE REFERENCE FALLBACK.')
            CALL RUN_DENSE_BACKEND(ASTD, EVAL_ALL, EVEC_STD, INFO)
            IF (INFO /= 0) THEN
               CALL EIG_LANCZOS_ARPACK
               RETURN
            ENDIF
            CALL STANDARDIZED_TO_ACTIVE(MASS_A, EVEC_STD, EVEC_ACTIVE)
         ELSE
            CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - FEAST NATIVE (CONDENSED GENERALIZED)')
            CALL RUN_FEAST_BACKEND(KCOND, MASS_A, EMIN, MAX(EMIN + EPSIL(1), EMAX*MAX(EIG_FEAST_SEARCH_SCALE,ONE())),              &
                                   EVAL_ALL, EVEC_ACTIVE, INFO)
            IF (INFO /= 0) THEN
               CALL WRITE_FALLBACK_WARNING(METHOD, 4952, 'FEAST NATIVE FAILED. USING DENSE REFERENCE FALLBACK.')
               CALL RUN_DENSE_BACKEND(ASTD, EVAL_ALL, EVEC_STD, INFO)
               IF (INFO /= 0) THEN
                  CALL EIG_LANCZOS_ARPACK
                  RETURN
               ENDIF
               CALL STANDARDIZED_TO_ACTIVE(MASS_A, EVEC_STD, EVEC_ACTIVE)
            ENDIF
         ENDIF
#else
         CALL WRITE_FALLBACK_WARNING(METHOD, 4953, 'FEAST EXTERNAL BACKEND NOT LINKED. USING DENSE REFERENCE FALLBACK.')
         CALL RUN_DENSE_BACKEND(ASTD, EVAL_ALL, EVEC_STD, INFO)
         IF (INFO /= 0) THEN
            CALL EIG_LANCZOS_ARPACK
            RETURN
         ENDIF
         CALL STANDARDIZED_TO_ACTIVE(MASS_A, EVEC_STD, EVEC_ACTIVE)
#endif

      CASE DEFAULT
         CALL WRITE_FALLBACK_WARNING(METHOD, 4954, 'UNKNOWN EXTRACT METHOD. USING ARPACK LANCZOS.')
         CALL EIG_LANCZOS_ARPACK
         RETURN
      END SELECT

      CALL SELECT_REQUESTED_MODES(EVAL_ALL, EMIN, EMAX, KEEP_IDX, KEEP_COUNT)
      IF (KEEP_COUNT <= 0) THEN
         CALL WRITE_FALLBACK_WARNING(METHOD, 4955, 'NO MODES SATISFIED THE REQUESTED DENSE/CONDENSED FILTER. USING ARPACK LANCZOS.')
         CALL EIG_LANCZOS_ARPACK
         RETURN
      ENDIF

      CALL EXPAND_FULL_EIGENVECTORS(KEEP_IDX, KEEP_COUNT, ACTIVE_POS, ZERO_POS, MASS_A, EVEC_ACTIVE, KZZ, KZA_COLPTR, KZA_ROW,     &
                                    KZA_VAL, EVEC_FULL, INFO)
      IF (INFO /= 0) THEN
         CALL WRITE_FALLBACK_WARNING(METHOD, 4956, 'FAILED TO EXPAND CONDENSED EIGENVECTORS. USING ARPACK LANCZOS.')
         CALL EIG_LANCZOS_ARPACK
         RETURN
      ENDIF

      CALL ALLOCATE_EIGEN1_MAT('EIGEN_VEC', NDOFL, KEEP_COUNT, SUBR_NAME)
      CALL ALLOCATE_EIGEN1_MAT('MODE_NUM' , NDOFL, 1,          SUBR_NAME)
      CALL ALLOCATE_EIGEN1_MAT('EIGEN_VAL', NDOFL, 1,          SUBR_NAME)

      DO I=1,KEEP_COUNT
         EIGEN_VAL(I) = EVAL_ALL(KEEP_IDX(I))
         MODE_NUM(I)  = I
         EIGEN_VEC(1:NDOFL,I) = EVEC_FULL(1:NDOFL,I)
      ENDDO

      NUM_EIGENS = KEEP_COUNT
      NVEC       = KEEP_COUNT

      END SUBROUTINE SOLVE_CONDENSED_MODAL

!***********************************************************************************************************************************
      SUBROUTINE BUILD_CONDENSED_MODAL_CORE ( NACTIVE, NZERO, ACTIVE_POS, ZERO_POS, ACTIVE_MAP, ZERO_MAP, MASS_A, KAA, KZZ,        &
                                              KZA_COLPTR, KZA_ROW, KZA_VAL, MDIAG, INFO )

      INTEGER(LONG), INTENT(OUT)      :: INFO
      INTEGER(LONG), INTENT(OUT)      :: NACTIVE
      INTEGER(LONG), INTENT(OUT)      :: NZERO
      INTEGER(LONG), ALLOCATABLE      :: ACTIVE_MAP(:)
      INTEGER(LONG), ALLOCATABLE      :: ACTIVE_POS(:)
      INTEGER(LONG), ALLOCATABLE      :: KZA_COLPTR(:)
      INTEGER(LONG), ALLOCATABLE      :: KZA_ROW(:)
      INTEGER(LONG), ALLOCATABLE      :: ZERO_MAP(:)
      INTEGER(LONG), ALLOCATABLE      :: ZERO_POS(:)
      REAL(DOUBLE), ALLOCATABLE       :: KAA(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: KZA_VAL(:)
      REAL(DOUBLE), ALLOCATABLE       :: KZZ(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: MASS_A(:)
      REAL(DOUBLE), ALLOCATABLE       :: MDIAG(:)

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: J
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: P
      INTEGER(LONG), ALLOCATABLE      :: KZA_COUNT(:)
      INTEGER(LONG), ALLOCATABLE      :: KZA_NEXT(:)
      REAL(DOUBLE)                    :: VAL

      INFO = 0

      ALLOCATE(MDIAG(NDOFL), ACTIVE_MAP(NDOFL), ZERO_MAP(NDOFL))
      MDIAG      = ZERO
      ACTIVE_MAP = 0
      ZERO_MAP   = 0

      DO I=1,NDOFL
         DO K=I_MLL(I),I_MLL(I+1)-1
            IF (J_MLL(K) == I) THEN
               MDIAG(I) = MLL(K)
               EXIT
            ENDIF
         ENDDO
      ENDDO

      NACTIVE = 0
      NZERO   = 0
      DO I=1,NDOFL
         IF (MDIAG(I) > EPSIL(1)) THEN
            NACTIVE = NACTIVE + 1
            ACTIVE_MAP(I) = NACTIVE
         ELSE
            NZERO = NZERO + 1
            ZERO_MAP(I) = NZERO
         ENDIF
      ENDDO

      IF (NACTIVE <= 0) THEN
         INFO = 1
         RETURN
      ENDIF

      ALLOCATE(ACTIVE_POS(NACTIVE), MASS_A(NACTIVE), KAA(NACTIVE,NACTIVE), KZA_COUNT(MAX(1,NACTIVE)))
      KAA       = ZERO
      KZA_COUNT = 0
      DO I=1,NDOFL
         IF (ACTIVE_MAP(I) > 0) THEN
            ACTIVE_POS(ACTIVE_MAP(I)) = I
            MASS_A(ACTIVE_MAP(I)) = MDIAG(I)
         ENDIF
      ENDDO

      IF (NZERO > 0) THEN
         ALLOCATE(ZERO_POS(NZERO), KZZ(NZERO,NZERO))
         KZZ = ZERO
         DO I=1,NDOFL
            IF (ZERO_MAP(I) > 0) ZERO_POS(ZERO_MAP(I)) = I
         ENDDO
      ELSE
         ALLOCATE(ZERO_POS(1), KZZ(1,1))
         ZERO_POS = 0
         KZZ = ZERO
      ENDIF

      DO I=1,NDOFL
         DO K=I_KLL(I),I_KLL(I+1)-1
            J = J_KLL(K)
            VAL = KLL(K)
            CALL ACCUMULATE_PARTITION(I, J, VAL, ACTIVE_MAP, ZERO_MAP, KAA, KZZ, KZA_COUNT)
         ENDDO
      ENDDO

      ALLOCATE(KZA_COLPTR(NACTIVE+1))
      KZA_COLPTR(1) = 1
      DO I=1,NACTIVE
         KZA_COLPTR(I+1) = KZA_COLPTR(I) + KZA_COUNT(I)
      ENDDO

      ALLOCATE(KZA_ROW(MAX(1,KZA_COLPTR(NACTIVE+1)-1)), KZA_VAL(MAX(1,KZA_COLPTR(NACTIVE+1)-1)), KZA_NEXT(MAX(1,NACTIVE)))
      KZA_ROW  = 0
      KZA_VAL  = ZERO
      KZA_NEXT = KZA_COLPTR(1:MAX(1,NACTIVE))

      DO I=1,NDOFL
         DO K=I_KLL(I),I_KLL(I+1)-1
            J = J_KLL(K)
            VAL = KLL(K)
            CALL FILL_KZA(I, J, VAL, ACTIVE_MAP, ZERO_MAP, KZA_NEXT, KZA_ROW, KZA_VAL)
         ENDDO
      ENDDO

      KAA = HALF_SYM(KAA)
      IF (NZERO > 0) KZZ(1:NZERO,1:NZERO) = HALF_SYM(KZZ(1:NZERO,1:NZERO))

      END SUBROUTINE BUILD_CONDENSED_MODAL_CORE

!***********************************************************************************************************************************
      SUBROUTINE ACCUMULATE_PARTITION ( IROW, JCOL, VAL, ACTIVE_MAP, ZERO_MAP, KAA, KZZ, KZA_COUNT )

      INTEGER(LONG), INTENT(IN)       :: IROW
      INTEGER(LONG), INTENT(IN)       :: JCOL
      INTEGER(LONG), INTENT(IN)       :: ACTIVE_MAP(:)
      INTEGER(LONG), INTENT(IN)       :: ZERO_MAP(:)
      INTEGER(LONG), INTENT(INOUT)    :: KZA_COUNT(:)
      REAL(DOUBLE), INTENT(IN)        :: VAL
      REAL(DOUBLE), INTENT(INOUT)     :: KAA(:,:)
      REAL(DOUBLE), INTENT(INOUT)     :: KZZ(:,:)

      INTEGER(LONG)                   :: IA
      INTEGER(LONG)                   :: IZ
      INTEGER(LONG)                   :: JA
      INTEGER(LONG)                   :: JZ

      IA = ACTIVE_MAP(IROW)
      JA = ACTIVE_MAP(JCOL)
      IZ = ZERO_MAP(IROW)
      JZ = ZERO_MAP(JCOL)

      IF (SPARSTOR == 'SYM   ') THEN
         IF ((IA > 0) .AND. (JA > 0)) THEN
            KAA(IA,JA) = KAA(IA,JA) + VAL
            IF (IA /= JA) KAA(JA,IA) = KAA(JA,IA) + VAL
         ELSE IF ((IZ > 0) .AND. (JZ > 0) .AND. (SIZE(KZZ,1) > 0)) THEN
            KZZ(IZ,JZ) = KZZ(IZ,JZ) + VAL
            IF (IZ /= JZ) KZZ(JZ,IZ) = KZZ(JZ,IZ) + VAL
         ELSE IF ((IA > 0) .AND. (JZ > 0)) THEN
            KZA_COUNT(IA) = KZA_COUNT(IA) + 1
         ELSE IF ((IZ > 0) .AND. (JA > 0)) THEN
            KZA_COUNT(JA) = KZA_COUNT(JA) + 1
         ENDIF
      ELSE
         IF ((IA > 0) .AND. (JA > 0)) THEN
            KAA(IA,JA) = KAA(IA,JA) + VAL
         ELSE IF ((IZ > 0) .AND. (JZ > 0) .AND. (SIZE(KZZ,1) > 0)) THEN
            KZZ(IZ,JZ) = KZZ(IZ,JZ) + VAL
         ELSE IF ((IZ > 0) .AND. (JA > 0)) THEN
            KZA_COUNT(JA) = KZA_COUNT(JA) + 1
         ENDIF
      ENDIF

      END SUBROUTINE ACCUMULATE_PARTITION

!***********************************************************************************************************************************
      SUBROUTINE FILL_KZA ( IROW, JCOL, VAL, ACTIVE_MAP, ZERO_MAP, KZA_NEXT, KZA_ROW, KZA_VAL )

      INTEGER(LONG), INTENT(IN)       :: IROW
      INTEGER(LONG), INTENT(IN)       :: JCOL
      INTEGER(LONG), INTENT(IN)       :: ACTIVE_MAP(:)
      INTEGER(LONG), INTENT(IN)       :: ZERO_MAP(:)
      INTEGER(LONG), INTENT(INOUT)    :: KZA_NEXT(:)
      INTEGER(LONG), INTENT(INOUT)    :: KZA_ROW(:)
      REAL(DOUBLE), INTENT(IN)        :: VAL
      REAL(DOUBLE), INTENT(INOUT)     :: KZA_VAL(:)

      INTEGER(LONG)                   :: IA
      INTEGER(LONG)                   :: IZ
      INTEGER(LONG)                   :: JA
      INTEGER(LONG)                   :: JZ
      INTEGER(LONG)                   :: P

      IA = ACTIVE_MAP(IROW)
      JA = ACTIVE_MAP(JCOL)
      IZ = ZERO_MAP(IROW)
      JZ = ZERO_MAP(JCOL)

      IF (SPARSTOR == 'SYM   ') THEN
         IF ((IA > 0) .AND. (JZ > 0)) THEN
            P = KZA_NEXT(IA)
            KZA_ROW(P) = JZ
            KZA_VAL(P) = VAL
            KZA_NEXT(IA) = P + 1
         ELSE IF ((IZ > 0) .AND. (JA > 0)) THEN
            P = KZA_NEXT(JA)
            KZA_ROW(P) = IZ
            KZA_VAL(P) = VAL
            KZA_NEXT(JA) = P + 1
         ENDIF
      ELSE
         IF ((IZ > 0) .AND. (JA > 0)) THEN
            P = KZA_NEXT(JA)
            KZA_ROW(P) = IZ
            KZA_VAL(P) = KZA_VAL(P) + VAL
            KZA_NEXT(JA) = P + 1
         ENDIF
      ENDIF

      END SUBROUTINE FILL_KZA

!***********************************************************************************************************************************
      SUBROUTINE FORM_CONDENSED_STIFFNESS ( NACTIVE, NZERO, KAA, KZZ, KZA_COLPTR, KZA_ROW, KZA_VAL, KCOND, INFO )

      INTEGER(LONG), INTENT(IN)       :: NACTIVE
      INTEGER(LONG), INTENT(IN)       :: NZERO
      INTEGER(LONG), INTENT(IN)       :: KZA_COLPTR(:)
      INTEGER(LONG), INTENT(IN)       :: KZA_ROW(:)
      INTEGER(LONG), INTENT(OUT)      :: INFO
      REAL(DOUBLE), INTENT(IN)        :: KAA(:,:)
      REAL(DOUBLE), INTENT(INOUT)     :: KZZ(:,:)
      REAL(DOUBLE), INTENT(IN)        :: KZA_VAL(:)
      REAL(DOUBLE), ALLOCATABLE       :: KCOND(:,:)

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: J
      REAL(DOUBLE), ALLOCATABLE       :: RHS(:)
      REAL(DOUBLE), ALLOCATABLE       :: Y(:)

      INFO = 0
      ALLOCATE(KCOND(NACTIVE,NACTIVE))
      KCOND = KAA

      IF (NZERO <= 0) RETURN

      CALL DPOTRF('U', NZERO, KZZ, NZERO, INFO)
      IF (INFO /= 0) RETURN

      ALLOCATE(RHS(NZERO), Y(NACTIVE))
      DO I=1,NACTIVE
         RHS = ZERO
         IF (KZA_COLPTR(I) < KZA_COLPTR(I+1)) THEN
            DO J=KZA_COLPTR(I),KZA_COLPTR(I+1)-1
               RHS(KZA_ROW(J)) = RHS(KZA_ROW(J)) + KZA_VAL(J)
            ENDDO
         ENDIF
         CALL DPOTRS('U', NZERO, 1, KZZ, NZERO, RHS, NZERO, INFO)
         IF (INFO /= 0) RETURN
         CALL CSC_TRANSPOSE_MATVEC(NACTIVE, KZA_COLPTR, KZA_ROW, KZA_VAL, RHS, Y)
         KCOND(:,I) = KCOND(:,I) - Y
      ENDDO

      KCOND = HALF_SYM(KCOND)

      END SUBROUTINE FORM_CONDENSED_STIFFNESS

!***********************************************************************************************************************************
      SUBROUTINE BUILD_STANDARDIZED_OPERATOR ( KCOND, MASS_A, ASTD )

      REAL(DOUBLE), INTENT(IN)        :: KCOND(:,:)
      REAL(DOUBLE), INTENT(IN)        :: MASS_A(:)
      REAL(DOUBLE), ALLOCATABLE       :: ASTD(:,:)

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: J

      ALLOCATE(ASTD(SIZE(KCOND,1),SIZE(KCOND,2)))
      DO I=1,SIZE(KCOND,1)
         DO J=1,SIZE(KCOND,2)
            ASTD(I,J) = KCOND(I,J)/SQRT(MASS_A(I)*MASS_A(J))
         ENDDO
      ENDDO

      END SUBROUTINE BUILD_STANDARDIZED_OPERATOR

!***********************************************************************************************************************************
      SUBROUTINE RUN_DENSE_BACKEND ( ASTD, EVAL_ALL, EVEC_STD, INFO )

      REAL(DOUBLE), INTENT(IN)        :: ASTD(:,:)
      INTEGER(LONG), INTENT(OUT)      :: INFO
      REAL(DOUBLE), ALLOCATABLE       :: EVAL_ALL(:)
      REAL(DOUBLE), ALLOCATABLE       :: EVEC_STD(:,:)

      INTEGER(LONG)                   :: LWORK
      REAL(DOUBLE), ALLOCATABLE       :: WORK(:)

      ALLOCATE(EVEC_STD(SIZE(ASTD,1),SIZE(ASTD,2)), EVAL_ALL(SIZE(ASTD,1)))
      EVEC_STD = ASTD

      LWORK = -1
      ALLOCATE(WORK(1))
      CALL DSYEV('V', 'U', SIZE(ASTD,1), EVEC_STD, SIZE(ASTD,1), EVAL_ALL, WORK, LWORK, INFO)
      IF (INFO /= 0) RETURN
      LWORK = MAX(1, NINT(WORK(1)))
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK))
      CALL DSYEV('V', 'U', SIZE(ASTD,1), EVEC_STD, SIZE(ASTD,1), EVAL_ALL, WORK, LWORK, INFO)

      END SUBROUTINE RUN_DENSE_BACKEND

!***********************************************************************************************************************************
      SUBROUTINE RUN_SUBSPACE_BACKEND ( ASTD, NEV, NSUB, TOL, MAXIT, EVAL_ALL, EVEC_STD, INFO )

      REAL(DOUBLE), INTENT(IN)        :: ASTD(:,:)
      REAL(DOUBLE), INTENT(IN)        :: TOL
      INTEGER(LONG), INTENT(IN)       :: MAXIT
      INTEGER(LONG), INTENT(IN)       :: NEV
      INTEGER(LONG), INTENT(IN)       :: NSUB
      INTEGER(LONG), INTENT(OUT)      :: INFO
      REAL(DOUBLE), ALLOCATABLE       :: EVAL_ALL(:)
      REAL(DOUBLE), ALLOCATABLE       :: EVEC_STD(:,:)

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: ITER
      INTEGER(LONG)                   :: LWORK
      INTEGER(LONG)                   :: NSUB_EFF
      REAL(DOUBLE)                    :: MAX_RESID
      REAL(DOUBLE), ALLOCATABLE       :: AFAC(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: AX(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: RESID(:)
      REAL(DOUBLE), ALLOCATABLE       :: T(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: W(:)
      REAL(DOUBLE), ALLOCATABLE       :: WORK(:)
      REAL(DOUBLE), ALLOCATABLE       :: X(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: Y(:,:)

      INFO = 0
      NSUB_EFF = MIN(SIZE(ASTD,1), MAX(NEV, NSUB))
      ALLOCATE(AFAC(SIZE(ASTD,1),SIZE(ASTD,2)), X(SIZE(ASTD,1),NSUB_EFF), Y(SIZE(ASTD,1),NSUB_EFF),                      &
               AX(SIZE(ASTD,1),NEV), T(NSUB_EFF,NSUB_EFF), W(NSUB_EFF), RESID(NEV))
      AFAC = ASTD
      CALL DPOTRF('U', SIZE(ASTD,1), AFAC, SIZE(ASTD,1), INFO)
      IF (INFO /= 0) RETURN

      CALL RANDOM_NUMBER(X)
      CALL ORTHONORMALIZE_MGS(X)

      LWORK = MAX(1, 8*NSUB_EFF)
      ALLOCATE(WORK(LWORK))
      DO ITER=1,MAX(1,MAXIT)
         Y = X
         CALL DPOTRS('U', SIZE(ASTD,1), NSUB_EFF, AFAC, SIZE(ASTD,1), Y, SIZE(ASTD,1), INFO)
         IF (INFO /= 0) RETURN
         CALL ORTHONORMALIZE_MGS(Y)
         T = MATMUL(TRANSPOSE(Y), MATMUL(ASTD, Y))
         CALL DSYEV('V', 'U', NSUB_EFF, T, NSUB_EFF, W, WORK, LWORK, INFO)
         IF (INFO /= 0) RETURN
         X = MATMUL(Y, T)
         AX = MATMUL(ASTD, X(:,1:NEV))
         MAX_RESID = ZERO
         DO I=1,NEV
            RESID(I) = SQRT(DOT_PRODUCT(AX(:,I) - W(I)*X(:,I), AX(:,I) - W(I)*X(:,I)))
            IF (RESID(I) > MAX_RESID) MAX_RESID = RESID(I)
         ENDDO
         IF (MAX_RESID < TOL) EXIT
      ENDDO

      ALLOCATE(EVAL_ALL(NEV), EVEC_STD(SIZE(ASTD,1),NEV))
      EVAL_ALL = W(1:NEV)
      EVEC_STD = X(:,1:NEV)

      END SUBROUTINE RUN_SUBSPACE_BACKEND

!***********************************************************************************************************************************
      SUBROUTINE RUN_CHASE_BACKEND ( ASTD, NEV, NEX, TOL, MAXIT, EVAL_ALL, EVEC_STD, INFO )

      REAL(DOUBLE), INTENT(IN)        :: ASTD(:,:)
      REAL(DOUBLE), INTENT(IN)        :: TOL
      INTEGER(LONG), INTENT(IN)       :: MAXIT
      INTEGER(LONG), INTENT(IN)       :: NEV
      INTEGER(LONG), INTENT(IN)       :: NEX
      INTEGER(LONG), INTENT(OUT)      :: INFO
      REAL(DOUBLE), ALLOCATABLE       :: EVAL_ALL(:)
      REAL(DOUBLE), ALLOCATABLE       :: EVEC_STD(:,:)

#ifdef MYSTRAN_HAVE_EXTERNAL_CHASE
      INTEGER(LONG)                   :: CHASE_IDO
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: J
      INTEGER(LONG)                   :: NEX_EFF
      REAL(DOUBLE), ALLOCATABLE       :: EVAL_TMP(:)
      REAL(DOUBLE), ALLOCATABLE       :: EVEC_TMP(:,:)
      REAL(DOUBLE)                    :: TMPVAL

      INFO = 0
      NEX_EFF = MIN(MAX(1,NEX), MAX(1,SIZE(ASTD,1) - MAX(1,NEV)))
      ALLOCATE(EVAL_TMP(NEV+NEX_EFF), EVEC_TMP(SIZE(ASTD,1),NEV+NEX_EFF))
      EVAL_TMP = ZERO
      EVEC_TMP = ZERO
      CHASE_IDO = 0

      CALL DCHASE_INIT(SIZE(ASTD,1), NEV, NEX_EFF, ASTD, SIZE(ASTD,1), EVEC_TMP, EVAL_TMP, CHASE_IDO)
      CALL DCHASE(MAXIT, TOL, 'R', 'S', 'C')
      CALL DCHASE_FINALIZE(CHASE_IDO)
      CALL SORT_EIGENPAIRS(EVAL_TMP, EVEC_TMP)

      ALLOCATE(EVAL_ALL(NEV), EVEC_STD(SIZE(ASTD,1),NEV))
      EVAL_ALL = EVAL_TMP(1:NEV)
      EVEC_STD = EVEC_TMP(:,1:NEV)
#else
      INFO = 1
#endif

      END SUBROUTINE RUN_CHASE_BACKEND

!***********************************************************************************************************************************
      SUBROUTINE RUN_FEAST_BACKEND ( KCOND, MASS_A, EMIN, EMAX, EVAL_ALL, EVEC_ACTIVE, INFO )

      REAL(DOUBLE), INTENT(IN)        :: EMAX
      REAL(DOUBLE), INTENT(IN)        :: EMIN
      REAL(DOUBLE), INTENT(IN)        :: KCOND(:,:)
      REAL(DOUBLE), INTENT(IN)        :: MASS_A(:)
      INTEGER(LONG), INTENT(OUT)      :: INFO
      REAL(DOUBLE), ALLOCATABLE       :: EVAL_ALL(:)
      REAL(DOUBLE), ALLOCATABLE       :: EVEC_ACTIVE(:,:)

#ifdef MYSTRAN_HAVE_EXTERNAL_FEAST
      INTEGER(LONG)                   :: FEAST_FPM(128)
      INTEGER(LONG)                   :: FEAST_INFO
      INTEGER(LONG)                   :: FEAST_LOOP
      INTEGER(LONG)                   :: FEAST_M0
      INTEGER(LONG)                   :: FEAST_MODE
      REAL(DOUBLE)                    :: EPSOUT
      REAL(DOUBLE), ALLOCATABLE       :: A(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: B(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: LAMBDA(:)
      REAL(DOUBLE), ALLOCATABLE       :: Q(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: RES(:)

      INTERFACE
         SUBROUTINE FEASTINIT(FPM)
            INTEGER :: FPM(*)
         END SUBROUTINE FEASTINIT
         SUBROUTINE DFEAST_SYGV(UPLO, N, A, LDA, B, LDB, FPM, EPSOUT, LOOP, EMIN, EMAX, M0, LAMBDA, Q, MODE, RES, INFO)
            CHARACTER(1), INTENT(IN) :: UPLO
            INTEGER, INTENT(IN) :: N, LDA, LDB, FPM(*), M0
            DOUBLE PRECISION, INTENT(INOUT) :: A(LDA,*), B(LDB,*)
            DOUBLE PRECISION, INTENT(IN) :: EMIN, EMAX
            DOUBLE PRECISION, INTENT(OUT) :: EPSOUT, LAMBDA(*), Q(N,*), RES(*)
            INTEGER, INTENT(OUT) :: LOOP, MODE, INFO
         END SUBROUTINE DFEAST_SYGV
      END INTERFACE

      FEAST_M0 = MAX(1, MIN(SIZE(KCOND,1), EIG_FEAST_M0))
      ALLOCATE(A(SIZE(KCOND,1),SIZE(KCOND,2)), B(SIZE(KCOND,1),SIZE(KCOND,2)), LAMBDA(FEAST_M0), Q(SIZE(KCOND,1),FEAST_M0),      &
               RES(FEAST_M0))
      A = KCOND
      B = ZERO
      B = ZERO
      DO INFO=1,SIZE(MASS_A)
         B(INFO,INFO) = MASS_A(INFO)
      ENDDO
      CALL FEASTINIT(FEAST_FPM)
      FEAST_FPM(1) = MERGE(1,0,SUPINFO == 'N')
      FEAST_FPM(3) = EIG_FEAST_N_CONTOUR
      FEAST_FPM(4) = EIG_FEAST_MAX_LOOP
      CALL DFEAST_SYGV('U', SIZE(KCOND,1), A, SIZE(KCOND,1), B, SIZE(KCOND,1), FEAST_FPM, EPSOUT, FEAST_LOOP, EMIN, EMAX,         &
                       FEAST_M0, LAMBDA, Q, FEAST_MODE, RES, FEAST_INFO)
      INFO = FEAST_INFO
      IF ((INFO == 0) .AND. (FEAST_MODE > 0)) THEN
         ALLOCATE(EVAL_ALL(FEAST_MODE), EVEC_ACTIVE(SIZE(KCOND,1),FEAST_MODE))
         EVAL_ALL = LAMBDA(1:FEAST_MODE)
         EVEC_ACTIVE = Q(:,1:FEAST_MODE)
         CALL SORT_EIGENPAIRS(EVAL_ALL, EVEC_ACTIVE)
      ENDIF
#else
      INFO = 1
#endif

      END SUBROUTINE RUN_FEAST_BACKEND

!***********************************************************************************************************************************
      SUBROUTINE STANDARDIZED_TO_ACTIVE ( MASS_A, EVEC_STD, EVEC_ACTIVE )

      REAL(DOUBLE), INTENT(IN)        :: MASS_A(:)
      REAL(DOUBLE), INTENT(IN)        :: EVEC_STD(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: EVEC_ACTIVE(:,:)

      INTEGER(LONG)                   :: I

      ALLOCATE(EVEC_ACTIVE(SIZE(EVEC_STD,1),SIZE(EVEC_STD,2)))
      DO I=1,SIZE(MASS_A)
         EVEC_ACTIVE(I,:) = EVEC_STD(I,:)/SQRT(MASS_A(I))
      ENDDO

      END SUBROUTINE STANDARDIZED_TO_ACTIVE

!***********************************************************************************************************************************
      SUBROUTINE SELECT_REQUESTED_MODES ( EVAL_ALL, EMIN, EMAX, KEEP_IDX, KEEP_COUNT )

      REAL(DOUBLE), INTENT(IN)        :: EMAX
      REAL(DOUBLE), INTENT(IN)        :: EMIN
      REAL(DOUBLE), INTENT(IN)        :: EVAL_ALL(:)
      INTEGER(LONG), INTENT(OUT)      :: KEEP_COUNT
      INTEGER(LONG), ALLOCATABLE      :: KEEP_IDX(:)

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: LIMIT

      LIMIT = MAX(1,EIG_N2)
      ALLOCATE(KEEP_IDX(MAX(1,MIN(LIMIT,SIZE(EVAL_ALL)))))
      KEEP_COUNT = 0
      DO I=1,SIZE(EVAL_ALL)
         IF (EVAL_ALL(I) < EMIN) CYCLE
         IF (EMAX > EMIN) THEN
            IF (EVAL_ALL(I) > EMAX) CYCLE
         ENDIF
         KEEP_COUNT = KEEP_COUNT + 1
         IF (KEEP_COUNT <= SIZE(KEEP_IDX)) KEEP_IDX(KEEP_COUNT) = I
         IF (KEEP_COUNT == LIMIT) EXIT
      ENDDO

      END SUBROUTINE SELECT_REQUESTED_MODES

!***********************************************************************************************************************************
      SUBROUTINE EXPAND_FULL_EIGENVECTORS ( KEEP_IDX, KEEP_COUNT, ACTIVE_POS, ZERO_POS, MASS_A, EVEC_ACTIVE, KZZ, KZA_COLPTR,      &
                                            KZA_ROW, KZA_VAL, EVEC_FULL, INFO )

      INTEGER(LONG), INTENT(IN)       :: KEEP_COUNT
      INTEGER(LONG), INTENT(IN)       :: ACTIVE_POS(:)
      INTEGER(LONG), INTENT(IN)       :: KEEP_IDX(:)
      INTEGER(LONG), INTENT(IN)       :: KZA_COLPTR(:)
      INTEGER(LONG), INTENT(IN)       :: KZA_ROW(:)
      INTEGER(LONG), INTENT(IN)       :: ZERO_POS(:)
      INTEGER(LONG), INTENT(OUT)      :: INFO
      REAL(DOUBLE), INTENT(IN)        :: EVEC_ACTIVE(:,:)
      REAL(DOUBLE), INTENT(INOUT)     :: KZZ(:,:)
      REAL(DOUBLE), INTENT(IN)        :: KZA_VAL(:)
      REAL(DOUBLE), INTENT(IN)        :: MASS_A(:)
      REAL(DOUBLE), ALLOCATABLE       :: EVEC_FULL(:,:)

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: J
      REAL(DOUBLE), ALLOCATABLE       :: RHS(:)

      INFO = 0
      ALLOCATE(EVEC_FULL(NDOFL,KEEP_COUNT))
      EVEC_FULL = ZERO

      DO J=1,KEEP_COUNT
         DO I=1,SIZE(ACTIVE_POS)
            EVEC_FULL(ACTIVE_POS(I),J) = EVEC_ACTIVE(I,KEEP_IDX(J))
         ENDDO
      ENDDO

      IF (SIZE(ZERO_POS) <= 1) RETURN

      ALLOCATE(RHS(SIZE(ZERO_POS)))
      DO J=1,KEEP_COUNT
         CALL CSC_MATVEC(SIZE(ACTIVE_POS), KZA_COLPTR, KZA_ROW, KZA_VAL, EVEC_ACTIVE(:,KEEP_IDX(J)), RHS)
         RHS = -RHS
         CALL DPOTRS('U', SIZE(ZERO_POS), 1, KZZ, SIZE(ZERO_POS), RHS, SIZE(ZERO_POS), INFO)
         IF (INFO /= 0) RETURN
         DO I=1,SIZE(ZERO_POS)
            IF (ZERO_POS(I) > 0) EVEC_FULL(ZERO_POS(I),J) = RHS(I)
         ENDDO
      ENDDO

      END SUBROUTINE EXPAND_FULL_EIGENVECTORS

!***********************************************************************************************************************************
      SUBROUTINE CSC_TRANSPOSE_MATVEC ( NCOL, COLPTR, ROWIND, VALUES, X, Y )

      INTEGER(LONG), INTENT(IN)       :: NCOL
      INTEGER(LONG), INTENT(IN)       :: COLPTR(:)
      INTEGER(LONG), INTENT(IN)       :: ROWIND(:)
      REAL(DOUBLE), INTENT(IN)        :: VALUES(:)
      REAL(DOUBLE), INTENT(IN)        :: X(:)
      REAL(DOUBLE), INTENT(OUT)       :: Y(:)

      INTEGER(LONG)                   :: J
      INTEGER(LONG)                   :: P

      Y = ZERO
      DO J=1,NCOL
         DO P=COLPTR(J),COLPTR(J+1)-1
            Y(J) = Y(J) + VALUES(P)*X(ROWIND(P))
         ENDDO
      ENDDO

      END SUBROUTINE CSC_TRANSPOSE_MATVEC

!***********************************************************************************************************************************
      SUBROUTINE CSC_MATVEC ( NCOL, COLPTR, ROWIND, VALUES, X, Y )

      INTEGER(LONG), INTENT(IN)       :: NCOL
      INTEGER(LONG), INTENT(IN)       :: COLPTR(:)
      INTEGER(LONG), INTENT(IN)       :: ROWIND(:)
      REAL(DOUBLE), INTENT(IN)        :: VALUES(:)
      REAL(DOUBLE), INTENT(IN)        :: X(:)
      REAL(DOUBLE), INTENT(OUT)       :: Y(:)

      INTEGER(LONG)                   :: J
      INTEGER(LONG)                   :: P

      Y = ZERO
      DO J=1,NCOL
         DO P=COLPTR(J),COLPTR(J+1)-1
            Y(ROWIND(P)) = Y(ROWIND(P)) + VALUES(P)*X(J)
         ENDDO
      ENDDO

      END SUBROUTINE CSC_MATVEC

!***********************************************************************************************************************************
      SUBROUTINE ORTHONORMALIZE_MGS ( X )

      REAL(DOUBLE), INTENT(INOUT)     :: X(:,:)

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: J
      REAL(DOUBLE)                    :: NRM

      DO J=1,SIZE(X,2)
         DO I=1,J-1
            X(:,J) = X(:,J) - DOT_PRODUCT(X(:,I), X(:,J))*X(:,I)
         ENDDO
         NRM = SQRT(MAX(DOT_PRODUCT(X(:,J), X(:,J)), ZERO))
         IF (NRM > ZERO) X(:,J) = X(:,J)/NRM
      ENDDO

      END SUBROUTINE ORTHONORMALIZE_MGS

!***********************************************************************************************************************************
      SUBROUTINE SORT_EIGENPAIRS ( EVALS, EVECS )

      REAL(DOUBLE), INTENT(INOUT)     :: EVALS(:)
      REAL(DOUBLE), INTENT(INOUT)     :: EVECS(:,:)

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: J
      REAL(DOUBLE)                    :: TMP
      REAL(DOUBLE), ALLOCATABLE       :: TMPV(:)

      ALLOCATE(TMPV(SIZE(EVECS,1)))
      DO I=1,SIZE(EVALS)-1
         DO J=I+1,SIZE(EVALS)
            IF (EVALS(J) < EVALS(I)) THEN
               TMP      = EVALS(I)
               EVALS(I) = EVALS(J)
               EVALS(J) = TMP
               TMPV     = EVECS(:,I)
               EVECS(:,I) = EVECS(:,J)
               EVECS(:,J) = TMPV
            ENDIF
         ENDDO
      ENDDO

      END SUBROUTINE SORT_EIGENPAIRS

!***********************************************************************************************************************************
      SUBROUTINE WRITE_FALLBACK_WARNING ( METHOD, CODE, MESSAGE )

      CHARACTER(LEN=*), INTENT(IN)    :: MESSAGE
      CHARACTER(LEN=*), INTENT(IN)    :: METHOD
      INTEGER(LONG), INTENT(IN)       :: CODE

      WARN_ERR = WARN_ERR + 1
      WRITE(ERR,'(A,I4,2A)') ' *WARNING ', CODE, ': ', TRIM(METHOD)//' - '//TRIM(MESSAGE)
      IF (SUPINFO == 'N') THEN
         WRITE(F06,'(A,I4,2A)') ' *WARNING ', CODE, ': ', TRIM(METHOD)//' - '//TRIM(MESSAGE)
      ENDIF

      END SUBROUTINE WRITE_FALLBACK_WARNING

!***********************************************************************************************************************************
      REAL(DOUBLE) FUNCTION FREQ_TO_LAMBDA ( FREQ )

      REAL(DOUBLE), INTENT(IN)        :: FREQ

      FREQ_TO_LAMBDA = (TWO_PI()*FREQ)**2

      END FUNCTION FREQ_TO_LAMBDA

!***********************************************************************************************************************************
      REAL(DOUBLE) FUNCTION ONE ()

      ONE = 1.0D0

      END FUNCTION ONE

!***********************************************************************************************************************************
      REAL(DOUBLE) FUNCTION TWO_PI ()

      TWO_PI = 2.0D0*ACOS(-1.0D0)

      END FUNCTION TWO_PI

!***********************************************************************************************************************************
      FUNCTION HALF_SYM ( A ) RESULT ( B )

      REAL(DOUBLE), INTENT(IN)        :: A(:,:)
      REAL(DOUBLE)                    :: B(SIZE(A,1),SIZE(A,2))

      B = 0.5D0*(A + TRANSPOSE(A))

      END FUNCTION HALF_SYM

      END MODULE EIGRL_EXTRACT_SOLVERS
! --- chase_feast_add --- end !
