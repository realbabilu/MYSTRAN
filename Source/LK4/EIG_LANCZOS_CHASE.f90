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

! !--- CHASE and FEAST --- begin!
      SUBROUTINE EIG_LANCZOS_CHASE

! CHASE backend for Lanczos entry point.
! Native path is optional (external wrapper/library). If unavailable or
! if native solve fails, this routine falls back to ARPACK Lanczos
! (not MGIV), so behavior stays consistent with Lanczos workflows.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, WARN_ERR, NDOFL, NUM_EIGENS, NVEC, NTERM_KLL, NTERM_MLL, SOL_NAME
      USE CONSTANTS_1, ONLY           :  ZERO, PI
      USE PARAMS, ONLY                :  SUPINFO
      USE MODEL_STUF, ONLY            :  EIG_FRQ1, EIG_FRQ2, EIG_N2, EIG_METH, EIG_MSGLVL
      USE SPARSE_MATRICES, ONLY       :  I_KLL, J_KLL, KLL, I_MLL, J_MLL, MLL
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL, EIGEN_VEC, MODE_NUM

      USE EIG_LANCZOS_CHASE_USE_IFs
      USE LINK_MESSAGE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EIG_LANCZOS_CHASE'
      INTEGER(LONG)                   :: I, J
      INTEGER(LONG)                   :: INFO
      INTEGER(LONG)                   :: NFOUND
      INTEGER(LONG)                   :: NEV_REQ
      INTEGER(LONG)                   :: MAXSUB
      REAL(DOUBLE)                    :: EMIN
      REAL(DOUBLE)                    :: EMAX
      REAL(DOUBLE), ALLOCATABLE       :: EVAL_TMP(:)
      REAL(DOUBLE), ALLOCATABLE       :: EVEC_TMP(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: RES_TMP(:)

#ifdef MYSTRAN_HAVE_EXTERNAL_CHASE
      INTERFACE
         SUBROUTINE MYSTRAN_CHASE_DSYGV_CRS(N,NNZA,IA,JA,A,NNZB,IB,JB,B,NEV,MAXSUB,TOL,MAXIT,EMIN,EMAX,EVAL,EVEC,RES,NFOUND,INFO)
            USE PENTIUM_II_KIND, ONLY : LONG, DOUBLE
            INTEGER(LONG) :: N, NNZA, IA(*), JA(*), NNZB, IB(*), JB(*), NEV, MAXSUB, MAXIT, NFOUND, INFO
            REAL(DOUBLE)  :: A(*), B(*), TOL, EMIN, EMAX, EVAL(*), EVEC(N,*), RES(*)
         END SUBROUTINE MYSTRAN_CHASE_DSYGV_CRS
      END INTERFACE
#endif

#ifdef MYSTRAN_HAVE_EXTERNAL_CHASE
      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,4915)
         IF (SUPINFO == 'N') WRITE(F06,4915)
         GOTO 900
      ENDIF

      IF (EIG_FRQ2 <= ZERO) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,4916)
         IF (SUPINFO == 'N') WRITE(F06,4916)
         GOTO 900
      ENDIF

      CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - CHASE NATIVE (GENERALIZED SYMMETRIC SPARSE)')

      WARN_ERR = WARN_ERR + 1
      WRITE(ERR,4918)
      IF (SUPINFO == 'N') WRITE(F06,4918)

      NEV_REQ = EIG_N2
      IF (NEV_REQ < 1) NEV_REQ = 1
      IF (NEV_REQ > NDOFL) NEV_REQ = NDOFL

      MAXSUB = MAX(2*NEV_REQ, NEV_REQ + 8)
      IF (MAXSUB > NDOFL) MAXSUB = NDOFL
      IF (MAXSUB < 8) MAXSUB = MIN(8, NDOFL)

      EMIN = ZERO
      IF (EIG_FRQ1 > ZERO) EMIN = (2.0D0*PI*EIG_FRQ1)**2
      EMAX = (2.0D0*PI*EIG_FRQ2)**2

      ALLOCATE(EVAL_TMP(MAXSUB), EVEC_TMP(NDOFL,MAXSUB), RES_TMP(MAXSUB))
      INFO = -999
      NFOUND = 0

      CALL MYSTRAN_CHASE_DSYGV_CRS( NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, &
                                    NTERM_MLL, I_MLL, J_MLL, MLL,         &
                                    NEV_REQ, MAXSUB, 1.0D-10, 300,        &
                                    EMIN, EMAX, EVAL_TMP, EVEC_TMP, RES_TMP, NFOUND, INFO )

      IF ((INFO == 0) .AND. (NFOUND > 0)) THEN
         IF (NFOUND > NEV_REQ) NFOUND = NEV_REQ

         CALL ALLOCATE_EIGEN1_MAT('EIGEN_VEC', NDOFL, NFOUND, SUBR_NAME)
         CALL ALLOCATE_EIGEN1_MAT('MODE_NUM' , NDOFL, 1,      SUBR_NAME)
         CALL ALLOCATE_EIGEN1_MAT('EIGEN_VAL', NDOFL, 1,      SUBR_NAME)

         DO I=1,NFOUND
            EIGEN_VAL(I) = EVAL_TMP(I)
            MODE_NUM(I)  = I
            DO J=1,NDOFL
               EIGEN_VEC(J,I) = EVEC_TMP(J,I)
            ENDDO
         ENDDO

         NUM_EIGENS = NFOUND
         NVEC       = NFOUND

         DEALLOCATE(EVAL_TMP, EVEC_TMP, RES_TMP)
         RETURN
      ENDIF

      WARN_ERR = WARN_ERR + 1
      WRITE(ERR,4917) INFO, NFOUND
      IF (SUPINFO == 'N') WRITE(F06,4917) INFO, NFOUND

      DEALLOCATE(EVAL_TMP, EVEC_TMP, RES_TMP)
#else
      WARN_ERR = WARN_ERR + 1
      WRITE(ERR,4912)
      IF (SUPINFO == 'N') WRITE(F06,4912)
#endif

 900  CONTINUE
      CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - CHASE FALLBACK (ARPACK LANCZOS)')
      CALL EIG_LANCZOS_ARPACK

      RETURN

 4912 FORMAT(' *WARNING 4912: CHASE EXTERNAL BACKEND NOT LINKED. USING ARPACK LANCZOS FALLBACK.')
 4915 FORMAT(' *WARNING 4915: CHASE NATIVE FOR BUCKLING (KLL,KLLD) NOT WIRED YET. USING ARPACK LANCZOS FALLBACK.')
 4916 FORMAT(' *WARNING 4916: CHASE NATIVE REQUIRES EIGRL FREQUENCY INTERVAL (V1/V2). USING ARPACK LANCZOS FALLBACK.')
 4917 FORMAT(' *WARNING 4917: CHASE NATIVE FAILED (INFO=',I8,', NFOUND=',I8,'). USING ARPACK LANCZOS FALLBACK.')
 4918 FORMAT(' *WARNING 4918: CHASE EXTERNAL BACKEND LINKED. RUNNING NATIVE CHASE PATH FOR MODES.')

      END SUBROUTINE EIG_LANCZOS_CHASE
! !--- CHASE and FEAST --- end!
