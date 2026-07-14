! ##################################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________

! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)

! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:

! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________

! End MIT license text.

      SUBROUTINE LINK4

! Calculates system eigenvalues, eigenvectors. There are 4 eigenvalue extraction methods in MYSTRAN, none of which seem suited to
! very large eigenvalue problems for one reason or another:

!   (1) LANCZOS method:

!          calculates some eigenvalues and some eigenvectors of KLL, MLL (or KLL, KLLD for BUCKLING). This is the most widely used
!          eigenvalue extraction method for large problems. The LANCZOS method in MYSTRAN uses 1 of 2 algorithms:
!             (a) LAPACK: requires KLL, MLL (or KLL, KLLD for buckling) to be in band storage - sparse storage can NOT be used

!          Tis Lanczos algorithms is not practical for very large eigenvalue problems since LAPACK/ARPACK will require
!          large amounts of memory to store the banded KLL, MLL matrices.

!   (2) GIV (Givens) method:
!          calculates all eigenvalues and some eigenvactors of KLL, MLL (or KLL, KLLD for buckling). This method is only practical
!          for relatively small problems. It requires MLL (or KLLD for buckling) to be a positive definite matrix. The algorithm
!          performs a Cholesky decomp of matrix MLL (or KLLD for buckling) which can be time consuming for large problems

!   (3) MGIV (modified Givens) method:
!          calculates all eigenvalues and some eigenvactors of KLL, MLL (or KLL, KLLD for buckling). This method is only practical
!          for relatively small problems. It requires KLL to be a positive definite matrix. The algorithm performs a Cholesky
!          decomp of matrix KLL which can be time consuming for large problems

!   (4) INV (Inverse Power) method:
!          calculates only the lowest eigenvalue and its eigenvector of KLL, MLL (or KLL, KLLD for buckling).


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR, ERR, ERRSTAT, F06, L1M, L3A, SC1
      USE IOUNT1, ONLY                :  LINK1M,  LINK2I,  LINK3A, L1M_MSG, L3A_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, COMM, FATAL_ERR, INT_SC_NUM, LINKNO, MBUG, NDOFG, NDOFL, NSUB,            &
                                         NMPC, NRIGEL, NTERM_RMG,                                                                  &
                                         NTERM_KLL, NTERM_KLLD, NTERM_KLLDn,                                                       &
                                         NTERM_MLL, NTERM_MLLn,                                                                    &
                                         NVEC, NUM_EIGENS, NUM_KLLD_DIAG_ZEROS, NUM_MLL_DIAG_ZEROS, SOL_NAME, WARN_ERR,           &
                                         NUM_MODES_SUBS, NUM_BUCKLING_SUBS, TOTAL_MODES, MODE_SUBCASE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  EPSIL, SCRSPEC, SOLLIB, SPARSTOR, SUPINFO
      USE MODEL_STUF, ONLY            :  CC_EIGR_SID, EIG_PARAMS, IS_MODES_SUBCASE, IS_BUCKLING_SUBCASE, NUM_EIGENS_SUB,          &
                                         EIG_COMP, EIG_CRIT, EIG_FRQ1, EIG_FRQ2, EIG_GRID, EIG_METH, EIG_MSGLVL,                  &
                                         EIG_LANCZOS_NEV_DELT, EIG_LAP_MAT_TYPE, EIG_MODE, EIG_N1, EIG_N2, EIG_NCVFACL, EIG_NORM, &
                                         EIG_SID,                                                                                    &
                                         EIG_SIGMA, EIG_VECS, MAXMIJ, MIJ_COL, MIJ_ROW, NUM_FAIL_CRIT,                            &
                                         EIG_EXTRACT_METHOD, EIG_EXTRACT_MODE, EIG_EXTRACT_SOURCE, EIG_FEAST_M0,                  &
                                         EIG_FEAST_TOL_DIGITS, EIG_FEAST_MAX_LOOP, EIG_FEAST_N_CONTOUR,                           &
                                         EIG_SUBSPACE_NSUB, EIG_SUBSPACE_MAX_ITER, EIG_DENSE_NEX, EIG_FEAST_SEARCH_SCALE,         &
                                         EIG_SUBSPACE_TOL, SCNUM

      USE SPARSE_MATRICES, ONLY       :  I_KLL, J_KLL, KLL, I_KLLD, J_KLLD, KLLD, I_KLLDn, J_KLLDn, KLLDn,                         &
                                         I_MLL, J_MLL, MLL, I_MLLn, J_MLLn, MLLn
      USE OUTPUT4_MATRICES, ONLY      :  NUM_OU4_REQUESTS
      USE EIGEN_MATRICES_1, ONLY      :  GEN_MASS, MODE_NUM, EIGEN_VAL, EIGEN_VEC
      USE LAPACK_DPB_MATRICES, ONLY   :  ABAND, BBAND
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE EIGRL_EXTRACT_SOLVERS, ONLY :  EIG_LANCZOS_FEAST, EIG_LANCZOS_SUBSPACE, EIG_LANCZOS_DENSE
      USE RESPONSE_SPECTRA_STUF, ONLY :  RS_NUM_SUPORT

      USE LINK4_USE_IFs
      USE LINK_MESSAGE_Interface

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LINK4'

      INTEGER(LONG)                   :: I,J                 ! DO loop indices or counters.
      INTEGER(LONG)                   :: CANONICAL_ISUB      ! First modes subcase, or one matching CC_EIGR_SID
      INTEGER(LONG)                   :: CUR_ISUB            ! Subcase solved in current iteration
      INTEGER(LONG)                   :: IERROR              ! Error count when reading records from a file.
      INTEGER(LONG)                   :: IDX                 ! Running index when concatenating modal results
      INTEGER(LONG)                   :: ITER                ! Modes-subcase iteration counter
      INTEGER(LONG)                   :: KCNT                ! Modes-subcase counter
      INTEGER(LONG)                   :: KMODE               ! Per-subcase mode loop index
      INTEGER(LONG)                   :: N_MODES_ITER        ! Number of modal solves to perform
      INTEGER(LONG)                   :: NTERM_KLL_BAK       ! Snapshot of KLL nonzero count for multi-iter modal solves
      INTEGER(LONG)                   :: RSA_ARPACK_DOF_AVAIL ! Effective modal DOF count available to adaptive ARPACK
      INTEGER(LONG)                   :: TOTAL_MODES_LOCAL   ! Sum of NUM_EIGENS_SUB across resolved modes subcases
      INTEGER(LONG)                   :: NVEC_USED           ! Safe vector count limited by allocated EIGEN_VEC columns.
      INTEGER(LONG)                   :: OUNT(2)             ! File units to write messages to. Input to subr UNFORMATTED_OPEN.
      INTEGER(LONG), PARAMETER        :: P_LINKNO = 2        ! Prior LINK no's that should have run before this LINK can execute.
      INTEGER(LONG), ALLOCATABLE      :: I_KLL_BAK(:)        ! KLL row-pointer backup for multi-subcase modal solves
      INTEGER(LONG), ALLOCATABLE      :: J_KLL_BAK(:)        ! KLL column backup for multi-subcase modal solves
      INTEGER(LONG)                   :: TARGET_PRELOAD
      INTEGER(LONG)                   :: IERR_RELOAD
      INTEGER(LONG)                   :: CURRENT_PRELOAD_ISUB
      LOGICAL                         :: IS_BUCK_MULTI

      REAL(DOUBLE)                    :: EPS1                ! Small number to compare variables against zero.
      REAL(DOUBLE)                    :: EIGEN_VEC_COL(NDOFL)! One eigenvector put into a 1-D array.
      REAL(DOUBLE), ALLOCATABLE       :: KLL_BAK(:)          ! KLL value backup for multi-subcase modal solves
      LOGICAL                         :: WRITE_MLL           ! write the MLL matrix
! **********************************************************************************************************************************
      LINKNO = 4

      EPS1   = EPSIL(1)
      WRITE_MLL = (DEBUG(42) == 2)


      ! Set time initializing parameters
      CALL TIME_INIT

      ! Initialize WRT_BUG
      DO I=0,MBUG-1
         WRT_BUG(I) = 0
      ENDDO

      ! Get date and time, write to screen
      CALL OURDAT
      CALL OURTIM
      WRITE(SC1,152) LINKNO

      ! Make units for writing errors the screen until we open output files
      OUNT(1) = SC1
      OUNT(2) = SC1

      ! Make units for writing errors the error file and output file
      OUNT(1) = ERR
      OUNT(2) = F06

      ! Write info to text files
      WRITE(F06,150) LINKNO
      WRITE(ERR,150) LINKNO

      ! Read LINK1A file
      CALL READ_L1A ( 'KEEP' )

      ! Check COMM for successful completion of prior LINKs
      IF (COMM(P_LINKNO) /= 'C') THEN
         WRITE(ERR,9998) P_LINKNO,P_LINKNO,LINKNO
         WRITE(F06,9998) P_LINKNO,P_LINKNO,LINKNO
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Prior LINK's didn't complete, so quit
      ENDIF

      ! Make sure we have correct SOL
      IF ((SOL_NAME(1:5) /= 'MODES') .AND. (SOL_NAME(1:12) /= 'GEN CB MODEL') .AND. (SOL_NAME(1:8) /= 'BUCKLING')) THEN
         WRITE(ERR,999) 'MODES or BUCKLING or GEN CB MODEL', SOL_NAME
         WRITE(F06,999) 'MODES or BUCKLING or GEN CB MODEL', SOL_NAME
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
      ! Read data from file LINK1M
      CALL READ_L1M ( IERROR )

      IF (DEBUG(184) > 0) THEN
         WRITE(F06,*   ) ' Data written to file L1M'
         WRITE(F06,9102) '   EIG_SID         ', EIG_SID
         WRITE(F06,9101) '   EIG_METH        ', EIG_METH
         WRITE(F06,9103) '   EIG_FRQ1        ', EIG_FRQ1
         WRITE(F06,9103) '   EIG_FRQ2        ', EIG_FRQ2
         WRITE(F06,9102) '   EIG_N1          ', EIG_N1
         WRITE(F06,9102) '   EIG_N2          ', EIG_N2
         WRITE(F06,9101) '   EIG_VECS        ', EIG_VECS
         WRITE(F06,9103) '   EIG_CRIT        ', EIG_CRIT
         WRITE(F06,9101) '   EIG_NORM        ', EIG_NORM
         WRITE(F06,9102) '   EIG_GRID        ', EIG_GRID
         WRITE(F06,9102) '   EIG_COMP        ', EIG_COMP
         WRITE(F06,9102) '   EIG_MODE        ', EIG_MODE
         WRITE(F06,9103) '   EIG_SIGMA       ', EIG_SIGMA
         WRITE(F06,9101) '   EIG_LAP_MAT_TYPE', EIG_LAP_MAT_TYPE
         WRITE(F06,9102) '   EIG_MSGLVL      ', EIG_MSGLVL
         WRITE(F06,9102) '   EIG_NCVFACL     ', EIG_NCVFACL
         WRITE(F06,9101) '   EIG_EXTRACT_METH', EIG_EXTRACT_METHOD
         WRITE(F06,9101) '   EIG_EXTRACT_MODE', EIG_EXTRACT_MODE
         WRITE(F06,9101) '   EIG_EXTRACT_SRC ', EIG_EXTRACT_SOURCE
         WRITE(F06,9102) '   EIG_FEAST_M0    ', EIG_FEAST_M0
         WRITE(F06,9102) '   EIG_FEAST_TDIG  ', EIG_FEAST_TOL_DIGITS
         WRITE(F06,9102) '   EIG_FEAST_LOOP  ', EIG_FEAST_MAX_LOOP
         WRITE(F06,9102) '   EIG_FEAST_CNTR  ', EIG_FEAST_N_CONTOUR
         WRITE(F06,9103) '   EIG_FEAST_SCALE ', EIG_FEAST_SEARCH_SCALE
         WRITE(F06,9102) '   EIG_SUBSP_NSUB  ', EIG_SUBSPACE_NSUB
         WRITE(F06,9103) '   EIG_SUBSP_TOL   ', EIG_SUBSPACE_TOL
         WRITE(F06,9102) '   EIG_SUBSP_ITR   ', EIG_SUBSPACE_MAX_ITER
         WRITE(F06,9102) '   EIG_DENSE_NEX   ', EIG_DENSE_NEX
         WRITE(F06,9102) '   NUM_FAIL_CRIT   ', NUM_FAIL_CRIT
         WRITE(F06,9103) '   MAXMIJ          ', MAXMIJ
         WRITE(F06,9102) '   MIJ_ROW         ', MIJ_ROW
         WRITE(F06,9102) '   MIJ_COL         ', MIJ_COL
         WRITE(F06,*)
      ENDIF

      IF (IERROR > 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      ! NUM_MLL_DIAG_ZEROS will be used for a message written
      ! when the eigen summary is printed in subr EIG_SUMMARY
      ! (if more than this number of eigens are requested)
      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CONTINUE
      ELSE
         CALL SPARSE_MAT_DIAG_ZEROS ( 'MLL', NDOFL, NTERM_MLL, I_MLL, J_MLL, NUM_MLL_DIAG_ZEROS )
      ENDIF

! Generate nonsymmetric storage form for KLLD (if BUCKLING soln) or MLL. This is done since subr MATMULT_SFF, used when subr DSBAND
! is called herein, will run faster. MATMULT_SFF is called in each "Reverse commumication loop" in DSBAND for the LANCZOS method.

      IF      (SPARSTOR == 'SYM   ') THEN

         IF (SOL_NAME(1:8) == 'BUCKLING') THEN

            CALL SPARSE_MAT_DIAG_ZEROS ( 'KLLD', NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, NUM_KLLD_DIAG_ZEROS )
            NTERM_KLLDn = 2*NTERM_KLLD  - (NDOFL - NUM_KLLD_DIAG_ZEROS)

            CALL LINK_MESSAGE('ALLOCATE SPARSE KLLDn ARRAYS')
            CALL ALLOCATE_SPARSE_MAT ( 'KLLDn', NDOFL, NTERM_KLLDn, SUBR_NAME )

            CALL LINK_MESSAGE('CONVERT SYM CRS KLLD TO NONSYM CRS KLLDn')
            CALL CRS_SYM_TO_CRS_NONSYM ( 'KLLD', NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, KLLD, 'KLLDn', NTERM_KLLDn,                    &
                                         I_KLLDn, J_KLLDn, KLLDn, 'Y' )

         ELSE

            CALL SPARSE_MAT_DIAG_ZEROS ( 'MLL', NDOFL, NTERM_MLL, I_MLL, J_MLL, NUM_MLL_DIAG_ZEROS )
            NTERM_MLLn = 2*NTERM_MLL  - (NDOFL - NUM_MLL_DIAG_ZEROS)

            CALL LINK_MESSAGE('ALLOCATE SPARSE MLLn ARRAYS')
            CALL ALLOCATE_SPARSE_MAT ( 'MLLn', NDOFL, NTERM_MLLn, SUBR_NAME )

            CALL LINK_MESSAGE('CONVERT SYM CRS MLL TO NONSYM CRS MLLn')
            CALL CRS_SYM_TO_CRS_NONSYM ( 'MLL', NDOFL, NTERM_MLL, I_MLL, J_MLL, MLL, 'MLLn', NTERM_MLLn, I_MLLn, J_MLLn, MLLn, 'Y' )

         ENDIF

      ELSE IF (SPARSTOR == 'NONSYM') THEN

         IF (SOL_NAME(1:8) == 'BUCKLING') THEN

            CALL LINK_MESSAGE('ALLOCATE ARRAYS FOR NONSYM STORAGE OF KLLD')
            NTERM_KLLDn = NTERM_KLLD
            CALL ALLOCATE_SPARSE_MAT ( 'KLLDn', NDOFL, NTERM_KLLDn, SUBR_NAME )

            CALL LINK_MESSAGE('GET VALUES FOR NONSYM FORM OF KLLD')
            DO I=1,NDOFL+1
               I_KLLDn(I) = I_KLLD(I)
            ENDDO
            DO J=1,NTERM_KLLDn
               J_KLLDn(J) = J_KLLD(J)
                 KLLDn(J) =   KLLD(J)
            ENDDO

         ELSE

            CALL LINK_MESSAGE('ALLOCATE ARRAYS FOR NONSYM STORAGE OF MLL')
            NTERM_MLLn = NTERM_MLL
            CALL ALLOCATE_SPARSE_MAT ( 'MLLn', NDOFL, NTERM_MLLn, SUBR_NAME )

            CALL LINK_MESSAGE('GET VALUES FOR NONSYM FORM OF MLL')
            DO I=1,NDOFL+1
               I_MLLn(I) = I_MLL(I)
            ENDDO
            DO J=1,NTERM_MLLn
               J_MLLn(J) = J_MLL(J)
                 MLLn(J) =   MLL(J)
            ENDDO

         ENDIF

      ELSE
         !      Error - incorrect SPARSTOR
         WRITE(ERR,932) SUBR_NAME, SPARSTOR
         WRITE(F06,932) SUBR_NAME, SPARSTOR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF

      IF (WRITE_MLL) THEN
         CALL WRITE_SPARSE_CRS ( ' MLLn', 'A ', 'A ', NTERM_MLLn, NDOFL, I_MLLn, J_MLLn, MLLn )
      ENDIF


! **********************************************************************************************************************************
      NUM_MODES_SUBS = 0
      CANONICAL_ISUB = 0
      IS_BUCK_MULTI  = .FALSE.
      IF (((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:8) == 'BUCKLING')) .AND. ALLOCATED(IS_MODES_SUBCASE)) THEN
         DO I=1,NSUB
            IF (IS_MODES_SUBCASE(I) == 'Y') THEN
               NUM_MODES_SUBS = NUM_MODES_SUBS + 1
               IF ((CANONICAL_ISUB == 0) .AND. ALLOCATED(EIG_PARAMS)) THEN
                  IF (EIG_PARAMS(I)%SID == CC_EIGR_SID) CANONICAL_ISUB = I
               ENDIF
            ENDIF
         ENDDO
         IF (CANONICAL_ISUB == 0) THEN
            DO I=1,NSUB
               IF (IS_MODES_SUBCASE(I) == 'Y') THEN
                  CANONICAL_ISUB = I
                  EXIT
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      IF (CANONICAL_ISUB == 0) CANONICAL_ISUB = 1
      N_MODES_ITER = MAX(1, NUM_MODES_SUBS)
      IS_BUCK_MULTI = ((SOL_NAME(1:8) == 'BUCKLING') .AND. (NUM_BUCKLING_SUBS > 1))
      ! Start unknown so ITER=1 also corrects UG_COL/KLLD if a prior link left
      ! the preload state pointing at a different static subcase.
      CURRENT_PRELOAD_ISUB = 0

      IF (N_MODES_ITER > 1) THEN
         NTERM_KLL_BAK = NTERM_KLL
         ALLOCATE(I_KLL_BAK(NDOFL+1))
         ALLOCATE(J_KLL_BAK(NTERM_KLL))
         ALLOCATE(KLL_BAK(NTERM_KLL))
         I_KLL_BAK = I_KLL
         J_KLL_BAK = J_KLL
         KLL_BAK   = KLL
      ENDIF

m_lp: DO ITER = 1, N_MODES_ITER

         IF (ITER > 1) THEN
            NTERM_KLL = NTERM_KLL_BAK
            IF (.NOT. ALLOCATED(I_KLL)) ALLOCATE(I_KLL(NDOFL+1))
            IF (.NOT. ALLOCATED(J_KLL)) ALLOCATE(J_KLL(NTERM_KLL))
            IF (.NOT. ALLOCATED(KLL))   ALLOCATE(KLL(NTERM_KLL))
            I_KLL = I_KLL_BAK
            J_KLL = J_KLL_BAK
            KLL   = KLL_BAK
         ENDIF

         IF (IS_BUCK_MULTI) THEN
            KCNT = 0
            CUR_ISUB = CANONICAL_ISUB
            DO I=1,NSUB
               IF (IS_MODES_SUBCASE(I) == 'Y') THEN
                  KCNT = KCNT + 1
                  IF (KCNT == ITER) THEN
                     CUR_ISUB = I
                     EXIT
                  ENDIF
               ENDIF
            ENDDO
            TARGET_PRELOAD = 0
            IF (ALLOCATED(EIG_PARAMS)) TARGET_PRELOAD = EIG_PARAMS(CUR_ISUB)%STATSUB_REF
            IF ((TARGET_PRELOAD > 0) .AND. (TARGET_PRELOAD /= CURRENT_PRELOAD_ISUB)) THEN
               CALL LINK_MESSAGE('RELOAD UG_COL FROM L5A FOR NEXT PRELOAD')
               CALL DEALLOCATE_COL_VEC ( 'UG_COL' )
               CALL ALLOCATE_COL_VEC ( 'UG_COL', NDOFG, SUBR_NAME )
               IERR_RELOAD = 0
               CALL READ_L5A_UG_FOR_SUBCASE ( TARGET_PRELOAD, IERR_RELOAD )
               IF (IERR_RELOAD /= 0) THEN
                  WRITE(ERR,9994) SUBR_NAME, TARGET_PRELOAD, IERR_RELOAD
                  WRITE(F06,9994) SUBR_NAME, TARGET_PRELOAD, IERR_RELOAD
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
               CURRENT_PRELOAD_ISUB = TARGET_PRELOAD
               ! ITER=1 can also arrive here with a stale preload/KLLD pair from an
               ! earlier link pass, so rebuild immediately after correcting UG_COL.
               IF (ITER == 1) THEN
                  CALL LINK_MESSAGE('REBUILD KLLD FROM CORRECTED UG_COL (STATSUB ITER=1)')
                  IF (ALLOCATED(KLLDn) .OR. ALLOCATED(I_KLLDn) .OR. ALLOCATED(J_KLLDn)) THEN
                     CALL DEALLOCATE_SPARSE_MAT ( 'KLLDn' )
                  ENDIF
                  CALL REBUILD_KLLD_FROM_KGGD
                  IF      (SPARSTOR == 'SYM   ') THEN
                     CALL SPARSE_MAT_DIAG_ZEROS ( 'KLLD', NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, NUM_KLLD_DIAG_ZEROS )
                     NTERM_KLLDn = 2*NTERM_KLLD - (NDOFL - NUM_KLLD_DIAG_ZEROS)
                     CALL ALLOCATE_SPARSE_MAT ( 'KLLDn', NDOFL, NTERM_KLLDn, SUBR_NAME )
                     CALL CRS_SYM_TO_CRS_NONSYM ( 'KLLD', NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, KLLD, 'KLLDn', NTERM_KLLDn,          &
                                                  I_KLLDn, J_KLLDn, KLLDn, 'Y' )
                  ELSE IF (SPARSTOR == 'NONSYM') THEN
                     NTERM_KLLDn = NTERM_KLLD
                     CALL ALLOCATE_SPARSE_MAT ( 'KLLDn', NDOFL, NTERM_KLLDn, SUBR_NAME )
                     DO I=1,NDOFL+1
                        I_KLLDn(I) = I_KLLD(I)
                     ENDDO
                     DO J=1,NTERM_KLLDn
                        J_KLLDn(J) = J_KLLD(J)
                        KLLDn(J)   = KLLD(J)
                     ENDDO
                  ENDIF
               ENDIF
            ENDIF
            IF (ITER > 1) THEN
               CALL LINK_MESSAGE('REBUILD KLLD FROM CURRENT UG_COL (STATSUB)')
               CALL REBUILD_KLLD_FROM_KGGD
               IF      (SPARSTOR == 'SYM   ') THEN
                  CALL SPARSE_MAT_DIAG_ZEROS ( 'KLLD', NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, NUM_KLLD_DIAG_ZEROS )
                  NTERM_KLLDn = 2*NTERM_KLLD - (NDOFL - NUM_KLLD_DIAG_ZEROS)
                  CALL ALLOCATE_SPARSE_MAT ( 'KLLDn', NDOFL, NTERM_KLLDn, SUBR_NAME )
                  CALL CRS_SYM_TO_CRS_NONSYM ( 'KLLD', NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, KLLD, 'KLLDn', NTERM_KLLDn,             &
                                               I_KLLDn, J_KLLDn, KLLDn, 'Y' )
               ELSE IF (SPARSTOR == 'NONSYM') THEN
                  NTERM_KLLDn = NTERM_KLLD
                  CALL ALLOCATE_SPARSE_MAT ( 'KLLDn', NDOFL, NTERM_KLLDn, SUBR_NAME )
                  DO I=1,NDOFL+1
                     I_KLLDn(I) = I_KLLD(I)
                  ENDDO
                  DO J=1,NTERM_KLLDn
                     J_KLLDn(J) = J_KLLD(J)
                     KLLDn(J)   = KLLD(J)
                  ENDDO
               ENDIF
            ENDIF
         ENDIF

         IF (NUM_MODES_SUBS == 0) THEN
            CUR_ISUB = CANONICAL_ISUB
         ELSE
            KCNT = 0
            CUR_ISUB = CANONICAL_ISUB
            DO I=1,NSUB
               IF (IS_MODES_SUBCASE(I) == 'Y') THEN
                  KCNT = KCNT + 1
                  IF (KCNT == ITER) THEN
                     CUR_ISUB = I
                     EXIT
                  ENDIF
               ENDIF
            ENDDO

            EIG_SID               = EIG_PARAMS(CUR_ISUB)%SID
            EIG_METH              = EIG_PARAMS(CUR_ISUB)%METHOD
            EIG_NORM              = EIG_PARAMS(CUR_ISUB)%NORM
            EIG_GRID              = EIG_PARAMS(CUR_ISUB)%GRID
            EIG_COMP              = EIG_PARAMS(CUR_ISUB)%COMP
            EIG_FRQ1              = EIG_PARAMS(CUR_ISUB)%FRQ1
            EIG_FRQ2              = EIG_PARAMS(CUR_ISUB)%FRQ2
            EIG_N1                = EIG_PARAMS(CUR_ISUB)%N1
            EIG_N2                = EIG_PARAMS(CUR_ISUB)%N2
            EIG_NCVFACL           = EIG_PARAMS(CUR_ISUB)%NCVFACL
            EIG_MSGLVL            = EIG_PARAMS(CUR_ISUB)%MSGLVL
            EIG_MODE              = EIG_PARAMS(CUR_ISUB)%MODE
            EIG_VECS              = EIG_PARAMS(CUR_ISUB)%VECS
            EIG_CRIT              = EIG_PARAMS(CUR_ISUB)%CRIT
            EIG_SIGMA             = EIG_PARAMS(CUR_ISUB)%SIGMA
            EIG_LAP_MAT_TYPE      = EIG_PARAMS(CUR_ISUB)%LAP_MAT_TYPE
            EIG_LANCZOS_NEV_DELT  = EIG_PARAMS(CUR_ISUB)%LANCZOS_NEV_DELT
            EIG_EXTRACT_METHOD    = EIG_PARAMS(CUR_ISUB)%EXTRACT_METHOD
            EIG_EXTRACT_MODE      = EIG_PARAMS(CUR_ISUB)%EXTRACT_MODE
            EIG_EXTRACT_SOURCE    = EIG_PARAMS(CUR_ISUB)%EXTRACT_SOURCE
            EIG_FEAST_M0          = EIG_PARAMS(CUR_ISUB)%FEAST_M0
            EIG_FEAST_TOL_DIGITS  = EIG_PARAMS(CUR_ISUB)%FEAST_TOL_DIGITS
            EIG_FEAST_MAX_LOOP    = EIG_PARAMS(CUR_ISUB)%FEAST_MAX_LOOP
            EIG_FEAST_N_CONTOUR   = EIG_PARAMS(CUR_ISUB)%FEAST_N_CONTOUR
            EIG_FEAST_SEARCH_SCALE= EIG_PARAMS(CUR_ISUB)%FEAST_SEARCH_SCALE
            EIG_SUBSPACE_NSUB     = EIG_PARAMS(CUR_ISUB)%SUBSPACE_NSUB
            EIG_SUBSPACE_TOL      = EIG_PARAMS(CUR_ISUB)%SUBSPACE_TOL
            EIG_SUBSPACE_MAX_ITER = EIG_PARAMS(CUR_ISUB)%SUBSPACE_MAX_ITER
            EIG_DENSE_NEX         = EIG_PARAMS(CUR_ISUB)%DENSE_NEX
         ENDIF
         INT_SC_NUM = CUR_ISUB

         ! Solve eigenvalue problem
         IF ((EIG_METH(1:3) == 'GIV') .OR. (EIG_METH(1:4) == 'MGIV')) THEN
            CALL EIG_GIV_MGIV

         ELSE IF (EIG_METH(1:3) == 'INV') THEN
            CALL EIG_INV_PWR

         ELSE IF (EIG_METH(1:7) == 'LANCZOS') THEN
            IF (EIG_EXTRACT_METHOD(1:5) == 'FEAST') THEN
               CALL EIG_LANCZOS_FEAST
            ELSE IF (EIG_EXTRACT_METHOD(1:5) == 'SUBSP') THEN
               CALL EIG_LANCZOS_SUBSPACE
            ELSE IF (EIG_EXTRACT_METHOD(1:5) == 'DENSE') THEN
               CALL EIG_LANCZOS_DENSE
            ELSE
               IF ((EIG_FRQ2 > EPS1) .AND. (SOL_NAME(1:8) /= 'BUCKLING') .AND. (SOL_NAME(1:12) /= 'GEN CB MODEL')) THEN
                  RSA_ARPACK_DOF_AVAIL = NDOFL - NUM_MLL_DIAG_ZEROS
                  IF (RSA_ARPACK_DOF_AVAIL <= 4) THEN
                     IF ((SCRSPEC == 'Y') .AND. (RS_NUM_SUPORT == 1)) THEN
                        CALL EIG_LANCZOS_ARPACK
                     ELSE
                        WARN_ERR = WARN_ERR + 1
                        WRITE(ERR,4971) EIG_N2, NDOFL
                        WRITE(F06,4971) EIG_N2, NDOFL
                        CALL EIG_LANCZOS_DENSE
                     ENDIF
                  ELSE
                     CALL EIG_LANCZOS_ARPACK_ADAPTIVE
                  ENDIF
               ELSE
                  CALL EIG_LANCZOS_ARPACK
               ENDIF
            ENDIF

         ELSE

            WRITE(ERR,4005) SUBR_NAME, EIG_METH
            WRITE(F06,4005) SUBR_NAME, EIG_METH
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )

         ENDIF

         IF (SOL_NAME(1:12) /= 'GEN CB MODEL') THEN
            IF (SOL_NAME(1:8) == 'BUCKLING') THEN
               WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KLLD', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KLLD' )
            ENDIF
         ENDIF

         NUM_FAIL_CRIT = 0
         MAXMIJ        = 0
         MIJ_ROW       = 0
         MIJ_COL       = 0

         CALL ALLOCATE_EIGEN1_MAT ( 'GEN_MASS', NUM_EIGENS, 1, SUBR_NAME )

         IF (NVEC > 0) THEN
            CALL LINK_MESSAGE('CALCULATE GENERALIZED MASS')
            CALL CALC_GEN_MASS

            IF (EIG_NORM == 'MASS') THEN
               CALL LINK_MESSAGE('RENORMALIZE EIGENVECTORS TO UNIT GEN MASS')
               CALL RENORM_ON_MASS ( NVEC, EPS1 )
            ENDIF
         ELSE
            DO I=1,NUM_EIGENS
               GEN_MASS(I) = ZERO
            ENDDO
         ENDIF

         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KLLDn', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KLLDn' )
         ENDIF

         IF (NUM_MODES_SUBS > 0) THEN
            IF (ALLOCATED(EIG_PARAMS(CUR_ISUB)%EIGEN_VAL)) DEALLOCATE(EIG_PARAMS(CUR_ISUB)%EIGEN_VAL)
            IF (ALLOCATED(EIG_PARAMS(CUR_ISUB)%MODE_NUM )) DEALLOCATE(EIG_PARAMS(CUR_ISUB)%MODE_NUM)
            IF (ALLOCATED(EIG_PARAMS(CUR_ISUB)%GEN_MASS )) DEALLOCATE(EIG_PARAMS(CUR_ISUB)%GEN_MASS)
            IF (ALLOCATED(EIG_PARAMS(CUR_ISUB)%EIGEN_VEC)) DEALLOCATE(EIG_PARAMS(CUR_ISUB)%EIGEN_VEC)
            ALLOCATE(EIG_PARAMS(CUR_ISUB)%EIGEN_VAL(NUM_EIGENS))
            ALLOCATE(EIG_PARAMS(CUR_ISUB)%MODE_NUM(NUM_EIGENS))
            ALLOCATE(EIG_PARAMS(CUR_ISUB)%GEN_MASS(NUM_EIGENS))
            ALLOCATE(EIG_PARAMS(CUR_ISUB)%EIGEN_VEC(NDOFL,MAX(1,NVEC)))
            EIG_PARAMS(CUR_ISUB)%EIGEN_VAL(1:NUM_EIGENS) = EIGEN_VAL(1:NUM_EIGENS)
            EIG_PARAMS(CUR_ISUB)%MODE_NUM(1:NUM_EIGENS)  = MODE_NUM(1:NUM_EIGENS)
            EIG_PARAMS(CUR_ISUB)%GEN_MASS(1:NUM_EIGENS)  = GEN_MASS(1:NUM_EIGENS)
            IF (NVEC > 0) THEN
               EIG_PARAMS(CUR_ISUB)%EIGEN_VEC(1:NDOFL,1:NVEC) = EIGEN_VEC(1:NDOFL,1:NVEC)
            ENDIF
            EIG_PARAMS(CUR_ISUB)%NUM_EIGENS    = NUM_EIGENS
            EIG_PARAMS(CUR_ISUB)%NVEC          = NVEC
            EIG_PARAMS(CUR_ISUB)%NUM_FAIL_CRIT = NUM_FAIL_CRIT
            EIG_PARAMS(CUR_ISUB)%MAXMIJ        = MAXMIJ
            EIG_PARAMS(CUR_ISUB)%MIJ_ROW       = MIJ_ROW
            EIG_PARAMS(CUR_ISUB)%MIJ_COL       = MIJ_COL
            NUM_EIGENS_SUB(CUR_ISUB)           = NUM_EIGENS
         ENDIF

         IF ((EIG_NORM == 'MASS    ') .OR. (EIG_NORM == 'NONE')) THEN
            CALL LINK_MESSAGE('WRITE EIGENVALUE SUMMARY TO OUTFIL')
            CALL EIG_SUMMARY ( CUR_ISUB )
         ENDIF

         IF (ITER < N_MODES_ITER) THEN
            CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL' )
            CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC' )
            CALL DEALLOCATE_EIGEN1_MAT ( 'MODE_NUM'  )
            CALL DEALLOCATE_EIGEN1_MAT ( 'GEN_MASS'  )
            CALL DEALLOCATE_LAPACK_MAT ( 'ABAND' )
            CALL DEALLOCATE_LAPACK_MAT ( 'BBAND' )
            CALL DEALLOCATE_LAPACK_MAT ( 'RFAC'  )
         ENDIF

      ENDDO m_lp

      IF ((SOL_NAME(1:5) == 'MODES')) THEN
         IF (ALLOCATED(MLL))  THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MLL ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MLL' )
         ENDIF
         IF (ALLOCATED(MLLn)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MLLn ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MLLn' )
         ENDIF
      ENDIF

      IF (ALLOCATED(I_KLL_BAK)) DEALLOCATE(I_KLL_BAK)
      IF (ALLOCATED(J_KLL_BAK)) DEALLOCATE(J_KLL_BAK)
      IF (ALLOCATED(KLL_BAK  )) DEALLOCATE(KLL_BAK)

      IF (IS_BUCK_MULTI) THEN
         CALL DEALLOCATE_MODEL_STUF ( 'MPC_IND_GRIDS' )
         CALL DEALLOCATE_MODEL_STUF ( 'SINGLE ELEMENT ARRAYS' )
         CALL DEALLOCATE_MODEL_STUF ( 'SUBLOD' )
         CALL DEALLOCATE_COL_VEC ( 'UG_COL' )
      ENDIF

      IF (NUM_MODES_SUBS > 1) THEN
         CALL LINK_MESSAGE('CONCATENATE MULTI-SUBCASE MODE RESULTS   ')
         TOTAL_MODES_LOCAL = 0
         DO I=1,NSUB
            IF (IS_MODES_SUBCASE(I) == 'Y') TOTAL_MODES_LOCAL = TOTAL_MODES_LOCAL + NUM_EIGENS_SUB(I)
         ENDDO
         TOTAL_MODES = TOTAL_MODES_LOCAL

         CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL' )
         CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC' )
         CALL DEALLOCATE_EIGEN1_MAT ( 'MODE_NUM'  )
         CALL DEALLOCATE_EIGEN1_MAT ( 'GEN_MASS'  )
         CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL', TOTAL_MODES_LOCAL, 1, SUBR_NAME )
         CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC', NDOFL, TOTAL_MODES_LOCAL, SUBR_NAME )
         CALL ALLOCATE_EIGEN1_MAT ( 'MODE_NUM' , TOTAL_MODES_LOCAL, 1, SUBR_NAME )
         CALL ALLOCATE_EIGEN1_MAT ( 'GEN_MASS' , TOTAL_MODES_LOCAL, 1, SUBR_NAME )

         IF (ALLOCATED(MODE_SUBCASE)) DEALLOCATE(MODE_SUBCASE)
         ALLOCATE(MODE_SUBCASE(TOTAL_MODES_LOCAL))

         IDX = 0
         DO I=1,NSUB
            IF (IS_MODES_SUBCASE(I) /= 'Y') CYCLE
            DO KMODE = 1, NUM_EIGENS_SUB(I)
               IDX = IDX + 1
               EIGEN_VAL(IDX)          = EIG_PARAMS(I)%EIGEN_VAL(KMODE)
               MODE_NUM(IDX)           = IDX
               GEN_MASS(IDX)           = EIG_PARAMS(I)%GEN_MASS(KMODE)
               EIGEN_VEC(1:NDOFL,IDX)  = EIG_PARAMS(I)%EIGEN_VEC(1:NDOFL,KMODE)
               MODE_SUBCASE(IDX)       = I
            ENDDO
         ENDDO
         NUM_EIGENS = TOTAL_MODES_LOCAL
         NVEC       = TOTAL_MODES_LOCAL

      ELSE
         IF (ALLOCATED(MODE_SUBCASE)) DEALLOCATE(MODE_SUBCASE)
         ALLOCATE(MODE_SUBCASE(MAX(1,NUM_EIGENS)))
         MODE_SUBCASE = CANONICAL_ISUB
         TOTAL_MODES  = NUM_EIGENS
      ENDIF

      IF (NUM_MODES_SUBS > 0) THEN
         EIG_SID          = EIG_PARAMS(CANONICAL_ISUB)%SID
         EIG_METH         = EIG_PARAMS(CANONICAL_ISUB)%METHOD
         EIG_NORM         = EIG_PARAMS(CANONICAL_ISUB)%NORM
         EIG_GRID         = EIG_PARAMS(CANONICAL_ISUB)%GRID
         EIG_COMP         = EIG_PARAMS(CANONICAL_ISUB)%COMP
         EIG_FRQ1         = EIG_PARAMS(CANONICAL_ISUB)%FRQ1
         EIG_FRQ2         = EIG_PARAMS(CANONICAL_ISUB)%FRQ2
         EIG_N1           = EIG_PARAMS(CANONICAL_ISUB)%N1
         EIG_N2           = EIG_PARAMS(CANONICAL_ISUB)%N2
         EIG_NCVFACL      = EIG_PARAMS(CANONICAL_ISUB)%NCVFACL
         EIG_MSGLVL       = EIG_PARAMS(CANONICAL_ISUB)%MSGLVL
         EIG_MODE         = EIG_PARAMS(CANONICAL_ISUB)%MODE
         EIG_VECS         = EIG_PARAMS(CANONICAL_ISUB)%VECS
         EIG_CRIT         = EIG_PARAMS(CANONICAL_ISUB)%CRIT
         EIG_SIGMA        = EIG_PARAMS(CANONICAL_ISUB)%SIGMA
         EIG_LAP_MAT_TYPE = EIG_PARAMS(CANONICAL_ISUB)%LAP_MAT_TYPE
         EIG_EXTRACT_METHOD     = EIG_PARAMS(CANONICAL_ISUB)%EXTRACT_METHOD
         EIG_EXTRACT_MODE       = EIG_PARAMS(CANONICAL_ISUB)%EXTRACT_MODE
         EIG_EXTRACT_SOURCE     = EIG_PARAMS(CANONICAL_ISUB)%EXTRACT_SOURCE
         EIG_FEAST_M0           = EIG_PARAMS(CANONICAL_ISUB)%FEAST_M0
         EIG_FEAST_TOL_DIGITS   = EIG_PARAMS(CANONICAL_ISUB)%FEAST_TOL_DIGITS
         EIG_FEAST_MAX_LOOP     = EIG_PARAMS(CANONICAL_ISUB)%FEAST_MAX_LOOP
         EIG_FEAST_N_CONTOUR    = EIG_PARAMS(CANONICAL_ISUB)%FEAST_N_CONTOUR
         EIG_FEAST_SEARCH_SCALE = EIG_PARAMS(CANONICAL_ISUB)%FEAST_SEARCH_SCALE
         EIG_SUBSPACE_NSUB      = EIG_PARAMS(CANONICAL_ISUB)%SUBSPACE_NSUB
         EIG_SUBSPACE_TOL       = EIG_PARAMS(CANONICAL_ISUB)%SUBSPACE_TOL
         EIG_SUBSPACE_MAX_ITER  = EIG_PARAMS(CANONICAL_ISUB)%SUBSPACE_MAX_ITER
         EIG_DENSE_NEX          = EIG_PARAMS(CANONICAL_ISUB)%DENSE_NEX
         NUM_FAIL_CRIT    = EIG_PARAMS(CANONICAL_ISUB)%NUM_FAIL_CRIT
         MAXMIJ           = EIG_PARAMS(CANONICAL_ISUB)%MAXMIJ
         MIJ_ROW          = EIG_PARAMS(CANONICAL_ISUB)%MIJ_ROW
         MIJ_COL          = EIG_PARAMS(CANONICAL_ISUB)%MIJ_COL
      ENDIF

      CALL LINK_MESSAGE('WRITE MERGED EIGENVALUE DATA TO L1M        ')
      CALL WRITE_L1M

      ! Open and set up file L3A (used to hold eigenvectors)
      CALL FILE_OPEN ( L3A, LINK3A, OUNT, 'REPLACE', L3A_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

! --- arpack_surgery begin --- !
      NVEC_USED = MIN( NVEC, NUM_EIGENS )
      IF (ALLOCATED(EIGEN_VEC)) THEN
         NVEC_USED = MIN( NVEC_USED, SIZE(EIGEN_VEC,2) )
      ELSE
         NVEC_USED = 0
      ENDIF
! --- arpack_surgery end --- !

      ! Write out computed eigenvectors to L3A
      CALL LINK_MESSAGE('WRITE EIGENVECTORS TO DISK FILE')
      DO J=1,NVEC_USED
         DO I=1,NDOFL
           WRITE(L3A) EIGEN_VEC(I,J)
         ENDDO
      ENDDO
      CALL FILE_CLOSE ( L3A, LINK3A, 'KEEP' )

      ! Optional eigenvector debug output
      IF (DEBUG(43) == 1) THEN
         DO J=1,NVEC_USED
            DO I=1,NDOFL
               EIGEN_VEC_COL(I) = EIGEN_VEC(I,J)
            ENDDO
            WRITE(F06,'(//,1X,''EIGENVECTOR'',I8/)') J
            CALL WRITE_VECTOR ('    A-SET EIGENVECTOR   ','DISPL',NDOFL, EIGEN_VEC_COL )
         ENDDO
      ENDIF

      ! Call OUTPUT4 processor to process output requests for OUTPUT4 matrices generated in this link
      IF (NUM_OU4_REQUESTS > 0) THEN
         CALL LINK_MESSAGE('WRITE OUTPUT4 NATRICES      ')
         WRITE(F06,*)
         CALL OUTPUT4_PROC ( SUBR_NAME )
      ENDIF

      ! Deallocate arrays
      CALL DEALLOCATE_LAPACK_MAT ( 'RFAC' )

      ! leave EIGEN_VAL until LINK9 since it may be needed there
!xx   CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL' )
      CALL DEALLOCATE_EIGEN1_MAT ( 'GEN_MASS' )
      CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC' )
      CALL DEALLOCATE_EIGEN1_MAT ( 'MODE_NUM' )

      CALL DEALLOCATE_LAPACK_MAT ( 'ABAND' )
      CALL DEALLOCATE_LAPACK_MAT ( 'BBAND' )

      ! Process is now complete so set COMM(LINKNO)
      COMM(LINKNO) = 'C'

      ! Write data to L1A
      CALL WRITE_L1A ( 'KEEP', 'Y' )

      ! Check allocation status of allocatable arrays, if requested
      IF (DEBUG(100) > 0) THEN
         CALL CHK_ARRAY_ALLOC_STAT
         IF (DEBUG(100) > 1) THEN
            CALL WRITE_ALLOC_MEM_TABLE ( 'at the end of '//SUBR_NAME )
         ENDIF
      ENDIF

      ! Write LINK4 end to F06
      CALL OURTIM
      WRITE(F06,151) LINKNO

      ! Close files
      IF (( DEBUG(193) == 4) .OR. (DEBUG(193) == 999)) THEN
         CALL FILE_INQUIRE ( 'near end of LINK4' )
      ENDIF

      ! Write LINK4 end to screen
      WRITE(SC1,153) LINKNO

      RETURN

! **********************************************************************************************************************************
  150 FORMAT(/,' >> LINK',I3,' BEGIN',/)

  151 FORMAT(/,' >> LINK',I3,' END',/)

  152 FORMAT(/,' >> LINK',I3,' BEGIN')

  153 FORMAT(  ' >> LINK',I3,' END')

  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

  999 FORMAT(' *ERROR   999: INCORRECT SOLUTION IN EXEC CONTROL. SHOULD BE "',A,'", BUT IS "',A,'"')

 4005 FORMAT(' *ERROR  4005: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CODE ONLY WRITTEN FOR METHOD = GIV, MGIV, OR LANCZOS BUT METHOD IS = ',A8)

 9101 FORMAT(1X,A,' =  ','"',A,'"')

  9102 FORMAT(1X,A,' =  ',I13)

  9103 FORMAT(1X,A,' =  ',1ES13.6)

 4971 FORMAT(' *WARNING    : SMALL MODAL MODEL DETECTED (REQUESTED MODES =',I8,', NDOFL =',I8,').',                           &
                    /,15X,'LINK4 WILL USE THE CONDENSED DENSE REFERENCE SOLVER INSTEAD OF ARPACK LANCZOS.')

 9994 FORMAT(' *ERROR  9994: SUBROUTINE ',A,' FAILED TO RELOAD UG_COL FROM FILE L5A FOR PRELOAD SUBCASE ',I8,                 &
                    ' (IOSTAT = ',I8,').')

 9998 FORMAT(' *ERROR  9998: COMM ',I3,' INDICATES UNSUCCESSFUL LINK ',I2,' COMPLETION.'                                           &
                    ,/,14X,' FATAL ERROR - CANNOT START LINK ',I2)

12345 FORMAT(A,10X,A)

99001 FORMAT(1X,6(1ES14.6))

! **********************************************************************************************************************************

      END SUBROUTINE LINK4
