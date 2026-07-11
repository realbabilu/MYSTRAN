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

      SUBROUTINE EIG_LANCZOS_ARPACK

! Solves for eigenvalues and eigenvectors when the Lanczos method is requested (on Bulk Data EIGRL entry)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, KMSM_SDIA, LINKNO, NDOFL, NTERM_KLL, NTERM_KLLD, NTERM_KMSM,     &
                                         NTERM_KMSMn, NTERM_KMSMs, NTERM_MLL, NTERM_ULL, NUM_EIGENS, NUM_KLLD_DIAG_ZEROS,          &
                                         NUM_MLL_DIAG_ZEROS, NVEC, SOL_NAME, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, PI
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  ARP_TOL, BAILOUT, DARPACK, EIGESTL, EPSIL, MXITERL, SOLLIB, SPARSE_FLAVOR, SPARSTOR,     &
                                         SUPINFO, SUPWARN
      USE DOF_TABLES, ONLY            :  TDOFI
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL, EIGEN_VEC, MODE_NUM
      USE MODEL_STUF, ONLY            :  EIG_FRQ1, EIG_FRQ2, EIG_LANCZOS_NEV_DELT, EIG_LAP_MAT_TYPE, EIG_MODE, EIG_N1, EIG_N2,     &
                                         EIG_NCVFACL, EIG_SIGMA
      USE ARPACK_MATRICES_1, ONLY     :  IWORK, RESID, RFAC, SELECT, VBAS, WORKD, WORKL
      USE SPARSE_MATRICES, ONLY       :  I_KLL, J_KLL, KLL, I_KLLD, J_KLLD, KLLD, I_MLL, J_MLL, MLL, SYM_KLL, SYM_KLLD, SYM_MLL,   &
                                         I_KMSM, J_KMSM, KMSM, I_KMSMn, J_KMSMn, KMSMn, I_KMSMs, J_KMSMs, KMSMs

      USE ARPACK_LANCZOS_EIG
! --- MUMPS_COO add begin --- !
      USE DMUMPS_STUF, ONLY           :  DMUMPS_FREE_FACTORS
! --- MUMPS_COO add end --- !

      USE EIG_LANCZOS_ARPACK_USE_IFs
      USE LINK_MESSAGE_Interface

      IMPLICIT NONE

      LOGICAL                         :: RVEC              ! = .TRUE. or .FALSE. Specifies whether eigenvectors are to be calculated
! --- arpack_surgery begin --- !
      LOGICAL                         :: EST_NUM_EIGENS_VALID
      LOGICAL                         :: EST_NUM_EIGENS_VALID1
      LOGICAL                         :: EST_NUM_EIGENS_VALID2
! --- arpack_surgery end --- !

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EIG_LANCZOS_ARPACK'
      CHARACTER( 1*BYTE)              :: BMAT              ! BMAT specifies the type of the matrix B that defines the semi-inner
!                                                            product for the ARPACK operator OP (see ARPACK subr DSBAND):
!                                                            = 'I' -> standard eigenvalue problem:    A*x = lambda*x
!                                                            = 'G' -> generalized eigenvalue problem: A*x = lambda*M*x
      CHARACTER( 1*BYTE)              :: HOWMNY            ! How many eigenvecs are wanted and the form of the eigenvec matrix Z
!                                                            = 'A': compute all eigenvecs
!                                                            = 'S': compute some of the eigenvecs, specified by logical array SELECT
!                                                                   (see documentation for ARPACK subr DSBAND)
      CHARACTER( 2*BYTE)              :: WHICH             ! When IPARAM(7)= 1 or 2,  WHICH can be set to any one of the following:
!                                                            = 'LM' -> want the NEV eigenvalues of largest magnitude.
!                                                            = 'SM' -> want the NEV eigenvalues of smallest magnitude.
!                                                            = 'LA' -> want the NEV eigenvalues of largest REAL part.
!                                                            = 'SA' -> want the NEV eigenvalues of smallest REAL part.
!                                                            = 'BE' -> Compute  NEV eigenvalues, half from each end of the spectrum
!                                                                      If NEV is odd, compute 1 more from high end than from low end
!                                                            When IPARAM(7) = 3, 4, or 5,  WHICH should be set to 'LM' only.

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: CALLED_SUBR = ' ' ! Name of a called subr (for output error purposes)

      INTEGER(LONG)                   :: COMPV             ! Component number (1-6) of a grid DOF
      INTEGER(LONG)                   :: GRIDV             ! Grid number
      INTEGER(LONG)                   :: I,J               ! DO loop index.
      INTEGER(LONG)                   :: INFO_ARPACK = 0   ! Error msg number from ARPACK routine to do eigenvalue/vec calc
      INTEGER(LONG)                   :: INFO_LAPACK = 0   ! Error msg number from ARPACK routine to do eigenvalue/vec calc
      INTEGER(LONG)                   :: IPARAM(11)        ! Integer array input to subr DSBAND (in module ARPACK_LANCZOS_1)
      INTEGER(LONG)                   :: KL                ! Max(number of subdiags of KLL, number of subdiags of MLL or KLLD)
      INTEGER(LONG)                   :: KU                ! Max(number of superdiags of KLL, number of superdiags of MLL or KLLD)
      INTEGER(LONG)                   :: LDRFAC            ! Leading dimension of RFAC (band form of KLL - EIG_SIGMA*MLL, or KLLD)
      INTEGER(LONG)                   :: LWORKL            ! Used to dimension a work array
      INTEGER(LONG)                   :: NEV               ! Number of eigenvalues to find
      INTEGER(LONG)                   :: NCV               ! Number of columns of the matrix V. Represents the dimension of the
!                                                            Lanczos basis constructed by ARPACK subr dsaupd
      INTEGER(LONG)                   :: NUM1              ! Number to use for max no. of eigens to find. Must be NUM1 <= NDOFL
      INTEGER(LONG)                   :: NUM_NEG_TERMS1    ! Number of negative terms on the diagonal of RFAC for EIG_FRQ1
      INTEGER(LONG)                   :: NUM_NEG_TERMS2    ! Number of negative terms on the diagonal of RFAC for EIG_FRQ2
      INTEGER(LONG)                   :: NUM_EST_EIGENS    ! Number of estimated eigens in the freq interval (EIG_FRQ2 - EIG_FRQ1)
      INTEGER(LONG)                   :: NUM_KMSM_DIAG_0   ! Number of zero diagonal terms in KMSM
      INTEGER(LONG)                   :: NVEC_SAFE         ! Safe converged-vector count bounded by allocated eigen arrays

! --- arpack_surgery begin --- !
      INTEGER(LONG)                   :: MIN_NCV, MAX_NCV, LNONZEROS, NDOFL_EFFECTIVE, NEV_EXTRA
! --- arpack_surgery end --- !

      REAL(DOUBLE)                    :: EPS1              ! A small number to compare zero to

      INTRINSIC                       :: MIN



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

      NUM_EST_EIGENS = 0
      NUM_NEG_TERMS1 = 0
      NUM_NEG_TERMS2 = 0
! --- arpack_surgery begin --- !
      EST_NUM_EIGENS_VALID  = .FALSE.
      EST_NUM_EIGENS_VALID1 = .TRUE.
      EST_NUM_EIGENS_VALID2 = .FALSE.
! --- arpack_surgery end --- !

!xx   IF      (SOLLIB == 'BANDED  ') THEN

         IF (NDOFL < EIGESTL) THEN
            IF (EIG_FRQ2 > EPS1) THEN
               IF (EIG_FRQ1 > EPS1) THEN
! --- arpack_surgery begin --- !
                  CALL EST_NUM_EIGENS_BANDED ( EIG_FRQ1, NUM_NEG_TERMS1, EST_NUM_EIGENS_VALID1 )
! --- arpack_surgery end --- !
               ENDIF
! --- arpack_surgery begin --- !
               CALL EST_NUM_EIGENS_BANDED ( EIG_FRQ2, NUM_NEG_TERMS2, EST_NUM_EIGENS_VALID2 )
               EST_NUM_EIGENS_VALID = EST_NUM_EIGENS_VALID1 .AND. EST_NUM_EIGENS_VALID2
               IF (EST_NUM_EIGENS_VALID .AND. (NUM_NEG_TERMS2 > 0)) THEN
! --- arpack_surgery end --- !
                  NUM_EST_EIGENS = NUM_NEG_TERMS2 - NUM_NEG_TERMS1
               ELSE
                  ! ?
               ENDIF
            ENDIF
         ENDIF

!xx   ELSE

!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,9991) SUBR_NAME, 'SOLLIB'
!xx      WRITE(F06,9991) SUBR_NAME, 'SOLLIB'
!xx      CALL OUTA_HERE ( 'Y' )

!xx   ENDIF

      IF (EIG_FRQ2 > EPS1) THEN
         WRITE(SC1,101) EIG_FRQ2, NUM_EST_EIGENS
         WRITE(F06,101) EIG_FRQ2, NUM_EST_EIGENS
! --- arpack_surgery begin --- !
         IF (.NOT. EST_NUM_EIGENS_VALID) THEN
            WRITE(SC1,'(A)') '       Lanczos eigen-count estimate unavailable; using fallback NEV logic'
            WRITE(F06,'(A)') '       Lanczos eigen-count estimate unavailable; using fallback NEV logic'
         ENDIF
! --- arpack_surgery end --- !
      ENDIF

! Calc KMSM = KLL - EIG_SIGMA*MLL (or + EIG_SIGMA*KLLD for BUCKLING) where EIG_SIGMA = shift freq

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CALL MATADD_SSS_NTERM ( NDOFL, 'KLL' , NTERM_KLL , I_KLL , J_KLL , SYM_KLL ,  'eig_sigma*KLLD',                           &
                                                NTERM_KLLD, I_KLLD, J_KLLD, SYM_KLLD, 'KMSM', NTERM_KMSM )
         CALL ALLOCATE_SPARSE_MAT ( 'KMSM', NDOFL, NTERM_KMSM, SUBR_NAME )
         CALL MATADD_SSS       ( NDOFL, 'KLL' , NTERM_KLL , I_KLL , J_KLL , KLL , ONE, 'eig_sigma*KLLD',                           &
                                                NTERM_KLLD, I_KLLD, J_KLLD, KLLD, EIG_SIGMA,                                       &
                                        'KMSM', NTERM_KMSM, I_KMSM, J_KMSM, KMSM )
      ELSE
         CALL MATADD_SSS_NTERM ( NDOFL, 'KLL',  NTERM_KLL , I_KLL , J_KLL , SYM_KLL ,  '-eig_sigma*MLL',                           &
                                                NTERM_MLL , I_MLL , J_MLL , SYM_MLL , 'KMSM', NTERM_KMSM )
         CALL ALLOCATE_SPARSE_MAT ( 'KMSM', NDOFL, NTERM_KMSM, SUBR_NAME )
         CALL MATADD_SSS       ( NDOFL, 'KLL' , NTERM_KLL , I_KLL , J_KLL , KLL , ONE, '-eig_sigma*MLL',                           &
                                                NTERM_MLL , I_MLL , J_MLL , MLL, -EIG_SIGMA,                                       &
                                        'KMSM', NTERM_KMSM, I_KMSM, J_KMSM, KMSM )
      ENDIF

! --- BANDED_optimizisation -begin-- !
      IF (DEBUG(49) > 0) THEN
         CALL WRITE_SPARSE_CRS ( ' KLL input to Lanczos/ARPACK', 'A ', 'A ', NTERM_KLL, NDOFL, I_KLL, J_KLL, KLL )
         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            CALL WRITE_SPARSE_CRS ( ' KLLD input to Lanczos/ARPACK', 'A ', 'A ', NTERM_KLLD, NDOFL, I_KLLD, J_KLLD, KLLD )
            CALL WRITE_SPARSE_CRS ( ' KMSM = KLL + sigma*KLLD input to Lanczos/ARPACK', 'A ', 'A ', NTERM_KMSM, NDOFL, I_KMSM,   &
                                     J_KMSM, KMSM )
         ELSE
            CALL WRITE_SPARSE_CRS ( ' MLL input to Lanczos/ARPACK', 'A ', 'A ', NTERM_MLL, NDOFL, I_MLL, J_MLL, MLL )
            CALL WRITE_SPARSE_CRS ( ' KMSM = KLL - sigma*MLL input to Lanczos/ARPACK', 'A ', 'A ', NTERM_KMSM, NDOFL, I_KMSM,    &
                                     J_KMSM, KMSM )
         ENDIF
      ENDIF
! --- BANDED_optimizisation -end-- !


! Det bandwidth of KMSM so BANDGEN can put it in LAPACK band form. KMSM_SDIA is the number of super-diags in the band form of KMSM

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CALL LINK_MESSAGE('CALCULATE BANDWIDTH OF [KLL + sigma*KLLD]')
      ELSE
         CALL LINK_MESSAGE('CALCULATE BANDWIDTH OF [KLL - sigma*MLL]')
      ENDIF
      CALL BANDSIZ ( NDOFL, NTERM_KMSM, I_KMSM, J_KMSM, KMSM_SDIA )
      WRITE(ERR,4905) KMSM_SDIA
      IF (SUPINFO == 'N') THEN
         WRITE(F06,4905) KMSM_SDIA
      ENDIF
! --- BANDED_optimizisation -begin-- !
      CALL REPORT_BANDED_STORAGE_ESTIMATE ( 'KMSM', NDOFL, NTERM_KMSM, I_KMSM, J_KMSM, KMSM_SDIA+1, SUPINFO )
! --- BANDED_optimizisation -end-- !

! EIG_LAP_MAT_TYPE was checked in BD_EIGRL for correctness, but make sure, here, that it is correct

      IF(SOLLIB(1:6) == 'SPARSE') THEN
         LDRFAC = 1 ! make it small because it's not used but parameter has to be allocated.
      ELSEIF      (EIG_LAP_MAT_TYPE(1:3) == 'DPB') THEN
         LDRFAC = KMSM_SDIA + 1
      ELSE IF (EIG_LAP_MAT_TYPE(1:3) == 'DGB') THEN
         LDRFAC = 3*KMSM_SDIA + 1
      ELSE
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,4003) SUBR_NAME, EIG_LAP_MAT_TYPE
         WRITE(F06,4003) SUBR_NAME, EIG_LAP_MAT_TYPE
         CALL OUTA_HERE ( 'Y' )
      ENDIF


! Allocate array RFAC = (KLL - EIG_SIGMA*MLL, or KLL + EIG_SIGMA*KLLD) for ARACK

! --- BANDED_optimizisation -begin-- !
      CALL REPORT_SOLVER_DISPATCH_POLICY ( 'KMSM', SUBR_NAME )
! --- BANDED_optimizisation -end-- !
      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CALL LINK_MESSAGE('ALLOCATE ARPACK BAND MAT: RFAC = KLL + sigma*KLLD')
      ELSE
         CALL LINK_MESSAGE('ALLOCATE ARPACK BAND MAT: RFAC = KLL - sigma*MLL')
      ENDIF
      CALL ALLOCATE_LAPACK_MAT ( 'RFAC', LDRFAC, NDOFL, SUBR_NAME )


      IF(SOLLIB(1:6) == 'BANDED') THEN

! Put KMSM in form required by LAPACK band matrix. Call result array RFAC

         CALL LINK_MESSAGE('PUT RFAC MATRIX IN ARPACK BAND FORM')
         IF      (EIG_LAP_MAT_TYPE(1:3) == 'DPB') THEN
            CALL BANDGEN_LAPACK_DPB ( 'KMSM', NDOFL, KMSM_SDIA, NTERM_KMSM, I_KMSM, J_KMSM, KMSM, RFAC, SUBR_NAME )
         ELSE IF (EIG_LAP_MAT_TYPE(1:3) == 'DGB') THEN
            CALL BANDGEN_LAPACK_DGB ( 'KMSM', NDOFL, KMSM_SDIA, NTERM_KMSM, I_KMSM, J_KMSM, KMSM, RFAC, SUBR_NAME )
         ENDIF

! Write RFAC, if requested

         IF (DEBUG(40) == 1) THEN
            IF      (EIG_LAP_MAT_TYPE(1:3) == 'DPB') THEN
               CALL WRITE_MATRIX_BY_ROWS ( 'MATRIX [KLL - sigma*KLLD] IN LAPACK BAND FORM', RFAC, LDRFAC, NDOFL, F06 )
            ELSE
               CALL WRITE_MATRIX_BY_ROWS ( 'MATRIX [KLL - sigma*MLL ] IN LAPACK BAND FORM', RFAC, LDRFAC, NDOFL, F06 )
            ENDIF
         ENDIF

      ENDIF

! If this is not a CB or BUCKLING soln, dellocate arrays for KLL.      ! Keep arrays MLL, KLLD. Need them later to calc gen mass

      IF ((SOL_NAME(1:12) /= 'GEN CB MODEL' ) .AND. (SOL_NAME(1:8) /= 'BUCKLING')) THEN
         CALL LINK_MESSAGE('DEALLOCATE SPARSE KLL ARRAYS')
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KLL', CR13
         CALL DEALLOCATE_SPARSE_MAT ( 'KLL' )
      ENDIF

! Generate nonsymmetric storage form for KMSM. This is done since subr MATMULT_SFF, used when subr DSBAND is called herein,
! will run faster. MATMULT_SFF is called in each "Reverse commumication loop" in DSBAND.

      IF      (SPARSTOR == 'SYM   ') THEN

         CALL SPARSE_MAT_DIAG_ZEROS ( 'KMSM', NDOFL, NTERM_KMSM, I_KMSM, J_KMSM, NUM_KMSM_DIAG_0 )
         NTERM_KMSMn = 2*NTERM_KMSM  - (NDOFL - NUM_KMSM_DIAG_0)
         CALL LINK_MESSAGE('ALLOCATE SPARSE KMSMn ARRAYS')
         CALL ALLOCATE_SPARSE_MAT ( 'KMSMn', NDOFL, NTERM_KMSMn, SUBR_NAME )

         CALL LINK_MESSAGE('CONVERT SYM CRS KMSM TO NONSYM CRS KMSMn')
         CALL CRS_SYM_TO_CRS_NONSYM ( 'KMSM', NDOFL, NTERM_KMSM, I_KMSM, J_KMSM, KMSM, 'KMSMn', NTERM_KMSMn, I_KMSMn, J_KMSMn,     &
                                       KMSMn, 'Y' )

      ELSE IF (SPARSTOR == 'NONSYM') THEN

         NTERM_KMSMn = NTERM_KMSM
         CALL ALLOCATE_SPARSE_MAT ( 'KMSMn', NDOFL, NTERM_KMSMn, SUBR_NAME )
         DO I=1,NDOFL+1
            I_KMSMn(I) = I_KMSM(I)
         ENDDO
         DO J=1,NTERM_KMSMn
            J_KMSMn(J) = J_KMSM(J)
              KMSMn(J) =   KMSM(J)
         ENDDO

      ELSE                                                 !      Error - incorrect SPARSTOR

         WRITE(ERR,932) SUBR_NAME, SPARSTOR
         WRITE(F06,932) SUBR_NAME, SPARSTOR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF

      IF (DEBUG(49) > 0) THEN
! --- BANDED_optimizisation -begin-- !
         CALL WRITE_SPARSE_CRS ( ' KMSMn nonsymmetric operator used by ARPACK', 'A ', 'A ', NTERM_KMSMn, NDOFL, I_KMSMn, J_KMSMn,&
                                 KMSMn )
! --- BANDED_optimizisation -end-- !
      ENDIF

! Now we can deallocate KMSM (since KMSMn will be used in subr DSBAND)

!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KMSM', CR13
      CALL DEALLOCATE_SPARSE_MAT ( 'KMSM' )

! Check that user did not ask for more than NDOFL-1 eigens (Lanczos can't be used to find all). If request is > NDOFL-1, decrease
! request to NDOFL-1 and give warning
      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         LNONZEROS = NDOFL - NUM_KLLD_DIAG_ZEROS
      ELSE
         LNONZEROS = NDOFL - NUM_MLL_DIAG_ZEROS
      ENDIF

! --- arpack_surgery begin --- !
      IF (DEBUG(185) == 0) THEN
         NDOFL_EFFECTIVE = LNONZEROS
      ELSE
         NDOFL_EFFECTIVE = NDOFL
      ENDIF
! --- arpack_surgery end --- !


      NEV_EXTRA = DARPACK

      IF (NUM_EST_EIGENS > 0) THEN
         NEV = NUM_EST_EIGENS + EIG_LANCZOS_NEV_DELT
      ELSE
         IF (EIG_N2 <= 2) THEN
            NEV_EXTRA = MAX(NEV_EXTRA, 4)
         ENDIF
         IF ((EIG_N2 <= 2) .AND. (EIG_LANCZOS_NEV_DELT > NEV_EXTRA)) THEN
            NEV_EXTRA = EIG_LANCZOS_NEV_DELT
         ENDIF
         ! prevent the addition of DARPACK from crashing ARPACK
! --- arpack_surgery begin --- !
         IF ((EIG_N2 + NEV_EXTRA) > (NDOFL_EFFECTIVE - 4)) THEN
! --- arpack_surgery end --- !
            WRITE(ERR,9775) NDOFL, DARPACK
            IF (SUPWARN == 'N') THEN
               WRITE(F06,9775) NDOFL, DARPACK
            ENDIF
            NEV_EXTRA = 0
         END IF
         NEV = EIG_N2 + NEV_EXTRA
      ENDIF
      IF (DEBUG(185) == 0) THEN                            ! If 0, only find finite eigens within the range requested
         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            NUM1 = NDOFL - NUM_KLLD_DIAG_ZEROS
         ELSE
            NUM1 = NDOFL - NUM_MLL_DIAG_ZEROS
         ENDIF
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,48006) NUM1
         IF (SUPWARN == 'N') THEN
            WRITE(F06,48006) NUM1
         ENDIF
      ELSE
         NUM1 = NDOFL
      ENDIF
      IF (NEV > (NUM1-1)) THEN
!xx      WARN_ERR = WARN_ERR + 1
!xx      WRITE(ERR,4901) EIG_N2, NDOFL-1-DARPACK, NDOFL-1-DARPACK
!xx      IF (SUPWARN == 'N') THEN
!xx         WRITE(F06,4901) EIG_N2, NDOFL-1-DARPACK, NDOFL-1-DARPACK
!xx      ENDIF
         NEV = MAX(1, NUM1 - 1)
      ENDIF

      MIN_NCV = NEV + 2
! --- arpack_surgery begin --- !
      MAX_NCV = NDOFL_EFFECTIVE - 2
! --- arpack_surgery end --- !

! --- validation_fix7 begin --- !
! Small mode-2 decks can have an infeasible effective NCV window when
! LNONZEROS excludes released/zero-mass rotational DOF's. In that narrow
! case, allow an NDOFL-based ceiling before declaring the request invalid.
      IF (MIN_NCV > MAX_NCV) THEN
         MAX_NCV = NDOFL - 2
      ENDIF
! --- validation_fix7 end --- !

      ! sanity check on feasible NCV range
      IF (MIN_NCV > MAX_NCV) THEN
         WRITE(ERR,9776) NDOFL, NEV
         WRITE(F06,9776) NDOFL, NEV
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      END IF

      ! NCV determination step (try user factor first)
      NCV = MIN_NCV
      DO I=EIG_NCVFACL,2,-1
         NCV = I*NEV
         IF (NCV >= MIN_NCV .AND. NCV <= MAX_NCV) THEN
            EXIT
         END IF
      END DO

      ! still invalid? (too large OR too small)
      IF (NCV > MAX_NCV .OR. NCV < MIN_NCV) THEN
         DO I=5,2,-1
            NCV = NEV+I
            IF (NCV >= MIN_NCV .AND. NCV <= MAX_NCV) THEN
               EXIT
            END IF
         END DO
      END IF

! Small requests such as NEV=2 can converge poorly with the minimum Krylov
! subspace (e.g. NCV=4). Give ARPACK a slightly wider window when available.
      IF (NEV <= 2) THEN
         NCV = MAX(NCV, MIN(MAX_NCV, NEV + 4))
      ENDIF

      ! no valid NCV.
      IF (NCV > MAX_NCV .OR. NCV < MIN_NCV) THEN
         WRITE(ERR,9777) NDOFL, NEV
         WRITE(F06,9777) NDOFL, NEV
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      END IF

      ! inform the user if we didn't use EIG_NCFACL*NEV
      IF (NCV < (EIG_NCVFACL*NEV)) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,9778) NEV, NDOFL, NCV, (EIG_NCVFACL*NEV)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,9778) NEV, NDOFL, NCV, (EIG_NCVFACL*NEV)
         ENDIF
      END IF


      LWORKL = NCV*(NCV + 8)                               ! Size of WORKL workspace used by ARPACK

      KL = KMSM_SDIA
      KU = KL

! Set IPARAM, used in the ARPACK code

      DO I=1,11
         IPARAM(I) = 0                                     ! IPARAM(1) = 1 is set in subr DSBAND
      ENDDO
      IPARAM(3) = MXITERL                                  ! Max number of Arnoldi update iterations allowed
      IPARAM(4) = 1                                        ! NB, the blocksize used in the recurrence. ARPACK only works for NB = 1
      IPARAM(7) = EIG_MODE                                 ! If EIG_MODE = 2 solve (1/lambda)*K*x = M*x, M symmetric, K sym pos def
!                                                            If EIG_MODE = 3 solve (K - EIG_SIGMA*M)x = lbar*Mx (shift invert mode)
!                                                            where, for EIG_MODE = 3 the eigenvalues are lambda = EIG_SIGMA + lbar
!                                                            NOTE: EIG_MODE was checked in subr BD_EIGRL for correct value
! Set other inputs to DSBAND (see above for definitions)

      RVEC   = .TRUE.
      HOWMNY = 'A'
      IF (EIG_MODE == 1) THEN
         BMAT   = 'I'
      ELSE
         BMAT   = 'G'
      ENDIF

      ! DSEUPD may use a Ritz-vector guard column while extracting a small
      ! number of modes. Keep the full NCV workspace so tiny beam models do
      ! not overwrite an NDOFL-by-NEV heap allocation before LINK9.
      CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC', NDOFL, NCV, SUBR_NAME )
      CALL ALLOCATE_EIGEN1_MAT ( 'MODE_NUM' , NDOFL, 1, SUBR_NAME )
      CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL', NDOFL, 1, SUBR_NAME )

      CALL ALLOCATE_LAPACK_MAT ( 'IWORK' , NDOFL  , 1, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'RESID' , NDOFL  , 1, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'SELECT', NCV    , 1, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'VBAS'  , NDOFL  , NCV, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'WORKD' , 3*NDOFL, 1, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'WORKL' , LWORKL , 1, SUBR_NAME )

      ! Use a deterministic nonzero starting vector so repeated Lanczos
      ! solves over multiple METHOD requests do not depend on ARPACK's
      ! internal random start-vector choice.
      DO I=1,NDOFL
         RESID(I) = ONE + REAL(I,DOUBLE)/REAL(MAX(1,NDOFL),DOUBLE)
      ENDDO
      INFO_ARPACK = 1

      DO I=1,NCV                                           ! With HOWMNY = 'A' we are calc'ing eigenvecs for all eigenvalues found
         SELECT(I) = .FALSE.                               ! so all members of SELECT are .FALSE.
      ENDDO

      CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - LANCZOS METH')
!xx   WRITE(SC1, * )                                       ! Make new line for DTBSV pass messages that overwrite each other

      IF (DEBUG(50) == 1) CALL DEBUG_EIG_LANCZOS
      WHICH = 'LM'
      CALL DSBAND ( RVEC, HOWMNY, SELECT, EIGEN_VAL, EIGEN_VEC, NDOFL, EIG_SIGMA, NDOFL, LDRFAC, RFAC, KL, KU, WHICH, BMAT,        &
                    NEV, ARP_TOL, RESID, NCV, VBAS, NDOFL, IPARAM, WORKD, WORKL, LWORKL, IWORK, INFO_ARPACK,                       &
                    INFO_LAPACK, 'Y', DEBUG(47) )

      NVEC_SAFE  = MIN( IPARAM(5), NEV, SIZE(EIGEN_VEC,2), SIZE(EIGEN_VAL), SIZE(MODE_NUM) )
      NVEC       = NVEC_SAFE                               ! With HOWMNY = 'A' we are calc'ing eigenvecs for all eigenvalues found
      NUM_EIGENS = NVEC_SAFE

!xx   WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KMSMn ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KMSMn' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate IWORK ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'IWORK' )
! --- validation_fix5 begin --- !
! RFAC ownership is in LINK4 cleanup; avoid duplicate deallocation here.
!      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate RFAC  ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'RFAC' )
! --- validation_fix5 end --- !
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate RESID ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'RESID' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate SELECT', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'SELECT' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate VBAS  ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'VBAS' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate WORKD ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'WORKD' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate WORKL ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'WORKL' )
! --- MUMPS_COO add begin --- !
      IF ((SOLLIB(1:6) == 'SPARSE') .AND. (SPARSE_FLAVOR(1:5) == 'MUMPS')) THEN
         CALL DMUMPS_FREE_FACTORS()
      ENDIF
! --- MUMPS_COO add end --- !

      WRITE(ERR,4002) IPARAM(3)                            ! DSBAND returns the actual number of Arnoldi update iterations taken
      IF (SUPINFO == 'N') THEN
         WRITE(F06,4002) IPARAM(3)
      ENDIF

      IF ((EIG_MODE == 1) .OR. (EIG_MODE == 2)) THEN       ! Modes came from subr DSBAND in reverse order since 1, 2 means we solved
         DO I=1,NUM_EIGENS                                 ! for 1/lambda and when we invert the eigens will go from high to low
            MODE_NUM(I) = NUM_EIGENS + 1 - I
         ENDDO
      ELSE IF (EIG_MODE == 3) THEN                         ! Modes came from subr DSBAND in correct order since we solved for lambda
         DO I=1,NUM_EIGENS                                 ! not 1/lambda
            MODE_NUM(I) = I
         ENDDO
      ENDIF

      IF (INFO_ARPACK < 0) THEN
         WRITE(ERR,9996)
         WRITE(F06,9996)
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALLED_SUBR = 'DPBTRF'
      IF      (INFO_LAPACK < 0) THEN                       ! LAPACK subr XERBLA should have reported error on an illegal argument
!                                                            in calling a LAPACK subr, so we should not have gotten here
         WRITE(ERR,993) SUBR_NAME, CALLED_SUBR
         WRITE(F06,993) SUBR_NAME, CALLED_SUBR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ELSE IF (INFO_LAPACK > 0) THEN                       ! NOTE: this will only occur from call to LAPACK dgbtrs in ARPACK

         CALL GET_GRID_AND_COMP ( 'A ', INFO_LAPACK, GRIDV, COMPV  )

         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            WRITE(ERR,989) '[KLL + SIGMA*KLLD]', CALLED_SUBR, INFO_LAPACK
            WRITE(F06,989) '[KLL + SIGMA*KLLD]', CALLED_SUBR, INFO_LAPACK
         ELSE
            WRITE(ERR,989) '[KLL - SIGMA*MLL]' , CALLED_SUBR, INFO_LAPACK
            WRITE(F06,989) '[KLL - SIGMA*MLL]' , CALLED_SUBR, INFO_LAPACK
         ENDIF
         FATAL_ERR = FATAL_ERR + 1
         IF ((GRIDV > 0) .AND. (COMPV > 0)) THEN
            WRITE(ERR,9892) GRIDV, COMPV
            WRITE(F06,9892) GRIDV, COMPV
         ENDIF
         IF (BAILOUT >= 0) THEN                            ! If BAILOUT >= 0 then quit (leading minor <= 0)
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

      ENDIF

      IF ((EIG_MODE == 1) .OR. (EIG_MODE == 2)) THEN       ! We found inverses of eigens, so invert them
         CALL INVERT_EIGENS ( NUM_EIGENS, NDOFL, EIGEN_VAL, EIGEN_VEC, MODE_NUM )
      ENDIF

! With HOWMNY = 'A' we are calculating eigenvecs for all eigenvalues found

      IF (NUM_EST_EIGENS > 0) THEN
         NVEC = NUM_EIGENS
      ELSE
         ! Keep at most the user-requested modes, but do not discard
         ! converged modes just because the over-requested DARPACK margin
         ! was larger than the actual number returned by ARPACK.
         NVEC = MIN(EIG_N2, NUM_EIGENS)
      ENDIF
      NUM_EIGENS = NVEC




      RETURN

! **********************************************************************************************************************************
  101 FORMAT(6X,'Estimated number of eigens below ',F8.0,' Hz is ',I7)

  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

  989 FORMAT(' *ERROR   989: THE FACTORIZATION OF THE MATRIX ',A,' COULD NOT BE COMPLETED BY LAPACK SUBR ',A                       &
                    ,/,14X,' THE LEADING MINOR OF ORDER ',I12,' IS NOT POSITIVE DEFINITE')

  993 FORMAT(' *ERROR   993: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' LAPACK SUBR XERBLA SHOULD HAVE REPORTED AN ERROR ON AN ILLEGAL ARGUMENT IN A CALL TO LAPACK SUBR '    &
                    ,/,15X,A,' (OR A SUBR CALLED BY IT) AND THEN ABORTED',/)

 4002 FORMAT(' *INFORMATION: THE ACTUAL NUMBER OF ARNOLDI UPDATE ITERATIONS TAKEN IN THE SOLUTION FOR THE EIGENVALUES WAS = ',I8)

 4003 FORMAT(' *ERROR  4003: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' EIG_LAP_MAT_TYPE MUST BE EITHER "DGB" OR "DPB" BUT IS = ',A)

 4901 FORMAT(' *WARNING    : REQUEST FOR ',I8,' EIGENVALUES CANNOT BE HONORED. LANCZOS CAN BE USED TO FIND NO MORE THAN',          &
                           ' NEV = '                                                                                               &
                 ,I8,/,14X,' ATTEMPT WILL BE MADE TO FIND ',I8,' EIGENVALUES (SEE MANUAL FOR DEFINITION OF PARAMETER DARPACK)')

 4905 FORMAT(' *INFORMATION: NUMBER OF SUPERDIAGONALS IN THE KMSM = [KLL - sigma*MLL] MATRIX UPPER TRIANGLE IS      = ',I12,/)

 9892 FORMAT('               THIS IS FOR ROW AND COL IN THE MATRIX FOR GRID POINT ',I8,' COMPONENT ',I3)

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

 9996 FORMAT('  PROCESSING STOPPED DUE TO ABOVE ARPACK ERRORS')

12345 FORMAT(A,10X,A)

48006 FORMAT(' *WARNING    : THE L-SET MASS MATRIX HAS ONLY ',I8,' NONZEROS ON ITS DIAGONAL. THERE ARE NO MORE FINITE EIGENVALUES',&
                           ' BEYOND THIS NUMBER'                                                                                   &
                    ,/,14x,' (NOTE: USE OF BULK DATA PARAM ART_MASS WITH SMALL VALUE MAY ALLOW MGIV TO FIND MORE EIGENVALUES)')

 9775 FORMAT(' *WARNING    : DUE TO THE SIZE OF THE PROBLEM (NDOFL=',I0,'), PARAMETER DARPACK',/,15X,&
             'WAS SET TO ZERO TO PREVENT ARPACK FROM CRASHING. ORIGINAL VALUE: ',I0)

 9776 FORMAT(' *ERROR  9776: TOO MANY EIGENVALUES REQUESTED FOR THIS PROBLEM SIZE: NDOFL=',I0,', NEV=',I0,'.')

 9777 FORMAT(' *ERROR  9777: UNABLE TO DETERMINE KRYLOV SUBSPACE SIZE NCV. NDOFL=',I0,', NEV=',I0,'.',/,15X,&
             'USER HAS LIKELY REQUESTED TOO MANY EIGENVALUES FOR THIS MODEL.',/,15X,&
             'NCV MUST BE AT LEAST NEV+2, BUT NO MORE THAN NDOFL-2.')

 9778 FORMAT(' *WARNING    : DUE TO THE REQUESTED EIGENVALUES (NEV=',I0,') RELATIVE TO MODEL SIZE (NDOFL=',I0,'),',/,15X,&
             'NCVFACL WAS NOT USED FOR DETERMINATION OF THE KRYLOV SUBSPACE SIZE NCV.',/,15X,&
             'NCV WAS SET TO ',I0,' INSTEAD OF NEV*NCVFACL=',I0,'.')


! **********************************************************************************************************************************

      CONTAINS

! ##################################################################################################################################

! --- arpack_surgery begin --- !
      SUBROUTINE EST_NUM_EIGENS_BANDED ( FREQ, NUM_NEG_TERMS, COUNT_VALID )
! --- arpack_surgery end --- !

      USE LAPACK_GIV_MGIV_EIG
      USE LAPACK_LIN_EQN_DGE
      USE LAPACK_BLAS_AUX

      IMPLICIT NONE

! --- arpack_surgery begin --- !
      LOGICAL                         :: COUNT_VALID       !
      INTEGER(LONG)                   :: INFO_DPTTRF       !
      INTEGER(LONG)                   :: INFO_DSBTRD       !
! --- arpack_surgery end --- !
      INTEGER(LONG)                   :: NUM_NEG_TERMS     ! Number of negative terms on the diagonal of RFAC

      REAL(DOUBLE)                    :: DIAG(NDOFL)       !
      REAL(DOUBLE)                    :: FREQ              !
      REAL(DOUBLE)                    :: OFF_DIAG(NDOFL-1) !
      REAL(DOUBLE)                    :: QMAT(1,NDOFL)     !
      REAL(DOUBLE)                    :: SIGMA             !
      REAL(DOUBLE)                    :: WORK(NDOFL)       !

! **********************************************************************************************************************************
      NUM_NEG_TERMS = 0
! --- arpack_surgery begin --- !
      COUNT_VALID = .FALSE.
! --- arpack_surgery end --- !
      SIGMA = (TWO*PI*FREQ)**2

! Add KLL minus SIGMA*MLL used in Lanczos

      CALL MATADD_SSS_NTERM ( NDOFL, 'KLL', NTERM_KLL,  I_KLL,  J_KLL, SYM_KLL, '-SIGMA*MLL', NTERM_MLL, I_MLL, J_MLL, SYM_MLL,&
                                    'KMSM', NTERM_KMSM )
      CALL ALLOCATE_SPARSE_MAT ( 'KMSM', NDOFL, NTERM_KMSM, SUBR_NAME )
      CALL MATADD_SSS      ( NDOFL, 'KLL',  NTERM_KLL,  I_KLL,  J_KLL,  KLL,  ONE, '-SIGMA*MLL', NTERM_MLL, I_MLL, J_MLL, MLL, &
                            -SIGMA, 'KMSM', NTERM_KMSM, I_KMSM, J_KMSM, KMSM )

! Det bandwidth of KMSM so BANDGEN can put it in LAPACK band form. KMSM_SDIA is the number of super-diags in the band form of KMSM

      CALL BANDSIZ ( NDOFL, NTERM_KMSM, I_KMSM, J_KMSM, KMSM_SDIA )

      IF      (EIG_LAP_MAT_TYPE(1:3) == 'DPB') THEN
         LDRFAC = KMSM_SDIA + 1
      ELSE IF (EIG_LAP_MAT_TYPE(1:3) == 'DGB') THEN
         LDRFAC = 3*KMSM_SDIA + 1
      ENDIF

! Allocate array RFAC = (KLL - SIGMA*MLL) for ARACK

! --- BANDED_optimizisation -begin-- !
      CALL REPORT_SOLVER_DISPATCH_POLICY ( 'KMSM', SUBR_NAME )
! --- BANDED_optimizisation -end-- !
      CALL ALLOCATE_LAPACK_MAT ( 'RFAC', LDRFAC, NDOFL, SUBR_NAME )

! Put KMSM in form required by LAPACK band matrix. Call result array RFAC

      IF      (EIG_LAP_MAT_TYPE(1:3) == 'DPB') THEN
         CALL BANDGEN_LAPACK_DPB ( 'KMSM', NDOFL, KMSM_SDIA, NTERM_KMSM, I_KMSM, J_KMSM, KMSM, RFAC, SUBR_NAME )
      ELSE IF (EIG_LAP_MAT_TYPE(1:3) == 'DGB') THEN
         CALL BANDGEN_LAPACK_DGB ( 'KMSM', NDOFL, KMSM_SDIA, NTERM_KMSM, I_KMSM, J_KMSM, KMSM, RFAC, SUBR_NAME )
      ENDIF

! Reduce KMSM to tridiag form and do Cholesky L*D*L' decomp on it

! --- arpack_surgery begin --- !
      INFO_DSBTRD = 0
      INFO_DPTTRF = 0

      CALL DSBTRD ( 'N', 'U', NDOFL, KMSM_SDIA, RFAC, KMSM_SDIA+1, DIAG, OFF_DIAG, QMAT, 1, WORK, INFO_DSBTRD )
      IF (INFO_DSBTRD /= 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,'(A,1ES14.6,A,I0)') ' WARNING: EST_NUM_EIGENS_BANDED skipped at FREQ=', FREQ,                  &
                                       ' due to DSBTRD INFO=', INFO_DSBTRD
         IF (SUPWARN == 'N') THEN
            WRITE(F06,'(A,1ES14.6,A,I0)') ' WARNING: EST_NUM_EIGENS_BANDED skipped at FREQ=', FREQ,               &
                                           ' due to DSBTRD INFO=', INFO_DSBTRD
         ENDIF
         GO TO 900
      ENDIF

      CALL DPTTRF_MYSTRAN ( NDOFL, DIAG, OFF_DIAG, INFO_DPTTRF )
      IF (INFO_DPTTRF /= 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,'(A,1ES14.6,A,I0)') ' WARNING: EST_NUM_EIGENS_BANDED skipped at FREQ=', FREQ,                  &
                                       ' due to DPTTRF_MYSTRAN INFO=', INFO_DPTTRF
         IF (SUPWARN == 'N') THEN
            WRITE(F06,'(A,1ES14.6,A,I0)') ' WARNING: EST_NUM_EIGENS_BANDED skipped at FREQ=', FREQ,               &
                                           ' due to DPTTRF_MYSTRAN INFO=', INFO_DPTTRF
         ENDIF
         GO TO 900
      ENDIF
! --- arpack_surgery end --- !
      NUM_NEG_TERMS = 0
      DO I=1,NDOFL
         IF (DIAG(I) < ZERO) THEN
            NUM_NEG_TERMS = NUM_NEG_TERMS + 1
         ENDIF
      ENDDO
! --- arpack_surgery begin --- !
      COUNT_VALID = .TRUE.

  900 CONTINUE
! --- arpack_surgery end --- !
      CALL DEALLOCATE_SPARSE_MAT ( 'KMSM' )
      CALL DEALLOCATE_LAPACK_MAT ( 'RFAC' )

      RETURN

! **********************************************************************************************************************************
12345 FORMAT(A,10X)

! **********************************************************************************************************************************

      END SUBROUTINE EST_NUM_EIGENS_BANDED

! ##################################################################################################################################

      SUBROUTINE DEBUG_EIG_LANCZOS

      WRITE(F06,'(A,I8)')    ' IN EIG_LANCZOS_ARPACK: NUM_EST_EIGENS       = ', NUM_EST_EIGENS
! --- arpack_surgery begin --- !
      WRITE(F06,'(A,L8)')    '                      : EST_NUM_EIGENS_VALID = ', EST_NUM_EIGENS_VALID
! --- arpack_surgery end --- !
      WRITE(F06,'(A,I8)')    '                      : EIG_LANCZOS_NEV_DELT = ', EIG_LANCZOS_NEV_DELT
      WRITE(F06,'(A,I8)')    '                      : EIG_N2               = ', EIG_N2
      WRITE(F06,'(A,I8)')    '                      : DARPACK              = ', DARPACK
      WRITE(F06,'(A,I8)')    '                      : NDOFL                = ', NDOFL
      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         WRITE(F06,'(A,I8)') '                      : NUM_KLLD_DIAG_ZEROS  = ', NUM_KLLD_DIAG_ZEROS
      ELSE
         WRITE(F06,'(A,I8)') '                      : NUM_MLL_DIAG_ZEROS   = ', NUM_MLL_DIAG_ZEROS
      ENDIF
      WRITE(F06,'(A,I8)')    '                      : NUM1                 = ', NUM1
      WRITE(F06,'(A,I8)')    '                      : EIG_NCVFACL          = ', EIG_NCVFACL
      WRITE(F06,'(A,I8)')    '                      : NCV                  = ', NCV
      WRITE(F06,'(A,I8)')    '                      : NEV                  = ', NEV
      WRITE(F06,*)

! **********************************************************************************************************************************

      END SUBROUTINE DEBUG_EIG_LANCZOS

      END SUBROUTINE EIG_LANCZOS_ARPACK
