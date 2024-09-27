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

       SUBROUTINE EIG_INV_PWR
#ifdef MKLDSS
         use mkl_dss   
#endif MKLDSS  
! Solves for eigenvalues and eigenvectors when method is INV. Code is only valid for the 1st eigenval/vec. Inverse Power is an
! iterative method
 
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, WRT_LOG, ERR, F04, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, KMSM_SDIA, LINKNO, NDOFL, NTERM_KLL, NTERM_KLLD, NTERM_KMSM,     &
                                         NTERM_KMSMs, NTERM_MLL, NUM_EIGENS, NVEC, SOL_NAME, WARN_ERR
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  BAILOUT, EPSIL, KLLRAT, MXITERI, SOLLIB, SPARSE_FLAVOR, SPARSTOR, SUPINFO, SUPWARN !, CRS_CSS
      USE SUBR_BEGEND_LEVELS, ONLY    :  EIG_INV_PWR_BEGEND
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL, EIGEN_VEC, MODE_NUM
      USE MODEL_STUF, ONLY            :  EIG_N2, EIG_SIGMA
      USE SPARSE_MATRICES, ONLY       :  I_KLL, J_KLL, KLL, I_KLLD, J_KLLD, KLLD, I_MLL, J_MLL, MLL,                               &
                                         I_KMSM, I2_KMSM, J_KMSM, KMSM, I_KMSMs, I2_KMSMs, J_KMSMs, KMSMs 
      USE SPARSE_MATRICES, ONLY       :  SYM_KLL, SYM_KLLD, SYM_MLL
      USE LAPACK_LIN_EQN_DPB 
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE SuperLU_STUF, ONLY          :  SLU_FACTORS
 
      USE EIG_INV_PWR_USE_IFs

      IMPLICIT NONE
      
#ifdef MKLDSS
      include 'mkl_pardiso.fi'
#endif MKLDSS
  
      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EIG_INV_PWR'
      CHARACTER(  1*BYTE)             :: EQUED             ! 'Y' if KLL stiff matrix was equilibrated in subr EQUILIBRATE    
      CHARACTER(44*BYTE)              :: MODNAM            ! Name to write to screen to describe module being run.

      INTEGER(LONG)                   :: DEB_PRT(2)        ! Debug numbers to say whether to write ABAND and/or its decomp to output
!                                                            file in called subr SYM_MAT_DECOMP_LAPACK (ABAND = band form of KLL)

      INTEGER(LONG)                   :: I, memerror       ! DO loop index
      INTEGER(LONG)                   :: INFO        = 0   ! 
      INTEGER(LONG)                   :: ITER_NUM          ! Number of iterations in converging on eigenvalue 

      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = EIG_INV_PWR_BEGEND

      REAL(DOUBLE),allocatable        :: EIGEN_VAL_APPROX(:)!(0:MXITERI)
                                                           ! Eigenvalue at a given iteration number

      REAL(DOUBLE)                    :: K_INORM           ! Inf norm of KOO matrix
      REAL(DOUBLE),allocatable        :: MVEC(:,:)!(NDOFL,1)     ! MLL*EIGEN_VEC (or KLLD*EIGEN_VEC for BUCKLING)
      REAL(DOUBLE)                    :: MAX_VALUE         ! Max value from EIGEN_VEC(I,1)
      REAL(DOUBLE),allocatable        :: NULL_SCALE_FACS(:)!(NDOFL)
                                                           ! KMSM will not be equilibrated so set these to zero
      REAL(DOUBLE)                    :: PERCENT_CHANGE    ! % change in eigenvalue estimate between two successive iterations
      REAL(DOUBLE)                    :: RCOND             ! Recrip of cond no. of the KLL. Det in  subr COND_NUM

#ifdef MKLDSS
      !DSS REAL
      TYPE(MKL_DSS_HANDLE)            :: handle ! Allocate storage for the solver handle.      !DSS var
      INTEGER                         :: perm(1) ! DSS VAR   
      INTEGER                         :: dsserror
      REAL(DOUBLE),allocatable        :: SOLN(:),rhs(:)       ! Solution
            ! pardiso var
      INTEGER                         :: pardisoerror
      TYPE(MKL_PARDISO_HANDLE)           pt(64)
      !.. All other variables
      INTEGER                         :: maxfct, mnum, mtype, phase, msglvl
      INTEGER                         :: iparm(64)
      INTEGER                         :: idum(1)
      REAL*8                          :: ddum(1)
#endif

      INTRINSIC                       :: MIN

      allocate(EIGEN_VAL_APPROX(0:MXITERI),stat=memerror)
      allocate(MVEC(NDOFL,1),  NULL_SCALE_FACS(NDOFL), stat=memerror )
      if (memerror.ne.0) stop 'error allocating memory at eig_inv_pwr'

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************
! Check that user did not ask for more than 1 eigen (currently Inverse Power can't be used to find more than 1 eigen).
! If request is > NDOFL-1, decrease request to 1 and give warning

      IF (EIG_N2 > 1) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,4901) EIG_N2
         IF (SUPWARN == 'N') THEN
            WRITE(F06,4901) EIG_N2
         ENDIF
         EIG_N2 = 1
      ENDIF

      CALL OURTIM
      MODNAM = 'SOLVE FOR EIGENVALS/VECTORS - INV POWER METH'
      WRITE(SC1,4092) LINKNO,MODNAM,HOUR,MINUTE,SEC,SFRAC

! Calc KMSM = KLL - EIG_SIGMA*MLL (or - EIG_SIGMA*KLLD for BUCKLING) where EIG_SIGMA = shift freq

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CALL MATADD_SSS_NTERM ( NDOFL, 'KLL',  NTERM_KLL , I_KLL , J_KLL , SYM_KLL ,  'eig_sigma*KLLD',                           &
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

! Allocate arrays for eigen matrices

      CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC', NDOFL, EIG_N2, SUBR_NAME )
      CALL ALLOCATE_EIGEN1_MAT ( 'MODE_NUM' , NDOFL, 1, SUBR_NAME )
      CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL', NDOFL, 1, SUBR_NAME )

!***********************************************************************************************************************************
! Factor KMSM

      DEB_PRT(1) = 44
      DEB_PRT(2) = 45

      EQUED = 'N'
      IF (SOLLIB == 'BANDED  ') THEN

         INFO = 0
         CALL SYM_MAT_DECOMP_LAPACK ( SUBR_NAME, 'KMSM', 'L ', NDOFL, NTERM_KMSM, I_KMSM, J_KMSM, KMSM, 'Y', KLLRAT, 'N', 'N',     &
                                      DEB_PRT, EQUED, KMSM_SDIA, K_INORM, RCOND, NULL_SCALE_FACS, INFO )

      ELSE IF (SOLLIB == 'SPARSE  ') THEN

         IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

            INFO = 0
            CALL SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME, 'KMSM', NDOFL, NTERM_KMSM, I_KMSM, J_KMSM, KMSM, INFO )

#ifdef MKLDSS
         ELSEIF  (SPARSE_FLAVOR(1:3) == 'DSS') THEN  !DSS STARTED           
            
            ! IF (CRS_CCS == 'CCS') STOP 'CCS NOT YET'
                
            WRITE(*,*) "Intel MKL Direct Sparse Solver Factoring"
                       
            ! Initialize the solver.
                allocate(SOLN(ndofl),rhs(ndofl) )
                dsserror = DSS_CREATE(handle, MKL_DSS_DEFAULTS)
                
                IF (dsserror /= MKL_DSS_SUCCESS)  stop 'DSS error in initializing :' 
                iNFO = dsserror
            
            ! Define the non-zero structure of the matrix.
                IF (SPARSTOR == 'SYM   ') THEN
                    dsserror =  DSS_DEFINE_STRUCTURE  (handle, MKL_DSS_SYMMETRIC, I_KMSM  , NDOFL, NDOFL, J_KMSM, NTERM_KMSM) !using KMSM  
                ELSE
                    dsserror =  DSS_DEFINE_STRUCTURE  (handle, MKL_DSS_NON_SYMMETRIC, I_KMSM  , NDOFL, NDOFL, J_KMSM, NTERM_KMSM) !using KMSM    
                ENDIF    
                
                        
                INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in non zero defining:'
                perm(1) = 0
                
                !reorder
                dsserror = DSS_REORDER(handle, MKL_DSS_DEFAULTS, perm)
                INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in reordering :' 
                
                
                ! Factor the matrix. 
                dsserror = DSS_FACTOR_REAL(handle, MKL_DSS_DEFAULTS, KMSM) !using KMSM
                INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in Factoring:' 
                
          
         ELSEIF  (SPARSE_FLAVOR(1:7) == 'PARDISO') THEN
             WRITE(*,*) "Intel MKL Pardiso"
             !IF (CRS_CCS == 'CCS') STOP 'CCS NOT YET'                       
                      
                    DO i = 1, 64
                        iparm(i) = 0
                    END DO!pardiso STARTED !iparm[64] This array is used to pass various parameters 
              
             
                    iparm(1) = 1 ! no solver default
                    iparm(2) = 2 ! fill-in reordering from METIS
                    iparm(3) = 1 ! numbers of processors       
                    iparm(4) = 0 ! no iterative-direct algorithm
                    iparm(5) = 0 ! no user fill-in reducing permutation
                    iparm(6) = 0 ! solution on the first n components of x
                    iparm(7) = 0 ! not in use
                    iparm(8) = 9 ! numbers of iterative refinement steps
                    iparm(9) = 0 ! not in use
                    iparm(10) = 13 ! perturb the pivot elements with 1E-13
                    iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
                    iparm(12) = 0 ! not in use
                    iparm(13) = 1 ! maximum weighted matching algorithm is ON
                    iparm(14) = 0 ! Output: number of perturbed pivots
                    iparm(15) = 0 ! not in use
                    iparm(16) = 0 ! not in use
                    iparm(17) = 0 ! not in use
                    iparm(18) = -1 ! Output: number of nonzeros in the factor LU
                    iparm(19) = -1 ! Output: Mflops for LU factorization
                    iparm(20) = 0 ! Output: Numbers of CG Iterations
                    pardisoerror = 0 ! initialize error flag
             
                    msglvl = 0 ! print statistical information 1=on 0=off        
                    maxfct = 1 ! Maximal number of factors in memory >0  Generally used value is 1 
                    mnum = 1 !The number of matrix (from 1 to maxfct) to solve; 
             
                    IF      (SPARSTOR == 'SYM   ') THEN
                        mtype = 1 ! real symmetric
                    ELSE
                        mtype = 11 ! real unsymmetric 
                    ENDIF
             
        
                    DO i = 1, 64         
                        pt(i)%DUMMY = 0
                    END DO
             
                    perm(1) = 0             !perm[n]
             
                    phase = 11 ! only reordering and symbolic factorization         
                    CALL pardiso (pt, maxfct, mnum, mtype, phase,NDOFL ,KMSM, I_KMSM, J_KMSM , idum, 1, iparm, msglvl, ddum, ddum, pardisoerror)
        
             
                    IF (pardisoerror .NE. 0) THEN
                 
                        WRITE(*,*) 'The following Pardiso ERROR was detected: ', pardisoerror
                        WRITE(err,*) 'The following Pardiso ERROR was detected: ', pardisoerror
                        STOP 1
             
                    END IF
    
                    !.. Factorization.
        
             
                    phase = 22 ! only factorization
             
                    CALL pardiso (pt, maxfct, mnum, mtype, phase, NDOFL, KMSM,I_KMSM,J_KMSM,idum, 1, iparm, msglvl, ddum, ddum, pardisoerror)
                    IF (pardisoerror .NE. 0) THEN
                        WRITE(*,*) 'The following Pardiso  ERROR was detected: ', pardisoerror
                        WRITE(err,*) 'The following Pardiso ERROR was detected: ', pardisoerror
                        STOP 1
                    END IF       
       
                
                
                
#endif

         ELSE

            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,9991) SUBR_NAME, 'SPARSE_FLAVOR'
            WRITE(F06,9991) SUBR_NAME, 'SPARSE_FLAVOR'
            CALL OUTA_HERE ( 'Y' )

         ENDIF

      ELSE

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9991) SUBR_NAME, 'SOLLIB'
         WRITE(F06,9991) SUBR_NAME, 'SOLLIB'
         CALL OUTA_HERE ( 'Y' )

      ENDIF

      IF (EQUED == 'Y') THEN                               ! If EQUED == 'Y' then error. We don't want KMSM equilibrated from the
         WRITE(ERR,4001) SUBR_NAME, EQUED                  ! call (above) to SYM_MAT_DECOMP_LAPACK 
         WRITE(F06,4001) SUBR_NAME, EQUED
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      DO I=1,NDOFL
         NULL_SCALE_FACS(I) = ZERO
      ENDDO

!***********************************************************************************************************************************
! Get starting value for eigenvector

      DO I=1,NDOFL
         EIGEN_VEC(I,1) = ONE
      ENDDO

! Loop until convergence

      EIGEN_VAL_APPROX(0) = ONE
!xx   WRITE(SC1, * ) '    Iterating on eigenvector number   1:'
!xx   WRITE(SC1, * ) '       Iter No.  Approx Eigenvalue  % Change from last'
!xx   WRITE(SC1, * )

      IF (DEBUG(46) == 1) THEN
         WRITE(F06,4902) EIGEN_VAL_APPROX(0)
      ENDIF

      ITER_NUM = 0
iters:DO

         ITER_NUM = ITER_NUM + 1

! Mult MLL*EIGEN_VEC (or KLLD*EIGEN_VEC for BUCKLING) to get MVEC. This is the "RHS" in the solution
! [KLL - sigma*MLL]*Vec = alpha*MLL*Vec (EIGEN_VAL = sigma + 1/alpha)

         IF (SOL_NAME(1:8) == 'BUCKLING') THEN

            CALL MATMULT_SFF ('KLLD', NDOFL, NDOFL, NTERM_KLLD, SYM_KLLD, I_KLLD, J_KLLD, KLLD, 'EIGENVEC', NDOFL, 1,              &
                               EIGEN_VEC, 'N', 'MVEC',-ONE, MVEC)
         ELSE

            CALL MATMULT_SFF ('MLL' , NDOFL, NDOFL, NTERM_MLL , SYM_MLL , I_MLL , J_MLL , MLL , 'EIGENVEC', NDOFL, 1,              &
                               EIGEN_VEC, 'N', 'MVEC', ONE, MVEC)
         ENDIF

         IF      (SOLLIB == 'BANDED  ') THEN

            CALL FBS_LAPACK ( 'N', NDOFL, KMSM_SDIA, NULL_SCALE_FACS, MVEC )

         ELSE IF (SOLLIB == 'SPARSE  ') THEN

            IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

               INFO = 0

               IF (SOL_NAME(1:8) == 'BUCKLING') THEN
                  CALL FBS_SUPRLU ( SUBR_NAME, 'KLLD', NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, KLLD, ITER_NUM, MVEC, INFO )
               ELSE
                  CALL FBS_SUPRLU ( SUBR_NAME, 'KLL' , NDOFL, NTERM_KLL , I_KLL , J_KLL , KLL , ITER_NUM, MVEC, INFO )
               ENDIF
#ifdef MKLDSS
            ELSEIF  (SPARSE_FLAVOR(1:3) == 'DSS') THEN
                   
               IF (SOL_NAME(1:8) == 'BUCKLING') THEN 
                IF (SPARSTOR == 'SYM   ') THEN
                    dsserror =  DSS_DEFINE_STRUCTURE  (handle, MKL_DSS_SYMMETRIC, I_KLLD  , NDOFL, NDOFL,J_KLLD, NTERM_KLLD) !using KMSM  
                ELSE
                    dsserror =  DSS_DEFINE_STRUCTURE  (handle, MKL_DSS_NON_SYMMETRIC, I_KLLD  , NDOFL, NDOFL, J_KLLD, NTERM_KLLD) !using KMSM    
                ENDIF    
                INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in non zero defining:'
                perm(1) = 0
                !reorder
                dsserror = DSS_REORDER(handle, MKL_DSS_DEFAULTS, perm)
                INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in reordering :' 
            
                ! Factor the matrix. 
                dsserror = DSS_FACTOR_REAL(handle, MKL_DSS_DEFAULTS, KLLD) !using KMSM
                INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in Factoring:'
                
                
                    !IF (CRS_CCS == 'CCS') STOP 'CCS NOT YET'
                    DO I=1,NDOFL
                         RHS(i) = MVEC(I,1)
                    ENDDO 
                    
                    dsserror = DSS_SOLVE_REAL (handle, MKL_DSS_DEFAULTS ,rhs , 1, SOLN)
                    INFO = dsserror
                    IF (dsserror /= MKL_DSS_SUCCESS) then 
                        stop 'DSS error in Solving :'
                    else
                      DO I=1,NDOFL
                         MVEC(I,1) = SOLN(I)
                      ENDDO 
                    write (F06,9902)  'KLLD','EIG_INV_PWR'
9902                FORMAT(' DSS FACTORIZATION OF MATRIX ', A, ' SUCCEEDED IN SUBR ', A)
                    endif
               else
                IF (SPARSTOR == 'SYM   ') THEN
                    dsserror =  DSS_DEFINE_STRUCTURE  (handle, MKL_DSS_SYMMETRIC, I_KLL  , NDOFL, NDOFL,J_KLL, NTERM_KLL) !using KMSM  
                ELSE
                    dsserror =  DSS_DEFINE_STRUCTURE  (handle, MKL_DSS_NON_SYMMETRIC, I_KLL  , NDOFL, NDOFL, J_KLL, NTERM_KLL) !using KMSM    
                ENDIF    
                INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in non zero defining:'
                perm(1) = 0
                !reorder
                dsserror = DSS_REORDER(handle, MKL_DSS_DEFAULTS, perm)
                INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in reordering :' 
            
                ! Factor the matrix. 
                dsserror = DSS_FACTOR_REAL(handle, MKL_DSS_DEFAULTS, KLL)
                INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in Factoring:'
                
                
                    ! IF (CRS_CCS == 'CCS') STOP 'CCS NOT YET'
                    DO I=1,NDOFL
                         RHS(i) = MVEC(I,1)
                    ENDDO 
                    
                    dsserror = DSS_SOLVE_REAL (handle, MKL_DSS_DEFAULTS ,rhs , 1, SOLN)
                    INFO = dsserror
                    IF (dsserror /= MKL_DSS_SUCCESS) then 
                        stop 'DSS error in Solving :'
                    else
                      DO I=1,NDOFL
                         MVEC(I,1) = SOLN(I)
                      ENDDO 
                    endif                   
               endif
                    
            ELSEIF  (SPARSE_FLAVOR(1:7) == 'PARDISO') THEN  !pardiso STARTED
                
                IF (SOL_NAME(1:8) == 'BUCKLING') THEN
                    DO i = 1, 64         
                        pt(i)%DUMMY = 0
                    END DO
             
                    DO I=1,NDOFL
                         RHS(i) = MVEC(I,1)
                    ENDDO 
                    
                    perm(1) = 0             !perm[n]
             
                    phase = 11 ! only reordering and symbolic factorization         
                    CALL pardiso (pt, maxfct, mnum, mtype, phase,NDOFL,KLLD, I_KLLD, J_KLLD , idum, 1, iparm, msglvl, ddum, ddum, pardisoerror)
        
             
                    IF (pardisoerror .NE. 0) THEN
                 
                        WRITE(*,*) 'The following Pardiso ERROR was detected: ', pardisoerror
                        WRITE(err,*) 'The following Pardiso ERROR was detected: ', pardisoerror
                        STOP 1
             
                    END IF
    
                    !.. Factorization.
        
             
                    phase = 22 ! only factorization
             
                    CALL pardiso (pt, maxfct, mnum, mtype, phase, NDOFL, KLLD,I_KLLD,J_KLLD,idum, 1, iparm, msglvl, ddum, ddum, pardisoerror)
                    IF (pardisoerror .NE. 0) THEN
                        WRITE(*,*) 'The following Pardiso  ERROR was detected: ', pardisoerror
                        WRITE(err,*) 'The following Pardiso ERROR was detected: ', pardisoerror
                        STOP 1
                    END IF       
     
                    
                    
                    
                    !.. Back substitution and iterative refinement
                iparm(8) = 2 ! max numbers of iterative refinement steps
                phase = 33 ! only factorization
                
                soln = 0.d0
                CALL pardiso (pt, maxfct, mnum, mtype, phase, ndofl, KLLD, I_KLLD, J_KLLD, idum, 1, iparm, msglvl,rhs, SOLN, pardisoerror)        
                    
                IF (pardisoerror /= 0) then 
                        
                    stop 'Pardiso error in Solving : '
                    
                else
                      
                      DO I=1,NDOFL
                         MVEC(I,1) = SOLN(I)
                      ENDDO 
                  
                    write (F06,9903)  'KLL','LINK3'
9903                FORMAT(' PARDISO FACTORIZATION OF MATRIX ', A, ' SUCCEEDED IN SUBR ', A)
                    
                endif
                    
                else
!----
                  
                    DO i = 1, 64         
                        pt(i)%DUMMY = 0
                    END DO
                    DO I=1,NDOFL
                         RHS(i) = MVEC(I,1)
                    ENDDO 
                    
                    perm(1) = 0             !perm[n]
             
                    phase = 11 ! only reordering and symbolic factorization         
                    CALL pardiso (pt, maxfct, mnum, mtype, phase,NDOFL,KLL, I_KLL, J_KLL , idum, 1, iparm, msglvl, ddum, ddum, pardisoerror)
       
                    IF (pardisoerror .NE. 0) THEN
                 
                        WRITE(*,*) 'The following Pardiso ERROR was detected: ', pardisoerror
                        WRITE(err,*) 'The following Pardiso ERROR was detected: ', pardisoerror
                        STOP 1
             
                    END IF
    
                    !.. Factorization.
        
             
                    phase = 22 ! only factorization
             
                    CALL pardiso (pt, maxfct, mnum, mtype, phase, NDOFL, KLL,I_KLL,J_KLL,idum, 1, iparm, msglvl, ddum, ddum, pardisoerror)
                    IF (pardisoerror .NE. 0) THEN
                        WRITE(*,*) 'The following Pardiso  ERROR was detected: ', pardisoerror
                        WRITE(err,*) 'The following Pardiso ERROR was detected: ', pardisoerror
                        STOP 1
                    END IF       
     
                    
                    
                    
                    !.. Back substitution and iterative refinement
                iparm(8) = 2 ! max numbers of iterative refinement steps
                phase = 33 ! only factorization
                
                soln = 0.d0
                CALL pardiso (pt, maxfct, mnum, mtype, phase, ndofl, KLL, I_KLL, J_KLL, idum, 1, iparm, msglvl,rhs, SOLN, pardisoerror)        
                    
                IF (pardisoerror /= 0) then 
                        
                    stop 'Pardiso error in Solving : '
                    
                else
                      DO I=1,NDOFL
                         MVEC(I,1) = SOLN(I)
                      ENDDO 
                                          
                endif
!----  
                endif
                
#endif

            ELSE

               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,9991) SUBR_NAME, 'SPARSE_FLAVOR'
               WRITE(F06,9991) SUBR_NAME, 'SPARSE_FLAVOR'
               CALL OUTA_HERE ( 'Y' )

            ENDIF

         ELSE

            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,9991) SUBR_NAME, 'SOLLIB'
            WRITE(F06,9991) SUBR_NAME, 'SOLLIB'
            CALL OUTA_HERE ( 'Y' )

         ENDIF

         MAX_VALUE = ZERO                                  ! Determine max term in approximate MVEC and check for > 0
         DO I=1,NDOFL
            IF (DABS(MVEC(I,1)) > DABS(MAX_VALUE)) THEN
               MAX_VALUE = MVEC(I,1)
            ENDIF
         ENDDO

         IF (DABS(MAX_VALUE) > EPSIL(1)) THEN              ! If max value in eigenvector > 0 then get next eigenvalue approx

!zz         EIGEN_VAL_APPROX(ITER_NUM) = ONE/MAX_VALUE     ! 11/25/11: WRONG !!    
                                                           ! If MVEC had its max numerical value repeated with a different sign,
!                                                            the algorithm could converge to the negative of the actual eigenvalue.
!                                                            This algorithm is valid only for positive eigens, so use ABS(MAX_VALUE)
            EIGEN_VAL_APPROX(ITER_NUM) = ONE/DABS(MAX_VALUE)

         ELSE                                              ! If all values in eigenvector are zero, write messages and quit

            FATAL_ERR = FATAL_ERR + 1
            IF (ITER_NUM > 1) THEN
               IF (DEBUG(46) == 1) THEN
                  WRITE(F06,*)
               ELSE
                  CALL WRITE_ITER_RESULTS
               ENDIF
            ENDIF
            WRITE(ERR,4007) ITER_NUM, MAX_VALUE
            WRITE(F06,4007) ITER_NUM, MAX_VALUE
            CALL OUTA_HERE ( 'Y' )

         ENDIF
                                                           ! Calc % change in 2 successive estimates (if den /= 0)
         IF (DABS(EIGEN_VAL_APPROX(ITER_NUM-1)) > EPSIL(1)) THEN
            PERCENT_CHANGE = 1.0D2*(EIGEN_VAL_APPROX(ITER_NUM) - EIGEN_VAL_APPROX(ITER_NUM-1))/EIGEN_VAL_APPROX(ITER_NUM-1)
            WRITE(SC1,12345,ADVANCE='NO') ITER_NUM, EIG_SIGMA+EIGEN_VAL_APPROX(ITER_NUM), PERCENT_CHANGE, CR13
            IF (DEBUG(46) == 1) THEN
               WRITE(F06,4903) ITER_NUM, EIG_SIGMA+EIGEN_VAL_APPROX(ITER_NUM), PERCENT_CHANGE
            ENDIF
         ELSE
            WRITE(SC1,22345,ADVANCE='NO') ITER_NUM, EIG_SIGMA+EIGEN_VAL_APPROX(ITER_NUM), CR13
            IF (DEBUG(46) == 1) THEN
               WRITE(F06,4904) ITER_NUM, EIG_SIGMA+EIGEN_VAL_APPROX(ITER_NUM)
            ENDIF
         ENDIF

         DO I=1,NDOFL                                      ! Normalize eigenvector to max value
            EIGEN_VEC(I,1) = MVEC(I,1)/MAX_VALUE
         ENDDO

         IF (DABS(PERCENT_CHANGE) < EPSIL(3)) THEN         ! Test for convergence
            IF (DEBUG(46) == 1) THEN
               WRITE(F06,*)
            ENDIF
            EXIT iters
         ELSE
            IF (ITER_NUM < MXITERI) THEN
               CYCLE iters
            ELSE
               IF (ITER_NUM > 1) THEN
                  IF (DEBUG(46) == 1) THEN
                     WRITE(F06,*)
                  ELSE
                     CALL WRITE_ITER_RESULTS
                  ENDIF
               ENDIF
               WRITE(ERR,4006) MXITERI, EPSIL(3)
               WRITE(F06,4006) MXITERI, EPSIL(3)
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF
         ENDIF

      ENDDO iters

      NVEC = 1
      NUM_EIGENS = 1
      EIGEN_VAL(1) = EIG_SIGMA + EIGEN_VAL_APPROX(ITER_NUM)
      WRITE(ERR,4008) ITER_NUM, EPSIL(3)
      IF (SUPINFO == 'N') THEN
         WRITE(F06,4008) ITER_NUM, EPSIL(3)
      ENDIF

      DO I=1,NUM_EIGENS
         MODE_NUM(I) = I
      ENDDO

!added

       FreeS:IF (SOLLIB == 'SPARSE  ') THEN  

         IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN
             INFO = 0
             CALL C_FORTRAN_DGSSV( 3, NDOFL, NTERM_KMSM, 1, KMSM , I_KMSM , J_KMSM , MVEC, NDOFL, SLU_FACTORS, INFO )

            IF (INFO .EQ. 0) THEN
               WRITE (*,*) 'SUPERLU STORAGE FREED'
            ELSE
               WRITE(*,*) 'SUPERLU STORAGE NOT FREED. INFO FROM SUPERLU FREE STORAGE ROUTINE = ', INFO
            ENDIF

#ifdef MKLDSS
         ELSEIF  (SPARSE_FLAVOR(1:3) == 'DSS') THEN  !DSS STARTED

                ! Deallocate solver storage and various local arrays.
                dsserror = DSS_DELETE(handle, MKL_DSS_DEFAULTS)
                deallocate(SOLN,rhs)
                IF (dsserror /= MKL_DSS_SUCCESS) STOP 'DSS error in CLEARING :' 
         ELSEIF  (SPARSE_FLAVOR(1:7) == 'PARDISO') THEN  !DSS STARTED
             !.. Termination and release of memory
             phase = -1 ! release internal memory
             CALL pardiso (pt, maxfct, mnum, mtype, phase, NDOFL, ddum, idum,  idum, idum, 1, iparm, msglvl, ddum, ddum, pardisoerror)
             deallocate(  SOLN,rhs)       ! Solution 
         
#endif
        ENDIF
      ENDIF FreeS        
      


!added
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages         
      WRITE(SC1,32345,ADVANCE='NO') '       Deallocate KMSM'
      CALL DEALLOCATE_SPARSE_MAT ( 'KMSM' )

! If this is not a CB or BUCKLING soln, dellocate arrays for KLL.

      IF ((SOL_NAME(1:12) /= 'GEN CB MODEL' ) .AND. (SOL_NAME(1:8) /= 'BUCKLING')) THEN
         CALL OURTIM
         MODNAM = 'DEALLOCATE SPARSE KLL ARRAYS'
         WRITE(SC1,4092) LINKNO,MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages         
         WRITE(SC1,32345,ADVANCE='NO') '       Deallocate KLL', CR13
         CALL DEALLOCATE_SPARSE_MAT ( 'KLL' )
      ENDIF

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF
      deallocate(EIGEN_VAL_APPROX)
      deallocate(MVEC, NULL_SCALE_FACS )

      RETURN

! **********************************************************************************************************************************
  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

 9892 FORMAT('               THIS IS FOR ROW AND COL IN THE MATRIX FOR GRID POINT ',I8,' COMPONENT ',I3)

 4001 FORMAT(' *ERROR  4001: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' MATRIX KMSM WAS EQUILIBRATED: EQUED = ',A,'. CODE NOT WRITTEN TO ALLOW THIS AS YET')

 4006 FORMAT(' *ERROR  4006: MAXIMUM ITERATIONS = ',I8,' HAVE BEEN TAKEN WITHOUT CONVERGING ON AN EIGENVALUE.'                     &
                    ,/,14X,' BULK DATA PARAM MXITERI CAN BE USED TO INCREASE THE NUMBER OF ITERATIONS OR PARAM EPSIL 3 = ',1ES14.6 &
                    ,/,14X,' (THE ITERATION % CHANGE CRITERIA) CAN BE CHANGED')

 4007 FORMAT(' *ERROR  4007: CANNOT CONTINUE EIGENVALUE ITERATION. ALL OF THE TERMS IN THE EIGENVECTOR FOR ITERATION NUMBER ',I4,  &
                           ' ARE LESS THAN ',1ES15.6)

 4008 FORMAT(' *INFORMATION: THE INVERSE POWER METHOD PERFORMED ',I8,' ITERATIONS TO CONVERGE TO THE FIRST EIGENVALUE WITHIN '     &
                            ,1ES9.1,'%')

 4092 FORMAT(1X,I2,'/',A44,18X,2X,I2,':',I2,':',I2,'.',I3)

 4901 FORMAT(' *WARNING    : REQUEST FOR ',I8,' EIGENVALUES CANNOT BE HONORED. INVERSE POWER CAN BE USED TO FIND NO MORE THAN ONE' &
                 ,I8,/,14X,' ATTEMPT WILL BE MADE TO FIND ONE EIGENVALUE')

 4902 FORMAT(43X,'Results of Inverse Power iteration on eigenvalue    1',//,42X,'Iter No.       Approx Eigenvalue     ',           &
                 '% Change from last',//,52X,1ES23.14)
 
 4903 FORMAT(45X,I4,3X,1ES23.14,8X,1ES9.2)

 4904 FORMAT(45X,I4,3X,1ES23.14)

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

12345 FORMAT(10X,I4,3X,1ES15.6,2X,1ES15.2,A)

22345 FORMAT(10X,I4,3X,1ES15.6)

32345 FORMAT(A,10X)

! ##################################################################################################################################
 
      CONTAINS
 
! ##################################################################################################################################

      SUBROUTINE WRITE_ITER_RESULTS

      IMPLICIT NONE

      INTEGER(LONG)                   :: II                ! DO loop index

! **********************************************************************************************************************************
      WRITE(f06,4912) EIGEN_VAL_APPROX(0)

      DO II=1,ITER_NUM
         WRITE(F06,4913) II, EIG_SIGMA+EIGEN_VAL_APPROX(II)
      ENDDO
      WRITE(F06,*)

! **********************************************************************************************************************************
 4912 FORMAT(39X,'Results of Inverse Power iteration on eigenvalue    1',//,52X,'Iter No.  Approx Eigenvalue',//,62X,1ES15.6)
 
 4913 FORMAT(55X,I4,3X,1ES15.6)

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_ITER_RESULTS

      END SUBROUTINE EIG_INV_PWR
