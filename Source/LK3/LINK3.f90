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
 
      SUBROUTINE LINK3
      #ifdef MKLDSS
      use mkl_dss  
      #endif
! LINK 3 solves the equation KLL*UL = PL where KLL, UL, PL are the L-set stiffness matrix, displs and loads. It solves the equation
! using one of three methods. For each method the solution is obtained in a 2 step process: (1) the KLL matrix is decomposed into
! triangular factors and (2) UL is solved for by forward-backward substitution (FBS). The 3 methods are:

!   a) The LAPACK freeware code. This code requires KLL to be in banded (NOT sparse) form. LAPACH has the advantage that
!      MYSTRAN contains the LAPACK source code so debugging is easy. Its disadvantage is that banded matrices require much more 
!      memory than sparse storage for large stiffness matrices.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_LOG, ERR, F04, F06, L3A, SC1, LINK3A, L3A_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, COMM, FATAL_ERR, KLL_SDIA, LINKNO, MBUG, NDOFL, NSUB,                       &
                                         NTERM_KLL, NTERM_PL, RESTART,  SOL_NAME, WARN_ERR, NTERM_KLLs
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC       
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, TEN
      USE PARAMS, ONLY                :  CRS_CCS, EPSERR, EPSIL, KLLRAT, RELINK3, RCONDK, SOLLIB, SUPWARN, SPARSE_FLAVOR, SPARSTOR
      USE SPARSE_MATRICES, ONLY       :  I_KLL, J_KLL, KLL, I_PL, J_PL, PL, I_KLLs, I2_KLLs, J_KLLs, KLLs
      USE LAPACK_DPB_MATRICES, ONLY   :  RES
      USE COL_VECS, ONLY              :  UL_COL, PL_COL
      USE MACHINE_PARAMS, ONLY        :  MACH_EPS, MACH_SFMIN
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE LAPACK_BLAS_AUX
      USE LAPACK_LIN_EQN_DPB
      USE SCRATCH_MATRICES, ONLY      :  I_CCS1, J_CCS1, CCS1
      USE SuperLU_STUF, ONLY          :  SLU_FACTORS, SLU_INFO

! Interface module not needed for subr's DPBTRF and DPBTRS. These are "CONTAIN'ed" in module LAPACK_LIN_EQN_DPB,
! which is "USE'd" above

!     USE LINK3_USE_IFs
                      
      IMPLICIT NONE
      
      #ifdef MKLDSS
      include 'mkl_pardiso.fi'
      #endif MKLDSS

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LINK3'
      CHARACTER(  2*BYTE)             :: L_SET    = 'L '   ! L-set designator
      CHARACTER(  1*BYTE)             :: EQUED             ! 'Y' if the stiff matrix was equilibrated in subr EQUILIBRATE    
      CHARACTER(  1*BYTE)             :: NULL_COL          ! 'Y' if a col of KAO(transpose) is null 
      CHARACTER( 54*BYTE)             :: MODNAM            ! Name to write to screen
 
      INTEGER(LONG)                   :: DEB_PRT(2)        ! Debug numbers to say whether to write ABAND and/or its decomp to output
!                                                            file in called subr SYM_MAT_DECOMP_LAPACK (ABAND = band form of KLL)

      INTEGER(LONG)                   :: IER_DECOMP        ! Overall error indicator
      INTEGER(LONG)                   :: ISUB              ! DO loop index for subcases 
      INTEGER(LONG)                   :: INFO     = 0      ! Info output from some routine that has been called
      INTEGER(LONG)                   :: I,J,memerror      ! DO loop indices            
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN  
      INTEGER(LONG), PARAMETER        :: P_LINKNO = 2      ! Prior LINK no's that should have run before this LINK can execute

      REAL(DOUBLE)                    :: BETA              ! Multiple for rhs for use in subr FBS
      REAL(DOUBLE)                    :: DEN               ! K_INORM*UL_INORM + PL_INORM
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero

      REAL(DOUBLE), allocatable       :: EQUIL_SCALE_FACS(:)!(NDOFL)
                                                           ! LAPACK_S values returned from subr SYM_MAT_DECOMP_LAPACK
      REAL(DOUBLE), allocatable       :: DUM_COL(:)!(NDOFL)    ! Temp variable used in SuperLU
      REAL(DOUBLE)                    :: K_INORM           ! Inf norm of KLL matrix (det in  subr COND_NUM)
      REAL(DOUBLE)                    :: LAP_ERR1          ! Bound on displ error = 2*OMEGAI/RCOND
      REAL(DOUBLE)                    :: OMEGAI            ! RES_INORM/DEN (similar to EPSILON)
      REAL(DOUBLE)                    :: OMEGAI0           ! Upper bound on OMEGAI. OMEGAI0 = 10*NDOFL*MACH_EPS
      REAL(DOUBLE)                    :: PL_INORM          ! Inf norm of load vector
      REAL(DOUBLE)                    :: RES_INORM         ! Inf norm of residual vector R = K*UL - PL 
      REAL(DOUBLE)                    :: RCOND             ! Recrip of cond no. of the KLL. Det in  subr COND_NUM
      REAL(DOUBLE)                    :: UL_INORM          ! Inf norm of displacement vector
 
      INTRINSIC                       :: DABS

      #ifdef MKLDSS 
      !DSS REAL
      INTEGER(LONG)                   :: NUM_KLL_DIAG_ZEROS  ! Number of zeros on the diag of KLL
      TYPE(MKL_DSS_HANDLE)            :: handle ! Allocate storage for the solver handle.      !DSS var
      INTEGER                         :: perm(1) ! DSS VAR
      integer                         :: KLLSused
      
      INTEGER                         :: dsserror
      REAL(DOUBLE),allocatable        :: SOLN(:)       ! Solution
      ! pardiso var
      INTEGER                         :: pardisoerror
      TYPE(MKL_PARDISO_HANDLE)           pt(64)
      !.. All other variables
      INTEGER                         :: maxfct, mnum, mtype, phase, msglvl
      INTEGER                         :: iparm(64)
      INTEGER                         :: idum(1)
      REAL*8                          :: ddum(1)
      #endif MKLDSS


!***********************************************************************************************************************************
      LINKNO = 3

      EPS1 = EPSIL(1)

      allocate(EQUIL_SCALE_FACS(NDOFL),DUM_COL(NDOFL),stat=memerror)
      if (memerror.ne.0) stop 'error in allocating EQUIL_SCALE_FACS in link3'

      #ifdef MKLDSS
      allocate( SOLN(NDOFL),stat=memerror)       ! Solution  
      if (memerror.ne.0) stop 'ERROR in allocating solution in link3'
      #endif MKLDSS

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
      IF (WRT_LOG > 0) THEN
         WRITE(F04,150) LINKNO
      ENDIF
      WRITE(ERR,150) LINKNO

! Read LINK1A file
 
      CALL READ_L1A ( 'KEEP', 'Y' )

! Check COMM for successful completion of prior LINKs

      IF (COMM(P_LINKNO) /= 'C') THEN
         WRITE(ERR,9998) P_LINKNO,P_LINKNO,LINKNO
         WRITE(F06,9998) P_LINKNO,P_LINKNO,LINKNO
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Prior LINK's didn't complete, so quit
      ENDIF

! Make sure SOL is STATICS, BUCKLING or NLSTATIC

      IF ((SOL_NAME(1:7) /= 'STATICS') .AND. (SOL_NAME(1:8) /= 'BUCKLING') .AND. (SOL_NAME(1:8) /= 'NLSTATIC')) THEN
         WRITE(ERR,999) SOL_NAME, 'STATICS or BUCKLING or NLSTATIC'
         WRITE(F06,999) SOL_NAME, 'STATICS or BUCKLING or NLSTATIC'
         CALL OUTA_HERE ( 'Y' )
      ENDIF

!***********************************************************************************************************************************
! Factor KLL

      DEB_PRT(1) = 34
      DEB_PRT(2) = 35
      IER_DECOMP = 0

      DO J=1,NDOFL                                         ! Need a null col of loads when SuperLU is called to factor KLL
         DUM_COL(J) = ZERO                                 ! (only because it appears in the calling list)
      ENDDO

      IF ((RESTART == 'Y') .AND. (RELINK3 == 'Y')) THEN     
sol_do:  DO
            WRITE(SC1,*) ' Input the value of SOLLIB (8 characters) to use in this restart:'
            READ (*,*) SOLLIB
            IF ((SOLLIB /= 'BANDED  ') .AND. (SOLLIB /= 'SPARSE  ')) THEN
               WRITE(SC1,*) '  Incorrect SOLLIB. Value must be BANDED or SPARSE'
               WRITE(SC1,*)
               CYCLE sol_do
            ELSE
               EXIT sol_do
            ENDIF
         ENDDO sol_do
      ENDIF

Factr:IF (SOLLIB == 'BANDED  ') THEN                       ! Use LAPACK

         INFO = 0
         CALL SYM_MAT_DECOMP_LAPACK ( SUBR_NAME, 'KLL', L_SET, NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, 'Y', KLLRAT, 'Y', RCONDK,      &
                                      DEB_PRT, EQUED, KLL_SDIA, K_INORM, RCOND, EQUIL_SCALE_FACS, INFO )

      ELSE IF (SOLLIB == 'SPARSE  ') THEN

         IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

            SLU_INFO = 0
            CALL SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME, 'KLL', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, SLU_INFO )
      #ifdef MKLDSS

         ELSEIF  (SPARSE_FLAVOR(1:3) == 'DSS') THEN  !DSS STARTED
            
            
            DO I=1,NDOFL                                         ! Need a null col of loads when SuperLU is called to factor KLL
                DUM_COL(I) = ZERO                                ! (only because it appears in the calling list)
            ENDDO 
            IF (CRS_CCS == 'CRS') THEN 
                
                
            WRITE(*,*) "Intel MKL Direct Sparse Solver Factoring"
                       
            IF      (SPARSTOR == 'SYM   ') THEN
            KLLSused = 0
            ! Initialize the solver.
                dsserror = DSS_CREATE(handle, MKL_DSS_DEFAULTS)
                
                IF (dsserror /= MKL_DSS_SUCCESS)  stop 'DSS error in initializing :' 
                SLU_INFO = dsserror
                dsserror =  DSS_DEFINE_STRUCTURE  (handle, MKL_DSS_SYMMETRIC, I_KLL,  NDOFL, NDOFL, J_KLL , NTERM_KLL)  
                SLU_INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in non zero defining:'
                perm(1) = 0
                dsserror = DSS_REORDER(handle, MKL_DSS_DEFAULTS, perm)
                SLU_INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in reordering :' 
                dsserror = DSS_FACTOR_REAL(handle, MKL_DSS_DEFAULTS, KLL)
                SLU_INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in Factoring:' 
            
            ELSE    
             !ALT1 convert to KLLs symmetry
             !KLLSused = 1
             !CALL SPARSE_MAT_DIAG_ZEROS ( 'KLL', NDOFL, NTERM_KLL, I_KLL, J_KLL, NUM_KLL_DIAG_ZEROS ) !get NUM_KLL_DIAG_ZEROS 
             !NTERM_KLLs = (NTERM_KLL  + (NDOFL - NUM_KLL_DIAG_ZEROS))/2
             !CALL ALLOCATE_SPARSE_MAT ( 'KLLs', NDOFL, NTERM_KLLs, SUBR_NAME )    
             !CALL CRS_NONSYM_TO_CRS_SYM ( 'KLL', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, 'KLLs', NTERM_KLLs, I_KLLs, J_KLLs, KLLs )
             
             !ALT2 just use nonsymmetry KLL
              KLLSused = 0
             
            ! Initialize the solver.
                dsserror = DSS_CREATE(handle, MKL_DSS_DEFAULTS)
                
                IF (dsserror /= MKL_DSS_SUCCESS)  stop 'DSS error in initializing :' 
                SLU_INFO = dsserror
            
            ! Define the non-zero structure of the matrix.
               !alt1
               ! dsserror =  DSS_DEFINE_STRUCTURE  (handle, MKL_DSS_SYMMETRIC, I_KLLs  , NDOFL, NDOFL, J_KLLs , NTERM_KLLs) !using KLLS
                                ! DSS_DEFINE_STRUCTURE ( HANDLE, MKL_DSS__SYMMETRIC, I_MAT    , NROWS, NROWS, J_MAT  , NTERMS ) 
                
               !alt2
                dsserror =  DSS_DEFINE_STRUCTURE  (handle, MKL_DSS_NON_SYMMETRIC, I_KLL  , NDOFL, NDOFL, J_KLL , NTERM_KLL) !using KLL
                !         error = DSS_DEFINE_STRUCTURE(handle, MKL_DSS_NON_SYMMETRIC, rowIndex=I_KLL,  nRows=NDOFL, nCols=NDOFL, columns=J_KLL , nNonZeros=NTERM_KLL)
                        
                SLU_INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in non zero defining:'
                perm(1) = 0
                
                
                !reorder
                dsserror = DSS_REORDER(handle, MKL_DSS_DEFAULTS, perm)
                SLU_INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in reordering :' 
                
                
                ! Factor the matrix. 
                !alt1
                !dsserror = DSS_FACTOR_REAL(handle, MKL_DSS_DEFAULTS, KLLs) !using KLLs
                
                !alt2
                dsserror = DSS_FACTOR_REAL(handle, MKL_DSS_DEFAULTS, KLL) !using KLL
                SLU_INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in Factoring:' 

            ENDIF

            ELSE IF (CRS_CCS == 'CCS') THEN    
                STOP 'CCS NOT YET'
            ENDIF    !DSS END
            
         ELSEIF  (SPARSE_FLAVOR(1:7) == 'PARDISO') THEN
             WRITE(*,*) "Intel MKL Pardiso"     
             
             DO i = 1, 64
                 iparm(i) = 0
             END DO!pardiso STARTED
             !iparm[64] This array is used to pass various parameters 
             !to Intel  oneAPI Math Kernel Library
             !PARDISO and to return some useful information after execution of the 
             !solver (see pardiso iparm Parameter for more details) 
             iparm(1) = 1 ! no solver default
             !If iparm[0] =0 Intel  oneAPI Math Kernel Library PARDISO fillsiparm [1]
             ! through iparm [63] with default values and uses them. 
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
             !n = bigsize !Number of equations in the sparse linear system A*X= B n>0         
             maxfct = 1 ! Maximal number of factors in memory >0  Generally used value is 1 
             mnum = 1 !The number of matrix (from 1 to maxfct) to solve; 
             ! Generally used value is 1                    
             IF      (SPARSTOR == 'SYM   ') THEN 
                 mtype = 1 ! real symmetric
             ELSE
                 mtype = 11 ! real unsymmetric 
             ENDIF
             
             !Matrix type 1 Real and structurally symmetric
             !            2 Real and symmetric positive definite
             !           -2 Real and symmetric indefinite
             !            3 Complex and structurally symmetric
             !            4 Complex and Hermitian positive definite
             !           -4 Complex and Hermitian indefinite
             !            6 Complex and symmetric matrix
             !           11 Real and nonsymmetric matrix
             !           13 Complex and nonsymmetric matrix 
        
        
             ! pt Solver internal data address pointer 0 
        
             ! Must be initialized with zeros and never be modified later 
        
             !.. Initialize the internal solver memory pointer. 
        
             !This is only necessary for the FIRST call of the PARDISO solver.
        
             DO i = 1, 64
            
                 pt(i)%DUMMY = 0
        
             END DO
        
        
             perm(1) = 0             !perm[n]
	    
             !Holds the permutation vector of size n , specifies elements used for 
             !computing a partial solution, or specifies differing values of the 
             !input matrices for low rank update >=0 
             !rhs and solution matrix allocated already
             ! rhs or B
             ! b[n*nrhs]
             ! Right-hand side vectors 
             ! On entry, contains the right-hand side vector/matrix B , which is placed contiguously in memory. The b[i+k*n]
             ! element must hold the i-th component of k-th right-hand side vector. Note that b is only accessed in the solution phase.
             ! On output, the array is replaced with the solution if iparm [5] =1.     
        
             !solution or x
        
             ! x [n*nrhs] Solution vectors 
             ! On output, if iparm [5] =0, contains solution vector/matrix X which is placed contiguously in memory. The 
             ! x[i+k*n] element must hold the i-th component of k-th solution vector. Note that x is only accessed in the solution phase.
             !.. Reordering and Symbolic Factorization, This step also allocates all memory that is necessary for the factorization
             
        
             phase = 11 ! only reordering and symbolic factorization
        
                
             !error
        
             !	0 No error
             !   -1 Input inconsistent
             ! -2 Not enough memory
             ! -3 Reordering problem        
             ! -4 Zero pivot, numerical factorization or iterative refinement problem
             ! -5 Unclassified (internal) error
             ! -6 Reordering failed (matrix types 11 and 13 only)
             ! -7 Diagonal matrix is singular
             ! -8 32-bit integer overflow problem
             ! -9 Not enough memory for OOC
             ! -10 Problems with opening OOC temporary files
             ! -11 Read/write problems with the OOC data file     
        
        
             !    pardiso (pt, maxfct, mnum, mtype, phase,n,a     ,ia       ,ja     & 
             !,perm, nrhs, iparm, msglvl, b   , x   , error)
        
         
             !error = DSS_DEFINE_STRUCTURE(handle, MKL_DSS_NON_SYMMETRIC, rowIndex=I_KLL,  nRows=NDOFL, nCols=NDOFL, columns=J_KLL , nNonZeros=NTERM_KLL)
       
             CALL pardiso (pt, maxfct, mnum, mtype, phase,NDOFL,KLL, I_KLL,J_KLL , idum, 1, iparm, msglvl, ddum, ddum, pardisoerror)
        
        
             ! WRITE(*,*) 'Reordering completed ... '
        
             IF (pardisoerror .NE. 0) THEN
                 WRITE(*,*) 'The following Pardiso ERROR was detected: ', pardisoerror
                 WRITE(err,*) 'The following Pardiso ERROR was detected: ', pardisoerror
                 STOP 1
             END IF
        
             !WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
             !WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)

        
             !.. Factorization.
        
             phase = 22 ! only factorization
          
             !error = DSS_FACTOR_REAL(handle, MKL_DSS_DEFAULTS, values)
        
             CALL pardiso (pt, maxfct, mnum, mtype, phase, NDOFL, KLL,I_KLL,J_KLL,idum, 1, iparm, msglvl, ddum, ddum, pardisoerror)
        
             !WRITE(*,*) 'Factorization completed ... '
        
             IF (pardisoerror .NE. 0) THEN
            
                 WRITE(*,*) 'The following Pardiso  ERROR was detected: ', pardisoerror
                 WRITE(err,*) 'The following Pardiso ERROR was detected: ', pardisoerror
                 STOP 1
        
             END IF
      #endif MKLDSS


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

      ENDIF Factr

!***********************************************************************************************************************************
!  Allocate col vector arrays for loads, displs and res vector

!xx   CALL ALLOCATE_COL_VEC ( 'UL_COL', NDOFL, SUBR_NAME )
!xx   CALL ALLOCATE_COL_VEC ( 'PL_COL', NDOFL, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'RES', NDOFL, 1, SUBR_NAME )

! Open file for writing displs to.
 
      CALL FILE_OPEN ( L3A, LINK3A, OUNT, 'REPLACE', L3A_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N', 'Y' )
 
! Loop on subcases

      WRITE(F06,*)
      BETA = ONE
Solve:DO ISUB = 1,NSUB

         SLU_INFO = 0
         CALL ALLOCATE_COL_VEC ( 'UL_COL', NDOFL, SUBR_NAME )
         CALL ALLOCATE_COL_VEC ( 'PL_COL', NDOFL, SUBR_NAME )

         CALL OURTIM                                       ! Get the loads for this subcase from I_PL, J_PL, PL and put into PL_COL
         MODNAM = 'GET COL OF PL LOADS FOR                        Subcase'
         WRITE(SC1,3093) LINKNO,MODNAM,ISUB,HOUR,MINUTE,SEC,SFRAC
         DO J=1,NDOFL
            PL_COL(J)  = ZERO
            DUM_COL(J) = ZERO
         ENDDO
         CALL GET_SPARSE_CRS_COL ( 'PL        ', ISUB, NTERM_PL, NDOFL, NSUB, I_PL, J_PL, PL, BETA, PL_COL, NULL_COL )
         DO J=1,NDOFL
            DUM_COL(J) = PL_COL(J)
         ENDDO
 
         IF (DEBUG(32) == 1) THEN                          ! DEBUG output of load vector for this subcase, if requested
            WRITE(F06,3020) ISUB
            CALL WRITE_VECTOR ( '      L-SET LOADS      ',' LOAD', NDOFL, PL_COL )
            WRITE(F06,*)
         ENDIF
 
         CALL OURTIM                                       ! Call FBS to solve for displacements for this subcase
         MODNAM = 'FBS - SOLVE FOR RHS ANSWERS FOR                   "'
         WRITE(SC1,3093) LINKNO,MODNAM,ISUB,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages

         IF      (SOLLIB == 'BANDED  ') THEN

            CALL FBS_LAPACK ( EQUED, NDOFL, KLL_SDIA, EQUIL_SCALE_FACS, DUM_COL )

         ELSE IF (SOLLIB == 'SPARSE  ') THEN

            IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

               SLU_INFO = 0
               CALL FBS_SUPRLU ( SUBR_NAME, 'KLL', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, ISUB, DUM_COL, SLU_INFO )

      #ifdef MKLDSS
            ELSEIF  (SPARSE_FLAVOR(1:3) == 'DSS') THEN  !DSS STARTED
                
               SLU_INFO = 0       
            
               iF      (CRS_CCS == 'CRS') THEN                      ! Use KLL stored in Compressed Row Storage (CRS) format
                    
                    dsserror = DSS_SOLVE_REAL(handle, MKL_DSS_DEFAULTS ,DUM_COL ,1, SOLN)
                    SLU_INFO = dsserror
                    IF (dsserror /= MKL_DSS_SUCCESS) then 
                        stop 'DSS error in Solving :'
                    else
                      DO I=1,NDOFL
                         DUM_COL(I) = SOLN(I)
                      ENDDO 
                    write (F06,9902)  'KLL','LINK3'
9902                FORMAT(' DSS FACTORIZATION OF MATRIX ', A, ' SUCCEEDED IN SUBR ', A)
                    endif    

                ELSE IF (CRS_CCS == 'CCS') THEN    
                    STOP 'CCS NOT YET'
                ENDIF    !DSS END
                
            ELSEIF  (SPARSE_FLAVOR(1:7) == 'PARDISO') THEN  !pardiso STARTED
                
                !.. Back substitution and iterative refinement
                iparm(8) = 2 ! max numbers of iterative refinement steps
                phase = 33 ! only factorization
                
                soln = 0.d0
                CALL pardiso (pt, maxfct, mnum, mtype, phase, ndofl, KLL, I_KLL, J_KLL, idum, 1, iparm, msglvl,DUM_COL, SOLN, pardisoerror)    
                
                    
                IF (pardisoerror /= 0) then 
                        
                    stop 'Pardiso error in Solving : '
                    
                else
                      
                    DO I=1,NDOFL
                         
                        DUM_COL(I) = SOLN(I)
                      
                    ENDDO
                  
                    write (F06,9903)  'KLL','LINK3'
9903                FORMAT(' PARDISO FACTORIZATION OF MATRIX ', A, ' SUCCEEDED IN SUBR ', A)
                    
                endif
                
                
      #endif MKLDSS

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

         DO J=1,NDOFL
            UL_COL(J) = DUM_COL(J)
         ENDDO

         IF (DEBUG(33) == 1) THEN                          ! DEBUG output of displs
            WRITE(F06,3022) ISUB
            CALL WRITE_VECTOR ( '      A-SET DISPL      ','DISPL', NDOFL, UL_COL )
            WRITE(F06,*)
         ENDIF
 
         IF (EPSERR == 'Y') THEN                           ! Calculate residual vector, R. Use RES to calculate EPSILON
            CALL OURTIM
            MODNAM = 'CALC  EPSILON ERROR ESTIMATE                      "'
            WRITE(SC1,3093) LINKNO,MODNAM,ISUB,HOUR,MINUTE,SEC,SFRAC
            CALL EPSCALC ( ISUB )
         ENDIF
                                                           ! Calculate the LAPACK error bounds
         IF ((RCONDK == 'Y') .AND. (SOLLIB == 'BANDED')) THEN 
            IF (DABS(RCOND) > MACH_SFMIN) THEN
               CALL OURTIM
               MODNAM = 'CALC LAPACK ERROR ESTIMATE                        "'
               WRITE(SC1,3093) LINKNO,MODNAM,ISUB,HOUR,MINUTE,SEC,SFRAC
               CALL VECINORM ( UL_COL, NDOFL,  UL_INORM )
               CALL VECINORM ( PL_COL, NDOFL,  PL_INORM )
               CALL VECINORM ( RES   , NDOFL, RES_INORM )
               DEN = K_INORM*UL_INORM + PL_INORM 
               IF (DABS(DEN) > EPS1) THEN
                  OMEGAI = (RES_INORM)/(DEN)
                  OMEGAI0 = TEN*NDOFL*MACH_EPS
                  LAP_ERR1 = TWO*OMEGAI/RCOND
                  WRITE(F06,3024) ISUB, LAP_ERR1, OMEGAI, RCOND, DEN, RES_INORM, K_INORM, UL_INORM, PL_INORM, OMEGAI0, MACH_EPS
               ELSE
                  WRITE(F06,3026)
               ENDIF
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,3025) ISUB, RCOND, MACH_SFMIN
               IF (SUPWARN == 'N') THEN 
                  WRITE(F06,3025) ISUB, RCOND, MACH_SFMIN
               ENDIF 
            ENDIF
         ENDIF

         DO J=1,NDOFL                                      ! Write UL to file L3A for this subcase
            WRITE(L3A) UL_COL(J)
         ENDDO

         CALL DEALLOCATE_COL_VEC  ( 'UL_COL' )
         CALL DEALLOCATE_COL_VEC  ( 'PL_COL' )
         IF (kllsused.eq.1) CALL DEALLOCATE_SPARSE_MAT ( 'KLLs' )

      ENDDO Solve

FreeS:IF (SOLLIB == 'SPARSE  ') THEN                       ! Last, free the storage allocated inside SuperLU

         IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

            DO J=1,NDOFL                                         ! Need a null col of loads when SuperLU is called to factor KLL
               DUM_COL(J) = ZERO                                  ! (only because it appears in the calling list)
            ENDDO

            CALL C_FORTRAN_DGSSV( 3, NDOFL, NTERM_KLL, 1, KLL , I_KLL , J_KLL , DUM_COL, NDOFL, SLU_FACTORS, SLU_INFO )

            IF (SLU_INFO .EQ. 0) THEN
               WRITE (*,*) 'SUPERLU STORAGE FREED'
            ELSE
               WRITE(*,*) 'SUPERLU STORAGE NOT FREED. INFO FROM SUPERLU FREE STORAGE ROUTINE = ', SLU_INFO
            ENDIF
      #ifdef MKLDSS
         ELSEIF  (SPARSE_FLAVOR(1:3) == 'DSS') THEN  !DSS STARTED

            DO J=1,NDOFL                                         ! Need a null col of loads when SuperLU is called to factor KLL
               DUM_COL(J) = ZERO                                  ! (only because it appears in the calling list)
            ENDDO

                ! Deallocate solver storage and various local arrays.
                dsserror = DSS_DELETE(handle, MKL_DSS_DEFAULTS)
                deallocate(SOLN)  
                IF (dsserror /= MKL_DSS_SUCCESS) STOP 'DSS error in CLEARING :' 
                
         ELSEIF  (SPARSE_FLAVOR(1:7) == 'PARDISO') THEN  !DSS STARTED
            DO J=1,NDOFL                                         ! Need a null col of loads when SuperLU is called to factor KLL
               DUM_COL(J) = ZERO                                  ! (only because it appears in the calling list)
            ENDDO
             !.. Termination and release of memory
            phase = -1 ! release internal memory
            CALL pardiso (pt, maxfct, mnum, mtype, phase, NDOFL, ddum, idum,  idum, idum, 1, iparm, msglvl, ddum, ddum, pardisoerror)
            deallocate(  SOLN)       ! Solution   
      #endif MKLDSS

         ENDIF

      ENDIF FreeS
 
! Dellocate arrays

      CALL OURTIM
      MODNAM = 'DEALLOCATE ARRAYS'
      WRITE(SC1,3092) LINKNO,MODNAM,HOUR,MINUTE,SEC,SFRAC
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages         

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CONTINUE
      ELSE
         IF (SOL_NAME(1:12) /= 'GEN CB MODEL' ) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KLL  ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'KLL' )
         ENDIF
      ENDIF

      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate ABAND ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'ABAND' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate RES   ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'RES' )
!xx   WRITE(SC1,12345,ADVANCE='NO') '       Deallocate UL_COL', CR13   ;   CALL DEALLOCATE_COL_VEC  ( 'UL_COL' )
!xx   WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PL_COL', CR13   ;   CALL DEALLOCATE_COL_VEC  ( 'PL_COL' )
!xx   WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PL    ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PL' )

      CALL FILE_CLOSE ( L3A, LINK3A, 'KEEP', 'Y' )

! Process is now complete so set COMM(LINKNO)
  
      COMM(LINKNO) = 'C'

! Write data to L1A

      CALL WRITE_L1A ( 'KEEP', 'Y', 'Y' )
  
! Check allocation status of allocatable arrays, if requested

      IF (DEBUG(100) > 0) THEN
         CALL CHK_ARRAY_ALLOC_STAT
         IF (DEBUG(100) > 1) THEN
            CALL WRITE_ALLOC_MEM_TABLE ( 'at the end of '//SUBR_NAME )
         ENDIF
      ENDIF

! Write LINK3 end to F04, F06

      CALL OURTIM
      IF (WRT_LOG > 0) THEN
         WRITE(F04,151) LINKNO
      ENDIF
      WRITE(F06,151) LINKNO

! Close files
  
      IF (( DEBUG(193) == 3) .OR. (DEBUG(193) == 999)) THEN
         CALL FILE_INQUIRE ( 'near end of LINK3' )
      ENDIF
      deallocate(EQUIL_SCALE_FACS,DUM_COL)
! Write LINK3 end to screen

      WRITE(SC1,153) LINKNO
!***********************************************************************************************************************************
  150 FORMAT(/,' >> LINK',I3,' BEGIN',/)

  151 FORMAT(/,' >> LINK',I3,' END',/)

  152 FORMAT(/,' >> LINK',I3,' BEGIN')

  153 FORMAT(  ' >> LINK',I3,' END')

  933 FORMAT(' *ERROR   933: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CRS_CCS  MUST BE EITHER "CRS" OR "CCS" BUT VALUE IS ',A)

  999 FORMAT(' *ERROR   999: INCORRECT SOLUTION IN EXEC CONTROL. SHOULD BE ',A,', BUT IS SOL = ',A)

 3020 FORMAT(//,18X,'LOAD VECTOR FOR SUBCASE ',I8)

 3022 FORMAT(//,18X,'DISPLACEMENTS FOR SUBCASE ',I8,/23X,'LSET DOF',10X,'DISP',14X,'S(J)')

 3024 FORMAT(' *INFORMATION: FOR INTERNAL SUBCASE NUMBER ',I8,' LAPACK ERROR EST (2*OMEGAI/RCOND) = ',1ES13.6,                     &
             ' Gen, slightly > than true err'                                                                                   ,/,&
                                          52X,'................................................................................',/,&
             '                                                    ... OMEGAI                        = ',1ES13.6,                   &
             ' (RES_INORM/DEN)              .'                                                                                  ,/,&
             '                                                    ... RCOND                         = ',1ES13.6,                   &
             ' (Recriprocal of KLL cond num).'                                                                                  ,/,&
             '                                                    ... DEN                           = ',1ES13.6,                   &
             ' (K_INORM*UL_INORM + PL_INORM).'                                                                                  ,/,&
             '                                                    ... RES_INORM                     = ',1ES13.6,                   &
             ' (Inf norm of KLL*UL - PL)    .'                                                                                  ,/,&
             '                                                    ... K_INORM                       = ',1ES13.6,                   &
             ' (Infinity norm of KLL)       .'                                                                                  ,/,&
             '                                                    ... UL_INORM                      = ',1ES13.6,                   &
             ' (Infinity norm of UL displs) .'                                                                                  ,/,&
             '                                                    ... PL_INORM                      = ',1ES13.6,                   &
             ' (Infinity norm of PL loads)  .'                                                                                  ,/,&
             '                                                    ... OMEGAI0 (OMEGAI upper bound)  = ',1ES13.6,                   &
             ' (10*NDOFL*MACH_EPS)          .'                                                                                  ,/,&
             '                                                    ... MACH_EPS                      = ',1ES13.6,                   &
             ' (Machine precision)          .'                                                                                  ,/,&
                                          52X,'................................................................................',/)

 3025 FORMAT(' *WARNING    : CANNOT CALCULATE LAPACK ERROR ESTIMATE FOR INTERNAL SUBCASE NUMBER ',I8                               &
                    ,/,14X,' THE RECIPROCAL OF THE CONDITION NUMBER OF KLL, RCOND         = ',1ES15.6,' CANNOT BE INVERTED.'       &
                    ,/,14X,' IT IS TOO SMALL COMPARED TO MACHINE SAFE MINIMUN (MACH_SFMIN) = ',1ES15.6,/)

 3026 FORMAT(' *INFORMATION: CANNOT CALCULATE OMEGAI. DEN = 0',/)

 3092 FORMAT(1X,I2,'/',A54,8X,2X,I2,':',I2,':',I2,'.',I3)

 3093 FORMAT(1X,I2,'/',A54,I8,2X,I2,':',I2,':',I2,'.',I3)

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

 9998 FORMAT(' *ERROR  9998: COMM ',I3,' INDICATES UNSUCCESSFUL LINK ',I2,' COMPLETION.'                                           &
                    ,/,14X,' FATAL ERROR - CANNOT START LINK ',I2)

12345 FORMAT(A,10X,A)

!***********************************************************************************************************************************

      END SUBROUTINE LINK3






