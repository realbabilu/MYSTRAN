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
 
      SUBROUTINE SOLVE_DLR
      #ifdef MKLDSS
         use mkl_dss   
      #endif MKLDSS
! Solves KLL*DLR = -KLR for matrix DLR. However, we will use rows of KRL instead of cols of KLR in the solution.
 
! For a description of Craig-Bamptom analyses, see Appendix D to the MYSTRAN User's Referance Manual


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, WRT_ERR, WRT_LOG, ERR, F04, F06, SCR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FACTORED_MATRIX, FATAL_ERR, KLL_SDIA, NDOFR, NDOFL, NTERM_DLR, NTERM_KLL,   &
                                         NTERM_KRL
      USE PARAMS, ONLY                :  EPSIL, PRTDLR, SOLLIB, SPARSE_FLAVOR, SPARSTOR, NOCOUNTS, CRS_CCS
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, TSEC
      USE SUBR_BEGEND_LEVELS, ONLY    :  SOLVE_DLR_BEGEND
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE SPARSE_MATRICES, ONLY       :  I2_DLR, I_DLR, J_DLR, DLR, I_DLRt, I2_DLRt, J_DLRt, DLRt, I_KRL, J_KRL, KRL,              &
                                         I_KLL, I2_KLL, J_KLL, KLL
      USE SuperLU_STUF, ONLY          :  SLU_FACTORS
 
                                         
      USE LAPACK_LIN_EQN_DPB

      USE SOLVE_DLR_USE_IFs

      IMPLICIT NONE

      #ifdef MKLDSS
      include 'mkl_pardiso.fi'
      #endif MKLDSS
      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SOLVE_DLR'
      CHARACTER(  1*BYTE)             :: CLOSE_IT          ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not 
      CHARACTER(  8*BYTE)             :: CLOSE_STAT        ! What to do with file when it is closed
      CHARACTER(  1*BYTE)             :: EQUED             ! 'Y' if KLL stiff matrix was equilibrated in subr EQUILIBRATE    
      CHARACTER( 24*BYTE)             :: MESSAG            ! File description. Input to subr UNFORMATTED_OPEN 
      CHARACTER( 22*BYTE)             :: MODNAM1           ! Name to write to screen to describe module being run
      CHARACTER(  1*BYTE)             :: READ_NTERM        ! 'Y' or 'N' Input to subr READ_MATRIX_1 
      CHARACTER(  1*BYTE)             :: NULL_COL          ! 'Y' if a col of KLR(transpose) is null 
      CHARACTER(  1*BYTE)             :: OPND              ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to open  a file or not 
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: SCRFIL            ! File name
 
      INTEGER(LONG)                   :: DEB_PRT(2)        ! Debug numbers to say whether to write ABAND and/or its decomp to output
!                                                            file in called subr SYM_MAT_DECOMP_LAPACK
      INTEGER(LONG)                   :: I,J,memerror      ! DO loop indices or counters
      INTEGER(LONG)                   :: INFO        = 0   ! Info on success of factorization or solve
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN   
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = SOLVE_DLR_BEGEND

      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero

      REAL(DOUBLE),allocatable        :: EQUIL_SCALE_FACS(:)!(NDOFL)
                                                           ! LAPACK_S values returned from subr SYM_MAT_DECOMP_LAPACK

      REAL(DOUBLE),allocatable        :: DLR_COL(:)!(NDOFL)    ! A column of DLR solved for herein
      REAL(DOUBLE),allocatable        :: INOUT_COL(:)!(NDOFL)  ! Temp variable for subr FBS
      REAL(DOUBLE)                    :: K_INORM           ! Inf norm of KLL matrix (det in  subr COND_NUM)
      REAL(DOUBLE)                    :: RCOND             ! Recrip of cond no. of the KLL. Det in  subr COND_NUM
 
      #ifdef MKLDSS
      !DSS REAL
      TYPE(MKL_DSS_HANDLE)            :: handle ! Allocate storage for the solver handle.      !DSS var
      INTEGER                         :: perm(1) ! DSS VAR     
      !INTEGER buff(bufLen) ! DSS VAR
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
      

      INTRINSIC                       :: DABS
      
      allocate(EQUIL_SCALE_FACS(NDOFL) , DLR_COL(NDOFL) , INOUT_COL(NDOFL) , stat=memerror )
      if (memerror.ne.0) stop 'error allocating memory at solve_dlr'


! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Decomp KLL

      DEB_PRT(1) = 64
      DEB_PRT(2) = 65

      IF (SOLLIB == 'BANDED') THEN

         INFO  = 0
         EQUED = 'N'
         CALL SYM_MAT_DECOMP_LAPACK ( SUBR_NAME, 'KLL', 'L ', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, 'Y', 'N', 'N', 'N', DEB_PRT, &
                                      EQUED, KLL_SDIA, K_INORM, RCOND, EQUIL_SCALE_FACS, INFO )
         IF (EQUED == 'Y') THEN                         ! If EQUED == 'Y' then error. We don't want KLL equilibrated
            WRITE(ERR,6001) SUBR_NAME, EQUED
            WRITE(F06,6001) SUBR_NAME, EQUED
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

      ELSE IF (SOLLIB == 'SPARSE  ') THEN

         IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

            INFO = 0
            CALL SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME, 'KLL', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, INFO )

      #ifdef MKLDSS
         ELSEIF (SPARSE_FLAVOR(1:3) == 'DSS') THEN
                IF (CRS_CCS == 'CCS') STOP 'CCS NOT YET in DSS'
                allocate( SOLN(NDOFL))       ! Solution
                INFO = 0
                ! Initialize the solver.
                dsserror = DSS_CREATE(handle, MKL_DSS_DEFAULTS)
                
                
                IF (dsserror /= MKL_DSS_SUCCESS)  stop 'DSS error in initializing :' 
                INFO = dsserror
                IF      (SPARSTOR == 'SYM   ') THEN
                    dsserror =  DSS_DEFINE_STRUCTURE  (handle, MKL_DSS_SYMMETRIC, I_KLL  , NDOFL, NDOFL, J_KLL , NTERM_KLL) !using KLL 
                ELSE
                    dsserror =  DSS_DEFINE_STRUCTURE  (handle, MKL_DSS_NON_SYMMETRIC, I_KLL  , NDOFL, NDOFL, J_KLL , NTERM_KLL) !using KLL           
                ENDIF

                INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in non zero defining:'
                perm(1) = 0
             
                                
                !reorder
                dsserror = DSS_REORDER(handle, MKL_DSS_DEFAULTS, perm)
                INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in reordering :' 
                
                
                !factor
                dsserror = DSS_FACTOR_REAL(handle, MKL_DSS_DEFAULTS, KLL) !using KLL
                INFO = dsserror
                IF (dsserror /= MKL_DSS_SUCCESS) stop 'DSS error in Factoring:' 
                
                                         
         ELSEIF  (SPARSE_FLAVOR(1:7) == 'PARDISO') THEN
             WRITE(*,*) "Intel MKL Pardiso" 
                             
             IF (CRS_CCS == 'CCS') STOP 'CCS NOT YET in DSS'
             allocate( SOLN(NDOFL))       ! Solution
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

      ENDIF
   
      DO I=1,NDOFL                                         ! Make sure that scale factors are one. We don't want any equil scaling
         EQUIL_SCALE_FACS(I) = ONE                         ! of KLL in this subr. FBS below has EQUIL_SCALE_FACS as an input but
      ENDDO                                                ! they shouldn't be used as EQUED = 'N' is also input there (1st arg)

!***********************************************************************************************************************************
! Solve for DLR

! Open a scratch file that will be used to write DLR nonzero terms to as we solve for columns of DLR. After all col's
! of DLR have been solved for, and we have a count on NTERM_DLR, we will allocate memory to the DLR arrays and read
! the scratch file values into those arrays. Then, in the calling subroutine, we will write NTERM_DLR, followed by
! DLR  row/col/value to a permanent file

      OUNT(1) = ERR
      OUNT(2) = F06
      SCRFIL(1:)  = ' '
      SCRFIL(1:9) = 'SCRATCH-991'
      OPEN (SCR(1),STATUS='SCRATCH',FORM='UNFORMATTED',ACTION='READWRITE',IOSTAT=IOCHK)
      IF (IOCHK /= 0) THEN
         CALL OPNERR ( IOCHK, SCRFIL, OUNT, 'Y' )
         CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE', 'Y' )
         CALL OUTA_HERE ( 'Y' )                            ! Can't open scratch file, so quit
      ENDIF
 
! Loop on columns of KLR using rows of KRL
 
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      NTERM_DLR = 0
      DO J = 1,NDOFR

         CALL OURTIM
         MODNAM1 = '    Solve for DLR col '
         IF (NOCOUNTS /= 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') MODNAM1, J, NDOFR, CR13
         ENDIF

! To solve for the j-th col of DLR, use the j-th row of KRL (j-th col of KLR) as a rhs vector. Get the j-th row of KRL and put
! the negative of it into array INOUT_COL:

         NULL_COL = 'Y'
         DO I=1,NDOFL
            INOUT_COL(I) = ZERO
            DLR_COL(I)   = ZERO
         ENDDO 
         CALL GET_SPARSE_CRS_ROW ( 'KRL', J,  NTERM_KRL, NDOFR, NDOFL, I_KRL, J_KRL, KRL, -ONE, INOUT_COL, NULL_COL )

! Calculate DLR_COL via forward/backward substitution.

         IF (NULL_COL == 'N') THEN                         ! FBS will solve for DLR_COL & load it into DLR array
                                                           ! DPBTRS will return DLR_COL = -KLL(-1)*RHS_col
!                                                            Note 1st arg = 'N' assures that EQUIL_SCAL_FACS will not be used
            IF      (SOLLIB == 'BANDED') THEN

               CALL FBS_LAPACK ( 'N', NDOFL, KLL_SDIA, EQUIL_SCALE_FACS, INOUT_COL )

            ELSE IF (SOLLIB == 'SPARSE  ') THEN

               IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

                  INFO = 0
                  CALL FBS_SUPRLU ( SUBR_NAME, 'KLL', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, J, INOUT_COL, INFO )
      #ifdef MKLDSS
               ELSEIF  (SPARSE_FLAVOR(1:3) == 'DSS') THEN  !DSS STARTED
                  INFO = 0       
                    
                  dsserror = DSS_SOLVE_REAL(handle, MKL_DSS_DEFAULTS ,INOUT_COL ,1, SOLN)
                  INFO = dsserror
                  IF (dsserror /= MKL_DSS_SUCCESS) then 
                     stop 'DSS error in Solving :'
                  else
                    DO I=1,NDOFL
                        INOUT_COL(I) = SOLN(I)
                    ENDDO 
                    write (F06,9902)  'KLL','LINK3'
9902                FORMAT(' DSS FACTORIZATION OF MATRIX ', A, ' SUCCEEDED IN SUBR ', A)
                  endif 
                  
                              
               ELSEIF  (SPARSE_FLAVOR(1:7) == 'PARDISO') THEN  !pardiso STARTED
                
                !.. Back substitution and iterative refinement
                iparm(8) = 2 ! max numbers of iterative refinement steps
                phase = 33 ! only factorization
                
                soln = 0.d0
                CALL pardiso (pt, maxfct, mnum, mtype, phase, ndofl, KLL, I_KLL, J_KLL, idum, 1, iparm, msglvl,INOUT_COL, SOLN, pardisoerror)    
                
                    
                IF (pardisoerror /= 0) then 
                        
                    stop 'Pardiso error in Solving : '
                    
                else
                      
                    DO I=1,NDOFL
                         
                        INOUT_COL(I) = SOLN(I)
                      
                    ENDDO
                  
                    write (F06,9903)  'KLL','solve_dlr'
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

            DO I=1,NDOFL
               DLR_COL(I) = INOUT_COL(I)
            ENDDO
            DO I=1,NDOFL                                   ! Count NTERM_DLR and write nonzero DLR to scratch file
               IF (DABS(DLR_COL(I)) > EPS1) THEN
                  NTERM_DLR = NTERM_DLR + 1
                  WRITE(SCR(1)) I, J, DLR_COL(I)
               ENDIF
            ENDDO 
         ENDIF
  
      ENDDO
  
! The DLR data in SCRATCH-991 is written one col at a time for DLR. Therefore it is rows of DLRt

      REWIND (SCR(1))
      MESSAG = 'SCRATCH: DLR ROW/COL/VAL'
      READ_NTERM = 'N'
      OPND       = 'Y'
      CLOSE_IT   = 'N'
      CLOSE_STAT = 'KEEP    '

      CALL ALLOCATE_SPARSE_MAT ( 'DLRt', NDOFR, NTERM_DLR, SUBR_NAME )
      CALL ALLOCATE_L6_2 ( 'DLRt', SUBR_NAME )

      CALL ALLOCATE_SPARSE_MAT ( 'DLR', NDOFL, NTERM_DLR, SUBR_NAME )
      CALL ALLOCATE_L6_2 ( 'DLR', SUBR_NAME )
                                                           ! J_DLRt is same as I2_DLR and I2_DLRt is same as J_DLR
      CALL READ_MATRIX_2 ( SCRFIL, SCR(1), OPND, CLOSE_IT, CLOSE_STAT, MESSAG, 'DLRt', NDOFR, NTERM_DLR, READ_NTERM,               &
                           J_DLRt, I2_DLRt, DLRt)

! Now get DLR from DLRt

      CALL GET_I_MAT_FROM_I2_MAT ( 'DLRt', NDOFR, NTERM_DLR, I2_DLRt, I_DLRt )

      CALL MATTRNSP_SS ( NDOFR, NDOFL, NTERM_DLR, 'DLRt', I_DLRt, J_DLRt, DLRt, 'DLR', I_DLR, J_DLR, DLR )

      CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE', 'Y' )

! Print out constraint matrix DLR, if requested

      IF ( PRTDLR == 1) THEN
         IF (NTERM_DLR > 0) THEN
            CALL WRITE_SPARSE_CRS ( 'CB BOUNDARY MODE MATRIX DLR', 'L ', 'R ', NTERM_DLR, NDOFL, I_DLR, J_DLR, DLR )
         ENDIF
      ENDIF

! clearing---added

      IF (SOLLIB == 'SPARSE  ') THEN
               IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN
                   INFO = 0
                   CALL C_FORTRAN_DGSSV( 3, NDOFL, NTERM_KLL, 1, KLL , I_KLL , J_KLL , INOUT_COL, NDOFL, SLU_FACTORS, INFO )
                   IF (INFO .EQ. 0) THEN
                       WRITE (*,*) 'SUPERLU STORAGE FREED'
                   ELSE
                       WRITE(*,*) 'SUPERLU STORAGE NOT FREED. INFO FROM SUPERLU FREE STORAGE ROUTINE = ', INFO
                   ENDIF
      #ifdef MKLDSS                    
               ELSEIF  (SPARSE_FLAVOR(1:3) == 'DSS') THEN 
                    dsserror = DSS_DELETE(handle, MKL_DSS_DEFAULTS) !clearing
                    deallocate( SOLN)       ! Solution
                    IF (dsserror /= MKL_DSS_SUCCESS) STOP 'DSS error in CLEARING :' 
                                        
               ELSEIF  (SPARSE_FLAVOR(1:7) == 'PARDISO') THEN  !DSS STARTED
                   !.. Termination and release of memory
                   phase = -1 ! release internal memory          
                   CALL pardiso (pt, maxfct, mnum, mtype, phase, NDOFL, ddum, idum,  idum, idum, 1, iparm, msglvl, ddum, ddum, pardisoerror)        
                   deallocate(  SOLN)       ! SolutionENDIF
               endif     
            
      #endif MKLDSS
      ENDIF
      
               
                   
                   
                   

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF
      deallocate(EQUIL_SCALE_FACS , DLR_COL , INOUT_COL  )
      RETURN

! **********************************************************************************************************************************
  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

 2092 FORMAT(4X,A44,20X,I2,':',I2,':',I2,'.',I3)

 6001 FORMAT(' *ERROR  6001: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' MATRIX KLL WAS EQUILIBRATED: EQUED = ',A,'. CODE NOT WRITTEN TO ALLOW THIS AS YET')

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

12345 FORMAT(3X,A,I8,' of ',I8,A) 






! **********************************************************************************************************************************
 
      END SUBROUTINE SOLVE_DLR        
