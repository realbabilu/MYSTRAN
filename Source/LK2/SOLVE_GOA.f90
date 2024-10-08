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

 
      SUBROUTINE SOLVE_GOA
#ifdef MKLDSS
         use mkl_dss   
#endif

! Solves the sustem of equations: KOO*GOA = -KAO' for matrix GOA which is used in the reduction of the F set stiffness, mass and
! load matrices from the F-set to the A, O_sets
 
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, WRT_LOG, ERR, F04, F06, SCR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FACTORED_MATRIX, FATAL_ERR, KOO_SDIA, NDOFA, NDOFO, NTERM_GOA, NTERM_KOO,   &
                                         NTERM_KAO
      USE PARAMS, ONLY                :  EPSIL, PRTGOA, CRS_CCS
      USE TIMDAT, ONLY                :  TSEC
      USE SUBR_BEGEND_LEVELS, ONLY    :  SOLVE_GOA_BEGEND
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  SOLLIB, SPARSE_FLAVOR
      USE SPARSE_MATRICES, ONLY       :  I2_GOA, I_GOA, J_GOA, GOA, I_KOO, J_KOO, KOO, I_KAO, J_KAO, KAO
     USE LAPACK_LIN_EQN_DPB

! Interface module not needed for subr's DPBTRF and DPBTRS. These are "CONTAIN'ed" in module LAPACK_LIN_EQN_DPB, which
! is "USE'd" above

      USE SOLVE_GOA_USE_IFs

      IMPLICIT NONE

#ifdef MKLDSS
      include 'mkl_pardiso.fi'
#endif

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SOLVE_GOA'
      CHARACTER(  1*BYTE)             :: CLOSE_IT          ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not 
      CHARACTER(  8*BYTE)             :: CLOSE_STAT        ! What to do with file when it is closed
      CHARACTER( 24*BYTE)             :: MESSAG            ! File description. Input to subr UNFORMATTED_OPEN 
      CHARACTER( 22*BYTE)             :: MODNAM1           ! Name to write to screen to describe module being run
      CHARACTER(  1*BYTE)             :: READ_NTERM        ! 'Y' or 'N' Input to subr READ_MATRIX_1 
      CHARACTER(  1*BYTE)             :: NULL_COL          ! 'Y' if a col of KAO(transpose) is null 
      CHARACTER(  1*BYTE)             :: OPND              ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to open  a file or not 
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: SCRFIL            ! File name
 
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: INFO              ! Info on success of factorization or solve
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN   
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = SOLVE_GOA_BEGEND

      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE),allocatable                    :: GOA_COL(:)!(NDOFO)    ! A column of GOA solved for herein
      REAL(DOUBLE),allocatable                    :: INOUT_COL(:)!(NDOFO)  ! A column of KAO'
      REAL(DOUBLE),allocatable                    :: KOO_SCALE_FACS(:)!(NDOFO)
                                                           ! KOO scale facs. KOO will not be equilibrated so these are set to 1.0

#ifdef MKLDSS
      TYPE(MKL_DSS_HANDLE)            :: handle ! Allocate storage for the solver handle.      !DSS var
      INTEGER                         :: perm(1) ! DSS VAR
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
      
#endif
 
      INTRINSIC                       :: DABS

      allocate(GOA_COL(NDOFO) , INOUT_COL(NDOFO)  , KOO_SCALE_FACS(NDOFO))

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************
! Make units for writing errors the screen and output file
 
      OUNT(1) = ERR
      OUNT(2) = F06
 
! Make sure that ABAND (KOO triangular factor) was successfully generated when SYM_MAT_DECOMP_LAPACK was run.

      IF (SOLLIB == 'BANDED  ') THEN
         IF (FACTORED_MATRIX(1:3) /= 'KOO') THEN
            WRITE(ERR,2504) SUBR_NAME, FACTORED_MATRIX
            WRITE(F06,2504) SUBR_NAME, FACTORED_MATRIX
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF
      ENDIF

! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Solve for GOA

! Open a scratch file that will be used to write GOA nonzero terms to as we solve for columns of GOA. After all col's
! of GOA have been solved for, and we have a count on NTERM_GOA, we will allocate memory to the GOA arrays and read
! the scratch file values into those arrays. Then, in the calling subroutine, we will write NTERM_GOA, followed by
! GOA  row/col/value to a permanent file

      SCRFIL(1:)  = ' '
      SCRFIL(1:9) = 'SCRATCH-991'
      OPEN (SCR(1),STATUS='SCRATCH',FORM='UNFORMATTED',ACTION='READWRITE',IOSTAT=IOCHK)
      IF (IOCHK /= 0) THEN
         CALL OPNERR ( IOCHK, SCRFIL, OUNT, 'Y' )
         CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE', 'Y' )
         CALL OUTA_HERE ( 'Y' )                            ! Can't open scratch file, so quit
      ENDIF
 
! Loop on columns of KAO(transpose)
 
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      NTERM_GOA = 0
      DO J=1,NDOFA

         DO I=1,NDOFO
            INOUT_COL(I)       = ZERO                      ! Initialize INOUT_COL since GET_SPARSE_CRS_COL won't
            KOO_SCALE_FACS(I)  = ONE                       ! Initialize these
         ENDDO 
         NULL_COL = 'Y'
         CALL GET_SPARSE_CRS_ROW ( 'KAO',J, NTERM_KAO, NDOFA, NDOFO, I_KAO, J_KAO, KAO, -ONE, INOUT_COL, NULL_COL )

         IF (NULL_COL == 'N') THEN                         ! FBS will solve for GOA_COL & load it into GOA array

            CALL OURTIM
            MODNAM1 = '    Solve for GOA col '
            WRITE(SC1,22345) MODNAM1, J, NDOFA
                                                           ! FBS should not equilibrate since KOO was prevented from equilibrating
            IF      (SOLLIB == 'BANDED  ') THEN

               CALL FBS_LAPACK ( 'N', NDOFO, KOO_SDIA, KOO_SCALE_FACS, INOUT_COL )

            ELSE IF (SOLLIB == 'SPARSE  ') THEN
         
               IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

                  INFO = 0
                  CALL FBS_SUPRLU ( SUBR_NAME, 'KOO', NDOFO, NTERM_KOO, I_KOO, J_KOO, KOO, J, INOUT_COL, INFO )

#ifdef MKLDSS
               elseIF (SPARSE_FLAVOR(1:3) == 'DSS') THEN
                    IF (CRS_CCS == 'CCS') STOP 'CCS NOT YET'
                    allocate( SOLN(NDOFO))       ! Solution
                    dsserror = DSS_SOLVE_REAL(handle, MKL_DSS_DEFAULTS ,INOUT_COL ,1, SOLN)
                    INFO = dsserror
                    IF (dsserror /= MKL_DSS_SUCCESS) then 
                        stop 'DSS error in Solving :'
                    else
                      DO I=1,NDOFO
                         INOUT_COL(I) = SOLN(I)
                      ENDDO 
                    deallocate( SOLN)       ! Solution  
                    write (F06,9902)  'KOO','solve_Goa'
9902                FORMAT(' DSS SOLVING OF MATRIX ', A, ' SUCCEEDED IN SUBR ', A)
                    endif 
                    
               ELSEIF  (SPARSE_FLAVOR(1:7) == 'PARDISO') THEN
             
                   WRITE(*,*) "Intel MKL Pardiso"
                   IF (CRS_CCS == 'CCS') STOP 'CCS NOT YET'
                   allocate( SOLN(NDOFO))       ! Solution 
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

                   
                   !.. Back substitution and iterative refinement
                   iparm(8) = 2 ! max numbers of iterative refinement steps
                   phase = 33 ! only factorization
                   soln = 0.d0
                   CALL pardiso (pt, maxfct, mnum, mtype, phase, NDOFO, KOO, I_KOO, J_KOO, idum, 1, iparm, msglvl,INOUT_COL, SOLN, pardisoerror)    
    
                
                   IF (pardisoerror /= 0) then  
                    stop 'Pardiso error in Solving : '
                   else
                    DO I=1,NDOFO
                        INOUT_COL(I) = SOLN(I)   
                    ENDDO
                    deallocate( SOLN)       ! Solution   
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

            DO I=1,NDOFO
               GOA_COL(I) = INOUT_COL(I)
            ENDDO

            DO I=1,NDOFO                                   ! Count NTERM_GOA and write nonzero GOA to scratch file
               IF (DABS(GOA_COL(I)) > EPS1) THEN
                  NTERM_GOA = NTERM_GOA + 1
                  WRITE(SCR(1)) I,J,GOA_COL(I)
               ENDIF
            ENDDO 

         ELSE

            DO I=1,NDOFO
               GOA_COL(I) = ZERO
            ENDDO 

         ENDIF
  
      ENDDO

! The GOA data in SCRATCH-991 is written one col at a time. We need it to be  written for one row at a time with all
! rows in numerical order.

      CALL ALLOCATE_SPARSE_MAT ( 'GOA', NDOFO, NTERM_GOA, SUBR_NAME )
      CALL ALLOCATE_L2_GOA_2 ( SUBR_NAME )
      REWIND (SCR(1))
      MESSAG = 'SCRATCH: GOA ROW/COL/VAL'
      READ_NTERM = 'N'
      OPND       = 'Y'
      CLOSE_IT   = 'N'
      CLOSE_STAT = 'KEEP    '
      CALL READ_MATRIX_2 (SCRFIL, SCR(1), OPND, CLOSE_IT, CLOSE_STAT, MESSAG,'GOA',NDOFO, NTERM_GOA, READ_NTERM, I2_GOA, J_GOA, GOA)
      CALL SORT_INT2_REAL1 (SUBR_NAME, 'I2_GOA, J_GOA, GOA', NTERM_GOA, I2_GOA, J_GOA, GOA )
      REWIND (SCR(1))
      WRITE(SCR(1)) NTERM_GOA
      DO K=1,NTERM_GOA
         WRITE(SCR(1)) I2_GOA(K),J_GOA(K),GOA(K)
      ENDDO  

! Reallocate memory to GOA based on the NTERM_GOA counted above and read values from scratch file into GOA arrays

      WRITE(SC1, * ) '    Reallocate GOA'
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages         
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GOA'     , CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'GOA' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate L2_GOA_2', CR13   ;   CALL DEALLOCATE_L2_GOA_2
      WRITE(SC1,12345,ADVANCE='NO') '       Allocate   GOA'     , CR13
      CALL ALLOCATE_SPARSE_MAT ('GOA', NDOFO, NTERM_GOA, SUBR_NAME)

      REWIND (SCR(1))
      SCRFIL(1:)  = ' '
      SCRFIL(1:9) = 'SCRATCH-991'
      MESSAG = 'SCRATCH: GOA ROW/COL/VAL'
      READ_NTERM = 'Y'
      OPND       = 'Y'
      CLOSE_IT   = 'Y'
      CLOSE_STAT = 'DELETE  '
      CALL READ_MATRIX_1 ( SCRFIL, SCR(1), OPND, CLOSE_IT, CLOSE_STAT, MESSAG, 'GOA', NTERM_GOA, READ_NTERM, NDOFO,                &
                           I_GOA, J_GOA, GOA)

      CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE', 'Y' )

! Print out constraint matrix GOA, if requested

      IF ( PRTGOA == 1) THEN
         IF (NTERM_GOA > 0) THEN
            CALL WRITE_SPARSE_CRS ( 'CONSTRAINT MATRIX GOA', 'O ', 'A ', NTERM_GOA, NDOFO, I_GOA, J_GOA, GOA )
         ENDIF
      ENDIF

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      deallocate   (GOA_COL , INOUT_COL  , KOO_SCALE_FACS)

      RETURN

! **********************************************************************************************************************************
 2504 FORMAT(' *ERROR  2504: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NAME OF MATRIX THAT HAS BEEN DECOMPOSED INTO TRIANGULAR FACTORS SHOULD BE "KOO".'                 &
                    ,/,14X,',HOWEVER, IT IS NAMED "',A,'". CANNOT CONTINUE')

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

12345 FORMAT(A,10X,A)

22345 FORMAT(3X,A,I8,' of ',I8,A) 

! **********************************************************************************************************************************
 
      END SUBROUTINE SOLVE_GOA        
