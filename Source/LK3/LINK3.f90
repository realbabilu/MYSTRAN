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

! LINK 3 solves the equation KLL*UL = PL where KLL, UL, PL are the L-set stiffness matrix, displs and loads. It solves the equation
! using one of three methods. For each method the solution is obtained in a 2 step process: (1) the KLL matrix is decomposed into
! triangular factors and (2) UL is solved for by forward-backward substitution (FBS). The 3 methods are:

!   a) The LAPACK freeware code. This code requires KLL to be in banded (NOT sparse) form. LAPACH has the advantage that
!      MYSTRAN contains the LAPACK source code so debugging is easy. Its disadvantage is that banded matrices require much more
!      memory than sparse storage for large stiffness matrices.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, ERR, F06, INFILE, L3A, SC1, LINK3A, L3A_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, COMM, FATAL_ERR, KLL_SDIA, LINKNO, MBUG, NDOFL, NSUB,                       &
                                         NTERM_KLL, NTERM_PL, RESTART,  SOL_NAME, WARN_ERR
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, TEN
      USE PARAMS, ONLY                :  CRS_CCS, EPSERR, EPSIL, KLLRAT, RELINK3, RCONDK, SOLLIB, SUPWARN, SPARSE_FLAVOR
      USE SPARSE_MATRICES, ONLY       :  I_KLL, J_KLL, KLL, I_PL, J_PL, PL, SYM_KLL
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
      USE LINK_MESSAGE_Interface

      IMPLICIT NONE

#ifdef DMUMPS_Solver
      INCLUDE 'dmumps_struc.h'
      INCLUDE 'mpif.h'
#endif

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LINK3'
      CHARACTER(  2*BYTE)             :: L_SET    = 'L '   ! L-set designator
      CHARACTER(  1*BYTE)             :: EQUED             ! 'Y' if the stiff matrix was equilibrated in subr EQUILIBRATE
      CHARACTER(  1*BYTE)             :: NULL_COL          ! 'Y' if a col of KAO(transpose) is null

      INTEGER(LONG)                   :: DEB_PRT(2)        ! Debug numbers to say whether to write ABAND and/or its decomp to output
!                                                            file in called subr SYM_MAT_DECOMP_LAPACK (ABAND = band form of KLL)

      INTEGER(LONG)                   :: IER_DECOMP        ! Overall error indicator
      INTEGER(LONG)                   :: ISUB              ! DO loop index for subcases
      INTEGER(LONG)                   :: INFO     = 0      ! Info output from some routine that has been called
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG), PARAMETER        :: P_LINKNO = 2      ! Prior LINK no's that should have run before this LINK can execute

      REAL(DOUBLE)                    :: BETA              ! Multiple for rhs for use in subr FBS
      REAL(DOUBLE)                    :: DEN               ! K_INORM*UL_INORM + PL_INORM
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero

      REAL(DOUBLE)                    :: EQUIL_SCALE_FACS(NDOFL)
                                                           ! LAPACK_S values returned from subr SYM_MAT_DECOMP_LAPACK
      REAL(DOUBLE)                    :: DUM_COL(NDOFL)    ! Temp variable used in SuperLU
      REAL(DOUBLE)                    :: K_INORM           ! Inf norm of KLL matrix (det in  subr COND_NUM)
      REAL(DOUBLE)                    :: LAP_ERR1          ! Bound on displ error = 2*OMEGAI/RCOND
      REAL(DOUBLE)                    :: OMEGAI            ! RES_INORM/DEN (similar to EPSILON)
      REAL(DOUBLE)                    :: OMEGAI0           ! Upper bound on OMEGAI. OMEGAI0 = 10*NDOFL*MACH_EPS
      REAL(DOUBLE)                    :: PL_INORM          ! Inf norm of load vector
      REAL(DOUBLE)                    :: RES_INORM         ! Inf norm of residual vector R = K*UL - PL
      REAL(DOUBLE)                    :: RCOND             ! Recrip of cond no. of the KLL. Det in  subr COND_NUM
      REAL(DOUBLE)                    :: UL_INORM          ! Inf norm of displacement vector

#ifdef DMUMPS_Solver
      TYPE(DMUMPS_STRUC)              :: MUMPS_PAR
      INTEGER                         :: MUMPS_COMM
      INTEGER                         :: MUMPS_INFOG1
      INTEGER                         :: MUMPS_MPI_IERR
      INTEGER, ALLOCATABLE, TARGET    :: MUMPS_IRN(:), MUMPS_JCN(:)
      REAL(DOUBLE), ALLOCATABLE, TARGET :: MUMPS_A(:)
      LOGICAL                         :: MUMPS_ACTIVE
      LOGICAL                         :: MUMPS_MPI_ACTIVE
#endif

      INTRINSIC                       :: DABS

! --- solverbattle_bridge begin --- !
      CHARACTER(1*BYTE)              :: BENCH_EXPORT = 'N'
      CHARACTER(LEN=256)             :: BENCH_UL_FILE = ''
      INTEGER(LONG)                  :: BENCH_UL_UNT = 0
! --- solverbattle_bridge end --- !

!***********************************************************************************************************************************
      LINKNO = 3

      EPS1 = EPSIL(1)

#ifdef DMUMPS_Solver
      MUMPS_COMM       = 0
      MUMPS_INFOG1     = 0
      MUMPS_MPI_IERR   = 0
      MUMPS_ACTIVE     = .FALSE.
      MUMPS_MPI_ACTIVE = .FALSE.
#endif

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

! --- BANDED_optimizisation -begin-- !
      CALL REPORT_SOLVER_DISPATCH_POLICY ( 'KLL', SUBR_NAME )
! --- BANDED_optimizisation -end-- !

Factr:IF (SOLLIB == 'BANDED  ') THEN                       ! Use LAPACK

         INFO = 0
         CALL SYM_MAT_DECOMP_LAPACK ( SUBR_NAME, 'KLL', L_SET, NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, 'Y', KLLRAT, 'Y', RCONDK,      &
                                      DEB_PRT, EQUED, KLL_SDIA, K_INORM, RCOND, EQUIL_SCALE_FACS, INFO )

      ELSE IF (SOLLIB == 'SPARSE  ') THEN

         IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

            SLU_INFO = 0
            CALL SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME, 'KLL', L_SET, NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, SLU_INFO )

         ELSE IF (SPARSE_FLAVOR(1:5) == 'MUMPS') THEN

#ifdef DMUMPS_Solver
            CALL SYM_MAT_DECOMP_MUMPS ( INFO )
#else
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,9992) SUBR_NAME, 'SPARSE_FLAVOR', 'MUMPS'
            WRITE(F06,9992) SUBR_NAME, 'SPARSE_FLAVOR', 'MUMPS'
            CALL OUTA_HERE ( 'Y' )
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

      ENDIF Factr

!***********************************************************************************************************************************
!  Allocate col vector arrays for loads, displs and res vector

!xx   CALL ALLOCATE_COL_VEC ( 'UL_COL', NDOFL, SUBR_NAME )
!xx   CALL ALLOCATE_COL_VEC ( 'PL_COL', NDOFL, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'RES', NDOFL, 1, SUBR_NAME )

! Open file for writing displs to.

      CALL FILE_OPEN ( L3A, LINK3A, OUNT, 'REPLACE', L3A_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

! --- solverbattle_bridge begin --- !
      IF (DEBUG(211) > 0) THEN
         BENCH_EXPORT = 'Y'
         CALL LINK_MESSAGE('EXPORT KLL/PL BENCHMARK INPUTS                         ')
         CALL EXPORT_LINK3_BENCHMARK_INPUTS
      ENDIF
! --- solverbattle_bridge end --- !

! Loop on subcases

      WRITE(F06,*)
      BETA = ONE
Solve:DO ISUB = 1,NSUB

         SLU_INFO = 0
         CALL ALLOCATE_COL_VEC ( 'UL_COL', NDOFL, SUBR_NAME )
         CALL ALLOCATE_COL_VEC ( 'PL_COL', NDOFL, SUBR_NAME )

                                                           ! Get the loads for this subcase from I_PL, J_PL, PL and put into PL_COL
         CALL LINK_MESSAGE_I('GET COL OF PL LOADS FOR                        Subcase', ISUB)
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

                                                           ! Call FBS to solve for displacements for this subcase
         CALL LINK_MESSAGE_I('FBS - SOLVE FOR RHS ANSWERS FOR                   "', ISUB)
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages

         IF      (SOLLIB == 'BANDED  ') THEN

            CALL FBS_LAPACK ( EQUED, NDOFL, KLL_SDIA, EQUIL_SCALE_FACS, DUM_COL )

         ELSE IF (SOLLIB == 'SPARSE  ') THEN

            IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

               SLU_INFO = 0
               CALL FBS_SUPRLU ( SUBR_NAME, 'KLL', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, ISUB, DUM_COL, SLU_INFO )

            ELSE IF (SPARSE_FLAVOR(1:5) == 'MUMPS') THEN

#ifdef DMUMPS_Solver
               CALL FBS_MUMPS ( DUM_COL, INFO )
#else
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,9992) SUBR_NAME, 'SPARSE_FLAVOR', 'MUMPS'
               WRITE(F06,9992) SUBR_NAME, 'SPARSE_FLAVOR', 'MUMPS'
               CALL OUTA_HERE ( 'Y' )
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

         DO J=1,NDOFL
            UL_COL(J) = DUM_COL(J)
         ENDDO

         IF (DEBUG(33) == 1) THEN                          ! DEBUG output of displs
            WRITE(F06,3022) ISUB
            CALL WRITE_VECTOR ( '      A-SET DISPL      ','DISPL', NDOFL, UL_COL )
            WRITE(F06,*)
         ENDIF

         IF (EPSERR == 'Y') THEN                           ! Calculate residual vector, R. Use RES to calculate EPSILON
            CALL LINK_MESSAGE_I('CALC  EPSILON ERROR ESTIMATE                      "', ISUB)
            CALL EPSCALC ( ISUB )
         ENDIF
                                                           ! Calculate the LAPACK error bounds
         IF ((RCONDK == 'Y') .AND. (SOLLIB == 'BANDED')) THEN
            IF (DABS(RCOND) > MACH_SFMIN) THEN
               CALL LINK_MESSAGE_I('CALC LAPACK ERROR ESTIMATE                        "', ISUB)
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

! --- solverbattle_bridge begin --- !
         IF (BENCH_EXPORT == 'Y') THEN
            CALL WRITE_LINK3_BENCHMARK_UL_COL ( UL_COL )
         ENDIF
! --- solverbattle_bridge end --- !

         CALL DEALLOCATE_COL_VEC  ( 'UL_COL' )
         CALL DEALLOCATE_COL_VEC  ( 'PL_COL' )


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

         ENDIF

      ENDIF FreeS

#ifdef DMUMPS_Solver
      IF ((SOLLIB == 'SPARSE  ') .AND. (SPARSE_FLAVOR(1:5) == 'MUMPS')) THEN
         CALL FREE_MUMPS_FACTORS
      ENDIF
#endif

! Dellocate arrays

      CALL LINK_MESSAGE('DEALLOCATE ARRAYS')
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

! --- solverbattle_bridge begin --- !
      IF ((BENCH_EXPORT == 'Y') .AND. (BENCH_UL_UNT > 0)) THEN
         CLOSE(BENCH_UL_UNT)
         BENCH_UL_UNT = 0
      ENDIF
! --- solverbattle_bridge end --- !

      CALL FILE_CLOSE ( L3A, LINK3A, 'KEEP' )

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

! Write LINK3 end to F06

      CALL OURTIM
      WRITE(F06,151) LINKNO

! Close files

      IF (( DEBUG(193) == 3) .OR. (DEBUG(193) == 999)) THEN
         CALL FILE_INQUIRE ( 'near end of LINK3' )
      ENDIF

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

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

 9992 FORMAT(' *ERROR  9992: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A,' = ',A,' WAS REQUESTED BUT THIS BUILD WAS NOT COMPILED WITH DMUMPS_Solver.')

 9998 FORMAT(' *ERROR  9998: COMM ',I3,' INDICATES UNSUCCESSFUL LINK ',I2,' COMPLETION.'                                           &
                    ,/,14X,' FATAL ERROR - CANNOT START LINK ',I2)

12345 FORMAT(A,10X,A)

!***********************************************************************************************************************************

      CONTAINS

! --- solverbattle_bridge begin --- !
! **********************************************************************************************************************************
      SUBROUTINE EXPORT_LINK3_BENCHMARK_INPUTS

      CHARACTER(LEN=256)             :: BENCH_PREFIX
      CHARACTER(LEN=256)             :: KLL_FILE
      CHARACTER(LEN=256)             :: META_FILE
      CHARACTER(LEN=256)             :: PL_FILE
      CHARACTER(LEN=256)             :: UL_FILE
      INTEGER(LONG)                  :: UNT

      BENCH_PREFIX = GET_BENCHMARK_PREFIX()
      KLL_FILE     = TRIM(BENCH_PREFIX) // '_kll.mtx'
      PL_FILE      = TRIM(BENCH_PREFIX) // '_pl.mtx'
      UL_FILE      = TRIM(BENCH_PREFIX) // '_ul.mtx'
      META_FILE    = TRIM(BENCH_PREFIX) // '_benchmark_meta.txt'

      CALL WRITE_MATRIX_MARKET_CRS ( KLL_FILE, NDOFL, NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, SYM_KLL )
      CALL WRITE_MATRIX_MARKET_CRS ( PL_FILE , NDOFL, NSUB , NTERM_PL , I_PL , J_PL , PL  )
      CALL OPEN_MATRIX_MARKET_ARRAY ( UL_FILE, NDOFL, NSUB, BENCH_UL_UNT )
      BENCH_UL_FILE = UL_FILE

      OPEN(NEWUNIT=UNT, FILE=META_FILE, STATUS='REPLACE', ACTION='WRITE')
      WRITE(UNT,'(A)') '# MYSTRAN LINK3 static benchmark export'
      WRITE(UNT,'(A)') 'input_file=' // TRIM(INFILE)
      WRITE(UNT,'(A)') 'matrix_file=' // TRIM(KLL_FILE)
      WRITE(UNT,'(A)') 'rhs_file=' // TRIM(PL_FILE)
      WRITE(UNT,'(A)') 'solution_file=' // TRIM(UL_FILE)
      WRITE(UNT,'(A)') 'matrix_set=KLL'
      WRITE(UNT,'(A)') 'rhs_set=PL'
      WRITE(UNT,'(A)') 'solution_set=UL'
      WRITE(UNT,'(A,I0)') 'ndofl=', NDOFL
      WRITE(UNT,'(A,I0)') 'nsub=', NSUB
      WRITE(UNT,'(A,A)') 'matrix_sparse_storage_symmetry=', SYM_KLL
      WRITE(UNT,'(A)') 'matrix_market_kll_export=expanded_full_general_when_sym_storage'
      WRITE(UNT,'(A)') 'matrix_market_ul_export=array_dense_column_major'
      WRITE(UNT,'(A)') 'matrix_market_rhs_columns_are_internal_subcase_numbers'
      CLOSE(UNT)

      WRITE(F06,'(/,A)') ' LINK3 benchmark export: Matrix Market files written for reduced static system.'
      WRITE(F06,'(A)')   '   KLL : ' // TRIM(KLL_FILE)
      WRITE(F06,'(A)')   '   PL  : ' // TRIM(PL_FILE)
      WRITE(F06,'(A)')   '   UL  : ' // TRIM(UL_FILE)
      WRITE(F06,'(A)')   '   META: ' // TRIM(META_FILE)
      WRITE(F06,*)

      END SUBROUTINE EXPORT_LINK3_BENCHMARK_INPUTS

! **********************************************************************************************************************************
      SUBROUTINE WRITE_MATRIX_MARKET_CRS ( FILNAM, NROWS, NCOLS, NTERM, I_MAT, J_MAT, MAT, SYM_STORAGE )

      CHARACTER(LEN=*), INTENT(IN)   :: FILNAM
      INTEGER(LONG), INTENT(IN)      :: NCOLS
      INTEGER(LONG), INTENT(IN)      :: NROWS
      INTEGER(LONG), INTENT(IN)      :: NTERM
      INTEGER(LONG), INTENT(IN)      :: I_MAT(NROWS+1)
      INTEGER(LONG), INTENT(IN)      :: J_MAT(NTERM)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: SYM_STORAGE
      CHARACTER(1*BYTE)              :: SYM_STOR
      INTEGER(LONG)                  :: NDIAG
      INTEGER(LONG)                  :: NOUT
      INTEGER(LONG)                  :: IROW
      INTEGER(LONG)                  :: K
      INTEGER(LONG)                  :: UNT
      REAL(DOUBLE), INTENT(IN)       :: MAT(NTERM)

      SYM_STOR = 'N'
      IF (PRESENT(SYM_STORAGE)) SYM_STOR = SYM_STORAGE(1:1)

      OPEN(NEWUNIT=UNT, FILE=FILNAM, STATUS='REPLACE', ACTION='WRITE')
      WRITE(UNT,'(A)') '%%MatrixMarket matrix coordinate real general'
      WRITE(UNT,'(A)') '% exported by MYSTRAN LINK3 benchmark bridge'
      IF ((SYM_STOR == 'Y') .AND. (NROWS == NCOLS)) THEN
         NDIAG = 0
         DO IROW=1,NROWS
            DO K=I_MAT(IROW),I_MAT(IROW+1)-1
               IF (J_MAT(K) == IROW) NDIAG = NDIAG + 1
            ENDDO
         ENDDO
         NOUT = 2*NTERM - NDIAG
         WRITE(UNT,'(I0,1X,I0,1X,I0)') NROWS, NCOLS, NOUT
         DO IROW=1,NROWS
            DO K=I_MAT(IROW),I_MAT(IROW+1)-1
               WRITE(UNT,'(I0,1X,I0,1X,ES25.16)') IROW, J_MAT(K), MAT(K)
               IF (J_MAT(K) /= IROW) THEN
                  WRITE(UNT,'(I0,1X,I0,1X,ES25.16)') J_MAT(K), IROW, MAT(K)
               ENDIF
            ENDDO
         ENDDO
      ELSE
         WRITE(UNT,'(I0,1X,I0,1X,I0)') NROWS, NCOLS, NTERM
         DO IROW=1,NROWS
            DO K=I_MAT(IROW),I_MAT(IROW+1)-1
               WRITE(UNT,'(I0,1X,I0,1X,ES25.16)') IROW, J_MAT(K), MAT(K)
            ENDDO
         ENDDO
      ENDIF
      CLOSE(UNT)

      END SUBROUTINE WRITE_MATRIX_MARKET_CRS

! **********************************************************************************************************************************
      SUBROUTINE OPEN_MATRIX_MARKET_ARRAY ( FILNAM, NROWS, NCOLS, UNT )

      CHARACTER(LEN=*), INTENT(IN)   :: FILNAM
      INTEGER(LONG), INTENT(IN)      :: NROWS
      INTEGER(LONG), INTENT(IN)      :: NCOLS
      INTEGER(LONG), INTENT(OUT)     :: UNT

      OPEN(NEWUNIT=UNT, FILE=FILNAM, STATUS='REPLACE', ACTION='WRITE')
      WRITE(UNT,'(A)') '%%MatrixMarket matrix array real general'
      WRITE(UNT,'(A)') '% exported by MYSTRAN LINK3 benchmark bridge'
      WRITE(UNT,'(I0,1X,I0)') NROWS, NCOLS

      END SUBROUTINE OPEN_MATRIX_MARKET_ARRAY

! **********************************************************************************************************************************
      SUBROUTINE WRITE_LINK3_BENCHMARK_UL_COL ( ULVEC )

      REAL(DOUBLE), INTENT(IN)       :: ULVEC(NDOFL)
      INTEGER(LONG)                  :: I

      DO I=1,NDOFL
         WRITE(BENCH_UL_UNT,'(ES25.16)') ULVEC(I)
      ENDDO

      END SUBROUTINE WRITE_LINK3_BENCHMARK_UL_COL
! --- solverbattle_bridge end --- !

#ifdef DMUMPS_Solver
! **********************************************************************************************************************************
      SUBROUTINE SYM_MAT_DECOMP_MUMPS ( INFO_OUT )

      INTEGER(LONG), INTENT(OUT)     :: INFO_OUT
      INTEGER(LONG)                  :: K
      INTEGER(LONG)                  :: IROW

      CALL INIT_MUMPS_RUNTIME
      CALL BUILD_COO_FROM_CRS

      MUMPS_PAR%JOB = -1
      MUMPS_PAR%PAR = 1
      MUMPS_PAR%SYM = 0
      IF (SYM_KLL == 'Y') MUMPS_PAR%SYM = 2
      MUMPS_PAR%COMM = MUMPS_COMM
      CALL DMUMPS(MUMPS_PAR)

      MUMPS_PAR%ICNTL(1) = -1
      MUMPS_PAR%ICNTL(2) = -1
      MUMPS_PAR%ICNTL(3) = -1
      MUMPS_PAR%ICNTL(4) = 0
      MUMPS_PAR%N   = NDOFL
      MUMPS_PAR%NZ  = SIZE(MUMPS_A)
      MUMPS_PAR%IRN => MUMPS_IRN
      MUMPS_PAR%JCN => MUMPS_JCN
      MUMPS_PAR%A   => MUMPS_A
      MUMPS_PAR%JOB = 4
      CALL DMUMPS(MUMPS_PAR)
      MUMPS_INFOG1 = MUMPS_PAR%INFOG(1)
      INFO_OUT = MUMPS_INFOG1

      IF (INFO_OUT < 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A,A,I12,A,A,A)') ' *ERROR  9811: THE FACTORIZATION OF THE MATRIX ', 'KLL',                                &
                                        ' BY MUMPS HAD ERROR WITH INFOG(1) = ', INFO_OUT, ' IN SUBR ', SUBR_NAME, '.'
         WRITE(F06,'(A,A,A,I12,A,A,A)') ' *ERROR  9811: THE FACTORIZATION OF THE MATRIX ', 'KLL',                                &
                                        ' BY MUMPS HAD ERROR WITH INFOG(1) = ', INFO_OUT, ' IN SUBR ', SUBR_NAME, '.'
         CALL OUTA_HERE ( 'Y' )
      ELSE
         MUMPS_ACTIVE = .TRUE.
         WRITE(F06,'(A,A,A,A)') ' MUMPS FACTORIZATION OF MATRIX ', 'KLL', ' SUCCEEDED IN SUBR ', SUBR_NAME
      ENDIF

      END SUBROUTINE SYM_MAT_DECOMP_MUMPS

! **********************************************************************************************************************************
      SUBROUTINE FBS_MUMPS ( RHS_COL, INFO_OUT )

      REAL(DOUBLE), TARGET, INTENT(INOUT) :: RHS_COL(NDOFL)
      INTEGER(LONG), INTENT(OUT)          :: INFO_OUT

      IF (.NOT. MUMPS_ACTIVE) THEN
         INFO_OUT = -999
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A,A)') ' *ERROR  9812: MUMPS SOLVE WAS REQUESTED IN SUBR ', SUBR_NAME,                                    &
                              ' BEFORE FACTORIZATION WAS AVAILABLE.'
         WRITE(F06,'(A,A,A)') ' *ERROR  9812: MUMPS SOLVE WAS REQUESTED IN SUBR ', SUBR_NAME,                                    &
                              ' BEFORE FACTORIZATION WAS AVAILABLE.'
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      MUMPS_PAR%NRHS = 1
      MUMPS_PAR%LRHS = NDOFL
      MUMPS_PAR%RHS  => RHS_COL
      MUMPS_PAR%JOB  = 3
      CALL DMUMPS(MUMPS_PAR)
      INFO_OUT = MUMPS_PAR%INFOG(1)

      IF (INFO_OUT < 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,I12,A,A,A)') ' *ERROR  9813: MUMPS SOLVE FAILED WITH INFOG(1) = ', INFO_OUT,                              &
                                     ' IN SUBR ', SUBR_NAME, '.'
         WRITE(F06,'(A,I12,A,A,A)') ' *ERROR  9813: MUMPS SOLVE FAILED WITH INFOG(1) = ', INFO_OUT,                              &
                                     ' IN SUBR ', SUBR_NAME, '.'
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      END SUBROUTINE FBS_MUMPS

! **********************************************************************************************************************************
      SUBROUTINE FREE_MUMPS_FACTORS

      IF (MUMPS_ACTIVE) THEN
         MUMPS_PAR%JOB = -2
         CALL DMUMPS(MUMPS_PAR)
         MUMPS_ACTIVE = .FALSE.
      ENDIF

      IF (ALLOCATED(MUMPS_IRN)) DEALLOCATE(MUMPS_IRN)
      IF (ALLOCATED(MUMPS_JCN)) DEALLOCATE(MUMPS_JCN)
      IF (ALLOCATED(MUMPS_A  )) DEALLOCATE(MUMPS_A  )

      IF (MUMPS_MPI_ACTIVE) THEN
         CALL MPI_FINALIZE(MUMPS_MPI_IERR)
         MUMPS_MPI_ACTIVE = .FALSE.
      ENDIF

      END SUBROUTINE FREE_MUMPS_FACTORS

! **********************************************************************************************************************************
      SUBROUTINE INIT_MUMPS_RUNTIME

      CALL MPI_INIT(MUMPS_MPI_IERR)
      IF (MUMPS_MPI_IERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,I12,A,A,A)') ' *ERROR  9814: MPI_INIT FAILED WITH IERR = ', MUMPS_MPI_IERR,                               &
                                     ' IN SUBR ', SUBR_NAME, '.'
         WRITE(F06,'(A,I12,A,A,A)') ' *ERROR  9814: MPI_INIT FAILED WITH IERR = ', MUMPS_MPI_IERR,                               &
                                     ' IN SUBR ', SUBR_NAME, '.'
         CALL OUTA_HERE ( 'Y' )
      ENDIF
      MUMPS_MPI_ACTIVE = .TRUE.
      MUMPS_COMM = MPI_COMM_WORLD

      END SUBROUTINE INIT_MUMPS_RUNTIME

! **********************************************************************************************************************************
      SUBROUTINE BUILD_COO_FROM_CRS

      INTEGER(LONG)                  :: IDX
      INTEGER(LONG)                  :: IROW
      INTEGER(LONG)                  :: K

      IF (ALLOCATED(MUMPS_IRN)) DEALLOCATE(MUMPS_IRN)
      IF (ALLOCATED(MUMPS_JCN)) DEALLOCATE(MUMPS_JCN)
      IF (ALLOCATED(MUMPS_A  )) DEALLOCATE(MUMPS_A  )

      ALLOCATE(MUMPS_IRN(NTERM_KLL))
      ALLOCATE(MUMPS_JCN(NTERM_KLL))
      ALLOCATE(MUMPS_A  (NTERM_KLL))

      IDX = 0
      DO IROW=1,NDOFL
         DO K=I_KLL(IROW),I_KLL(IROW+1)-1
            IDX = IDX + 1
            MUMPS_IRN(IDX) = IROW
            MUMPS_JCN(IDX) = J_KLL(K)
            MUMPS_A(IDX)   = KLL(K)
         ENDDO
      ENDDO

      END SUBROUTINE BUILD_COO_FROM_CRS
#endif

! **********************************************************************************************************************************
      FUNCTION GET_BENCHMARK_PREFIX() RESULT(PREFIX)

      CHARACTER(LEN=256)             :: PREFIX
      INTEGER(LONG)                  :: I
      INTEGER(LONG)                  :: LAST_DOT
      INTEGER(LONG)                  :: LAST_SEP
      INTEGER(LONG)                  :: NAME_END

      PREFIX   = 'mystran_link3'
      LAST_SEP = 0
      LAST_DOT = 0
      NAME_END = LEN_TRIM(INFILE)

      IF (NAME_END <= 0) RETURN

      DO I=1,NAME_END
         IF ((INFILE(I:I) == '\') .OR. (INFILE(I:I) == '/')) THEN
            LAST_SEP = I
            LAST_DOT = 0
         ELSE IF (INFILE(I:I) == '.') THEN
            LAST_DOT = I
         ENDIF
      ENDDO

      IF ((LAST_DOT > LAST_SEP + 1) .AND. (LAST_DOT <= NAME_END)) THEN
         PREFIX = ADJUSTL(INFILE(LAST_SEP+1:LAST_DOT-1))
      ELSE
         PREFIX = ADJUSTL(INFILE(LAST_SEP+1:NAME_END))
      ENDIF

      IF (LEN_TRIM(PREFIX) <= 0) PREFIX = 'mystran_link3'

      END FUNCTION GET_BENCHMARK_PREFIX

! **********************************************************************************************************************************

      END SUBROUTINE LINK3
