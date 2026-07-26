! --- MUMPS_COO add begin --- !
      MODULE DMUMPS_STUF

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE

      IMPLICIT NONE

#ifdef DMUMPS_Solver
      INCLUDE 'dmumps_struc.h'
      INCLUDE 'mpif.h'

      TYPE(DMUMPS_STRUC)              :: DMUMPS_PAR
      INTEGER                         :: DMUMPS_COMM = 0
      INTEGER                         :: DMUMPS_MPI_IERR = 0
      INTEGER, ALLOCATABLE, TARGET    :: DMUMPS_IRN(:), DMUMPS_JCN(:)
      REAL(DOUBLE), ALLOCATABLE, TARGET :: DMUMPS_A(:)
      LOGICAL                         :: DMUMPS_ACTIVE = .FALSE.
      LOGICAL                         :: DMUMPS_MPI_ACTIVE = .FALSE.
#endif

      CONTAINS

      LOGICAL FUNCTION DMUMPS_COMPILED_IN()

#ifdef DMUMPS_Solver
      DMUMPS_COMPILED_IN = .TRUE.
#else
      DMUMPS_COMPILED_IN = .FALSE.
#endif

      END FUNCTION DMUMPS_COMPILED_IN

      LOGICAL FUNCTION DMUMPS_CRS_IS_NUMERICALLY_SYMMETRIC ( N, NTERM, I_CRS, J_CRS, A_CRS )

      INTEGER(LONG), INTENT(IN)       :: N
      INTEGER(LONG), INTENT(IN)       :: NTERM
      INTEGER(LONG), INTENT(IN)       :: I_CRS(N+1)
      INTEGER(LONG), INTENT(IN)       :: J_CRS(NTERM)
      REAL(DOUBLE), INTENT(IN)        :: A_CRS(NTERM)

      INTEGER(LONG)                   :: IROW
      INTEGER(LONG)                   :: JCOL
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: KK
      LOGICAL                         :: FOUND_PAIR
      REAL(DOUBLE)                    :: SCALE
      REAL(DOUBLE), PARAMETER         :: SYM_TOL = 1.0D-10

      DMUMPS_CRS_IS_NUMERICALLY_SYMMETRIC = .TRUE.

      DO IROW=1,N
         DO K=I_CRS(IROW),I_CRS(IROW+1)-1
            JCOL = J_CRS(K)
            IF (JCOL == IROW) CYCLE

            FOUND_PAIR = .FALSE.
            DO KK=I_CRS(JCOL),I_CRS(JCOL+1)-1
               IF (J_CRS(KK) == IROW) THEN
                  FOUND_PAIR = .TRUE.
                  SCALE = MAX(1.0D0, ABS(A_CRS(K)), ABS(A_CRS(KK)))
                  IF (ABS(A_CRS(K) - A_CRS(KK)) > SYM_TOL*SCALE) THEN
                     DMUMPS_CRS_IS_NUMERICALLY_SYMMETRIC = .FALSE.
                     RETURN
                  ENDIF
                  EXIT
               ENDIF
            ENDDO

            IF (.NOT. FOUND_PAIR) THEN
               DMUMPS_CRS_IS_NUMERICALLY_SYMMETRIC = .FALSE.
               RETURN
            ENDIF
         ENDDO
      ENDDO

      END FUNCTION DMUMPS_CRS_IS_NUMERICALLY_SYMMETRIC

      SUBROUTINE DMUMPS_FACTOR_CRS ( N, NTERM, I_CRS, J_CRS, A_CRS, SYM_FLAG, INFO_OUT )

      INTEGER(LONG), INTENT(IN)       :: N
      INTEGER(LONG), INTENT(IN)       :: NTERM
      INTEGER(LONG), INTENT(IN)       :: I_CRS(N+1)
      INTEGER(LONG), INTENT(IN)       :: J_CRS(NTERM)
      REAL(DOUBLE), INTENT(IN)        :: A_CRS(NTERM)
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_FLAG
      INTEGER(LONG), INTENT(OUT)      :: INFO_OUT

      INTEGER(LONG)                   :: IROW
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: IDX
      INTEGER(LONG)                   :: NZ_KEEP
      LOGICAL                         :: KEEP_LOWER_TRI

#ifdef DMUMPS_Solver
      INFO_OUT = 0

      CALL DMUMPS_FREE_FACTORS()
      CALL DMUMPS_INIT_RUNTIME(INFO_OUT)
      IF (INFO_OUT /= 0) RETURN

      KEEP_LOWER_TRI = (SYM_FLAG(1:1) == 'Y')

      IF (KEEP_LOWER_TRI) THEN
         NZ_KEEP = 0
         DO IROW=1,N
            DO K=I_CRS(IROW),I_CRS(IROW+1)-1
               IF (J_CRS(K) <= IROW) NZ_KEEP = NZ_KEEP + 1
            ENDDO
         ENDDO
      ELSE
         NZ_KEEP = NTERM
      ENDIF

      ALLOCATE(DMUMPS_IRN(NZ_KEEP))
      ALLOCATE(DMUMPS_JCN(NZ_KEEP))
      ALLOCATE(DMUMPS_A(NZ_KEEP))

      IDX = 0
      DO IROW=1,N
         DO K=I_CRS(IROW),I_CRS(IROW+1)-1
            IF (KEEP_LOWER_TRI) THEN
               IF (J_CRS(K) > IROW) CYCLE
            ENDIF
            IDX = IDX + 1
            DMUMPS_IRN(IDX) = IROW
            DMUMPS_JCN(IDX) = J_CRS(K)
            DMUMPS_A(IDX)   = A_CRS(K)
         ENDDO
      ENDDO

      DMUMPS_PAR%JOB = -1
      DMUMPS_PAR%PAR = 1
      DMUMPS_PAR%SYM = 0
      IF (SYM_FLAG(1:1) == 'Y') DMUMPS_PAR%SYM = 2
      DMUMPS_PAR%COMM = DMUMPS_COMM
      CALL DMUMPS(DMUMPS_PAR)

      DMUMPS_PAR%ICNTL(1) = -1
      DMUMPS_PAR%ICNTL(2) = -1
      DMUMPS_PAR%ICNTL(3) = -1
      DMUMPS_PAR%ICNTL(4) = 0
      DMUMPS_PAR%N   = N
      DMUMPS_PAR%NZ  = NZ_KEEP
      DMUMPS_PAR%IRN => DMUMPS_IRN
      DMUMPS_PAR%JCN => DMUMPS_JCN
      DMUMPS_PAR%A   => DMUMPS_A
      DMUMPS_PAR%JOB = 4
      CALL DMUMPS(DMUMPS_PAR)
      INFO_OUT = DMUMPS_PAR%INFOG(1)
      IF (INFO_OUT == 0) DMUMPS_ACTIVE = .TRUE.
#else
      INFO_OUT = -999
#endif

      END SUBROUTINE DMUMPS_FACTOR_CRS

      SUBROUTINE DMUMPS_SOLVE_VECTOR ( N, RHS_COL, INFO_OUT )

      INTEGER(LONG), INTENT(IN)          :: N
      REAL(DOUBLE), TARGET, INTENT(INOUT):: RHS_COL(N)
      INTEGER(LONG), INTENT(OUT)         :: INFO_OUT

#ifdef DMUMPS_Solver
      IF (.NOT. DMUMPS_ACTIVE) THEN
         INFO_OUT = -998
         RETURN
      ENDIF

      DMUMPS_PAR%NRHS = 1
      DMUMPS_PAR%LRHS = N
      DMUMPS_PAR%RHS  => RHS_COL
      DMUMPS_PAR%JOB  = 3
      CALL DMUMPS(DMUMPS_PAR)
      INFO_OUT = DMUMPS_PAR%INFOG(1)
#else
      INFO_OUT = -999
#endif

      END SUBROUTINE DMUMPS_SOLVE_VECTOR

      SUBROUTINE DMUMPS_FREE_FACTORS()

#ifdef DMUMPS_Solver
      IF (DMUMPS_ACTIVE) THEN
         DMUMPS_PAR%JOB = -2
         CALL DMUMPS(DMUMPS_PAR)
         DMUMPS_ACTIVE = .FALSE.
      ENDIF

      IF (ALLOCATED(DMUMPS_IRN)) DEALLOCATE(DMUMPS_IRN)
      IF (ALLOCATED(DMUMPS_JCN)) DEALLOCATE(DMUMPS_JCN)
      IF (ALLOCATED(DMUMPS_A  )) DEALLOCATE(DMUMPS_A  )

      IF (DMUMPS_MPI_ACTIVE) THEN
         CALL MPI_FINALIZE(DMUMPS_MPI_IERR)
         DMUMPS_MPI_ACTIVE = .FALSE.
      ENDIF
#endif

      END SUBROUTINE DMUMPS_FREE_FACTORS

      SUBROUTINE DMUMPS_INIT_RUNTIME ( INFO_OUT )

      INTEGER(LONG), INTENT(OUT)      :: INFO_OUT

#ifdef DMUMPS_Solver
      INFO_OUT = 0
      IF (DMUMPS_MPI_ACTIVE) RETURN

      CALL MPI_INIT(DMUMPS_MPI_IERR)
      IF (DMUMPS_MPI_IERR /= 0) THEN
         INFO_OUT = DMUMPS_MPI_IERR
         RETURN
      ENDIF
      DMUMPS_MPI_ACTIVE = .TRUE.
      DMUMPS_COMM = MPI_COMM_WORLD
#else
      INFO_OUT = -999
#endif

      END SUBROUTINE DMUMPS_INIT_RUNTIME

      END MODULE DMUMPS_STUF
! --- MUMPS_COO add end --- !
