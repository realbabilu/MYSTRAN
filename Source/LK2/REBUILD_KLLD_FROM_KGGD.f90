      SUBROUTINE REBUILD_KLLD_FROM_KGGD

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, SOL_NAME,                                                                   &
                                         NTERM_KNND, NTERM_KNMD, NTERM_KMMD,                                                       &
                                         NTERM_KFFD, NTERM_KFSD, NTERM_KSSD,                                                       &
                                         NTERM_KAAD, NTERM_KAOD, NTERM_KOOD,                                                       &
                                         NTERM_KLLD, NTERM_KRLD, NTERM_KRRD
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE SPARSE_MATRICES, ONLY       :  I_KNND, J_KNND, KNND, I_KNMD, J_KNMD, KNMD, I_KMMD, J_KMMD, KMMD,                         &
                                         I_KFFD, J_KFFD, KFFD, I_KFSD, J_KFSD, KFSD, I_KSSD, J_KSSD, KSSD,                         &
                                         I_KAAD, J_KAAD, KAAD, I_KAOD, J_KAOD, KAOD, I_KOOD, J_KOOD, KOOD,                         &
                                         I_KLLD, J_KLLD, KLLD, I_KRLD, J_KRLD, KRLD, I_KRRD, J_KRRD, KRRD,                         &
                                         I_GMN , J_GMN , GMN , I_GMNt, J_GMNt, GMNt, I_KMND, J_KMND, KMND,                         &
                                         I_GOA , J_GOA , GOA , I_GOAt, J_GOAt, GOAt

      USE REBUILD_KLLD_FROM_KGGD_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REBUILD_KLLD_FROM_KGGD'

      IF ((SOL_NAME(1:8) /= 'BUCKLING') .OR. (LOAD_ISTEP /= 2)) THEN
         WRITE(ERR,9101) SUBR_NAME, SOL_NAME, LOAD_ISTEP
         WRITE(F06,9101) SUBR_NAME, SOL_NAME, LOAD_ISTEP
         RETURN
      ENDIF

      IF (ALLOCATED(KNND) .OR. ALLOCATED(I_KNND) .OR. ALLOCATED(J_KNND)) CALL DEALLOCATE_SPARSE_MAT ( 'KNND' )
      IF (ALLOCATED(KNMD) .OR. ALLOCATED(I_KNMD) .OR. ALLOCATED(J_KNMD)) CALL DEALLOCATE_SPARSE_MAT ( 'KNMD' )
      IF (ALLOCATED(KMMD) .OR. ALLOCATED(I_KMMD) .OR. ALLOCATED(J_KMMD)) CALL DEALLOCATE_SPARSE_MAT ( 'KMMD' )
      IF (ALLOCATED(KFFD) .OR. ALLOCATED(I_KFFD) .OR. ALLOCATED(J_KFFD)) CALL DEALLOCATE_SPARSE_MAT ( 'KFFD' )
      IF (ALLOCATED(KFSD) .OR. ALLOCATED(I_KFSD) .OR. ALLOCATED(J_KFSD)) CALL DEALLOCATE_SPARSE_MAT ( 'KFSD' )
      IF (ALLOCATED(KSSD) .OR. ALLOCATED(I_KSSD) .OR. ALLOCATED(J_KSSD)) CALL DEALLOCATE_SPARSE_MAT ( 'KSSD' )
      IF (ALLOCATED(KAAD) .OR. ALLOCATED(I_KAAD) .OR. ALLOCATED(J_KAAD)) CALL DEALLOCATE_SPARSE_MAT ( 'KAAD' )
      IF (ALLOCATED(KAOD) .OR. ALLOCATED(I_KAOD) .OR. ALLOCATED(J_KAOD)) CALL DEALLOCATE_SPARSE_MAT ( 'KAOD' )
      IF (ALLOCATED(KOOD) .OR. ALLOCATED(I_KOOD) .OR. ALLOCATED(J_KOOD)) CALL DEALLOCATE_SPARSE_MAT ( 'KOOD' )
      IF (ALLOCATED(KLLD) .OR. ALLOCATED(I_KLLD) .OR. ALLOCATED(J_KLLD)) CALL DEALLOCATE_SPARSE_MAT ( 'KLLD' )
      IF (ALLOCATED(KRLD) .OR. ALLOCATED(I_KRLD) .OR. ALLOCATED(J_KRLD)) CALL DEALLOCATE_SPARSE_MAT ( 'KRLD' )
      IF (ALLOCATED(KRRD) .OR. ALLOCATED(I_KRRD) .OR. ALLOCATED(J_KRRD)) CALL DEALLOCATE_SPARSE_MAT ( 'KRRD' )

      IF (ALLOCATED(GMNt) .OR. ALLOCATED(I_GMNt) .OR. ALLOCATED(J_GMNt)) CALL DEALLOCATE_SPARSE_MAT ( 'GMNt' )
      IF (ALLOCATED(KMND) .OR. ALLOCATED(I_KMND) .OR. ALLOCATED(J_KMND)) CALL DEALLOCATE_SPARSE_MAT ( 'KMND' )
      IF (ALLOCATED(GOAt) .OR. ALLOCATED(I_GOAt) .OR. ALLOCATED(J_GOAt)) CALL DEALLOCATE_SPARSE_MAT ( 'GOAt' )
      IF (ALLOCATED(GMN ) .OR. ALLOCATED(I_GMN ) .OR. ALLOCATED(J_GMN )) CALL DEALLOCATE_SPARSE_MAT ( 'GMN'  )
      IF (ALLOCATED(GOA ) .OR. ALLOCATED(I_GOA ) .OR. ALLOCATED(J_GOA )) CALL DEALLOCATE_SPARSE_MAT ( 'GOA'  )

      CALL BUILD_KGGD_FROM_UG

      NTERM_KNND = 0
      NTERM_KNMD = 0
      NTERM_KMMD = 0
      NTERM_KFFD = 0
      NTERM_KFSD = 0
      NTERM_KSSD = 0
      NTERM_KAAD = 0
      NTERM_KAOD = 0
      NTERM_KOOD = 0
      NTERM_KLLD = 0
      NTERM_KRLD = 0
      NTERM_KRRD = 0

      CALL REDUCE_G_NM
      CALL REDUCE_N_FS
      CALL REDUCE_F_AO
      CALL REDUCE_A_LR

      RETURN

 9101 FORMAT(' *ERROR  9101: ',A,' was called with SOL_NAME = "',A,'" and LOAD_ISTEP = ',I0,                                       &
             '. This routine is only valid for BUCKLING step 2; ignoring call.')

      END SUBROUTINE REBUILD_KLLD_FROM_KGGD
