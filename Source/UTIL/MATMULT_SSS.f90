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
 
      SUBROUTINE MATMULT_SSS ( MAT_A_NAME, NROW_A, NTERM_A, SYM_A, I_A, J_A, A,                                                    &
                               MAT_B_NAME, NCOL_B, NTERM_B, SYM_B, J_B, I_B, B, AROW_MAX_TERMS, MAT_C_NAME, CONS,                  &
                                                   NTERM_C,        I_C, J_C, C )
 
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! Subroutine MATMULT_SSS_NTERM must be run before this subroutine to calculate NTERM_C, an input to this subroutine, that is the
! number of nonzero terms in C. Then memory can be allocated to C before this subroutine is called
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

! Sparse matrix multiply to obtain C = cons*A*B with A, B and C in sparse format. A and C are stored in compressed row storage
! (CRS) format. In addition, if A is symmetric it can be stored  with only the terms on, and above, the diagonal.
! Matrix B must be stored in compressed col storage (CCS) and cannot be stored symmetric (i.e. all terms in B must be stored)

! Input  matrix A is stored in compressed row storage (CRS) format using arrays I_A(NROW_A+1), J_A(NTERM_A) A(NTERM_A) where NROW_A
! is the number of rows in matrix A and NTERM_A are the number of nonzero terms in matrix A:

!      I_A is an array of NROW_A+1 integers that is used to specify the number of nonzero terms in rows of matrix A. That is:
!          I_A(I+1) - I_A(I) are the number of nonzero terms in row I of matrix A

!      J_A is an integer array giving the col numbers of the NTERM_A nonzero terms in matrix A

!        A is a real array of the nonzero terms in matrix A. If SYM_A='Y' then only the terms on, and above, the diag are stored.

! Input  matrix B is stored in compressed col storage (CCS) format using arrays J_B(NCOL_B+1), I_B(NTERM_B), B(NTERM_B) where NCOL_B
! is the number of columns in matrix B and NTERM_B are the number of nonzero terms in matrix B:

!      J_B is an array of NCOL_B+1 integers that is used to specify the number of nonzero terms in rows of matrix A. That is:
!          J_B(I+1) - J_B(I) are the number of nonzero terms in col I of matrix B

!      I_B is an integer array giving the row numbers of the NTERM_B nonzero terms in matrix B

!        B is a real array of the nonzero terms in matrix B. All of the nonzero terms of B must be included (i.e., no SYM option)

! Output matrix C, which must have the same number of rows as matrix A and the same number of columns as matrix B
! is stored in compressed row storage (CRS) format using arrays I_C(NROW_A+1), J_C(NTERM_C) C(NTERM_C) where NROW_A is
! the number of rows in matrix A and NTERM_C are the number of nonzero terms in matrix C:

!      I_C is an array of NROW_A+1 integers that is used to specify the number of nonzero terms in rows of matrix C. That is:
!          I_C(I+1) - I_C(I) are the number of nonzero terms in row I of matrix C

!      J_C is an integer array giving the col numbers of the NTERM_C nonzero terms in matrix C

!        C is a real array of the nonzero terms in matrix C. All of the nonzero terms of C are stored (i.e., no SYM option for C)

! This subr determines integer arrays I_C and J_C and real array C.

! In order to handle symmetric A matrices, which do not have all terms in a row (only those from the diagonal out), array AROW is
! used. AROW is a 1D array that contains all nonzero terms in one row of A (including those that are not explicitly in array A due
! to symmetry). The multiplication of matrix A times matrix B is then accomplished (row by row of result matrix C) by multiplying
! AROW times matrix B. AROW is a compact array containing only the nonzero terms from one row of matrix A. Thus, integer array
! J_AROW is needed to give the column numbers, from matrix A (for one row), that the terms in array AROW are for.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F04, F06, WRT_ERR, WRT_LOG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE SUBR_BEGEND_LEVELS, ONLY    :  MATMULT_SSS_BEGEND
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
 
      USE MATMULT_SSS_USE_IFs

      IMPLICIT NONE
 
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATMULT_SSS'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_A_NAME            ! Name of matrix A
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_B_NAME            ! Name of matrix B
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_C_NAME            ! Name of matrix C
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_A                 ! ='Y' if matrix A is input symmetric (terms on and above diag only)
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_B                  ! ='Y' if matrix B is input sym (terms on and above diag only)

      INTEGER(LONG), INTENT(IN )      :: AROW_MAX_TERMS        ! Max number of terms in any row of A
      INTEGER(LONG), INTENT(IN )      :: NCOL_B                ! Number of cols in input matrix B
      INTEGER(LONG), INTENT(IN )      :: NROW_A                ! Num rows in input matrix A
      INTEGER(LONG), INTENT(IN )      :: NTERM_A               ! Num non 0's in input  matrix A
      INTEGER(LONG), INTENT(IN )      :: NTERM_B               ! Num non 0's in input  matrix B
      INTEGER(LONG), INTENT(IN )      :: NTERM_C               ! Size of arrays J_C and C (MUST be determined by subr MATMULT_SSS)
      INTEGER(LONG), INTENT(IN )      :: I_A(NROW_A+1)         ! I_A(I+1) - I_A(I) = num nonzeros in row I of matrix A (CRS format)
      INTEGER(LONG), INTENT(IN )      :: J_A(NTERM_A)          ! Col no's for nonzero terms in matrix A
      INTEGER(LONG), INTENT(IN )      :: J_B(NCOL_B+1)         ! J_B(I+1) - J_B(I) = num nonzeros in col I of matrix B (CCS format)
      INTEGER(LONG), INTENT(IN )      :: I_B(NTERM_B)          ! Row no's for nonzero terms in matrix B
      INTEGER(LONG), INTENT(OUT)      :: I_C(NROW_A+1)         ! I_C(I+1) - I_C(I) = num nonzeros in row I of matrix C (CRS format)
      INTEGER(LONG), INTENT(OUT)      :: J_C(NTERM_C)          ! Col no's for nonzero terms in matrix C
      INTEGER(LONG)                   :: A_ROW_BEG             ! Index into array A where a row of matrix A begins
      INTEGER(LONG)                   :: A_ROW_END             ! Index into array A where a row of matrix A ends
      INTEGER(LONG)                   :: A_COL_NUM             ! A col number in matrix A
      INTEGER(LONG)                   :: B_COL_BEG             ! Index into array B where a col of matrix B begins
      INTEGER(LONG)                   :: B_COL_END             ! Index into array B where a col of matrix B ends
      INTEGER(LONG)                   :: B_ROW_NUM             ! A row number in matrix B
      INTEGER(LONG)                   :: A_NTERM_ROW_I         ! Number of terms in row I of matrix A
      INTEGER(LONG)                   :: B_NTERM_COL_J         ! Number of terms in col J of matrix B
      INTEGER(LONG)                   :: DELTA_KTERM_C         ! Incr in KTERM_C (0 or 1) when mult row I of A times col J of B
      INTEGER(LONG)                   :: I,J,K,L,II,memerror   ! DO loop indices
      INTEGER(LONG)                   :: I1,I2                 ! DO loop range
      INTEGER(LONG), allocatable      :: J_AROW(:)!(AROW_MAX_TERMS)! Col numbers of terms in real array AROW (see below)
      INTEGER(LONG)                   :: KTERM_C               ! Count of number of nonzero terms put into output matrix C
      INTEGER(LONG)                   :: NHITS                 ! Number of "hits" of terms in a row of A existing where terms in a
!                                                                col of B exist when a row of A is multiplied by a col of B
      INTEGER(LONG)                   :: NTERM_AROW            ! Number of nonzero terms in AROW (one row of A)
      INTEGER(LONG), allocatable      :: A_ROW_COLJ_BEG(:)!(NROW_A)! jth term is row number in array A where col j nonzeros begin 
      INTEGER(LONG), allocatable      :: A_ROW_COLJ_END(:)!(NROW_A)! jth term is row number in MATIN where col j nonzeros end
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = MATMULT_SSS_BEGEND
       
      REAL(DOUBLE) , INTENT(IN )      :: CONS                  ! Constant multiplier in cons*A*B to get C
      REAL(DOUBLE) , INTENT(IN )      :: A(NTERM_A)            ! Nonzero values in matrix A
      REAL(DOUBLE) , INTENT(IN )      :: B(NTERM_B)            ! Nonzero values in matrix B
      REAL(DOUBLE) , INTENT(OUT)      :: C(NTERM_C)            ! Nonzero values in matrix C
      REAL(DOUBLE)                    :: CTEMP                 ! A value accumulated as the nonzero terms from one row of A are
!                                                                multiplied by the corresponding nonzero terms from one col of B
      REAL(DOUBLE), allocatable       :: AROW(:)!(AROW_MAX_TERMS)  ! Array containing the nonzero terms from one row of A

      INTRINSIC                       :: MAX

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************
! Make sure B is stored as nonsym (all terms)

      IF (SYM_B /= 'N') THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,948) SUBR_NAME, MAT_B_NAME, SYM_B
         WRITE(F06,948) SUBR_NAME, MAT_B_NAME, SYM_B
         CALL OUTA_HERE ( 'Y' )
      ENDIF 

! Initialize outputs

      DO I=1,NROW_A+1
         I_C(I) = 0
      ENDDO

      DO I=1,NTERM_C
         J_C(I) = 0
           C(I) = ZERO
      ENDDO

      IF ((DEBUG(84) == 2) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_DEB ( '1', '   ' )

! Create arrays that give the row numbers at which col j begins and ends. When SYM_A = 'Y'', we have to get
! terms for the A matrix that are not explicitly in A. This is done by getting A terms in the column (above the
! diagonal of a row) as well as the explicit terms from A that are there from the diagonal out to the end of the row. The two
! arrays: A_ROW_COLJ_BEG and A_ROW_COLJ_END are used to aid in getting the terms in the column above the diagonal.
! A_ROW_COLJ_BEG is an array that gives, for each col of A, the starting row number of nonzero terms in that column.  
! A_ROW_COLJ_END is an array that gives, for each col of A, the ending   row number of nonzero terms in that column.
! The span: A_ROW_COLJ_BEG to A_ROW_COLJ_END is used when we search for terms in the columns.
! We only need A_ROW_COLJ_BEG and A_ROW_COLJ_END when A is input symmetric and MATOUT is not to be output as symmetric  
      allocate (A_ROW_COLJ_BEG(NROW_A),A_ROW_COLJ_END(NROW_A),stat=memerror)
      if (memerror.ne.0) stop 'Failed in allocating A_ROW_COLJ_BEGw(NROW_A) in MATMULT_SSS'
      if(.not.allocated(J_AROW)) allocate(J_AROW(AROW_MAX_TERMS),stat=memerror)
      if (memerror.ne.0) stop 'Fail in allocating J_AROW in MATMULT_SSS'
      if(.not.allocated(AROW)) allocate(AROw(AROW_MAX_TERMS),stat=memerror)
      if (memerror.ne.0) stop 'Fail in allocating AROW in MATMULT_SSS'
         

      IF (SYM_A == 'Y') THEN                              ! The number of cols in A is NROW_A (due to SYM_A = 'Y')

         CALL ROW_AT_COLJ_BEGEND ( MAT_A_NAME, NROW_A, NROW_A, NTERM_A, I_A, J_A, A_ROW_COLJ_BEG, A_ROW_COLJ_END )

      ENDIF

! Do the multiply, using values put into AROW and J_AROW for each row of A. This is done to facilitate the SYM option for matrix A

      DELTA_KTERM_C = 0                                    ! Initialize variables used in the matrix multiplication in the DO I loop
      KTERM_C       = 0
      NHITS         = 0
      A_ROW_BEG        = 1
      CTEMP         = ZERO
      I_C(1)        = 1
i_do: DO I=1,NROW_A                                        ! Matrix multiply loop. Range over the rows in A

         I_C(I+1) = I_C(I)                                 ! End value in I_C for next this row is initially set at beginning value
         A_NTERM_ROW_I = I_A(I+1) - I_A(I)                 ! Number of terms in matrix A in row I 
         IF (A_NTERM_ROW_I == 0) CYCLE i_do 
         A_ROW_END = A_ROW_BEG + A_NTERM_ROW_I - 1         ! A_ROW_BEG to A_ROW_END is range of indices of terms in A for row I of A
         IF ((DEBUG(84) == 2) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_DEB ( '2', '   ' )

         DO K=1,AROW_MAX_TERMS                             ! Null J_AROW and AROW each time we begin a new row of A
            AROW(K)   = ZERO
            J_AROW(K) = 0
         ENDDO 
                                                           ! Build AROW for one row of matrix A
         NTERM_AROW = 0                                    ! 1st, look for terms that would be in this row, but are not, due to SYM
         IF (SYM_A == 'Y') THEN                            
            DO K=1,A_ROW_BEG-1
               IF (J_A(K) == I) THEN
                  NTERM_AROW         = NTERM_AROW + 1
                  AROW(NTERM_AROW)   = A(K)
                  I1 = A_ROW_COLJ_BEG(I)
                  I2 = A_ROW_COLJ_END(I)
                  DO II=I1,I2
                     IF ((K >= I_A(II)) .AND. (K < I_A(II+1))) THEN
                        J_AROW(NTERM_AROW) = II
                        IF ((DEBUG(84) == 2) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_DEB ( '5', ' #2' )
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF

         DO K=A_ROW_BEG,A_ROW_END                          ! 2nd, get terms from this row of A from the diagonal out
            NTERM_AROW = NTERM_AROW + 1
            AROW(NTERM_AROW)   = A(K)
            J_AROW(NTERM_AROW) = J_A(K)
            IF ((DEBUG(84) == 2) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_DEB ( '6', ' #1' )
         ENDDO

         B_COL_BEG = 1
j_do:    DO J=1,NCOL_B                                     ! J loops over the number of columns in B
            B_NTERM_COL_J = J_B(J+1) - J_B(J)
            IF (B_NTERM_COL_J == 0) CYCLE j_do
            B_COL_END = B_COL_BEG + B_NTERM_COL_J - 1      ! B_COL_BEG to B_COL_END is range of indices of terms in B for col J of B

k_do:       DO K=1,NTERM_AROW                              ! The following 2 loops produce the ij-th term of C
               A_COL_NUM = J_AROW(K)
l_do:          DO L=B_COL_BEG,B_COL_END
                  B_ROW_NUM = I_B(L)
                  IF (A_COL_NUM == B_ROW_NUM) THEN
                     NHITS = NHITS + 1
                     DELTA_KTERM_C = 1
                     CTEMP = CTEMP + CONS*AROW(K)*B(L)     ! This is the nonzero ij-th term in matrix C. 
                  ENDIF
               ENDDO l_do
            ENDDO k_do
            B_COL_BEG   = B_COL_END + 1

            KTERM_C  = KTERM_C + DELTA_KTERM_C             ! Now update sparse CRS representation of C
            IF (KTERM_C > NTERM_C) CALL ARRAY_SIZE_ERROR_1( SUBR_NAME, NTERM_C, MAT_C_NAME ) 
             IF (NHITS > 0) THEN
               I_C(I+1) = I_C(I+1) + 1
               J_C(KTERM_C) = J
                 C(KTERM_C) = CTEMP
               IF ((DEBUG(84) == 2) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_DEB ( '7', '   ' )
            ENDIF
            DELTA_KTERM_C = 0
            CTEMP = ZERO
            NHITS = 0

         ENDDO j_do
         A_ROW_BEG = A_ROW_END + 1
         IF ((DEBUG(84) == 2) .OR. (DEBUG(84) == 3)) THEN
            WRITE(F06,*)
         ENDIF

      ENDDO i_do

      IF ((DEBUG(84) == 2) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_DEB ( '9', '   ' )

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF
      if (allocated(A_ROW_COLJ_BEG)) deallocate (A_ROW_COLJ_BEG)
      if (allocated(A_ROW_COLJ_END)) deallocate (A_ROW_COLJ_END)
      if(allocated(J_AROW)) deallocate(J_AROW)
      if(allocated(AROW)) deallocate(AROW)
      RETURN

! **********************************************************************************************************************************
  948 FORMAT(' *ERROR   948: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' MATRIX B = ',A,' MUST BE STORED AS SYM_B = ','''N''',' FOR THIS SUBR TO WORK BUT IS STORED AS',   &
                           ' SYM_B = ',''',A,''')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE MATMULT_SSS_DEB ( WHICH, ALG )

      CHARACTER(LEN=*), INTENT(IN)    :: ALG                    ! Which algorithm is used (#1 for terms above diag when SYM_A='Y'
!                                                                 or #2 for terms in row from diag out)
      CHARACTER( 1*BYTE)              :: WHICH                  ! Decides what to print out for this call to this subr
      INTEGER(LONG)                   :: I,J,K                  ! Local loop indices

! **********************************************************************************************************************************
      IF      (WHICH == '1') THEN

         WRITE(F06,*)
         WRITE(F06,1011)
         WRITE(F06,1012)
         WRITE(F06,1013)
         WRITE(F06,1014) MAT_A_NAME, MAT_B_NAME, MAT_C_NAME
         WRITE(F06,1015) CONS
         WRITE(F06,1016) NROW_A, NTERM_A, NCOL_B, NTERM_B, NTERM_C
         IF (SYM_A == 'Y') THEN
            WRITE(F06,1017)
         ELSE
            WRITE(F06,1018)
         ENDIF
         WRITE(F06,1019)
         WRITE(F06,*)

      ELSE IF (WHICH == '2') THEN

         WRITE(F06,1021)
         WRITE(F06,1022) I
         WRITE(F06,1023) I, I, A_ROW_BEG, A_ROW_END
         WRITE(F06,1024)
         WRITE(F06,*)

      ELSE IF (WHICH == '3') THEN

      ELSE IF (WHICH == '4') THEN

      ELSE IF (WHICH == '5') THEN

         WRITE(F06,1051) ALG,                     K , II, J_A(K) , A(K), NTERM_AROW, J_AROW(NTERM_AROW), AROW(NTERM_AROW)

      ELSE IF (WHICH == '6') THEN

         WRITE(F06,1061) ALG, K, I, J_A(K), A(K),                         NTERM_AROW, J_AROW(NTERM_AROW), AROW(NTERM_AROW)

      ELSE IF (WHICH == '7') THEN

         IF (NHITS > 0) THEN
            IF (J == 1) THEN
               WRITE(F06,*)
               WRITE(F06,1071) I
            ENDIF
            WRITE(F06,1072) I, J, NHITS, KTERM_C, I, J_C(KTERM_C), C(KTERM_C)
         ENDIF

      ELSE IF (WHICH == '8') THEN

      ELSE IF (WHICH == '9') THEN

         WRITE(F06,*)
         WRITE(F06,1091) MAT_C_NAME
         DO I=1,NROW_A+1                                   ! The number of rows in C is the same as that in A
            WRITE(F06,9192) I,I_C(I)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,1093)
         DO K=1,NTERM_C
            WRITE(F06,1094) K, J_C(K), C(K)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,1095)

      ENDIF

! **********************************************************************************************************************************
 1011 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' :::::::::::::::::::::::::::::::::::::::START DEBUG(84) OUTPUT FROM SUBROUTINE MATMULT_SSS:::::::::::::::::::::::::',&
              ':::::::::::::::::',/)

 1012 FORMAT(' SSS SPARSE MATRIX MULTIPLY ROUTINE: Multiply matrix A, stored in sparse Compressed Row Storage (CRS) format, times',&
' matrix B, stored in',/,' -----------------------------------',/,' sparse Compressed Column Storage (CCS) format, to obtain'     ,&
' sparse CRS matrix C',/)

 1013 FORMAT(' A may be stored as symmetric (only terms on and above the diagonal) or with all nonzero terms included.'         ,/,&
' Matrix B must be stored with all nonzero terms. Result matrix C will also be stored with all nonzero terms',/)

 1014 FORMAT(40X,' The name of CRS formatted matrix A is: ',A                                                                   ,/,&
             40x,' The name of CRS formatted matrix B is: ',A                                                                   ,/,&
             40x,' The name of CRS formatted matrix C is: ',A,/)

 1015 FORMAT(' Multiply ',1ES14.6,' times the product of matrix A and matrix B to obtain matrix C',/)

 1016 FORMAT(36X,' Matrix A has ',I8,' rows and '  ,I12,' nonzero terms'                                                        ,/,&
             36X,' Matrix B has ',I8,' cols and '  ,I12,' nonzero terms'                                                        ,/,&
             36X,' Matrix C will have             ',I12,' nonzero terms*'                                                       ,/,&
             22X,'*(as detrmined by subr MATMULT_SSS_NTERM which had to have been run prior to this subr)'/)

 1017 FORMAT(' Matrix A was input as a symmetric CRS array (only those terms on and above the diagonal are stored in array A)',/)

 1018 FORMAT(' Matrix A was input as a CRS array with all nonzero terms stored in array A',/)

 1019 FORMAT(                                                                                                                      &
' In order to handle symmetric A matrices, which do not have all terms in a row (only those from the diagonal out), arrays AROW',/,&
' and J_AROW are used. AROW is a 1D array that contains all nonzero terms from one row of A (including those that are not'      ,/,&
' explicitly in array A due to symmetry). The multiplication of matrix A times matrix B is then accomplished (row by row of'    ,/,&
' result matrix C) by multiplying AROW times matrix B. Since AROW is a compact array containing only the nonzero terms from one',/,&
' row of matrix A, integer array J_AROW is also needed to give the col numbers, from matrix A (for one row), for the terms in'  ,/,&
' array AROW.'                                                                                                                  ,//&
' Alg #1 (below) gets data for arrays J_AROW and AROW directly from the Compressed Row Storage (CRS) format of array A'         ,/,&
' Alg #2 (below) is only needed if matrix A is input as symmetric (only terms on and above the diagonal) and gets terms for'    ,/,&
'         J_AROW and AROW from column I of matrix A while working on row I of matrix A. These are the terms that would be below',/,&
'         the diag in matrix A but are not explicitly in the array due to symmetry storage'                                    ,//,&
' For each row of matrix A, the following shows the development of arrays J_AROW and AROW and the result of mult AROW times'   ,   &
' matrix B to get one row of result matrix C. Output is only given for non null rows of matrix A and non null cols of matrix B',/)

 1021 FORMAT(' ******************************************************************************************************************',&
              '*****************')

 1022 FORMAT(30X,' W O R K I N G   O N   R O W ',I8,'   O F   O U T P U T   M A T R I X   C',/)

 1023 FORMAT(' Multiply row ',I8,' of matrix A times all columns of matrix B to get row ',I8,' of matrix C'                     ,/,&
             ' This row of A begins in array A(K) at index K  = ',I8,' and ends at index K = ',I8,//)


 1024 FORMAT(16X,'Data from input array A             Data below diag of matrix A not in array A             Data for array AROW',/&
,9X,'----------------------------------------   ------------------------------------------      --------------------------------',/&
,' Alg     Index       Row     Col        Value         Index       Row     Col        Value          Index    Col No        Value'&
,/,13X,'K         I    J_A(K)       A(K)             K        II    J_A(K)       A(K)              M  J_AROW(M)      AROW(M)')

 1051 FORMAT(1X,A3,10X,10X,10X,15X    ,I10,I10,I10,1ES15.6,I11,I11,1ES16.6)

 1061 FORMAT(1X,A3,I10,I10,I10,1ES15.6,10X,10X,10X,15X    ,I11,I11,1ES16.6)

 1071 FORMAT('                                                                          Data for row',I8,' of output matrix C'  ,/,&
             '                                                                       --------------------------------------------' &
          ,/,'                                                                          Index       Row     Col        Value'   ,/,&
             '                                                                              K         I    J_C(K)       C(K)',/)

 1072 FORMAT(' Row ',I8,' of A times col ',I8,' of B gets ',I8,' hits and:   ',I10,I10,i10,1ES16.6)


 1091 FORMAT(' ******************************************************************************************************************',&
              '*****************'                                                                                               ,/,&
             ' SUMMARY: Compressed Row Storage (CRS) format of matrix C = ',A,':',/,' -------'                                 ,//,&
             ' 1) Index, L, and array I_C(L) for matrix C, where I_C(L+1) - I_C(L) is the number of nonzero terms in row L of',    &
             ' matrix C.',/,'    (also, I_C(L) is the index, K, in array C(K) where row L begins - up to, but not including, the', &
             ' last entry in I_C(L)).',/)

 9192 FORMAT('    L, I_C(L)       = ',2I12)

 1093 FORMAT(' 2) Index, K, and arrays J_C(K) and C(K). C(K) are the nonzeros in matrix C and J_C(K) is the col number in matrix',&
             ' C for term C(K).',/)

 1094 FORMAT('    K, J_C(K), C(K) = ',2I12,1ES15.6)

 1095 FORMAT(' ::::::::::::::::::::::::::::::::::::::END DEBUG(84) OUTPUT FROM SUBROUTINE MATMULT_SSS::::::::::::::::::::::::::::',&
              ':::::::::::::::::'                                                                                               ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE MATMULT_SSS_DEB

      END SUBROUTINE MATMULT_SSS
