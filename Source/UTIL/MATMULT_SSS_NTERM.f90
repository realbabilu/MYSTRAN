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
 
      SUBROUTINE MATMULT_SSS_NTERM ( MAT_A_NAME, NROW_A, NTERM_A, SYM_A, I_A, J_A,                                                 &
                                     MAT_B_NAME, NCOL_B, NTERM_B, SYM_B, J_B, I_B, AROW_MAX_TERMS, MAT_C_NAME, NTERM_C )

! Setup routine for subr MATMULT_SSS (which performs the matrix operation C = cons*A*B With A and B in sparse format
! Matrices A and C are stored in compressed row storage (CRS) format. In addition, if A is symmetric it can be stored  with only the
! terms on, and above, the diagonal. Matrix B must be stored in compressed col storage (CCS) format.

! Input  matrix A is stored in compressed row storage (CRS) format using arrays I_A(NROW_A+1), J_A(NTERM_A) A(NTERM_A) where NROW_A
! is the number of rows in matrix A and NTERM_A are the number of nonzero terms in matrix A:

!      I_A is an array of NROW_A+1 integers that is used to specify the number of nonzero terms in rows of matrix A. That is:
!          I_A(I+1) - I_A(I) are the number of nonzero terms in row I of matrix A

!      J_A is an integer array giving the col numbers of the NTERM_A nonzero terms in matrix A

!        A is a real array of the nonzero terms in matrix A. If SYM_A='Y' then only the terms on, and above, the diag are stored.
!          array A is not used in this subr - it is used in the actual matrix multiplication routine.

! Input  matrix B is stored in compressed col storage (CCS) format using arrays J_B(NCOL_B+1), I_B(NTERM_B), B(NTERM_B) where NCOL_B
! is the number of columns in matrix B and NTERM_B are the number of nonzero terms in matrix B:

!      J_B is an array of NCOL_B+1 integers that is used to specify the number of nonzero terms in rows of matrix A. That is:
!          J_B(I+1) - J_B(I) are the number of nonzero terms in col I of matrix B

!      I_B is an integer array giving the row numbers of the NTERM_B nonzero terms in matrix B

!        B is a real array of the nonzero terms in matrix B. All of the nonzero terms of B must be included (i.e., no SYM option)
!          array B is not used in this subr - it is used in the actual matrix multiplication routine.

! This subr determines the storage required (NTERM_C) for array J_C (and C) from the data in arrays I_A, J_A, I_B and J_B. This
! information is used to allocate memory for arrays J_C and C prior to calling subr MATMULT_SSS so that it can it can do the sparse
! matrix multiply

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, WRT_LOG, ERR, F04, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE SUBR_BEGEND_LEVELS, ONLY    :  MATMULT_SSS_NTERM_BEGEND
      USE SPARSE_ALG_ARRAYS, ONLY     :  J_AROW
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
 
      USE MATMULT_SSS_NTERM_USE_IFs

      IMPLICIT NONE
 
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATMULT_SSS_NTERM'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_A_NAME             ! Name of matrix A
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_B_NAME             ! Name of matrix B
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_C_NAME             ! Name of matrix C
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_A                  ! ='Y' if matrix A is input sym (terms on and above diag only)
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_B                  ! ='Y' if matrix B is input sym (terms on and above diag only)

      INTEGER(LONG), INTENT(IN )      :: NCOL_B                 ! Number of cols in input matrix B
      INTEGER(LONG), INTENT(IN )      :: NROW_A                 ! Number of rows in input matrix A
      INTEGER(LONG), INTENT(IN )      :: NTERM_A                ! Number of nonzero terms in input matrix A
      INTEGER(LONG), INTENT(IN )      :: NTERM_B                ! Number of nonzero terms in input matrix B
      INTEGER(LONG), INTENT(IN )      :: I_A(NROW_A+1)          ! I_A(I+1) - I_A(I) = no. terms in row I of matrix A
      INTEGER(LONG), INTENT(IN )      :: J_A(NTERM_A)           ! Col no's for nonzero terms in matrix A
      INTEGER(LONG), INTENT(IN )      :: J_B(NCOL_B+1)          ! J_B(I+1) - J_B(I) = no. terms in row I of matrix B
      INTEGER(LONG), INTENT(IN )      :: I_B(NTERM_B)           ! Row no's for nonzero terms in matrix B
      INTEGER(LONG), INTENT(OUT)      :: AROW_MAX_TERMS         ! Max number of terms in any row of A
      INTEGER(LONG), INTENT(OUT)      :: NTERM_C                ! Number of nonzero terms in output matrix C
      INTEGER(LONG), allocatable      :: A_ROW_COLJ_BEG(:)!(NROW_A) ! jth term is row number in A where col j nonzeros begin 
      INTEGER(LONG), allocatable      :: A_ROW_COLJ_END(:)!(NROW_A) ! jth term is row number in A where col j nonzeros end
      INTEGER(LONG)                   :: A_NTERM_ROW_I          ! Number of terms in row I of matrix A
      INTEGER(LONG)                   :: A_COL_NUM              ! A col number in matrix A
      INTEGER(LONG)                   :: A_ROW_BEG              ! Index into array I_A where a row of matrix A begins
      INTEGER(LONG)                   :: A_ROW_END              ! Index into array I_A where a row of matrix A ends
      INTEGER(LONG)                   :: B_NTERM_COL_J          ! Number of terms in col J of matrix B
      INTEGER(LONG)                   :: B_ROW_NUM              ! A row number in matrix B
      INTEGER(LONG)                   :: B_COL_BEG              ! Index into array J_B where a col of matrix B begins
      INTEGER(LONG)                   :: B_COL_END              ! Index into array J_B where a col of matrix B ends
      INTEGER(LONG)                   :: DELTA_NTERM_C = 0      ! Incr in NTERM_C (0,1) resulting from mult row of A times col of B
      INTEGER(LONG)                   :: I,J,K,L,II,memerror    ! DO loop indices
      INTEGER(LONG)                   :: I1,I2                  ! DO loop range
      INTEGER(LONG)                   :: NHITS                  ! Num of "hits" of terms in a col of A existing where terms in a
!                                                                 row of B exist when a row of A is multiplied by a col of B
      INTEGER(LONG)                   :: NHITS_TOT_FOR_ROW_OF_A ! Num of "hits" of terms in a col of A existing where terms in any
!                                                                 row of B exist when a row of A is multiplied by al cols of B
      INTEGER(LONG)                   :: NTERM_AROW             ! Max number of nonzero terms in one row of A
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = MATMULT_SSS_NTERM_BEGEND
       
      INTRINSIC                       :: MAX

      allocate(A_ROW_COLJ_END(NROW_A),A_ROW_COLJ_BEG(NROW_A),stat=memerror)
      if (memerror.ne.0) STOP 'A_ROW_COLJ_END in MATMULT_SSS_NTERM cant allocated'
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

      AROW_MAX_TERMS = 0
      NTERM_C        = 0

! First, set up array J_AROW so it can handle the row from A that has the most terms in it

      A_ROW_BEG         = 1
      DO I=1,NROW_A
         A_NTERM_ROW_I = I_A(I+1) - I_A(I)
         A_ROW_END = A_ROW_BEG + A_NTERM_ROW_I - 1
         NTERM_AROW = 0
         IF (SYM_A == 'Y') THEN
            DO K=1,A_ROW_BEG-1
               IF (J_A(K) == I) THEN
                  NTERM_AROW = NTERM_AROW + 1
               ENDIF
            ENDDO
         ENDIF
         DO K=A_ROW_BEG,A_ROW_END
            NTERM_AROW = NTERM_AROW + 1
         ENDDO
         IF (NTERM_AROW > AROW_MAX_TERMS) THEN
            AROW_MAX_TERMS = NTERM_AROW
         ENDIF
         A_ROW_BEG = A_ROW_END + 1
      ENDDO
      NTERM_AROW = 0
                                                           ! Allocate integer vector to hold the column no's of the terms in the row
      CALL ALLOCATE_SPARSE_ALG ( 'J_AROW', AROW_MAX_TERMS, 0, SUBR_NAME)

      IF ((DEBUG(84) == 1) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_NTERM_DEB ( '1', '   ' )

! Create arrays that give the row numbers at which col j begins and ends. When SYM_A = 'Y'', we have to get
! terms for the A matrix that are not explicitly in A. This is done by getting A terms in the column (above the
! diagonal of a row) as well as the explicit terms from A that are there from the diagonal out to the end of the row. The two
! arrays: A_ROW_COLJ_BEG and A_ROW_COLJ_END are used to aid in getting the terms in the column above the diagonal.
! A_ROW_COLJ_BEG is an array that gives, for each col of A, the starting row number of nonzero terms in that column.  
! A_ROW_COLJ_END is an array that gives, for each col of A, the ending   row number of nonzero terms in that column.
! The span: A_ROW_COLJ_BEG to A_ROW_COLJ_END is used when we search for terms in the columns.
! We only need A_ROW_COLJ_BEG and A_ROW_COLJ_END when A is input symmetric and MATOUT is not to be output as symmetric  

      IF (SYM_A == 'Y') THEN

         CALL ROW_AT_COLJ_BEGEND ( MAT_A_NAME, NROW_A, NROW_A, NTERM_A, I_A, J_A, A_ROW_COLJ_BEG, A_ROW_COLJ_END )

      ENDIF

! Now count the terms that go into C, using values put into AROW and J_AROW for each row of A. This is done to facilitate the
! SYM option for matrix A.

      NHITS   = 0                                          ! Initialize variables used in the matrix multiplication in the DO I loop
      NHITS_TOT_FOR_ROW_OF_A = 0
      A_ROW_BEG  = 1      
i_do: DO I=1,NROW_A                                        ! Matrix multiply loop. Range over the rows in A

         A_NTERM_ROW_I = I_A(I+1) - I_A(I)                 ! Number of terms in matrix A in row I
         IF (A_NTERM_ROW_I == 0) CYCLE i_do 
         A_ROW_END = A_ROW_BEG + A_NTERM_ROW_I - 1         ! A_ROW_BEG to A_ROW_END is the range of indices of terms in A for row I
         IF ((DEBUG(84) == 1) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_NTERM_DEB ( '2', '   ' )

         DO K=1,AROW_MAX_TERMS                             ! Null J_AROW and AROW each time we begin a new row of A
            J_AROW(K) = 0
         ENDDO 
                                                           ! Formulate J_AROW
         NTERM_AROW = 0                                    ! 1st, look for terms that would be in this row, but are not, due to SYM
         IF (SYM_A == 'Y') THEN
            DO K=1,A_ROW_BEG-1
               IF (J_A(K) == I) THEN
                  NTERM_AROW = NTERM_AROW + 1
                  I1 = A_ROW_COLJ_BEG(I)
                  I2 = A_ROW_COLJ_END(I)
                  DO II=I1,I2
                     IF ((K >= I_A(II)) .AND. (K < I_A(II+1))) THEN
                        J_AROW(NTERM_AROW) = II
                        IF ((DEBUG(84) == 1) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_NTERM_DEB ( '5', ' #2' )
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
         DO K=A_ROW_BEG,A_ROW_END                          ! 2nd, get terms from this row of A from the diagonal out
            NTERM_AROW = NTERM_AROW + 1
            J_AROW(NTERM_AROW) = J_A(K)
            IF ((DEBUG(84) == 1) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_NTERM_DEB ( '6', ' #1' )
         ENDDO

         B_COL_BEG = 1
j_do:    DO J=1,NCOL_B                                     ! J loops over the number of columns in B
            B_NTERM_COL_J = J_B(J+1) - J_B(J)
            IF (B_NTERM_COL_J == 0) CYCLE j_do
            B_COL_END = B_COL_BEG + B_NTERM_COL_J - 1

k_do:       DO K=1,NTERM_AROW                              ! The following 2 loops det  the ij-th term of C 
               A_COL_NUM = J_AROW(K)
l_do:          DO L=B_COL_BEG,B_COL_END
                  B_ROW_NUM = I_B(L)
                  IF (A_COL_NUM == B_ROW_NUM) THEN
                     NHITS = NHITS + 1
                     NHITS_TOT_FOR_ROW_OF_A = NHITS_TOT_FOR_ROW_OF_A + 1
                     DELTA_NTERM_C = 1
                  ENDIF
               ENDDO l_do
            ENDDO k_do
            B_COL_BEG = B_COL_END + 1
 
            NTERM_C = NTERM_C + DELTA_NTERM_C
            IF ((DEBUG(84) == 1) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_NTERM_DEB ( '7', '   ' )
            DELTA_NTERM_C = 0
            NHITS = 0

         ENDDO j_do

         IF ((DEBUG(84) == 1) .OR. (DEBUG(84) == 3)) THEN
            WRITE(F06,*)
         ENDIF

         IF (NHITS_TOT_FOR_ROW_OF_A == 0) THEN
             IF ((DEBUG(83) == 1) .OR. (DEBUG(83) == 3)) CALL MATMULT_SSS_NTERM_DEB ( '8', '   ' )
         ENDIF
         NHITS_TOT_FOR_ROW_OF_A = 0
         A_ROW_BEG = A_ROW_END + 1

      ENDDO i_do

      IF ((DEBUG(84) == 1) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_NTERM_DEB ( '9', '   ' )
            
      CALL DEALLOCATE_SPARSE_ALG ( 'J_AROW' )

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF
      if (allocated(A_ROW_COLJ_END)) deallocate(A_ROW_COLJ_END)
      if (allocated(A_ROW_COLJ_BEG)) deallocate(A_ROW_COLJ_BEG) 
      RETURN

! **********************************************************************************************************************************
  948 FORMAT(' *ERROR   948: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' MATRIX B WHICH IS "',A,'" MUST BE STORED AS SYM_B = "N" FOR THIS SUBR TO WORK BUT IS STORED AS',      &
                           ' SYM_B = "',A,'"')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE MATMULT_SSS_NTERM_DEB ( WHICH, ALG )

      CHARACTER( 3*BYTE), INTENT(IN)  :: ALG                    ! Which algorithm is used (#1 for terms above diag when SYM_A='Y'
!                                                                 or #2 for terms in row from diag out)
      CHARACTER( 1*BYTE)              :: WHICH                  ! Decides what to print out for this call to this subr

! **********************************************************************************************************************************
      IF      (WHICH == '1') THEN

         WRITE(F06,*)
         WRITE(F06,2011)
         WRITE(F06,2012)
         WRITE(F06,2013)
         WRITE(F06,2014) MAT_A_NAME, MAT_B_NAME, MAT_C_NAME
         WRITE(F06,2016) NROW_A, NTERM_A, NCOL_B, NTERM_B
         IF (SYM_A == 'Y') THEN
            WRITE(F06,2017)
         ELSE
            WRITE(F06,2018)
         ENDIF
         WRITE(F06,2019)
         WRITE(F06,*)

      ELSE IF (WHICH == '2') THEN

         WRITE(F06,2021)
         WRITE(F06,2022) I
         WRITE(F06,2023) I, I, A_ROW_BEG, A_ROW_END
         WRITE(F06,2024)
         WRITE(F06,*)

      ELSE IF (WHICH == '3') THEN

      ELSE IF (WHICH == '4') THEN

      ELSE IF (WHICH == '5') THEN

         WRITE(F06,2051) ALG,               K , II, J_A(K) , NTERM_AROW, J_AROW(NTERM_AROW)

      ELSE IF (WHICH == '6') THEN

         WRITE(F06,2061) ALG, K, I, J_A(K),                  NTERM_AROW, J_AROW(NTERM_AROW) 

      ELSE IF (WHICH == '7') THEN

         IF (NHITS > 0) THEN
            IF (J == 1) THEN
               WRITE(F06,*)
            ENDIF
            WRITE(F06,2071) I, J, NHITS, NTERM_C
         ENDIF

      ELSE IF (WHICH == '8') THEN

         WRITE(F06,2081) I
         WRITE(F06,*)

      ELSE IF (WHICH == '9') THEN

         WRITE(F06,*)
         WRITE(F06,2021)
         WRITE(F06,2091) MAT_C_NAME, NTERM_C
         WRITE(F06,9100)
         WRITE(F06,*)

      ENDIF

! **********************************************************************************************************************************
 2011 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' ::::::::::::::::::::::::::::::::::::START DEBUG(84) OUTPUT FROM SUBROUTINE MATMULT_SSS_NTERM::::::::::::::::::::::',&
              ':::::::::::::::::',/)

 2012 FORMAT(' SSS SETUP FOR SPARSE MATRIX MULTIPLY ROUTINE: Determine memory required for sparse Compressed Row Storage (CRS)'   ,&
' formatted matrix C',/,' --------------------------------------------',/,' resulting from the multiplication of sparse CRS'      ,&
' matrix A and sparse Compressed Column Storage (CCS) formatted matrix B',/)

 2013 FORMAT(' A may be stored as symmetric (only terms on and above the diagonal) or with all nonzero terms included.'         ,/,&
' Matrix B must be stored with all nonzero terms.',/)

 2014 FORMAT(40X,' The name of CRS formatted matrix A is: ',A                                                                   ,/,&
             40x,' The name of CRS formatted matrix B is: ',A                                                                   ,/,&
             40x,' The name of CRS formatted matrix C is: ',A,/)

 2016 FORMAT(36X,' Matrix A has ',I8,' rows and '  ,I12,' nonzero terms'                                                        ,/,&
             36X,' Matrix B has ',I8,' cols and '  ,I12,' nonzero terms'                                                        ,/,&
             36X,' Matrix C will have                       20 nonzero terms*'                                                  ,/,&
             22X,' *(as detrmined by subr MATMULT_SSS_NTERM which had to have been run prior to this subr',/)

 2017 FORMAT(' Matrix A was flagged as being a symmetric CRS array (only nonzero terms on and above the diagonal)',/)

 2018 FORMAT(' Matrix A was flagged as a CRS array that contains all nonzero terms',/)

 2019 FORMAT(                                                                                                                      &
' In order to handle symmetric A matrices, which do not have all terms in a row (only those from the diag out), array J_AROW is',/,&
' used. J_AROW is a 1D array that contains the column numbers for all nonzero terms in one row of matrix A including those that',/,&
' may not actually be in array A due to symmetry storage. This array is used in a simulated multiply of one row of matrix A',/    ,&
' times columns of matrix B to determine whether there are any "hits" of terms in the row of A with terms in a column of B. In' ,/,&
' this manner, the number of nonzero terms (NTERM_C) that will be in array C can be estimated. The estimate is exact except for',/,&
' the case where a row of A multiplies a column of B and has "hits" that turn out to sum to zero when subr MATMULT_SSS is run' ,//,&
' Alg #1 (below) gets data for array J_AROW directly from the Compressed Row Storage (CRS) format of array A'                   ,/,&
' Alg #2 (below) is only needed if matrix A is input as symmetric (only terms on and above the diagonal) and gets terms for'    ,/,&
'         J_AROW from column I of matrix A while working on row I of matrix A. These are the terms that would be below the diag',/,&
'         in matrix A but are not explicitly in the array due to symmetry storage',//,                                             &
' For each row of matrix A, the following shows the development of array AROW and the result of multiplying AROW times matrix B',/,&
' to get one row of result matrix C. Output is only given for non null rows of matrix A and non null cols of matrix B',/)

 2021 FORMAT(' ******************************************************************************************************************',&
              '*****************')

 2022 FORMAT(30X,' W O R K I N G   O N   R O W ',I8,'   O F   O U T P U T   M A T R I X   C',/)

 2023 FORMAT(' Simulate the multiplication of row ',I8,' of matrix A times all cols of matrix B to get row ',I8,' of matrix C'  ,/,&
             ' This row of A begins in array A(K) at index K  = ',I8,' and ends at index K = ',I8,//)


 2024 FORMAT(11X,'Data from input array A           Data below diag of matrix A not in array A             Data for array AROW' ,/,&
 9X,'--------------------------          ------------------------------------------             -------------------'            ,/,&
 ' Alg     Index       Row     Col                      Index       Row     Col                         Index    Col No'           &
,/,13X,'K         I    J_A(K)                        K        II    J_A(K)                           M  J_AROW(M)')

 2051 FORMAT(1X,A3,10X,10X,10X,I25,I10,I10,I28,I11)

 2061 FORMAT(1X,A3,I10,I10,I10,25X,10X,10X,I28,I11)

 2071 FORMAT(' Row ',I8,' of A times col ',I8,' of B results in ',I8,' hits and accum number of nonzero terms in C = ',I9)

 2081 FORMAT(' There were no "hits" of any terms in row ',I8,' of matrix A multiplying terms from any col in matrix B')

 2091 FORMAT(' Total number of nonzero terms needed for allocating arrays J_C and C for matrix C = ',A,' is NTERM_C = ',I8,//)

 9100 FORMAT(' :::::::::::::::::::::::::::::::::::::END DEBUG(84) OUTPUT FROM SUBROUTINE MATMULT_SSS_NTERM:::::::::::::::::::::::',&
              ':::::::::::::::::'                                                                                               ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE MATMULT_SSS_NTERM_DEB

      END SUBROUTINE MATMULT_SSS_NTERM
