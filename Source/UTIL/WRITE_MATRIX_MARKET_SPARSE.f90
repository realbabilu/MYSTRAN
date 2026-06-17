! ##################################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________
!
! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
! LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
! EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
! THE USE OR OTHER DEALINGS IN THE SOFTWARE.
! _______________________________________________________________________________________________________
!
! End MIT license text.

      SUBROUTINE WRITE_MATRIX_MARKET_SPARSE ( MAT_NAME, NROWS, NCOLS, SYM, NTERM_MAT, I_MAT, J_MAT, MAT )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, INFILE, LEN_INPUT_FNAME
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, WARN_ERR

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_MATRIX_MARKET_SPARSE'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_NAME
      CHARACTER(1*BYTE), INTENT(IN)   :: SYM
      INTEGER(LONG), INTENT(IN)       :: NROWS
      INTEGER(LONG), INTENT(IN)       :: NCOLS
      INTEGER(LONG), INTENT(IN)       :: NTERM_MAT
      INTEGER(LONG), INTENT(IN)       :: I_MAT(NROWS+1)
      INTEGER(LONG), INTENT(IN)       :: J_MAT(NTERM_MAT)
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: J
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: LBASE
      INTEGER(LONG)                   :: LNAME
      INTEGER(LONG)                   :: UNT
      INTEGER(LONG)                   :: IOCHK
      REAL(DOUBLE), INTENT(IN)        :: MAT(NTERM_MAT)
      CHARACTER(16*BYTE)              :: STORAGE_KIND
      CHARACTER(256*BYTE)             :: FILNAM
      CHARACTER(32*BYTE)              :: MAT_TAG

      UNT = 91
      FILNAM = ' '
      MAT_TAG = ' '
      STORAGE_KIND = 'general'

      LBASE = LEN_INPUT_FNAME - 1
      IF (LBASE < 1) THEN
         LBASE = LEN_TRIM(INFILE)
      ENDIF

      LNAME = MIN(LEN_TRIM(MAT_NAME), LEN(MAT_TAG))
      MAT_TAG(1:LNAME) = MAT_NAME(1:LNAME)

      DO I=1,LNAME
         IF (MAT_TAG(I:I) == ' ') MAT_TAG(I:I) = '_'
      ENDDO

      FILNAM(1:LBASE) = INFILE(1:LBASE)
      FILNAM(LBASE+1:LBASE+1) = '_'
      FILNAM(LBASE+2:LBASE+1+LNAME) = MAT_TAG(1:LNAME)
      FILNAM(LBASE+2+LNAME:LBASE+5+LNAME) = '.mtx'

      IF (SYM == 'Y') STORAGE_KIND = 'symmetric'

      OPEN(UNIT=UNT,FILE=FILNAM,STATUS='REPLACE',FORM='FORMATTED',ACTION='WRITE',IOSTAT=IOCHK)
      IF (IOCHK /= 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,1001) SUBR_NAME, MAT_NAME, FILNAM(1:LEN_TRIM(FILNAM))
         WRITE(F06,1001) SUBR_NAME, MAT_NAME, FILNAM(1:LEN_TRIM(FILNAM))
         RETURN
      ENDIF

      WRITE(UNT,'(A,A)') '%%MatrixMarket matrix coordinate real ', TRIM(STORAGE_KIND)
      WRITE(UNT,'(A)') '% MYSTRAN sparse matrix export'
      WRITE(UNT,'(I0,1X,I0,1X,I0)') NROWS, NCOLS, NTERM_MAT

      DO I=1,NROWS
         DO K=I_MAT(I),I_MAT(I+1)-1
            J = J_MAT(K)
            WRITE(UNT,'(I0,1X,I0,1X,ES24.16E3)') I, J, MAT(K)
         ENDDO
      ENDDO

      CLOSE(UNT)

      WRITE(ERR,1002) MAT_NAME, FILNAM(1:LEN_TRIM(FILNAM))
      WRITE(F06,1002) MAT_NAME, FILNAM(1:LEN_TRIM(FILNAM))

 1001 FORMAT(' *WARNING    : ',A,' COULD NOT OPEN MATRIX MARKET FILE FOR ',A,' : ',A)
 1002 FORMAT(' *INFORMATION: MATRIX MARKET EXPORT WRITTEN FOR ',A,' TO FILE ',A)

      END SUBROUTINE WRITE_MATRIX_MARKET_SPARSE
