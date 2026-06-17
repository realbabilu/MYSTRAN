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
! LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
! _______________________________________________________________________________________________________
!
! End MIT license text.

      SUBROUTINE WRITE_CBEAM_STRESS (NUM, WRITE_F06)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  BARTOR, BLNK_SUB_NAM, MOGEL
      USE LINK9_STUFF, ONLY           :  CBEAM_XL_OUT, EID_OUT_ARRAY, MAXREQ, MSPRNT, OGEL

      USE WRITE_CBEAM_STRESS_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_CBEAM_STRESS'

      CHARACTER(170*BYTE)             :: BLINE1A
      CHARACTER(170*BYTE)             :: BLINE1B
      CHARACTER(170*BYTE)             :: BLINE2A
      CHARACTER(170*BYTE)             :: BLINE2B
      CHARACTER(130*BYTE)             :: BOUT1
      CHARACTER(130*BYTE)             :: BOUT2
      CHARACTER( 14*BYTE)             :: BSTA
      CHARACTER( 10*BYTE)             :: BMS1
      CHARACTER( 10*BYTE)             :: BMS2
      CHARACTER( 14*BYTE)             :: BMS3
      CHARACTER(  1*BYTE)             :: BMSF1
      CHARACTER(  1*BYTE)             :: BMSF2
      CHARACTER(  1*BYTE)             :: BMSF3
      CHARACTER( 14*BYTE)             :: BTOR
      CHARACTER(14*BYTE)              :: OGEL_CHAR(MOGEL)
      CHARACTER(  1*BYTE)             :: MSFLAG

      INTEGER(LONG), INTENT(IN)       :: NUM
      LOGICAL,       INTENT(IN)       :: WRITE_F06
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: J
      INTEGER(LONG)                   :: K

! --- cbeam_stations begin --- !
      IF (.NOT. WRITE_F06) THEN
         RETURN
      ENDIF

      K = 0
      DO I=1,NUM
         BLINE1A(1:) = ' '
         BLINE1B(1:) = ' '
         BLINE2A(1:) = ' '
         BLINE2B(1:) = ' '

         BOUT1(1:)   = ' '
         BOUT2(1:)   = ' '
         BSTA(1:)    = ' '
         BMS1(1:)    = ' '
         BMS2(1:)    = ' '
         BMS3(1:)    = ' '
         BMSF1       = ' '
         BMSF2       = ' '
         BMSF3       = ' '
         BTOR(1:)    = ' '

         K = K + 1
         CALL WRT_REAL_TO_CHAR_VAR ( OGEL, MAXREQ, MOGEL, K, OGEL_CHAR )
         WRITE(BSTA,'(F8.4)') CBEAM_XL_OUT(I)
         WRITE(BOUT1,9011) EID_OUT_ARRAY(I,1), BSTA, (OGEL_CHAR(J),J=1,7)

         MSFLAG = ' '
         IF (MSPRNT(K,1) /= '0') THEN
            IF (OGEL(K,8) < 0.0D0) MSFLAG = '*'
            WRITE(BMS1, 9021) OGEL(K,8)
            WRITE(BMSF1,9031) MSFLAG
         ELSE
            WRITE(BMS1, 9022)
            WRITE(BMSF1,9032)
         ENDIF
         IF (BARTOR == 'Y') THEN
            WRITE(BTOR, 9041) OGEL(K,9)
         ENDIF

         K = K + 1
         CALL WRT_REAL_TO_CHAR_VAR ( OGEL, MAXREQ, MOGEL, K, OGEL_CHAR )
         WRITE(BOUT2,9012) BSTA, (OGEL_CHAR(J),J=1,4), (OGEL_CHAR(J),J=6,7)

         MSFLAG = ' '
         IF (MSPRNT(K,2) /= '0') THEN
            IF (OGEL(K,8) < 0.0D0) MSFLAG = '*'
            WRITE(BMS2, 9021) OGEL(K,8)
            WRITE(BMSF2,9031) MSFLAG
         ELSE
            WRITE(BMS2, 9022)
            WRITE(BMSF2,9032)
         ENDIF

         MSFLAG = ' '
         IF (MSPRNT(K,3) /= '0') THEN
            IF (OGEL(K,9) < 0.0D0) MSFLAG = '*'
            WRITE(BMS3, 9023) OGEL(K,9)
            WRITE(BMSF3,9031) MSFLAG
         ELSE
            WRITE(BMS3, 9024)
            WRITE(BMSF3,9032)
         ENDIF

         WRITE(F06,*)
         IF (BARTOR == 'Y') THEN
            BLINE1A = BOUT1//BMS1//BMSF1//BTOR
            BLINE2A = BOUT2//BMS2//BMSF2//BMS3//BMSF3
            WRITE(F06,9031) BLINE1A
            WRITE(F06,9031) BLINE2A
         ELSE
            BLINE1B = BOUT1//BMS1//BMSF1
            BLINE2B = BOUT2//BMS2//BMSF2
            WRITE(F06,9031) BLINE1B
            WRITE(F06,9031) BLINE2B
         ENDIF
      ENDDO
! --- cbeam_stations end --- !

      RETURN

 9011 FORMAT(1X,I8,1X,A,7A)
 9012 FORMAT(1X,8X,1X,A,4A,14X,2A)
 9021 FORMAT(ES10.2)
 9022 FORMAT('          ')
 9023 FORMAT(4X,ES10.2)
 9024 FORMAT('              ')
 9031 FORMAT(A)
 9032 FORMAT(' ')
 9041 FORMAT(ES14.6)

      END SUBROUTINE WRITE_CBEAM_STRESS
