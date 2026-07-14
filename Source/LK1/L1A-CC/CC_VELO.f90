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

      SUBROUTINE CC_VELO ( CARD )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  PCHSTAT
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_CMD_DESCRIBERS, LSUB, NSUB, NCCCD
      USE TIMDAT, ONLY                :  TSEC
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  VELO_OUT
      USE MODEL_STUF, ONLY            :  SC_VELO

      USE CC_VELO_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_VELO'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD
      CHARACTER( 1*BYTE)              :: FOUND_PRINT
      CHARACTER( 1*BYTE)              :: FOUND_PLOT
      CHARACTER( 1*BYTE)              :: FOUND_PUNCH
      CHARACTER( 1*BYTE)              :: FOUND_NEU
      CHARACTER( 1*BYTE)              :: FOUND_CSV

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: SETID

      CALL CC_OUTPUTS ( CARD, 'VELO', SETID )

      FOUND_PRINT = 'N'
      FOUND_PLOT  = 'N'
      FOUND_PUNCH = 'N'
      FOUND_NEU   = 'N'
      FOUND_CSV   = 'N'
      DO I=1,NCCCD
         IF (CC_CMD_DESCRIBERS(I)(1:5) == 'PRINT') FOUND_PRINT = 'Y'
         IF ((CC_CMD_DESCRIBERS(I)(1:4) == 'PLOT') .OR. (CC_CMD_DESCRIBERS(I)(1:4) == 'POST')) FOUND_PLOT = 'Y'
         IF (CC_CMD_DESCRIBERS(I)(1:5) == 'PUNCH') FOUND_PUNCH = 'Y'
         IF (CC_CMD_DESCRIBERS(I)(1:3) == 'NEU')   FOUND_NEU   = 'Y'
         IF (CC_CMD_DESCRIBERS(I)(1:3) == 'CSV')   FOUND_CSV   = 'Y'
      ENDDO

      VELO_OUT = TRIM(FOUND_PRINT) // TRIM(FOUND_PLOT) // TRIM(FOUND_PUNCH) // TRIM(FOUND_NEU) // TRIM(FOUND_CSV)
      IF (VELO_OUT(1:5) == 'NNNNN') THEN
         VELO_OUT = 'YYNNN'
      ENDIF
      IF (FOUND_PUNCH == 'Y') PCHSTAT = 'KEEP    '

      IF (NSUB == 0) THEN
         DO I = 1,LSUB
            SC_VELO(I) = SETID
         ENDDO
      ELSE
         SC_VELO(NSUB) = SETID
      ENDIF

      END SUBROUTINE CC_VELO
