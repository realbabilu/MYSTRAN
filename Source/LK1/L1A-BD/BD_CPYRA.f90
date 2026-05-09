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

      SUBROUTINE BD_CPYRA ( CARD, LARGE_FLD_INP, NUM_GRD )

! --- solids_add begin --- !
! Processes CPYRA Bulk Data Cards for 5-node linear and 14-node Liu quadratic
! pyramid solid elements.
! --- solids_add end --- !

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, JCARD_LEN, NCPYRA5, NCPYRA14, NELE
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  ETYPE

      USE BD_CPYRA_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CPYRA'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP
      CHARACTER(LEN(CARD))            :: CHILD
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)
      CHARACTER(LEN(JCARD))           :: ID
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)
      CHARACTER(LEN(JCARD))           :: NAME

      INTEGER(LONG), INTENT(OUT)      :: NUM_GRD
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: ICONT = 0
      INTEGER(LONG)                   :: IERR  = 0

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      NAME = JCARD(1)
      ID   = JCARD(2)

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

! Parent: EID, PID, G1-G5. Field 9 must remain blank for the 5-node parent.
      CALL ELEPRO ( 'Y', JCARD_EDAT, 7, 7, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'N' )
      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,6,7,8,0 )
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )
      CALL CRDERR ( CARD )

      ETYPE(NELE)(1:) = ' '
      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      IF (ICONT == 0) THEN
         NCPYRA5    = NCPYRA5 + 1
         ETYPE(NELE) = 'PYRA5   '
         NUM_GRD     = 5
      ELSE
         IF (CARD(1:) /= ' ') THEN
            NCPYRA14   = NCPYRA14 + 1
            ETYPE(NELE) = 'PYRA14  '
            NUM_GRD     = 14

            DO I=1,10
               JCARD_EDAT(I) = JCARD(I)
            ENDDO
            CALL ELEPRO ( 'N', JCARD_EDAT, 8, 8, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y' )
            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )
            CALL CRDERR ( CARD )

            IF (LARGE_FLD_INP == 'N') THEN
               CALL NEXTC  ( CARD, ICONT, IERR )
            ELSE
               CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
               CARD = CHILD
            ENDIF
            CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
            IF (ICONT == 1) THEN
               DO I=1,10
                  JCARD_EDAT(I) = JCARD(I)
               ENDDO
               CALL ELEPRO ( 'N', JCARD_EDAT, 1, 1, 'Y', 'N', 'N', 'N', 'N', 'N', 'N', 'N' )
               CALL BD_IMBEDDED_BLANK( JCARD,2,0,0,0,0,0,0,0 )
               CALL CARD_FLDS_NOT_BLANK(JCARD,0,3,4,5,6,7,8,9)
               CALL CRDERR ( CARD )
            ELSE
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1136) NAME, ID
               WRITE(F06,1136) NAME, ID
            ENDIF
         ELSE
            NCPYRA5    = NCPYRA5 + 1
            ETYPE(NELE) = 'PYRA5   '
            NUM_GRD     = 5
         ENDIF
      ENDIF

      RETURN

 1136 FORMAT(' *ERROR  1136: REQUIRED CONTINUATION FOR ',A,' ID = ',A,' MISSING')

      END SUBROUTINE BD_CPYRA
