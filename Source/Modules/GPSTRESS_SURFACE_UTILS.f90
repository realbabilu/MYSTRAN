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
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________
!
! End MIT license text.

      MODULE GPSTRESS_SURFACE_UTILS

      USE PENTIUM_II_KIND, ONLY         :  BYTE, LONG
      USE IOUNT1, ONLY                  :  ERR, F06
      USE SCONTR, ONLY                  :  BLNK_SUB_NAM, CC_ENTRY_LEN, FATAL_ERR, LSETLN, MAX_TOKEN_LEN, SETLEN
      USE MODEL_STUF, ONLY              :  ALL_SETS_ARRAY, EPNT, ESORT1, ETYPE, EDAT
      USE CC_OUTPUT_DESCRIBERS, ONLY    :  GP_SURFACE_IDS, GP_SURFACE_NORMAL_MODE, GP_SURFACE_SETIDS, NUM_GP_SURFACE, &
                                           GP_POST_SET_IDS, GP_POST_SET_TEXT, NUM_GP_POST_SET, REGISTER_GP_SURFACE,   &
                                           REGISTER_GP_POST_SET

      USE GET_ARRAY_ROW_NUM_Interface
      USE STOKEN_Interface
      USE TOKCHK_Interface

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE PARSE_GPSTRESS_SET_CARD ( CARD )

      CHARACTER(LEN=*), INTENT(IN)    :: CARD

      CHARACTER(LEN=CC_ENTRY_LEN)     :: CARD1
      CHARACTER(LEN=LSETLN)           :: SET_TEXT
      CHARACTER(16*BYTE)              :: KEY1
      INTEGER(LONG)                   :: IERR
      INTEGER(LONG)                   :: IOS
      INTEGER(LONG)                   :: ICMT
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: SETID
      INTEGER(LONG)                   :: TEXT_BEG
      INTEGER(LONG)                   :: TEXT_END
      LOGICAL                         :: SEEN_SETID

      CARD1 = CARD
      CALL GPSTRESS_TO_UPPER_LINE ( CARD1 )
      ICMT = INDEX(CARD1,'$')
      IF (ICMT > 0) CARD1(ICMT:) = ' '

      KEY1 = ' '
      SETID = 0
      READ(CARD1,*,IOSTAT=IOS) KEY1, SETID
      IF ((IOS /= 0) .OR. (TRIM(KEY1) /= 'SET') .OR. (SETID <= 0)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR    : CANNOT PARSE GPSTRESS OUTPUT(POST) SET ENTRY: ', TRIM(CARD)
         WRITE(F06,'(A,A)') ' *ERROR    : CANNOT PARSE GPSTRESS OUTPUT(POST) SET ENTRY: ', TRIM(CARD)
         RETURN
      ENDIF

      SEEN_SETID = .FALSE.
      TEXT_BEG = 0
      TEXT_END = 0
      K = 0
      DO K=5,CC_ENTRY_LEN
         IF (CARD1(K:K) /= ' ') THEN
            SEEN_SETID = .TRUE.
         ELSE IF (SEEN_SETID) THEN
            TEXT_BEG = K + 1
            EXIT
         ENDIF
      ENDDO
      IF (TEXT_BEG == 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR    : GPSTRESS OUTPUT(POST) SET ENTRY HAS NO SET DATA: ', TRIM(CARD)
         WRITE(F06,'(A,A)') ' *ERROR    : GPSTRESS OUTPUT(POST) SET ENTRY HAS NO SET DATA: ', TRIM(CARD)
         RETURN
      ENDIF
      IF (CARD1(TEXT_BEG:TEXT_BEG) == '=') TEXT_BEG = TEXT_BEG + 1
      DO WHILE ((TEXT_BEG <= CC_ENTRY_LEN) .AND. (CARD1(TEXT_BEG:TEXT_BEG) == ' '))
         TEXT_BEG = TEXT_BEG + 1
      ENDDO
      TEXT_END = CC_ENTRY_LEN
      DO K=CC_ENTRY_LEN,TEXT_BEG,-1
         IF (CARD1(K:K) /= ' ') THEN
            TEXT_END = K
            EXIT
         ENDIF
      ENDDO

      SET_TEXT = ' '
      IF (TEXT_END >= TEXT_BEG) THEN
         SET_TEXT(1:MIN(LSETLN,TEXT_END-TEXT_BEG+1)) = CARD1(TEXT_BEG:MIN(TEXT_END,TEXT_BEG+LSETLN-1))
      ENDIF
      CALL REGISTER_GP_POST_SET ( SETID, SET_TEXT, IERR )
      IF (IERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,I8)') ' *ERROR    : TOO MANY GPSTRESS OUTPUT(POST) SET DEFINITIONS. LAST SET ID = ', SETID
         WRITE(F06,'(A,I8)') ' *ERROR    : TOO MANY GPSTRESS OUTPUT(POST) SET DEFINITIONS. LAST SET ID = ', SETID
      ENDIF

      END SUBROUTINE PARSE_GPSTRESS_SET_CARD

      SUBROUTINE PARSE_GPSTRESS_SURFACE_CARD ( CARD )

      CHARACTER(LEN=*), INTENT(IN)    :: CARD

      CHARACTER(LEN=CC_ENTRY_LEN)     :: CARD1
      CHARACTER(16*BYTE)              :: KEY1, KEY2, KEY3, KEY4
      INTEGER(LONG)                   :: IERR
      INTEGER(LONG)                   :: IOS
      INTEGER(LONG)                   :: ICMT
      INTEGER(LONG)                   :: SETID
      INTEGER(LONG)                   :: SURFACE_ID

      CARD1 = CARD
      CALL GPSTRESS_TO_UPPER_LINE ( CARD1 )
      ICMT = INDEX(CARD1,'$')
      IF (ICMT > 0) CARD1(ICMT:) = ' '

      KEY1 = ' '
      KEY2 = ' '
      KEY3 = ' '
      KEY4 = ' '
      SURFACE_ID = 0
      SETID = 0

      READ(CARD1,*,IOSTAT=IOS) KEY1, SURFACE_ID, KEY2, SETID, KEY3, KEY4
      IF (IOS /= 0) THEN
         KEY1 = ' '
         KEY2 = ' '
         KEY3 = ' '
         KEY4 = ' '
         SURFACE_ID = 0
         SETID = 0
         READ(CARD1,*,IOSTAT=IOS) KEY1, SURFACE_ID, KEY2, SETID
      ENDIF

      IF (IOS /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR    : CANNOT PARSE GPSTRESS SURFACE ENTRY: ', TRIM(CARD)
         WRITE(F06,'(A,A)') ' *ERROR    : CANNOT PARSE GPSTRESS SURFACE ENTRY: ', TRIM(CARD)
         RETURN
      ENDIF

      IF ((TRIM(KEY1) /= 'SURFACE') .OR. (TRIM(KEY2) /= 'SET')) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR    : GPSTRESS SURFACE ENTRY MUST LOOK LIKE "SURFACE id SET setid [NORMAL dir]": ', TRIM(CARD)
         WRITE(F06,'(A,A)') ' *ERROR    : GPSTRESS SURFACE ENTRY MUST LOOK LIKE "SURFACE id SET setid [NORMAL dir]": ', TRIM(CARD)
         RETURN
      ENDIF

      IF ((SURFACE_ID <= 0) .OR. (SETID <= 0)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR    : GPSTRESS SURFACE ENTRY HAS NONPOSITIVE SURFACE/SET ID: ', TRIM(CARD)
         WRITE(F06,'(A,A)') ' *ERROR    : GPSTRESS SURFACE ENTRY HAS NONPOSITIVE SURFACE/SET ID: ', TRIM(CARD)
         RETURN
      ENDIF

      IF (TRIM(KEY3) == 'NORMAL') THEN
         KEY3 = KEY4
      ENDIF
      IF (LEN_TRIM(KEY3) == 0) KEY3 = 'BASIC'

      CALL REGISTER_GP_SURFACE ( SURFACE_ID, SETID, KEY3, IERR )
      IF (IERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,I8)') ' *ERROR    : TOO MANY GPSTRESS SURFACE DEFINITIONS. LAST SURFACE ID = ', SURFACE_ID
         WRITE(F06,'(A,I8)') ' *ERROR    : TOO MANY GPSTRESS SURFACE DEFINITIONS. LAST SURFACE ID = ', SURFACE_ID
      ENDIF

      END SUBROUTINE PARSE_GPSTRESS_SURFACE_CARD

      SUBROUTINE GPSTRESS_GET_SET_STRING ( SETID_IN, NULSET, SLEN, TOKSTR_OUT )

      INTEGER(LONG), INTENT(IN)       :: SETID_IN
      INTEGER(LONG), INTENT(OUT)      :: NULSET
      INTEGER(LONG), INTENT(OUT)      :: SLEN
      CHARACTER(LEN=*), INTENT(OUT)   :: TOKSTR_OUT

      CHARACTER(LSETLN*BYTE)          :: SETCHR
      CHARACTER(8*BYTE)               :: TOKENI
      CHARACTER(LSETLN*BYTE)          :: TOKSTR
      CHARACTER(8*BYTE)               :: TOKTY1
      INTEGER(LONG)                   :: DATA_BEG, DATA_END, ECOL, I, IOCHK, K, POSN, SETID
      INTEGER(LONG)                   :: SET_1_BEG, SET_2_BEG, SID_BEG, SID_END

      TOKSTR_OUT = ' '
      NULSET = 0
      SLEN   = 0

      DO I=1,NUM_GP_POST_SET
         IF (GP_POST_SET_IDS(I) == SETID_IN) THEN
            TOKSTR_OUT = ' '
            SLEN = LEN_TRIM(GP_POST_SET_TEXT(I))
            IF (SLEN > 0) THEN
               TOKSTR_OUT(1:MIN(SLEN,LEN(TOKSTR_OUT))) = GP_POST_SET_TEXT(I)(1:MIN(SLEN,LEN(TOKSTR_OUT)))
               NULSET = 1
            ENDIF
            RETURN
         ENDIF
      ENDDO

      DO I=1,LSETLN
         SETCHR(I:I) = ALL_SETS_ARRAY(I)
      ENDDO

      SET_1_BEG = INDEX(SETCHR(1:),'SET')
      IF (SET_1_BEG == 0) THEN
         WRITE(ERR,'(A,I8,A)') ' *ERROR  1405: SET ID ', SETID_IN, ' NOT FOUND'
         WRITE(F06,'(A,I8,A)') ' *ERROR  1405: SET ID ', SETID_IN, ' NOT FOUND'
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

outer:DO
         TOKSTR(1:) = ' '
         POSN = INDEX(SETCHR(SET_1_BEG+1:),'SET')
         IF (POSN == 0) THEN
            SET_2_BEG = SETLEN + 1
         ELSE
            SET_2_BEG = SET_1_BEG + POSN
         ENDIF

         POSN = INDEX(SETCHR(SET_1_BEG:),'=')
         ECOL = SET_1_BEG + POSN - 1
         DATA_BEG = ECOL + 1
         DATA_END = SET_2_BEG - 1

         DO I=DATA_BEG,DATA_END
            IF ((SETCHR(I:I) == ' ') .OR. (SETCHR(I:I) == ',')) THEN
               DATA_BEG = DATA_BEG + 1
            ELSE
               EXIT
            ENDIF
         ENDDO

         TOKENI(1:) = ' '
         SID_BEG = SET_1_BEG + 3
         SID_END = ECOL - 1
         SETID = 0
         K = 0
         DO I=SID_BEG,SID_END
            IF ((SETCHR(I:I) == ' ') .OR. (SETCHR(I:I) == ',')) CYCLE
            K = K + 1
            IF (K > MAX_TOKEN_LEN) THEN
               NULSET = 0
               SET_1_BEG = SET_2_BEG
               CYCLE outer
            ENDIF
            TOKENI(K:K) = SETCHR(I:I)
         ENDDO

         CALL TOKCHK ( TOKENI, TOKTY1 )
         IF (TOKTY1 == 'INTEGER ') THEN
            READ(TOKENI,'(I8)',IOSTAT=IOCHK) SETID
            IF (IOCHK /= 0) THEN
               NULSET = 0
               SET_1_BEG = SET_2_BEG
               CYCLE outer
            ENDIF
         ELSE
            NULSET = 0
            SET_1_BEG = SET_2_BEG
            CYCLE outer
         ENDIF

         IF (SETID == SETID_IN) THEN
            NULSET = 1
            SLEN = DATA_END - DATA_BEG + 1
            TOKSTR(1:SLEN) = SETCHR(DATA_BEG:DATA_END)
            EXIT
         ELSE
            NULSET = 0
            SLEN = 0
            IF (SET_2_BEG > SETLEN) THEN
               WRITE(ERR,'(A,I8,A)') ' *ERROR  1405: SET ID ', SETID_IN, ' NOT FOUND'
               WRITE(F06,'(A,I8,A)') ' *ERROR  1405: SET ID ', SETID_IN, ' NOT FOUND'
               FATAL_ERR = FATAL_ERR + 1
               RETURN
            ELSE
               SET_1_BEG = SET_2_BEG
               CYCLE outer
            ENDIF
         ENDIF
      ENDDO outer

      IF (NULSET /= 0) THEN
         TOKSTR_OUT(1:MIN(SLEN,LEN(TOKSTR_OUT))) = TOKSTR(1:MIN(SLEN,LEN(TOKSTR_OUT)))
      ENDIF

      END SUBROUTINE GPSTRESS_GET_SET_STRING

      SUBROUTINE GPSTRESS_COLLECT_SURFACE_PATCH ( SURF_INDEX, ELEM_IDS, NUM_ELEMS, GRID_IDS, NUM_GRIDS, IERR )

      INTEGER(LONG), INTENT(IN)       :: SURF_INDEX
      INTEGER(LONG), INTENT(OUT)      :: ELEM_IDS(:)
      INTEGER(LONG), INTENT(OUT)      :: NUM_ELEMS
      INTEGER(LONG), INTENT(OUT)      :: GRID_IDS(:)
      INTEGER(LONG), INTENT(OUT)      :: NUM_GRIDS
      INTEGER(LONG), INTENT(OUT)      :: IERR

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GPSTRESS_COLLECT_SURFACE_PATCH'
      CHARACTER(LSETLN*BYTE)          :: ERRTOK
      CHARACTER(3*BYTE)               :: EXCEPT
      CHARACTER(3*BYTE)               :: THRU
      CHARACTER(8*BYTE)               :: TOKEN(3)
      CHARACTER(8*BYTE)               :: TOKTYP(3)
      CHARACTER(LSETLN*BYTE)          :: TOKSTR
      INTEGER(LONG)                   :: ACT_ELEM, AELEM_HI, AELEM_LO, ESORT1_ROW_NUM, I, IERROR, ISTART, K
      INTEGER(LONG)                   :: NTOKEN, NULSET, SLEN, TOKLEN
      INTEGER(LONG)                   :: SHELL_GRIDS(8), NGR
      LOGICAL                         :: KEEP

      ELEM_IDS = 0
      GRID_IDS = 0
      NUM_ELEMS = 0
      NUM_GRIDS = 0
      IERR = 0

      IF ((SURF_INDEX < 1) .OR. (SURF_INDEX > NUM_GP_SURFACE)) THEN
         IERR = 1
         RETURN
      ENDIF

      CALL GPSTRESS_GET_SET_STRING ( GP_SURFACE_SETIDS(SURF_INDEX), NULSET, SLEN, TOKSTR )
      IF (NULSET == 0) THEN
         IERR = 2
         RETURN
      ENDIF

      IF (TRIM(TOKSTR(1:SLEN)) == 'ALL') THEN
         DO I=1,SIZE(ESORT1)
            IF ((ETYPE(I)(1:5) == 'TRIA3') .OR. (ETYPE(I)(1:5) == 'QUAD4') .OR. (ETYPE(I) == 'QUADR   ')) THEN
               CALL ADD_ELEM_IF_SHELL ( ESORT1(I) )
            ENDIF
         ENDDO
         GOTO 900
      ENDIF

      ISTART = 1
      THRU   = 'OFF'
      EXCEPT = 'OFF'
      TOKLEN = SLEN

token_loop: DO
         CALL STOKEN ( SUBR_NAME, TOKSTR, ISTART, TOKLEN, NTOKEN, IERROR, TOKTYP, TOKEN, ERRTOK, THRU, EXCEPT )
         IF (IERROR /= 0) THEN
            IERR = 3
            RETURN
         ENDIF
         IF (NTOKEN <= 0) EXIT token_loop

         IF (NTOKEN == 1) THEN
            IF (TOKTYP(1) /= 'INTEGER ') CYCLE
            READ(TOKEN(1),'(I8)') ACT_ELEM
            IF (EXCEPT == 'OFF') CALL ADD_ELEM_IF_SHELL ( ACT_ELEM )
            IF (EXCEPT == 'ON ') THEN
               IF ((ACT_ELEM >= AELEM_LO) .AND. (ACT_ELEM <= AELEM_HI)) THEN
                  CALL REMOVE_ELEM ( ACT_ELEM )
               ELSE
                  CALL ADD_ELEM_IF_SHELL ( ACT_ELEM )
                  EXCEPT = 'OFF'
                  THRU   = 'OFF'
               ENDIF
            ENDIF

         ELSE IF (NTOKEN == 3) THEN
            IF ((TOKTYP(1) == 'INTEGER ') .AND. (TOKTYP(2) == 'THRU    ') .AND. (TOKTYP(3) == 'INTEGER ')) THEN
               READ(TOKEN(1),'(I8)') AELEM_LO
               READ(TOKEN(3),'(I8)') AELEM_HI
               IF (AELEM_LO <= AELEM_HI) THEN
                  DO ACT_ELEM=AELEM_LO,AELEM_HI
                     CALL ADD_ELEM_IF_SHELL ( ACT_ELEM )
                  ENDDO
               ENDIF
            ENDIF
         ENDIF

         IF (ISTART > TOKLEN) EXIT token_loop
      ENDDO token_loop

      DO I=1,NUM_ELEMS
         CALL GET_SHELL_CORNER_GRIDS ( ELEM_IDS(I), SHELL_GRIDS, NGR, KEEP )
         IF (.NOT. KEEP) CYCLE
         DO K=1,NGR
            CALL ADD_UNIQUE_INT ( SHELL_GRIDS(K), GRID_IDS, NUM_GRIDS )
         ENDDO
      ENDDO

  900 CONTINUE
      DO I=1,NUM_ELEMS
         CALL GET_SHELL_CORNER_GRIDS ( ELEM_IDS(I), SHELL_GRIDS, NGR, KEEP )
         IF (.NOT. KEEP) CYCLE
         DO K=1,NGR
            CALL ADD_UNIQUE_INT ( SHELL_GRIDS(K), GRID_IDS, NUM_GRIDS )
         ENDDO
      ENDDO

      CONTAINS

      SUBROUTINE ADD_ELEM_IF_SHELL ( ACTUAL_ELEM_ID )
      INTEGER(LONG), INTENT(IN)       :: ACTUAL_ELEM_ID
      INTEGER(LONG)                   :: J

      CALL GET_ARRAY_ROW_NUM ( 'ESORT1', SUBR_NAME, SIZE(ESORT1), ESORT1, ACTUAL_ELEM_ID, ESORT1_ROW_NUM )
      IF (ESORT1_ROW_NUM == -1) RETURN

      IF ((ETYPE(ESORT1_ROW_NUM)(1:5) == 'TRIA3') .OR. (ETYPE(ESORT1_ROW_NUM)(1:5) == 'QUAD4') .OR.                     &
          (ETYPE(ESORT1_ROW_NUM)      == 'QUADR   ')) THEN
         DO J=1,NUM_ELEMS
            IF (ELEM_IDS(J) == ACTUAL_ELEM_ID) RETURN
         ENDDO
         IF (NUM_ELEMS < SIZE(ELEM_IDS)) THEN
            NUM_ELEMS = NUM_ELEMS + 1
            ELEM_IDS(NUM_ELEMS) = ACTUAL_ELEM_ID
         ELSE
            IERR = 4
         ENDIF
      ENDIF
      END SUBROUTINE ADD_ELEM_IF_SHELL

      SUBROUTINE REMOVE_ELEM ( ACTUAL_ELEM_ID )
      INTEGER(LONG), INTENT(IN)       :: ACTUAL_ELEM_ID
      INTEGER(LONG)                   :: J
      DO J=1,NUM_ELEMS
         IF (ELEM_IDS(J) == ACTUAL_ELEM_ID) THEN
            IF (J < NUM_ELEMS) THEN
               ELEM_IDS(J:NUM_ELEMS-1) = ELEM_IDS(J+1:NUM_ELEMS)
            ENDIF
            ELEM_IDS(NUM_ELEMS) = 0
            NUM_ELEMS = NUM_ELEMS - 1
            EXIT
         ENDIF
      ENDDO
      END SUBROUTINE REMOVE_ELEM

      END SUBROUTINE GPSTRESS_COLLECT_SURFACE_PATCH

      SUBROUTINE GET_SHELL_CORNER_GRIDS ( ACTUAL_ELEM_ID, SHELL_GRIDS, NGR, KEEP )

      INTEGER(LONG), INTENT(IN)       :: ACTUAL_ELEM_ID
      INTEGER(LONG), INTENT(OUT)      :: SHELL_GRIDS(8)
      INTEGER(LONG), INTENT(OUT)      :: NGR
      LOGICAL, INTENT(OUT)            :: KEEP

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_SHELL_CORNER_GRIDS'
      INTEGER(LONG)                   :: EP, I, ROW

      SHELL_GRIDS = 0
      NGR = 0
      KEEP = .FALSE.

      CALL GET_ARRAY_ROW_NUM ( 'ESORT1', SUBR_NAME, SIZE(ESORT1), ESORT1, ACTUAL_ELEM_ID, ROW )
      IF (ROW == -1) RETURN

      EP = EPNT(ROW)
      IF (ETYPE(ROW)(1:5) == 'TRIA3') THEN
         NGR = 3
         DO I=1,NGR
            SHELL_GRIDS(I) = EDAT(EP + 1 + I)
         ENDDO
         KEEP = .TRUE.
      ELSE IF ((ETYPE(ROW)(1:5) == 'QUAD4') .OR. (ETYPE(ROW) == 'QUADR   ')) THEN
         NGR = 4
         DO I=1,NGR
            SHELL_GRIDS(I) = EDAT(EP + 1 + I)
         ENDDO
         KEEP = .TRUE.
      ENDIF

      END SUBROUTINE GET_SHELL_CORNER_GRIDS

      SUBROUTINE ADD_UNIQUE_INT ( VALUE, ARRAY, NUSED )

      INTEGER(LONG), INTENT(IN)       :: VALUE
      INTEGER(LONG), INTENT(INOUT)    :: ARRAY(:)
      INTEGER(LONG), INTENT(INOUT)    :: NUSED

      INTEGER(LONG)                   :: I

      IF (VALUE <= 0) RETURN
      DO I=1,NUSED
         IF (ARRAY(I) == VALUE) RETURN
      ENDDO
      IF (NUSED < SIZE(ARRAY)) THEN
         NUSED = NUSED + 1
         ARRAY(NUSED) = VALUE
      ENDIF

      END SUBROUTINE ADD_UNIQUE_INT

      SUBROUTINE GPSTRESS_TO_UPPER_LINE ( LINE )

      CHARACTER(LEN=*), INTENT(INOUT) :: LINE
      INTEGER(LONG)                   :: I, ICODE

      DO I=1,LEN(LINE)
         ICODE = ICHAR(LINE(I:I))
         IF ((ICODE >= ICHAR('a')) .AND. (ICODE <= ICHAR('z'))) THEN
            LINE(I:I) = CHAR(ICODE - 32)
         ENDIF
      ENDDO

      END SUBROUTINE GPSTRESS_TO_UPPER_LINE

      END MODULE GPSTRESS_SURFACE_UTILS
