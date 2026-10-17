      SUBROUTINE BD_PYRAM ( CARD, LARGE_FLD_INP, NUM_GRD )

! Standard Nastran 5/13-node connectivity; reject nonstandard 14-node input.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, JCARD_LEN, NPYRAM5, NPYRAM13, NELE
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  ETYPE

      USE BD_PYRAM_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PYRAM'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP
      CHARACTER(LEN(CARD))            :: CHILD
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)
      CHARACTER(LEN(JCARD))           :: ID
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)
      CHARACTER(LEN(JCARD))           :: NAME

      INTEGER(LONG), INTENT(OUT)      :: NUM_GRD
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: NF
      INTEGER(LONG)                   :: ICONT = 0
      INTEGER(LONG)                   :: IERR  = 0

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      NAME = JCARD(1)
      ID   = JCARD(2)
      IF ((TRIM(NAME) /= 'CPYRAM') .AND. (TRIM(NAME) /= 'CPYRA') .AND. &
          (TRIM(NAME) /= 'CPYRAM*') .AND. (TRIM(NAME) /= 'CPYRA*')) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR: NONSTANDARD PYRAMID CARD NAME ',TRIM(NAME)
         WRITE(F06,'(A,A)') ' *ERROR: NONSTANDARD PYRAMID CARD NAME ',TRIM(NAME)
         NUM_GRD = 0
         RETURN
      ENDIF

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

      NF = 7
      IF (JCARD(9) /= ' ') NF = 8
      CALL ELEPRO ( 'Y', JCARD_EDAT, NF, NF, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y' )
      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )
      CALL CRDERR ( CARD )
      NUM_GRD = NF - 2
      DO
         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         IF (ICONT == 0) EXIT
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         NF = 0
         DO I=2,9
            IF (JCARD(I) /= ' ') NF = I - 1
         ENDDO
         IF (NF == 0) EXIT
         IF (NUM_GRD + NF > 13) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1137) NAME, ID, NUM_GRD + NF
            WRITE(F06,1137) NAME, ID, NUM_GRD + NF
            EXIT
         ENDIF
         CALL ELEPRO ( 'N', JCARD, NF, NF, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y' )
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )
         CALL CRDERR ( CARD )
         NUM_GRD = NUM_GRD + NF
      ENDDO
      SELECT CASE (NUM_GRD)
      CASE (5)
         NPYRAM5 = NPYRAM5 + 1
         ETYPE(NELE) = 'PYRAM5  '
      CASE (13)
         NPYRAM13 = NPYRAM13 + 1
         ETYPE(NELE) = 'PYRAM13 '
      CASE DEFAULT
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1137) NAME, ID, NUM_GRD
         WRITE(F06,1137) NAME, ID, NUM_GRD
      END SELECT

      RETURN

 1137 FORMAT(' *ERROR  1137: ',A,' ID = ',A,' HAS ',I4,' GRIDS; EXPECTED 5 OR 13')

      END SUBROUTINE BD_PYRAM
