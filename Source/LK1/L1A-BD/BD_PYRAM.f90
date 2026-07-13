      SUBROUTINE BD_PYRAM ( CARD, LARGE_FLD_INP, NUM_GRD )

! Processes pyramid solid Bulk Data Cards for 5-node linear and 14-node
! quadratic elements. Legacy CPYRAM/CPYRA card names remain accepted on input.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, JCARD_LEN, NPYRAM5, NPYRAM14, NELE
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
      INTEGER(LONG)                   :: ICONT = 0
      INTEGER(LONG)                   :: IERR  = 0

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      NAME = JCARD(1)
      ID   = JCARD(2)

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

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
         NPYRAM5     = NPYRAM5 + 1
         ETYPE(NELE) = 'PYRAM5  '
         NUM_GRD     = 5
      ELSE
         IF (CARD(1:) /= ' ') THEN
            NPYRAM14    = NPYRAM14 + 1
            ETYPE(NELE) = 'PYRAM14 '
            NUM_GRD      = 14

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
            NPYRAM5     = NPYRAM5 + 1
            ETYPE(NELE) = 'PYRAM5  '
            NUM_GRD     = 5
         ENDIF
      ENDIF

      RETURN

 1136 FORMAT(' *ERROR  1136: REQUIRED CONTINUATION FOR ',A,' ID = ',A,' MISSING')

      END SUBROUTINE BD_PYRAM
