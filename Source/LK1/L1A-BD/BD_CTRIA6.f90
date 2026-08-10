! ##################################################################################################################################
      SUBROUTINE BD_CTRIA6 ( CARD, LARGE_FLD_INP, NUM_GRD )

! Processes CTRIA6 Bulk Data Cards.
! Parent card layout: EID, PID, G1, G2, G3, G4, G5, G6.
! Optional continuation: THETA/MCID, ZOFFS, T1, T2, T3.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, IERRFL, FATAL_ERR, JCARD_LEN, JF, LMATANGLE, LPLATEOFF, LPLATETHICK,       &
                                         MEDAT_CTRIA6, NCTRIA6, NEDAT, NELE, NMATANGLE, NPLATEOFF, NPLATETHICK
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  EDAT, ETYPE, MATANGLE, PLATEOFF, PLATETHICK

      USE MKJCARD_Interface
      USE ELEPRO_Interface
      USE TOKCHK_Interface
      USE OUTA_HERE_Interface
      USE R8FLD_Interface
      USE I4FLD_Interface
      USE BD_IMBEDDED_BLANK_Interface
      USE CRDERR_Interface
      USE NEXTC_Interface
      USE NEXTC2_Interface
      USE CARD_FLDS_NOT_BLANK_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CTRIA6'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP
      CHARACTER(LEN(CARD))            :: CHILD
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)
      CHARACTER( 8*BYTE)              :: TOKEN
      CHARACTER( 8*BYTE)              :: TOKTYP

      INTEGER(LONG), INTENT(OUT)      :: NUM_GRD
      INTEGER(LONG)                   :: I, J
      INTEGER(LONG)                   :: I4INP
      INTEGER(LONG)                   :: INT41, INT42
      INTEGER(LONG)                   :: ICONT = 0
      INTEGER(LONG)                   :: IERR  = 0

      REAL(DOUBLE)                    :: R8INP = ZERO

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

      IF (JCARD(3)(1:) == ' ') THEN
         JCARD_EDAT(3) = JCARD(2)
      ENDIF

      CALL ELEPRO ( 'Y', JCARD_EDAT, 8, MEDAT_CTRIA6, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y' )

      NUM_GRD = 6
      NCTRIA6 = NCTRIA6 + 1
      ETYPE(NELE) = 'TRIA6   '

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,6,7,8,9 )
      CALL CRDERR ( CARD )

      INT41 = 0
      INT42 = 0

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      IF (ICONT == 1) THEN
         IF (JCARD(2)(1:) /= ' ') THEN
            TOKEN = JCARD(2)(1:8)
            CALL TOKCHK ( TOKEN, TOKTYP )
            IF      ((TOKTYP /= 'INTEGER ') .AND. (TOKTYP /= 'FL PT   ')) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE (ERR,1196) JCARD(1), JCARD(2)
               WRITE (F06,1196) JCARD(1), JCARD(2)
            ELSE IF (TOKTYP == 'FL PT   ') THEN
               NMATANGLE = NMATANGLE + 1
               IF (NMATANGLE > LMATANGLE) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1141) SUBR_NAME,LMATANGLE
                  WRITE(F06,1141) SUBR_NAME,LMATANGLE
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
               INT41 = NMATANGLE
               INT42 = 0
               CALL R8FLD ( JCARD(2), JF(2), R8INP )
               IF (IERRFL(2) == 'N') MATANGLE(NMATANGLE) = R8INP
            ELSE IF (TOKTYP == 'INTEGER ') THEN
               CALL I4FLD ( JCARD(2), JF(2), I4INP )
               INT41 = -I4INP
               IF      (I4INP > 0) THEN
                  INT42 = 2
               ELSE IF (I4INP == 0) THEN
                  INT42 = 1
               ELSE
                  INT42 = 0
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1197) JCARD(1), JF(2), JCARD(2)
                  WRITE(F06,1197) JCARD(1), JF(2), JCARD(2)
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      NEDAT = NEDAT + 1
      EDAT(NEDAT) = INT41
      NEDAT = NEDAT + 1
      EDAT(NEDAT) = INT42

      IF ((ICONT == 1) .AND. (JCARD(3)(1:) /= ' ')) THEN
         NPLATEOFF = NPLATEOFF + 1
         IF (NPLATEOFF > LPLATEOFF) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1144) SUBR_NAME,' TOO MANY PLATE OFFSETS. LIMIT IS LPLATEOFF =  ',LPLATEOFF
            WRITE(F06,1144) SUBR_NAME,' TOO MANY PLATE OFFSETS. LIMIT IS LPLATEOFF =  ',LPLATEOFF
            CALL OUTA_HERE ( 'Y' )
         ENDIF
         NEDAT = NEDAT + 1
         EDAT(NEDAT) = NPLATEOFF
         CALL R8FLD ( JCARD(3), JF(3), R8INP )
         IF (IERRFL(3) == 'N') PLATEOFF(NPLATEOFF) = R8INP
      ELSE
         NEDAT = NEDAT + 1
         EDAT(NEDAT) = 0
      ENDIF

      NEDAT = NEDAT + 1
      EDAT(NEDAT) = 0

      NEDAT = NEDAT + 1
      EDAT(NEDAT) = 0
      IF ((ICONT == 1) .AND. ((JCARD(4)(1:) /= ' ') .OR. (JCARD(5)(1:) /= ' ') .OR. (JCARD(6)(1:) /= ' '))) THEN
         EDAT(NEDAT) = NPLATETHICK + 1
         DO J=4,6
            NPLATETHICK = NPLATETHICK + 1
            IF (NPLATETHICK > LPLATETHICK) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1144) SUBR_NAME,' TOO MANY PLATE THICKNESSES. LIMIT IS LPLATETHICK = ',LPLATETHICK
               WRITE(F06,1144) SUBR_NAME,' TOO MANY PLATE THICKNESSES. LIMIT IS LPLATETHICK = ',LPLATETHICK
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            CALL R8FLD ( JCARD(J), JF(J), R8INP )
            IF (IERRFL(J) == 'N') PLATETHICK(NPLATETHICK) = R8INP
         ENDDO
      ENDIF

      IF (ICONT == 1) THEN
         CALL BD_IMBEDDED_BLANK   ( JCARD,0,2,3,4,5,6,0,0 )
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,7,8 )
         CALL CRDERR ( CARD )
      ENDIF

      RETURN

 1141 FORMAT(' *ERROR  1141: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY PLATE ELEMENT MATERIAL PROPERTY ANGLES. LIMIT IS NMATANGLE =  ',I8)
 1144 FORMAT(' *ERROR  1144: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A,I8)
 1196 FORMAT(' *ERROR  1196: VALUE FOR MATERIAL ANGLE ON ',A,A,' MUST BE AN INTEGER OR REAL NUMBER')
 1197 FORMAT(' *ERROR  1197: FOR ',A,' THE COORD SYS ID IN FIELD ',I2,' MUST BE >= 0. HOWEVER, THE VALUE INPUT WAS ',A)

      END SUBROUTINE BD_CTRIA6
