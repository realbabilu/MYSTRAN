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

      SUBROUTINE CC_STRE ( CARD, IS_GPSTRESS_ALIAS )

      ! Processes Case Control STRE cards for element stress output requests.
      ! GPSTRESS/GSTRESS are accepted as a compatibility alias in this branch.
      ! They request corner-style shell stress recovery for the separate
      ! GPSTRESS surface/volume writer without changing ordinary STRESS output.
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_CMD_DESCRIBERS, ECHO, LSUB, NSUB, NCCCD, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRE_LOC, STRE_OUT, GPSTRESS_OUT, GPSTRESS_REQ, GPSTRESS_SETID, STRESS_USER_REQ
      USE MODEL_STUF, ONLY            :  SC_STRE
      USE PARAMS, ONLY                :  SUPWARN

      USE CC_STRE_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_STRE'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      LOGICAL, OPTIONAL, INTENT(IN)    :: IS_GPSTRESS_ALIAS ! =.TRUE. when called for GPSTRESS/GSTRESS
      CHARACTER( 1*BYTE)              :: FOUND_PRINT       ! CC_CMD_DESCRIBERS has request for "PRINT"
      CHARACTER( 1*BYTE)              :: FOUND_PLOT        ! CC_CMD_DESCRIBERS has request for "PLOT"
      CHARACTER( 1*BYTE)              :: FOUND_PUNCH       ! CC_CMD_DESCRIBERS has request for "PUNCH"
      CHARACTER( 1*BYTE)              :: FOUND_NEU         ! CC_CMD_DESCRIBERS has request for "NEU"
      CHARACTER( 1*BYTE)              :: FOUND_CSV         ! CC_CMD_DESCRIBERS has request for "CSV"
      CHARACTER(LEN(CC_CMD_DESCRIBERS)):: REQUEST_OUT
      LOGICAL                          :: GPSTRESS_ALIAS

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: RECOVERY_SETID    ! Element set used for hidden recovery
      INTEGER(LONG)                   :: SETID             ! Set ID on this Case Control card




! **********************************************************************************************************************************
      ! CC_OUTPUTS processes all output type Case Control entries (they all
      ! have some common code so it is put there)

      GPSTRESS_ALIAS = .FALSE.
      IF (PRESENT(IS_GPSTRESS_ALIAS)) GPSTRESS_ALIAS = IS_GPSTRESS_ALIAS

      CALL CC_OUTPUTS ( CARD, 'STRE', SETID )
      RECOVERY_SETID = SETID

      ! Check to see if PLOT, PRINT, PUNCH, NEU, CSV were in the request
      FOUND_PRINT = 'N'
      FOUND_PLOT  = 'N'
      FOUND_PUNCH = 'N'
      FOUND_NEU   = 'N'
      FOUND_CSV   = 'N'
      DO I=1,NCCCD
         IF (CC_CMD_DESCRIBERS(I)(1:5) == 'PRINT') FOUND_PRINT = 'Y'
         IF (CC_CMD_DESCRIBERS(I)(1:4) == 'PLOT')  FOUND_PLOT  = 'Y'
         IF (CC_CMD_DESCRIBERS(I)(1:5) == 'PUNCH') FOUND_PUNCH = 'Y'
         IF (CC_CMD_DESCRIBERS(I)(1:3) == 'NEU')   FOUND_NEU   = 'Y'
         IF (CC_CMD_DESCRIBERS(I)(1:3) == 'CSV')   FOUND_CSV   = 'Y'
      ENDDO
      ! concatenate the strings
      REQUEST_OUT = TRIM(FOUND_PRINT) // TRIM(FOUND_PLOT) // TRIM(FOUND_PUNCH) // TRIM(FOUND_NEU) // TRIM(FOUND_CSV)

      ! For bare "STRE = ALL" requests, default to PRINT+PLOT so classic
      ! OP2 stress tables are emitted without requiring explicit "(PLOT)"
      ! qualifiers on the Case Control entry.
      IF (REQUEST_OUT(1:5) == 'NNNNN') THEN
        REQUEST_OUT = 'YYNNN'
      ENDIF

      IF (GPSTRESS_ALIAS) THEN
         GPSTRESS_REQ = .TRUE.
         GPSTRESS_SETID = SETID
         GPSTRESS_OUT = REQUEST_OUT
         STRE_LOC = 'CORNER'
         RECOVERY_SETID = -1
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,901) CARD(1:LEN_TRIM(CARD))
         IF (SUPWARN == 'N') THEN
            IF (ECHO == 'NONE  ') THEN
               WRITE(F06,901) CARD(1:LEN_TRIM(CARD))
            ENDIF
         ENDIF
      ELSE
         STRESS_USER_REQ = .TRUE.
         STRE_OUT = REQUEST_OUT
         IF (GPSTRESS_REQ) THEN
            ! GPSTRESS needs grid/corner stress recovery even when a later
            ! STRESS(CENTER) request is present in the same deck.
            STRE_LOC = 'CORNER'
         ENDIF
      ENDIF

      ! Set CASE CONTROL output request variable to SETID
      IF (NSUB == 0) THEN
         DO I = 1,LSUB
            SC_STRE(I) = RECOVERY_SETID
         ENDDO
      ELSE
         SC_STRE(NSUB) = RECOVERY_SETID
      ENDIF



      RETURN

! **********************************************************************************************************************************
  901 FORMAT(' *WARNING    : ',A,/,14X,' IS MAPPED TO STRESS(CORNER) FOR GPSTRESS RECOVERY IN THIS MYSTRAN BUILD.',/, &
             14X,' OUTPUT(POST) SURFACE/VOLUME CARDS ARE ACCEPTED; OGS1 OUTPUT IS A BASELINE COMPATIBILITY WRITER.')

      END SUBROUTINE CC_STRE
