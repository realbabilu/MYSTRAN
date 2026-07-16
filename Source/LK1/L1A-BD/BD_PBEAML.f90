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

      SUBROUTINE BD_PBEAML ( CARD, LARGE_FLD_INP )

! Processes a minimal phase-1 PBEAML Bulk Data card.
! Current support is intentionally narrow:
!   - selected section types with direct property mapping into internal PBEAM
!   - small/free-field continuation chain only
! The parsed section is mapped into the existing internal PBEAM storage so the
! current CBEAM station runtime can be reused without a separate formulation path.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, IN1
      USE PARAMS, ONLY                :  BEAMAMO, BEAMAMO_PID, BEAMAMO_VAL, BEAMM1MO, BEAMM1MO_PID, BEAMM1MO_VAL,      &
                                         BEAMM2MO, BEAMM2MO_PID, BEAMM2MO_VAL, BEAMTMO, BEAMTMO_PID, BEAMTMO_VAL,      &
                                         BEAMV1MO, BEAMV1MO_PID, BEAMV1MO_VAL, BEAMV2MO, BEAMV2MO_PID, BEAMV2MO_VAL,    &
                                         EPSIL, MBEAMAMO_PID, MBEAMM1MO_PID, MBEAMM2MO_PID, MBEAMTMO_PID,              &
                                         MBEAMV1MO_PID, MBEAMV2MO_PID, NBEAMAMO_PID, NBEAMM1MO_PID, NBEAMM2MO_PID,      &
                                         NBEAMTMO_PID, NBEAMV1MO_PID, NBEAMV2MO_PID, SUPINFO
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LPBEAM, NPBEAM,             &
                                          MPBEAM_STATIONS, WARN_ERR
      USE CONSTANTS_1, ONLY           :  ONE, ZERO
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  PBEAM, PBEAM_NSTATIONS, PBEAM_XL, PBEAM_RPROPS, RPBEAM

      USE BD_PBEAML_USE_IFs

      IMPLICIT NONE

      INTEGER(LONG), PARAMETER        :: MAXTOK = 512
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PBEAML'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP
      CHARACTER(LEN=LEN(CARD))        :: CARD_WORK
      CHARACTER(LEN=LEN(CARD))        :: RAW_LINE
      CHARACTER(LEN=LEN(CARD))        :: PARSE_LINE
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)
      CHARACTER(LEN=JCARD_LEN)        :: TOKENS(MAXTOK)
      CHARACTER(LEN=JCARD_LEN)        :: WORDS(32)
      CHARACTER(LEN=JCARD_LEN)        :: CARD_NAME
      CHARACTER(LEN=JCARD_LEN)        :: ID
      CHARACTER(LEN=JCARD_LEN)        :: SEC_TYPE
      CHARACTER(LEN=JCARD_LEN)        :: SOFLAG
      CHARACTER(LEN=JCARD_LEN)        :: TOKEN_WORK

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: IERR
      INTEGER(LONG)                   :: ICOM
      INTEGER(LONG)                   :: NDIM_SEC
      INTEGER(LONG)                   :: IFIRST
      INTEGER(LONG)                   :: IOCHK
      INTEGER(LONG)                   :: ISTA
      INTEGER(LONG)                   :: ITOK
      INTEGER(LONG)                   :: MATERIAL_ID
      INTEGER(LONG)                   :: NWORDS
      INTEGER(LONG)                   :: NTOK
      INTEGER(LONG)                   :: PROPERTY_ID
      INTEGER(LONG)                   :: STATION_COUNT
      INTEGER(LONG)                   :: NSEG_STATIONS
      INTEGER(LONG)                   :: NSTATION_EXTRA
      INTEGER(LONG)                   :: TAPER_MODE

      REAL(DOUBLE)                    :: DIMS_A(10)
      REAL(DOUBLE)                    :: DIMS_B(10)
      REAL(DOUBLE)                    :: DIMS_CUR(10)
      REAL(DOUBLE)                    :: NSM_A
      REAL(DOUBLE)                    :: NSM_B
      REAL(DOUBLE)                    :: NSM_CUR
      REAL(DOUBLE)                    :: AREA_MOD
      REAL(DOUBLE)                    :: I1_MOD
      REAL(DOUBLE)                    :: I2_MOD
      REAL(DOUBLE)                    :: K1_MOD
      REAL(DOUBLE)                    :: K2_MOD
      REAL(DOUBLE)                    :: J_MOD
      REAL(DOUBLE)                    :: ROFSET
      REAL(DOUBLE)                    :: STIFFMOD
      REAL(DOUBLE)                    :: STATION_XL
      REAL(DOUBLE)                    :: STATION_EXTRA(3)
      LOGICAL                         :: IS_PBEAMZ
      LOGICAL                         :: DIM0A_SEEN
      LOGICAL                         :: DIM1A_SEEN
      LOGICAL                         :: PBEAMZ_B_SECTION_READ

! **********************************************************************************************************************************
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      CARD_NAME = JCARD(1)
      CALL TO_UPPER ( CARD_NAME )
      IS_PBEAMZ = (CARD_NAME(1:6) == 'PBEAMZ')
      ID = JCARD(2)

      IF (LARGE_FLD_INP == 'Y') THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1301) TRIM(CARD_NAME)
         WRITE(F06,1301) TRIM(CARD_NAME)
         RETURN
      ENDIF

      CALL I4FLD ( JCARD(2), JF(2), PROPERTY_ID )
      CALL I4FLD ( JCARD(3), JF(3), MATERIAL_ID )

      SEC_TYPE = ADJUSTL(JCARD(5))
      IF (SEC_TYPE(1:1) == ' ') SEC_TYPE = ADJUSTL(JCARD(4))
      IF (SEC_TYPE(1:1) == ' ') THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1302) ID
         WRITE(F06,1302) ID
         RETURN
      ENDIF
      CALL TO_UPPER ( SEC_TYPE )
      NDIM_SEC = GET_SECTION_NDIMS ( SEC_TYPE )
      IF (NDIM_SEC <= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1303) ID, SEC_TYPE
         WRITE(F06,1303) ID, SEC_TYPE
         RETURN
      ENDIF

      NPBEAM = NPBEAM + 1
      PBEAM_NSTATIONS(NPBEAM) = 1
      PBEAM_XL(NPBEAM,1) = ZERO
      PBEAM(NPBEAM,1) = PROPERTY_ID
      PBEAM(NPBEAM,2) = MATERIAL_ID
      PBEAM(NPBEAM,3) = 1

      DO I=1,NPBEAM-1
         IF (PROPERTY_ID == PBEAM(I,1)) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1145) TRIM(CARD_NAME), PROPERTY_ID
            WRITE(F06,1145) TRIM(CARD_NAME), PROPERTY_ID
            RETURN
         ENDIF
      ENDDO

      DO I=1,MAXTOK
         TOKENS(I)(1:) = ' '
      ENDDO
      NTOK = 0
      CARD_WORK = CARD
      DIM0A_SEEN = .FALSE.
      DIM1A_SEEN = .FALSE.
      PBEAMZ_B_SECTION_READ = .FALSE.

collect_tokens: DO
         READ(IN1,'(A)',IOSTAT=IOCHK) RAW_LINE
         IF (IOCHK < 0) EXIT collect_tokens
         IF (IOCHK > 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1305) ID
            WRITE(F06,1305) ID
            RETURN
         ENDIF
         IFIRST = 0
         DO I=1,LEN(RAW_LINE)
            IF (RAW_LINE(I:I) /= ' ') THEN
               IFIRST = I
               EXIT
            ENDIF
         ENDDO
         IF (IFIRST == 0) CYCLE collect_tokens
         IF (RAW_LINE(IFIRST:IFIRST) == '$') CYCLE collect_tokens
         IF ((RAW_LINE(IFIRST:IFIRST) /= '+') .AND. (RAW_LINE(IFIRST:IFIRST) /= ',')) THEN
            BACKSPACE(IN1)
            EXIT collect_tokens
         ENDIF

         PARSE_LINE = RAW_LINE
         ICOM = INDEX(PARSE_LINE, '$')
         IF (ICOM > 0) PARSE_LINE(ICOM:) = ' '
         PARSE_LINE(IFIRST:IFIRST) = ' '

         CALL PARSE_CHAR_STRING ( PARSE_LINE, LEN_TRIM(PARSE_LINE), 32, JCARD_LEN, NWORDS, WORDS, IERR )
         DO I=1,NWORDS
            TOKEN_WORK = ADJUSTL(WORDS(I))
            IF (TOKEN_WORK(1:1) == ' ') CYCLE
            IF (IS_PBEAMZ .AND. IS_PBEAMZ_DIM_LABEL(TOKEN_WORK)) THEN
               CALL SET_PBEAMZ_DIM_LABEL_FLAGS ( TOKEN_WORK, DIM0A_SEEN, DIM1A_SEEN )
               CYCLE
            ENDIF
            IF (NTOK >= MAXTOK) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1304) ID, MAXTOK
               WRITE(F06,1304) ID, MAXTOK
               RETURN
            ENDIF
            NTOK = NTOK + 1
            TOKENS(NTOK) = TOKEN_WORK
         ENDDO
      ENDDO collect_tokens

      IF (IS_PBEAMZ .AND. (.NOT. DIM0A_SEEN)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR  1316: PBEAMZ ENTRY DOES NOT DEFINE DIM0A'
         WRITE(F06,'(A,A)') ' *ERROR  1316: PBEAMZ ENTRY DOES NOT DEFINE DIM0A'
         RETURN
      ENDIF

      IF (IS_PBEAMZ .AND. DIM1A_SEEN .AND. (.NOT. DIM0A_SEEN)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR  1316: PBEAMZ ENTRY HAS DIM1A BUT DOES NOT DEFINE DIM0A'
         WRITE(F06,'(A,A)') ' *ERROR  1316: PBEAMZ ENTRY HAS DIM1A BUT DOES NOT DEFINE DIM0A'
         RETURN
      ENDIF

      IF (IS_PBEAMZ) THEN
         CALL PARSE_PBEAMZ_TOKENS ( NDIM_SEC, NTOK, TOKENS, DIMS_A, DIMS_B, NSM_A, NSM_B, DIMS_CUR, STIFFMOD, ROFSET, TAPER_MODE, &
                                    AREA_MOD, I1_MOD, I2_MOD, K1_MOD, K2_MOD, J_MOD, NSEG_STATIONS, NSTATION_EXTRA, STATION_EXTRA )
         RETURN
      ENDIF

      IF (NTOK < NDIM_SEC) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR  1305: PBEAMZ ID = ', ID
         WRITE(F06,'(A,A)') ' *ERROR  1305: PBEAMZ ID = ', ID
         RETURN
      ENDIF

      DIMS_A = ZERO
      DIMS_B = ZERO
      DIMS_CUR = ZERO
      STATION_EXTRA = ZERO

      DO I=1,NDIM_SEC
         READ(TOKENS(I),*,ERR=900) DIMS_A(I)
      ENDDO
      ITOK = NDIM_SEC + 1
      NSM_A = ZERO
      NSEG_STATIONS = 10
      NSTATION_EXTRA = 0
      IF (ITOK <= NTOK) THEN
         IF (LEN_TRIM(TOKENS(ITOK)) > 0) THEN
            IF ((.NOT. IS_SO_TOKEN(TOKENS(ITOK))) .AND. (.NOT. (IS_PBEAMZ .AND. IS_PBEAMZ_OPTION(TOKENS(ITOK))))) THEN
               READ(TOKENS(ITOK),*,ERR=900) NSM_A
               ITOK = ITOK + 1
            ENDIF
         ENDIF
      ENDIF

      CALL LOAD_SECTION_A ( SEC_TYPE, NDIM_SEC, DIMS_A, NSM_A )
      DIMS_B = DIMS_A
      NSM_B  = NSM_A
      STATION_COUNT = 1
      STIFFMOD   = 1.0D0
      ROFSET     = 0.0D0
      TAPER_MODE = 0
      AREA_MOD   = 1.0D0
      I1_MOD     = 1.0D0
      I2_MOD     = 1.0D0
      K1_MOD     = 1.0D0
      K2_MOD     = 1.0D0
      J_MOD      = 1.0D0

! --- cbeam_pbeaml_constant begin --- !
! Support the standard constant-section PBEAML form where the only continuation
! fields are the section dimensions (and optional NSM) with no SO/XL station
! data. In that case the section is uniform from end A to end B.
      IF (ITOK > NTOK) THEN
         STATION_COUNT = 2
         PBEAM_NSTATIONS(NPBEAM) = STATION_COUNT
         PBEAM_XL(NPBEAM,STATION_COUNT) = 1.0D0
         CALL LOAD_SECTION_B ( SEC_TYPE, NDIM_SEC, DIMS_B, NSM_B )
         CALL STORE_PBEAMZ_META ( STIFFMOD, ROFSET, TAPER_MODE, AREA_MOD, I1_MOD, I2_MOD, K1_MOD, K2_MOD, J_MOD )
         IF (IS_PBEAMZ) CALL FINALIZE_PBEAMZ_GENERATED_PBEAM ( TAPER_MODE, NSEG_STATIONS, AREA_MOD, I1_MOD, I2_MOD, K1_MOD, K2_MOD, J_MOD, ROFSET, NSTATION_EXTRA, STATION_EXTRA )
         RETURN
      ENDIF
! --- cbeam_pbeaml_constant end --- !

station_parse: DO WHILE (ITOK <= NTOK)
         SOFLAG = TOKENS(ITOK)
         IF (.NOT. IS_SO_TOKEN(SOFLAG)) THEN
            IF (IS_PBEAMZ .AND. IS_PBEAMZ_OPTION(SOFLAG)) THEN
               CALL READ_PBEAMZ_OPTION ( SOFLAG, ITOK, STIFFMOD, ROFSET, TAPER_MODE, AREA_MOD, I1_MOD, I2_MOD, K1_MOD, K2_MOD, J_MOD, NSM_A, NSM_B, NSEG_STATIONS, NSTATION_EXTRA, STATION_EXTRA )
               CYCLE station_parse
            ENDIF
            IF (IS_PBEAMZ .AND. DIM1A_SEEN .AND. (.NOT. PBEAMZ_B_SECTION_READ)) THEN
               IF (ITOK + NDIM_SEC - 1 > NTOK) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1307) ID
                  WRITE(F06,1307) ID
                  RETURN
               ENDIF
               DO I=1,NDIM_SEC
                  READ(TOKENS(ITOK),*,ERR=900) DIMS_B(I)
                  ITOK = ITOK + 1
               ENDDO
               NSM_B = NSM_A
               IF (ITOK <= NTOK) THEN
                  IF ((.NOT. IS_SO_TOKEN(TOKENS(ITOK))) .AND. (.NOT. (IS_PBEAMZ .AND. IS_PBEAMZ_OPTION(TOKENS(ITOK))))) THEN
                     READ(TOKENS(ITOK),*,ERR=900) NSM_B
                     ITOK = ITOK + 1
                  ENDIF
               ENDIF
               PBEAMZ_B_SECTION_READ = .TRUE.
               CYCLE station_parse
            ENDIF
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1306) ID, TOKENS(ITOK)
            WRITE(F06,1306) ID, TOKENS(ITOK)
            RETURN
         ENDIF
         ITOK = ITOK + 1
         IF (ITOK + NDIM_SEC > NTOK + 1) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1307) ID
            WRITE(F06,1307) ID
            RETURN
         ENDIF

         READ(TOKENS(ITOK),*,ERR=900) STATION_XL
         ITOK = ITOK + 1
         DO I=1,NDIM_SEC
            READ(TOKENS(ITOK),*,ERR=900) DIMS_CUR(I)
            ITOK = ITOK + 1
         ENDDO
         NSM_CUR = ZERO
         IF (ITOK <= NTOK) THEN
            IF ((.NOT. IS_SO_TOKEN(TOKENS(ITOK))) .AND. (.NOT. (IS_PBEAMZ .AND. IS_PBEAMZ_OPTION(TOKENS(ITOK))))) THEN
               READ(TOKENS(ITOK),*,ERR=900) NSM_CUR
               ITOK = ITOK + 1
            ENDIF
         ENDIF

         IF (STATION_COUNT < MPBEAM_STATIONS) THEN
            STATION_COUNT = STATION_COUNT + 1
            PBEAM_XL(NPBEAM,STATION_COUNT) = STATION_XL
            PBEAM_NSTATIONS(NPBEAM) = STATION_COUNT
            CALL STORE_SECTION_PROPS ( STATION_COUNT, SEC_TYPE, NDIM_SEC, DIMS_CUR, NSM_CUR )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,1197) 'PBEAML', ID, MPBEAM_STATIONS
            WRITE(F06,1197) 'PBEAML', ID, MPBEAM_STATIONS
         ENDIF

         DIMS_B = DIMS_CUR
         NSM_B  = NSM_CUR
      ENDDO station_parse

      IF (PBEAM_NSTATIONS(NPBEAM) <= 1) THEN
         IF (IS_PBEAMZ) THEN
            PBEAM_NSTATIONS(NPBEAM) = 2
            PBEAM_XL(NPBEAM,2) = 1.0D0
            CALL LOAD_SECTION_B ( SEC_TYPE, NDIM_SEC, DIMS_B, NSM_B )
            CALL STORE_PBEAMZ_META ( STIFFMOD, ROFSET, TAPER_MODE, AREA_MOD, I1_MOD, I2_MOD, K1_MOD, K2_MOD, J_MOD )
            CALL FINALIZE_PBEAMZ_GENERATED_PBEAM ( TAPER_MODE, NSEG_STATIONS, AREA_MOD, I1_MOD, I2_MOD, K1_MOD, K2_MOD, J_MOD, ROFSET, NSTATION_EXTRA, STATION_EXTRA )
            RETURN
         ENDIF
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1308) ID
         WRITE(F06,1308) ID
         RETURN
      ENDIF

      CALL LOAD_SECTION_B ( SEC_TYPE, NDIM_SEC, DIMS_B, NSM_B )
      CALL STORE_PBEAMZ_META ( STIFFMOD, ROFSET, TAPER_MODE, AREA_MOD, I1_MOD, I2_MOD, K1_MOD, K2_MOD, J_MOD )
      IF (IS_PBEAMZ) CALL FINALIZE_PBEAMZ_GENERATED_PBEAM ( TAPER_MODE, NSEG_STATIONS, AREA_MOD, I1_MOD, I2_MOD, K1_MOD, K2_MOD, J_MOD, ROFSET, NSTATION_EXTRA, STATION_EXTRA )

      RETURN

  900 CONTINUE
      FATAL_ERR = FATAL_ERR + 1
      WRITE(ERR,1309) ID
      WRITE(F06,1309) ID
      RETURN

! **********************************************************************************************************************************
 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)
 1197 FORMAT(' *WARNING 1197: ',A,' ID = ',A,' HAS MORE THAN ',I8,' STORED STATIONS. EXTRA x/L VALUES WILL BE IGNORED')
 1301 FORMAT(' *ERROR  1301: ',A,' IS NOT SUPPORTED IN THIS PHASE-1 PBEAML IMPLEMENTATION')
 1302 FORMAT(' *ERROR  1302: PBEAML ID = ',A,' IS MISSING A SECTION TYPE')
 1303 FORMAT(' *ERROR  1303: PBEAML ID = ',A,' HAS UNSUPPORTED SECTION TYPE "',A,'" IN THIS PHASE-1 IMPLEMENTATION')
 1304 FORMAT(' *ERROR  1304: PBEAML ID = ',A,' EXCEEDED TOKEN BUFFER OF ',I8,' CONTINUATION TOKENS')
 1305 FORMAT(' *ERROR  1305: PBEAML ID = ',A,' DOES NOT CONTAIN ENOUGH DIMENSION DATA')
 1306 FORMAT(' *ERROR  1306: PBEAML ID = ',A,' EXPECTED SO TOKEN YES/YESA/NO BUT FOUND "',A,'"')
 1307 FORMAT(' *ERROR  1307: PBEAML ID = ',A,' HAS AN INCOMPLETE STATION DEFINITION')
 1308 FORMAT(' *ERROR  1308: PBEAML ID = ',A,' DID NOT DEFINE ANY STATION BEYOND END A')
 1309 FORMAT(' *ERROR  1309: PBEAML ID = ',A,' HAS NONNUMERIC TYPE=I DIMENSION DATA')
 1310 FORMAT(' *ERROR  1310: PBEAML ID = ',A,' CURRENTLY REQUIRES FREE-FIELD COMMA INPUT IN THIS PHASE-1 IMPLEMENTATION')
 1311 FORMAT(' *ERROR  1311: ',A,' ENTRY HAS OPTION "',A,'" WITHOUT A FOLLOWING VALUE')
 1312 FORMAT(' *ERROR  1312: ',A,' ENTRY HAS INVALID VALUE "',A,'" FOR OPTION "',A,'"')
 1313 FORMAT(' *ERROR  1313: ',A,' ENTRY HAS INVALID TAPER VALUE "',A,'" IN TOKEN "',A,'"')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      LOGICAL FUNCTION IS_SO_TOKEN ( TOKEN )

      CHARACTER(LEN=*), INTENT(IN) :: TOKEN
      CHARACTER(LEN=JCARD_LEN)     :: TOKEN_UP

      TOKEN_UP = TOKEN
      CALL TO_UPPER ( TOKEN_UP )
      IF ((TOKEN_UP(1:3) == 'YES') .OR. (TOKEN_UP(1:4) == 'YESA') .OR. (TOKEN_UP(1:2) == 'NO')) THEN
         IS_SO_TOKEN = .TRUE.
      ELSE
         IS_SO_TOKEN = .FALSE.
      ENDIF

      END FUNCTION IS_SO_TOKEN

! ##################################################################################################################################

      INTEGER(LONG) FUNCTION GET_SECTION_NDIMS ( SEC )

      CHARACTER(LEN=*), INTENT(IN) :: SEC

      GET_SECTION_NDIMS = -1
      IF (SEC(1:5) == 'TUBE2') THEN
         GET_SECTION_NDIMS = 2
      ELSE IF (SEC(1:4) == 'ROD ') THEN
         GET_SECTION_NDIMS = 1
      ELSE IF (SEC(1:4) == 'TUBE') THEN
         GET_SECTION_NDIMS = 2
      ELSE IF (SEC(1:4) == 'BAR ') THEN
         GET_SECTION_NDIMS = 2
      ELSE IF (SEC(1:4) == 'BOX ') THEN
         GET_SECTION_NDIMS = 4
      ELSE IF (SEC(1:4) == 'H   ') THEN
         GET_SECTION_NDIMS = 4
      ELSE IF (SEC(1:4) == 'CHAN') THEN
         GET_SECTION_NDIMS = 4
      ELSE IF (SEC(1:4) == 'T   ') THEN
         GET_SECTION_NDIMS = 4
      ELSE IF (SEC(1:4) == 'L   ') THEN
         GET_SECTION_NDIMS = 4
      ELSE IF (SEC(1:4) == 'I   ') THEN
         GET_SECTION_NDIMS = 6
      ENDIF

      END FUNCTION GET_SECTION_NDIMS

! ##################################################################################################################################

      SUBROUTINE LOAD_SECTION_A ( SEC_TYPE_IN, NDIM_IN, DIMS, NSM )

      CHARACTER(LEN=*), INTENT(IN) :: SEC_TYPE_IN
      INTEGER(LONG), INTENT(IN)    :: NDIM_IN
      REAL(DOUBLE), INTENT(IN)     :: DIMS(10), NSM
      REAL(DOUBLE)             :: AREA, I1, I2, I12, JTOR
      REAL(DOUBLE)             :: STRE(8)
      REAL(DOUBLE)             :: K1, K2, YC, ZC, YS, ZS, IWARP
      INTEGER(LONG)            :: IERR_LOC

      CALL CALC_SECTION ( SEC_TYPE_IN, NDIM_IN, DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )
      RPBEAM(NPBEAM, 1) = AREA
      RPBEAM(NPBEAM, 2) = I1
      RPBEAM(NPBEAM, 3) = I2
      RPBEAM(NPBEAM, 4) = I12
      RPBEAM(NPBEAM, 5) = JTOR
      RPBEAM(NPBEAM, 6) = NSM
      DO I=1,8
         RPBEAM(NPBEAM,6+I) = STRE(I)
      ENDDO
      RPBEAM(NPBEAM,30) = K1
      RPBEAM(NPBEAM,31) = K2
      RPBEAM(NPBEAM,36) = IWARP
      RPBEAM(NPBEAM,38) = YS
      RPBEAM(NPBEAM,39) = ZS
      RPBEAM(NPBEAM,42) = YC
      RPBEAM(NPBEAM,43) = ZC
      CALL STORE_SECTION_PROPS ( 1, SEC_TYPE_IN, NDIM_IN, DIMS, NSM )
      CALL CHECK_BAR_MOIs ( 'PBEAML', ID, RPBEAM(NPBEAM,2), RPBEAM(NPBEAM,3), RPBEAM(NPBEAM,4), IERR_LOC )
      IF (IERR_LOC /= 0) FATAL_ERR = FATAL_ERR + 1

      END SUBROUTINE LOAD_SECTION_A

! ##################################################################################################################################

      SUBROUTINE LOAD_SECTION_B ( SEC_TYPE_IN, NDIM_IN, DIMS, NSM )

      CHARACTER(LEN=*), INTENT(IN) :: SEC_TYPE_IN
      INTEGER(LONG), INTENT(IN)    :: NDIM_IN
      REAL(DOUBLE), INTENT(IN)     :: DIMS(10), NSM
      REAL(DOUBLE)             :: AREA, I1, I2, I12, JTOR
      REAL(DOUBLE)             :: STRE(8)
      REAL(DOUBLE)             :: K1, K2, YC, ZC, YS, ZS, IWARP
      INTEGER(LONG)            :: IERR_LOC

      CALL CALC_SECTION ( SEC_TYPE_IN, NDIM_IN, DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )
      RPBEAM(NPBEAM,15) = PBEAM_XL(NPBEAM,PBEAM_NSTATIONS(NPBEAM))
      RPBEAM(NPBEAM,16) = AREA
      RPBEAM(NPBEAM,17) = I1
      RPBEAM(NPBEAM,18) = I2
      RPBEAM(NPBEAM,19) = I12
      RPBEAM(NPBEAM,20) = JTOR
      RPBEAM(NPBEAM,21) = NSM
      DO I=1,8
         RPBEAM(NPBEAM,21+I) = STRE(I)
      ENDDO
      RPBEAM(NPBEAM,30) = K1
      RPBEAM(NPBEAM,31) = K2
      RPBEAM(NPBEAM,37) = IWARP
      RPBEAM(NPBEAM,40) = YS
      RPBEAM(NPBEAM,41) = ZS
      RPBEAM(NPBEAM,44) = YC
      RPBEAM(NPBEAM,45) = ZC
      CALL STORE_SECTION_PROPS ( PBEAM_NSTATIONS(NPBEAM), SEC_TYPE_IN, NDIM_IN, DIMS, NSM )
      CALL CHECK_BAR_MOIs ( 'PBEAML', ID, RPBEAM(NPBEAM,17), RPBEAM(NPBEAM,18), RPBEAM(NPBEAM,19), IERR_LOC )
      IF (IERR_LOC /= 0) FATAL_ERR = FATAL_ERR + 1

      END SUBROUTINE LOAD_SECTION_B

! ##################################################################################################################################

      SUBROUTINE STORE_SECTION_PROPS ( ISTA_IN, SEC_TYPE_IN, NDIM_IN, DIMS, NSM )

      INTEGER(LONG), INTENT(IN)    :: ISTA_IN
      INTEGER(LONG), INTENT(IN)    :: NDIM_IN
      CHARACTER(LEN=*), INTENT(IN) :: SEC_TYPE_IN
      REAL(DOUBLE), INTENT(IN)     :: DIMS(10), NSM
      REAL(DOUBLE)                 :: AREA, I1, I2, I12, JTOR
      REAL(DOUBLE)                 :: K1, K2, YC, ZC, YS, ZS, IWARP, STRE(8)

      IF ((ISTA_IN < 1) .OR. (ISTA_IN > MPBEAM_STATIONS)) RETURN

      CALL CALC_SECTION ( SEC_TYPE_IN, NDIM_IN, DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )
      PBEAM_RPROPS(NPBEAM,ISTA_IN,1) = AREA
      PBEAM_RPROPS(NPBEAM,ISTA_IN,2) = I1
      PBEAM_RPROPS(NPBEAM,ISTA_IN,3) = I2
      PBEAM_RPROPS(NPBEAM,ISTA_IN,4) = I12
      PBEAM_RPROPS(NPBEAM,ISTA_IN,5) = JTOR
      PBEAM_RPROPS(NPBEAM,ISTA_IN,6) = NSM

! Default runtime modifiers for all beam properties. PBEAMZ can overwrite them later.
      RPBEAM(NPBEAM,46) = ONE
      RPBEAM(NPBEAM,47) = ZERO
      RPBEAM(NPBEAM,48) = ZERO
      RPBEAM(NPBEAM,49) = ONE
      RPBEAM(NPBEAM,50) = ONE
      RPBEAM(NPBEAM,51) = ONE
      RPBEAM(NPBEAM,52) = ONE
      RPBEAM(NPBEAM,53) = ONE
      RPBEAM(NPBEAM,54) = ONE

      CALL WRITE_PBEAML_CONVERTED_PBEAM_DEBUG ( 'S', PBEAM_XL(NPBEAM,ISTA_IN), SEC_TYPE_IN, NDIM_IN, DIMS, AREA, I1, I2, I12,   &
                                                JTOR, NSM, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )

      END SUBROUTINE STORE_SECTION_PROPS

! ##################################################################################################################################

      SUBROUTINE STORE_PBEAMZ_META ( STIFFMOD, ROFSET, TAPER_MODE, AREA_MOD, I1_MOD, I2_MOD, K1_MOD, K2_MOD, J_MOD )

      REAL(DOUBLE), INTENT(IN)     :: STIFFMOD
      REAL(DOUBLE), INTENT(IN)     :: ROFSET
      INTEGER(LONG), INTENT(IN)    :: TAPER_MODE
      REAL(DOUBLE), INTENT(IN)     :: AREA_MOD
      REAL(DOUBLE), INTENT(IN)     :: I1_MOD
      REAL(DOUBLE), INTENT(IN)     :: I2_MOD
      REAL(DOUBLE), INTENT(IN)     :: K1_MOD
      REAL(DOUBLE), INTENT(IN)     :: K2_MOD
      REAL(DOUBLE), INTENT(IN)     :: J_MOD

      IF (.NOT. IS_PBEAMZ) RETURN

      ! Do not use STIFFMOD as a global scalar multiplier for PBEAMZ.
      ! Keep the legacy column neutral; per-component modifiers live in cols 49-54.
      RPBEAM(NPBEAM,46) = ONE
      RPBEAM(NPBEAM,47) = ROFSET
      RPBEAM(NPBEAM,48) = DBLE(TAPER_MODE)
      RPBEAM(NPBEAM,49) = AREA_MOD
      RPBEAM(NPBEAM,50) = I1_MOD
      RPBEAM(NPBEAM,51) = I2_MOD
      RPBEAM(NPBEAM,52) = K1_MOD
      RPBEAM(NPBEAM,53) = K2_MOD
      RPBEAM(NPBEAM,54) = J_MOD

      CALL STORE_PBEAMZ_BEAM_PARAM_OVERRIDE ( 'BEAMAMO', BEAMAMO, BEAMAMO_PID, BEAMAMO_VAL, NBEAMAMO_PID, MBEAMAMO_PID, &
                                              INT(PBEAM(NPBEAM,1), KIND=LONG), AREA_MOD )
      CALL STORE_PBEAMZ_BEAM_PARAM_OVERRIDE ( 'BEAMV1MO', BEAMV1MO, BEAMV1MO_PID, BEAMV1MO_VAL, NBEAMV1MO_PID, MBEAMV1MO_PID, &
                                              INT(PBEAM(NPBEAM,1), KIND=LONG), K1_MOD )
      CALL STORE_PBEAMZ_BEAM_PARAM_OVERRIDE ( 'BEAMV2MO', BEAMV2MO, BEAMV2MO_PID, BEAMV2MO_VAL, NBEAMV2MO_PID, MBEAMV2MO_PID, &
                                              INT(PBEAM(NPBEAM,1), KIND=LONG), K2_MOD )
      CALL STORE_PBEAMZ_BEAM_PARAM_OVERRIDE ( 'BEAMM1MO', BEAMM1MO, BEAMM1MO_PID, BEAMM1MO_VAL, NBEAMM1MO_PID, MBEAMM1MO_PID, &
                                              INT(PBEAM(NPBEAM,1), KIND=LONG), I1_MOD )
      CALL STORE_PBEAMZ_BEAM_PARAM_OVERRIDE ( 'BEAMM2MO', BEAMM2MO, BEAMM2MO_PID, BEAMM2MO_VAL, NBEAMM2MO_PID, MBEAMM2MO_PID, &
                                              INT(PBEAM(NPBEAM,1), KIND=LONG), I2_MOD )
      CALL STORE_PBEAMZ_BEAM_PARAM_OVERRIDE ( 'BEAMTMO', BEAMTMO, BEAMTMO_PID, BEAMTMO_VAL, NBEAMTMO_PID, MBEAMTMO_PID, &
                                              INT(PBEAM(NPBEAM,1), KIND=LONG), J_MOD )

      END SUBROUTINE STORE_PBEAMZ_META

! ##################################################################################################################################

      SUBROUTINE STORE_PBEAMZ_BEAM_PARAM_OVERRIDE ( PARNAM, GLOBAL_VAL, PID_ARR, VAL_ARR, NPID, MPID, PID_IN, VALUE_IN )

      CHARACTER(LEN=*), INTENT(IN)    :: PARNAM
      REAL(DOUBLE), INTENT(INOUT)     :: GLOBAL_VAL
      INTEGER(LONG), INTENT(INOUT)    :: PID_ARR(:)
      REAL(DOUBLE), INTENT(INOUT)     :: VAL_ARR(:)
      INTEGER(LONG), INTENT(INOUT)    :: NPID
      INTEGER(LONG), INTENT(IN)       :: MPID
      INTEGER(LONG), INTENT(IN)       :: PID_IN
      REAL(DOUBLE), INTENT(IN)        :: VALUE_IN
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: ISLOT

      ISLOT = 0
      DO I=1,NPID
         IF (PID_ARR(I) == PID_IN) THEN
            ISLOT = I
            EXIT
         ENDIF
      ENDDO

      IF (ISLOT == 0) THEN
         IF (NPID < MIN(MPID,SIZE(PID_ARR))) THEN
            NPID = NPID + 1
            ISLOT = NPID
         ENDIF
      ENDIF

      IF (ISLOT > 0) THEN
         PID_ARR(ISLOT) = PID_IN
         VAL_ARR(ISLOT) = VALUE_IN
      ENDIF

      END SUBROUTINE STORE_PBEAMZ_BEAM_PARAM_OVERRIDE

! ##################################################################################################################################

      SUBROUTINE FINALIZE_PBEAMZ_GENERATED_PBEAM ( TAPER_MODE, NSEG_STATIONS, AREA_MOD, I1_MOD, I2_MOD, K1_MOD, K2_MOD, J_MOD, ROFSET, NSTATION_EXTRA, STATION_EXTRA )

      INTEGER(LONG), INTENT(IN)    :: TAPER_MODE
      INTEGER(LONG), INTENT(IN)    :: NSEG_STATIONS
      REAL(DOUBLE), INTENT(IN)     :: AREA_MOD
      REAL(DOUBLE), INTENT(IN)     :: I1_MOD
      REAL(DOUBLE), INTENT(IN)     :: I2_MOD
      REAL(DOUBLE), INTENT(IN)     :: K1_MOD
      REAL(DOUBLE), INTENT(IN)     :: K2_MOD
      REAL(DOUBLE), INTENT(IN)     :: J_MOD
      REAL(DOUBLE), INTENT(IN)     :: ROFSET
      INTEGER(LONG), INTENT(IN)    :: NSTATION_EXTRA
      REAL(DOUBLE), INTENT(IN)     :: STATION_EXTRA(3)

      IF (.NOT. IS_PBEAMZ) RETURN

      IF (PBEAM_NSTATIONS(NPBEAM) <= 2) THEN
         CALL EXPAND_PBEAMZ_DEFAULT_STATIONS ( TAPER_MODE, NSEG_STATIONS )
      ENDIF

      IF (NSTATION_EXTRA > 0) THEN
         CALL INSERT_PBEAMZ_EXTRA_STATIONS ( NSTATION_EXTRA, STATION_EXTRA )
      ENDIF

      CALL REFRESH_PBEAMZ_RPBEAM_SNAPSHOTS ()

      END SUBROUTINE FINALIZE_PBEAMZ_GENERATED_PBEAM

! ##################################################################################################################################

      SUBROUTINE REFRESH_PBEAMZ_RPBEAM_SNAPSHOTS ()

      INTEGER(LONG)                :: NSTA_LOC

      NSTA_LOC = PBEAM_NSTATIONS(NPBEAM)
      IF (NSTA_LOC < 1) RETURN

      RPBEAM(NPBEAM, 1) = PBEAM_RPROPS(NPBEAM,1,1)
      RPBEAM(NPBEAM, 2) = PBEAM_RPROPS(NPBEAM,1,2)
      RPBEAM(NPBEAM, 3) = PBEAM_RPROPS(NPBEAM,1,3)
      RPBEAM(NPBEAM, 4) = PBEAM_RPROPS(NPBEAM,1,4)
      RPBEAM(NPBEAM, 5) = PBEAM_RPROPS(NPBEAM,1,5)
      RPBEAM(NPBEAM, 6) = PBEAM_RPROPS(NPBEAM,1,6)
      RPBEAM(NPBEAM,15) = PBEAM_XL(NPBEAM,NSTA_LOC)
      RPBEAM(NPBEAM,16) = PBEAM_RPROPS(NPBEAM,NSTA_LOC,1)
      RPBEAM(NPBEAM,17) = PBEAM_RPROPS(NPBEAM,NSTA_LOC,2)
      RPBEAM(NPBEAM,18) = PBEAM_RPROPS(NPBEAM,NSTA_LOC,3)
      RPBEAM(NPBEAM,19) = PBEAM_RPROPS(NPBEAM,NSTA_LOC,4)
      RPBEAM(NPBEAM,20) = PBEAM_RPROPS(NPBEAM,NSTA_LOC,5)
      RPBEAM(NPBEAM,21) = PBEAM_RPROPS(NPBEAM,NSTA_LOC,6)

      END SUBROUTINE REFRESH_PBEAMZ_RPBEAM_SNAPSHOTS

! ##################################################################################################################################

      SUBROUTINE EXPAND_PBEAMZ_DEFAULT_STATIONS ( TAPER_MODE, NSEG_DEFAULT )

      INTEGER(LONG), INTENT(IN)    :: TAPER_MODE
      INTEGER(LONG), INTENT(IN)    :: NSEG_DEFAULT
      INTEGER(LONG)                :: ISTA
      INTEGER(LONG)                :: TAPER_EXP
      REAL(DOUBLE)                 :: FRAC
      REAL(DOUBLE)                 :: A1, A2, I1A, I1B, I2A, I2B, I12A, I12B, J1, J2, NSM1, NSM2
      REAL(DOUBLE)                 :: AREA_OUT, I1_OUT, I2_OUT, I12_OUT, JTOR_OUT, NSM_OUT

      IF (.NOT. IS_PBEAMZ) RETURN
      IF (PBEAM_NSTATIONS(NPBEAM) > 2) RETURN
      IF (NSEG_DEFAULT <= 0) RETURN

      A1   = PBEAM_RPROPS(NPBEAM,1,1)
      I1A  = PBEAM_RPROPS(NPBEAM,1,2)
      I2A  = PBEAM_RPROPS(NPBEAM,1,3)
      I12A = PBEAM_RPROPS(NPBEAM,1,4)
      J1   = PBEAM_RPROPS(NPBEAM,1,5)
      NSM1 = PBEAM_RPROPS(NPBEAM,1,6)
      A2   = PBEAM_RPROPS(NPBEAM,2,1)
      I1B  = PBEAM_RPROPS(NPBEAM,2,2)
      I2B  = PBEAM_RPROPS(NPBEAM,2,3)
      I12B = PBEAM_RPROPS(NPBEAM,2,4)
      J2   = PBEAM_RPROPS(NPBEAM,2,5)
      NSM2 = PBEAM_RPROPS(NPBEAM,2,6)

      PBEAM_NSTATIONS(NPBEAM) = NSEG_DEFAULT + 1
      TAPER_EXP = MAX(0, MIN(3, TAPER_MODE))

      DO ISTA=1,NSEG_DEFAULT+1
         FRAC = DBLE(ISTA-1)/DBLE(NSEG_DEFAULT)
         PBEAM_XL(NPBEAM,ISTA) = FRAC
         CALL INTERP_PBEAMZ_STATION ( FRAC, TAPER_EXP, A1, A2, I1A, I1B, I2A, I2B, I12A, I12B, J1, J2, NSM1, NSM2,   &
                                      AREA_OUT, I1_OUT, I2_OUT, I12_OUT, JTOR_OUT, NSM_OUT )
         PBEAM_RPROPS(NPBEAM,ISTA,1) = AREA_OUT
         PBEAM_RPROPS(NPBEAM,ISTA,2) = I1_OUT
         PBEAM_RPROPS(NPBEAM,ISTA,3) = I2_OUT
         PBEAM_RPROPS(NPBEAM,ISTA,4) = I12_OUT
         PBEAM_RPROPS(NPBEAM,ISTA,5) = JTOR_OUT
         PBEAM_RPROPS(NPBEAM,ISTA,6) = NSM_OUT
      ENDDO

      PBEAM_XL(NPBEAM,1) = ZERO
      PBEAM_XL(NPBEAM,NSEG_DEFAULT+1) = ONE
      RPBEAM(NPBEAM,15) = ONE
      RPBEAM(NPBEAM,16) = PBEAM_RPROPS(NPBEAM,NSEG_DEFAULT+1,1)
      RPBEAM(NPBEAM,17) = PBEAM_RPROPS(NPBEAM,NSEG_DEFAULT+1,2)
      RPBEAM(NPBEAM,18) = PBEAM_RPROPS(NPBEAM,NSEG_DEFAULT+1,3)
      RPBEAM(NPBEAM,19) = PBEAM_RPROPS(NPBEAM,NSEG_DEFAULT+1,4)
      RPBEAM(NPBEAM,20) = PBEAM_RPROPS(NPBEAM,NSEG_DEFAULT+1,5)
      RPBEAM(NPBEAM,21) = PBEAM_RPROPS(NPBEAM,NSEG_DEFAULT+1,6)

      END SUBROUTINE EXPAND_PBEAMZ_DEFAULT_STATIONS

! ##################################################################################################################################

      SUBROUTINE INSERT_PBEAMZ_EXTRA_STATIONS ( NSTATION_EXTRA, STATION_EXTRA )

      INTEGER(LONG), INTENT(IN)    :: NSTATION_EXTRA
      REAL(DOUBLE), INTENT(IN)     :: STATION_EXTRA(3)

      INTEGER(LONG), PARAMETER     :: NEXTRA_MAX = 3
      INTEGER(LONG)                :: I, J, K, NSTA_OLD, NSTA_NEW, INSERT_POS, NEXTRA, IDX(3)
      REAL(DOUBLE)                 :: FRAC, EPSF, F1, F2, W, V1, V2
      REAL(DOUBLE)                 :: ROW(6)

      IF (.NOT. IS_PBEAMZ) RETURN
      NEXTRA = MAX(0, MIN(NEXTRA_MAX, NSTATION_EXTRA))
      IF (NEXTRA <= 0) RETURN

      IDX = (/1, 2, 3/)
      DO I=1,NEXTRA-1
         DO J=I+1,NEXTRA
            IF (STATION_EXTRA(IDX(J)) < STATION_EXTRA(IDX(I))) THEN
               K = IDX(I)
               IDX(I) = IDX(J)
               IDX(J) = K
            ENDIF
         ENDDO
      ENDDO

      EPSF = 1.0D-10
      DO I=1,NEXTRA
         FRAC = STATION_EXTRA(IDX(I))
         IF ((FRAC <= ZERO) .OR. (FRAC >= ONE)) CYCLE

         NSTA_OLD = PBEAM_NSTATIONS(NPBEAM)
         INSERT_POS = 1
         DO WHILE ((INSERT_POS <= NSTA_OLD) .AND. (PBEAM_XL(NPBEAM,INSERT_POS) < FRAC))
            INSERT_POS = INSERT_POS + 1
         ENDDO

         IF ((INSERT_POS <= NSTA_OLD) .AND. (ABS(PBEAM_XL(NPBEAM,INSERT_POS) - FRAC) <= EPSF)) CYCLE

         IF (NSTA_OLD >= MPBEAM_STATIONS) THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,'(A,A,A)') ' *WARNING 1197: PBEAMZ ID = ', ID, ' cannot add extra stations because the station buffer is full'
            WRITE(F06,'(A,A,A)') ' *WARNING 1197: PBEAMZ ID = ', ID, ' cannot add extra stations because the station buffer is full'
            RETURN
         ENDIF

         IF (INSERT_POS > NSTA_OLD) CYCLE

         F1 = PBEAM_XL(NPBEAM,MAX(1,INSERT_POS-1))
         F2 = PBEAM_XL(NPBEAM,INSERT_POS)
         IF (ABS(F2 - F1) <= EPSF) CYCLE
         W = (FRAC - F1)/(F2 - F1)

         DO J=NSTA_OLD,INSERT_POS,-1
            PBEAM_XL(NPBEAM,J+1) = PBEAM_XL(NPBEAM,J)
            DO K=1,6
               PBEAM_RPROPS(NPBEAM,J+1,K) = PBEAM_RPROPS(NPBEAM,J,K)
            ENDDO
         ENDDO

         PBEAM_XL(NPBEAM,INSERT_POS) = FRAC
         DO K=1,6
            V1 = PBEAM_RPROPS(NPBEAM,INSERT_POS-1,K)
            V2 = PBEAM_RPROPS(NPBEAM,INSERT_POS+1,K)
            ROW(K) = (ONE - W)*V1 + W*V2
         ENDDO
         DO K=1,6
            PBEAM_RPROPS(NPBEAM,INSERT_POS,K) = ROW(K)
         ENDDO
         PBEAM_NSTATIONS(NPBEAM) = NSTA_OLD + 1
      ENDDO

      END SUBROUTINE INSERT_PBEAMZ_EXTRA_STATIONS

! ##################################################################################################################################

      SUBROUTINE INTERP_PBEAMZ_STATION ( FRAC, TAPER_MODE, A1, A2, I1A, I1B, I2A, I2B, I12A, I12B, J1, J2, NSM1, NSM2, &
                                         AREA_OUT, I1_OUT, I2_OUT, I12_OUT, JTOR_OUT, NSM_OUT )

      REAL(DOUBLE), INTENT(IN)     :: FRAC
      INTEGER(LONG), INTENT(IN)    :: TAPER_MODE
      REAL(DOUBLE), INTENT(IN)     :: A1, A2, I1A, I1B, I2A, I2B, I12A, I12B, J1, J2, NSM1, NSM2
      REAL(DOUBLE), INTENT(OUT)    :: AREA_OUT, I1_OUT, I2_OUT, I12_OUT, JTOR_OUT, NSM_OUT

      INTEGER(LONG)                :: NEXP
      REAL(DOUBLE)                 :: INVN

      IF (TAPER_MODE <= 0) THEN
         AREA_OUT = A1
         I1_OUT   = I1A
         I2_OUT   = I2A
         I12_OUT  = I12A
         JTOR_OUT = J1
         NSM_OUT  = NSM1
         RETURN
      ENDIF

      NEXP = MAX(1, MIN(3, TAPER_MODE))
      INVN = ONE/DBLE(NEXP)

      AREA_OUT = (ONE - FRAC)*A1 + FRAC*A2
      IF ((I1A > ZERO) .AND. (I1B > ZERO)) THEN
         I1_OUT = ((I1A**INVN)*(ONE - FRAC) + (I1B**INVN)*FRAC)**DBLE(NEXP)
      ELSE
         I1_OUT = (ONE - FRAC)*I1A + FRAC*I1B
      ENDIF
      IF ((I2A > ZERO) .AND. (I2B > ZERO)) THEN
         I2_OUT = ((I2A**INVN)*(ONE - FRAC) + (I2B**INVN)*FRAC)**DBLE(NEXP)
      ELSE
         I2_OUT = (ONE - FRAC)*I2A + FRAC*I2B
      ENDIF
      ! Keep torsion linear by default; the taper mode only governs section inertia interpolation.
      JTOR_OUT = (ONE - FRAC)*J1 + FRAC*J2
      I12_OUT = (ONE - FRAC)*I12A + FRAC*I12B
      NSM_OUT  = (ONE - FRAC)*NSM1  + FRAC*NSM2

      END SUBROUTINE INTERP_PBEAMZ_STATION

! ##################################################################################################################################

      LOGICAL FUNCTION IS_PBEAMZ_OPTION ( TOKEN )

      CHARACTER(LEN=*), INTENT(IN) :: TOKEN
      CHARACTER(LEN=JCARD_LEN)     :: TOKEN_UP

      TOKEN_UP = TOKEN
      CALL TO_UPPER ( TOKEN_UP )

      IF ((TRIM(TOKEN_UP) == 'STIFFMOD') .OR. (TRIM(TOKEN_UP) == 'ROFSET') .OR. (TRIM(TOKEN_UP) == 'RIOFFSET') .OR. (TRIM(TOKEN_UP) == 'TAPER') .OR. &
          (TRIM(TOKEN_UP) == 'NSM') .OR. (TRIM(TOKEN_UP) == 'END') .OR. &
          (TRIM(TOKEN_UP) == 'AREAMOD') .OR. (TRIM(TOKEN_UP) == 'I1MOD') .OR. (TRIM(TOKEN_UP) == 'I2MOD') .OR. (TRIM(TOKEN_UP) == 'K1MOD') .OR. &
          (TRIM(TOKEN_UP) == 'K2MOD') .OR. (TRIM(TOKEN_UP) == 'JMOD') .OR. (TRIM(TOKEN_UP) == 'STATIONS')) THEN
         IS_PBEAMZ_OPTION = .TRUE.
      ELSE
         IS_PBEAMZ_OPTION = .FALSE.
      ENDIF

      END FUNCTION IS_PBEAMZ_OPTION

! ##################################################################################################################################

      LOGICAL FUNCTION IS_PBEAMZ_DIM_LABEL ( TOKEN )

      CHARACTER(LEN=*), INTENT(IN) :: TOKEN
      CHARACTER(LEN=JCARD_LEN)     :: TOKEN_UP

      TOKEN_UP = TOKEN
      CALL TO_UPPER ( TOKEN_UP )

      IF (LEN_TRIM(TOKEN_UP) == 5) THEN
         IF ((TOKEN_UP(1:3) == 'DIM') .AND. ((TOKEN_UP(4:4) == '0') .OR. (TOKEN_UP(4:4) == '1'))) THEN
            IF (((TOKEN_UP(5:5) >= 'A') .AND. (TOKEN_UP(5:5) <= 'Z'))) THEN
               IS_PBEAMZ_DIM_LABEL = .TRUE.
               RETURN
            ENDIF
         ENDIF
      ENDIF

      IS_PBEAMZ_DIM_LABEL = .FALSE.

      END FUNCTION IS_PBEAMZ_DIM_LABEL

! ##################################################################################################################################

      SUBROUTINE SET_PBEAMZ_DIM_LABEL_FLAGS ( TOKEN, DIM0A_SEEN, DIM1A_SEEN )

      CHARACTER(LEN=*), INTENT(IN) :: TOKEN
      LOGICAL, INTENT(INOUT)       :: DIM0A_SEEN
      LOGICAL, INTENT(INOUT)       :: DIM1A_SEEN
      CHARACTER(LEN=JCARD_LEN)     :: TOKEN_UP

      TOKEN_UP = TOKEN
      CALL TO_UPPER ( TOKEN_UP )
      IF (TOKEN_UP(1:5) == 'DIM0A') THEN
         DIM0A_SEEN = .TRUE.
      ELSE IF (TOKEN_UP(1:5) == 'DIM1A') THEN
         DIM1A_SEEN = .TRUE.
      ENDIF

      END SUBROUTINE SET_PBEAMZ_DIM_LABEL_FLAGS

! ##################################################################################################################################

      SUBROUTINE PARSE_PBEAMZ_TOKENS ( NDIM_SEC, NTOK, TOKENS, DIMS_A, DIMS_B, NSM_A, NSM_B, DIMS_CUR, STIFFMOD, ROFSET, TAPER_MODE, &
                                       AREA_MOD, I1_MOD, I2_MOD, K1_MOD, K2_MOD, J_MOD, NSEG_STATIONS, NSTATION_EXTRA, STATION_EXTRA )

      INTEGER(LONG), INTENT(IN)       :: NDIM_SEC
      INTEGER(LONG), INTENT(IN)       :: NTOK
      CHARACTER(LEN=JCARD_LEN), INTENT(IN) :: TOKENS(MAXTOK)
      REAL(DOUBLE), INTENT(INOUT)     :: DIMS_A(10)
      REAL(DOUBLE), INTENT(INOUT)     :: DIMS_B(10)
      REAL(DOUBLE), INTENT(INOUT)     :: DIMS_CUR(10)
      REAL(DOUBLE), INTENT(INOUT)     :: NSM_A
      REAL(DOUBLE), INTENT(INOUT)     :: NSM_B
      REAL(DOUBLE), INTENT(INOUT)     :: STIFFMOD
      REAL(DOUBLE), INTENT(INOUT)     :: ROFSET
      INTEGER(LONG), INTENT(INOUT)    :: TAPER_MODE
      REAL(DOUBLE), INTENT(INOUT)     :: AREA_MOD
      REAL(DOUBLE), INTENT(INOUT)     :: I1_MOD
      REAL(DOUBLE), INTENT(INOUT)     :: I2_MOD
      REAL(DOUBLE), INTENT(INOUT)     :: K1_MOD
      REAL(DOUBLE), INTENT(INOUT)     :: K2_MOD
      REAL(DOUBLE), INTENT(INOUT)     :: J_MOD
      INTEGER(LONG), INTENT(INOUT)    :: NSEG_STATIONS
      INTEGER(LONG), INTENT(INOUT)    :: NSTATION_EXTRA
      REAL(DOUBLE), INTENT(INOUT)     :: STATION_EXTRA(3)

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: ITOK_LOC
      LOGICAL                         :: HAVE_END_SECTION

      IF (NTOK < NDIM_SEC) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR  1305: PBEAMZ ID = ', ID
         WRITE(F06,'(A,A)') ' *ERROR  1305: PBEAMZ ID = ', ID
         RETURN
      ENDIF

      DIMS_A = ZERO
      DIMS_B = ZERO
      DIMS_CUR = ZERO
      STATION_EXTRA = ZERO
      NSM_A = ZERO
      NSM_B = ZERO
      STIFFMOD   = ONE
      ROFSET     = ZERO
      TAPER_MODE = 0
      AREA_MOD   = ONE
      I1_MOD     = ONE
      I2_MOD     = ONE
      K1_MOD     = ONE
      K2_MOD     = ONE
      J_MOD      = ONE
      NSEG_STATIONS = 10
      NSTATION_EXTRA = 0
      HAVE_END_SECTION = .FALSE.

      DO I=1,NDIM_SEC
         READ(TOKENS(I),*,ERR=900) DIMS_A(I)
      ENDDO
      ITOK_LOC = NDIM_SEC + 1
      IF (ITOK_LOC <= NTOK) THEN
         IF (.NOT. IS_PBEAMZ_OPTION(TOKENS(ITOK_LOC))) THEN
            READ(TOKENS(ITOK_LOC),*,ERR=900) NSM_A
            ITOK_LOC = ITOK_LOC + 1
         ENDIF
      ENDIF

      CALL LOAD_SECTION_A ( SEC_TYPE, NDIM_SEC, DIMS_A, NSM_A )
      DIMS_B = DIMS_A
      NSM_B  = NSM_A

      DO WHILE (ITOK_LOC <= NTOK)
         IF (IS_PBEAMZ_OPTION(TOKENS(ITOK_LOC))) THEN
            CALL READ_PBEAMZ_OPTION ( TOKENS(ITOK_LOC), ITOK_LOC, STIFFMOD, ROFSET, TAPER_MODE, AREA_MOD, I1_MOD, I2_MOD, K1_MOD, K2_MOD, &
                                      J_MOD, NSM_A, NSM_B, NSEG_STATIONS, NSTATION_EXTRA, STATION_EXTRA )
            CYCLE
         ENDIF

         IF ((TAPER_MODE > 0) .AND. (.NOT. HAVE_END_SECTION)) THEN
            IF (ITOK_LOC + NDIM_SEC - 1 > NTOK) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,'(A,A)') ' *ERROR  1307: PBEAMZ ID = ', ID
               WRITE(F06,'(A,A)') ' *ERROR  1307: PBEAMZ ID = ', ID
               RETURN
            ENDIF
            DO I=1,NDIM_SEC
               READ(TOKENS(ITOK_LOC),*,ERR=900) DIMS_B(I)
               ITOK_LOC = ITOK_LOC + 1
            ENDDO
            NSM_B = NSM_A
            IF (ITOK_LOC <= NTOK) THEN
               IF (.NOT. IS_PBEAMZ_OPTION(TOKENS(ITOK_LOC))) THEN
                  READ(TOKENS(ITOK_LOC),*,ERR=900) NSM_B
                  ITOK_LOC = ITOK_LOC + 1
               ENDIF
            ENDIF
            HAVE_END_SECTION = .TRUE.
            CYCLE
         ENDIF

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A,A)') ' *ERROR  1317: PBEAMZ ID = ', ID, ' HAS UNEXPECTED RAW TOKEN OUTSIDE DIM0A/DIM1A OR OPTION BLOCKS'
         WRITE(F06,'(A,A,A)') ' *ERROR  1317: PBEAMZ ID = ', ID, ' HAS UNEXPECTED RAW TOKEN OUTSIDE DIM0A/DIM1A OR OPTION BLOCKS'
         RETURN
      ENDDO

      IF ((TAPER_MODE > 0) .AND. (.NOT. HAVE_END_SECTION)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR  1318: PBEAMZ ENTRY DEFINES TAPER BUT DOES NOT SUPPLY A DIM1A SECTION FOR ID = ', ID
         WRITE(F06,'(A,A)') ' *ERROR  1318: PBEAMZ ENTRY DEFINES TAPER BUT DOES NOT SUPPLY A DIM1A SECTION FOR ID = ', ID
         RETURN
      ENDIF

      PBEAM_NSTATIONS(NPBEAM) = 2
      PBEAM_XL(NPBEAM,1) = ZERO
      PBEAM_XL(NPBEAM,2) = ONE
      CALL LOAD_SECTION_B ( SEC_TYPE, NDIM_SEC, DIMS_B, NSM_B )
      CALL STORE_PBEAMZ_META ( STIFFMOD, ROFSET, TAPER_MODE, AREA_MOD, I1_MOD, I2_MOD, K1_MOD, K2_MOD, J_MOD )
      CALL FINALIZE_PBEAMZ_GENERATED_PBEAM ( TAPER_MODE, NSEG_STATIONS, AREA_MOD, I1_MOD, I2_MOD, K1_MOD, K2_MOD, J_MOD, ROFSET, &
                                             NSTATION_EXTRA, STATION_EXTRA )

      RETURN

  900 CONTINUE
      FATAL_ERR = FATAL_ERR + 1
      WRITE(ERR,'(A,A)') ' *ERROR  1309: PBEAMZ ID = ', ID
      WRITE(F06,'(A,A)') ' *ERROR  1309: PBEAMZ ID = ', ID
      RETURN

      END SUBROUTINE PARSE_PBEAMZ_TOKENS

! ##################################################################################################################################

      SUBROUTINE READ_PBEAMZ_OPTION ( OPTION_TOKEN, ITOK, STIFFMOD, ROFSET, TAPER_MODE, AREA_MOD, I1_MOD, I2_MOD, K1_MOD, K2_MOD, J_MOD, NSM_A, NSM_B, NSEG_STATIONS, NSTATION_EXTRA, STATION_EXTRA )

      CHARACTER(LEN=*), INTENT(IN)    :: OPTION_TOKEN
      INTEGER(LONG), INTENT(INOUT)    :: ITOK
      REAL(DOUBLE), INTENT(INOUT)     :: STIFFMOD
      REAL(DOUBLE), INTENT(INOUT)     :: ROFSET
      INTEGER(LONG), INTENT(INOUT)    :: TAPER_MODE
      REAL(DOUBLE), INTENT(INOUT)     :: AREA_MOD
      REAL(DOUBLE), INTENT(INOUT)     :: I1_MOD
      REAL(DOUBLE), INTENT(INOUT)     :: I2_MOD
      REAL(DOUBLE), INTENT(INOUT)     :: K1_MOD
      REAL(DOUBLE), INTENT(INOUT)     :: K2_MOD
      REAL(DOUBLE), INTENT(INOUT)     :: J_MOD
      REAL(DOUBLE), INTENT(INOUT)     :: NSM_A
      REAL(DOUBLE), INTENT(INOUT)     :: NSM_B
      INTEGER(LONG), INTENT(INOUT)    :: NSEG_STATIONS
      INTEGER(LONG), INTENT(INOUT)    :: NSTATION_EXTRA
      REAL(DOUBLE), INTENT(INOUT)     :: STATION_EXTRA(3)

      CHARACTER(LEN=JCARD_LEN)        :: OPTION_UP
      CHARACTER(LEN=JCARD_LEN)        :: VALUE_TOKEN
      CHARACTER(LEN=JCARD_LEN)        :: NEXT_TOKEN
      INTEGER(LONG)                   :: IOCHK
      INTEGER(LONG)                   :: JVAL
      INTEGER(LONG)                   :: NUM_COMPACT_MODS
      LOGICAL                         :: HIT_PBEAMZ_KEYWORD
      REAL(DOUBLE)                    :: MODVALS(6)
      REAL(DOUBLE)                    :: TMP_MODVALS(6)

      OPTION_UP = OPTION_TOKEN
      CALL TO_UPPER ( OPTION_UP )

      IF (TRIM(OPTION_UP) == 'END') THEN
         ITOK = NTOK + 1
         RETURN
      ENDIF

      IF (TRIM(OPTION_UP) == 'STIFFMOD') THEN
         IF (ITOK >= NTOK) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            WRITE(F06,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            RETURN
         ENDIF

         ! PBEAMZ compact syntax:
         !   STIFFMOD,area,i1,i2,k1,k2,j
         ! The numeric tail is interpreted as per-component modifiers in the
         ! order axial area, I1, I2, K1, K2, J. Missing trailing values default
         ! to 1.0.
         MODVALS(1) = ONE
         MODVALS(2) = ONE
         MODVALS(3) = ONE
         MODVALS(4) = ONE
         MODVALS(5) = ONE
         MODVALS(6) = ONE
         TMP_MODVALS = ONE
         JVAL = 1
         NUM_COMPACT_MODS = 0
         HIT_PBEAMZ_KEYWORD = .FALSE.
         DO WHILE ((ITOK < NTOK) .AND. (JVAL <= 6))
            ITOK = ITOK + 1
            NEXT_TOKEN = TOKENS(ITOK)
            CALL TO_UPPER ( NEXT_TOKEN )
            IF (IS_PBEAMZ_OPTION ( NEXT_TOKEN )) THEN
               HIT_PBEAMZ_KEYWORD = .TRUE.
               EXIT
            ENDIF
            READ(TOKENS(ITOK),*,IOSTAT=IOCHK) TMP_MODVALS(JVAL)
            IF (IOCHK /= 0) EXIT
            NUM_COMPACT_MODS = NUM_COMPACT_MODS + 1
            JVAL = JVAL + 1
         ENDDO

         IF (.NOT. HIT_PBEAMZ_KEYWORD) ITOK = NTOK + 1

         IF (NUM_COMPACT_MODS > 0) THEN
            AREA_MOD = TMP_MODVALS(1)
            I1_MOD   = TMP_MODVALS(2)
            I2_MOD   = TMP_MODVALS(3)
            K1_MOD   = TMP_MODVALS(4)
            K2_MOD   = TMP_MODVALS(5)
            J_MOD    = TMP_MODVALS(6)
         ELSE
            AREA_MOD = ONE
         ENDIF

      ELSE IF ((TRIM(OPTION_UP) == 'ROFSET') .OR. (TRIM(OPTION_UP) == 'RIOFFSET')) THEN
         IF (ITOK >= NTOK) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            WRITE(F06,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            RETURN
         ENDIF
         ITOK = ITOK + 1
         READ(TOKENS(ITOK),*,IOSTAT=IOCHK) ROFSET
         IF (IOCHK /= 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            WRITE(F06,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            RETURN
         ENDIF
         ITOK = ITOK + 1

      ELSE IF (TRIM(OPTION_UP) == 'NSM') THEN
         IF (ITOK >= NTOK) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            WRITE(F06,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            RETURN
         ENDIF
         ITOK = ITOK + 1
         READ(TOKENS(ITOK),*,IOSTAT=IOCHK) NSM_A
         IF (IOCHK /= 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            WRITE(F06,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            RETURN
         ENDIF
         NSM_B = NSM_A

         IF (ITOK < NTOK) THEN
            NEXT_TOKEN = TOKENS(ITOK+1)
            IF (.NOT. IS_PBEAMZ_OPTION ( NEXT_TOKEN )) THEN
               ITOK = ITOK + 1
               READ(TOKENS(ITOK),*,IOSTAT=IOCHK) NSM_B
               IF (IOCHK /= 0) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
                  WRITE(F06,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
                  RETURN
               ENDIF
            ENDIF
         ENDIF
         ITOK = ITOK + 1

      ELSE IF (TRIM(OPTION_UP) == 'STATIONS') THEN
         IF (ITOK >= NTOK) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            WRITE(F06,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            RETURN
         ENDIF
         ITOK = ITOK + 1
         READ(TOKENS(ITOK),*,IOSTAT=IOCHK) NSEG_STATIONS
         IF (IOCHK /= 0) THEN
            WARN_ERR = WARN_ERR + 1
            NSEG_STATIONS = 10
            WRITE(ERR,'(A,A,A,A,A)') ' *WARNING 1315: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID STATIONS COUNT "', TRIM(TOKENS(ITOK)), '" - DEFAULTING TO 10'
            WRITE(F06,'(A,A,A,A,A)') ' *WARNING 1315: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID STATIONS COUNT "', TRIM(TOKENS(ITOK)), '" - DEFAULTING TO 10'
         ELSE IF (NSEG_STATIONS <= 0) THEN
            WARN_ERR = WARN_ERR + 1
            NSEG_STATIONS = 10
            WRITE(ERR,'(A,A,A,I0)') ' *WARNING 1315: ', TRIM(CARD_NAME), ' ENTRY HAS STATIONS COUNT OUT OF RANGE, DEFAULTING TO 10. VALUE = ', NSEG_STATIONS
            WRITE(F06,'(A,A,A,I0)') ' *WARNING 1315: ', TRIM(CARD_NAME), ' ENTRY HAS STATIONS COUNT OUT OF RANGE, DEFAULTING TO 10. VALUE = ', NSEG_STATIONS
         ENDIF
         STATION_EXTRA = ZERO
         NSTATION_EXTRA = 0
         DO JVAL=1,3
            IF (ITOK >= NTOK) EXIT
            NEXT_TOKEN = TOKENS(ITOK+1)
            IF (IS_PBEAMZ_OPTION ( NEXT_TOKEN )) EXIT
            ITOK = ITOK + 1
            READ(TOKENS(ITOK),*,IOSTAT=IOCHK) STATION_EXTRA(JVAL)
            IF (IOCHK /= 0) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
               WRITE(F06,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
               RETURN
            ENDIF
            NSTATION_EXTRA = NSTATION_EXTRA + 1
         ENDDO
         ITOK = ITOK + 1

      ELSE IF (TRIM(OPTION_UP) == 'AREAMOD') THEN
         IF (ITOK >= NTOK) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            WRITE(F06,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            RETURN
         ENDIF
         ITOK = ITOK + 1
         READ(TOKENS(ITOK),*,IOSTAT=IOCHK) AREA_MOD
         IF (IOCHK /= 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            WRITE(F06,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            RETURN
         ENDIF
         ITOK = ITOK + 1

      ELSE IF (TRIM(OPTION_UP) == 'I1MOD') THEN
         IF (ITOK >= NTOK) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            WRITE(F06,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            RETURN
         ENDIF
         ITOK = ITOK + 1
         READ(TOKENS(ITOK),*,IOSTAT=IOCHK) I1_MOD
         IF (IOCHK /= 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            WRITE(F06,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            RETURN
         ENDIF
         ITOK = ITOK + 1

      ELSE IF (TRIM(OPTION_UP) == 'I2MOD') THEN
         IF (ITOK >= NTOK) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            WRITE(F06,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            RETURN
         ENDIF
         ITOK = ITOK + 1
         READ(TOKENS(ITOK),*,IOSTAT=IOCHK) I2_MOD
         IF (IOCHK /= 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            WRITE(F06,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            RETURN
         ENDIF
         ITOK = ITOK + 1

      ELSE IF (TRIM(OPTION_UP) == 'K1MOD') THEN
         IF (ITOK >= NTOK) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            WRITE(F06,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            RETURN
         ENDIF
         ITOK = ITOK + 1
         READ(TOKENS(ITOK),*,IOSTAT=IOCHK) K1_MOD
         IF (IOCHK /= 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            WRITE(F06,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            RETURN
         ENDIF
         ITOK = ITOK + 1

      ELSE IF (TRIM(OPTION_UP) == 'K2MOD') THEN
         IF (ITOK >= NTOK) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            WRITE(F06,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            RETURN
         ENDIF
         ITOK = ITOK + 1
         READ(TOKENS(ITOK),*,IOSTAT=IOCHK) K2_MOD
         IF (IOCHK /= 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            WRITE(F06,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            RETURN
         ENDIF
         ITOK = ITOK + 1

      ELSE IF (TRIM(OPTION_UP) == 'JMOD') THEN
         IF (ITOK >= NTOK) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            WRITE(F06,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            RETURN
         ENDIF
         ITOK = ITOK + 1
         READ(TOKENS(ITOK),*,IOSTAT=IOCHK) J_MOD
         IF (IOCHK /= 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            WRITE(F06,'(A,A,A,A,A,A,A)') ' *ERROR  1312: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID VALUE "', TRIM(TOKENS(ITOK)), '" FOR OPTION "', TRIM(OPTION_UP), '"'
            RETURN
         ENDIF
         ITOK = ITOK + 1

      ELSE IF (TRIM(OPTION_UP) == 'TAPER') THEN
         IF (ITOK >= NTOK) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            WRITE(F06,'(A,A,A,A,A)') ' *ERROR  1311: ', TRIM(CARD_NAME), ' ENTRY HAS OPTION "', TRIM(OPTION_UP), '" WITHOUT A FOLLOWING VALUE'
            RETURN
         ENDIF
         ITOK = ITOK + 1
         VALUE_TOKEN = TOKENS(ITOK)
         READ(VALUE_TOKEN,*,IOSTAT=IOCHK) TAPER_MODE
      IF (IOCHK /= 0) THEN
         CALL TO_UPPER ( VALUE_TOKEN )
            IF ((TRIM(VALUE_TOKEN) == 'LINEAR') .OR. (TRIM(VALUE_TOKEN) == 'DEFAULT')) THEN
               TAPER_MODE = 1
            ELSE IF (TRIM(VALUE_TOKEN) == 'PARABOLIC') THEN
               TAPER_MODE = 2
            ELSE IF (TRIM(VALUE_TOKEN) == 'CUBIC') THEN
               TAPER_MODE = 3
            ELSE IF ((TRIM(VALUE_TOKEN) == 'NONE') .OR. (TRIM(VALUE_TOKEN) == 'NONTAPER') .OR. &
                     (TRIM(VALUE_TOKEN) == 'NON-TAPER') .OR. (TRIM(VALUE_TOKEN) == 'CONSTANT') .OR. &
                     (TRIM(VALUE_TOKEN) == 'PRISMATIC') .OR. (TRIM(VALUE_TOKEN) == 'OFF')) THEN
               TAPER_MODE = 0
            ELSE
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,'(A,A,A,A,A,A,A)') ' *ERROR  1313: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID TAPER VALUE "', TRIM(VALUE_TOKEN), '" IN TOKEN "', TRIM(OPTION_UP), '"'
               WRITE(F06,'(A,A,A,A,A,A,A)') ' *ERROR  1313: ', TRIM(CARD_NAME), ' ENTRY HAS INVALID TAPER VALUE "', TRIM(VALUE_TOKEN), '" IN TOKEN "', TRIM(OPTION_UP), '"'
               RETURN
            ENDIF
         ENDIF
         ITOK = ITOK + 1

      ENDIF

      END SUBROUTINE READ_PBEAMZ_OPTION

! ##################################################################################################################################

      SUBROUTINE WRITE_PBEAML_CONVERTED_PBEAM_DEBUG ( WHICH, XL, SEC_TYPE_IN, NDIM_IN, DIMS, AREA, I1, I2, I12, JTOR, NSM, K1, K2, &
                                                      YC, ZC, YS, ZS, IWARP, STRE )

      CHARACTER(LEN=*), INTENT(IN) :: WHICH
      CHARACTER(LEN=*), INTENT(IN) :: SEC_TYPE_IN
      INTEGER(LONG), INTENT(IN)    :: NDIM_IN
      REAL(DOUBLE), INTENT(IN)     :: XL, DIMS(10), AREA, I1, I2, I12, JTOR, NSM, K1, K2, YC, ZC, YS, ZS, IWARP, STRE(8)
      INTEGER(LONG)                :: IDIM
      CHARACTER(LEN=256)           :: LINE

      WRITE(F06,'(A)') ' '
      WRITE(F06,'(A)') '*** PBEAML CONVERTED PBEAM DEBUG *********************************************'
      WRITE(F06,'(A,I0)') '  Property ID      : ', PROPERTY_ID
      WRITE(F06,'(A,I0)') '  Material ID      : ', PBEAM(NPBEAM,2)
      WRITE(F06,'(A,A)')  '  PBEAML section   : ', TRIM(SEC_TYPE_IN)
      IF (WHICH == 'A') THEN
         WRITE(F06,'(A)') '  Snapshot         : End A / station 1'
      ELSE IF (WHICH == 'B') THEN
         CALL WRITE_LABEL_VALUE ( '  Snapshot         : End B at x/L = ', XL )
      ELSE
         CALL WRITE_LABEL_VALUE ( '  Snapshot         : Intermediate station at x/L = ', XL )
      ENDIF
      WRITE(F06,'(A)') '  Input dimensions :'
      DO IDIM=1,NDIM_IN
         WRITE(LINE,'(A,I0,A,A)') '    DIM(', IDIM, ') = ', TRIM(FMT_REAL_SHORT(DIMS(IDIM)))
         WRITE(F06,'(A)') TRIM(LINE)
      ENDDO
      WRITE(F06,'(A)') '  Converted PBEAM fields :'
      CALL WRITE_LABEL_VALUE ( '    Area           = ', AREA )
      CALL WRITE_LABEL_VALUE ( '    I1             = ', I1 )
      CALL WRITE_LABEL_VALUE ( '    I2             = ', I2 )
      CALL WRITE_LABEL_VALUE ( '    I12            = ', I12 )
      CALL WRITE_LABEL_VALUE ( '    J              = ', JTOR )
      CALL WRITE_LABEL_VALUE ( '    NSM            = ', NSM )
      CALL WRITE_LABEL_VALUE ( '    K1 shear       = ', K1 )
      CALL WRITE_LABEL_VALUE ( '    K2 shear       = ', K2 )
      CALL WRITE_LABEL_VALUE ( '    CW             = ', IWARP )
      CALL WRITE_LABEL_VALUE ( '    Neutral axis Y = ', YC )
      CALL WRITE_LABEL_VALUE ( '    Neutral axis Z = ', ZC )
      CALL WRITE_LABEL_VALUE ( '    Shear center Y = ', YS )
      CALL WRITE_LABEL_VALUE ( '    Shear center Z = ', ZS )
      WRITE(F06,'(A)') '  Stress recovery points :'
      CALL WRITE_POINT_VALUE ( '    C = ', STRE(1), STRE(2) )
      CALL WRITE_POINT_VALUE ( '    D = ', STRE(3), STRE(4) )
      CALL WRITE_POINT_VALUE ( '    E = ', STRE(5), STRE(6) )
      CALL WRITE_POINT_VALUE ( '    F = ', STRE(7), STRE(8) )
      WRITE(F06,'(A)') '  NX-style converted PBEAM card image:'
      WRITE(F06,'(A)') '     THE USER SUPPLIED PBEAML BULK DATA ENTRY IS REPRESENTED INTERNALLY AS:'
      WRITE(LINE,'(A,I0,1X,I0)') '  PBEAM      ', PROPERTY_ID, PBEAM(NPBEAM,2)
      CALL WRITE_VALUE_LIST_6 ( TRIM(LINE), AREA, I1, I2, I12, JTOR, NSM )
      CALL WRITE_VALUE_LIST_8 ( '              ', STRE(1), STRE(2), STRE(3), STRE(4), STRE(5), STRE(6), STRE(7), STRE(8) )
      CALL WRITE_VALUE_LIST_8 ( '              ', K1, K2, 0.0D0, 0.0D0, 0.0D0, 0.0D0, IWARP, IWARP )
      CALL WRITE_VALUE_LIST_8 ( '              ', YS, ZS, YS, ZS, YC, ZC, YC, ZC )
      WRITE(F06,'(A)') '***************************************************************************'

      END SUBROUTINE WRITE_PBEAML_CONVERTED_PBEAM_DEBUG

! ##################################################################################################################################

      CHARACTER(LEN=24) FUNCTION FMT_REAL_SHORT ( VALUE )

      REAL(DOUBLE), INTENT(IN) :: VALUE
      CHARACTER(LEN=32)        :: BUFFER
      CHARACTER(LEN=24)        :: MANT
      CHARACTER(LEN=8)         :: EXPSTR
      INTEGER(LONG)            :: EPOS, IEND

      IF (DABS(VALUE) <= 1.0D-12) THEN
         FMT_REAL_SHORT = '0.0'
         RETURN
      ENDIF

      WRITE(BUFFER,'(ES16.8E2)') VALUE
      BUFFER = ADJUSTL(BUFFER)
      EPOS = INDEX(BUFFER,'E')
      IF (EPOS <= 0) THEN
         FMT_REAL_SHORT = TRIM(BUFFER)
         RETURN
      ENDIF

      MANT = BUFFER(1:EPOS-1)
      EXPSTR = BUFFER(EPOS:)
      IEND = LEN_TRIM(MANT)

      DO WHILE (IEND > 1)
         IF (MANT(IEND:IEND) /= '0') EXIT
         IEND = IEND - 1
      ENDDO

      IF (MANT(IEND:IEND) == '.') THEN
         MANT(IEND+1:IEND+1) = '0'
         IEND = IEND + 1
      ENDIF
      MANT = MANT(1:IEND)

      IF ((TRIM(EXPSTR) == 'E+00') .OR. (TRIM(EXPSTR) == 'E-00')) THEN
         FMT_REAL_SHORT = TRIM(MANT)
      ELSE
         FMT_REAL_SHORT = TRIM(MANT)//TRIM(EXPSTR)
      ENDIF

      END FUNCTION FMT_REAL_SHORT

! ##################################################################################################################################

      SUBROUTINE WRITE_LABEL_VALUE ( LABEL, VALUE )

      CHARACTER(LEN=*), INTENT(IN) :: LABEL
      REAL(DOUBLE), INTENT(IN)     :: VALUE

      WRITE(F06,'(A,1X,A)') TRIM(LABEL), TRIM(FMT_REAL_SHORT(VALUE))

      END SUBROUTINE WRITE_LABEL_VALUE

! ##################################################################################################################################

      SUBROUTINE WRITE_POINT_VALUE ( LABEL, VALUE1, VALUE2 )

      CHARACTER(LEN=*), INTENT(IN) :: LABEL
      REAL(DOUBLE), INTENT(IN)     :: VALUE1, VALUE2
      CHARACTER(LEN=128)           :: LINE

      LINE = LABEL//'('//TRIM(FMT_REAL_SHORT(VALUE1))//', '//TRIM(FMT_REAL_SHORT(VALUE2))//')'
      WRITE(F06,'(A)') TRIM(LINE)

      END SUBROUTINE WRITE_POINT_VALUE

! ##################################################################################################################################

      SUBROUTINE WRITE_VALUE_LIST_6 ( PREFIX, V1, V2, V3, V4, V5, V6 )

      CHARACTER(LEN=*), INTENT(IN) :: PREFIX
      REAL(DOUBLE), INTENT(IN)     :: V1, V2, V3, V4, V5, V6
      CHARACTER(LEN=256)           :: LINE

      LINE = PREFIX//' '//TRIM(FMT_REAL_SHORT(V1))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V2))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V3))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V4))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V5))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V6))
      WRITE(F06,'(A)') TRIM(LINE)

      END SUBROUTINE WRITE_VALUE_LIST_6

! ##################################################################################################################################

      SUBROUTINE WRITE_VALUE_LIST_8 ( PREFIX, V1, V2, V3, V4, V5, V6, V7, V8 )

      CHARACTER(LEN=*), INTENT(IN) :: PREFIX
      REAL(DOUBLE), INTENT(IN)     :: V1, V2, V3, V4, V5, V6, V7, V8
      CHARACTER(LEN=256)           :: LINE

      LINE = PREFIX//' '//TRIM(FMT_REAL_SHORT(V1))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V2))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V3))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V4))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V5))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V6))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V7))
      LINE = TRIM(LINE)//' '//TRIM(FMT_REAL_SHORT(V8))
      WRITE(F06,'(A)') TRIM(LINE)

      END SUBROUTINE WRITE_VALUE_LIST_8

! ##################################################################################################################################

      SUBROUTINE CALC_SECTION ( SEC_TYPE_IN, NDIM_IN, DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )

      CHARACTER(LEN=*), INTENT(IN)  :: SEC_TYPE_IN
      INTEGER(LONG), INTENT(IN)     :: NDIM_IN
      REAL(DOUBLE), INTENT(IN)      :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT)     :: AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE(8)

      CHARACTER(LEN=JCARD_LEN) :: SEC

      SEC = SEC_TYPE_IN
      CALL TO_UPPER ( SEC )

      AREA  = ZERO ; I1    = ZERO ; I2    = ZERO ; I12  = ZERO ; JTOR = ZERO
      K1    = 1.D0 ; K2    = 1.D0 ; YC    = ZERO ; ZC   = ZERO ; YS   = ZERO ; ZS = ZERO ; IWARP = ZERO
      STRE  = ZERO

      IF (SEC(1:5) == 'TUBE2') THEN
         CALL CALC_TUBE_SECTION ( DIMS, .TRUE., AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      ELSE IF (SEC(1:4) == 'I   ') THEN
         CALL CALC_I_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )
      ELSE IF (SEC(1:4) == 'ROD ') THEN
         CALL CALC_ROD_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      ELSE IF (SEC(1:4) == 'TUBE') THEN
         CALL CALC_TUBE_SECTION ( DIMS, .FALSE., AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      ELSE IF (SEC(1:4) == 'BAR ') THEN
         CALL CALC_BAR_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      ELSE IF (SEC(1:4) == 'BOX ') THEN
         CALL CALC_BOX_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      ELSE IF (SEC(1:4) == 'H   ') THEN
         CALL CALC_H_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )
      ELSE IF (SEC(1:4) == 'CHAN') THEN
         CALL CALC_CHAN_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, STRE )
      ELSE IF (SEC(1:4) == 'T   ') THEN
         CALL CALC_T_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, STRE )
      ELSE IF (SEC(1:4) == 'L   ') THEN
         CALL CALC_L_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, STRE )
      ELSE
         FATAL_ERR = FATAL_ERR + 1
      ENDIF

      END SUBROUTINE CALC_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_I_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )

      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE(8)

      REAL(DOUBLE) :: H, B1, B2, TW, T1, T2
      REAL(DOUBLE) :: A1, A2, A3
      REAL(DOUBLE) :: Z1, Z2, Z3, ZBAR, HWEB

      H  = DIMS(1)
      B1 = DIMS(2)
      B2 = DIMS(3)
      TW = DIMS(4)
      T1 = DIMS(5)
      T2 = DIMS(6)

      HWEB = H - T1 - T2
      IF (HWEB < ZERO) HWEB = ZERO

      A1 = B1*T1
      A2 = HWEB*TW
      A3 = B2*T2
      AREA = A1 + A2 + A3

      Z1 = 0.5D0*T2
      Z2 = T2 + 0.5D0*HWEB
      Z3 = H - 0.5D0*T1
      IF (AREA > ZERO) THEN
         ZBAR = (A1*Z3 + A2*Z2 + A3*Z1)/AREA
      ELSE
         ZBAR = ZERO
      ENDIF

      YC = ZERO
      ZC = ZBAR - 0.5D0*H
      YS = ZERO
      ZS = ZERO
      IWARP = ZERO

      I1 = (B1*T1**3)/12.D0 + A1*(Z3-ZBAR)**2                                                  &
         + (TW*HWEB**3)/12.D0 + A2*(Z2-ZBAR)**2                                                &
         + (B2*T2**3)/12.D0 + A3*(Z1-ZBAR)**2

      I2 = (T1*B1**3)/12.D0 + (HWEB*TW**3)/12.D0 + (T2*B2**3)/12.D0
      I12 = ZERO
      JTOR = (B1*T1**3 + B2*T2**3 + HWEB*TW**3)/3.D0
      K1 = 0.5D0
      K2 = 0.5D0

      STRE(1) =  0.5D0*B1
      STRE(2) =  0.5D0*H - ZBAR
      STRE(3) = -0.5D0*B1
      STRE(4) =  0.5D0*H - ZBAR
      STRE(5) = -0.5D0*B2
      STRE(6) = -0.5D0*H - ZBAR
      STRE(7) =  0.5D0*B2
      STRE(8) = -0.5D0*H - ZBAR

      END SUBROUTINE CALC_I_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_H_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE(8)

! Temporary phase-1 support: treat H as an I-section rotated 90 degrees with equal flange thicknesses.
! NX reverse-engineering for H is still thin, so keep this mapping conservative and explicit.
      REAL(DOUBLE) :: IDIMS(10)

      IDIMS = ZERO
      IDIMS(1) = DIMS(1)
      IDIMS(2) = DIMS(3)
      IDIMS(3) = DIMS(3)
      IDIMS(4) = DIMS(4)
      IDIMS(5) = DIMS(2)
      IDIMS(6) = DIMS(2)

      CALL CALC_I_SECTION ( IDIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, IWARP, STRE )

      END SUBROUTINE CALC_H_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_ROD_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, STRE(8)
      REAL(DOUBLE) :: R, PI
      PI = 4.D0*DATAN(1.D0)
      R = DIMS(1)
      AREA = PI*R**2
      I1 = PI*R**4/4.D0
      I2 = I1
      I12 = ZERO
      JTOR = PI*R**4/2.D0
      K1 = 0.9D0 ; K2 = 0.9D0
      STRE = (/ R, ZERO, ZERO, R, -R, ZERO, ZERO, -R /)
      END SUBROUTINE CALC_ROD_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_TUBE_SECTION ( DIMS, USE_THICKNESS, AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      LOGICAL, INTENT(IN)       :: USE_THICKNESS
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, STRE(8)
      REAL(DOUBLE) :: RO, RI, T, PI
      PI = 4.D0*DATAN(1.D0)
      IF (USE_THICKNESS) THEN
         RO = 0.5D0*DIMS(1)
         T  = DIMS(2)
         RI = RO - T
      ELSE
         RO = DIMS(1)
         RI = DIMS(2)
      ENDIF
      IF (RI < ZERO) RI = ZERO
      AREA = PI*(RO**2 - RI**2)
      I1 = PI*(RO**4 - RI**4)/4.D0
      I2 = I1
      I12 = ZERO
      JTOR = PI*(RO**4 - RI**4)/2.D0
      K1 = 0.5D0 ; K2 = 0.5D0
      STRE = (/ RO, ZERO, ZERO, RO, -RO, ZERO, ZERO, -RO /)
      END SUBROUTINE CALC_TUBE_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_BAR_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, STRE(8)
      REAL(DOUBLE) :: H, B, A, BT, RATIO
      B = DIMS(1)
      H = DIMS(2)
      AREA = B*H
      I1 = B*H**3/12.D0
      I2 = H*B**3/12.D0
      I12 = ZERO
      A = DMAX1(H,B)/2.D0
      BT = DMIN1(H,B)/2.D0
      ! Use the same closed-form rectangle expression for the full aspect-ratio range,
      ! including squares, so generated BAR torsion constants stay consistent with the
      ! explicit PBEAM station values used in validation decks.
      RATIO = BT/A
      JTOR = A*BT**3*(16.D0/3.D0 - 3.36D0*RATIO*(1.D0 - RATIO**4/12.D0))
      K1 = 5.D0/6.D0 ; K2 = 5.D0/6.D0
      STRE = (/ 0.5D0*B, 0.5D0*H, -0.5D0*B, 0.5D0*H, -0.5D0*B, -0.5D0*H, 0.5D0*B, -0.5D0*H /)
      END SUBROUTINE CALC_BAR_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_BOX_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, STRE(8)
      REAL(DOUBLE) :: H, B, T1, T2, HI, BI, AM, PERIM, TAVG
      H = DIMS(1); B = DIMS(2); T1 = DIMS(3); T2 = DIMS(4)
      AREA = 2.D0*(B*T1 + (H - 2.D0*T1)*T2)
      HI = H - 2.D0*T1
      BI = B - 2.D0*T2
      IF (HI < ZERO) HI = ZERO
      IF (BI < ZERO) BI = ZERO
      I1 = B*H**3/12.D0 - BI*HI**3/12.D0
      I2 = H*B**3/12.D0 - HI*BI**3/12.D0
      I12 = ZERO
      AM = (B - T2)*(H - T1)
      PERIM = 2.D0*((B - T2) + (H - T1))
      TAVG = ZERO
      IF ((B/T2 + H/T1) > ZERO) TAVG = PERIM/(2.D0*(B/T2 + H/T1))
      IF ((PERIM > ZERO) .AND. (TAVG > ZERO)) THEN
         JTOR = 4.D0*AM**2/(PERIM/TAVG)
      ELSE
         JTOR = ZERO
      ENDIF
      K1 = 0.5D0 ; K2 = 0.5D0
      STRE = (/ 0.5D0*B, 0.5D0*H, -0.5D0*B, 0.5D0*H, -0.5D0*B, -0.5D0*H, 0.5D0*B, -0.5D0*H /)
      END SUBROUTINE CALC_BOX_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_CHAN_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, YS, ZS, STRE(8)
      REAL(DOUBLE) :: B, H, TW, TF, AWEB, AFLG, AWEB_EFF, AFLG_EFF, YWEB, YFLG, ZBAR, E
      B  = DIMS(1)
      H  = DIMS(2)
      TW = DIMS(3)
      TF = DIMS(4)
      AWEB = TW*H
      AFLG = 2.D0*B*TF
      AREA = AWEB + AFLG
      YWEB = H/2.D0
      YFLG = H - TF/2.D0
      ZBAR = ZERO
      IF (AREA > ZERO) ZBAR = (AWEB*YWEB + AFLG*YFLG)/AREA
      YC = ZERO
      ZC = ZBAR - H/2.D0
      I1 = TW*H**3/12.D0 + AWEB*(YWEB - ZBAR)**2 + 2.D0*(B*TF**3/12.D0 + B*TF*(YFLG - ZBAR)**2)
      I2 = H*TW**3/12.D0 + 2.D0*TF*B**3/12.D0
      I12 = ZERO
      JTOR = (2.D0*B*TF**3 + H*TW**3)/3.D0
      AWEB_EFF = TW*(H - 2.D0*TF)
      IF (AWEB_EFF < ZERO) AWEB_EFF = ZERO
      AFLG_EFF = ZERO
      IF (B > ZERO) AFLG_EFF = AFLG*(1.D0 - TF/(H + B))
      K1 = ZERO
      K2 = ZERO
      IF (AREA > ZERO) THEN
         K1 = AWEB_EFF/AREA
         K2 = AFLG_EFF/AREA
      ENDIF
      E = ZERO
      IF (I2 > ZERO) E = B**2*TF*H**2*TW/(4.D0*I2)
      YS = E
      ZS = ZERO
      STRE = (/ B, 0.5D0*H, ZERO, 0.5D0*H, ZERO, -0.5D0*H, B, -0.5D0*H /)
      END SUBROUTINE CALC_CHAN_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_T_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, STRE(8)
      REAL(DOUBLE) :: B, H, TF, TW, AWEB, AFLG, YWEB, YFLG, ZBAR
      B  = DIMS(1)
      H  = DIMS(2)
      TF = DIMS(3)
      TW = DIMS(4)
      AWEB = TW*(H - TF)
      AFLG = B*TF
      AREA = AWEB + AFLG
      YWEB = 0.5D0*(H - TF)
      YFLG = H - 0.5D0*TF
      ZBAR = ZERO
      IF (AREA > ZERO) ZBAR = (AWEB*YWEB + AFLG*YFLG)/AREA
      YC = ZERO
      ZC = ZBAR - H/2.D0
      I1 = TW*(H - TF)**3/12.D0 + AWEB*(YWEB - ZBAR)**2 + B*TF**3/12.D0 + AFLG*(YFLG - ZBAR)**2
      I2 = (H - TF)*TW**3/12.D0 + TF*B**3/12.D0
      I12 = ZERO
      JTOR = (B*TF**3 + (H - TF)*TW**3)/3.D0
      K1 = ZERO
      K2 = ZERO
      IF (AREA > ZERO) THEN
         K1 = AWEB/AREA
         K2 = AFLG/AREA
      ENDIF
      STRE = (/ 0.5D0*TW, ZERO, 0.5D0*TW, 0.5D0*H, -B + 0.5D0*TW, ZERO, 0.5D0*TW, -0.5D0*H /)
      END SUBROUTINE CALC_T_SECTION

! ##################################################################################################################################

      SUBROUTINE CALC_L_SECTION ( DIMS, AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, STRE )
      REAL(DOUBLE), INTENT(IN)  :: DIMS(10)
      REAL(DOUBLE), INTENT(OUT) :: AREA, I1, I2, I12, JTOR, K1, K2, YC, ZC, STRE(8)
      REAL(DOUBLE) :: B, H, TF, TW, A1, A2, Y1, Z1, Y2, Z2, YBAR, ZBAR
      B  = DIMS(1)
      H  = DIMS(2)
      TF = DIMS(3)
      TW = DIMS(4)
      A1 = (H - TF)*TW
      A2 = (B - TW)*TF
      AREA = A1 + A2
      Y1 = 0.5D0*(H - TF) ; Z1 = 0.5D0*TW
      Y2 = 0.5D0*TF ; Z2 = TW + 0.5D0*(B - TW)
      YBAR = ZERO ; ZBAR = ZERO
      IF (AREA > ZERO) THEN
         YBAR = (A1*Y1 + A2*Y2)/AREA
         ZBAR = (A1*Z1 + A2*Z2)/AREA
      ENDIF
      YC = YBAR - H/2.D0
      ZC = ZBAR - B/2.D0
      I1 = TW*(H - TF)**3/12.D0 + A1*(Y1 - YBAR)**2 + (B - TW)*TF**3/12.D0 + A2*(Y2 - YBAR)**2
      I2 = (H - TF)*TW**3/12.D0 + A1*(Z1 - ZBAR)**2 + TF*(B - TW)**3/12.D0 + A2*(Z2 - ZBAR)**2
      I12 = A1*(Y1 - YBAR)*(Z1 - ZBAR) + A2*(Y2 - YBAR)*(Z2 - ZBAR)
      JTOR = ((H - TF)*TW**3 + (B - TW)*TF**3)/3.D0
      K1 = ZERO
      K2 = ZERO
      IF (AREA > ZERO) THEN
         K1 = A1/AREA
         K2 = A2/AREA
      ENDIF
      STRE = (/ B - ZBAR, 0.5D0*TF, -ZBAR, H - YBAR, -ZBAR, -YBAR, B - ZBAR, -YBAR /)
      END SUBROUTINE CALC_L_SECTION

! ##################################################################################################################################

      SUBROUTINE TO_UPPER ( STR )

      CHARACTER(LEN=*), INTENT(INOUT) :: STR
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: ICODE

      DO K=1,LEN(STR)
         ICODE = IACHAR(STR(K:K))
         IF ((ICODE >= IACHAR('a')) .AND. (ICODE <= IACHAR('z'))) THEN
            STR(K:K) = ACHAR(ICODE - 32)
         ENDIF
      ENDDO

      END SUBROUTINE TO_UPPER

! ##################################################################################################################################

      END SUBROUTINE BD_PBEAML

