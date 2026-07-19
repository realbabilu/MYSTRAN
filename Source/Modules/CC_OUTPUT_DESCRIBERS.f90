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

      MODULE CC_OUTPUT_DESCRIBERS

      USE PENTIUM_II_KIND, ONLY         :  BYTE, LONG
      USE SCONTR, ONLY                  :  CC_CMD_DESCRIBERS, LSETLN

      IMPLICIT NONE

      SAVE

      ! The following are the default values for the Case Control command
      ! describers.  These are used to check the values in parens () in Case
      ! Control entries to see if they are valid for MYSTRAN.  The ones below
      ! are ones MYSTRAN honors.  If other values are encountered in the user's
      !  Case Control, warning messages are written to the MYSTRAN output file.

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: ACCE_SORT = 'SORT1   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: ACCE_OUT  = 'YNNNN   '  ! print, plot, punch, neu, csv
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: ACCE_MAG  = 'MAG     '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: DISP_SORT = 'SORT1   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: DISP_OUT  = 'YNNNN   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: DISP_MAG  = 'MAG     '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: VELO_SORT = 'SORT1   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: VELO_OUT  = 'YNNNN   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: VELO_MAG  = 'MAG     '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: GPFO_SORT = 'SORT1   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: GPFO_OUT  = 'YNNNN   '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: MPCF_SORT = 'SORT1   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: MPCF_OUT  = 'YNNNN   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: MPCF_MAG  = 'MAG     '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: OLOA_SORT = 'SORT1   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: OLOA_OUT  = 'YNNNN   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: OLOA_MAG  = 'MAG     '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: SPCF_SORT = 'SORT1   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: SPCF_OUT  = 'YNNNN   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: SPCF_MAG  = 'MAG     '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: FORC_SORT = 'SORT1   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: FORC_OUT  = 'YNNNN   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: FORC_MAG  = 'MAG     '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: FORC_LOC  = 'CENTER  '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRN_SORT = 'SORT1   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRN_OUT  = 'YNNNN   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRN_MAG  = 'MAG     '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRN_OPT  = 'VONMISES'
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRN_CUR  = 'STRCUR  '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRN_LOC  = 'CENTER  '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRE_SORT = 'SORT1   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRE_OUT  = 'YNNNN   '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRE_MAG  = 'MAG     '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRE_OPT  = 'VONMISES'
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRE_LOC  = 'CENTER  '

      INTEGER(LONG), PARAMETER          :: MAX_GP_SURFACES = 256
      INTEGER(LONG), PARAMETER          :: MAX_GP_POST_SETS = 256
      INTEGER(LONG), PARAMETER          :: GP_POST_SET_TEXT_LEN = 1024

      LOGICAL                           :: GPSTRESS_REQ      = .FALSE.
      LOGICAL                           :: STRFIELD_REQ      = .FALSE.
      LOGICAL                           :: OUTPUT_POST_REQ   = .FALSE.
      INTEGER(LONG)                     :: GPSTRESS_SETID    = 0
      INTEGER(LONG)                     :: NUM_GP_SURFACE    = 0
      INTEGER(LONG)                     :: NUM_GP_VOLUME     = 0
      INTEGER(LONG)                     :: GP_SURFACE_IDS(MAX_GP_SURFACES)       = 0
      INTEGER(LONG)                     :: GP_SURFACE_SETIDS(MAX_GP_SURFACES)    = 0
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: GP_SURFACE_NORMAL_MODE(MAX_GP_SURFACES)= '        '
      INTEGER(LONG)                     :: NUM_GP_POST_SET = 0
      INTEGER(LONG)                     :: GP_POST_SET_IDS(MAX_GP_POST_SETS) = 0
      CHARACTER(LEN=GP_POST_SET_TEXT_LEN):: GP_POST_SET_TEXT(MAX_GP_POST_SETS) = ' '

      CONTAINS

      SUBROUTINE RESET_GPSTRESS_POST_STATE

      INTEGER(LONG)                   :: I

      GPSTRESS_REQ    = .FALSE.
      STRFIELD_REQ    = .FALSE.
      OUTPUT_POST_REQ = .FALSE.
      GPSTRESS_SETID  = 0
      NUM_GP_SURFACE  = 0
      NUM_GP_VOLUME   = 0
      NUM_GP_POST_SET = 0

      DO I=1,MAX_GP_SURFACES
         GP_SURFACE_IDS(I) = 0
         GP_SURFACE_SETIDS(I) = 0
         GP_SURFACE_NORMAL_MODE(I) = '        '
      ENDDO
      DO I=1,MAX_GP_POST_SETS
         GP_POST_SET_IDS(I) = 0
         GP_POST_SET_TEXT(I) = ' '
      ENDDO

      END SUBROUTINE RESET_GPSTRESS_POST_STATE

      SUBROUTINE REGISTER_GP_SURFACE ( SURFACE_ID, SETID, NORMAL_MODE, IERR )

      INTEGER(LONG), INTENT(IN)       :: SURFACE_ID
      INTEGER(LONG), INTENT(IN)       :: SETID
      CHARACTER(LEN=*), INTENT(IN)    :: NORMAL_MODE
      INTEGER(LONG), INTENT(OUT)      :: IERR

      INTEGER(LONG)                   :: I

      IERR = 0
      IF (NUM_GP_SURFACE >= MAX_GP_SURFACES) THEN
         IERR = 1
         RETURN
      ENDIF

      DO I=1,NUM_GP_SURFACE
         IF (GP_SURFACE_IDS(I) == SURFACE_ID) THEN
            GP_SURFACE_SETIDS(I) = SETID
            GP_SURFACE_NORMAL_MODE(I) = NORMAL_MODE
            RETURN
         ENDIF
      ENDDO

      NUM_GP_SURFACE = NUM_GP_SURFACE + 1
      GP_SURFACE_IDS(NUM_GP_SURFACE) = SURFACE_ID
      GP_SURFACE_SETIDS(NUM_GP_SURFACE) = SETID
      GP_SURFACE_NORMAL_MODE(NUM_GP_SURFACE) = NORMAL_MODE

      END SUBROUTINE REGISTER_GP_SURFACE

      SUBROUTINE REGISTER_GP_POST_SET ( SETID, SET_TEXT, IERR )

      INTEGER(LONG), INTENT(IN)       :: SETID
      CHARACTER(LEN=*), INTENT(IN)    :: SET_TEXT
      INTEGER(LONG), INTENT(OUT)      :: IERR

      INTEGER(LONG)                   :: I

      IERR = 0
      DO I=1,NUM_GP_POST_SET
         IF (GP_POST_SET_IDS(I) == SETID) THEN
            GP_POST_SET_TEXT(I) = ' '
            GP_POST_SET_TEXT(I)(1:MIN(LEN(GP_POST_SET_TEXT(I)),LEN_TRIM(SET_TEXT))) = SET_TEXT(1:MIN(LEN(GP_POST_SET_TEXT(I)),LEN_TRIM(SET_TEXT)))
            RETURN
         ENDIF
      ENDDO

      IF (NUM_GP_POST_SET >= MAX_GP_POST_SETS) THEN
         IERR = 1
         RETURN
      ENDIF

      NUM_GP_POST_SET = NUM_GP_POST_SET + 1
      GP_POST_SET_IDS(NUM_GP_POST_SET) = SETID
      GP_POST_SET_TEXT(NUM_GP_POST_SET) = ' '
      GP_POST_SET_TEXT(NUM_GP_POST_SET)(1:MIN(LEN(GP_POST_SET_TEXT(NUM_GP_POST_SET)),LEN_TRIM(SET_TEXT))) = &
         SET_TEXT(1:MIN(LEN(GP_POST_SET_TEXT(NUM_GP_POST_SET)),LEN_TRIM(SET_TEXT)))

      END SUBROUTINE REGISTER_GP_POST_SET

      END MODULE CC_OUTPUT_DESCRIBERS
