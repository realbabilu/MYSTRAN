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

      SUBROUTINE CC_MPF_MEFM ( CARD, WHAT, CALC_FLAG )

! Processes Case Control MEFFMASS/MPFACTOR entries.
!
! Current compatibility policy:
!   1) Legacy MYSTRAN forms such as "MEFFMASS = ALL" still work as before.
!   2) Nastran/Altair-style forms such as "MEFFMASS(ALL)=YES" are accepted.
!   3) Core parenthesized descriptors SUMMARY, PARTFAC, MEFFM, MEFFW, FRACSUM, and GRID=
!      are mapped onto current MYSTRAN controls.
!   4) This is still a compatibility bridge. Descriptor parsing is Nastran-style, but unsupported
!      extended options still fall back to current MYSTRAN behavior.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  WARN_ERR, BLNK_SUB_NAM, LSUB, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN, MEFMLOC, MEFMGRID
      USE MODEL_STUF, ONLY            :  MEFFMASS_CALC, MPFACTOR_CALC, MEFFMASS_REQ_SUMMARY, MEFFMASS_REQ_MEFFM,                &
                                         MEFFMASS_REQ_MEFFW, MEFFMASS_REQ_FRACSUM, MPFACTOR_REQ_PARTFAC,                        &
                                         MEFFMASS_CALC_SUB, MPFACTOR_CALC_SUB, MEFFMASS_REQ_SUMMARY_SUB,                        &
                                         MEFFMASS_REQ_MEFFM_SUB, MEFFMASS_REQ_MEFFW_SUB, MEFFMASS_REQ_FRACSUM_SUB,              &
                                         MPFACTOR_REQ_PARTFAC_SUB, MEFMLOC_SUB, MEFMGRID_SUB

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_MPF_MEFM'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD
      CHARACTER(LEN=*), INTENT(IN)    :: WHAT
      CHARACTER( 1*BYTE), INTENT(INOUT):: CALC_FLAG

      CHARACTER(LEN=LEN(CARD))        :: REST
      CHARACTER(LEN=LEN(CARD))        :: RHS
      CHARACTER(LEN=LEN(CARD))        :: DESCR

      INTEGER(LONG)                   :: EQ
      INTEGER(LONG)                   :: GRID_EQ
      INTEGER(LONG)                   :: GRID_END
      INTEGER(LONG)                   :: GRID_ID
      INTEGER(LONG)                   :: IOS
      INTEGER(LONG)                   :: LP
      INTEGER(LONG)                   :: RP
      INTEGER(LONG)                   :: LWHAT
      INTEGER(LONG)                   :: I

      LOGICAL                         :: COMPAT_STYLE
      LOGICAL                         :: ENABLE_REQ
      LOGICAL                         :: HAS_ALL
      LOGICAL                         :: HAS_SUMMARY
      LOGICAL                         :: HAS_PARTFAC
      LOGICAL                         :: HAS_MEFFM
      LOGICAL                         :: HAS_MEFFW
      LOGICAL                         :: HAS_FRACSUM
      LOGICAL                         :: HAS_GRID
      LOGICAL                         :: UNKNOWN_RHS
      LOGICAL                         :: UNKNOWN_DESCR

      CHARACTER(LEN=32)               :: GRID_TOKEN




! **********************************************************************************************************************************
! Process MEFFMASS/MPFACTOR cards

      LWHAT        = LEN_TRIM(WHAT)
      REST         = ' '
      RHS          = ' '
      DESCR        = ' '
      COMPAT_STYLE = .FALSE.
      ENABLE_REQ   = .TRUE.
      HAS_ALL      = .FALSE.
      HAS_SUMMARY  = .FALSE.
      HAS_PARTFAC  = .FALSE.
      HAS_MEFFM    = .FALSE.
      HAS_MEFFW    = .FALSE.
      HAS_FRACSUM  = .FALSE.
      HAS_GRID     = .FALSE.
      UNKNOWN_RHS  = .FALSE.
      UNKNOWN_DESCR = .FALSE.
      GRID_ID = -1

      IF (ALLOCATED(MEFMLOC_SUB)) THEN
         IF (NSUB == 0) THEN
            DO I=1,LSUB
               MEFMLOC_SUB(I) = MEFMLOC
               MEFMGRID_SUB(I) = MEFMGRID
            ENDDO
         ELSE
            MEFMLOC_SUB(NSUB) = MEFMLOC
            MEFMGRID_SUB(NSUB) = MEFMGRID
         ENDIF
      ENDIF

      IF (LEN_TRIM(CARD) > LWHAT) THEN
         REST = ADJUSTL(CARD(LWHAT+1:))
      ENDIF

      EQ = INDEX(REST,'=')
      LP = INDEX(REST,'(')
      RP = INDEX(REST,')')

      IF ((LP > 0) .AND. (RP > LP) .AND. ((EQ == 0) .OR. (LP < EQ))) THEN
         COMPAT_STYLE = .TRUE.
         DESCR = ADJUSTL(REST(LP+1:RP-1))
         IF (EQ > 0) THEN
            RHS = ADJUSTL(REST(EQ+1:))
         ENDIF
      ELSE IF (EQ > 0) THEN
         RHS = ADJUSTL(REST(EQ+1:))
      ELSE
         RHS = ADJUSTL(REST)
      ENDIF

      IF (LEN_TRIM(RHS) > 0) THEN
         IF ((RHS(1:2) == 'NO') .OR. (RHS(1:4) == 'NONE')) THEN
            ENABLE_REQ = .FALSE.
         ELSE IF ((RHS(1:3) == 'YES') .OR. (RHS(1:3) == 'ALL')) THEN
            ENABLE_REQ = .TRUE.
         ELSE
            ENABLE_REQ = .TRUE.
            UNKNOWN_RHS = .TRUE.
         ENDIF
      ENDIF

      IF (LEN_TRIM(DESCR) > 0) THEN
         HAS_ALL     = (INDEX(DESCR,'ALL')     > 0)
         HAS_SUMMARY = (INDEX(DESCR,'SUMMARY') > 0)
         HAS_PARTFAC = (INDEX(DESCR,'PARTFAC') > 0)
         HAS_MEFFM   = (INDEX(DESCR,'MEFFM')   > 0)
         HAS_MEFFW   = (INDEX(DESCR,'MEFFW')   > 0)
         HAS_FRACSUM = (INDEX(DESCR,'FRACSUM') > 0)
         GRID_EQ     = INDEX(DESCR,'GRID=')
         IF (GRID_EQ > 0) THEN
            HAS_GRID = .TRUE.
            GRID_TOKEN = ' '
            GRID_END = INDEX(DESCR(GRID_EQ:),',')
            IF (GRID_END > 0) THEN
               GRID_TOKEN = ADJUSTL(DESCR(GRID_EQ+5:GRID_EQ+GRID_END-2))
            ELSE
               GRID_TOKEN = ADJUSTL(DESCR(GRID_EQ+5:LEN_TRIM(DESCR)))
            ENDIF
            READ(GRID_TOKEN,*,IOSTAT=IOS) GRID_ID
            IF ((IOS == 0) .AND. (GRID_ID >= 0)) THEN
               MEFMLOC  = 'GRID  '
               MEFMGRID = GRID_ID
               IF (ALLOCATED(MEFMLOC_SUB)) THEN
                  IF (NSUB == 0) THEN
                     DO I=1,LSUB
                        MEFMLOC_SUB(I) = 'GRID  '
                        MEFMGRID_SUB(I) = GRID_ID
                     ENDDO
                  ELSE
                     MEFMLOC_SUB(NSUB) = 'GRID  '
                     MEFMGRID_SUB(NSUB) = GRID_ID
                  ENDIF
               ENDIF
            ELSE
               GRID_ID = -1
            ENDIF
         ENDIF
      ENDIF

      IF (ENABLE_REQ) THEN

         IF (.NOT. COMPAT_STYLE) THEN
            IF (WHAT == 'MEFFMASS') THEN
               MEFFMASS_CALC = 'Y'
               MEFFMASS_REQ_SUMMARY = 'Y'
               MEFFMASS_REQ_MEFFM   = 'Y'
               MEFFMASS_REQ_FRACSUM = 'Y'
               IF (ALLOCATED(MEFFMASS_CALC_SUB)) THEN
                  IF (NSUB == 0) THEN
                     DO I=1,LSUB
                        MEFFMASS_CALC_SUB(I) = 'Y'
                        MEFFMASS_REQ_SUMMARY_SUB(I) = 'Y'
                        MEFFMASS_REQ_MEFFM_SUB(I) = 'Y'
                        MEFFMASS_REQ_FRACSUM_SUB(I) = 'Y'
                     ENDDO
                  ELSE
                     MEFFMASS_CALC_SUB(NSUB) = 'Y'
                     MEFFMASS_REQ_SUMMARY_SUB(NSUB) = 'Y'
                     MEFFMASS_REQ_MEFFM_SUB(NSUB) = 'Y'
                     MEFFMASS_REQ_FRACSUM_SUB(NSUB) = 'Y'
                  ENDIF
               ENDIF
            ELSE IF (WHAT == 'MPFACTOR') THEN
               MPFACTOR_CALC = 'Y'
               MPFACTOR_REQ_PARTFAC = 'Y'
               IF (ALLOCATED(MPFACTOR_CALC_SUB)) THEN
                  IF (NSUB == 0) THEN
                     DO I=1,LSUB
                        MPFACTOR_CALC_SUB(I) = 'Y'
                        MPFACTOR_REQ_PARTFAC_SUB(I) = 'Y'
                     ENDDO
                  ELSE
                     MPFACTOR_CALC_SUB(NSUB) = 'Y'
                     MPFACTOR_REQ_PARTFAC_SUB(NSUB) = 'Y'
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            IF ((LEN_TRIM(DESCR) == 0) .OR. HAS_ALL) THEN
               IF (WHAT == 'MEFFMASS') THEN
                  MEFFMASS_CALC = 'Y'
                  MEFFMASS_REQ_SUMMARY = 'Y'
                  MEFFMASS_REQ_MEFFM   = 'Y'
                  MEFFMASS_REQ_MEFFW   = 'Y'
                  MEFFMASS_REQ_FRACSUM = 'Y'
                  MPFACTOR_CALC        = 'Y'
                  MPFACTOR_REQ_PARTFAC = 'Y'
                  IF (ALLOCATED(MEFFMASS_CALC_SUB)) THEN
                     IF (NSUB == 0) THEN
                        DO I=1,LSUB
                           MEFFMASS_CALC_SUB(I) = 'Y'
                           MEFFMASS_REQ_SUMMARY_SUB(I) = 'Y'
                           MEFFMASS_REQ_MEFFM_SUB(I) = 'Y'
                           MEFFMASS_REQ_MEFFW_SUB(I) = 'Y'
                           MEFFMASS_REQ_FRACSUM_SUB(I) = 'Y'
                           MPFACTOR_CALC_SUB(I) = 'Y'
                           MPFACTOR_REQ_PARTFAC_SUB(I) = 'Y'
                        ENDDO
                     ELSE
                        MEFFMASS_CALC_SUB(NSUB) = 'Y'
                        MEFFMASS_REQ_SUMMARY_SUB(NSUB) = 'Y'
                        MEFFMASS_REQ_MEFFM_SUB(NSUB) = 'Y'
                        MEFFMASS_REQ_MEFFW_SUB(NSUB) = 'Y'
                        MEFFMASS_REQ_FRACSUM_SUB(NSUB) = 'Y'
                        MPFACTOR_CALC_SUB(NSUB) = 'Y'
                        MPFACTOR_REQ_PARTFAC_SUB(NSUB) = 'Y'
                     ENDIF
                  ENDIF
               ELSE IF (WHAT == 'MPFACTOR') THEN
                  MPFACTOR_CALC        = 'Y'
                  MPFACTOR_REQ_PARTFAC = 'Y'
                  IF (ALLOCATED(MPFACTOR_CALC_SUB)) THEN
                     IF (NSUB == 0) THEN
                        DO I=1,LSUB
                           MPFACTOR_CALC_SUB(I) = 'Y'
                           MPFACTOR_REQ_PARTFAC_SUB(I) = 'Y'
                        ENDDO
                     ELSE
                        MPFACTOR_CALC_SUB(NSUB) = 'Y'
                        MPFACTOR_REQ_PARTFAC_SUB(NSUB) = 'Y'
                     ENDIF
                  ENDIF
               ENDIF
            ELSE
               IF (WHAT == 'MEFFMASS') THEN
                  IF (HAS_SUMMARY) THEN
                     MEFFMASS_CALC = 'Y'
                     MEFFMASS_REQ_SUMMARY = 'Y'
                     IF (ALLOCATED(MEFFMASS_CALC_SUB)) THEN
                        IF (NSUB == 0) THEN
                           DO I=1,LSUB
                              MEFFMASS_CALC_SUB(I) = 'Y'
                              MEFFMASS_REQ_SUMMARY_SUB(I) = 'Y'
                           ENDDO
                        ELSE
                           MEFFMASS_CALC_SUB(NSUB) = 'Y'
                           MEFFMASS_REQ_SUMMARY_SUB(NSUB) = 'Y'
                        ENDIF
                     ENDIF
                  ENDIF
                  IF (HAS_MEFFM) THEN
                     MEFFMASS_CALC = 'Y'
                     MEFFMASS_REQ_MEFFM = 'Y'
                     IF (ALLOCATED(MEFFMASS_CALC_SUB)) THEN
                        IF (NSUB == 0) THEN
                           DO I=1,LSUB
                              MEFFMASS_CALC_SUB(I) = 'Y'
                              MEFFMASS_REQ_MEFFM_SUB(I) = 'Y'
                           ENDDO
                        ELSE
                           MEFFMASS_CALC_SUB(NSUB) = 'Y'
                           MEFFMASS_REQ_MEFFM_SUB(NSUB) = 'Y'
                        ENDIF
                     ENDIF
                  ENDIF
                  IF (HAS_MEFFW) THEN
                     MEFFMASS_CALC = 'Y'
                     MEFFMASS_REQ_MEFFW = 'Y'
                     IF (ALLOCATED(MEFFMASS_CALC_SUB)) THEN
                        IF (NSUB == 0) THEN
                           DO I=1,LSUB
                              MEFFMASS_CALC_SUB(I) = 'Y'
                              MEFFMASS_REQ_MEFFW_SUB(I) = 'Y'
                           ENDDO
                        ELSE
                           MEFFMASS_CALC_SUB(NSUB) = 'Y'
                           MEFFMASS_REQ_MEFFW_SUB(NSUB) = 'Y'
                        ENDIF
                     ENDIF
                  ENDIF
                  IF (HAS_FRACSUM) THEN
                     MEFFMASS_CALC = 'Y'
                     MEFFMASS_REQ_FRACSUM = 'Y'
                     IF (ALLOCATED(MEFFMASS_CALC_SUB)) THEN
                        IF (NSUB == 0) THEN
                           DO I=1,LSUB
                              MEFFMASS_CALC_SUB(I) = 'Y'
                              MEFFMASS_REQ_FRACSUM_SUB(I) = 'Y'
                           ENDDO
                        ELSE
                           MEFFMASS_CALC_SUB(NSUB) = 'Y'
                           MEFFMASS_REQ_FRACSUM_SUB(NSUB) = 'Y'
                        ENDIF
                     ENDIF
                  ENDIF
                  IF (HAS_PARTFAC) THEN
                     MPFACTOR_CALC = 'Y'
                     MPFACTOR_REQ_PARTFAC = 'Y'
                     IF (ALLOCATED(MPFACTOR_CALC_SUB)) THEN
                        IF (NSUB == 0) THEN
                           DO I=1,LSUB
                              MPFACTOR_CALC_SUB(I) = 'Y'
                              MPFACTOR_REQ_PARTFAC_SUB(I) = 'Y'
                           ENDDO
                        ELSE
                           MPFACTOR_CALC_SUB(NSUB) = 'Y'
                           MPFACTOR_REQ_PARTFAC_SUB(NSUB) = 'Y'
                        ENDIF
                     ENDIF
                  ENDIF
               ELSE IF (WHAT == 'MPFACTOR') THEN
                  MPFACTOR_CALC = 'Y'
                  MPFACTOR_REQ_PARTFAC = 'Y'
                  IF (ALLOCATED(MPFACTOR_CALC_SUB)) THEN
                     IF (NSUB == 0) THEN
                        DO I=1,LSUB
                           MPFACTOR_CALC_SUB(I) = 'Y'
                           MPFACTOR_REQ_PARTFAC_SUB(I) = 'Y'
                        ENDDO
                     ELSE
                        MPFACTOR_CALC_SUB(NSUB) = 'Y'
                        MPFACTOR_REQ_PARTFAC_SUB(NSUB) = 'Y'
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

      ENDIF

      CALC_FLAG = MERGE('Y','N',CALC_FLAG == 'Y')
      IF (WHAT == 'MEFFMASS') THEN
         IF ((MEFFMASS_REQ_SUMMARY == 'Y') .OR. (MEFFMASS_REQ_MEFFM == 'Y') .OR. (MEFFMASS_REQ_MEFFW == 'Y') .OR.                 &
             (MEFFMASS_REQ_FRACSUM == 'Y')) THEN
            CALC_FLAG = 'Y'
         ENDIF
      ELSE IF (WHAT == 'MPFACTOR') THEN
         IF (MPFACTOR_REQ_PARTFAC == 'Y') THEN
            CALC_FLAG = 'Y'
         ENDIF
      ENDIF

      IF (COMPAT_STYLE .AND. (LEN_TRIM(DESCR) > 0) .AND. (.NOT. HAS_ALL)) THEN
         UNKNOWN_DESCR = .TRUE.
         IF (WHAT == 'MEFFMASS') THEN
            IF (HAS_SUMMARY .OR. HAS_PARTFAC .OR. HAS_MEFFM .OR. HAS_MEFFW .OR. HAS_FRACSUM .OR. HAS_GRID) THEN
               UNKNOWN_DESCR = .FALSE.
            ENDIF
         ELSE IF (WHAT == 'MPFACTOR') THEN
            IF (HAS_PARTFAC .OR. HAS_GRID) THEN
               UNKNOWN_DESCR = .FALSE.
            ENDIF
         ENDIF
      ENDIF

      IF (COMPAT_STYLE) THEN
         WARN_ERR = WARN_ERR + 1
         IF (WHAT == 'MEFFMASS') THEN
            WRITE(ERR,9885) CARD(1:LEN_TRIM(CARD))
         ELSE
            WRITE(ERR,9889) CARD(1:LEN_TRIM(CARD))
         ENDIF
         IF (SUPWARN == 'N') THEN
            IF (WHAT == 'MEFFMASS') THEN
               WRITE(F06,9885) CARD(1:LEN_TRIM(CARD))
            ELSE
               WRITE(F06,9889) CARD(1:LEN_TRIM(CARD))
            ENDIF
         ENDIF
         IF (UNKNOWN_DESCR) THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,9887) DESCR(1:LEN_TRIM(DESCR)), WHAT
            IF (SUPWARN == 'N') THEN
               WRITE(F06,9887) DESCR(1:LEN_TRIM(DESCR)), WHAT
            ENDIF
         ENDIF
         IF (HAS_GRID .AND. (GRID_ID < 0)) THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,9888) WHAT, DESCR(1:LEN_TRIM(DESCR))
            IF (SUPWARN == 'N') THEN
               WRITE(F06,9888) WHAT, DESCR(1:LEN_TRIM(DESCR))
            ENDIF
         ENDIF
      ELSE IF (UNKNOWN_RHS) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,9886) WHAT, RHS(1:LEN_TRIM(RHS)), WHAT
         IF (SUPWARN == 'N') THEN
            WRITE(F06,9886) WHAT, RHS(1:LEN_TRIM(RHS)), WHAT
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 9885 FORMAT(' *WARNING    : NASTRAN-STYLE MEFFMASS CASE CONTROL ENTRY "',A,'" WAS RECOGNIZED. MYSTRAN 18a CURRENTLY MAPS',       &
            ' ALL, SUMMARY, PARTFAC, MEFFM, MEFFW, FRACSUM, AND GRID= ONTO THE CURRENT OUTPUT CONTROLS.',                         &
      14X,  ' OTHER NASTRAN-STYLE OPTIONS STILL FALL BACK TO CURRENT MYSTRAN BEHAVIOR OR ARE IGNORED IN THIS RELEASE.',           &
      14X,  ' A FUTURE RELEASE IS EXPECTED TO MOVE CLOSER TO FULL NASTRAN SEMANTICS.')
 9886 FORMAT(' *WARNING    : UNRECOGNIZED OPTION "',A,'" ON ',A,' CASE CONTROL ENTRY. CURRENT MYSTRAN BEHAVIOR WILL TREAT THIS',   &
            ' AS AN ENABLED ',A,' REQUEST.')
 9887 FORMAT(' *WARNING    : DESCRIPTOR LIST "',A,'" ON ',A,' CASE CONTROL ENTRY WAS ONLY PARTIALLY MAPPED. UNSUPPORTED',          &
            ' DESCRIPTORS ARE BEING IGNORED IN THIS RELEASE.')
 9888 FORMAT(' *WARNING    : COULD NOT PARSE A NONNEGATIVE GRID ID FROM ',A,' CASE CONTROL DESCRIPTORS "',A,'". EXISTING',         &
            ' PARAM MEFMLOC/MEFMGRID SETTINGS WILL BE RETAINED.')
 9889 FORMAT(' *WARNING    : NASTRAN-STYLE MPFACTOR CASE CONTROL ENTRY "',A,'" WAS RECOGNIZED. MYSTRAN 18a CURRENTLY MAPS',       &
            ' ALL, PARTFAC, AND GRID= ONTO THE CURRENT OUTPUT CONTROLS.',                                                           &
      14X,  ' OTHER NASTRAN-STYLE OPTIONS STILL FALL BACK TO CURRENT MYSTRAN BEHAVIOR OR ARE IGNORED IN THIS RELEASE.',           &
      14X,  ' A FUTURE RELEASE IS EXPECTED TO MOVE CLOSER TO FULL NASTRAN SEMANTICS.')

! **********************************************************************************************************************************

      END SUBROUTINE CC_MPF_MEFM
