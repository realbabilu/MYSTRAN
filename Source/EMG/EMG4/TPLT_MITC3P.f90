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

      SUBROUTINE TPLT_MITC3P ( OPT, AREA, X2E, X3E, Y3E, CALC_EMATS, IERROR, KV, PTV, PPV, B2V, B3V, S2V, S3V, BIG_BB,            &
                               MN4T_QD, TRIA_NUM, PSI )

! --- MITC3+_add begin --- !
! First MYSTRAN integration point for PARAM,TRIA3TYP,MITC3+.
!
! This keeps the legacy CTRIA3 data contract intact and routes the plate branch
! through the existing MIN3 recovery/load machinery while adding explicit MITC3+
! identity and drilling rotational stiffness. The full bubble-condensed MITC3+
! plate kernel can replace the TPLT2 delegate here without changing parser or
! TREL1 routing.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NSUB, NTSUB
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  KE

      USE TPLT2_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TPLT_MITC3P'
      CHARACTER(1*BYTE), INTENT(IN)   :: CALC_EMATS
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      CHARACTER(LEN=*) , INTENT(IN)   :: MN4T_QD

      INTEGER(LONG), INTENT(OUT)      :: IERROR
      INTEGER(LONG), INTENT(IN)       :: TRIA_NUM
      INTEGER(LONG)                   :: I
      INTEGER(LONG), PARAMETER        :: RZ_DOF(3) = (/ 6, 12, 18 /)

      REAL(DOUBLE) , INTENT(IN)       :: AREA
      REAL(DOUBLE) , INTENT(IN)       :: PSI
      REAL(DOUBLE) , INTENT(IN)       :: X2E
      REAL(DOUBLE) , INTENT(IN)       :: X3E
      REAL(DOUBLE) , INTENT(IN)       :: Y3E
      REAL(DOUBLE) , INTENT(OUT)      :: BIG_BB(3,18,1)
      REAL(DOUBLE) , INTENT(OUT)      :: B2V(3,9)
      REAL(DOUBLE) , INTENT(OUT)      :: B3V(3,9)
      REAL(DOUBLE) , INTENT(OUT)      :: KV(9,9)
      REAL(DOUBLE) , INTENT(OUT)      :: PPV(9,NSUB)
      REAL(DOUBLE) , INTENT(OUT)      :: PTV(9,NTSUB)
      REAL(DOUBLE) , INTENT(OUT)      :: S2V(3,9)
      REAL(DOUBLE) , INTENT(OUT)      :: S3V(3,9)

      REAL(DOUBLE)                    :: DRILL_PEN
      REAL(DOUBLE)                    :: KTRACE

      CALL TPLT2 ( OPT, AREA, X2E, X3E, Y3E, CALC_EMATS, IERROR, KV, PTV, PPV, B2V, B3V, S2V, S3V, BIG_BB, MN4T_QD, TRIA_NUM, PSI )

      IF (OPT(4) == 'Y') THEN
         KTRACE = ZERO
         DO I=1,18
            KTRACE = KTRACE + DABS(KE(I,I))
         ENDDO
         DRILL_PEN = 1.0D-8*KTRACE/3.0D0
         DO I=1,3
            KE(RZ_DOF(I),RZ_DOF(I)) = KE(RZ_DOF(I),RZ_DOF(I)) + DRILL_PEN
         ENDDO
      ENDIF

      RETURN
! --- MITC3+_add end --- !

      END SUBROUTINE TPLT_MITC3P
