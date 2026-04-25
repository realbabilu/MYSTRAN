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

! !--- CHASE and FEAST --- begin!
      SUBROUTINE EIG_LANCZOS_FEAST

! Stage-1 FEAST entry point.
! The method selector is integrated, but if the native backend is not linked
! this wrapper falls back to ARPACK Lanczos.

      USE PENTIUM_II_KIND, ONLY       :  BYTE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, WARN_ERR
      USE PARAMS, ONLY                :  SUPINFO

      USE EIG_LANCZOS_FEAST_USE_IFs
      USE LINK_MESSAGE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EIG_LANCZOS_FEAST'

      WARN_ERR = WARN_ERR + 1
      WRITE(ERR,4911)
      IF (SUPINFO == 'N') WRITE(F06,4911)

      CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - FEAST FALLBACK (ARPACK LANCZOS)')
      CALL EIG_LANCZOS_ARPACK

      RETURN

 4911 FORMAT(' *WARNING 4911: FEAST NATIVE BACKEND NOT LINKED IN THIS BUILD. USING ARPACK LANCZOS FALLBACK.')

      END SUBROUTINE EIG_LANCZOS_FEAST
! !--- CHASE and FEAST --- end!
