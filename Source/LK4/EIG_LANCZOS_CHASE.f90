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

      SUBROUTINE EIG_LANCZOS_CHASE

! CHASE backend surrogate for Lanczos entry point.
! Current behavior: run non-ARPACK surrogate via MGIV core so we can
! do head-to-head validation before external ChASE library integration.

      USE PENTIUM_II_KIND, ONLY       :  BYTE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, WARN_ERR
      USE PARAMS, ONLY                :  SUPINFO
      USE MODEL_STUF, ONLY            :  EIG_METH

      USE EIG_LANCZOS_CHASE_USE_IFs
      USE LINK_MESSAGE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EIG_LANCZOS_CHASE'
      CHARACTER(8*BYTE)               :: EIG_METH_SAVE = ' '

#ifdef MYSTRAN_HAVE_EXTERNAL_CHASE
      CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - CHASE NATIVE HOOK (TEMP SURROGATE CORE)')

      WARN_ERR = WARN_ERR + 1
      WRITE(ERR,4914)
      IF (SUPINFO == 'N') THEN
         WRITE(F06,4914)
      ENDIF
#else
      CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - CHASE SURROGATE (MGIV CORE)')

      WARN_ERR = WARN_ERR + 1
      WRITE(ERR,4912)
      IF (SUPINFO == 'N') THEN
         WRITE(F06,4912)
      ENDIF
#endif

      EIG_METH_SAVE = EIG_METH
      EIG_METH      = 'MGIV    '

      CALL EIG_GIV_MGIV

      EIG_METH = EIG_METH_SAVE

      RETURN

 4912 FORMAT(' *WARNING 4912: CHASE EXTERNAL BACKEND NOT LINKED. USING MGIV SURROGATE (NON-ARPACK).')
 4914 FORMAT(' *WARNING 4914: CHASE EXTERNAL BACKEND IS LINKED, BUT GENERALIZED BRIDGE IS NOT FINAL. USING MGIV SURROGATE.')

      END SUBROUTINE EIG_LANCZOS_CHASE
