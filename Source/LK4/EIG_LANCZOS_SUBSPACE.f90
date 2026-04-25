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

      SUBROUTINE EIG_LANCZOS_SUBSPACE

! Experimental subspace entry point.
! Stage-1 integration maps Subspace controls to ARPACK Lanczos backend.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, WARN_ERR
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  ARP_TOL, MXITERL, SUBSPITR, SUBSPTOL, SUBSPMAX, SUPINFO
      USE MODEL_STUF, ONLY            :  EIG_FRQ1, EIG_FRQ2, EIG_N2

      USE EIG_LANCZOS_SUBSPACE_USE_IFs
      USE LINK_MESSAGE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EIG_LANCZOS_SUBSPACE'

      INTEGER(LONG)                   :: MXITERL_SAVE
      INTEGER(LONG)                   :: NEV_REQ
      REAL(DOUBLE)                    :: ARP_TOL_SAVE

      MXITERL_SAVE = MXITERL
      ARP_TOL_SAVE = ARP_TOL

      MXITERL = SUBSPITR
      ARP_TOL = SUBSPTOL

      NEV_REQ = EIG_N2
      IF ((EIG_N2 <= 1) .AND. (EIG_FRQ1 <= ZERO) .AND. (EIG_FRQ2 <= ZERO)) THEN
         NEV_REQ = SUBSPMAX
      ENDIF
      IF (NEV_REQ < 1) NEV_REQ = SUBSPMAX
      EIG_N2 = NEV_REQ

      WARN_ERR = WARN_ERR + 1
      WRITE(ERR,4921)
      IF (SUPINFO == 'N') WRITE(F06,4921)

      WRITE(ERR,4922) SUBSPITR, SUBSPTOL, SUBSPMAX, EIG_N2
      IF (SUPINFO == 'N') WRITE(F06,4922) SUBSPITR, SUBSPTOL, SUBSPMAX, EIG_N2

      CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - SUBSPACE (ARPACK-MAPPED STAGE-1)')
      CALL EIG_LANCZOS_ARPACK

      MXITERL = MXITERL_SAVE
      ARP_TOL = ARP_TOL_SAVE

      RETURN

 4921 FORMAT(' *WARNING 4921: LANCMETH=SUBSP IS CURRENTLY STAGE-1 MAPPED TO ARPACK BACKEND. ',                                &
                    'NATIVE SUBSPACE ITERATION SOLVER IS NOT WIRED YET.')

 4922 FORMAT(' *INFORMATION: SUBSPACE CONTROL ACTIVE: SUBSPITR=',I8,', SUBSPTOL=',1ES11.4,', SUBSPMAX=',I8,                 &
                    ', EFFECTIVE_ND=',I8)

      END SUBROUTINE EIG_LANCZOS_SUBSPACE
