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

      MODULE TREL1_USE_IFs

! USE Interface statements for all subroutines called by SUBROUTINE TREL1

      USE OURTIM_Interface
      USE TMEM1_Interface
      USE TPLT1_Interface
      USE TPLT2_Interface
      USE TPLT_DKMT_Interface
      USE outa_here_Interface
      USE MATMULT_FFF_Interface
      USE MATMULT_FFF_T_Interface

! --- cquadr/ctriar begin --- !
      INTERFACE
         SUBROUTINE TMEM_ALLT3 ( OPT, AREA, WRT_BUG_THIS_TIME, BIG_BM )
            USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
            USE SCONTR, ONLY                :  NSUB, NTSUB
            USE MODEL_STUF, ONLY            :  ELDOF
            IMPLICIT NONE
            CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
            CHARACTER(1*BYTE), INTENT(IN)   :: WRT_BUG_THIS_TIME
            REAL(DOUBLE)   , INTENT(IN)     :: AREA
            REAL(DOUBLE)   , INTENT(OUT)    :: BIG_BM(3,ELDOF,1)
         END SUBROUTINE TMEM_ALLT3
         SUBROUTINE TMEM_T3FREE ( OPT, AREA, BIG_BM )
            USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
            USE MODEL_STUF, ONLY            :  ELDOF
            IMPLICIT NONE
            CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
            REAL(DOUBLE), INTENT(IN)        :: AREA
            REAL(DOUBLE), INTENT(OUT)       :: BIG_BM(3,ELDOF,1)
         END SUBROUTINE TMEM_T3FREE
! --- MITC3+_add begin --- !
         SUBROUTINE TPLT_MITC3P ( OPT, AREA, X2E, X3E, Y3E, BIG_BB )
            USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
            USE MODEL_STUF, ONLY            :  ELDOF
            IMPLICIT NONE
            CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
            REAL(DOUBLE), INTENT(IN)        :: AREA
            REAL(DOUBLE), INTENT(IN)        :: X2E
            REAL(DOUBLE), INTENT(IN)        :: X3E
            REAL(DOUBLE), INTENT(IN)        :: Y3E
            REAL(DOUBLE), INTENT(OUT)       :: BIG_BB(3,ELDOF,1)
         END SUBROUTINE TPLT_MITC3P
! --- MITC3+_add end --- !
      END INTERFACE
! --- cquadr/ctriar end --- !

      END MODULE TREL1_USE_IFs
