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

      MODULE EMG_USE_IFs

! USE Interface statements for all subroutines called by SUBROUTINE EMG

      USE IS_ELEM_PCOMP_PROPS_Interface
      USE OURTIM_Interface
      USE ELMDAT1_Interface
      USE OUTA_HERE_Interface
      USE ELMGM1_Interface
      USE ELMGM2_Interface
      USE ELMGM3_Interface
      USE GET_MATANGLE_FROM_CID_Interface
      USE MATERIAL_PROPS_2D_Interface
      USE ROT_AXES_MATL_TO_LOC_Interface
      USE MATERIAL_PROPS_3D_Interface
      USE ELMOUT_Interface
      USE SHELL_ABD_MATRICES_Interface
      USE ELMDAT2_Interface
      USE ELAS1_Interface
      USE BREL1_Interface
      USE BUSH_Interface
      USE TREL1_Interface
      USE CTRIA3_T3FF_Interface
      USE CTRIA6_MITC6_Interface
      USE CTRIA6_MH6T_Interface
      USE CTRIA6_REZAIEE_Interface
      USE CTRIA6_SIMO1993_Interface
      USE CTRIAR_DKMT18_Interface
      USE CTRIAR_MITC3PHB_Interface
      USE CTRIAR_T3FFD_Interface
       USE QDEL1_Interface
       USE CQUAD4_DSQK_RHR_Interface
       USE CQUAD4_DKMQ20_RHR_Interface
       USE CQUAD4_SIMO1989_Interface
       USE CQUADR_DKM24AU_Interface
       USE CQUADR_DKM24EA_Interface
       USE CQUADR_DKMQ24R_Interface
       USE CQUADR_DKMQ24_Interface
       USE CQUADR_DKMQ24N_Interface
       USE CQUADR_HW20_Interface
       USE CQUADR_MITC4PHB_B_Interface
       USE CQUADR_MITC4PHB_Interface
       USE CQUADR_MBP1C0_Interface
       USE CQUADR_Q4EASANS_Interface
       USE CQUADR_Q4RS_Interface
       USE CQUADR_SIMO1993_Interface
       USE HEXA_Interface
      USE PENTA_Interface
      USE PYRAM_Interface
      USE TETRA_Interface
      USE KUSER1_Interface
      USE USERIN_Interface
      USE ELMOFF_Interface
      USE MITC8_Interface
      USE CQUAD8_MITC8D_Interface
      USE CQUAD8_SIMOQ8_Interface
      USE CQUAD8_HBQ8_Interface
      USE CQUAD8_ANS8BDG6_Interface
      USE CQUAD8_MACQ8D_Interface

      END MODULE EMG_USE_IFs
