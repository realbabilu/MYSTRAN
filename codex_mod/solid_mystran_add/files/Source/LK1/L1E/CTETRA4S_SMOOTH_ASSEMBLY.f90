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

      SUBROUTINE CTETRA4S_SMOOTH_ASSEMBLY ( ACTION, LTERM, NTERM )

! --- newsolid_add begin --- !
! Entry point for the CTETRA4 smooth-blend assembly path.
!
! Keep legacy CTETRA/TETRA untouched. The smooth Tet4 candidate from Python is
! a nodal-patch/global assembly formulation, not a local element matrix tweak,
! so it belongs beside the stiffness assembly processor instead of in TETRA.f90.
!
! ACTION = 'COUNT' estimates/records the additional STF topology needed by the
!          nodal-patch smooth terms.
! ACTION = 'ADD  ' adds the alpha*(K_smooth - K_standard) correction into STF.
! --- newsolid_add end --- !

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NELE, NGRID
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO
      USE PARAMS, ONLY                :  EPSIL, SOLIDTYP, SPARSTOR
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  EDAT, EPNT, ETYPE, GRID_ID, MATL, PSOLID, RGRID, RMATL
      USE STF_ARRAYS, ONLY            :  STFKEY, STF3

      USE GET_ARRAY_ROW_NUM_Interface
      USE TDOF_COL_NUM_Interface

      IMPLICIT NONE

      CHARACTER( 5*BYTE), INTENT(IN)  :: ACTION
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CTETRA4S_SMOOTH_ASSEMBLY'

      INTEGER(LONG), INTENT(IN)       :: LTERM
      INTEGER(LONG), INTENT(INOUT)    :: NTERM

      INTEGER(LONG)                   :: CENTER_GRID
      INTEGER(LONG)                   :: I
      INTEGER(LONG), ALLOCATABLE      :: ADJ_ELEMS(:)

      REAL(DOUBLE)                    :: EPS1

      INTRINSIC                       :: DABS

      SUBR_NAME = SUBR_NAME
      IF ((ACTION /= 'COUNT') .AND. (ACTION /= 'ADD  ')) THEN
         RETURN
      ENDIF
      IF (LTERM < 0) THEN
         RETURN
      ENDIF
      IF (SOLIDTYP /= 'NEWSOLID') THEN
         RETURN
      ENDIF
      EPS1 = EPSIL(1)

      ALLOCATE(ADJ_ELEMS(NELE))

      DO I=1,NGRID
         CENTER_GRID = GRID_ID(I)
         CALL PROCESS_GRID_PATCH ( CENTER_GRID, ADJ_ELEMS )
      ENDDO

      DEALLOCATE(ADJ_ELEMS)
      RETURN

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE ADD_STF_TERM ( KGG_ROW_IN, KGG_COL_IN, VALUE )

! Adds or counts one term in the LINK1 STF linked-list matrix.

      INTEGER(LONG), INTENT(IN)       :: KGG_ROW_IN
      INTEGER(LONG), INTENT(IN)       :: KGG_COL_IN
      REAL(DOUBLE),  INTENT(IN)       :: VALUE

      INTEGER(LONG)                   :: IDUM
      INTEGER(LONG)                   :: IS
      INTEGER(LONG)                   :: ISS
      INTEGER(LONG)                   :: KGG_COL
      INTEGER(LONG)                   :: KGG_ROW

      IF ((ACTION == 'ADD  ') .AND. (DABS(VALUE) < EPS1)) THEN
         RETURN
      ENDIF

      KGG_ROW = KGG_ROW_IN
      KGG_COL = KGG_COL_IN

      IF (SPARSTOR == 'SYM') THEN
         IF (KGG_COL < KGG_ROW) THEN
            IDUM    = KGG_ROW
            KGG_ROW = KGG_COL
            KGG_COL = IDUM
         ENDIF
      ENDIF

      IS = STFKEY(KGG_ROW)

      IF (IS == 0) THEN
         NTERM = NTERM + 1
         IF (NTERM > LTERM) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1624) SUBR_NAME, 'STIFFNESS','LTERM', LTERM
            WRITE(F06,1624) SUBR_NAME, 'STIFFNESS','LTERM', LTERM
            RETURN
         ENDIF
         STFKEY(KGG_ROW)   = NTERM
         STF3(NTERM)%Col_1 = KGG_COL
         STF3(NTERM)%Col_2 = 0
         IF (ACTION == 'ADD  ') THEN
            STF3(NTERM)%Col_3 = VALUE
         ELSE
            STF3(NTERM)%Col_3 = ZERO
         ENDIF
         RETURN
      ENDIF

stfpnt0: DO
         IF (KGG_COL == STF3(IS)%Col_1) THEN
            IF (ACTION == 'ADD  ') THEN
               STF3(IS)%Col_3 = STF3(IS)%Col_3 + VALUE
            ENDIF
            RETURN
         ENDIF

         ISS = IS
         IS  = STF3(IS)%Col_2
         IF (IS == 0) THEN
            NTERM = NTERM + 1
            IF (NTERM > LTERM) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1624) SUBR_NAME, 'STIFFNESS','LTERM', LTERM
               WRITE(F06,1624) SUBR_NAME, 'STIFFNESS','LTERM', LTERM
               RETURN
            ENDIF
            STF3(ISS)%Col_2   = NTERM
            STF3(NTERM)%Col_1 = KGG_COL
            STF3(NTERM)%Col_2 = 0
            IF (ACTION == 'ADD  ') THEN
               STF3(NTERM)%Col_3 = VALUE
            ELSE
               STF3(NTERM)%Col_3 = ZERO
            ENDIF
            RETURN
         ENDIF
      ENDDO stfpnt0

 1624 FORMAT(' *ERROR  1624: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A,' MATRIX ARRAY SIZE BASED ON ',A,' = ',I12,' IS TOO SMALL')

      END SUBROUTINE ADD_STF_TERM

! ##################################################################################################################################

      SUBROUTINE PROCESS_GRID_PATCH ( CENTER_GRID, ADJ_ELEMS )

! Assembles the alpha*(K_smooth - K_standard) correction for one CTETRA4 nodal patch.

      INTEGER(LONG), INTENT(IN)       :: CENTER_GRID
      INTEGER(LONG), INTENT(INOUT)    :: ADJ_ELEMS(NELE)

      INTEGER(LONG)                   :: ADJ_COUNT
      INTEGER(LONG)                   :: DOF_COUNT
      INTEGER(LONG)                   :: E
      INTEGER(LONG)                   :: GDOF(12)
      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: J
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: KSTART
      INTEGER(LONG)                   :: L
      INTEGER(LONG)                   :: MAT_ROW
      INTEGER(LONG)                   :: PATCH_COL
      INTEGER(LONG), ALLOCATABLE      :: PATCH_GDOF(:)

      REAL(DOUBLE), PARAMETER         :: ALPHA = 0.9D0
      REAL(DOUBLE)                    :: B(6,12)
      REAL(DOUBLE), ALLOCATABLE       :: BBAR(:,:)
      REAL(DOUBLE)                    :: D(6,6)
      REAL(DOUBLE), ALLOCATABLE       :: DBBAR(:,:)
      REAL(DOUBLE)                    :: DBLOC(6,12)
      REAL(DOUBLE)                    :: KVAL
      REAL(DOUBLE)                    :: TOTAL_VOLUME
      REAL(DOUBLE)                    :: VOLUME

      CALL FIND_CTETRA4_PATCH ( CENTER_GRID, ADJ_ELEMS, ADJ_COUNT )
      IF (ADJ_COUNT <= 0) THEN
         RETURN
      ENDIF

      CALL GET_ELEM_MAT_ROW ( ADJ_ELEMS(1), MAT_ROW )
      DO I=2,ADJ_COUNT
         CALL GET_ELEM_MAT_ROW ( ADJ_ELEMS(I), J )
         IF (J /= MAT_ROW) THEN
            RETURN
         ENDIF
      ENDDO

      IF (.NOT. BUILD_MATERIAL_D ( MAT_ROW, D )) THEN
         RETURN
      ENDIF

      ALLOCATE(PATCH_GDOF(12*ADJ_COUNT))
      ALLOCATE(BBAR(6,12*ADJ_COUNT))
      ALLOCATE(DBBAR(6,12*ADJ_COUNT))

      DOF_COUNT    = 0
      TOTAL_VOLUME = ZERO
      BBAR(:,:)    = ZERO

      DO E=1,ADJ_COUNT
         CALL CTETRA4_B_GLOBAL_FROM_ELEM ( ADJ_ELEMS(E), VOLUME, B, GDOF )
         IF (VOLUME <= EPS1) THEN
            CYCLE
         ENDIF

         TOTAL_VOLUME = TOTAL_VOLUME + VOLUME

         DO J=1,12
            PATCH_COL = PATCH_DOF_COL ( GDOF(J), PATCH_GDOF, DOF_COUNT )
            DO I=1,6
               BBAR(I,PATCH_COL) = BBAR(I,PATCH_COL) + VOLUME*B(I,J)
            ENDDO
         ENDDO

         DO J=1,12
            DO I=1,6
               DBLOC(I,J) = ZERO
               DO K=1,6
                  DBLOC(I,J) = DBLOC(I,J) + D(I,K)*B(K,J)
               ENDDO
            ENDDO
         ENDDO

         DO J=1,12
            IF (SPARSTOR == 'SYM') THEN
               KSTART = J
            ELSE
               KSTART = 1
            ENDIF
            DO K=KSTART,12
               KVAL = ZERO
               DO L=1,6
                  KVAL = KVAL + B(L,J)*DBLOC(L,K)
               ENDDO
               KVAL = -ALPHA*(VOLUME/4.0D0)*KVAL
               CALL ADD_STF_TERM ( GDOF(J), GDOF(K), KVAL )
            ENDDO
         ENDDO
      ENDDO

      IF (TOTAL_VOLUME > EPS1) THEN
         DO J=1,DOF_COUNT
            DO I=1,6
               BBAR(I,J) = BBAR(I,J)/TOTAL_VOLUME
            ENDDO
         ENDDO

         DO J=1,DOF_COUNT
            DO I=1,6
               DBBAR(I,J) = ZERO
               DO K=1,6
                  DBBAR(I,J) = DBBAR(I,J) + D(I,K)*BBAR(K,J)
               ENDDO
            ENDDO
         ENDDO

         DO J=1,DOF_COUNT
            IF (SPARSTOR == 'SYM') THEN
               KSTART = J
            ELSE
               KSTART = 1
            ENDIF
            DO K=KSTART,DOF_COUNT
               KVAL = ZERO
               DO L=1,6
                  KVAL = KVAL + BBAR(L,J)*DBBAR(L,K)
               ENDDO
               KVAL = ALPHA*(TOTAL_VOLUME/4.0D0)*KVAL
               CALL ADD_STF_TERM ( PATCH_GDOF(J), PATCH_GDOF(K), KVAL )
            ENDDO
         ENDDO
      ENDIF

      DEALLOCATE(BBAR)
      DEALLOCATE(DBBAR)
      DEALLOCATE(PATCH_GDOF)

      END SUBROUTINE PROCESS_GRID_PATCH

! ##################################################################################################################################

      SUBROUTINE FIND_CTETRA4_PATCH ( CENTER_GRID, ADJ_ELEMS, ADJ_COUNT )

      INTEGER(LONG), INTENT(IN)       :: CENTER_GRID
      INTEGER(LONG), INTENT(OUT)      :: ADJ_COUNT
      INTEGER(LONG), INTENT(OUT)      :: ADJ_ELEMS(NELE)

      INTEGER(LONG)                   :: E
      INTEGER(LONG)                   :: G
      INTEGER(LONG)                   :: GRIDS(4)

      ADJ_COUNT = 0
      DO E=1,NELE
         IF (ETYPE(E) /= 'TETRA4  ') THEN
            CYCLE
         ENDIF
         CALL GET_CTETRA4_GRIDS ( E, GRIDS )
         DO G=1,4
            IF (GRIDS(G) == CENTER_GRID) THEN
               ADJ_COUNT = ADJ_COUNT + 1
               ADJ_ELEMS(ADJ_COUNT) = E
               EXIT
            ENDIF
         ENDDO
      ENDDO

      END SUBROUTINE FIND_CTETRA4_PATCH

! ##################################################################################################################################

      INTEGER(LONG) FUNCTION PATCH_DOF_COL ( GDOF_IN, PATCH_GDOF, DOF_COUNT )

      INTEGER(LONG), INTENT(IN)       :: GDOF_IN
      INTEGER(LONG), INTENT(INOUT)    :: DOF_COUNT
      INTEGER(LONG), INTENT(INOUT)    :: PATCH_GDOF(:)

      INTEGER(LONG)                   :: I

      DO I=1,DOF_COUNT
         IF (PATCH_GDOF(I) == GDOF_IN) THEN
            PATCH_DOF_COL = I
            RETURN
         ENDIF
      ENDDO

      DOF_COUNT = DOF_COUNT + 1
      PATCH_GDOF(DOF_COUNT) = GDOF_IN
      PATCH_DOF_COL = DOF_COUNT

      END FUNCTION PATCH_DOF_COL

! ##################################################################################################################################

      SUBROUTINE GET_CTETRA4_GRIDS ( INT_ELEM_ID, ACTUAL_GRIDS )

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID
      INTEGER(LONG), INTENT(OUT)      :: ACTUAL_GRIDS(4)

      INTEGER(LONG)                   :: EPNTK

      EPNTK = EPNT(INT_ELEM_ID)
      ACTUAL_GRIDS(1) = EDAT(EPNTK+2)
      ACTUAL_GRIDS(2) = EDAT(EPNTK+3)
      ACTUAL_GRIDS(3) = EDAT(EPNTK+4)
      ACTUAL_GRIDS(4) = EDAT(EPNTK+5)

      END SUBROUTINE GET_CTETRA4_GRIDS

! ##################################################################################################################################

      SUBROUTINE GET_ELEM_MAT_ROW ( INT_ELEM_ID, MAT_ROW )

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID
      INTEGER(LONG), INTENT(OUT)      :: MAT_ROW

      INTEGER(LONG)                   :: EPNTK
      INTEGER(LONG)                   :: PID_ROW

      EPNTK = EPNT(INT_ELEM_ID)
      PID_ROW = EDAT(EPNTK+1)
      MAT_ROW = PSOLID(PID_ROW,2)

      END SUBROUTINE GET_ELEM_MAT_ROW

! ##################################################################################################################################

      LOGICAL FUNCTION BUILD_MATERIAL_D ( MAT_ROW, D )

      INTEGER(LONG), INTENT(IN)       :: MAT_ROW
      REAL(DOUBLE),  INTENT(OUT)      :: D(6,6)

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: J
      INTEGER(LONG)                   :: K
      REAL(DOUBLE)                    :: DEN1
      REAL(DOUBLE)                    :: DEN2
      REAL(DOUBLE)                    :: E
      REAL(DOUBLE)                    :: E0
      REAL(DOUBLE)                    :: G
      REAL(DOUBLE)                    :: NU
      REAL(DOUBLE)                    :: NU0

      D(:,:) = ZERO
      BUILD_MATERIAL_D = .FALSE.

      IF (MATL(MAT_ROW,2) == 1) THEN
         E  = RMATL(MAT_ROW,1)
         G  = RMATL(MAT_ROW,2)
         NU = RMATL(MAT_ROW,3)

         DEN1 = ONE + NU
         DEN2 = ONE - TWO*NU
         IF ((DABS(DEN1) <= EPS1) .OR. (DABS(DEN2) <= EPS1)) THEN
            RETURN
         ENDIF

         E0  = E/DEN1
         NU0 = (ONE - NU)/DEN2

         D(1,1) = E0*NU0
         D(1,2) = E0*NU/DEN2
         D(1,3) = D(1,2)
         D(2,1) = D(1,2)
         D(2,2) = D(1,1)
         D(2,3) = D(1,2)
         D(3,1) = D(1,2)
         D(3,2) = D(1,2)
         D(3,3) = D(1,1)
         D(4,4) = G
         D(5,5) = G
         D(6,6) = G

         BUILD_MATERIAL_D = .TRUE.

      ELSE IF (MATL(MAT_ROW,2) == 9) THEN
         K = 0
         DO I=1,6
            DO J=I,6
               K = K + 1
               D(I,J) = RMATL(MAT_ROW,K)
               D(J,I) = D(I,J)
            ENDDO
         ENDDO
         BUILD_MATERIAL_D = .TRUE.
      ENDIF

      END FUNCTION BUILD_MATERIAL_D

! ##################################################################################################################################

      SUBROUTINE CTETRA4_B_GLOBAL_FROM_ELEM ( INT_ELEM_ID, VOLUME, BGLOBAL, GDOF )

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID
      INTEGER(LONG), INTENT(OUT)      :: GDOF(12)
      REAL(DOUBLE),  INTENT(OUT)      :: BGLOBAL(6,12)
      REAL(DOUBLE),  INTENT(OUT)      :: VOLUME

      INTEGER(LONG)                   :: ACTUAL_GRIDS(4)

      CALL GET_CTETRA4_GRIDS ( INT_ELEM_ID, ACTUAL_GRIDS )
      CALL CTETRA4_B_GLOBAL ( ACTUAL_GRIDS, VOLUME, BGLOBAL, GDOF )

      END SUBROUTINE CTETRA4_B_GLOBAL_FROM_ELEM

! ##################################################################################################################################

      SUBROUTINE GET_GRID_TRANSLATION_GDOFS ( ACTUAL_GRID, GDOF )

! Returns G-set DOF numbers for components 1, 2, 3 of one physical grid.

      INTEGER(LONG), INTENT(IN)       :: ACTUAL_GRID
      INTEGER(LONG), INTENT(OUT)      :: GDOF(3)

      INTEGER(LONG)                   :: G_SET_COL_NUM
      INTEGER(LONG)                   :: IGRID
      INTEGER(LONG)                   :: ROW_NUM_START

      CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, ACTUAL_GRID, IGRID )
      ROW_NUM_START = TDOF_ROW_START(IGRID)
      CALL TDOF_COL_NUM ( 'G ',  G_SET_COL_NUM )

      GDOF(1) = TDOF(ROW_NUM_START    , G_SET_COL_NUM)
      GDOF(2) = TDOF(ROW_NUM_START + 1, G_SET_COL_NUM)
      GDOF(3) = TDOF(ROW_NUM_START + 2, G_SET_COL_NUM)

      END SUBROUTINE GET_GRID_TRANSLATION_GDOFS

! ##################################################################################################################################

      SUBROUTINE GET_GRID_BASIC_COORDS ( ACTUAL_GRID, XYZ )

! Returns basic-system coordinates for one grid.

      INTEGER(LONG), INTENT(IN)       :: ACTUAL_GRID
      REAL(DOUBLE),  INTENT(OUT)      :: XYZ(3)

      INTEGER(LONG)                   :: IGRID

      CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, ACTUAL_GRID, IGRID )

      XYZ(1) = RGRID(IGRID,1)
      XYZ(2) = RGRID(IGRID,2)
      XYZ(3) = RGRID(IGRID,3)

      END SUBROUTINE GET_GRID_BASIC_COORDS

! ##################################################################################################################################

      SUBROUTINE CTETRA4_B_GLOBAL ( ACTUAL_GRIDS, VOLUME, BGLOBAL, GDOF )

! Builds the constant-strain CTETRA4 B matrix in the basic/global Cartesian basis.

      INTEGER(LONG), INTENT(IN)       :: ACTUAL_GRIDS(4)
      INTEGER(LONG), INTENT(OUT)      :: GDOF(12)
      REAL(DOUBLE),  INTENT(OUT)      :: BGLOBAL(6,12)
      REAL(DOUBLE),  INTENT(OUT)      :: VOLUME

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: J
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: NODE_GDOF(3)
      REAL(DOUBLE)                    :: A(3,3)
      REAL(DOUBLE)                    :: DET
      REAL(DOUBLE)                    :: GRAD(3,4)
      REAL(DOUBLE)                    :: INV_A(3,3)
      REAL(DOUBLE)                    :: XYZ(3,4)

      DO I=1,4
         CALL GET_GRID_BASIC_COORDS ( ACTUAL_GRIDS(I), XYZ(1:3,I) )
         CALL GET_GRID_TRANSLATION_GDOFS ( ACTUAL_GRIDS(I), NODE_GDOF )
         GDOF(3*(I-1)+1) = NODE_GDOF(1)
         GDOF(3*(I-1)+2) = NODE_GDOF(2)
         GDOF(3*(I-1)+3) = NODE_GDOF(3)
      ENDDO

      DO I=1,3
         A(I,1) = XYZ(I,2) - XYZ(I,1)
         A(I,2) = XYZ(I,3) - XYZ(I,1)
         A(I,3) = XYZ(I,4) - XYZ(I,1)
      ENDDO

      CALL INVERT_3X3 ( A, DET, INV_A )
      VOLUME = DABS(DET)/6.0D0

      IF (VOLUME <= EPS1) THEN
         BGLOBAL(:,:) = ZERO
         RETURN
      ENDIF

      DO I=1,3
         GRAD(I,2) = INV_A(1,I)
         GRAD(I,3) = INV_A(2,I)
         GRAD(I,4) = INV_A(3,I)
         GRAD(I,1) = -GRAD(I,2) - GRAD(I,3) - GRAD(I,4)
      ENDDO

      BGLOBAL(:,:) = ZERO
      DO I=1,4
         J = 3*(I-1)
         BGLOBAL(1,J+1) = GRAD(1,I)
         BGLOBAL(2,J+2) = GRAD(2,I)
         BGLOBAL(3,J+3) = GRAD(3,I)
         BGLOBAL(4,J+1) = GRAD(2,I)
         BGLOBAL(4,J+2) = GRAD(1,I)
         BGLOBAL(5,J+2) = GRAD(3,I)
         BGLOBAL(5,J+3) = GRAD(2,I)
         BGLOBAL(6,J+1) = GRAD(3,I)
         BGLOBAL(6,J+3) = GRAD(1,I)
      ENDDO

      END SUBROUTINE CTETRA4_B_GLOBAL

! ##################################################################################################################################

      SUBROUTINE INVERT_3X3 ( A, DET, AINV )

      REAL(DOUBLE), INTENT(IN)        :: A(3,3)
      REAL(DOUBLE), INTENT(OUT)       :: AINV(3,3)
      REAL(DOUBLE), INTENT(OUT)       :: DET

      DET = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2))                                                                       &
          - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1))                                                                       &
          + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))

      IF (DABS(DET) <= EPS1) THEN
         AINV(:,:) = ZERO
         RETURN
      ENDIF

      AINV(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))/DET
      AINV(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))/DET
      AINV(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2))/DET
      AINV(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))/DET
      AINV(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))/DET
      AINV(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))/DET
      AINV(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1))/DET
      AINV(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))/DET
      AINV(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))/DET

      END SUBROUTINE INVERT_3X3

! ##################################################################################################################################

      END SUBROUTINE CTETRA4S_SMOOTH_ASSEMBLY
