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

      MODULE GPSTRESS_SURFACE_UTILS

      USE PENTIUM_II_KIND, ONLY         :  BYTE, LONG
      USE IOUNT1, ONLY                  :  ERR, F06
      USE SCONTR, ONLY                  :  BLNK_SUB_NAM, CC_ENTRY_LEN, FATAL_ERR, LSETLN, MAX_TOKEN_LEN, SETLEN
      USE MODEL_STUF, ONLY              :  ALL_SETS_ARRAY, EPNT, ESORT1, ETYPE, EDAT
      USE CC_OUTPUT_DESCRIBERS, ONLY    :  GP_SURFACE_IDS, GP_SURFACE_NORMAL_MODE, GP_SURFACE_SETIDS, NUM_GP_SURFACE, &
                                           GP_POST_SET_IDS, GP_POST_SET_TEXT, NUM_GP_POST_SET, REGISTER_GP_SURFACE,   &
                                           REGISTER_GP_POST_SET, GP_VOLUME_DEFINITION, GP_VOLUMES, NUM_GP_VOLUME

      USE GET_ARRAY_ROW_NUM_Interface
      USE STOKEN_Interface
      USE TOKCHK_Interface

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE PARSE_GPSTRESS_VOLUME_CARD(CARD)
      CHARACTER(LEN=*), INTENT(IN) :: CARD
      CHARACTER(LEN=CC_ENTRY_LEN) :: LINE
      CHARACTER(32) :: TOK(12)
      TYPE(GP_VOLUME_DEFINITION), ALLOCATABLE :: TEMP(:)
      INTEGER(LONG) :: I, J, N, IOS, VID, SID, CID, START, FINISH
      LOGICAL :: DIRECT, PRINCIPAL, SPECIFIED
      LINE = ADJUSTL(CARD)
      CALL GPSTRESS_TO_UPPER_LINE(LINE)
      I = INDEX(LINE,'$')
      IF (I > 0) LINE(I:) = ' '
      DO I=1,LEN_TRIM(LINE)
         IF (LINE(I:I) == ',') LINE(I:I) = ' '
      ENDDO
      TOK = ' '
      N = 0
      START = 1
      DO WHILE (START <= LEN_TRIM(LINE))
         IF (LINE(START:START) == ' ') THEN
            START = START + 1
            CYCLE
         ENDIF
         FINISH = INDEX(LINE(START:),' ')
         IF (FINISH == 0) FINISH = LEN(LINE)-START+2
         FINISH = START+FINISH-2
         N = N+1
         IF (N > SIZE(TOK)) GOTO 900
         TOK(N) = LINE(START:FINISH)
         START = FINISH+2
      ENDDO
      IF (N < 4) GOTO 900
      IF (TOK(1) /= 'VOLUME' .OR. TOK(3) /= 'SET') GOTO 900
      READ(TOK(2),*,IOSTAT=IOS) VID
      IF (IOS /= 0) GOTO 900
      READ(TOK(4),*,IOSTAT=IOS) SID
      IF (IOS /= 0) GOTO 900
      IF (VID <= 0 .OR. SID <= 0) GOTO 900
      CID = 0
      DIRECT=.FALSE.; PRINCIPAL=.FALSE.; SPECIFIED=.FALSE.
      I = 5
      DO WHILE (I <= N)
         SELECT CASE (TRIM(TOK(I)))
         CASE ('DIRECT')
            DIRECT=.TRUE.; SPECIFIED=.TRUE.
            I = I+1
         CASE ('PRINCIPAL')
            PRINCIPAL=.TRUE.; SPECIFIED=.TRUE.
            I = I+1
         CASE ('SYSTEM')
            I = I+1
            IF (I > N) GOTO 900
            IF (TOK(I) == 'BASIC') THEN
               CID = 0
            ELSE IF (TOK(I) == 'CORD') THEN
               I = I+1
               IF (I > N) GOTO 900
               READ(TOK(I),*,IOSTAT=IOS) CID
               IF (IOS /= 0 .OR. CID < 0) GOTO 900
            ELSE
               GOTO 900
            ENDIF
            I = I+1
         CASE DEFAULT
            GOTO 900
         END SELECT
      ENDDO
      DO J=1,NUM_GP_VOLUME
         IF (GP_VOLUMES(J)%ID == VID) GOTO 900
      ENDDO
      ALLOCATE(TEMP(NUM_GP_VOLUME+1))
      IF (NUM_GP_VOLUME > 0) TEMP(1:NUM_GP_VOLUME) = GP_VOLUMES
      TEMP(NUM_GP_VOLUME+1) = GP_VOLUME_DEFINITION(VID,SID,CID,DIRECT .OR. .NOT. SPECIFIED,PRINCIPAL .OR. .NOT. SPECIFIED)
      CALL MOVE_ALLOC(TEMP,GP_VOLUMES)
      NUM_GP_VOLUME = NUM_GP_VOLUME+1
      RETURN
  900 FATAL_ERR = FATAL_ERR+1
      WRITE(ERR,'(A,A)') ' *ERROR: INVALID OR DUPLICATE GPSTRESS VOLUME: ',TRIM(CARD)
      WRITE(F06,'(A,A)') ' *ERROR: INVALID OR DUPLICATE GPSTRESS VOLUME: ',TRIM(CARD)
      END SUBROUTINE PARSE_GPSTRESS_VOLUME_CARD

! Resolve element sets against sorted model IDs. Storage scales with the model,
! never with the largest external ID. Explicitly reject missing/unsupported IDs.
      SUBROUTINE GPSTRESS_VOLUME_ELEMENTS(V, ELEMENTS)
      INTEGER(LONG), INTENT(IN) :: V
      INTEGER(LONG), ALLOCATABLE, INTENT(OUT) :: ELEMENTS(:)
      LOGICAL :: SELECTED(SIZE(ESORT1))
      CHARACTER(LSETLN) :: TEXT, ERRTOK
      CHARACTER(8) :: TOKEN(3), TOKTYP(3)
      CHARACTER(3) :: THRU, EXCEPT
      INTEGER(LONG) :: NULSET, SLEN, POS, NTOKEN, IERROR, TOKLEN, LO, HI, ID, ROW, I, N
      SELECTED = .FALSE.
      LO=0; HI=-1
      CALL GPSTRESS_GET_SET_STRING(GP_VOLUMES(V)%SETID,NULSET,SLEN,TEXT)
      IF (NULSET == 0) THEN
         CALL GP_VOLUME_ERROR(V,'UNDEFINED OR EMPTY ELEMENT SET')
         RETURN
      ENDIF
      IF (TRIM(TEXT) == 'ALL') THEN
         SELECTED = .TRUE.
      ELSE
         POS = 1
         TOKLEN = SLEN
         THRU = 'OFF'
         EXCEPT = 'OFF'
         DO
            CALL STOKEN('GPSTRESS_VOLUME_ELEMENTS',TEXT,POS,TOKLEN,NTOKEN,IERROR,TOKTYP,TOKEN,ERRTOK,THRU,EXCEPT)
            IF (IERROR /= 0) THEN
               CALL GP_VOLUME_ERROR(V,'INVALID ELEMENT SET EXPRESSION')
               RETURN
            ENDIF
            IF (NTOKEN <= 0) EXIT
            IF (NTOKEN == 1 .AND. TOKTYP(1) == 'EXCEPT  ') CYCLE
            IF (NTOKEN == 1 .AND. TOKTYP(1) == 'INTEGER ') THEN
               READ(TOKEN(1),*) ID
               ROW = GP_FIND_ID(ID,ESORT1)
               IF (ROW == 0) THEN
                  CALL GP_VOLUME_ERROR(V,'MISSING ELEMENT ID',ID)
                  RETURN
               ENDIF
               IF (EXCEPT == 'ON ' .AND. ID >= LO .AND. ID <= HI) THEN
                  SELECTED(ROW) = .FALSE.
               ELSE
                  SELECTED(ROW) = .TRUE.
                  EXCEPT = 'OFF'
               ENDIF
            ELSE IF (NTOKEN == 3 .AND. TOKTYP(2) == 'THRU    ') THEN
               READ(TOKEN(1),*) LO
               READ(TOKEN(3),*) HI
               IF (HI < LO) THEN
                  CALL GP_VOLUME_ERROR(V,'REVERSED ELEMENT RANGE')
                  RETURN
               ENDIF
               N = COUNT(ESORT1 >= LO .AND. ESORT1 <= HI)
               IF (INT(N,8) /= INT(HI,8)-INT(LO,8)+1) THEN
                  CALL GP_VOLUME_ERROR(V,'ELEMENT RANGE CONTAINS MISSING IDS',LO)
                  RETURN
               ENDIF
               WHERE (ESORT1 >= LO .AND. ESORT1 <= HI) SELECTED = .TRUE.
            ELSE
               CALL GP_VOLUME_ERROR(V,'UNSUPPORTED ELEMENT SET EXPRESSION')
               RETURN
            ENDIF
            IF (POS > TOKLEN) EXIT
         ENDDO
      ENDIF
      IF (.NOT. ANY(SELECTED)) THEN
         CALL GP_VOLUME_ERROR(V,'EMPTY VOLUME')
         RETURN
      ENDIF
      DO I=1,SIZE(SELECTED)
         IF (.NOT. SELECTED(I)) CYCLE
         IF (GP_SOLID_NODES(ETYPE(I)) == 0) THEN
            CALL GP_VOLUME_ERROR(V,'UNSUPPORTED ELEMENT FAMILY',ESORT1(I))
            RETURN
         ENDIF
      ENDDO
      ELEMENTS = PACK([(I,I=1,SIZE(SELECTED))],SELECTED)
      END SUBROUTINE GPSTRESS_VOLUME_ELEMENTS

      INTEGER(LONG) FUNCTION GP_SOLID_NODES(NAME) RESULT(N)
      CHARACTER(LEN=*), INTENT(IN) :: NAME
      N = 0
      SELECT CASE (TRIM(NAME))
      CASE ('TETRA4'); N=4
      CASE ('TETRA10'); N=10
      CASE ('PENTA6'); N=6
      CASE ('PENTA15'); N=15
      CASE ('HEXA8'); N=8
      CASE ('HEXA20'); N=20
      CASE ('PYRAM5'); N=5
      CASE ('PYRAM13'); N=13
      END SELECT
      END FUNCTION GP_SOLID_NODES

! VOLUME results are recovered at corner grids, including quadratic solids.
      INTEGER(LONG) FUNCTION GP_SOLID_CORNERS(NAME) RESULT(N)
      CHARACTER(LEN=*), INTENT(IN) :: NAME
      N=0
      SELECT CASE (NAME(1:4))
      CASE ('TETR'); N=4
      CASE ('PENT'); N=6
      CASE ('HEXA'); N=8
      CASE ('PYRA'); N=5
      END SELECT
      END FUNCTION GP_SOLID_CORNERS

      INTEGER(LONG) FUNCTION GP_FIND_ID(ID,IDS) RESULT(ROW)
      INTEGER(LONG), INTENT(IN) :: ID, IDS(:)
      INTEGER(LONG) :: LO, HI, MID
      LO=1; HI=SIZE(IDS); ROW=0
      DO WHILE (LO <= HI)
         MID=LO+(HI-LO)/2
         IF (IDS(MID) == ID) THEN
            ROW=MID
            RETURN
         ELSE IF (IDS(MID) < ID) THEN
            LO=MID+1
         ELSE
            HI=MID-1
         ENDIF
      ENDDO
      END FUNCTION GP_FIND_ID

      SUBROUTINE GP_VOLUME_ERROR(V,MESSAGE,ID)
      INTEGER(LONG), INTENT(IN) :: V
      CHARACTER(LEN=*), INTENT(IN) :: MESSAGE
      INTEGER(LONG), OPTIONAL, INTENT(IN) :: ID
      INTEGER(LONG) :: VALUE
      VALUE=0
      IF (PRESENT(ID)) VALUE=ID
      FATAL_ERR=FATAL_ERR+1
      WRITE(ERR,'(A,I0,2A,I0)') ' *ERROR: GPSTRESS VOLUME ',GP_VOLUMES(V)%ID,': ',TRIM(MESSAGE)//' ID=',VALUE
      WRITE(F06,'(A,I0,2A,I0)') ' *ERROR: GPSTRESS VOLUME ',GP_VOLUMES(V)%ID,': ',TRIM(MESSAGE)//' ID=',VALUE
      END SUBROUTINE GP_VOLUME_ERROR

! Initial recovery mode: interior-point tensor extrapolation; center projection
! for pyramids; equal element weights. MSC numerical parity remains unvalidated.
      SUBROUTINE GP_RECOVER_VOLUME(V,GRIDS,VALUES)
      USE PENTIUM_II_KIND, ONLY : DOUBLE
      USE IEEE_ARITHMETIC, ONLY : IEEE_IS_FINITE, IEEE_VALUE, IEEE_QUIET_NAN
      USE MODEL_STUF, ONLY : GRID_ID, AGRID, ELGP, TYPE, STRESS, CORD, RCORD
      USE PARAMS, ONLY : STR_CID
      USE EMG_Interface
      USE ELMDIS_Interface
      USE ELEM_STRE_STRN_ARRAYS_Interface
      USE MODEL_STUF, ONLY : NUM_EMG_FATAL_ERRS
      USE MODEL_STUF, ONLY : SUBLOD
      USE SCONTR, ONLY : INT_SC_NUM
      INTEGER(LONG), INTENT(IN) :: V
      INTEGER(LONG), ALLOCATABLE, INTENT(OUT) :: GRIDS(:)
      REAL(DOUBLE), ALLOCATABLE, INTENT(OUT) :: VALUES(:,:)
      INTEGER(LONG), ALLOCATABLE :: ELEMENTS(:), COUNTS(:)
      LOGICAL :: USED(SIZE(GRID_ID))
      INTEGER(LONG) :: I,J,K,N,ROW,OLD_CID,BEFORE,ICORD,NG
      REAL(DOUBLE) :: NODAL(6,20)
      REAL(DOUBLE) :: R(3,3),TENSOR(3,3),ROTATED(3,3)
      CHARACTER(1) :: OPT(6)

      CALL GPSTRESS_VOLUME_ELEMENTS(V,ELEMENTS)
      IF (.NOT. ALLOCATED(ELEMENTS)) RETURN
      ICORD=0
      IF (GP_VOLUMES(V)%CID /= 0) THEN
         DO I=1,SIZE(CORD,1)
            IF (CORD(I,2) == GP_VOLUMES(V)%CID) ICORD=I
         ENDDO
         IF (ICORD == 0) THEN
            CALL GP_VOLUME_ERROR(V,'UNDEFINED OUTPUT COORDINATE SYSTEM',GP_VOLUMES(V)%CID)
            RETURN
         ENDIF
! Curvilinear systems require a grid-dependent basis and are not supported yet.
         IF (MOD(CORD(ICORD,1),10) /= 1) THEN
            CALL GP_VOLUME_ERROR(V,'ONLY RECTANGULAR OUTPUT COORDINATES SUPPORTED',GP_VOLUMES(V)%CID)
            RETURN
         ENDIF
      ENDIF
      USED=.FALSE.
      DO I=1,SIZE(ELEMENTS)
         ROW=ELEMENTS(I)
         N=GP_SOLID_CORNERS(ETYPE(ROW))
         DO J=1,N
            K=GP_FIND_ID(EDAT(EPNT(ROW)+1+J),GRID_ID)
            IF (K == 0) THEN
               CALL GP_VOLUME_ERROR(V,'MISSING CONNECTED GRID',EDAT(EPNT(ROW)+1+J))
               RETURN
            ENDIF
            USED(K)=.TRUE.
         ENDDO
      ENDDO
      GRIDS=PACK(GRID_ID,USED)
      NG=SIZE(GRIDS)
      ALLOCATE(VALUES(6,NG),COUNTS(NG))
      VALUES=0D0; COUNTS=0
      OLD_CID=STR_CID
! Recover in BASIC through the shared solid pipeline. Apply the volume basis
! here: RCORD stores the coordinate-to-basic rotation, so its transpose is R.
      STR_CID=0
      R=0D0
      DO I=1,3
         R(I,I)=1D0
      ENDDO
      IF (ICORD /= 0) THEN
         DO I=1,3
            DO J=1,3
               R(J,I)=RCORD(ICORD,3+3*(I-1)+J)
            ENDDO
         ENDDO
      ENDIF
      OPT=['N','N','N','N','N','N']
! With every option N, EMG returns before ELMDAT2 loads temperatures. Ask
! for the thermal load operator in thermal subcases so recovery sees DT.
      IF (SUBLOD(INT_SC_NUM,2)>0) OPT(2)='Y'
      DO I=1,SIZE(ELEMENTS)
         BEFORE=FATAL_ERR
         OPT(3)='N'
         IF (ETYPE(ELEMENTS(I))(1:4) == 'PYRA') OPT(3)='Y'
         CALL EMG(ELEMENTS(I),OPT,'N','GP_RECOVER_VOLUME','N')
         IF (NUM_EMG_FATAL_ERRS > 0 .OR. FATAL_ERR > BEFORE) THEN
            CALL GP_VOLUME_ERROR(V,'ELEMENT RECOVERY/JACOBIAN REJECTED',ESORT1(ELEMENTS(I)))
            EXIT
         ENDIF
         CALL ELMDIS
         N=ELGP
         NODAL=IEEE_VALUE(0D0,IEEE_QUIET_NAN)
         IF (TYPE(1:4) == 'PYRA') THEN
            CALL ELEM_STRE_STRN_ARRAYS(1_LONG)
            DO J=1,N
               NODAL(:,J)=STRESS(1:6)
            ENDDO
         ELSE
            CALL GP_SOLID_INTERIOR_RECOVERY(NODAL(:,1:N))
         ENDIF
         DO J=1,N
            TENSOR=RESHAPE([NODAL(1,J),NODAL(4,J),NODAL(6,J), &
                            NODAL(4,J),NODAL(2,J),NODAL(5,J), &
                            NODAL(6,J),NODAL(5,J),NODAL(3,J)],[3,3])
            ROTATED=MATMUL(R,MATMUL(TENSOR,TRANSPOSE(R)))
            NODAL(:,J)=[ROTATED(1,1),ROTATED(2,2),ROTATED(3,3),ROTATED(1,2),ROTATED(2,3),ROTATED(3,1)]
         ENDDO
         IF (FATAL_ERR > BEFORE .OR. .NOT. ALL(IEEE_IS_FINITE(NODAL(:,1:N)))) THEN
            CALL GP_VOLUME_ERROR(V,'INVALID/NONFINITE STRESS CONTRIBUTION',ESORT1(ELEMENTS(I)))
            EXIT
         ENDIF
         DO J=1,GP_SOLID_CORNERS(TYPE)
            K=GP_FIND_ID(AGRID(J),GRIDS)
            IF (K == 0) THEN
               CALL GP_VOLUME_ERROR(V,'INCONSISTENT RECOVERY CONNECTIVITY',AGRID(J))
               EXIT
            ENDIF
            VALUES(:,K)=VALUES(:,K)+NODAL(:,J)
            COUNTS(K)=COUNTS(K)+1
         ENDDO
      ENDDO
      STR_CID=OLD_CID
      IF (FATAL_ERR > BEFORE) THEN
         DEALLOCATE(GRIDS,VALUES)
         RETURN
      ENDIF
      IF (ANY(COUNTS == 0)) THEN
         CALL GP_VOLUME_ERROR(V,'GRID HAS NO VALID CONTRIBUTIONS')
         DEALLOCATE(GRIDS,VALUES)
         RETURN
      ENDIF
      DO K=1,NG
         VALUES(:,K)=VALUES(:,K)/COUNTS(K)
      ENDDO
      IF (.NOT. ALL(IEEE_IS_FINITE(VALUES))) THEN
         CALL GP_VOLUME_ERROR(V,'NONFINITE ACCUMULATED GRID STRESS')
         DEALLOCATE(GRIDS,VALUES)
      ENDIF
      END SUBROUTINE GP_RECOVER_VOLUME

! Reuse the solver's shape, Jacobian, B and constitutive recovery routines.
! Sample within the element; extrapolate a trilinear HEXA, linear tetrahedral,
! or wedge tensor field. Midside targets use the same interpolated field.
      SUBROUTINE GP_SOLID_INTERIOR_RECOVERY(NODAL)
      USE PENTIUM_II_KIND, ONLY : DOUBLE
      USE MODEL_STUF, ONLY : ELGP, TYPE, BE1, BE2, STRESS, NUM_EMG_FATAL_ERRS
      USE SHP3DH_Interface
      USE SHP3DT_Interface
      USE SHP3DP_Interface
      USE JAC3D_Interface
      USE B3D_ISOPARAMETRIC_Interface
      USE ELEM_STRE_STRN_ARRAYS_Interface
      REAL(DOUBLE), INTENT(OUT) :: NODAL(6,ELGP)
      REAL(DOUBLE) :: TARGET(3,20),POINTS(3,8),SAMPLES(6,8),WEIGHTS(8)
      REAL(DOUBLE) :: PSH(ELGP),DN(3,ELGP),DX(3,ELGP),JAC(3,3),INVJ(3,3),DETJ,BMAT(6,3*ELGP)
      REAL(DOUBLE) :: G,A,B,R,S,T,LI(4),TW(3)
      INTEGER(LONG) :: NC,NP,EDGES(2,12),I,J,K,Q,BEFORE
      TARGET=0D0; POINTS=0D0; EDGES=0
      G=1D0/SQRT(3D0)
      A=(5D0+3D0*SQRT(5D0))/20D0
      B=(5D0-SQRT(5D0))/20D0
      SELECT CASE (TYPE(1:4))
      CASE ('HEXA')
         NC=8; NP=8
         TARGET(:,1:8)=RESHAPE([-1D0,-1D0,-1D0,1D0,-1D0,-1D0,1D0,1D0,-1D0,-1D0,1D0,-1D0, &
                               -1D0,-1D0,1D0,1D0,-1D0,1D0,1D0,1D0,1D0,-1D0,1D0,1D0],[3,8])
         POINTS=G*TARGET(:,1:8)
         EDGES=RESHAPE([1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5],[2,12])
      CASE ('PENT')
         NC=6; NP=6
         TARGET(:,1:6)=RESHAPE([0D0,0D0,-1D0,1D0,0D0,-1D0,0D0,1D0,-1D0, &
                               0D0,0D0,1D0,1D0,0D0,1D0,0D0,1D0,1D0],[3,6])
         POINTS(:,1:6)=RESHAPE([1D0/6,1D0/6,-G,2D0/3,1D0/6,-G,1D0/6,2D0/3,-G, &
                               1D0/6,1D0/6,G,2D0/3,1D0/6,G,1D0/6,2D0/3,G],[3,6])
         EDGES(:,1:9)=RESHAPE([1,2,2,3,3,1,1,4,2,5,3,6,4,5,5,6,6,4],[2,9])
      CASE ('TETR')
         NC=4; NP=4
         TARGET(:,1:4)=RESHAPE([0D0,0D0,0D0,1D0,0D0,0D0,0D0,1D0,0D0,0D0,0D0,1D0],[3,4])
         POINTS(:,1:4)=RESHAPE([B,B,B,A,B,B,B,A,B,B,B,A],[3,4])
         EDGES(:,1:6)=RESHAPE([1,2,2,3,3,1,1,4,2,4,3,4],[2,6])
      CASE DEFAULT
         FATAL_ERR=FATAL_ERR+1
         RETURN
      END SELECT
      DO I=NC+1,ELGP
         K=I-NC
         TARGET(:,I)=0.5D0*(TARGET(:,EDGES(1,K))+TARGET(:,EDGES(2,K)))
      ENDDO
      BEFORE=FATAL_ERR
      DO Q=1,NP
         R=POINTS(1,Q); S=POINTS(2,Q); T=POINTS(3,Q)
         SELECT CASE (TYPE(1:4))
         CASE ('HEXA')
            CALL SHP3DH(1_LONG,1_LONG,1_LONG,ELGP,'GPSTRESS','INTERIOR RECOVERY',2_LONG,R,S,T,'N',PSH,DN)
         CASE ('PENT')
            CALL SHP3DP(1_LONG,1_LONG,1_LONG,ELGP,'GPSTRESS','INTERIOR RECOVERY',2_LONG,2_LONG,R,S,T,'N',PSH,DN)
         CASE ('TETR')
            CALL SHP3DT(1_LONG,ELGP,'GPSTRESS','INTERIOR RECOVERY',2_LONG,R,S,T,'N',PSH,DN)
         END SELECT
         CALL JAC3D(R,S,T,DN,'N',JAC,INVJ,DETJ)
         IF (FATAL_ERR > BEFORE .OR. NUM_EMG_FATAL_ERRS > 0 .OR. DETJ <= 0D0) RETURN
         DX=MATMUL(INVJ,DN)
         CALL B3D_ISOPARAMETRIC(DX,0_LONG,1_LONG,1_LONG,1_LONG,'all strains','N',BMAT)
         BE1(:,:,1)=0D0; BE2(:,:,1)=0D0
         DO I=1,ELGP
            DO J=1,3
               K=6*(I-1)+J
               BE1(:,K,1)=BMAT(1:3,3*(I-1)+J)
               BE2(:,K,1)=BMAT(4:6,3*(I-1)+J)
            ENDDO
         ENDDO
         CALL ELEM_STRE_STRN_ARRAYS(1_LONG,POINTS(:,Q))
         IF (FATAL_ERR > BEFORE) RETURN
         SAMPLES(:,Q)=STRESS(1:6)
      ENDDO
      DO I=1,ELGP
         R=TARGET(1,I); S=TARGET(2,I); T=TARGET(3,I)
         WEIGHTS=0D0
         SELECT CASE (TYPE(1:4))
         CASE ('HEXA')
            DO Q=1,NP
               WEIGHTS(Q)=PRODUCT(1D0+TARGET(:,I)*POINTS(:,Q)/(G*G))/8D0
            ENDDO
         CASE ('PENT')
            TW=[5D0/3D0-2D0*(R+S),2D0*R-1D0/3D0,2D0*S-1D0/3D0]
            WEIGHTS(1:3)=TW*(1D0-T/G)*0.5D0
            WEIGHTS(4:6)=TW*(1D0+T/G)*0.5D0
         CASE ('TETR')
            LI=[1D0-R-S-T,R,S,T]
            WEIGHTS(1:4)=(LI-B)/(A-B)
         END SELECT
         NODAL(:,I)=MATMUL(SAMPLES(:,1:NP),WEIGHTS(1:NP))
      ENDDO
      END SUBROUTINE GP_SOLID_INTERIOR_RECOVERY

      SUBROUTINE WRITE_GPSTRESS_VOLUMES
      USE PENTIUM_II_KIND, ONLY : DOUBLE
      USE CC_OUTPUT_DESCRIBERS, ONLY : GPSTRESS_REQ, GP_VOLUME_REQUESTS
      USE MODEL_STUF, ONLY : SCNUM
      USE IOUNT1, ONLY : OP2
      USE SCONTR, ONLY : INT_SC_NUM, SOL_NAME
      INTEGER(LONG) :: V,I,REQUEST_INDEX,SETID
      CHARACTER(8) :: REQUEST_OUT
      INTEGER(LONG), ALLOCATABLE :: GRIDS(:)
      REAL(DOUBLE), ALLOCATABLE :: VALUES(:,:)
      LOGICAL :: OP2_OPENED
      IF (.NOT. GPSTRESS_REQ .OR. NUM_GP_VOLUME == 0) RETURN
      REQUEST_INDEX=0
      IF (.NOT. ALLOCATED(GP_VOLUME_REQUESTS)) RETURN
      DO I=1,SIZE(GP_VOLUME_REQUESTS)
         IF (GP_VOLUME_REQUESTS(I)%SUBCASE_INDEX == 0) REQUEST_INDEX=I
      ENDDO
      DO I=1,SIZE(GP_VOLUME_REQUESTS)
         IF (GP_VOLUME_REQUESTS(I)%SUBCASE_INDEX == INT_SC_NUM) REQUEST_INDEX=I
      ENDDO
      IF (REQUEST_INDEX == 0) RETURN
      SETID=GP_VOLUME_REQUESTS(REQUEST_INDEX)%SETID
      REQUEST_OUT=GP_VOLUME_REQUESTS(REQUEST_INDEX)%OUTPUT
      IF (SETID == 0) RETURN
      IF (SOL_NAME /= 'STATICS') THEN
         CALL GP_VOLUME_ERROR(1_LONG,'VOLUME RECOVERY CURRENTLY REQUIRES SOL 101')
         RETURN
      ENDIF
      DO V=1,NUM_GP_VOLUME
         IF (.NOT. GP_VOLUME_SELECTED(V,SETID)) CYCLE
         CALL GP_RECOVER_VOLUME(V,GRIDS,VALUES)
         IF (.NOT. ALLOCATED(VALUES)) RETURN
         IF (REQUEST_OUT(1:1) == 'Y') THEN
            WRITE(F06,'(/,A)') ' G R I D   P O I N T   S T R E S S E S'
            WRITE(F06,'(A,I0)') ' SUBCASE = ',SCNUM(INT_SC_NUM)
            WRITE(F06,'(A,I0)') ' VOLUME  = ',GP_VOLUMES(V)%ID
            IF (GP_VOLUMES(V)%CID == 0) THEN
               WRITE(F06,'(A)') ' SYSTEM  = BASIC'
            ELSE
               WRITE(F06,'(A,I0)') ' SYSTEM  = CORD ',GP_VOLUMES(V)%CID
            ENDIF
            WRITE(F06,'(A)') ' RECOVERY = INTERIOR EXTRAPOLATION / PYRAMID CENTER, EQUAL WEIGHTS'
            WRITE(F06,'(A)') ' GRID                 SXX                SYY                SZZ'// &
                           '                SXY                SYZ                SZX'
            DO I=1,SIZE(GRIDS)
               WRITE(F06,'(I10,6(1X,ES18.10))') GRIDS(I),VALUES(:,I)
            ENDDO
         ENDIF
         INQUIRE(UNIT=OP2,OPENED=OP2_OPENED)
         IF (REQUEST_OUT(2:2) == 'Y' .AND. OP2_OPENED) CALL GP_WRITE_OGS_VOLUME(V,GRIDS,VALUES)
         DEALLOCATE(GRIDS,VALUES)
      ENDDO
      END SUBROUTINE WRITE_GPSTRESS_VOLUMES

      LOGICAL FUNCTION GP_VOLUME_SELECTED(V,SETID) RESULT(SELECTED)
      INTEGER(LONG), INTENT(IN) :: V,SETID
      CHARACTER(LSETLN) :: TEXT,ERRTOK
      CHARACTER(8) :: TOKEN(3),TOKTYP(3)
      CHARACTER(3) :: THRU,EXCEPT
      INTEGER(LONG) :: NULSET,SLEN,POS,NTOKEN,IERROR,LO,HI,ID
      SELECTED=SETID == -1
      IF (SETID == -1) RETURN
      CALL GPSTRESS_GET_SET_STRING(SETID,NULSET,SLEN,TEXT)
      IF (NULSET == 0) THEN
         CALL GP_VOLUME_ERROR(V,'UNDEFINED GPSTRESS DOMAIN SELECTION SET',SETID)
         RETURN
      ENDIF
      IF (TRIM(TEXT) == 'ALL') THEN
         SELECTED=.TRUE.
         RETURN
      ENDIF
      POS=1; THRU='OFF'; EXCEPT='OFF'; LO=0; HI=-1
      DO
         CALL STOKEN('GP_VOLUME_SELECTED',TEXT,POS,SLEN,NTOKEN,IERROR,TOKTYP,TOKEN,ERRTOK,THRU,EXCEPT)
         IF (IERROR /= 0) THEN
            CALL GP_VOLUME_ERROR(V,'INVALID GPSTRESS DOMAIN SELECTION SET',SETID)
            RETURN
         ENDIF
         IF (NTOKEN == 0) EXIT
         IF (NTOKEN == 1 .AND. TOKTYP(1) == 'EXCEPT  ') CYCLE
         IF (NTOKEN == 1 .AND. TOKTYP(1) == 'INTEGER ') THEN
            READ(TOKEN(1),*) ID
            IF (EXCEPT == 'ON ' .AND. ID >= LO .AND. ID <= HI) THEN
               IF (GP_VOLUMES(V)%ID == ID) SELECTED=.FALSE.
            ELSE
               EXCEPT='OFF'
               IF (GP_VOLUMES(V)%ID == ID) SELECTED=.TRUE.
            ENDIF
         ELSE IF (NTOKEN == 3 .AND. TOKTYP(2) == 'THRU    ') THEN
            READ(TOKEN(1),*) LO
            READ(TOKEN(3),*) HI
            IF (GP_VOLUMES(V)%ID >= LO .AND. GP_VOLUMES(V)%ID <= HI) SELECTED=.TRUE.
         ELSE
            CALL GP_VOLUME_ERROR(V,'UNSUPPORTED GPSTRESS DOMAIN SELECTION SET',SETID)
            RETURN
         ENDIF
         IF (POS > SLEN) EXIT
      ENDDO
      END FUNCTION GP_VOLUME_SELECTED

! OGS1 real SORT1 volume direct: table code 27, width 9 (i8f).
! Header word 3 is volume ID; word 8 is reference CID.
! Confirmed against pyNastran OGS._read_ogs1_table27_numwide9.
      SUBROUTINE GP_WRITE_OGS_VOLUME(V,GRIDS,VALUES)
      USE PENTIUM_II_KIND, ONLY : DOUBLE
      USE IOUNT1, ONLY : OP2
      USE LAPACK_STD_EIG_1
      USE MODEL_STUF, ONLY : TITLE, STITLE, LABEL, SCNUM
      USE SCONTR, ONLY : INT_SC_NUM
      INTEGER(LONG), INTENT(IN) :: V,GRIDS(:)
      REAL(DOUBLE), INTENT(IN) :: VALUES(:,:)
      INTEGER(LONG) :: HEADER(50),ITABLE,I,J,K,INFO,ORDER(3)
      REAL(DOUBLE) :: AXES(3,3),EIG(3),WORK(32),PRINC(14,SIZE(GRIDS)),SCALE
      CHARACTER(128) :: T,S,L
      REAL(DOUBLE) :: PRESSURE(SIZE(GRIDS)),VM(SIZE(GRIDS))
      HEADER=0
      HEADER(1)=13
      HEADER(2)=27
      HEADER(3)=GP_VOLUMES(V)%ID
      HEADER(4)=SCNUM(INT_SC_NUM)
      HEADER(5)=SCNUM(INT_SC_NUM)
      HEADER(8)=GP_VOLUMES(V)%CID
      HEADER(9)=1
      HEADER(10)=9
      HEADER(11)=1
      HEADER(12)=3
      T=TITLE(INT_SC_NUM); S=STITLE(INT_SC_NUM); L=LABEL(INT_SC_NUM)
      DO I=1,SIZE(GRIDS)
         PRESSURE(I)=-SUM(VALUES(1:3,I))/3D0
         VM(I)=SQRT(0.5D0*((VALUES(1,I)-VALUES(2,I))**2+(VALUES(2,I)-VALUES(3,I))**2+ &
                       (VALUES(3,I)-VALUES(1,I))**2)+3D0*SUM(VALUES(4:6,I)**2))
      ENDDO
      IF (GP_VOLUMES(V)%DIRECT) THEN
         CALL WRITE_TABLE_HEADER('OGS1    ')
         ITABLE=-3
         CALL WRITE_ITABLE(ITABLE)
         WRITE(OP2) 146_LONG
         WRITE(OP2) HEADER,T,S,L
         ITABLE=ITABLE-1
         CALL WRITE_ITABLE(ITABLE)
         ITABLE=ITABLE-1
         WRITE(OP2) 9*SIZE(GRIDS)
         WRITE(OP2) (GRIDS(I)*10+3,(REAL(VALUES(J,I),4),J=1,6),REAL(PRESSURE(I),4),REAL(VM(I),4),I=1,SIZE(GRIDS))
         CALL END_OP2_TABLE(ITABLE)
      ENDIF
      IF (.NOT. GP_VOLUMES(V)%PRINCIPAL) RETURN
! Observed 5-001 table 28 layout: three principal values, then the
! direction-cosine matrix by rows, pressure and von Mises. Axis order A,B,C
! corresponds to descending maximum, minimum, intermediate eigenvalues.
      ORDER=[3,1,2]
      DO I=1,SIZE(GRIDS)
         AXES=RESHAPE([VALUES(1,I),VALUES(4,I),VALUES(6,I), &
                       VALUES(4,I),VALUES(2,I),VALUES(5,I), &
                       VALUES(6,I),VALUES(5,I),VALUES(3,I)],[3,3])
         CALL DSYEV('V','U',3_LONG,AXES,3_LONG,EIG,WORK,32_LONG,INFO)
         IF (INFO /= 0) THEN
            CALL GP_VOLUME_ERROR(V,'PRINCIPAL STRESS EIGENSOLUTION FAILED',GRIDS(I))
            RETURN
         ENDIF
         SCALE=MAX(1D0,MAXVAL(ABS(EIG)))
         DO J=1,3
            K=MAXLOC(ABS(AXES(:,J)),DIM=1)
            IF (AXES(K,J)<0D0) AXES(:,J)=-AXES(:,J)
         ENDDO
! Repeated principal stresses have no unique axis. The reference marks
! their direction cosines zero; use single-precision output resolution.
         DO J=1,2
            IF (ABS(EIG(J+1)-EIG(J))<=1D-6*SCALE) AXES(:,J:J+1)=0D0
         ENDDO
         PRINC(1:3,I)=EIG(ORDER)
         DO J=1,3
            PRINC(3+3*(J-1)+1:3+3*J,I)=AXES(J,ORDER)
         ENDDO
         PRINC(13,I)=PRESSURE(I); PRINC(14,I)=VM(I)
      ENDDO
      HEADER(2)=28; HEADER(10)=15
      CALL WRITE_TABLE_HEADER('OGS1    ')
      ITABLE=-3
      CALL WRITE_ITABLE(ITABLE)
      WRITE(OP2) 146_LONG
      WRITE(OP2) HEADER,T,S,L
      ITABLE=ITABLE-1
      CALL WRITE_ITABLE(ITABLE)
      ITABLE=ITABLE-1
      WRITE(OP2) 15*SIZE(GRIDS)
      WRITE(OP2) (GRIDS(I)*10+3,(REAL(PRINC(J,I),4),J=1,14),I=1,SIZE(GRIDS))
      CALL END_OP2_TABLE(ITABLE)
      END SUBROUTINE GP_WRITE_OGS_VOLUME

      SUBROUTINE PARSE_GPSTRESS_SET_CARD ( CARD )

      CHARACTER(LEN=*), INTENT(IN)    :: CARD

      CHARACTER(LEN=CC_ENTRY_LEN)     :: CARD1
      CHARACTER(LEN=LSETLN)           :: SET_TEXT
      CHARACTER(16*BYTE)              :: KEY1
      INTEGER(LONG)                   :: IERR
      INTEGER(LONG)                   :: IOS
      INTEGER(LONG)                   :: ICMT
      INTEGER(LONG)                   :: K
      INTEGER(LONG)                   :: SETID
      INTEGER(LONG)                   :: TEXT_BEG
      INTEGER(LONG)                   :: TEXT_END
      LOGICAL                         :: SEEN_SETID

      CARD1 = CARD
      CALL GPSTRESS_TO_UPPER_LINE ( CARD1 )
      ICMT = INDEX(CARD1,'$')
      IF (ICMT > 0) CARD1(ICMT:) = ' '

      KEY1 = ' '
      SETID = 0
      READ(CARD1,*,IOSTAT=IOS) KEY1, SETID
      IF ((IOS /= 0) .OR. (TRIM(KEY1) /= 'SET') .OR. (SETID <= 0)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR    : CANNOT PARSE GPSTRESS OUTPUT(POST) SET ENTRY: ', TRIM(CARD)
         WRITE(F06,'(A,A)') ' *ERROR    : CANNOT PARSE GPSTRESS OUTPUT(POST) SET ENTRY: ', TRIM(CARD)
         RETURN
      ENDIF

      SEEN_SETID = .FALSE.
      TEXT_BEG = 0
      TEXT_END = 0
      K = 0
      DO K=5,CC_ENTRY_LEN
         IF (CARD1(K:K) /= ' ') THEN
            SEEN_SETID = .TRUE.
         ELSE IF (SEEN_SETID) THEN
            TEXT_BEG = K + 1
            EXIT
         ENDIF
      ENDDO
      IF (TEXT_BEG == 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR    : GPSTRESS OUTPUT(POST) SET ENTRY HAS NO SET DATA: ', TRIM(CARD)
         WRITE(F06,'(A,A)') ' *ERROR    : GPSTRESS OUTPUT(POST) SET ENTRY HAS NO SET DATA: ', TRIM(CARD)
         RETURN
      ENDIF
      IF (CARD1(TEXT_BEG:TEXT_BEG) == '=') TEXT_BEG = TEXT_BEG + 1
      DO WHILE ((TEXT_BEG <= CC_ENTRY_LEN) .AND. (CARD1(TEXT_BEG:TEXT_BEG) == ' '))
         TEXT_BEG = TEXT_BEG + 1
      ENDDO
      TEXT_END = CC_ENTRY_LEN
      DO K=CC_ENTRY_LEN,TEXT_BEG,-1
         IF (CARD1(K:K) /= ' ') THEN
            TEXT_END = K
            EXIT
         ENDIF
      ENDDO

      SET_TEXT = ' '
      IF (TEXT_END >= TEXT_BEG) THEN
         SET_TEXT(1:MIN(LSETLN,TEXT_END-TEXT_BEG+1)) = CARD1(TEXT_BEG:MIN(TEXT_END,TEXT_BEG+LSETLN-1))
      ENDIF
      CALL REGISTER_GP_POST_SET ( SETID, SET_TEXT, IERR )
      IF (IERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,I8)') ' *ERROR    : TOO MANY GPSTRESS OUTPUT(POST) SET DEFINITIONS. LAST SET ID = ', SETID
         WRITE(F06,'(A,I8)') ' *ERROR    : TOO MANY GPSTRESS OUTPUT(POST) SET DEFINITIONS. LAST SET ID = ', SETID
      ENDIF

      END SUBROUTINE PARSE_GPSTRESS_SET_CARD

      SUBROUTINE PARSE_GPSTRESS_SURFACE_CARD ( CARD )

      CHARACTER(LEN=*), INTENT(IN)    :: CARD

      CHARACTER(LEN=CC_ENTRY_LEN)     :: CARD1
      CHARACTER(16*BYTE)              :: KEY1, KEY2, KEY3, KEY4
      INTEGER(LONG)                   :: IERR
      INTEGER(LONG)                   :: IOS
      INTEGER(LONG)                   :: ICMT
      INTEGER(LONG)                   :: SETID
      INTEGER(LONG)                   :: SURFACE_ID

      CARD1 = CARD
      CALL GPSTRESS_TO_UPPER_LINE ( CARD1 )
      ICMT = INDEX(CARD1,'$')
      IF (ICMT > 0) CARD1(ICMT:) = ' '

      KEY1 = ' '
      KEY2 = ' '
      KEY3 = ' '
      KEY4 = ' '
      SURFACE_ID = 0
      SETID = 0

      READ(CARD1,*,IOSTAT=IOS) KEY1, SURFACE_ID, KEY2, SETID, KEY3, KEY4
      IF (IOS /= 0) THEN
         KEY1 = ' '
         KEY2 = ' '
         KEY3 = ' '
         KEY4 = ' '
         SURFACE_ID = 0
         SETID = 0
         READ(CARD1,*,IOSTAT=IOS) KEY1, SURFACE_ID, KEY2, SETID
      ENDIF

      IF (IOS /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR    : CANNOT PARSE GPSTRESS SURFACE ENTRY: ', TRIM(CARD)
         WRITE(F06,'(A,A)') ' *ERROR    : CANNOT PARSE GPSTRESS SURFACE ENTRY: ', TRIM(CARD)
         RETURN
      ENDIF

      IF ((TRIM(KEY1) /= 'SURFACE') .OR. (TRIM(KEY2) /= 'SET')) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR    : GPSTRESS SURFACE ENTRY MUST LOOK LIKE "SURFACE id SET setid [NORMAL dir]": ', TRIM(CARD)
         WRITE(F06,'(A,A)') ' *ERROR    : GPSTRESS SURFACE ENTRY MUST LOOK LIKE "SURFACE id SET setid [NORMAL dir]": ', TRIM(CARD)
         RETURN
      ENDIF

      IF ((SURFACE_ID <= 0) .OR. (SETID <= 0)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,A)') ' *ERROR    : GPSTRESS SURFACE ENTRY HAS NONPOSITIVE SURFACE/SET ID: ', TRIM(CARD)
         WRITE(F06,'(A,A)') ' *ERROR    : GPSTRESS SURFACE ENTRY HAS NONPOSITIVE SURFACE/SET ID: ', TRIM(CARD)
         RETURN
      ENDIF

      IF (TRIM(KEY3) == 'NORMAL') THEN
         KEY3 = KEY4
      ENDIF
      IF (LEN_TRIM(KEY3) == 0) KEY3 = 'BASIC'

      CALL REGISTER_GP_SURFACE ( SURFACE_ID, SETID, KEY3, IERR )
      IF (IERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,'(A,I8)') ' *ERROR    : TOO MANY GPSTRESS SURFACE DEFINITIONS. LAST SURFACE ID = ', SURFACE_ID
         WRITE(F06,'(A,I8)') ' *ERROR    : TOO MANY GPSTRESS SURFACE DEFINITIONS. LAST SURFACE ID = ', SURFACE_ID
      ENDIF

      END SUBROUTINE PARSE_GPSTRESS_SURFACE_CARD

      SUBROUTINE GPSTRESS_GET_SET_STRING ( SETID_IN, NULSET, SLEN, TOKSTR_OUT )

      INTEGER(LONG), INTENT(IN)       :: SETID_IN
      INTEGER(LONG), INTENT(OUT)      :: NULSET
      INTEGER(LONG), INTENT(OUT)      :: SLEN
      CHARACTER(LEN=*), INTENT(OUT)   :: TOKSTR_OUT

      CHARACTER(LSETLN*BYTE)          :: SETCHR
      CHARACTER(8*BYTE)               :: TOKENI
      CHARACTER(LSETLN*BYTE)          :: TOKSTR
      CHARACTER(8*BYTE)               :: TOKTY1
      INTEGER(LONG)                   :: DATA_BEG, DATA_END, ECOL, I, IOCHK, K, POSN, SETID
      INTEGER(LONG)                   :: SET_1_BEG, SET_2_BEG, SID_BEG, SID_END

      TOKSTR_OUT = ' '
      NULSET = 0
      SLEN   = 0

      DO I=1,NUM_GP_POST_SET
         IF (GP_POST_SET_IDS(I) == SETID_IN) THEN
            TOKSTR_OUT = ' '
            SLEN = LEN_TRIM(GP_POST_SET_TEXT(I))
            IF (SLEN > 0) THEN
               TOKSTR_OUT(1:MIN(SLEN,LEN(TOKSTR_OUT))) = GP_POST_SET_TEXT(I)(1:MIN(SLEN,LEN(TOKSTR_OUT)))
               NULSET = 1
            ENDIF
            RETURN
         ENDIF
      ENDDO

      IF (.NOT. ALLOCATED(ALL_SETS_ARRAY)) RETURN
      IF (SETLEN <= 0) RETURN
      DO I=1,LSETLN
         SETCHR(I:I) = ALL_SETS_ARRAY(I)
      ENDDO

      SET_1_BEG = INDEX(SETCHR(1:),'SET')
      IF (SET_1_BEG == 0) THEN
         WRITE(ERR,'(A,I8,A)') ' *ERROR  1405: SET ID ', SETID_IN, ' NOT FOUND'
         WRITE(F06,'(A,I8,A)') ' *ERROR  1405: SET ID ', SETID_IN, ' NOT FOUND'
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

outer:DO
         TOKSTR(1:) = ' '
         POSN = INDEX(SETCHR(SET_1_BEG+1:),'SET')
         IF (POSN == 0) THEN
            SET_2_BEG = SETLEN + 1
         ELSE
            SET_2_BEG = SET_1_BEG + POSN
         ENDIF

         POSN = INDEX(SETCHR(SET_1_BEG:),'=')
         ECOL = SET_1_BEG + POSN - 1
         DATA_BEG = ECOL + 1
         DATA_END = SET_2_BEG - 1

         DO I=DATA_BEG,DATA_END
            IF ((SETCHR(I:I) == ' ') .OR. (SETCHR(I:I) == ',')) THEN
               DATA_BEG = DATA_BEG + 1
            ELSE
               EXIT
            ENDIF
         ENDDO

         TOKENI(1:) = ' '
         SID_BEG = SET_1_BEG + 3
         SID_END = ECOL - 1
         SETID = 0
         K = 0
         DO I=SID_BEG,SID_END
            IF ((SETCHR(I:I) == ' ') .OR. (SETCHR(I:I) == ',')) CYCLE
            K = K + 1
            IF (K > MAX_TOKEN_LEN) THEN
               NULSET = 0
               SET_1_BEG = SET_2_BEG
               CYCLE outer
            ENDIF
            TOKENI(K:K) = SETCHR(I:I)
         ENDDO

         CALL TOKCHK ( TOKENI, TOKTY1 )
         IF (TOKTY1 == 'INTEGER ') THEN
            READ(TOKENI,'(I8)',IOSTAT=IOCHK) SETID
            IF (IOCHK /= 0) THEN
               NULSET = 0
               SET_1_BEG = SET_2_BEG
               CYCLE outer
            ENDIF
         ELSE
            NULSET = 0
            SET_1_BEG = SET_2_BEG
            CYCLE outer
         ENDIF

         IF (SETID == SETID_IN) THEN
            NULSET = 1
            SLEN = DATA_END - DATA_BEG + 1
            TOKSTR(1:SLEN) = SETCHR(DATA_BEG:DATA_END)
            EXIT
         ELSE
            NULSET = 0
            SLEN = 0
            IF (SET_2_BEG > SETLEN) THEN
               WRITE(ERR,'(A,I8,A)') ' *ERROR  1405: SET ID ', SETID_IN, ' NOT FOUND'
               WRITE(F06,'(A,I8,A)') ' *ERROR  1405: SET ID ', SETID_IN, ' NOT FOUND'
               FATAL_ERR = FATAL_ERR + 1
               RETURN
            ELSE
               SET_1_BEG = SET_2_BEG
               CYCLE outer
            ENDIF
         ENDIF
      ENDDO outer

      IF (NULSET /= 0) THEN
         TOKSTR_OUT(1:MIN(SLEN,LEN(TOKSTR_OUT))) = TOKSTR(1:MIN(SLEN,LEN(TOKSTR_OUT)))
      ENDIF

      END SUBROUTINE GPSTRESS_GET_SET_STRING

      SUBROUTINE GPSTRESS_COLLECT_SURFACE_PATCH ( SURF_INDEX, ELEM_IDS, NUM_ELEMS, GRID_IDS, NUM_GRIDS, IERR )

      INTEGER(LONG), INTENT(IN)       :: SURF_INDEX
      INTEGER(LONG), INTENT(OUT)      :: ELEM_IDS(:)
      INTEGER(LONG), INTENT(OUT)      :: NUM_ELEMS
      INTEGER(LONG), INTENT(OUT)      :: GRID_IDS(:)
      INTEGER(LONG), INTENT(OUT)      :: NUM_GRIDS
      INTEGER(LONG), INTENT(OUT)      :: IERR

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GPSTRESS_COLLECT_SURFACE_PATCH'
      CHARACTER(LSETLN*BYTE)          :: ERRTOK
      CHARACTER(3*BYTE)               :: EXCEPT
      CHARACTER(3*BYTE)               :: THRU
      CHARACTER(8*BYTE)               :: TOKEN(3)
      CHARACTER(8*BYTE)               :: TOKTYP(3)
      CHARACTER(LSETLN*BYTE)          :: TOKSTR
      INTEGER(LONG)                   :: ACT_ELEM, AELEM_HI, AELEM_LO, ESORT1_ROW_NUM, I, IERROR, ISTART, K
      INTEGER(LONG)                   :: NTOKEN, NULSET, SLEN, TOKLEN
      INTEGER(LONG)                   :: SHELL_GRIDS(8), NGR
      LOGICAL                         :: KEEP

      ELEM_IDS = 0
      GRID_IDS = 0
      NUM_ELEMS = 0
      NUM_GRIDS = 0
      IERR = 0

      IF ((SURF_INDEX < 1) .OR. (SURF_INDEX > NUM_GP_SURFACE)) THEN
         IERR = 1
         RETURN
      ENDIF

      CALL GPSTRESS_GET_SET_STRING ( GP_SURFACE_SETIDS(SURF_INDEX), NULSET, SLEN, TOKSTR )
      IF (NULSET == 0) THEN
         IERR = 2
         RETURN
      ENDIF

      IF (TRIM(TOKSTR(1:SLEN)) == 'ALL') THEN
         DO I=1,SIZE(ESORT1)
            IF ((ETYPE(I)(1:5) == 'TRIA3') .OR. (ETYPE(I)(1:5) == 'TRIA6') .OR.                              &
                (ETYPE(I)(1:5) == 'QUAD4') .OR. (ETYPE(I) == 'QUADR   ') .OR. (ETYPE(I)(1:5) == 'QUAD8')) THEN
               CALL ADD_ELEM_IF_SHELL ( ESORT1(I) )
            ENDIF
         ENDDO
         GOTO 900
      ENDIF

      ISTART = 1
      THRU   = 'OFF'
      EXCEPT = 'OFF'
      TOKLEN = SLEN

token_loop: DO
         CALL STOKEN ( SUBR_NAME, TOKSTR, ISTART, TOKLEN, NTOKEN, IERROR, TOKTYP, TOKEN, ERRTOK, THRU, EXCEPT )
         IF (IERROR /= 0) THEN
            IERR = 3
            RETURN
         ENDIF
         IF (NTOKEN <= 0) EXIT token_loop

         IF (NTOKEN == 1) THEN
            IF (TOKTYP(1) /= 'INTEGER ') CYCLE
            READ(TOKEN(1),'(I8)') ACT_ELEM
            IF (EXCEPT == 'OFF') CALL ADD_ELEM_IF_SHELL ( ACT_ELEM )
            IF (EXCEPT == 'ON ') THEN
               IF ((ACT_ELEM >= AELEM_LO) .AND. (ACT_ELEM <= AELEM_HI)) THEN
                  CALL REMOVE_ELEM ( ACT_ELEM )
               ELSE
                  CALL ADD_ELEM_IF_SHELL ( ACT_ELEM )
                  EXCEPT = 'OFF'
                  THRU   = 'OFF'
               ENDIF
            ENDIF

         ELSE IF (NTOKEN == 3) THEN
            IF ((TOKTYP(1) == 'INTEGER ') .AND. (TOKTYP(2) == 'THRU    ') .AND. (TOKTYP(3) == 'INTEGER ')) THEN
               READ(TOKEN(1),'(I8)') AELEM_LO
               READ(TOKEN(3),'(I8)') AELEM_HI
               IF (AELEM_LO <= AELEM_HI) THEN
                  DO ACT_ELEM=AELEM_LO,AELEM_HI
                     CALL ADD_ELEM_IF_SHELL ( ACT_ELEM )
                  ENDDO
               ENDIF
            ENDIF
         ENDIF

         IF (ISTART > TOKLEN) EXIT token_loop
      ENDDO token_loop

      DO I=1,NUM_ELEMS
         CALL GET_SHELL_CORNER_GRIDS ( ELEM_IDS(I), SHELL_GRIDS, NGR, KEEP )
         IF (.NOT. KEEP) CYCLE
         DO K=1,NGR
            CALL ADD_UNIQUE_INT ( SHELL_GRIDS(K), GRID_IDS, NUM_GRIDS )
         ENDDO
      ENDDO

  900 CONTINUE
      DO I=1,NUM_ELEMS
         CALL GET_SHELL_CORNER_GRIDS ( ELEM_IDS(I), SHELL_GRIDS, NGR, KEEP )
         IF (.NOT. KEEP) CYCLE
         DO K=1,NGR
            CALL ADD_UNIQUE_INT ( SHELL_GRIDS(K), GRID_IDS, NUM_GRIDS )
         ENDDO
      ENDDO

      CONTAINS

      SUBROUTINE ADD_ELEM_IF_SHELL ( ACTUAL_ELEM_ID )
      INTEGER(LONG), INTENT(IN)       :: ACTUAL_ELEM_ID
      INTEGER(LONG)                   :: J

      CALL GET_ARRAY_ROW_NUM ( 'ESORT1', SUBR_NAME, SIZE(ESORT1), ESORT1, ACTUAL_ELEM_ID, ESORT1_ROW_NUM )
      IF (ESORT1_ROW_NUM == -1) RETURN

      IF ((ETYPE(ESORT1_ROW_NUM)(1:5) == 'TRIA3') .OR. (ETYPE(ESORT1_ROW_NUM)(1:5) == 'TRIA6') .OR.                     &
          (ETYPE(ESORT1_ROW_NUM)(1:5) == 'QUAD4') .OR. (ETYPE(ESORT1_ROW_NUM)(1:5) == 'QUAD8') .OR.                     &
          (ETYPE(ESORT1_ROW_NUM)      == 'QUADR   ')) THEN
         DO J=1,NUM_ELEMS
            IF (ELEM_IDS(J) == ACTUAL_ELEM_ID) RETURN
         ENDDO
         IF (NUM_ELEMS < SIZE(ELEM_IDS)) THEN
            NUM_ELEMS = NUM_ELEMS + 1
            ELEM_IDS(NUM_ELEMS) = ACTUAL_ELEM_ID
         ELSE
            IERR = 4
         ENDIF
      ENDIF
      END SUBROUTINE ADD_ELEM_IF_SHELL

      SUBROUTINE REMOVE_ELEM ( ACTUAL_ELEM_ID )
      INTEGER(LONG), INTENT(IN)       :: ACTUAL_ELEM_ID
      INTEGER(LONG)                   :: J
      DO J=1,NUM_ELEMS
         IF (ELEM_IDS(J) == ACTUAL_ELEM_ID) THEN
            IF (J < NUM_ELEMS) THEN
               ELEM_IDS(J:NUM_ELEMS-1) = ELEM_IDS(J+1:NUM_ELEMS)
            ENDIF
            ELEM_IDS(NUM_ELEMS) = 0
            NUM_ELEMS = NUM_ELEMS - 1
            EXIT
         ENDIF
      ENDDO
      END SUBROUTINE REMOVE_ELEM

      END SUBROUTINE GPSTRESS_COLLECT_SURFACE_PATCH

      SUBROUTINE GET_SHELL_CORNER_GRIDS ( ACTUAL_ELEM_ID, SHELL_GRIDS, NGR, KEEP )

      INTEGER(LONG), INTENT(IN)       :: ACTUAL_ELEM_ID
      INTEGER(LONG), INTENT(OUT)      :: SHELL_GRIDS(8)
      INTEGER(LONG), INTENT(OUT)      :: NGR
      LOGICAL, INTENT(OUT)            :: KEEP

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_SHELL_CORNER_GRIDS'
      INTEGER(LONG)                   :: EP, I, ROW

      SHELL_GRIDS = 0
      NGR = 0
      KEEP = .FALSE.

      CALL GET_ARRAY_ROW_NUM ( 'ESORT1', SUBR_NAME, SIZE(ESORT1), ESORT1, ACTUAL_ELEM_ID, ROW )
      IF (ROW == -1) RETURN

      EP = EPNT(ROW)
      IF (ETYPE(ROW)(1:5) == 'TRIA3') THEN
         NGR = 3
         DO I=1,NGR
            SHELL_GRIDS(I) = EDAT(EP + 1 + I)
         ENDDO
         KEEP = .TRUE.
      ELSE IF ((ETYPE(ROW)(1:5) == 'QUAD4') .OR. (ETYPE(ROW) == 'QUADR   ')) THEN
         NGR = 4
         DO I=1,NGR
            SHELL_GRIDS(I) = EDAT(EP + 1 + I)
         ENDDO
         KEEP = .TRUE.
      ELSE IF (ETYPE(ROW)(1:5) == 'TRIA6') THEN
         NGR = 6
         DO I=1,NGR
            SHELL_GRIDS(I) = EDAT(EP + 1 + I)
         ENDDO
         KEEP = .TRUE.
      ELSE IF (ETYPE(ROW)(1:5) == 'QUAD8') THEN
         NGR = 8
         DO I=1,NGR
            SHELL_GRIDS(I) = EDAT(EP + 1 + I)
         ENDDO
         KEEP = .TRUE.
      ENDIF

      END SUBROUTINE GET_SHELL_CORNER_GRIDS

      SUBROUTINE ADD_UNIQUE_INT ( VALUE, ARRAY, NUSED )

      INTEGER(LONG), INTENT(IN)       :: VALUE
      INTEGER(LONG), INTENT(INOUT)    :: ARRAY(:)
      INTEGER(LONG), INTENT(INOUT)    :: NUSED

      INTEGER(LONG)                   :: I

      IF (VALUE <= 0) RETURN
      DO I=1,NUSED
         IF (ARRAY(I) == VALUE) RETURN
      ENDDO
      IF (NUSED < SIZE(ARRAY)) THEN
         NUSED = NUSED + 1
         ARRAY(NUSED) = VALUE
      ENDIF

      END SUBROUTINE ADD_UNIQUE_INT

      SUBROUTINE GPSTRESS_TO_UPPER_LINE ( LINE )

      CHARACTER(LEN=*), INTENT(INOUT) :: LINE
      INTEGER(LONG)                   :: I, ICODE

      DO I=1,LEN(LINE)
         ICODE = ICHAR(LINE(I:I))
         IF ((ICODE >= ICHAR('a')) .AND. (ICODE <= ICHAR('z'))) THEN
            LINE(I:I) = CHAR(ICODE - 32)
         ENDIF
      ENDDO

      END SUBROUTINE GPSTRESS_TO_UPPER_LINE

      END MODULE GPSTRESS_SURFACE_UTILS
