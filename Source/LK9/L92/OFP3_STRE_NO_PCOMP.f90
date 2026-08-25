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

      SUBROUTINE OFP3_STRE_NO_PCOMP ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )

! Processes element stress output requests for non PCOMP elements for one subcase. Also write Output Transformation Matrices (OTM's)
! for stresses for Craig-Bampton models)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, ERR, F06, NEU
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELOUT_STRE_BIT, FATAL_ERR, IBIT, INT_SC_NUM,                                &
                                         MAX_STRESS_POINTS, MBUG, MOGEL,                                                           &
                                         NELE, NGRID, NCBAR, NCBEAM, NCBUSH, NCELAS1, NCELAS2, NCELAS3, NCELAS4, NCHEXA8, NCHEXA20,&
                                         NCPENTA6,                                                                                   &
                                         NCPENTA15, NPYRAM5, NPYRAM14, NCTETRA4, NCTETRA10, NCQUAD4, NCQUAD4K, NCQUADR, NCROD,    &
                                         NCSHEAR,                                                                                     &
                                         NCTRIA3, NCTRIA3K, NCTRIA6,                                                                &
                                         SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, HALF, ONE, THREE, FOUR
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_NUMS, FEMAP_EL_VECS
      USE PARAMS, ONLY                :  OTMSKIP, QUAD4TYP, QUADRTYP, TRIA3TYP, TRIARTYP
      USE LINK9_STUFF, ONLY           :  WRITE_NEU_STRE
      USE MODEL_STUF, ONLY            :  AGRID, ANY_STRE_OUTPUT, CBEAM_ACTIVE_NSTATIONS, CBEAM_ACTIVE_XL, EDAT, EPNT, ETYPE, EID, &
                                         ELGP, ELMTYP, ELOUT, METYPE, NUM_SEi, NUM_EMG_FATAL_ERRS, OGROUT, PCOMP_PROPS, PLY_NUM,   &
                                         STRESS, PBEAM_NSTATIONS, TE, TYPE, SHELL_STR_ANGLE, ZS, GRID_ID, XEB
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRE_LOC, STRE_OPT, GPSTRESS_REQ
      USE LINK9_STUFF, ONLY           :  CBEAM_XL_OUT, EID_OUT_ARRAY, GID_OUT_ARRAY, MAXREQ, OGEL, SHELL_OUT_TE,                 &
                                         SHELL_STRESS_IN_LOCAL, POLY_FIT_ERR, POLY_FIT_ERR_INDEX
      USE OUTPUT4_MATRICES, ONLY      :  OTM_STRE, TXT_STRE

      USE PLANE_COORD_TRANS_21_Interface
      USE TRANSFORM_SHELL_STR_Interface
      USE OFP3_STRE_NO_PCOMP_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'OFP3_STRE_NO_PCOMP'
      CHARACTER( 1*BYTE), PARAMETER   :: IHDR      = 'Y'   ! An input to subr WRITE_GRID_OUTPUTS, called herein
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option indicators for subr EMG, called herein
      CHARACTER(31*BYTE)              :: OT4_DESCRIPTOR    ! Descriptor for rows of OT4 file
      CHARACTER(30*BYTE)              :: REQUEST           ! Text for error message
      CHARACTER(20*BYTE)              :: STRESS_ITEM(20)   ! Char description of element stresses

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID      ! Set ID for FEMAP output
      INTEGER(LONG), INTENT(IN)       :: ITE               ! Unit number for text files for OTM row descriptors
      INTEGER(LONG), INTENT(IN)       :: JVEC              ! Solution vector number
      INTEGER(LONG), INTENT(INOUT)    :: OT4_EROW          ! Row number in OT4 file for elem related OTM descriptors
      INTEGER(LONG)                   :: ELOUT_STRE        ! If > 0, there are STRESS   requests for some elems
      INTEGER(LONG)                   :: I,J,K,L,M         ! DO loop indices
      INTEGER(LONG)                   :: IERROR    = 0     ! Local error count
!xx   INTEGER(LONG)                   :: IROW_MAT          ! Row number in OTM's
!xx   INTEGER(LONG)                   :: IROW_TXT          ! Row number in OTM text file
      INTEGER(LONG)                   :: NDUM              ! Value initialized to zero and used in call to CALC_ELEM_STRESSES
      INTEGER(LONG)                   :: NELREQ(METYPE)    ! Count of the no. of requests for ELFORCE(NODE or ENGR) or STRESS
      INTEGER(LONG)                   :: NUM_OGEL_ROWS     ! No. elems processed prior to writing results to F06 file
      INTEGER(LONG)                   :: NUM_FROWS         ! No. elems processed for FEMAP
      INTEGER(LONG)                   :: NUM_OGEL          ! No. rows written to array OGEL prior to writing results to F06 file
!                                                            (this can be > NUM_OGEL_ROWS since more than 1 row is written to OGEL
!                                                            for ELFORCE(NODE) - elem nodal forces)

      INTEGER(LONG)                   :: NUM_OTM_ENTRIES   ! Number of entries in OGEL for a particular element type
      INTEGER(LONG)                   :: NUM_PTS(METYPE)   ! Num diff stress points for one element (3rd dim in arrays SEi, STEi)
      INTEGER(LONG)                   :: NUM_PTS_CUR       ! Actual number of stress points for the current element
      INTEGER(LONG)                   :: NUM_PTS_ELEM      ! Actual number of stress points for the current element in request counting
      INTEGER(LONG)                   :: RECOVERY_POINT    ! Actual SEi/STEi recovery point used for this output point

                                                           ! Stress index (1 through 9) where poly fit err is max
      INTEGER(LONG)                   :: STRESS_OUT_ERR_INDEX(MAX_STRESS_POINTS+1)



                                                           ! Array of %errs from subr POLYNOM_FIT_STRE_STRN (only NUM_PTS vals used)
      REAL(DOUBLE)                    :: STRESS_OUT_PCT_ERR(MAX_STRESS_POINTS+1)

      REAL(DOUBLE)                    :: PCT_ERR_MAX       ! Max value from array STRESS_OUT_PCT_ERR
      REAL(DOUBLE)                    :: C1,C2,D1,D2,E1,E2,F1,F2
      REAL(DOUBLE)                    :: EA0,EA1,EA2,EA3,EA4,EAMAX,EAMIN
      REAL(DOUBLE)                    :: EB0,EB1,EB2,EB3,EB4,EBMAX,EBMIN
                                                           ! Array of values from array STRESS for all stress points
      REAL(DOUBLE)                    :: STRESS_RAW(9,MAX_STRESS_POINTS+1)

                                                           ! Array of output stress values after surface fit
      REAL(DOUBLE)                    :: STRESS_OUT(9,MAX_STRESS_POINTS+1)
      REAL(DOUBLE)                    :: TEL(3,3)          ! Transformation matrix from cartesian local (L) to element (E) coordinates.

      ! OP2 stuff
      CHARACTER(8*BYTE)               :: TABLE_NAME   ! name of the op2 table name
      INTEGER(LONG)                   :: ITABLE       ! the subtable
      LOGICAL                         :: WRITE_NEU
      LOGICAL                         :: HAVE_SECTION_POINTS
      LOGICAL                         :: SHELL_GPSTRESS_RECOVERY

      INTRINSIC DABS, DMAX1, DMIN1, IAND
      ITABLE = 0
      TABLE_NAME = "OES ERR "

      WRITE_NEU = WRITE_NEU_STRE

! **********************************************************************************************************************************
! Process element stress output (STRESS) requests for all elems except composite shells

      OPT(1) = 'N'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'Y'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'N'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff


! Find out how many output requests were made for each element type.

      DO I=1,METYPE                                        ! Initialize the array containing the no. requests/elem.
         NELREQ(I) = 0
      ENDDO

      DO I=1,METYPE
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               IF (ETYPE(J) == ELMTYP(I)) THEN
! --- cbeam_stations begin --- !
                  IF (ETYPE(J) == 'BEAM    ') THEN
                     NUM_PTS_ELEM = PBEAM_NSTATIONS(EDAT(EPNT(J)+1))
                     IF (NUM_PTS_ELEM <= 0) NUM_PTS_ELEM = 5
                     IF (NUM_PTS_ELEM > NUM_PTS(I)) NUM_PTS(I) = NUM_PTS_ELEM
                  ELSE
                     SHELL_GPSTRESS_RECOVERY = GPSTRESS_REQ .AND.                                                                 &
                        ((ETYPE(J)(1:5) == 'TRIA3') .OR. (ETYPE(J)(1:5) == 'QUAD4') .OR. (ETYPE(J) == 'QUADR   '))
                     IF (SHELL_GPSTRESS_RECOVERY .AND. (ETYPE(J)(1:5) == 'TRIA3') .AND. (NUM_SEi(I) == 1)) THEN
                        NUM_PTS_ELEM = ELGP + 1
                     ELSE IF ((STRE_LOC == 'CORNER  ') .OR.                                                                       &
                         (STRE_LOC == 'GAUSS   ') .OR.                                                                            &
                         SHELL_GPSTRESS_RECOVERY .OR.                                                                             &
                         (ETYPE(J)(1:4) == 'HEXA') .OR.                                                                           &
                         (ETYPE(J)(1:5) == 'PYRAM') .OR.                                                                          &
                         (ETYPE(J)(1:5) == 'PENTA') .OR.                                                                          &
                         (ETYPE(J)(1:5) == 'TETRA') .OR.                                                                          &
                         (ETYPE(J)(1:5) == 'QUAD8') .OR. (ETYPE(J)(1:5) == 'TRIA6')) THEN
                        NUM_PTS_ELEM = NUM_SEi(I)
                     ELSE
                        NUM_PTS_ELEM = 1
                     ENDIF
                     NUM_PTS(I) = NUM_PTS_ELEM
                  ENDIF
! --- cbeam_stations end --- !
                  ELOUT_STRE = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRE_BIT))
                  IF (ELOUT_STRE > 0) THEN
                     NELREQ(I) = NELREQ(I) + NUM_PTS_ELEM
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      DO I=1,MAXREQ
         DO J=1,MOGEL
            OGEL(I,J) = ZERO
         ENDDO
         DO J=1,3
            SHELL_OUT_TE(J,1,I) = ZERO
            SHELL_OUT_TE(J,2,I) = ZERO
            SHELL_OUT_TE(J,3,I) = ZERO
         ENDDO
      ENDDO

! 101  FORMAT("*DEBUG:      ",A,"; ELEMENT_TYPE_INT=",I8,"; TABLE_NAME=",A)
!xx   IROW_MAT = 0
!xx   IROW_TXT = 0
      OT4_DESCRIPTOR = 'Element stress'
reqs5:DO I=1,METYPE
         IF (NELREQ(I) == 0) CYCLE reqs5
         NUM_OGEL_ROWS = 0
         NUM_OGEL      = 0

elems_5: DO J = 1,NELE

            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J) == ELMTYP(I)) THEN
               ELOUT_STRE = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRE_BIT))
               IF (ELOUT_STRE > 0) THEN
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE elems_5
                  ENDIF
                  CALL ELMDIS
                  IF (TYPE == 'BEAM    ') THEN
                     CALL CALC_ELEM_NODE_FORCES
                  ENDIF

! --- CBEAM_standard begin --- !
                   NUM_PTS_CUR = NUM_PTS(I)
                   IF (TYPE == 'BEAM    ') THEN
                      NUM_PTS_CUR = CBEAM_ACTIVE_NSTATIONS
                      IF (NUM_PTS_CUR <= 0) NUM_PTS_CUR = 1
                      IF ((DEBUG(233) > 0) .AND. ((EID == 14) .OR. (EID == 15))) THEN
                         WRITE(ERR,9233) EID, NUM_PTS(I), CBEAM_ACTIVE_NSTATIONS, NUM_PTS_CUR,                                     &
                                          (CBEAM_ACTIVE_XL(K), K=1,NUM_PTS_CUR)
                      ENDIF
                   ENDIF
! --- CBEAM_standard end --- !
                  SHELL_GPSTRESS_RECOVERY = GPSTRESS_REQ .AND.                                                                    &
                     ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4') .OR. (TYPE == 'QUADR   '))
                   DO M=1,NUM_PTS_CUR
                      RECOVERY_POINT = M
                      IF (SHELL_GPSTRESS_RECOVERY .AND. (TYPE(1:5) == 'TRIA3') .AND. (NUM_SEi(I) == 1)) THEN
                         RECOVERY_POINT = 1
                      ENDIF
                      CALL ELEM_STRE_STRN_ARRAYS ( RECOVERY_POINT )
                      STRESS_RAW(:,M) = STRESS(:)
                   ENDDO

! --- cbeam_stations begin --- !
                  IF (TYPE == 'BEAM    ') THEN
                     STRESS_OUT(:,:) = STRESS_RAW(:,:)
                  ELSE
                     STRESS_OUT(:,1) = STRESS(:)         ! Set STRESS_OUT for NUM_PTS(I) = 1
                  ENDIF
! --- cbeam_stations end --- !

                  IF ((STRE_LOC == 'CORNER  ') .OR.                                                                                &
                      (STRE_LOC == 'GAUSS   ') .OR.                                                                                &
                      SHELL_GPSTRESS_RECOVERY .OR.                                                                                 &
                      (TYPE(1:4) == 'HEXA') .OR.                                                                                   &
                      (TYPE(1:5) == 'PYRAM') .OR.                                                                                  &
                      (TYPE(1:5) == 'PENTA') .OR.                                                                                  &
                      (TYPE(1:5) == 'TETRA') .OR.                                                                                  &
                      (TYPE(1:5) == 'QUAD8')) THEN

                     IF (((TYPE == 'QUADR   ') .AND. (QUADRTYP == 'Q4RS    ')) .OR.                                               &
                         ((TYPE(1:5) == 'QUAD4') .AND. (QUAD4TYP == 'DSQK  '))) THEN
! Q4RS and DSQK recovery rows are evaluated directly at the requested output points.
                        STRESS_OUT(:,:) = STRESS_RAW(:,:)

                     ELSE IF ((TYPE(1:5) == 'QUAD4') .OR. (TYPE == 'QUADR   ')) THEN
                         CALL POLYNOM_FIT_STRE_STRN ( STRESS_RAW, 9, NUM_PTS_CUR, STRESS_OUT, STRESS_OUT_PCT_ERR,                  &
                                                      STRESS_OUT_ERR_INDEX, PCT_ERR_MAX )

                     ELSE IF (TYPE(1:5) == 'QUAD8') THEN
                         CALL POLYNOM_FIT_STRE_STRN ( STRESS_RAW, 9, NUM_PTS_CUR, STRESS_OUT, STRESS_OUT_PCT_ERR,                  &
                                                      STRESS_OUT_ERR_INDEX, PCT_ERR_MAX )

                                                           ! Transform stress from the cartesian local coordinate system to
                                                           ! the element coordinate system
                         DO M=2,NUM_PTS_CUR
                            CALL PLANE_COORD_TRANS_21( SHELL_STR_ANGLE( M ), TEL, '')
                            CALL TRANSFORM_SHELL_STR( TEL, STRESS_OUT(:,M), ONE)
                         ENDDO
                                                           ! Center stress is the average of corner stress in element coordinates.
                                                           ! This is how MSC does it.
                        STRESS_OUT(:,1) = (STRESS_OUT(:,2) + STRESS_OUT(:,3) + STRESS_OUT(:,4) + STRESS_OUT(:,5)) / FOUR

                     ELSE IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'TRIA6')) THEN
! Stresses are directly recovered at the triangular element output points.
                        STRESS_OUT(:,:) = STRESS_RAW(:,:)

                     ELSE IF ((TYPE(1:4) == 'HEXA') .OR.                                                                           &
                              (TYPE(1:5) == 'PYRAM') .OR.                                                                          &
                              (TYPE(1:5) == 'PENTA') .OR.                                                                          &
                              (TYPE(1:5) == 'TETRA')) THEN
! Stresses are directly evaluated at the corner grid points. If they are going to be evaluated at Gauss points
! then extrapolated to grid points, that should be done here, in POLYNOM_FIT_STRE_STRN, or in an equivalent subroutine.
                        STRESS_OUT(:,:) = STRESS_RAW(:,:)

                     ENDIF

                  ENDIF

 do_stress_pts:    DO M=1,NUM_PTS_CUR

                     DO K=1,9
                        STRESS(K) = STRESS_OUT(K,M)
                     ENDDO

                     CALL CALC_ELEM_STRESSES ( MAXREQ, NUM_OGEL, J, 'Y', 'N' )
                                                           ! If CB soln, write rows of OGEL, from CALC_ELEM_STRESSES, to OTM_STRE
                     IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN

                        CALL GET_STRESS_ITEM_DATA

                        IF ((TYPE == 'BAR     ') .OR. (TYPE == 'TRIA3   ') .OR. (TYPE == 'TRIA6   ') .OR.                          &
                            ((TYPE == 'QUAD4   ') .OR. (TYPE == 'QUADR   ')) .OR. (TYPE == 'SHEAR   ')) THEN
                           DO L=1,2
                              DO K=1,NUM_OTM_ENTRIES
                                 OT4_EROW = OT4_EROW + 1
                                 OTM_STRE(OT4_EROW,JVEC) = OGEL(NUM_OGEL-2+L,K)
                                 IF (JVEC == 1) THEN
                                    IF ((STRE_LOC == 'CORNER  ') .OR. (STRE_LOC == 'GAUSS   ')) THEN
                                       IF (M == 1) THEN
                                          IF ((TYPE(1:5) == 'QUAD4') .OR. (TYPE == 'QUADR   ')) THEN
                                             WRITE(TXT_STRE(OT4_EROW),9190) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID,                   &
                                                                           STRESS_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                          ELSE
                                             WRITE(TXT_STRE(OT4_EROW), 9193) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID,                  &
                                                                             STRESS_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                          ENDIF
                                       ELSE
                                          WRITE(TXT_STRE(OT4_EROW),9191) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, AGRID(M-1),          &
                                                                         STRESS_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                       ENDIF
                                    ELSE
                                       WRITE(TXT_STRE(OT4_EROW), 9192) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID,                        &
                                                                       STRESS_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                    ENDIF
                                 ENDIF
                              ENDDO
                           ENDDO
                        ELSE
                           DO K=1,NUM_OTM_ENTRIES
                              OT4_EROW = OT4_EROW + 1
                              OTM_STRE(OT4_EROW,JVEC) = OGEL(NUM_OGEL,K)
                              IF (JVEC == 1) THEN
                                 WRITE(TXT_STRE(OT4_EROW), 9193) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, STRESS_ITEM(K)
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDIF

                     IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (JVEC == 1) .AND. (OT4_EROW >= 1)) THEN
                        DO K=1,OTMSKIP                        ! Write OTMSKIP blank separator lines
                           OT4_EROW = OT4_EROW + 1
                           WRITE(TXT_STRE(OT4_EROW), 9199)
                        ENDDO
                     ENDIF

                     NUM_OGEL_ROWS = NUM_OGEL_ROWS + 1
                     EID_OUT_ARRAY(NUM_OGEL_ROWS,1) = EID
! --- cbeam_stations begin --- !
                     IF (TYPE == 'BEAM    ') THEN
                        CBEAM_XL_OUT(NUM_OGEL_ROWS) = CBEAM_ACTIVE_XL(M)
                     ELSE
                        CBEAM_XL_OUT(NUM_OGEL_ROWS) = ZERO
                     ENDIF
! --- cbeam_stations end --- !
                     SHELL_OUT_TE(1:3,1:3,NUM_OGEL_ROWS) = ZERO
                     SHELL_STRESS_IN_LOCAL(NUM_OGEL_ROWS) = .FALSE.
                     IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'TRIA6') .OR. (TYPE(1:5) == 'QUAD4') .OR.                     &
                         (TYPE == 'QUADR   ') .OR.                                                                                 &
                         (TYPE(1:5) == 'QUAD8')) THEN
                        SHELL_OUT_TE(1:3,1:3,NUM_OGEL_ROWS) = TE(1:3,1:3)
                        CALL SET_SHELL_STRESS_BASIS_FOR_OUTPUT ( NUM_OGEL_ROWS, M )
                     ENDIF
                     GID_OUT_ARRAY(NUM_OGEL_ROWS,1) = 0
                     IF ((STRE_LOC == 'CORNER  ') .OR. (STRE_LOC == 'GAUSS   ') .OR. SHELL_GPSTRESS_RECOVERY) THEN
                        IF ((TYPE(1:5) == 'QUAD4') .OR. (TYPE == 'QUADR   ')) THEN
                           POLY_FIT_ERR(NUM_OGEL_ROWS)       = STRESS_OUT_PCT_ERR(M)
                           POLY_FIT_ERR_INDEX(NUM_OGEL_ROWS) = STRESS_OUT_ERR_INDEX(M)
                        ENDIF
                     ENDIF
                     DO K=1,ELGP
                        GID_OUT_ARRAY(NUM_OGEL_ROWS,K+1) = AGRID(K)
                     ENDDO

                  ENDDO do_stress_pts

                  IF (ETYPE(J)(1:5) /='USER1') THEN
                     IF (NUM_OGEL_ROWS == NELREQ(I)) THEN
                        CALL CHK_OGEL_ZEROS ( NUM_OGEL )
                        CALL SET_OES_TABLE_NAME(TYPE, TABLE_NAME, ITABLE)
                        CALL WRITE_ELEM_STRESSES ( JVEC, NUM_OGEL_ROWS, IHDR, NUM_PTS_CUR, ITABLE )
                        EXIT
                     ENDIF
                  ENDIF

               ENDIF

            ENDIF

         ENDDO elems_5

      ENDDO reqs5

      IF ((TABLE_NAME .NE. "OES ERR ") .AND. (ITABLE < 0)) THEN
        CALL END_OP2_TABLE(ITABLE)
      ENDIF
!===========================
      IF (WRITE_NEU .AND. (ANY_STRE_OUTPUT > 0)) THEN

         NDUM = 0
! --- neu_upgrade begin --- !
         NUM_FROWS= 0                                      ! Write out BEAM stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCBEAM, 14, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:4) == 'BEAM') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=1,14
                     FEMAP_EL_VECS(NUM_FROWS,K) = ZERO
                  ENDDO
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  IF (TYPE == 'BEAM    ') THEN
                     OPT(5) = 'Y'
                  ELSE
                     OPT(5) = 'N'
                  ENDIF
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
! --- CBEAM_standard begin --- !
! Beam stress recovery in FEMAP/NEU needs the same element nodal force reconstruction
! used by the standard F06 path; otherwise the beam stress vectors stay effectively zero.
                  CALL CALC_ELEM_NODE_FORCES
                  HAVE_SECTION_POINTS = .FALSE.
                  DO K=1,8
                     IF (DABS(ZS(K)) > ZERO) THEN
                        HAVE_SECTION_POINTS = .TRUE.
                        EXIT
                     ENDIF
                  ENDDO
                  IF (HAVE_SECTION_POINTS) THEN
                     C1 = ZS(1); C2 = ZS(2)
                     D1 = ZS(3); D2 = ZS(4)
                     E1 = ZS(5); E2 = ZS(6)
                     F1 = ZS(7); F2 = ZS(8)
                  ELSE
                     C1 = ZERO; C2 = ZERO
                     D1 = ZERO; D2 = ZERO
                     E1 = ZERO; E2 = ZERO
                     F1 = ZERO; F2 = ZERO
                  ENDIF
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  EA0   = STRESS(1)
                  EA1   = EA0 - (C1*STRESS(2) + C2*STRESS(3))
                  EA2   = EA0 - (D1*STRESS(2) + D2*STRESS(3))
                  EA3   = EA0 - (E1*STRESS(2) + E2*STRESS(3))
                  EA4   = EA0 - (F1*STRESS(2) + F2*STRESS(3))
                  EAMAX = DMAX1(EA1,EA2,EA3,EA4)
                  EAMIN = DMIN1(EA1,EA2,EA3,EA4)
                  FEMAP_EL_VECS(NUM_FROWS, 1) = EA1
                  FEMAP_EL_VECS(NUM_FROWS, 3) = EA2
                  FEMAP_EL_VECS(NUM_FROWS, 5) = EA3
                  FEMAP_EL_VECS(NUM_FROWS, 7) = EA4
                  FEMAP_EL_VECS(NUM_FROWS, 9) = EAMAX
                  FEMAP_EL_VECS(NUM_FROWS,11) = EAMIN
                  IF (CBEAM_ACTIVE_NSTATIONS > 1) THEN
                     CALL ELEM_STRE_STRN_ARRAYS ( CBEAM_ACTIVE_NSTATIONS )
                  ELSE
                     CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  ENDIF
                  EB0   = STRESS(1)
                  EB1   = EB0 - (C1*STRESS(2) + C2*STRESS(3))
                  EB2   = EB0 - (D1*STRESS(2) + D2*STRESS(3))
                  EB3   = EB0 - (E1*STRESS(2) + E2*STRESS(3))
                  EB4   = EB0 - (F1*STRESS(2) + F2*STRESS(3))
                  EBMAX = DMAX1(EB1,EB2,EB3,EB4)
                  EBMIN = DMIN1(EB1,EB2,EB3,EB4)
                  FEMAP_EL_VECS(NUM_FROWS, 2) = EB1
                  FEMAP_EL_VECS(NUM_FROWS, 4) = EB2
                  FEMAP_EL_VECS(NUM_FROWS, 6) = EB3
                  FEMAP_EL_VECS(NUM_FROWS, 8) = EB4
                  FEMAP_EL_VECS(NUM_FROWS,10) = EBMAX
                  FEMAP_EL_VECS(NUM_FROWS,12) = EBMIN
! --- CBEAM_standard end --- !
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'BEAM    ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA
! --- neu_upgrade end --- !

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out BUSH stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCBUSH, 6, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:4) == 'BUSH') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCBUSH, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'BUSH   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out ELAS1 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCELAS1, 2, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'ELAS1') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCELAS1, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'ELAS1   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out ELAS2 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCELAS2, 2, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'ELAS2') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCELAS2, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'ELAS2   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out ELAS3 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCELAS3, 2, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'ELAS3') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCELAS3, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'ELAS3   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out ELAS4 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCELAS4, 2, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'ELAS4') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCELAS4, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'ELAS4   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out ROD stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCROD, 4, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:3) == 'ROD') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCROD, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'ROD     ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out BAR stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCBAR, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:3) == 'BAR') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCBAR, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'BAR     ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TRIA3K stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTRIA3K, 24, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'TRIA3K') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCTRIA3K, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'TRIA3K  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TRIA3 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTRIA3, 24, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'TRIA3 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCTRIA3, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'TRIA3   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out QUAD4K stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCQUAD4K, 24, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'QUAD4K') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCQUAD4K, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'QUAD4K  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out QUAD4/CQUADR stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCQUAD4 + NCQUADR, 24, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF ((ETYPE(J)(1:6) == 'QUAD4 ') .OR. (ETYPE(J) == 'QUADR   ')) THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCQUAD4 + NCQUADR, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'QUADR   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out PYRAM5 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NPYRAM5, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'PYRAM5 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NPYRAM5, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'PYRAM5  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out PYRAM14 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NPYRAM14, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'PYRAM14') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NPYRAM14, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'PYRAM14 ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out HEXA8 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCHEXA8, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'HEXA8 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCHEXA8, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'HEXA8   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out HEXA20 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCHEXA20, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'HEXA20') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCHEXA20, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'HEXA20  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out PENTA6 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCPENTA6, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'PENTA6 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCPENTA6, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'PENTA6  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out PENTA15 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCPENTA15, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'PENTA15') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCPENTA15, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'PENTA15 ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TETRA4 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTETRA4, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'TETRA4 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCTETRA4, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'TETRA4  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TETRA10 stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTETRA10, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'TETRA10') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCTETRA10, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'TETRA10 ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out SHEAR stresses
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCSHEAR, 24, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'SHEAR') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRESSES ( NCSHEAR, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRE_VECS ( 'SHEAR   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

      ENDIF

      IF (IERROR > 0) THEN
         REQUEST = 'ELEMENT STRESS'
         WRITE(ERR,9201) TYPE, REQUEST, EID
         WRITE(F06,9201) TYPE, REQUEST, EID
      ENDIF



      RETURN

! **********************************************************************************************************************************
 9190 FORMAT(I8,1X,A,A8,I8,4X,'CENTER',9X,A20)

 9191 FORMAT(I8,1X,A,A8,I8,4X,'GRID',I8,3X,A20)

 9192 FORMAT(I8,1X,A,A8,I8,4X,A20)

 9193 FORMAT(I8,1X,A,A8,I8,19X,A20)

 9199 FORMAT(' ')

 9201 FORMAT(' *ERROR  9201: DUE TO ABOVE LISTED ERRORS, CANNOT CALCULATE ',A,' REQUESTS FOR ',A,' ELEMENT ID = ',I8)

 9233 FORMAT(' *CBEAM STRE DEBUG: EID=',I8,' NUM_PTS(I)=',I8,' ACTIVE=',I8,' CUR=',I8,' XL=',21(1X,ES12.5))

 1001 FORMAT(2(I8,','),'       1,')
 1002 FORMAT(A)
 1003 FORMAT(3(1ES17.6,','))
 1004 FORMAT(10(I8,','))
 1005 FORMAT(2(I8,','),'       1,       7,',/,'       1,       1,       1')
 1006 FORMAT(I8,',',1ES17.6,',')
 1007 FORMAT('      -1,     0.          ,')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE SET_SHELL_STRESS_BASIS_FOR_OUTPUT ( ROW_NUM, POINT_NUM )

      INTEGER(LONG), INTENT(IN)       :: ROW_NUM
      INTEGER(LONG), INTENT(IN)       :: POINT_NUM

      LOGICAL                         :: OK
      REAL(DOUBLE)                    :: BASIS(3,3)

      IF (ROW_NUM <= 0) RETURN

      IF (TYPE(1:5) == 'TRIA3') THEN
         IF ((TRIA3TYP == 'T3FF  ') .OR. (TRIA3TYP == 'MITC3+') .OR.                                                               &
             (TRIARTYP == 'T3FFD   ') .OR. (TRIARTYP == 'MITC3+HB')) THEN
            CALL BUILD_TRIA_STRESS_BASIS ( BASIS, OK )
            IF (OK) SHELL_OUT_TE(1:3,1:3,ROW_NUM) = BASIS(1:3,1:3)
            SHELL_STRESS_IN_LOCAL(ROW_NUM) = OK
         ENDIF
      ELSE IF ((TYPE(1:5) == 'QUAD4') .AND. (QUAD4TYP == 'DSQK  ')) THEN
         SHELL_STRESS_IN_LOCAL(ROW_NUM) = .TRUE.
      ELSE IF ((TYPE == 'QUADR   ') .AND. (QUADRTYP == 'MITC4PD ')) THEN
         SHELL_STRESS_IN_LOCAL(ROW_NUM) = .TRUE.
      ELSE IF ((TYPE == 'QUADR   ') .AND. (QUADRTYP == 'Q4RS    ')) THEN
         CALL BUILD_QUAD_STRESS_BASIS ( POINT_NUM, BASIS, OK )
         IF (OK) SHELL_OUT_TE(1:3,1:3,ROW_NUM) = BASIS(1:3,1:3)
         SHELL_STRESS_IN_LOCAL(ROW_NUM) = OK
      ENDIF

      END SUBROUTINE SET_SHELL_STRESS_BASIS_FOR_OUTPUT

! ##################################################################################################################################

      SUBROUTINE BUILD_TRIA_STRESS_BASIS ( BASIS, OK )

      REAL(DOUBLE), INTENT(OUT)       :: BASIS(3,3)
      LOGICAL, INTENT(OUT)            :: OK

      INTEGER(LONG)                   :: II
      REAL(DOUBLE)                    :: E1(3), E2(3), E3(3), V12(3), V13(3)

      BASIS = ZERO
      OK = .FALSE.

      DO II=1,3
         V12(II) = XEB(2,II) - XEB(1,II)
         V13(II) = XEB(3,II) - XEB(1,II)
      ENDDO
      IF (VEC_NORM(V12) <= 1.0D-14) RETURN
      E1 = V12 / VEC_NORM(V12)
      CALL CROSS3 ( E1, V13, E3 )
      IF (VEC_NORM(E3) <= 1.0D-14) RETURN
      E3 = E3 / VEC_NORM(E3)
      CALL CROSS3 ( E3, E1, E2 )

      BASIS(1,1:3) = E1
      BASIS(2,1:3) = E2
      BASIS(3,1:3) = E3
      OK = .TRUE.

      END SUBROUTINE BUILD_TRIA_STRESS_BASIS

! ##################################################################################################################################

      SUBROUTINE BUILD_QUAD_STRESS_BASIS ( POINT_NUM, BASIS, OK )

      INTEGER(LONG), INTENT(IN)       :: POINT_NUM
      REAL(DOUBLE), INTENT(OUT)       :: BASIS(3,3)
      LOGICAL, INTENT(OUT)            :: OK

      INTEGER(LONG)                   :: II
      REAL(DOUBLE)                    :: DNXI(4), DNETA(4), XI, ETA
      REAL(DOUBLE)                    :: E1(3), E2(3), E3(3), G1(3), G2(3)

      BASIS = ZERO
      OK = .FALSE.

      IF (POINT_NUM == 2) THEN
         XI = -ONE
         ETA = -ONE
      ELSE IF (POINT_NUM == 3) THEN
         XI = ONE
         ETA = -ONE
      ELSE IF (POINT_NUM == 4) THEN
         XI = ONE
         ETA = ONE
      ELSE IF (POINT_NUM == 5) THEN
         XI = -ONE
         ETA = ONE
      ELSE
         XI = ZERO
         ETA = ZERO
      ENDIF

      DNXI(1) = -0.25D0*(ONE - ETA)
      DNXI(2) =  0.25D0*(ONE - ETA)
      DNXI(3) =  0.25D0*(ONE + ETA)
      DNXI(4) = -0.25D0*(ONE + ETA)
      DNETA(1) = -0.25D0*(ONE - XI)
      DNETA(2) = -0.25D0*(ONE + XI)
      DNETA(3) =  0.25D0*(ONE + XI)
      DNETA(4) =  0.25D0*(ONE - XI)

      G1 = ZERO
      G2 = ZERO
      DO II=1,4
         G1(1:3) = G1(1:3) + DNXI(II)*XEB(II,1:3)
         G2(1:3) = G2(1:3) + DNETA(II)*XEB(II,1:3)
      ENDDO
      IF (VEC_NORM(G1) <= 1.0D-14) RETURN
      E1 = G1 / VEC_NORM(G1)
      CALL CROSS3 ( E1, G2, E3 )
      IF (VEC_NORM(E3) <= 1.0D-14) CALL CROSS3 ( G1, G2, E3 )
      IF (VEC_NORM(E3) <= 1.0D-14) RETURN
      E3 = E3 / VEC_NORM(E3)
      CALL CROSS3 ( E3, E1, E2 )

      BASIS(1,1:3) = E1
      BASIS(2,1:3) = E2
      BASIS(3,1:3) = E3
      OK = .TRUE.

      END SUBROUTINE BUILD_QUAD_STRESS_BASIS

! ##################################################################################################################################

      SUBROUTINE CROSS3 ( A, B, C )

      REAL(DOUBLE), INTENT(IN)        :: A(3), B(3)
      REAL(DOUBLE), INTENT(OUT)       :: C(3)

      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)

      END SUBROUTINE CROSS3

! ##################################################################################################################################

      REAL(DOUBLE) FUNCTION VEC_NORM ( V )

      REAL(DOUBLE), INTENT(IN)        :: V(3)

      VEC_NORM = DSQRT(DOT_PRODUCT(V,V))

      END FUNCTION VEC_NORM

! ##################################################################################################################################

      SUBROUTINE GET_STRESS_ITEM_DATA

      IMPLICIT NONE

      INTEGER(LONG)                   :: II               ! DO loop index

! **********************************************************************************************************************************
      DO II=1,18
         STRESS_ITEM(II)(1:) = ' '
      ENDDO

      IF       (TYPE(1:5) == 'ELAS1') THEN
         NUM_OTM_ENTRIES = 1
         STRESS_ITEM( 1) = 'Spring elem stress  '

      ELSE IF  (TYPE(1:3) == 'ROD'  ) THEN
         NUM_OTM_ENTRIES = 4
         STRESS_ITEM( 1) = 'Axial Stress        '
         STRESS_ITEM( 2) = 'MS - Axial          '
         STRESS_ITEM( 3) = 'Torsional Stress    '
         STRESS_ITEM( 4) = 'MS - Torsion        '

      ELSE IF  (TYPE(1:3) == 'BAR'  ) THEN
         NUM_OTM_ENTRIES = 9
         STRESS_ITEM( 1) = 'SA1: Stress Pt1 EndA'  ;  STRESS_ITEM(10) = 'SB1: Stress Pt1 EndB'
         STRESS_ITEM( 2) = 'SA2: Stress Pt2 EndA'  ;  STRESS_ITEM(11) = 'SB2: Stress Pt2 EndB'
         STRESS_ITEM( 3) = 'SA3: Stress Pt3 EndA'  ;  STRESS_ITEM(12) = 'SB3: Stress Pt3 EndB'
         STRESS_ITEM( 4) = 'SA4: Stress Pt4 EndA'  ;  STRESS_ITEM(13) = 'SB4: Stress Pt4 EndB'
         STRESS_ITEM( 5) = 'Axial Stress        '  ;  STRESS_ITEM(14) = 'Axial stress        '
         STRESS_ITEM( 6) = 'SA-Max              '  ;  STRESS_ITEM(15) = 'SB-Max              '
         STRESS_ITEM( 7) = 'SA-Min              '  ;  STRESS_ITEM(16) = 'SB-Min              '
         STRESS_ITEM( 8) = 'MS-Tension          '  ;  STRESS_ITEM(17) = 'MS-Compression      '
         STRESS_ITEM( 9) = 'Torsional Stress    '  ;  STRESS_ITEM(18) = 'MS-Torsion          '

      ELSE IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'TRIA6') .OR. (TYPE(1:5) == 'QUAD4') .OR. (TYPE == 'QUADR   ')) THEN
         NUM_OTM_ENTRIES = 10
         STRESS_ITEM( 1) = 'Fibre Dist      -Z1 '  ;  STRESS_ITEM(11) = 'Fibre Dist      +Z1 '
         STRESS_ITEM( 2) = 'Normal X Stress -Z1 '  ;  STRESS_ITEM(12) = 'Normal X Stress +Z1 '
         STRESS_ITEM( 3) = 'Normal Y Stress -Z1 '  ;  STRESS_ITEM(13) = 'Normal Y Stress +Z1 '
         STRESS_ITEM( 4) = 'Shear XY Stress -Z1 '  ;  STRESS_ITEM(14) = 'Shear XY Stress +Z1 '
         STRESS_ITEM( 5) = 'Princ Angle  at -Z1 '  ;  STRESS_ITEM(15) = 'Princ Angle  at +Z1 '
         STRESS_ITEM( 6) = 'Major Stress at -Z1 '  ;  STRESS_ITEM(16) = 'Major Stress at +Z1 '
         STRESS_ITEM( 7) = 'Minor Stress at -Z1 '  ;  STRESS_ITEM(17) = 'Minor Stress at +Z1 '
         IF      (STRE_OPT(1:8) == 'VONMISES') THEN
            STRESS_ITEM( 8) = 'von Mises XY at -Z1 '  ;  STRESS_ITEM(18) = 'von Mises XY at +Z1 '
         ELSE IF (STRE_OPT(1:4) == 'MAXS'    ) THEN
            STRESS_ITEM( 8) = 'Max Shear XY at -Z1 '  ;  STRESS_ITEM(18) = 'Max Shear XY at +Z1 '
         ELSE
            STRESS_ITEM( 8) = '**** undefined **** '  ;  STRESS_ITEM(18) = '**** undefined **** '
         ENDIF
         STRESS_ITEM( 9) = 'Shear XZ Stress avg '  ;  STRESS_ITEM(19) = 'Shear XZ Stress avg '
         STRESS_ITEM(10) = 'Shear YZ Stress avg '  ;  STRESS_ITEM(20) = 'Shear XZ Stress avg '

      ELSE IF  (TYPE(1:6) == 'USERIN') THEN
         NUM_OTM_ENTRIES = 9
         STRESS_ITEM( 1) = '*** Not Defined ****'
         STRESS_ITEM( 2) = '*** Not Defined ****'
         STRESS_ITEM( 3) = '*** Not Defined ****'
         STRESS_ITEM( 4) = '*** Not Defined ****'
         STRESS_ITEM( 5) = '*** Not Defined ****'
         STRESS_ITEM( 6) = '*** Not Defined ****'
         STRESS_ITEM( 7) = '*** Not Defined ****'
         STRESS_ITEM( 8) = '*** Not Defined ****'
         STRESS_ITEM( 9) = '*** Not Defined ****'

      ELSE IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PYRAM') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
         NUM_OTM_ENTRIES = 8
         STRESS_ITEM( 1) = 'Normal x Stress     '
         STRESS_ITEM( 2) = 'Normal y Stress     '
         STRESS_ITEM( 3) = 'Normal z Stress     '
         STRESS_ITEM( 4) = 'Shear xy Stress     '
         STRESS_ITEM( 5) = 'Shear yz Stress     '
         STRESS_ITEM( 6) = 'Shear zx Stress     '
         STRESS_ITEM( 7) = 'Oct Direct Stress   '
         STRESS_ITEM( 8) = 'Oct Shear Stress    '

      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE GET_STRESS_ITEM_DATA

!====================================================================================================

      END SUBROUTINE OFP3_STRE_NO_PCOMP

