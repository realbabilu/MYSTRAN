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

      SUBROUTINE OFP3_ELFE_2D ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )

! Processes element engr force output requests for 2D (TRIA3, QUAD4, SHEAR) elements for one subcase. Results go into array OGEL
! for later output in LINK9

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELOUT_ELFE_BIT, FATAL_ERR, IBIT, INT_SC_NUM, MBUG, MOGEL,                   &
                                         WARN_ERR, NELE, NCQUAD4, NCQUAD4K, NCQUADR, NCSHEAR, NCTRIA3, NCTRIA3K, SOL_NAME,         &
                                         MAX_STRESS_POINTS
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, FOUR
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_NUMS, FEMAP_EL_VECS
      USE PARAMS, ONLY                :  OTMSKIP, QUAD4TYP, QUADRTYP
      USE LINK9_STUFF, ONLY           :  WRITE_NEU_ELFO
      use model_stuf, only            :  pcomp_props
      USE MODEL_STUF, ONLY            :  ANY_ELFE_OUTPUT, EDAT, EPNT, ETYPE, FCONV, EID, ELMTYP, ELOUT, METYPE, NUM_EMG_FATAL_ERRS,&
                                         PLY_NUM, TYPE, STRESS, SHELL_STR_ANGLE, NUM_SEi, ELGP, AGRID, GRID_ID, RGRID
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  FORC_LOC, GPSTRESS_REQ, NUM_GP_SURFACE, MAX_GP_SURFACES, GP_SURFACE_IDS,                 &
                                         GP_SURFACE_NORMAL_MODE
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, GID_OUT_ARRAY, MAXREQ, OGEL
      USE OUTPUT4_MATRICES, ONLY      :  OTM_ELFE, TXT_ELFE
      USE GPSTRESS_SURFACE_UTILS, ONLY:  GPSTRESS_COLLECT_SURFACE_PATCH
      USE FAST_OUTPUT_FORMATTERS, ONLY:  FAST_FMT_I8_RJ, FAST_FMT_F06_E14_6

      USE GET_ARRAY_ROW_NUM_Interface
      USE PLANE_COORD_TRANS_21_Interface
      USE TRANSFORM_SHELL_STR_Interface
      USE OFP3_ELFE_2D_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'OFP3_ELFE_2D'
      CHARACTER( 1*BYTE), PARAMETER   :: IHDR      = 'Y'   ! An input to subr WRITE_GRID_OUTPUTS, called herein
      CHARACTER(20*BYTE)              :: FORCE_ITEM(8)     ! Char description of element engineering forces
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option indicators for subr EMG, called herein
      CHARACTER(31*BYTE)              :: OT4_DESCRIPTOR    ! Descriptor for rows of OT4 file
      CHARACTER(30*BYTE)              :: REQUEST           ! Text for error message

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID      ! Set ID for FEMAP output
      INTEGER(LONG), INTENT(IN)       :: ITE               ! Unit number for text files for OTM row descriptors
      INTEGER(LONG), INTENT(IN)       :: JVEC              ! Solution vector number
     !INTEGER(LONG), INTENT(INOUT)    :: ITABLE            ! the op2 subtable number, should be -3, -5, ...
      INTEGER(LONG), INTENT(INOUT)    :: OT4_EROW          ! Row number in OT4 file for elem related OTM descriptors
      INTEGER(LONG)                   :: ELOUT_ELFE        ! If > 0, there are ELFORCE(ENGR) requests for some elems
      INTEGER(LONG)                   :: I,J,K,M           ! DO loop indices
      INTEGER(LONG)                   :: IERROR       = 0  ! Local error count
!xx   INTEGER(LONG)                   :: IROW_MAT          ! Row number in OTM's
!xx   INTEGER(LONG)                   :: IROW_TXT          ! Row number in OTM text file
      INTEGER(LONG)                   :: NELREQ(METYPE)    ! Count of the no. of requests for ELFORCE(NODE or ENGR) or STRESS
      INTEGER(LONG)                   :: NUM_OGEL_ROWS     ! No. elem points processed prior to writing results to F06 file
      INTEGER(LONG)                   :: NUM_FROWS         ! No. elems processed for FEMAP
      INTEGER(LONG)                   :: NUM_OGEL          ! No. rows written to array OGEL prior to writing results to F06 file
!                                                            (this can be > NUM_OGEL_ROWS since more than 1 row is written to OGEL
!                                                            for ELFORCE(NODE) - elem nodal forces)
      INTEGER(LONG)                   :: NUM_PTS(METYPE)   ! Num diff force points for one element
      integer(long)                   :: num_pcomp_elems   ! number of elements that are composites (used to prevent output of engr
!                                                            forces for PCOMP elems until I fix that output)

                                                           ! Stress index (1 through 9) where poly fit err is max
      INTEGER(LONG)                   :: STRESS_OUT_ERR_INDEX(MAX_STRESS_POINTS)

                                                           ! Array of %errs from subr POLYNOM_FIT_STRE_STRN (only NUM_PTS vals used)
      REAL(DOUBLE)                    :: STRESS_OUT_PCT_ERR(MAX_STRESS_POINTS)

      REAL(DOUBLE)                    :: PCT_ERR_MAX       ! Max value from array STRESS_OUT_PCT_ERR

      REAL(DOUBLE)                    :: TEL(3,3)          ! Transformation matrix from cartesian local (L) to element (E) coordinates.
                                                           ! Array of values from array STRESS for all stress points
      REAL(DOUBLE)                    :: STRESS_RAW(9,MAX_STRESS_POINTS)
                                                           ! Array of output stress values after surface fit
      REAL(DOUBLE)                    :: STRESS_OUT(9,MAX_STRESS_POINTS)


      ! OP2 parameters
      INTEGER(LONG)                   :: ITABLE            ! the op2 subtable number
      CHARACTER(8*BYTE)               :: TABLE_NAME        ! the op2 table name

      LOGICAL                         :: WRITE_NEU
      LOGICAL                         :: DIRECT_SHELL_RECOVERY
      INTRINSIC IAND

! **********************************************************************************************************************************
!     Initialize
      TABLE_NAME = "OEF ERR "
      ITABLE = 0

      WRITE_NEU = WRITE_NEU_ELFO

! **********************************************************************************************************************************
! Process element engineering force requests for plate and USERIN elements.
! For the MIN3,4 elements, the stiffness matrix has to be generated to get FCONV(3). Therefore, set OPT(4) to 'N' initially, and if
! the element being processed is a MIN3 or MIN4, reset OPT(4) to 'Y' in the calculation loop prior to EMG call.

      OPT(1) = 'N'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'Y'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'N'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff

      FORCE_ITEM(1) = 'Nxx: Normal x Force '
      FORCE_ITEM(2) = 'Nyy: Normal y Force '
      FORCE_ITEM(3) = 'Nxy: Shear xy Force '
      FORCE_ITEM(4) = 'Mxx: Moment x Plane '
      FORCE_ITEM(5) = 'Myy: Moment y Plane '
      FORCE_ITEM(6) = 'Mxy: Twist Mom xy   '
      FORCE_ITEM(7) = 'Qx : Transv Shear x '
      FORCE_ITEM(8) = 'Qy : Transv Shear y '

! Find out how many output requests were made for each element type.

      DO I=1,METYPE                                        ! Initialize the array containing the number of requests per element
         NELREQ(I) = 0
      ENDDO

      num_pcomp_elems = 0                                  ! Remove lower case code when I fix engr force output for PCOMP's
      DO I=1,METYPE
         DO J=1,NELE
            IF((ETYPE(J)(1:5) == 'TRIA3') .OR. (ETYPE(J)(1:5) == 'QUAD4') .OR. (ETYPE(J) == 'QUADR   ') .OR.                      &
               (ETYPE(J)(1:5) == 'QUAD8') .OR.                                                                                      &
               (ETYPE(J)(1:5) == 'SHEAR') .OR. (ETYPE(J)(1:6) == 'USERIN')) THEN
               IF (ETYPE(J) == ELMTYP(I)) THEN
                  call is_elem_pcomp_props ( j )
                  if (pcomp_props == 'N') then
                     IF ((FORC_LOC == 'CORNER  ') .OR.                                                                             &
                         (ETYPE(J)(1:5) == 'QUAD8')) THEN
                        NUM_PTS(I) = NUM_SEi(I)
                     ELSE
                        NUM_PTS(I) = 1
                     ENDIF
                     ELOUT_ELFE = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_ELFE_BIT))
                     IF (ELOUT_ELFE > 0) THEN
                        NELREQ(I) = NELREQ(I) + NUM_PTS(I)
                     ENDIF
                  else
                     num_pcomp_elems = num_pcomp_elems + 1
                  endif
               ENDIF
            ENDIF
         ENDDO
      ENDDO

!xx   DO I=1,METYPE                                        ! Engr force requests for ELAS elems not honored. Give message
!xx      IF ((ELMTYP(I)(1:4) == 'ELAS') .AND. (NELREQ(I) > 0)) THEN
!xx         WARN_ERR = WARN_ERR + 1
!xx         WRITE(ERR,9204) NELREQ(I), ELMTYP(I)
!xx         IF (SUPWARN == 'N') THEN
!xx            WRITE(F06,9204) NELREQ(I), ELMTYP(I)
!xx         ENDIF
!xx      ENDIF
!xx   ENDDO

      DO I=1,MAXREQ
         DO J=1,MOGEL
            OGEL(I,J) = ZERO
         ENDDO
      ENDDO

!xx   IROW_MAT = 0
!xx   IROW_TXT = 0
      OT4_DESCRIPTOR = 'Element engineering force'
reqs3:DO I=1,METYPE
         IF (NELREQ(I) == 0) CYCLE reqs3
         NUM_OGEL_ROWS = 0
         NUM_OGEL = 0

elems_3: DO J = 1,NELE
            call is_elem_pcomp_props ( j )                 ! Remove lower case code when I fix engr force output for PCOMP's
            if (pcomp_props == 'N') then
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF((ETYPE(J)(1:5) == 'TRIA3') .OR. (ETYPE(J)(1:5) == 'QUAD4') .OR. (ETYPE(J) == 'QUADR   ') .OR.                   &
                  (ETYPE(J)(1:5) == 'QUAD8') .OR.                                                                                   &
                  (ETYPE(J)(1:5) == 'SHEAR') .OR. (ETYPE(J)(1:6) == 'USERIN')) THEN
                  IF (ETYPE(J) == ELMTYP(I)) THEN
                     ELOUT_ELFE = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_ELFE_BIT))
                     IF (ELOUT_ELFE > 0) THEN
                        IF((ETYPE(J)(1:5) == 'TRIA3') .OR. (ETYPE(J)(1:5) == 'QUAD4') .OR. (ETYPE(J) == 'QUADR   ') .OR.          &
                           (ETYPE(J)(1:5) == 'QUAD8') .OR.                                                                          &
                           (ETYPE(J)(1:5) == 'SHEAR')) THEN
                           OPT(4) = 'Y'
                        ENDIF
                        DO K=0,MBUG-1
                           WRT_BUG(K) = 0
                        ENDDO
                        PLY_NUM = 0                        ! 'N' in call to EMG means do not write to BUG file
                        CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                        IF (NUM_EMG_FATAL_ERRS > 0) THEN
                           IERROR = IERROR + 1
                           CYCLE elems_3
                        ENDIF
                        OPT(4) = 'N'
                        CALL ELMDIS


                        DO M=1,NUM_PTS(I)                  ! Gauss point stress
                           CALL ELEM_STRE_STRN_ARRAYS ( M )
                           STRESS_RAW(:,M) = STRESS(:)
                        ENDDO

                        STRESS_OUT(:,1) = STRESS_RAW(:,1)  ! Set STRAIN_OUT for NUM_PTS(I) = 1
                        DIRECT_SHELL_RECOVERY = ((TYPE == 'QUADR   ') .AND. ((QUADRTYP == 'DKM24EA ') .OR.                       &
                                                                              (QUADRTYP == 'DKM24AU ') .OR.                       &
                                                                              (QUADRTYP == 'SIMO    ') .OR.                       &
                                                                              (QUADRTYP == 'MBP1C0  ') .OR.                       &
                                                                              (QUADRTYP == 'Q4EASANS'))) .OR.                    &
                                               ((TYPE == 'QUAD4   ') .AND. (QUAD4TYP == 'DKMQ20'))

                        IF ((FORC_LOC == 'CORNER  ') .OR.                                                                          &
                            (ETYPE(J)(1:5) == 'QUAD8')) THEN

                           IF ((TYPE(1:5) == 'QUAD4') .OR. (TYPE == 'QUADR   ')) THEN

                                                           ! Extrapolate stress to corners
                              IF (DIRECT_SHELL_RECOVERY) THEN
! Rows 1:5 and 7:9 remain Gauss samples. Row 6 is supplied by the
! DKMQ/Simo kernels as direct center/nodal twisting recovery.
                                 CALL POLYNOM_FIT_STRE_STRN ( STRESS_RAW, 9, NUM_PTS(I), STRESS_OUT, STRESS_OUT_PCT_ERR,           &
                                       STRESS_OUT_ERR_INDEX, PCT_ERR_MAX )
                                 STRESS_OUT(6,1:NUM_PTS(I)) = STRESS_RAW(6,1:NUM_PTS(I))
                              ELSE
                                 CALL POLYNOM_FIT_STRE_STRN ( STRESS_RAW, 9, NUM_PTS(I), STRESS_OUT, STRESS_OUT_PCT_ERR,           &
                                       STRESS_OUT_ERR_INDEX, PCT_ERR_MAX )
                              ENDIF

                           ELSEIF (ETYPE(J)(1:5) == 'QUAD8') THEN

                                                           ! Extrapolate stress to corners
                              CALL POLYNOM_FIT_STRE_STRN ( STRESS_RAW, 9, NUM_PTS(I), STRESS_OUT, STRESS_OUT_PCT_ERR,              &
                                STRESS_OUT_ERR_INDEX, PCT_ERR_MAX )

                                                           ! Transfrom stress to element coordinate system
                              DO M=2,NUM_PTS(I)
                                 CALL PLANE_COORD_TRANS_21( SHELL_STR_ANGLE( M ), TEL, '')
                                 CALL TRANSFORM_SHELL_STR( TEL, STRESS_OUT(:,M), ONE)
                              ENDDO

                                                           ! Center stress is the average of corner stress in element coordinates.
                                                           ! In MSC, the center force and moment resultants are the average but
                                                           ! this is equivalent.
                              STRESS_OUT(:,1) = (STRESS_OUT(:,2) + STRESS_OUT(:,3) + STRESS_OUT(:,4) + STRESS_OUT(:,5)) / FOUR

                           ENDIF

                        ENDIF

                        DO M=1,NUM_PTS(I)                  ! Calculate forces and moments from stresses
                           STRESS(:) = STRESS_OUT(:,M)
                           CALL SHELL_ENGR_FORCE_OGEL ( NUM_OGEL )

                           NUM_OGEL_ROWS = NUM_OGEL_ROWS + 1
                           EID_OUT_ARRAY(NUM_OGEL_ROWS,1) = EID
                           GID_OUT_ARRAY(NUM_OGEL_ROWS,1) = 0
                           DO K=1,ELGP
                              GID_OUT_ARRAY(NUM_OGEL_ROWS,K+1) = AGRID(K)
                           ENDDO

                        ENDDO


                        IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                           DO K=1,8
                              OT4_EROW = OT4_EROW + 1
                              OTM_ELFE(OT4_EROW,JVEC) = OGEL(NUM_OGEL,K)
                              IF (JVEC == 1) THEN
                                 WRITE(TXT_ELFE(OT4_EROW), 9192) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, FORCE_ITEM(K)
                              ENDIF
                           ENDDO
                        ENDIF

                        IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (JVEC == 1) .AND. (OT4_EROW >= 1)) THEN
                           DO K=1,OTMSKIP                     ! Write OTMSKIP blank separator lines
                              OT4_EROW = OT4_EROW + 1
                              WRITE(TXT_ELFE(OT4_EROW), 9199)
                           ENDDO
                        ENDIF

                        IF (NUM_OGEL_ROWS == NELREQ(I)) THEN
                           CALL CHK_OGEL_ZEROS ( NUM_OGEL )

 100                       FORMAT("*DEBUG:      ",A,"; ELEMENT_TYPE=",A,"; TABLE_NAME=",A,"; ITABLE=",I8)
                           WRITE(ERR,100) "F3",ETYPE(J),TABLE_NAME,ITABLE
                           CALL SET_OEF_TABLE_NAME(ETYPE(J), TABLE_NAME, ITABLE)
                           WRITE(ERR,100) "F4",ETYPE(J),TABLE_NAME,ITABLE
                           CALL WRITE_ELEM_ENGR_FORCE ( JVEC, NUM_OGEL_ROWS, IHDR, NUM_PTS(I), ITABLE )
                           CALL WRITE_SURFACE_AVERAGED_FORCES_F06 ( NUM_OGEL_ROWS, NUM_PTS(I), ETYPE(J) )
                           EXIT
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO elems_3
      ENDDO reqs3

 10   FORMAT("*DEBUG:      OEF_END 2D:    TABLE_NAME",A)
      WRITE(ERR,10) TABLE_NAME
      IF ((TABLE_NAME .NE. "OEF ERR ") .AND. (ITABLE < 0)) THEN
        CALL END_OP2_TABLE(ITABLE)
      ENDIF

      IF (WRITE_NEU .AND. (ANY_ELFE_OUTPUT > 0)) THEN

         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTRIA3K, 8, SUBR_NAME )
         DO J=1,NELE                                       ! Write out TRIA3K engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:6) == 'TRIA3K') THEN
               NUM_FROWS= NUM_FROWS+ 1
               OPT(4) = 'Y'
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               OPT(4) = 'N'
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )
               FEMAP_EL_VECS(NUM_FROWS,1) = FCONV(1)*STRESS(1)       ! X  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,2) = FCONV(1)*STRESS(2)       ! Y  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,3) = FCONV(1)*STRESS(3)       ! XY Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,4) = FCONV(2)*STRESS(4)       ! X  Moment
               FEMAP_EL_VECS(NUM_FROWS,5) = FCONV(2)*STRESS(5)       ! Y  Moment
               FEMAP_EL_VECS(NUM_FROWS,6) = FCONV(2)*STRESS(6)       ! XY Moment
               FEMAP_EL_VECS(NUM_FROWS,7) = FCONV(3)*STRESS(7)       ! X  Transverse Shear
               FEMAP_EL_VECS(NUM_FROWS,8) = FCONV(3)*STRESS(8)       ! Y  Transverse Shear
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'TRIA3K  ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTRIA3, 8, SUBR_NAME )
         DO J=1,NELE                                       ! Write out TRIA3 engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:6) == 'TRIA3 ') THEN
               NUM_FROWS= NUM_FROWS+ 1
               OPT(4) = 'Y'
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               OPT(4) = 'N'
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )
               FEMAP_EL_VECS(NUM_FROWS,1) = FCONV(1)*STRESS(1)       ! X  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,2) = FCONV(1)*STRESS(2)       ! Y  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,3) = FCONV(1)*STRESS(3)       ! XY Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,4) = FCONV(2)*STRESS(4)       ! X  Moment
               FEMAP_EL_VECS(NUM_FROWS,5) = FCONV(2)*STRESS(5)       ! Y  Moment
               FEMAP_EL_VECS(NUM_FROWS,6) = FCONV(2)*STRESS(6)       ! XY Moment
               FEMAP_EL_VECS(NUM_FROWS,7) = FCONV(3)*STRESS(7)       ! X  Transverse Shear
               FEMAP_EL_VECS(NUM_FROWS,8) = FCONV(3)*STRESS(8)       ! Y  Transverse Shear
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'TRIA3   ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCQUAD4K, 8, SUBR_NAME )
         DO J=1,NELE                                       ! Write out QUAD4K engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:6) == 'QUAD4K') THEN
               NUM_FROWS= NUM_FROWS+ 1
               OPT(4) = 'Y'
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               OPT(4) = 'N'
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )
               FEMAP_EL_VECS(NUM_FROWS,1) = FCONV(1)*STRESS(1)       ! X  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,2) = FCONV(1)*STRESS(2)       ! Y  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,3) = FCONV(1)*STRESS(3)       ! XY Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,4) = FCONV(2)*STRESS(4)       ! X  Moment
               FEMAP_EL_VECS(NUM_FROWS,5) = FCONV(2)*STRESS(5)       ! Y  Moment
               FEMAP_EL_VECS(NUM_FROWS,6) = FCONV(2)*STRESS(6)       ! XY Moment
               FEMAP_EL_VECS(NUM_FROWS,7) = FCONV(3)*STRESS(7)       ! X  Transverse Shear
               FEMAP_EL_VECS(NUM_FROWS,8) = FCONV(3)*STRESS(8)       ! Y  Transverse Shear
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'QUAD4K  ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCQUAD4 + NCQUADR, 8, SUBR_NAME )
         DO J=1,NELE                                       ! Write out QUAD4/CQUADR engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF ((ETYPE(J)(1:6) == 'QUAD4 ') .OR. (ETYPE(J) == 'QUADR   ')) THEN
               NUM_FROWS= NUM_FROWS+ 1
               OPT(4) = 'Y'
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               OPT(4) = 'N'
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )
               FEMAP_EL_VECS(NUM_FROWS,1) = FCONV(1)*STRESS(1)       ! X  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,2) = FCONV(1)*STRESS(2)       ! Y  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,3) = FCONV(1)*STRESS(3)       ! XY Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,4) = FCONV(2)*STRESS(4)       ! X  Moment
               FEMAP_EL_VECS(NUM_FROWS,5) = FCONV(2)*STRESS(5)       ! Y  Moment
               FEMAP_EL_VECS(NUM_FROWS,6) = FCONV(2)*STRESS(6)       ! XY Moment
               FEMAP_EL_VECS(NUM_FROWS,7) = FCONV(3)*STRESS(7)       ! X  Transverse Shear
               FEMAP_EL_VECS(NUM_FROWS,8) = FCONV(3)*STRESS(8)       ! Y  Transverse Shear
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'QUADR   ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NUM_FROWS= 0
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCSHEAR, 8, SUBR_NAME )
         DO J=1,NELE                                       ! Write out SHEAR engineering forces
            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J)(1:6) == 'SHEAR') THEN
               NUM_FROWS= NUM_FROWS+ 1
               OPT(4) = 'Y'
               DO K=0,MBUG-1
                  WRT_BUG(K) = 0
               ENDDO
               PLY_NUM = 0
               CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' ) ! 'N' in call to EMG means do not write to BUG file
               OPT(4) = 'N'
               FEMAP_EL_NUMS(NUM_FROWS,1) = EID
               IF (NUM_EMG_FATAL_ERRS > 0) THEN
                  IERROR = IERROR + 1
                  CYCLE
               ENDIF
               CALL ELMDIS
               CALL ELEM_STRE_STRN_ARRAYS ( 1 )
               FEMAP_EL_VECS(NUM_FROWS,1) = FCONV(1)*STRESS(1)       ! X  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,2) = FCONV(1)*STRESS(2)       ! Y  Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,3) = FCONV(1)*STRESS(3)       ! XY Membrane Force
               FEMAP_EL_VECS(NUM_FROWS,4) = FCONV(2)*STRESS(4)       ! X  Moment
               FEMAP_EL_VECS(NUM_FROWS,5) = FCONV(2)*STRESS(5)       ! Y  Moment
               FEMAP_EL_VECS(NUM_FROWS,6) = FCONV(2)*STRESS(6)       ! XY Moment
               FEMAP_EL_VECS(NUM_FROWS,7) = FCONV(3)*STRESS(7)       ! X  Transverse Shear
               FEMAP_EL_VECS(NUM_FROWS,8) = FCONV(3)*STRESS(8)       ! Y  Transverse Shear
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_ELFO_VECS ( 'SHEAR   ', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

      ENDIF

      IF (IERROR > 0) THEN
         REQUEST = 'ELEMENT ENGINEERING FORCE'
         WRITE(ERR,9201) TYPE, REQUEST, EID
         WRITE(F06,9201) TYPE, REQUEST, EID
      ENDIF



      RETURN

! **********************************************************************************************************************************
 9192 FORMAT(I8,1X,A,A8,I8,4X,A20)

 9199 FORMAT(' ')

 9201 FORMAT(' *ERROR  9201: DUE TO ABOVE LISTED ERRORS, CANNOT CALCULATE ',A,' REQUESTS FOR ',A,' ELEMENT ID = ',I8)

 9204 FORMAT(' *WARNING    : ELEMENT ENGINEERING FORCE OUTPUT REQUESTS FOR ',I8,1X,A,' NOT ALLOWED'                                &
                    ,/,14X,' REQUEST STRESS OUTPUT IN CASE CONTROL FOR THESE ELEMENTS')


! **********************************************************************************************************************************

      CONTAINS

!----------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE WRITE_SURFACE_AVERAGED_FORCES_F06 ( NUM, NUM_PTS_PER_ELEM, FAMILY )

      INTEGER(LONG), INTENT(IN)       :: NUM
      INTEGER(LONG), INTENT(IN)       :: NUM_PTS_PER_ELEM
      CHARACTER(LEN=*), INTENT(IN)    :: FAMILY

      INTEGER(LONG)                   :: IERR
      INTEGER(LONG)                   :: SURF
      INTEGER(LONG)                   :: NUM_ELEMS
      INTEGER(LONG)                   :: NUM_GRIDS
      INTEGER(LONG)                   :: I, K, GPOS, ELEM_POS
      INTEGER(LONG)                   :: ISTART
      INTEGER(LONG)                   :: NELGP
      INTEGER(LONG)                   :: GID
      INTEGER(LONG)                   :: POINT_ROW
      INTEGER(LONG), ALLOCATABLE      :: PATCH_ELEMS(:)
      INTEGER(LONG), ALLOCATABLE      :: PATCH_GRIDS(:)
      INTEGER(LONG), ALLOCATABLE      :: OUT_EIDS(:)
      REAL(DOUBLE), ALLOCATABLE       :: SUM_FORCE(:,:)
      REAL(DOUBLE), ALLOCATABLE       :: SUM_WT(:)
      REAL(DOUBLE)                    :: WT
      REAL(DOUBLE)                    :: FORCE_LOCAL(8)
      REAL(DOUBLE)                    :: FORCE_SURF(8)

      IF (.NOT. GPSTRESS_REQ) RETURN
      IF (NUM_GP_SURFACE <= 0) RETURN
      IF (NUM <= 0) RETURN
      IF (.NOT. ((FAMILY(1:4) == 'TRIA') .OR. (FAMILY(1:4) == 'QUAD'))) RETURN

      ALLOCATE(PATCH_ELEMS(MAX(1_LONG,NUM)))
      ALLOCATE(PATCH_GRIDS(MAX(1_LONG,4*NUM)))

      DO SURF=1,NUM_GP_SURFACE
         CALL GPSTRESS_COLLECT_SURFACE_PATCH ( SURF, PATCH_ELEMS, NUM_ELEMS, PATCH_GRIDS, NUM_GRIDS, IERR )
         IF ((IERR /= 0) .OR. (NUM_ELEMS <= 0) .OR. (NUM_GRIDS <= 0)) CYCLE

         ALLOCATE(SUM_FORCE(8,NUM_GRIDS))
         ALLOCATE(SUM_WT(NUM_GRIDS))
         ALLOCATE(OUT_EIDS(NUM_GRIDS))
         SUM_FORCE = ZERO
         SUM_WT = ZERO
         OUT_EIDS = 0

         DO ISTART=1,NUM,MAX(1_LONG,NUM_PTS_PER_ELEM)
            IF (FIND_INT(EID_OUT_ARRAY(ISTART,1), PATCH_ELEMS, NUM_ELEMS) == 0) CYCLE

            NELGP = GET_ELEM_NUM_CORNERS ( ISTART, FAMILY )
            IF (NELGP <= 0) CYCLE
            WT = GET_ELEM_AREA ( ISTART, NELGP ) / REAL(NELGP,DOUBLE)
            IF (WT <= ZERO) WT = ONE / REAL(NELGP,DOUBLE)

            DO K=1,NELGP
               GID = GID_OUT_ARRAY(ISTART,K+1)
               GPOS = FIND_INT(GID, PATCH_GRIDS, NUM_GRIDS)
               IF (GPOS == 0) CYCLE

               IF (NUM_PTS_PER_ELEM > 1) THEN
                  POINT_ROW = ISTART + K
                  IF (POINT_ROW > NUM) POINT_ROW = ISTART
               ELSE
                  POINT_ROW = ISTART
               ENDIF

               FORCE_LOCAL(1:8) = OGEL(POINT_ROW,1:8)
               CALL TRANSFORM_SURFACE_FORCE8 ( SURF, POINT_ROW, FORCE_LOCAL, FORCE_SURF )
               SUM_FORCE(1:8,GPOS) = SUM_FORCE(1:8,GPOS) + WT * FORCE_SURF(1:8)
               SUM_WT(GPOS) = SUM_WT(GPOS) + WT
               IF (OUT_EIDS(GPOS) == 0) OUT_EIDS(GPOS) = EID_OUT_ARRAY(ISTART,1)
            ENDDO
         ENDDO

         WRITE(F06,'(//,33X,A,I8)') 'F O R C E S   A T   G R I D   P O I N T S   - -   S U R F A C E', GP_SURFACE_IDS(SURF)
         WRITE(F06,'(A,22X,A,A1,8X,A)') '0', 'SURFACE X-AXIS X  NORMAL(Z-AXIS)  ', GP_SURFACE_NORMAL_MODE(SURF)(1:1),              &
                                        'REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID        0'
         WRITE(F06,'(/,8X,A,7X,A,13X,A,11X,A,11X,A,11X,A,11X,A,11X,A,12X,A,12X,A)')                                      &
                    'GRID', 'ELEM', 'NXX', 'NYY', 'NXY', 'MXX', 'MYY', 'MXY', 'QX', 'QY'
         DO I=1,NUM_GRIDS
            IF (SUM_WT(I) > ZERO) THEN
               FORCE_SURF(1:8) = SUM_FORCE(1:8,I) / SUM_WT(I)
               CALL WRITE_SURFACE_FORCE_ROW ( PATCH_GRIDS(I), OUT_EIDS(I), FORCE_SURF )
            ENDIF
         ENDDO
         WRITE(F06,*)

         DEALLOCATE(SUM_FORCE)
         DEALLOCATE(SUM_WT)
         DEALLOCATE(OUT_EIDS)
      ENDDO

      DEALLOCATE(PATCH_ELEMS)
      DEALLOCATE(PATCH_GRIDS)

      END SUBROUTINE WRITE_SURFACE_AVERAGED_FORCES_F06

!----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG) FUNCTION GET_ELEM_NUM_CORNERS ( ISTART, FAMILY )

      INTEGER(LONG), INTENT(IN)       :: ISTART
      CHARACTER(LEN=*), INTENT(IN)    :: FAMILY

      INTEGER(LONG)                   :: K

      GET_ELEM_NUM_CORNERS = 0
      IF (FAMILY(1:4) == 'TRIA') THEN
         GET_ELEM_NUM_CORNERS = 3
      ELSE IF (FAMILY(1:4) == 'QUAD') THEN
         GET_ELEM_NUM_CORNERS = 4
      ENDIF

      DO K=GET_ELEM_NUM_CORNERS,1,-1
         IF (GID_OUT_ARRAY(ISTART,K+1) > 0) RETURN
      ENDDO
      GET_ELEM_NUM_CORNERS = 0

      END FUNCTION GET_ELEM_NUM_CORNERS

!----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE) FUNCTION GET_ELEM_AREA ( ISTART, NELGP )

      INTEGER(LONG), INTENT(IN)       :: ISTART
      INTEGER(LONG), INTENT(IN)       :: NELGP

      REAL(DOUBLE)                    :: X1(3), X2(3), X3(3), X4(3)

      IF (NELGP == 3) THEN
         CALL GET_GRID_BASIC_COORDS ( GID_OUT_ARRAY(ISTART,2), X1 )
         CALL GET_GRID_BASIC_COORDS ( GID_OUT_ARRAY(ISTART,3), X2 )
         CALL GET_GRID_BASIC_COORDS ( GID_OUT_ARRAY(ISTART,4), X3 )
         GET_ELEM_AREA = TRI_AREA_FROM_XYZ ( X1, X2, X3 )
      ELSE IF (NELGP == 4) THEN
         CALL GET_GRID_BASIC_COORDS ( GID_OUT_ARRAY(ISTART,2), X1 )
         CALL GET_GRID_BASIC_COORDS ( GID_OUT_ARRAY(ISTART,3), X2 )
         CALL GET_GRID_BASIC_COORDS ( GID_OUT_ARRAY(ISTART,4), X3 )
         CALL GET_GRID_BASIC_COORDS ( GID_OUT_ARRAY(ISTART,5), X4 )
         GET_ELEM_AREA = TRI_AREA_FROM_XYZ ( X1, X2, X3 ) + TRI_AREA_FROM_XYZ ( X1, X3, X4 )
      ELSE
         GET_ELEM_AREA = ZERO
      ENDIF

      END FUNCTION GET_ELEM_AREA

!----------------------------------------------------------------------------------------------------------------------------------
      REAL(DOUBLE) FUNCTION TRI_AREA_FROM_XYZ ( X1, X2, X3 )

      REAL(DOUBLE), INTENT(IN)        :: X1(3), X2(3), X3(3)
      REAL(DOUBLE)                    :: C(3)

      C(1) = (X2(2)-X1(2))*(X3(3)-X1(3)) - (X2(3)-X1(3))*(X3(2)-X1(2))
      C(2) = (X2(3)-X1(3))*(X3(1)-X1(1)) - (X2(1)-X1(1))*(X3(3)-X1(3))
      C(3) = (X2(1)-X1(1))*(X3(2)-X1(2)) - (X2(2)-X1(2))*(X3(1)-X1(1))
      TRI_AREA_FROM_XYZ = 0.5D0 * DSQRT(C(1)*C(1) + C(2)*C(2) + C(3)*C(3))

      END FUNCTION TRI_AREA_FROM_XYZ

!----------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE GET_GRID_BASIC_COORDS ( GRID_NUM, XYZ )

      INTEGER(LONG), INTENT(IN)       :: GRID_NUM
      REAL(DOUBLE), INTENT(OUT)       :: XYZ(3)

      CHARACTER(32*BYTE)              :: LOCAL_SUBR_NAME = 'OFP3_ELFE_2D'
      INTEGER(LONG)                   :: IGRID

      CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', LOCAL_SUBR_NAME, SIZE(GRID_ID), GRID_ID, GRID_NUM, IGRID )
      IF (IGRID > 0) THEN
         XYZ(1) = RGRID(IGRID,1)
         XYZ(2) = RGRID(IGRID,2)
         XYZ(3) = RGRID(IGRID,3)
      ELSE
         XYZ = ZERO
      ENDIF

      END SUBROUTINE GET_GRID_BASIC_COORDS

!----------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE TRANSFORM_SURFACE_FORCE8 ( SURF_INDEX, POINT_INDEX, FORCE_LOCAL, FORCE_SURF )

      INTEGER(LONG), INTENT(IN)       :: SURF_INDEX
      INTEGER(LONG), INTENT(IN)       :: POINT_INDEX
      REAL(DOUBLE), INTENT(IN)        :: FORCE_LOCAL(8)
      REAL(DOUBLE), INTENT(OUT)       :: FORCE_SURF(8)

      REAL(DOUBLE)                    :: SURF_BASIS(3,3)
      REAL(DOUBLE)                    :: LOCAL_TENSOR(3,3)
      REAL(DOUBLE)                    :: SURF_TENSOR(3,3)
      REAL(DOUBLE)                    :: LOCAL_VEC(3)
      REAL(DOUBLE)                    :: SURF_VEC(3)

      FORCE_SURF = FORCE_LOCAL

      CALL GET_SURFACE_BASIS ( SURF_INDEX, SURF_BASIS )

      LOCAL_TENSOR = ZERO
      LOCAL_TENSOR(1,1) = FORCE_LOCAL(1)
      LOCAL_TENSOR(2,2) = FORCE_LOCAL(2)
      LOCAL_TENSOR(1,2) = FORCE_LOCAL(3)
      LOCAL_TENSOR(2,1) = FORCE_LOCAL(3)
      SURF_TENSOR = MATMUL(SURF_BASIS, MATMUL(LOCAL_TENSOR, TRANSPOSE(SURF_BASIS)))
      FORCE_SURF(1) = SURF_TENSOR(1,1)
      FORCE_SURF(2) = SURF_TENSOR(2,2)
      FORCE_SURF(3) = SURF_TENSOR(1,2)

      LOCAL_TENSOR = ZERO
      LOCAL_TENSOR(1,1) = FORCE_LOCAL(4)
      LOCAL_TENSOR(2,2) = FORCE_LOCAL(5)
      LOCAL_TENSOR(1,2) = FORCE_LOCAL(6)
      LOCAL_TENSOR(2,1) = FORCE_LOCAL(6)
      SURF_TENSOR = MATMUL(SURF_BASIS, MATMUL(LOCAL_TENSOR, TRANSPOSE(SURF_BASIS)))
      FORCE_SURF(4) = SURF_TENSOR(1,1)
      FORCE_SURF(5) = SURF_TENSOR(2,2)
      FORCE_SURF(6) = SURF_TENSOR(1,2)

      LOCAL_VEC = ZERO
      LOCAL_VEC(1) = FORCE_LOCAL(7)
      LOCAL_VEC(2) = FORCE_LOCAL(8)
      SURF_VEC = MATMUL(SURF_BASIS, LOCAL_VEC)
      FORCE_SURF(7) = SURF_VEC(1)
      FORCE_SURF(8) = SURF_VEC(2)

      END SUBROUTINE TRANSFORM_SURFACE_FORCE8

!----------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE GET_SURFACE_BASIS ( SURF_INDEX, SURF_BASIS )

      INTEGER(LONG), INTENT(IN)       :: SURF_INDEX
      REAL(DOUBLE), INTENT(OUT)       :: SURF_BASIS(3,3)

      CHARACTER(8*BYTE)               :: NMODE

      SURF_BASIS = ZERO
      NMODE = GP_SURFACE_NORMAL_MODE(SURF_INDEX)

      IF (NMODE(1:1) == 'X') THEN
         SURF_BASIS(1,2) = ONE
         SURF_BASIS(2,3) = ONE
         SURF_BASIS(3,1) = ONE
      ELSE IF (NMODE(1:1) == 'Y') THEN
         SURF_BASIS(1,3) = ONE
         SURF_BASIS(2,1) = ONE
         SURF_BASIS(3,2) = ONE
      ELSE
         SURF_BASIS(1,1) = ONE
         SURF_BASIS(2,2) = ONE
         SURF_BASIS(3,3) = ONE
      ENDIF

      END SUBROUTINE GET_SURFACE_BASIS

!----------------------------------------------------------------------------------------------------------------------------------
      INTEGER(LONG) FUNCTION FIND_INT ( VALUE, ARRAY, NUSED )

      INTEGER(LONG), INTENT(IN)       :: VALUE
      INTEGER(LONG), INTENT(IN)       :: ARRAY(:)
      INTEGER(LONG), INTENT(IN)       :: NUSED

      INTEGER(LONG)                   :: I

      FIND_INT = 0
      DO I=1,NUSED
         IF (ARRAY(I) == VALUE) THEN
            FIND_INT = I
            RETURN
         ENDIF
      ENDDO

      END FUNCTION FIND_INT

!----------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE WRITE_SURFACE_FORCE_ROW ( GRID_NUM, ELEM_NUM, VALUES )

      INTEGER(LONG), INTENT(IN)       :: GRID_NUM
      INTEGER(LONG), INTENT(IN)       :: ELEM_NUM
      REAL(DOUBLE), INTENT(IN)        :: VALUES(8)

      CHARACTER(160*BYTE)             :: LINE_BUF
      INTEGER(LONG)                   :: POS, J

      LINE_BUF = ' '
      CALL FAST_FMT_I8_RJ ( GRID_NUM, LINE_BUF(2:9) )
      CALL FAST_FMT_I8_RJ ( ELEM_NUM, LINE_BUF(12:19) )
      POS = 25
      DO J=1,8
         CALL FAST_FMT_F06_E14_6 ( VALUES(J), LINE_BUF(POS:POS+13) )
         POS = POS + 14
      ENDDO
      WRITE(F06,'(A)') TRIM(LINE_BUF)

      END SUBROUTINE WRITE_SURFACE_FORCE_ROW

      END SUBROUTINE OFP3_ELFE_2D
