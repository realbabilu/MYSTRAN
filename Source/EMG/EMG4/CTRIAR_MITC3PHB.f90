! ##################################################################################################################################
! CTRIAR MITC3+HB selector for PARAM,TRIARTYP,MITC3+HB.

      SUBROUTINE CTRIAR_MITC3PHB ( OPT, INT_ELEM_ID )

! Port target: D:\18a\bending_only\Shell\gemini2\shit\MITC3pD_HughesBrezzi_ShellElement.py
! MITC3+ triangular shell with two internal bubble bending/shear DOF statically
! condensed, plus Hughes-Brezzi drilling coupling.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NSUB, SOL_NAME
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, THREE, TWELVE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  COUPMASS
      USE MODEL_STUF, ONLY            :  EID, ELGP, KE, KED, ME, BE1, BE2, BE3, EPROP, MASS_PER_UNIT_AREA, PRESS, PPE,             &
                                         TE, NUM_EMG_FATAL_ERRS, SHELL_A, SHELL_D, SHELL_T, XEB
      USE CTRIAR_DKMT18_Interface
      USE OUTA_HERE_Interface

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CTRIAR_MITC3PHB'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID
      CHARACTER(1*BYTE)                :: REC_OPT(6)

      INTEGER(LONG), PARAMETER        :: NNODE = 3
      INTEGER(LONG), PARAMETER        :: NDOF  = 18
      INTEGER(LONG), PARAMETER        :: NALL  = 20
      REAL(DOUBLE), PARAMETER         :: KAPPA = 5.0D0/6.0D0

      INTEGER(LONG)                   :: I, J, K, GP, JSUB
      REAL(DOUBLE)                    :: XYZ(NNODE,3), EG(3,3), XY(NNODE,2), DNX(NNODE), DNY(NNODE)
      REAL(DOUBLE)                    :: AREA, THICK, GVAL, CDRILL, MASS_NODE
      REAL(DOUBLE)                    :: KFULL(NALL,NALL), KCOND(NDOF,NDOF), KBB(2,2), KBBI(2,2), KBC(2,NDOF), BMAP(2,NDOF)
      REAL(DOUBLE)                    :: BM(3,NALL), BB(3,NALL), BS(2,NALL), BD(1,NALL)
      REAL(DOUBLE)                    :: BM18(3,NDOF), BB18(3,NDOF), BS18(2,NDOF), T18(NDOF,NDOF)
      REAL(DOUBLE)                    :: RGP(3), SGP(3), WT, UNIT_PPE(NDOF)

      IF (ELGP /= 3) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9001) SUBR_NAME, EID, ELGP
         WRITE(F06,9001) SUBR_NAME, EID, ELGP
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      THICK = EPROP(1)
      IF (THICK <= ZERO) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9002) SUBR_NAME, EID, THICK
         WRITE(F06,9002) SUBR_NAME, EID, THICK
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL LOAD_BASIC_COORDS ( XYZ )
      CALL MITC3PHB_GEOMETRY ( XYZ, EG, XY, DNX, DNY, AREA )
      CALL BUILD_T18 ( EG, T18 )

      TE = ZERO
      TE(1,1) = ONE
      TE(2,2) = ONE
      TE(3,3) = ONE

      RGP = (/ ONE/6.0D0, TWO/THREE, ONE/6.0D0 /)
      SGP = (/ ONE/6.0D0, ONE/6.0D0, TWO/THREE /)
      WT  = ONE/THREE

      IF ((OPT(3) == 'Y') .OR. (OPT(4) == 'Y')) THEN
         KFULL = ZERO
         CALL MITC3PHB_BM ( DNX, DNY, BM )
         KFULL = KFULL + MATMUL(TRANSPOSE(BM), MATMUL(SHELL_A, BM)) * AREA

         GVAL = ZERO
         IF (DABS(SHELL_T(1,1)) > 1.0D-20) GVAL = SHELL_T(1,1)/(KAPPA*THICK)
         CDRILL = 1.0D-4 * GVAL * THICK

         DO GP=1,3
            CALL MITC3PHB_BCURV ( RGP(GP), SGP(GP), DNX, DNY, BB )
            CALL MITC3PHB_BSHEAR ( RGP(GP), SGP(GP), XY, AREA, BS )
            CALL MITC3PHB_BDRILL ( RGP(GP), SGP(GP), DNX, DNY, BD )
            KFULL = KFULL + MATMUL(TRANSPOSE(BB), MATMUL(SHELL_D, BB)) * (WT*AREA)
            KFULL = KFULL + MATMUL(TRANSPOSE(BS), MATMUL(SHELL_T, BS)) * (WT*AREA)
            KFULL = KFULL + MATMUL(TRANSPOSE(BD), CDRILL*BD) * (WT*THICK*AREA)
         ENDDO

         KBB = KFULL(19:20,19:20)
         KBC = KFULL(19:20,1:NDOF)
         CALL INV2 ( KBB, KBBI )
         KCOND = KFULL(1:NDOF,1:NDOF) - MATMUL(TRANSPOSE(KBC), MATMUL(KBBI, KBC))
         KCOND = 0.5D0*(KCOND + TRANSPOSE(KCOND))

         IF (OPT(4) == 'Y') THEN
            KE(1:NDOF,1:NDOF) = MATMUL(TRANSPOSE(T18), MATMUL(KCOND, T18))
         ENDIF

         IF (OPT(3) == 'Y') THEN
            REC_OPT = 'N'
            REC_OPT(3) = 'Y'
            CALL CTRIAR_DKMT18 ( REC_OPT, INT_ELEM_ID )
         ENDIF

         IF ((DEBUG(233) > 0) .AND. (OPT(4) == 'Y')) THEN
            WRITE(F06,'(A,I8,A,ES15.7)') 'CTRIAR_MITC3PHB EID=', EID, ' KE_NORM=', DSQRT(SUM(KE(1:NDOF,1:NDOF)*KE(1:NDOF,1:NDOF)))
         ENDIF
      ENDIF

      IF (OPT(1) == 'Y') THEN
         ME(1:NDOF,1:NDOF) = ZERO
         IF ((SOL_NAME(1:5) == 'MODES') .AND. (COUPMASS > 0)) THEN
            DO I=1,NNODE
               DO J=1,NNODE
                  IF (I == J) THEN
                     MASS_NODE = MASS_PER_UNIT_AREA * AREA / 6.0D0
                  ELSE
                     MASS_NODE = MASS_PER_UNIT_AREA * AREA / 12.0D0
                  ENDIF
                  DO K=1,3
                     ME(6*(I-1)+K,6*(J-1)+K) = MASS_NODE
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            MASS_NODE = MASS_PER_UNIT_AREA * AREA / THREE
            DO I=1,NNODE
               DO J=1,3
                  ME(6*(I-1)+J,6*(I-1)+J) = MASS_NODE
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      IF (OPT(5) == 'Y') THEN
         UNIT_PPE = ZERO
         DO I=1,NNODE
            UNIT_PPE(6*(I-1)+1) = AREA*EG(1,3)/THREE
            UNIT_PPE(6*(I-1)+2) = AREA*EG(2,3)/THREE
            UNIT_PPE(6*(I-1)+3) = AREA*EG(3,3)/THREE
         ENDDO
         DO JSUB=1,NSUB
            PPE(1:NDOF,JSUB) = PPE(1:NDOF,JSUB) + UNIT_PPE(1:NDOF) * PRESS(3,JSUB)
         ENDDO
      ENDIF

      IF (OPT(6) == 'Y') THEN
         KED(1:NDOF,1:NDOF) = ZERO
      ENDIF

      RETURN

 9001 FORMAT(' *ERROR: ',A,' expects ELGP=3 for element ',I8,' but got ',I8)
 9002 FORMAT(' *ERROR: ',A,' element ',I8,' has nonpositive thickness ',ES15.7)

      CONTAINS

      SUBROUTINE MITC3PHB_BUILD_RECOVERY ( BM18_IN, BB18_IN, BS18_IN, T18_IN )
      REAL(DOUBLE), INTENT(IN) :: BM18_IN(3,NDOF), BB18_IN(3,NDOF), BS18_IN(2,NDOF), T18_IN(NDOF,NDOF)

      BE1(1:3,1:NDOF,1) = MATMUL(BM18_IN, T18_IN)
      BE2(1:3,1:NDOF,1) = MATMUL(BB18_IN, T18_IN)
      BE3(1:2,1:NDOF,1) = MATMUL(BS18_IN, T18_IN)

      IF (DEBUG(233) > 0) THEN
         WRITE(F06,'(A,I8)') 'CTRIAR_MITC3PHB RECOVERY EID=', EID
      ENDIF
      END SUBROUTINE MITC3PHB_BUILD_RECOVERY

      SUBROUTINE LOAD_BASIC_COORDS ( XYZOUT )
      REAL(DOUBLE), INTENT(OUT) :: XYZOUT(NNODE,3)
      INTEGER(LONG) :: II, JJ
      DO II=1,NNODE
         DO JJ=1,3
            XYZOUT(II,JJ) = XEB(II,JJ)
         ENDDO
      ENDDO
      END SUBROUTINE LOAD_BASIC_COORDS

      SUBROUTINE MITC3PHB_GEOMETRY ( XYZN, EG, XY, DNX, DNY, AREA )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(NNODE,3)
      REAL(DOUBLE), INTENT(OUT) :: EG(3,3), XY(NNODE,2), DNX(NNODE), DNY(NNODE), AREA
      REAL(DOUBLE) :: V12(3), V13(3), NM, X2, X3, Y3
      V12 = XYZN(2,:) - XYZN(1,:)
      V13 = XYZN(3,:) - XYZN(1,:)
      NM = VNORM(V12)
      IF (NM <= 1.0D-14) CALL GEOM_FATAL
      EG(:,1) = V12/NM
      CALL CROSSV(V12, V13, EG(:,3))
      NM = VNORM(EG(:,3))
      IF (NM <= 1.0D-14) CALL GEOM_FATAL
      EG(:,3) = EG(:,3)/NM
      CALL CROSSV(EG(:,3), EG(:,1), EG(:,2))
      X2 = DOT_PRODUCT(V12, EG(:,1))
      X3 = DOT_PRODUCT(V13, EG(:,1))
      Y3 = DOT_PRODUCT(V13, EG(:,2))
      XY = ZERO
      XY(2,1) = X2
      XY(3,1) = X3
      XY(3,2) = Y3
      AREA = 0.5D0*X2*Y3
      IF (AREA <= 1.0D-14) CALL GEOM_FATAL
      DNX(2) = ONE/X2
      DNY(2) = -X3/(X2*Y3)
      DNX(3) = ZERO
      DNY(3) = ONE/Y3
      DNX(1) = -DNX(2) - DNX(3)
      DNY(1) = -DNY(2) - DNY(3)
      END SUBROUTINE MITC3PHB_GEOMETRY

      SUBROUTINE GEOM_FATAL
      NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
      FATAL_ERR = FATAL_ERR + 1
      WRITE(ERR,*) ' *ERROR: CTRIAR_MITC3PHB degenerate geometry on element ', EID
      WRITE(F06,*) ' *ERROR: CTRIAR_MITC3PHB degenerate geometry on element ', EID
      CALL OUTA_HERE ( 'Y' )
      END SUBROUTINE GEOM_FATAL

      SUBROUTINE MITC3PHB_BM ( DNX, DNY, BM )
      REAL(DOUBLE), INTENT(IN)  :: DNX(NNODE), DNY(NNODE)
      REAL(DOUBLE), INTENT(OUT) :: BM(3,NALL)
      INTEGER(LONG) :: II, C
      BM = ZERO
      DO II=1,NNODE
         C = (II-1)*6
         BM(1,C+1) = DNX(II)
         BM(2,C+2) = DNY(II)
         BM(3,C+1) = DNY(II)
         BM(3,C+2) = DNX(II)
      ENDDO
      END SUBROUTINE MITC3PHB_BM

      SUBROUTINE MITC3PHB_BCURV ( R1, R2, DNX, DNY, BB )
      REAL(DOUBLE), INTENT(IN)  :: R1, R2, DNX(NNODE), DNY(NNODE)
      REAL(DOUBLE), INTENT(OUT) :: BB(3,NALL)
      INTEGER(LONG) :: II, C
      REAL(DOUBLE) :: F4, DF1, DF2, DFX, DFY
      BB = ZERO
      DO II=1,NNODE
         C = (II-1)*6
         BB(1,C+5) =  DNX(II)
         BB(2,C+4) = -DNY(II)
         BB(3,C+4) = -DNX(II)
         BB(3,C+5) =  DNY(II)
      ENDDO
      CALL BUBBLE ( R1, R2, F4, DF1, DF2 )
      DFX = DF1*DNX(2)
      DFY = DF1*DNY(2) + DF2*DNY(3)
      BB(1,20) = DFX
      BB(2,19) = -DFY
      BB(3,19) = -DFX
      BB(3,20) = DFY
      END SUBROUTINE MITC3PHB_BCURV

      SUBROUTINE MITC3PHB_BSHEAR ( R1, R2, XY, AREA, BS )
      REAL(DOUBLE), INTENT(IN)  :: R1, R2, XY(NNODE,2), AREA
      REAL(DOUBLE), INTENT(OUT) :: BS(2,NALL)
      REAL(DOUBLE) :: BMID(3,NALL), NG(2,3), JMAT(2,2), TMP(2,NALL), F4, DF1, DF2
      REAL(DOUBLE) :: X21, Y21, X31, Y31, X32, Y32
      BMID = ZERO
      X21 = XY(2,1) - XY(1,1); Y21 = XY(2,2) - XY(1,2)
      X31 = XY(3,1) - XY(1,1); Y31 = XY(3,2) - XY(1,2)
      X32 = XY(3,1) - XY(2,1); Y32 = XY(3,2) - XY(2,2)
      BMID(1,3)  = -ONE; BMID(1,4)  = -Y21/TWO; BMID(1,5)  = X21/TWO
      BMID(1,9)  =  ONE; BMID(1,10) = -Y21/TWO; BMID(1,11) = X21/TWO
      BMID(2,3)  = -ONE; BMID(2,4)  = -Y31/TWO; BMID(2,5)  = X31/TWO
      BMID(2,15) =  ONE; BMID(2,16) = -Y31/TWO; BMID(2,17) = X31/TWO
      BMID(3,9)  = -ONE; BMID(3,10) = -Y32/TWO; BMID(3,11) = X32/TWO
      BMID(3,15) =  ONE; BMID(3,16) = -Y32/TWO; BMID(3,17) = X32/TWO
      JMAT(1,1) = Y31/(TWO*AREA)
      JMAT(1,2) = -Y21/(TWO*AREA)
      JMAT(2,1) = -X31/(TWO*AREA)
      JMAT(2,2) = X21/(TWO*AREA)
      NG(1,:) = (/ ONE - R2, -R2, R2 /)
      NG(2,:) = (/ -R1, ONE - R1, R1 /)
      TMP = MATMUL(NG, BMID)
      BS = MATMUL(JMAT, TMP)
      CALL BUBBLE ( R1, R2, F4, DF1, DF2 )
      BS(1,20) = BS(1,20) - F4
      BS(2,19) = BS(2,19) + F4
      END SUBROUTINE MITC3PHB_BSHEAR

      SUBROUTINE MITC3PHB_BDRILL ( R1, R2, DNX, DNY, BD )
      REAL(DOUBLE), INTENT(IN)  :: R1, R2, DNX(NNODE), DNY(NNODE)
      REAL(DOUBLE), INTENT(OUT) :: BD(1,NALL)
      REAL(DOUBLE) :: NVAL(NNODE)
      INTEGER(LONG) :: II, C
      BD = ZERO
      NVAL = (/ ONE - R1 - R2, R1, R2 /)
      DO II=1,NNODE
         C = (II-1)*6
         BD(1,C+1) = -0.5D0*DNY(II)
         BD(1,C+2) =  0.5D0*DNX(II)
         BD(1,C+6) = -NVAL(II)
      ENDDO
      END SUBROUTINE MITC3PHB_BDRILL

      SUBROUTINE BUBBLE ( R1, R2, F4, DF1, DF2 )
      REAL(DOUBLE), INTENT(IN)  :: R1, R2
      REAL(DOUBLE), INTENT(OUT) :: F4, DF1, DF2
      REAL(DOUBLE) :: LAM
      LAM = ONE - R1 - R2
      F4 = 27.0D0*R1*R2*LAM
      DF1 = 27.0D0*R2*(LAM - R1)
      DF2 = 27.0D0*R1*(LAM - R2)
      END SUBROUTINE BUBBLE

      SUBROUTINE INV2 ( A, AINV )
      REAL(DOUBLE), INTENT(IN)  :: A(2,2)
      REAL(DOUBLE), INTENT(OUT) :: AINV(2,2)
      REAL(DOUBLE) :: DET
      DET = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      IF (DABS(DET) <= 1.0D-20) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,*) ' *ERROR: CTRIAR_MITC3PHB singular bubble block on element ', EID
         WRITE(F06,*) ' *ERROR: CTRIAR_MITC3PHB singular bubble block on element ', EID
         CALL OUTA_HERE ( 'Y' )
      ENDIF
      AINV(1,1) =  A(2,2)/DET
      AINV(1,2) = -A(1,2)/DET
      AINV(2,1) = -A(2,1)/DET
      AINV(2,2) =  A(1,1)/DET
      END SUBROUTINE INV2

      SUBROUTINE BUILD_T18 ( EG, T18 )
      REAL(DOUBLE), INTENT(IN)  :: EG(3,3)
      REAL(DOUBLE), INTENT(OUT) :: T18(NDOF,NDOF)
      INTEGER(LONG) :: BI, R, C, OFF
      T18 = ZERO
      DO BI=1,NNODE
         OFF = 6*(BI-1)
         DO R=1,3
            DO C=1,3
               T18(OFF+R,OFF+C) = EG(C,R)
               T18(OFF+3+R,OFF+3+C) = EG(C,R)
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE BUILD_T18

      SUBROUTINE CROSSV(A,B,C)
      REAL(DOUBLE), INTENT(IN) :: A(3), B(3)
      REAL(DOUBLE), INTENT(OUT) :: C(3)
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      END SUBROUTINE CROSSV

      FUNCTION VNORM ( V ) RESULT(NM)
      REAL(DOUBLE), INTENT(IN) :: V(3)
      REAL(DOUBLE) :: NM
      NM = DSQRT(DOT_PRODUCT(V,V))
      END FUNCTION VNORM

      END SUBROUTINE CTRIAR_MITC3PHB
