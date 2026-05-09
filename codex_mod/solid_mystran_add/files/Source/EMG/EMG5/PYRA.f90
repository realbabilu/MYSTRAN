      SUBROUTINE PYRA ( OPT, INT_ELEM_ID, RED_INT_SHEAR, WRITE_WARN )

! --- newsolid_add begin --- !
! CPYRA solid element kernel.
!
! PYRA5  : linear five-node pyramid, 2x2x2 integration, matching
!          solid3d_cpyra5.py.
! PYRA14 : Liu symmetric composite 14-node quadratic pyramid trial, matching
!          solid3d_cpyra14_liu.py including finite-difference reference
!          derivatives.
! --- newsolid_add end --- !

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE CONSTANTS_1, ONLY           :  HALF, ZERO, ONE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE SCONTR, ONLY                :  NTSUB
      USE MODEL_STUF, ONLY            :  ALPVEC, BE1, BE2, DT, EID, ELGP, NUM_EMG_FATAL_ERRS, ES, KE, KED, ME, PTE, RHO, SE1, SE2, &
                                         STE1, STRESS, TREF, TYPE, XEL
      USE PARAMS, ONLY                :  EPSIL, SOLIDTYP

      USE PYRA_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'PYRA'
      CHARACTER( 1*BYTE), INTENT(IN)  :: RED_INT_SHEAR
      CHARACTER( 1*BYTE), INTENT(IN)  :: OPT(6)
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID
      INTEGER(LONG)                   :: I,J,K,L
      INTEGER(LONG)                   :: ID(3*ELGP)
      INTEGER(LONG)                   :: II,JJ
      INTEGER(LONG)                   :: NGP
      INTEGER(LONG), PARAMETER        :: EAS54_N = 54

      REAL(DOUBLE)                    :: B(6,3*ELGP)
      REAL(DOUBLE)                    :: CBAR(3,3*ELGP)
      REAL(DOUBLE)                    :: DETJ
      REAL(DOUBLE)                    :: DNX(3,ELGP)
      REAL(DOUBLE)                    :: DUM0(3*ELGP)
      REAL(DOUBLE)                    :: DUM3(3*ELGP,3*ELGP)
      REAL(DOUBLE)                    :: DUM4(6,3*ELGP)
      REAL(DOUBLE)                    :: DUM5(3*ELGP,3*ELGP)
      REAL(DOUBLE)                    :: DUM6(3,3*ELGP)
      REAL(DOUBLE)                    :: EALP(6)
      REAL(DOUBLE)                    :: EAS_DM(6,EAS54_N)
      REAL(DOUBLE)                    :: EAS_KAA(EAS54_N,EAS54_N)
      REAL(DOUBLE)                    :: EAS_KUA(3*ELGP,EAS54_N)
      REAL(DOUBLE)                    :: EAS_M(6,EAS54_N)
      REAL(DOUBLE)                    :: EAS_X(EAS54_N,3*ELGP)
      REAL(DOUBLE)                    :: EPS1
      REAL(DOUBLE)                    :: INT_R(64)
      REAL(DOUBLE)                    :: INT_S(64)
      REAL(DOUBLE)                    :: INT_T(64)
      REAL(DOUBLE)                    :: INT_W(64)
      REAL(DOUBLE)                    :: M_1DOF(ELGP,ELGP)
      REAL(DOUBLE)                    :: N(ELGP)
      REAL(DOUBLE)                    :: KWW(3,3)
      REAL(DOUBLE)                    :: RHO1
      REAL(DOUBLE)                    :: ALP(6)
      REAL(DOUBLE)                    :: SIGxx
      REAL(DOUBLE)                    :: SIGxy
      REAL(DOUBLE)                    :: SIGyy
      REAL(DOUBLE)                    :: SIGyz
      REAL(DOUBLE)                    :: SIGzx
      REAL(DOUBLE)                    :: SIGzz
      REAL(DOUBLE)                    :: TBAR(NTSUB)
      REAL(DOUBLE)                    :: TEMP
      REAL(DOUBLE)                    :: TREF1
      REAL(DOUBLE)                    :: VOLUME

      IF ((RED_INT_SHEAR == 'Y') .AND. (WRITE_WARN == 'Y') .AND. (INT_ELEM_ID < 0)) THEN
         RETURN
      ENDIF

      EPS1 = EPSIL(1)
      TREF1 = TREF(1)

      DO I=1,6
         ALP(I) = ALPVEC(I,1)
      ENDDO
      CALL MATMULT_FFF ( ES, ALP, 6, 6, 1, EALP )

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         DO J=1,NTSUB
            TBAR(J) = ZERO
            DO K=1,ELGP
               TBAR(J) = TBAR(J) + DT(K,J)
            ENDDO
            TBAR(J) = TBAR(J)/ELGP - TREF1
         ENDDO
      ENDIF

      JJ = 0
      DO I=1,ELGP
         II = 6*(I - 1)
         DO J=1,3
            II = II + 1
            JJ = JJ + 1
            ID(JJ) = II
         ENDDO
      ENDDO

      IF (TYPE == 'PYRA5   ') THEN
         CALL PYRA5_POINTS ( NGP, INT_R, INT_S, INT_T, INT_W )
      ELSE IF (TYPE == 'PYRA14  ') THEN
         CALL PYRA14_POINTS ( NGP, INT_R, INT_S, INT_T, INT_W )
      ELSE
         RETURN
      ENDIF

      VOLUME = ZERO
      DO I=1,NGP
         CALL PYRA_BMAT ( INT_R(I), INT_S(I), INT_T(I), DETJ, B, N )
         VOLUME = VOLUME + INT_W(I)*DETJ
      ENDDO

      IF (VOLUME < EPS1) THEN
         WRITE(ERR,1925) EID, TYPE, 'VOLUME', VOLUME
         WRITE(F06,1925) EID, TYPE, 'VOLUME', VOLUME
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

      IF (OPT(1) == 'Y') THEN
         RHO1 = RHO(1)
         IF (RHO1 /= ZERO) THEN
            M_1DOF(:,:) = ZERO
            DO I=1,NGP
               CALL PYRA_BMAT ( INT_R(I), INT_S(I), INT_T(I), DETJ, B, N )
               DO J=1,ELGP
                  DO K=1,ELGP
                     M_1DOF(J,K) = M_1DOF(J,K) + RHO1*INT_W(I)*DETJ*N(J)*N(K)
                  ENDDO
               ENDDO
            ENDDO
            CALL EXPAND_MASS_DOFS ( M_1DOF )
         ENDIF
      ENDIF

      IF (OPT(2) == 'Y') THEN
         DO L=1,NTSUB
            DO J=1,3*ELGP
               DUM0(J) = ZERO
            ENDDO
            DO I=1,NGP
               CALL PYRA_BMAT ( INT_R(I), INT_S(I), INT_T(I), DETJ, B, N )
               CALL MATMULT_FFF_T ( B, EALP, 6, 3*ELGP, 1, DUM3 )
               TEMP = ZERO
               DO J=1,ELGP
                  TEMP = TEMP + N(J)*DT(J,L)
               ENDDO
               TEMP = TEMP - TREF1
               DO J=1,3*ELGP
                  DUM0(J) = DUM0(J) + DUM3(J,1)*TEMP*INT_W(I)*DETJ
               ENDDO
            ENDDO
            DO J=1,3*ELGP
               PTE(ID(J),L) = DUM0(J)
            ENDDO
         ENDDO
      ENDIF

      IF (OPT(4) == 'Y') THEN
         IF ((TYPE == 'PYRA5   ') .AND. (SOLIDTYP == 'NEWSOLID')) THEN
            CALL PYRA5_EAS54_STIFFNESS ( ID, DUM3, EAS_KUA, EAS_KAA, EAS_X, EAS_M, EAS_DM )
         ELSE
            DO I=1,NGP
               CALL PYRA_BMAT ( INT_R(I), INT_S(I), INT_T(I), DETJ, B, N )
               CALL MATMULT_FFF ( ES, B, 6, 6, 3*ELGP, DUM4 )
               CALL MATMULT_FFF_T ( B, DUM4, 6, 3*ELGP, 3*ELGP, DUM5 )
               DO J=1,3*ELGP
                  DO K=1,3*ELGP
                     KE(ID(J),ID(K)) = KE(ID(J),ID(K)) + DUM5(J,K)*INT_W(I)*DETJ
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

         DO I=2,6*ELGP
            DO J=1,I-1
               KE(I,J) = KE(J,I)
            ENDDO
         ENDDO
      ENDIF

      IF ((OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         CALL PYRA_BMAT ( ZERO, ZERO, 0.25D0, DETJ, B, N )
         CALL MATMULT_FFF ( ES, B, 6, 6, 3*ELGP, DUM4 )
         DO I=1,3
            DO J=1,3*ELGP
               BE1(I,ID(J),1) = B(I  ,J)
               BE2(I,ID(J),1) = B(I+3,J)
               SE1(I,ID(J),1) = DUM4(I  ,J)
               SE2(I,ID(J),1) = DUM4(I+3,J)
            ENDDO
         ENDDO
         DO L=1,NTSUB
            DO I=1,3
               STE1(I,L,1) = EALP(I)*TBAR(L)
            ENDDO
         ENDDO
      ENDIF

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN

         CALL ELMDIS
         CALL ELEM_STRE_STRN_ARRAYS ( 1 )

         SIGxx = STRESS(1)
         SIGyy = STRESS(2)
         SIGzz = STRESS(3)
         SIGxy = STRESS(4)
         SIGyz = STRESS(5)
         SIGzx = STRESS(6)

         KWW(1,1) = SIGyy + SIGzz   ;   KWW(1,2) = -SIGxy   ;   KWW(1,3) = -SIGzx
         KWW(2,2) = SIGxx + SIGzz   ;   KWW(2,3) = -SIGyz
         KWW(3,3) = SIGxx + SIGyy
         KWW(2,1) = KWW(1,2)
         KWW(3,1) = KWW(1,3)
         KWW(3,2) = KWW(2,3)

         DUM3(:,:) = ZERO
         DO I=1,NGP
            CALL PYRA_BMAT ( INT_R(I), INT_S(I), INT_T(I), DETJ, B, N, DNX )
            CBAR(:,:) = ZERO
            DO L=1,ELGP
               CBAR(1,3*(L-1)+1) =  ZERO              ;  CBAR(1,3*(L-1)+2) = -HALF*DNX(3,L)  ;  CBAR(1,3*(L-1)+3) =  HALF*DNX(2,L)
               CBAR(2,3*(L-1)+1) =  HALF*DNX(3,L)    ;  CBAR(2,3*(L-1)+2) =  ZERO           ;  CBAR(2,3*(L-1)+3) = -HALF*DNX(1,L)
               CBAR(3,3*(L-1)+1) = -HALF*DNX(2,L)    ;  CBAR(3,3*(L-1)+2) =  HALF*DNX(1,L)  ;  CBAR(3,3*(L-1)+3) =  ZERO
            ENDDO
            CALL MATMULT_FFF ( KWW, CBAR, 3, 3, 3*ELGP, DUM6 )
            CALL MATMULT_FFF_T ( CBAR, DUM6, 3, 3*ELGP, 3*ELGP, DUM5 )
            DO J=1,3*ELGP
               DO K=1,3*ELGP
                  DUM3(J,K) = DUM3(J,K) + DUM5(J,K)*INT_W(I)*DETJ
               ENDDO
            ENDDO
         ENDDO

         DO I=1,3*ELGP
            DO J=1,3*ELGP
               KED(ID(I),ID(J)) = DUM3(I,J)
            ENDDO
         ENDDO

         DO I=2,6*ELGP
            DO J=1,I-1
               KED(I,J) = KED(J,I)
            ENDDO
         ENDDO

         IF (DEBUG(178) >= 1) THEN
            WRITE(F06,2001) TRIM(SUBR_NAME), OPT(6), LOAD_ISTEP, TYPE, EID
         ENDIF

      ENDIF

      RETURN

 1925 FORMAT(' *ERROR  1925: ELEMENT ',I8,', TYPE ',A,' HAS NONPOSITIVE ',A,' = ',1ES14.6)
 2001 FORMAT(' In ', A, ' with OPT(6) = ', A, ', and LOAD_ISTEP = ', I8, ': KED Differential Stiffness Matrix for ', A,            &
                    ' element ', I8)

      CONTAINS

      SUBROUTINE PYRA5_POINTS ( NPT, RR, SS, TT, WW )

      INTEGER(LONG), INTENT(OUT)      :: NPT
      REAL(DOUBLE),  INTENT(OUT)      :: RR(:),SS(:),TT(:),WW(:)

      INTEGER(LONG)                   :: A,B,C
      REAL(DOUBLE)                    :: G(2),GT(2),W(2),WT(2)

      G(1)  = -ONE/DSQRT(3.0D0)
      G(2)  =  ONE/DSQRT(3.0D0)
      W(1)  =  ONE
      W(2)  =  ONE
      GT(1) =  0.5D0*(ONE - ONE/DSQRT(3.0D0))
      GT(2) =  0.5D0*(ONE + ONE/DSQRT(3.0D0))
      WT(1) =  0.5D0
      WT(2) =  0.5D0

      NPT = 0
      DO A=1,2
         DO B=1,2
            DO C=1,2
               NPT = NPT + 1
               RR(NPT) = G(A)
               SS(NPT) = G(B)
               TT(NPT) = GT(C)
               WW(NPT) = W(A)*W(B)*WT(C)
            ENDDO
         ENDDO
      ENDDO

      END SUBROUTINE PYRA5_POINTS

      SUBROUTINE PYRA14_POINTS ( NPT, RR, SS, TT, WW )

      INTEGER(LONG), INTENT(OUT)      :: NPT
      REAL(DOUBLE),  INTENT(OUT)      :: RR(:),SS(:),TT(:),WW(:)

      INTEGER(LONG)                   :: A,B
      REAL(DOUBLE)                    :: BARY(4,4)
      REAL(DOUBLE)                    :: TET(4,4,3)
      REAL(DOUBLE)                    :: VOL

      BARY(1,1) = 0.5854101966249685D0; BARY(1,2) = 0.1381966011250105D0; BARY(1,3) = 0.1381966011250105D0; BARY(1,4) = 0.1381966011250105D0
      BARY(2,1) = 0.1381966011250105D0; BARY(2,2) = 0.5854101966249685D0; BARY(2,3) = 0.1381966011250105D0; BARY(2,4) = 0.1381966011250105D0
      BARY(3,1) = 0.1381966011250105D0; BARY(3,2) = 0.1381966011250105D0; BARY(3,3) = 0.5854101966249685D0; BARY(3,4) = 0.1381966011250105D0
      BARY(4,1) = 0.1381966011250105D0; BARY(4,2) = 0.1381966011250105D0; BARY(4,3) = 0.1381966011250105D0; BARY(4,4) = 0.5854101966249685D0

      TET(:,:,:) = ZERO
      TET(1,2,1) = -ONE; TET(1,2,2) = -ONE; TET(1,3,1) =  ONE; TET(1,3,2) = -ONE; TET(1,4,3) = ONE
      TET(2,2,1) =  ONE; TET(2,2,2) = -ONE; TET(2,3,1) =  ONE; TET(2,3,2) =  ONE; TET(2,4,3) = ONE
      TET(3,2,1) =  ONE; TET(3,2,2) =  ONE; TET(3,3,1) = -ONE; TET(3,3,2) =  ONE; TET(3,4,3) = ONE
      TET(4,2,1) = -ONE; TET(4,2,2) =  ONE; TET(4,3,1) = -ONE; TET(4,3,2) = -ONE; TET(4,4,3) = ONE

      VOL = 1.0D0/3.0D0
      NPT = 0
      DO A=1,4
         DO B=1,4
            NPT = NPT + 1
            RR(NPT) = BARY(B,1)*TET(A,1,1) + BARY(B,2)*TET(A,2,1) + BARY(B,3)*TET(A,3,1) + BARY(B,4)*TET(A,4,1)
            SS(NPT) = BARY(B,1)*TET(A,1,2) + BARY(B,2)*TET(A,2,2) + BARY(B,3)*TET(A,3,2) + BARY(B,4)*TET(A,4,2)
            TT(NPT) = BARY(B,1)*TET(A,1,3) + BARY(B,2)*TET(A,2,3) + BARY(B,3)*TET(A,3,3) + BARY(B,4)*TET(A,4,3)
            WW(NPT) = VOL/4.0D0
         ENDDO
      ENDDO

      END SUBROUTINE PYRA14_POINTS

      SUBROUTINE PYRA5_EAS54_POINTS ( NPT, RR, SS, TT, WW )

      INTEGER(LONG), INTENT(OUT)      :: NPT
      REAL(DOUBLE),  INTENT(OUT)      :: RR(:),SS(:),TT(:),WW(:)

      INTEGER(LONG)                   :: A,B,C
      REAL(DOUBLE)                    :: G(3),GT(3),W(3),WT(3)

      G(1)  = -DSQRT(3.0D0/5.0D0)
      G(2)  =  ZERO
      G(3)  =  DSQRT(3.0D0/5.0D0)
      W(1)  =  5.0D0/9.0D0
      W(2)  =  8.0D0/9.0D0
      W(3)  =  5.0D0/9.0D0
      GT(1) =  0.5D0*(ONE - DSQRT(3.0D0/5.0D0))
      GT(2) =  0.5D0
      GT(3) =  0.5D0*(ONE + DSQRT(3.0D0/5.0D0))
      WT(1) =  0.5D0*5.0D0/9.0D0
      WT(2) =  0.5D0*8.0D0/9.0D0
      WT(3) =  0.5D0*5.0D0/9.0D0

      NPT = 0
      DO A=1,3
         DO B=1,3
            DO C=1,3
               NPT = NPT + 1
               RR(NPT) = G(A)
               SS(NPT) = G(B)
               TT(NPT) = GT(C)
               WW(NPT) = W(A)*W(B)*WT(C)
            ENDDO
         ENDDO
      ENDDO

      END SUBROUTINE PYRA5_EAS54_POINTS

      SUBROUTINE PYRA5_EAS54_STIFFNESS ( ID, KUU, KUA, KAA, X, M, DM )

      INTEGER(LONG), INTENT(IN)       :: ID(3*ELGP)
      REAL(DOUBLE), INTENT(OUT)       :: KUU(3*ELGP,3*ELGP)
      REAL(DOUBLE), INTENT(OUT)       :: KUA(3*ELGP,EAS54_N)
      REAL(DOUBLE), INTENT(OUT)       :: KAA(EAS54_N,EAS54_N)
      REAL(DOUBLE), INTENT(OUT)       :: X(EAS54_N,3*ELGP)
      REAL(DOUBLE), INTENT(OUT)       :: M(6,EAS54_N)
      REAL(DOUBLE), INTENT(OUT)       :: DM(6,EAS54_N)

      INTEGER(LONG)                   :: A,BMODE,GP,II,JJ,PP,QQ
      INTEGER(LONG)                   :: NPT
      INTEGER(LONG)                   :: IERR
      REAL(DOUBLE)                    :: DETJ0,DETJ1
      REAL(DOUBLE)                    :: RR(27),SS(27),TT(27),WW(27)
      REAL(DOUBLE)                    :: SCALE
      REAL(DOUBLE)                    :: TRACE,REG

      KUU(:,:) = ZERO
      KUA(:,:) = ZERO
      KAA(:,:) = ZERO

      CALL PYRA_BMAT ( ZERO, ZERO, 0.25D0, DETJ0, B, N )
      CALL PYRA5_EAS54_POINTS ( NPT, RR, SS, TT, WW )

      DO GP=1,NPT
         CALL PYRA_BMAT ( RR(GP), SS(GP), TT(GP), DETJ1, B, N )
         CALL MATMULT_FFF ( ES, B, 6, 6, 3*ELGP, DUM4 )
         CALL MATMULT_FFF_T ( B, DUM4, 6, 3*ELGP, 3*ELGP, DUM5 )
         DO II=1,3*ELGP
            DO JJ=1,3*ELGP
               KUU(II,JJ) = KUU(II,JJ) + DUM5(II,JJ)*WW(GP)*DETJ1
            ENDDO
         ENDDO

         SCALE = DETJ0/DETJ1
         CALL PYRA5_EAS54_MODES ( RR(GP), SS(GP), TT(GP), SCALE, M )
         DM(:,:) = ZERO
         DO A=1,6
            DO BMODE=1,EAS54_N
               DO II=1,6
                  DM(A,BMODE) = DM(A,BMODE) + ES(A,II)*M(II,BMODE)
               ENDDO
            ENDDO
         ENDDO
         DO II=1,3*ELGP
            DO BMODE=1,EAS54_N
               DO A=1,6
                  KUA(II,BMODE) = KUA(II,BMODE) + B(A,II)*DM(A,BMODE)*WW(GP)*DETJ1
               ENDDO
            ENDDO
         ENDDO
         DO A=1,EAS54_N
            DO BMODE=1,EAS54_N
               DO II=1,6
                  KAA(A,BMODE) = KAA(A,BMODE) + M(II,A)*DM(II,BMODE)*WW(GP)*DETJ1
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      TRACE = ZERO
      DO A=1,EAS54_N
         TRACE = TRACE + KAA(A,A)
      ENDDO
      REG = 1.0D-14*DABS(TRACE)
      IF (REG < 1.0D-14) THEN
         REG = 1.0D-14
      ENDIF
      DO A=1,EAS54_N
         KAA(A,A) = KAA(A,A) + REG
      ENDDO

      DO A=1,EAS54_N
         DO II=1,3*ELGP
            X(A,II) = KUA(II,A)
         ENDDO
      ENDDO
      CALL SOLVE_EAS54_SYSTEM ( KAA, X, IERR )
      IF (IERR == 0) THEN
         DO II=1,3*ELGP
            DO JJ=1,3*ELGP
               DO A=1,EAS54_N
                  KUU(II,JJ) = KUU(II,JJ) - KUA(II,A)*X(A,JJ)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      DO II=1,3*ELGP
         DO JJ=1,3*ELGP
            KE(ID(II),ID(JJ)) = HALF*(KUU(II,JJ) + KUU(JJ,II))
         ENDDO
      ENDDO

      END SUBROUTINE PYRA5_EAS54_STIFFNESS

      SUBROUTINE PYRA5_EAS54_MODES ( R, S, T, SCALE, M )

      REAL(DOUBLE), INTENT(IN)        :: R,S,T,SCALE
      REAL(DOUBLE), INTENT(OUT)       :: M(6,EAS54_N)

      INTEGER(LONG)                   :: A,STRAIN
      REAL(DOUBLE)                    :: TB
      REAL(DOUBLE)                    :: MODES(9)

      TB = T - 0.25D0
      MODES(1) = R
      MODES(2) = S
      MODES(3) = TB
      MODES(4) = R*S
      MODES(5) = R*TB
      MODES(6) = S*TB
      MODES(7) = R*R - ONE/3.0D0
      MODES(8) = S*S - ONE/3.0D0
      MODES(9) = TB*TB - 0.0375D0

      M(:,:) = ZERO
      DO STRAIN=1,6
         DO A=1,9
            M(STRAIN,9*(STRAIN-1)+A) = SCALE*MODES(A)
         ENDDO
      ENDDO

      END SUBROUTINE PYRA5_EAS54_MODES

      SUBROUTINE SOLVE_EAS54_SYSTEM ( A_INOUT, B_INOUT, IERR )

      REAL(DOUBLE), INTENT(INOUT)     :: A_INOUT(EAS54_N,EAS54_N)
      REAL(DOUBLE), INTENT(INOUT)     :: B_INOUT(EAS54_N,3*ELGP)
      INTEGER(LONG), INTENT(OUT)      :: IERR

      INTEGER(LONG)                   :: I,J,KROW,KCOL,PIV
      REAL(DOUBLE)                    :: FACTOR,MAXVAL,TMP

      IERR = 0
      DO KCOL=1,EAS54_N
         PIV = KCOL
         MAXVAL = DABS(A_INOUT(KCOL,KCOL))
         DO I=KCOL+1,EAS54_N
            IF (DABS(A_INOUT(I,KCOL)) > MAXVAL) THEN
               MAXVAL = DABS(A_INOUT(I,KCOL))
               PIV = I
            ENDIF
         ENDDO
         IF (MAXVAL <= EPS1) THEN
            IERR = 1
            RETURN
         ENDIF
         IF (PIV /= KCOL) THEN
            DO J=KCOL,EAS54_N
               TMP = A_INOUT(KCOL,J)
               A_INOUT(KCOL,J) = A_INOUT(PIV,J)
               A_INOUT(PIV,J) = TMP
            ENDDO
            DO J=1,3*ELGP
               TMP = B_INOUT(KCOL,J)
               B_INOUT(KCOL,J) = B_INOUT(PIV,J)
               B_INOUT(PIV,J) = TMP
            ENDDO
         ENDIF
         DO I=KCOL+1,EAS54_N
            FACTOR = A_INOUT(I,KCOL)/A_INOUT(KCOL,KCOL)
            A_INOUT(I,KCOL) = ZERO
            DO J=KCOL+1,EAS54_N
               A_INOUT(I,J) = A_INOUT(I,J) - FACTOR*A_INOUT(KCOL,J)
            ENDDO
            DO J=1,3*ELGP
               B_INOUT(I,J) = B_INOUT(I,J) - FACTOR*B_INOUT(KCOL,J)
            ENDDO
         ENDDO
      ENDDO

      DO KROW=EAS54_N,1,-1
         DO J=1,3*ELGP
            TMP = B_INOUT(KROW,J)
            DO I=KROW+1,EAS54_N
               TMP = TMP - A_INOUT(KROW,I)*B_INOUT(I,J)
            ENDDO
            B_INOUT(KROW,J) = TMP/A_INOUT(KROW,KROW)
         ENDDO
      ENDDO

      END SUBROUTINE SOLVE_EAS54_SYSTEM

      SUBROUTINE PYRA_SHAPE ( R, S, T, N, DN )

      REAL(DOUBLE), INTENT(IN)        :: R,S,T
      REAL(DOUBLE), INTENT(OUT)       :: N(ELGP), DN(3,ELGP)

      REAL(DOUBLE)                    :: H

      IF (TYPE == 'PYRA5   ') THEN
         CALL PYRA5_SHAPE ( R, S, T, N, DN )
      ELSE
         CALL LIU_SHAPE ( R, S, T, N )
         H = 1.0D-6
         CALL LIU_DERIV ( R, S, T, H, DN )
      ENDIF

      END SUBROUTINE PYRA_SHAPE

      SUBROUTINE PYRA5_SHAPE ( R, S, T, N, DN )

      REAL(DOUBLE), INTENT(IN)        :: R,S,T
      REAL(DOUBLE), INTENT(OUT)       :: N(ELGP), DN(3,ELGP)
      REAL(DOUBLE)                    :: Q

      Q = ONE - T
      N(1) = 0.25D0*(ONE - R)*(ONE - S)*Q
      N(2) = 0.25D0*(ONE + R)*(ONE - S)*Q
      N(3) = 0.25D0*(ONE + R)*(ONE + S)*Q
      N(4) = 0.25D0*(ONE - R)*(ONE + S)*Q
      N(5) = T

      DN(:,:) = ZERO
      DN(1,1) = -0.25D0*(ONE - S)*Q
      DN(1,2) =  0.25D0*(ONE - S)*Q
      DN(1,3) =  0.25D0*(ONE + S)*Q
      DN(1,4) = -0.25D0*(ONE + S)*Q
      DN(2,1) = -0.25D0*(ONE - R)*Q
      DN(2,2) = -0.25D0*(ONE + R)*Q
      DN(2,3) =  0.25D0*(ONE + R)*Q
      DN(2,4) =  0.25D0*(ONE - R)*Q
      DN(3,1) = -0.25D0*(ONE - R)*(ONE - S)
      DN(3,2) = -0.25D0*(ONE + R)*(ONE - S)
      DN(3,3) = -0.25D0*(ONE + R)*(ONE + S)
      DN(3,4) = -0.25D0*(ONE - R)*(ONE + S)
      DN(3,5) =  ONE

      END SUBROUTINE PYRA5_SHAPE

      SUBROUTINE PYRA_BMAT ( R, S, T, DETJ, B, N, DNX_OUT )

      REAL(DOUBLE), INTENT(IN)        :: R,S,T
      REAL(DOUBLE), INTENT(OUT)       :: DETJ
      REAL(DOUBLE), INTENT(OUT)       :: B(6,3*ELGP)
      REAL(DOUBLE), INTENT(OUT)       :: N(ELGP)
      REAL(DOUBLE), INTENT(OUT), OPTIONAL :: DNX_OUT(3,ELGP)

      INTEGER(LONG)                   :: A,C
      INTEGER(LONG)                   :: P,Q
      REAL(DOUBLE)                    :: DN(3,ELGP)
      REAL(DOUBLE)                    :: DNX(3,ELGP)
      REAL(DOUBLE)                    :: INVJ(3,3)
      REAL(DOUBLE)                    :: JAC(3,3)

      CALL PYRA_SHAPE ( R, S, T, N, DN )

      JAC(:,:) = ZERO
      DO A=1,ELGP
         DO P=1,3
            DO Q=1,3
               JAC(P,Q) = JAC(P,Q) + DN(P,A)*XEL(A,Q)
            ENDDO
         ENDDO
      ENDDO

      CALL INV3 ( JAC, DETJ, INVJ )
      IF (DETJ <= EPS1) THEN
         WRITE(ERR,2925) EID, TYPE, 'JACOBIAN', DETJ
         WRITE(F06,2925) EID, TYPE, 'JACOBIAN', DETJ
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         B(:,:) = ZERO
         RETURN
      ENDIF

      DO A=1,ELGP
         DO P=1,3
            DNX(P,A) = ZERO
            DO Q=1,3
               DNX(P,A) = DNX(P,A) + INVJ(P,Q)*DN(Q,A)
            ENDDO
         ENDDO
      ENDDO

      IF (PRESENT(DNX_OUT)) THEN
         DNX_OUT(:,:) = DNX(:,:)
      ENDIF

      B(:,:) = ZERO
      DO A=1,ELGP
         C = 3*(A-1)
         B(1,C+1) = DNX(1,A)
         B(2,C+2) = DNX(2,A)
         B(3,C+3) = DNX(3,A)
         B(4,C+1) = DNX(2,A)
         B(4,C+2) = DNX(1,A)
         B(5,C+2) = DNX(3,A)
         B(5,C+3) = DNX(2,A)
         B(6,C+1) = DNX(3,A)
         B(6,C+3) = DNX(1,A)
      ENDDO

 2925 FORMAT(' *ERROR  1925: ELEMENT ',I8,', TYPE ',A,' HAS NONPOSITIVE ',A,' = ',1ES14.6)

      END SUBROUTINE PYRA_BMAT

      SUBROUTINE LIU_DERIV ( R, S, T, H, DN )

      REAL(DOUBLE), INTENT(IN)        :: R,S,T,H
      REAL(DOUBLE), INTENT(OUT)       :: DN(3,ELGP)
      REAL(DOUBLE)                    :: NP(ELGP),NM(ELGP)

      CALL LIU_SHAPE ( R+H, S, T, NP ); CALL LIU_SHAPE ( R-H, S, T, NM ); DN(1,:) = (NP(:) - NM(:))/(2.0D0*H)
      CALL LIU_SHAPE ( R, S+H, T, NP ); CALL LIU_SHAPE ( R, S-H, T, NM ); DN(2,:) = (NP(:) - NM(:))/(2.0D0*H)
      CALL LIU_SHAPE ( R, S, T+H, NP ); CALL LIU_SHAPE ( R, S, T-H, NM ); DN(3,:) = (NP(:) - NM(:))/(2.0D0*H)

      END SUBROUTINE LIU_DERIV

      SUBROUTINE LIU_SHAPE ( X, Y, Z, N )

      REAL(DOUBLE), INTENT(IN)        :: X,Y,Z
      REAL(DOUBLE), INTENT(OUT)       :: N(ELGP)
      INTEGER(LONG)                   :: A
      INTEGER(LONG)                   :: MM(14)
      REAL(DOUBLE)                    :: Q(14),QM(14)

      MM = (/2,1,4,3,5,6,9,8,7,10,12,11,14,13/)
      CALL LIU_RAW ( X, Y, Z, Q )
      CALL LIU_RAW ( -X, Y, Z, QM )
      DO A=1,14
         N(A) = 0.5D0*(Q(A) + QM(MM(A)))
      ENDDO
      N(5) = Q(5)

      END SUBROUTINE LIU_SHAPE

      SUBROUTINE LIU_RAW ( X, Y, Z, Q )

      REAL(DOUBLE), INTENT(IN)        :: X,Y,Z
      REAL(DOUBLE), INTENT(OUT)       :: Q(14)

      Q(:) = ZERO
      IF (X > Y) THEN
         Q(1)  =  0.25D0*(X + Z)*(X + Z - ONE)*(Y - Z - ONE)*(Y - Z)
         Q(2)  = -0.25D0*(X + Z)*(Y - Z)*((X + Z + ONE)*(-Y + Z + ONE) - 4.0D0*Z) - Z*(X - Y)
         Q(3)  =  0.25D0*(Y - Z)*(X + Z)*(Y - Z + ONE)*(X + Z + ONE)
         Q(4)  =  0.25D0*(X + Z)*(Y - Z)*(Y - Z + ONE)*(X + Z - ONE)
         Q(6)  = -0.5D0*(X + Z - ONE)*(((Y - Z - ONE)*(X + ONE)*Y - Z) + Z*(2.0D0*X + ONE))
         Q(7)  = -0.5D0*(Y - Z + ONE)*(((X + Z + ONE)*(Y - ONE)*X - Z) + Z*(2.0D0*Y + ONE))
         Q(8)  = -0.5D0*(Y - Z + ONE)*(X + Z - ONE)*(X + ONE)*Y
         Q(9)  = -0.5D0*(Y - Z + ONE)*(X + Z - ONE)*(Y - ONE)*X
         Q(10) = (Y - Z + ONE)*(X + Z - ONE)*((Y - ONE)*(X + ONE) + Z*(X - Y + Z + ONE))
         Q(11) = Z*(X + Z - ONE)*(Y - Z - ONE)
         Q(12) = -Z*((X + Z + ONE)*(Y - Z - ONE) + 4.0D0*Z)
         Q(13) = Z*(Y - Z + ONE)*(X + Z + ONE)
         Q(14) = -Z*(Y - Z + ONE)*(X + Z - ONE)
      ELSE
         Q(1)  =  0.25D0*(Y + Z)*(Y + Z - ONE)*(X - Z - ONE)*(X - Z)
         Q(2)  = -0.25D0*(X - Z)*(Y + Z)*(X - Z + ONE)*(-Y - Z + ONE)
         Q(3)  =  0.25D0*(Y + Z)*(X - Z)*(X - Z + ONE)*(Y + Z + ONE)
         Q(4)  =  0.25D0*(X - Z)*(Y + Z)*((X - Z - ONE)*(Y + Z + ONE) + 4.0D0*Z) + Z*(X - Y)
         Q(6)  = -0.5D0*(X - Z + ONE)*(Y + Z - ONE)*(X - ONE)*Y
         Q(7)  = -0.5D0*(X - Z + ONE)*(Y + Z - ONE)*(Y + ONE)*X
         Q(8)  = -0.5D0*(X - Z + ONE)*(((Y + Z + ONE)*(X - ONE)*Y - Z) + Z*(2.0D0*X + ONE))
         Q(9)  = -0.5D0*(Y + Z - ONE)*(((X - Z - ONE)*(Y + ONE)*X - Z) + Z*(2.0D0*Y + ONE))
         Q(10) = (X - Z + ONE)*(Y + Z - ONE)*((Y + ONE)*(X - ONE) - Z*(X - Y - Z - ONE))
         Q(11) = Z*(Y + Z - ONE)*(X - Z - ONE)
         Q(12) = -Z*(X - Z + ONE)*(Y + Z - ONE)
         Q(13) = Z*(X - Z + ONE)*(Y + Z + ONE)
         Q(14) = -Z*((Y + Z + ONE)*(X - Z - ONE) + 4.0D0*Z)
      ENDIF
      Q(5) = Z*(2.0D0*Z - ONE)

      END SUBROUTINE LIU_RAW

      SUBROUTINE INV3 ( A, DET, AI )

      REAL(DOUBLE), INTENT(IN)        :: A(3,3)
      REAL(DOUBLE), INTENT(OUT)       :: DET,AI(3,3)

      DET = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) +                   &
            A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
      IF (DABS(DET) <= EPS1) THEN
         AI(:,:) = ZERO
         RETURN
      ENDIF
      AI(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))/DET
      AI(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))/DET
      AI(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2))/DET
      AI(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))/DET
      AI(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))/DET
      AI(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))/DET
      AI(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1))/DET
      AI(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))/DET
      AI(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))/DET

      END SUBROUTINE INV3

      END SUBROUTINE PYRA
