! Translational midsurface mass in basic DOFs; no rotary/drilling inertia.
! COUPMASS>0: consistent; otherwise positive scaled diagonal (not row sums).
! Density includes rho*t and PSHELL NSM via MASS_PER_UNIT_AREA.
      SUBROUTINE QUADRATIC_SURFACE_MASS(NN,XYZ,CORR)
      USE PENTIUM_II_KIND, ONLY: LONG,DOUBLE
      USE MODEL_STUF, ONLY: ME,MASS_PER_UNIT_AREA,EID
      USE PARAMS, ONLY: COUPMASS
      USE IOUNT1, ONLY: ERR,F06
      USE SCONTR, ONLY: FATAL_ERR
      USE OUTA_HERE_Interface
      IMPLICIT NONE
      INTEGER(LONG),INTENT(IN) :: NN
      REAL(DOUBLE),INTENT(IN) :: XYZ(NN,3),CORR(8)
      REAL(DOUBLE) :: N(NN),DR(NN),DS(NN),NF(NN),G1(3),G2(3),A(3),M(NN,NN)
      REAL(DOUBLE) :: R,S,L,W,JAC,TOTAL,TRACE
      REAL(DOUBLE),PARAMETER :: GP(5)=(/-0.906179845938664D0,-0.538469310105683D0,0D0, &
                                       0.538469310105683D0,0.906179845938664D0/)
      REAL(DOUBLE),PARAMETER :: GW(5)=(/0.236926885056189D0,0.478628670499366D0,0.568888888888889D0, &
                                       0.478628670499366D0,0.236926885056189D0/)
      REAL(DOUBLE),PARAMETER :: AR(4)=(/-1D0,1D0,1D0,-1D0/),AS(4)=(/-1D0,-1D0,1D0,1D0/)
      INTEGER(LONG) :: I,J,K,H,D
      ME=0D0
      IF (MASS_PER_UNIT_AREA==0D0) RETURN
      M=0D0; TOTAL=0D0
      DO I=1,5
         DO J=1,5
            R=GP(I); S=GP(J); W=GW(I)*GW(J)
               IF (NN==8) THEN
                  N(1:4)=0.25D0*(1D0+AR*R)*(1D0+AS*S)*(AR*R+AS*S-1D0)
                  DR(1:4)=0.25D0*AR*(1D0+AS*S)*(2D0*AR*R+AS*S)
                  DS(1:4)=0.25D0*AS*(1D0+AR*R)*(AR*R+2D0*AS*S)
                  N(5:8)=(/0.5D0*(1D0-R*R)*(1D0-S),0.5D0*(1D0+R)*(1D0-S*S), &
                            0.5D0*(1D0-R*R)*(1D0+S),0.5D0*(1D0-R)*(1D0-S*S)/)
                  DR(5:8)=(/-R*(1D0-S),0.5D0*(1D0-S*S),-R*(1D0+S),-0.5D0*(1D0-S*S)/)
                  DS(5:8)=(/-0.5D0*(1D0-R*R),-(1D0+R)*S,0.5D0*(1D0-R*R),-(1D0-R)*S/)
                  NF=N+(1D0-R*R)*(1D0-S*S)*CORR
               ELSE
                  R=(GP(I)+1D0)/2D0; S=(1D0-R)*(GP(J)+1D0)/2D0
                  W=W*(1D0-R)/4D0; L=1D0-R-S
                  N=(/L*(2D0*L-1D0),R*(2D0*R-1D0),S*(2D0*S-1D0),4D0*L*R,4D0*R*S,4D0*S*L/)
                  DR=(/1D0-4D0*L,4D0*R-1D0,0D0,4D0*(L-R),4D0*S,-4D0*S/)
                  DS=(/1D0-4D0*L,0D0,4D0*S-1D0,-4D0*R,4D0*R,4D0*(L-S)/)
                  NF=N
               ENDIF
               G1=MATMUL(DR,XYZ); G2=MATMUL(DS,XYZ)
               A=(/G1(2)*G2(3)-G1(3)*G2(2),G1(3)*G2(1)-G1(1)*G2(3),G1(1)*G2(2)-G1(2)*G2(1)/)
               JAC=SQRT(DOT_PRODUCT(A,A))
            IF (JAC<=TINY(JAC)) THEN
               FATAL_ERR=FATAL_ERR+1
               WRITE(ERR,*) ' *ERROR: Degenerate mass surface on element ',EID
               WRITE(F06,*) ' *ERROR: Degenerate mass surface on element ',EID
               CALL OUTA_HERE('Y')
            ENDIF
            W=W*JAC*MASS_PER_UNIT_AREA
            TOTAL=TOTAL+W
            DO K=1,NN
               DO H=1,NN
                  M(K,H)=M(K,H)+W*NF(K)*NF(H)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      IF (COUPMASS<=0) THEN
         TRACE=0D0
         DO K=1,NN
            TRACE=TRACE+M(K,K)
         ENDDO
         DO K=1,NN
            W=M(K,K)*TOTAL/TRACE
            M(K,:)=0D0
            M(K,K)=W
         ENDDO
      ENDIF
      DO K=1,NN
         DO H=1,NN
            DO D=1,3
               ME(6*(K-1)+D,6*(H-1)+D)=M(K,H)
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE QUADRATIC_SURFACE_MASS
