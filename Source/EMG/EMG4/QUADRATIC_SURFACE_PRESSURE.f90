! Shared linear midsurface traction integration. Output uses basic DOFs.
! Geometry: standard Q8/T6. Q8 field: standard plus center-bubble correction.
      SUBROUTINE QUADRATIC_SURFACE_PRESSURE(INT_ELEM_ID,NN,XYZ,CORR)
      USE PENTIUM_II_KIND, ONLY: LONG,DOUBLE
      USE MODEL_STUF, ONLY: PPE,PDATA,PPNT,PTYPE,PRESS,EID
      USE IOUNT1, ONLY: ERR,F06
      USE SCONTR, ONLY: FATAL_ERR,NSUB
      USE OUTA_HERE_Interface
      IMPLICIT NONE
      INTEGER(LONG),INTENT(IN) :: INT_ELEM_ID,NN
      REAL(DOUBLE),INTENT(IN) :: XYZ(NN,3),CORR(8)
      REAL(DOUBLE) :: N(NN),DR(NN),DS(NN),NF(NN),PC(4),V(3),G1(3),G2(3),A(3)
      REAL(DOUBLE) :: R,S,L,W,P,JAC,VMAG
      REAL(DOUBLE),PARAMETER :: GP(5)=(/-0.906179845938664D0,-0.538469310105683D0,0D0, &
                                       0.538469310105683D0,0.906179845938664D0/)
      REAL(DOUBLE),PARAMETER :: GW(5)=(/0.236926885056189D0,0.478628670499366D0,0.568888888888889D0, &
                                       0.478628670499366D0,0.236926885056189D0/)
      INTEGER(LONG) :: I,J,K,JS,IP
      REAL(DOUBLE),PARAMETER :: AR(4)=(/-1D0,1D0,1D0,-1D0/),AS(4)=(/-1D0,-1D0,1D0,1D0/)

      PPE(1:6*NN,:)=0D0
      DO JS=1,NSUB
         IP=PPNT(INT_ELEM_ID,JS)
         IF (IP==0) CYCLE
         V=0D0
         IF (PTYPE(INT_ELEM_ID)=='1') THEN
            PC=PRESS(3,JS)
         ELSE
            PC=PDATA(IP:IP+3)
            V=PDATA(IP+5:IP+7)
         ENDIF
         VMAG=SQRT(DOT_PRODUCT(V,V))
         IF (VMAG>0D0) V=V/VMAG
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
                  P=0.25D0*SUM((1D0+AR*R)*(1D0+AS*S)*PC)
               ELSE
                  R=(GP(I)+1D0)/2D0; S=(1D0-R)*(GP(J)+1D0)/2D0
                  W=W*(1D0-R)/4D0; L=1D0-R-S
                  N=(/L*(2D0*L-1D0),R*(2D0*R-1D0),S*(2D0*S-1D0),4D0*L*R,4D0*R*S,4D0*S*L/)
                  DR=(/1D0-4D0*L,4D0*R-1D0,0D0,4D0*(L-R),4D0*S,-4D0*S/)
                  DS=(/1D0-4D0*L,0D0,4D0*S-1D0,-4D0*R,4D0*R,4D0*(L-S)/)
                  NF=N; P=L*PC(1)+R*PC(2)+S*PC(3)
               ENDIF
               G1=MATMUL(DR,XYZ); G2=MATMUL(DS,XYZ)
               A=(/G1(2)*G2(3)-G1(3)*G2(2),G1(3)*G2(1)-G1(1)*G2(3),G1(1)*G2(2)-G1(2)*G2(1)/)
               JAC=SQRT(DOT_PRODUCT(A,A))
               IF (JAC<=TINY(JAC)) THEN
                  FATAL_ERR=FATAL_ERR+1
                  WRITE(ERR,*) ' *ERROR: Degenerate pressure surface on element ',EID
                  WRITE(F06,*) ' *ERROR: Degenerate pressure surface on element ',EID
                  CALL OUTA_HERE('Y')
               ENDIF
               IF (VMAG>0D0) A=JAC*V
               DO K=1,NN
                  PPE(6*(K-1)+1:6*(K-1)+3,JS)=PPE(6*(K-1)+1:6*(K-1)+3,JS)+NF(K)*P*A*W
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE QUADRATIC_SURFACE_PRESSURE
