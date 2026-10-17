! Exercise real PYRAM recovery operators on nonconstant fields and temperatures.
program solid_recovery_probe
use PENTIUM_II_KIND
use SCONTR, only: NTSUB,MELGP,MELDOF,FATAL_ERR,SOL_NAME
use MODEL_STUF
use PARAMS, only: SOLIDTYP
use NONLINEAR_PARAMS, only: LOAD_ISTEP
implicit none
integer :: n,i,j,p,c,id
real(double),allocatable :: u(:)
real(double) :: strain_out(6),stress_out(6)
character(1) :: opt(6)
read(*,*) TYPE,n
ELGP=n;EID=1;NTSUB=1;MELGP=n;MELDOF=6*n;LOAD_ISTEP=1
SOLIDTYP='NEWSOLID';SOL_NAME='STATICS'
allocate(XEL(n,3),DT(n,1),BE1(3,6*n,n+1),BE2(3,6*n,n+1))
allocate(SE1(3,6*n,n+1),SE2(3,6*n,n+1),STE1(3,1,n+1),STE2(3,1,n+1),u(3*n))
do i=1,n
 read(*,*)XEL(i,:)
end do
do i=1,6
 read(*,*)ES(i,:)
end do
read(*,*)u
read(*,*)DT(:,1)
read(*,*)ALPVEC(:,1)
TREF=0;BE1=0;BE2=0;SE1=0;SE2=0;STE1=0;STE2=0
opt=['N','N','Y','N','N','N']
call PYRAM(opt,1,'N','N')
if(FATAL_ERR>0)stop 1
do p=1,6
 strain_out=0;stress_out=0
 do j=1,3*n
  id=6*((j-1)/3)+mod(j-1,3)+1
  do c=1,3
   strain_out(c)=strain_out(c)+BE1(c,id,p)*u(j)
   strain_out(c+3)=strain_out(c+3)+BE2(c,id,p)*u(j)
   stress_out(c)=stress_out(c)+SE1(c,id,p)*u(j)
   stress_out(c+3)=stress_out(c+3)+SE2(c,id,p)*u(j)
  end do
 end do
 stress_out(1:3)=stress_out(1:3)-STE1(:,1,p)
 stress_out(4:6)=stress_out(4:6)-STE2(:,1,p)
 write(*,'(*(ES25.16,1X))')strain_out,stress_out
end do
end program
