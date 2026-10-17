program kernel_probe
use PENTIUM_II_KIND
use SCONTR, only: NTSUB, MELGP, MELDOF, FATAL_ERR, SOL_NAME
use MODEL_STUF
use PARAMS, only: SOLIDTYP, COUPMASS, PENTATYP
use NONLINEAR_PARAMS, only: LOAD_ISTEP
implicit none
integer :: i,j,k,n,ii,jj,station
real(DOUBLE), allocatable :: sn(:),dn(:,:),dx(:,:),bm(:,:)
real(DOUBLE) :: point(3),jac(3,3),ji(3,3),det,gg(10),gw(10),rr(10),ss(10),tw(10)
character(1) :: opt(6), reduced
character(16) :: requested
read(*,*) TYPE,SOLIDTYP,n,reduced
PENTATYP='EAS9' ! Historical EAS9 regression; ABE tested separately.
ELGP=n; EID=1; NTSUB=1; MELGP=n; MELDOF=6*n; LOAD_ISTEP=1
allocate(XEL(n,3),KE(6*n,6*n),ME(6*n,6*n),KED(6*n,6*n),DT(n,1),PTE(6*n,1))
allocate(BE1(3,6*n,n+1),BE2(3,6*n,n+1),SE1(3,6*n,n+1),SE2(3,6*n,n+1),STE1(3,1,n+1))
do i=1,n
 read(*,*) XEL(i,:)
end do
do i=1,6
 read(*,*) ES(i,:)
end do
TE=0;do i=1,3; TE(i,i)=1;end do
ALPVEC=0; ALPVEC(1:3,1)=1d-5; TREF=0;RHO=2
COUPMASS=1; SOL_NAME='MODES'; DT=10
KE=0;ME=0;KED=0;PTE=0;BE1=0;BE2=0;SE1=0;SE2=0;STE1=0
opt=['Y','Y','Y','Y','N','N']
call get_command_argument(1,requested)
if(trim(requested)=='weights') then
 call ORDER_GAUSS(3,gg,gw)
 call ORDER_TRIA(7,rr,ss,tw)
 do i=1,3
  write(*,'(*(ES25.16,1X))') gg(i),gw(i)
 end do
 do i=1,7
  write(*,'(*(ES25.16,1X))') rr(i),ss(i),tw(i)
 end do
 stop
end if
if(trim(requested)=='shape') then
 read(*,*) point
 allocate(sn(n),dn(3,n),dx(3,n),bm(6,3*n))
 if(n==20) then
  call SHP3DH(1,1,1,n,'quadratic probe','shape',3,point(1),point(2),point(3),'N',sn,dn)
 else
  call SHP3DP(1,1,1,n,'quadratic probe','shape',7,3,point(1),point(2),point(3),'N',sn,dn)
 end if
 call JAC3D(point(1),point(2),point(3),dn,'N',jac,ji,det)
 dx=matmul(ji,dn)
 call B3D_ISOPARAMETRIC(dx,0,1,1,1,'all strains','N',bm)
 write(*,'(*(ES25.16,1X))') sn
 do i=1,3
  write(*,'(*(ES25.16,1X))') dn(i,:)
 end do
 write(*,'(ES25.16)') det
 do i=1,6
  write(*,'(*(ES25.16,1X))') bm(i,:)
 end do
 stop
end if
if(trim(requested)=='stiffness') opt=['N','N','N','Y','N','N']
if(trim(requested)=='modal_lumped') COUPMASS=0
station=1
if(trim(requested)=='gradient_corner') then
 DT(:,1)=10d0+4d0*XEL(:,1)
 station=2
end if
select case(TYPE)
case('PYRAM5  ','PYRAM13 ')
 call PYRAM(opt,1,'N','N')
case('PENTA6  ','PENTA15 ')
 call PENTA(opt,1,7,3,reduced,'N')
case('HEXA8   ','HEXA20  ')
 call HEXA(opt,1,3,reduced,'N')
end select
if(FATAL_ERR>0) stop 1
do i=1,3*n
 ii=6*((i-1)/3)+mod(i-1,3)+1
 write(*,'(*(ES25.16,1X))') (KE(ii,6*((j-1)/3)+mod(j-1,3)+1),j=1,3*n)
end do
do i=1,3*n
 ii=6*((i-1)/3)+mod(i-1,3)+1
 write(*,'(*(ES25.16,1X))') (ME(ii,6*((j-1)/3)+mod(j-1,3)+1),j=1,3*n)
end do
do i=1,3*n
 ii=6*((i-1)/3)+mod(i-1,3)+1
 write(*,'(ES25.16)') PTE(ii,1)
end do
do k=1,6
 do j=1,3*n
  jj=6*((j-1)/3)+mod(j-1,3)+1
  if(k<=3) then
   write(*,'(ES25.16)') BE1(k,jj,station)
  else
   write(*,'(ES25.16)') BE2(k-3,jj,station)
  end if
 end do
end do
end program

