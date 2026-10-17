program tied_probe
use PARAMS, only: HEXATYP,HEXATAX
use EAS9_BASE, only: HEX8_EAS9_STATE
implicit none
real(8)::x(8,3),d(6,6),alpha(6),temp(8),u(24),point(3),k(24,24),f(24),b(6,24),e(6),s(6),kg(24,24),young,nu,mu,lam
integer::i,j,err
character(1024)::output
read(*,*)HEXATYP,HEXATAX,young,nu
read(*,*)((x(i,j),j=1,3),i=1,8)
read(*,*)u,temp,point
mu=young/(2*(1+nu));lam=young*nu/((1+nu)*(1-2*nu));d=0;d(1:3,1:3)=lam
alpha=[1d-5,1d-5,1d-5,0d0,0d0,0d0]
do i=1,3
d(i,i)=lam+2*mu;d(i+3,i+3)=mu
enddo
call HEX8_EAS9_STATE(x,d,alpha,temp,0d0,u,point,k,f,b,e,s,kg,err)
if(err/=0)stop 2
call get_command_argument(1,output)
open(20,file=trim(output),access='stream',form='unformatted',status='replace')
write(20)k,f,b,e,s,kg
close(20)
end program
