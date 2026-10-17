module EAS9_BASE
   use PARAMS, only: PENTATYP
   use iso_fortran_env, only: real64
   implicit none
   private
   public :: MAP_NATURAL_EAS9, HEX8_EAS9_STATE, HEX8_POINT, HEX8_EAS9_RECOVERY
contains
   subroutine MAP_NATURAL_EAS9(invj0,mhat,m)
      ! v2: covariant natural strain -> local Cartesian strain at center.
      ! J0(i,a)=dx_a/dxi_i; Ecart=inv(J0)*Enat*transpose(inv(J0)).
      ! Voigt ordering xx,yy,zz,xy,yz,zx uses engineering shear strains.
      real(real64), intent(in) :: invj0(3,3),mhat(6,9)
      real(real64), intent(out) :: m(6,9)
      real(real64) :: e(3,3),c(3,3)
      integer :: k
      do k=1,9
         e=0
         e(1,1)=mhat(1,k); e(2,2)=mhat(2,k); e(3,3)=mhat(3,k)
         e(1,2)=0.5_real64*mhat(4,k); e(2,1)=e(1,2)
         e(2,3)=0.5_real64*mhat(5,k); e(3,2)=e(2,3)
         e(1,3)=0.5_real64*mhat(6,k); e(3,1)=e(1,3)
         c=matmul(invj0,matmul(e,transpose(invj0)))
         m(:,k)=[c(1,1),c(2,2),c(3,3),2*c(1,2),2*c(2,3),2*c(1,3)]
      end do
   end subroutine

   subroutine HEX8_POINT(x,p,n,g,b,invj,det,ierr)
      real(real64),intent(in)::x(8,3),p(3)
      real(real64),intent(out)::n(8),g(3,8),b(6,24),invj(3,3),det
      integer,intent(out)::ierr
      integer,parameter::sx(8)=[-1,1,1,-1,-1,1,1,-1],sy(8)=[-1,-1,1,1,-1,-1,1,1],sz(8)=[-1,-1,-1,-1,1,1,1,1]
      real(real64)::dn(3,8),a(3,3)
      integer::i,c
      do i=1,8
         n(i)=(1+sx(i)*p(1))*(1+sy(i)*p(2))*(1+sz(i)*p(3))/8.0_real64
         dn(1,i)=sx(i)*(1+sy(i)*p(2))*(1+sz(i)*p(3))/8.0_real64
         dn(2,i)=sy(i)*(1+sx(i)*p(1))*(1+sz(i)*p(3))/8.0_real64
         dn(3,i)=sz(i)*(1+sx(i)*p(1))*(1+sy(i)*p(2))/8.0_real64
      end do
      a=matmul(dn,x);b=0;g=0;invj=0;ierr=0
      det=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
         +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
      if(det<=0) then
         ierr=1;return
      end if
      invj(1,:)=[a(2,2)*a(3,3)-a(2,3)*a(3,2),a(1,3)*a(3,2)-a(1,2)*a(3,3),a(1,2)*a(2,3)-a(1,3)*a(2,2)]/det
      invj(2,:)=[a(2,3)*a(3,1)-a(2,1)*a(3,3),a(1,1)*a(3,3)-a(1,3)*a(3,1),a(1,3)*a(2,1)-a(1,1)*a(2,3)]/det
      invj(3,:)=[a(2,1)*a(3,2)-a(2,2)*a(3,1),a(1,2)*a(3,1)-a(1,1)*a(3,2),a(1,1)*a(2,2)-a(1,2)*a(2,1)]/det
      g=matmul(invj,dn)
      do i=1,8
         c=3*(i-1)
         b(1,c+1)=g(1,i);b(2,c+2)=g(2,i);b(3,c+3)=g(3,i)
         b(4,c+1)=g(2,i);b(4,c+2)=g(1,i);b(5,c+2)=g(3,i);b(5,c+3)=g(2,i)
         b(6,c+1)=g(3,i);b(6,c+3)=g(1,i)
      end do
   end subroutine

   subroutine HEX8_MODES(p,invj0,det0,det,m)
      real(real64),intent(in)::p(3),invj0(3,3),det0,det
      real(real64),intent(out)::m(6,9)
      real(real64)::mh(6,9)
      mh=0
      mh(1,1)=p(1);mh(2,2)=p(2);mh(3,3)=p(3)
      mh(4,4)=p(1);mh(4,5)=p(2);mh(5,6)=p(2);mh(5,7)=p(3);mh(6,8)=p(3);mh(6,9)=p(1)
      call MAP_NATURAL_EAS9(invj0,mh*det0/det,m)
   end subroutine

   subroutine HEX8_EAS9_STATE(x,d,alpha,temp,tref,u,point,k,f,beff,eps,sigma,kg,ierr)
      ! Same Kuu/Kua/Kaa and covariant transformation as Python EAS9 v2.
      ! Local alpha = inv(Kaa)*(Fa-Kau*u); no global enrichment DOFs.
      ! Kg is the physical initial-stress operator, not a finite-strain EAS tangent.
      real(real64),intent(in)::x(8,3),d(6,6),alpha(6),temp(8),tref,u(24),point(3)
      real(real64),intent(out)::k(24,24),f(24),beff(6,24),eps(6),sigma(6),kg(24,24)
      integer,intent(out)::ierr
      real(real64)::p(3,8),n(8),g(3,8),b(6,24),iv(3,3),iv0(3,3),det,det0,m(6,9)
      real(real64)::ku(24,9),ka(9,9),fa(9),rhs(9,25),a(9,9),factor,tmp(25),sol(9,25)
      real(real64)::eta(9),thermal(6),bg(6,24,8),mg(6,9,8),gg(3,8,8),dv(8),nt(8)
      real(real64)::sg(3,3),h(8,8),v(6),q
      integer::i,j,l,c,r,piv,t
      k=0;f=0;beff=0;eps=0;sigma=0;kg=0;ku=0;ka=0;fa=0;ierr=0
      call HEX8_POINT(x,[0.0_real64,0.0_real64,0.0_real64],n,g,b,iv0,det0,ierr)
      if(ierr/=0)return
      q=1/sqrt(3.0_real64);t=0
      do l=1,2
         do j=1,2
            do i=1,2
               t=t+1;p(:,t)=[(2*i-3)*q,(2*j-3)*q,(2*l-3)*q]
               call HEX8_POINT(x,p(:,t),n,g,b,iv,det,ierr)
               if(ierr/=0)return
               call HEX8_MODES(p(:,t),iv0,det0,det,m)
               bg(:,:,t)=b;mg(:,:,t)=m;gg(:,:,t)=g;dv(t)=det;nt(t)=dot_product(n,temp)-tref
               thermal=matmul(d,alpha)*nt(t)
               k=k+matmul(transpose(b),matmul(d,b))*det
               ku=ku+matmul(transpose(b),matmul(d,m))*det
               ka=ka+matmul(transpose(m),matmul(d,m))*det
               f=f+matmul(transpose(b),thermal)*det
               fa=fa+matmul(transpose(m),thermal)*det
            end do
         end do
      end do
      a=ka;rhs(:,1:24)=transpose(ku);rhs(:,25)=fa
      ! Partial pivoting, solve all displacement and thermal RHS together.
      do i=1,9
         piv=i-1+maxloc(abs(a(i:9,i)),dim=1)
         if(abs(a(piv,i))<=epsilon(1.0_real64)*maxval(abs(ka)))then
            ierr=2;return
         end if
         do j=1,9
            factor=a(i,j);a(i,j)=a(piv,j);a(piv,j)=factor
         end do
         tmp=rhs(i,:);rhs(i,:)=rhs(piv,:);rhs(piv,:)=tmp
         do r=i+1,9
            factor=a(r,i)/a(i,i);a(r,i:9)=a(r,i:9)-factor*a(i,i:9);rhs(r,:)=rhs(r,:)-factor*rhs(i,:)
         end do
      end do
      sol=0
      do i=9,1,-1
         sol(i,:)=(rhs(i,:)-matmul(a(i,i+1:9),sol(i+1:9,:)))/a(i,i)
      end do
      k=k-matmul(ku,sol(:,1:24));k=(k+transpose(k))*0.5_real64
      f=f-matmul(ku,sol(:,25));eta=sol(:,25)-matmul(sol(:,1:24),u)
      call HEX8_POINT(x,point,n,g,b,iv,det,ierr)
      if(ierr/=0)return
      call HEX8_MODES(point,iv0,det0,det,m)
      beff=b-matmul(m,sol(:,1:24))
      eps=matmul(b,u)+matmul(m,eta)
      sigma=matmul(d,eps-alpha*(dot_product(n,temp)-tref))
      do t=1,8
         v=matmul(d,matmul(bg(:,:,t),u)+matmul(mg(:,:,t),eta)-alpha*nt(t))
         sg(1,:)=[v(1),v(4),v(6)];sg(2,:)=[v(4),v(2),v(5)];sg(3,:)=[v(6),v(5),v(3)]
         h=matmul(transpose(gg(:,:,t)),matmul(sg,gg(:,:,t)))*dv(t)
         do i=1,8
            do j=1,8
               do c=1,3
                  kg(3*(i-1)+c,3*(j-1)+c)=kg(3*(i-1)+c,3*(j-1)+c)+h(i,j)
               end do
            end do
         end do
      end do
      kg=(kg+transpose(kg))*0.5_real64
   end subroutine

   subroutine HEX8_EAS9_RECOVERY(x,d,alpha,temp,tref,u,point,k,f,beff,eps,sigma,kg,ierr)
      ! Ordinary corner output is a trilinear extrapolation of the eight
      ! valid interior EAS states. Never invert a corner mapping for output.
      ! Interior/center requests retain the direct enhanced state.
      real(real64),intent(in)::x(8,3),d(6,6),alpha(6),temp(8),tref,u(24),point(3)
      real(real64),intent(out)::k(24,24),f(24),beff(6,24),eps(6),sigma(6),kg(24,24)
      integer,intent(out)::ierr
      real(real64)::q,p(3),w,bs(6,24),es(6),ss(6)
      integer::i,j,l
      if(maxval(abs(point))<1.0_real64-1.0e-12_real64)then
         call HEX8_EAS9_STATE(x,d,alpha,temp,tref,u,point,k,f,beff,eps,sigma,kg,ierr)
         return
      end if
      q=1/sqrt(3.0_real64);beff=0;eps=0;sigma=0;ierr=0
      do l=1,2
         do j=1,2
            do i=1,2
               p=[(2*i-3)*q,(2*j-3)*q,(2*l-3)*q]
               w=product((1+point*p/(q*q))/2)
               call HEX8_EAS9_STATE(x,d,alpha,temp,tref,u,p,k,f,bs,es,ss,kg,ierr)
               if(ierr/=0)return
               beff=beff+w*bs;eps=eps+w*es;sigma=sigma+w*ss
            end do
         end do
      end do
   end subroutine


end module
