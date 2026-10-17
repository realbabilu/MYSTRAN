module EAS9_FRAME
   use PARAMS, only: PENTATYP, HEXATYP, HEXATAX
   use iso_fortran_env, only: real64
   implicit none
   private
   public :: MAP_NATURAL_EAS9, HEX8_EAS9_STATE, HEX8_POINT, PENTA6_ANSEAS_STATE, HEX8_EAS9_RECOVERY, HEX8_AXIS
contains
   subroutine HEX8_AXIS(x,axis,lengths,fallback)
      real(real64),intent(in)::x(8,3)
      integer,intent(out)::axis
      real(real64),intent(out)::lengths(3)
      logical,intent(out)::fallback
      integer,parameter::edges(2,4,3)=reshape([1,2,4,3,5,6,8,7,1,4,2,3,5,8,6,7,1,5,2,6,3,7,4,8],[2,4,3])
      real(real64)::len(4,3),v(3),jac(3,3),spread,skew
      integer::a,i,j
      do a=1,3
         jac(a,:)=0
         do i=1,4
            v=x(edges(2,i,a),:)-x(edges(1,i,a),:)
            len(i,a)=sqrt(sum(v*v));jac(a,:)=jac(a,:)+v/4
         enddo
         lengths(a)=sum(len(:,a))/4
      enddo
      axis=3;fallback=.false.
      select case(trim(HEXATAX))
      case('XI');axis=1
      case('ETA');axis=2
      case('AUTO')
         axis=minloc(lengths,dim=1)
         spread=maxval(maxval(len,dim=1)/max(minval(len,dim=1),tiny(1.0_real64)))
         skew=0
         do i=1,3
            do j=i+1,3
               skew=max(skew,abs(dot_product(jac(i,:),jac(j,:)))/ &
                  max(sqrt(sum(jac(i,:)**2)*sum(jac(j,:)**2)),tiny(1.0_real64)))
            enddo
         enddo
         v=lengths;v(axis)=huge(1.0_real64)
         fallback=lengths(axis)>0.5_real64*minval(v).or.spread>1.5_real64.or.skew>0.3_real64
         if(fallback)axis=3
      end select
   end subroutine

   subroutine HEX8_PROJECTORS(x,axis,px,py,ierr)
      real(real64),intent(in)::x(8,3)
      integer,intent(in)::axis
      real(real64),intent(out)::px(6,6),py(6,6)
      integer,intent(out)::ierr
      integer,parameter::signs(8,3)=reshape([-1,1,1,-1,-1,1,1,-1, &
         -1,-1,1,1,-1,-1,1,1,-1,-1,-1,-1,1,1,1,1],[8,3])
      real(real64)::jac(3,3),rot(3,3),t(3),n(3),v(6),e(3,3),local(3,3),part(3,3),out(3,3),norm
      integer::a,j,h
      jac=matmul(transpose(real(signs,real64)),x)/8
      n=jac(axis,:);norm=sqrt(sum(n*n));ierr=1
      if(norm<=tiny(norm))return
      n=n/norm;t=jac(mod(axis,3)+1,:);t=t-n*dot_product(n,t);norm=sqrt(sum(t*t))
      if(norm<=epsilon(norm)*maxval(abs(jac)))return
      t=t/norm;rot(:,1)=t;rot(:,3)=n
      rot(:,2)=[n(2)*t(3)-n(3)*t(2),n(3)*t(1)-n(1)*t(3),n(1)*t(2)-n(2)*t(1)]
      do j=1,6
         v=0;v(j)=1;e=0;e(1,1)=v(1);e(2,2)=v(2);e(3,3)=v(3)
         e(1,2)=v(4)/2;e(2,1)=e(1,2);e(2,3)=v(5)/2;e(3,2)=e(2,3)
         e(1,3)=v(6)/2;e(3,1)=e(1,3)
         local=matmul(transpose(rot),matmul(e,rot))
         do h=1,2
            part=0;part(h,3)=local(h,3);part(3,h)=local(3,h)
            out=matmul(rot,matmul(part,transpose(rot)))
            v=[out(1,1),out(2,2),out(3,3),2*out(1,2),2*out(2,3),2*out(1,3)]
            if(h==1)px(:,j)=v
            if(h==2)py(:,j)=v
         enddo
      enddo
      ierr=0
   end subroutine

   subroutine HEX8_TIED_B(x,p,axis,px,py,b,ierr)
      real(real64),intent(in)::x(8,3),p(3),px(6,6),py(6,6)
      integer,intent(in)::axis
      real(real64),intent(out)::b(6,24)
      integer,intent(out)::ierr
      real(real64)::n(8),g(3,8),iv(3,3),det,bx(6,24),by(6,24),q(3)
      call HEX8_POINT(x,p,n,g,b,iv,det,ierr)
      if(ierr/=0)return
      q=p;q(mod(axis,3)+1)=0
      call HEX8_POINT(x,q,n,g,bx,iv,det,ierr)
      if(ierr/=0)return
      q=p;q(mod(axis+1,3)+1)=0
      call HEX8_POINT(x,q,n,g,by,iv,det,ierr)
      if(ierr/=0)return
      b=b+matmul(px,bx-b)+matmul(py,by-b)
   end subroutine

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
      real(real64),intent(in)::x(8,3),d(6,6),alpha(6),temp(8),tref,u(24),point(3)
      real(real64),intent(out)::k(24,24),f(24),beff(6,24),eps(6),sigma(6),kg(24,24)
      integer,intent(out)::ierr
      if(HEXATYP=='TIED')then
         call HEXA8_TIED_EAS9(x,d,alpha,temp,tref,u,point,k,f,beff,eps,sigma,kg,ierr)
      else
         call HEX8_EAS9_LEGACY(x,d,alpha,temp,tref,u,point,k,f,beff,eps,sigma,kg,ierr)
      endif
   end subroutine

   subroutine HEX8_EAS9_LEGACY(x,d,alpha,temp,tref,u,point,k,f,beff,eps,sigma,kg,ierr)
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

   subroutine HEXA8_TIED_EAS9(x,d,alpha,temp,tref,u,point,k,f,beff,eps,sigma,kg,ierr)
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
      integer::i,j,l,c,r,piv,t,axis,info
      real(real64)::px(6,6),py(6,6),meanb(6,24),bt(6,24),lengths(3)
      logical::fallback
      k=0;f=0;beff=0;eps=0;sigma=0;kg=0;ku=0;ka=0;fa=0;ierr=0
      call HEX8_POINT(x,[0.0_real64,0.0_real64,0.0_real64],n,g,b,iv0,det0,ierr)
      if(ierr/=0)return
      call HEX8_AXIS(x,axis,lengths,fallback)
      call HEX8_PROJECTORS(x,axis,px,py,ierr)
      if(ierr/=0)return
      meanb=0
      q=1/sqrt(3.0_real64);t=0
      do l=1,2
         do j=1,2
            do i=1,2
               t=t+1;p(:,t)=[(2*i-3)*q,(2*j-3)*q,(2*l-3)*q]
               call HEX8_POINT(x,p(:,t),n,g,b,iv,det,ierr)
               if(ierr/=0)return
               call HEX8_MODES(p(:,t),iv0,det0,det,m)
               bg(:,:,t)=b;mg(:,:,t)=m;gg(:,:,t)=g;dv(t)=det;nt(t)=dot_product(n,temp)-tref
               call HEX8_TIED_B(x,p(:,t),axis,px,py,bt,ierr)
               if(ierr/=0)return
               meanb=meanb+(bt-b)*det
               bg(:,:,t)=bt
            end do
         end do
      end do
      meanb=meanb/sum(dv)
      do t=1,8
               b=bg(:,:,t)-meanb;bg(:,:,t)=b;m=mg(:,:,t);det=dv(t)
               thermal=matmul(d,alpha)*nt(t)
               k=k+matmul(transpose(b),matmul(d,b))*det
               ku=ku+matmul(transpose(b),matmul(d,m))*det
               ka=ka+matmul(transpose(m),matmul(d,m))*det
               f=f+matmul(transpose(b),thermal)*det
               fa=fa+matmul(transpose(m),thermal)*det
      end do
      ! Symmetric positive definite factorization with reciprocal condition check.
      a=ka;rhs(:,1:24)=transpose(ku);rhs(:,25)=fa
      call dpotrf('U',9,a,9,info)
      if(info/=0)then
         ierr=2;return
      endif
      block
         real(real64)::rcond,work(27),anorm
         integer::iw(9)
         anorm=maxval(sum(abs(ka),dim=1))
         call dpocon('U',9,a,9,anorm,rcond,work,iw,info)
         if(info/=0.or.rcond<1.0e-12_real64)then
            ierr=2;return
         endif
      end block
      sol=rhs
      call dpotrs('U',9,25,a,9,sol,9,info)
      if(info/=0)then
         ierr=2;return
      endif
      k=k-matmul(ku,sol(:,1:24));k=(k+transpose(k))*0.5_real64
      f=f-matmul(ku,sol(:,25));eta=sol(:,25)-matmul(sol(:,1:24),u)
      call HEX8_POINT(x,point,n,g,b,iv,det,ierr)
      if(ierr/=0)return
      call HEX8_MODES(point,iv0,det0,det,m)
      call HEX8_TIED_B(x,point,axis,px,py,b,ierr)
      if(ierr/=0)return
      b=b-meanb
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

   subroutine PENTA6_POINT(x,p,n,g,b,jac,iv,det,ierr)
      real(real64),intent(in)::x(6,3),p(3)
      real(real64),intent(out)::n(6),g(3,6),b(6,18),jac(3,3),iv(3,3),det
      integer,intent(out)::ierr
      real(real64)::dn(3,6),l(3),lo,hi,dr(3),ds(3)
      integer::i,c
      l=[1-p(1)-p(2),p(1),p(2)];lo=(1-p(3))/2;hi=(1+p(3))/2
      dr=[-1.0_real64,1.0_real64,0.0_real64];ds=[-1.0_real64,0.0_real64,1.0_real64]
      n(1:3)=l*lo;n(4:6)=l*hi
      dn(1,1:3)=dr*lo;dn(1,4:6)=dr*hi;dn(2,1:3)=ds*lo;dn(2,4:6)=ds*hi
      dn(3,1:3)=-l/2;dn(3,4:6)=l/2
      jac=matmul(dn,x);b=0;g=0;iv=0;ierr=0
      det=jac(1,1)*(jac(2,2)*jac(3,3)-jac(2,3)*jac(3,2))-jac(1,2)*(jac(2,1)*jac(3,3)-jac(2,3)*jac(3,1)) &
         +jac(1,3)*(jac(2,1)*jac(3,2)-jac(2,2)*jac(3,1))
      if(det<=0)then
         ierr=1;return
      end if
      iv(1,:)=[jac(2,2)*jac(3,3)-jac(2,3)*jac(3,2),jac(1,3)*jac(3,2)-jac(1,2)*jac(3,3), &
               jac(1,2)*jac(2,3)-jac(1,3)*jac(2,2)]/det
      iv(2,:)=[jac(2,3)*jac(3,1)-jac(2,1)*jac(3,3),jac(1,1)*jac(3,3)-jac(1,3)*jac(3,1), &
               jac(1,3)*jac(2,1)-jac(1,1)*jac(2,3)]/det
      iv(3,:)=[jac(2,1)*jac(3,2)-jac(2,2)*jac(3,1),jac(1,2)*jac(3,1)-jac(1,1)*jac(3,2), &
               jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)]/det
      g=matmul(iv,dn)
      do i=1,6
         c=3*(i-1)
         b(1,c+1)=g(1,i);b(2,c+2)=g(2,i);b(3,c+3)=g(3,i)
         b(4,c+1)=g(2,i);b(4,c+2)=g(1,i);b(5,c+2)=g(3,i);b(5,c+3)=g(2,i)
         b(6,c+1)=g(3,i);b(6,c+3)=g(1,i)
      end do
   end subroutine

   subroutine PENTA6_ANS_POINT(x,p,ti,t0,det0,n,g,b,m,det,ierr)
      ! Claude original edge-tying operator. No affine projector/correction.
      real(real64),intent(in)::x(6,3),p(3),ti(6,6),t0(6,6),det0
      real(real64),intent(out)::n(6),g(3,6),b(6,18),m(6,9),det
      integer,intent(out)::ierr
      real(real64)::jac(3,3),iv(3,3),bn(6,18),bt(6,18),e1(18),e2(18),corr(18),nt(6),gt(3,6),dtmp
      call PENTA6_POINT(x,p,n,g,b,jac,iv,det,ierr)
      if(ierr/=0)return
      if(PENTATYP=='EAS9    ')then
         m(:,1)=ti(:,1)*(p(1)-1.0_real64/3);m(:,2)=ti(:,2)*(p(2)-1.0_real64/3)
         m(:,3)=ti(:,3)*p(3);m(:,4)=ti(:,4)*(p(1)-1.0_real64/3)
         m(:,5)=ti(:,4)*(p(2)-1.0_real64/3);m(:,6)=ti(:,5)*(p(2)-1.0_real64/3)
         m(:,7)=ti(:,5)*p(3);m(:,8)=ti(:,6)*p(3);m(:,9)=ti(:,6)*(p(1)-1.0_real64/3)
         m=m*(det0/det)
         return
      end if
      bn=matmul(t0,b)
      call PENTA6_POINT(x,[0.5_real64,0.0_real64,p(3)],nt,gt,bt,jac,iv,dtmp,ierr)
      if(ierr/=0)return
      bt=matmul(t0,bt);e1=bt(6,:)
      call PENTA6_POINT(x,[0.0_real64,0.5_real64,p(3)],nt,gt,bt,jac,iv,dtmp,ierr)
      if(ierr/=0)return
      bt=matmul(t0,bt);e2=bt(5,:)
      call PENTA6_POINT(x,[0.5_real64,0.5_real64,p(3)],nt,gt,bt,jac,iv,dtmp,ierr)
      if(ierr/=0)return
      bt=matmul(t0,bt);corr=e2-e1-(bt(5,:)-bt(6,:))
      bn(6,:)=e1+corr*p(2);bn(5,:)=e2-corr*p(1)
      b=matmul(ti,bn)
      m(:,1)=ti(:,3)*p(3);m(:,2)=ti(:,5)*p(3);m(:,3)=ti(:,6)*p(3)
      m(:,4)=ti(:,1)*(p(1)-1.0_real64/3);m(:,5)=ti(:,1)*(p(2)-1.0_real64/3)
      m(:,6)=ti(:,2)*(p(1)-1.0_real64/3);m(:,7)=ti(:,2)*(p(2)-1.0_real64/3)
      m(:,8)=ti(:,4)*(p(1)-1.0_real64/3);m(:,9)=ti(:,4)*(p(2)-1.0_real64/3)
      m=m*(det0/det)
   end subroutine

   subroutine PENTA6_ANSEAS_STATE(x,d,alpha,temp,tref,u,point,k,f,beff,eps,sigma,kg,ierr,cache_allowed)
      ! ABE edge-tied ANS plus nine enhanced modes.
      ! Kg below is physical initial stress with standard displacement gradients;
      ! it is not a derived finite-strain ANS/EAS consistent tangent.
      real(real64),intent(in)::x(6,3),d(6,6),alpha(6),temp(6),tref,u(18),point(3)
      real(real64),intent(out)::k(18,18),f(18),beff(6,18),eps(6),sigma(6),kg(18,18)
      integer,intent(out)::ierr
      real(real64)::n(6),g(3,6),b(6,18),jac(3,3),iv(3,3),j0(3,3),iv0(3,3),det0,det
      real(real64)::basis(6,9),mapped(6,9),ti(6,6),t0(6,6),m(6,9),ku(18,9),ka(9,9),fa(9),eta(9),thermal(6)
      real(real64)::pts(3,6),bb(6,18,6),mm(6,9,6),gg(3,6,6),dv(6),nt(6),tr(2,3),q,z,p(3),v(6),sg(3,3),h(6,6)
      real(real64)::rhs(9,19),sol(9,19),a(9,9),rcond,anorm,work(27)
      integer::i,j,c,t,l,iwork(9),info
      logical,intent(in),optional::cache_allowed
      ! One-entry cache of point-independent state; exact input key, no approximation.
      logical,save::valid=.false.
      logical::want_cache,hit
      character(8),save::cached_variant=''
      real(real64),save::cx(6,3),cd(6,6),calpha(6),ctemp(6),ctref,cu(18)
      real(real64),save::cti(6,6),ct0(6,6),cdet0,csol(9,19),ceta(9),ck(18,18),cf(18),ckg(18,18)
!$omp threadprivate(valid,cached_variant,cx,cd,calpha,ctemp,ctref,cu,cti,ct0,cdet0,csol,ceta,ck,cf,ckg)
      k=0;f=0;beff=0;eps=0;sigma=0;kg=0;ku=0;ka=0;fa=0;ierr=0
      want_cache=.true.
      if(present(cache_allowed))want_cache=cache_allowed
      hit=.false.
      if(want_cache.and.valid)then
         hit=cached_variant==PENTATYP.and.all(cx==x).and.all(cd==d).and.all(calpha==alpha).and. &
             all(ctemp==temp).and.ctref==tref.and.all(cu==u)
      end if
      if(hit)then
         k=ck;f=cf;kg=ckg;ti=cti;t0=ct0;det0=cdet0;sol=csol;eta=ceta
         go to 100
      end if
      call PENTA6_POINT(x,[1.0_real64/3,1.0_real64/3,0.0_real64],n,g,b,j0,iv0,det0,ierr)
      if(ierr/=0)return
      basis=0
      do i=1,6
         basis(i,i)=1
      end do
      call MAP_NATURAL_EAS9(iv0,basis,mapped);ti=mapped(:,1:6)
      call MAP_NATURAL_EAS9(j0,basis,mapped);t0=mapped(:,1:6)
      tr(:,1)=[1.0_real64/6,1.0_real64/6];tr(:,2)=[2.0_real64/3,1.0_real64/6]
      tr(:,3)=[1.0_real64/6,2.0_real64/3];q=1/sqrt(3.0_real64);t=0
      do i=1,3
         do j=1,2
            t=t+1;p=[tr(1,i),tr(2,i),(2*j-3)*q];pts(:,t)=p
            call PENTA6_ANS_POINT(x,p,ti,t0,det0,n,g,b,m,det,ierr)
            if(ierr/=0)return
            dv(t)=det/6;nt(t)=dot_product(n,temp)-tref;bb(:,:,t)=b;mm(:,:,t)=m;gg(:,:,t)=g
            thermal=matmul(d,alpha)*nt(t)
            k=k+matmul(transpose(b),matmul(d,b))*dv(t)
            ku=ku+matmul(transpose(b),matmul(d,m))*dv(t)
            ka=ka+matmul(transpose(m),matmul(d,m))*dv(t)
            f=f+matmul(transpose(b),thermal)*dv(t);fa=fa+matmul(transpose(m),thermal)*dv(t)
         end do
      end do
      ! Positive definite local solve, condition gate, no regularization.
      a=ka;rhs(:,1:18)=transpose(ku);rhs(:,19)=fa
      anorm=maxval(sum(abs(ka),dim=1))
      call dpotrf('U',9,a,9,info)
      if(info/=0)then
         ierr=2;return
      end if
      call dpocon('U',9,a,9,anorm,rcond,work,iwork,info)
      if(info/=0.or.rcond<1.0e-12_real64)then
         ierr=3;return
      end if
      sol=rhs
      call dpotrs('U',9,19,a,9,sol,9,info)
      if(info/=0)then
         ierr=2;return
      end if
      k=k-matmul(ku,sol(:,1:18));k=(k+transpose(k))/2
      f=f-matmul(ku,sol(:,19));eta=sol(:,19)-matmul(sol(:,1:18),u)
100   call PENTA6_ANS_POINT(x,point,ti,t0,det0,n,g,b,m,det,ierr)
      if(ierr/=0)return
      beff=b-matmul(m,sol(:,1:18))
      eps=matmul(b,u)+matmul(m,eta);sigma=matmul(d,eps-alpha*(dot_product(n,temp)-tref))
      if(hit)return
      do t=1,6
         v=matmul(d,matmul(bb(:,:,t),u)+matmul(mm(:,:,t),eta)-alpha*nt(t))
         sg(1,:)=[v(1),v(4),v(6)];sg(2,:)=[v(4),v(2),v(5)];sg(3,:)=[v(6),v(5),v(3)]
         h=matmul(transpose(gg(:,:,t)),matmul(sg,gg(:,:,t)))*dv(t)
         do i=1,6
            do j=1,6
               do c=1,3
                  kg(3*(i-1)+c,3*(j-1)+c)=kg(3*(i-1)+c,3*(j-1)+c)+h(i,j)
               end do
            end do
         end do
      end do
      kg=(kg+transpose(kg))/2
      if(want_cache)then
         cx=x;cd=d;calpha=alpha;ctemp=temp;ctref=tref;cu=u;cached_variant=PENTATYP
         cti=ti;ct0=t0;cdet0=det0;csol=sol;ceta=eta;ck=k;cf=f;ckg=kg;valid=.true.
      end if
   end subroutine
end module
