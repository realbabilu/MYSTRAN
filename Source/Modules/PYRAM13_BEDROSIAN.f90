module PYRAM13_BEDROSIAN
   use iso_fortran_env, only: real64
   implicit none
   private
   public :: BEDROSIAN13_SHAPE, BEDROSIAN13_QUADRATURE, BEDROSIAN13_POINT, BEDROSIAN13_RECOVER, BEDROSIAN13_KG
contains
   subroutine BEDROSIAN13_SHAPE(a,b,z,n,dn)
      ! Bedrosian (1992), algebraically cancelled in collapsed coordinates.
      ! Matches handoff/python/cpyram13/bedrosian13.py node/edge order.
      ! dn differentiates (a,b,z), NOT physical pyramid (xi,eta,zeta).
      ! Integrate with the collapsed Jacobian; do not add (1-z)**2 again.
      ! Values and collapsed derivatives are finite at the apex. The physical
      ! Jacobian there is singular, so never use the apex as a volume point.
      real(real64), intent(in) :: a,b,z
      real(real64), intent(out) :: n(13),dn(3,13)
      integer :: i
      integer, parameter :: sx(4)=[-1,1,1,-1], sy(4)=[-1,-1,1,1]
      real(real64) :: q,u,v,h,c
      q=1.0_real64-z
      do i=1,4
         u=1.0_real64+sx(i)*a
         v=1.0_real64+sy(i)*b
         h=sx(i)*a+sy(i)*b
         c=q*h-1.0_real64
         n(i)=0.25_real64*q*u*v*c
         dn(1,i)=0.25_real64*q*sx(i)*v*(c+q*u)
         dn(2,i)=0.25_real64*q*sy(i)*u*(c+q*v)
         dn(3,i)=-0.25_real64*u*v*(2.0_real64*q*h-1.0_real64)
         n(i+9)=z*q*u*v
         dn(1,i+9)=z*q*sx(i)*v
         dn(2,i+9)=z*q*sy(i)*u
         dn(3,i+9)=(1.0_real64-2.0_real64*z)*u*v
      end do
      n(5)=z*(2.0_real64*z-1.0_real64)
      dn(:,5)=[0.0_real64,0.0_real64,4.0_real64*z-1.0_real64]
      n(6)=0.5_real64*q*q*(1-a*a)*(1-b)
      n(7)=0.5_real64*q*q*(1-b*b)*(1+a)
      n(8)=0.5_real64*q*q*(1-a*a)*(1+b)
      n(9)=0.5_real64*q*q*(1-b*b)*(1-a)
      dn(:,6)=[-q*q*a*(1-b),-0.5_real64*q*q*(1-a*a),-q*(1-a*a)*(1-b)]
      dn(:,7)=[0.5_real64*q*q*(1-b*b),-q*q*b*(1+a),-q*(1-b*b)*(1+a)]
      dn(:,8)=[-q*q*a*(1+b),0.5_real64*q*q*(1-a*a),-q*(1-a*a)*(1+b)]
      dn(:,9)=[-0.5_real64*q*q*(1-b*b),-q*q*b*(1-a),-q*(1-b*b)*(1-a)]
   end subroutine

   subroutine BEDROSIAN13_QUADRATURE(order,p,w,npt,ierr)
      integer, intent(in) :: order
      real(real64), intent(out) :: p(3,27),w(27)
      integer, intent(out) :: npt,ierr
      real(real64) :: g(3),wg(3),z(3),wz(3),s
      integer :: i,j,k
      p=0; w=0; npt=0; ierr=0
      select case(order)
      case(2)
         g(1:2)=[-1.0_real64,1.0_real64]/sqrt(3.0_real64); wg=1
         s=sqrt(10.0_real64)
         z(1:2)=[5.0_real64-s,5.0_real64+s]/15.0_real64
         wz(1:2)=[1.0_real64/6+s/48,1.0_real64/6-s/48]
      case(3)
         g=[-sqrt(0.6_real64),0.0_real64,sqrt(0.6_real64)]
         wg=[5.0_real64/9,8.0_real64/9,5.0_real64/9]
         z=[0.072994024073149699_real64,0.34700376603835181_real64,0.70500220988849838_real64]
         wz=[0.15713636106488646_real64,0.14624626925986611_real64,0.029950703008580715_real64]
      case default
         ierr=1; return
      end select
      do i=1,order
         do j=1,order
            do k=1,order
               npt=npt+1; p(:,npt)=[g(i),g(j),z(k)]
               ! Collapsed Jacobian already contains (1-z)^2.
               w(npt)=wg(i)*wg(j)*wz(k)/(1-z(k))**2
            end do
         end do
      end do
   end subroutine

   subroutine BEDROSIAN13_POINT(x,p,n,grad,b,det,ierr)
      real(real64), intent(in) :: x(13,3),p(3)
      real(real64), intent(out) :: n(13),grad(3,13),b(6,39),det
      integer, intent(out) :: ierr
      real(real64) :: dn(3,13),a(3,3),ai(3,3)
      integer :: i,c
      call BEDROSIAN13_SHAPE(p(1),p(2),p(3),n,dn)
      a=matmul(dn,x); grad=0; b=0; ierr=0
      det=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
         -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
      if(det<=0) then
         ierr=2; return
      end if
      ai(1,:)=[a(2,2)*a(3,3)-a(2,3)*a(3,2),a(1,3)*a(3,2)-a(1,2)*a(3,3),a(1,2)*a(2,3)-a(1,3)*a(2,2)]/det
      ai(2,:)=[a(2,3)*a(3,1)-a(2,1)*a(3,3),a(1,1)*a(3,3)-a(1,3)*a(3,1),a(1,3)*a(2,1)-a(1,1)*a(2,3)]/det
      ai(3,:)=[a(2,1)*a(3,2)-a(2,2)*a(3,1),a(1,2)*a(3,1)-a(1,1)*a(3,2),a(1,1)*a(2,2)-a(1,2)*a(2,1)]/det
      grad=matmul(ai,dn)
      do i=1,13
         c=3*(i-1)
         b(1,c+1)=grad(1,i); b(2,c+2)=grad(2,i); b(3,c+3)=grad(3,i)
         b(4,c+1)=grad(2,i); b(4,c+2)=grad(1,i)
         b(5,c+2)=grad(3,i); b(5,c+3)=grad(2,i)
         b(6,c+1)=grad(3,i); b(6,c+3)=grad(1,i)
      end do
   end subroutine

   subroutine BEDROSIAN13_RECOVER(x,d,u,alpha,temp,tref,order,stress,strain,npt,ierr)
      real(real64), intent(in) :: x(13,3),d(6,6),u(39),alpha(6),temp(13),tref
      integer, intent(in) :: order
      real(real64), intent(out) :: stress(6,27),strain(6,27)
      integer, intent(out) :: npt,ierr
      real(real64) :: p(3,27),w(27),n(13),grad(3,13),b(6,39),det
      integer :: k
      stress=0; strain=0
      call BEDROSIAN13_QUADRATURE(order,p,w,npt,ierr)
      if(ierr/=0) return
      do k=1,npt
         call BEDROSIAN13_POINT(x,p(:,k),n,grad,b,det,ierr)
         if(ierr/=0) return
         strain(:,k)=matmul(b,u)
         stress(:,k)=matmul(d,strain(:,k)-alpha*(dot_product(n,temp)-tref))
      end do
   end subroutine

   subroutine BEDROSIAN13_KG(x,stress,order,kg,ierr)
      ! Physical Cauchy prestress in local element axes; compression negative.
      ! No rotation-only surrogate and no global enrichment DOFs.
      real(real64), intent(in) :: x(13,3),stress(6,27)
      integer, intent(in) :: order
      real(real64), intent(out) :: kg(39,39)
      integer, intent(out) :: ierr
      real(real64) :: p(3,27),w(27),n(13),grad(3,13),b(6,39),det,s(3,3),h(13,13)
      integer :: npt,k,i,j,c
      kg=0
      call BEDROSIAN13_QUADRATURE(order,p,w,npt,ierr)
      if(ierr/=0) return
      do k=1,npt
         call BEDROSIAN13_POINT(x,p(:,k),n,grad,b,det,ierr)
         if(ierr/=0) then
            kg=0; return
         end if
         s(1,:)=[stress(1,k),stress(4,k),stress(6,k)]
         s(2,:)=[stress(4,k),stress(2,k),stress(5,k)]
         s(3,:)=[stress(6,k),stress(5,k),stress(3,k)]
         h=matmul(transpose(grad),matmul(s,grad))*w(k)*det
         do i=1,13
            do j=1,13
               do c=1,3
                  kg(3*(i-1)+c,3*(j-1)+c)=kg(3*(i-1)+c,3*(j-1)+c)+h(i,j)
               end do
            end do
         end do
      end do
      kg=0.5_real64*(kg+transpose(kg))
   end subroutine
end module
