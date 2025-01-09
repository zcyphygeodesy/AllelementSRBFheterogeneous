      subroutine Pndpn_dt(p,dp,n,t)
	implicit none
	integer::n,m,k
	real*8 ::p(n+5),dp(n+5),t
!---------------------------------------------------------------------------
!	计算全部Legendre函数pn及其导数dpn.p(1)=p0,dp(1)=dp0
      p(1)=1.d0;p(2)=t;p(3)=1.5d0*t**2-0.5d0;p(4)=2.5d0*t**3-1.5d0*t
      dp(1)=0.d0;dp(2)=1.d0;dp(3)=3.d0*t;dp(4)=7.5d0*t**2-1.5d0
      do m=3,n+1
         k=m-1
	   p(m)=real(2*k-1)/real(k)*t*p(m-1)-real(k-1)/real(k)*p(m-2)
	   dp(m)=t*dp(k)+real(m)*p(k)
      enddo
      end
!
!******************************************************************************
!
      subroutine LegPn_dt1(pn,dp1,n,t)
      !计算Pn(t)及其对θ导数t=cosθ
	implicit none
	integer::n,k
	real*8::pn(n),dp1(n),t,u
!---------------------------------------------------------------------------
      u=dsqrt(1-t**2)
      pn(1)=t;pn(2)=1.5d0*t**2-0.5d0
      dp1(1)=-u;dp1(2)=-3.d0*t*u
      do k=3,n
        pn(k)=dble(2*k-1)/dble(k)*t*pn(k-1)-dble(k-1)*pn(k-2)/dble(k)
        dp1(k)=dble(2*k-1)*(t*dp1(k-1)-u*pn(k-1))-dble(k-1)*dp1(k-2)
        dp1(k)=dp1(k)/dble(k)
      enddo
      end
!
!******************************************************************************
!
      subroutine LegPn_dt2(pn,dp1,dp2,n,t)
      !Legendre function and its first and second derivatives to θ. 
	implicit none
	integer::n,k
	real*8::pn(n),dp1(n),dp2(n),t,u
!---------------------------------------------------------------------------
      u=dsqrt(1-t**2)
      pn(1)=t;pn(2)=1.5d0*t**2-0.5d0
      dp1(1)=-u;dp1(2)=-3.d0*t*u
      dp2(1)=-t;dp2(2)=3.d0*(1.d0-2.d0*t**2)
      do k=3,n
        pn(k)=dble(2*k-1)/dble(k)*t*pn(k-1)-dble(k-1)*pn(k-2)/dble(k)
        dp1(k)=dble(2*k-1)*(t*dp1(k-1)-u*pn(k-1))-dble(k-1)*dp1(k-2)
        dp1(k)=dp1(k)/dble(k)
        dp2(k)=(t*dp2(k-1)-2.d0*u*dp1(k-1)-t*pn(k-1))*dble(2*k-1)
        dp2(k)=(dp2(k)-dble(k-1)*dp2(k-2))/dble(k)
      enddo
      end
