      subroutine RBFvalue(RBF,NF,dr,xx,val)
      !!由夹角xx°内插RBF值val
	implicit none
	integer::NF,i,k
	real*8::xx,dr,val,dt,RBF(NF+1),px,tmp
!---------------------------------------------------------------------------
      val=0.d0;px=0.d0;dt=dr/dble(NF)!间隔°
      k=nint(xx/dt)!与xx最近的数组序号
      do i=k-2,k+2
        if(i<1.or.i>NF+1)goto 1001
        tmp=(i-1.d0)*dt
        val=val+RBF(i)*tmp;px=px+tmp
1001    continue
      enddo
      if(px>1.d-12)val=val/px
      end
!
!******************************************************************
!
      subroutine SRBF5all(RBF,order,krbf,mpn,mdp,minN,maxN,NF,nta)
      !5种的RBF（krbf）曲线（NF+1）0~dr°
      !1扰动重力2高程异常3空间异常4重力梯度5总垂线偏差
      !krbf-order次0径向多极子1Possion小波
      implicit none
      integer::order,krbf,minN,maxN,NF
      real*8::RBF(NF+1,5),mpn(maxN-minN+1,NF+1),mdp(maxN-minN+1,NF+1),nta
      integer::i,j,n,m,k
      real*8::onei,Bn,pi,RAD,CnmCalc
!---------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      m=order;RBF=0.d0
      onei=0.d0
      do i=1,NF+1!采用扰动位RFB系数归一化
        do n=minN,maxN
          k=n-minN+1
          if(m==0)then
              Bn=1.d0
          else
            if(krbf==0)Bn=CnmCalc(n,m)*dexp(-dlog(nta)*m)
            if(krbf==1)Bn=(2.d0*n+1.d0)*(-dble(n)*dlog(nta))**m
          endif
          RBF(i,1)=RBF(i,1)+(n+1.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
          RBF(i,2)=RBF(i,2)+Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBF(i,3)=RBF(i,3)+(n-1.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
          RBF(i,4)=RBF(i,4)+(n+1.d0)*(n+2.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
          RBF(i,5)=RBF(i,5)+Bn*dexp(dlog(nta)*n)*mdp(k,i)
          onei=onei+Bn*dexp(dlog(nta)*n)*mpn(k,i)
        enddo
      enddo
      RBF(1:NF+1,1:5)=RBF(1:NF+1,1:5)/onei
      end
!
!******************************************************************
!
      subroutine SRBF4one(RBF4,RBFn,order,krbf,mpn,mdp,minN,maxN,NF,nta)
      !RBF4-4种扰动场元的RBF（krbf）空域曲线（NF+1）0~dr°
      !RBF4-4种扰动场元的RBF谱域曲线
      !1扰动重力2高程异常3重力梯度4总垂线偏差
      !krbf-order次0径向多极子1Possion小波
      implicit none
      integer::order,krbf,fknd,minN,maxN,NF
      real*8::RBF4(NF+1,4),RBFn(maxN-minN+1,4),mpn(maxN-minN+1,NF+1),mdp(maxN-minN+1,NF+1),nta
      integer::i,j,n,m,k
      real*8::onei(4),Bn,pi,RAD,CnmCalc,tmp
!---------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      m=order;RBF4=0.d0;RBFn=0.d0
      tmp=0.d0;  onei=0.d0
      do i=1,NF+1
        do n=minN,maxN
          k=n-minN+1
          if(krbf==0)Bn=CnmCalc(n,m)*dexp(-dlog(nta)*m)
          if(krbf==1)then
             if(m>0) Bn=(2.d0*n+1.d0)*(-dble(n)*dlog(nta))**m
             if(m==0)Bn=1.d0
          endif
          RBF4(i,1)=RBF4(i,1)+(n+1.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
          RBF4(i,2)=RBF4(i,2)+Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBF4(i,3)=RBF4(i,3)+(n+1.d0)*(n+2.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
          RBF4(i,4)=RBF4(i,4)+Bn*dexp(dlog(nta)*n)*mdp(k,i)
          RBFn(k,1)=RBFn(k,1)+((n+1.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i))**2
          RBFn(k,2)=RBFn(k,2)+(Bn*dexp(dlog(nta)*n)*mpn(k,i))**2
          RBFn(k,3)=RBFn(k,3)+((n+1.d0)*(n+2.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i))**2
          RBFn(k,4)=RBFn(k,4)+(Bn*dexp(dlog(nta)*n)*mdp(k,i))**2
        enddo
      enddo
      do i=1,4
        tmp=maxval(RBF4(1:NF+1,i))-minval(RBF4(1:NF+1,i))
        RBF4(1:NF+1,i)=RBF4(1:NF+1,i)/tmp
        tmp=maxval(RBFn(1:maxN-minN+1,i))-minval(RBFn(1:maxN-minN+1,i))
        RBFn(1:maxN-minN+1,i)=RBFn(1:maxN-minN+1,i)/tmp
      enddo
      end
!
!******************************************************************
!
      subroutine SRBF4kgrav(RBF4,RBFn,order,knd,mpn,mdp,minN,maxN,NF,nta)
      !knd=0扰动重力1高程异常2总垂线偏差3重力梯度
      implicit none
      integer::order,knd,minN,maxN,NF
      real*8::RBF4(NF+1,4),RBFn(maxN-minN+1,4),mpn(maxN-minN+1,NF+1),mdp(maxN-minN+1,NF+1),nta
      integer::i,j,n,m,k
      real*8::onei(4),Bn(4),pi,RAD,CnmCalc,tmp
!---------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      m=order;RBF4=0.d0;RBFn=0.d0
      do i=1,NF+1
        do n=minN,maxN
          k=n-minN+1
          Bn(1)=1.d0
          Bn(2)=2.d0*n+1.d0
          Bn(3)=CnmCalc(n,m)*dexp(-dlog(nta)*m)
          Bn(4)=(2.d0*n+1.d0)*(-dble(n)*dlog(nta))**m!dexp(dlog((-dlog(nta)*n)*m))
          if(knd==0)then
             RBF4(i,1:4)=RBF4(i,1:4)+(n+1.d0)*Bn(1:4)*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
             RBFn(k,1:4)=RBFn(k,1:4)+((n+1.d0)*Bn(1:4)*dexp(dlog(nta)*(n-1.d0)))**2
          endif
          if(knd==1)then
             RBF4(i,1:4)=RBF4(i,1:4)+Bn(1:4)*dexp(dlog(nta)*n)*mpn(k,i)
             RBFn(k,1:4)=RBFn(k,1:4)+(Bn(1:4)*dexp(dlog(nta)*n))**2
          endif
          if(knd==2)then
             RBF4(i,1:4)=RBF4(i,1:4)+Bn(1:4)*dexp(dlog(nta)*(n-1.d0))*mdp(k,i)
             RBFn(k,1:4)=RBFn(k,1:4)+(Bn(1:4)*dexp(dlog(nta)*(n-1.d0)))**2
          endif
          if(knd==3)then
             RBF4(i,1:4)=RBF4(i,1:4)+(n+1.d0)*(n+2.d0)*Bn(1:4)*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
             RBFn(k,1:4)=RBFn(k,1:4)+((n+1.d0)*(n+2.d0)*Bn(1:4)*dexp(dlog(nta)*(n-1.d0)))**2 
          endif
        enddo
      enddo
      do i=1,4
        tmp=maxval(RBF4(1:NF+1,i))-minval(RBF4(1:NF+1,i))
        RBF4(1:NF+1,i)=RBF4(1:NF+1,i)/tmp
        tmp=maxval(RBFn(1:maxN-minN+1,i))-minval(RBFn(1:maxN-minN+1,i))
        RBFn(1:maxN-minN+1,i)=RBFn(1:maxN-minN+1,i)/tmp
      enddo
      end
!
!*********************************************************
!
      real*8 function CnmCalc(n,m) 
      implicit none
      integer :: i,k,n,m 
      real*8 :: Cnm(10),tmp
      if(m==0)then
        CnmCalc=1.d0;return
      endif
      if(m==1)then
        CnmCalc=n;return
      endif
      if(m==2)then
        CnmCalc=n*(n-1)/2;return
      endif
      Cnm(1)=n;Cnm(2)=n*(n-1)/2  !n=m+1->C(m+1,m)   
      do i=3,m
        Cnm(3)=Cnm(1)+Cnm(2)
        Cnm(1)=Cnm(2);Cnm(2)=Cnm(1)
      enddo
      CnmCalc=Cnm(3)
      return
      end
