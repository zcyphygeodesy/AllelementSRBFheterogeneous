      subroutine SRBFheterogeneous(observationfl,surfhgtgrdfl,para)
      implicit none
	character*800::observationfl,surfhgtgrdfl,outfl,nstr(80)
	character*800::rgafl,ksifl,grafl,grrfl,dftfl
	character*80000::line,line0,str,astr,ldstr
      character(len=15)::znm(900000)
	integer sn,kln,astat(20),i,j,k,NF,NF4,Kt,nn,mm,lvl,order,kobs,mn,ks,kp,ka,aa,k0
      real*8 rec(8000),BLH(3),rln(3),rlnk(3),gr,NFD(5),dln(2),GMr,rhd(4),rr,ab(6),cntwgh(6),kgm(5)
	real*8 para(20),GRS(6),hd(6),GM,ae,pi,RAD,wgh,r0,dpth,nta,br,sta(4),dout(200)!br-Bjerhammar
	real*8 val,unit(6),dlat,dr,tmp,mr,st(6,4),st0(6,4),bf(8),blat,sinf,tt,cosa,sina,lmt!第一行平行圈格网中心地心纬度°
	integer nlon,nlat,minN,maxN,typrw,wghrw,nk,krbf
	integer ki,kj,nd,mk,ni,nj,mthd,edgn,obsn(6)
	integer::status=0
	real*8,allocatable::mpn(:,:),mdp(:,:),mpn4(:,:),mdp4(:,:),obs(:,:,:),obsp(:,:),RBF4(:,:),RBFn(:,:)
      real*8,allocatable::BPB(:,:,:),BPL(:,:),BB(:),B2(:),xx(:),chs(:,:),sr(:),rlatlon(:,:),dl(:),lon(:,:)
	real*8,allocatable::APA(:,:),APL(:),hgt(:,:),rst(:,:),rst2(:,:),RBF(:,:)!目标场元
	integer,allocatable::nln(:),nrd(:,:),node(:),enode(:),gpnt(:,:)!格网、未知数序号,每个观测量有效节点序号
	real*8,allocatable::rga(:),ksi(:),gra(:),grr(:),vms(:),vmw(:,:)
!---------------------------------------------------------------------
      GRS(1)= 3.986004415d14; GRS(2)=6378136.3d0; GRS(3)=1.082636277388d-3
      GRS(4) = 7.292115d-5; GRS(5) = 1.d0/298.25641153d0
      GM=GRS(1)*1.d-7;ae=GRS(2)
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0;mr=36.d2/RAD
      unit(1)=1.d5;unit(2)=1.d0;unit(3)=1.d5;unit(4)=1.d9;unit(5)=mr;unit(6)=mr
      minN=nint(para(1));maxN=nint(para(2));order=nint(para(3));krbf=nint(para(4))
      lvl=nint(para(5));dpth=para(7)*1.d3;dr=para(6)*1.d3/ae/RAD!球面角距°
      typrw=nint(para(8));wghrw=nint(para(9))
      cntwgh=1.d0;if(para(10)>0.5)cntwgh(nint(para(10)))=para(11)
      cntwgh(6)=cntwgh(5);mthd=nint(para(12))
      !打开计算面大地高格网文件
      open(unit=8,file=surfhgtgrdfl,status="old",iostat=status)
      if(status/=0) goto 902
      read(8,'(a)') line0
      call PickReclong(line0,kln,rec,sn)
      if(sn<6)then
         close(8);goto 902
      endif
      hd(1:6)=rec(1:6)
	nlat=nint((hd(4)-hd(3))/hd(6))
	nlon=nint((hd(2)-hd(1))/hd(5))
	hd(5)=(hd(2)-hd(1))/dble(nlon)
	hd(6)=(hd(4)-hd(3))/dble(nlat)
 	allocate(hgt(nlat,nlon), stat=astat(1))
 	allocate(rst(nlat,nlon), stat=astat(2))
 	allocate(rst2(nlat,nlon), stat=astat(3))
	if (sum(astat(1:3)) /= 0) then
          close(8);goto 902
      endif
 	do i=1,nlat
	   read(8,*,end=905)(hgt(i,j),j=1,nlon)
      enddo
905   close(8)
      if(minN<order)minN=order; k0=0; k=0
      open(unit=8,file=observationfl,status="old",iostat=status)
      if(status/=0) goto 901
      read(8,'(a)') line
      do while(.not.eof(8))
        read(8,'(a)') line
        call PickRecstr(line,kln,nstr,sn)
        if(sn>3)then
          k=k+1;znm(k)=nstr(1)
        endif
        call PickReclong(line,kln,rec,sn)
        if(sn>3)then
            if(rec(2)<hd(2).and.rec(2)>hd(1).and.rec(3)<hd(4).and.rec(3)>hd(3))k0=k0+1
        endif
      enddo
      close(8);mn=nlat*nlon
      if(k0<6) goto 901
      NF=nint(dr*3600)!影响半径等分,间隔1″,NF+1→[0,dr]
      NF4=nint(dr*600)!间隔6″,NF4+1→[0,dr]
 	allocate(mpn(maxN-minN+1,NF+1), stat=astat(1))
 	allocate(mdp(maxN-minN+1,NF+1), stat=astat(2))
 	allocate(mpn4(maxN-minN+1,NF4+1), stat=astat(3))
 	allocate(mdp4(maxN-minN+1,NF4+1), stat=astat(4))
 	allocate(obs(3*k,6,6), stat=astat(5))!6种观测场元-球坐标，观测量，权值, 站点名对应的序号
 	allocate(chs(mn*2+3*k,6), stat=astat(6))!6种观测场元-球坐标，残差量，权值
 	allocate(obsp(mn+k*2,3), stat=astat(7))!全部观测点球坐标
	if (sum(astat(1:7)) /= 0) goto 901
      mn=k;obsn=0;r0=0.d0;k=0!重新读取观测文件，转换为观测点球面坐标，计算平均地心距r0
      open(unit=8,file=observationfl,status="old",iostat=status)
      read(8,'(a)') line
      do while(.not.eof(8))
        read(8,'(a)') line
        call PickReclong(line,kln,rec,sn)
        if(typrw>sn)goto 306
        if(rec(typrw)<-0.1.or.rec(typrw)>5.5)goto 306
        BLH(1)=rec(3);BLH(2)=rec(2);BLH(3)=rec(4);wgh=1.d0
        if(wghrw>0.and.wghrw<sn+1)wgh=rec(wghrw);if(wgh<-1.d-5)wgh=1.d0
        call BLH_RLAT(GRS,BLH,rln)
        k=k+1;r0=r0+rln(1)/ae;obsp(k,1:3)=rln(1:3)
        nk=nint(rec(typrw)+1);obsn(nk)=obsn(nk)+1  !!6种观测场元-球坐标，观测量，权值
        obs(obsn(nk),1:3,nk)=rln(1:3);obs(obsn(nk),4,nk)=rec(5);obs(obsn(nk),5,nk)=wgh;obs(obsn(nk),6,nk)=k!站点名对应的序号
306     continue
      enddo
      close(8);mn=k;r0=r0/dble(k)*ae!平均地心距
   !计算格网节点/未知数个数Kt
      !nln(nn)-平行圈方向格网数nn=maxi-mini+1
      !sr(nn)平行圈方向格网面积与赤道格网面积之差的百分比
      !dl(nn)平行圈方向经度间隔°Kt格网总点数-节点数、未知数
      !lon(nn,mm)格网中心经度,mm为平行圈方向最多格网数
      dlat=180.d0/lvl;nd=nint(dr/dlat+0.5d0)!dlat格网间隔,积分半径对应的格网数
      rhd(1:4)=hd(1:4);BLH(2)=(hd(1)+hd(2))/2.d0;BLH(3)=0.d0!!!!!!!目标格网范围用球坐标表示
      BLH(1)=hd(3);call BLH_RLAT(GRS,BLH,rln);rhd(3)=rln(2)
      BLH(1)=hd(4);call BLH_RLAT(GRS,BLH,rln);rhd(4)=rln(2)!!!!!!!!!
      nn=nint((rhd(4)-rhd(3))/dlat+0.5);mm=nint((rhd(2)-rhd(1))/dlat+0.5)!mm平行圈方向最大格网数
      allocate(nln(nn),sr(nn),dl(nn),nrd(nn,mm),gpnt(nn,mm),lon(nn,mm),rlatlon(2*(nn+mm),2),enode(2*(nn+mm)))
      call ReuterGrid(rhd,lvl,Kt,blat,nn,mm,nln,sr,dl,nrd,lon)!Kt节点数/未知数个数
      gpnt=0!计算格网中测点数，修正Reuter格网节点数Kt,序号nrd
      call Edgnode(enode,rlatlon,lvl,edgn,lon,blat,nln,gpnt,nn,mm)
      chs(1:edgn,1)=r0;chs(1:edgn,2:3)=rlatlon(1:edgn,1:2)
      chs(edgn+1:edgn+mn,1:3)=obsp(1:mn,1:3);ks=edgn+mn
      call AdjReuterGrd(chs(1:ks,1:3),ks,Kt,blat,lvl,nn,mm,nln,dl,lon,nrd,gpnt)
      allocate(RBF(NF+1,5),BPL(Kt,6),APL(Kt),BB(Kt),B2(Kt),xx(Kt),node(Kt),RBF4(NF4+1,4),RBFn(maxN-minN+1,4))
 	allocate(BPB(Kt,Kt,6), stat=astat(1))
 	allocate(APA(Kt,Kt), stat=astat(2))
	if (sum(astat(1:2)) /= 0) goto 903
   !计算minN~maxN阶NF+1组勒让德函数[minN~maxN]×[0,dr]
      call LegPn01(mpn,mdp,minN,maxN,NF,dr)
      call LegPn01(mpn4,mdp4,minN,maxN,NF4,dr)
   !由初始补偿深度dpth,计算SRFB曲线
      br=r0-dpth;nta=br/r0;rlnk(1)=br!初始补偿深度dpth和宽度参数nta
      call SRBF5all(RBF,order,krbf,mpn,mdp,minN,maxN,NF,nta)
      st0=0.d0;BPL=0.d0;BPB=0.d0!构造观测方程和法方程
      do kobs=0,5!!!!!***
        do k=1,obsn(kobs+1)!6种观测场元-球坐标，观测量，权值obs(k,1:5,kobs+1)
          rln(1:3)=obs(k,1:3,kobs+1);wgh=obs(k,5,kobs+1);val=obs(k,4,kobs+1)/unit(kobs+1);rr=rln(1)
          if(kobs==0.or.kobs==2)GMr=GM/rr/rr
          if(kobs==3)GMr=GM/rr/rr/rr
          if(kobs==1.or.kobs==4)then
            call RLAT_BLH(GRS,rln,BLH);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
            if(kobs==1)GMr=GM/rr/gr
            if(kobs==4.or.kobs==5)GMr=GM/rr/rr/gr
          endif
          call RtGridij(rln,ki,kj,blat,lvl,nn,mm,nln,dl,lon)!计算rln在格网中的序号ki>0,kj>0
          mk=0;BB=0.d0;node=0!BB-RBF系数,node(Kt)观测量作用半径范围内有效节点序号
          do i=ki-nd,ki+nd!增加实现节点序号与格网行列号对应，即节点序号二维整数数组
            if(i<1.or.i>nn)goto 1001!lon(i,j),blat第一行平行圈格网中心地心纬度°
            rlnk(2)=blat+(i-1.d0)*dlat
            do j=kj-nd,kj+nd
              if(j<1.or.j>nln(i)) goto 1002 
              if(nrd(i,j)<1) goto 1002 !nrd(i,j)未知数序号
              rlnk(3)=lon(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
              if(dln(2)>dr)goto 1002
              if(kobs<5)call RBFvalue(RBF(:,kobs+1),NF,dr,dln(2),tmp)!由球面角距°dln(2)内插RBF值tmp
              if(kobs==5)call RBFvalue(RBF(:,5),NF,dr,dln(2),tmp)!由球面角距°dln(2)内插RBF值tmp
              mk=mk+1;node(mk)=nrd(i,j)
              BB(nrd(i,j))=GMr*tmp*dexp((minN-1.0)*dlog(r0/rr))
              if(kobs>3)then
                tt=dcos(dln(2)*RAD);sinf=dsqrt(1.d0-tt**2);if(sinf<1.d-12)sinf=1.d-12
                if(kobs==4)then
	            cosa=(dcos(rln(2)*RAD)*dsin(rlnk(2)*RAD)-dsin(rln(2)*RAD)*dcos(rlnk(2)*RAD)*dcos(rlnk(3)*RAD-rln(3)*RAD))/sinf
                  BB(nrd(i,j))= BB(nrd(i,j))*cosa
                else
	            sina=dcos(rlnk(2)*RAD)*dsin(rlnk(3)*RAD-rln(3)*RAD)/sinf
                  BB(nrd(i,j))= BB(nrd(i,j))*sina
                endif
              endif
1002          continue
            enddo
1001        continue    
          enddo
          do i=1,mk
            ki=node(i)
            BPL(ki,kobs+1)=BPL(ki,kobs+1)+BB(ki)*val*wgh
            do j=1,i
              kj=node(j)
              BPB(ki,kj,kobs+1)=BPB(ki,kj,kobs+1)+BB(ki)*BB(kj)*wgh
	      enddo
          enddo
        enddo
        if(obsn(kobs+1)>0)call Stat1d(obs(1:obsn(kobs+1),4,kobs+1),obsn(kobs+1),st0(kobs+1,1:4))
      enddo!!!!***
      APA=0.d0;APL=0.d0;ab=0.d0
      do k=1,6!计算法方程对角线标准差
        if(obsn(k)<1)goto 3131
        do i=1,Kt
           B2(i)=BPB(i,i,k)
        enddo
        ab(k)=maxval(B2(1:Kt))-minval(B2(1:Kt))
3131    continue
      enddo
      do k=1,6
        if(obsn(k)<1)goto 3232
        tmp=cntwgh(k)/ab(k)!法方程对角线标准差法组合,调控场元贡献
        do i=1,Kt
          APL(i)=APL(i)+tmp*BPL(i,k)
          do j=1,i
            APA(i,j)=APA(i,j)+tmp*BPB(i,j,k)
	    enddo
        enddo
3232    continue
      enddo
      tmp=0.d0
	do i=1,Kt
	   do j=1,i-1
	      APA(j,i)=APA(i,j)
         enddo
         tmp=tmp+APA(i,i)**2/dble(Kt)
      enddo
      tmp=dsqrt(tmp)
      !以Reuter格网四周节点未知数为零组成观测方程，抑制边缘效应。
      !节点序号数组enode
      do i=1,edgn!edgn-Reuter格网四周节点数
         ki=enode(i); APA(ki,ki)=APA(ki,ki)+tmp/dsqrt(dble(mn))
      enddo
	do i=1,Kt
         APA(i,i)=APA(i,i)+tmp*1.d-4
      enddo
5001  xx=0.d0
      !mthd=1 LU分解,2 Cholesky分解,3 最小二乘QR分解,4最小范数奇异值分解,5-岭估计
      if(mthd<5)call Equsolve(APA,xx,Kt,APL,mthd,bf)
      if(mthd==5)call RidgeEstimate(APA,xx,Kt,APL)
      chs=0.d0!计算残差并统计
      do kobs=0,5!!!!!***
        do k=1,obsn(kobs+1)!6种观测场元-球坐标，观测量，权值obs(k,1:5,kobs+1)
          rln(1:3)=obs(k,1:3,kobs+1);rr=rln(1)
          if(kobs==0.or.kobs==2)GMr=GM/rr/rr
          if(kobs==3)GMr=GM/rr/rr/rr
          if(kobs==1.or.kobs==4)then
            call RLAT_BLH(GRS,rln,BLH);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
            if(kobs==1)GMr=GM/rr/gr
            if(kobs==4.or.kobs==5)GMr=GM/rr/rr/gr
          endif
          call RtGridij(rln,ki,kj,blat,lvl,nn,mm,nln,dl,lon)!计算rln在格网中的序号ki>0,kj>0
          mk=0;BB=0.d0;node=0!s-RBF系数,node(Kt)观测量作用半径范围内有效节点序号
          do i=ki-nd,ki+nd!增加实现节点序号与格网行列号对应，即节点序号二维整数数组，0表示
            if(i<1.or.i>nn)goto 2001!lon(i,j),blat第一行平行圈格网中心地心纬度°
            rlnk(2)=blat+(i-1.d0)*dlat
            do j=kj-nd,kj+nd
              if(j<1.or.j>nln(i)) goto 2002 
              if(nrd(i,j)<1) goto 2002 !nrd(i,j)未知数序号
              rlnk(3)=lon(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
              if(dln(2)>dr)goto 2002
              if(kobs<5)call RBFvalue(RBF(:,kobs+1),NF,dr,dln(2),tmp)!由球面角距°dln(2)内插RBF值tmp
              if(kobs==5)call RBFvalue(RBF(:,5),NF,dr,dln(2),tmp)!由球面角距°dln(2)内插RBF值tmp
              mk=mk+1;node(mk)=nrd(i,j)
              BB(nrd(i,j))=GMr*tmp*dexp((minN-1.0)*dlog(r0/rr))
              if(kobs>3)then
                tt=dcos(dln(2)*RAD);sinf=dsqrt(1.d0-tt**2);if(sinf<1.d-12)sinf=1.d-12
                if(kobs==4)then
	            cosa=(dcos(rln(2)*RAD)*dsin(rlnk(2)*RAD)-dsin(rln(2)*RAD)*dcos(rlnk(2)*RAD)*dcos(rlnk(3)*RAD-rln(3)*RAD))/sinf
                  BB(nrd(i,j))= BB(nrd(i,j))*cosa
                else
	            sina=dcos(rlnk(2)*RAD)*dsin(rlnk(3)*RAD-rln(3)*RAD)/sinf
                  BB(nrd(i,j))= BB(nrd(i,j))*sina
                endif
              endif
2002          continue
            enddo
2001        continue    
          enddo
          val=0.d0
          do i=1,mk
            ki=node(i)
            val=val+BB(ki)*xx(ki)
          enddo
          chs(k,kobs+1)=obs(k,4,kobs+1)-val*unit(kobs+1)
        enddo
        if(obsn(kobs+1)>0)call Stat1d(chs(1:obsn(kobs+1),kobs+1),obsn(kobs+1),st(kobs+1,1:4))
      enddo!!!!***
      open(unit=40,file="residuals.txt",status="replace")
      dout(41:100)=0.d0
      do kobs=0,5
        dout(40+kobs*10+1)=kobs
	  if(obsn(kobs+1)>0)then
          write(40,'(I3,2F10.4,2F11.4,a,2F10.4,40F11.4)')kobs,(st0(kobs+1,i),i=1,4),"  residuals:",(st(kobs+1,i),i=1,4)
          dout(40+kobs*10+2:40+kobs*10+5)=st0(kobs+1,1:4);dout(40+kobs*10+6:40+kobs*10+9)=st(kobs+1,1:4)
        endif
      enddo
      ki=0
      do kobs=0,5
        do k=1,obsn(kobs+1)!6种观测场元-球坐标，观测量，权值obs(k,1:5,kobs+1)
          rln(1:3)=obs(k,1:3,kobs+1);call RLAT_BLH(GRS,rln,BLH);ki=ki+1;kp=nint(obs(k,6,kobs+1))
	    write(40,'(a8,2F12.5,F11.3,2F12.4,I4,F8.3)')trim(znm(kp)),BLH(2),BLH(1),BLH(3),chs(k,kobs+1),
     *          obs(k,4,kobs+1),kobs,obs(k,5,kobs+1)
        enddo
      enddo
      close(40)
      allocate(rga(nlon),ksi(nlon),gra(nlon),grr(nlon),vms(nlon),vmw(nlat,nlon))
	write(rgafl,*) "SRBFhetero.rga"
	write(ksifl,*) "SRBFhetero.ksi"
	write(grafl,*) "SRBFhetero.gra"
	write(grrfl,*) "SRBFhetero.grr"
	write(dftfl,*) "SRBFhetero.dft"
      aa=0;open(unit=10,file="SRBFhetero.txt",status="replace")
      write(ldstr,*)"ID lon lat ellipshgt gravity disturbance(mGal) height anomaly(m) gravity anomaly(mGal),"
      ldstr=trim(ldstr)//"gravity gradient(E) vertical deflection(S,W)"
      write(10,'(a)')trim(ldstr)
      open(unit=20,file=rgafl,status="replace");write(20,'(a)')line0
      open(unit=22,file=ksifl,status="replace");write(22,'(a)')line0
      open(unit=24,file=grafl,status="replace");write(24,'(a)')line0
      open(unit=26,file=grrfl,status="replace");write(26,'(a)')line0
      open(unit=28,file=dftfl,status="replace");write(28,'(a)')line0
 	do ni=1,nlat
        BLH(1)=hd(3)+(real(ni)-0.5d0)*hd(6)
        do nj=1,nlon
	    BLH(2)=hd(1)+(real(nj)-0.5d0)*hd(5);aa=aa+1
          BLH(3)=hgt(ni,nj);call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
          kgm(1)=GM/rr/rr;kgm(3)=kgm(1); kgm(4)=kgm(1)/rr
          call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
          kgm(2)=GM/rr/gr;kgm(5)=kgm(2)/rr
          call RtGridij(rln,ki,kj,blat,lvl,nn,mm,nln,dl,lon)!计算rln在格网中的序号ki>0,kj>0
          do kp=1,5!0扰动重力，1高程异常，2空间异常，3扰动重力梯度，4垂线偏差
            mk=0;BB=0.d0;B2=0.d0;node=0!s-RBF系数,node(Kt)观测量作用半径范围内有效节点序号
            do i=ki-nd,ki+nd!增加实现节点序号与格网行列号对应，即节点序号二维整数数组，0表示
              if(i<1.or.i>nn)goto 3001!lon(i,j),blat第一行平行圈格网中心地心纬度°
              rlnk(2)=blat+(i-1.d0)*dlat
              do j=kj-nd,kj+nd
                if(j<1.or.j>nln(i)) goto 3002 
                if(nrd(i,j)<1) goto 3002 !nrd(i,j)未知数序号
                rlnk(3)=lon(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
                if(dln(2)>dr)goto 3002
                call RBFvalue(RBF(:,kp),NF,dr,dln(2),tmp)!由球面角距°dln(2)内插RBF值val
                mk=mk+1;node(mk)=nrd(i,j)
                BB(nrd(i,j))=kgm(kp)*tmp*dexp((minN-1.0)*dlog(r0/rr))
                if(kp==5)then
                  tt=dcos(dln(2)*RAD);sinf=dsqrt(1.d0-tt**2);if(sinf<1.d-12)sinf=1.d-12
	            cosa=(dcos(rln(2)*RAD)*dsin(rlnk(2)*RAD)-dsin(rln(2)*RAD)*dcos(rlnk(2)*RAD)*dcos(rlnk(3)*RAD-rln(3)*RAD))/sinf
	            sina=dcos(rlnk(2)*RAD)*dsin(rlnk(3)*RAD-rln(3)*RAD)/sinf
                  BB(nrd(i,j))= BB(nrd(i,j))*cosa
                  B2(nrd(i,j))= BB(nrd(i,j))*sina
                endif
3002            continue
              enddo
3001          continue    
            enddo
            val=0.d0;tmp=0.d0
            do i=1,mk
              ka=node(i)
              val=val+BB(ka)*xx(ka)
              if(kp==5)tmp=tmp+B2(ka)*xx(ka)
            enddo
            if(kp==1)rga(nj)=val*unit(kp);
            if(kp==2)ksi(nj)=val*unit(kp);
            if(kp==3)gra(nj)=val*unit(kp);
            if(kp==4)grr(nj)=val*unit(kp);
            if(kp==5)then
                vms(nj)=val*unit(kp);vmw(ni,nj)=tmp*unit(kp)
            endif
          enddo!kp
          write(10,'(I7,2F11.5,F10.3,15F11.4)')aa,BLH(2),BLH(1),BLH(3),rga(nj),ksi(nj),gra(nj),grr(nj),vms(nj),vmw(ni,nj)
        enddo
        write(20,'(15F12.4)')(rga(nj),nj=1,nlon)
        write(22,'(15F12.4)')(ksi(nj),nj=1,nlon)
        write(24,'(15F12.4)')(gra(nj),nj=1,nlon)
        write(26,'(15F12.4)')(grr(nj),nj=1,nlon)
        write(28,'(15F12.4)')(vms(nj),nj=1,nlon)
      enddo
      do ni=1,nlat
	  write(28,'(15F12.4)')(vmw(ni,nj),nj=1,nlon)
      enddo
      close(10);close(20);close(22);close(24);close(26);close(28)
      deallocate(ksi,gra,rga,grr,vms,vmw)
      call SRBF4one(RBF4,RBFn,order,krbf,mpn4,mdp4,minN,maxN,NF4,nta)
      open(unit=10,file="SRBFspc.txt",status="replace")
	write(10,'(2I3,2I6,F8.2,a)')krbf,order,minN,maxN,dpth*1.d-3,
     *  " gravity disturbance, height anomaly, gravity gradient, vertical deflection" 
      tmp=dr/dble(NF4)*RAD*ae*1.d-3
 	do i=1,NF4
	   write(10,'(F12.3,8F13.5)')-(NF4-i+1.0)*tmp,(RBF4(NF4-i+2,j),j=1,4)
      enddo
 	do i=1,NF4+1
	   write(10,'(F12.3,8F13.5)')(i-1.0)*tmp,(RBF4(i,j),j=1,4)
      enddo
      close(10)
      open(unit=10,file="SRBFdgr.txt",status="replace")
	write(10,'(2I3,2I6,F8.2,a)')krbf,order,minN,maxN,dpth*1.d-3,
     * " gravity disturbance, height anomaly, gravity gradient, vertical deflection" 
 	do i=minN,maxN
	   write(10,'(I8,8F13.5)')i,(RBFn(i-minN+1,j),j=1,4)
      enddo
      close(10)
      open(unit=10,file="SRBFcenter.txt",status="replace")
	write(10,'(4I6,F8.3)')lvl,Kt,nn,nln(1),dlat*60.d0
      rln(1)=rr;k=0
      do i=1,nn!增加实现节点序号与格网行列号对应，即节点序号二维整数数组
          rln(2)=blat+(i-0.5d0)*dlat
          do j=1,nln(i)
            if(gpnt(i,j)>0)then
              k=k+1;rln(3)=lon(i,j);call RLAT_BLH(GRS,rln,BLH)
	        write(10,'(I6,2F12.6,2F8.3)')k,BLH(2),BLH(1),sr(i),dl(i)*60.d0
            endif
          enddo
      enddo
      close(10)
909   continue
904   deallocate(BPB,APA)
903   deallocate(nln,sr,dl,nrd,lon,node,enode,RBF,RBF4,RBFn,gpnt,rlatlon)
      deallocate(mpn,mdp,mpn4,mdp4,obs,obsp,BB,B2,BPL,APL,xx,chs)
901   deallocate(hgt,rst,rst2)
902	continue
101   format(a,40F12.4)
      end
!
!******************************************************************
!
      subroutine drln(rln1,rln2,dln)
      !由两点球坐标计算距离与夹角°
      implicit none
      real*8::rln1(3),rln2(3),dln(2),XYZ1(3),XYZ2(3),L2
      integer::i,j,n,m,k
      real*8::onei,Bn,pi,RAD,CnmCalc
!---------------------------------------------------------------
      RAD=datan(1.d0)/45.d0
      call RLAT_XYZ(rln1,XYZ1)
      call RLAT_XYZ(rln2,XYZ2)
      L2=(XYZ1(1)-XYZ2(1))**2+(XYZ1(2)-XYZ2(2))**2+(XYZ1(3)-XYZ2(3))**2
      dln(1)=dsqrt(L2)
      dln(2)=dabs(dacos((rln1(1)**2+rln2(1)**2-L2)/2.d0/rln1(1)/rln2(1)))/RAD
      end
!
!******************************************************************************
!
      subroutine LegPn01(mpn,mdp,minN,maxN,NF,dr)
      !计算minN~maxN阶NF+1组勒让德函数[minN~maxN]×[0,dr]
	implicit none
	integer::minN,maxN,NF,nn,i,k,n
	real*8::dr,t,dt,pi,RAD
	real*8::mpn(maxN-minN+1,NF+1),mdp(maxN-minN+1,NF+1)
	real*8,allocatable::pn(:),dp1(:)
!---------------------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      dt=dr/dble(NF)*RAD
      allocate(pn(maxN+2),dp1(maxN+2))
      do k=1,NF+1
        t=dcos(dble(k-1)*dt)
        call LegPn_dt1(pn,dp1,maxN+1,t)
        do n=minN,maxN
          i=n-minN+1
          mpn(i,k)=pn(n+1)
          mdp(i,k)=dp1(n+1)
        enddo
      enddo
      deallocate(pn,dp1)
      end
!
!******************************************************************************
!
      subroutine Stat1d(dt,nn,rst)
      implicit none
	integer::nn,i,kk
      real*8::dt(nn),rst(4),pv,std,maxt,mint
!---------------------------------------------------------------------
	pv=0.d0;std=0.d0;maxt=-9.d28;mint=9.d28
      kk=0
      do i=1,nn
        if(dt(i)>9900.d0)goto 1001
        kk=kk+1;pv=pv+dt(i)
	  if(maxt<dt(i))maxt=dt(i)
        if(mint>dt(i))mint=dt(i)
1001    continue
      enddo
      pv=pv/dble(kk)
      do i=1,nn
        if(dt(i)>9900.d0)goto 1002
        std=std+(dt(i)-pv)**2
1002    continue
      enddo
      std=dsqrt(std/dble(kk))
      rst(1)=pv;rst(2)=std;rst(3)=mint;rst(4)=maxt
      end
