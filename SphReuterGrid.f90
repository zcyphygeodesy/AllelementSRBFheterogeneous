      subroutine ReuterGrid(hd,KK,Pn,blat,nn,mm,nln,ss,dl,nrd,lon)
      !nln(nn)-平行圈方向格网数nn=maxi-mini+1
      !ss(nn)平行圈方向格网面积与赤道格网面积之差的百分比
      !dl(nn)平行圈方向经度间隔°,pn格网总点数
      !lon(nn,mm)格网中心经度,mm为平行圈方向最多格网数
      !nrd(nn,mm)节点序号
      implicit none
      integer::i,j,ki,kj,KK,nn,mm,mini,maxi,minj,maxj
      integer::nln(nn),nrd(nn,mm),ML,Pn,k
	real*8::hd(4),ss(nn),dl(nn),lon(nn,mm),blat!第一行平行圈格网中心地心纬度°
 	real*8::pi,RAD,dlat,dlon,lat
!---------------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0;k=0
      nln=0;nrd=-1;dlat=pi/KK!格网纬度间隔
      mini=nint((hd(3)+90.d0)*RAD/dlat+0.5d0)
      maxi=nint((hd(4)+90.d0)*RAD/dlat-0.5d0)
      do i=mini,maxi!南北极需要单独处理
        ki=i-mini+1
        lat=-pi/2.d0+(dble(i)-0.5d0)*dlat
        if(ki==1)blat=lat/RAD
        ML=nint(2.d0*pi*dcos(lat)/dlat)!平行圈方向格网数近似
        dlon=2.d0*pi/dble(ML)!经度间隔
        !计算hd(1),hd(2)的格网序号minj,maxj
        minj=nint(hd(1)*RAD/dlon+0.5d0)
        maxj=nint(hd(2)*RAD/dlon-0.5d0)
        nln(ki)=maxj-minj+1!平行圈方向格网数
        dl(ki)=dlon/RAD
        ss(ki)=(dlon*dcos(lat)/dlat-1.d0)*1.d2!dl平行圈方向格网长度
        do j=minj,maxj
          kj=j-minj+1;k=k+1;nrd(ki,kj)=k
          lon(ki,kj)=(dble(j)-0.5d0)*dlon/RAD
        enddo
      enddo
      Pn=k
902   continue
      end
!
!************************************************************************************
!
      subroutine RtGridij(rln,ki,kj,blat,KK,nn,mm,nln,dl,lon)
      !已知计算点rln(3)，返回其落在Reuter格网中行列号ki,kj
      !结果用于缩小搜索范围，无需精准
      !nln(nn)-平行圈方向格网数
      !dl(nn)平行圈方向经度间隔°,pn格网总点数
      !当ki=0,kj=0时，表示rln(3)不在格网内
      implicit none
      integer::ki,kj,KK,nn,mm
      integer::nln(nn)
	real*8::rln(3),blat,dl(nn),lon(nn,mm),dlat,dlon
!---------------------------------------------------------------------
      ki=0;kj=0
      dlat=180.d0/KK!格网纬度间隔°
      ki=nint((rln(2)-blat)/dlat+0.5d0)
      if(ki<1.or.ki>nn)then
        ki=0;kj=0;return
      endif
      if(nln(ki)<1)then
        ki=0;kj=0;return
      endif
      dlon=dl(ki)!平行圈方向经度间隔°
      kj=nint((rln(3)-lon(ki,1))/dlon+0.5001d0)
      if(kj<1.or.kj>nln(ki))then
        ki=0;kj=0
      endif
      return
      end
!
!************************************************************************************
!
      subroutine AdjReuterGrd(pos,mn,Mt,blat,KK,nn,mm,nln,dl,lon,nrd,gpnt)
      implicit none
      integer::mn,Mt,KK,nn,mm,npv,k,ki,kj,i,j
      integer::nln(nn),nrd(nn,mm),gpnt(nn,mm)
	real*8::pos(mn,3),rln(3),blat,dl(nn),lon(nn,mm),dlat,dlon
!---------------------------------------------------------------------
      nrd=0;dlat=180.d0/KK!格网纬度间隔°kk单数检查!!!!!!!!!!!!
      do k=1,mn
         rln=pos(k,1:3)
         call RtGridij(rln,ki,kj,blat,KK,nn,mm,nln,dl,lon)
         if(ki>0.and.ki<=nn)then
           if(kj>0.and.kj<=nln(ki))gpnt(ki,kj)=gpnt(ki,kj)+1
         endif
      enddo
      k=0
      do i=1,nn
        do j=1,nln(i)
          if(gpnt(i,j)<1)then
              nrd(i,j)=0
          else
              k=k+1;nrd(i,j)=k
          endif
        enddo 
      enddo
      Mt=k
      return
      end
!
!******************************************************************
!
      subroutine Edgnode(enode,rlatlon,lvl,edgn,lon,blat,nln,gpnt,nn,mm)
      !返回Reuter格网周边节点序号enode(edgn)和球坐标rlatlon(edgn,2)
      implicit none
      integer::lvl,edgn,nk,nn,mm,npv,k,ki,kj,i,j
      integer::nln(nn),gpnt(nn,mm),enode(2*(nn+mm))
	real*8::blat,rlatlon(2*(nn+mm),2),lon(nn,mm),dlat
!---------------------------------------------------------------------
      nk=0;dlat=180.d0/lvl!格网纬度间隔nln(nn)
      do j=1,nln(1)
         nk=nk+1;gpnt(1,j)=1;enode(nk)=nk
         rlatlon(nk,2)=lon(1,j);rlatlon(nk,1)=blat
      enddo
      do j=1,nln(nn)
         nk=nk+1;gpnt(nn,j)=1;enode(nk)=nk
         rlatlon(nk,2)=lon(nn,j);rlatlon(nk,1)=blat+(nn-1.d0)*dlat
      enddo
      do i=2,nn-1
          nk=nk+1;gpnt(i,1)=1;enode(nk)=nk
          rlatlon(nk,2)=lon(i,1);rlatlon(nk,1)=blat+(i-1.d0)*dlat
          nk=nk+1;gpnt(i,nln(i))=1;enode(nk)=nk
          rlatlon(nk,2)=lon(i,nln(i));rlatlon(nk,1)=blat+(i-1.d0)*dlat
      enddo
      edgn=nk
      return
      end
