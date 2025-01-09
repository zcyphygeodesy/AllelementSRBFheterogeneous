      SUBROUTINE rnl_XYZd(rln,drln,dxyz)
      implicit none
      real*8::XYZ(3),rln(3),drln(3),dxyz(3),RAD
      real*8::sinu,cosu,sinl,cosl,rr,dr,du,dl
!-----------------------------------------------------------------------
	RAD=datan(1.d0)/45.d0;rr=rln(3) !令rr=1
      sinu=dsin(rln(2)*RAD);cosu=dcos(rln(2)*RAD)
      sinl=dsin(rln(3)*RAD);cosl=dcos(rln(3)*RAD)
      dl=drln(1);du=drln(2);dr=drln(3)
      dxyz(1)=dr*cosu*cosl-du*sinu*cosl-dl*sinl
      dxyz(2)=dr*cosu*sinl-du*sinu*sinl+dl*cosl
      dxyz(3)=dr*sinu+du*cosu
      return
      end
!
!********************************************************************************
!
      SUBROUTINE  BLH_XYZ(GRS,BLH,XYZ)
!输入BLH－纬度、经度(度)、大地高
      implicit none
      real*8::GRS(6),a,f,BLH(3),XYZ(3)
      real*8 slat,clat,slon,clon,v,h,e2,RAD
!-----------------------------------------------------------------------
      a=GRS(2); f=GRS(5) !1/f
	e2=(2.d0-f)*f
	h=BLH(3)
	RAD=datan(1.d0)/45.d0
      slat = dsin(BLH(1)*RAD)
      clat = dcos(BLH(1)*RAD)
      slon = dsin(BLH(2)*RAD)
      clon = dcos(BLH(2)*RAD)
      v = a / dsqrt(1.d0-e2*slat*slat)

      XYZ(1) = (v+h) * clat * clon
      XYZ(2) = (v+h) * clat * slon
      XYZ(3) = (v*(1.d0-e2)+h) * slat
      return
      end
!
!**************************************************************************
!
      SUBROUTINE BLH_RLAT(GRS,BLH,RLAT)
!2005年3月27日,章传银
       IMPLICIT REAL*8(A-H,O-Z)
	 real*8 GRS(6),a,f,BLH(3),RLAT(3)
       real*8 N,K,flat,RLATB
!----------------------------------------------------------
       a=GRS(2); f=GRS(5) !1/f
       AE=a
	 E2=(2.d0-f)*f
	 pi=datan(1.d0)*4.d0
	 RAD=datan(1.d0)/45.d0
	 FLAT=BLH(1)
	 FLON=BLH(2)
	 HT=BLH(3)

       FLATR=FLAT*RAD
       FLONR=FLON*RAD
       T1=DSIN(FLATR)**2
       N=AE/DSQRT(1.d0-E2*T1)
       T2=(N+HT)*DCOS(FLATR)
       X=T2*DCOS(FLONR)
       Y=T2*DSIN(FLONR)
       Z=(N*(1.d0-E2)+HT)*DSIN(FLATR)
       N=AE/DSQRT(1.d0-E2*T1)
! COMPUTE THE GEOCENTRIC RADIUS
       RE=DSQRT(X**2+Y**2+Z**2)
! COMPUTE THE GEOCENTRIC LATITUDE
       RLATB=DATAN(Z/DSQRT(X**2+Y**2))
! COMPUTE NORMAL GRAVITY:UNITS ARE M/SEC**2
!      GR=GEQT*(1.+K*T1)/DSQRT(1.-E2*T1)
	   RLAT(1)=RE
	   RLAT(2)=RLATB/RAD
	   RLAT(3)=FLON
       RETURN
      END
!
!*********************************************************************** xyz2ell
!
      subroutine XYZ_BLH(GRS,XYZ,BLH)
      implicit none
      real*8 GRS(6),a,f,XYZ(3),BLH(3)
      real*8 ae, e2, lat, lon, h, x, y, z
      real*8 pi, dlat, dh, lat0, h0, elat, eh, v, p
      data   elat/1.d-12/, eh/1.d-5/
!-----------------------------------------------------------------------
      a=GRS(2); f=GRS(5) !1/f
      pi=datan(1.d0)*4.d0
	  ae=a
	  e2=(2.d0-f)*f
	  x=XYZ(1)
	  y=XYZ(2)
	  z=XYZ(3)
      p   = sqrt(x*x + y*y)
      lat = datan2(z, p/(1.d0-e2))
      h   = 0.d0

!  Iterate until height and latitude converge

 10   continue
      lat0 = lat
      h0   = h
      v    = ae / sqrt( 1.d0-e2*sin(lat)*sin(lat) )
      h    = p / cos(lat) - v
      lat  = atan2( z, p*(1.d0-e2*v/(v+h)) )
      dlat = abs(lat - lat0)
      dh   = abs(h - h0)
      if (dlat.gt.elat .or. dh.gt.eh) goto 10

!  Compute longitude

      lon = datan2(y,x)
	  lon=dmod(lon,2.d0*pi)
      if (lon.lt.0.d0) lon = lon + 2.d0*pi
	  lat=lat+pi/2.d0
	  lat=dmod(lat,pi)
	  lat=lat-pi/2.d0
	  BLH(1)=lat*180.d0/pi
	  BLH(2)=lon*180.d0/pi
	  BLH(3)=h

      return
      end
!
!********************************************************************************
!
      SUBROUTINE RLAT_XYZ(RLAT,XYZ)
      implicit none
      real*8::XYZ(3),RLAT(3),RAD
	RAD=datan(1.d0)/45.d0
	XYZ(1)=dcos(RLAT(2)*RAD)*dcos(RLAT(3)*RAD)
	XYZ(2)=dcos(RLAT(2)*RAD)*dsin(RLAT(3)*RAD)
	XYZ(3)=dsin(RLAT(2)*RAD)
	XYZ=XYZ*RLAT(1)
      return
      end
!
!********************************************************************************
!
      SUBROUTINE XYZ_RLAT(GRS,XYZ,RLAT)
      implicit none
      real*8::GRS(6),XYZ(3),BLH(3),RLAT(3)
	  call XYZ_BLH(GRS,XYZ,BLH)
	  call BLH_RLAT(GRS,BLH,RLAT)
      return
      end
!
!********************************************************************************
!
      SUBROUTINE RLAT_BLH(GRS,RLAT,BLH)
      implicit none
      real*8::GRS(6),XYZ(3),BLH(3),RLAT(3)
	  call RLAT_XYZ(RLAT,XYZ)
	  call XYZ_BLH(GRS,XYZ,BLH)
      return
      end
