      subroutine PickRecord(str0,kln,rec,nn)
      implicit none
	character*800::str,str0,b
	integer::kln,nn,i,j,k,m,err
	real*8::rec(800)
!---------------------------------------------------------------------
      str=trim(AdjustL(str0))
      kln=len_trim(str)
      if(kln<2)then
         rec=0.d0;nn=0;return
      endif
      k=1; m=1
      do i=1,kln
         if(str(i:i)==','.or.str(i:i)==' '.or.iachar(str(i:i))==9.or.iachar(str(i:i))==11) then
            if(i-1>=m) then
               do j=1,kln
                  b(j:j)=' '
               enddo
               b(1:i-m+1)=str(m:i-1)
               read(b,*,iostat=err) rec(k)
               if(abs(err)>0)rec(k)=-1.d0
               k=k+1
            endif
            m=i+1
         endif
      enddo
      if(i-m<1)then
         k=k-1; goto 1101
      endif   
      b(1:i-m+1)=str(m:kln)
      read(b,*,iostat=err) rec(k)
      if(abs(err)>0)rec(k)=-1.d0
1101  nn=k
      end
!
!********************************************************
!
      subroutine PickRecstr(line0,kln,str,nn)
      implicit none
	character*800::line,line0,b,str(80)
	integer::kln,nn,i,j,k,m,err,nval,sgn
!---------------------------------------------------------------------
      line=trim(AdjustL(line0))
      kln=len_trim(line)
      if(kln<2)then
         nn=0;return
      endif
      k=1; m=1
      do i=1,kln
         nval=iachar(line(i:i))
         if(nval==32.or.nval==44.or.nval==59.or.nval==9)line(i:i)=' '
         nval=iachar(line(i:i));sgn=1
         if(nval==32)sgn=0
         if(nval>42.and.nval<47)sgn=0
         if(nval>47.and.nval<58)sgn=0
         if(nval>64.and.nval<91)sgn=0
         if(nval>96.and.nval<123)sgn=0
         if(nval==95)sgn=0
         if(sgn>0)line(i:i)=''
         if(line(i:i)==' ') then
            if(i-1>=m) then
               do j=1,kln
                  b(j:j)=' '
               enddo
               b(1:i-m+1)=line(m:i-1)
               read(b,*,iostat=err)str(k)
               if(abs(err)>0)str(k)=''
               k=k+1
            endif
            m=i+1
         endif
      enddo
      if(i-m<1)then
         k=k-1; goto 1101
      endif   
      b(1:i-m+1)=line(m:kln)
      read(b,*,iostat=err)str(k)
      if(abs(err)>0)str(k)=''
1101  nn=k
      do i=1,nn
        str(i)=trim(str(i))
      enddo
      end
!
!*******************************************************************************************
!
      subroutine PickReclong(str0,kln,rec,nn)
      implicit none
	character*80000::str,str0,b
	integer::kln,nn,i,j,k,m,err
	real*8::rec(8000)
!---------------------------------------------------------------------
      str=trim(AdjustL(str0))
      kln=len_trim(str)
      if(kln<2)then
         rec=0.d0;nn=0;return
      endif
      k=1; m=1
      do i=1,kln
         if(str(i:i)==','.or.str(i:i)==' '.or.iachar(str(i:i))==9.or.iachar(str(i:i))==11) then
            if(i-1>=m) then
               do j=1,kln
                  b(j:j)=' '
               enddo
               b(1:i-m+1)=str(m:i-1)
               read(b,*,iostat=err) rec(k)
               if(abs(err)>0)rec(k)=-1.d0
               k=k+1
            endif
            m=i+1
         endif
      enddo
      if(i-m<1)then
         k=k-1; goto 1101
      endif   
      b(1:i-m+1)=str(m:kln)
      read(b,*,iostat=err) rec(k)
      if(abs(err)>0)rec(k)=-1.d0
1101  nn=k
      end
!
!********************************************************
!
      subroutine PickRstrlg(line0,kln,str,nn)
      implicit none
	character*80000::line,line0,b
      character(len=25)::str(8000)
	integer::kln,nn,i,j,k,m,err,nval,sgn
!---------------------------------------------------------------------
      line=trim(AdjustL(line0))
      kln=len_trim(line)
      if(kln<2)then
         nn=0;return
      endif
      k=1; m=1
      do i=1,kln
         nval=iachar(line(i:i))
         if(nval==32.or.nval==44.or.nval==59.or.nval==9)line(i:i)=' '
         nval=iachar(line(i:i));sgn=1
         if(nval==32)sgn=0
         if(nval>42.and.nval<47)sgn=0
         if(nval>47.and.nval<58)sgn=0
         if(nval>64.and.nval<91)sgn=0
         if(nval>96.and.nval<123)sgn=0
         if(nval==95)sgn=0
         if(sgn>0)line(i:i)=''
         if(line(i:i)==' ') then
            if(i-1>=m) then
               do j=1,kln
                  b(j:j)=' '
               enddo
               b(1:i-m+1)=line(m:i-1)
               read(b,*,iostat=err)str(k)
               if(abs(err)>0)str(k)=''
               k=k+1
            endif
            m=i+1
         endif
      enddo
      if(i-m<1)then
         k=k-1; goto 1101
      endif   
      b(1:i-m+1)=line(m:kln)
      read(b,*,iostat=err)str(k)
      if(abs(err)>0)str(k)=''
1101  nn=k
      do i=1,nn
        str(i)=trim(str(i))
      enddo
      end
