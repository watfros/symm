c   Copyright (C) 2016 Laszlo Gyevi-Nagy, Gyula Tasi
c
c   This file is part of SYVA.
c
c   SYVA is free software; you can redistribute it and/or
c   modify it under the terms of the GNU Lesser General Public
c   License as published by the Free Software Foundation; either
c   version 2.1 of the License, or (at your option) any later version.
c
c   SYVA is distributed in the hope that it will be useful,
c   but WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c   Lesser General Public License for more details.
c
c   You should have received a copy of the GNU Lesser General Public
c   License along with SYVA; if not, write to the Free Software
c   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

************************************************************************
      integer function igcd(i1,i2)
************************************************************************
c
c Calculation of the greatest common divisor of two integers
c
      integer i1,i2,i,j,k
      i=i1
      j=i2
      do
         k=j
         j=mod(i,k)
         i=k
         if(j.eq.0) then
            igcd=i
            return
         end if
      end do
      end
************************************************************************
      double precision function dot(x,y,n)                         
************************************************************************
c
c Scalar (dot or inner) product of two vectors: <x|y>=x.y
c
      implicit double precision(a-h,o-z)                               
      implicit integer(i-n)
      dimension x(*),y(*)                                         
      dot=0.d0                                                          
      do i=1,n                                                          
         dot=dot+x(i)*y(i)                                              
      end do                                                            
      return                                                            
      end                                                               
************************************************************************
      subroutine crossp(x,y,z)
************************************************************************
c
c Vector (cross or outer) product of two vectors: z=[x,y] 
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      dimension x(3),y(3),z(3)
      z(1)=x(2)*y(3)-x(3)*y(2)
      z(2)=-x(1)*y(3)+x(3)*y(1)
      z(3)=x(1)*y(2)-x(2)*y(1)
      return
      end
************************************************************************
      subroutine timestamp(string)
************************************************************************
c
      implicit none
      character ampm*8,date*8
      character month(12)*9
      character time*10
      character string*14
      integer d,h,m,mm,n,s,y
      save month
      data month /
     &  'January  ', 'February ', 'March    ', 'April    ', 
     &  'May      ', 'June     ', 'July     ', 'August   ', 
     &  'September', 'October  ', 'November ', 'December ' /
      call date_and_time(date,time)
      read(date,'(i4,i2,i2)') y,m,d
      read(time,'(i2,i2,i2,1x,i3)') h,n,s,mm
      if(h.lt.12) then
        ampm='AM'
      elseif(h.eq.12) then
        if (n.eq.0.and.s.eq.0 ) then
          ampm='Noon'
        else
          ampm='PM'
        endif
      else
        h=h-12
        if(h.lt.12) then
          ampm='PM'
        elseif(h.eq.12) then
          if(n.eq.0.and.s.eq.0 ) then
            ampm='Midnight'
          else
            ampm='AM'
          endif
        endif
      endif
      write(*,'(/a,i3,1x,a,i5,i3,a1,i2.2,a1,i2.2,a1,i3.3,1x,a/)' )
     &  trim(string),d,trim(month(m)),y,h,':',n,':',s,'.',mm,ampm
      return
      end
************************************************************************
      subroutine swap(a,b)
************************************************************************
      double precision a,b,temp
      temp=a
      a=b
      b=temp
      end
************************************************************************
      subroutine swapi(a,b)
************************************************************************
      integer a,b,temp
      temp=a
      a=b
      b=temp
      end

************************************************************************
      logical function issubgroup(pg1, pg2)
************************************************************************
      implicit none
      character pg1*3, pg2*3

      integer nir, nsymop, nrotharm
      double precision chtab
      character pgsymb*3, irsymb*4
      common/chartab/ nir(2,55),chtab(14,322),
     &   nsymop(14,4,55),nrotharm(3,322),pgsymb(57),irsymb(322)
      integer nsgb, nsgr
      common/subgroups/ nsgb(2,57),nsgr(406)

      integer i, npg

      do npg=1,57
         if(adjustl(pg2).eq.adjustl(pgsymb(npg))) exit
      end do
      do i=nsgb(1,npg),nsgb(1,npg)+nsgb(2,npg)-1
        if(adjustl(pgsymb(nsgr(i))).eq.adjustl(pg1)) then
          issubgroup=.true.
          return
        end if
      end do
      issubgroup=.false.
      end function issubgroup
