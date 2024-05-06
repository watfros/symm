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
      subroutine inversion(natoms,nat,coord,delta,nc,ntrans,delta3)
************************************************************************
c
c Performs an inversion to the origin
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      parameter(mat=250)
      dimension nat(mat),coord(3,mat),cord(3,mat),ntrans(mat)
      do i=1,3
         do j=1,natoms
            cord(i,j)=-coord(i,j)
         end do
      end do
      call check(natoms,delta,nat,coord,cord,nc,ntrans,delta3)
      if(nc.ne.natoms) then
         do i=1,natoms
            ntrans(i)=i
         end do
      end if
      return
      end
************************************************************************
      subroutine rotate(natoms,nat,coord,v1,sina,cosa,delta,nc,ntrans,
     &   delta3)
************************************************************************
c
c Performs a rotation around axis v1
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      parameter(mat=250)
      dimension nat(mat),coord(3,mat),cord(3,mat),ntrans(mat)
      dimension v1(3),v2(3),v3(3),p0(3),p(3)
      do j=1,natoms
         p0(1)=coord(1,j)
         p0(2)=coord(2,j)
         p0(3)=coord(3,j)
         call crossp(v1,p0,v2)
         call crossp(v1,v2,v3)
         p=p0+sina*v2+(1.d0-cosa)*v3
         cord(1,j)=p(1)
         cord(2,j)=p(2)
         cord(3,j)=p(3)
      end do
      call check(natoms,delta,nat,coord,cord,nc,ntrans,delta3)
      return
      end
************************************************************************
      subroutine reflect(natoms,nat,coord,v,p0,delta,nc,ntrans,delta3)
************************************************************************
c
c Performs a reflection to the plane v
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      parameter(mat=250)
      dimension nat(mat),coord(3,mat),cord(3,mat),ntrans(mat)
      dimension v(3),p0(3),p(3)
      do i=1,natoms
         p(1)=coord(1,i)
         p(2)=coord(2,i)
         p(3)=coord(3,i)
         vk=-dot(v,p-p0,3)
         cord(1,i)=coord(1,i)+2.d0*vk*v(1)
         cord(2,i)=coord(2,i)+2.d0*vk*v(2)
         cord(3,i)=coord(3,i)+2.d0*vk*v(3)
      end do
      call check(natoms,delta,nat,coord,cord,nc,ntrans,delta3)
      return
      end
************************************************************************
      subroutine srotate(natoms,nat,coord,v1,sina,cosa,delta,nc,ntrans,
     &   delta3)
************************************************************************
c
c Improper rotation around axis v1
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      parameter(mat=250)
      dimension nat(mat),coord(3,mat),cord(3,mat),cc(3,mat),ntrans(mat)
      dimension v1(3),v2(3),v3(3),p0(3),p(3)
      do l=1,natoms
         p0(1)=coord(1,l)
         p0(2)=coord(2,l)
         p0(3)=coord(3,l)
         call crossp(v1,p0,v2)
         call crossp(v1,v2,v3)
         p=p0+sina*v2+(1.d0-cosa)*v3
         cord(1,l)=p(1)
         cord(2,l)=p(2)
         cord(3,l)=p(3)
      end do
      do j=1,natoms
         p0(1)=cord(1,j)
         p0(2)=cord(2,j)
         p0(3)=cord(3,j)
         vk=-dot(v1,p0,3)
         cc(1,j)=cord(1,j)+2.d0*vk*v1(1)
         cc(2,j)=cord(2,j)+2.d0*vk*v1(2)
         cc(3,j)=cord(3,j)+2.d0*vk*v1(3)
      end do
      call check(natoms,delta,nat,coord,cc,nc,ntrans,delta3)
      return
      end

