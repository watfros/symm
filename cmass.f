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
      subroutine cmass(natoms,nat,wt,coord,wmol,cmx,cmy,cmz)
************************************************************************
* Subroutine CMASS calculates the centre of mass of a molecule         *
************************************************************************
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      parameter(mat=250)
      dimension wt(90),coord(3,mat),nat(mat)
      sumwx=0.d0
      sumwy=0.d0
      sumwz=0.d0
      wmol=0.d0
      do i=1,natoms
         nati=nat(i)
         wmol=wmol+wt(nati)
         sumwx=sumwx+wt(nati)*coord(1,i)
         sumwy=sumwy+wt(nati)*coord(2,i)
         sumwz=sumwz+wt(nati)*coord(3,i)
      end do
      cmx=sumwx/wmol
      cmy=sumwy/wmol
      cmz=sumwz/wmol
      return
      end
************************************************************************
      subroutine cshift(natoms,coord,pc)
************************************************************************
* Subroutine cshift shifts the centre of mass of a molecule into the
* origin
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      parameter(mat=250)
      dimension coord(3,mat),pc(3)
      do i=1,natoms
         coord(1,i)=coord(1,i)-pc(1)
         coord(2,i)=coord(2,i)-pc(2)
         coord(3,i)=coord(3,i)-pc(3)
      end do
      return
      end
