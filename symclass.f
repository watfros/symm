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
      subroutine symclass(natoms,nprm,nper,nseq,nccl,nscl,nat,symb,nout)
************************************************************************
c
c Subroutine symclass detemines the equivalence classes defined by the
c symmetry operations
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      character symb*2
      parameter(mat=250)
      dimension nper(mat,mat),nscl(mat,mat),nccl(mat),nat(mat),symb(90)
      nseq=0
      outer: do i=1,natoms
         do j=1,nseq
            do k=1,nccl(j)
               if(nscl(k,j).eq.i) cycle outer
            end do
         end do
         nseq=nseq+1
         nccl(nseq)=0
         inner: do j=1,nprm
            do k=1,nccl(nseq)
               if(nper(i,j).eq.nscl(k,nseq)) cycle inner
            end do
            nccl(nseq)=nccl(nseq)+1
            nscl(nccl(nseq),nseq)=nper(i,j)
         end do inner
      end do outer
      do i=1,nseq
         do j=1,nccl(i)-1
            ii=j
            do k=j+1,nccl(i)
               if(nscl(k,i).lt.nscl(ii,i)) ii=k
            end do
            if(ii.ne.j) then
               itemp=nscl(j,i)
               nscl(j,i)=nscl(ii,i)
               nscl(ii,i)=itemp
            end if
         end do
      end do
      if(nout.eq.2) then
         write(*,'(/a,i3)')
     &      '-- Symmetry-equivalence classes of atoms: ',nseq
         do i=1,nseq
            write(*,'(/5x,a,i3,a7,a2,a1)') '#',i,
     &         ' (atom ',symb(nat(nscl(1,i))),')'
            write(*,'(5x,15i4)') (nscl(j,i),j=1,nccl(i))
         end do
      end if
      end
