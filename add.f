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
c   modified by Jinglin Mu(mjlink@126.com), to align with the VMD.
************************************************************************
      subroutine check(natoms,delta,nat,coord,c,nc,ntrans,delta3)
************************************************************************
c
c Subroutine check counts the number of atoms unchanged by an operation
c and generates the permutation describing it
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      parameter(mat=250)
      dimension coord(3,mat),c(3,mat),nat(mat),ntrans(mat),diff(3)
      nc=0
      delta3=0.d0
      outer: do i=1,natoms
         do j=1,natoms
c find closest atom
            if(nat(i).ne.nat(j)) cycle
            diff(1)=coord(1,i)-c(1,j)
            diff(2)=coord(2,i)-c(2,j)
            diff(3)=coord(3,i)-c(3,j)
            vn=dsqrt(dot(diff,diff,3))
            if(vn.le.delta) then
               nc=nc+1
c permutation
               ntrans(i)=j
               if(vn.gt.delta3) delta3=vn
               cycle outer
            end if
         end do
         return
      end do outer
      end
************************************************************************
      subroutine add_perm(natoms,ntrans,nprm,nper)
************************************************************************
c
c Subroutine add_perm stores new permutations
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      parameter(mat=250)
      dimension ntrans(mat),nper(mat,mat)
c check already stored permutations
      outer: do i=1,nprm
         do j=1,natoms
            if(ntrans(j).ne.nper(j,i)) cycle outer
         end do
         return
      end do outer
c add a new permutation
      nprm=nprm+1
      if(nprm.gt.mat) then
         write(*,'(a)')
     &    'ERROR: Too many symmetry operations. Try a lower tolerance.'
         stop
      end if
      do i=1,natoms
         nper(i,nprm)=ntrans(i)
      end do
      end
************************************************************************
      subroutine add_SG(nsg,sigman,v,p,delta)
************************************************************************
c
c Subroutine add_SG stores new reflection planes
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      parameter(nmax=150)
      dimension sigman(nmax,3)
      dimension v(3),p(3)
c check found planes
      do k=1,nsg
         p(1)=sigman(k,1)
         p(2)=sigman(k,2)
         p(3)=sigman(k,3)
         vk=dot(v,p,3)
         if(dabs(vk).gt.(1.d0-delta*delta).and.
     &      dabs(vk).lt.(1.d0+delta*delta)) return
      end do
c add a new plane
      nsg=nsg+1
      if(nsg.gt.nmax) then
         write(*,'(a)')
     &    'ERROR: Too many symmetry operations. Try a lower tolerance.'
         stop
      end if
      sigman(nsg,1)=v(1)
      sigman(nsg,2)=v(2)
      sigman(nsg,3)=v(3)
      end
************************************************************************
      subroutine add_Cn(nrot,rotn,rota,v,p,alpha,delta)
************************************************************************
c
c Subroutine add_Cn stores new rotation axes
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      parameter(nmax=150)
      dimension rotn(nmax,3),rota(nmax)
      dimension v(3),p(3)
c check found axes
      do k=1,nrot
         p(1)=rotn(k,1)
         p(2)=rotn(k,2)
         p(3)=rotn(k,3)
         vk=dot(v,p,3)
         if(dabs(vk).gt.(1.d0-delta*delta).and.
     &      dabs(vk).lt.(1.d0+delta*delta)) then 
            if(rota(k).gt.alpha) rota(k)=alpha
            return
         endif
      end do
c add a new axis
      nrot=nrot+1
      if(nrot.gt.nmax) then
         write(*,'(a)')
     &    'ERROR: Too many symmetry operations. Try a lower tolerance.'
         stop
      end if
      rotn(nrot,1)=v(1)
      rotn(nrot,2)=v(2)
      rotn(nrot,3)=v(3)
      rota(nrot)=alpha
      end
