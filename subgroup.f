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

c Selects symmetry operations of a chosen point group(subgroup). 
      subroutine select_subgroup(point_group,nsg,ncr,nsr,nsym,symn,
     &         nsg_ret,ncr_ret,nsr_ret,nsym_ret,symn_ret,ng_ret,np_ret,
     &         ni_ret,nmax)

      implicit none

      integer nmax   ! maximal dimension
      character point_group*3   ! selected point group
      integer nsg,ncr,nsr  ! number of different symmetry opertaions in the original point group
      integer nsg_ret,ncr_ret,nsr_ret  ! number of different symmetry opertaions in the selected point group
      integer ng_ret,np_ret,ni_ret     ! order of the group,the principal axis and inversion
      integer nsym   ! description of the symmetry operations of the original point group
      integer nsym_ret   ! description of the symmetry operations of the selected point group
      double precision symn
      double precision symn_ret

      character pgsymb*3,irsymb*3   ! point groups and irreducible representations
      integer nir,nsymop,nrotharm,nsgb,nsgroup
      double precision chtab

      double precision PI
      double precision TA
      integer nsg1,nsg2    ! first and last symmetry planes
      integer ncr1,ncr2    ! first and last proper axes
      integer nsr1,nsr2    ! first and last improper axes
      integer i,j,k   ! cycle variables
      integer c2,c2_atoms
      integer order,order2     ! order of the principal axis,two times the order
      integer pr_axis      ! which symmetry operation is the principal axis
      integer prev_axis    ! type of the previous axis
      double precision praxis_dv     ! direction vector of the principal axis
      double precision TOL     ! tolerance
      double precision theta     ! rotation angle
      double precision vec,vec_rot    ! vector to rotate,rotated vector
      double precision vec2,vec3

      double precision dot
      integer mm,mm2
      double precision temp1

      dimension nsym(nmax,5),nsym_ret(nmax,5),symn(3,nmax),
     &      symn_ret(3,nmax),praxis_dv(3),vec(3),vec_rot(3),vec2(3),
     &      vec3(3)
       
      common /chartab/ nir(2,55),chtab(14,322),
     &      nsymop(14,4,55),nrotharm(3,322),pgsymb(57),irsymb(322)
      common /subgroups/ nsgb(2,55),nsgroup(402)

c initialization
      TOL=2.d-2
      PI=4*datan(1.d0)
      TA=1./3.
      nsg1=2
      nsg2=nsg+1
      ncr1=nsg2+1
      ncr2=nsg2+ncr
      nsr1=ncr2+1
      nsr2=ncr2+nsr
      nsg_ret=0
      ncr_ret=0
      nsr_ret=0

      if(trim(point_group).eq.'Civ'.or.trim(point_group).eq.'Dih') then
c full point group copied
            ng_ret=-1
            np_ret=-1
            nsg_ret=nsg
            ncr_ret=ncr
            nsr_ret=nsr
            do mm2=1,5
               do mm=1,5
                  nsym_ret(mm2,mm)=nsym(mm2,mm)
               end do
               do mm=1,3
                  symn_ret(mm,mm2)=symn(mm,mm2)
               end do
            end do
            if(nsym_ret(1,2).eq.1) then
               ni_ret=1
            else
               ni_ret=0
            end if

      elseif(trim(point_group).eq.'C1') then
c only symmetry operation is E
            ng_ret=1
            np_ret=0
            ni_ret=0
            nsg_ret=0
            ncr_ret=0
            nsr_ret=0
            do mm=1,5
               nsym_ret(1,mm)=0
            end do
            do mm=1,3
               symn_ret(mm,1)=0.d0
            end do

      elseif(trim(point_group).eq.'Ci') then
c only symmetry operations are E and i
            ng_ret=2
            np_ret=0
            ni_ret=1
            nsg_ret=0
            ncr_ret=0
            nsr_ret=0
            do mm=1,5
               nsym_ret(1,mm)=nsym(1,mm)
            end do
            do mm=1,3
               symn_ret(mm,1)=0.d0
            end do

      elseif(trim(point_group).eq.'Cs') then
c only symmetry operations are E and sigma_h
            ng_ret=2
            np_ret=0
            ni_ret=0
            nsg_ret=1
            ncr_ret=0
            nsr_ret=0
            do mm=1,5
               nsym_ret(1,mm)=0
            end do
            do mm=1,3
               symn_ret(mm,1)=0.d0
            end do
            do mm=1,5
               nsym_ret(2,mm)=nsym(nsg1,mm)
            end do
            nsym_ret(2,4)=1
            do mm=1,3
               symn_ret(mm,2)=symn(mm,nsg1)
            end do

      elseif(trim(point_group).eq.'Ih') then
c all elements
            ng_ret=120
            np_ret=5
            ni_ret=1
            do mm2=1,120
               do mm=1,5
                  nsym_ret(mm2,mm)=nsym(mm2,mm)
               end do
               do mm=1,3
                  symn_ret(mm,mm2)=symn(mm,mm2)
               end do
            end do
            nsg_ret=nsg
            ncr_ret=ncr
            nsr_ret=nsr

      elseif(trim(point_group).eq.'I') then
c only the proper rotations
            ng_ret=60
            np_ret=5
            ni_ret=0
            do mm=1,5
               nsym_ret(1,mm)=0
            end do
            nsym_ret(1,3)=nsym(1,3)
            nsym_ret(1,5)=nsym(1,5)
            do mm=1,3
               symn_ret(mm,1)=0.d0
            end do
            do mm2=2,60
               do mm=1,5
                  nsym_ret(mm2,mm)=nsym(ncr1-2+mm2,mm)
               end do
               do mm=1,3
                  symn_ret(mm,mm2)=symn(mm,ncr1-2+mm2)
               end do
            end do
            nsg_ret=0
            ncr_ret=59
            nsr_ret=0

      elseif(trim(point_group).eq.'Oh') then
c all elements
            ng_ret=48
            np_ret=4
            ni_ret=1
            do mm2=1,48
               do mm=1,5
                  nsym_ret(mm2,mm)=nsym(mm2,mm)
               end do
               do mm=1,3
                  symn_ret(mm,mm2)=symn(mm,mm2)
               end do
            end do
            nsg_ret=nsg
            ncr_ret=ncr
            nsr_ret=nsr

      elseif(trim(point_group).eq.'O') then
c only the proper rotations
            ng_ret=24
            np_ret=4
            ni_ret=0
            do mm=1,5
               nsym_ret(1,mm)=0
            end do
            nsym_ret(1,3)=nsym(1,3)
            nsym_ret(1,5)=nsym(1,5)
            do mm=1,3
               symn_ret(mm,1)=0.d0
            end do
            do mm=1,5
               do mm2=2,24
                  nsym_ret(mm2,mm)=nsym(ncr1-2+mm2,mm)
               end do
            end do
            do mm=1,3
               do mm2=2,24
                  symn_ret(mm,mm2)=symn(mm,ncr1-2+mm2)
               end do
            end do
            nsg_ret=0
            ncr_ret=23
            nsr_ret=0

      elseif(trim(point_group).eq.'Td') then
c all C3s,S4s and sigma_ds
            np_ret=3
            ni_ret=0
            do mm=1,5
               nsym_ret(1,5)=0
            end do
            nsym_ret(1,3)=nsym(1,3)
            nsym_ret(1,5)=nsym(1,5)
            do mm=1,3
               symn_ret(mm,1)=0.d0
            end do
            j=1
            do i=nsg1,nsg2
               if(nsym(i,4).eq.3) then
                  j=j+1
                  do mm=1,5
                     nsym_ret(j,mm)=nsym(i,mm)
                  end do
                  do mm=1,3
                     symn_ret(mm,j)=symn(mm,i)
                  end do
               end if
            end do
            nsg_ret=j-1
            do i=ncr1,ncr2
               if(nsym(i,2).eq.3) then
                  j=j+1
                  do mm=1,5
                     nsym_ret(j,mm)=nsym(i,mm)
                  end do
                  do mm=1,3
                     symn_ret(mm,j)=symn(mm,i)
                  end do
               end if
            end do
            ncr_ret=j-1-nsg_ret
c find the C2 axes
            do i=nsr1,nsr2
               if(nsym(i,2).eq.4) then
                  j=j+1
                  do mm=1,5
                     nsym_ret(j,mm)=nsym(i,mm)
                  end do
                  do mm=1,3
                     symn_ret(mm,j)=symn(mm,i)
                  end do
               end if
            end do
            nsr_ret=j-1-nsg_ret-ncr_ret
            do i=ncr1,ncr2
               if(nsym(i,2).ne.2) cycle
               do k=2+nsg_ret+ncr_ret,1+nsg_ret+ncr_ret+nsr_ret
                  if(dabs(dabs(dot(symn_ret(1,k),symn(1,i),3))-1)
     &                  .le.TOL) then
                     do j=1+nsg_ret+ncr_ret+nsr_ret,2+nsg_ret+ncr_ret,-1
                        do mm=1,3
                           symn_ret(mm,j+1)=symn_ret(mm,j)
                        end do
                        do mm=1,5
                           nsym_ret(j+1,mm)=nsym_ret(j,mm)
                        end do
                     end do
                     ncr_ret=ncr_ret+1
                     do mm=1,5
                        nsym_ret(1+nsg_ret+ncr_ret,mm)=nsym(i,mm)
                     end do
                     do mm=1,3
                        symn_ret(mm,1+nsg_ret+ncr_ret)=symn(mm,i)
                     end do
                     exit
                  end if
               end do
            end do

            ng_ret=1+nsg_ret+ncr_ret+nsr_ret

      elseif(trim(point_group).eq.'T'.or.trim(point_group).eq.'Th') then
         if(point_group.eq.'T') then
c E
            np_ret=3
            ni_ret=0
            do mm=1,5
               nsym_ret(1,mm)=0
            end do
            do mm=1,3
               symn_ret(mm,1)=0.d0
            end do
         elseif(point_group.eq.'Th') then
c i
            np_ret=3
            ni_ret=1
            do mm=1,5
               nsym_ret(1,mm)=nsym(1,mm)
            end do
            do mm=1,3
               symn_ret(mm,1)=symn(mm,1)
            end do
         end if
         nsym_ret(1,3)=nsym(1,3)
         nsym_ret(1,5)=nsym(1,5)

c find a C3 axis
            do i=ncr1,ncr2
               if(nsym(i,2).eq.3) then
                  do mm=1,3
                     praxis_dv(mm)=symn(mm,i)
                  end do
                  exit
               end if
            end do
            do i=ncr1,ncr2
               if(nsym(i,2).ne.3) cycle
               if(dabs(dabs(dot(symn(1,i),praxis_dv,3))-1).le.TOL) then
                  ncr_ret=ncr_ret+1
                  do mm=1,5
                     nsym_ret(1+ncr_ret,mm)=nsym(i,mm)
                  end do
                  do mm=1,3
                     symn_ret(mm,1+ncr_ret)=symn(mm,i)
                  end do
               end if
            end do
c find the other C3 axes
            do i=ncr1,ncr2
               if(nsym(i,2).ne.3) cycle
               if(dabs(dabs(dot(symn(1,i),praxis_dv,3))-TA).le.TOL) then
                  do mm=1,3
                     vec(mm)=symn(mm,i)
                  end do
                  exit
               end if
            end do
            do i=1,3
               call subgroup_rotate(vec,praxis_dv,i*2*PI/3,vec_rot)
               do j=ncr1,ncr2
                  if(dabs(dabs(dot(symn(1,j),vec_rot,3))-1).le.TOL) then
                     ncr_ret=ncr_ret+1
                     do mm=1,5
                        nsym_ret(1+ncr_ret,mm)=nsym(j,mm)
                     end do
                     do mm=1,3
                        symn_ret(mm,1+ncr_ret)=symn(mm,j)
                     end do
                  end if
               end do
            end do

c find the C2 axes
            do i=3,1+ncr_ret
               if(dabs(dabs(dot(symn_ret(1,i),symn_ret(1,2),3))-1)
     &               .le.TOL) cycle
               k=floor(dot(symn_ret(1,i),symn_ret(1,2),3))
               if(k.lt.0) then
                  do mm=1,3
                     vec(mm)=symn_ret(mm,2)+symn_ret(mm,i)
                  end do
               else
                  do mm=1,3
                     vec(mm)=symn_ret(mm,2)-symn_ret(mm,i)
                  end do
               end if
               exit
            end do
            temp1=dsqrt(dot(vec,vec,3))
            do mm=1,3
               vec(mm)=vec(mm)/temp1
            end do
            do i=3,1+ncr_ret
               if(dabs(dabs(dot(symn_ret(1,i),symn_ret(1,2),3))-1)
     &               .le.TOL) cycle
               k=floor(dot(symn_ret(1,i),symn_ret(1,2),3))
               if(k.lt.0) then
                  do mm=1,3
                     vec2(mm)=symn_ret(mm,2)+symn_ret(mm,i)
                  end do
               else
                  do mm=1,3
                     vec2(mm)=symn_ret(mm,2)-symn_ret(mm,i)
                  end do
               end if
               temp1=dsqrt(dot(vec2,vec2,3))
               do mm=1,3
                  vec2(mm)=vec2(mm)/temp1
               end do
               if(dabs(dot(vec,vec2,3)).le.TOL) exit
            end do
            call crossp(vec,vec2,vec3)
            do i=ncr1,ncr2
               if(nsym(i,2).ne.2) cycle
               if(dabs(dabs(dot(vec,symn(1,i),3))-1).le.TOL.or.
     &               dabs(dabs(dot(vec2,symn(1,i),3))-1).le.TOL.or.
     &               dabs(dabs(dot(vec3,symn(1,i),3))-1).le.TOL) then
                  ncr_ret=ncr_ret+1
                  do mm=1,5
                     nsym_ret(1+ncr_ret,mm)=nsym(i,mm)
                  end do
                  do mm=1,3
                     symn_ret(mm,1+ncr_ret)=symn(mm,i)
                  end do
               end if
            end do

            if(point_group.eq.'Th') then
c find S6s
               do i=nsr1,nsr2
                  if(nsym(i,2).ne.6) cycle
                  do j=2,1+ncr_ret
                     if(nsym_ret(j,2).ne.3) cycle
                     if(dabs(dabs(dot(symn(1,i),symn_ret(1,j),3))-1)
     &                     .le.TOL) then
                        nsr_ret=nsr_ret+1
                        do mm=1,5
                           nsym_ret(1+ncr_ret+nsr_ret,mm)=nsym(i,mm)
                        end do
                        do mm=1,3
                           symn_ret(mm,1+ncr_ret+nsr_ret)=symn(mm,i)
                        end do
                        exit
                     end if
                  end do
               end do

c find sigma_hs
               do i=nsg1,nsg2
                  do j=2+nsg_ret,1+nsg_ret+ncr_ret
                     if(nsym_ret(j,2).ne.2) cycle
                     if(dabs(dabs(dot(symn(1,i),symn_ret(1,j),3))-1)
     &                     .le.TOL) then
                        do k=1+nsg_ret+ncr_ret+nsr_ret,2+nsg_ret,-1
                           do mm=1,5
                              nsym_ret(k+1,mm)=nsym_ret(k,mm)
                           end do
                           do mm=1,3
                              symn_ret(mm,k+1)=symn_ret(mm,k)
                           end do
                        end do
                        nsg_ret=nsg_ret+1
                        do mm=1,5
                           nsym_ret(1+nsg_ret,mm)=nsym(i,mm)
                        end do
                        nsym_ret(1+nsg_ret,4)=1
                        do mm=1,3
                           symn_ret(mm,1+nsg_ret)=symn(mm,i)
                        end do
                        exit
                     end if
                  end do
               end do

            end if

            ng_ret=1+ni_ret+nsg_ret+ncr_ret+nsr_ret

      else
c determine order of the principal axis
            read(point_group(2:2),'(i1)') order
            nsg_ret=0
            ncr_ret=0
            nsr_ret=0
c E
            do mm=1,5
               nsym_ret(1,mm)=0
            end do
            do mm=1,3
               symn_ret(mm,1)=0.d0
            end do
c O subspace
            if(point_group(1:1).eq.'D'.or.point_group(3:3).eq.'h') then
               nsym_ret(1,3)=nsym(1,3)
               nsym_ret(1,5)=nsym(1,5)
            end if
c find principal axis
            do i=ncr1,nsr2
               if(nsym(i,2).eq.order) then
                  pr_axis=i
                  do mm=1,3
                     praxis_dv(mm)=symn(mm,pr_axis)
                  end do
                  exit
               end if
            end do
            theta=PI/order

            if(point_group(1:1).eq.'S') then
c S2n point groups
c i
               if(mod(order,4).eq.2) then
                  do mm=1,2
                     nsym_ret(1,mm)=nsym(1,mm)
                  end do
               end if
               do mm=3,5
                  nsym_ret(1,mm)=nsym(1,mm)
               end do
c copy appropriate Cn axis
               do i=ncr1,ncr2
                  if(dabs(dabs(dot(praxis_dv,symn(1,i),3))-1).le.TOL
     &                  .and.mod(order/2,nsym(i,2)).eq.0) then
                     ncr_ret=ncr_ret+1
                     do mm=1,3
                        symn_ret(mm,1+ncr_ret)=symn(mm,i)
                     end do
                     do mm=1,5
                        nsym_ret(1+ncr_ret,mm)=nsym(i,mm)
                     end do
                     if(nsym_ret(1+ncr_ret,2).eq.order/2) then
                        nsym_ret(1+ncr_ret,4)=1
                     else
                        nsym_ret(1+ncr_ret,4)=0
                     end if
                  end if
               end do
c copy appropriate Sn axis
               do i=nsr1,nsr2
                  if(dabs(dabs(dot(praxis_dv,symn(1,i),3))-1).le.TOL
     &                  .and.nsym(i,2).eq.order) then
                     nsr_ret=nsr_ret+1
                     do mm=1,3
                        symn_ret(mm,1+ncr_ret+nsr_ret)=symn(mm,i)
                     end do
                     do mm=1,5
                        nsym_ret(1+ncr_ret+nsr_ret,mm)=nsym(i,mm)
                     end do
                  end if
               end do

            else
c copy appropriate Cn axis
               if(order.ne.2) then
                  do i=ncr1,ncr2
                     if(dabs(dabs(dot(praxis_dv,symn(1,i),3))-1).le.TOL
     &                     .and.mod(order,nsym(i,2)).eq.0) then
                        ncr_ret=ncr_ret+1
                        do mm=1,3
                           symn_ret(mm,1+ncr_ret)=symn(mm,i)
                        end do
                        do mm=1,5
                           nsym_ret(1+ncr_ret,mm)=nsym(i,mm)
                        end do
                        if(nsym_ret(1+ncr_ret,2).eq.order) then
                           nsym_ret(1+ncr_ret,4)=1
                        else
                           nsym_ret(1+ncr_ret,4)=0
                        end if
                     end if
                  end do

               else
                  if(point_group.eq.'D2d') then
c find the S4 axis
                     do i=nsr1,nsr2
                        if(nsym(i,2).eq.4) then
                           pr_axis=i
                           do mm=1,3
                              praxis_dv(mm)=symn(mm,i)
                           end do
                           exit
                        end if
                     end do
                     do i=ncr1,ncr2
                        if(nsym(i,2).eq.2.and.dabs(dabs(dot(praxis_dv,
     &                        symn(1,i),3))-1).le.TOL) then
                           pr_axis=i
                           do mm=1,3
                              praxis_dv(mm)=symn(mm,i)
                           end do
                           ncr_ret=ncr_ret+1
                           do mm=1,3
                              symn_ret(mm,1+ncr_ret)=symn(mm,i)
                           end do
                           do mm=1,5
                              nsym_ret(1+ncr_ret,mm)=nsym(i,mm)
                           end do
                           nsym_ret(1+ncr_ret,4)=1
                           exit
                        end if
                     end do
c select the C2 axis with the most atoms on it
                  else
                     c2=0
                     c2_atoms=-1
                     do i=ncr1,ncr2
                        if(nsym(i,2).eq.2) then
                           if(nsym(i,5).gt.c2_atoms) then
                              c2=i
                              c2_atoms=nsym(i,5)
                           end if
                        end if
                     end do
                     pr_axis=c2
                     do mm=1,3
                        praxis_dv(mm)=symn(mm,c2)
                     end do
                     ncr_ret=ncr_ret+1
                     do mm=1,3
                        symn_ret(mm,1+ncr_ret)=symn(mm,c2)
                     end do
                     do mm=1,5
                        nsym_ret(1+ncr_ret,mm)=nsym(c2,mm)
                     end do
                     nsym_ret(1+ncr_ret,4)=1
                  end if
               end if

               if(point_group(1:1).eq.'D') then
c select perpendicular C2 axes with the most atoms
                  c2=0
                  c2_atoms=-1
                  do i=ncr1,ncr2
                     if(nsym(i,2).eq.2.and.
     &                     dabs(dot(praxis_dv,symn(1,i),3)).le.TOL) then
                        if(nsym(i,5).gt.c2_atoms) then
                           c2=i
                           c2_atoms=nsym(i,5)
                        end if
                     end if
                  end do
                  ncr_ret=ncr_ret+1
                  do mm=1,3
                     symn_ret(mm,1+ncr_ret)=symn(mm,c2)
                     vec(mm)=symn(mm,c2)
                  end do
                  do mm=1,5
                     nsym_ret(1+ncr_ret,mm)=nsym(c2,mm)
                  end do
                  nsym_ret(1+ncr_ret,4)=2
                  prev_axis=2
                  do i=1,order-1
                     call subgroup_rotate(vec,praxis_dv,i*theta,vec_rot)
                     do j=ncr1,ncr2
                        if(nsym(j,2).ne.2.or.dabs(dabs(dot(vec_rot,
     &                        symn(1,j),3))-1).gt.TOL) cycle
                        ncr_ret=ncr_ret+1
                        do mm=1,3
                           symn_ret(mm,1+ncr_ret)=symn(mm,j)
                        end do
                        do mm=1,5
                           nsym_ret(1+ncr_ret,mm)=nsym(j,mm)
                        end do
                        if(point_group(3:3).eq.'d'.or.mod(order,2).eq.1
     &                        .or.prev_axis.eq.3) then
                           nsym_ret(1+ncr_ret,4)=2
                           prev_axis=2
                        else
                           nsym_ret(1+ncr_ret,4)=3
                           prev_axis=3
                        end if
                     end do
                  end do
               end if

c copy Sn axes
               if((point_group(3:3).eq.'d').or.(point_group(3:3).eq.'h'
     &               .and.order.ne.2)) then
                  if(point_group(3:3).eq.'d') then
                     order2=2*order
                  else
                     order2=order
                  end if
                  do i=nsr1,nsr2
                     if(dabs(dabs(dot(praxis_dv,symn(1,i),3))-1).le.TOL)
     &                     then
                        if(point_group(3:3).eq.'h') then
                           if(order.eq.6.or.order.eq.8) then
                              if(mod(order2,nsym(i,2)).eq.0) then
                                 nsr_ret=nsr_ret+1
                                 do mm=1,3
                                    symn_ret(mm,1+ncr_ret+nsr_ret)=
     &                                    symn(mm,i)
                                 end do
                                 do mm=1,5
                                    nsym_ret(1+ncr_ret+nsr_ret,mm)=
     &                                    nsym(i,mm)
                                 end do
                              end if
                           else
                              if(nsym(i,2).eq.order2) then
                                 nsr_ret=nsr_ret+1
                                 do mm=1,3
                                    symn_ret(mm,1+ncr_ret+nsr_ret)=
     &                                    symn(mm,i)
                                 end do
                                 do mm=1,5
                                    nsym_ret(1+ncr_ret+nsr_ret,mm)=
     &                                    nsym(i,mm)
                                 end do
                              end if
                           end if
                        else
                           if(order.eq.6) then
                              if(nsym(i,2).eq.12.or.nsym(i,2).eq.4) then
                                 nsr_ret=nsr_ret+1
                                 do mm=1,3
                                    symn_ret(mm,1+ncr_ret+nsr_ret)=
     &                                    symn(mm,i)
                                 end do
                                 do mm=1,5
                                    nsym_ret(1+ncr_ret+nsr_ret,mm)=
     &                                    nsym(i,mm)
                                 end do
                              end if
                           else
                              if(nsym(i,2).eq.order2) then
                                 nsr_ret=nsr_ret+1
                                 do mm=1,3
                                    symn_ret(mm,1+ncr_ret+nsr_ret)=
     &                                    symn(mm,i)
                                 end do
                                 do mm=1,5
                                    nsym_ret(1+ncr_ret+nsr_ret,mm)=
     &                                    nsym(i,mm)
                                 end do
                              end if
                           end if
                        end if
                     end if
                  end do
               end if

c copy inversion
               if((point_group(3:3).eq.'h'.and.mod(order,2).eq.0).or.
     &             (point_group(3:3).eq.'d'.and.mod(order,2).eq.1)) then
                  do mm=1,5
                     nsym_ret(1,mm)=nsym(1,mm)
                  end do
                  do mm=1,3
                     symn_ret(mm,1)=0.d0
                  end do
               end if

c copy sigma_h
               if(point_group(3:3).eq.'h') then
                  do i=nsg1,nsg2
                     if(dabs(dabs(dot(praxis_dv,symn(1,i),3))-1).le.TOL)
     &                     then
                        do j=ncr_ret+nsr_ret+1,2,-1
                           do mm=1,3
                              symn_ret(mm,j+1)=symn_ret(mm,j)
                           end do
                           do mm=1,5
                              nsym_ret(j+1,mm)=nsym_ret(j,mm)
                           end do
                        end do
                        do mm=1,3
                           symn_ret(mm,2)=symn(mm,i)
                        end do
                        do mm=1,5
                           nsym_ret(2,mm)=nsym(i,mm)
                        end do
                        nsym_ret(2,4)=1
                        nsg_ret=1
                        exit
                     end if
                  end do
               end if

c copy sigma_d,Dnd point groups
               if(point_group(3:3).eq.'d') then
                  if(mod(order,2).eq.0) then
                     do mm=1,3
                        vec(mm)=symn_ret(mm,1+ncr_ret)+
     &                        symn_ret(mm,ncr_ret)
                     end do
                     temp1=dsqrt(dot(vec,vec,3))
                     do mm=1,3
                        vec(mm)=vec(mm)/temp1
                     end do
                  else
                     do mm=1,3
                        vec(mm)=symn_ret(mm,ncr_ret)
                     end do
                  end if

                  do i=1,order
                     call subgroup_rotate(vec,praxis_dv,i*theta,vec_rot)
                     do j=nsg1,nsg2
                        if(dabs(dabs(dot(vec_rot,symn(1,j),3))-1).gt.
     &                        TOL) cycle
                        do k=nsg_ret+ncr_ret+nsr_ret+1,nsg_ret+2,-1
                           do mm=1,3
                              symn_ret(mm,k+1)=symn_ret(mm,k)
                           end do
                           do mm=1,5
                              nsym_ret(k+1,mm)=nsym_ret(k,mm)
                           end do
                        end do
                        nsg_ret=nsg_ret+1
                        do mm=1,3
                           symn_ret(mm,1+nsg_ret)=symn(mm,j)
                        end do
                        do mm=1,5
                           nsym_ret(1+nsg_ret,mm)=nsym(j,mm)
                        end do
                        nsym_ret(1+nsg_ret,4)=3
                     end do
                  end do
               end if

c copy remaining sigmas
               if((point_group(3:3).eq.'v').or.
     &                (point_group(1:1).eq.'D'.and.
     &                  point_group(3:3).eq.'h')) then
                  do i=nsg1,nsg2
c  find a plane with orthogonal normal vector
                     if(dabs(dot(praxis_dv,symn(1,i),3))
     &                     .le.TOL) then
                        do j=nsg_ret+ncr_ret+nsr_ret+1,nsg_ret+2,-1
                           do mm=1,3
                              symn_ret(mm,j+1)=symn_ret(mm,j)
                           end do
                           do mm=1,5
                              nsym_ret(j+1,mm)=nsym_ret(j,mm)
                           end do
                        end do
                        nsg_ret=nsg_ret+1
                        do mm=1,3
                           symn_ret(mm,1+nsg_ret)=symn(mm,i)
                        end do
                        do mm=1,3
                           vec(mm)=symn(mm,i)
                        end do
                        do mm=1,5
                           nsym_ret(1+nsg_ret,mm)=nsym(i,mm)
                        end do
                        if(mod(order,2).eq.1) then
                           nsym_ret(1+nsg_ret,4)=2
                        else
                           if(nsym(i,4).ne.1) then
                              prev_axis=nsym(i,4)
                           else
                              prev_axis=2
                           end if
                           nsym_ret(1+nsg_ret,4)=prev_axis
                        end if
                        exit
                     end if
                  end do

                  do i=1,order-1
                     call subgroup_rotate(vec,praxis_dv,i*theta,vec_rot)
                     do j=nsg1,nsg2
                        if(dabs(dabs(dot(vec_rot,symn(1,j),3))-1).gt.
     &                        TOL) cycle
                        do k=nsg_ret+ncr_ret+nsr_ret+1,nsg_ret+2,-1
                           do mm=1,3
                              symn_ret(mm,k+1)=symn_ret(mm,k)
                           end do
                           do mm=1,5
                              nsym_ret(k+1,mm)=nsym_ret(k,mm)
                           end do
                        end do
                        nsg_ret=nsg_ret+1
                        do mm=1,3
                           symn_ret(mm,1+nsg_ret)=symn(mm,j)
                        end do
                        do mm=1,5
                           nsym_ret(1+nsg_ret,mm)=nsym(j,mm)
                        end do
                        if(mod(order,2).eq.1) then
                           nsym_ret(1+nsg_ret,4)=2
                        else
                           if(prev_axis.eq.2) then
                              nsym_ret(1+nsg_ret,4)=3
                              prev_axis=3
                           else
                              nsym_ret(1+nsg_ret,4)=2
                              prev_axis=2
                           end if
                        end if
                     end do
                  end do
               end if

c C2v subgroup of Dnd point groups when there's no atom on the principal axis
               if(point_group.eq.'C2v'.and.nsg_ret.ne.2) then
                  do i=ncr1,ncr2
                     if(nsym(i,4).eq.1) then
                        nsym_ret(4,1)=2
                        nsym_ret(4,2)=2
                        nsym_ret(4,3)=1
                        nsym_ret(4,4)=1
                        nsym_ret(4,5)=nsym(i,5)
                        do mm=1,3
                           symn_ret(mm,4)=symn(mm,i)
                        end do
                        exit
                     end if
                  end do

                  nsym_ret(3,1)=1
                  nsym_ret(3,2)=0
                  nsym_ret(3,3)=0
                  nsym_ret(3,4)=2
                  nsym_ret(3,5)=nsym(nsg1,5)
                  do mm=1,3
                     symn_ret(mm,3)=symn(mm,nsg1)
                  end do

                  do i=nsg1+1,nsg2
                     if(dabs(dot(symn(1,nsg1),symn(1,i),3)).le.TOL) then
                        nsym_ret(2,1)=1
                        nsym_ret(2,2)=0
                        nsym_ret(2,3)=0
                        nsym_ret(2,4)=3
                        nsym_ret(2,5)=nsym(i,5)
                        do mm=1,3
                           symn_ret(mm,2)=symn(mm,i)
                        end do
                        exit
                     end if
                  end do

                  nsg_ret=2
                  ncr_ret=1
                  nsr_ret=0
               end if

            end if

            if(nsym_ret(1,2).eq.1) then
               ni_ret=1
            else
               ni_ret=0
            end if

            ng_ret=1+ni_ret+nsg_ret+ncr_ret+nsr_ret
            np_ret=order
            if(point_group(1:1).eq.'S') np_ret=np_ret/2

      end if

      contains
      subroutine subgroup_rotate(vector,direction,angle,rotated)
      implicit none
      double precision vector,rotated,direction
      double precision angle
      double precision vec1,vec2
      integer i
      dimension vector(3),rotated(3),direction(3),vec1(3),vec2(3)
      call crossp(direction,vector,vec1)
      call crossp(direction,vec1,vec2)
      do i=1,3
         rotated(i)=vector(i)+dsin(angle)*vec1(i)
     &      +(1-dcos(angle))*vec2(i)
      end do
      end subroutine subgroup_rotate

      end subroutine select_subgroup
