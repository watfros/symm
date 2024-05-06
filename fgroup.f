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
      subroutine symrank(i,k,natsym,nsym)
************************************************************************
c
c Subroutine symrank whether atom k lies on symmetry element i
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      parameter(mat=250,nmax=150)
      dimension natsym(mat),nsym(nmax,5)
      j=natsym(k)
      if(j.eq.0) then
         natsym(k)=i
      elseif(j.eq.1.or.nsym(j,2).eq.0) then
         return
      else
         i1=nsym(i,1)
         i2=nsym(i,2)
         if(i1.eq.1) i2=nsym(i,4)
         j1=nsym(j,1)
         j2=nsym(j,2)
         if(j1.eq.1) j2=nsym(j,4)
         if(i1.lt.j1) return 
         if((i1.eq.1).and.(i2.gt.j2)) return
         if((i1.eq.2).and.(i2.lt.j2)) return
         natsym(k)=i
      endif
      return
      end
************************************************************************
      subroutine fw_group(natoms,nat,coord,symb,nsg,ncr,nsr,symn,nsym,
     &                    natsym,delta,pgrp,nout)
************************************************************************
c
c Subroutine fw_group determines the type of symmetry operations based on
c the point group and the framework group
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      character symb*2,symel*3,pgrp*3,fwgrp*90,ch*3
      logical se1,se2,center
      parameter(mat=250,nmax=150)
      dimension coord(3,mat),nat(mat),symb(90)
      dimension symn(3,nmax),nsym(nmax,5),natsym(mat)
      dimension v0(3),v1(3),p1(3),p2(3),p3(3),p(3)
      dimension msym(nmax,2),matsym(nmax,10,2)
      dimension a1(3),a2(3),d(10)
      msym = 0
      matsym = 0
      pi=4.d0*datan(1.d0)
c
c Set the exact type of symmetry elements
c
c Planes of symmetry
c
      if(pgrp(len_trim(pgrp):len_trim(pgrp)).eq.'d') then
         do i=2,nsg+1
            nsym(i,4)=3
         end do
      elseif(trim(pgrp).eq.'C2v'.or.trim(pgrp).eq.'C4v'.or.
     &       trim(pgrp).eq.'C6v'.or.trim(pgrp).eq.'C8v') then
         se1=.false.
         do i=2,nsg+1
            if(nsym(i,4).eq.3) se1=.true.
         end do
         if(.not.se1) then
            do i1=nsg+2,nsg+ncr+1
               if(nsym(i1,4).eq.1) then
                  v1=symn(:,i1)
                  alpha=2.d0*pi/dble(nsym(i1,2))
                  exit
               end if
            end do
            p1=symn(:,2)
            call crossp(v1,p1,p2)
            call crossp(v1,p2,v0)
            do i=0,nsym(i1,2)-1
               sina=dsin(dble(i)*alpha)
               cosa=dcos(dble(i)*alpha)
               p=p1+sina*p2+(1.d0-cosa)*v0
               do k=2,nsg+1
                  p3=symn(:,k)
                  if(dabs(dabs(dot(p3,p,3))-1.d0).lt.delta) nsym(k,4)=3
               end do
            end do
         end if
      elseif(trim(pgrp).eq.'Oh') then
         do i=2,nsg+1
            if(nsym(i,4).eq.2) nsym(i,4)=3
         end do
      elseif(trim(pgrp).eq.'Ih'.or.trim(pgrp).eq.'Th') then
         do i=2,nsg+1
            nsym(i,4)=1
         end do
      elseif(trim(pgrp).eq.'Cs') then
         do i=2,nsg+1
            nsym(i,4)=1
         end do
      end if
c
c Symmetry axes
c
      if(trim(pgrp).eq.'Ih'.or.trim(pgrp).eq.'I') then
         do i=nsg+2,nsg+ncr+1
            nsym(i,4)=0
         end do
      elseif(trim(pgrp).eq.'Oh'.or.trim(pgrp).eq.'O') then
         do i=nsg+2,nsg+ncr+1
            if(nsym(i,2).eq.2) then
               p1=symn(:,i)
               do j=nsg+2,nsg+ncr+1
                  if(nsym(j,2).eq.4) then
                     p2=symn(:,j)
                     if(dabs(dabs(dot(p1,p2,3))-1.d0).lt.delta) then
                        nsym(i,4)=0
                        goto 5
                     end if
                  end if
               end do
               nsym(i,4)=2
            end if
 5          continue
         end do
      elseif(trim(pgrp).eq.'D2') then
         if(nsym(2,5).ge.nsym(3,5).and.nsym(2,5).ge.nsym(4,5)) then
            if(nsym(3,5).ge.nsym(4,5)) then
               nsym(2,4)=1
               nsym(3,4)=2
               nsym(4,4)=3
            else
               nsym(2,4)=1
               nsym(3,4)=3
               nsym(4,4)=2
            end if
         elseif(nsym(3,5).ge.nsym(2,5).and.nsym(3,5).ge.nsym(4,5)) then
            if(nsym(2,5).ge.nsym(4,5)) then
               nsym(2,4)=2
               nsym(3,4)=1
               nsym(4,4)=3
            else
               nsym(2,4)=3
               nsym(3,4)=1
               nsym(4,4)=2
            end if
         elseif(nsym(4,5).ge.nsym(2,5).and.nsym(4,5).ge.nsym(3,5)) then
            if(nsym(2,5).ge.nsym(3,5)) then
               nsym(2,4)=2
               nsym(3,4)=3
               nsym(4,4)=1
            else
               nsym(2,4)=3
               nsym(3,4)=2
               nsym(4,4)=1
            end if
         end if
      elseif(trim(pgrp).eq.'D4'.or.trim(pgrp).eq.'D6'.or.
     &   trim(pgrp).eq.'D8') then
         se1=.false.
         do i=nsg+2,nsg+ncr+1
            if(nsym(i,4).eq.3) se1=.true.
         end do
         if(.not.se1) then
            nc2=0
            nc2_at=-1
            do i1=nsg+2,nsg+ncr+1
               if(nsym(i1,2).eq.2.and.nsym(i1,4).eq.2) then
                  if(nsym(i1,5).gt.nc2_at) then
                     nc2=i1
                     nc2_at=nsym(i1,5)
                  end if
                  nsym(i1,4)=3
               end if
               if(nsym(i1,4).eq.1) then
                  v1=symn(:,i1)
                  alpha=2.d0*pi/dble(nsym(i1,2))
                  ii=i1
               end if
            end do
            p1=symn(:,nc2)
            call crossp(v1,p1,p2)
            call crossp(v1,p2,v0)
            do i=0,nsym(ii,2)-1
               sina=dsin(dble(i)*alpha)
               cosa=dcos(dble(i)*alpha)
               p=p1+sina*p2+(1.d0-cosa)*v0
               do k=nsg+2,nsg+ncr+1
                  p3=symn(:,k)
                  if(dabs(dabs(dot(p3,p,3))-1.d0).lt.delta) nsym(k,4)=2
               end do
            end do
         end if
      end if
c
c Dnh point groups
c
      if(trim(pgrp).eq.'D2h') then
         if(nsym(5,4).eq.1) then
            i1=6
            i2=7
         elseif(nsym(6,4).eq.1) then
            i1=5
            i2=7
         elseif(nsym(7,4).eq.1) then
            i1=5
            i2=6
         end if
         if(nsym(i1,5).ge.nsym(i2,5)) then
            nsym(i1,4)=2
            nsym(i2,4)=3
            p1=symn(:,i1)
         elseif(nsym(i2,5).ge.nsym(i1,5)) then
            nsym(i2,4)=2
            nsym(i1,4)=3
            p1=symn(:,i2)
         end if
         do i=2,nsg+1
            if(nsym(i,4).eq.1) cycle
            p2=symn(:,i)
            if(dabs(dabs(dot(p1,p2,3))-1.d0).lt.delta) then
               nsym(i,4)=2
            else
               nsym(i,4)=3
            end if
         end do
      elseif(trim(pgrp).eq.'D4h'.or.
     &   trim(pgrp).eq.'D6h'.or.trim(pgrp).eq.'D8h') then
         se1=.false.
         se2=.false.
         do i=2,nsg+1
            if(nsym(i,4).eq.3) se1=.true.
         end do
         do i=nsg+2,nsg+ncr+1
            if(nsym(i,2).eq.2.and.nsym(i,4).eq.3) se2=.true.
         end do
         if(.not.se1.and..not.se2) then
            do i1=nsg+2,nsg+ncr+1
               if(nsym(i1,2).eq.2.and.nsym(i1,4).eq.2) p1=symn(:,i1)
               if(nsym(i1,4).eq.1) then
                  v1=symn(:,i1)
                  alpha=2.d0*pi/dble(nsym(i1,2))
                  ii=i1
               end if
            end do
            call crossp(v1,p1,p2)
            call crossp(v1,p2,v0)
            do i=0,nsym(ii,2)-1
               sina=dsin(dble(i)*alpha)
               cosa=dcos(dble(i)*alpha)
               p=p1+sina*p2+(1.d0-cosa)*v0
               do k=nsg+2,nsg+ncr+1
                  p3=symn(:,k)
                  if(dabs(dabs(dot(p3,p,3))-1.d0).lt.delta) nsym(k,4)=3
               end do
            end do
            do i=nsg+2,nsg+nsr+1
               if(nsym(i,2).ne.2) cycle
               if(nsym(i,4).eq.3) then
                  p1=symn(:,i)
                  do j=2,nsg+1
                     if(nsym(j,4).ne.2) cycle
                     p2=symn(:,j)
                     if(dabs(dabs(dot(p1,p2,3))-1.d0).lt.delta)
     &                  nsym(j,4)=3
                  end do
               end if
            end do
            se1=.true.
            se2=.true.
         elseif(.not.se1.and.se2) then
            do i=nsg+2,nsg+nsr+1
               if(nsym(i,2).ne.2) cycle
               if(nsym(i,4).eq.3) then
                  p1=symn(:,i)
                  do j=2,nsg+1
                     if(nsym(j,4).ne.2) cycle
                     p2=symn(:,j)
                     if(dabs(dabs(dot(p1,p2,3))-1.d0).lt.delta)
     &                  nsym(j,4)=3
                  end do
               end if
            end do
            se1=.true.
         elseif(se1.and..not.se2) then
            do i=2,nsg+1
               if(nsym(i,4).eq.3) then
                  p1=symn(:,i)
                  do j=nsg+2,nsg+ncr+1
                     if(nsym(j,2).ne.2.or.nsym(j,4).ne.2) cycle
                     p2=symn(:,j)
                     if(dabs(dabs(dot(p1,p2,3))-1.d0).lt.delta)
     &                  nsym(j,4)=3
                  end do
               end if
            end do
            se2=.true.
         end if
      end if
c Determine framework group
c
c X subspace
      do i=1,natoms
         natsym(i)=0
      end do
c O subspace
      if(nsym(1,3).gt.0) natsym(nsym(1,3))=1
c Cn subspace
      do i=nsg+2,nsg+ncr+1
         if(nsym(i,3).ne.1) cycle
         v1=symn(:,i)
         do j=1,natoms
            p1(1)=coord(1,j)
            p1(2)=coord(2,j)
            p1(3)=coord(3,j)                                  
            call crossp(v1,p1,v0)
            vn=dsqrt(dot(v0,v0,3))
            if(dabs(vn).lt.delta) then
               call symrank(i,j,natsym,nsym)
            endif
         end do
      end do
c SGH, SGV, SGD subspaces
      do i=2,nsg+1
         v1=symn(:,i)
         do j=1,natoms
            p1(1)=coord(1,j)
            p1(2)=coord(2,j)
            p1(3)=coord(3,j)
            sp=dot(v1,p1,3)
            if(dabs(sp).lt.delta) then
               call symrank(i,j,natsym,nsym)
            endif         
         end do
      end do
c
c Types of atoms and their number on each symmetry element
c
      do i=1,nsg+ncr+1
         a1=0.d0
         nm=0
         center=.false.
         do j=1,natoms
            if(natsym(j).ne.i) cycle
            if(i.gt.nsg+1) then
               if(dsqrt(dabs(dot(a1,a1,3))).lt.delta) then
                  a1=coord(:,j)
               end if
               a2=coord(:,j)
               dist=dot(a1,a2,3)
            end if
            natj=nat(j)
            do k=1,10
               if(i.le.nsg+1) then
                  if(natj.eq.matsym(i,k,1)) then
                     matsym(i,k,2)=matsym(i,k,2)+1
                     exit
                  end if
               end if
               if(matsym(i,k,1).eq.0) then
                  nm=k
                  matsym(i,k,1)=natj
                  matsym(i,k,2)=1
                  if(i.gt.nsg+1) d(k)=dist
                  if(dsqrt(dabs(dot(a2,a2,3))).lt.delta) then
                     d(k)=0.d0
                     center=.true.
                  end if
                  exit
               end if
            end do
         end do
         if(i.gt.nsg+1) then
            nm=nm+1
            matsym(i,nm,1)=-1
            matsym(i,nm,2)=1
            d(nm)=0.d0
            if(center) then
               d(nm)=-delta
               nm=nm+1
               matsym(i,nm,1)=-1
               matsym(i,nm,2)=1
               d(nm)=delta
            end if
         end if
c Ordering in decreasing atomic number
         do j=1,nm-1
            maxnat=j
            do k=j+1,nm
               if(i.le.nsg+1) then
                  if(matsym(i,k,1).gt.matsym(i,maxnat,1)) maxnat=k
               else
                  if(d(k).gt.d(maxnat)) maxnat=k
               end if
            end do
            call swapi(matsym(i,j,1),matsym(i,maxnat,1))
            call swapi(matsym(i,j,2),matsym(i,maxnat,2))
            if(i.gt.nsg+1) call swap(d(j),d(maxnat))
         end do
         if(i.gt.nsg+1) then
            do k=1,nm
               if(matsym(i,k,1).eq.-1) ncm=k
            end do
            if(ncm.eq.nm.and..not.center) then
               matsym(i,nm,1)=0
               matsym(i,nm,2)=0
            else
               if(nm-ncm+1.lt.ncm) then
                  do k=1,nm/2
                     call swapi(matsym(i,k,1),matsym(i,nm-k+1,1))
                     call swapi(matsym(i,k,2),matsym(i,nm-k+1,2))
                  end do
               end if
            end if
         end if
         msym(i,1)=0
         msym(i,2)=i
      end do
c X subspace
      ii=nsg+ncr+2
      do j=1,natoms
         if(natsym(j).eq.0) then
            do k=1,10
               if(nat(j).eq.matsym(ii,k,1)) then
                  matsym(ii,k,2)=matsym(ii,k,2)+1
                  exit
               end if
               if(matsym(ii,k,1).eq.0) then
                  matsym(ii,k,1)=nat(j)
                  matsym(ii,k,2)=1
                  exit
               end if
            end do
         end if
      end do
      do j=1,9
         maxnat=j
         do k=j+1,10
            if(matsym(ii,k,1).gt.matsym(ii,maxnat,1)) maxnat=k
         end do
         call swapi(matsym(ii,j,1),matsym(ii,maxnat,1))
         call swapi(matsym(ii,j,2),matsym(ii,maxnat,2))
      end do
c
c Counting equivalent symmetry elements
c
      msym(1,1)=1
      msym(ii,1)=1
      do i=2,nsg+ncr+1
         if(matsym(i,1,1).le.0.and.matsym(i,2,1).le.0) cycle
         do j=1,i-1
            if(nsym(i,1).eq.nsym(j,1).and.nsym(i,2).eq.nsym(j,2).and.
     &         nsym(i,4).eq.nsym(j,4).and.nsym(i,5).eq.nsym(j,5)) then
               do k=1,10
                  if(matsym(i,k,1).ne.matsym(j,k,1).or.
     &               matsym(i,k,2).ne.matsym(j,k,2)) then
                     goto 10
                  end if
               end do
               msym(j,1)=msym(j,1)+1
               goto 15
            end if
 10         continue
         end do
         msym(i,1)=msym(i,1)+1
 15      continue
      end do
c
c Ordering: 
c     1. O  subspace
c     2. Cn subspaces, increasing order
c     3. SG subspaces (SGH, SHV, SHD)
c     4. X  subspace (ii=nsg+ncr+2)
c
      do i=2,nsg+ncr+1
         next=i
         do j=i+1,nsg+ncr+1
            if(matsym(j,1,1).eq.0) cycle
            if(nsym(msym(j,2),1).gt.nsym(msym(next,2),1)) then
               next=j
            elseif(nsym(msym(j,2),1).eq.nsym(msym(next,2),1)) then
               if(nsym(msym(j,2),1).eq.1) then
                  if(nsym(msym(j,2),4).lt.nsym(msym(next,2),4)) then
                     next=j
                  elseif(nsym(msym(j,2),4).eq.
     &               nsym(msym(next,2),4)) then
                     if(nsym(msym(j,2),5).gt.nsym(msym(next,2),5))
     &                  next=j
                  end if
               else
                  if(nsym(msym(j,2),2).gt.nsym(msym(next,2),2)) then
                     next=j
                  elseif(nsym(msym(j,2),2).eq.nsym(msym(next,2),2)) then
                     if(nsym(msym(j,2),4).lt.nsym(msym(next,2),4)) then
                        next=j
                     elseif(nsym(msym(j,2),4).eq.
     &                  nsym(msym(next,2),4)) then
                        if(nsym(msym(j,2),5).gt.
     &                     nsym(msym(next,2),5)) next=j
                     end if
                  end if
               end if
            end if
         end do
         call swapi(msym(i,1),msym(next,1))
         call swapi(msym(i,2),msym(next,2))
         do j=1,10
            call swapi(matsym(i,j,1),matsym(next,j,1))
            call swapi(matsym(i,j,2),matsym(next,j,2))
         end do
      end do
c
c Generating framework group string
c
      fwgrp=trim(pgrp)//'['
      do i=1,nsg+ncr+1
         if(nsym(msym(i,2),5).eq.0.or.msym(i,1).eq.0) cycle
         if(fwgrp(len_trim(fwgrp):len_trim(fwgrp)).ne.'[')
     &      fwgrp=trim(fwgrp)//','
         write(ch,'(i3)') msym(i,1)
         if(msym(i,1).gt.1) fwgrp=trim(fwgrp)//adjustl(ch)
         if(nsym(msym(i,2),1).eq.0) then
c O subspace
            fwgrp=trim(fwgrp)//'O('//
     &         trim(adjustl(symb(matsym(i,1,1))))//')'
         elseif(nsym(msym(i,2),1).eq.2) then
c C subspaces
            fwgrp=trim(fwgrp)//'C'
            if(nsym(msym(i,2),4).eq.2) fwgrp=trim(fwgrp)//''''
            if(nsym(msym(i,2),4).eq.3) fwgrp=trim(fwgrp)//'"'
            if(nsym(msym(i,2),2).gt.0) then
               write(ch,'(i3)') nsym(msym(i,2),2)
               fwgrp=trim(fwgrp)//adjustl(ch)
            else
               fwgrp=trim(fwgrp)//'inf'
            end if
            fwgrp=trim(fwgrp)//'('
            do j=1,10
               if(matsym(i,j,1).eq.0) cycle
               if(matsym(i,j,1).eq.-1) then
                  fwgrp=trim(fwgrp)//'.'
                  cycle
               end if
               if(fwgrp(len_trim(fwgrp):len_trim(fwgrp)).ne.'('.and.
     &            fwgrp(len_trim(fwgrp):len_trim(fwgrp)).ne.'.')
     &            fwgrp=trim(fwgrp)//','
               write(ch,'(i3)') matsym(i,j,2)
               fwgrp=trim(fwgrp)//
     &            trim(adjustl(symb(matsym(i,j,1))))
               if(matsym(i,j,2).gt.1) 
     &            fwgrp=trim(fwgrp)//trim(adjustl(ch))
            end do
            fwgrp=trim(fwgrp)//')'
         elseif(nsym(msym(i,2),1).eq.1) then
c SG subspaces
            fwgrp=trim(fwgrp)//'SG'
            if(nsym(msym(i,2),4).eq.1) then
               fwgrp=trim(fwgrp)//'H'
            elseif(nsym(msym(i,2),4).eq.2) then
               fwgrp=trim(fwgrp)//'V'
            elseif(nsym(msym(i,2),4).eq.3) then
               fwgrp=trim(fwgrp)//'D'
            end if
            fwgrp=trim(fwgrp)//'('
            do j=1,10
               if(matsym(i,j,1).eq.0) exit
               if(fwgrp(len_trim(fwgrp):len_trim(fwgrp)).ne.'(')
     &            fwgrp=trim(fwgrp)//','
               write(ch,'(i3)') matsym(i,j,2)
               fwgrp=trim(fwgrp)//
     &            trim(adjustl(symb(matsym(i,j,1))))
               if(matsym(i,j,2).gt.1) 
     &            fwgrp=trim(fwgrp)//trim(adjustl(ch))
            end do
            fwgrp=trim(fwgrp)//')'
         end if
      end do
      if(matsym(nsg+ncr+2,1,1).ne.0) then
c X subspace
         if(fwgrp(len_trim(fwgrp):len_trim(fwgrp)).ne.'[')
     &      fwgrp=trim(fwgrp)//','
         fwgrp=trim(fwgrp)//'X('
         do j=1,10
            if(matsym(i,j,1).eq.0) exit
            if(fwgrp(len_trim(fwgrp):len_trim(fwgrp)).ne.'(')
     &         fwgrp=trim(fwgrp)//','
               write(ch,'(i3)') matsym(i,j,2)
               fwgrp=trim(fwgrp)//
     &            trim(adjustl(symb(matsym(i,j,1))))
               if(matsym(i,j,2).gt.1) 
     &            fwgrp=trim(fwgrp)//trim(adjustl(ch))
         end do
         fwgrp=trim(fwgrp)//')'
      end if
      fwgrp=trim(fwgrp)//']'
c
      if(nout.ge.1) write(*,'(/a)') '-- FRAMEWORK GROUP --'
c
      if(nout.ge.1) then
         write(*,'(/)')
         do i=1,natoms
            k=natsym(i)
            symel='   '
            if(k.eq.0) then
               symel='X  '
               write(*,'(15x,a1,i3,5x,a2,5x,a3)') 
     &            '#',i,symb(nat(i)),symel
               cycle
            endif
            i1=nsym(k,1)
            i2=nsym(k,2)
            i3=nsym(k,3)
            i4=nsym(k,4)
            i5=nsym(k,5)
            if(nsym(k,1).eq.0) then
               symel='O  '
            elseif(nsym(k,1).eq.1.and.nsym(k,4).eq.0) then
               symel='SG '
            elseif(nsym(k,1).eq.1.and.nsym(k,4).eq.1) then
               symel='SGH'
            elseif(nsym(k,1).eq.1.and.nsym(k,4).eq.2) then
               symel='SGV'
            elseif(nsym(k,1).eq.1.and.nsym(k,4).eq.3) then
               symel='SGD'
            elseif(nsym(k,1).eq.2) then
               symel='C'//char(48+nsym(k,2))
               if(nsym(k,4).eq.2) then
                  symel=trim(symel)//"'"
               elseif(nsym(k,4).eq.3) then
                  symel=trim(symel)//'"'
               end if
            endif
            if(nsym(k,1).eq.1) then
               if(k.gt.0) then
                  write(*,'(15x,a1,i3,5x,a2,5x,a3,a,i3,a)') 
     &               '#',i,symb(nat(i)),symel,' (#',k,')'
               else
                  write(*,'(15x,a1,i3,5x,a2,5x,a3)') 
     &               '#',i,symb(nat(i)),symel
               endif
            else
               if(k.gt.0) then
                  write(*,'(15x,a1,i3,5x,a2,5x,a3,a,i3,a)') 
     &              '#',i,symb(nat(i)),symel,' (#',k,')'
               else
                  write(*,'(15x,a1,i3,5x,a2,5x,a3)') 
     &              '#',i,symb(nat(i)),symel
               endif
            endif
         end do
      end if
      if(nout.ge.1) write(*,'(/5x,2a)') 'Framework group: ',trim(fwgrp)
      return
      end
