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
      subroutine permutations(nat,coord,natoms,nsym,symn,nsg,ncr,nsr,
     &      nperm,nmax,nout)
************************************************************************
c
c Subroutine permutations determines permutations according to symmetry
c elements stored in nsym and symn arrays
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      character symb*2,op*7,nu*2
      logical symcen
      parameter(tol=1.d-2)
      dimension coord(3,nmax),nsym(150,5),symn(3,150),nperm(nmax,150),
     &      nat(nmax),o(3)
      common/data/ wt(90),symb(90)
      data o/0.d0,0.d0,0.d0/
      pi=4.d0*datan(1.d0)

c E and i
      if(nsym(1,2).eq.1) then
         do i=1,natoms
            nperm(i,1)=i
         end do
         ii=2
      else
         ii=1
      end if
      call inversion(natoms,nat,coord,tol,nc,nperm(1,ii),delta3)

c sigmas
      do i=2,nsg+1
         ii=ii+1
         call reflect(natoms,nat,coord,symn(1,i),o,tol,nc,nperm(1,ii),
     &         delta3)
      end do

c Cn axes
      do i=nsg+2,nsg+ncr+1
         ii=ii+1
         if(nsym(i,2).ne.0) then
            alpha=2.d0*pi/dble(nsym(i,2))
         else
            alpha=0.d0
         end if
         call rotate(natoms,nat,coord,symn(1,i),
     &         dsin(nsym(i,3)*alpha),dcos(nsym(i,3)*alpha),
     &         tol,nv,nperm(1,ii),delta3)
      end do

c Sn axes
      do i=nsg+ncr+2,nsg+ncr+nsr+1
         ii=ii+1
         if(nsym(i,2).ne.0) then
            alpha=2.d0*pi/dble(nsym(i,2))
         else
            alpha=0.d0
         end if
         call srotate(natoms,nat,coord,symn(1,i),
     &         dsin(nsym(i,3)*alpha),dcos(nsym(i,3)*alpha),
     &         tol,nv,nperm(1,ii),delta3)
      end do

c
      if(nout.lt.2) return
      symcen=nsym(1,2).eq.1
      write(*,'(/a)') '-- PERMUTATIONS --'
      write(*,'(/5x,a/)') '  # Operation   Permutation '
      write(*,'(5x,i3,2x,a1,6x,60(5(i6,a3,a3,i3,a3),:,/,17x))')
     &      1,adjustl('E'),
     &      (i,symb(nat(i)),' ->',i,symb(nat(i)),i=1,natoms)
      write(*,'()')
      if(symcen) then
         ii=1
      else
         ii=0
      end if
      do i=1,nsg+ncr+nsr+1
         if(i.eq.1.and.nsym(1,2).ne.1) cycle
         if(nsym(i,1).eq.1) then
            op='Sigma_'
            if(nsym(i,4).eq.1) then
               op=trim(op)//'h'
            elseif(nsym(i,4).eq.2) then
               op=trim(op)//'v'
            elseif(nsym(i,4).eq.3) then
               op=trim(op)//'d'
            end if
         elseif(nsym(i,1).eq.2) then
            write(nu,'(i2)') nsym(i,2)
            op='C'//adjustl(nu)
            if(nsym(i,3).gt.1) then
               write(nu,'(i2)') nsym(i,3)
               op=trim(op)//'^'//adjustl(nu)
            end if
            if(nsym(i,2).eq.2) then
               if(nsym(i,4).eq.2) then
                  op=trim(op)//"'"
               elseif(nsym(i,4).eq.3) then
                  op=trim(op)//'"'
               end if
            end if
         elseif(nsym(i,1).eq.3) then
            write(nu,'(i2)') nsym(i,2)
            op='S'//adjustl(nu)
            if(nsym(i,3).gt.1) then
               write(nu,'(i2)') nsym(i,3)
               op=trim(op)//'^'//adjustl(nu)
            end if
         elseif(nsym(i,1).eq.0.and.nsym(1,2).ne.1) then
            op='E'
         elseif(nsym(i,1).eq.0.and.nsym(1,2).eq.1) then
            op='i'
         end if
         write(*,'(5x,i3,2x,a7,60(5(i6,a3,a3,i3,a3),:,/,17x))') 
     &         i+ii,op,(j,symb(nat(j)),' ->',
     &         nperm(j,i+ii),symb(nat(nperm(j,i+ii))),j=1,natoms)
         write(*,'()')
      end do

      end
