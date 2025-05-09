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
      integer function element_number(element_name,symb)
       implicit none
       character(len=2)             element_name,symb
       dimension                    symb(90)
      !=====
       integer :: ielement
      !=====
      
       ielement=1
       do while( ADJUSTL(element_name) /= ADJUSTL(symb(ielement)) )
	     !write(*,*)  "element_name  ",ADJUSTL(element_name)
		 !write(*,*)  "symb  ",ADJUSTL(symb(ielement))
         if( ielement == 91 ) then
           write(*,*)    'element name error, element number > 90'
           stop
         endif
         ielement = ielement + 1
       enddo
      
       element_number = ielement
	   !write(*,*)    element_number
       
      
      end function element_number

************************************************************************
      subroutine sym_elements(natoms,nat,coord,symb,delta,norder,ni,nsg,
     &   ncr,nsr,np,symn,nsym,nout,nprm,nper,del_opt,nseq,nccl,nscl)
************************************************************************
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      character symb*2,symel*3
      logical symcen,linear,planar
      parameter(mat=250,nmax=150)
      dimension coord(3,mat),nat(mat),symb(90),ntrans(mat)
      dimension meq(mat),ieq(10,mat),sigman(nmax,3)
      dimension p0(3),p1(3),p2(3),p3(3)
      dimension v1(3),v2(3),v3(3),v0(3)
      dimension rotn(nmax,3),rota(nmax)
      dimension a(3),b(3),c(3)
      dimension symn(3,nmax),nsym(nmax,5)
      dimension nper(mat,mat),nscl(mat,mat),nccl(mat)
c
      delta2=0.d0
      delta3=0.d0
c
c Identical permutation (E operation)
c
      nprm=1
      do i=1,natoms
         nper(i,1)=i
      end do
c
c Generating the value of PI
c
      pi=4.d0*datan(1.d0)
c
c Partitioning the set of atoms into classes based on their 
c atomic numbers	
c
      neq=1
      meq(neq)=1
      ieq(neq,1)=1
      do i=2,natoms
         nati=nat(i)
         do j=1,neq
            natj=nat(ieq(j,1))
            if(nati.eq.natj) then
               meq(j)=meq(j)+1
               ieq(j,meq(j))=i
               goto 10
            endif
         end do
         neq=neq+1
         meq(neq)=1
         ieq(neq,1)=i
 10      continue
      end do
      if(nout.eq.2) then
         write(*,'(/a,i3)') '-- Equivalence classes of atoms: ',neq
         do i=1,neq
            write(*,'(/5x,a,i3,a7,a2,a1)') '#',i,
     &           ' (atom ',symb(nat(ieq(i,1))),')'
            write(*,'(5x,15i4)') (ieq(i,j),j=1,meq(i))
         end do
      end if
c
      symcen=.false.
      linear=.false.
      planar=.false.
      nsg=0
      nrot=0
c
c Centre of symmetry
c
      call inversion(natoms,nat,coord,delta,nc,nper(1,2),del)
      if(nc.eq.natoms) then
         symcen=.true.
         if(del.gt.delta3) delta3=del
         nprm=2
      end if
      icent=0
c
      symn(:,1)=0.d0
      nsym(1,:)=0
c
      if(symcen) then
         if(nout.ge.1) write(*,'(/a)') '-- CENTRE OF SYMMETRY: {i}'
         nsym(1,2)=1
      endif
      do i=1,natoms
         p0(1)=coord(1,i)
         p0(2)=coord(2,i)
         p0(3)=coord(3,i)
         sp=dsqrt(dot(p0,p0,3))
         if(sp.le.delta) then
            icent=i
            if(sp.gt.delta3) delta3=sp
         end if
      end do
      if(icent.ne.0.and.nout.eq.2) write(*,'(/a,a2,a,i3,a,a)') 
     &     '-- Atom ',symb(nat(icent)),' (',icent,')',
     &     ' in the COM'
      nsym(1,3)=icent
      if(icent.gt.0) nsym(1,5)=1
c
      do i=1,natoms-1
         if(i.eq.icent) cycle
         p1(1)=coord(1,i)
         p1(2)=coord(2,i)
         p1(3)=coord(3,i)
         do j=i+1,natoms
            if(j.eq.icent) cycle
            p2(1)=coord(1,j)
            p2(2)=coord(2,j)         
            p2(3)=coord(3,j)
            call crossp(p1,p2,p0)
            vn=dsqrt(dot(p0,p0,3))
            if(vn.gt.delta) goto 20
            if(vn.gt.delta3) delta3=vn
         end do
      end do
      linear=.true.
c
      if(nout.ge.1) write(*,'(/a)') '-- LINEAR MOLECULE'
c
      if(symcen) then
         if(nout.ge.1) write(*,'(/a,a6,a)') 
     &       '-- The structure should belong to the ','Dinf_h',
     &       ' point group.'
         if(nout.ge.1) write(*,'(/a)') '-- PLANES OF SYMMETRY --'
         nsg = 1
         nsym(2, 1) = 1
         nsym(2, 2) = 0
         nsym(2, 3) = 0
         nsym(2, 4) = 2
         nsym(2, 5) = natoms
         symn(:, 2) = 0.d0
         if(nout.ge.1) write(*,'(/a)') '-- Infinite planes' 
         if(nout.eq.2) write(*,'(5x,a)') 'All atoms included.'
         if(nout.ge.1) write(*,'(/a)')
     &      '-- Distinct PROPER ROTATIONAL AXES --'
         ncr = 2
         nsym(3, 1) = 2
         nsym(3, 2) = 0
         nsym(3, 3) = 1
         nsym(3, 4) = 1
         nsym(3, 5) = natoms
         if(icent.ne.1) then
            symn(:, 3) = coord(:, 1)
         else
            symn(:, 3) = coord(:, 2)
         end if
         symn(:, 3) = symn(:, 3) / dsqrt(dot(symn(:, 3), symn(:, 3), 3))
         nsym(4, 1) = 2
         nsym(4, 2) = 2
         nsym(4, 3) = 1
         nsym(4, 4) = 2
         nsym(4, 5) = 0
         symn(:, 4) = 0.d0
         if(nout.ge.1)
     &      write(*,'(/a,i3,a,f8.2,a)') '-- Axis #',1,': C(',0.d0,')'
         if(nout.eq.2) then
            write(*,'(2x,a,3f12.5)') ' d: ',(symn(k, 3),k=1,3)
            write(*,'(5x,a)') 'All atoms included.'
         end if
c write to a.tcl
	  write(1001,*) " variable axes { \"
	  write(1001,'(a,3f8.4,a,a)') 
     &            "{{",(symn(k, 3),k=1,3)," }"," -1 {} }} "
	  write(1001,*) " variable  planes { \"
	  write(1001,'(a,3f8.4,a,a)') 
     &            "{{",(symn(k, 3),k=1,3)," }"," horizontal }} "
	  write(1001,*) " variable rraxes  {}"
         if(nout.ge.1)
     &      write(*,'(/a,i3,a,f8.2,a)') '-- Axis #',2,': C(',180.d0,')'
         if(nout.eq.2) then
            write(*,'(2x,a,3f12.5)') ' d: ',(symn(k, 4),k=1,3)
            write(*,'(5x,a)') 'Atoms included: '
            if(icent.ne.0) write(*,'(10x,2a2,i3,a2)') 
     &      symb(nat(nsym(1,3))),' (',nsym(1,3),')'
         end if

         nsr = 1
         nsym(5, 1) = 3
         nsym(5, 2) = 0
         nsym(5, 3) = 1
         nsym(5, 4) = 0
         nsym(5, 5) = 1
         if(icent.ne.1) then
            symn(:, 5) = coord(:, 1)
         else
            symn(:, 5) = coord(:, 2)
         end if
         symn(:, 5) = symn(:, 5) / dsqrt(dot(symn(:, 5), symn(:, 5), 3))
         ni = 1
         norder = -1
         np = -1
         if(nout.ge.1) write(*,'(/a)') 
     &      '-- Number of symmetry operations = infinite' 
      else
         if(nout.ge.1) write(*,'(/a,a6,a)') 
     &       '-- The structure should belong to the ','Cinf_v',
     &       ' point group.'
         if(nout.ge.1) write(*,'(/a)') '-- PLANES OF SYMMETRY --'
         nsg = 1
         nsym(2, 1) = 1
         nsym(2, 2) = 0
         nsym(2, 3) = 0
         nsym(2, 4) = 2
         nsym(2, 5) = natoms
         symn(:, 2) = 0.d0
         if(nout.ge.1) write(*,'(/a)') '-- Infinite planes' 
         if(nout.eq.2) write(*,'(5x,a)') 'All atoms included.'
         if(nout.ge.1) write(*,'(/a)')
     &      '-- Distinct PROPER ROTATIONAL AXES --'
         ncr = 1
         nsym(3, 1) = 2
         nsym(3, 2) = 0
         nsym(3, 3) = 1
         nsym(3, 4) = 1
         nsym(3, 5) = natoms
         if(icent.ne.1) then
            symn(:, 3) = coord(:, 1)
         else
            symn(:, 3) = coord(:, 2)
         end if
         symn(:, 3) = symn(:, 3) / dsqrt(dot(symn(:, 3), symn(:, 3), 3))
         if(nout.ge.1)
     &      write(*,'(/a,i3,a,f8.2,a)') '-- Axis #',1,': C(',0.d0,')'
         if(nout.eq.2) then
            write(*,'(2x,a,3f12.5)') ' d: ',(symn(k, 3),k=1,3)
            write(*,'(5x,a)') 'All atoms included.'
         end if
c write to a.tcl
c write to a.tcl
	  write(1001,*) " variable axes { \"
	  write(1001,'(a,3f8.4,a,a)') 
     &            "{{",(symn(k, 3),k=1,3)," }"," -1 {} }} "
	  write(1001,*) " variable  planes {} "
	  write(1001,*) " variable rraxes  {}"
         nsr = 0
         norder = -1
         ni = 0
         np = -1
         if(nout.ge.1) write(*,'(/a)') 
     &      '-- Number of symmetry operations = infinite' 
      endif
      goto 100
 20   continue
      v1=p0/vn
      do i=1,natoms
         p3(1)=coord(1,i)
         p3(2)=coord(2,i)
         p3(3)=coord(3,i)
         sp=dot(v1,p3,3)
         if(dabs(sp).gt.delta) goto 30
         if(dabs(sp).gt.delta3) delta3=dabs(sp)
      end do
      planar=.true.
      nsg=nsg+1
      sigman(nsg,1)=v1(1)
      sigman(nsg,2)=v1(2)
      sigman(nsg,3)=v1(3)
      if(nout.ge.1) write(*,'(/a)') '-- PLANAR MOLECULE'
      if(nout.eq.2) write(*,'(2x,a,3f12.5)') ' n: ',(v1(k),k=1,3)
      if(symcen.and.planar) then
         do i=12,2,-1
            alpha=2.d0*pi/dble(i)
            sp=alpha*180.d0/pi
            sina=dsin(alpha)
            cosa=dcos(alpha)
            call rotate(natoms,nat,coord,v1,sina,cosa,delta,nc,ntrans,
     &         del)
            if(nc.eq.natoms) then
               call add_Cn(nrot,rotn,rota,v1,p3,sp,delta)
               call add_perm(natoms,ntrans,nprm,nper)
               if(del.gt.delta3) delta3=del
            end if
         end do
      endif         
 30   continue
c
c Planes of symmetry
c
      do i=1,neq
         meqi=meq(i)
         if(meqi.eq.1) cycle
         do j1=1,meqi-1
            i1=ieq(i,j1)
            p1(1)=coord(1,i1)
            p1(2)=coord(2,i1)
            p1(3)=coord(3,i1)
            do j2=j1+1,meqi
               i2=ieq(i,j2)
               p2(1)=coord(1,i2)
               p2(2)=coord(2,i2)
               p2(3)=coord(3,i2)
               p0=(p1+p2)/2.d0
               v1=p2-p0
               vn=dsqrt(dot(v1,v1,3))
               if(vn.gt.delta) then
                  v1=v1/vn
                  sp=dot(v1,-p0,3)
                  if(dabs(sp).lt.delta) then
                    call reflect(natoms,nat,coord,v1,p0,delta,nc,ntrans,
     &                  del)
                     if(nc.eq.natoms) then
                        if(del.gt.delta3) delta3=del
                        call add_SG(nsg,sigman,v1,p3,delta)
                        call add_perm(natoms,ntrans,nprm,nper)
                     end if
                  endif
               endif
               call crossp(p1,p2,v2)
               vn=dsqrt(dot(v2,v2,3))
               if(vn.gt.delta) then
                  v2=v2/vn
                  sp=dot(v2,-p0,3)
                  if(dabs(sp).lt.delta) then
                    call reflect(natoms,nat,coord,v2,p0,delta,nc,ntrans,
     &                  del)
                     if(nc.eq.natoms) then
                        if(del.gt.delta3) delta3=del
                        call add_SG(nsg,sigman,v2,p3,delta)
                        call add_perm(natoms,ntrans,nprm,nper)
                     end if
                  endif
               endif
               call crossp(v1,v2,v3)
               vn=dsqrt(dot(v3,v3,3))
               if(vn.gt.delta) then
                  v3=v3/vn
                  sp=dot(v3,-p0,3)
                  if(dabs(sp).lt.delta) then
                    call reflect(natoms,nat,coord,v3,p0,delta,nc,ntrans,
     &                  del)
                     if(nc.eq.natoms) then
                        if(del.gt.delta3) delta3=del
                        call add_SG(nsg,sigman,v3,p3,delta)
                        call add_perm(natoms,ntrans,nprm,nper)
                     end if
                  endif
               endif
            end do 
         end do
      end do
c
      if(nout.ge.1) write(*,'(/a)') '-- PLANES OF SYMMETRY --'
c
c write to a.tcl
	  write(1001,*) "variable  planes {  \"
      do i=1,nsg
         if(nout.ge.1) write(*,'(/a,i3)') '-- Plane #',i
         v1(1)=sigman(i,1)
         v1(2)=sigman(i,2)
         v1(3)=sigman(i,3)
         if(nout.eq.2) then
            write(*,'(2x,a,3f12.5)') ' n: ',(v1(k),k=1,3)
            write(*,'(5x,a)') 'Atoms included: '
         end if
c write to a.tcl
	    write(1001,'(a,3f8.4,a)') " {{",v1," } plane } \"
         m=0
         do j=1,natoms
            p3(1)=coord(1,j)
            p3(2)=coord(2,j)
            p3(3)=coord(3,j)
            sp=dot(v1,p3,3)
            if(dabs(sp).le.delta) then
               if(nout.eq.2)
     &            write(*,'(10x,2a2,i3,a2)') symb(nat(j)),' (',j,')'
               m=m+1
               if(dabs(sp).gt.delta3) delta3=dabs(sp)
            endif         
         end do
c Vector n
         symn(:,i+1)=v1
c SG
         nsym(i+1,1)=1
c Not used
         nsym(i+1,2)=0
         nsym(i+1,3)=0
c SGH=1, SGV=2, SGD=3 
         nsym(i+1,4)=0
c Number of atoms included 
         nsym(i+1,5)=m
      end do
c write to a.tcl
	  write(1001,*) "} "
c
      if(nout.ge.1) write(*,'(/a)') 
     &'-- Proper rotations due to the centre and planes of symmetry --'
c
      do i=1,nsg-1
         v1(1)=sigman(i,1)
         v1(2)=sigman(i,2)
         v1(3)=sigman(i,3)
         do j=i+1,nsg
            v2(1)=sigman(j,1)
            v2(2)=sigman(j,2)
            v2(3)=sigman(j,3)
            call crossp(v1,v2,v3)
            vn=dsqrt(dot(v3,v3,3))
            v3=v3/vn
            sp=dot(v1,v2,3)
            sp=dacos(sp)*180.d0/pi
            if(sp.gt.90.d0) sp=180.d0-sp
            sp=2.d0*sp
            call add_Cn(nrot,rotn,rota,v3,p3,sp,delta)
         end do
      end do
c
      do i=1,nrot
         m=0
         sp=rota(i)
         if(nout.ge.1) write(*,'(/a,i3,a,f8.2,a)')
     &      '-- Rotation #',i,': C(',sp,')'
         v1(1)=rotn(i,1)
         v1(2)=rotn(i,2)
         v1(3)=rotn(i,3)
         if(nout.eq.2) then
            write(*,'(2x,a,3f12.5)') ' d: ',(v1(k),k=1,3)
            write(*,'(5x,a)') 'Atoms included: '
         end if
         do j=1,natoms
            p3(1)=coord(1,j)
            p3(2)=coord(2,j)
            p3(3)=coord(3,j)                                  
            v2=p0+v1
            call crossp(v2,p3,v0)
            vn=dsqrt(dot(v0,v0,3))
            if(dabs(vn).le.delta) then
               if(nout.eq.2)
     &            write(*,'(10x,2a2,i3,a2)') symb(nat(j)),' (',j,')'
               m=m+1
               if(dabs(vn).gt.delta3) delta3=dabs(vn)
            end if
         end do
      end do
c
c
c Proper rotational axes
c
      if(nout.ge.1) write(*,'(/a)')
     &   '-- Distinct PROPER ROTATIONAL AXES --'
c
c Cn (for each atom)
c
      do i=1,neq
         meqi=meq(i)
         do j=1,meqi
            i1=ieq(i,j)
            p0(1)=coord(1,i1)
            p0(2)=coord(2,i1)
            p0(3)=coord(3,i1)
            vn=dsqrt(dot(p0,p0,3))
            if(vn.lt.delta) cycle
            v1=p0/vn
            do k=12,2,-1
               alpha=2.d0*pi/dble(k)
               sp=alpha*180.d0/pi
               sina=dsin(alpha)
               cosa=dcos(alpha)
c
            call rotate(natoms,nat,coord,v1,sina,cosa,delta,nc,ntrans,
     &         del)
               if(nc.eq.natoms)  then
                  call add_Cn(nrot,rotn,rota,v1,p3,sp,delta)
                  call add_perm(natoms,ntrans,nprm,nper)
                  if(del.gt.delta3) delta3=del
               end if
            end do
         end do
      end do
c
c Cn (for each pair of atoms)
c
      do i=1,neq
         meqi=meq(i)
         if(meqi.lt.2) cycle
         do j1=1,meqi-1
            i1=ieq(i,j1)
            p1(1)=coord(1,i1)
            p1(2)=coord(2,i1)
            p1(3)=coord(3,i1)
            do j2=j1+1,meqi
               i2=ieq(i,j2)
               p2(1)=coord(1,i2)
               p2(2)=coord(2,i2)
               p2(3)=coord(3,i2)
               p0=(p1+p2)/2.d0
               vn=dsqrt(dot(p0,p0,3))
               if(vn.lt.delta) cycle
               v1=p0/vn
               alpha=pi
               sp=180.d0
               sina=dsin(alpha)
               cosa=dcos(alpha)
              call rotate(natoms,nat,coord,v1,sina,cosa,delta,nc,ntrans,
     &            del)
               if(nc.eq.natoms) then
                  call add_Cn(nrot,rotn,rota,v1,p3,sp,delta)
                  call add_perm(natoms,ntrans,nprm,nper)
                  if(del.gt.delta3) delta3=del
               end if
            end do
         end do
      end do
c
c Cn (n > 2)
c
      do i=1,neq
         meqi=meq(i)
         if(meqi.lt.3) cycle
         do j1=1,meqi-2
            i1=ieq(i,j1)
            a(1)=coord(1,i1)
            a(2)=coord(2,i1)
            a(3)=coord(3,i1)
            do j2=j1+1,meqi-1
               i2=ieq(i,j2)
               b(1)=coord(1,i2)
               b(2)=coord(2,i2)
               b(3)=coord(3,i2)
               p1=(a+b)/2.d0
               p3=b-p1
               vn=dsqrt(dot(p3,p3,3))
               v1=p3/vn
               do j3=j2+1,meqi
                  i3=ieq(i,j3)
                  c(1)=coord(1,i3)
                  c(2)=coord(2,i3)
                  c(3)=coord(3,i3)
                  p2=(b+c)/2.d0
                  p3=c-p2
                  vn=dsqrt(dot(p3,p3,3))
                  v2=p3/vn
                  call crossp(v1,v2,v3)
                  vn=dsqrt(dot(v3,v3,3))
                  if(vn.lt.delta) cycle
                  v3=v3/vn
                  sp=dacos(dot(v1,v2,3))
                  if(dabs(sp).lt.delta) cycle 
                  m=idint(2.d0*pi/sp+delta)
                  if((m*sp).lt.(2.d0*pi-delta)) cycle 
                  if((m.lt.3).or.(m.gt.meqi)) cycle
                  alpha=sp
                  sp=alpha*180.d0/pi
                  sina=dsin(alpha)
                  cosa=dcos(alpha)
                  call rotate(natoms,nat,coord,v3,sina,cosa,delta,nc,
     &               ntrans,del)
                  if(nc.eq.natoms) then
                     call add_Cn(nrot,rotn,rota,v3,p3,sp,delta)
                     call add_perm(natoms,ntrans,nprm,nper)
                     if(del.gt.delta3) delta3=del
                  end if
               end do
            end do
         end do
      end do
c
      do i=1,nrot
         m=0
         sp=rota(i)
         if(nout.ge.1)
     &      write(*,'(/a,i3,a,f8.2,a)') '-- Axis #',i,': C(',sp,')'
         v1(1)=rotn(i,1)
         v1(2)=rotn(i,2)
         v1(3)=rotn(i,3)
         if(nout.eq.2) then
            write(*,'(2x,a,3f12.5)') ' d: ',(v1(k),k=1,3)
            write(*,'(5x,a)') 'Atoms included: '
         end if
         do j=1,natoms
            p3(1)=coord(1,j)
            p3(2)=coord(2,j)
            p3(3)=coord(3,j)                                  
            call crossp(v1,p3,v0)
            vn=dsqrt(dot(v0,v0,3))
            if(dabs(vn).le.delta) then
               if(nout.eq.2)
     &            write(*,'(10x,2a2,i3,a2)') symb(nat(j)),' (',j,')'
               m=m+1
               if(dabs(vn).gt.delta3) delta3=dabs(vn)
            end if
         end do
      end do
c
      if(nout.ge.1)
     &   write(*,'(/a/)') '-- PROPER ROTATIONAL AXES & ROTATIONS --'
c write to a.tcl
	  write(1001,*) "variable axes {  \"
c
      nsgi=nsg+1
      ii=0
      do i=1,nrot
         v1(1)=rotn(i,1)
         v1(2)=rotn(i,2)
         v1(3)=rotn(i,3)
         do k=12,2,-1
            alpha=2.d0*pi/dble(k)
            sina=dsin(alpha)
            cosa=dcos(alpha)
            call rotate(natoms,nat,coord,v1,sina,cosa,delta,nc,ntrans,
     &         del)
            if(nc.eq.natoms) then
               if(del.gt.delta3) delta3=del
               call add_perm(natoms,ntrans,nprm,nper)
               ii=ii+1
               m=0
               do j=1,natoms
                  p3(1)=coord(1,j)
                  p3(2)=coord(2,j)
                  p3(3)=coord(3,j)                                  
                  call crossp(v1,p3,v0)
                  vn=dsqrt(dot(v0,v0,3))
                  if(dabs(vn).le.delta) then
                     m=m+1 
                     if(dabs(vn).gt.delta3) delta3=dabs(vn)
                  end if
               end do
c Vector n
               symn(:,nsgi+ii)=v1
c Cn
               nsym(nsgi+ii,1)=2
c n-fold (n ^ 1) 
               nsym(nsgi+ii,2)=k
               nsym(nsgi+ii,3)=1
c Principal axis
               nsym(nsgi+ii,4)=0
c Number of atoms included (unmoved atoms) 
               nsym(nsgi+ii,5)=m
c
               if(nout.ge.1) write(*,'(a,i3,a,i3,a,i2,a)')
     &              '-- #',i,'-',ii,': C(',k,')'
c write to a.tcl
	           write(1001,'(a,3f8.4,a,i4,a)') 
     &            "{{",v1," }",k," axis } \"
			   
               if(k.gt.2) then
                  do kk=2,k-1
                     alpha2=dble(kk)*alpha
                     sina2=dsin(alpha2)
                     cosa2=dcos(alpha2)
                     call rotate(natoms,nat,coord,v1,sina2,cosa2,delta,
     &                  nc2,ntrans,del)
                     if(nc2.eq.natoms) then
                        if(del.gt.delta3) delta3=del
                        call add_perm(natoms,ntrans,nprm,nper)
                        ii=ii+1
c Vector n
                        symn(:,nsgi+ii)=v1
c Cn
                        nsym(nsgi+ii,1)=2
c Principal axis
                        nsym(nsgi+ii,4)=0
c Number of atoms included 
                        nsym(nsgi+ii,5)=m
c
                        ngcd=igcd(k,kk)
                        nsym(nsgi+ii,2)=k/ngcd
                        nsym(nsgi+ii,3)=kk/ngcd
                        if(nout.ge.1) write(*,'(a,i3,a,i3,a,i2,a,i2,a)') 
     &                    '-- #',i,'-',ii,': C(',k,' ^',kk,')'
                     endif
                  end do
                  exit
               endif
            endif
         end do
      end do
c write to a.tcl
	  write(1001,*) "} "
      ncr=ii
c
      if(nout.ge.1)
     &   write(*,'(/a/)') '-- IMPROPER ROTATIONAL AXES & ROTATIONS --'
c write to a.tcl
	  write(1001,*) "variable rraxes  {  \"
      nsgicn=nsgi+ncr
      ii=0
      do i=1,nrot
         v1(1)=rotn(i,1)
         v1(2)=rotn(i,2)
         v1(3)=rotn(i,3)
         do k=24,2,-1
            alpha=2.d0*pi/dble(k)
            sina=dsin(alpha)
            cosa=dcos(alpha)
            call srotate(natoms,nat,coord,v1,sina,cosa,delta,nc,ntrans,
     &         del)
            if(nc.eq.natoms.and.k.gt.2) then
               if(del.gt.delta3) delta3=del
               call add_perm(natoms,ntrans,nprm,nper)
               ii=ii+1
               m=0
               if(icent.gt.0) m=1
c Vector n
               symn(:,nsgicn+ii)=v1
c Sn
               nsym(nsgicn+ii,1)=3
c n-fold (n ^ 1) 
               nsym(nsgicn+ii,2)=k
               nsym(nsgicn+ii,3)=1
c
               nsym(nsgicn+ii,4)=0
c Number of atoms included 
               nsym(nsgicn+ii,5)=m
c
               if(nout.ge.1) write(*,'(a,i3,a,i3,a,i2,a)') 
     &                 '-- #',i,'-',ii,': S(',k,')'
c write to a.tcl
	           write(1001,'(a,3f8.4,a,i4,a)') "{{",v1," }",k,
     &			   " r_axis } \"
               kv=k-1
               if(mod(k,2).ne.0) kv=2*k-1
               do kk=2,kv
                  alpha2=dble(kk)*alpha
                  sina2=dsin(alpha2)
                  cosa2=dcos(alpha2)
                  call srotate(natoms,nat,coord,v1,sina2,cosa2,delta,
     &               nc2,ntrans,del)
                  if((nc2.eq.natoms).and.(mod(kk,2).ne.0).and.
     &               (kk.ne.k).and.(igcd(k,kk).eq.1)) then
                     if(del.gt.delta3) delta3=del
                     call add_perm(natoms,ntrans,nprm,nper)
                     ii=ii+1
                     symn(:,nsgicn+ii)=v1
                     nsym(nsgicn+ii,1)=3
                     nsym(nsgicn+ii,5)=m
                     nsym(nsgicn+ii,2)=k
                     nsym(nsgicn+ii,3)=kk
                     if(nout.ge.1) write(*,'(a,i3,a,i3,a,i2,a,i2,a)') 
     &                   '-- #',i,'-',ii,': S(',k,'^',kk,')'
                  endif
               end do
            endif
         end do
      end do
c write to a.tcl
	  write(1001,*) "} "
      nsr=ii
      ni=0
      if(symcen) ni=1 
      norder=1+ni+nsg+ncr+nsr 
      if(nout.ge.1) write(*,'(/a,i5)') 
     &'-- Number of symmetry operations (including E) = ',norder 
c
c Determination of the principal axis
c
      np=0
      do i=nsg+2,nsg+ncr+1 
         if((nsym(i,2).gt.np).and.(nsym(i,3).eq.1)) then
            np=nsym(i,2)
         end if
      end do
      do i=nsg+2,nsg+ncr+1 
          if((nsym(i,2).eq.np).and.(nsym(i,3).eq.1)) nsym(i,4)=1
      end do
c
c Rotation axes: principal axes & orthogonal C2 axes
c
      nnp=0
      do i=nsg+2,nsg+ncr+1 
         if((nsym(i,4).eq.1)) nnp=nnp+1
      end do
      if(nnp.eq.3.and.np.eq.2) then
         do i=nsg+2,nsg+ncr+1 
            nsym(i,4)=2
         end do
c D2, D2d, D2h
         if(nsg.eq.3) then
c D2h
            nm=2
            do i=2,nsg+1 
               if((nsym(i,5).gt.nsym(nm,5))) nm=i
            end do
            nsym(nm,4)=1
            v2=symn(:,nm)          
            do i=nsg+2,nsg+ncr+1 
               v1=symn(:,i)
               vk=dot(v1,v2,3)
               if(dabs(dabs(vk)-1.d0).le.delta) then
                  nsym(i,4)=1
                  if(dabs(dabs(vk)-1.d0).gt.delta2)
     &               delta2=dabs(dabs(vk)-1.d0)
               end if
            end do
         elseif(nsg.eq.2) then
c D2d
            do i=nsg+2,nsg+ncr+1 
               v1=symn(:,i)
               do j=nsg+ncr+2,nsg+ncr+nsr+1 
                  v2=symn(:,j)           
                  vk=dot(v1,v2,3)
               if(dabs(dabs(vk)-1.d0).le.delta) then
                  nsym(i,4)=1
                  if(dabs(dabs(vk)-1.d0).gt.delta2)
     &               delta2=dabs(dabs(vk)-1.d0)
               end if
               end do
            end do
         endif
      endif
c
      do i=nsg+2,nsg+ncr+1 
         if(nsym(i,4).ne.1) cycle
         v1=symn(:,i)
         do j=nsg+2,nsg+ncr+1 
            if((nsym(j,2).eq.2).and.(nsym(j,3).eq.1)) then
               v2=symn(:,j)           
               vk=dot(v1,v2,3)
               if(dabs(vk).lt.delta) then
                  nsym(j,4)=2
                  if(dabs(vk).gt.delta2) delta2=dabs(vk)
               end if
            endif
         end do
      end do
c
c Perpendicular C2 axes
c
      maxcn=0
      if(ncr.gt.0) then
         do i=nsg+2,nsg+ncr+1 
            if(nsym(i,4).ne.2) cycle
            if(nsym(i,5).gt.maxcn) maxcn=nsym(i,5)
         end do
         do i=nsg+2,nsg+ncr+1 
            if(nsym(i,4).ne.2) cycle
            if(nsym(i,5).lt.maxcn) nsym(i,4)=3
         end do
      endif
c
c Planes of symmetry: SGH, SGV, SGD
c
      if(nsg.gt.0) then
c SGH
         do i=nsg+2,nsg+ncr+1 
            if(nsym(i,4).ne.1) cycle
            v1=symn(:,i)
            do j=2,nsg+1
               v2=symn(:,j)           
               vk=dot(v1,v2,3)
               if(dabs(dabs(vk)-1.d0).le.delta) then
                  nsym(j,4)=1
                  if(dabs(dabs(vk)-1.d0).gt.delta2)
     &               delta2=dabs(dabs(vk)-1.d0)
               end if
            end do
         end do
c SGV
         do i=2,nsg+1 
            if(nsym(i,4).eq.1) cycle
            v1=symn(:,i)
            do j=nsg+2,nsg+ncr+1 
               if(nsym(j,4).ne.1) cycle
               v2=symn(:,j)           
               vk=dot(v1,v2,3)
               if(dabs(vk).lt.delta) then
                  nsym(i,4)=2
                  if(dabs(vk).gt.delta2) delta2=dabs(vk)
               end if
            end do
         end do
         maxsg=0
         do i=2,nsg+1
            if(nsym(i,4).ne.2) cycle
            if(nsym(i,5).gt.maxsg) maxsg=nsym(i,5)
         end do
         do i=2,nsg+1
            if(nsym(i,4).ne.2) cycle
            if(nsym(i,5).lt.maxsg) nsym(i,4)=3
         end do
      endif
 100  continue
c
      if(nout.ge.1) write(*,'(/a/)') '-- SYMMETRY OPERATIONS --'
c
c COM and Inversion Center
c
      if(nout.eq.2) then
         if(nsym(1,2).eq.0) then
            if(nsym(1,3).gt.0) then
               write(*,'(15x,a,i3,a,a2,a,i3,a)') 
     &              '#',1,': COM    -- with atom ',symb(nat(nsym(1,3))),
     &              ' (#',nsym(1,3),')'
            else
               write(*,'(15x,a,i3,a)') 
     &              '#',1,': COM'
            endif
         elseif(nsym(1,2).eq.1) then     
            if(nsym(1,3).gt.0) then
               write(*,'(15x,a,i3,a,a2,a,i3,a)') 
     &              '#',1,': INVERSION CENTER  -- with atom ',
     &              symb(nat(nsym(1,3))),' (#',nsym(1,3),')'
            else
               write(*,'(15x,a,i3,a)') 
     &              '#',1,': INVERSION CENTER '
            endif
         endif
c
c SG
c
         do k=2,nsg+1
            if(nsym(k,1).eq.1.and.nsym(k,4).eq.0) then
               symel='SG '
            elseif(nsym(k,1).eq.1.and.nsym(k,4).eq.1) then
               symel='SGH'
            elseif(nsym(k,1).eq.1.and.nsym(k,4).eq.2) then
               symel='SGV'
            elseif(nsym(k,1).eq.1.and.nsym(k,4).eq.3) then
               symel='SGD'
            endif
            if(nsym(k,5).eq.0) then
               write(*,'(15x,a,i3,a,a3)') 
     &              '#',k,': ',symel
            else
               write(*,'(15x,a,i3,a,a3,a,i3,a)') 
     &       '#',k,': ',symel,'     -- with ',nsym(k,5),' unmoved atoms'
      
            endif
         end do
c
c C(n^k)
c 
         do k=nsg+2,nsg+ncr+1
            if(nsym(k,3).gt.1) then
               write(*,'(15x,a,i3,a,i2,a,i2,a)') 
     &              '#',k,': C(',nsym(k,2),'^',
     &              nsym(k,3),')'
            else
               if(nsym(k,5).gt.0) then
                  if(nsym(k,4).eq.2) then
                     write(*,'(15x,a,i3,a,i2,a,a,i3,a)') 
     &                    '#',k,': C''(',nsym(k,2),')',
     &                    '  -- with ',nsym(k,5),' unmoved atoms'
                  elseif(nsym(k,4).eq.3) then
                     write(*,'(15x,a,i3,a,i2,a,a,i3,a)') 
     &                    '#',k,': C"(',nsym(k,2),')',
     &                    '  -- with ',nsym(k,5),' unmoved atoms'
                  else
                     write(*,'(15x,a,i3,a,i2,a,a,i3,a)') 
     &                    '#',k,': C(',nsym(k,2),')',
     &                    '   -- with ',nsym(k,5),' unmoved atoms'
                  endif
               else
                  if(nsym(k,4).eq.2) then
                     write(*,'(15x,a,i3,a,i2,a)') 
     &                    '#',k,': C''(',nsym(k,2),')'
                  elseif(nsym(k,4).eq.3) then
                     write(*,'(15x,a,i3,a,i2,a,a,i3,a)') 
     &                    '#',k,': C"(',nsym(k,2),')'
                  else
                     write(*,'(15x,a,i3,a,i2,a)') 
     &                    '#',k,': C(',nsym(k,2),')'
                  endif
               endif
            endif
         end do
c
c S(n^k)
c 
         do k=nsg+ncr+2,nsg+ncr+nsr+1
            if(nsym(k,3).gt.1) then
               write(*,'(15x,a,i3,a,i2,a,i2,a)') 
     &              '#',k,': S(',nsym(k,2),'^',
     &              nsym(k,3),')'
            else
               if(nsym(k,5).gt.0) then
                  write(*,'(15x,a,i3,a,i2,a,a,i3,a)') 
     &                 '#',k,': S(',nsym(k,2),')',
     &                 '   -- with ',nsym(k,5),' unmoved atoms'
               else
                  write(*,'(15x,a,i3,a,i2,a)') 
     &                 '#',k,': S(',nsym(k,2),')'
               endif
            endif
         end do
      end if
c
c Number of symmetry operations
c
      if(nout.ge.1) then
         write(*,'(/5x,a)') 
     &   '       g       E       i      SG      Cn      Sn' 
         write(*,'(8x,a)') 
     &   '-------------------------------------------------'
         write(*,'(10x,i3,5x,i3,4(5x,i3))')
     &        norder,1,ni,nsg,ncr,nsr
c
c Determination of symmetry equivalence classes
c
      end if
      call symclass(natoms,nprm,nper,nseq,nccl,nscl,nat,symb,nout)
c
      if(nout.ge.1) then
         write(*,'(//5x,a,f18.8)') 
     &      'Distorsion of the found symmetry elements: ',delta2
         write(*,'(5x,a,f12.8)') 
     &      'Distorsion of geometry due to symmetry elements: ',delta3
      end if
c write to a.tcl
      write(1001,'(a,F8.4)')  " variable rmsd ",delta3
      del_opt=dmax1(delta2,delta3)
      return
      end

