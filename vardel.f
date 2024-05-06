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
      subroutine var_delta(natoms,nat,coord,symb,
     &   del_from,del_to,del_opt1,del_opt2, ii, pgrps, range,nout)
************************************************************************
c
c Subroutine var_delta determines largest symmetry
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      character symb*2
      character pgroup*3,pgmax*3,pgprev*3,pgrps*3
      character pgsymb*3, irsymb*4
      common/chartab/ nir(2,55),chtab(14,322),
     &   nsymop(14,4,55),nrotharm(3,322),pgsymb(57),irsymb(322)
      common/subgroups/ nsgb(2,57),nsgr(406)
      parameter(mat=250,nmax=150)
      dimension nat(mat),coord(3,mat),symb(90)
      dimension symn(3,nmax),nsym(nmax,5)
      dimension pgrps(20),range(2,20)
      dimension nper(mat,mat),nscl(mat,mat),nccl(mat)
      ngmax=0
      pgprev='   '
      ii=0
      if(del_from.lt.del_to) call swap(del_from,del_to)
c
      if(nout.gt.0) write(*,'(/a)') '-- VARYING TOLERANCE --'
      if(nout.gt.1) write(*,'(/,5x,a,f12.8,a,f12.8/)') 
     &   'Scanning range: ',del_to,' < delta < ',del_from
c
      do while(del_from.ge.del_to)
         call sym_elements(natoms,nat,coord,symb,del_from,ng,ni,nsg,
     &      ncr,nsr,np,symn,nsym,0,nprm,nper,delta,nseq,nccl,nscl)
        call point_group(ng,ni,nsg,ncr,nsr,np,pgroup,0)
         if(len_trim(pgroup).eq.0) then
            del_from=0.97d0*del_from
            pgprev='   '
            cycle
         end if
         if(pgroup.ne.pgprev) then
            ii=ii+1
            pgrps(ii)=pgroup
            range(1,ii)=del_from
            range(2,ii)=delta
            pgprev=pgroup
         elseif(pgroup.eq.pgprev) then
            range(2,ii)=delta
         end if
         if(iabs(ng).gt.ngmax) then
            ngmax=iabs(ng)
            pgmax=pgroup
         end if
         del_from=delta*0.98d0
      end do
c
      do i=1,ii
         if(trim(pgrps(i)).eq.trim(pgmax)) then
            del_opt1=range(2,i)
            del_opt2=range(1,i)
         end if
      end do
c
      do i=1,ii
        if(nout.gt.0) write(*,'(5x,2a,f12.8,a,f12.8)') pgrps(i),
     &      ' point group found in range ',range(2,i),' < delta < ',
     &      range(1,i)
      end do
      if(nout.gt.0) write(*,'(/5x,3a)') 'Largest symmetry should be ',
     &      trim(pgmax),'.'
c
c Write subgroups
c
      do npg = 1, 57
         if(pgmax .eq. adjustl(pgsymb(npg))) exit
      end do
      if(nout.gt.1) then
          write(*, '(/5x,a)') 'Available subgroups:'
          write(*,'(5(10x,20a,/))',advance='no') (pgsymb(nsgr(i)),
     &      ',  ', i = nsgb(1, npg), nsgb(1, npg) + nsgb(2, npg) - 2)
          write(*,'(a)')
     &      pgsymb(nsgr(nsgb(1,npg)+nsgb(2,npg)-1))
      end if
c
      end
