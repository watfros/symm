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
      subroutine redrep(pgrp,nsymu,rep,rirrep)
************************************************************************
c
c Subroutine redrep reduces a representation
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      character pgsymb*3,irsymb*4,pgrp*3,repres*120,ch,ch2*2,
     &   trans*50,rotat*50,raman*50
      parameter(nmax=150)
      dimension rep(nmax),nsymu(nmax,5),rirrep(nmax)
      common/chartab/ nir(2,55),chtab(14,322),
     &   nsymop(14,4,55),nrotharm(3,322),pgsymb(57),irsymb(322)
c
      write(*,'(/a/)') '-- VIBRATIONS --'
      if(pgrp.eq.'Civ'.or.pgrp.eq.'Dih') then
         write(*,'(5x,a)') 'Not available.'
         return
      end if
      do ipg=1,55
         if(trim(pgrp).eq.trim(adjustl(pgsymb(ipg)))) exit
      end do
      if(ipg.eq.56) then
         write(*,'(a)') '   ERROR: Wrong point group.'
         return
      end if
      nord=0
      do i=1,14
         nord=nord+nsymop(i,4,ipg)
      end do
      do i=1,nir(2,ipg)
         do j=i,nir(2,ipg)
            if(nsymu(j,1).eq.nsymop(i,1,ipg)) then
               if(nsymu(j,1).eq.2.and.
     &            nsymu(j,2).eq.nsymop(i,2,ipg)) then
                  if(nsymu(j,2).eq.2) then
                     if(nsymu(j,4).eq.nsymop(i,3,ipg)) then
                        nch=j
                        exit
                     end if
                  else
                     if(nsymu(j,3).eq.nsymop(i,3,ipg)) then
                        nch=j
                        exit
                     end if
                  end if
               elseif(nsymu(j,1).eq.3.and.
     &            nsymu(j,2).eq.nsymop(i,2,ipg)) then
                  if(nsymu(j,3).eq.nsymop(i,3,ipg)) then
                     nch=j
                     exit
                  end if
               elseif(nsymu(j,1).eq.1) then
                  if(nsymu(j,4).eq.nsymop(i,2,ipg)) then
                     nch=j
                     exit
                  end if
               elseif(nsymu(j,1).eq.0.or.nsymu(j,1).eq.4) then
                  nch=j
                  exit
               end if
            end if
         end do
         if(j.eq.nir(2,ipg)+1) then
            write(*,'(/a)') '   ERROR: Couldn''t identify classes.'
            return
         end if
         if(nch.ne.i) then
            do j=1,5
               ntemp=nsymu(i,j)
               nsymu(i,j)=nsymu(nch,j)
               nsymu(nch,j)=ntemp
            end do
            tmp=rep(i)
            rep(i)=rep(nch)
            rep(nch)=tmp
         end if
      end do
      nf=0
      do i=1,nir(2,ipg)
         rirrep(i)=0.d0
         do j=1,nir(2,ipg)
            rirrep(i)=rirrep(i)+
     &         nsymop(j,4,ipg)*rep(j)*chtab(j,nir(1,ipg)+i-1)
         end do
         rirrep(i)=rirrep(i)/dble(nord)
         if((pgrp(1:1).eq.'C'.and.
     &         (pgrp(3:3).eq.' '.or.pgrp(3:3).eq.'h')).or.
     &         pgrp.eq.'T'.or.pgrp.eq.'Th'.or.pgrp(1:1).eq.'S') then
            if(scan(irsymb(nir(1,ipg)+i-1),'E').ne.0) 
     &            rirrep(i)=rirrep(i)/2.d0
         end if
         if(nint(rirrep(i)).ne.0.and.nf.eq.0) nf=i
      end do
c
      repres='        g'
      do i=1,nir(2,ipg)
         if(nsymop(i,1,ipg).eq.4) then
            repres=trim(repres)//'     E'
         elseif(nsymop(i,1,ipg).eq.0) then
            repres=trim(repres)//'      i'
         elseif(nsymop(i,1,ipg).eq.1) then
            repres=trim(repres)//'    SG'
            if(nsymop(i,2,ipg).eq.1) then
               repres=trim(repres)//'H'
            elseif(nsymop(i,2,ipg).eq.2) then
               repres=trim(repres)//'V'
            elseif(nsymop(i,2,ipg).eq.3) then
               repres=trim(repres)//'D'
            end if
         elseif(nsymop(i,1,ipg).eq.2) then
            if(nsymop(i,2,ipg).eq.2.and.nsymop(i,3,ipg).gt.1) then
               repres=trim(repres)//'    C'
            elseif(nsymop(i,3,ipg).gt.1) then
               repres=trim(repres)//'   C'
            else
               repres=trim(repres)//'     C'
            end if
            write(ch,'(i1)') nsymop(i,2,ipg)
            repres=trim(repres)//trim(ch)
            if(nsymop(i,2,ipg).eq.2) then
               if(nsymop(i,3,ipg).eq.2) then
                  repres=trim(repres)//"'"
               elseif(nsymop(i,3,ipg).eq.3) then
                  repres=trim(repres)//'"'
               end if
            elseif(nsymop(i,3,ipg).gt.1) then
               repres=trim(repres)//'^'
               write(ch,'(i1)') nsymop(i,3,ipg)
               repres=trim(repres)//trim(ch)
            end if
         elseif(nsymop(i,1,ipg).eq.3) then
            if(nsymop(i,3,ipg).eq.1) then
               if(nsymop(i,2,ipg).lt.10) then
                  repres=trim(repres)//'     S'
               else
                  repres=trim(repres)//'    S'
               end if
            elseif(nsymop(i,3,ipg).lt.10) then
               if(nsymop(i,2,ipg).lt.10) then
                  repres=trim(repres)//'   S'
               else
                  repres=trim(repres)//'  S'
               end if
            else
               if(nsymop(i,2,ipg).lt.10) then
                  repres=trim(repres)//'  S'
               else
                  repres=trim(repres)//' S'
               end if
            end if
            if(nsymop(i,2,ipg).lt.10) then
               write(ch,'(i1)') nsymop(i,2,ipg)
               repres=trim(repres)//trim(ch)
            else
               write(ch2,'(i2)') nsymop(i,2,ipg)
               repres=trim(repres)//trim(ch2)
            end if
            if(nsymop(i,3,ipg).gt.1) then
               repres=trim(repres)//'^'
               write(ch,'(i1)') nsymop(i,3,ipg)
               repres=trim(repres)//trim(ch)
            end if
         end if
      end do
c
      write(*,'(5x,a,/)') 
     &  'Reducible representation generated by the normal coordinates:'
      write(*,'(7x,a)') trim(repres)
      write(*,'(10x,a2,i4)',advance='no') 'gi',nord
      write(*,'(i6)',advance='no') nsymop(1,4,ipg)
      do i=2,nir(2,ipg)
         write(*,'(14(i7))',advance='no') nsymop(i,4,ipg)
      end do
      write(*,'()')
      write(*,'(7x)',advance='no')
      do i=0,nir(2,ipg)
         write(*,'(a7)',advance='no') '-------'
      end do
      write(*,'(a4)') '----'
      write(*,'(15x)',advance='no')
      do i=1,nir(2,ipg)
         write(*,'(14(f7.1))',advance='no') rep(i)
      end do
      write(*,'()')
      write(*,'(/,5x,a//,7x,i3,1x,a)',advance='no') 
     &   'Decomposed into irreps:',nint(rirrep(nf)),
     &   trim(adjustl(irsymb(nir(1,ipg)-1+nf)))
      do i=nf+1,nir(2,ipg)
         if(nint(rirrep(i)).ne.0)
     &      write(*,'(a,i3,1x,a)',advance='no') ' + ',nint(rirrep(i)),
     &   trim(adjustl(irsymb(nir(1,ipg)+i-1)))
      end do
      write(*,'()')
c
      trans=''
      rotat=''
      raman=''
      do i=1,nir(2,ipg)
         if(nrotharm(1,nir(1,ipg)-1+i).gt.0) then
            rirrep(i)=rirrep(i)-1.d0
            rotat=trim(rotat)//', '//
     &         trim(adjustl(irsymb(nir(1,ipg)-1+i)))
         end if
         if(nrotharm(2,nir(1,ipg)-1+i).gt.0) then
            rirrep(i)=rirrep(i)-1.d0
            trans=trim(trans)//', '//
     &         trim(adjustl(irsymb(nir(1,ipg)-1+i)))
         end if
         if(nrotharm(3,nir(1,ipg)-1+i).gt.0) then
            raman=trim(raman)//', '//
     &         trim(adjustl(irsymb(nir(1,ipg)-1+i)))
         end if
      end do
      write(*,'(/,5x,2a)') 'Irreps of translation: ',
     &   trans(3:len_trim(trans))
      write(*,'(5x,2a)') 'Irreps of  rotation:   ',
     &   rotat(3:len_trim(rotat))
      write(*,'(/,5x,a//,7x,i3,1x,a)',advance='no') 
     &   'Irreps of the vibrational normal coordinates:',
     &   nint(rirrep(nf)),trim(adjustl(irsymb(nir(1,ipg)-1+nf)))
      do i=nf+1,nir(2,ipg)
         if(nint(rirrep(i)).ne.0)
     &      write(*,'(a,i3,1x,a)',advance='no') ' + ',nint(rirrep(i)),
     &   trim(adjustl(irsymb(nir(1,ipg)+i-1)))
      end do
      write(*,'()')
c
      nvib=0
      do i=1,nir(2,ipg)
         nvib=nvib+nint(rirrep(i))
      end do
      ninfra=0
      nraman=0
      do i=1,nir(2,ipg)
         if(nrotharm(2,nir(1,ipg)-1+i).gt.0)
     &      ninfra=ninfra+nint(rirrep(i))
         if(nrotharm(3,nir(1,ipg)-1+i).gt.0)
     &      nraman=nraman+nint(rirrep(i))
      end do
      write(*,'(/,5x,a,i3,a)') 'There are ',nvib,
     &   ' different vibrational frequencies:'
      if(ninfra.eq.0) then
         write(*,'(/,10x,a)') '-  no  IR   active vibrations'
      elseif(ninfra.eq.1) then
         write(*,'(/,10x,a,i3,2a)') '- ',ninfra,
     &      '  IR   active vibration, with symmetry: ',
     &      trans(2:len_trim(trans))
      elseif(ninfra.gt.1) then
         write(*,'(/,10x,a,i3,2a)') '- ',ninfra,
     &      '  IR   active vibrations, with symmetry: ',
     &      trans(2:len_trim(trans))
      end if
      if(nraman.eq.0) then
         write(*,'(10x,a)') '-  no Raman active vibrations'
      elseif(nraman.eq.1) then
         write(*,'(10x,a,i3,2a)') '- ',nraman,
     &      ' Raman active vibration, with symmetry: ',
     &      raman(2:len_trim(raman))
      elseif(nraman.gt.1) then
         write(*,'(10x,a,i3,2a)') '- ',nraman,
     &      ' Raman active vibrations, with symmetry: ',
     &      raman(2:len_trim(raman))
      end if
      end
