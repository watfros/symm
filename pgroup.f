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
      subroutine point_group(ngp,ni,nsg,ncr,nsr,np,pgrp,nout)
************************************************************************
c
c Data for 55 Point Groups
c
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      character sg*3,cg*56,pgrp*3
      character S1, S2, Sn*2, tot*4
      logical yesno
      dimension sg(57),ng(57,6),cg(57)
      pgrp='   '
      nf=0
c 使用假循环避免多次嵌套

      Do i=1,1
c   Dih or Civ
        if(np==-1) then 
            if(ni==1) then 
                pgrp='Dih'
              exit
              else 
                pgrp='Civ'
              exit
	        endif
        endif
        
c   Ci or Cs
        if(np==0) then 
            if(nsg==1) then
              pgrp='Cs'
              exit
            endif
            if(ni==1) then
              pgrp='Ci'
              exit
            endif
        endif
    
c   Ih
        if(ngp==24*np) then 
            pgrp='Ih'
            exit
        endif
        
c   I or Oh
        if(ngp==12*np) then 
            if(np==5) then
              pgrp='I'
              exit
            endif
            if(np==4) then
              pgrp='Oh'
              exit
            endif
        endif

c   O
        if(ngp==6*np) then 
            pgrp='O'
            exit
        endif

c   Th or Td
        if(ngp==8*np) then 
            if(ni==1) then 
                pgrp='Th'
              exit
              else 
                pgrp='Td'
              exit
	    endif
        endif

c   T
        if(ngp==4*np .and. ni==0 .and. nsg==0) then 
            pgrp='T'
            exit
        endif

c   Dnh or Dnd
        if(ngp==4*np) then 
            if(nsg==np) then 
                write(pgrp,'(I2)') np
                pgrp="D"//trim(adjustl((pgrp)))//"d"
                exit
              else 
                write(pgrp,'(I2)') np
                pgrp='D'//trim(adjustl((pgrp)))//'h'
                exit
	      endif
        endif
c       另一个程序中有点群命名的方法，大致如下
c       write(Sn,'(I2)') np
c       tot=S1//Sn//S2
c       call compatta(tot,4,k)
c       pgrp = tot(1:3)

c   Dn
        if(ngp==2*np .and. ni==0 .and. nsg==0 .and. nsr==0) then 
            write(pgrp,'(I2)') np
            pgrp='D'//trim(adjustl((pgrp)))
            exit
        endif

c   Sn, n=np*2
        if(ngp==2*np .and. nsg==0 .and. ncr==np-1) then 
            write(pgrp,'(I2)') 2*np
            pgrp='S'//trim(adjustl((pgrp)))
            exit
        endif

c   Cnh
        if(ngp==2*np .and. nsg==1) then 
            write(pgrp,'(I2)') np
            pgrp='C'//trim(adjustl((pgrp)))//'h'
            exit
        endif
        
c   Cnv
        if(ngp==2*np .and. nsg==np) then 
            write(pgrp,'(I2)') np
            pgrp='C'//trim(adjustl((pgrp)))//'v'
            exit
        endif
        
c   Cn
        if(ngp==np) then 
            write(pgrp,'(I2)') np
            pgrp='C'//trim(adjustl((pgrp)))
            exit
        endif
		
c   C1
        if(ngp==1) then 
            pgrp='C1'
            exit
        endif
		
        write(*,'(/a)') 'Failed to determine the point group'
        exit
      enddo
      
	  write(*,*) '  ngp,ni ,nsg,ncr,nsr,np'
	  write(*,'(6I4)')  ngp,ni,nsg,ncr,nsr,np

c write to a.tcl
      write(1001,'(a,a)')  " variable pointgroup  ",pgrp
      write(1001,'(a,I4)')  " variable inversion",ni
      write(1001,'(a,I3,a,I3,a,I3,a,I3,a)')  
     & " variable elements {",ni,"*(i)",nsg,"*(sg)",ncr,"*(C)"
     & ,nsr,"*(S)}"
      write(*,'(/a)') '---- POINT GROUP ----'
      write(*,*) pgrp
      return
      end
