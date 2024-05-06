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

      subroutine gaussian_output(natoms,coord,nat,nmax,fname)
c Writes a Gaussian input file
      implicit double precision(a-h, o-z)
      implicit integer(i-n)
      character fname*100,lot*30
      dimension coord(3,nmax),nat(nmax)

c default level of theory, charge, multiplicity
      data lot/'RHF/6-31G(d)'/,kharge/0/,mult/1/

c open Gaussian input file
      open(15,file=trim(fname)//'.gjf')

      write(15,'(a)') '#'//trim(lot)
      write(15,'(/a)') fname
      write(15,'(/2i3)') kharge,mult
      do i=1,natoms
         write(15,'(i3,3f12.6)') nat(i),(coord(j,i),j=1,3)
      end do
      write(15,*)

      close(15)

      end
