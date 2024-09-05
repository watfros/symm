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
c   
c   modified by Jinglin Mu(mjlink@126.com), to align with the VMD.
************************************************************************
      program main
************************************************************************
*                                                                      *
*   Program SYMMETRY determines all the symmetry elements & symmetry   *
*                    operations of rigid molecules.                    *
*           (C) G. Tasi, L. Gyevi-Nagy, R. Tobias, T.S. Tasi           *
*               J. Math. Chem., 51, 2187-2195 (2013)                   *
*         Department of Applied and Environmental Chemistry            *
*                     University of Szeged                             *
*                       Rerrich B. ter 1.                              *
*                         H-6720 Szeged                                *
*                            Hungary                                   *
*                                                                      *
************************************************************************
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      integer element_number
      character title*70,symb*2,string*14,pgrp*3,pgrps*3,
     &      sym_opt*3,filename*128,argtemp*50,fname*100,atom_symbol*2
      logical maxsym,error,fromfile,all,lperm,gaussian,issubgroup,
     &    subset
c
c   mat: maximum number of atoms
c
      parameter(mat=250,nmax=150)
      common/data/ wt(90),symb(90)
      character pgsymb*3,irsymb*4
      common/chartab/ nir(2,55),chtab(14,322),
     &   nsymop(14,4,55),nrotharm(3,322),pgsymb(57),irsymb(322)
      common/subgroups/ nsgb(2,57),nsgr(406)
      dimension coord(3,mat),nat(mat),pc(3), coord_ret(3,mat)
c
      dimension symn(3,nmax),nsym(nmax,5),natsym(mat)
      dimension symn2(3,nmax),nsym2(nmax,5)
c
      dimension rep(nmax),nsymu(nmax,5),rirrep(nmax)
      dimension nper(mat,mat),nscl(mat,mat),nccl(mat)
      dimension pgrps(mat),range(2,20)
      dimension nperm(mat,nmax)

      dimension nsubatoms(mat),coord_all(3,mat),nat_all(mat)
c
c delta: a positive number representing the distortion of the geometry
c
      data delta/0.010000000d0/
c
      data nout/2/,maxsym/.false./,all/.false./
      data del_from,del_to/5.d-2,5.d-3/
      data lperm/.false./
      data ncr,nsr,nsg/0,0,0/
      data filename/''/,fromfile/.false./
      data gaussian/.false./
      data subset/.false./
c
      do i=1,iargc()
         call getarg(i,argtemp)
         if(argtemp.eq.'all') then
c optimize all subgroups
            all=.true.
            maxsym=.true.
         elseif(index(argtemp,'tol=').ne.0) then
c set tolerance
            read(argtemp(scan(argtemp,'=')+1:),*) delta
         elseif(index(argtemp,'tolsup=').ne.0) then
c set interval for determination of the highest order point group
            read(argtemp(scan(argtemp,'=')+1:),*) del_from
         elseif(index(argtemp,'tolinf=').ne.0) then
            read(argtemp(scan(argtemp,'=')+1:),*) del_to
         elseif(index(argtemp,'out=').ne.0) then
c controls output
            read(argtemp(scan(argtemp,'=')+1:),*) nout
         elseif(index(argtemp,'permutations').ne.0) then
c determination of permutations
            lperm=.true.
         elseif(index(argtemp,'maxsym').ne.0) then
c search for highest symmetry between tol_high and tol_low
            maxsym=.true.
         elseif(index(argtemp,'subset').ne.0) then
c search for highest symmetry between tol_high and tol_low
            subset=.true.
         elseif(index(argtemp,'gaussian').ne.0) then
c write Gaussian input
            gaussian=.true.
         else
c any other argument is threated as an input file
            fromfile=.true.
            filename=argtemp
         end if
      end do

c decide source of input
      if(fromfile) then
         open(1,file=filename,status='old',iostat=nerr)
         if(nerr.ne.0) then
            write(*,'(//a)')
     &         'ERROR: Couldn''t open file '//trim(filename)//'.'
            stop
         end if
         in=1
      else
         in=5
      end if

      string='Job starts at '
      call timestamp(string)
c
      write(*,'(a)') repeat('*',71)
      write(*,'(10x,a)')
     &'Program SYMMETRY determines all the symmetry elements'
      write(*,'(15x,a)') 'and symmetry operations of rigid molecules.'
      write(*,'(12x,a)') 
     &'(C) G. Tasi, L. Gyevi-Nagy, R. Tobias, T.S. Tasi'    
      write(*,'(12x,a)') 
     &'Department of Applied and Environmental Chemistry'     
      write(*,'(26x,a)') 'University of Szeged'
      write(*,'(29x,a)') 'H-6720 Szeged'
      write(*,'(32x,a)') 'Hungary'
      write(*,'(a)') repeat('*',71)
      write(*,'(/a,f16.8)') '-- Delta (distortion parameter)=',delta
c
c output file
      open(1001,file="a.tcl")

c Read in the atomic number and Cartesian coordinates of atoms
c
      title = "Molecule"
      read(in,*) natoms
      if(natoms.gt.mat) then
         write(*,'(a)') 'ERROR: Too many atoms.'
         stop
      end if
      read(in,*) 
      do i=1,natoms
         read(in,*) atom_symbol,(coord(j,i),j=1,3)
		 !write(*,*) element_number(atom_symbol,symb)
         nat(i)=element_number(atom_symbol,symb)
      end do
      if(subset) then
         read(in,*) nsub
         read(in,*) (nsubatoms(j),j=1,nsub)
      end if
c
      write(*,'(//a,a/)') '-- Title: ',title
      write(*,'(a)') '-- Molecular geometry (Angstroms) --'
      write(*,'(/1x,a)') 
     &     'Atomic'
      write(*,'(1x,a,10x,a1,16x,a1,15x,a1,15x,a1)') 
     &     'number','x','y','z','w'
      do i=1,natoms
         write(*,'(1x,i3,2x,4f16.6)') 
     &        nat(i),(coord(j,i),j=1,3),wt(nat(i))
      end do
c
c Calculation of the COM (centre of mass) of the molecule
c
      call cmass(natoms,nat,wt,coord,wmol,cmx,cmy,cmz)
      pc(1)=cmx
      pc(2)=cmy
      pc(3)=cmz
      write(*,'(/a,f16.6)') '-- Molecular weight= ',wmol
      write(*,'(/a/)') '-- Centre of mass'
      write(*,'(3(3x,a5,f13.6))') 
     &      'cmx= ',pc(1),' cmy= ',pc(2),' cmz= ',pc(3)
c write to a.tcl
      write(1001,'(a,f9.4,f9.4,f9.4,a)') 
     &" variable com {",pc(1),pc(2),pc(3)," }"
      write(*,'(3(3x,a5,f13.6))') 
      if(subset) then
         write(*,'(/a)') 
     &      '-- Using only the requested subset from this point'
         do i=1,natoms
            do j=1,3
               coord_all(j,i)=coord(j,i)
            end do
            nat_all(i)=nat(i)
         end do
         do i=1,nsub
            do j=1,3
               coord(j,i)=coord_all(j,nsubatoms(i))
            end do
            nat(i)=nat_all(nsubatoms(i))
         end do
         natoms_all=natoms
         natoms=nsub
      end if
c
c Shift the origin of the Cartesian system to COM
c
      call cshift(natoms,coord,pc)
c
      if(nout.eq.2) then
         write(*,'(/a/)')
     &     '-- Cartesian coordinates related to centre of mass --'
         do i=1,natoms
            write(*,'(i3,2x,3f16.6)') nat(i),(coord(j,i),j=1,3)
         end do
      end if
c
c Find symmetry operations
c
      call sym_elements(natoms,nat,coord,symb,delta,ng,ni,nsg,ncr,nsr,
     &  np,symn,nsym,nout,nprm,nper,del_opt,nseq,nccl,nscl)
c
c Determine point group and framework group
c
      call point_group(ng,ni,nsg,ncr,nsr,np,pgrp,nout)
      string='Job ends at '
      call timestamp(string)
c
      end
