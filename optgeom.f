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

c Subroutine opt_geom optimizes the geometry by placing the atoms exactly 
c   on the symmetry elements.

      subroutine opt_geom(natoms,coord,nat,pgroup,ng,ni,nsg,ncr,nsr,np,
     &    nsym,symn,natsym,coord_ret,nout,mat,nmax)

      implicit none

      integer mat,nmax   ! dimensions
      integer natoms    ! number of atoms
      double precision coord    ! coordinates of the atoms
      integer nsg,ncr,nsr   ! number of symmetry operations (sigmas,C's and S's)
      integer nsym    ! types of the symmetry operations
      double precision symn   ! direction and normal vectors of the symmetry elements
      integer natsym    ! atoms are on these elements
      double precision coord_ret   ! symmetized geometry
      character pgroup*3
      integer nat
      integer nout    ! output control

      double precision PI
      integer nsg1,nsg2    ! first and last symmetry planes
      integer ncr1,ncr2    ! first and last proper axes
      integer nsr1,nsr2    ! first and last improper axes
      integer i,j,mm   ! cycle variables
      integer elem
      double precision atom_transformed    ! coordinates of an atom after transformation
      logical optimized
      double precision datan,dot
      integer close,nopt

      integer np,ni,ng,nprm,nseq,nper,nccl,nscl
      double precision del_opt
      character pgrp*3,symb*2
      dimension symb(90),nper(mat,mat),nccl(mat,mat),nscl(mat,mat)

      dimension coord(3,mat),nsym(nmax,5),symn(3,nmax),
     &      natsym(mat),coord_ret(3,mat),nat(mat),
     &      atom_transformed(3),optimized(mat)

      integer npg
      logical check, issubgroup
      integer nsgb, nsgr
      common/subgroups/ nsgb(2,57),nsgr(406)
      character pgsymb*3,irsymb*4
      integer nir, nsymop, nrotharm
      double precision chtab
      common/chartab/ nir(2,55),chtab(14,322),
     &   nsymop(14,4,55),nrotharm(3,322),pgsymb(57),irsymb(322)
      double precision del_opt1, del_opt2, range, del_from, del_to
      integer ngroups
      character pgrps*3
      dimension pgrps(20), range(2,20)

      del_from=1.d-3
      del_to=1.d-4

c initialization
      do mm=1,natoms
         optimized(mm)=.false.
      end do
      PI=4.d0*datan(1.d0)
      nsg1=2
      nsg2=nsg+1
      ncr1=nsg2+1
      ncr2=nsg2+ncr
      nsr1=ncr2+1
      nsr2=ncr2+nsr
      nopt=0

      if(nout.gt.0) write(*,'(/a)') '-- OPTIMIZING GEOMETRY --'

c optimize every symmetry-equivalence class
      do i=1,natoms
c ith element of each class on symmetry element elem
        if(optimized(i)) cycle
        elem=natsym(i)

c ith element is in the center of mass
        if(elem.eq.1) then
          do mm=1,3
            coord_ret(mm,i)=0.d0
          end do
          optimized(i)=.true.
          cycle

c ith element is on a symmetry axis
        elseif(ncr1.le.elem.and.elem.le.ncr2) then
c project point onto the axis
          do mm=1,3
            coord_ret(mm,i)=dot(coord(1,i),symn(1,elem),3)*
     &         symn(mm,elem)
          end do
          optimized(i)=.true.

c ith element is on a symmetry plane
        elseif(nsg1.le.elem.and.elem.le.nsg2) then
c project point onto the plane
          do mm=1,3
            coord_ret(mm,i)=coord(mm,i)-
     &         dot(coord(1,i),symn(1,elem),3)*symn(mm,elem)
          end do
          optimized(i)=.true.

c ith element is not on any symmetry element
        elseif(elem.eq.0) then
          do mm=1,3
            coord_ret(mm,i)=coord(mm,i)
          end do
          optimized(i)=.true.

        end if

c perform all symmetry operations on the ith element
        do j=1,nsr2+1

          if(nsym(j,1).eq.0) then
c inversion
            if(nsym(j,2).ne.1) cycle
            call optgeom_invert(coord_ret(1,i),atom_transformed)

c reflection
          elseif(nsym(j,1).eq.1) then
             if(pgroup.eq.'Civ'.or.pgroup.eq.'Dih') exit
             call optgeom_reflect(coord_ret(1,i),
     &            atom_transformed,symn(1,j))

c proper rotation
          elseif(nsym(j,1).eq.2) then
             call optgeom_rotate(coord_ret(1,i),
     &          atom_transformed,symn(1,j),2*PI/nsym(j,2)*nsym(j,3))

c improper rotation
          elseif(nsym(j,1).eq.3) then
             call optgeom_srotate(coord_ret(1,i),atom_transformed,
     &            symn(1,j),2*PI/nsym(j,2)*nsym(j,3))

          end if

          call find_close_atom(atom_transformed,natoms,coord,mat,
     &         close)
          if(close.ne.0.and..not.optimized(close)) then
            do mm=1,3
              coord_ret(mm,close)=atom_transformed(mm)
            end do
            optimized(close)=.true.
          else if(close.eq.0) then
             write(*,'(a)') 'Couldn''t find close atom.'
          end if

        end do

      end do

      nopt=0
      do i=1,natoms
        if(optimized(i)) nopt=nopt+1
      end do

      call var_delta(natoms,nat,coord_ret,symb,
     &      del_from,del_to,del_opt1,del_opt2, ngroups, pgrps, range,0)
      call sym_elements(natoms,nat,coord_ret,symb,del_opt1,ng,ni,nsg,
     &    ncr,nsr,np,symn,nsym,0,nprm,nper,del_opt,nseq,nccl,nscl)
      call point_group(ng,ni,nsg,ncr,nsr,np,pgrp,0)
      do npg=1,57
         if(pgrp.eq.adjustl(pgsymb(npg))) exit
      end do
      check=issubgroup(pgroup, pgrp)

      if(nout.eq.2) write(*,'(/5x,a,i4)') 'Number of optimized atoms: ',
     &   nopt

      if(check) then
        if(nout.gt.0) then
          write(*,'(/5x,a)') 'Geometry optimization succesful.'
        end if
        write(*,'(/5x,3a/)') 'Optimized ',trim(pgroup),' geometry:'
        do i=1,natoms
           write(*,'(i3,2x,3f16.8)') nat(i),(coord_ret(j,i),j=1,3)
        end do
      else
        write(*,'(/5x,a)') 'Unsuccesful geometry optimization of ' 
     &         // trim(pgroup) // ' symmetry.'
        if(nout.gt.0) then
          if(nopt.ne.natoms) then
            write(*,'(5x,a,i4)') 'Number of unoptimized atoms: ',
     &         natoms-nopt
          else
            write(*,'(5x,a)') 'Wrong point group after symmetrization.'
          end if
        end if
      end if

      contains
      subroutine optgeom_invert(atom,transformed)
      implicit none
      double precision atom
      double precision transformed
      integer mm
      dimension atom(3),transformed(3)
      do mm=1,3
         transformed(mm)=-atom(mm)
      end do
      end subroutine optgeom_invert

      subroutine optgeom_reflect(atom,transformed,plane_normal)
      implicit none
      double precision atom
      double precision transformed
      double precision plane_normal
      integer i
      double precision dot
      dimension atom(3),transformed(3),plane_normal(3)
      do i=1,3
         transformed(i)=atom(i)-2*dot(plane_normal,atom,3)*
     &         plane_normal(i)
      end do
      end subroutine optgeom_reflect

      subroutine optgeom_rotate(atom,transformed,axis_direction,angle)
      implicit none
      double precision atom
      double precision transformed
      double precision axis_direction
      double precision angle
      integer i
      double precision vec1,vec2
      intrinsic dsin,dcos
      dimension atom(3),transformed(3),axis_direction(3),vec1(3),vec2(3)
      call crossp(axis_direction,atom,vec1)
      call crossp(axis_direction,vec1,vec2)
      do i=1,3
         transformed(i)=atom(i)+dsin(angle)*vec1(i)+(1-dcos(angle))*
     &         vec2(i)
      end do
      end subroutine optgeom_rotate

      subroutine optgeom_srotate(atom,transformed,axis_direction,angle)
      implicit none
      double precision atom
      double precision transformed
      double precision axis_direction
      double precision angle
      double precision rotated
      dimension atom(3),transformed(3),axis_direction(3),rotated(3)
      call optgeom_rotate(atom,rotated,axis_direction,angle)
      call optgeom_reflect(rotated,transformed,axis_direction)
      end subroutine optgeom_srotate

      subroutine find_close_atom(transformed,natoms,coord,mat,
     &      close)
      implicit none
      integer mat
      double precision transformed
      integer natoms
      double precision coord
      integer close
      integer i,j
      double precision distance,dist,diff
      dimension transformed(3),coord(3,mat),diff(3)
      close=0
      distance=1.d10
      do i=1,natoms
        do j=1,3
          diff(j)=coord(j,i)-transformed(j)
        end do
        dist=dot(diff,diff,3)
        if(dist.lt.distance) then
          close=i
          distance=dist
        end if
      end do
      end subroutine find_close_atom

      end subroutine opt_geom
