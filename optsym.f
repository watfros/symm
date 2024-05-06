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

c Subroutine opt_symmetry optimizes symmetry operations stored in symn.

      subroutine opt_symmetry(nsg,ncr,nsr,symn,nsym,pgroup,nout,error)

      implicit none

c arguments
      double precision  symn    ! directions of the symmetry elements
      integer  nsg,ncr,nsr,nsym,nout    ! number and type of symmetry elements,output control
      character pgroup*3    ! point group
      logical error

c local variables
      integer nsg1,nsg2    ! first and last symmetry planes
      integer ncr1,ncr2    ! first and last proper axes
      integer nsr1,nsr2    ! first and last improper axes
      integer i,j,k    ! cycle variables
      integer c3s    ! C3 axes (cubic point groups)
      integer praxis,ord_prax,sig_h
      integer nopt   ! number of optimized symmetry operations
      integer ng   ! total number of symmetry operations to optimize
      double precision datan,PI,dot
      double precision TOL
      double precision TOL_I
      double precision praxis_dv   ! direction vector of the principal axis
      double precision alpha
      double precision s,ri,dsqrt
      double precision vec,vec_rot   ! a vector and that rotated around the principal axis
      double precision v1,v2
      double precision vec_comp  ! vector to compare with
      logical optimized
      integer mm
      double precision temp1,temp2

      parameter(TOL=1.d-1,TOL_I=4.d-2)

      dimension symn(3,150),nsym(150,5),c3s(4),praxis_dv(3),vec(3),
     &      vec_rot(3),v1(3),v2(3),vec_comp(3),optimized(2:151)

c initialization
      praxis=0
      sig_h=0
      nopt=0
      PI=4*datan(1.d0)
      nsg1=2
      nsg2=nsg+1
      ncr1=nsg2+1
      ncr2=nsg2+ncr
      nsr1=ncr2+1
      nsr2=ncr2+nsr
      ng=nsg+ncr+nsr
      nopt=0

      do mm=2,ng+1
         optimized(mm)=.false.
      end do

      if(nout.ge.1) write(*,'(/a)') '-- OPTIMIZING SYMMETRY ELEMENTS --'

      if(trim(pgroup).eq.'') then

c missing point group
          write(*,'(/5x,a/)') 'ERROR: No point group.'
          error=.true.
          return

c nothing to optimize
       elseif(trim(pgroup).eq.'C1'.or.trim(pgroup).eq.'Ci'.or.
     &      trim(pgroup).eq.'Cs'.or.trim(pgroup).eq.'Civ'.or.
     &      trim(pgroup).eq.'Dih'.or.trim(pgroup).eq.'C2'.or.
     &      trim(pgroup).eq.'C3'.or.trim(pgroup).eq.'C4'.or.
     &      trim(pgroup).eq.'C5'.or.trim(pgroup).eq.'C6'.or.
     &      trim(pgroup).eq.'C7'.or.trim(pgroup).eq.'C8'.or.
     &      trim(pgroup).eq.'S4'.or.trim(pgroup).eq.'S6'.or.
     &      trim(pgroup).eq.'S8') then
          if(nout.ge.1) write(*,'(/5x,a)') 'Nothing to optimize.'
          error=.false.
          return

c tetrahedral and octahedral symmetry
       elseif(trim(pgroup).eq.'T'.or.trim(pgroup).eq.'Th'.or.
     &      trim(pgroup).eq.'Td'.or.trim(pgroup).eq.'O'.or.
     &      trim(pgroup).eq.'Oh') then
c C2 axis are perpendicular
          do i=ncr1,ncr2
            if(nsym(i,2).ne.2.or.nsym(i,4).eq.2) cycle
            if(praxis.eq.0) then
              praxis=i
              do mm=1,3
                praxis_dv(mm)=symn(mm,i)
              end do
              optimized(i)=.true.
              cycle
            else
c find a perpendicular C2 axis
              do mm=1,3
                vec(mm)=symn(mm,i)-dot(symn(1,i),praxis_dv,3)*
     &               praxis_dv(mm)
              end do
              temp1=dsqrt(dot(vec,vec,3))
              do mm=1,3
                vec(mm)=vec(mm)/temp1
                symn(mm,i)=vec(mm)
              end do
              optimized(i)=.true.
c optimize 3rd C2 axis
              call crossp(vec,praxis_dv,vec_rot)
              do j=ncr2,ncr1,-1
                if(nsym(j,2).eq.2.and.nsym(j,4).ne.2) then
                  do mm=1,3
                    symn(mm,j)=vec_rot(mm)
                  end do
                  optimized(j)=.true.
                  exit
                end if
              end do
              exit
            end if
          end do

c optimize parallel vectors (sigma_h,S4)
          do i=nsg1,nsr2
            if(optimized(i)) cycle
            if(dabs(dabs(dot(symn(1,i),praxis_dv,3))-1).le.TOL) then
              do mm=1,3
                symn(mm,i)=praxis_dv(mm)
              end do
              optimized(i)=.true.
            else if(dabs(dabs(dot(symn(1,i),vec,3))-1).le.TOL) then
              do mm=1,3
                symn(mm,i)=vec(mm)
              end do
              optimized(i)=.true.
            else if(dabs(dabs(dot(symn(1,i),vec_rot,3))-1).le.TOL) then
              do mm=1,3
                symn(mm,i)=vec_rot(mm)
              end do
              optimized(i)=.true.
            end if
          end do

c construct C3 axes out of C2 axes
          do mm=1,3
            v1(mm)=praxis_dv(mm)+vec(mm)+vec_rot(mm)
            v2(mm)=praxis_dv(mm)+vec(mm)-vec_rot(mm)
          end do
          temp1=dsqrt(dot(v1,v1,3))
          temp2=dsqrt(dot(v2,v2,3))
          do mm=1,3
            v1(mm)=v1(mm)/temp1
            v2(mm)=v2(mm)/temp2
          end do
          do i=nsg1,nsr2
            if(optimized(i)) cycle
            if(dabs(dabs(dot(symn(1,i),v1,3))-1).le.TOL) then
              do mm=1,3
                symn(mm,i)=v1(mm)
              end do
              optimized(i)=.true.
            else if(dabs(dabs(dot(symn(1,i),v2,3))-1).le.TOL) then
              do mm=1,3
                symn(mm,i)=v2(mm)
              end do
              optimized(i)=.true.
            end if
          end do
          do mm=1,3
            v1(mm)=praxis_dv(mm)-vec(mm)+vec_rot(mm)
            v2(mm)=praxis_dv(mm)-vec(mm)-vec_rot(mm)
          end do
          temp1=dsqrt(dot(v1,v1,3))
          temp2=dsqrt(dot(v2,v2,3))
          do mm=1,3
            v1(mm)=v1(mm)/temp1
            v2(mm)=v2(mm)/temp2
          end do
          do i=nsg1,nsr2
            if(optimized(i)) cycle
            if(dabs(dabs(dot(symn(1,i),v1,3))-1).le.TOL) then
              do mm=1,3
                symn(mm,i)=v1(mm)
              end do
              optimized(i)=.true.
            else if(dabs(dabs(dot(symn(1,i),v2,3))-1).le.TOL) then
              do mm=1,3
                symn(mm,i)=v2(mm)
              end do
              optimized(i)=.true.
            end if
          end do

c optimize sigma_d planes and C2' axes
c find C3 axes
          j=0
          do i=ncr1,ncr2
            if(nsym(i,2).eq.3.and.nsym(i,3).eq.1) then
              j=j+1
              c3s(j)=i
            end if
          end do
c construct vectors perpendicular to two C3 axes
          do i=1,3
            do j=i+1,4
              call crossp(symn(1,c3s(i)),symn(1,c3s(j)),vec_comp)
              temp1=dsqrt(dot(vec_comp,vec_comp,3))
              do mm=1,3
                vec_comp(mm)=vec_comp(mm)/temp1
              end do
c find non-optimized pairs
              do k=nsg1,nsr2
                if(optimized(k)) cycle
                if(dabs(dabs(dot(vec_comp,symn(1,k),3))-1).le.TOL) then
                  do mm=1,3
                     symn(mm,k)=vec_comp(mm)
                  end do
                  optimized(k)=.true.
                end if
              end do
            end do
          end do

c icosahedral symmetry
      elseif(trim(pgroup).eq.'I'.or.trim(pgroup).eq.'Ih') then
          s=dtan(54.d0/180.d0*PI)/2*dsqrt(2*(1-dcos(PI-datan(2.d0))))
          ri=1.d0/2.d0*dsqrt(5.d0/2.d0+11.d0/10.d0*dsqrt(5.d0))
          alpha=dacos(1-s**2/2/ri**2)

c optimize C5 axes
          outer: do i=ncr1,ncr2
            if(nsym(i,2).eq.5.and.nsym(i,3).eq.1) then
              do mm=1,3
                vec(mm)=symn(mm,i)
              end do
              optimized(i)=.true.
              do j=nsg1,nsr2
                if(optimized(j)) cycle
                if(dabs(dabs(dot(vec,symn(1,j),3))-1).le.TOL_I) then
                  do mm=1,3
                    symn(mm,j)=vec(mm)
                  end do
                  optimized(j)=.true.
                end if
              end do
              do j=i+1,ncr2
                if(nsym(j,2).eq.5.and.nsym(j,3).eq.1) then
                  call crossp(vec,symn(1,j),praxis_dv)
                  temp1=dsqrt(dot(praxis_dv,praxis_dv,3))
                  do mm=1,3
                    praxis_dv(mm)=praxis_dv(mm)/temp1
                  end do
                  exit outer
                end if
              end do
            end if
          end do outer

c rotation around the perpendicular axis
          call optsym_rotate(vec,praxis_dv,
     &        sign(alpha,dot(vec,symn(1,j),3)),vec_rot)

          alpha=2*PI/5
          do i=1,5
c rotation around a C5 axis
            call optsym_rotate(vec_rot,vec,i*alpha,vec_comp)

c find non-optimized pairs
            do j=nsg1,nsr2
              if(optimized(j)) cycle
              if(dabs(dabs(dot(vec_comp,symn(1,j),3))-1).le.TOL_I) then
                do mm=1,3
                  symn(mm,j)=vec_comp(mm)
                end do
                optimized(j)=.true.
              end if
            end do
          end do

c optimize C3 axes
          call optsym_rotate(vec_rot,vec,alpha,praxis_dv)
          do mm=1,3
            vec_comp(mm)=vec(mm)+vec_rot(mm)+praxis_dv(mm)
          end do
          temp1=dsqrt(dot(vec_comp,vec_comp,3))
          do mm=1,3
            vec_comp(mm)=vec_comp(mm)/temp1
            praxis_dv(mm)=vec_comp(mm)
          end do
          do i=1,5
            call optsym_rotate(praxis_dv,vec,i*alpha,vec_comp)
            do j=nsg1,nsr2
              if(optimized(j)) cycle
              if(dabs(dabs(dot(vec_comp,symn(1,j),3))-1).le.TOL_I) then
                do mm=1,3
                   symn(mm,j)=vec_comp(mm)
                end do
                optimized(j)=.true.
              end if
            end do
          end do
          call optsym_rotate(vec_rot,vec,2*alpha,praxis_dv)
          do mm=1,3
            vec_comp(mm)=vec_rot(mm)-praxis_dv(mm)
          end do
          call optsym_rotate(vec_rot,vec,3*alpha,praxis_dv)
          do mm=1,3
            vec_comp(mm)=vec_comp(mm)-praxis_dv(mm)
          end do
          temp1=dsqrt(dot(vec_comp,vec_comp,3))
          do mm=1,3
            vec_comp(mm)=vec_comp(mm)/temp1
            praxis_dv(mm)=vec_comp(mm)
          end do
          do i=1,5
            call optsym_rotate(praxis_dv,vec,i*alpha,vec_comp)
            do j=nsg1,nsr2
              if(optimized(j)) cycle
              if(dabs(dabs(dot(vec_comp,symn(1,j),3))-1).le.TOL_I) then
                do mm=1,3
                   symn(mm,j)=vec_comp(mm)
                end do
                optimized(j)=.true.
              end if
            end do
          end do

c optimize C2 axes
          do mm=1,3
             vec_comp(mm)=vec(mm)+vec_rot(mm)
          end do
          temp1=dsqrt(dot(vec_comp,vec_comp,3))
          do mm=1,3
             vec_comp(mm)=vec_comp(mm)/temp1
             praxis_dv(mm)=vec_comp(mm)
          end do
          do i=1,5
            call optsym_rotate(praxis_dv,vec,i*alpha,vec_comp)
            do j=nsg1,nsr2
              if(optimized(j)) cycle
              if(dabs(dabs(dot(vec_comp,symn(1,j),3))-1).le.TOL_I) then
                do mm=1,3
                  symn(mm,j)=vec_comp(mm)
                end do
                optimized(j)=.true.
              end if
            end do
          end do
          call optsym_rotate(vec_rot,vec,alpha,praxis_dv)
          do mm=1,3
            vec_comp(mm)=vec_rot(mm)+praxis_dv(mm)
          end do
          temp1=dsqrt(dot(vec_comp,vec_comp,3))
          do mm=1,3
            vec_comp(mm)=vec_comp(mm)/temp1
            praxis_dv(mm)=vec_comp(mm)
          end do
          do i=1,5
            call optsym_rotate(praxis_dv,vec,i*alpha,vec_comp)
            do j=nsg1,nsr2
              if(optimized(j)) cycle
              if(dabs(dabs(dot(vec_comp,symn(1,j),3))-1).le.TOL_I) then
                do mm=1,3
                  symn(mm,j)=vec_comp(mm)
                end do
                optimized(j)=.true.
              end if
            end do
          end do
          call optsym_rotate(vec_rot,vec,2*alpha,praxis_dv)
          do mm=1,3
            vec_comp(mm)=vec_rot(mm)-praxis_dv(mm)
          end do
          temp1=dsqrt(dot(vec_comp,vec_comp,3))
          do mm=1,3
            vec_comp(mm)=vec_comp(mm)/temp1
            praxis_dv(mm)=vec_comp(mm)
          end do
          do i=1,5
            call optsym_rotate(praxis_dv,vec,i*alpha,vec_comp)
            do j=nsg1,nsr2
              if(optimized(j)) cycle
              if(dabs(dabs(dot(vec_comp,symn(1,j),3))-1).le.TOL_I) then
                do mm=1,3
                  symn(mm,j)=vec_comp(mm)
                end do
                optimized(j)=.true.
              end if
            end do
          end do

c all other point groups
      else
c find principal axis
          do i=ncr1,ncr2
            if(nsym(i,4).eq.1) then
              praxis=i
              ord_prax=nsym(i,2)
              alpha=PI/ord_prax/2
              optimized(i)=.true.
            end if
          end do
          do mm=1,3
            praxis_dv(mm)=symn(mm,praxis)
          end do

c optimize parallel vectors
          do i=nsg1,nsr2
            if(optimized(i)) cycle
            if(dabs(dabs(dot(praxis_dv,symn(1,i),3))-1).le.TOL) then
              do mm=1,3
                symn(mm,i)=praxis_dv(mm)
              end do
              optimized(i)=.true.
            end if
          end do

c optimize perpendicular vectors
c find vector to rotate
          do i=nsg1,nsr2
            if(nsym(i,4).lt.2) cycle
            do mm=1,3
              vec(mm)=symn(mm,i)-dot(symn(1,i),praxis_dv,3)*
     &            praxis_dv(mm)
            end do
            temp1=dsqrt(dot(vec,vec,3))
            do mm=1,3
              vec(mm)=vec(mm)/temp1
            end do
            exit
          end do

          do i=1,2*ord_prax
c rotation around the principal axis
            call optsym_rotate(vec,praxis_dv,i*alpha,vec_rot)

c find non-optimized pairs
            do j=nsg1,nsr2
              if(optimized(j).or.nsym(j,4).lt.2) cycle
              do mm=1,3
                vec_comp(mm)=symn(mm,j)
              end do
              if(dabs(dabs(dot(vec_rot,vec_comp,3))-1).le.TOL/ord_prax)
     &            then
                do mm=1,3
                  symn(mm,j)=vec_rot(mm)
                end do
                optimized(j)=.true.
              end if
            end do
          end do

       end if

      do i=2,ng+1
        if(optimized(i)) nopt=nopt+1
      end do

      if(nout.eq.2) write(*,'(/5x,a,i3)') 
     &      'Number of optimized symmetry operations: ',nopt

      if(nopt.eq.ng) then
        if(nout.ge.1) write(*,'(/5x,3a)') 'Succesful optimization of ',
     &      trim(pgroup),' symmetry.'
        error=.false.
      else
        if(nout.ge.1) then
          write(*,'(/5x,a)') 'Unsuccesful optimization of ' // 
     &         trim(pgroup) // ' symmetry.'
          write(*,'(5x,a,i3)') 'Number of unoptimized operations: ',
     &         ng-nopt
        end if
        error=.true.
      end if

      contains
      subroutine optsym_rotate(vector,direction,angle,rotated)
      implicit none
      double precision vector,rotated,direction
      double precision angle
      double precision vec1,vec2
      integer i
      dimension vector(3),rotated(3),direction(3),vec1(3),vec2(3)
      call crossp(direction,vector,vec1)
      call crossp(direction,vec1,vec2)
      do i=1,3
        rotated(i)=vector(i)+dsin(angle)*vec1(i)+(1-dcos(angle))*vec2(i)
      end do
      end subroutine optsym_rotate

      end subroutine opt_symmetry
