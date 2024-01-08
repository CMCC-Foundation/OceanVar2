subroutine ini_bmd

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2007 Srdjan Dobricic, CMCC, Bologna                         !
!                                                                          !
!    This file is part of OceanVar.                                        !
!                                                                          !
!    OceanVar is free software: you can redistribute it and/or modify.     !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!    OceanVar is distributed in the hope that it will be useful,           !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
!    GNU General Public License for more details.                          !
!                                                                          !
!    You should have received a copy of the GNU General Public License     !
!    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.     !
!                                                                          !
!---------------------------------------------------------------------------

!-----------------------------------------------------------------------
!                                                                      !
! Initialise the barotropic model                                      !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
! Version 1.1: P.Oddo 2014                                             !
!-----------------------------------------------------------------------


  use set_knd
  use drv_str
  use grd_str
  use bmd_str
  use mpi_str

  implicit none

  include 'mpif.h'

  INTEGER(i4)    :: i, j
  REAL(r8)       :: bnmp, bnma
  INTEGER        :: ierr
      
       bmd%g     = 9.8066499999999994_r8
       bmd%nstp  = int( 24._r8 * 3600._r8 / bmd%dt )
       bmd%nstps = int( bmd%ndy * bmd%nstp )
       bmd%nstpa = int( bmd%ady * bmd%nstp )
       bmd%alp2  = 1.0_r8 - bmd%alp1


       bmd%df1 = bmd%fc1 * grd%adxdy**2
       bmd%df2 = bmd%fc2 * grd%adxdy**2

     ALLOCATE ( bmd%itr(bmd%nstps) )
     ALLOCATE ( bmd%mst(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%msu(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%msv(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%hgt(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%hgu(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%hgv(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%dxu(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%dxv(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%dyu(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%dyv(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%a1(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%a2(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%a3(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%a4(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%a0(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%a00(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%bx(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%by(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%b_x(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
     ALLOCATE ( bmd%b_y(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
     ALLOCATE ( bmd%dns(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
     ALLOCATE ( bmd%bxby(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%rgh(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%etb(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%ub(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%vb(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%etn(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%un(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%vn(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%eta(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%ua(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%va(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%etm(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%um(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%vm(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%div(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%cu(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%cv(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%dux(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%duy(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%dvx(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%dvy(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%etx(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%ety(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%bfu(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( bmd%bfv(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )

     bmd%mst(:,:) = 0.0_r8
     bmd%msu(:,:) = 0.0_r8
     bmd%msv(:,:) = 0.0_r8
     bmd%hgt(:,:) = 0.0_r8
     bmd%hgu(:,:) = 0.0_r8
     bmd%hgv(:,:) = 0.0_r8
     bmd%dxu(:,:) = 0.0_r8
     bmd%dyv(:,:) = 0.0_r8
     bmd%dyu(:,:) = 0.0_r8
     bmd%dxv(:,:) = 0.0_r8
     bmd%a1 (:,:) = 0.0_r8
     bmd%a2 (:,:) = 0.0_r8
     bmd%a3 (:,:) = 0.0_r8
     bmd%a4 (:,:) = 0.0_r8
     bmd%a0 (:,:) = 0.0_r8
     bmd%a00(:,:) = 0.0_r8

     bmd%itr(:) = 0


     do j=1-grd%jas,grd%jm+grd%jae
      do i=1-grd%ias,grd%im+grd%iae
       if(grd%hgt(i,j).gt.0.0_r8) then
!         bmd%hgt(i,j) = max(15.,grd%hgt(i,j))
          bmd%hgt(i,j) = grd%hgt(i,j)
!sd
!         bmd%hgt(i,j) = min(4000.,max(15.,grd%hgt(i,j)))
!sd
       else
          bmd%hgt(i,j) = 0.0_r8
       endif
      enddo
     enddo

     if(mpi%ir.eq.1      ) bmd%hgt(1,:)      = 0.0_r8
     if(mpi%jr.eq.1      ) bmd%hgt(:,1)      = 0.0_r8
     if(mpi%ir.eq.mpi%irm) bmd%hgt(grd%im,:) = 0.0_r8
     if(mpi%jr.eq.mpi%jrm) bmd%hgt(:,grd%jm) = 0.0_r8

!-------------------------------------------------

     do j=1-grd%jas,grd%jm+grd%jae
      do i=1-grd%ias,grd%im+grd%iae
       if(bmd%hgt(i,j).gt.0.0_r8) then
         bmd%mst(i,j) = 1.0_r8
       else
         bmd%mst(i,j) = 0.0_r8
       endif
      enddo
     enddo

     bnma = 0.0_r8
     do j=2-grd%jas,grd%jm-1+grd%jae
      do i=2-grd%ias,grd%im-1+grd%iae
       if(bmd%mst(i,j).eq.1.) bnma = bnma + 1._r8
      enddo
     enddo

 if(mpi%nproc.gt.1)then

  call mpi_reduce( bnma, bnmp, 1, mpi%r8, mpi_sum, 0, mpi%comm, ierr)
  if(mpi%myrank.eq.0)then
     bnma = bnmp
  endif
  call mpi_bcast( bnma, 1, mpi%r8, 0, mpi%comm, ierr)

 endif

  bmd%bnm = bnma


     do j=1-grd%jas,grd%jm+grd%jae
      do i=2-grd%ias,grd%im+grd%iae
!       bmd%hgu(i,j) = (bmd%hgt(i,j)+bmd%hgt(i-1,j))*0.5
       bmd%hgu(i,j) = min(bmd%hgt(i,j),bmd%hgt(i-1,j))
       bmd%dxu(i,j) = (grd%dx(i,j)+grd%dx(i-1,j))*0.5_r8
       bmd%dyu(i,j) = (grd%dy(i,j)+grd%dy(i-1,j))*0.5_r8
       bmd%msu(i,j) =  bmd%mst(i,j)*bmd%mst(i-1,j)
      enddo
     enddo
       if(mpi%ir.eq.1 ) bmd%dxu(1,:) = bmd%dxu(2,:)
       if(mpi%ir.eq.1 ) bmd%dyu(1,:) = bmd%dyu(2,:)
       if(mpi%ir.eq.1 ) bmd%msu(1,:) = 0.0_r8

 if(mpi%nproc.gt.1) then
  call eao_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%dyu)
 endif

     do j=2-grd%jas,grd%jm+grd%jae
      do i=1-grd%ias,grd%im+grd%iae
!       bmd%hgv(i,j) = (bmd%hgt(i,j)+bmd%hgt(i,j-1))*0.5
       bmd%hgv(i,j) = min(bmd%hgt(i,j),bmd%hgt(i,j-1))
       bmd%dxv(i,j) = (grd%dx(i,j)+grd%dx(i,j-1))*0.5_r8
       bmd%dyv(i,j) = (grd%dy(i,j)+grd%dy(i,j-1))*0.5_r8
       bmd%msv(i,j) =  bmd%mst(i,j)*bmd%mst(i,j-1)
      enddo
     enddo
       if(mpi%jr.eq.1 ) bmd%dxv(:,1) = bmd%dxv(:,2)
       if(mpi%jr.eq.1 ) bmd%dyv(:,1) = bmd%dyv(:,2)
       if(mpi%jr.eq.1 ) bmd%msv(:,1) = 0.0_r8

 if(mpi%nproc.gt.1) then
  call eao_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%dxv)
 endif


     do j= 2-grd%jas,grd%jm-1+grd%jae
      do i= 2-grd%ias,grd%im-1+grd%iae
        bmd%a1(i,j) = bmd%alp1**2 *(bmd%dt**2)*bmd%g*bmd%hgu(i+1,j)/grd%dx(i,j)**2*bmd%msu(i+1,j)
        bmd%a2(i,j) = bmd%alp1**2 *(bmd%dt**2)*bmd%g*bmd%hgu(i  ,j)/grd%dx(i,j)**2*bmd%msu(i  ,j)
        bmd%a3(i,j) = bmd%alp1**2 *(bmd%dt**2)*bmd%g*bmd%hgv(i,j+1)/grd%dy(i,j)**2*bmd%msv(i,j+1)
        bmd%a4(i,j) = bmd%alp1**2 *(bmd%dt**2)*bmd%g*bmd%hgv(i,j  )/grd%dy(i,j)**2*bmd%msv(i,j  )
        bmd%a0(i,j) = (bmd%a1(i,j)+bmd%a2(i,j)+bmd%a3(i,j)+bmd%a4(i,j) +1.0_r8) 
       bmd%a00(i,j) = (bmd%a1(i,j)+bmd%a2(i,j)+bmd%a3(i,j)+bmd%a4(i,j)) *bmd%mst(i,j)
      enddo
     enddo

!ADANI  Bottom Friction not used
     do j=1-grd%jas,grd%jm+grd%jae
      do i=2-grd%ias,grd%im+grd%iae
       bmd%bfu(i,j) = 0.1_r8/max(0.1_r8,bmd%hgu(i,j)) * bmd%msu(i,j)
       bmd%bfu(i,j) = 0.01_r8* bmd%msu(i,j)
      enddo
     enddo
     do j=2-grd%jas,grd%jm+grd%jae
      do i=1-grd%ias,grd%im+grd%iae
       bmd%bfv(i,j) = 0.1_r8/max(0.1_r8,bmd%hgv(i,j)) * bmd%msv(i,j)
!ADANI to make it similar to OceanVar       bmd%bfu(i,j) = 0.01* bmd%msu(i,j)
       bmd%bfv(i,j) = 0.01_r8* bmd%msv(i,j)
      enddo
     enddo

end subroutine ini_bmd
