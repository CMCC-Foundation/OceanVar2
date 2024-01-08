subroutine bar_mod_ad
  
!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2007 Srdjan Dobricic, CMCC, Bologna                         !
!                                                                          !
!    This file is part of OceanVar.                                          !
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
!    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.       !
!                                                                          !
!---------------------------------------------------------------------------

!-----------------------------------------------------------------------
!                                                                      !
! Barotropic model  (adjoint)                                          !
!                                                                      !
! Version 1:   S.Dobricic 2007                                           !
! Version 1.1: M.Adani    2023 
! 
!-----------------------------------------------------------------------
  
  
  use set_knd
  use drv_str
  use grd_str
  use bmd_str
  use mpi_str

  implicit none

  INTEGER(i4)    :: i, j, kstp, ierr
  REAL(r8)       :: tempb,tempb0,tempb1,tempb2,tempb3
  
  
  bmd%bxby(:,:) = 0._r8
  bmd%cu(:,:) = 0._r8
  bmd%cv(:,:) = 0._r8
  bmd%div(:,:) = 0._r8
  bmd%dux(:,:) = 0._r8
  bmd%duy(:,:) = 0._r8
  bmd%dvx(:,:) = 0._r8
  bmd%dvy(:,:) = 0._r8
  bmd%eta(:,:) = 0._r8
  bmd%etb(:,:) = 0._r8
  bmd%etx(:,:) = 0._r8
  bmd%ety(:,:) = 0._r8
  bmd%rgh(:,:) = 0._r8
  bmd%ua(:,:) = 0._r8
  bmd%ub(:,:) = 0._r8
  bmd%un(:,:) = 0._r8
  bmd%va(:,:) = 0._r8
  bmd%vb(:,:) = 0._r8
  bmd%vn(:,:) = 0._r8
  grd%bx(:,:) = 0._r8
  grd%by(:,:) = 0._r8
  !bmd%etm(:,:) = 0._r8
  !bmd%um(:,:) = 0._r8
  !bmd%vm(:,:) = 0._r8
  
  !---
  
  do kstp = bmd%nstps, 1, -1
  
  
  
     if (mod(kstp,bmd%nstpa) .eq. 0) then
  
        bmd%vm(:,:) = 0.0_r8
        bmd%um(:,:) = 0.0_r8
  
        bmd%etm(:,:) = grd%eta_ad(:,:)
        grd%eta_ad(:,:) = 0.0_r8
  
        bmd%vm(:,:) =  bmd%vm(:,:)/bmd%nstpa
        bmd%um(:,:) =  bmd%um(:,:)/bmd%nstpa
        bmd%etm(:,:) = bmd%etm(:,:)/bmd%nstpa*bmd%mst
  
     endif
     bmd%vb(:,:)  = bmd%vb(:,:) + bmd%vm(:,:)
     bmd%ub(:,:)  = bmd%ub(:,:) + bmd%um(:,:)
     bmd%etb(:,:) = bmd%etb(:,:) + bmd%etm(:,:)
     bmd%va(:,:)  = bmd%va(:,:) + bmd%vn(:,:)
     bmd%vn(:, :) = bmd%vb(:, :)
     bmd%ua(:,:)  = bmd%ua(:,:) + bmd%un(:,:)
     bmd%un(:, :) = bmd%ub(:, :)
     bmd%vb(:, :) = 0.05_r8*bmd%vn(:, :)
     bmd%ub(:, :) = 0.05_r8*bmd%un(:, :)
     bmd%eta(:,:) = bmd%eta(:,:) + bmd%etb(:,:)
     bmd%etb(:,:) = 0.0_r8
     bmd%va(:,:)  = bmd%va(:,:) + 0.05_r8*bmd%vn(:,:)
     bmd%vn(:, :) = (1.0_r8-2.0_r8*0.05_r8)*bmd%vn(:, :)
     bmd%ua(:,:)  = bmd%ua(:,:) + 0.05_r8*bmd%un(:,:)
     bmd%un(:, :) = (1.0_r8-2.0_r8*0.05_r8)*bmd%un(:, :)
  
  
   if(mpi%nproc.gt.1)then
  !  call exa_mpi( 0_i4, 0_i4, grd%im+1, 1_i4, grd%im, 0_i4, grd%jm+1, 1_i4, grd%jm,     &
  !                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%eta)
  !  call eao_mpi( 0_i4, 0_i4, grd%im+1, 1_i4, grd%im, 0_i4, grd%jm+1, 1_i4, grd%jm,     &
  !                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%ua )
  !  call eao_mpi( 0_i4, 0_i4, grd%im+1, 1_i4, grd%im, 0_i4, grd%jm+1, 1_i4, grd%jm,     &
  !                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%va )
   endif
     if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm) bmd%dvx(grd%im+1,:)  = 0.0_r8
     if(mpi%nproc.gt.1 .and. mpi%jr.lt.mpi%jrm) bmd%dvy(:,grd%jm+1)  = 0.0_r8
     do j=grd%jm+grd%jae-1,2-grd%jas,-1
        do i=grd%im+grd%iae-1,2-grd%ias,-1
           tempb3 = bmd%df2*bmd%msv(i, j)*bmd%va(i, j)
           tempb2 = tempb3/bmd%dxv(i, j)
           tempb1 = tempb3/bmd%dyv(i, j)
           bmd%dvy(i, j+1) = bmd%dvy(i, j+1) + tempb1
           bmd%dvy(i, j)   = bmd%dvy(i, j) - tempb1
           bmd%dvx(i+1, j) = bmd%dvx(i+1, j) + tempb2
           bmd%dvx(i, j)   = bmd%dvx(i, j) - tempb2
        end do
     end do
     if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm) bmd%dux(grd%im+1,:)  = 0.0_r8
     if(mpi%nproc.gt.1 .and. mpi%jr.lt.mpi%jrm) bmd%duy(:,grd%jm+1)  = 0.0_r8
     do j=grd%jm+grd%jae-1,2-grd%jas,-1
        do i=grd%im+grd%iae-1,2-grd%ias,-1
           tempb3 = bmd%df2*bmd%msu(i, j)*bmd%ua(i, j)
           tempb2 = tempb3/bmd%dxu(i, j)
           tempb1 = tempb3/bmd%dyu(i, j)
           bmd%duy(i, j+1) = bmd%duy(i, j+1) + tempb1
           bmd%duy(i, j) = bmd%duy(i, j) - tempb1
           bmd%dux(i+1, j) = bmd%dux(i+1, j) + tempb2
           bmd%dux(i, j) = bmd%dux(i, j) - tempb2
        end do
     end do
     if(mpi%nproc.gt.1)then
        call exo_mpi( 1_i4, 0_i4,                                       &
                      1_i4, grd%im+1, 1_i4,                             &
                      0_i4, grd%jm, 0_i4,                               &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dux)
        call exo_mpi( 1_i4, 0_i4,                                       &
                      1_i4, grd%im+1, 1_i4,                             &
                      0_i4, grd%jm, 0_i4,                               &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dvx)
        call exo_mpi( 1_i4, 0_i4,                                       &
                      0_i4, grd%im, 0_i4,                               &
                      1_i4, grd%jm+1, 1_i4,                             &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%duy)
        call exo_mpi( 1_i4, 0_i4,                                       &
                      0_i4, grd%im, 0_i4,                               &
                      1_i4, grd%jm+1, 1_i4,                             &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dvy)
     endif
  
     if(mpi%nproc.gt.1 .and. mpi%jr.gt.1) bmd%ua(:,0)  = 0.0_r8
     if(mpi%nproc.gt.1 .and. mpi%jr.gt.1) bmd%va(:,0)  = 0.0_r8
     do j=grd%jm,2-grd%jas,-1
        do i=grd%im,1,-1
           tempb3 = bmd%dvy(i, j)/bmd%dyv(i, j)
           bmd%dvy(i, j) = 0.0_r8
           bmd%va(i, j) = bmd%va(i, j) + tempb3
           bmd%va(i, j-1) = bmd%va(i, j-1) - tempb3
           tempb3 = bmd%duy(i, j)/bmd%dyu(i, j)
           bmd%duy(i, j) = 0.0_r8
           bmd%ua(i, j) = bmd%ua(i, j) + tempb3
           bmd%ua(i, j-1) = bmd%ua(i, j-1) - tempb3
        end do
     end do
     if(mpi%nproc.gt.1 .and. mpi%ir.gt.1) bmd%ua(0,:)  = 0.0_r8
     if(mpi%nproc.gt.1 .and. mpi%ir.gt.1) bmd%va(0,:)  = 0.0_r8
     do j=grd%jm,1,-1
        do i=grd%im,2-grd%ias,-1
           tempb3 = bmd%dvx(i, j)/bmd%dxv(i, j)
           bmd%dvx(i, j) = 0.0_r8
           bmd%va(i, j) = bmd%va(i, j) + tempb3
           bmd%va(i-1, j) = bmd%va(i-1, j) - tempb3
           tempb3 = bmd%dux(i, j)/bmd%dxu(i, j)
           bmd%dux(i, j) = 0.0_r8
           bmd%ua(i, j) = bmd%ua(i, j) + tempb3
           bmd%ua(i-1, j) = bmd%ua(i-1, j) - tempb3
        end do
     end do
  
     if(mpi%nproc.gt.1)then
        call exo_mpi( 1_i4, 0_i4,                                       &
                     -1_i4, 0_i4, grd%im,                               &
                     -1_i4, 0_i4, grd%jm,                               &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%ua)
        call exo_mpi( 1_i4, 0_i4,                                       &
                     -1_i4, 0_i4, grd%im,                               &
                     -1_i4, 0_i4, grd%jm,                               &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%va)
     endif
  
     if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm) bmd%dvx(grd%im+1,:)  = 0.0_r8
     if(mpi%nproc.gt.1 .and. mpi%jr.lt.mpi%jrm) bmd%dvy(:,grd%jm+1)  = 0.0_r8
     if(mpi%nproc.gt.1 .and. mpi%jr.gt.1      ) bmd%eta(:,0)         = 0.0_r8
     do j=grd%jm+grd%jae-1,2-grd%jas,-1
        do i=grd%im+grd%iae-1,1,-1
           bmd%vb(i, j) = bmd%vb(i, j) + bmd%va(i, j)
           tempb3 = -(bmd%dt*bmd%msv(i, j)*bmd%va(i, j))
           tempb1 = bmd%df1*bmd%msv(i, j)*bmd%va(i, j)
           bmd%ety(i, j) = bmd%ety(i, j) - bmd%va(i, j)
           bmd%va(i, j) = 0.0_r8
           tempb0 = tempb1/bmd%dxv(i, j)
           tempb = tempb1/bmd%dyv(i, j)
           bmd%dvy(i, j+1) = bmd%dvy(i, j+1) + tempb
           bmd%dvy(i, j) = bmd%dvy(i, j) - tempb
           bmd%dvx(i+1, j) = bmd%dvx(i+1, j) + tempb0
           bmd%dvx(i, j) = bmd%dvx(i, j) - tempb0
           bmd%cv(i, j) = bmd%cv(i, j) + tempb3
           tempb2 = bmd%alp1*bmd%g*bmd%hgv(i, j)*tempb3/bmd%dyv(i, j)
           grd%by(i, j) = grd%by(i, j) + tempb3
           bmd%eta(i, j) = bmd%eta(i, j) + tempb2
           bmd%eta(i, j-1) = bmd%eta(i, j-1) - tempb2
        end do
     end do
     if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm) bmd%dux(grd%im+1,:)  = 0.0_r8
     if(mpi%nproc.gt.1 .and. mpi%jr.lt.mpi%jrm) bmd%duy(:,grd%jm+1)  = 0.0_r8
     if(mpi%nproc.gt.1 .and. mpi%ir.gt.1      ) bmd%eta(0,:)         = 0.0_r8
     do j=grd%jm+grd%jae-1,1,-1
        do i=grd%im+grd%iae-1,2-grd%ias,-1
           bmd%ub(i, j) = bmd%ub(i, j) + bmd%ua(i, j)
           tempb1 = -(bmd%dt*bmd%msu(i, j)*bmd%ua(i, j))
           tempb = bmd%df1*bmd%msu(i, j)*bmd%ua(i, j)
           bmd%etx(i, j) = bmd%etx(i, j) - bmd%ua(i, j)
           bmd%ua(i, j) = 0.0_r8
           tempb2 = tempb/bmd%dxu(i, j)
           tempb3 = tempb/bmd%dyu(i, j)
           bmd%duy(i, j+1) = bmd%duy(i, j+1) + tempb3
           bmd%duy(i, j) = bmd%duy(i, j) - tempb3
           bmd%dux(i+1, j) = bmd%dux(i+1, j) + tempb2
           bmd%dux(i, j) = bmd%dux(i, j) - tempb2
           bmd%cu(i, j) = bmd%cu(i, j) + tempb1
           tempb0 = bmd%alp1*bmd%g*bmd%hgu(i, j)*tempb1/bmd%dxu(i, j)
           grd%bx(i, j) = grd%bx(i, j) + tempb1
           bmd%eta(i, j) = bmd%eta(i, j) + tempb0
           bmd%eta(i-1, j) = bmd%eta(i-1, j) - tempb0
        end do
     end do
  
     if(mpi%nproc.gt.1)then
  !  call exa_mpi( 1_i4, 0_i4, grd%im+1, 1_i4, grd%im, 0_i4, grd%jm+1, 1_i4, grd%jm,     &
  !                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%eta)
        call exo_mpi( 1_i4, 0_i4,                                       &
                     -1_i4, 0_i4, grd%im,                               &
                     -1_i4, 0_i4, grd%jm,                               &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%eta)
        call exo_mpi( 1_i4, 0_i4,                                       &
                      1_i4, grd%im+1, 1_i4,                             &
                      0_i4, grd%jm, 0_i4,                               &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dux)
        call exo_mpi( 1_i4, 0_i4,                                       &
                      1_i4, grd%im+1, 1_i4,                             &
                      0_i4, grd%jm, 0_i4,                               &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dvx)
        call exo_mpi( 1_i4, 0_i4,                                       &
                      0_i4, grd%im, 0_i4,                               &
                      1_i4, grd%jm+1, 1_i4,                             &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%duy)
        call exo_mpi( 1_i4, 0_i4,                                       &
                      0_i4, grd%im, 0_i4,                               &
                      1_i4, grd%jm+1, 1_i4,                             &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dvy)
     endif
  
     if(mpi%nproc.gt.1 .and. mpi%jr.gt.1) bmd%ub (:,0)  = 0.0_r8
     if(mpi%nproc.gt.1 .and. mpi%jr.gt.1) bmd%vb (:,0)  = 0.0_r8
     do j=grd%jm,2-grd%jas,-1
        do i=grd%im,1,-1
           tempb1 = bmd%dvy(i, j)/bmd%dyv(i, j)
           bmd%dvy(i, j) = 0.0_r8
           bmd%vb(i, j) = bmd%vb(i, j) + tempb1
           bmd%vb(i, j-1) = bmd%vb(i, j-1) - tempb1
           tempb1 = bmd%duy(i, j)/bmd%dyv(i, j)
           bmd%duy(i, j) = 0.0_r8
           bmd%ub(i, j) = bmd%ub(i, j) + tempb1
           bmd%ub(i, j-1) = bmd%ub(i, j-1) - tempb1
        end do
     end do
     if(mpi%nproc.gt.1 .and. mpi%ir.gt.1) bmd%ub (0,:)  = 0.0_r8
     if(mpi%nproc.gt.1 .and. mpi%ir.gt.1) bmd%vb (0,:)  = 0.0_r8
     do j=grd%jm,1,-1
        do i=grd%im,2-grd%ias,-1
           tempb1 = bmd%dvx(i, j)/bmd%dxu(i, j)
           bmd%dvx(i, j) = 0.0_r8
           bmd%vb(i, j) = bmd%vb(i, j) + tempb1
           bmd%vb(i-1, j) = bmd%vb(i-1, j) - tempb1
           tempb1 = bmd%dux(i, j)/bmd%dxu(i, j)
           bmd%dux(i, j) = 0.0_r8
           bmd%ub(i, j) = bmd%ub(i, j) + tempb1
           bmd%ub(i-1, j) = bmd%ub(i-1, j) - tempb1
        end do
     end do
  
  
     if(mpi%nproc.gt.1)then
        call exo_mpi( 1_i4, 0_i4,                                       &
                     -1_i4, 0_i4, grd%im,                               &
                     -1_i4, 0_i4, grd%jm,                               &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%ub)
        call exo_mpi( 1_i4, 0_i4,                                       &
                     -1_i4, 0_i4, grd%im,                               &
                     -1_i4, 0_i4, grd%jm,                               &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%vb)
     endif
  
  
     if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm)  bmd%un (grd%im+1,:)  = 0.0_r8
     if(mpi%nproc.gt.1 .and. mpi%jr.gt.1      )  bmd%un (:,0)         = 0.0_r8
     do j=grd%jm+grd%jae-1,2-grd%jas,-1
        do i=grd%im+grd%iae-1,1,-1
           tempb1 = 0.25_r8*bmd%cv(i, j)/bmd%dyv(i, j)
           bmd%cv(i, j) = 0.0_r8
           tempb0 = grd%f(i, j)*tempb1
           tempb = grd%f(i+1, j)*tempb1
           bmd%un(i+1, j) = bmd%un(i+1, j) + bmd%dyu(i+1, j)*tempb
           bmd%un(i+1, j-1) = bmd%un(i+1, j-1) + bmd%dyu(i+1, j-1)*tempb
           bmd%un(i, j) = bmd%un(i, j) + bmd%dyu(i, j)*tempb0
           bmd%un(i, j-1) = bmd%un(i, j-1) + bmd%dyu(i, j-1)*tempb0
        end do
     end do
     if(mpi%nproc.gt.1 .and. mpi%ir.gt.1      )  bmd%vn (0,:)         = 0.0_r8
     if(mpi%nproc.gt.1 .and. mpi%jr.lt.mpi%jrm)  bmd%vn (:,grd%jm+1)  = 0.0_r8
     do j=grd%jm+grd%jae-1,1,-1
        do i=grd%im+grd%iae-1,2-grd%ias,-1
           tempb1 = -(0.25_r8*bmd%cu(i, j)/bmd%dxu(i, j))
           bmd%cu(i, j) = 0.0_r8
           tempb0 = grd%f(i, j)*tempb1
           tempb = grd%f(i, j+1)*tempb1
           bmd%vn(i, j+1) = bmd%vn(i, j+1) + bmd%dxv(i, j+1)*tempb
           bmd%vn(i-1, j+1) = bmd%vn(i-1, j+1) + bmd%dxv(i-1, j+1)*tempb
           bmd%vn(i, j) = bmd%vn(i, j) + bmd%dxv(i, j)*tempb0
           bmd%vn(i-1, j) = bmd%vn(i-1, j) + bmd%dxv(i-1, j)*tempb0
        end do
     end do
  
  
     if(mpi%nproc.gt.1)then
        call exo_mpi( 1_i4, 0_i4,                                       &
                      1_i4, grd%im+1, 1_i4,                             &
                     -1_i4, 0_i4, grd%jm,                               &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%un)
        call exo_mpi( 1_i4, 0_i4,                                       &
                     -1_i4, 0_i4, grd%im,                               &
                      1_i4, grd%jm+1, 1_i4,                             &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%vn)
     endif
  
     call invrt_ad( kstp )
  
     if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm)  bmd%etx(grd%im+1,:)  = 0.0_r8
     if(mpi%nproc.gt.1 .and. mpi%jr.lt.mpi%jrm)  bmd%ety(:,grd%jm+1)  = 0.0_r8
     do j=grd%jm+grd%jae-1,2-grd%jas,-1
        do i=grd%im+grd%iae-1,2-grd%ias,-1
           bmd%bxby(i, j) = bmd%bxby(i, j) + bmd%rgh(i, j)
           bmd%etb(i, j) = bmd%etb(i, j) - bmd%rgh(i, j)
           bmd%div(i, j) = bmd%div(i, j) + bmd%dt*bmd%rgh(i, j)
           tempb1 = -(bmd%alp1*bmd%dt*bmd%rgh(i, j)/grd%dx(i, j))
           tempb0 = -(bmd%alp1*bmd%dt*bmd%rgh(i, j)/grd%dy(i, j))
           bmd%rgh(i, j) = 0.0_r8
           bmd%ety(i, j+1) = bmd%ety(i, j+1) + tempb0
           bmd%ety(i, j) = bmd%ety(i, j) - tempb0
           bmd%etx(i+1, j) = bmd%etx(i+1, j) + tempb1
           bmd%etx(i, j) = bmd%etx(i, j) - tempb1
        end do
     end do
  
  
     if(mpi%nproc.gt.1)then
        call exo_mpi( 1_i4, 0_i4,                                         &
                      1_i4, grd%im+1, 1_i4,                               &
                      0_i4, 1_i4, grd%jm+1,                               &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%etx)
        call exo_mpi( 1_i4, 0_i4,                                         &
                      0_i4, 1_i4, grd%im+1,                               &
                      1_i4, grd%jm+1, 1_i4,                               &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%ety)
     endif
  
     if(mpi%nproc.gt.1 .and. mpi%jr.gt.1) bmd%etb(:,0)  = 0.0_r8
     do j=grd%jm,2-grd%jas,-1
        do i=grd%im+grd%iae-1,2-grd%ias,-1
           tempb1 = bmd%alp2*bmd%dt*bmd%g*bmd%hgv(i, j)*bmd%msv(i, j)*   &
                    bmd%ety(i, j)/bmd%dyv(i, j)
           bmd%ety(i, j) = 0.0_r8
           bmd%etb(i, j) = bmd%etb(i, j) + tempb1
           bmd%etb(i, j-1) = bmd%etb(i, j-1) - tempb1
        end do
     end do
     if(mpi%nproc.gt.1 .and. mpi%ir.gt.1) bmd%etb(0,:)  = 0.0_r8
     do j=grd%jm+grd%jae-1,2-grd%jas,-1
        do i=grd%im,2-grd%ias,-1
           tempb1 = bmd%alp2*bmd%dt*bmd%g*bmd%hgu(i, j)*bmd%msu(i, j)*    &
                    bmd%etx(i, j)/bmd%dxu(i, j)
           bmd%etx(i, j) = 0.0_r8
           bmd%etb(i, j) = bmd%etb(i, j) + tempb1
           bmd%etb(i-1, j) = bmd%etb(i-1, j) - tempb1
        end do
     end do
     if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm) bmd%ub(grd%im+1,:)  = 0.0_r8
     if(mpi%nproc.gt.1 .and. mpi%jr.lt.mpi%jrm) bmd%vb(:,grd%jm+1)  = 0.0_r8
     do j=grd%jm+grd%jae-1,2-grd%jas,-1
        do i=grd%im+grd%iae-1,2-grd%ias,-1
           tempb1 = bmd%mst(i, j)*bmd%div(i, j)
           bmd%div(i, j) = 0.0_r8
           tempb0 = tempb1/grd%dx(i, j)
           tempb = tempb1/grd%dy(i, j)
           bmd%vb(i, j+1) = bmd%vb(i, j+1) + tempb
           bmd%vb(i, j) = bmd%vb(i, j) - tempb
           bmd%ub(i+1, j) = bmd%ub(i+1, j) + tempb0
           bmd%ub(i, j) = bmd%ub(i, j) - tempb0
        end do
     end do
  
  
     if(mpi%nproc.gt.1)then
        call exo_mpi( 1_i4, 0_i4,                                      &
                     -1_i4, 0_i4, grd%im,                              &
                     -1_i4, 0_i4, grd%jm,                              &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%etb)
        call exo_mpi( 1_i4, 0_i4,                                      &
                      1_i4, grd%im+1, 1_i4,                            &
                      0_i4, grd%jm+1, 1_i4,                            &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%ub )
        call exo_mpi( 1_i4, 0_i4,                                      &
                      0_i4, grd%im+1, 1_i4,                            &
                      1_i4, grd%jm+1, 1_i4,                            &
                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%vb )
     endif
  enddo   ! kstp
  bmd%dvy = 0.0_r8
  bmd%duy = 0.0_r8
  bmd%dvx = 0.0_r8
  bmd%dux = 0.0_r8
  bmd%cv = 0.0_r8
  bmd%cu = 0.0_r8
  bmd%div = 0.0_r8
  bmd%ety = 0.0_r8
  bmd%etx = 0.0_r8
  bmd%va = 0.0_r8
  bmd%ua = 0.0_r8
  bmd%eta = 0.0_r8
  bmd%vn = 0.0_r8
  bmd%un = 0.0_r8
  bmd%vb = 0.0_r8
  bmd%ub = 0.0_r8
  bmd%etb = 0.0_r8
  
  
  if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm) grd%bx(grd%im+1,:)  = 0.0_r8
  if(mpi%nproc.gt.1 .and. mpi%jr.lt.mpi%jrm) grd%by(:,grd%jm+1)  = 0.0_r8
  do j=grd%jm-1+grd%jae,2-grd%jas,-1        ! 2,grd%jm-1
     do i=grd%im-1+grd%iae,2-grd%ias,-1       ! 2,grd%im-1
        bmd%bxby(i, j) = bmd%bxby(i, j) + bmd%alp1**2*bmd%rgh(i, j)
        bmd%rgh(i, j) = 0.0_r8
        tempb = -(bmd%dt**2*bmd%mst(i, j)*bmd%bxby(i, j))
        bmd%bxby(i, j) = 0.0_r8
        tempb0 = tempb/grd%dx(i, j)
        tempb1 = tempb/grd%dy(i, j)
        grd%by(i, j+1) = grd%by(i, j+1) + tempb1
        grd%by(i, j) = grd%by(i, j) - tempb1
        grd%bx(i+1, j) = grd%bx(i+1, j) + tempb0
        grd%bx(i, j) = grd%bx(i, j) - tempb0
     enddo
  enddo
  
  if(mpi%nproc.gt.1)then
  
    call exo_mpi( 1_i4, 0_i4,                                         &
                  1_i4, grd%im+1, 1_i4,                               &
                  0_i4, 1_i4, grd%jm+1,                               &
                  1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, grd%bx )
    call exo_mpi( 1_i4, 0_i4,                                         &
                  0_i4, 1_i4, grd%im+1,                               &
                  1_i4, grd%jm+1, 1_i4,                               &
                  1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, grd%by )
  endif
      
end subroutine bar_mod_ad
