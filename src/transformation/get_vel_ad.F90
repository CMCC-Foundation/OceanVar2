subroutine get_vel_ad

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
! Calculate horizontal velocity from geostrophic formula (adjoint)     !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
! Bug correction 21.04.2009  thanks to Andrea Storto                   !
!-----------------------------------------------------------------------


  use set_knd
  use grd_str
  use bmd_str
  use mpi_str

  implicit none

  REAL(r8), ALLOCATABLE, DIMENSION (:,:,:)  :: ud, vd
  integer(i4)    :: k, i, j
  integer ierr

  ALLOCATE ( ud(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
  ALLOCATE ( vd(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))


      ud(:,:,:) = 0.0
      vd(:,:,:) = 0.0

     do k=grd%km,1,-1

        if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm) vd(grd%im+1,:,k) = 0.0
        if(mpi%nproc.gt.1 .and. mpi%jr.gt.1      ) vd(:,0,k)        = 0.0
      do j=2-grd%jas,grd%jm      ! 2:grd%jm
       do i=1,grd%im-1+grd%iae   ! 1:grd%im-1
        vd(i  ,j  ,k) = vd(i  ,j  ,k) + grd%vvl_ad(i,j,k)*0.25 * grd%msk(i,j,k)*grd%msk(i,j-1,k)
        vd(i  ,j-1,k) = vd(i  ,j-1,k) + grd%vvl_ad(i,j,k)*0.25 * grd%msk(i,j,k)*grd%msk(i,j-1,k)
        vd(i+1,j  ,k) = vd(i+1,j  ,k) + grd%vvl_ad(i,j,k)*0.25 * grd%msk(i,j,k)*grd%msk(i,j-1,k)
        vd(i+1,j-1,k) = vd(i+1,j-1,k) + grd%vvl_ad(i,j,k)*0.25 * grd%msk(i,j,k)*grd%msk(i,j-1,k)
        grd%vvl_ad(i,j,k) = 0.0
       enddo
      enddo
        if(mpi%nproc.gt.1 .and. mpi%ir.gt.1      ) ud(0,:       ,k) = 0.0
        if(mpi%nproc.gt.1 .and. mpi%jr.lt.mpi%jrm) ud(:,grd%jm+1,k) = 0.0
      do j=1,grd%jm-1+grd%jae
       do i=2-grd%ias,grd%im
        ud(i  ,j  ,k) = ud(i  ,j  ,k) + grd%uvl_ad(i,j,k)*0.25 * grd%msk(i,j,k)*grd%msk(i-1,j,k)
        ud(i-1,j  ,k) = ud(i-1,j  ,k) + grd%uvl_ad(i,j,k)*0.25 * grd%msk(i,j,k)*grd%msk(i-1,j,k)
        ud(i  ,j+1,k) = ud(i  ,j+1,k) + grd%uvl_ad(i,j,k)*0.25 * grd%msk(i,j,k)*grd%msk(i-1,j,k)
        ud(i-1,j+1,k) = ud(i-1,j+1,k) + grd%uvl_ad(i,j,k)*0.25 * grd%msk(i,j,k)*grd%msk(i-1,j,k)
        grd%uvl_ad(i,j,k) = 0.0
       enddo
      enddo

     enddo

 if(mpi%nproc.gt.1)then
  call exo_mpi( 1_i4, 0_i4,                                       &
                -1_i4, 0_i4,   grd%im,                               &
                 1_i4, grd%jm+1, 1_i4,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, grd%km, ud )
  call exo_mpi( 1_i4, 0_i4,                                       &
                 1_i4, grd%im+1, 1_i4,                               &
                -1_i4, 0_i4,   grd%jm,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, grd%km, vd )
 endif


      if(mpi%nproc.gt.1 .and. mpi%jr.gt.1) grd%eta_ad(:,0) = 0.0
      if(mpi%nproc.gt.1 .and. mpi%ir.gt.1) grd%eta_ad(0,:) = 0.0
    do k=1,grd%km
      do j=2-grd%jas,grd%jm      ! 2:grd%jm
       do i=1,grd%im               ! 1:grd%im
        grd%b_y(i,j,k)    = grd%b_y(i,j,k)    - ud(i,j,k)      / grd%dy(i,j) * grd%msk(i,j,k)*grd%msk(i,j-1,k) / grd%f(i,j)
        grd%eta_ad(i,j  ) = grd%eta_ad(i,j  ) - ud(i,j,k)*9.81 / grd%dy(i,j) * grd%msk(i,j,k)*grd%msk(i,j-1,k) / grd%f(i,j)
        grd%eta_ad(i,j-1) = grd%eta_ad(i,j-1) + ud(i,j,k)*9.81 / grd%dy(i,j) * grd%msk(i,j,k)*grd%msk(i,j-1,k) / grd%f(i,j)
       enddo
      enddo
      do j=1,grd%jm
       do i=2-grd%ias,grd%im
        grd%b_x(i,j,k)    = grd%b_x(i,j,k)    + vd(i,j,k)      / grd%dx(i,j) * grd%msk(i,j,k)*grd%msk(i-1,j,k) / grd%f(i,j)
        grd%eta_ad(i-1,j) = grd%eta_ad(i-1,j) - vd(i,j,k)*9.81 / grd%dx(i,j) * grd%msk(i,j,k)*grd%msk(i-1,j,k) / grd%f(i,j)
        grd%eta_ad(i  ,j) = grd%eta_ad(i  ,j) + vd(i,j,k)*9.81 / grd%dx(i,j) * grd%msk(i,j,k)*grd%msk(i-1,j,k) / grd%f(i,j)
       enddo
      enddo
      ud(:,:,k) = 0.0
      vd(:,:,k) = 0.0
     enddo

   call exo_mpi( 0_i4, 0_i4,                                         &
                -1_i4, 0_i4,   grd%im,                               &
                -1_i4, 0_i4,   grd%jm,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, grd%eta_ad )

     grd%uvl_ad(:,:,:) = 0.0
     grd%vvl_ad(:,:,:) = 0.0

      grd%b_y(:,2-grd%jas:grd%jm,1) = grd%b_y(:,2-grd%jas:grd%jm,1) * 0.5
      grd%b_x(2-grd%ias:grd%im,:,1) = grd%b_x(2-grd%ias:grd%im,:,1) * 0.5

     do k=2,grd%km
      grd%b_y(:,2-grd%jas:grd%jm,k-1) = grd%b_y(:,2-grd%jas:grd%jm,k-1) + grd%b_y(:,2-grd%jas:grd%jm,k) * 0.5
      grd%b_y(:,2-grd%jas:grd%jm,k)   = grd%b_y(:,2-grd%jas:grd%jm,k) * 0.5
      grd%b_x(2-grd%ias:grd%im,:,k-1) = grd%b_x(2-grd%ias:grd%im,:,k-1) + grd%b_x(2-grd%ias:grd%im,:,k) * 0.5
      grd%b_x(2-grd%ias:grd%im,:,k)   = grd%b_x(2-grd%ias:grd%im,:,k) * 0.5
     enddo


  DEALLOCATE ( ud, vd )



end subroutine get_vel_ad
