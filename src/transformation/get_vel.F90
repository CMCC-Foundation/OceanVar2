subroutine get_vel

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
! Calculate horizontal velocity from geostrophic formula               !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
! Bug correction 21.04.2009  thanks to Andrea Storto                   !
!-----------------------------------------------------------------------


  use set_knd
  use drv_str
  use grd_str
  use mpi_str

  implicit none

  include 'mpif.h'

  REAL(r8), ALLOCATABLE, DIMENSION (:,:,:)  :: ud, vd

  integer(i4)    :: k, i, j, lr


  ALLOCATE ( ud(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
  ALLOCATE ( vd(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))

 if(mpi%nproc.gt.1)then
  call exo_mpi( 0_i4, 1_i4,                                       &
                1_i4, grd%im, 0_i4,                               &
                1_i4, grd%jm, 0_i4,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, grd%eta)
 endif

     do k=grd%km,2,-1
      grd%b_x(:,:,k) = ( grd%b_x(:,:,k) + grd%b_x(:,:,k-1) ) * 0.5
      grd%b_y(:,:,k) = ( grd%b_y(:,:,k) + grd%b_y(:,:,k-1) ) * 0.5
     enddo
      grd%b_x(:,:,1) = grd%b_x(:,:,1) * 0.5
      grd%b_y(:,:,1) = grd%b_y(:,:,1) * 0.5

     grd%uvl(:,:,:) = 0.0
     grd%vvl(:,:,:) = 0.0

      ud(:,:,:) = 0.0
      vd(:,:,:) = 0.0

    do k=1,grd%km

      do j=1,grd%jm
       do i=2-grd%ias,grd%im
        vd(i,j,k) =   ( (grd%eta(i,j)-grd%eta(i-1,j))*9.81 + grd%b_x(i,j,k) )               &
                             / grd%dx(i,j) * grd%msk(i,j,k)*grd%msk(i-1,j,k) / grd%f(i,j)
       enddo
      enddo
      do j=2-grd%jas,grd%jm
       do i=1,grd%im
        ud(i,j,k) = - ( (grd%eta(i,j)-grd%eta(i,j-1))*9.81 + grd%b_y(i,j,k) )               &
                             / grd%dy(i,j) * grd%msk(i,j,k)*grd%msk(i,j-1,k) / grd%f(i,j)
       enddo
      enddo

    enddo

 if(mpi%nproc.gt.1)then
  call exo_mpi( 1_i4, 1_i4,                                       &
                1_i4, grd%im, 0_i4,                               &
               -1_i4, 1_i4, grd%jm+1,                             &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, ud)
  call exo_mpi( 1_i4, 1_i4,                                       &
               -1_i4, 1_i4, grd%im+1,                             &
                1_i4, grd%jm, 0_i4,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, vd)
 endif

    do k=1,grd%km
      do j=1,grd%jm-1+grd%jae
       do i=2-grd%ias,grd%im
        grd%uvl(i,j,k) = ( ud(i,j,k)+ud(i-1,j,k)+ud(i,j+1,k)+ud(i-1,j+1,k) )*0.25 * grd%msk(i,j,k)*grd%msk(i-1,j,k)
       enddo
      enddo
      do j=2-grd%jas,grd%jm      ! 2:grd%jm
       do i=1,grd%im-1+grd%iae   ! 1:grd%im-1
        grd%vvl(i,j,k) = ( vd(i,j,k)+vd(i,j-1,k)+vd(i+1,j,k)+vd(i+1,j-1,k) )*0.25  * grd%msk(i,j,k)*grd%msk(i,j-1,k)
       enddo
      enddo
    enddo

  DEALLOCATE ( ud, vd )


end subroutine get_vel
