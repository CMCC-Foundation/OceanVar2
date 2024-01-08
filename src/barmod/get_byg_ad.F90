subroutine get_byg_ad

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
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
! Calculate vertical integral of bouyancy gradient (adjoint)           !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
!-----------------------------------------------------------------------


  use set_knd
  use drv_str
  use grd_str
  use bmd_str
  use mpi_str

  implicit none

  integer(i4)    :: i, j, k

      grd%dns(:,:,:) = 0.0

    do j=2-grd%jas,grd%jm
     do i=2-grd%ias,grd%im
      grd%bx(i,j) = grd%bx(i,j)/grd%dx(i,j)
      grd%by(i,j) = grd%by(i,j)/grd%dy(i,j)
     enddo
    enddo

   do k=1,grd%km
    do j=2-grd%jas,grd%jm
     do i=2-grd%ias,grd%im
      grd%b_x(i,j,k) = grd%b_x(i,j,k) + grd%bx(i,j)*grd%dz(k) * grd%msk(i,j,k)*grd%msk(i-1,j,k)
      grd%b_y(i,j,k) = grd%b_y(i,j,k) + grd%by(i,j)*grd%dz(k) * grd%msk(i,j,k)*grd%msk(i,j-1,k)
     enddo
    enddo
   enddo

   do k=grd%km,2,-1
    do j=2-grd%jas,grd%jm
     do i=2-grd%ias,grd%im
      grd%dns(i,j  ,k) = grd%dns(i,j  ,k) +               &
                                grd%b_y(i,j,k)*grd%dz(k)*grd%msk(i,j,k)*grd%msk(i,j-1,k)
      grd%dns(i,j-1,k) = grd%dns(i,j-1,k) -               &
                                grd%b_y(i,j,k)*grd%dz(k)*grd%msk(i,j,k)*grd%msk(i,j-1,k)
      grd%b_y(i,j,k-1) = grd%b_y(i,j,k-1) + grd%b_y(i,j,k)
      grd%b_y(i,j,k)   = 0.0
     enddo
    enddo
   enddo
    do j=2-grd%jas,grd%jm
     do i=2-grd%ias,grd%im
      grd%dns(i,j  ,1) = grd%dns(i,j  ,1) +               &
                                grd%b_y(i,j,1)*grd%dz(1)*grd%msk(i,j,1)*grd%msk(i,j-1,1)
      grd%dns(i,j-1,1) = grd%dns(i,j-1,1) -               &
                                grd%b_y(i,j,1)*grd%dz(1)*grd%msk(i,j,1)*grd%msk(i,j-1,1)
     enddo
    enddo

   do k=grd%km,2,-1
    do j=2-grd%jas,grd%jm
     do i=2-grd%ias,grd%im
      grd%dns(i  ,j,k) = grd%dns(i  ,j,k) +               &
                                grd%b_x(i  ,j,k)*grd%dz(k)*grd%msk(i  ,j,k)*grd%msk(i-1,j,k)
      grd%dns(i-1,j,k) = grd%dns(i-1,j,k) -               &
                                grd%b_x(i  ,j,k)*grd%dz(k)*grd%msk(i  ,j,k)*grd%msk(i-1,j,k)
      grd%b_x(i  ,j,k-1) = grd%b_x(i  ,j,k-1) + grd%b_x(i  ,j,k)
      grd%b_x(i  ,j,k)   = 0.0
     enddo
    enddo
   enddo
    do j=2-grd%jas,grd%jm
     do i=2-grd%ias,grd%im
      grd%dns(i  ,j,1) = grd%dns(i  ,j,1) +               &
                                grd%b_x(i  ,j,1)*grd%dz(k)*grd%msk(i  ,j,1)*grd%msk(i-1,j,1)
      grd%dns(i-1,j,1) = grd%dns(i-1,j,1) -               &
                                grd%b_x(i  ,j,1)*grd%dz(k)*grd%msk(i  ,j,1)*grd%msk(i-1,j,1)
     enddo
    enddo

 if(mpi%nproc.gt.1)then

  call exo_mpi( 0_i4, 0_i4,                                        &
                -1_i4, 0_i4, grd%im,                               &
                -1_i4, 0_i4, grd%jm,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%dns)

 endif

      grd%tem_ad(:,:,:) = grd%tem_ad(:,:,:) - 0.24*grd%dns(:,:,:) * bmd%g/1025. * grd%msk(:,:,:)
      grd%sal_ad(:,:,:) = grd%sal_ad(:,:,:) + 0.74*grd%dns(:,:,:) * bmd%g/1025. * grd%msk(:,:,:)

end subroutine get_byg_ad
