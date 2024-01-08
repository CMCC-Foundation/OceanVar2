subroutine div_dmp_ad

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
! Divergence damping (adjoint)
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
!-----------------------------------------------------------------------


  use set_knd
  use drv_str
  use grd_str
  use mpi_str

  implicit none

  INTEGER(i4)    :: k, kdiv, i, j
  REAL(r8), ALLOCATABLE   :: div(:,:,:)

  ALLOCATE ( div(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))


   do kdiv = 1,100

    div(:,:,:) = 0.0

    do k=grd%km,1,-1
      do j=2-grd%jas,grd%jm-1+grd%jae
       do i=1,grd%im-1+grd%iae
        div(i,j  ,k) = div(i,j  ,k) + grd%vvl(i,j,k)*0.2 * grd%adxdy**2 / grd%dy(i,j) * grd%msk(i,j,k)*grd%msk(i,j-1,k)
        div(i,j-1,k) = div(i,j-1,k) + grd%vvl(i,j,k)*0.2 * grd%adxdy**2 / grd%dy(i,j) * grd%msk(i,j,k)*grd%msk(i,j-1,k)
       enddo
      enddo
      do j=1,grd%jm-1+grd%jae
       do i=2-grd%ias,grd%im-1+grd%iae
        div(i  ,j,k) = div(i  ,j,k) + grd%uvl(i,j,k)*0.2 * grd%adxdy**2 / grd%dx(i,j) * grd%msk(i,j,k)*grd%msk(i-1,j,k)
        div(i-1,j,k) = div(i-1,j,k) + grd%uvl(i,j,k)*0.2 * grd%adxdy**2 / grd%dx(i,j) * grd%msk(i,j,k)*grd%msk(i-1,j,k)
       enddo
      enddo
    enddo

 if(mpi%nproc.gt.1)then
  call exo_mpi( 0_i4, 0_i4,                                       &
               -1_i4, 0_i4, grd%im,                               &
               -1_i4, 0_i4, grd%jm,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, grd%km, div)
 endif

    do k=grd%km,1,-1
        if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm) grd%uvl(grd%im+1,:,:) = 0.0
        if(mpi%nproc.gt.1 .and. mpi%jr.lt.mpi%jrm) grd%vvl(:,grd%jm+1,:) = 0.0
      do j=1,grd%jm-1+grd%jae
       do i=1,grd%im-1+grd%iae
        grd%uvl(i  ,j  ,k) = grd%uvl(i  ,j  ,k) - div(i,j,k) * grd%dy(i  ,j  ) / grd%dx(i,j) / grd%dy(i,j)
        grd%uvl(i+1,j  ,k) = grd%uvl(i+1,j  ,k) + div(i,j,k) * grd%dy(i+1,j  ) / grd%dx(i,j) / grd%dy(i,j)
        grd%vvl(i  ,j  ,k) = grd%vvl(i  ,j  ,k) - div(i,j,k) * grd%dx(i  ,j  ) / grd%dy(i,j) / grd%dx(i,j)
        grd%vvl(i  ,j+1,k) = grd%vvl(i  ,j+1,k) + div(i,j,k) * grd%dx(i  ,j+1) / grd%dy(i,j) / grd%dx(i,j)
       enddo
      enddo
    enddo

 if(mpi%nproc.gt.1)then
  call exo_mpi( 0_i4, 1_i4,                                       &
                1_i4, grd%im+1, 1_i4,                             &
                0_i4, 1_i4, grd%jm+1,                             &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, grd%km, grd%uvl)
  call exo_mpi( 0_i4, 1_i4,                                       &
                0_i4, 1_i4, grd%im+1,                             &
                1_i4, grd%jm+1, 1_i4,                             &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, grd%km, grd%vvl)
 endif

  enddo

   DEALLOCATE (div)

end subroutine div_dmp_ad
