subroutine div_dmp

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
! Divergence damping of velocity fields
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
!-----------------------------------------------------------------------


  use set_knd
  use drv_str
  use grd_str
  use mpi_str

  implicit none

  REAL(r8), ALLOCATABLE   :: div(:,:,:)
  INTEGER(i4)    :: k, kdiv, i, j


  ALLOCATE ( div(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))

   do kdiv = 1,100

 if(mpi%nproc.gt.1)then
  call exo_mpi( 0_i4, 1_i4,                                       &
               -1_i4, 1_i4, grd%im+1,                             &
                0_i4, 1_i4, grd%jm+1,                             &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, grd%km, grd%uvl)
  call exo_mpi( 0_i4, 1_i4,                                       &
                0_i4, 1_i4, grd%im+1,                             &
               -1_i4, 1_i4, grd%jm+1,                             &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, grd%km, grd%vvl)
 endif


    do k=1,grd%km
      do j=1,grd%jm-1+grd%jae
       do i=1,grd%im-1+grd%iae
        div(i,j,k) = ( grd%uvl(i+1,j  ,k)*grd%dy(i+1,j  ) - grd%uvl(i  ,j  ,k)*grd%dy(i  ,j  ) +       &
                       grd%vvl(i  ,j+1,k)*grd%dx(i  ,j+1) - grd%vvl(i  ,j  ,k)*grd%dx(i  ,j  ) )       &
                     / grd%dx(i,j) / grd%dy(i,j)
       enddo
      enddo
    enddo

 if(mpi%nproc.gt.1)then
  call exo_mpi( 0_i4, 1_i4,                                       &
                1_i4, grd%im, 0_i4,                               &
                1_i4, grd%jm, 0_i4,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, grd%km, div)
 endif

    do k=1,grd%km
      do j=1,grd%jm-1+grd%jae
       do i=2-grd%ias,grd%im-1+grd%iae
        grd%uvl(i,j,k) = grd%uvl(i,j,k)  + 0.2  * grd%adxdy**2 * (div(i  ,j  ,k) - div(i-1,j  ,k))    &
                         / grd%dx(i,j) * grd%msk(i,j,k)*grd%msk(i-1,j  ,k)
       enddo
      enddo
      do j=2-grd%jas,grd%jm-1+grd%jae
       do i=1,grd%im-1+grd%iae
       grd%vvl(i,j,k) = grd%vvl(i,j,k)  + 0.2  * grd%adxdy**2 * (div(i  ,j  ,k) - div(i  ,j-1,k))    &
                        / grd%dy(i,j) * grd%msk(i,j,k)*grd%msk(i  ,j-1,k)
       enddo
      enddo
    enddo


   enddo

  DEALLOCATE ( div )

end subroutine div_dmp
