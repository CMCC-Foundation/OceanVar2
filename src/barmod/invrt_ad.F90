subroutine invrt_ad(kstp)

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
! Implicit solver - overrelaxation (adjoint)                           !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
!-----------------------------------------------------------------------


  use set_knd
  use grd_str
  use bmd_str
  use mpi_str

  implicit none

 INTEGER(i4)    :: kstp

 REAL(r8), allocatable  :: res(:,:)

 INTEGER(i4)    :: i, j, icnt

!  ALLOCATE (res(grd%im,grd%jm))
  ALLOCATE (res(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )

   res(:,:) = 0.0_r8

      do icnt=bmd%itr(kstp),1,-1

! MPI exchange
!   if(mpi%nproc.gt.1 .and. icnt.lt.bmd%itr(kstp))                                                                  &
!       call exa_mpi( 1_i4, 0_i4, grd%im+1, 1_i4, grd%im, 0_i4, grd%jm+1, 1_i4, grd%jm,     &
!                      1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%eta)

          if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm) bmd%eta(grd%im+1,:) = 0.0_r8
          if(mpi%nproc.gt.1 .and. mpi%jr.lt.mpi%jrm) bmd%eta(:,grd%jm+1) = 0.0_r8
          if(mpi%nproc.gt.1 .and. mpi%jr.gt.1      ) bmd%eta(:,0)        = 0.0_r8
          if(mpi%nproc.gt.1 .and. mpi%ir.gt.1      ) bmd%eta(0,:)        = 0.0_r8

         do j=grd%js2,grd%jm-1+grd%jae,2          ! 2,jm-1,2
          do i=grd%is3,grd%im-1+grd%iae,2           ! 3,im-1,2
            res(i,j) = res(i,j) + bmd%eta(i,j)*bmd%ovr/bmd%a0(i,j)*bmd%mst(i,j)
            bmd%eta(i,j-1) = bmd%eta(i,j-1) + res(i,j)*bmd%a4(i,j)
            bmd%eta(i,j+1) = bmd%eta(i,j+1) + res(i,j)*bmd%a3(i,j)
            bmd%eta(i-1,j) = bmd%eta(i-1,j) + res(i,j)*bmd%a2(i,j)
            bmd%eta(i+1,j) = bmd%eta(i+1,j) + res(i,j)*bmd%a1(i,j)
            bmd%eta(i,j) = bmd%eta(i,j) - res(i,j)*bmd%a0(i,j)
            bmd%rgh(i,j) = bmd%rgh(i,j) - res(i,j)
            res(i,j) = 0.0_r8
          enddo
         enddo
         do j=grd%js3,grd%jm-1+grd%jae,2          ! 3,jm-1,2
          do i=grd%is2,grd%im-1+grd%iae,2           ! 2,im-1,2
            res(i,j) = res(i,j) + bmd%eta(i,j)*bmd%ovr/bmd%a0(i,j)*bmd%mst(i,j)
            bmd%eta(i,j-1) = bmd%eta(i,j-1) + res(i,j)*bmd%a4(i,j)
            bmd%eta(i,j+1) = bmd%eta(i,j+1) + res(i,j)*bmd%a3(i,j)
            bmd%eta(i-1,j) = bmd%eta(i-1,j) + res(i,j)*bmd%a2(i,j)
            bmd%eta(i+1,j) = bmd%eta(i+1,j) + res(i,j)*bmd%a1(i,j)
            bmd%eta(i,j) = bmd%eta(i,j) - res(i,j)*bmd%a0(i,j)
            bmd%rgh(i,j) = bmd%rgh(i,j) - res(i,j)
            res(i,j) = 0.0_r8
          enddo
         enddo

! MPI exchange
    if(mpi%nproc.gt.1)                                                                  &
        call exa_mpi( 0_i4, 0_i4, grd%im+1, 1_i4, grd%im, 0_i4, grd%jm+1, 1_i4, grd%jm,     &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%eta)

          if(mpi%nproc.gt.1 .and. mpi%ir.lt.mpi%irm) bmd%eta(grd%im+1,:) = 0.0_r8
          if(mpi%nproc.gt.1 .and. mpi%jr.lt.mpi%jrm) bmd%eta(:,grd%jm+1) = 0.0_r8
          if(mpi%nproc.gt.1 .and. mpi%jr.gt.1      ) bmd%eta(:,0)        = 0.0_r8
          if(mpi%nproc.gt.1 .and. mpi%ir.gt.1      ) bmd%eta(0,:)        = 0.0_r8

         do j=grd%js3,grd%jm-1+grd%jae,2          ! 3,jm-1,2
          do i=grd%is3,grd%im-1+grd%iae,2         ! 3,im-1,2
            res(i,j) = res(i,j) + bmd%eta(i,j)*bmd%ovr/bmd%a0(i,j)*bmd%mst(i,j)
            bmd%eta(i,j-1) = bmd%eta(i,j-1) + res(i,j)*bmd%a4(i,j)
            bmd%eta(i,j+1) = bmd%eta(i,j+1) + res(i,j)*bmd%a3(i,j)
            bmd%eta(i-1,j) = bmd%eta(i-1,j) + res(i,j)*bmd%a2(i,j)
            bmd%eta(i+1,j) = bmd%eta(i+1,j) + res(i,j)*bmd%a1(i,j)
            bmd%eta(i,j) = bmd%eta(i,j) - res(i,j)*bmd%a0(i,j)
            bmd%rgh(i,j) = bmd%rgh(i,j) - res(i,j)
            res(i,j) = 0.0_r8
          enddo
         enddo
         do j=grd%js2,grd%jm-1+grd%jae,2          ! 2,jm-1,2
          do i=grd%is2,grd%im-1+grd%iae,2         ! 2,jm-1,2
            res(i,j) = res(i,j) + bmd%eta(i,j)*bmd%ovr/bmd%a0(i,j)*bmd%mst(i,j)
            bmd%eta(i,j-1) = bmd%eta(i,j-1) + res(i,j)*bmd%a4(i,j)
            bmd%eta(i,j+1) = bmd%eta(i,j+1) + res(i,j)*bmd%a3(i,j)
            bmd%eta(i-1,j) = bmd%eta(i-1,j) + res(i,j)*bmd%a2(i,j)
            bmd%eta(i+1,j) = bmd%eta(i+1,j) + res(i,j)*bmd%a1(i,j)
            bmd%eta(i,j) = bmd%eta(i,j) - res(i,j)*bmd%a0(i,j)
            bmd%rgh(i,j) = bmd%rgh(i,j) - res(i,j)
            res(i,j) = 0.0_r8
          enddo
         enddo

            bmd%eta(:,:) = bmd%eta(:,:) * bmd%mst(:,:)

    if(mpi%nproc.gt.1)                                                                  &
        call exa_mpi( 0_i4, 0_i4, grd%im+1, 1_i4, grd%im, 0_i4, grd%jm+1, 1_i4, grd%jm,     &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%eta)

      enddo


  DEALLOCATE (res)

end subroutine invrt_ad
