subroutine invrt( kstp )

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
! Implicit solver - overrelaxation                                     !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
! Version 1.1:  P. Oddo 2014                                           !
!-----------------------------------------------------------------------


  use set_knd
  use bmd_str
  use grd_str
  use mpi_str

  implicit none

  include 'mpif.h'

  REAL(r8), ALLOCATABLE  :: res(:,:)

  REAL(r8)       :: reser, resem, reserp
  INTEGER(i4)    :: i, j, icnt, kstp
  INTEGER        :: ierr


  ALLOCATE (res(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )

! MPI exchange
   if(mpi%nproc.gt.1)                                                                  &
   call exa_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%eta)

   reser = 1.e20_r8

   bmd%itr(kstp) = 0_r8

      do icnt=1,bmd%ncnt      ! DO CONVERGENCE LOOP

       if(reser.gt.bmd%resem)then  ! IF CHECK RESIDUAL

         bmd%itr(kstp) = bmd%itr(kstp) + 1
         ! ------------------------------------------------------------------------------
         do j=grd%js2,grd%jm-1+grd%jae,2          ! 2,jm-1,2
          do i=grd%is2,grd%im-1+grd%iae,2         ! 2,jm-1,2
            res(i,j) = bmd%a1(i,j)*bmd%eta(i+1,j  )+bmd%a2(i,j)*bmd%eta(i-1,j  )+        &
                       bmd%a3(i,j)*bmd%eta(i  ,j+1)+bmd%a4(i,j)*bmd%eta(i  ,j-1)-        &
                       bmd%a0(i,j)*bmd%eta(i  ,j  ) - bmd%rgh(i,j)
            bmd%eta(i,j) = bmd%eta(i,j)+bmd%ovr*res(i,j)/bmd%a0(i,j)*bmd%mst(i,j)
          enddo
         enddo

         do j=grd%js3,grd%jm-1+grd%jae,2          ! 3,jm-1,2
          do i=grd%is3,grd%im-1+grd%iae,2         ! 3,im-1,2
            res(i,j) = bmd%a1(i,j)*bmd%eta(i+1,j  )+bmd%a2(i,j)*bmd%eta(i-1,j  )+        &
                       bmd%a3(i,j)*bmd%eta(i  ,j+1)+bmd%a4(i,j)*bmd%eta(i  ,j-1)-        &
                       bmd%a0(i,j)*bmd%eta(i  ,j  ) - bmd%rgh(i,j)
             bmd%eta(i,j) = bmd%eta(i,j)+bmd%ovr*res(i,j)/bmd%a0(i,j)*bmd%mst(i,j)
          enddo
         enddo
         ! ------------------------------------------------------------------------------

! MPI exchange
   if(mpi%nproc.gt.1)                                                                  &
   call exa_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%eta)

         ! ------------------------------------------------------------------------------
         do j=grd%js3,grd%jm-1+grd%jae,2          ! 3,jm-1,2
          do i=grd%is2,grd%im-1+grd%iae,2         ! 2,im-1,2
            res(i,j) = bmd%a1(i,j)*bmd%eta(i+1,j  )+bmd%a2(i,j)*bmd%eta(i-1,j  )+        &
                       bmd%a3(i,j)*bmd%eta(i  ,j+1)+bmd%a4(i,j)*bmd%eta(i  ,j-1)-        &
                       bmd%a0(i,j)*bmd%eta(i  ,j  ) - bmd%rgh(i,j)
            bmd%eta(i,j) = bmd%eta(i,j)+bmd%ovr*res(i,j)/bmd%a0(i,j)*bmd%mst(i,j)
          enddo
         enddo

         do j=grd%js2,grd%jm-1+grd%jae,2          ! 2,jm-1,2
          do i=grd%is3,grd%im-1+grd%iae,2         ! 3,im-1,2
            res(i,j) = bmd%a1(i,j)*bmd%eta(i+1,j  )+bmd%a2(i,j)*bmd%eta(i-1,j  )+        &
                       bmd%a3(i,j)*bmd%eta(i  ,j+1)+bmd%a4(i,j)*bmd%eta(i  ,j-1)-        &
                      bmd%a0(i,j)*bmd%eta(i  ,j  ) - bmd%rgh(i,j)
             bmd%eta(i,j) = bmd%eta(i,j)+bmd%ovr*res(i,j)/bmd%a0(i,j)*bmd%mst(i,j)
          enddo
         enddo
         ! ------------------------------------------------------------------------------

! MPI exchange
   if(mpi%nproc.gt.1)                                                                  &
   call exa_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%eta)

       reser = 0.0_r8
       do j=2-grd%jas,grd%jm-1+grd%jae                  ! 2,jm-1
        do i=2-grd%ias,grd%im-1+grd%iae                 ! 2,im-1
           reser = reser + abs(res(i,j))
        enddo
       enddo
       if(bmd%bnm.gt.0.0_r8) reser = reser/bmd%bnm

 if(mpi%nproc.gt.1)then
  reserp = reser
  call mpi_allreduce( reserp, reser, 1, mpi%r8, mpi_sum, mpi%comm, ierr)
 endif

       ! write(*,*) ' Barotropic Model Diagnostic '
       ! write(*,*) ' At Iteraction n ', bmd%itr(kstp), icnt
       ! write(*,*) ' Residual is ', reser, ' limit=', bmd%resem

       endif                  ! IF CHECK RESIDUAL

      enddo                   ! DO CONVERGENCE LOOP

  DEALLOCATE (res)

end subroutine invrt
