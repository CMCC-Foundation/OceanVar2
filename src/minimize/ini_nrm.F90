subroutine ini_nrm


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
! Minimise the cost function                                           !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use drv_str
 use obs_str
 use grd_str
 use eof_str
 use ctl_str
 use mpi_str

 implicit none

 include 'mpif.h'
 INTEGER               :: ierr
 REAL(r8)              :: mpism

 INTEGER(i4)           :: k
 REAL(r8)              :: maxpg
 REAL(r8), allocatable, dimension (:) :: x_s, g_s


 ALLOCATE ( x_s(ctl%n), g_s(ctl%n) )

! Save in temprary arrays
          x_s(:) = ctl%x_c(:)
          g_s(:) = ctl%g_c(:)

! Initialize the control vector to zero
          ctl%x_c(:) =  0.0

! Calculate the cost function and its gradient
          call costf

! Calculate the norm of the gradient
            maxpg = 0.0
           do k=1,ctl%n
!sd            maxpg = max(maxpg,abs(ctl%g_c(k)))
            maxpg = maxpg+ctl%g_c(k)**2
           enddo


    if(mpi%nproc.gt.1) then
      call mpi_reduce( maxpg, mpism, 1, mpi%r8,       &
                      mpi_sum, 0, mpi%comm, ierr)
!sd                      mpi_max, 0, mpi%comm, ierr)
      if(mpi%myrank==0)then
           maxpg = mpism
      endif
      call mpi_bcast(maxpg, 1, mpi%r8, 0, mpi%comm, ierr)
    endif
!sd
           maxpg = sqrt(maxpg)
!sd

            ctl%pgtol = maxpg * ctl%pgper

     if(mpi%myrank==0) print*,' maxpg  is: ',maxpg, ctl%pgtol

          ctl%x_c(:) = x_s(:)
          ctl%g_c(:) = g_s(:)

 DEALLOCATE ( x_s, g_s )

end subroutine ini_nrm
