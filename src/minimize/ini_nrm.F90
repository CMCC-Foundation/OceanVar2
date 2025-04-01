!======================================================================
!
! This file is part of Oceanvar.
!
!  Copyright (C) 2025 OceanVar System Team ( oceanvar@cmcc.it )
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! any later version (GPL-3.0-or-later).
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!======================================================================
!-----------------------------------------------------------------------
!                                                                      !
! Inizialization cost function                                         !
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE ini_nrm

   USE set_knd
   USE drv_str
   USE obs_str
   USE grd_str
   USE eof_str
   USE ctl_str
   USE mpi_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'
   INTEGER                             :: ierr
   REAL(r8)                            :: mpism

   INTEGER(i4)                         :: k
   REAL(r8)                            :: maxpg
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: x_s, g_s

   ALLOCATE ( x_s(ctl%n), g_s(ctl%n) )

! Save in temprary arrays
   x_s(:) = ctl%x_c(:)
   g_s(:) = ctl%g_c(:)

! Initialize the control vector to zero
   ctl%x_c(:) =  0.0_r8

! Calculate the cost function and its gradient
   CALL costf

! Calculate the norm of the gradient
   maxpg = 0.0_r8
   DO k=1,ctl%n
      maxpg = maxpg+ctl%g_c(k)**2
   ENDDO

   IF ( mpi%nproc .GT. 1 ) THEN

      CALL MPI_REDUCE( maxpg, mpism, 1, mpi%r8,       &
                       MPI_SUM, 0, mpi%comm, ierr)
      IF ( mpi%myrank==0 ) THEN
         maxpg = mpism
      ENDIF
      CALL mpi_bcast(maxpg, 1, mpi%r8, 0, mpi%comm, ierr)

   ENDIF
   maxpg = sqrt(maxpg)

   ctl%pgtol = maxpg * ctl%pgper

   ctl%x_c(:) = x_s(:)
   ctl%g_c(:) = g_s(:)

   DEALLOCATE ( x_s, g_s )

END SUBROUTINE ini_nrm
