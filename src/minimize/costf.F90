subroutine costf

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
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
! Calclate the cost function and its gradient                          !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use obs_str
 use grd_str
 use eof_str
 use ctl_str
 use mpi_str
 use drv_str

 implicit none

 include 'mpif.h'

  integer       :: ierr, i, j
  real(r8)      :: mpism

! -------------------------------------------------------
! calculate backgorund cost term
! -------------------------------------------------------

    ctl%f_b = 0.5_r8 * dot_product( ctl%x_c, ctl%x_c)

! -------------------------------------------------------
! calculate observational cost term
! -------------------------------------------------------
! --------
! Convert the control vector to v
   call cnv_ctv

! --------
! Control to physical space 
   call ver_hor

! --------
! Apply observational operators
   call obsop

! --------
! Calculate residuals
   call resid

! --------
! calculate cost
    ctl%f_o = 0.5_r8 * dot_product( obs%amo, obs%amo)

! -------------------------------------------------------
! Cost function
! -------------------------------------------------------

    ctl%f_c = ctl%f_b + ctl%f_o

  if(mpi%nproc.gt.1) then
    call mpi_reduce( ctl%f_c, mpism, 1, mpi%r8,      &
                      mpi_sum, 0, mpi%comm, ierr)
    if(mpi%myrank==0)then
       ctl%f_c = mpism
    endif
    call mpi_bcast( ctl%f_c, 1, mpi%r8, 0, mpi%comm, ierr)
  endif
    if(mpi%myrank==0) write(drv%dia,*)' -------------------------------------'
    if(mpi%myrank==0) write(drv%dia,*)' Cost function ',ctl%f_c
    if(mpi%myrank==0) write(drv%dia,*)' Modl Background ',ctl%f_b
    if(mpi%myrank==0) write(drv%dia,*)' Obsr Background ',ctl%f_o

! -------------------------------------------------------
! calculate the cost function gradient
! -------------------------------------------------------

! --------
! Reset the increments
   call res_inc

! --------
! Observational operators
   call obsop_ad

! --------
! Control to physical space 
   call ver_hor_ad

! --------
! Convert the control vector 
   call cnv_ctv_ad

! -------------------------------------------------------
! Cost function gradient
! -------------------------------------------------------

   ctl%g_c(:) = ctl%x_c(:) + ctl%g_c(:)

         
end subroutine costf
