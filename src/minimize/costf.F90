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
!> Calclate the cost function and its gradient                          
!!
!! It computes:
!! 1) background cost term
!! 2) observational term
!! 3) total cost function
!! 4) adjoint routines to get the gradient
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!            Mario Adani and Francesco Carere 2024                     !
!-----------------------------------------------------------------------
SUBROUTINE costf

   USE set_knd
   USE obs_str
   USE grd_str
   USE eof_str
   USE ctl_str
   USE mpi_str
   USE drv_str

#ifdef REPRO
   USE rpr_str
#endif
   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER                                :: ierr, i, j
   REAL(r8)                               :: mpism,mpismo,mpismb
   REAL(r8),    ALLOCATABLE, DIMENSION(:) :: all_vals
   INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: all_nvals, disp_nvals
   INTEGER(i4)                            :: tot_nvals
#ifdef REPRO
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sort_array
#endif

! -------------------------------------------------------
! 1) calculate backgorund cost term
!-------------------------------------------------------
#ifndef REPRO

   ctl%f_b = 0.5_r8 * DOT_PRODUCT( ctl%x_c, ctl%x_c)

#else

   IF ( mpi%nproc .GT. 1 ) THEN            ! Parallel

      ! allocate recvcounts and displs arrays
      IF ( mpi%myrank==0 ) THEN
         ALLOCATE ( all_nvals(mpi%nproc) )
         ALLOCATE ( disp_nvals(mpi%nproc) )
      ELSE
         ALLOCATE ( all_nvals(1) )
         ALLOCATE ( disp_nvals(1) )
      ENDIF
      ! Gather the number of data to receive from each process
      CALL MPI_GATHER(ctl%n, 1, mpi%i4, all_nvals, 1, mpi%i4, 0, mpi%comm, ierr)
      ! Calculate displacements
      IF ( mpi%myrank==0 ) THEN
         tot_nvals = 0
         DO i = 1, mpi%nproc
            disp_nvals(i) = tot_nvals
            tot_nvals     = tot_nvals + all_nvals(i)
         ENDDO
         ! allocate receive buffer
         ALLOCATE ( all_vals(tot_nvals) )
         ALLOCATE ( sort_array(tot_nvals) )
      ELSE
         ALLOCATE( all_vals(1) )
         ALLOCATE( sort_array(1) )
      ENDIF
      ! Gather the final data
      all_vals = 0.0_r8
      CALL MPI_GATHERV(ctl%x_c, ctl%n, mpi%r8, all_vals, all_nvals, disp_nvals, mpi%r8, 0, mpi%comm, ierr)
      ! Perform final dot product to compute the bkg component of the cost function
      IF ( mpi%myrank==0 ) THEN
         CALL gth_to_seq(all_vals,sort_array)
         ctl%f_b = 0.5_r8 * DOT_PRODUCT( sort_array, sort_array)
      ELSE
         ctl%f_b = 0.0_r8
      ENDIF
      DEALLOCATE ( all_vals,all_nvals,disp_nvals,sort_array )

   ELSE                                      ! Single processor

      ctl%f_b= 0.5_r8 * DOT_PRODUCT( ctl%x_c, ctl%x_c)

   ENDIF

#endif

! -------------------------------------------------------
! 2) calculate observational cost term
! -------------------------------------------------------
! --------
! Convert the control vector to v
   CALL cnv_ctv
! --------
! Control to physical space
   CALL ver_hor
! --------
! Apply observational operators
   CALL obsop
! --------
! Calculate residuals
   CALL resid
! --------
!  Perform final dot product to compute the obs component of the cost function
#ifndef REPRO
   IF ( (ctl%isave(30) .GT. huberqc%iter) .AND.  &
      (huberqc%any) ) THEN
      CALL huber_costf(ctl%f_o)
   ELSE
      ctl%f_o = 0.5_r8 * DOT_PRODUCT( obs%amo, obs%amo)
   ENDIF

! -------------------------------------------------------
! 3) calculate total cost function: bkg + obs
! -------------------------------------------------------
   ctl%f_c = ctl%f_b + ctl%f_o

   IF ( mpi%nproc .GT. 1 ) THEN
      CALL MPI_REDUCE( ctl%f_c, mpism, 1, mpi%r8,      &
         MPI_SUM, 0, mpi%comm, ierr)
      CALL MPI_REDUCE( ctl%f_o, mpismo, 1, mpi%r8,      &
         MPI_SUM, 0, mpi%comm, ierr)
      CALL MPI_REDUCE( ctl%f_b, mpismb, 1, mpi%r8,      &
         MPI_SUM, 0, mpi%comm, ierr)
      IF ( mpi%myrank==0 ) THEN
         ctl%f_c = mpism
         ctl%f_o = mpismo
         ctl%f_b = mpismb
      ENDIF
      CALL mpi_bcast( ctl%f_c, 1, mpi%r8, 0, mpi%comm, ierr)
   ENDIF

#else

   ! Perform final dot product to compute the obs component of the cost function
   IF ( mpi%nproc .GT. 1 ) THEN
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, rpr%obs_amo,rpr%obs_nog, mpi%r8, MPI_SUM, mpi%comm, ierr)
      ctl%f_o= 0.5_r8 * DOT_PRODUCT(rpr%obs_amo, rpr%obs_amo)
   ELSE
      ctl%f_o = 0.5_r8 * DOT_PRODUCT( obs%amo, obs%amo)
   ENDIF

   IF ( mpi%nproc .GT. 1 ) THEN
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, ctl%f_b, 1, mpi%r8, MPI_SUM, mpi%comm, ierr)
   ENDIF
   ctl%f_c = ctl%f_b + ctl%f_o

#endif
   IF ( mpi%myrank==0 ) WRITE (drv%dia,*)' -------------------------------------'
   IF ( mpi%myrank==0 ) WRITE (drv%dia,*)' Cost Function   ',ctl%f_c
   IF ( mpi%myrank==0 ) WRITE (drv%dia,*)' Modl Background ',ctl%f_b
   IF ( mpi%myrank==0 ) WRITE (drv%dia,*)' Obsr Background ',ctl%f_o

! -------------------------------------------------------
! calculate the cost function gradient
! -------------------------------------------------------

! --------
! Reset the increments
   CALL res_inc

! --------
! Observational OPERATORs
   CALL obsop_ad

! --------
! Control to physical space
   CALL ver_hor_ad
! --------
! Convert the control vector
   CALL cnv_ctv_ad

! -------------------------------------------------------
! Cost function gradient
! -------------------------------------------------------
   ctl%g_c(:) = ctl%x_c(:) + ctl%g_c(:)

END SUBROUTINE costf
