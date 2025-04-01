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
!> Apply observational operators                                       
!!
!! 
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE obsop

   USE set_knd
   USE obs_str
   USE grd_str
   USE mpi_str

   IMPLICIT NONE

! ---
   IF ( mpi%nproc .GT. 1 ) THEN
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                   1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae ,   1_i4, grd%eta)
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                   1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%tem)
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                   1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%sal)
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                   1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%uvl)
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                   1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%vvl)
   ENDIF

! ---
! Satellite observations of SLA
   CALL obs_sla

! ---
! Observations by ARGO floats
   CALL obs_arg

! ---
! Observations by XBT profiles
   CALL obs_xbt

! ---
! Observations by gliders
   CALL obs_gld

! ---
! Observations of Argo trajectories
   IF ( tra%no .GT. 0 ) CALL obs_tra

! ---
! Observations of trajectories of surface drIFters
   IF ( trd%no .GT. 0 ) CALL obs_trd

! ---
! Observations of velocities by drIFters
   CALL obs_vdr

! ---
! Observations of velocities by gliders
   CALL obs_gvl

! ---
! Observations of SST
   CALL obs_sst

END SUBROUTINE obsop
