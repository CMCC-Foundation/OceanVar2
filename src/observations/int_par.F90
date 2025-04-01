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
!> Calculate interpolation parameters                                   
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE int_par

   USE set_knd
   USE obs_str

   IMPLICIT NONE

! ----
! Load SLA observations
   CALL int_par_sla

! ----
! Load ARGO observations
   CALL int_par_arg

! ----
! Load XBT observations
   CALL int_par_xbt

! ----
! Load glider observations
   CALL int_par_gld

! ----
! Load observations of Argo trajectories
   CALL int_par_tra

! ----
! Load observations of trajectories of surface drifters
   CALL int_par_trd

! ----
! Load observations of drifter velocities
   CALL int_par_vdr

! ----
! Load observations of glider velocities
   CALL int_par_gvl

! ----
! Load observations of SST
   CALL int_par_sst

END SUBROUTINE int_par
