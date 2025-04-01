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
!> Load observations                                                   
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE get_obs

   USE set_knd
   USE obs_str

   IMPLICIT NONE

! ----
! Load SLA observations
   CALL get_obs_sla

! ----
! Load ARGO observations
   CALL get_obs_arg

! ----
! Load XBT observations
   CALL get_obs_xbt

! ----
! Load glider observations
   CALL get_obs_gld

! ----
! Load Argo trajectory observations
   CALL get_obs_tra

! ----
! Load trajectory observations of surface drifters
   CALL get_obs_trd

! ----
! Load observations of velocity by drifters
   CALL get_obs_vdr

! ----
! Load observations of velocity by drifters
   CALL get_obs_gvl

! ----
! Load observations of SST
   CALL get_obs_sst

END SUBROUTINE get_obs
