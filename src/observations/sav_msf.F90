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
!> Save the observational vector                                     
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE sav_msf

   USE set_knd
   USE drv_str
   USE obs_str

   IMPLICIT NONE

! SLA observations
   IF ( sla%no .GT. 0 ) THEN
      sla%rss = sla%res
      sla%fls = sla%flc
      sla%flc = 0
      sla%ns = sla%nc
      sla%nc = 0
   ENDIF

! ARGO observations
   IF ( arg%no .GT. 0 ) THEN
      arg%rss = arg%res
      arg%fls = arg%flc
      arg%flc = 0
      arg%ns = arg%nc
      arg%nc = 0
   ENDIF

! XBT observations
   IF ( xbt%no .GT. 0 ) THEN
      xbt%rss = xbt%res
      xbt%fls = xbt%flc
      xbt%flc = 0
      xbt%ns = xbt%nc
      xbt%nc = 0
   ENDIF

! Glider observations
   IF ( gld%no .GT. 0 ) THEN
      gld%rss = gld%res
      gld%fls = gld%flc
      gld%flc = 0
      gld%ns = gld%nc
      gld%nc = 0
   ENDIF

! Argo trajectory observations
   IF ( tra%no .GT. 0 ) THEN
      tra%rsx = tra%rex
      tra%rsy = tra%rey
      tra%fls = tra%flc
      tra%flc = 0
      tra%ns = tra%nc
      tra%nc = 0
      tra%ncs = tra%ncc
      tra%ncc = 0
   ENDIF

! Trajectory observations of surface drIFters
   IF ( trd%no .GT. 0 ) THEN
      trd%rsx = trd%rex
      trd%rsy = trd%rey
      trd%fls = trd%flc
      trd%flc = 0
      trd%ns = trd%nc
      trd%nc = 0
      trd%ncs = trd%ncc
      trd%ncc = 0
   ENDIF

! Observations of drIFter velocities
   IF ( vdr%no .GT. 0 ) THEN
      vdr%rss = vdr%res
      vdr%fls = vdr%flc
      vdr%flc = 0
      vdr%ns = vdr%nc
      vdr%nc = 0
   ENDIF

! Observations of glider velocities
   IF ( gvl%no .GT. 0 ) THEN
      gvl%rss = gvl%res
      gvl%fls = gvl%flc
      gvl%flc = 0
      gvl%ns = gvl%nc
      gvl%nc = 0
   ENDIF

END SUBROUTINE sav_msf
