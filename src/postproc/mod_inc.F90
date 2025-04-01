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
!> Create the observational vector                                      
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE mod_inc

   USE set_knd
   USE drv_str
   USE obs_str
   USE grd_str

   IMPLICIT NONE

! SLA observations
   IF ( sla%nc .GT. 0 ) THEN
      sla%res(:) = sla%rss(:)
      sla%inc(:) = sla%inc(:) + sla%ins(:)
   ENDIF

! ARGO observations
   IF ( arg%nc .GT. 0 ) THEN
      arg%res(:) = arg%rss(:)
      arg%inc(:) = arg%inc(:) + arg%ins(:)
   ENDIF

! XBT observations
   IF ( xbt%nc .GT. 0 ) THEN
      xbt%res(:) = xbt%rss(:)
      xbt%inc(:) = xbt%inc(:) + xbt%ins(:)
   ENDIF

! Glider observations
   IF ( gld%nc .GT. 0 ) THEN
      gld%res(:) = gld%rss(:)
      gld%inc(:) = gld%inc(:) + gld%ins(:)
   ENDIF

! Argo trajectory observations
   IF ( tra%nc .GT. 0 ) THEN
      tra%rex(:) = tra%rsx(:)
      tra%inx(:) = tra%inx(:) + tra%isx(:)
      tra%rey(:) = tra%rsy(:)
      tra%iny(:) = tra%iny(:) + tra%isy(:)
   ENDIF

! Trajectory observations of surface drIFters
   IF ( trd%nc .GT. 0 ) THEN
      trd%rex(:) = trd%rsx(:)
      trd%inx(:) = trd%inx(:) + trd%isx(:)
      trd%rey(:) = trd%rsy(:)
      trd%iny(:) = trd%iny(:) + trd%isy(:)
   ENDIF

! Observations of drIFter velocities
   IF ( vdr%nc .GT. 0 ) THEN
      vdr%res(:) = vdr%rss(:)
      vdr%inc(:) = vdr%inc(:) + vdr%ins(:)
   ENDIF

! Observations of glider velocities
   IF ( gvl%nc .GT. 0 ) THEN
      gvl%res(:) = gvl%rss(:)
      gvl%inc(:) = gvl%inc(:) + gvl%ins(:)
   ENDIF

! Increments
   grd%tem(:,:,:) = grd%tem(:,:,:) + grd%tes(:,:,:)
   grd%sal(:,:,:) = grd%sal(:,:,:) + grd%sas(:,:,:)
   grd%uvl(:,:,:) = grd%uvl(:,:,:) + grd%uvs(:,:,:)
   grd%vvl(:,:,:) = grd%vvl(:,:,:) + grd%vvs(:,:,:)
   grd%eta(:,:)   = grd%eta(:,:)   + grd%ets(:,:)

! Observations of SST
   sst%flc(:) = sst%fls(:)
   sst%nc     = sst%ns

END SUBROUTINE mod_inc
