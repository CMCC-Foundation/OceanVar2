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
!> Create the observational vector                                      !
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE mod_msf

   USE set_knd
   USE drv_str
   USE obs_str
   USE grd_str

   IMPLICIT NONE

! SLA observations
   IF ( sla%no.GT.0 ) THEN
      sla%res(:) = sla%res(:) - sla%inc(:)
      sla%ins(:) = sla%inc(:)
      sla%flc(:) = sla%fls(:)
      sla%fls(:) = 0
      sla%nc = sla%ns
   ENDIF

! ARGO observations
   IF ( arg%no.GT.0 ) THEN
      arg%res(:) = arg%res(:) - arg%inc(:)
      arg%ins(:) = arg%inc(:)
      arg%flc(:) = arg%fls(:)
      arg%fls(:) = 0
      arg%nc = arg%ns
   ENDIF

! XBT observations
   IF ( xbt%no.GT.0 ) THEN
      xbt%res(:) = xbt%res(:) - xbt%inc(:)
      xbt%ins(:) = xbt%inc(:)
      xbt%flc(:) = xbt%fls(:)
      xbt%fls(:) = 0
      xbt%nc = xbt%ns
   ENDIF

! Glider observations
   IF ( gld%no.GT.0 ) THEN
      gld%res(:) = gld%res(:) - gld%inc(:)
      gld%ins(:) = gld%inc(:)
      gld%flc(:) = gld%fls(:)
      gld%fls(:) = 0
      gld%nc = gld%ns
   ENDIF

! Argo trajectory observations
   IF ( tra%no.GT.0 ) THEN
      tra%rex(:) = tra%rex(:) - tra%inx(:)
      tra%isx(:) = tra%inx(:)
      tra%rey(:) = tra%rey(:) - tra%iny(:)
      tra%isy(:) = tra%iny(:)
      tra%flc(:) = tra%fls(:)
      tra%fls(:) = 0
      tra%nc = tra%ns
      tra%ncc = tra%ncs
      tra%ncs = 0
   ENDIF

! Trajectory observations of surface drIFters
   IF ( trd%no.GT.0 ) THEN
      trd%rex(:) = trd%rex(:) - trd%inx(:)
      trd%isx(:) = trd%inx(:)
      trd%rey(:) = trd%rey(:) - trd%iny(:)
      trd%isy(:) = trd%iny(:)
      trd%flc(:) = trd%fls(:)
      trd%fls(:) = 0
      trd%nc = trd%ns
      trd%ncc = trd%ncs
      trd%ncs = 0
   ENDIF

! Observations of drIFter velocities
   IF ( vdr%no.GT.0 ) THEN
      vdr%res(:) = vdr%res(:) - vdr%inc(:)
      vdr%ins(:) = vdr%inc(:)
      vdr%flc(:) = vdr%fls(:)
      vdr%fls(:) = 0
      vdr%nc = vdr%ns
   ENDIF

! Observations of glider velocities
   IF ( gvl%no.GT.0 ) THEN
      gvl%res(:) = gvl%res(:) - gvl%inc(:)
      gvl%ins(:) = gvl%inc(:)
      gvl%flc(:) = gvl%fls(:)
      gvl%fls(:) = 0
      gvl%nc = gvl%ns
   ENDIF

   ALLOCATE ( grd%tes(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%sas(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%uvs(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%vvs(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%ets(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )

   grd%tes(:,:,:) = grd%tem(:,:,:)
   grd%sas(:,:,:) = grd%sal(:,:,:)
   grd%uvs(:,:,:) = grd%uvl(:,:,:)
   grd%vvs(:,:,:) = grd%vvl(:,:,:)
   grd%ets(:,:)   = grd%eta(:,:)

! No SST assimilation
   IF ( sst%no .GT. 0 ) THEN
      sst%fls(:) = sst%flc(:)
      sst%flc(:) = 0
      sst%ns = sst%nc
      sst%nc = 0
   ENDIF

END SUBROUTINE mod_msf
