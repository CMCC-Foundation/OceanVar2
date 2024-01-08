subroutine mod_msf


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
! Create the observational vector                                      !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use drv_str
 use obs_str
 use grd_str

 implicit none

! SLA observations
  if( sla%no.gt.0 )then
   sla%res(:) = sla%res(:) - sla%inc(:)
   sla%ins(:) = sla%inc(:)
   sla%flc(:) = sla%fls(:)
   sla%fls(:) = 0
   sla%nc = sla%ns
  endif

! ARGO observations
  if( arg%no.gt.0 )then
   arg%res(:) = arg%res(:) - arg%inc(:)
   arg%ins(:) = arg%inc(:)
   arg%flc(:) = arg%fls(:)
   arg%fls(:) = 0
   arg%nc = arg%ns
  endif

! XBT observations
  if( xbt%no.gt.0 )then
   xbt%res(:) = xbt%res(:) - xbt%inc(:)
   xbt%ins(:) = xbt%inc(:)
   xbt%flc(:) = xbt%fls(:)
   xbt%fls(:) = 0
   xbt%nc = xbt%ns
  endif

! Glider observations
  if( gld%no.gt.0 )then
   gld%res(:) = gld%res(:) - gld%inc(:)
   gld%ins(:) = gld%inc(:)
   gld%flc(:) = gld%fls(:)
   gld%fls(:) = 0
   gld%nc = gld%ns
  endif

! Argo trajectory observations
  if( tra%no.gt.0 )then
   tra%rex(:) = tra%rex(:) - tra%inx(:)
   tra%isx(:) = tra%inx(:)
   tra%rey(:) = tra%rey(:) - tra%iny(:)
   tra%isy(:) = tra%iny(:)
   tra%flc(:) = tra%fls(:)
   tra%fls(:) = 0
   tra%nc = tra%ns
   tra%ncc = tra%ncs
   tra%ncs = 0
  endif

! Trajectory observations of surface drifters
  if( trd%no.gt.0 )then
   trd%rex(:) = trd%rex(:) - trd%inx(:)
   trd%isx(:) = trd%inx(:)
   trd%rey(:) = trd%rey(:) - trd%iny(:)
   trd%isy(:) = trd%iny(:)
   trd%flc(:) = trd%fls(:)
   trd%fls(:) = 0
   trd%nc = trd%ns
   trd%ncc = trd%ncs
   trd%ncs = 0
  endif

! Observations of drifter velocities
  if( vdr%no.gt.0 )then
   vdr%res(:) = vdr%res(:) - vdr%inc(:)
   vdr%ins(:) = vdr%inc(:)
   vdr%flc(:) = vdr%fls(:)
   vdr%fls(:) = 0
   vdr%nc = vdr%ns
  endif

! Observations of glider velocities
  if( gvl%no.gt.0 )then
   gvl%res(:) = gvl%res(:) - gvl%inc(:)
   gvl%ins(:) = gvl%inc(:)
   gvl%flc(:) = gvl%fls(:)
   gvl%fls(:) = 0
   gvl%nc = gvl%ns
  endif

     ALLOCATE ( grd%tes(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%sas(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%uvs(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%vvs(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%ets(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae))

     grd%tes(:,:,:) = grd%tem(:,:,:)
     grd%sas(:,:,:) = grd%sal(:,:,:)
     grd%uvs(:,:,:) = grd%uvl(:,:,:)
     grd%vvs(:,:,:) = grd%vvl(:,:,:)
     grd%ets(:,:)   = grd%eta(:,:)  


! No SST assimilation
  if( sst%no.gt.0 )then
     sst%fls(:) = sst%flg(:)
     sst%flg(:) = 0
     sst%ns = sst%nc
     sst%nc = 0
  endif

end subroutine mod_msf
