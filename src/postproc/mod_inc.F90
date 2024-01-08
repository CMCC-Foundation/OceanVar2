subroutine mod_inc


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
  if( sla%nc.gt.0 )then
   sla%res(:) = sla%rss(:) 
   sla%inc(:) = sla%inc(:) + sla%ins(:)
  endif

! ARGO observations
  if( arg%nc.gt.0 )then
   arg%res(:) = arg%rss(:)
   arg%inc(:) = arg%inc(:) + arg%ins(:)
  endif

! XBT observations
  if( xbt%nc.gt.0 )then
   xbt%res(:) = xbt%rss(:)
   xbt%inc(:) = xbt%inc(:) + xbt%ins(:)
  endif

! Glider observations
  if( gld%nc.gt.0 )then
   gld%res(:) = gld%rss(:)
   gld%inc(:) = gld%inc(:) + gld%ins(:)
  endif

! Argo trajectory observations
  if( tra%nc.gt.0 )then
   tra%rex(:) = tra%rsx(:)
   tra%inx(:) = tra%inx(:) + tra%isx(:)
   tra%rey(:) = tra%rsy(:)
   tra%iny(:) = tra%iny(:) + tra%isy(:)
  endif

! Trajectory observations of surface drifters
  if( trd%nc.gt.0 )then
   trd%rex(:) = trd%rsx(:)
   trd%inx(:) = trd%inx(:) + trd%isx(:)
   trd%rey(:) = trd%rsy(:)
   trd%iny(:) = trd%iny(:) + trd%isy(:)
  endif

! Observations of drifter velocities
  if( vdr%nc.gt.0 )then
   vdr%res(:) = vdr%rss(:)
   vdr%inc(:) = vdr%inc(:) + vdr%ins(:)
  endif

! Observations of glider velocities
  if( gvl%nc.gt.0 )then
   gvl%res(:) = gvl%rss(:)
   gvl%inc(:) = gvl%inc(:) + gvl%ins(:)
  endif


! Increments
     grd%tem(:,:,:) = grd%tem(:,:,:) + grd%tes(:,:,:)
     grd%sal(:,:,:) = grd%sal(:,:,:) + grd%sas(:,:,:)
     grd%uvl(:,:,:) = grd%uvl(:,:,:) + grd%uvs(:,:,:)
     grd%vvl(:,:,:) = grd%vvl(:,:,:) + grd%vvs(:,:,:)
     grd%eta(:,:)   = grd%eta(:,:)   + grd%ets(:,:)  

! Observations of SST
   sst%flg(:) = sst%fls(:)
   sst%nc = sst%ns

end subroutine mod_inc
