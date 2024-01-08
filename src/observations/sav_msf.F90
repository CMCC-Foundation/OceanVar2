subroutine sav_msf


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

 implicit none


! SLA observations
  if( sla%no.gt.0 )then
   sla%rss(:) = sla%res(:)
   sla%fls(:) = sla%flc(:)
   sla%flc(:) = 0
   sla%ns = sla%nc
   sla%nc = 0
  endif

! ARGO observations
  if( arg%no.gt.0 )then
   arg%rss(:) = arg%res(:)
   arg%fls(:) = arg%flc(:)
   arg%flc(:) = 0
   arg%ns = arg%nc
   arg%nc = 0
  endif

! XBT observations
  if( xbt%no.gt.0 )then
   xbt%rss(:) = xbt%res(:)
   xbt%fls(:) = xbt%flc(:)
   xbt%flc(:) = 0
   xbt%ns = xbt%nc
   xbt%nc = 0
  endif

! Glider observations
  if( gld%no.gt.0 )then
   gld%rss(:) = gld%res(:)
   gld%fls(:) = gld%flc(:)
   gld%flc(:) = 0
   gld%ns = gld%nc
   gld%nc = 0
  endif

! Argo trajectory observations
  if( tra%no.gt.0 )then
   tra%rsx(:) = tra%rex(:)
   tra%rsy(:) = tra%rey(:)
   tra%fls(:) = tra%flc(:)
   tra%flc(:) = 0
   tra%ns = tra%nc
   tra%nc = 0
   tra%ncs = tra%ncc
   tra%ncc = 0
  endif

! Trajectory observations of surface drifters
  if( trd%no.gt.0 )then
   trd%rsx(:) = trd%rex(:)
   trd%rsy(:) = trd%rey(:)
   trd%fls(:) = trd%flc(:)
   trd%flc(:) = 0
   trd%ns = trd%nc
   trd%nc = 0
   trd%ncs = trd%ncc
   trd%ncc = 0
  endif

! Observations of drifter velocities
  if( vdr%no.gt.0 )then
   vdr%rss(:) = vdr%res(:)
   vdr%fls(:) = vdr%flc(:)
   vdr%flc(:) = 0
   vdr%ns = vdr%nc
   vdr%nc = 0
  endif

! Observations of glider velocities
  if( gvl%no.gt.0 )then
   gvl%rss(:) = gvl%res(:)
   gvl%fls(:) = gvl%flc(:)
   gvl%flc(:) = 0
   gvl%ns = gvl%nc
   gvl%nc = 0
  endif


end subroutine sav_msf
