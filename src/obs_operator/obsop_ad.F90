subroutine obsop_ad

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
! Apply observational operators - adjoint
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use obs_str
 use grd_str
 use mpi_str

 implicit none


! ---

  obs%k = 0

! ---
! Satellite observations of SLA
  call obs_sla_ad

! ---
! ARGO observations 
  call obs_arg_ad

! ---
! XBT observations 
  call obs_xbt_ad

! ---
! Glider observations 
  call obs_gld_ad

! ---
! Observations of Argo trajectories
  if(tra%no.gt.0) call obs_tra_ad

! ---
! Observations of trajectories of surface drifters
  if(trd%no.gt.0) call obs_trd_ad

! ---
! Observations of velocity from drifters
  call obs_vdr_ad

! ---
! Observations of velocity from gliders
  call obs_gvl_ad

! ---
! Observations of SST
  call obs_sst_ad

! ---

! ---
 if(mpi%nproc.gt.1) then
  call exo_mpi( 1_i4, 0_i4,                                         &
                1_i4, grd%im+1, 1_i4,                               &
                1_i4, grd%jm+1, 1_i4,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae ,   1_i4, grd%eta_ad)
  call exo_mpi( 1_i4, 0_i4,                                         &
                1_i4, grd%im+1, 1_i4,                               &
                1_i4, grd%jm+1, 1_i4,                              &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%tem_ad)
  call exo_mpi( 1_i4, 0_i4,                                         &
                1_i4, grd%im+1, 1_i4,                               &
                1_i4, grd%jm+1, 1_i4,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%sal_ad)
  call exo_mpi( 1_i4, 0_i4,                                         &
                1_i4, grd%im+1, 1_i4,                               &
                1_i4, grd%jm+1, 1_i4,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%uvl_ad)
  call exo_mpi( 1_i4, 0_i4,                                         &
                1_i4, grd%im+1, 1_i4,                               &
                1_i4, grd%jm+1, 1_i4,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%vvl_ad)
 endif


end subroutine obsop_ad
